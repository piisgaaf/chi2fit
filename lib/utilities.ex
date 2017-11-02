defmodule Chi2fit.Utilities do

  # Copyright 2015-2017 Pieter Rijken
  #
  # Licensed under the Apache License, Version 2.0 (the "License");
  # you may not use this file except in compliance with the License.
  # You may obtain a copy of the License at
  #
  #     http://www.apache.org/licenses/LICENSE-2.0
  #
  # Unless required by applicable law or agreed to in writing, software
  # distributed under the License is distributed on an "AS IS" BASIS,
  # WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  # See the License for the specific language governing permissions and
  # limitations under the License.

  @moduledoc """
  Provides various utilities:
  
    * Bootstrapping
    * Derivatives
    * Creating Cumulative Distribution Functions / Histograms from sample data
    * Solving linear, quadratic, and cubic equations
    * Autocorrelation coefficients
  
  """

  import Chi2fit.FFT
  
  @typedoc "Cumulative Distribution Function"
  @type cdf :: ((number)->{number,number,number})

  @typedoc "Algorithm used to assign errors to frequencey data: Wald score and Wilson score."
  @type algorithm :: :wilson | :wald

  @doc """
  Converts a list of number to frequency data.
  
  The data is divived into bins of size `binsize` and the number of data points inside a bin are counted. A map
  is returned with the bin's index as a key and value the number of data points in that bin.
  
  ## Examples
  
      iex> make_histogram [1,2,3]
      [{1, 1}, {2, 1}, {3, 1}]
      
      iex> make_histogram [1,2,3,4,5,6,5,4,3,4,5,6,7,8,9]
      [{1, 1}, {2, 1}, {3, 2}, {4, 3}, {5, 3}, {6, 2}, {7, 1}, {8, 1}, {9, 1}]
      
      iex> make_histogram [1,2,3,4,5,6,5,4,3,4,5,6,7,8,9], 3
      [{0, 2}, {1, 8}, {2, 4}, {3, 1}]

  """
  @spec make_histogram([number],number) :: %{required(number) => pos_integer}
  def make_histogram(list,binsize \\ 1) do
    Enum.reduce(list, %{}, fn
      (number,acc) ->
        acc |> Map.update(trunc(number/binsize),1,&(1+&1))
    end) |> Enum.reduce([], fn (pair,acc)->[pair|acc] end) |> Enum.sort_by(fn ({k,_v})->k end)
  end

  @doc """
  Returns a cumulative distribution corresponding to the input data.
  
  ## Example
  
      iex> to_cdf [1,2,3,4,5], 0.5, 1
      [{0.5, 0.0}, {1.5, 1.0}, {2.5, 2.0}, {3.5, 3.0}, {4.5, 4.0}, {5.5, 5.0}]
      
      iex> to_cdf [1,2,3,4,5,6,5,4,3,4,5,6,7,8,9], 0.5, 2
      [{0.5, 0.0}, {2.5, 2.0}, {4.5, 7.0}, {6.5, 12.0}, {8.5, 14.0}, {10.5, 15.0}]
  """
  @spec to_cdf([number],number,number) :: [ {float,float} ]
  def to_cdf(list, bin, interval \\ 0), do: to_cdf(list, bin, interval, 0.0, [])

  defp to_cdf([], _bin, _interval, _sum, result), do: Enum.reverse(result)
  defp to_cdf(list, bin, interval, sum, result) do
    {in_bin, out_bin} = list |> Enum.partition(fn (x)->x<=bin end)
    to_cdf(out_bin, bin+interval, interval, sum+length(in_bin), [{bin,sum+length(in_bin)}|result])
  end

  defmodule UnknownSampleErrorAlgorithmError do
    defexception message: "unknown sample error algorithm"
  end

  @doc """
  Converts a list of x,y data into a Cumulative Distribution function.
  
  Supports two ways of assigning errors: Wald score or Wilson score. See [1]. Valie values for the `algorithm`
  argument are `:wald` or `:wilson`.
  
  The second argument `numpoints` specifies the size of the original sample.
  
  The returned function returns tuples for its argument where the first element is the actual value of the
  function and the second and third elements gice the minimum and maximum confidence bounds.

  ## References
  
      [1] See https://en.wikipedia.org/wiki/Cumulative_frequency_analysis
      [2] https://arxiv.org/pdf/1112.2593v3.pdf
      [3] See https://en.wikipedia.org/wiki/Student%27s_t-distribution:
          90% confidence ==> t = 1.645 for many data points (> 120)
          70% confidence ==> t = 1.000
          
  ## Example
  
      iex(1)> fun = [1,2,3,4,5]
      ...> |> to_cdf(0.5, 1)
      ...> |> Enum.map(fn {x,y}->{x,y/5} end)
      ...> |> to_cdf_fun(5,:wilson)
      iex(2)> fun.(2.2)
      {0.2, 0.027223328910987405, 0.5233625708498564}
  
  """
  @spec to_cdf_fun([{x::number,y::number}],pos_integer,algorithm) :: cdf
  def to_cdf_fun(data,numpoints,algorithm \\ :wilson) do
    fn (x) ->
      y = data |> Enum.reverse |> Enum.find({nil,0.0}, fn ({xx,_})-> xx<=x end) |> elem(1)
      # t = 1.96
      t = 1.00

      case algorithm do
        :wald ->
          sd = :math.sqrt(y*(1.0-y)/numpoints)
          ylow = y - 2*y*t*sd
          yhigh = y + 2*(1.0-y)*t*sd
          {y,ylow,yhigh}

        :wilson ->
          ylow = if y > 0 do
            splus = t*t - 1/numpoints + 4*numpoints*y*(1-y) + (4*y - 2)
            srtplus = 1.0 + t*:math.sqrt(splus)
            max(0.0, (2*numpoints*y + t*t - srtplus)/2/(numpoints + t*t))
          else
            0.0
          end

          yhigh = if y < 1 do
            smin =  t*t - 1/numpoints + 4*numpoints*y*(1-y) - (4*y - 2)
            srtmin =  1.0 + t*:math.sqrt(smin)
            min(1.0, (2*numpoints*y + t*t + srtmin )/2/(numpoints + t*t))
          else
            1.0
          end

          {y,ylow,yhigh}

        other ->
          raise UnknownSampleErrorAlgorithmError, message: "unknown algorithm '#{inspect other}'"
      end
    end
  end

  @doc """
  Generates an empirical Cumulative Distribution Function from sample data.
  
  Three parameters determine the resulting empirical distribution:

    1) algorithm for assigning errors,
    
    2) the size of the bins,
    
    3) a correction for limiting the bounds on the 'y' values
  
  When e.g. task effort/duration is modeled, some tasks measured have 0 time. In practice
  what is actually is meant, is that the task effort is between 0 and 1 hour. This is where
  binning of the data happens. Specify a size of the bins to control how this is done. A bin
  size of 1 means that 0 effort will be mapped to 1/2 effort (at the middle of the bin).
  This also prevents problems when the fited distribution cannot cope with an effort os zero.
  
  In the handbook of MCMC [1] a cumulative distribution is constructed. For the largest 'x' value
  in the sample, the 'y' value is exactly one (1). In combination with the Wald score this
  gives zero errors on the value '1'. If the resulting distribution is used to fit a curve
  this may give an infinite contribution to the maximum likelihood function.
  Use the correction number to have a 'y' value of slightly less than 1 to prevent this from
  happening.
  Especially the combination of 0 correction, algorithm `:wald`, and 'linear' model for
  handling asymmetric errors gives problems.
  
  The algorithm parameter determines how the errors onthe 'y' value are determined. Currently
  supported values include `:wald` and `:wilson`.
  
  ## References
  
  [1] "Handbook of Monte Carlo Methods" by Kroese, Taimre, and Botev, section 8.4
  """
  @correction 0.01
  @spec empirical_cdf([{float,number}],integer,algorithm) :: {cdf,bins :: [float], numbins :: pos_integer, sum :: float}
  def empirical_cdf(data,binsize \\ 1,algorithm \\ :wilson) do
    {bins,sum} = data
      |> Enum.sort(fn ({x1,_},{x2,_})->x1<x2 end)
      |> Enum.reduce({[],0}, fn ({x,y},{acc,sum}) -> {[{binsize*(x+1/2),y+sum}|acc],sum+y} end)

    normbins = bins
      |> Enum.reverse
      |> Enum.map(fn ({x,y})->{x,y/(sum+trunc(Float.ceil(sum*@correction)))} end)

    {normbins |> to_cdf_fun(length(bins),algorithm),
     normbins,
     length(bins),
     sum}
  end

  @doc """
  Calculates the empirical CDF from a sample.
  
  Convenience function that chains `make_histogram/2` and `empirical_cdf/3`.
  """
  @spec get_cdf([number], number, algorithm) :: {cdf,bins :: [float], numbins :: pos_integer, sum :: float}
  def get_cdf(data, binsize \\ 1,algorithm \\ :wilson) do
    data
    |> make_histogram(binsize)
    |> empirical_cdf(binsize,algorithm)
  end

  @doc """
  Converts a CDF function to a list of data points.
  
  ## Example
  
      iex> convert_cdf {fn x->{:math.exp(-x),:math.exp(-x)/16,:math.exp(-x)/4} end, [1,4]}
      [{1, 0.6321205588285577, 0.9080301397071394, 0.9770075349267848},
                  {2, 0.8646647167633873, 0.9661661791908468, 0.9915415447977117},
                  {3, 0.950212931632136, 0.987553232908034, 0.9968883082270085},
                  {4, 0.9816843611112658, 0.9954210902778164, 0.9988552725694542}]

  """
  @spec convert_cdf({cdf,range :: [float,...]}) :: [{float,float,float,float}]
  def convert_cdf({cdf,[mindur,maxdur]}) do
    round(mindur)..round(maxdur)
    |> Stream.map(fn (x)->
        {y,y1,y2} = cdf.(x)
        {x,y,y1,y2}
      end)
    |> Stream.map(fn ({x,y,y1,y2})->{x,1.0-y,1.0-y2,1.0-y1} end)
    |> Enum.to_list
  end

  defp expand_pars(list,h) do
    list |> Enum.map(
          fn
            ({{x,0,factor}}) -> {{x,0,factor}}
            ({{x,0}}) -> {{x,0,1.0}}
            ({{x,n,factor}}) when n>0 ->
              xplus = x*(1.0+h)
              xmin = x*(1.0-h)
              dx = xplus-xmin
              [{{xplus,n-1,factor*dx}},{xmin,n-1,factor*dx}] |> expand_pars(h) |> List.flatten 
            ({{x,n}}) when n>0 ->
              xplus = x*(1.0+h)
              xmin = x*(1.0-h)
              dx = xplus-xmin
              [{{xplus,n-1,dx}},{xmin,n-1,dx}] |> expand_pars(h) |> List.flatten 
            ({x,0,factor}) -> {x,0,factor}
            ({x,0}) -> {x,0,1.0}
            ({x,n,factor}) when n>0 ->
              xplus = x*(1.0+h)
              xmin = x*(1.0-h)
              dx = xplus-xmin
              [{xplus,n-1,factor*dx},{{xmin,n-1,factor*dx}}] |> expand_pars(h) |> List.flatten 
            ({x,n}) when n>0 ->
              xplus = x*(1.0+h)
              xmin = x*(1.0-h)
              dx = xplus-xmin
              [{xplus,n-1,dx},{{xmin,n-1,dx}}] |> expand_pars(h) |> List.flatten 
            (x) when is_number(x) -> {x,0,1.0}
          end)
  end

  defp reduce_pars(list) do
    list |> Enum.reduce([{[],1,1.0}],
      fn
        (list,acc) when is_list(list) ->
          Enum.flat_map(list,
            fn
              ({{x,0,dx1}}) -> Enum.map(acc, fn ({y,n,dx2})->{[x|y],-n,dx1*dx2} end)
              ({x,0,dx1}) -> Enum.map(acc, fn ({y,n,dx2})->{[x|y],n,dx1*dx2} end)
            end)
        ({x,0,dx1},acc) -> Enum.map(acc, fn ({y,n,dx2})->{[x|y],n,dx1*dx2} end)
      end)
      |> Enum.map(fn ({l,n,dx}) -> {Enum.reverse(l),n,dx} end)
  end

  @doc """
  Calculates the partial derivative of a function and returns the value.
  
  ## Examples

      The function value at a point:
      iex> der([3.0], fn [x]-> x*x end) |> Float.round(3)
      9.0

      The first derivative of a function at a point:
      iex> der([{3.0,1}], fn [x]-> x*x end) |> Float.round(3)
      6.0

      The second derivative of a function at a point:
      iex> der([{3.0,2}], fn [x]-> x*x end) |> Float.round(3)
      2.0

      Partial derivatives with respect to two variables:
      iex> der([{2.0,1},{3.0,1}], fn [x,y] -> 3*x*x*y end) |> Float.round(3)
      12.0

  """
  @default_h 1.5e-3
  @spec der([float|{float,integer}], (([float])->float), Keyword.t) :: float
  def der(parameters, fun, options \\ []) do
    h = options[:h] || @default_h

    parameters
    |> expand_pars(h)
    |> reduce_pars
    |> Enum.reduce(0.0, fn ({x,n,dx},sum) when is_list(x) -> sum+n*fun.(x)/dx end)
  end

  defp jacobian(x=[_|_], k, fun) when k>0 and k<=length(x) and is_function(fun,1) do
    x |> List.update_at(k-1, fn (val) -> {val,1} end) |> der(fun,h: 1.0e-6)
  end

  @doc """
  Calculates the jacobian of the function at the point `x`.
  
  ## Examples
  
      iex> jacobian([2.0,3.0], fn [x,y] -> x*y end) |> Enum.map(&Float.round(&1))
      [3.0, 2.0]

  """
  @spec jacobian(x :: [float], (([float])->float)) :: [float]
  def jacobian(x, fun) do
    jacfun = &(jacobian(x, &1, fun))
    Enum.reduce(length(x)..1, [], fn (k,acc) -> [jacfun.(k)|acc] end)
  end

  @doc """
  Converts the input so that the result is a Puiseaux diagram, that is a strict convex shape.
  
  ## Examples
  
      iex> puiseaux [1]
      [1]
      
      iex> puiseaux [5,3,3,2]
      [5, 3, 2.5, 2]

  """
  @h 1.0e-10
  @spec puiseaux([number],[number],boolean) :: [number]
  def puiseaux(list,result \\ [],flag \\ false)
  def puiseaux([x],result,false), do: Enum.reverse [x|result]
  def puiseaux([x,y],result,false), do: Enum.reverse [y,x|result]
  def puiseaux([x,y],result,true), do: Enum.reverse([y,x|result]) |> puiseaux
  def puiseaux([x,y,z|rest],result,flag) do
    if y>(x+z)/2+@h do
      [(x+z)/2,z|rest] |> puiseaux([x|result],true)
    else
      [y,z|rest] |> puiseaux([x|result],flag)
    end
  end

  @doc """
  Calculates the autocorrelation coefficient of a list of observations.
  
  The implementation uses the discrete Fast Fourier Transform to calculate the autocorrelation.
  For available options see `Chi2fit.FFT.fft/2`. Returns a list of the autocorrelation coefficients.
  
  ## Example
  
      iex> auto [1,2,3]
      [14.0, 7.999999999999999, 2.999999999999997]

  """
  @spec auto([number],Keyword.t) :: [number]
  def auto(list,opts \\ [nproc: 1])
  
  def auto([],_opts), do: []
  def auto([x],_opts), do: [x*x]
  def auto(list,opts) do
    n = length(list)
    List.duplicate(0,n)
    |> Enum.concat(list)
    |> fft(opts) |> normv |> ifft(opts)
    |> Stream.take(n)
    |> Stream.map(&(elem(&1,0)))
    |> Enum.to_list
  end

  @doc """
  Calculates and returns the error associated with a list of observables.
  
  Usually these are the result of a Markov Chain Monte Carlo simulation run.
  
  The only supported method is the so-called `Initial Sequence Method`. See section 1.10.2 (Initial sequence method)
  of [1].
  
  Input is a list of autocorrelation coefficients. This may be the output of `auto/2`.
  
  ## References
  
    [1] 'Handbook of Markov Chain Monte Carlo'
  """
  @spec error([{gamma :: number,k :: pos_integer}], :initial_sequence_method) :: {var :: number, lag :: number}
  def error(nauto, :initial_sequence_method) do
    ## For reversible Markov Chains
    gamma = nauto |> Stream.chunk(2) |> Stream.map(fn ([{x,k},{y,_}])->{k/2,x+y} end) |> Enum.to_list
    gamma0 = nauto |> Stream.take(1) |> Enum.to_list |> (&(elem(hd(&1),0))).()
    m = gamma |> Stream.take_while(fn ({_k,x})->x>0 end) |> Enum.count
    gammap = gamma |> Stream.take_while(fn ({_k,x})->x>0 end) |> Stream.map(fn {_,x}->x end) |> Stream.concat([0.0]) |> Enum.to_list
    gammap = gammap |> puiseaux
    var = -gamma0 + 2.0*(gammap |> Enum.sum)

    if var < 0, do: throw {:negative_variance, var, 2*m}
    {var,2*m}
  end

  @doc """
  Implements bootstrapping procedure as resampling with replacement.
  
  It supports saving intermediate results to a file using `:dets`. Use the options `:safe` and `:filename` (see below)
  
  ## Arguments:
  
      `total` - Total number resmaplings to perform
      `data` - The sample data
      `fun` - The function to evaluate
      `options` - A keyword list of options, see below.
      
  ## Options
  
      `:safe` - Whether to safe intermediate results to a file, so as to support continuation when it is interrupted.
            Valid values are `:safe` and `:cont`.
      `:filename` - The filename to use for storing intermediate results
  
  """
  @spec bootstrap(total :: integer, data :: [number], fun :: (([number],integer)->number), options :: Keyword.t) :: [number]
  def bootstrap(total, data, fun, options) do
    safe = options |> Keyword.get(:safe, false)

    {start,continuation} = case safe do
      :safe ->
        file = options |> Keyword.fetch!(:filename)
        {:ok,:storage} = :dets.open_file :storage, type: :set, file: file, auto_save: 1000, estimated_no_objects: total
        :ok = :dets.delete_all_objects :storage
        {1,[]}        
      :cont ->
        file = options |> Keyword.fetch!(:filename)
        {:ok,:storage} = :dets.open_file :storage, type: :set, file: file, auto_save: 1000, estimated_no_objects: total
        objects = :dets.select(:storage, [{{:'_',:'$1'},[],[:'$1']}])
        {length(objects)+1,objects}
      _ ->
        {1,[]}
    end

    if start>total, do: raise ArgumentError, message: "start cannot be larger than the total"

    1..total |> Enum.reduce(continuation, fn (k,acc) ->
      try do
        ## Evaluate the function
        result = data |> Enum.map(fn _ -> Enum.random(data) end) |> fun.(k)

        if safe, do: true = :dets.insert_new :storage, {k,result}

        [result|acc]
      rescue
        _error ->
          [nil|acc]
      end
    end)
  end
  
  @doc """
  Reads data from a file specified by `filename` and returns a stream with the data parsed as floats.

  It expects a single data point on a separate line and removes entries that:
  
    * are not floats, and
    * smaller than zero (0)

  """
  @spec read_data(filename::String.t) :: Stream.t
  def read_data(filename) do
    filename
    |> File.stream!([],:line)
    |> Stream.flat_map(&String.split(&1,"\r",trim: true))
    |> Stream.filter(&is_tuple(Float.parse(&1)))
    |> Stream.map(&elem(Float.parse(&1),0))
    |> Stream.filter(&(&1 >= 0.0))
  end

end
