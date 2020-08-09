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
  
  alias Chi2fit.Distribution, as: D
  alias Chi2fit.Fit, as: F
  alias Chi2fit.Matrix, as: M
  
  @typedoc "Cumulative Distribution Function"
  @type cdf :: ((number)->{number,number,number})

  @typedoc "Binned data with error bounds specified through low and high values"
  @type ecdf :: [{float,float,float,float}]

  @typedoc "Algorithm used to assign errors to frequencey data: Wald score and Wilson score."
  @type algorithm :: :wilson | :wald

  @typedoc "Supported numerical integration methods"
  @type method :: :gauss | :gauss2 | :gauss3 | :romberg | :romberg2 | :romberg3

  @typedoc "Average and standard deviationm (error)"
  @type avgsd :: {avg :: float, sd :: float}

  @doc """
  Converts a list of numbers to frequency data.
  
  The data is divided into bins of size `binsize` and the number of data points inside a bin are counted. A map
  is returned with the bin's index as a key and as value the number of data points in that bin.
  
  ## Examples
  
      iex> make_histogram [1,2,3]
      [{1, 1}, {2, 1}, {3, 1}]
      
      iex> make_histogram [1,2,3], 1.0, 0
      [{1, 1}, {2, 1}, {3, 1}]
      
      iex> make_histogram [1,2,3,4,5,6,5,4,3,4,5,6,7,8,9]
      [{1, 1}, {2, 1}, {3, 2}, {4, 3}, {5, 3}, {6, 2}, {7, 1}, {8, 1}, {9, 1}]
      
      iex> make_histogram [1,2,3,4,5,6,5,4,3,4,5,6,7,8,9], 3, 1.5
      [{0, 1}, {1, 6}, {2, 6}, {3, 2}]

  """
  @spec make_histogram([number],number,number) :: [{non_neg_integer,pos_integer}]
  def make_histogram(list,binsize \\ 1.0,offset \\ 0.5)
  def make_histogram(list,binsize,offset) when binsize>offset do
    Enum.reduce(list, %{}, fn
      (number,acc) ->
        acc |> Map.update(if(number<=offset,do: 0, else: trunc(Float.ceil((number-offset)/binsize))),1,&(1+&1))
    end) |> Enum.reduce([], fn (pair,acc)->[pair|acc] end) |> Enum.sort_by(fn ({k,_v})->k end)
  end
  def make_histogram(_list,_binsize,_offset), do: raise ArgumentError, message: "binsize must be larger than bin offset"

  defmodule UnknownSampleErrorAlgorithmError do
    defexception message: "unknown sample error algorithm"
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
  
  Supports two ways of assigning errors: Wald score or Wilson score. See [1]. Valie values for the `algorithm`
  argument are `:wald` or `:wilson`.

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
      [2] See https://en.wikipedia.org/wiki/Cumulative_frequency_analysis
      [3] https://arxiv.org/pdf/1112.2593v3.pdf
      [4] See https://en.wikipedia.org/wiki/Student%27s_t-distribution:
          90% confidence ==> t = 1.645 for many data points (> 120)
          70% confidence ==> t = 1.000
  """
  @spec empirical_cdf([{float,number}],{number,number},algorithm,integer) :: {cdf,bins :: [float], numbins :: pos_integer, sum :: float}
  def empirical_cdf(data,bin \\ {1.0,0.5},algorithm \\ :wilson,correction \\ 0)
  def empirical_cdf(data,{binsize,offset},algorithm,correction) do
    {bins,sum} = data
      |> Enum.sort(fn ({x1,_},{x2,_})->x1<x2 end)
      |> Enum.reduce({[],0}, fn ({n,y},{acc,sum}) -> {[{offset+binsize*n,y+sum}|acc],sum+y} end)

    normbins = bins
      |> Enum.reverse
      |> Enum.map(fn ({x,y})->{x,y/(sum+correction),y} end)

    {normbins |> to_cdf_fun(sum,algorithm),
     normbins,
     length(bins),
     sum}
  end

  @doc """
  Calculates the empirical CDF from a sample.
  
  Convenience function that chains `make_histogram/2` and `empirical_cdf/3`.
  """
  @spec get_cdf([number], number|{number,number}, algorithm, integer) :: {cdf,bins :: [float], numbins :: pos_integer, sum :: float}
  def get_cdf(data, binsize \\ {1.0,0.5},algorithm \\ :wilson, correction \\ 0)
  def get_cdf(data, {binsize,offset},algorithm, correction) do
    data
    |> make_histogram(binsize,offset)
    |> empirical_cdf({binsize,offset},algorithm,correction)
  end

  @doc """
  Converts a CDF function to a list of data points.
  
  ## Example
  
      iex> convert_cdf {fn x->{:math.exp(-x),:math.exp(-x)/16,:math.exp(-x)/4} end, {1,4}}
      [{1, 0.36787944117144233, 0.022992465073215146, 0.09196986029286058},
       {2, 0.1353352832366127, 0.008458455202288294, 0.033833820809153176},
       {3, 0.049787068367863944, 0.0031116917729914965, 0.012446767091965986},
       {4, 0.01831563888873418, 0.0011447274305458862, 0.004578909722183545}]

  """
  @type range :: {float,float} | [float,...]
  @spec convert_cdf({cdf,range}) :: [{float,float,float,float}]
  def convert_cdf({cdf,{mindur,maxdur}}), do: round(mindur)..round(maxdur) |> y_with_errors(cdf)
  def convert_cdf({cdf,categories}) when is_list(categories), do: categories |> y_with_errors(cdf)
  defp y_with_errors(list,cdf), do: list |> Enum.map(&Tuple.insert_at(cdf.(&1),0,&1)) 

  @doc """
  Converts raw data to binned data with (asymmetrical) errors.
  """
  @spec to_bins(data :: [number], binsize :: {number,number}) :: ecdf()
  def to_bins(data,binsize \\ {1.0,0.5}) do
    # Convert the raw data to binned data (histogram or frequency data):
    {cdf,bins,_,_} = get_cdf data, binsize

    # Add the errors based on the binomial distribution (Wilson score):
    convert_cdf {cdf,bins|>Enum.map(&elem(&1,0))}
  end
  
  @doc """
  Calculates the nth moment of the sample.
  
  ## Example
  
      iex> moment [1,2,3,4,5,6], 1
      3.5
  """
  @spec moment(sample::[number],n::pos_integer) :: float
  def moment(sample,n) when length(sample)>0 and is_integer(n) and n>0 do
    (sample |> Stream.map(fn x-> :math.pow(x,n) end) |> Enum.sum)/length(sample)
  end

  @doc """
  Calculates the nth centralized moment of the sample.
  
  ## Example
  
      iex> momentc [1,2,3,4,5,6], 1
      0.0
  
      iex> momentc [1,2,3,4,5,6], 2
      2.9166666666666665
  """
  @spec momentc(sample::[number],n::pos_integer) :: float
  def momentc(sample,n) when length(sample)>0 and is_integer(n) and n>0 do
    mean = sample |> moment(1)
    sample |> momentc(n,mean)
  end
  
  @doc """
  Calculates the nth centralized moment of the sample.
  
  ## Example
  
      iex> momentc [1,2,3,4,5,6], 2, 3.5
      2.9166666666666665
  """
  @spec momentc(sample::[number],n::pos_integer,mu::float) :: float
  def momentc(sample,n,mu) when length(sample)>0 and is_integer(n) and n>0 do
    (sample |> Stream.map(fn x-> :math.pow(x-mu,n) end) |> Enum.sum)/length(sample)
  end

  @doc """
  Calculates the nth normalized moment of the sample.
  
  ## Example
  
      iex> momentn [1,2,3,4,5,6], 1
      0.0
  
      iex> momentn [1,2,3,4,5,6], 2
      1.0
  
      iex> momentn [1,2,3,4,5,6], 4
      1.7314285714285718
  """
  @spec momentn(sample::[number],n::pos_integer) :: float
  def momentn(sample,n) when length(sample)>0 and is_integer(n) and n>0 do
    mean = sample |> moment(1)
    sample |> momentn(n,mean)
  end

  @doc """
  Calculates the nth normalized moment of the sample.
  
  ## Example
  
      iex> momentn [1,2,3,4,5,6], 4, 3.5
      1.7314285714285718
  """
  @spec momentn(sample::[number],n::pos_integer,mu::float) :: float
  def momentn(sample,n,mu) when length(sample)>0 and is_integer(n) and n>0 do
    sigma = :math.sqrt(sample |> momentc(2,mu))
    (sample |> momentc(n,mu))/:math.pow(sigma,n)
  end

  @doc """
  Calculates the nth normalized moment of the sample.
  """
  @spec momentn(sample::[number],n::pos_integer,mu::float,sigma::float) :: float
  def momentn(sample,n,mu,sigma) when length(sample)>0 and is_integer(n) and n>0 and sigma>0.0 do
    (sample |> momentc(n,mu))/:math.pow(sigma,n)
  end

  @type cullenfrey :: [{squared_skewness::float,kurtosis::float}|nil]
  
  @doc """
  Generates a Cullen & Frey plot for the sample data.
  
  The kurtosis returned is the 'excess kurtosis'.
  """
  @spec cullen_frey(sample::[number], n::integer) :: cullenfrey
  def cullen_frey(sample,n \\ 100) do
    bootstrap(n,sample,
      fn
        data,_i ->
          mean = data |> moment(1)
          sigma = :math.sqrt(data |> momentc(2))
          skewness = data |> momentn(3,mean,sigma)
          kurtosis = data |> momentn(4,mean,sigma)
          {skewness*skewness,kurtosis-3.0}
      end)
  end

  @doc """
  Extracts data point with standard deviation from Cullen & Frey plot data.
  """
  @spec cullen_frey_point(data::cullenfrey) :: {{x::float,dx::float},{y::float,dy::float}}
  def cullen_frey_point(data) do
    {skew,kurt} = data |> Stream.filter(fn x -> x end) |> Enum.unzip
    {
      {moment(skew,1),momentc(skew,2)},
      {moment(kurt,1),momentc(kurt,2)}
    }
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
  @default_h 0.001
  @spec der([float|{float,integer}], (([float])->float), Keyword.t) :: float
  def der(parameters, fun, options \\ []) do
    richardson(fn acc ->
        result = parameters
        |> expand_pars(acc)
        |> reduce_pars
        |> Enum.reduce(0.0, fn ({x,n,dx},sum) when is_list(x) -> sum+n*fun.(x)/dx end)
        {result,acc/2.0}
      end,
      @default_h,4.0,options)
  end

  @doc """
  Calculates the jacobian of the function at the point `x`.
  
  ## Examples
  
      iex> jacobian([2.0,3.0], fn [x,y] -> x*y end) |> Enum.map(&Float.round(&1))
      [3.0, 2.0]

  """
  @spec jacobian(x :: [float], (([float])->float)) :: [float]
  def jacobian(x, fun, options \\ []) do
    jacfun = &(jacobian(x, &1, fun, options))
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
    gamma = nauto |> Stream.chunk_every(2) |> Stream.map(fn ([{x,k},{y,_}])->{k/2,x+y} end) |> Enum.to_list
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
  
      `total` - Total number resamplings to perform
      `data` - The sample data
      `fun` - The function to evaluate
      `options` - A keyword list of options, see below.
      
  ## Options
  
      `:safe` - Whether to safe intermediate results to a file, so as to support continuation when it is interrupted.
            Valid values are `:safe` and `:cont`.
      `:filename` - The filename to use for storing intermediate results
  
  """
  @spec bootstrap(total :: integer, data :: [number], fun :: (([number],integer)->number), options :: Keyword.t) :: [any]
  def bootstrap(total, data, fun, options \\ []) do
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
        objects = :dets.select(:storage, [{{:_,:'$1'},[],[:'$1']}])
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
  Reamples the subsequences of numbers contained in the list as determined by `analyze/2`
  """
  @spec resample(data :: [number], options :: Keyword.t) :: [number]
  def resample(data,options) do
    data
    |> analyze(fn dat,opt ->
          F.find_all(dat,opt) |> Enum.flat_map(fn {_,_,d}->resample(d) end)
      end,
      options)
  end

  defp resample(data), do: Enum.map(data,fn _ -> Enum.random(data) end)
  
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

  ## TODO: implement gauss-kronrad integration (progressive gauss)
  @doc """
  Numerical integration providing Gauss and Romberg types.
  """
  @default_points 32
  @spec integrate(method, ((float)->float), a::float, b::float, options::Keyword.t) :: float
  def integrate(method, func, a, b, options \\ [])
  def integrate(:gauss, func, a, b, options) do
    npoints = options[:points] || @default_points

    factor_min = (b-a)/2.0
    factor_plus = (b+a)/2.0

    {weights,abscissa} = case npoints do
      4 ->
        {
          [ 0.6521451548625461,0.3478548451374538 ],
          [ 0.3399810435848563,0.8611363115940526 ]
         }
      8 ->
        {
          [ 0.3626837833783620,0.3137066458778873,0.2223810344533745,0.1012285362903763 ],
          [ 0.1834346424956498,0.5255324099163290,0.7966664774136267,0.9602898564975363 ]
        }
      32 ->
        {
          [ 0.0965400885147278,0.0956387200792749,0.0938443990808046,0.0911738786957639,0.0876520930044038,0.0833119242269467,0.0781938957870703,0.0723457941088485,0.0658222227763618,0.0586840934785355,0.0509980592623762,0.0428358980222267,0.0342738629130214,0.0253920653092621,0.0162743947309057,0.0070186100094701 ],
          [ 0.0483076656877383,0.1444719615827965,0.2392873622521371,0.3318686022821277,0.4213512761306353,0.5068999089322294,0.5877157572407623,0.6630442669302152,0.7321821187402897,0.7944837959679424,0.8493676137325700,0.8963211557660521,0.9349060759377397,0.9647622555875064,0.9856115115452684,0.9972638618494816 ]
        }
    end

    factor_min * (Enum.zip(abscissa,weights) |> Enum.map(fn {x,w} -> w*( func.(factor_min*x+factor_plus) + func.(-factor_min*x+factor_plus) ) end) |> Enum.sum)
  end
  def integrate(:gauss2, func, a, :infinity, options) do
    fac = 500.0 ## t = tanh(x/fac)
    fac*integrate(:gauss, fn t -> (func.(fac*:math.atanh(t)))/(1.0-t*t) end, :math.tanh(a/fac), 1.0, options)
  end
  def integrate(:gauss2, func, a, b, options) do
    fac = 500.0 ## t = tanh(x/fac)
    fac*integrate(:gauss, fn t -> (func.(fac*:math.atanh(t)))/(1.0-t*t) end, :math.tanh(a/fac), :math.tanh(b/fac), options)
  end
  def integrate(:gauss3, func, a, :infinity, options) do
    ## x = t/(1-t) = -1 + 1/(1-t), dx = dt/(1-t)^2
    integrate(:gauss, fn t -> (func.(t/(1.0-t)))/(1.0-t)/(1.0-t) end, a/(a+1.0), 1.0, options)
  end
  def integrate(:gauss3, func, a, b, options) do
    ## x = t/(1-t) = -1 + 1/(1-t), dx = dt/(1-t)^2
    integrate(:gauss, fn t -> (func.(t/(1.0-t)))/(1.0-t)/(1.0-t) end, a/(a+1.0), b/(b+1.0), options)
  end

  @default_tolerance 1.0e-6
  def integrate(:romberg, func, a, b, options) do
    richardson(fn acc ->
        case acc do
          [] ->
            f1 = try do func.(a) rescue _e -> 0.0 end
            f2 = try do func.(b) rescue _e -> 0.0 end
            result = (b-a) * ( f1 + f2 )/2.0
            {result,[{a,f1},{b,f2}]}

          values ->
            vals = values
            |> Stream.transform(nil, fn
              {x2,f},nil -> {[{x2,f}],x2}
              {x2,f},x1 -> {[{(x2+x1)/2.0,func.((x2+x1)/2.0)},{x2,f}],x2}
            end)
            |> Enum.to_list

            result = vals
            |> Stream.chunk_every(2,1,:discard)
            |> Stream.map(fn [{x1,f1},{x2,f2}] -> (x2-x1)*( f1 + f2 )/2.0 end)
            |> Enum.sum
            {result,vals}
        end
      end, [], 4.0, options)
  end
  def integrate(:romberg2, func, a, :infinity, options) do
    fac = 500.0 ## t = tanh(x/fac)
    integrate(:romberg, fn t -> (func.(fac*:math.atanh(t)))*fac/(1.0-t*t) end, :math.tanh(a/fac), 1.0, options)
  end
  def integrate(:romberg2, func, a, b, options) do
    fac = 500.0 ## t = tanh(x/fac)
    integrate(:romberg, fn t -> (func.(fac*:math.atanh(t)))*fac/(1.0-t*t) end, :math.tanh(a/fac), :math.tanh(b/fac), options)
  end
  def integrate(:romberg3, func, a, :infinity, options) do
    ## x = t/(1-t) = -1 + 1/(1-t), dx = dt/(1-t)^2
    integrate(:romberg, fn t -> (func.(t/(1.0-t)))/(1.0-t)/(1.0-t) end, a/(a+1.0), 1.0, options)
  end
  def integrate(:romberg3, func, a, b, options) do
    ## x = t/(1-t) = -1 + 1/(1-t), dx = dt/(1-t)^2
    integrate(:romberg, fn t -> (func.(t/(1.0-t)))/(1.0-t)/(1.0-t) end, a/(a+1.0), b/(b+1.0), options)
  end

  @doc """
  Richardson extrapolation.
  """
  @default_tolerance 1.0e-6
  @spec richardson(func::((term)->{float,term}), init::term, factor::float, results::[float], options::Keyword.t) :: float
  def richardson(func, init, factor, results \\ [], options)
  def richardson(func, init, factor, results, options) do
    tolerance = options[:tolerance] || @default_tolerance
    max = options[:itermax]

    {result,acc} = func.(init)
    {new,last,error,_} = results |> Enum.reduce({[],result,nil,factor}, fn
        _prev,{acc,item,0.0,order} ->
          {acc,item,0.0,order}
        prev,{acc,item,_,order} ->
          diff = (order*item - prev)/(order-1.0)
          {[diff|acc],diff,if(diff==0, do: 0.0, else: abs((diff-item)/diff)),order*factor}
        end)
    cond do
      max && (length(new) > max) ->
        last
      error < tolerance ->
        last
      true ->
        richardson(func, acc, factor, [result|Enum.reverse(new)], options)
    end
  end

  @doc """
  Newton-Fourier method for locating roots and returning the interval where the root is located.
  
  See [https://en.wikipedia.org/wiki/Newton%27s_method#Newton.E2.80.93Fourier_method]
  """
  @spec newton(a::float,b::float,func::((x::float)->float),maxiter::non_neg_integer,options::Keyword.t) :: {float, {float,float}, {float,float}}
  def newton(a,b,func,maxiter \\ 10, options), do: newton(a,b,func,maxiter,{(a+b)/2,{a,b},{nil,nil}},options)

  @doc """
  Unzips lists of 1-, 2-, 3-, 4-, 5-, 6-, 7-, and 8-tuples.
  """
  @spec unzip(list::[tuple]) :: tuple
  def unzip([]), do: {}
  def unzip(list=[{_}|_]), do: {Enum.map(list,fn {x}->x end)}
  def unzip(list=[{_,_}|_]), do: Enum.unzip(list)
  def unzip(list=[{_,_,_}|_]) do
    {
      list |> Enum.map(&elem(&1,0)),
      list |> Enum.map(&elem(&1,1)),
      list |> Enum.map(&elem(&1,2))
    }
  end
  def unzip(list=[{_,_,_,_}|_]) do
    {
      list |> Enum.map(&elem(&1,0)),
      list |> Enum.map(&elem(&1,1)),
      list |> Enum.map(&elem(&1,2)),
      list |> Enum.map(&elem(&1,3))
    }
  end
  def unzip(list=[{_,_,_,_,_}|_]) do
    {
      list |> Enum.map(&elem(&1,0)),
      list |> Enum.map(&elem(&1,1)),
      list |> Enum.map(&elem(&1,2)),
      list |> Enum.map(&elem(&1,3)),
      list |> Enum.map(&elem(&1,4))
    }
  end
  def unzip(list=[{_,_,_,_,_,_}|_]) do
    {
      list |> Enum.map(&elem(&1,0)),
      list |> Enum.map(&elem(&1,1)),
      list |> Enum.map(&elem(&1,2)),
      list |> Enum.map(&elem(&1,3)),
      list |> Enum.map(&elem(&1,4)),
      list |> Enum.map(&elem(&1,5))
    }
  end
  def unzip(list=[{_,_,_,_,_,_,_}|_]) do
    {
      list |> Enum.map(&elem(&1,0)),
      list |> Enum.map(&elem(&1,1)),
      list |> Enum.map(&elem(&1,2)),
      list |> Enum.map(&elem(&1,3)),
      list |> Enum.map(&elem(&1,4)),
      list |> Enum.map(&elem(&1,5)),
      list |> Enum.map(&elem(&1,6))
    }
  end
  def unzip(list=[_|_]) do
    0..tuple_size(hd(list))-1
    |> Enum.reduce({},fn i,tup -> Tuple.append(tup,list |> Enum.map(&elem(&1,i))) end)
  end

  ##
  ## Local functions
  ##
  
  @spec to_cdf_fun([{x::number,y::number,n::integer}],pos_integer,algorithm) :: cdf
  defp to_cdf_fun(data,numpoints,algorithm) do
    fn (x) ->
      y = data |> Enum.reverse |> Enum.find({nil,0.0}, fn ({xx,_,_})-> xx<=x end) |> elem(1)
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
            if splus < 0.0 do
              0.0
            else
              srtplus = 1.0 + t*:math.sqrt(splus)
              max(0.0, (2*numpoints*y + t*t - srtplus)/2/(numpoints + t*t))
            end
          else
            0.0
          end

          yhigh = if y < 1 do
            smin =  t*t - 1/numpoints + 4*numpoints*y*(1-y) - (4*y - 2)
            if smin < 0.0 do
              1.0
            else
              srtmin =  1.0 + t*:math.sqrt(smin)
              min(1.0, (2*numpoints*y + t*t + srtmin )/2/(numpoints + t*t))
            end
          else
            1.0
          end

          {y,ylow,yhigh}

        other ->
          raise UnknownSampleErrorAlgorithmError, message: "unknown algorithm '#{inspect other}'"
      end
    end
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

  defp jacobian(x=[_|_], k, fun, options) when k>0 and k<=length(x) and is_function(fun,1) do
    x |> List.update_at(k-1, fn (val) -> {val,1} end) |> der(fun,options)
  end

  @default_rel_tolerance 1.0e-6
  defp newton(_a,_b,func,0,{root,{l,r},_},_options), do: {root,{l,r},{func.(l),func.(r)}}
  defp newton(a,b,func,maxiter,{prev,{left,right},{vleft,vright}},options) do
    tolerance = options[:tolerance] || @default_rel_tolerance

    x0 = func.(right)
    z0 = func.(left)

    if x0*z0 > 0 do
      raise ArgumentError, message: "Interval does not contain root"
    end

    derx0 = der([{right,1}], fn [x]->func.(x) end, options)

    if derx0 == 0 do
      raise ArithmeticError,
        message: "Interval contains local minimum/maximum [left/z0=#{left}/#{z0}; right/x0=#{right}/#{x0}; der=#{derx0}]"
    end

    x1 = right - x0/derx0
    z1 = left - z0/derx0
    root = (x1+z1)/2.0

    cond do
      z1 < left ->
        newton(a,b,func,0,{prev,{left,right},{vleft,vright}},options)

      x1 > right ->
        newton(a,b,func,0,{prev,{left,right},{vleft,vright}},options)

      z1 < x1 and abs(x1-z1) < tolerance ->
        newton(a,b,func,0,{root,{z1,x1},{z0,x0}},options)

      z1 > x1 and abs(x1-z1) < tolerance ->
        newton(a,b,func,0,{root,{x1,z1},{z0,x0}},options)

      z1 > x1 ->
        newton(a,b,func,maxiter-1,{prev,{x1,z1},{z0,x0}},options)

      true ->
        newton(a,b,func,maxiter-1,{root,{z1,x1},{z0,x0}},options)
    end
  end

  @doc """
  Outputs and formats the errors that result from a call to `Chi2fit.Fit.chi2/4`
  
  Errors are tuples of length 2 and larger: `{[min1,max1], [min2,max2], ...}`.
  """
  @spec puts_errors(device :: IO.device(), errors :: tuple()) :: none()
  def puts_errors(device \\ :stdio, errors) do
    errors
    |> Tuple.to_list
    |> Enum.with_index
    |> Enum.each(fn
        {[mn,mx],0} -> IO.puts device, "\t\t\tchi2:\t\t#{mn}\t-\t#{mx}"
        {[mn,mx],_} -> IO.puts device, "\t\t\tparameter:\t#{mn}\t-\t#{mx}"
    end)
  end

  @doc """
  Forecasts how many time periods are needed to complete `size` items
  
  Related functions: `forecast_duration/2` and `forecast_items/2`.
  """
  @spec forecast(fun :: (() -> non_neg_integer),size :: pos_integer, tries :: pos_integer, update :: (() -> number)) :: number
  def forecast(fun, size, tries \\ 0,update \\ fn -> 1 end)
  def forecast(fun, size, tries, update) when size>0 do
      forecast(fun, size-fun.(),tries+update.(),update)
  end
  def forecast(_fun,_size,tries,_update), do: tries

  @doc """
  Returns a function for forecasting the duration to complete a number of items.
  
  This function is a wrapper for `forecast/4`.

  ## Arguments
  
      `data` - either a data set to base the forecasting on, or a function that returns (random) numbers
      `size` - the number of items to complete

  """
  @spec forecast_duration(data :: [number] | (()->number), size :: pos_integer) :: (() -> number)
  def forecast_duration(data, size) when is_list(data) do
    fn -> forecast(fn -> Enum.random(data) end, size) end
  end
  def forecast_duration(fun, size) when is_function(fun,0) do
    fn -> forecast(fun, size) end
  end

  @doc """
  Returns a function for forecasting the number of completed items in a number periods.
  
  This function is a wrapper for `forecast/4`.

  ## Arguments
  
      `data` - either a data set to base the forecasting on, or a function that returns (random) numbers
      `periods` - the number of periods to forecast the number of completed items for

  """
  @spec forecast_items(data :: [number] | (()->number), periods :: pos_integer) :: (() -> number)
  def forecast_items(data, periods) when is_list(data) do
    fn -> forecast(fn -> 1 end, periods, 0, fn -> Enum.random(data) end) end
  end
  def forecast_items(fun, periods) when is_function(fun,0) do
    fn -> forecast(fn -> 1 end, periods, 0, fun) end
  end

  @doc """
  Basic Monte Carlo simulation to repeatedly run a simulation multiple times.
  
  ## Options
  
      `:collect_all?` - If true, collects data from each individual simulation run and returns this an the third element of the result tuple
  
  """
  @spec mc(iterations :: pos_integer, fun :: ((pos_integer) -> float), options :: Keyword.t) :: {avg :: float, sd :: float, tries :: [float]} | {avg :: float, sd :: float}
  def mc(iterations,fun,options \\ []) do
    all? = options[:collect_all?] || false

    tries = 1..iterations |> Enum.map(fn _ -> fun.() end)
    avg = moment tries, 1
    sd = :math.sqrt momentc(tries,2,avg)
    if all?, do: {avg,sd,tries}, else: {avg,sd}
  end

  @hours 24.0
  @default_workday {8.0,18.0}
  @default_epoch ~D[1970-01-01]

  @doc ~s"""
  Adjusts the times to working hours and/or work days.
  
  ## Options
  
      `workhours` - a 2-tuple containing the starting and ending hours of the work day (defaults
          to #{inspect @default_workday})
      `epoch` - the epoch to which all data elements are relative (defaults to #{@default_epoch})
      `saturday` - number of days since the epoch that corresponds to a Saturday (defaults
          to #{13 - Date.day_of_week(@default_epoch)})
      `correct` - whether to correct the times for working hours and weekdays; possible values
          `:worktime`, `:weekday`, `:"weekday+worktime"` (defaults to `false`)

  """
  @spec adjust_times(Enumerable.t, options :: Keyword.t) :: Enumerable.t
  def adjust_times(data, options) do
    {startofday,endofday} = options[:workhours] || @default_workday
    correct = options[:correct] || false
    epoch = options[:epoch] || @default_epoch
    sat = 13 - Date.day_of_week(epoch)
    saturday = options[:saturday] || sat

    data
    |> Stream.map(fn x ->
        case correct do
          :worktime -> map2workhours(x, startofday, endofday)
          :weekday -> map2weekdays(x, saturday)
          :"weekday+worktime" -> x |> map2workhours(startofday, endofday) |> map2weekdays(saturday)
          _ -> x
        end
      end)
    |> Enum.sort(&(&1>&2)) # Sort on new delivery times
  end

  @default_cutoff 0.01

  @doc """
  Returns a list of time differences (assumes an ordered list as input)
  
  ## Options
  
      `cutoff` - time differences below the cutoff are changed to the cutoff value (defaults to `#{@default_cutoff}`)
      `drop?` - whether to drop time differences below the cutoff (defaults to `false`)

  """
  @spec time_diff(data :: Enumrable.t, options :: Keyword.t) :: Enumerable.t
  def time_diff(data, options) do
    cutoff = options[:cutoff] || @default_cutoff
    drop? = options[:drop] || false

    data
    |> Stream.chunk_every(2,1,:discard)
    |> Stream.map(fn [x,y]->x-y end)
    |> Stream.transform(nil,fn x,_acc ->
          {
            cond do
              x < cutoff and drop? -> []
              x < cutoff -> [cutoff]
              true -> [x]
            end,
            nil
          }
        end)
    |> (& if is_function(data, 2), do: &1, else: Enum.into(&1, [])).()
  end

  @doc """
  Calculates the systematic errors for bins due to uncertainty in assigning data to bins.
  
  ## Options
  
      `bin` - the size of bins to use (defaults to 1)
      `iterations` - the number of iterations to use to estimate the error due to noise (defatuls to 100)

  """
  @spec binerror(data :: [number], noise_fun :: ((Enumerable.t) -> Enumerable.t), options :: Keyword.t) :: [{bin :: number, avg :: number, error :: number}]
  def binerror(data, noise_fun, options \\ []) do
    binsize = options[:bin] || 1
    iterations = options[:iterations] || 100

    1..iterations
    |> Stream.map(fn _ ->
        data
        |> noise_fun.()
        |> to_bins({binsize,0})
        |> Stream.map(fn {x,y,low,high}->{x,[{y,low,high}]} end)
        |> Map.new()
      end)
    |> Enum.reduce(%{}, fn map,acc -> Map.merge(map,acc, fn _k, v1,v2 -> v1++v2 end) end)
    |> Stream.map(fn {k,list} ->
        {ys,lows,highs} = unzip list
        avg = moment ys,1
        avg_low = moment lows,1
        avg_high = moment highs,1
        sd = :math.sqrt momentc ys,2,avg
        {k,avg,avg_low,avg_high,sd}
      end)
    |> Stream.map(fn {x,y,ylow,yhigh,err} ->
        {
          x,
          y,
          max(0.0,y-:math.sqrt((y-ylow)*(y-ylow)+err*err)),
          min(1.0,y+:math.sqrt((yhigh-y)*(yhigh-y)+err*err))
        }
      end)
    |> Enum.sort(fn t1,t2 -> elem(t1,0)<elem(t2,0) end)
  end

  @doc """
  Displays results of the function `Chi2fit.Fit.chi2probe/4`
  """
  @spec display(device :: IO.device(), F.chi2probe() | avgsd()) :: none()
  def display(device \\ :stdio, results)
  def display(device,{chi2, parameters,errors,_saved}) do
      IO.puts device,"Initial guess:"
      IO.puts device,"    chi2:\t\t#{chi2}"
      IO.puts device,"    pars:\t\t#{inspect parameters}"
      IO.puts device,"    ranges:\t\t#{inspect errors}\n"
  end
  def display(device,{avg,sd,direction}) do
    op = case direction do
      :+ -> &Kernel.+/2
      :- -> &Kernel.-/2
    end
    IO.puts device,"50%    => #{:math.ceil(avg)} units"
    IO.puts device,"84%    => #{:math.ceil(op.(avg,sd))} units"
    IO.puts device,"97.5%  => #{:math.ceil(op.(avg,2*sd))} units"
    IO.puts device,"99.85% => #{:math.ceil(op.(avg,3*sd))} units"
  end
  
  @doc """
  Displays results of the function `Chi2fit.Fit.chi2fit/4`
  """
  @spec display(device :: IO.device(),hdata :: ecdf(),model :: D.model(),F.chi2fit(),options :: Keyword.t) :: none()
  def display(device \\ :stdio,hdata,model,{chi2, cov, parameters, errors},options) do
      param_errors = cov |> M.diagonal |> Enum.map(fn x->x|>abs|>:math.sqrt end)

      IO.puts device,"Final:"
      IO.puts device,"    chi2:\t\t#{chi2}"
      IO.puts device,"    Degrees of freedom:\t#{length(hdata)-Distribution.size(model)}"
      IO.puts device,"    gradient:\t\t#{inspect jacobian(parameters,&F.chi2(hdata,fn x->Distribution.cdf(model).(x,&1) end,fn _->0.0 end,options),options)}"
      IO.puts device,"    parameters:\t\t#{inspect parameters}"
      IO.puts device,"    errors:\t\t#{inspect param_errors}"
      IO.puts device,"    ranges:"
      puts_errors device,errors
  end

  @doc """
  Pretty prints subsequences.
  """
  @spec display_subsequences(device :: IO.device(), trends :: list(), intervals :: [NaiveDateTime.t]) :: none()
  def display_subsequences(device \\ :stdio, trends, intervals) do
    trends
    |> Stream.transform(1, fn arg={_,_,data}, index -> { [{arg, Enum.at(intervals,index)}], index+length(data)} end)
    |> Stream.each(fn
        {{chimin, [rate], subdata},date} ->
          IO.puts device, "Subsequence ending @#{Timex.format!(date,~S({Mshort}, {D} {YYYY}))}"
          IO.puts device, "----------------------------------"
          IO.puts device, "    chi2@minimum:  #{Float.round(chimin,1)}"
          IO.puts device, "    delivery rate: #{Float.round(rate,1)}"
          IO.puts device, "    subsequence:   #{inspect(subdata, charlists: :as_lists)}"
          IO.puts device, ""
      end)
    |> Stream.run()
  end

  @doc """
  Maps the time of a day into the working hour period
  
  Scales the resulting part of the day between 0..1.
  
  ## Arguments
  
      `t` - date and time of day as a float; the integer part specifies the day and the fractional part the hour of the day
      `startofday` - start of the work day in hours
      `endofday` - end of the working day in hours
  
  ## Example
  
      iex> map2workhours(43568.1, 8, 18)
      43568.0
  
      iex> map2workhours(43568.5, 8, 18)
      43568.4
  
  """
  @spec map2workhours(t :: number, startofday :: number, endofday :: number) :: number
  def map2workhours(t,startofday,endofday)
    when startofday>0 and startofday<endofday and endofday<@hours do
    frac = t - trunc(t)
    hours = endofday - startofday
    trunc(t) + min(max(0.0,frac - startofday/@hours),hours/@hours) * @hours/hours
  end

  @doc """
  Maps the date to weekdays such that weekends are eliminated; it does so with respect to a given Saturday
  
  ## Example
  
      iex> map2weekdays(43568.123,43566)
      43566.123
  
      iex> map2weekdays(43574.123,43566)
      43571.123
  
  """
  @spec map2weekdays(t :: number, sat :: pos_integer) :: number
  def map2weekdays(t, sat) when is_integer(sat) do
    offset = rem trunc(t)-sat, 7
    weeks = div trunc(t)-sat, 7
    
    part_of_day = t - trunc(t)
    sat + 5*weeks + max(0.0,offset-2.0) + part_of_day
  end

  @doc """
  Walks a map structure while applying the function `fun`.
  """
  @spec analyze(map :: %{}, fun :: (([number],Keyword.t) -> Keyword.t), options :: Keyword.t) :: Keyword.t
  def analyze(map = %{}, fun, options) do
      map |> Enum.reduce(%{}, fn {k,v},acc -> Map.put(acc,k,analyze(v,fun,options)) end)
  end
  def analyze(data, fun, options) when is_list(data) do
      cond do
          Keyword.keyword?(data) ->
              Keyword.merge(data, fun.(data,Keyword.put(options,:bin,data[:bin])))
          true ->
              analyze([throughput: data, bin: options[:bin]], fun, options)
      end
  end

  @doc """
  Pretty-prints a nested array-like structure (list or tuple) as a table.
  """
  @spec as_table(rows :: [any], header :: tuple()) :: list()
  def as_table(rows, header) do
    map = 1..tuple_size(header) |> Enum.map(&{&1,0}) |> Enum.into(%{})
    table = [header|rows] |> _to_string()
    map = Enum.reduce(table, map, fn row,acc ->
      row
      |> Enum.with_index(1)
      |> Enum.reduce(acc, fn {str,i},acc2 -> Map.update!(acc2, i, fn v -> max(v,String.length(str)) end) end)
    end)
    table
    |> Enum.with_index()
    |> Enum.each(fn
      {row, 0} ->
        IO.puts row |> Enum.with_index(1) |> Enum.map(fn {str,i} -> String.pad_trailing(str, map[i]) end) |> Enum.join("|")
        IO.puts row |> Enum.with_index(1) |> Enum.map(fn {_,i} -> String.duplicate("-",map[i]) end) |> Enum.join("|")
      {row, _} ->
        IO.puts row |> Enum.with_index(1) |> Enum.map(fn {str,i} -> String.pad_trailing(str, map[i]) end) |> Enum.join("|")
    end)
    
    rows
  end
  
  defp _to_string(list) when is_list(list), do: list |> Enum.map(&_to_string/1)
  defp _to_string(tuple) when is_tuple(tuple), do: tuple |> Tuple.to_list |> _to_string()
  defp _to_string(string) when is_binary(string), do: string
  defp _to_string(float) when is_float(float), do: "#{float}"
  defp _to_string(integer) when is_integer(integer), do: "#{integer}"

  @doc ~S"""
  Reads CSV data, extracts one column, and returns it as a list of `NaiveDateTime`.
  
  ## Examples
  
      iex> csv = ["Done","2019/05/01","2019/06/01"] |> Stream.map(& &1)
      ...> csv_to_list csv, "Done", header?: true
      [~N[2019-06-01 00:00:00], ~N[2019-05-01 00:00:00]]
  
      iex> csv = ["Done","2019/May/01","2019/Jun/01"] |> Stream.map(& &1)
      ...> csv_to_list csv, "Done", header?: true, format: "{YYYY}/{Mshort}/{0D}"
      [~N[2019-06-01 00:00:00], ~N[2019-05-01 00:00:00]]
  
      iex> csv = ["Done","2019/May/01","2019/06/01"] |> Stream.map(& &1)
      ...> csv_to_list csv, "Done", header?: true, format: "{YYYY}/{Mshort}/{0D}"
      [~N[2019-05-01 00:00:00]]
  
      iex> csv = ["Done","2019/May/01","2019/06/01"] |> Stream.map(& &1)
      ...> csv_to_list csv, "Done", header?: true, format: ["{YYYY}/{Mshort}/{0D}","{YYYY}/{0M}/{0D}"]
      [~N[2019-06-01 00:00:00], ~N[2019-05-01 00:00:00]]

      iex> csv = ["Done","2019/May/01","2019/Jun/01"] |> Stream.map(& &1)
      ...> csv_to_list csv, "Done", header?: true, format: ["%Y/%b/%d"], parser: :strftime
      [~N[2019-06-01 00:00:00], ~N[2019-05-01 00:00:00]]

  """
  @spec csv_to_list(csvcata :: Enumerable.t, key :: String.t, options :: Keyword.t) :: [NaiveDateTime.t]
  def csv_to_list(csvdata, key, options \\ []) do
    header? = options[:header?] || false
    format = options[:format] || "{YYYY}/{0M}/{0D}"
    separator = options[:separator] || ?,
    parser = case options[:parser] do
      :strftime -> Timex.Format.DateTime.Formatters.Strftime
      :default -> Timex.Format.DateTime.Formatters.Default
      nil -> Timex.Format.DateTime.Formatters.Default
      _ -> Timex.Format.DateTime.Formatters.Default
    end

    formats = if is_list(format), do: format, else: [format]

    csvdata
    |> CSV.decode!(headers: header?, separator: separator)
    |> Stream.filter(& Map.fetch!(&1, key) != "")
    |> Stream.map(& Map.fetch!(&1, key))
    |> Stream.map(fn t -> Enum.reduce_while(formats,nil,fn f,_acc ->
        case Timex.parse(t, f, parser) do
          {:ok,n} -> {:halt,n}
          {:error,_msg} -> {:cont,nil}
        end
      end)
    end)
    |> Stream.filter(& &1 != nil)
    |> Enum.sort(& NaiveDateTime.compare(&1,&2) === :gt)
  end

  @doc """
  Returns a `Stream` that generates a stream of dates.
  
  ## Examples
  
      iex> intervals(end: ~D[2019-06-01]) |> Enum.take(4)
      [~D[2019-06-01], ~D[2019-05-16], ~D[2019-05-01], ~D[2019-04-16]]
  
      iex> intervals(end: ~D[2019-06-01], type: :weekly) |> Enum.take(4)
      [~D[2019-06-01], ~D[2019-05-18], ~D[2019-05-04], ~D[2019-04-20]]
  
      iex> intervals(end: ~D[2019-06-01], type: :weekly, weeks: 1) |> Enum.take(4)
      [~D[2019-06-01], ~D[2019-05-25], ~D[2019-05-18], ~D[2019-05-11]]

      iex> intervals(end: ~D[2019-06-01], type: :weekly, weeks: [3,2]) |> Enum.take(4)
      [~D[2019-06-01], ~D[2019-05-11], ~D[2019-04-27], ~D[2019-04-13]]

  """
  @spec intervals(options :: Keyword.t) :: Stream.t
  def intervals(options \\ []) do
    type = options[:type] || :half_month
    periods = case options[:weeks] do
      nil -> [2]
      x when is_number(x) -> [x]
      list when is_list(list) -> list
    end
    last = options[:end] || Date.utc_today()

    case type do
      :half_month ->
        recent = case last do
          date = %Date{day: day} when day > 16 -> %Date{date | day: 1, month: date.month+1}
          date = %Date{day: day} when day > 1 -> %Date{date | day: 16}
          date = %Date{day: 1} -> date
        end

        recent |> Stream.iterate(fn
          previous = %Date{day: 16} -> Timex.shift previous, days: -15
          previous = %Date{day: 1} -> Timex.shift previous, days: +15, months: -1
        end)
      
      :weekly ->
        Stream.resource(
          fn -> {last,periods} end,
          fn
            {current,[p]} ->
              next = Timex.shift(current, weeks: -p)
              {[current], {next,[p]}}
            {current,[p|rest]} ->
              next = Timex.shift(current, weeks: -p)
              {[current], {next,rest}}
          end,
          fn _ -> [] end)
    end
  end

  @doc """
  Counts the number of dates (`datelist`) that is between consecutive dates in `intervals` and returns the result as a list of numbers.
  """
  @spec throughput(intervals :: Enumerable.t, datelist :: [NaiveDateTime.t]) :: [number]
  def throughput(intervals, datelist) do
    intervals
    |> Stream.chunk_every(2, 1, :discard)
    |> Stream.transform(datelist, fn
        _, [] ->
          {:halt, []}
        [d1,d2], acc ->
          {left,right} = Enum.split_with(acc, fn d -> Timex.between?(d,d2,d1) end)
          {[{d1,Enum.count(left)}],right}
      end)
    |> Enum.map(fn {_d,count} -> count end)
  end

  @doc ~S"""

  ## Examples

      iex> subsequences []
      []
  
      iex> subsequences [:a, :b]
      [[:a], [:a, :b]]
  
      iex> Stream.cycle([1,2,3]) |> subsequences |> Enum.take(4)
      [[1], [1, 2], [1, 2, 3], [1, 2, 3, 1]]

  """
  @spec subsequences(Enumerable.t) :: Enumerable.t
  def subsequences(stream) when is_function(stream, 2) do
    stream
    |> Stream.transform([], fn x,acc -> {[Enum.reverse([x|acc])], [x|acc]} end)
  end
  def subsequences(list) do
    {result, _} = list
    |> Enum.reduce({[],[]}, fn  x, {res,acc} -> {[Enum.reverse([x|acc])|res], [x|acc]} end)
    Enum.reverse(result)
  end

end
