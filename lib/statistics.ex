defmodule Chi2fit.Statistics do

  # Copyright 2015-2021 Pieter Rijken
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

  import Chi2fit.FFT
  import Chi2fit.Utilities

  alias Chi2fit.Fit, as: F
  
  @typedoc "Algorithm used to assign errors to frequencey data: Wald score and Wilson score."
  @type algorithm :: :wilson | :wald

  @typedoc "Cumulative Distribution Function"
  @type cdf :: ((number)->{number, number, number})

  @typedoc "Binned data with error bounds specified through low and high values"
  @type ecdf :: [{float, float, float, float}]

  @doc """
  Converts a list of numbers to frequency data.
  
  The data is divided into bins of size `binsize` and the number of data points inside a bin are counted. A map
  is returned with the bin's index as a key and as value the number of data points in that bin.

  The function returns a list of 2-tuples. Each tuple contains the index of the bin and the value of the count of the
  number of items in the bin. The index of the bins start at 1 in the following way:

    * [0..1) has index 1 (including 0 and excludes 1),
    * [1..2) has index 2,
    * etc.

  When an offset is used, the bin starting from the offset, i.e. [offset..offset+1) gets index 1. Values less than
  the offset are gathered in a bin with index 0.
  
  ## Examples
  
      iex> make_histogram [1,2,3]
      [{2, 1}, {3, 1}, {4, 1}]
      
      iex> make_histogram [1,2,3], 1.0, 0
      [{2, 1}, {3, 1}, {4, 1}]
      
      iex> make_histogram [1,2,3,4,5,6,5,4,3,4,5,6,7,8,9]
      [{2, 1}, {3, 1}, {4, 2}, {5, 3}, {6, 3}, {7, 2}, {8, 1}, {9, 1}, {10  , 1}]
      
      iex> make_histogram [1,2,3,4,5,6,5,4,3,4,5,6,7,8,9], 3, 1.5
      [{0, 1}, {1, 6}, {2, 6}, {3, 2}]

      iex> make_histogram [0,0,0,1,3,4,3,2,6,7],1
      [{1,3},{2,1},{3,1},{4,2},{5,1},{7,1},{8,1}]

      iex> make_histogram [0,0,0,1,3,4,3,2,6,7],1,0.5
      [{0,3},{1,1},{2,1},{3,2},{4,1},{6,1},{7,1}]


  """
  @spec make_histogram([number], number, number) :: [{non_neg_integer, pos_integer}]
  def make_histogram(list, binsize \\ 1.0, offset \\ 0.0)
  def make_histogram(list, binsize, offset) when binsize > offset do
    Enum.reduce(list, %{}, fn
      (number, acc) ->
        acc |> Map.update(if(number < offset, do: 0, else: trunc(Float.floor(1 + (number - offset) / binsize))), 1, &(1 + &1))
    end) |> Enum.reduce([], fn (pair, acc) -> [pair|acc] end) |> Enum.sort_by(fn {k, _v} -> k end)
  end
  def make_histogram(_list, _binsize, _offset), do: raise ArgumentError, message: "binsize must be larger than bin offset"


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
  @spec empirical_cdf([{float, number}], {number, number}, algorithm, integer) :: {cdf, bins :: [float], numbins :: pos_integer, sum :: float}
  def empirical_cdf(data, bin \\ {1.0, 0.5}, algorithm \\ :wilson, correction \\ 0)
  def empirical_cdf(data, {binsize, offset}, algorithm, correction) do
    {bins, sum} = data
      |> Enum.sort(fn ({x1, _}, {x2, _}) -> x1 < x2 end)
      |> Enum.reduce({[], 0}, fn ({n, y}, {acc, sum}) -> {[{offset + binsize * n, y + sum}|acc], sum + y} end)

    normbins = bins
      |> Enum.reverse
      |> Enum.map(fn ({x, y}) -> {x, y / (sum + correction), y} end)

    {normbins |> to_cdf_fun(sum, algorithm),
     normbins,
     length(bins),
     sum}
  end

  @doc """
  Calculates the empirical CDF from a sample.
  
  Convenience function that chains `make_histogram/2` and `empirical_cdf/3`.
  """
  @spec get_cdf([number], number|{number, number}, algorithm, integer) :: {cdf, bins :: [float], numbins :: pos_integer, sum :: float}
  def get_cdf(data, binsize \\ {1.0, 0.5}, algorithm \\ :wilson, correction \\ 0)
  def get_cdf(data, {binsize, offset}, algorithm, correction) do
    data
    |> make_histogram(binsize, offset)
    |> empirical_cdf({binsize, offset}, algorithm, correction)
  end

  defp get_add_cdf(data, n, binsize, algorithm \\ :wilson, correction \\ 0) when n > 1 do
    data
    |> List.duplicate(n)
    |> Enum.reduce([], fn
        list, [] -> list
        list, acc -> Enum.flat_map(acc, fn d1 -> Enum.map(list, & d1+&1) end)
    end)
    |> get_cdf(binsize, algorithm, correction)
  end

  defp get_max_cdf(data, n, binsize, algorithm \\ :wilson, correction \\ 0) when n > 1 do
    data
    |> List.duplicate(n)
    |> Enum.reduce([], fn
        list, [] -> list
        list, acc -> Enum.flat_map(acc, fn d1 -> Enum.map(list, & max(d1,&1)) end)
    end)
    |> get_cdf(binsize, algorithm, correction)
  end

  defp statistic(data, fun1, fun2) do
    data
    |> Enum.map(&elem(&1, 0))
    |> Enum.map(&{1 - elem(fun1.(&1), 0), 1 - elem(fun2.(&1), 0)})
    |> Enum.filter(fn {d1, d2} -> d1 > 0 and d2 > 0 end)
    |> Enum.map(fn {d1, d2} -> d1 / d2 end)
    |> Enum.take(-1)
    |> hd
  end

  @doc """
  Calculates the test statistic for subexponentiality of a sample.

  A value close to 0 is a strong indication that the sample shows subexponential behaviour (extremistan), i.e. is fat-tailed.
  """
  def subexponential_stat(data, test \\ :sum, n \\ 2, binsize \\ {1, 0})
  def subexponential_stat(data, :sum, n, binsize) do
    {cdf_fun, _, _, _} = data |> get_cdf(binsize)
    {cdf_add_fun, cdf_emp, _, _} = data |> get_add_cdf(n, binsize)

    test = statistic(cdf_emp, cdf_add_fun, cdf_fun)
    abs(test - n)
  end
  def subexponential_stat(data, :max, n, binsize) do
    {cdf_max_fun, _, _, _} = data |> get_max_cdf(n, binsize)
    {cdf_add_fun, cdf_emp, _, _} = data |> get_add_cdf(n, binsize)

    test = statistic(cdf_emp, cdf_add_fun, cdf_max_fun)
    abs(test - 1)
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
  @type range :: {float, float} | [float, ...]
  @spec convert_cdf({cdf, range}) :: [{float, float, float, float}]
  def convert_cdf({cdf, {mindur, maxdur}}), do: round(mindur)..round(maxdur) |> y_with_errors(cdf)
  def convert_cdf({cdf, categories}) when is_list(categories), do: categories |> y_with_errors(cdf)
  defp y_with_errors(list, cdf), do: list |> Enum.map(&Tuple.insert_at(cdf.(&1), 0,& 1)) 

  @doc """
  Converts raw data to binned data with (asymmetrical) errors.
  """
  @spec to_bins(data :: [number], binsize :: {number, number}) :: ecdf()
  def to_bins(data, binsize \\ {1.0, 0.5}) do
    # Convert the raw data to binned data (histogram or frequency data):
    {cdf, bins, _, _} = get_cdf data, binsize

    # Add the errors based on the binomial distribution (Wilson score):
    convert_cdf {cdf, bins|>Enum.map(&elem(&1, 0))}
  end
  
  @doc """
  Calculates the nth moment of the sample.
  
  ## Example
  
      iex> moment [1,2,3,4,5,6], 1
      3.5
  """
  @spec moment(sample::[number], n::pos_integer) :: float
  def moment(sample, n) when length(sample)>0 and is_integer(n) and n>0 do
    (sample |> Stream.map(fn x-> :math.pow(x, n) end) |> Enum.sum)/length(sample)
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
  def cullen_frey(sample, n \\ 100) do
    bootstrap(n, sample,
      fn
        data, _i ->
          mean = data |> moment(1)
          sigma = :math.sqrt(data |> momentc(2))
          skewness = data |> momentn(3, mean, sigma)
          kurtosis = data |> momentn(4, mean, sigma)
          {skewness*skewness, kurtosis-3.0}
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

end
