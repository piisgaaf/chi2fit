defmodule Chi2fit.Utilities do

  # Copyright 2017 Pieter Rijken
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

  @h 1.0e-10

  require Integer
  require Logger
  import Kernel, except: [*: 2, /: 2,+: 2, -: 2]

  alias Decimal, as: D

  @typedoc "A real number."
  @type real :: number
  
  @typedoc "A complex number with a real part and an imaginary part."
  @type complex :: {real,real}
  
  @typedoc "Cumulative Distribution Function"
  @type cdf :: ((number)->{number,number,number})

  @typedoc "Algorithm used to assign errors to frequencey data: Wald score and Wilson score."
  @type algorithm :: :wilson | :wald

  defp {x1,x2} * {y1,y2}, do: {x1*y1-x2*y2,x1*y2+x2*y1}
  defp x * {y1,y2}, do: {x*y1,x*y2}
  defp x * y, do: Kernel.*(x,y)

  defp {x1,x2} / y, do: {x1/y,x2/y}
  defp x / y, do: Kernel./(x,y)

  defp {x1,x2} + {y1,y2}, do: {x1+y1,x2+y2}
  defp x + {y1,y2}, do: {x+y1,y2}
  defp {x1,x2} + y, do: {x1+y,x2}
  defp x + y, do: Kernel.+(x,y)

  defp {x1,x2} - {y1,y2}, do: {x1-y1,x2-y2}
  defp x - {y1,y2}, do: {x-y1,-y2}
  defp x - y, do: Kernel.-(x,y)

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
          splus = t*t - 1/numpoints + 4*numpoints*y*(1-y) + (4*y - 2)
          smin =  t*t - 1/numpoints + 4*numpoints*y*(1-y) - (4*y - 2)
          srtplus = 1.0 + t*:math.sqrt(splus)
          srtmin =  1.0 + t*:math.sqrt(smin)
          ylow =  max(0.0, (2*numpoints*y + t*t - srtplus)/2/(numpoints + t*t))
          yhigh = min(1.0, (2*numpoints*y + t*t + srtmin )/2/(numpoints + t*t))
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
  @spec empirical_cdf([{float,number}],integer,algorithm) :: {cdf,[float],pos_integer,float}
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
  @spec get_cdf([number], number, algorithm) :: {cdf,[float],pos_integer,float,[number]}
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

  @doc """
  Returns the real roots of polynoms of order 1, 2 and 3 as a list.
  
  ## Examples
  
      Solve `2.0*x + 5.0 = 0`
      iex> solve [2.0,5.0]
      [-2.5]
      
      iex> solve [2.0,-14.0,24.0]
      [4.0,3.0]
      
      iex> solve [1.0,0.0,5.0,6.0]
      [-0.9999999999999999]
  """
  @spec solve([float]) :: [float]
  def solve([0.0|rest]), do: solve rest

  def solve([a1,a0]), do: [-a0/a1]

  def solve([a2,a1,a0]) do
    sqr = a1*a1-4*a2*a0
    cond do
      sqr == 0 -> -a1/2/a2
      sqr > 0 -> [(-a1+:math.sqrt(sqr))/2/a2,(-a1-:math.sqrt(sqr))/2/a2]
      true -> []
    end
  end

  def solve([1.0,0.0,p,q]) do
    ## For details see equations (83) and (84) in http://mathworld.wolfram.com/CubicFormula.html
    c = -0.5*q*:math.pow(3/abs(p),1.5)
    cond do
      p>0 ->
        [:math.sinh(1.0/3.0*:math.asinh(c))]
      c>=1 ->
        [:math.cosh(1.0/3.0*:math.acosh(c))]
      c<=-1 ->
        [-:math.cosh(1.0/3.0*:math.acosh(abs(c)))]
      true ->
        ## Three real solutions
        [:math.cos(1.0/3.0*:math.acos(c)),:math.cos(1.0/3.0*:math.acos(c) + 2*:math.pi()/3.0),:math.cos(1.0/3.0*:math.acos(c) + 4*:math.pi()/3.0)]
    end
    |> Enum.map(&(&1*2*:math.sqrt(abs(p)/3.0)))
  end
  def solve([1.0,a2,a1,a0]), do: solve([1.0,0.0,(3*a1-a2*a2)/3.0,(2*a2*a2*a2-9*a1*a2+27*a0)/27.0]) |> Enum.map(&(&1-a2/3.0))
  def solve([a3,a2,a1,a0]), do: solve([1.0,a2/a3,a1/a3,a0/a3])

  defp expand_pars(list,h) do
    list |> Enum.map(
          fn
            ({{x,0}}) -> {x}
            ({{x,n}}) -> List.flatten expand_pars([{{D.add(x,h),n-1}},{x,n-1}],h)
            ({x,0}) -> x
            ({x,n}) -> List.flatten expand_pars([{D.add(x,h),n-1},{{x,n-1}}],h)
            (x) -> x
          end)
  end

  defp reduce_pars(list) do
    list |> Enum.reduce([{[],1}],
      fn
        (list,acc) when is_list(list) ->
          Enum.flat_map(list,
            fn
              ({x}) -> Enum.map(acc, fn ({y,n})->{[x|y],-n} end)
              (x) -> Enum.map(acc, fn ({y,n})->{[x|y],n} end)
            end)
        (x,acc) -> Enum.map(acc, fn ({y,n})->{[x|y],n} end)
      end)
      |> Enum.map(fn ({l,n}) -> {Enum.reverse(l),n} end)
  end

  defp pow(_d,0), do: D.new(1)
  defp pow(d,n) when is_integer(n) and n > 0, do: D.mult(d,pow(d,n-1))

  defp to_decimals(list) do
    list |> Enum.map(
      fn
        {x,n} -> {D.new(x),n}
        x -> D.new(x)
      end)
  end

  @doc """
  Calculates the partial derivative of a function and returns the value.
  
  ## Examples

      The function value at a point:
      iex> der [3.0], fn [x]-> D.mult(x,x) end
      #Decimal<9.0>

      The first derivative of a function at a point:
      iex> der [{3.0,1}], fn [x]-> D.mult(x,x) end
      #Decimal<6.0000000001>

      The second derivative of a function at a point:
      iex> der [{3.0,2}], fn [x]-> D.mult(x,x) end
      #Decimal<2>

      Partial derivatives with respect to two variables:
      iex> der [{2.0,1},{3.0,1}], fn [x,y] -> D.mult(D.new(3),D.mult(x,D.mult(x,y))) end
      #Decimal<12.000>

  """
  @spec der([float|{float,integer}], (([D.t])->D.t), Keyword.t) :: float
  def der(parameters, fun, options \\ []) do
    h = D.new(options[:h] || @h)

    factor = Enum.reduce(parameters,D.new(1.0),
      fn
        ({_x,n},acc) -> D.mult(acc,pow(h,n))
        (_x,acc) -> acc
      end)
    parameters
    |> to_decimals
    |> expand_pars(h)
    |> reduce_pars
    |> Enum.reduce(D.new(0), fn ({x,n},sum) when is_list(x) -> D.add(sum,D.mult(D.new(n),fun.(x))) end)
    |> D.div(factor)
  end

  defp jacobian(x=[_|_], k, fun) when k>0 and k<=length(x) and is_function(fun,1) do
    x |> List.update_at(k-1, fn (val) -> {val,1} end) |> der(fun)
  end

  @doc """
  Calculates the jacobian of the function at the point `x`.
  
  ## Examples
  
      iex> jacobian [2.0,3.0], fn [x,y] -> D.mult(x,y) end
      [D.new(3), D.new(2)]

  """
  @spec jacobian(x :: [float], (([D.t])->D.t)) :: [D.t]
  def jacobian(x, fun) do
    jacfun = &(jacobian(x, &1, fun))
    Enum.reduce(length(x)..1, [], fn (k,acc) -> [jacfun.(k)|acc] end)
  end

  defp weight(r,m,n), do: weight(r*m,n)
  defp weight(rm,n), do: weight(rm/n)
  defp weight(x), do: {:math.cos(2*:math.pi()*x),-:math.sin(2*:math.pi()*x)}

  defp split_evenodds(list) when Integer.is_even(length(list)) do
    list
    |> List.foldr({[[],[]],false},
      fn
        (item,{[e,o],true}) -> {[[item|e],o],false}
        (item,{[e,o],false}) -> {[e,[item|o]],true}
      end)
    |> elem(0)
  end

  @doc """
  Calculates the discrete Fast Fourier Transform of a list of numbers.
  
  Provides a parallel version (see options below). See [1] for details of the algorithm implemented.
  
  ## Options
    * `:phase` - Correction factor to use in the weights of the FFT algorithm. Defaults to 1.
    * `:nproc` - Parellel version. Number of processes to use. See [2]. Defaults to 1.

  ## References
  
    [1] Zie: https://en.wikipedia.org/wiki/Cooley%E2%80%93Tukey_FFT_algorithm
  
    [2] Parallel version of FFT; see http://www.webabode.com/articles/Parallel%20FFT%20implementations.pdf
  
  ## Examples
    
      iex> fft [4]
      [{4.0, 0.0}]
    
      iex> fft [1,2,3,4,5,6]
      [{21.0, 0.0}, {-3.0000000000000053, 5.19615242270663},
                {-3.0000000000000036, 1.7320508075688736}, {-3.0, 0.0},
                {-2.9999999999999982, -1.7320508075688799},
                {-2.999999999999991, -5.196152422706634}]

  """
  @spec fft([real],opts :: Keyword.t) :: [complex]
  def fft(list,opts \\ [])

  def fft([],_opts), do: []
  def fft([x,y],opts) do
    fac = opts[:phase] || 1
    [x*weight(fac*0,0,2)+y*weight(fac*0,1,2),x*weight(fac*1,0,2)+y*weight(fac*1,1,2)]
  end
  def fft(list=[_|_],opts) do
    fac = opts[:phase] || 1
    nproc = opts[:nproc] || 1

    nn = length(list)
    cond do
      Integer.is_even(length(list)) ->
        zipped = cond do
          nproc == 2 or nproc == 4 ->
            list
            |> split_evenodds
            |> Enum.map(fn x-> Task.async(fn -> fft(x,Keyword.merge(opts,[nproc: nproc/2])) end) end)
            |> Task.yield_many(3_600_000)
            |> Enum.map(fn ({_task,{:ok,result}})->result end)
            |> (&(apply(fn x,y->Stream.zip(x,y) end,&1))).()
        
          nproc == 1 ->
            list
            |> split_evenodds
            |> Enum.map(fn arg->fft(arg,opts) end)
            |> (&(apply(fn x,y->Stream.zip(x,y) end,&1))).()
        end

        n = nn/2
        zipped
        |> Stream.concat(zipped)
        |> Stream.with_index(0)
        |> Stream.map(
          fn
            ({{x,y},m}) when m<n -> x + (weight(fac*1,m,2*n)*y)
            ({{x,y},m}) when m>=n -> x - (weight(fac*1,m-n,2*n)*y)
          end)
        |> Enum.to_list
        
      true ->
        0..nn-1 |> Enum.map(
          fn m ->
            list |> Stream.with_index(0) |> Stream.map(fn ({item,k})-> item*weight(fac*m,k,nn) end) |> Enum.reduce(0,fn (x,acc)->x+acc end)
          end)
    end
  end
  
  @doc """
  Calculates the inverse FFT.
  
  For available options see `fft/2`.

  ## Examples
  
      iex> ifft [4.0]
      [{4.0, 0.0}]
  
      iex> ifft [1.0,2.0,3.0]
      [{2.0, 0.0}, {-0.5000000000000003, -0.2886751345948125},
                  {-0.4999999999999995, 0.28867513459481353}]

      iex> [1.0,5.0] |> fft |> ifft
      [{1.0, -3.061616997868383e-16}, {5.0, 6.123233995736767e-17}]

  """
  @spec ifft([real],Keyword.t) :: [complex]
  def ifft(list,opts \\ [nproc: 1]) do
    n = length(list)
    list |> fft(Keyword.merge(opts,[phase: -1])) |> Enum.map(&(&1/n))
  end

  @doc """
  Calculates the norm of a complex number or list of complex numbers.
  
  ## Examples
  
      iex> normv []
      []
      
      iex> normv {2,3}
      13
      
      iex> normv [{2,3},{1,2}]
      [13,5]

  """
  @spec normv([complex]|complex) :: real
  def normv({x,y}), do: x*x+y*y
  def normv(list) when is_list(list), do: list |> Enum.map(&normv/1)

  @doc """
  Converts the input so that the result is a Puiseaux diagram, that is a strict convex shape.
  
  ## Examples
  
      iex> puiseaux [1]
      [1]
      
      iex> puiseaux [5,3,3,2]
      [5, 3, 2.5, 2]

  """
  @spec puiseaux([real],[real],boolean) :: [real]
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
  
  For available options see `fft/2`. Returns a list of the autocorrelation coefficients.
  
  ## Example
  
      iex> auto [1,2,3]
      [14.0, 7.999999999999999, 2.999999999999997]

  """
  @spec auto([real],Keyword.t) :: [real]
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
  @spec error([{gamma :: number,k :: pos_integer}], :initial_sequence_method) :: {number, number}
  def error(nauto, :initial_sequence_method) do
    ## For reversible Markov Chains
    gamma = nauto |> Stream.chunk(2) |> Stream.map(fn ([{x,k},{y,_}])->{k/2,x+y} end) |> Enum.to_list
    gamma0 = nauto |> Stream.take(1) |> Enum.to_list |> (&(elem(hd(&1),0))).()
    m = gamma |> Stream.take_while(fn ({_k,x})->x>0 end) |> Enum.count
    gammap = gamma |> Stream.take_while(fn ({_k,x})->x>0 end) |> Stream.map(fn {_,x}->x end) |> Stream.concat([0.0]) |> Enum.to_list
    gammap = gammap |> puiseaux
    cov = -gamma0 + 2.0*(gammap |> Enum.sum)

    if cov < 0, do: throw {:negative_covariance, cov, 2*m}
    {cov,2*m}
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
