defmodule Chi2fit.Utilities do

  @h 1.0e-10

  @type cdf :: ((number)->{number,number})

  require Integer
  require Logger
  import Kernel, except: [*: 2, /: 2,+: 2, -: 2]

  @type real :: number
  @type complex :: {real,real}
  
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

  @spec make_histogram([number],number) :: %{required(number) => pos_integer}
  def make_histogram(list,binsize \\ 1) do
    Enum.reduce(list, %{}, fn
      (number,acc) ->
        acc |> Map.update(trunc(number/binsize),1,&(1+&1))
    end) |> Enum.reduce([], fn (pair,acc)->[pair|acc] end) |> Enum.sort_by(fn ({k,_v})->k end)
  end

  @spec to_cdf([number],number,number) :: [{float,float}]
  def to_cdf(list, bin, interval \\ 0), do: to_cdf(list, bin, interval, 0.0, [])

  defp to_cdf([], _bin, _interval, _sum, result), do: Enum.reverse(result)
  defp to_cdf(list, bin, interval, sum, result) do
    {in_bin, out_bin} = list |> Enum.partition(fn (x)->x<=interval+bin end)
    to_cdf(out_bin, bin, interval+bin, sum+length(in_bin), [[interval+bin,sum+length(in_bin)]|result])
  end

  defmodule UnknownSampleErrorAlgorithmError do
    defexception message: "unknown sample error algorithm"
  end

  @doc """
  See https://en.wikipedia.org/wiki/Cumulative_frequency_analysis
  And: https://arxiv.org/pdf/1112.2593v3.pdf
  See https://en.wikipedia.org/wiki/Student%27s_t-distribution:
    90% confidence ==> t = 1.645 for many data points (> 120)
    70% confidence ==> t = 1.000
  """
  @type algorithm :: :wilson | :wald
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
          ## Wilson score:
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
  See section 8.4 in "Handbook of Monte Carlo Methods" by Kroese, Taimre, and Botev

  Three parameters determine the resulting empirical distribution:
  1) algorithm for assigning errors,
  2) the size of the bins,
  3) a correction for limiting the bounds on the 'y' values
  
  When e.g. task effort/duration is modeled, some tasks measured have 0 time. In practice
  what is actually is meant, is that the task effort is between 0 and 1 hour. This is where
  binning of the data happens. Specify a size of the bins to control how this is done. A bin
  size of 1 means that 0 effort will be mapped to 1/2 effort (at the middle of the bin).
  This also prevents problems when the fited distribution cannot cope with an effort os zero.
  
  In the handbook of MCMC a cumulative distribution is constructed. For the largest 'x' value
  in the sample, the 'y' value is exactly one (1). In combination with the Wald score this
  gives zero errors on the value '1'. If the resulting distribution is used to fit a curve
  this may give an infinite contribution to the maximum likelihood function.
  Use the correction number to have a 'y' value of slightly less than 1 to prevent this from
  happening.
  Especially the combination of 0 correction, algorithm ':wald', and 'linear' model for
  handling asymmetric errors gives problems.
  
  The algorithm parameter determines how the errors onthe 'y' value are determined. Currently
  supported values include ':wald' and 'wilson'.
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

  @spec get_cdf([number], number) :: {cdf,[float],pos_integer,float,[number]}
  def get_cdf(data, binsize \\ 1,algorithm \\ :wilson) do
    data
    |> make_histogram(binsize)
    |> empirical_cdf(binsize,algorithm)
  end

  def convert_cdf({cdf,[mindur,maxdur]}) do
    round(mindur)..round(maxdur)
    |> Stream.map(fn (x)->
        {y,y1,y2} = cdf.(x)
        {x,y,y1,y2}
      end)
    |> Stream.map(fn ({x,y,y1,y2})->{x,1.0-y,1.0-y2,1.0-y1} end)
    |> Enum.to_list
  end

  ## Solve polynomial equations
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
    import :math
    c = -0.5*q*pow(3/abs(p),1.5)
    cond do
      p>0 ->
        [sinh(1.0/3.0*asinh(c))]
      c>=1 ->
        [cosh(1.0/3.0*acosh(c))]
      c<=-1 ->
        [-cosh(1.0/3.0*acosh(abs(c)))]
      true ->
        ## Three real solutions
        [cos(1.0/3.0*acos(c)),cos(1.0/3.0*acos(c) + 2*pi()/3.0),cos(1.0/3.0*acos(c) + 4*pi()/3.0)]
    end
    |> Enum.map(&(&1*2*sqrt(abs(p)/3.0)))
  end
  def solve([1.0,a2,a1,a0]), do: solve([1.0,0.0,(3*a1-a2*a2)/3.0,(2*a2*a2*a2-9*a1*a2+27*a0)/27.0]) |> Enum.map(&(&1-a2/3.0))
  def solve([a3,a2,a1,a0]), do: solve([1.0,a2/a3,a1/a3,a0/a3])

  defmacro mapder(list, delta \\ 0.0) do
    quote do
      unquote(list) |> Enum.flat_map(fn
        ({x,1}) when is_number(x) -> (x/1.0 + unquote(delta))
        (x) when is_number(x) -> x/1.0
      end)
    end
  end

  defp expand_pars(list) do
    list |> Enum.map(
          fn
            ({x,0}) when is_number(x) -> x/1.0
            ({x,n}) when is_number(x) -> List.flatten expand_pars([{x/1.0 + @h,n-1},{{x/1.0,n-1}}])
            (x) when is_number(x) -> x/1.0
            ({{x,0}}) when is_number(x) -> {x/1.0}
            ({{x,n}}) when is_number(x) -> List.flatten expand_pars([{{x/1.0 + @h,n-1}},{x/1.0,n-1}])
          end)
  end

  defp reduce_pars(list) do
    list |> Enum.reduce([{[],1}],
      fn
        (x,acc) when is_number(x) -> Enum.map(acc, fn ({y,n})->{[x|y],n} end)
        (list,acc) when is_list(list) ->
          Enum.flat_map(list,
            fn
              (x) when is_number(x) -> Enum.map(acc, fn ({y,n})->{[x|y],n} end)
              ({x}) when is_number(x) -> Enum.map(acc, fn ({y,n})->{[x|y],-n} end)
            end)
      end)
      |> Enum.map(fn ({l,n}) -> {Enum.reverse(l),n} end)
  end

  @spec der([float|{float,integer}], (([float])->float)) :: float
  def der(parameters, fun, debug \\ false) do
    factor = Enum.reduce(parameters,1.0,
      fn
        (x,acc) when is_number(x) -> acc
        ({x,n},acc) when is_number(x) -> acc*:math.pow(@h,n)
      end)
    d = parameters |> expand_pars |> reduce_pars
    if debug, do: Logger.debug "===> #{inspect d}"
    d |> Enum.reduce(0.0, fn ({x,n},sum) when is_list(x) -> sum + n*fun.(x) end) |> Kernel./(factor)
  end

  defp jacobian(x=[_|_], k, fun) when k>0 and k<=length(x) and is_function(fun,1) do
    x |> List.update_at(k-1, fn (val) -> {val,1} end) |> der(fun)
  end

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

  ## Zie: https://en.wikipedia.org/wiki/Cooley%E2%80%93Tukey_FFT_algorithm
  # Parallel version of FFT; see http://www.webabode.com/articles/Parallel%20FFT%20implementations.pdf
  @spec fft([real],Keyword.t) :: [complex]
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
  
  @spec ifft([real],Keyword.t) :: [complex]
  def ifft(list,opts \\ [nproc: 1]) do
    n = length(list)
    list |> fft(Keyword.merge(opts,[phase: -1])) |> Enum.map(&(&1/n))
  end

  @spec normv([complex]|complex) :: real
  def normv({x,y}), do: x*x+y*y
  def normv(list) when is_list(list), do: list |> Enum.map(&normv/1)

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

  @spec auto([real],Keyword.t) :: [real]
  def auto(list,opts \\ [nproc: 1])
  
  def auto([],_opts), do: []
  def auto([x],_opts), do: [x*x]
  def auto(list,opts) do
    n = length(list)
    List.duplicate(0,n) |> Enum.concat(list) |> fft(opts) |> normv |> ifft(opts) |> Stream.take(n) |> Stream.map(&(elem(&1,0))) |> Enum.to_list
  end

  ##
  ## See section 1.10.2 (Initial sequence method) of 'Handbook of Markov Chain Monte Carlo'
  ## Input is a list of gamma_k
  ##
  @spec error([{gamma :: number,k :: pos_integer}], :initial_sequence_method) :: {number, number}
  def error(nauto, :initial_sequence_method) do
    ## For reversible Markov Chains
    gamma = nauto |> Stream.chunk(2) |> Stream.map(fn ([{x,k},{y,_}])->{k/2,x+y} end) |> Enum.to_list
    gamma0 = nauto |> Stream.take(1) |> Enum.to_list |> (&(elem(hd(&1),0))).()
    m = gamma |> Stream.take_while(fn ({_k,x})->x>0 end) |> Enum.count
    gammap = gamma |> Stream.take_while(fn ({_k,x})->x>0 end) |> Stream.map(fn {_,x}->x end) |> Stream.concat([0.0]) |> Enum.to_list
    gammap = gammap |> puiseaux
    cov = -gamma0 + 2.0*(gammap |> Enum.sum)

    if cov < 0, do: Logger.debug "WARNING: cov<0 [nauto=#{length nauto}::#{inspect nauto}]"
    {cov,2*m}
  end

  ##
  ## Bootstrapping
  ##
  def bootstrap(total, data, fun, options) do
    debug? = options |> Keyword.get(:debug, false)
    safe = options |> Keyword.get(:safe, false)

    {start,continuation} = case safe do
      :safe ->
        file = options |> Keyword.fetch!(:filename)
        {:ok,:storage} = :dets.open_file :storage, type: :set, file: file, auto_save: 1000, estimated_no_objects: total
        :ok = :dets.delete_all_objects :storage
        {1,[]}        
      :cont ->
        file = options |> Keyword.fetch!(:filename)
        if debug?, do: Logger.debug "Reading saved data from previous run..."
        {:ok,:storage} = :dets.open_file :storage, type: :set, file: file, auto_save: 1000, estimated_no_objects: total
        if debug?, do: Logger.debug "#{inspect :dets.info :storage}"
        objects = :dets.select(:storage, [{{:'_',:'$1'},[],[:'$1']}])
        {length(objects)+1,objects}
      _ ->
        {1,[]}
    end

    if start>total, do: raise ArgumentError, message: "start cannot be larger than the total"

    1..total |> Enum.reduce(continuation, fn (k,acc) ->
      try do
        ## Run Monte Carlo
        result = data |> Enum.map(fn _ -> Enum.random(data) end) |> fun.(k)

        if safe, do: true = :dets.insert_new :storage, {k,result}

        [result|acc]
      rescue
        _error ->
          stack = System.stacktrace
          Logger.debug "#{inspect stack}"
          [nil|acc]
      end
    end)
  end
  
  def read_data(filename) do
    filename
    |> File.stream!([],:line)
    |> Stream.flat_map(&String.split(&1,"\r",trim: true))
    |> Stream.filter(&is_tuple(Float.parse(&1)))
    |> Stream.map(&elem(Float.parse(&1),0))
    |> Stream.filter(&(&1 >= 0.0))
  end

end
