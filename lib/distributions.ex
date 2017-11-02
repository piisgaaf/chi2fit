defmodule Chi2fit.Distributions do

  # Copyright 2012-2017 Pieter Rijken
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

  @type distribution() :: ((...) :: number())

  ## Various distributions
  @spec uniform(min::integer(),max::integer()) :: distribution
  def uniform([]), do: uniform(0, 2.0)
  def uniform([avg: average]), do: uniform(0,2*average)
  def uniform(list) when is_list(list), do: fn () -> Enum.random(list) end
  def uniform(min,max) when max>=min, do: fn () -> random(min,max) end

  @spec constant(number | Keyword.t) :: distribution
  def constant([avg: average]), do: fn () -> average end
  def constant(average) when is_number(average), do: fn () -> average end

  @spec coin(integer) :: distribution
  def coin(value), do: uniform([0.0,value])

  @spec dice([] | number) :: distribution
  def dice([]), do: dice(1.0)
  def dice([avg: avg]), do: dice(avg)
  def dice(avg), do: uniform([avg*1,avg*2,avg*3,avg*4,avg*5,avg*6])

  @spec dice_gk4([] | number) :: distribution
  def dice_gk4([]), do: dice_gk4(1.0)
  def dice_gk4([avg: avg]), do: dice_gk4(avg)
  def dice_gk4(avg), do: uniform([avg*3,avg*4,avg*4,avg*5,avg*5,avg*6])

  @spec exponential(Keyword.t) :: distribution
  def exponential([avg: average]) do
    fn () ->
      u = :rand.uniform()
      -average*:math.log(u)
    end
  end
  def exponential([cdf: rate]), do: fn (t) -> 1.0 - :math.exp(-rate*t) end

# -spec erlang(Mean::number(),M::pos_integer()) -> distribution().
# erlang(Mean, M) when is_integer(M) andalso M>0 ->
#     List = lists:seq(1,M),
#     fun
#     () ->
#       U = random:uniform(),
#       -(Mean/M)*math:log(lists:foldl(fun (_E,Acc) -> U*Acc end, 1, List))
#   end.

  @gamma53 0.902745292950933611297
  @gamma32 0.886226925452758013649

  def weibull(1.0, [avg: average]), do: weibull(1.0, average)
  def weibull(1.5, [avg: average]), do: weibull(1.5, average/@gamma53)
  def weibull(2.0, [avg: average]), do: weibull(2.0, average/@gamma32)
  def weibull(alpha, beta) when is_number(alpha) and is_number(beta) do
    fn () ->
      u = :rand.uniform()
      beta*:math.pow(-:math.log(u),1.0/alpha)
    end
  end
  
  def weibullCDF(0,_,_), do: 0.0
  def weibullCDF(0.0,_,_), do: 0.0
  def weibullCDF(x,_,_) when x<0, do: 0.0
  def weibullCDF(_,k,_) when k<0, do: 0.0
  def weibullCDF(_,_,lambda) when lambda<0, do: 0.0
  def weibullCDF(x,k,lambda) when is_number(x) and is_number(k) and is_number(lambda) do
    require Logger
    try do
      if :math.log(x/lambda)*k > 100, do: 0.0, else: 1.0 - :math.exp -:math.pow(x/lambda,k)
    rescue
      e ->
        stack=System.stacktrace
        Logger.error "args=#{x},#{k},#{lambda}"
        Logger.error "ERROR: #{inspect e} #{inspect stack}"
        raise e
    end
  end

  # @spec normal(mean::number(),sigma::number()) :: distribution()
  # def normal(mean,sigma) when is_number(mean) and is_number(sigma) and sigma>=0 do
  #   fn () ->
  #     {w,v1,_} = polar()
  #     y = :math.sqrt(-2*:math.log(w)/w)
  #     mean + sigma*(v1*y)
  #   end
  # end

  defp bernoulli(value) when is_number(value) do
   fn () ->
       u = :rand.uniform()
       if u <= value, do: 1, else: 0
   end
  end

  @spec wald(mu::number(),lambda::number()) :: distribution
  def wald(mu,lambda) when is_number(mu) and is_number(lambda) do
   fn
  ##     (:average) -> mu
  ##     (:stddev) -> :math.sqrt(mu*mu*mu/lambda)
     () ->
       w = :rand.uniform()
       y = w*w
       z = mu + mu*mu*y/2/lambda + mu/2/lambda*:math.sqrt(4*mu*lambda*y+mu*mu*y*y)
       case (bernoulli(mu/(mu+z))).() do
         1 -> z
         _else -> mu*mu/z
       end
   end
  end
  def wald([avg: average],lambda), do: wald(average,lambda)

  def waldCDF(x,_,_) when x == 0, do: 0.0
  def waldCDF(x,_,_) when x < 0, do: 0.0
  def waldCDF(_,mu,_) when mu < 0, do: raise ArithmeticError, "Wald is only defined for positive average"
  def waldCDF(_,_,lambda) when lambda < 0, do: raise ArithmeticError, "Wald is only defined for positive shape"
  def waldCDF(x,mu,lambda) when x>0 and lambda>=0 do
   phi(:math.sqrt(lambda/x) * (x/mu-1.0)) + :math.exp(2.0*lambda/mu) * phi(-:math.sqrt(lambda/x) * (x/mu+1.0))
  end

  def poissonCDF(x,_) when x == 0, do: 0.0
  def poissonCDF(x,_) when x < 0, do: 0.0
  def poissonCDF(x,lambda) when is_float(x), do: poissonCDF Float.ceil(x),lambda
  def poissonCDF(x,lambda) when x>0 and is_integer(x) do
   :math.exp(-lambda)*(0..x-1 |> Enum.reduce({1.0,0.0},
     fn
       (0,{_,_})->{1.0,1.0}
       (k,{acc,sum})->
         delta=acc*lambda/k
         {delta,sum+delta}
     end) |> elem(1))
  end

  ##
  ## Local Functions
  ##
  @spec random(min::number(),max::number()) :: number()
  defp random(min,max) when max >= min do
    min + (max-min)*:rand.uniform()
  end
  
  @spec phi(x :: float) :: float
  defp phi(x) do
   (1.0 + :math.erf(x/:math.sqrt(2.0)))/2.0
  end

  defmodule UnsupportedDistributionError do
    defexception message: "Unsupported distribution function"
  end

  def model(name,ranges) do
    result = case name do
      "wald" -> [
        fun: fn (x,[mu,lambda]) -> 1.0-waldCDF(x,mu,lambda) end,
        curve: fn ([k,lambda]) -> fn x->waldCDF(x,k,lambda) end end,
        df: 2,
        init: [65.0,1.0],
        probe: [{10.0,80.0},{0.1,20.0}]
      ]
      "weibull" -> [
        fun: fn (x,[k,lambda]) -> 1.0-weibullCDF(x,k,lambda) end,
        curve: fn ([k,lambda]) -> fn x->weibullCDF(x,k,lambda) end end,
        df: 2,
        init: [1.0,1.0],
        probe: [{0.55,0.65},{26.0,27.0}]
      ]
      "cpoisson" -> [
        fun: fn (x,[lambda]) -> 1.0-poissonCDF(x,lambda) end,
        curve: fn ([lambda]) -> fn x->poissonCDF(x,lambda) end end,
        df: 1,
        init: [1.0],
        probe: [{0.01,9.9}]
      ]
      unknown ->
        raise UnsupportedDistributionError, message: "Unsupported cumulative distribution function '#{inspect unknown}'"
    end
    if ranges, do: Keyword.put(result,:probe,ranges), else: result
  end

  # @spec polar() :: {number(), number(), number()}
  # defp polar() do
  #   v1 = 2*:random.uniform()-1
  #   v2 = 2*:random.uniform()-1
  #   w = v1*v1 + v2*v2
  #
  #   cond do
  #     w > 1.0 -> polar()
  #     true -> {w,v1,v2}
  #   end
  # end

end
