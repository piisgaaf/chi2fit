defmodule Chi2fit.Distribution do

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

  @moduledoc """
  Provides various distributions.
  """

  @type distribution() :: ((...) :: number())
  @type cdf() :: ((number) :: number())

  defmodule UnsupportedDistributionError do
    defexception message: "Unsupported distribution function"
  end

  ###
  ### Standard distributions
  ###

  @doc """
  Uniform distribution.
  """
  @spec uniform(Keyword.t) :: distribution
  def uniform([]), do: uniform(0, 2.0)
  def uniform([avg: average]), do: uniform(0,2*average)
  def uniform(list) when is_list(list), do: fn () -> Enum.random(list) end

  @doc """
  Uniform distribution.
  """
  @spec uniform(min::integer(),max::integer()) :: distribution
  def uniform(min,max) when max>=min, do: fn () -> random(min,max) end

  @doc """
  Constant distribution.
  """
  @spec constant(number | Keyword.t) :: distribution
  def constant([avg: average]), do: fn () -> average end
  def constant(average) when is_number(average), do: fn () -> average end

  @doc """
  The exponential distribution.
  """
  @spec exponential(Keyword.t) :: distribution
  def exponential([avg: average]) do
    fn () ->
      u = :rand.uniform()
      -average*:math.log(u)
    end
  end
  def exponential([cdf: rate]), do: fn (t) -> 1.0 - :math.exp(-rate*t) end

  @doc """
  The Erlang distribution.
  """
  @spec erlang(mean::number(),m::pos_integer()) :: distribution
  def erlang(mean, m) when is_integer(m) and m>0 do
    list = 1..m
    fn () ->
      -(mean/m)*:math.log(list |> Enum.reduce(1.0, fn (_,acc) -> :rand.uniform()*acc end))
    end
  end

  @gamma53 0.902745292950933611297
  @gamma32 0.886226925452758013649

  @doc """
  The Weibull distribution.
  """
  @spec weibull(number, number|Keyword.t) :: distribution
  def weibull(1.0, [avg: average]), do: weibull(1.0, average)
  def weibull(1.5, [avg: average]), do: weibull(1.5, average/@gamma53)
  def weibull(2.0, [avg: average]), do: weibull(2.0, average/@gamma32)
  def weibull(alpha, beta) when is_number(alpha) and is_number(beta) do
    fn () ->
      u = :rand.uniform()
      beta*:math.pow(-:math.log(u),1.0/alpha)
    end
  end
  
  @doc """
  The Weibull cumulative distribution function.
  """
  @spec weibullCDF(number,number) :: cdf
  def weibullCDF(k,_) when k<0, do: raise ArithmeticError, "Weibull is only defined for positive shape"
  def weibullCDF(_,lambda) when lambda<0, do: raise ArithmeticError, "Weibull is only defined for positive scale"
  def weibullCDF(k,lambda) when is_number(k) and is_number(lambda) do
    fn
      0 -> 0.0
      0.0 -> 0.0
      x when x<0 -> 0.0
      x ->
        if :math.log(x/lambda)*k > 100 do
          0.0
        else
          1.0 - :math.exp -:math.pow(x/lambda,k)
        end
    end
  end

  @doc """
  The normal or Gauss distribution
  """
  @spec normal(mean::number(),sigma::number()) :: distribution()
  def normal(mean,sigma) when is_number(mean) and is_number(sigma) and sigma>=0 do
    fn () ->
      {w,v1,_} = polar()
      y = :math.sqrt(-2*:math.log(w)/w)
      mean + sigma*(v1*y)
    end
  end

  @doc """
  The Bernoulli distribution.
  """
  @spec bernoulli(value :: number) :: distribution
  def bernoulli(value) when is_number(value) do
   fn () ->
       u = :rand.uniform()
       if u <= value, do: 1, else: 0
   end
  end

  @doc """
  Wald or Inverse Gauss distribution.
  """
  @spec wald(mu::number(),lambda::number()) :: distribution
  def wald(mu,lambda) when is_number(mu) and is_number(lambda) do
   fn () ->
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

  @doc """
  The Wald cumulative distribution function.
  """
  @spec waldCDF(number,number) :: cdf
  def waldCDF(mu,_) when mu < 0, do: raise ArithmeticError, "Wald is only defined for positive average"
  def waldCDF(_,lambda) when lambda < 0, do: raise ArithmeticError, "Wald is only defined for positive shape"
  def waldCDF(mu,lambda) do
    fn
      x when x == 0 -> 0.0
      x when x < 0 -> 0.0
      x when x>0 ->
        phi(:math.sqrt(lambda/x) * (x/mu-1.0)) + :math.exp(2.0*lambda/mu) * phi(-:math.sqrt(lambda/x) * (x/mu+1.0))
    end
  end

  ###
  ### Special distributions
  ###

  @doc """
  Distribution for flipping coins.
  """
  @spec coin(integer) :: distribution
  def coin(value), do: uniform([0.0,value])

  @doc """
  Distribution simulating a dice (1..6)
  """
  @spec dice([] | number) :: distribution
  def dice([]), do: dice(1.0)
  def dice([avg: avg]), do: dice(avg)
  def dice(avg), do: uniform([avg*1,avg*2,avg*3,avg*4,avg*5,avg*6])

  @doc """
  Distribution simulating the dice in the GetKanban V4 simulation game.
  """
  @spec dice_gk4([] | number) :: distribution
  def dice_gk4([]), do: dice_gk4(1.0)
  def dice_gk4([avg: avg]), do: dice_gk4(avg)
  def dice_gk4(avg), do: uniform([avg*3,avg*4,avg*4,avg*5,avg*5,avg*6])

  @doc """
  Returns the model for a name.
  """
  @spec model(name::String.t) :: [fun: cdf, df: pos_integer()]
  def model(name) do
    case name do
      "wald" -> [
        fun: fn (x,[k,lambda]) -> waldCDF(k,lambda).(x) end,
        df: 2
      ]
      "weibull" -> [
        fun: fn (x,[k,lambda]) -> weibullCDF(k,lambda).(x) end,
        df: 2
      ]
      unknown ->
        raise UnsupportedDistributionError, message: "Unsupported cumulative distribution function '#{inspect unknown}'"
    end
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

  @spec polar() :: {number(), number(), number()}
  defp polar() do
    v1 = random(-1,1)
    v2 = random(-1,1)
    w = v1*v1 + v2*v2

   cond do
     w > 1.0 -> polar()
     true -> {w,v1,v2}
   end
 end

end
