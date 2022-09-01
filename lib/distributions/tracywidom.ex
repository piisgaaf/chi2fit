defmodule Chi2fit.Distribution.TracyWidom do

  # Copyright 2019 Pieter Rijken
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
  Tracy-Widom distribution.
  """

  defstruct [:pars, :type, name: "tracy-widom"]
  
  @type t() :: %__MODULE__{
    pars: [number()] | nil,
    type: 1|2|4,
    name: String.t
  }

end


defimpl Chi2fit.Distribution, for: Chi2fit.Distribution.TracyWidom do
  alias Chi2fit.Distribution, as: D

  import D.TracyWidom
  alias D.TracyWidom

  import Exboost.Math, only: [tgamma: 1, tgamma_lower: 2]

  @t1k 46.446
  @t1theta 0.186054
  @t1alpha 9.84801

  @t2k 79.6595
  @t2theta 0.101037
  @t2alpha 9.81961

  @t4k 146.021
  @t4theta 0.0595445
  @t4alpha 11.0016

  @e :math.exp(1.0)
  
  defp _gamma(k) when k > 1.0 do
    delta = k - :math.floor(k)
    
    u = :rand.uniform()
    v = :rand.uniform()
    w = :rand.uniform()
    
    if u <= @e/(@e+delta) do
      x = :math.pow(v, 1/delta)
      if w > :math.exp(-x), do: _gamma(k), else: x
    else
      x = 1.0 - :math.log(v)
      if w > :math.pow(x, delta-1), do: _gamma(k), else: x
    end
  end
  
  defp gamma(k,1.0) when k > 1.0 do
    _gamma(k) - (1..trunc(k) |> Enum.map(fn _ -> :math.log(:rand.uniform()) end) |> Enum.sum)
  end
  
  @spec tracywidom(number, number, number) :: ((...) -> number)
  defp tracywidom(k, theta, alpha) do
    fn -> theta*gamma(k,1.0) - alpha end
  end
  
  @spec tracywidomCDF(number,number,number) :: (number -> number)
  defp tracywidomCDF(k,theta,alpha) do
    fn
      x when x == alpha -> 0.0
      x when x < -alpha -> 0.0
      x when x > -alpha ->
        1/tgamma(k)*tgamma_lower(k, (x + alpha)/theta)
    end
  end

  @spec tracywidomPDF(number,number,number) :: (number -> number)
  defp tracywidomPDF(k,theta,alpha) do
    fn
      x when x == alpha -> 0.0
      x when x < -alpha -> 0.0
      x when x > -alpha ->
        1/tgamma(k)*:math.pow(theta,-k)*:math.pow(x+alpha,k-1)*:math.exp(-(x+alpha)/theta)
    end
  end
  
  def skewness(%TracyWidom{type: 1}), do: fn -> 2.0/:math.sqrt(@t1k) end
  def skewness(%TracyWidom{type: 2}), do: fn -> 2.0/:math.sqrt(@t2k) end
  def skewness(%TracyWidom{type: 4}), do: fn -> 2.0/:math.sqrt(@t4k) end

  def kurtosis(%TracyWidom{type: 1}), do: fn -> 0.1652429384 end
  def kurtosis(%TracyWidom{type: 2}), do: fn -> 0.0934480876 end
  def kurtosis(%TracyWidom{type: 4}), do: fn -> 0.0491951565 end

  def size(%TracyWidom{}), do: 0

  def cdf(%TracyWidom{type: 1}), do: fn x,[] -> tracywidomCDF(@t1k,@t1theta,@t1alpha).(x) end
  def cdf(%TracyWidom{type: 2}), do: fn x,[] -> tracywidomCDF(@t2k,@t2theta,@t2alpha).(x) end
  def cdf(%TracyWidom{type: 4}), do: fn x,[] -> tracywidomCDF(@t4k,@t4theta,@t4alpha).(x) end

  def pdf(%TracyWidom{type: 1}), do: fn x,[] -> tracywidomPDF(@t1k,@t1theta,@t1alpha).(x) end
  def pdf(%TracyWidom{type: 2}), do: fn x,[] -> tracywidomPDF(@t2k,@t2theta,@t2alpha).(x) end
  def pdf(%TracyWidom{type: 4}), do: fn x,[] -> tracywidomPDF(@t4k,@t4theta,@t4alpha).(x) end

  def random(%TracyWidom{type: 1}), do: fn -> tracywidom(@t1k,@t1theta,@t1alpha).() end
  def random(%TracyWidom{type: 2}), do: fn -> tracywidom(@t2k,@t2theta,@t2alpha).() end
  def random(%TracyWidom{type: 4}), do: fn -> tracywidom(@t4k,@t4theta,@t4alpha).() end

  def name(model), do: model.name
  
end

defimpl Inspect, for: Chi2fit.Distribution.TracyWidom do
  def inspect(dict, _opts) do
    type = cond do
      is_integer(dict.type) -> "_#{dict.type}"
      true -> ""
    end
        
    "#TracyWidom#{type}<>"
  end

end