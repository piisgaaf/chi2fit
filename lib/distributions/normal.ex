defmodule Distribution.Normal do

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
  The normal or Gauss distribution
  """

  defstruct [:pars, name: "normal"]
  
  @type t() :: %__MODULE__{
    pars: [number()] | nil,
    name: String.t
  }

end

defimpl Distribution, for: Distribution.Normal do
  import Distribution.Normal
  alias Distribution.Normal

  defp normal(mean,sigma) when is_number(mean) and is_number(sigma) and sigma>=0 do
    fn () ->
      {w,v1,_} = polar()
      y = v1*:math.sqrt(-2*:math.log(w)/w)
      mean + sigma*y
    end
  end

  defp normalCDF(_mean,sigma) when sigma<0, do: raise ArgumentError
  defp normalCDF(mean,sigma) when is_number(mean) and is_number(sigma) and sigma>=0 do
    fn
      x when (x-mean)/sigma < 4.0 ->
        0.5*:math.erfc(-(x-mean)/sigma/:math.sqrt(2.0)) 
      x ->
        0.5*( 1.0 + :math.erf((x-mean)/sigma/:math.sqrt(2.0)) )
    end
  end

  @spec polar() :: {number(), number(), number()}
  defp polar() do
    v1 = random(-1,1)
    v2 = random(-1,1)
    w = v1*v1 + v2*v2

    cond do
      w >= 1.0 -> polar()
      true -> {w,v1,v2}
    end
  end

  @spec random(min::number(),max::number()) :: number()
  defp random(min,max) when max >= min do
    min + (max-min)*:rand.uniform()
  end
  
  def skewness(%Normal{pars: nil}), do: fn _ -> 0 end
  def kurtosis(%Normal{pars: nil}), do: fn _ -> 0 end
  def size(%Normal{}), do: 2
  
  def cdf(%Normal{pars: nil}), do: fn x, [mu,sigma] -> normalCDF(mu,sigma).(x) end
  def cdf(%Normal{pars: [mu,sigma]}), do: normalCDF(mu,sigma)

  def pdf(%Normal{pars: nil}), do: fn x, [mu,sigma] -> 1/:math.sqrt(2*:math.pi)/sigma*:math.exp(-(x-mu)*(x-mu)/2/sigma/sigma) end

  def random(%Normal{pars: [mu,sigma]}), do: normal(mu,sigma).()
  def random(%Normal{pars: nil}), do: fn [mu,sigma] -> normal(mu,sigma).() end

  def name(model), do: model.name
  
end

defimpl Inspect, for: Distribution.Normal do
  import Inspect.Algebra
  
  def inspect(dict, opts) do
    case dict.pars do
      nil ->
        "#Normal<>"
      [mu,sigma] ->
        concat ["#Normal<", to_doc("mu=#{mu}, sigma=#{sigma}", opts), ">"]
      list ->
        concat ["#Normal<", to_doc(list, opts), ">"]
    end
  end

end