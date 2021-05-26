defmodule Chi2fit.Distribution.Frechet do

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
  The Fréchet distribution, also known inverse Weibull distribution.
  """

  defstruct [:pars, name: "frechet"]
  
  @type t() :: %__MODULE__{
    pars: [number()] | nil,
    name: String.t
  }

end

defimpl Chi2fit.Distribution, for: Chi2fit.Distribution.Frechet do
  alias Chi2fit.Distribution, as: D

  import D.Frechet
  alias D.Frechet

  import Exboost.Math, only: [tgamma: 1]

  @spec frechet(scale::number(),shape::number()) :: ((...) -> number)
  defp frechet(scale,shape) when is_number(scale) and is_number(shape) do
   fn ->
     u = :rand.uniform()
     scale * :math.pow(-:math.log(u),-1.0/shape)
   end
  end

  @spec frechetCDF(scale :: float,shape :: float) :: (number -> number)
  defp frechetCDF(scale,shape) when scale>0 and shape>0 do
    fn
      x when x==0.0 ->
        0.0
      x ->
        :math.exp(-:math.pow(x/scale,-shape))
    end
  end
  defp frechetCDF(_scale,_shape), do: raise ArithmeticError, "Fréchet is only defined for positive scale and shape"

  def skewness(%Frechet{pars: nil}) do
    fn [_scale,shape] ->
      g1 = tgamma(1.0-1.0/shape)
      g2 = tgamma(1.0-2.0/shape)
      g3 = tgamma(1.0-3.0/shape)
      (g3 - 3*g2*g1 + 2*g1*g1*g1)/:math.pow(g2 - g1*g1,1.5)
    end
  end
  def kurtosis(%Frechet{pars: nil}) do
    fn [_scale,shape] ->
      g1 = tgamma(1.0-1.0/shape)
      g2 = tgamma(1.0-2.0/shape)
      g3 = tgamma(1.0-3.0/shape)
      g4 = tgamma(1.0-4.0/shape)
      -6 + (g4 - 4*g3*g1 + 3*g2*g2)/:math.pow(g2 - g1*g1,2.0)
    end
  end
  def size(%Frechet{}), do: 2
  
  def cdf(%Frechet{pars: nil}), do: fn x, [scale,shape] -> frechetCDF(scale,shape).(x) end

  def pdf(%Frechet{pars: nil}), do: raise D.FunctionNotSupportedError, message: "pdf is not supported for the Frechet distribution"
  def random(%Frechet{pars: [scale,shape]}), do: frechet(scale, shape).()

  def name(model), do: model.name
  
end

defimpl Inspect, for: Chi2fit.Distribution.Frechet do
  import Inspect.Algebra
  
  def inspect(dict, opts) do
    case dict.pars do
      nil ->
        "#Frechet<>"
      [scale,shape] ->
        concat ["#Frechet<", to_doc("scale=#{scale}, shape=#{shape}", opts), ">"]
      list ->
        concat ["#Frechet<", to_doc(list, opts), ">"]
    end
  end

end