defmodule Chi2fit.Distribution.Nakagami do

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
  The Nakagami distribution.
  """

  defstruct [:pars, name: "nakagami"]
  
  @type t() :: %__MODULE__{
    pars: [number()] | nil,
    name: String.t
  }

end

defimpl Chi2fit.Distribution, for: Chi2fit.Distribution.Nakagami do
  alias Chi2fit.Distribution, as: D

  import D.Nakagami
  alias D.Nakagami

  import Exboost.Math, only: [tgamma: 1]

  @spec nakagami(scale::number(),shape::number()) :: ((...) -> number)
  defp nakagami(scale,shape) do
    fn ->
      u = :rand.uniform()
      scale*:math.sqrt(Exboost.Math.gamma_p_inv(shape,u)/shape)
    end
  end

  @spec nakagamiCDF(scale :: float,shape :: float) :: (number -> number)
  defp nakagamiCDF(scale,shape) when scale>0 and shape>0 do
    fn
      x ->
        Exboost.Math.tgamma_lower(shape,shape*(x/scale)*(x/scale))
    end
  end
  defp nakagamiCDF(_scale,_shape), do: raise ArithmeticError, "Nakagami is only defined for positive scale and shape"

  def skewness(%Nakagami{pars: nil}) do
    fn [_scale,shape] ->
      g = tgamma(shape)
      g1_2 = tgamma(shape+0.5)
      g1 = tgamma(shape+1.0)
      g3_2 = tgamma(shape+1.5)
      num = 2*g1_2*g1_2*g1_2 + g*g*( g3_2 - 3*shape*g1_2 )
      den = g*g*g*:math.pow(shape*(1.0-shape*g1_2*g1_2/g1/g1),1.5)
      num/den
    end
  end
  def kurtosis(%Nakagami{pars: nil}) do
    fn [_scale,shape] ->
      g = tgamma(shape)
      g1_2 = tgamma(shape+0.5)
      g2 = tgamma(shape+2.0)
      gdouble = tgamma(2*shape)
      num = -6*g1_2*g1_2*g1_2*g1_2 - 3*shape*shape*g*g*g*g + g*g*g*g2 + :math.pow(2,3-4*shape)*(4*shape-1)*:math.pi*gdouble*gdouble
      den = :math.pow(abs(g1_2*g1_2 - shape*g*g),2)
      num/den
    end
  end
  def size(%Nakagami{}), do: 2
  
  def cdf(%Nakagami{pars: nil}), do: fn x, [scale,shape] -> nakagamiCDF(scale,shape).(x) end

  def pdf(%Nakagami{pars: nil}), do: raise D.FunctionNotSupportedError, message: "pdf is not supported for the Nakagami distribution"
  def random(%Nakagami{pars: [scale,shape]}), do: nakagami(scale, shape).()

  def name(model), do: model.name
  
end

defimpl Inspect, for: Chi2fit.Distribution.Nakagami do
  import Inspect.Algebra
  
  def inspect(dict, opts) do
    case dict.pars do
      nil ->
        "#Nakagami<>"
      [scale,shape] ->
        concat ["#Nakagami<", to_doc("scale=#{scale}, shape=#{shape}", opts), ">"]
      list ->
        concat ["#Nakagami<", to_doc(list, opts), ">"]
    end
  end

end