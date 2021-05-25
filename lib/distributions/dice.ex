defmodule Chi2fit.Distribution.Dice do

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
  Provides the Dice distribution.
  """

  defstruct [mode: :regular, name: "dice"]

  @type t() :: %__MODULE__{
    mode: :regular | :gk4,
    name: String.t
  }

end

defimpl Chi2fit.Distribution, for: Chi2fit.Distribution.Dice do
  alias Chi2fit.Distribution, as: D

  import Chi2fit.Distribution.Dice
  alias Chi2fit.Distribution.Dice

  @spec dice([] | number) :: ((...) -> number)
  defp dice([]), do: dice(1)
  defp dice([avg: avg]), do: dice(avg)
  defp dice(avg), do: %D.Uniform{pars: [avg*1,avg*2,avg*3,avg*4,avg*5,avg*6]}

  @spec dice_gk4([] | number) :: ((...) -> number)
  defp dice_gk4([]), do: dice_gk4(1)
  defp dice_gk4([avg: avg]), do: dice_gk4(avg)
  defp dice_gk4(avg), do: %D.Uniform{pars: [avg*3,avg*4,avg*4,avg*5,avg*5,avg*6]}

  def skewness(%Dice{}), do: raise D.FunctionNotSupportedError, message: "skewness is not supported for the Dice distribution"
  def kurtosis(%Dice{}), do: raise D.FunctionNotSupportedError, message: "kurtosis is not supported for the Dice distribution"
  def size(%Dice{}), do: 0

  def cdf(%Dice{}), do: raise D.FunctionNotSupportedError, message: "cdf is not supported for the Dice distribution"
  def pdf(%Dice{}), do: raise D.FunctionNotSupportedError, message: "pdf is not supported for the Dice distribution"

  def random(%Dice{mode: :regular}), do: D.random(dice([]))
  def random(%Dice{mode: :gk4}), do: D.random(dice_gk4([]))
  
  def name(model), do: model.name
  
end

defimpl Inspect, for: Chi2fit.Distribution.Dice do
  import Inspect.Algebra
  
  def inspect(dict, opts) do
    case dict.mode do
      nil ->
        "#Dice<>"
      mode ->
        concat ["#Dice<", to_doc("mode=#{mode}", opts), ">"]
    end
  end

end