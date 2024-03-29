defmodule Chi2fit.Distribution.Constant do

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
  Distribution for constant numbers.
  """

  @enforce_keys [:pars]
  defstruct [:pars, name: "constant"]
  
  @type t() :: %__MODULE__{
    name: String.t
  }
  
end

defimpl Chi2fit.Distribution, for: Chi2fit.Distribution.Constant do
  alias Chi2fit.Distribution, as: D

  import Chi2fit.Distribution.Constant
  alias Chi2fit.Distribution.Constant

  @spec constant(number | Keyword.t) :: ((...) -> number)
  defp constant([avg: average]), do: fn -> average end
  defp constant(average) when is_number(average), do: fn -> average end
  
  def skewness(%Constant{}), do: fn _ -> 0 end
  def kurtosis(%Constant{}), do: fn _ -> 0 end
  def size(%Constant{}), do: 1
  
  def cdf(%Constant{}), do: raise D.FunctionNotSupportedError, message: "cdf is not supported for the Constant distribution"
  def pdf(%Constant{}), do: raise D.FunctionNotSupportedError, message: "pdf is not supported for the Constant distribution"

  @doc ~S"""
  ## Examples:
  
      iex> random %Distribution.Constant{pars: [5.0]}
      5.0
  """
  def random(%Constant{pars: [avg]}), do: constant(avg).()
  
  def name(model), do: model.name
  
end

defimpl Inspect, for: Chi2fit.Distribution.Constant do
  import Inspect.Algebra
  
  def inspect(dict, opts) do
    case dict.pars do
      nil ->
        "#Constant<>"
      [value] ->
        concat ["#Constant<", to_doc("value=#{value}", opts), ">"]
      pars ->
        concat ["#Constant<", to_doc("pars=#{inspect pars}", opts), ">"]
    end
  end

end