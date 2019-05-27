defmodule Distribution.Uniform do

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
  Provides the Uniform distribution.
  """

  defstruct [:list,:pars,:range]
  
  @type t() :: %__MODULE__{
    list: [number()],
    pars: [number()],
    range: {number,number}
  }

end

defimpl Distribution, for: Distribution.Uniform do
  import Distribution.Uniform
  alias Distribution.Uniform

  @spec uniform(Keyword.t) :: ((...) -> number)
  defp uniform([]), do: uniform(0, 2.0)
  defp uniform([avg: average]), do: uniform(0,2*average)
  defp uniform(list) when is_list(list), do: fn () -> Enum.random(list) end

  @spec uniform(min::integer(),max::integer()) :: ((...) -> number)
  defp uniform(min,max) when max>=min, do: fn () -> random(min,max) end

  @spec random(min::number(),max::number()) :: number()
  defp random(min,max) when max >= min do
    min + (max-min)*:rand.uniform()
  end

  def skewness(%Uniform{}), do: raise Distribution.FunctionNotSupportedError, message: "skewness is not supported for the Uniform distribution"
  def kurtosis(%Uniform{}), do: raise Distribution.FunctionNotSupportedError, message: "kurtosis is not supported for the Uniform distribution"
  def size(%Uniform{}), do: 1

  def cdf(%Uniform{}), do: raise Distribution.FunctionNotSupportedError, message: "cdf is not supported for the Uniform distribution"
  def pdf(%Uniform{}), do: raise Distribution.FunctionNotSupportedError, message: "pdf is not supported for the Uniform distribution"

  @doc """
  ## Examples:

      iex> %Distribution.Uniform{range: {0,0}} |> random().()

  """
  def random(%Uniform{range: {min,max}}), do: uniform(min,max).()
  def random(%Uniform{list: list}), do: uniform(list).()
  def random(%Uniform{}), do: uniform([]).()
  
end
