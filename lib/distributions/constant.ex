defmodule Distribution.Constant do

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
  defstruct [:pars]
  
  @type t() :: %__MODULE__{
  }
  
end

defimpl Distribution, for: Distribution.Constant do
  import Distribution.Constant
  alias Distribution.Constant

  @spec constant(number | Keyword.t) :: ((...) -> number)
  defp constant([avg: average]), do: fn () -> average end
  defp constant(average) when is_number(average), do: fn () -> average end
  
  def skewness(%Constant{}), do: fn _ -> 0 end
  def kurtosis(%Constant{}), do: fn _ -> 0 end
  def size(%Constant{}), do: 1
  
  def cdf(%Constant{}), do: raise Distribution.FunctionNotSupportedError, message: "cdf is not supported for the Constant distribution"
  def pdf(%Constant{}), do: raise Distribution.FunctionNotSupportedError, message: "pdf is not supported for the Constant distribution"

  @doc """
  ## Examples:
  
      iex> %Distribution.Constant{pars: [5.0]} |> Distribute.random()

  """
  def random(%Constant{pars: [avg]}), do: constant(avg).()
  
end
