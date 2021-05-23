defmodule Distribution.Coin do

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
  Distribution for flipping coins.
  """

  defstruct [name: "coin"]
  
  @type t() :: %__MODULE__{
    name: String.t
  }

end

defimpl Distribution, for: Distribution.Coin do
  import Distribution.Coin
  alias Distribution.Coin

  @spec coin() :: ((...) -> :heads|:tails)
  defp coin() do
    fn -> if(:rand.uniform()<0.5, do: :heads, else: :tails) end
  end
  
  def skewness(%Coin{}), do: fn _ -> 0 end
  def kurtosis(%Coin{}), do: fn _ -> -2 end
  def size(%Coin{}), do: 0
  
  def cdf(%Coin{}) do
    fn 
      x, _ when x < 0.0 -> 0.0
      x, _ when 0.0 <= x and x < 1.0 -> 0.5
      _, _ -> 1.0
    end
  end

  def pdf(%Coin{}), do: fn _, _ -> 0.5 end
  def random(%Coin{}), do: coin().()

  def name(model), do: model.name
  
end

defimpl Inspect, for: Distribution.Coin do
  def inspect(_dict, _opts), do: "#Coin<>"

end