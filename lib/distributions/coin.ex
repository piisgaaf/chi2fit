defmodule Distribution.Coin do

  @moduledoc """
  Distribution for flipping coins.
  """

  defstruct []
  
  @type t() :: %__MODULE__{
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
end
