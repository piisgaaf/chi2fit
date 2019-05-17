defmodule Distribution.Dice do

  @moduledoc """
  Provides the Dice distribution.
  """

  defstruct [mode: :regular]
  
  @type t() :: %__MODULE__{
    mode: :regular | :gk4
  }

end

defimpl Distribution, for: Distribution.Dice do
  import Distribution.Dice
  alias Distribution.Dice

  @spec dice([] | number) :: ((...) -> number)
  defp dice([]), do: dice(1.0)
  defp dice([avg: avg]), do: dice(avg)
  defp dice(avg), do: %Distribution.Uniform{list: [avg*1,avg*2,avg*3,avg*4,avg*5,avg*6]}

  @spec dice_gk4([] | number) :: ((...) -> number)
  defp dice_gk4([]), do: dice_gk4(1.0)
  defp dice_gk4([avg: avg]), do: dice_gk4(avg)
  defp dice_gk4(avg), do: %Distribution.Uniform{list: [avg*3,avg*4,avg*4,avg*5,avg*5,avg*6]}

  def skewness(%Dice{}), do: raise Distribution.FunctionNotSupportedError, message: "skewness is not supported for the Dice distribution"
  def kurtosis(%Dice{}), do: raise Distribution.FunctionNotSupportedError, message: "kurtosis is not supported for the Dice distribution"
  def size(%Dice{}), do: 0

  def cdf(%Dice{}), do: raise Distribution.FunctionNotSupportedError, message: "cdf is not supported for the Dice distribution"
  def pdf(%Dice{}), do: raise Distribution.FunctionNotSupportedError, message: "pdf is not supported for the Dice distribution"

  def random(%Dice{mode: :regular}), do: dice([]).()
  def random(%Dice{mode: :gk4}), do: dice_gk4([]).()
  
end
