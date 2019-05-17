defmodule Distribution.Uniform do

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

  def random(%Uniform{range: {min,max}}), do: uniform(min,max).()
  def random(%Uniform{list: list}), do: uniform(list).()
  def random(%Uniform{}), do: uniform([]).()
  
end
