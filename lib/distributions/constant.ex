defmodule Distribution.Constant do
  
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
  def random(%Constant{pars: [avg]}), do: constant(avg).()
  
end
