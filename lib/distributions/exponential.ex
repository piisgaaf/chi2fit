defmodule Distribution.Exponential do

  @moduledoc """
  The exponential distribution.
  """

  defstruct [:pars]
  
  @type t() :: %__MODULE__{
    pars: [number()] | nil
  }

end

defimpl Distribution, for: Distribution.Exponenbtial do
  import Distribution.Exponential
  alias Distribution.Exponential

  defp exponential([avg: average]) do
    fn ->
      u = :rand.uniform()
      -average*:math.log(u)
    end
  end
  defp exponential(rate), do: exponential([avg: 1.0/rate])
  defp exponentialCDF(rate) when rate > 0.0, do: fn t -> 1.0 - :math.exp(-rate*t) end

  def skewness(%Exponential{pars: nil}), do: fn _ -> 2 end
  def kurtosis(%Exponential{pars: nil}), do: fn _ -> 6 end
  def size(%Exponential{}), do: 1
  def cdf(%Exponential{pars: nil}), do: fn x, [lambda] -> exponentialCDF(lambda).(x) end
  def pdf(%Exponential{pars: nil}), do: fn x, [lambda] -> lambda * :math.exp( -lambda*x ) end
  def random(%Exponential{pars: [lambda]}), do: exponential(lambda).()
end
