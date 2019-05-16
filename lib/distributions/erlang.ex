defmodule Distribution.Erlang do
  
  @moduledoc """
  The Erlang distribution.
  """

  alias Exboost.Math, as: M
  
  defstruct [:pars]
  
  @type t() :: %__MODULE__{
    pars: [number()]
  }

end

defimpl Distribution, for: Distribution.Erlang do
  import Distribution.Erlang
  alias Distribution.Erlang
  alias Exboost.Math, as: M

  defp erlang(k,lambda) when is_integer(k) and k>0 do
    fn ->
      -1.0/lambda*:math.log(1..k |> Enum.reduce(1.0, fn _,acc -> :rand.uniform()*acc end))
    end
  end
  
  defp erlangCDF(k,lambda) when k<0 or lambda<0, do: raise ArithmeticError, "Erlang is only defined for positive shape and mode"
  defp erlangCDF(k,lambda) when k>0, do: &M.gamma_p(k,lambda*&1)

  def skewness(%Erlang{pars: nil}), do: fn [k, _] -> 2/:math.sqrt(k) end
  def kurtosis(%Erlang{pars: nil}), do: fn [k, _] -> 6/k end
  def size(%Erlang{}), do: 2
  def cdf(%Erlang{pars: nil}), do: fn x, [k, lambda] -> erlangCDF(k,lambda).(x) end
  def pdf(%Erlang{pars: nil}), do: fn x, [k, lambda] -> :math.exp( -lambda*x + k*:math.log(lambda) + (k-1)*:math.log(x) - M.lgamma(k) ) end
  def random(%Erlang{pars: [k, lambda]}), do: erlang(k,lambda).()
end