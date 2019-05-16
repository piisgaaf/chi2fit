defmodule Distribution.Wald do

  @moduledoc """
  Wald or Inverse Gauss distribution.
  """

  defstruct [:pars]
  
  @type t() :: %__MODULE__{
    pars: [number()] | nil,
  }

end

defimpl Distribution, for: Distribution.Wald do
  import Distribution.Wald
  alias Distribution.Wald

  @spec phi(x :: float) :: float
  defp phi(x) do
    d = %Distribution.Normal{pars: [0.0, 1.0]}
    Distribution.cdf(d).(x)
  end

  @spec wald(mu::number(),lambda::number()) :: ((...) -> number)
  defp wald(mu,lambda) when is_number(mu) and is_number(lambda) do
    d = %Distribution.Normal{pars: [0.0,1.0]}
    fn () ->
      w = Distribution.random(d)
      y = w*w
      x = mu + mu*mu*y/2/lambda - mu/2/lambda*:math.sqrt(4*mu*lambda*y + mu*mu*y*y)

      z = :rand.uniform()
      if z <= mu/(mu+x), do: x, else: mu*mu/x
    end
  end
  defp wald([avg: average],lambda), do: wald(average,lambda)

  @spec waldCDF(number,number) :: (number -> number)
  defp waldCDF(mu,_) when mu < 0, do: raise ArithmeticError, "Wald is only defined for positive average"
  defp waldCDF(_,lambda) when lambda < 0, do: raise ArithmeticError, "Wald is only defined for positive shape"
  defp waldCDF(mu,lambda) do
    fn
      x when x == 0 -> 0.0
      x when x < 0 -> 0.0
      x when x > 0 ->
        phi(:math.sqrt(lambda/x) * (x/mu-1.0)) + :math.exp(2.0*lambda/mu) * phi(-:math.sqrt(lambda/x) * (x/mu+1.0))
    end
  end

  def skewness(%Wald{pars: nil}), do: fn [k,lambda] -> 3*:math.sqrt(k/lambda) end
  def kurtosis(%Wald{pars: nil}), do: fn [k,lambda] -> 15*k/lambda end
  def size(%Wald{}), do: 2
  def cdf(%Wald{pars: nil}), do: fn x,[k,lambda] -> waldCDF(k,lambda).(x) end
  def pdf(%Wald{pars: nil}), do: fn x, [mu,lambda] -> :math.sqrt(lambda/2/:math.pi/x/x/x) * :math.exp( -lambda*(x-mu)*(x-mu)/2/x/mu/mu ) end
  def random(%Wald{pars: [k,lambda]}), do: wald(k,lambda).()

end
