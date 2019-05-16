defmodule Distribution.Weibull do

  @moduledoc """
  Weibull distribution.
  """

  defstruct [:pars]
  
  @type t() :: %__MODULE__{
    pars: [number()] | nil,
  }

end

defimpl Distribution, for: Distribution.Weibull do
  import Distribution.Weibull
  alias Distribution.Weibull

  import Exboost.Math, only: [tgamma: 1]

  @gamma53 0.902745292950933611297
  @gamma32 0.886226925452758013649

  @spec weibull(number, number|Keyword.t) :: ((...) -> number)
  defp weibull(1.0, [avg: average]), do: weibull(1.0, average)
  defp weibull(1.5, [avg: average]), do: weibull(1.5, average/@gamma53)
  defp weibull(2.0, [avg: average]), do: weibull(2.0, average/@gamma32)
  defp weibull(alpha, beta) when is_number(alpha) and is_number(beta) do
    fn ->
      u = :rand.uniform()
      beta*:math.pow(-:math.log(u),1.0/alpha)
    end
  end
  
  @spec weibullCDF(number,number) :: (number -> number)
  defp weibullCDF(k,_) when k<0, do: raise ArithmeticError, "Weibull is only defined for positive shape"
  defp weibullCDF(_,lambda) when lambda<0, do: raise ArithmeticError, "Weibull is only defined for positive scale"
  defp weibullCDF(k,lambda) when is_number(k) and is_number(lambda) do
    fn
      0 -> 0.0
      0.0 -> 0.0
      x when x<0 -> 0.0
      x ->
        lg = :math.log(x/lambda)*k
        cond do
          lg > 100.0 ->
            0.0
          lg < -18.0 ->
            ## With -18 (x/lambda)^2k < 10^(-16)
            t = :math.pow(x/lambda,k)
            t*(1 - 0.5*t)
          true ->
            1.0 - :math.exp -:math.pow(x/lambda,k)
        end
    end
  end

  def skewness(%Weibull{pars: nil}) do
    fn [k,lambda] ->
      mu = lambda*tgamma(1+1/k)
      sigma = lambda*:math.sqrt(tgamma(1+2/k) - tgamma(1+1/k)*tgamma(1+1/k))
      tgamma(1+3/k)*:math.pow(lambda/sigma,3) - 3*mu/sigma - :math.pow(mu/sigma,3)
    end
  end
  def kurtosis(%Weibull{pars: nil}) do
    fn [k,lambda] ->
      mu = lambda*tgamma(1+1/k)
      sigma = lambda*:math.sqrt(tgamma(1+2/k) - tgamma(1+1/k)*tgamma(1+1/k))
      skew = tgamma(1+3/k)*:math.pow(lambda/sigma,3) - 3*mu/sigma - :math.pow(mu/sigma,3)
      tgamma(1+4/k)*:math.pow(lambda/sigma,4) - 4*mu/sigma*skew - 6*:math.pow(mu/sigma,2) - :math.pow(mu/sigma,4) - 3.0
    end
  end
  def size(%Weibull{}), do: 2
  def cdf(%Weibull{pars: nil}), do: fn x,[k,lambda] -> weibullCDF(k,lambda).(x) end
  def pdf(%Weibull{pars: nil}), do: fn x, [k,lambda] -> k/lambda*:math.pow(x/lambda, k-1)*:math.exp( -:math.pow(x/lambda,k) ) end
  def random(%Weibull{pars: [k,lambda]}), do: weibull(k,lambda).()

end
