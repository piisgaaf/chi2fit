defmodule Chi2fit.Distribution.Weibull do

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
  Weibull distribution.
  """

  defstruct [:pars, name: "weibull"]

  @type t() :: %__MODULE__{
    pars: [number()] | nil,
    name: String.t
  }

end

defimpl Chi2fit.Distribution, for: Chi2fit.Distribution.Weibull do
  alias Chi2fit.Distribution, as: D

  import D.Weibull
  alias D.Weibull

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
  defp weibullCDF(k,_) when k<0, do: raise(ArithmeticError, "Weibull is only defined for positive shape")
  defp weibullCDF(_,lambda) when lambda<0, do: raise(ArithmeticError, "Weibull is only defined for positive scale")
  defp weibullCDF(k,lambda) when is_number(k) and is_number(lambda) do
    fn
      0 -> 0.0
      +0.0 -> 0.0
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
            1.0 - :math.exp(-:math.pow(x/lambda,k))
        end
    end
  end

  @spec weibullPDF(number,number) :: (number -> number)
  defp weibullPDF(k,_) when k<0, do: raise(ArithmeticError, "Weibull is only defined for positive shape")
  defp weibullPDF(_,lambda) when lambda<0, do: raise(ArithmeticError, "Weibull is only defined for positive scale")
  defp weibullPDF(k,lambda) when is_number(k) and is_number(lambda) do
    fn
      0 -> 0.0
      +0.0 -> 0.0
      x when x<0 -> 0.0
      x ->
        t = :math.pow(x/lambda, k)
        k*t/x*:math.exp( -t )
    end
  end

  @spec skewness(%Weibull{}) :: ([number] -> float)
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
  def size(%Weibull{pars: nil}), do: 2
  def size(%Weibull{pars: [k,nil]}) when is_number(k), do: 1
  def size(%Weibull{pars: [nil,lambda]}) when is_number(lambda), do: 1
  def size(%Weibull{pars: [k,lambda]}) when is_number(k) and is_number(lambda), do: 0

  def cdf(%Weibull{pars: nil}), do: fn x,[k,lambda] -> weibullCDF(k,lambda).(x) end
  def cdf(%Weibull{pars: [k,nil]}), do: fn x,[lambda] -> weibullCDF(k,lambda).(x) end
  def cdf(%Weibull{pars: [nil,lambda]}), do: fn x,[k] -> weibullCDF(k,lambda).(x) end
  def cdf(%Weibull{pars: [k,lambda]}), do: fn x -> weibullCDF(k,lambda).(x) end

  def pdf(%Weibull{pars: nil}), do: fn x, [k,lambda] -> weibullPDF(k,lambda).(x) end
  def pdf(%Weibull{pars: [k,nil]}), do: fn x, [lambda] -> weibullPDF(k,lambda).(x) end
  def pdf(%Weibull{pars: [nil,lambda]}), do: fn x, [k] -> weibullPDF(k,lambda).(x) end
  def pdf(%Weibull{pars: [k,lambda]}), do: fn x -> weibullPDF(k,lambda).(x) end

  def random(%Weibull{pars: [k,lambda]}), do: weibull(k,lambda).()
  def random(%Weibull{pars: nil}), do: fn [k,lambda] -> weibull(k,lambda).() end

  def name(model), do: model.name

end

defimpl Inspect, for: Chi2fit.Distribution.Weibull do
  import Inspect.Algebra

  def inspect(dict, opts) do
    case dict.pars do
      nil ->
        "#Weibull<>"
      [alpha,beta] ->
        concat ["#Weibull<", to_doc("alpha=#{alpha}, beta=#{beta}", opts), ">"]
      list ->
        concat ["#Weibull<", to_doc(list, opts), ">"]
    end
  end

end
