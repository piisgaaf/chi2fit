defprotocol Chi2fit.Distribution do

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

  defmodule FunctionNotSupportedError do
    defexception message: "Unsupported distribution function"
  end

  @spec skewness(t) :: number
  def skewness(distrib)
  
  @spec kurtosis(t) :: number
  def kurtosis(distrib)
  
  @spec size(t) :: number
  def size(distrib)
  
  @spec cdf(t) :: number
  def cdf(distrib)
  
  @spec pdf(t) :: number
  def pdf(distrib)
  
  @spec random(t) :: number
  def random(distrib)

  @spec name(t) :: String.t
  def name(distrib)

end

defimpl Enumerable, for: [
  Chi2fit.Distribution.Constant,
  Chi2fit.Distribution.Coin,
  Chi2fit.Distribution.Dice,
  Chi2fit.Distribution.Uniform,
  Chi2fit.Distribution.Bernoulli,
  Chi2fit.Distribution.Exponential,
  Chi2fit.Distribution.Poisson,
  Chi2fit.Distribution.Erlang,
  Chi2fit.Distribution.Normal,
  Chi2fit.Distribution.Wald,
  Chi2fit.Distribution.Weibull,
  Chi2fit.Distribution.Frechet,
  Chi2fit.Distribution.Nakagami,
  Chi2fit.Distribution.SEP,
  Chi2fit.Distribution.BiModal,
  Chi2fit.Distribution.TriModal] do

  alias Chi2fit.Distribution, as: D

  def count(_lazy), do: {:error, __MODULE__}

  def member?(_lazy, _value), do: {:error, __MODULE__}

  def slice(_lazy), do: {:error, __MODULE__}

  def reduce(_distrib, {:halt, acc}, _fun), do: {:halted, acc}
  def reduce(distrib, {:suspend, acc}, fun), do: {:suspended, acc, &reduce(distrib, &1, fun)}
  def reduce(distrib, {:cont, acc}, fun), do: reduce(distrib, fun.(D.random(distrib), acc), fun)

end
