defprotocol Distribution do

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

end

defimpl Enumerable, for: [
  Distribution.Constant,
  Distribution.Coin,
  Distribution.Dice,
  Distribution.Uniform,
  Distribution.Bernoulli,
  Distribution.Exponential,
  Distribution.Poisson,
  Distribution.Erlang,
  Distribution.Normal,
  Distribution.Wald,
  Distribution.Weibull,
  Distribution.Frechet,
  Distribution.Nakagami,
  Distribution.SEP,
  Distribution.BiModal,
  Distribution.TriModal] do

  alias Distribution, as: D

  def count(_lazy), do: {:error, __MODULE__}

  def member?(_lazy, _value), do: {:error, __MODULE__}

  def slice(_lazy), do: {:error, __MODULE__}

  def reduce(_distrib, {:halt, acc}, _fun), do: {:halted, acc}
  def reduce(distrib, {:suspend, acc}, fun), do: {:suspended, acc, &reduce(distrib, &1, fun)}
  def reduce(distrib, {:cont, acc}, fun), do: reduce(distrib, fun.(D.random(distrib), acc), fun)
end
