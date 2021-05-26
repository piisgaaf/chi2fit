defmodule Chi2fit.Distribution.Exponential do

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
  The exponential distribution.
  """

  defstruct [:pars, name: "exponential"]

  @type t() :: %__MODULE__{
    pars: [number()] | nil,
    name: String.t
  }

end

defimpl Chi2fit.Distribution, for: Chi2fit.Distribution.Exponential do
  alias Chi2fit.Distribution, as: D

  import D.Exponential
  alias D.Exponential

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
  def cdf(%Exponential{pars: lambda}), do: fn x -> exponentialCDF(lambda).(x) end

  def pdf(%Exponential{pars: nil}), do: fn x, [lambda] -> lambda * :math.exp( -lambda*x ) end
  def random(%Exponential{pars: [lambda]}), do: exponential(lambda).()
  def random(%Exponential{pars: nil}), do: fn [lambda] -> exponential(lambda).() end

  def name(model), do: model.name

end

defimpl Inspect, for: Chi2fit.Distribution.Exponential do
  import Inspect.Algebra

  def inspect(dict, opts) do
    case dict.pars do
      nil ->
        "#Exponential<>"
      [rate] ->
        concat ["#Exponential<", to_doc("rate=#{rate}", opts), ">"]
      list ->
        concat ["#Exponential<", to_doc(list, opts), ">"]
    end
  end

end