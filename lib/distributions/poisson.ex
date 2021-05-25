defmodule Chi2fit.Distribution.Poisson do

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
  The Poisson distribution.
  
  For the implementation, see https://en.wikipedia.org/wiki/Poisson_distribution, 'Generating Poisson-distributed random variables'
  """

  alias Exboost.Math, as: M
  
  defstruct [:pars, period: 1.0, name: "poisson"]
  
  @type t() :: %__MODULE__{
    pars: [number()] | nil,
    period: number(),
    name: String.t
  }

end

defimpl Chi2fit.Distribution, for: Distribution.Poisson do
  alias Chi2fit.Distribution, as: D

  import D.Poisson
  alias D.Poisson
  alias Exboost.Math, as: M

  @step 500
  @expstep :math.exp(@step)
  defp iterate(p,r) when p<1 and r>0 and r>@step, do: iterate(p*@expstep,r-@step)
  defp iterate(p,r) when p<1 and r>0, do: iterate(p*:math.exp(r),0)
  defp iterate(p,r), do: {p,r}

  defp _poisson(rate, k \\ 0, p \\ 1.0)
  defp _poisson(rate,k,p) do
    k = k+1
    p = p*:rand.uniform()
    {p,rate} = iterate p,rate
    if p>1, do: _poisson(rate,k,p), else: k-1
  end

  defp poisson(rate), do: fn -> _poisson(rate) end
  defp poissonCDF(rate) when rate > 0.0 do
    fn t -> 1.0 - M.gamma_p(Float.floor(t+1.0),rate) end
  end

  def skewness(%Poisson{pars: nil, period: factor}), do: fn [lambda] -> 1/:math.sqrt(lambda*factor) end
  def kurtosis(%Poisson{pars: nil, period: factor}), do: fn [lambda] -> 1/lambda/factor end
  def size(%Poisson{}), do: 1
  def cdf(%Poisson{pars: nil, period: factor}), do: fn x, [lambda] -> poissonCDF(lambda*factor).(x) end
  def cdf(%Poisson{pars: [lambda], period: factor}), do: fn x -> poissonCDF(lambda*factor).(x) end
  def pdf(%Poisson{pars: nil, period: factor}), do: fn x, [lambda] -> :math.exp( x*:math.log(lambda*factor) - lambda*factor - M.lgamma(x+1) ) end
  def pdf(%Poisson{pars: [lambda], period: factor}), do: fn x -> :math.exp( x*:math.log(lambda*factor) - lambda*factor - M.lgamma(x+1) ) end

  def random(%Poisson{pars: [lambda], period: factor}), do: poisson(lambda*factor).()
  def random(%Poisson{pars: nil, period: factor}), do: fn [lambda] -> poisson(lambda*factor).() end

  def name(model), do: model.name

end

defimpl Inspect, for: Chi2fit.Distribution.Poisson do
  import Inspect.Algebra
  
  def inspect(dict, opts) do
    case dict.pars do
      nil ->
        "#Poisson<>"
      [rate] ->
        concat ["#Poisson<", to_doc("rate=#{rate}", opts), ">"]
      list ->
        concat ["#Poisson<", to_doc(list, opts), ">"]
    end
  end

end
