defmodule Chi2fit.Distribution.Bernoulli do

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
  Provides the Bernoulli distribution.
  """

  @enforce_keys [:pars]
  defstruct [:pars, name: "bernoulli"]

  @type t() :: %__MODULE__{
    pars: [float],
    name: String.t
  }

end

defimpl Chi2fit.Distribution, for: Chi2fit.Distribution.Bernoulli do
  alias Chi2fit.Distribution, as: D

  import D.Bernoulli
  alias D.Bernoulli

  @spec bernoulli(value :: number) :: ((...) -> number)
  defp bernoulli(value) when is_number(value) do
   fn () ->
       u = :rand.uniform()
       if u <= value, do: 1, else: 0
   end
  end

  def skewness(%Bernoulli{pars: [p]}), do: fn _ -> (1-2*p)/:math.sqrt(p*(1.0-p)) end
  def kurtosis(%Bernoulli{pars: [p]}), do: fn _ -> (1-6*p*(1.0-p))/p/(1.0-p) end
  def size(%Bernoulli{}), do: 1

  def cdf(%Bernoulli{}), do: raise(D.FunctionNotSupportedError, message: "cdf is not supported for the Constant distribution")
  def pdf(%Bernoulli{}), do: raise(D.FunctionNotSupportedError, message: "pdf is not supported for the Constant distribution")
  def random(%Bernoulli{pars: [value]}), do: bernoulli(value).()

  def name(model), do: model.name

end

defimpl Inspect, for: Chi2fit.Distribution.Bernoulli do
  import Inspect.Algebra

  def inspect(dict, opts) do
    case dict.pars do
      nil ->
        "#Bernoulli<>"
      [p] ->
        concat ["#Bernoulli<", to_doc("p=#{p}", opts), ">"]
      list ->
        concat ["#Bernoulli<", to_doc(list, opts), ">"]
    end
  end

end
