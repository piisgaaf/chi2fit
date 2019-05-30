defmodule Distribution.Erlang do

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
  The Erlang distribution.
  """

  alias Exboost.Math, as: M
  
  defstruct [:pars,:batches]
  
  @type t() :: %__MODULE__{
    pars: [number()],
    batches: nil | pos_integer
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

  def skewness(%Erlang{pars: nil, batches: nil}), do: fn [k, _] -> 2/:math.sqrt(k) end
  def skewness(%Erlang{pars: nil, batches: k}) when k>0, do: fn [_] -> 2/:math.sqrt(k) end

  def kurtosis(%Erlang{pars: nil, batches: nil}), do: fn [k, _] -> 6/k end
  def kurtosis(%Erlang{pars: nil, batches: k}) when k>0, do: fn [_] -> 6/k end

  def size(%Erlang{batches: nil}), do: 2
  def size(%Erlang{batches: k}) when k>0, do: 1

  def cdf(%Erlang{pars: nil, batches: nil}), do: fn x, [k, lambda] -> erlangCDF(k,lambda).(x) end
  def cdf(%Erlang{pars: nil, batches: k}) when k>0, do: fn x, [lambda] -> erlangCDF(k,lambda).(x) end

  def pdf(%Erlang{pars: nil, batches: nil}), do: fn x, [k, lambda] -> :math.exp( -lambda*x + k*:math.log(lambda) + (k-1)*:math.log(x) - M.lgamma(k) ) end
  def pdf(%Erlang{pars: nil, batches: k}) when k>0, do: fn x, [lambda] -> :math.exp( -lambda*x + k*:math.log(lambda) + (k-1)*:math.log(x) - M.lgamma(k) ) end

  def random(%Erlang{pars: [k, lambda], batches: nil}), do: erlang(k,lambda).()
  def random(%Erlang{pars: [lambda], batches: k}), do: erlang(k,lambda).()
  def random(%Erlang{pars: nil, batches: nil}), do: fn [k,lambda] -> erlang(k,lambda).() end

end