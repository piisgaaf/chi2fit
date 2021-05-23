defmodule Distribution.TriModal do

  # Copyright 2020 Pieter Rijken
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
  Bimodal distribution.
  """

  defstruct [:weights,:distribs,name: "trimodal"]
  
  @type t() :: %__MODULE__{
    weights: [number()] | nil,
    distribs: [Distribution.t()] | nil,
    name: String.t
  }

end

defimpl Distribution, for: Distribution.TriModal do
  import Distribution.TriModal
  alias Distribution.TriModal

  def skewness(%TriModal{distribs: nil}), do: raise ArithmeticError, "Skewness not supported for TriModal distribution"
  def kurtosis(%TriModal{distribs: nil}), do: raise ArithmeticError, "Kurtosis not supported for TriModal distribution"
  
  def size(%TriModal{distribs: distribs}), do: 2 + (distribs|>Enum.map(&Distribution.size(&1))|>Enum.sum)
  
  def cdf(%TriModal{weights: nil, distribs: distribs}) do
    fn x,[w1,w2|parameters] ->
      distribs
      |> Enum.map(&{&1,Distribution.size(&1)})
      |> Enum.reduce({[],parameters},fn {d,size},{result,rest} -> {[{d,Enum.take(rest,size)}|result],Enum.drop(rest,size)} end)
      |> elem(0)
      |> Enum.reverse()
      |> Enum.zip([w1,(1-w1)*w2,(1-w1)*(1-w2)])
      |> Enum.map(fn {tup,p} -> Tuple.append(tup,p) end)
      |> Enum.map(fn {d,pars,p} -> p*Distribution.cdf(d).(x,pars) end)
      |> Enum.sum
    end
  end
  
  def pdf(%TriModal{weights: nil, distribs: distribs}) do
    fn x,[w1,w2|parameters] ->
      distribs
      |> Enum.map(&{&1,Distribution.size(&1)})
      |> Enum.reduce({[],parameters},fn {d,size},{result,rest} -> {[{d,Enum.take(rest,size)}|result],Enum.drop(rest,size)} end)
      |> elem(0)
      |> Enum.reverse()
      |> Enum.zip([w1,(1-w1)*w2,(1-w1)*(1-w2)])
      |> Enum.map(fn {tup,p} -> Tuple.append(tup,p) end)
      |> Enum.map(fn {d,pars,p} -> p*Distribution.pdf(d).(x,pars) end)
      |> Enum.sum
    end
  end

  def random(%TriModal{weights: [w1,w2], distribs: distribs}) do
    distribs
    |> Enum.zip([w1,(1-w1)*w2,(1-w1)*(1-w2)])
    |> Enum.map(fn {d,p} -> p*Distribution.random(d) end)
    |> Enum.sum
  end

  def name(model), do: model.name
  
end

defimpl Inspect, for: Distribution.TriModal do
  import Inspect.Algebra
  
  def inspect(dict, opts) do
    case {dict.weights,dict.distribs} do
      {_,nil} ->
        "#TriModal<>"
      {nil,[d1,d2,d3]} ->
        concat ["#TriModal<", to_doc([d1,d2,d3], opts), ">"]
      {[w1,w2],[d1,d2,d3]} ->
        concat ["#TriModal<", "weights=(#{w1},#{(1-w1)*w2},#{(1-w1)*(1-w2)})", "distribs=", to_doc([d1,d2,d3], opts), ">"]
    end
  end

end