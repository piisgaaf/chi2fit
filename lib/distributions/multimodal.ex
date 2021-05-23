defmodule Distribution.MultiModal do

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

  defstruct [:weights,:distribs,name: "multimodal"]
  
  @type t() :: %__MODULE__{
    weights: [number()] | nil,
    distribs: [Distribution.t()] | nil,
    name: String.t
  }

  @spec weights([number()]) :: [number()]
  def weights(list) do
    list ++ [1.0]
    |> Enum.reduce({[],1.0}, fn w,{result,last} -> {[last*w|result],last*(1-w)} end)
    |> elem(0)
    |> Enum.reverse
  end

end

defimpl Distribution, for: Distribution.MultiModal do
  require Distribution.MultiModal
  alias Distribution.MultiModal

  def skewness(%MultiModal{distribs: nil}), do: raise ArithmeticError, "Skewness not supported for MultiModal distribution"
  def kurtosis(%MultiModal{distribs: nil}), do: raise ArithmeticError, "Kurtosis not supported for MultiModal distribution"
  
  def size(%MultiModal{distribs: distribs}), do: length(distribs)-1 + (distribs|>Enum.map(&Distribution.size(&1))|>Enum.sum)
  
  def cdf(%MultiModal{weights: nil, distribs: distribs}) do
    fn x,list when is_list(list) ->
      {weights,parameters} = Enum.split(list,length(distribs))
      distribs
      |> Enum.map(&{&1,Distribution.size(&1)})
      |> Enum.reduce({[],parameters},fn {d,size},{result,rest} -> {[{d,Enum.take(rest,size)}|result],Enum.drop(rest,size)} end)
      |> elem(0)
      |> Enum.reverse()
      |> Enum.zip(MultiModal.weights(weights))
      |> Enum.map(fn {tup,p} -> Tuple.append(tup,p) end)
      |> Enum.map(fn {d,pars,p} -> p*Distribution.cdf(d).(x,pars) end)
      |> Enum.sum
    end
  end
  
  def pdf(%MultiModal{weights: nil, distribs: distribs}) do
    fn x,list when is_list(list) ->
      {weights,parameters} = Enum.split(list,length(distribs))
      distribs
      |> Enum.map(&{&1,Distribution.size(&1)})
      |> Enum.reduce({[],parameters},fn {d,size},{result,rest} -> {[{d,Enum.take(rest,size)}|result],Enum.drop(rest,size)} end)
      |> elem(0)
      |> Enum.reverse()
      |> Enum.zip(MultiModal.weights(weights))
      |> Enum.map(fn {tup,p} -> Tuple.append(tup,p) end)
      |> Enum.map(fn {d,pars,p} -> p*Distribution.pdf(d).(x,pars) end)
      |> Enum.sum
    end
  end

  def random(%MultiModal{weights: weights, distribs: distribs}) do
    distribs
    |> Enum.zip(MultiModal.weights(weights))
    |> Enum.map(fn {d,p} -> p*Distribution.random(d) end)
    |> Enum.sum
  end

  def name(model), do: model.name
  
end

defimpl Inspect, for: Distribution.MultiModal do
  import Inspect.Algebra
  alias Distribution.MultiModal
  
  def inspect(dict, opts) do
    case {dict.weights,dict.distribs} do
      {_,nil} ->
        "#MultiModal<>"
      {nil,list=[_|_]} ->
        concat ["#MultiModal<", to_doc(list, opts), ">"]
      {weights=[_|_],list=[_|_]} ->
        concat ["#MultiModal<", "weights=(",Enum.join(MultiModal.weights(weights),","),"),", "distribs=", to_doc(list, opts), ">"]
    end
  end

end