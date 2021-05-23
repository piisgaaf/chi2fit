defmodule Distribution.BiModal do

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

  defstruct [:weights,:distribs,name: "bimodal"]
  
  @type t() :: %__MODULE__{
    weights: [number()] | nil,
    distribs: [Distribution.t()] | nil,
    name: String.t
  }

end

defimpl Distribution, for: Distribution.BiModal do
  import Distribution.BiModal
  alias Distribution.BiModal

  def skewness(%BiModal{distribs: nil}), do: raise ArithmeticError, "Skewness not supported for BiModal distribution"
  def kurtosis(%BiModal{distribs: nil}), do: raise ArithmeticError, "Kurtosis not supported for BiModal distribution"
  
  def size(%BiModal{distribs: distribs}), do: 1 + (distribs|>Enum.map(&Distribution.size(&1))|>Enum.sum)
  
  def cdf(%BiModal{weights: nil, distribs: distribs}) do
    fn x,[w|parameters] ->
      distribs
      |> Enum.map(&{&1,Distribution.size(&1)})
      |> Enum.reduce({[],parameters},fn {d,size},{result,rest} -> {[{d,Enum.take(rest,size)}|result],Enum.drop(rest,size)} end)
      |> elem(0)
      |> Enum.reverse()
      |> Enum.zip([w,1-w])
      |> Enum.map(fn {tup,p} -> Tuple.append(tup,p) end)
      |> Enum.map(fn {d,pars,p} -> p*Distribution.cdf(d).(x,pars) end)
      |> Enum.sum
    end
  end
  
  def pdf(%BiModal{weights: nil, distribs: distribs}) do
    fn x,[w|parameters] ->
      distribs
      |> Enum.map(&{&1,Distribution.size(&1)})
      |> Enum.reduce({[],parameters},fn {d,size},{result,rest} -> {[{d,Enum.take(rest,size)}|result],Enum.drop(rest,size)} end)
      |> elem(0)
      |> Enum.reverse()
      |> Enum.zip([w,1-w])
      |> Enum.map(fn {tup,p} -> Tuple.append(tup,p) end)
      |> Enum.map(fn {d,pars,p} -> p*Distribution.pdf(d).(x,pars) end)
      |> Enum.sum
    end
  end

  def random(%BiModal{weights: nil, distribs: distribs}) do
	  fn [w|parameters] ->
			rnd = :rand.uniform()
			
      distribs
      |> Enum.map(&{&1,Distribution.size(&1)})
      |> Enum.reduce({[],parameters},fn {d,size},{result,rest} -> {[{d,Enum.take(rest,size)}|result],Enum.drop(rest,size)} end)
      |> elem(0)
      |> Enum.reverse()
    	|> Enum.zip([w,1])
      |> Enum.map(fn {tup,p} -> Tuple.append(tup,p) end)
    	|> Enum.map(fn {d,pars,p} -> {Distribution.random(d).(pars),p} end)
			|> Enum.reduce(nil, fn ({r,p},nil) -> if(rnd<p, do: r, else: nil); (_,acc) -> acc end)
		end
  end

  def name(model), do: model.name
  
end

defimpl Inspect, for: Distribution.BiModal do
  import Inspect.Algebra
  
  def inspect(dict, opts) do
    case {dict.weights,dict.distribs} do
      {_,nil} ->
        "#BiModal<>"
      {nil,[d1,d2]} ->
        concat ["#BiModal<", to_doc([d1,d2], opts), ">"]
      {[w],[d1,d2]} ->
        concat ["#BiModal<", "weights=(#{w},#{1-w})", "distribs=", to_doc([d1,d2], opts), ">"]
    end
  end

end