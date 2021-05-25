defmodule Chi2fit.Distribution.SEP do

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
  The Skew Exponential Power cumulative distribution (Azzalini).

  ## Options
      `:method` - the integration method to use, :gauss and :romberg types are supported, see below
      `:tolerance` - re-iterate until the tolerance is reached (only for :romberg)
      `:points` - the number of points to use in :gauss method

  ## Integration methods
      `:gauss` - n-point Gauss rule,
      `:gauss2` - n-point Guass rule with tanh transformation,
      `:gauss3` - n-point Gauss rule with linear transformstion,
      `:romberg` - Romberg integration,
      `:romberg2` - Romberg integration with tanh transformation,
      `:romberg3` - Romberg integration with linear transformstion.
  """

  defstruct [:pars, :offset, options: [], name: "sep"]
  
  @type t() :: %__MODULE__{
    pars: [number()] | nil,
    offset: number() | nil,
    options: Keyword.t,
    name: String.t
  }

end

defimpl Chi2fit.Distribution, for: Chi2fit.Distribution.SEP do
  alias Chi2fit.Distribution, as: D

  import D.SEP
  alias D.SEP

  import Chi2fit.Math, only: [integrate: 5]

  @pi :math.pi()
  
  @spec sepCDF(a :: float,b :: float,lambda :: float,alpha :: float, options :: Keyword.t) :: (number -> number)
  defp sepCDF(a,b,lambda,alpha,options) do
    method = options[:method] || :romberg2
    endpoint = if method in [:gauss2,:gauss3,:romberg2,:romberg3], do: :infinity, else: 1000.0
    fn
      x ->
        result2 = integrate(method, sepPDF(a,b,lambda,alpha), 0.0, x, options)
        result3 = integrate(method, sepPDF(a,b,lambda,alpha), 0.0, endpoint, options)
        result2/result3
    end
  end

  @spec sepPDF(a::float,b::float,lambda::float,alpha::float) :: (number -> number)
  defp sepPDF(a,b,lambda,alpha) do
    fn x ->
      z = (x-a)/b
      t = :math.pow(abs(z),alpha/2.0)
      w = lambda*:math.sqrt(2.0/alpha)*t

      if z > 0.0 do
        :math.exp(-t*t/alpha) * 0.5 * ( 1.0 + :math.erf(w/:math.sqrt(2.0)) )
      else
        :math.exp(-t*t/alpha) * 0.5 * ( :math.erfc(w/:math.sqrt(2.0)) )
      end
    end
  end

  def skewness(%SEP{pars: nil}) do
    fn [_a,_b,lambda,_alpha] ->
      delta = lambda/:math.sqrt(1+lambda*lambda)
      0.5*(4-@pi)*:math.pow(delta*:math.sqrt(2/@pi),3)/:math.pow(1-2*delta*delta/@pi,1.5)
    end
  end
  def kurtosis(%SEP{pars: nil}) do
    fn [_a,_b,lambda,_alpha] ->
        delta = lambda/:math.sqrt(1+lambda*lambda)
        2*(@pi-3)*:math.pow(delta*:math.sqrt(2/@pi),4)/:math.pow(1-2*delta*delta/@pi,2)
    end
  end

  def size(%SEP{offset: nil}), do: 4
  def size(%SEP{offset: offset}) when is_number(offset), do: 3
  
  def cdf(%SEP{pars: nil, options: options}), do: fn x,[a,b,lambda,alpha] -> sepCDF(a,b,lambda,alpha,options).(x) end

  def pdf(%SEP{pars: nil}), do: fn x,[a,b,lambda,alpha] -> sepPDF(a,b,lambda,alpha).(x) end
  def random(%SEP{}), do: raise D.FunctionNotSupportedError, message: "random is not supported for the SEP distribution"

  def name(model), do: model.name
  
end

defimpl Inspect, for: Chi2fit.Distribution.SEP do
  import Inspect.Algebra
  
  def inspect(dict, opts) do
    case {dict.pars,dict.offset} do
      {nil,nil} ->
        "#SEP<>"
      {nil,offset} ->
        concat ["#SEP<", to_doc("offset=#{offset}", opts), ">"]
      {[scale,lambda,alpha],nil} ->
        concat ["#SEP<", to_doc("scale=#{scale}, lambda=#{lambda}, alpha=#{alpha}", opts), ">"]
      {[scale,lambda,alpha],offset} ->
        concat ["#SEP<", to_doc("offset=#{offset}, scale=#{scale}, lambda=#{lambda}, alpha=#{alpha}", opts), ">"]
      {list,offset} ->
        concat ["#SEP<", "offset=#{offset}, ", to_doc(list, opts), ">"]
    end
  end

end