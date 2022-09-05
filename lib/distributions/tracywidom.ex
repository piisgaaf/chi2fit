defmodule Chi2fit.Distribution.TracyWidom do

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
  Tracy-Widom distribution.
  """

  defstruct [:pars, :type, name: "tracy-widom"]
  
  @type t() :: %__MODULE__{
    pars: [number()] | nil,
    type: 1|2|4,
    name: String.t
  }

end


defimpl Chi2fit.Distribution, for: Chi2fit.Distribution.TracyWidom do
  alias Chi2fit.Distribution, as: D

  import D.TracyWidom
  alias D.TracyWidom

  import Exboost.Math, only: [tgamma: 1, tgamma_lower: 2]

  # "Distribution of the largest eigenvalue for real Wishart and Gaussian random matrices and a simple approximation for the Tracy-Widom distribution", arXiv:1209.3394, Journal of Multivariate Analysis, Vol. 129, p. 69-81, 2014
  # See https://arxiv.org/pdf/1209.3394.pdf, Table 1
  @t1k 46.446
  @t1theta 0.186054
  @t1alpha 9.84801

  @t2k 79.6595
  @t2theta 0.101037
  @t2alpha 9.81961

  @t4k 146.021
  @t4theta 0.0595445
  @t4alpha 11.0016

  @e :math.exp(1.0)
  
  defp _gamma(k) when k > 1.0 do
    delta = k - :math.floor(k)
    
    u = :rand.uniform()
    v = :rand.uniform()
    w = :rand.uniform()
    
    if u <= @e/(@e+delta) do
      x = :math.pow(v, 1/delta)
      if w > :math.exp(-x), do: _gamma(k), else: x
    else
      x = 1.0 - :math.log(v)
      if w > :math.pow(x, delta-1), do: _gamma(k), else: x
    end
  end
  
  defp gamma(k,1.0) when k > 1.0 do
    _gamma(k) - (1..trunc(k) |> Enum.map(fn _ -> :math.log(:rand.uniform()) end) |> Enum.sum)
  end
  
  @spec tracywidom(number, number, number, number, number) :: ((...) -> number)
  defp tracywidom(mu, scale, k, theta, alpha) do
    fn -> mu + scale*( theta*gamma(k,1.0) - alpha ) end
  end
  
  @spec tracywidomCDF(number,number,number,number,number) :: (number -> number)
  defp tracywidomCDF(mu,scale,k,theta,alpha) when is_number(mu) and is_number(scale) do
    fn
      x when x == mu - scale*alpha -> 0.0
      x when x < mu - scale*alpha -> 0.0
      x when x > mu - scale*alpha ->
        1/tgamma(k)*tgamma_lower(k, (x - mu + scale*alpha)/theta/scale)
    end
  end

  @spec tracywidomPDF(number,number,number,number,number) :: (number -> number)
  defp tracywidomPDF(mu,scale,k,theta,alpha) when is_number(mu) and is_number(scale) do
    fn
      x when x == mu - scale*alpha -> 0.0
      x when x < mu - scale*alpha -> 0.0
      x when x > mu - scale*alpha ->
        1/tgamma(k)*:math.pow(theta*scale,-k)*:math.pow(x - mu + scale*alpha,k-1)*:math.exp(-(x - mu + scale*alpha)/theta/scale)
    end
  end

  defp mean(%TracyWidom{type: 1}), do: fn [mu,scale] -> (-1.2065335745820 - mu)/scale end
  defp mean(%TracyWidom{type: 2}), do: fn [mu,scale] -> (-1.7710868074110 - mu)/scale end
  defp mean(%TracyWidom{type: 4}), do: fn [mu,scale] -> (-2.3068848932410 - mu)/scale end
  
  defp variantie(d=%TracyWidom{type: 1}), do: fn [mu,scale] -> (1.6077810345810 - 2*mu*scale*mean(d).([mu,scale]) - mu*mu)/scale/scale end
  defp variantie(d=%TracyWidom{type: 2}), do: fn [mu,scale] -> (0.8131947928320 - 2*mu*scale*mean(d).([mu,scale]) - mu*mu)/scale/scale end
  defp variantie(d=%TracyWidom{type: 4}), do: fn [mu,scale] -> (0.5177237207726 - 2*mu*scale*mean(d).([mu,scale]) - mu*mu)/scale/scale end

  def skewness(d=%TracyWidom{type: 1}), do: fn [mu,scale] -> (2.0/:math.sqrt(@t1k) -3*mu*scale*scale*variantie(d).([mu,scale]) - 3*mu*mu*scale*mean(d).([mu,scale]) - mu*mu*mu)/scale/scale/scale end
  def skewness(d=%TracyWidom{type: 2}), do: fn [mu,scale] -> (2.0/:math.sqrt(@t2k) -3*mu*scale*scale*variantie(d).([mu,scale]) - 3*mu*mu*scale*mean(d).([mu,scale]) - mu*mu*mu)/scale/scale/scale end
  def skewness(d=%TracyWidom{type: 4}), do: fn [mu,scale] -> (2.0/:math.sqrt(@t4k) -3*mu*scale*scale*variantie(d).([mu,scale]) - 3*mu*mu*scale*mean(d).([mu,scale]) - mu*mu*mu)/scale/scale/scale end

  def kurtosis(d=%TracyWidom{type: 1}), do: fn [mu,scale] -> (0.1652429384 - 4*mu*scale*scale*scale*skewness(d).([mu,scale]) - 6*mu*mu*scale*scale*variantie(d).([mu,scale]) - 4*mu*mu*mu*scale*mean(d).([mu,scale]) - mu*mu*mu*mu)/scale/scale/scale/scale end
  def kurtosis(d=%TracyWidom{type: 2}), do: fn [mu,scale] -> (0.0934480876 - 4*mu*scale*scale*scale*skewness(d).([mu,scale]) - 6*mu*mu*scale*scale*variantie(d).([mu,scale]) - 4*mu*mu*mu*scale*mean(d).([mu,scale]) - mu*mu*mu*mu)/scale/scale/scale/scale end
  def kurtosis(d=%TracyWidom{type: 4}), do: fn [mu,scale] -> (0.0491951565 - 4*mu*scale*scale*scale*skewness(d).([mu,scale]) - 6*mu*mu*scale*scale*variantie(d).([mu,scale]) - 4*mu*mu*mu*scale*mean(d).([mu,scale]) - mu*mu*mu*mu)/scale/scale/scale/scale end

  def size(%TracyWidom{}), do: 2

  def cdf(%TracyWidom{pars: nil, type: 1}), do: fn x,[mu,scale] -> tracywidomCDF(mu,scale,@t1k,@t1theta,@t1alpha).(x) end
  def cdf(%TracyWidom{pars: nil, type: 2}), do: fn x,[mu,scale] -> tracywidomCDF(mu,scale,@t2k,@t2theta,@t2alpha).(x) end
  def cdf(%TracyWidom{pars: nil, type: 4}), do: fn x,[mu,scale] -> tracywidomCDF(mu,scale,@t4k,@t4theta,@t4alpha).(x) end

  def pdf(%TracyWidom{pars: nil, type: 1}), do: fn x,[mu,scale] -> tracywidomPDF(mu,scale,@t1k,@t1theta,@t1alpha).(x) end
  def pdf(%TracyWidom{pars: nil, type: 2}), do: fn x,[mu,scale] -> tracywidomPDF(mu,scale,@t2k,@t2theta,@t2alpha).(x) end
  def pdf(%TracyWidom{pars: nil, type: 4}), do: fn x,[mu,scale] -> tracywidomPDF(mu,scale,@t4k,@t4theta,@t4alpha).(x) end

  def random(%TracyWidom{pars: nil, type: 1}), do: fn [mu,scale] -> tracywidom(mu,scale,@t1k,@t1theta,@t1alpha).() end
  def random(%TracyWidom{pars: nil, type: 2}), do: fn [mu,scale] -> tracywidom(mu,scale,@t2k,@t2theta,@t2alpha).() end
  def random(%TracyWidom{pars: nil, type: 4}), do: fn [mu,scale] -> tracywidom(mu,scale,@t4k,@t4theta,@t4alpha).() end

  def name(model), do: model.name
  
end

defimpl Inspect, for: Chi2fit.Distribution.TracyWidom do
  def inspect(dict, opts) do
    import Inspect.Algebra
  
    type = cond do
      is_integer(dict.type) -> "_#{dict.type}"
      true -> ""
    end
        
    case dict.pars do
      nil ->
        "#TracyWidom#{type}<>"
      [mu,scale] ->
        concat ["#TracyWidom#{type}<", to_doc("mu=#{mu}, scale=#{scale}", opts), ">"]
      list ->
        concat ["#TracyWidom#{type}<", to_doc(list, opts), ">"]
    end
  end

end