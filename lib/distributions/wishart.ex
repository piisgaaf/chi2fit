defmodule Chi2fit.Distribution.Wishart do

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
  Wishart distribution.
  """

  defstruct [:pars, :dim, name: "wishart"]

  @type t() :: %__MODULE__{
    pars: [number()] | nil,
    dim: [integer()],
    name: String.t
  }

end


defimpl Chi2fit.Distribution, for: Chi2fit.Distribution.Wishart do
  alias Chi2fit.Distribution, as: D
  alias Chi2fit.Matrix, as: M

  import D.Wishart
  alias D.Wishart

  import Exboost.Math, only: [gamma_p: 2, tgamma: 1]

  defp makerow([_],_,_pars,row), do: Enum.reverse(row)
  defp makerow([{_,_,_,rprev},{i,p,g,r}|pgrrest],[{_,q}|qrest],{a,pi,gi},row) do
    a = a - pi*rprev + 2*q/gi/g
    makerow([{i,p,g,r}|pgrrest],qrest,{a,pi,gi},[a|row])
  end

  defp makemat(pgrlist, qlist), do: makemat(pgrlist,qlist,[])
  defp makemat([], _, mat), do: Enum.reverse(mat)
  defp makemat(pgrlist, qlist, mat) do
    [{i,p,g,_}|pgrrest] = pgrlist
    makemat(pgrrest,tl(qlist),[makerow(pgrlist,qlist,{0,p,g},List.duplicate(0,i))|mat])
  end

  defp gamma(m,x) do
    :math.pow(:math.pi(),m*(m-1)/4)*Enum.reduce(1..m, 1.0, fn i,acc -> acc*tgamma(x-(i-1)/2) end)
  end

  @spec wishartCDF(number,number,number,number) :: (number -> number)
  defp wishartCDF(mu,scale,m,n) when is_number(mu) and is_number(scale) and scale>0 do
    nmin = min(m,n)
    nmax = max(m,n)

    alpha = (nmax - nmin - 1)/2

    k = :math.pow(:math.pi(),nmin*nmin/2)/:math.pow(2,nmin*nmax/2)/gamma(nmin,nmax/2)/gamma(nmin,nmin/2)

    nmat = if rem(nmin,2) == 0, do: nmin, else: nmin + 1
    kprime = k*:math.pow(2,alpha*nmat + nmat*(nmat+1)/2)*Enum.reduce(1..nmat, 1, fn i,acc -> acc*tgamma(alpha+i) end)

    fn
      x when x == mu -> 0.0
      x when x < mu -> 0.0
      x when x > mu ->
        xprime = (x-mu)/scale

        pgrlist = 1..nmin |> Enum.map(& {&1, gamma_p(alpha+&1, xprime/2), tgamma(alpha+&1), :math.exp(-xprime/2)*:math.pow(xprime/2,alpha+&1)/tgamma(alpha+1+&1)})
        qlist = 2..2*nmin-1 |> Enum.map(& {&1, :math.pow(2, -2*alpha-&1)*tgamma(2*alpha+&1)*gamma_p(2*alpha+&1,xprime)})

        mata = makemat(pgrlist,qlist)
        matt = M.transpose(mata)
        M.subtract(mata, matt) |> M.det |> :math.sqrt |> then(& kprime*&1)
    end
  end
  defp wishartCDF(_mu,_scale,_m,_n), do: raise(ArithmeticError, "Wishart is only defined for positive scale")

  def skewness(%Wishart{}), do: raise(ArithmeticError, "Skewness not supported for Wishart distribution")
  def kurtosis(%Wishart{}), do: raise(ArithmeticError, "Kurtosis not supported for Wishart distribution")
  def pdf(%Wishart{pars: nil}), do: raise(ArithmeticError, "pdf not supported for Wishart distribution")
  def random(%Wishart{pars: nil}), do: raise(ArithmeticError, "random not supported for Wishart distribution")

  def size(%Wishart{pars: nil}), do: 2

  def cdf(%Wishart{pars: nil, dim: [m,n]}), do: fn x,[mu,scale] -> wishartCDF(mu,scale,m,n).(x) end

  def name(model), do: model.name

end

defimpl Inspect, for: Chi2fit.Distribution.Wishart do
  def inspect(dict, opts) do
    import Inspect.Algebra

    [m,n] = dict.dim

    case dict.pars do
      nil ->
        "#Wishart<#{m},#{n}|[]>"
      list ->
        concat ["#Wishart<#{m},#{n}|", to_doc(list, opts), ">"]
    end
  end

end
