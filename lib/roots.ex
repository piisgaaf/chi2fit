defmodule Chi2fit.Roots do

  # Copyright 2015-2017 Pieter Rijken
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
  Solves roots for linear, quadratic, and cubic equations.
  """

  @doc """
  Returns the real roots of polynoms of order 1, 2 and 3 as a list.

  ## Examples

      Solve `2.0*x + 5.0 = 0`
      iex> solve [2.0,5.0]
      [-2.5]

      iex> solve [2.0,-14.0,24.0]
      [4.0,3.0]

      iex> solve [1.0,0.0,5.0,6.0]
      [-0.9999999999999999]
  """
  @spec solve([float]) :: [float]
  def solve([+0.0|rest]), do: solve rest

  def solve([a1,a0]), do: [-a0/a1]

  def solve([a2,a1,a0]) do
    sqr = a1*a1-4*a2*a0
    cond do
      sqr == 0 -> -a1/2/a2
      sqr > 0 -> [(-a1+:math.sqrt(sqr))/2/a2,(-a1-:math.sqrt(sqr))/2/a2]
      true -> []
    end
  end

  def solve([1.0,+0.0,p,q]) do
    ## For details see equations (83) and (84) in http://mathworld.wolfram.com/CubicFormula.html
    c = -0.5*q*:math.pow(3/abs(p),1.5)
    cond do
      p>0 ->
        [:math.sinh(1.0/3.0*:math.asinh(c))]
      c>=1 ->
        [:math.cosh(1.0/3.0*:math.acosh(c))]
      c<=-1 ->
        [-:math.cosh(1.0/3.0*:math.acosh(abs(c)))]
      true ->
        ## Three real solutions
        [:math.cos(1.0/3.0*:math.acos(c)),:math.cos(1.0/3.0*:math.acos(c) + 2*:math.pi()/3.0),:math.cos(1.0/3.0*:math.acos(c) + 4*:math.pi()/3.0)]
    end
    |> Enum.map(&(&1*2*:math.sqrt(abs(p)/3.0)))
  end
  def solve([1.0,a2,a1,a0]), do: solve([1.0,0.0,(3*a1-a2*a2)/3.0,(2*a2*a2*a2-9*a1*a2+27*a0)/27.0]) |> Enum.map(&(&1-a2/3.0))
  def solve([a3,a2,a1,a0]), do: solve([1.0,a2/a3,a1/a3,a0/a3])

end
