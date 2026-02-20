defmodule UtilitiesTest do

  # Copyright 2017 Pieter Rijken
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

  use ExUnit.Case, async: true
  import Chi2fit.Math
  import Chi2fit.Times
  import Chi2fit.Utilities

  doctest Chi2fit.Math
  doctest Chi2fit.Times
  doctest Chi2fit.Utilities

  test "Third derivatives" do
    assert 6.0 == der([{3.0,3}], fn [x]-> x*x*x end) |> Float.round(3)
  end

  test "Newton-Fourier" do
    {root, {l, r}, _} = newton(2.0,3.0,fn x->x*x-6.25 end, tolerance: 1.0e-6)
    assert_in_delta 2.5, root, 1.0e-6
    assert_in_delta 2.5, l, 1.0e-6
    assert_in_delta 2.5, r, 1.0e-6

    {root, {l, r}, _} = newton(0.5,3.0,fn x->x*x-3.0 end, tolerance: 1.0e-4)
    assert_in_delta :math.sqrt(3), root, 1.0e-4
    assert_in_delta :math.sqrt(3), l, 1.0e-4
    assert_in_delta :math.sqrt(3), r, 1.0e-4
  end

  test "Newton-Fourier - no (unique) root in interval" do
    assert_raise ArgumentError, fn -> newton(3.0,3.5,fn x->x*x-6.25 end, []) end
  end

  test "Unzip of large tuples" do
    assert {[0.9047387007685814, 1.9047387007685814], [0.6192499721483649, 1.6192499721483649], [1.66747661178371, 2.66747661178371], [7.341904765315738, 8.341904765315737], [15.86694141620683, 16.86694141620683], [14.925856759588413, 15.925856759588413], [6.435352365675433, 7.435352365675433], [8.533850697685317, 9.533850697685317]} == unzip [
      {0.9047387007685814, 0.6192499721483649, 1.66747661178371, 7.341904765315738, 15.86694141620683, 14.925856759588413, 6.435352365675433, 8.533850697685317},
      {1.9047387007685814, 1.6192499721483649, 2.66747661178371, 8.341904765315738, 16.86694141620683, 15.925856759588413, 7.435352365675433, 9.533850697685317}
    ]
  end
end
