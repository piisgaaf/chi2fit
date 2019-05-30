defmodule Chi2fitDistributionTest do

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

  doctest Distribution.Distribution.Constant, import: true
  doctest Distribution.Distribution.Uniform, import: true
  doctest Distribution.Distribution.Exponential, import: true
  doctest Distribution.Distribution.Poisson, import: true
  doctest Distribution.Distribution.Erlang, import: true
  doctest Distribution.Distribution.Normal, import: true

  @tag long: true
  test "exponential distribution" do
    dist = %Distribution.Exponential{pars: [5.0]}

    total = 10_000_000

    data = 1..total |> Stream.map(fn _ -> Distribution.random(dist) end)
    avg = Enum.sum(data)/total
    sd = :math.sqrt(Enum.reduce(data, 0, fn (x,sum)->sum+(x-avg)*(x-avg) end)/total)
    assert_in_delta 5.0, 1/avg, 0.005
    assert_in_delta 5.0, 1/sd, 0.005
  end

  @tag long: true
  test "exponential distribution - convolution" do
    dist = %Distribution.Exponential{pars: [5.0]}

    total = 10_000_000
    chunk = 10

    data = 1..total |> Stream.map(fn _ ->Distribution.random(dist) end) |> Stream.chunk_every(chunk) |> Stream.map(&Enum.sum/1)
    avg = Enum.sum(data)/(total/chunk)
    sd = :math.sqrt(Enum.reduce(data, 0, fn (x,sum)->sum+(x-avg)*(x-avg) end)/(total/chunk))
    assert_in_delta 5.0, chunk/avg, 0.2
    assert_in_delta 5.0, :math.sqrt(chunk)/sd, 0.1
  end
end
