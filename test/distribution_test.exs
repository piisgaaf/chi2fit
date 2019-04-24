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
  import Chi2fit.Distribution

  doctest Chi2fit.Distribution

  test "uniform distribution" do
    assert uniform(0,0).() == 0
    assert uniform(1,1).() == 1
    
    numbers = for _n<-List.duplicate(1, 25), do: uniform(1,5).()
    assert Enum.all? numbers, fn (n) -> n>=1 and n<=5 end
  end
  
  test "constant distribution" do
    assert constant([avg: 5]).() == 5
    
    numbers = for _n<-List.duplicate(1, 25), do: constant([avg: 1]).()
    assert Enum.all? numbers, fn (n) -> n==1 end
  end

  @tag long: true
  test "exponential distribution" do
    total = 10_000_000

    data = 1..total |> Stream.map(fn (_)->exponential([avg: 5.0]).() end)
    avg = Enum.sum(data)/total
    sd = :math.sqrt(Enum.reduce(data, 0, fn (x,sum)->sum+(x-avg)*(x-avg) end)/total)
    assert_in_delta 5.0, avg, 0.005
    assert_in_delta 5.0, sd, 0.005
  end

  @tag long: true
  test "exponential distribution - convolution" do
    total = 10_000_000
    chunk = 10

    data = 1..total |> Stream.map(fn (_)->exponential([avg: 5.0]).() end) |> Stream.chunk_every(chunk) |> Stream.map(&Enum.sum/1)
    avg = Enum.sum(data)/(total/chunk)
    sd = :math.sqrt(Enum.reduce(data, 0, fn (x,sum)->sum+(x-avg)*(x-avg) end)/(total/chunk))
    assert_in_delta 50.0, avg, 0.2
    assert_in_delta 50.0, sd*:math.sqrt(chunk), 0.1
  end
end
