defmodule AutocorrelationTest do

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

  use ExUnit.Case

  import Chi2fit.Statistics
  import Chi2fit.FFT
    
  @h 1.0e-5

  defp [] <~> [], do: true
  defp [x|restx] <~> [y|resty], do: (x <~> y) and (restx <~> resty)
  defp [] <~> [_|_], do: false
  defp [_|_] <~> [], do: false
  defp {x1,y1} <~> {x2,y2}, do: (x1<~>x2) and (y1<~>y2)
  defp x1 <~> {x2,y2}, do: (x1<~>x2) and (0<~>y2)
  defp x <~> y when is_number(x) and is_number(y), do: abs(x-y) < @h

  defp toms(t=%Time{microsecond: {ms,_}}) do
    (t.hour*3600 + t.minute*60 + t.second)*1_000_000 + ms
  end

  test "FFT trivial" do
    assert [] == fft []
  end

  test "FFT length returned" do
    assert 0 == length(fft [])
    assert 1 == length(fft [1])
    assert 2 == length(fft [1, 2])
    assert 11 == length(fft List.duplicate(1.1,11))
    assert 1001 == length(fft List.duplicate(1.1,1001)) ## some larger vector
    assert 1024 == length(fft List.duplicate(1.1,1024)) ## power of 2
    assert 60 == length(fft List.duplicate(1.1,60)) ## 4 * 15
  end

  test "FFT Single" do
    assert [{3.1,0.0}] == fft [3.1]
  end

  test "FFT Double" do
    assert [{3.7,0.0},{-1.5,0.0}] <~> fft([1.1, 2.6])
  end
  
  test "FFT odd" do
    assert [{6.0,0.0},{-1.5,0.8660254037844375},{-1.5,-0.8660254037844375}] <~> fft([1,2,3])
    assert [{15.0,0.0},{-2.5,3.44095},{-2.5,0.812299},{-2.5,-0.812299},{-2.5,-3.44095}] <~> fft([1,2,3,4,5])
  end

  test "FFT even" do
    assert [4.0,0.0,0.0,0.0] <~> fft([1,1,1,1])
    assert [10.0,{-2.0,2.0},-2.0,{-2.0,-2.0}] <~> fft([1,2,3,4])
    assert [21.0,{-3.0,5.196152},{-3.0,1.73205},-3.0,{-3.0,-1.73205},{-3.0,-5.196152}] <~> fft([1,2,3,4,5,6])
  end
  
  test "FFT - multi-core" do
    assert [4.0,0.0,0.0,0.0] <~> fft([1,1,1,1],nproc: 2)
    assert [4.0,0.0,0.0,0.0] <~> fft([1,1,1,1],nproc: 4)
  end

  @tag performance: true
  test "FFT - performance" do
    list = Stream.repeatedly(&:rand.uniform/0) |> Enum.take(4096*2)
    t0 = Time.utc_now() |> toms
    fft(list,nproc: 1)
    t1 = Time.utc_now() |> toms
    fft(tl(list),nproc: 1)
    t2 = Time.utc_now() |> toms
    
    IO.puts "\nPower of 2 test:"
    IO.puts "t1 = #{(t1-t0)/1_000_000} seconds"
    IO.puts "t2 = #{(t2-t1)/1_000_000} seconds"
  end

  @tag performance: true
  test "FFT - multi-core - performance" do
    list = Stream.repeatedly(&:rand.uniform/0) |> Enum.take(65536*8)
    t0 = Time.utc_now() |> toms
    fft(list,nproc: 1)
    t1 = Time.utc_now() |> toms
    fft(list,nproc: 2)
    t2 = Time.utc_now() |> toms
    fft(list,nproc: 4)
    t4 = Time.utc_now() |> toms
    
    IO.puts "\nFFT - parallel (multi-core) test:"
    IO.puts "single (nproc=1) = #{(t1-t0)/1_000_000} seconds"
    IO.puts "dual   (nproc=2) = #{(t2-t1)/1_000_000} seconds"
    IO.puts "quad   (nproc=4) = #{(t4-t2)/1_000_000} seconds"
  end

  test "Inverse FFT - trivial" do
    assert [] <~> ([] |> fft |> ifft)
    assert [1.1] <~> ([1.1] |> fft |> ifft)
    assert [1.1,2.6] <~> ([1.1,2.6] |> fft |> ifft)
  end
  
  test "Inverse FFT - nontrivial" do
    assert [1,2,3] <~> ([1,2,3] |> fft |> ifft)
    assert [1,2,3,4] <~> ([1,2,3,4] |> fft |> ifft)
    assert [1,2,3,4,5] <~> ([1,2,3,4,5] |> fft |> ifft)
    assert [1,2,3,4,5,6] <~> ([1,2,3,4,5,6] |> fft |> ifft)
  end

  test "Autocorrelation - trivial" do
    assert [] == auto []
    assert [1.21] <~> auto([1.1])
  end

  test "Autocorrelation - simple" do
    assert [5.0,2.0] <~> auto([1,2])
    assert [38.0,16.0,15.0] <~> auto([5,2,3])
    assert [13.0,10.0,6.0,2.0] <~> auto([1,2,2,2])
  end

  @tag timeout: 3_600_000
  @tag performance: true
  test "Autocorrelation - multi-core - performance" do
    list = Stream.repeatedly(&:rand.uniform/0) |> Enum.take(65536*2)
    t0 = Time.utc_now() |> toms
    auto(list,nproc: 1)
    t1 = Time.utc_now() |> toms
    auto(list,nproc: 2)
    t2 = Time.utc_now() |> toms
    auto(list,nproc: 4)
    t4 = Time.utc_now() |> toms
    
    IO.puts "\nAutocorrelation - parallel (multi-core) test:"
    IO.puts "single (nproc=1) = #{(t1-t0)/1_000_000} seconds"
    IO.puts "dual   (nproc=2) = #{(t2-t1)/1_000_000} seconds"
    IO.puts "quad   (nproc=4) = #{(t4-t2)/1_000_000} seconds"
  end

end
