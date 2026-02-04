defmodule Chi2fit.FFT do

  # Copyright 2016-2017 Pieter Rijken
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
  Provides Fast Fourier Transform.
  """

  require Integer

  import Kernel, except: [*: 2, /: 2,+: 2, -: 2]

  alias Chi2fit.Utilities, as: U

  @typedoc "A real number."
  @opaque real :: number

  @typedoc "A complex number with a real part and an imaginary part."
  @opaque complex :: {real,real}

  defp {x1,x2} * {y1,y2}, do: {x1*y1-x2*y2,x1*y2+x2*y1}
  defp x * {y1,y2}, do: {x*y1,x*y2}
  defp x * y, do: Kernel.*(x,y)

  defp {x1,x2} / y, do: {x1/y,x2/y}
  defp x / y, do: Kernel./(x,y)

  defp {x1,x2} + {y1,y2}, do: {x1+y1,x2+y2}
  defp x + {y1,y2}, do: {x+y1,y2}
  defp {x1,x2} + y, do: {x1+y,x2}
  defp x + y, do: Kernel.+(x,y)

  defp {x1,x2} - {y1,y2}, do: {x1-y1,x2-y2}
  defp x - {y1,y2}, do: {x-y1,-y2}
  defp x - y, do: Kernel.-(x,y)

  @doc """
  Calculates the discrete Fast Fourier Transform of a list of numbers.

  Provides a parallel version (see options below). See [1] for details of the algorithm implemented.

  ## Options
      `:phase` - Correction factor to use in the weights of the FFT algorithm. Defaults to 1.
      `:nproc` - Parellel version. Number of processes to use. See [2]. Defaults to 1.

  ## References

      [1] Zie: https://en.wikipedia.org/wiki/Cooley%E2%80%93Tukey_FFT_algorithm
      [2] Parallel version of FFT; see http://www.webabode.com/articles/Parallel%20FFT%20implementations.pdf

  ## Examples

      iex> fft [4]
      [{4.0, 0.0}]

      iex> fft [1,2,3,4,5,6]
      [{21.0, 0.0}, {-3.000000000000005, 5.196152422706631},
                {-3.0000000000000027, 1.7320508075688736}, {-3.0, 0.0},
                {-2.9999999999999987, -1.7320508075688794},
                {-2.999999999999992, -5.196152422706634}]

  """
  @spec fft([real],opts :: Keyword.t) :: [complex]
  def fft(list,opts \\ [])

  def fft([],_opts), do: []
  def fft([x,y],opts) do
    fac = opts[:phase] || 1
    [x*weight(fac*0,0,2)+y*weight(fac*0,1,2), x*weight(fac*1,0,2)+y*weight(fac*1,1,2)]
  end
  def fft(list=[_|_],opts) do
    fac = opts[:phase] || 1
    nproc = opts[:nproc] || 1

    nn = length(list)
    cond do
      Integer.is_even(length(list)) ->
        zipped = cond do
          nproc == 2 or nproc == 4 ->
            list
            |> split_evenodds
            |> Enum.map(fn x-> Task.async(fn -> fft(x,Keyword.merge(opts,[nproc: nproc/2])) end) end)
            |> Task.yield_many(3_600_000)
            |> Enum.map(fn ({_task,{:ok,result}})->result end)
            |> (&(apply(fn x,y->Stream.zip(x,y) end,&1))).()

          nproc == 1 ->
            list
            |> split_evenodds
            |> Enum.map(fn arg->fft(arg,opts) end)
            |> (&(apply(fn x,y->Stream.zip(x,y) end,&1))).()
        end

        n = nn/2
        zipped
        |> Stream.concat(zipped)
        |> Stream.with_index(0)
        |> Stream.map(
          fn
            ({{x,y},m}) when m<n -> x + (weight(fac*1,m,2*n)*y)
            ({{x,y},m}) when m>=n -> x - (weight(fac*1,m-n,2*n)*y)
          end)
        |> Enum.to_list

      true ->
        0..nn-1 |> Enum.map(
          fn m ->
            list |> Stream.with_index(0) |> Stream.map(fn ({item,k})-> item*weight(fac*m,k,nn) end) |> Enum.reduce(0,fn (x,acc)->x+acc end)
          end)
    end
  end

  @doc """
  Calculates the inverse FFT.

  For available options see `fft/2`.

  ## Examples

      iex> ifft [4.0]
      [{4.0, 0.0}]

      iex> ifft [1.0,2.0,3.0]
      [{2.0, 0.0}, {-0.5000000000000003, -0.2886751345948126},
                  {-0.49999999999999956, 0.28867513459481353}]

      iex> [1.0,5.0] |> fft |> ifft
      [{1.0, -3.061616997868383e-16}, {5.0, 6.123233995736767e-17}]

  """
  @spec ifft([real],Keyword.t) :: [complex]
  def ifft(list,opts \\ [nproc: 1]) do
    n = length(list)
    list |> fft(Keyword.merge(opts,[phase: -1])) |> Enum.map(&(&1/n))
  end

  @doc """
  Calculates the norm of a complex number or list of complex numbers.

  ## Examples

      iex> normv []
      []

      iex> normv {2,3}
      13

      iex> normv [{2,3},{1,2}]
      [13,5]

  """
  @spec normv([complex]|complex) :: real
  def normv({x,y}), do: x*x+y*y
  def normv(list) when is_list(list), do: list |> Enum.map(&normv/1)

  defp weight(r,m,n), do: weight(r*m,n)
  defp weight(rm,n), do: weight(rm/n)
  defp weight(x), do: {:math.cos(2*:math.pi()*x),U.normalize_zero(-:math.sin(2*:math.pi()*x))}

  defp split_evenodds(list) when Integer.is_even(length(list)) do
    list
    |> List.foldr({[[],[]],false},
      fn
        (item,{[e,o],true}) -> {[[item|e],o],false}
        (item,{[e,o],false}) -> {[e,[item|o]],true}
      end)
    |> elem(0)
  end

end
