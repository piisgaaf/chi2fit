defmodule Chi2fit.Matrix do

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
  This module provides matrix inverse operations and supporting functions.

  It provides 2 types of matrix norms and an iterative approach to calculating the matrix inverse.
  The implementation is based on the work [1].
  
  ## References
  
  [1] F. Soleymani, A Rapid Numerical Algorithm to Compute Matrix Inversion, International Journal of Mathematics and
  Mathematical Sciences, Volume 2012, Article ID 134653, doi:10.1155/2012/134653

  """

  @inverse_tolerance 1.0e-8
  @default_inverse_iterations 500
  @default_range 100
  @default_size 100
  @default_algorithm :lie

  @typedoc """
  A list of numbers
  """
  @type vector :: [number]
  
  @typedoc """
  A list of vectors (list of lists of numbers)
  """
  @type matrix :: [vector]

  import ExAlgebra.Matrix
  import ExAlgebra.Vector, only: [dot: 2]

  @doc """
  Constructs a unit matrix of size n. All diagonal elements have value 1 and the rest has value 0.
  
  ## Examples
  
      iex> unit(3)
      [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
      
      iex> unit(0)
      ** (ArgumentError) Illegal argument '0'

      iex> unit -1
      ** (ArgumentError) Illegal argument '-1'

      iex> unit 1.3
      ** (ArgumentError) Illegal argument '1.3'

  """
  @spec unit(n :: pos_integer) :: [[0|1]]
  def unit(n) when is_integer(n) and n > 0 do
    {result,_} = List.duplicate(0,n) |> List.duplicate(n) |> Enum.reduce({[],-1}, fn (list,{acc,m}) -> {[list |> List.replace_at(m,1)|acc],m-1} end)
    result
  end
  def unit(n), do: raise ArgumentError, message: "Illegal argument '#{inspect n}'"
  
  defp abssum(list), do: list |> Enum.map(&abs/1) |> Enum.sum

  @doc """
  Calculates the norm of the matrix as the sum of the absolutes value of all elements.
  
  ## Example
  
      iex> norm [[1,2,3],[4,5,6],[7,8,9]]
      45

  """
  @spec norm(matrix) :: number
  def norm(matrix), do: matrix |> Enum.map(&abssum/1) |> Enum.sum

  @doc """
  Calculates the norm of the matrix. All absolute values of the elements of each row are summed. The maximum
  value is returned
  
  ## Example
  
      iex> norm_1 [[1,2,3],[4,5,6],[7,8,9]]
      24

  """
  @spec norm_1(matrix) :: number
  def norm_1(matrix), do: matrix |> Enum.map(&abssum/1) |> Enum.max

  @doc """
  Calculates the norm of the matrix as `norm_1/1` but over the columns instead of over the rows.
  
  ## Example
  
      iex> norm_inf [[1,2,3],[4,5,6],[7,8,9]]
      18

  """
  @spec norm_inf(matrix) :: number
  def norm_inf(matrix), do: matrix |> transpose |> norm_1

  defmacrop telescope(matrix, []), do: quote(do: unit(length unquote(matrix)))
  defmacrop telescope(matrix, [a|rest]) do
    quote do
      mat = unquote(matrix)
      subtract(scalar_multiply(unit(length mat),unquote(a)),multiply(mat,telescope(mat, unquote(rest))))
    end
  end

  defp findv0(:way2, matrix, _options) do
    v0 = matrix |> transpose |> scalar_multiply(1.0/norm_1(matrix)/norm_inf(matrix))
    test = matrix |> length |> unit |> subtract(multiply(matrix,v0)) |> norm_1
    {v0,test}
  end

  defp findv0(:way3, matrix, options) do
    range = options[:range] || @default_range
    size = options[:size] || @default_size

    List.duplicate(0,size)
    |> Enum.map(fn (_x)->range*(2*:rand.uniform() - 1) end)
    |> Enum.reduce({nil,:infinity},fn
      (factor,{_,:infinity}) ->
        v0 = matrix |> length |> unit |> scalar_multiply(factor)
        test = matrix |> length |> unit |> subtract(multiply(matrix,v0)) |> norm_1 |> abs
        {v0,test}
      (factor,{v,error}) ->
        v0 = matrix |> length |> unit |> scalar_multiply(factor)
        test = matrix |> length |> unit |> subtract(multiply(matrix,v0)) |> norm_1 |> abs
        if test < error, do: {v0,test}, else: {v,error}
      end)
  end

  @doc """
  Returns the matrix inverse of the argument.
  
  ## Options

  * `:tolerance` - Iterate until the `norm_1/1` of I-AV is less than this value
  
  * `:algorithm` - Four algorithms are supported: `:hotelling_bodewig` (second order), `:lie` (third order),
  `:krishnamurthy_sen` (sixth order), and `:soleymani` (seventh order); defaults to `#{inspect @default_algorithm}`
  
  * `:max_iterations` - Maximum number of iterations to perform; defaults to #{@default_inverse_iterations}
  
  * `:range` - Range of values from -range to +range as a multiple of the unit matrix to try as an estimate
  of the inverse matric; defaults to #{@default_range}
  
  * `:size` - Number of tries to estimate initial inverse; defautls to #{@default_size}
  
  ## Examples
  
      iex> inverse [[3]]
      {:ok,[[0.3333333333333333]]}

      iex> inverse [[1,2],[3,4]]
      {:ok,[[-2.0, 1.0], [1.5, -0.5]]}
      
      iex> inverse [[3,2,0],[0,0,1],[2,-2,1]]
      {:ok,
        [[0.2000000000000267, -0.19999999999873883, 0.1999999999998383],
         [0.19999999999995877, 0.29999999999804805, -0.29999999999974974],
         [-1.204377690852969e-13, 0.9999999999943042, 7.3034853221406e-13]]}
       
      iex> inverse [[3,2,0],[0,0,1],[2,-2,1]], algorithm: :soleymani
      {:ok,
        [[0.2000000000000003, -0.19999999999999973, 0.19999999999999926],
         [0.2000000000000003, 0.29999999999999966, -0.29999999999999893],
         [-6.617444900424221e-23, 0.9999999999999988, -1.1183551388713516e-21]]}
      
      iex> inverse [[3,2,0],[0,0,1],[2,-2,1]], tolerance: 1.0e-15
      {:ok,
        [[0.20000000000000004, -0.2, 0.2],
         [0.20000000000000004, 0.3, -0.3],
         [-5.048709793414476e-29, 0.9999999999999999, -1.6087621596054126e-28]]}
      
  For matrices that have no inverse:
  
      iex> try do inverse [[1,2,3],[4,5,6],[7,8,9]] catch x->x end
      :no_inverse
  
  """
  @spec inverse(matrix, options :: Keyword.t) :: {:ok,inverse :: matrix} | :failed_to_find_v0 | :no_inverse | {:failed_to_reach_tolerance,inverse :: matrix,error :: float}
  def inverse(matrix, options \\ [])
  def inverse([[x]], _options), do: {:ok,[[1.0/x]]}
  def inverse([[x1,x2],[y1,y2]], _options), do: {:ok,[[y2,-x2],[-y1,x1]] |> scalar_multiply(1.0/(x1*y2-x2*y1))}
  def inverse(matrix, options) do
    max_iter = options[:max_iterations] || @default_inverse_iterations

    {v0,test} = findv0(:way2,matrix, options)
    try do
      if test < 2.0 do
        iterate(matrix,v0,max_iter,options)
      else
        {v0,test} = findv0(:way3,matrix,options)
        if test < 1.0 do
          iterate(matrix,v0,max_iter,options)
        else
          throw :no_v0
        end
      end
    catch
      {:impossible_inverse,v,error} ->
        throw {:failed_to_reach_tolerance,v,error}
      :no_v0 ->
        throw :failed_to_find_v0
    rescue
      _e in ArithmeticError ->
        throw :no_inverse
    end
  end

  defp forward(:hotelling_bodewig, av), do: telescope(av,[2.0])
  defp forward(:lie, av), do: telescope(av,[3.0,3.0])
  defp forward(:krishnamurthy_sen, av), do: telescope(av,[2.0]) |> multiply(telescope(av,[3.0,3.0]) |> multiply(telescope(av,[1.0,1.0])))
  defp forward(:soleymani, av), do: telescope(av,[120.0,393.0,735.0,861.0,651.0,315.0,93.0,15.0]) |> scalar_multiply(1/16.0)

  defp iterate(matrix,v0,0,_options) do
    throw {:impossible_inverse,v0,subtract(unit(length matrix),multiply(matrix,v0)) |> norm_1}
  end
  defp iterate(matrix,v0,max,options) when is_integer(max) and max > 0 do
    tolerance = options[:tolerance] || @inverse_tolerance
    algorithm = options[:algorithm] || @default_algorithm

    u = unit(length matrix)
    av = multiply(matrix,v0)
    test = subtract(u,av) |> norm_1
    unless test < tolerance do
      matrix |> iterate(multiply(v0,forward(algorithm, av)),max-1,options)
    else
      {:ok,v0}
    end
  end

  @doc """
  Returns the diagonal elements of the matrix as a vector.
  
  ## Example
  
      iex> diagonal [[1,2,3],[4,5,6],[7,8,9]]
      [1, 5, 9]

  """
  @spec diagonal(matrix) :: vector
  def diagonal(matrix) do
    matrix |> Enum.reduce({[],0}, fn (row, {acc,index})->{[Enum.at(row,index)|acc],index+1} end) |> elem(0) |> Enum.reverse
  end

  @doc """
  Returns a matrix with the supplied vector as its diagonal elements.
  
  ## Examples
  
      iex> from_diagonal [1,5,9]
      [[1, 0, 0], [0, 5, 0], [0, 0, 9]]
  
  """
  @spec from_diagonal(vector) :: matrix
  def from_diagonal(vector) do
    vector |> Enum.reduce({[],0}, fn (elem, {acc,index})->{[List.duplicate(0,length(vector)) |> List.replace_at(index,elem)|acc],index+1} end) |> elem(0) |> Enum.reverse
  end

  @doc """
  Calculates the inner product of two vectors.
  
  ## Examples
  
      iex> dotproduct [1,2], [3,4]
      11

      iex> dotproduct [], []
      0
  
      iex> dotproduct [1,2], []
      ** (ArgumentError) Vectors of unequal length

      iex> dotproduct [1,2], [1]
      ** (ArgumentError) Vectors of unequal length

  """
  @spec dotproduct(vector,vector) :: number
  def dotproduct(vector1, vector2) when length(vector1) != length(vector2) do
    raise ArgumentError, message: "Vectors of unequal length"
  end
  def dotproduct(vector1, vector2), do: dot(vector1,vector2)

end
