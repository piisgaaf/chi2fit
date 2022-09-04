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

      `:tolerance` - Iterate until the `norm_1/1` of I-AV is less than this value
      `:algorithm` - Four algorithms are supported: `:hotelling_bodewig` (second order), `:lie` (third order),
          `:krishnamurthy_sen` (sixth order), and `:soleymani` (seventh order); defaults to `#{inspect @default_algorithm}`
      `:max_iterations` - Maximum number of iterations to perform; defaults to #{@default_inverse_iterations}
      `:range` - Range of values from -range to +range as a multiple of the unit matrix to try as an estimate
          of the inverse matric; defaults to #{@default_range}
      `:size` - Number of tries to estimate initial inverse; defautls to #{@default_size}
  
  ## Examples
  
      iex> inverse [[3]]
      {:ok,[[0.3333333333333333]]}

      iex> inverse [[1,2],[3,4]]
      {:ok,[[-2.0, 1.0], [1.5, -0.5]]}
      
      iex> inverse([[3,2,0],[0,0,1],[2,-2,1]]) |> elem(1) |> Enum.map(fn row -> Enum.map(row, & Float.round(&1,10)) end)
      [[0.2, -0.2, 0.2], [0.2, 0.3, -0.3], [0.0, 1.0, 0.0]]
       
      iex> inverse([[3,2,0],[0,0,1],[2,-2,1]], algorithm: :soleymani) |> elem(1) |> Enum.map(fn row -> Enum.map(row, & Float.round(&1,14)) end)
      [[0.2, -0.2, 0.2], [0.2, 0.3, -0.3], [0.0, 1.0, 0.0]]
      
      iex> inverse([[3,2,0],[0,0,1],[2,-2,1]], tolerance: 1.0e-15) |> elem(1) |> Enum.map(fn row -> Enum.map(row, & Float.round(&1,14)) end)
      [[0.2, -0.2, 0.2], [0.2, 0.3, -0.3], [0.0, 1.0, 0.0]]
      
  For matrices that have no inverse:
  
      iex> try do inverse [[1,2,3],[4,5,6],[7,8,9]] catch x->x end
      :no_inverse
  
  """
  @spec inverse(matrix, options :: Keyword.t) :: {:ok,inverse :: matrix} | :failed_to_find_v0 | :no_inverse | {:failed_to_reach_tolerance,inverse :: matrix,error :: float}
  def inverse(matrix, options \\ [])
  def inverse([[x]], _options), do: {:ok,[[1.0/x]]}
  def inverse([[x1,x2],[y1,y2]], _options) when x1*y2-x2*y1==0.0, do: throw {:impossible_inverse,"Zero determinant"}
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
    rescue
      _e in ArithmeticError ->
        throw :no_inverse
    catch
      {:impossible_inverse,v,error} ->
        throw {:failed_to_reach_tolerance,v,error}
      :no_v0 ->
        throw :failed_to_find_v0
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
  def dotproduct(vector1, vector2, sum \\ 0)
  def dotproduct(vector1, vector2, _sum) when length(vector1) != length(vector2) do
    raise ArgumentError, message: "Vectors of unequal length"
  end
  def dotproduct([], [], sum), do: sum
  def dotproduct([v1|rest1], [v2|rest2], sum), do: dotproduct(rest1, rest2, sum + v1*v2) 

  @doc """
  Subtracts two matrices and returns the result.
  """
  @spec subtract(matrix,matrix) :: matrix
  def subtract(matrix1, matrix2), do: subtract(matrix1,matrix2,[])
  defp subtract([], [], result), do: Enum.reverse(result)
  defp subtract([row1|rest1], [row2|rest2], result) do
    subtract(rest1, rest2, [subtractv(row1,row2)|result])
  end

  defp subtractv(vector1, vector2, result \\ [])
  defp subtractv([], [], result), do: Enum.reverse(result)
  defp subtractv([v1|rest1], [v2|rest2], result) do
    subtractv(rest1, rest2, [v1-v2|result])
  end

  defp scalar_multiply(matrix, scalar) when is_number(scalar) do
    matrix |> Enum.map(fn row -> Enum.map(row, & &1*scalar) end)
  end

  defp multiply(matrix1, matrix2, result \\ [])
  defp multiply([], _matrix2, result), do: Enum.reverse(result)
  defp multiply([row|rest1], matrix2, result) do
    multiply(rest1, matrix2, [matrix2
      |> Stream.zip
      |> Enum.reduce([], fn col, acc -> [dotproduct(row, Tuple.to_list(col))|acc] end)
      |> Enum.reverse | result])
  end

  @doc """
  Returns the tranpose of the matrix

  ## Examples:

      iex> transpose [ [1] ]
      [[1]]

      iex> transpose [ [1,2], [3,4] ]
      [[1, 3], [2, 4]]

  """
  @spec transpose(matrix) :: matrix
  def transpose(matrix), do: matrix |> Stream.zip |> Stream.map(&Tuple.to_list/1) |> Enum.into([])

  @doc """
  Adds two vectors.

  ## Examples

      iex> add [1,2], [3,4]
      [4,6]

      iex> add [], []
      []

      iex> add [1], [5]
      [6]

  """
  @spec add(vector, vector) :: vector
  def add(vector1, vector2), do: Stream.zip(vector1, vector2) |> Enum.map(fn {v1,v2} -> v1 + v2 end)

  @doc """
  Calculates the determinant of the matrix.
  """
  @spec det(matrix) :: number 
  def det([x]), do: x
  def det([ [a11,a12], [a21,a22]] ), do: a11*a22 - a12*a21
  
end
