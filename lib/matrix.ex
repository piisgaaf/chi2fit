defmodule Chi2fit.Matrix do

  @inverse_tolerance 1.0e-8
  @default_inverse_iterations 500

  @type matrix :: [[...]]
  @type vector :: [...] 

  import ExAlgebra.Matrix
  import ExAlgebra.Vector, only: [dot: 2]

  #######################################################################################################
  ## Inverse matrix stuff
  ##

  @spec unit(n :: pos_integer) :: [[0|1]]
  def unit(n) do
    {result,_} = List.duplicate(0,n) |> List.duplicate(n) |> Enum.reduce({[],-1}, fn (list,{acc,m}) -> {[list |> List.replace_at(m,1)|acc],m-1} end)
    result
  end

  defp abssum(list), do: list |> Enum.map(&abs/1) |> Enum.sum

  @spec norm(matrix) :: number
  def norm(matrix), do: matrix |> Enum.map(&abssum/1) |> Enum.sum

  @spec norm_1(matrix) :: number
  def norm_1(matrix), do: matrix |> Enum.map(&abssum/1) |> Enum.max

  @spec norm_inf(matrix) :: number
  def norm_inf(matrix), do: matrix |> transpose |> norm_1

  defmacrop telescope(matrix, []), do: quote(do: unit(length unquote(matrix)))
  defmacrop telescope(matrix, [a|rest]) do
    quote do
      mat = unquote(matrix)
      subtract(scalar_multiply(unit(length mat),unquote(a)),multiply(mat,telescope(mat, unquote(rest))))
    end
  end

  defp findv0(matrix, range \\ 100, size \\ 100) do
    {v, error} = List.duplicate(0,size)
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
    if error < 1.0, do: v, else: throw :no_v0
  end

  @spec inverse(matrix) :: matrix
  def inverse([[x]]), do: [[1.0/x]]
  def inverse([[x1,x2],[y1,y2]]), do: [[y2,-x2],[-y1,x1]] |> scalar_multiply(1.0/(x1*y2-x2*y1))
  def inverse(matrix) do
    require Logger

    v0 = matrix |> transpose |> scalar_multiply(1.0/norm_1(matrix)/norm_inf(matrix))
    test = matrix |> length |> unit |> subtract(multiply(matrix,v0)) |> norm_1
    if test < 2.0 do
      try do
        iterate(matrix,v0)
      catch
        {:impossible_inverse,v,_} ->
          Logger.warn "inverse: failed to reached tolerance"
          v
      end
    else
      v0 = findv0(matrix)
      iterate(matrix,v0)
    end
  end

  defp iterate(matrix,v0,maxn \\ @default_inverse_iterations)
  defp iterate(matrix,v0,0), do: throw {:impossible_inverse,v0,subtract(unit(length matrix),multiply(matrix,v0)) |> norm_1}
  defp iterate(matrix,v0,max) when is_integer(max) and max > 0 do
    u = unit(length matrix)
    test = subtract(u,multiply(matrix,v0)) |> norm_1
    unless test < @inverse_tolerance do
      matrix |> iterate(multiply(v0,telescope(multiply(matrix,v0),[2.0])),max-1)
    else
      v0
    end
  end

  @spec diagonal(matrix) :: vector
  def diagonal(matrix) do
    matrix |> Enum.reduce({[],0}, fn (row, {acc,index})->{[Enum.at(row,index)|acc],index+1} end) |> elem(0) |> Enum.reverse
  end

  @spec from_diagonal(vector) :: matrix
  def from_diagonal(vector) do
    vector |> Enum.reduce({[],0}, fn (elem, {acc,index})->{[List.duplicate(0,length(vector)) |> List.replace_at(index,elem)|acc],index+1} end) |> elem(0) |> Enum.reverse
  end

  @spec dotproduct(vector,vector) :: number
  def dotproduct(vector1, vector2), do: dot(vector1,vector2)

end
