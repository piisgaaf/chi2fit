defmodule Chi2fit.Math do

  # Copyright 2016-2021 Pieter Rijken
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

  @typedoc "Supported numerical integration methods"
  @type method :: :gauss | :gauss2 | :gauss3 | :romberg | :romberg2 | :romberg3


  @doc """
  Calculates the partial derivative of a function and returns the value.
  
  ## Examples

      The function value at a point:
      iex> der([3.0], fn [x]-> x*x end) |> Float.round(3)
      9.0

      The first derivative of a function at a point:
      iex> der([{3.0,1}], fn [x]-> x*x end) |> Float.round(3)
      6.0

      The second derivative of a function at a point:
      iex> der([{3.0,2}], fn [x]-> x*x end) |> Float.round(3)
      2.0

      Partial derivatives with respect to two variables:
      iex> der([{2.0,1},{3.0,1}], fn [x,y] -> 3*x*x*y end) |> Float.round(3)
      12.0

  """
  @default_h 0.001
  @spec der([float|{float,integer}], (([float])->float), Keyword.t) :: float
  def der(parameters, fun, options \\ []) do
    richardson(fn acc ->
        result = parameters
        |> expand_pars(acc)
        |> reduce_pars
        |> Enum.reduce(0.0, fn ({x,n,dx},sum) when is_list(x) -> sum+n*fun.(x)/dx end)
        {result,acc/2.0}
      end,
      @default_h,4.0,options)
  end

  @doc """
  Calculates the jacobian of the function at the point `x`.
  
  ## Examples
  
      iex> jacobian([2.0,3.0], fn [x,y] -> x*y end) |> Enum.map(&Float.round(&1))
      [3.0, 2.0]

  """
  @spec jacobian(x :: [float], (([float])->float)) :: [float]
  def jacobian(x, fun, options \\ []) do
    jacfun = &(jacobian(x, &1, fun, options))
    Enum.reduce(length(x)..1, [], fn (k,acc) -> [jacfun.(k)|acc] end)
  end

  ## TODO: implement gauss-kronrad integration (progressive gauss)
  @doc """
  Numerical integration providing Gauss and Romberg types.
  """
  @default_points 32
  @spec integrate(method, ((float)->float), a::float, b::float, options::Keyword.t) :: float
  def integrate(method, func, a, b, options \\ [])
  def integrate(:gauss, func, a, b, options) do
    npoints = options[:points] || @default_points

    factor_min = (b-a)/2.0
    factor_plus = (b+a)/2.0

    {weights,abscissa} = case npoints do
      4 ->
        {
          [ 0.6521451548625461,0.3478548451374538 ],
          [ 0.3399810435848563,0.8611363115940526 ]
         }
      8 ->
        {
          [ 0.3626837833783620,0.3137066458778873,0.2223810344533745,0.1012285362903763 ],
          [ 0.1834346424956498,0.5255324099163290,0.7966664774136267,0.9602898564975363 ]
        }
      32 ->
        {
          [ 0.0965400885147278,0.0956387200792749,0.0938443990808046,0.0911738786957639,0.0876520930044038,0.0833119242269467,0.0781938957870703,0.0723457941088485,0.0658222227763618,0.0586840934785355,0.0509980592623762,0.0428358980222267,0.0342738629130214,0.0253920653092621,0.0162743947309057,0.0070186100094701 ],
          [ 0.0483076656877383,0.1444719615827965,0.2392873622521371,0.3318686022821277,0.4213512761306353,0.5068999089322294,0.5877157572407623,0.6630442669302152,0.7321821187402897,0.7944837959679424,0.8493676137325700,0.8963211557660521,0.9349060759377397,0.9647622555875064,0.9856115115452684,0.9972638618494816 ]
        }
    end

    factor_min * (Enum.zip(abscissa,weights) |> Enum.map(fn {x,w} -> w*( func.(factor_min*x+factor_plus) + func.(-factor_min*x+factor_plus) ) end) |> Enum.sum)
  end
  def integrate(:gauss2, func, a, :infinity, options) do
    fac = 500.0 ## t = tanh(x/fac)
    fac*integrate(:gauss, fn t -> (func.(fac*:math.atanh(t)))/(1.0-t*t) end, :math.tanh(a/fac), 1.0, options)
  end
  def integrate(:gauss2, func, a, b, options) do
    fac = 500.0 ## t = tanh(x/fac)
    fac*integrate(:gauss, fn t -> (func.(fac*:math.atanh(t)))/(1.0-t*t) end, :math.tanh(a/fac), :math.tanh(b/fac), options)
  end
  def integrate(:gauss3, func, a, :infinity, options) do
    ## x = t/(1-t) = -1 + 1/(1-t), dx = dt/(1-t)^2
    integrate(:gauss, fn t -> (func.(t/(1.0-t)))/(1.0-t)/(1.0-t) end, a/(a+1.0), 1.0, options)
  end
  def integrate(:gauss3, func, a, b, options) do
    ## x = t/(1-t) = -1 + 1/(1-t), dx = dt/(1-t)^2
    integrate(:gauss, fn t -> (func.(t/(1.0-t)))/(1.0-t)/(1.0-t) end, a/(a+1.0), b/(b+1.0), options)
  end

  @default_tolerance 1.0e-6
  def integrate(:romberg, func, a, b, options) do
    richardson(fn acc ->
        case acc do
          [] ->
            f1 = try do func.(a) rescue _e -> 0.0 end
            f2 = try do func.(b) rescue _e -> 0.0 end
            result = (b-a) * ( f1 + f2 )/2.0
            {result,[{a,f1},{b,f2}]}

          values ->
            vals = values
            |> Stream.transform(nil, fn
              {x2,f},nil -> {[{x2,f}],x2}
              {x2,f},x1 -> {[{(x2+x1)/2.0,func.((x2+x1)/2.0)},{x2,f}],x2}
            end)
            |> Enum.to_list

            result = vals
            |> Stream.chunk_every(2,1,:discard)
            |> Stream.map(fn [{x1,f1},{x2,f2}] -> (x2-x1)*( f1 + f2 )/2.0 end)
            |> Enum.sum
            {result,vals}
        end
      end, [], 4.0, options)
  end
  def integrate(:romberg2, func, a, :infinity, options) do
    fac = 500.0 ## t = tanh(x/fac)
    integrate(:romberg, fn t -> (func.(fac*:math.atanh(t)))*fac/(1.0-t*t) end, :math.tanh(a/fac), 1.0, options)
  end
  def integrate(:romberg2, func, a, b, options) do
    fac = 500.0 ## t = tanh(x/fac)
    integrate(:romberg, fn t -> (func.(fac*:math.atanh(t)))*fac/(1.0-t*t) end, :math.tanh(a/fac), :math.tanh(b/fac), options)
  end
  def integrate(:romberg3, func, a, :infinity, options) do
    ## x = t/(1-t) = -1 + 1/(1-t), dx = dt/(1-t)^2
    integrate(:romberg, fn t -> (func.(t/(1.0-t)))/(1.0-t)/(1.0-t) end, a/(a+1.0), 1.0, options)
  end
  def integrate(:romberg3, func, a, b, options) do
    ## x = t/(1-t) = -1 + 1/(1-t), dx = dt/(1-t)^2
    integrate(:romberg, fn t -> (func.(t/(1.0-t)))/(1.0-t)/(1.0-t) end, a/(a+1.0), b/(b+1.0), options)
  end

  @doc """
  Richardson extrapolation.
  """
  @default_tolerance 1.0e-6
  @spec richardson(func::((term)->{float,term}), init::term, factor::float, results::[float], options::Keyword.t) :: float
  def richardson(func, init, factor, results \\ [], options)
  def richardson(func, init, factor, results, options) do
    tolerance = options[:tolerance] || @default_tolerance
    max = options[:itermax]

    {result,acc} = func.(init)
    {new,last,error,_} = results |> Enum.reduce({[],result,nil,factor}, fn
        _prev,{acc,item,0.0,order} ->
          {acc,item,0.0,order}
        prev,{acc,item,_,order} ->
          diff = (order*item - prev)/(order-1.0)
          {[diff|acc],diff,if(diff==0, do: 0.0, else: abs((diff-item)/diff)),order*factor}
        end)
    cond do
      max && (length(new) > max) ->
        last
      error < tolerance ->
        last
      true ->
        richardson(func, acc, factor, [result|Enum.reverse(new)], options)
    end
  end

  @doc """
  Newton-Fourier method for locating roots and returning the interval where the root is located.
  
  See [https://en.wikipedia.org/wiki/Newton%27s_method#Newton.E2.80.93Fourier_method]
  """
  @spec newton(a::float,b::float,func::((x::float)->float),maxiter::non_neg_integer,options::Keyword.t) :: {float, {float,float}, {float,float}}
  def newton(a,b,func,maxiter \\ 10, options), do: newton(a,b,func,maxiter,{(a+b)/2,{a,b},{nil,nil}},options)

  ##
  ## Local functions
  ##

  defp jacobian(x=[_|_], k, fun, options) when k>0 and k<=length(x) and is_function(fun,1) do
    x |> List.update_at(k-1, fn (val) -> {val,1} end) |> der(fun,options)
  end

  @default_rel_tolerance 1.0e-6
  defp newton(_a,_b,func,0,{root,{l,r},_},_options), do: {root,{l,r},{func.(l),func.(r)}}
  defp newton(a,b,func,maxiter,{prev,{left,right},{vleft,vright}},options) do
    tolerance = options[:tolerance] || @default_rel_tolerance

    x0 = func.(right)
    z0 = func.(left)

    if x0*z0 > 0 do
      raise ArgumentError, message: "Interval does not contain root"
    end

    derx0 = der([{right,1}], fn [x]->func.(x) end, options)

    if derx0 == 0 do
      raise ArithmeticError,
        message: "Interval contains local minimum/maximum [left/z0=#{left}/#{z0}; right/x0=#{right}/#{x0}; der=#{derx0}]"
    end

    x1 = right - x0/derx0
    z1 = left - z0/derx0
    root = (x1+z1)/2.0

    cond do
      z1 < left ->
        newton(a,b,func,0,{prev,{left,right},{vleft,vright}},options)

      x1 > right ->
        newton(a,b,func,0,{prev,{left,right},{vleft,vright}},options)

      z1 < x1 and abs(x1-z1) < tolerance ->
        newton(a,b,func,0,{root,{z1,x1},{z0,x0}},options)

      z1 > x1 and abs(x1-z1) < tolerance ->
        newton(a,b,func,0,{root,{x1,z1},{z0,x0}},options)

      z1 > x1 ->
        newton(a,b,func,maxiter-1,{prev,{x1,z1},{z0,x0}},options)

      true ->
        newton(a,b,func,maxiter-1,{root,{z1,x1},{z0,x0}},options)
    end
  end

  defp reduce_pars(list) do
    list |> Enum.reduce([{[],1,1.0}],
      fn
        (list,acc) when is_list(list) ->
          Enum.flat_map(list,
            fn
              ({{x,0,dx1}}) -> Enum.map(acc, fn ({y,n,dx2})->{[x|y],-n,dx1*dx2} end)
              ({x,0,dx1}) -> Enum.map(acc, fn ({y,n,dx2})->{[x|y],n,dx1*dx2} end)
            end)
        ({x,0,dx1},acc) -> Enum.map(acc, fn ({y,n,dx2})->{[x|y],n,dx1*dx2} end)
      end)
      |> Enum.map(fn ({l,n,dx}) -> {Enum.reverse(l),n,dx} end)
  end

  defp expand_pars(list,h) do
    list |> Enum.map(
          fn
            ({{x,0,factor}}) -> {{x,0,factor}}
            ({{x,0}}) -> {{x,0,1.0}}
            ({{x,n,factor}}) when n>0 ->
              xplus = x*(1.0+h)
              xmin = x*(1.0-h)
              dx = xplus-xmin
              [{{xplus,n-1,factor*dx}},{xmin,n-1,factor*dx}] |> expand_pars(h) |> List.flatten 
            ({{x,n}}) when n>0 ->
              xplus = x*(1.0+h)
              xmin = x*(1.0-h)
              dx = xplus-xmin
              [{{xplus,n-1,dx}},{xmin,n-1,dx}] |> expand_pars(h) |> List.flatten 
            ({x,0,factor}) -> {x,0,factor}
            ({x,0}) -> {x,0,1.0}
            ({x,n,factor}) when n>0 ->
              xplus = x*(1.0+h)
              xmin = x*(1.0-h)
              dx = xplus-xmin
              [{xplus,n-1,factor*dx},{{xmin,n-1,factor*dx}}] |> expand_pars(h) |> List.flatten 
            ({x,n}) when n>0 ->
              xplus = x*(1.0+h)
              xmin = x*(1.0-h)
              dx = xplus-xmin
              [{xplus,n-1,dx},{{xmin,n-1,dx}}] |> expand_pars(h) |> List.flatten 
            (x) when is_number(x) -> {x,0,1.0}
          end)
  end

end