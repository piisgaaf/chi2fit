defmodule Chi2fit do

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

  require Logger

  import Chi2fit.Matrix
  import Chi2fit.Utilities

  @type observable :: {x :: float, y :: float, dy :: float}
  @type observables :: [observable]
  @type model :: {[float], ((float,[float])->float)}

  @type chi2 :: float
  @type cov :: Matrix.matrix
  @type params :: [{float,float}]
  
  @arithmic_penalty 1_000_000_000
  
  defp nopenalties(_,_), do: 0.0

  defp dchi2_simple(y, y1, y2,f),              do: (f-y)/abs(y-(y1+y2)/2)
  defp dchi2_asimple(y, y1,_y2,f) when f<y,    do: (f-y)/(y1-y)
  defp dchi2_asimple(y,_y1,_y2,f) when y == f, do: 0.0
  defp dchi2_asimple(y,_y1, y2,f),             do: (f-y)/(y2-y)

  defp dchi2_linear(y,y1,y2,f) do
    delta = f-y
    splus = y2-y
    smin = y-y1

    cond do
      # Special cases
      f==1.0 and y2==1.0 -> 1.0
      f==0.0 and y1==0.0 -> 1.0
      f==y -> 0.0
      # Extreme punishment
      f==1.0 -> 1_000_000
      f==0.0 -> 1_000_000
      # Model
      # 
      # Linear transformation that:
      # - is continuous in u=0,
      # - passes through the point sigma+ at u=1,
      # - asymptotically reaches 1-y at u->infinity
      # - pass through the point -sigma- at u=-1,
      # - asymptotically reaches -y at u->-infinity
      #
      delta>0 -> (1.0-y2)/(1.0-f) * delta/splus
      true ->          y1/f       * delta/smin
    end
  end

  @spec chi2(observables, ((float)->float), ((float)->float), Keyword.t) :: float
  def chi2(observables, fun, penalties \\ fn (_)->0.0 end, options \\ [])
  def chi2(observables, fun, penalties, []), do: chi2(observables, fun, penalties, [model: :simple])
  def chi2(observables, fun, penalties, options) do
    observables
    |> Stream.map(
      fn
        ({x,y,dy}) ->
          # Symmetric errors
          tmp = (y-fun.(x))/dy
          tmp*tmp + penalties.(x)
        ({x,y,y1,y2}) ->
          ## Carefully handle asymmetric errors
          ## See Bohm (DESY), formula (8.5)
          ## See https://arxiv.org/pdf/physics/0401042v1.pdf
          try do
            tmp = case options[:model] do
              :linear -> dchi2_linear y,y1,y2,fun.(x)
              :simple -> dchi2_simple y,y1,y2,fun.(x)
              :asimple -> dchi2_asimple y,y1,y2,fun.(x)
            end
            tmp*tmp + penalties.(x)
          rescue
            ArithmeticError -> @arithmic_penalty
          end
      end)
    |> Enum.sum
  end

  defp beta(observables, {parameters, fun, penalties}) do
    betafun = &(beta({&1,&2}, observables, {parameters, fun, penalties}))
    Enum.reduce(length(parameters)..1, [], fn
      (k,acc) -> [
        Enum.reduce(length(parameters)..1, [], fn (j,acc)->[betafun.(k, j)|acc] end)
        |acc]
      end)
  end

  @spec beta({pos_integer,pos_integer}, observables, model) :: float
  def beta(index, observables, {parameters, fun}), do: beta(index, observables, {parameters, fun, &nopenalties/2})
  def beta({k,j}, observables, {parameters, fun, _penalties}) when k>0 and k<=length(parameters) and j>0 and j<=length(parameters) do
    params_k = parameters |> List.update_at(k-1, fn (val) -> {val,1} end)
    params_j = parameters |> List.update_at(j-1, fn (val) -> {val,1} end)

    observables
    |> Stream.map(
      fn
        ({x,_y,dy}) ->
          der(params_k,&fun.(x,&1))*der(params_j,&fun.(x,&1))/dy/dy
        ({x,_y,y1,y2}) ->
          dy = max(0.000001,y2-y1)/2
          der(params_k,&fun.(x,&1))*der(params_j,&fun.(x,&1))/dy/dy
      end)
    |> Enum.sum
  end

  defp gamma(observables, {parameters, fun, penalties, options}) do
    gammafun = &(gamma(&1,observables, {parameters, fun,penalties, options}))
    Enum.reduce(length(parameters)..1, [], fn (k,acc)->[gammafun.(k)|acc] end)
  end

  @spec gamma(pos_integer, observables, model) :: float
  def gamma(k, observables, {parameters, fun, penalties, options}) when k>0 and k<=length(parameters) do
    params_k = parameters |> List.update_at(k-1, fn (val) -> {val,1} end)
    -0.5*der(params_k, fn (pars)->chi2smooth(observables, pars, {fun,penalties},options[:smoothing],options) end)
  end

  defp alpha(observables, {parameters, fun, penalties, options}) do
    alphafun = &(alpha({&1,&2}, observables, {parameters, fun, penalties,options}))
    Enum.reduce(length(parameters)..1, [], fn
      (k,acc) -> [
        Enum.reduce(length(parameters)..1, [], fn (j,acc)->[alphafun.(k, j)|acc] end)
        |acc]
      end)
  end

  defp derive_par(list, index), do: list |> List.update_at(index-1, fn (val) when is_number(val) -> {val,1}; ({val,n}) -> {val,n+1} end)

  @spec alpha({pos_integer,pos_integer}, observables, model) :: float
  def alpha({k,j}, observables, {parameters, fun, penalties, options}) when k>0 and k<=length(parameters) and j>0 and j<=length(parameters) do
    params_kj = parameters |> derive_par(k-1) |> derive_par(j-1)
    0.5*der(params_kj,fn (pars)->chi2smooth(observables, pars, {fun,penalties},options[:smoothing],options) end)
  end

  #######################################################################################################
  ## Chi squared fit
  ##

  defp chi2smooth(observables,parameters,{fun,penalties},true,options) do
    rx = 5.0e-4
    ry = 5.0e-3
    n = 1
    (for dx<- -n..n, dy<- -n..n, do: {rx*dx,ry*dy})
    |> Stream.map(fn ({dx,dy})-> [p1,p2]=parameters; [p1+dx,p2+dy] end)
    |> Stream.map(fn (pars)-> chi2(observables, &(fun.(&1,pars)), &(penalties.(&1,pars)), options)/(2*n+1)/(2*n+1) end)
    |> Enum.sum
  end
  defp chi2smooth(observables,parameters,{fun,penalties},false,options) do
    chi2(observables, &(fun.(&1,parameters)), &(penalties.(&1,parameters)), options)
  end

  defp sample(list) do
    list |> Enum.map(fn
        ({low,high})->low + :rand.uniform()*(high-low)
        (x)->x
      end)
  end

  def chi2probe(observables, parranges, fun_penalties, options) do
    chi2probe(observables, parranges, fun_penalties, options[:num], nil, options)
  end

  defp chi2probe(_observables, _parranges, {_fun,_penalties}, 0, best, _options) do
    ## Refactor this!!!!!
    {chi2,parameters,saved} = best
    {_chis,plists} = saved |> Enum.unzip
    {plist1,plist2} = plists |> Stream.map(&List.to_tuple/1) |> Enum.unzip
    {chi2,parameters,{[Enum.min(plist1),Enum.max(plist1)],[Enum.min(plist2),Enum.max(plist2)]}}
  end
  defp chi2probe(observables, parranges, {fun,penalties}, num, best, options) do
    if options[:progress] do
      cond do
        options[:mark][:m] and rem(num,1000) == 0 -> IO.write "M"
        options[:mark][:c] and rem(num,100) == 0 -> IO.write "C"
        options[:mark][:x] and rem(num,10) == 0 -> IO.write "x"
        true -> :ok
      end
    end
    try do
      parameters = parranges |> sample
      chi2 = chi2smooth observables,parameters,{fun,penalties},options[:smoothing],options
      if options[:print?] do
        parameters |> Enum.each(fn (p)->IO.binwrite options[:print], "#{p} " end)
        IO.binwrite options[:print], "#{chi2}\n"
      end
      options[:save] && options[:save].(parameters,chi2)
      chi2probe(observables, parranges, {fun,penalties}, num-1,
        case best do
          nil ->
            if options[:debug], do: Logger.debug "debug: chi2 -> #{chi2} #{inspect parameters}"
            {chi2,parameters,[{chi2,parameters}]}
          {oldchi2,_,saved} when chi2<oldchi2 ->
            if options[:progress], do: IO.write "*"
            if options[:debug], do: Logger.debug "debug: chi2 -> #{inspect chi2} #{inspect parameters}"
            {chi2,parameters,[{chi2,parameters}|Enum.filter(saved,fn ({x,_})-> x < chi2+1.0 end)]}
          {oldchi2,oldpars,saved} when chi2<oldchi2+1.0 ->
            {oldchi2,oldpars,[{chi2,parameters}|saved]}
          _else ->
            best
        end,
        options)
    rescue
      ArithmeticError ->
        chi2probe(observables, parranges, {fun,penalties}, num-1, best, options)
      err ->
        Logger.debug "\nError: #{inspect err} #{inspect System.stacktrace}"
        reraise err, "Error!"
    end
  end

  defp vary_params(parameters, num_variations \\ 100) when is_list(parameters) do
    -1..length(parameters)
    |> Stream.map(&(List.duplicate(&1,num_variations)))
    |> Stream.concat
    |> Stream.flat_map(
      fn
        (-1) -> [List.duplicate(:rand.uniform(),length(parameters)), List.duplicate(:rand.uniform()/10_000,length(parameters))]
        (0) -> [List.duplicate(0.0,length(parameters)) |> Enum.map(fn (_)->:rand.uniform() end)]
        (n) when is_integer(n) and n>0 -> [List.duplicate(0.0,length(parameters)) |> List.replace_at(n-1, :rand.uniform()),List.duplicate(0.0,length(parameters)) |> List.replace_at(n-1, :rand.uniform()/10_000)]
      end)
  end

  @spec chi2fit(observables, model, pos_integer, Keyword.t) :: {chi2,cov,params}
  def chi2fit(observables, model, max \\ 100, error \\ nil, options \\ [debug: false])
  def chi2fit(observables, {parameters, fun}, max, error, options), do: chi2fit observables, {parameters, fun, &nopenalties/2}, max, error, options
  def chi2fit(observables, {parameters, fun, penalties}, 0, {cov,_error}, options) do
    {chi2(observables, &(fun.(&1,parameters)), &(penalties.(&1,parameters)), options), cov, parameters}
  end
  def chi2fit observables, {parameters, fun, penalties}, 0, nil, options do
    chi2 = chi2(observables, &(fun.(&1,parameters)), &(penalties.(&1,parameters)), options)
    alpha = alpha(observables, {parameters, fun, penalties, options})

    cov = try do
        alpha |> inverse
      catch
        {:impossible_inverse,error} ->
          throw {:inverse_error, error, chi2, parameters}
      rescue
        ArithmeticError ->
          throw {:inverse_error, ArithmeticError, chi2, parameters}
      end

    error = cov |> diagonal
    chi2fit observables, {parameters, fun, penalties}, 0, {cov,error}, options
  end
  def chi2fit observables, {parameters, fun, penalties}, max, preverror, options do
    matb = beta(observables, {parameters, fun, penalties})
    vecg = gamma(observables, {parameters, fun, penalties, options})
    chi2 = chi2(observables, &(fun.(&1,parameters)), &(penalties.(&1,parameters)),options)
    alpha = alpha(observables, {parameters, fun, penalties,options})

    try do
      cov = alpha |> inverse
      error = cov |> diagonal

      betainv = matb |> inverse
      delta = betainv |> Enum.map(&(dotproduct(&1,vecg)))
      {params,_chi2} = parameters
      |> vary_params
      |> Enum.reduce({parameters,chi2},
        fn
          (factor,{pars,oldchi}) ->
            dvec = factor |> from_diagonal |> Enum.map(&dotproduct(&1,delta))
            vec = ExAlgebra.Vector.add(dvec,parameters)
            try do
              newchi = chi2smooth observables,vec,{fun,penalties},options[:smoothing],options
              if newchi < oldchi do
                options[:onstep] && options[:onstep].(%{delta: dvec, chi2: newchi, params: vec})
                {vec,newchi}
              else
                {pars,oldchi}
              end
            rescue
              ArithmeticError ->
                Logger.debug "chi2fit: arithmetic error [#{inspect vec}] [#{inspect System.stacktrace}]"
                {pars,oldchi}
            end
        end)

      cond do
        Enum.all?(delta, &(&1 == 0)) -> 
          chi2fit observables, {params,fun,penalties}, 0, {cov,error}, options

        true ->
          chi2fit observables, {params,fun,penalties}, max-1, {cov,error}, options
      end
    catch
      {:impossible_inverse,error} ->
        Logger.debug "chi2: impossible inverse: #{error}"
        chi2fit observables, {parameters,fun,penalties}, 0, preverror, options
    rescue
      ArithmeticError ->
        Logger.debug "chi2: arithmetic error"
        chi2fit observables, {parameters,fun,penalties}, 0, preverror, options
    end

  end

end
