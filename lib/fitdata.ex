defmodule Chi2fit.Fit do

  # Copyright 2012-2017 Pieter Rijken
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
  Implements fitting a distribution function to sample data. It minimizes the liklihood function.
  
  ## Asymmetric Errors
  To handle asymmetric errors the module provides three ways of determining the contribution to the likelihood function:
      `simple` - value difference of the observable and model divided by the averaged error lower and upper bounds;
      `asimple` - value difference of the observable and model divided by the difference between upper/lower bound and the observed
        value depending on whether the model is larger or smaller than the observed value;
      `linear` - value difference of the observable and model divided by a linear tranformation (See below).
      
  ### 'linear': Linear transformation
  Linear transformation that:
      - is continuous in u=0,
      - passes through the point sigma+ at u=1,
      - asymptotically reaches 1-y at u->infinity
      - pass through the point -sigma- at u=-1,
      - asymptotically reaches -y at u->-infinity

  ## References
  [1] See https://arxiv.org/pdf/physics/0401042v1.pdf
  """

  require Logger

  import Chi2fit.Matrix
  import Chi2fit.Utilities

  @typedoc "Observation with symmetric errors 'dy'."
  @type observable_symm ::  {x :: float, y :: float, dy :: float}

  @typedoc "Observation with asymmetric bounds 'y1 < y < y2'."
  @type observable_asym ::  {x :: float, y :: float, y1 :: float, y2 :: float}

  @type observable :: observable_symm | observable_asym
  @type observables :: [observable]

  @typedoc "Cumulative distribution mapping 'x' and parameters to a float in the range [0,1]."
  @type distribution :: ((x::float,[parameter::float])->float)

  @typedoc "Tuple describing the parameter values and the distribution function."
  @type model :: {[float], distribution}

  @typedoc "Chi-squared statistic"
  @type chi2 :: float

  @typedoc "Covariance matrix"
  @type cov :: Chi2fit.Matrix.matrix

  @typedoc "List of parameter ranges"
  @type params :: [{float,float}]

  @arithmic_penalty 1_000_000_000

  defp nopenalties(_,_), do: 0.0

  defp dchi2_simple(y, y1, y2,f),            do: (f-y)/abs(y-(y1+y2)/2)
  defp dchi2_asimple(y, y1,_y2,f) when f<y,  do: (y-f)/(y-y1)
  defp dchi2_asimple(y,_y1,_y2,f) when y==f, do: 0.0
  defp dchi2_asimple(y,_y1, y2,f),           do: (f-y)/(y2-y)
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
      delta>0 -> (1.0-y2)/(1.0-f) * delta/splus
      true ->          y1/f       * delta/smin
    end
  end

  defp likelihood_contrib(:linear, y,y1,y2,f), do: dchi2_linear y,y1,y2,f
  defp likelihood_contrib(:simple, y,y1,y2,f), do: dchi2_simple y,y1,y2,f
  defp likelihood_contrib(:asimple, y,y1,y2,f), do: dchi2_asimple y,y1,y2,f

  @doc """
  Calculates the Chi-squared function for a list of observables.

  The `observables` are given as a list. Each observation has an error associated with it. The errors can be either
  symmetric or asymmetric.

  A 'penalties'-function is used to assign penalties and these contribute to the chi-squared function. It may be used
  to 'forbid' certain parameter, x combinations.

  ## Options

      `model` - Required. Determines the contribution to chi-squared taking the asymmetric errors into account.
              Vaid values are `:linear`, `:simple`, and `:asimple`. See Errors below

  ## Errors
      `simple` - Use for asymmetric errors when the sigma+ and sigma- are close to each other
      `asimple` - Use for asymmetric errors when y-values are not bound.
      `linear` - Use this model when the y-values ar bound between 0 and 1. Linear transformation that:
          - is continuous in u=0,
          - passes through the point sigma+ at u=1,
          - asymptotically reaches 1-y at u->infinity
          - pass through the point -sigma- at u=-1,
          - asymptotically reaches -y at u->-infinity
  """
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
          try do
            tmp = likelihood_contrib options[:model], y,y1,y2,fun.(x)
            tmp*tmp + penalties.(x)
          rescue
            ArithmeticError -> @arithmic_penalty
          end
      end)
    |> Enum.sum
  end

  defp gamma(observables, {parameters, fun, penalties, options}) do
    gammafun = &(gamma(&1,observables, {parameters, fun, penalties, options}))
    Enum.reduce(length(parameters)..1, [], fn (k,acc)->[gammafun.(k)|acc] end)
  end

  @spec gamma(pos_integer, observables, model) :: float
  defp gamma(k, observables, {parameters, fun, penalties, options}) when k>0 and k<=length(parameters) do
    params_k = parameters |> derive_par(k)
    -0.5*der params_k, fn (pars)->chi2smooth(observables, pars, {fun,penalties},options[:smoothing],options) end, options
  end

  defp alpha(observables, {parameters, fun, penalties, options}) do
    alphafun = &(alpha({&1,&2}, observables, {parameters, fun, penalties,options}))
    Enum.reduce(length(parameters)..1, [], fn
      (k,acc) -> [
        Enum.reduce(length(parameters)..1, [], fn
          (j,acc)->[alphafun.(k,j)|acc] end)
        |acc]
      end)
  end

  defp derive_par(list, index) do
    list |> List.update_at(index-1, fn
        (val) when is_number(val) ->
          {val,1}

        ({val,n}) ->
          {val,n+1}
      end)
  end

  @spec alpha({pos_integer,pos_integer}, observables, model) :: float
  defp alpha({k,j}, observables, {parameters, fun, penalties, options}) when k>0 and k<=length(parameters) and j>0 and j<=length(parameters) do
    params_kj = parameters |> derive_par(k) |> derive_par(j)
    0.5*der params_kj,fn (pars)->chi2smooth(observables, pars, {fun,penalties},options[:smoothing],options) end, options
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

  @doc """
  Probes the chi-squared surface within a certain range of the parameters.
  
  It does so by randomly selecting parameter value combinations and calculate the chi-squared for the list
  of observations based on the selected parameter values. This routine is used to roughly probe the chi-squared
  surface and perform more detailed and expensive calculations to precisely determine the minimum by `chi2fit/5`.
  
  Returns the minimum chi-squared found, the parameter values, and all probes that resulted in chi-squared difference
  less than 1 with the minimum. The parameter values found in this set correspond with the errors in determining
  the parameters.
  
  ## Options
      `num` or `probes` - the number of points to calculate,
      `mark` - progress indicator: a keyword list with keys `m`, `c`, `x`, and `*`; the value must be a call back
      function taking zero arguments. These are called when 1000, 100, 10, probes have been done. The value of
      key `*` is called when a new chi-squared minimum has been found,
      `smoothing` - boolean value indicating whether the chi-squared is smoothened using a Gauss distribution. This
      is used in case the surface is rough because of numerical instabilities to smoothen the surface,
      `model` - See chi2/3 and chi2/4

  """
  @spec chi2probe(observables, [float], (...->any), Keyword.t) :: {chi2::float,[parameters::float],{[float],[float]}}
  def chi2probe(observables, parranges, fun_penalties, options) do
    chi2probe(observables, parranges, fun_penalties, options[:num] || options[:probes], nil, options)
  end

  defp chi2probe(_observables, _parranges, {_fun,_penalties}, 0, best, _options) do
    ## Refactor this!!!!!
    {chi2,parameters,saved} = best
    {_chis,plists} = saved |> Enum.unzip
    {chi2,parameters,
      plists
      |> Enum.map(&List.to_tuple/1)
      |> unzip
      |> Tuple.to_list
      |> Enum.map(fn plist -> [ Enum.min(plist),Enum.max(plist) ] end)
      |> List.to_tuple}
  end
  defp chi2probe(observables, parranges, {fun,penalties}, num, best, options) do
    if options[:progress] do
      cond do
        rem(num,1000) == 0 -> options[:mark][:m].()
        rem(num,100) == 0 -> options[:mark][:c].()
        rem(num,10) == 0 -> options[:mark][:x].()
        true -> :ok
      end
    end
    try do
      parameters = parranges |> sample
      chi2 = chi2smooth observables,parameters,{fun,penalties},options[:smoothing],options
      chi2probe(observables, parranges, {fun,penalties}, num-1,
        case best do
          nil ->
            {chi2,parameters,[{chi2,parameters}]}
          {oldchi2,_,saved} when chi2<oldchi2 ->
            options[:mark][:*].()
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
        reraise err, "Error!"
    end
  end

  defp vary_params(parameters, num_variations \\ 100) when is_list(parameters) do
    -1..length(parameters)
    |> Stream.map(&(List.duplicate(&1,num_variations)))
    |> Stream.concat
    |> Stream.flat_map(
      fn
        (-1) -> [List.duplicate(:rand.uniform(),length(parameters)), List.duplicate(:rand.uniform()/100,length(parameters))]
        (0) -> [List.duplicate(0.0,length(parameters)) |> Enum.map(fn (_)->:rand.uniform() end)]
        (n) when is_integer(n) and n>0 -> [List.duplicate(0.0,length(parameters)) |> List.replace_at(n-1, :rand.uniform()),List.duplicate(0.0,length(parameters)) |> List.replace_at(n-1, :rand.uniform()/100)]
      end)
  end

  @doc """
  Fits observables to a known model.
  
  Returns the found minimum chi-squared value, covariance matrix, gradient at the minimum, and the corresponding parameter values including
  error estimates.
  For a good fit check the following:
      `chi2 per degree of freedom` - this should be about 1 or less,
      `gradient` - at the minimum the gradient should be zero at all directions.
  
  For asymmetric errors use the option `model` equal to `linear`.
  Rough chi-squared surfaces or if numerically unstable, use the option `smoothing` set to `true`.

  ## Arguments
      `observables` - list of measurements including errors,
      `model` - `{parameters, fun}`: set of initial parameter values and a function to fit against the measurements

  ## Options
      `onstep` - call back function; it is called with a map with keys `delta`, `chi2`, and `params`,
      `smoothing` - boolean value indicating whether the chi-squared is smoothened using a Gauss distribution. This
      is used in case the surface is rough because of numerical instabilities to smoothen the surface,
      `model` - The same values as in chi2/3 and chi2/4
  """
  @spec chi2fit(observables, model, iterations::pos_integer, options::Keyword.t) :: {chi2,cov,params}
  def chi2fit(observables, model, max \\ 100, error \\ nil, options \\ [])
  def chi2fit(observables, {parameters, fun}, max, error, options), do: chi2fit observables, {parameters, fun, &nopenalties/2}, max, error, options
  def chi2fit(observables, {parameters, fun, penalties}, 0, {cov,_error}, options) do
    {chi2(observables, &(fun.(&1,parameters)), &(penalties.(&1,parameters)), options), cov, parameters}
  end
  def chi2fit observables, {parameters, fun, penalties}, 0, nil, options do
    chi2 = chi2(observables, &(fun.(&1,parameters)), &(penalties.(&1,parameters)), options)
    alpha = alpha(observables, {parameters, fun, penalties, options})

    {:ok,cov} = try do
        IO.inspect alpha
        alpha |> inverse
      catch
        {:impossible_inverse,error} ->
          throw {:inverse_error, error, chi2, parameters}

        {:failed_to_reach_tolerance,_pars,error} ->
          throw {:failed_to_reach_tolerance, error, chi2, parameters}

      rescue
        ArithmeticError ->
          throw {:inverse_error, ArithmeticError, chi2, parameters}
      end

    error = cov |> diagonal
    chi2fit observables, {parameters, fun, penalties}, 0, {cov,error}, options
  end
  def chi2fit observables, {parameters, fun, penalties}, max, preverror, options do
    vecg = gamma(observables, {parameters, fun, penalties, options})
    chi2 = chi2(observables, &(fun.(&1,parameters)), &(penalties.(&1,parameters)),options)
    alpha = alpha(observables, {parameters, fun, penalties,options})

    try do
      {:ok,cov} = alpha |> inverse
      error = cov |> diagonal

      delta = cov |> Enum.map(&(dotproduct(&1,vecg)))
      {params,_chi2} = parameters
      |> vary_params
      |> Enum.reduce({parameters,chi2},
        fn
          (factor,{pars,oldchi}) ->
            dvec = factor |> from_diagonal |> Enum.map(&dotproduct(&1,delta))
            vec = ExAlgebra.Vector.add(parameters,dvec)
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
      {:failed_to_reach_tolerance,_pars,error} ->
        Logger.debug "chi2: failed to reach tolerance: #{error}"
        chi2fit observables, {parameters,fun,penalties}, 0, preverror, options
    rescue
      ArithmeticError ->
        Logger.debug "chi2: arithmetic error"
        chi2fit observables, {parameters,fun,penalties}, 0, preverror, options
    end

  end

end
