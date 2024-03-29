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

  alias Chi2fit.Distribution, as: D
  alias Chi2fit.Matrix, as: M
  alias Chi2fit.Math, as: Ma
  alias Chi2fit.Statistics, as: S
  alias Chi2fit.Utilities, as: U

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
  @type cov :: M.matrix()

  @typedoc "List of parameter ranges"
  @type params :: [{float,float}]

  @typedoc "Tuple with chi-squared, parameter values, and parameter errors at the found minimum (see `chi2probe/4`)"
  @type chi2probe_simple :: {chi2(),[parameters::float],{[float],[float]}}
  
  @typedoc "Tuple with chi-squared, parameter values, parameter errors, and list of intermediate fit results (see `chi2probe/4`)"
  @type chi2probe_saved :: {chi2(),[parameters::float],{[float],[float]},[{float,[float]}]}
  
  @typedoc "Result of chi-squared probe (see &chi2probe/4)"
  @type chi2probe :: chi2probe_simple() | chi2probe_saved()
  
  @typedoc "Tuple holding chi-squared value, covariance matrix, parameter values, and parameter errors at the minimum chi2fit(see `chi2fit/4`)"
  @type chi2fit :: {chi2(), cov(), parameters :: [float], errors :: [float]}
  
  @arithmic_penalty 1_000_000_000
  @extreme_punishment 1_000_000

  def nopenalties(_,_), do: 0.0

  @cutoff 0.0001
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
      #
      f<@cutoff and f<y1 -> dchi2_linear(y,y1,y2,min(@cutoff,y1))
      f>(1-@cutoff) and f>y2 -> dchi2_linear(y,y1,y2,max(1-@cutoff,y2))
      # Extreme punishment
      f==1.0 -> @extreme_punishment
      f==0.0 -> @extreme_punishment
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

  ## Examples
  
      iex> fun = &(&1)
      ...> chi2 [{0,3,1}], fun
      9.0

      iex> fun = &(&1)
      ...> chi2 [{0,3,1},{1,7,1},{2,3,1}], fun
      46.0

      iex> fun = &(&1)
      ...> chi2 [{0,3,3},{1,7,1},{2,3,1}], fun
      38.0

      iex> fun = &(&1-2)
      ...> chi2 [{0,3,1}], fun
      25.0

  end

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
            ArgumentError -> @arithmic_penalty
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
    smoothing? = options[:smoothing] || false
    params_k = parameters |> derive_par(k)
    -0.5*Ma.der params_k, fn (pars)->chi2smooth(observables, pars, {fun,penalties},smoothing?,options) end, options
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
    smoothing? = options[:smoothing] || false
    params_kj = parameters |> derive_par(k) |> derive_par(j)
    0.5*Ma.der params_kj,fn (pars)->chi2smooth(observables, pars, {fun,penalties},smoothing?,options) end, options
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
  surface and perform more detailed and expensive calculations to precisely determine the minimum by `chi2fit/4`.
  
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
      `model` - See `chi2/3` and `chi2/4`

  """
  @spec chi2probe(observables, [float], (...->any), Keyword.t) :: chi2probe()
  def chi2probe(observables, parranges, fun_penalties, options) do
    chi2probe(observables, parranges, fun_penalties, options[:num] || options[:probes], nil, options)
  end

  defp chi2probe(_observables, _parranges, {_fun,_penalties}, 0, best, options) do
    ## Refactor this!!!!!
    {chi2,parameters,saved} = best
    {_chis,plists} = saved |> Enum.unzip
    result = {chi2,parameters,
      plists
      |> Enum.map(&List.to_tuple/1)
      |> U.unzip
      |> Tuple.to_list
      |> Enum.map(fn plist -> [ Enum.min(plist),Enum.max(plist) ] end)
      |> List.to_tuple}
    if options[:saved?], do: Tuple.append(result,saved), else: result
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
      smoothing = options[:smoothing] || false
      parameters = parranges |> sample
      chi2 = chi2smooth observables,parameters,{fun,penalties},smoothing,options
      if options[:surfacefile], do: IO.puts options[:surfacefile], "#{Enum.join(parameters,",")},#{chi2}"
      if options[:save], do: options[:save].(parameters,chi2)
      chi2probe(observables, parranges, {fun,penalties}, num-1,
        case best do
          nil ->
            {chi2,parameters,[{chi2,parameters}]}
          {oldchi2,_,saved} when chi2<oldchi2 ->
            is_function(options[:mark][:*]) && options[:mark][:*].()
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
        reraise err, __STACKTRACE__
    end
  end

  defp filter_param(parameters, flags) do
    parameters |> Enum.zip(flags) |> Enum.map(fn {par, true}->par; {_par, false}->0.0 end)
  end

  defp vary_params(parameters, flags, _first, _second, num_variations) when is_list(parameters) do
    -1..length(parameters) |> Enum.to_list
    |> Stream.map(&(List.duplicate(&1,num_variations)))
    |> Stream.concat
    |> Stream.flat_map(
      fn
        (-1) ->
          [
            List.duplicate(:rand.uniform(),length(parameters))|>filter_param(flags),
            List.duplicate(:rand.uniform()/100,length(parameters))|>filter_param(flags)
          ]
        (0) ->
          [
            List.duplicate(0.0,length(parameters)) |> Enum.map(fn (_)-> 2* :rand.uniform() end)|>filter_param(flags),
            List.duplicate(0.0,length(parameters)) |> Enum.map(fn (_)->    :rand.uniform() end)|>filter_param(flags),
            List.duplicate(0.0,length(parameters)) |> Enum.map(fn (_)->  - :rand.uniform() end)|>filter_param(flags)
          ]
        (n) when is_integer(n) and n>0 ->
          [
            List.duplicate(0.0,length(parameters)) |> List.replace_at(n-1, :rand.uniform())|>filter_param(flags),
            List.duplicate(0.0,length(parameters)) |> List.replace_at(n-1, :rand.uniform()/100)|>filter_param(flags),
            List.duplicate(0.0,length(parameters)) |> List.replace_at(n-1, :rand.uniform()/1000)|>filter_param(flags),
            List.duplicate(0.0,length(parameters)) |> List.replace_at(n-1, :rand.uniform()/10000)|>filter_param(flags)
          ]
      end)
    |> Stream.filter(fn list -> Enum.any?(list, fn x->x != 0.0 end) end)
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
      `model` - The same values as in `chi2/3` and `chi2/4`
      `grid?` - Performs a grid search: per step, tries to fit only one parameter and keeps the others fixed; selects the parameter in
      a round-robin fashion
      `probes` -- a list of tuples containing the result of the `chi2probe/4` function. Each tuple contains the chi2 value and parameter list.
      Defaults to the empty list.
  """
  @spec chi2fit(observables, model, iterations::pos_integer, options::Keyword.t) :: chi2fit()
  def chi2fit(observables, model, max \\ 100, options \\ []) do
    probes = options[:probes] || []
    chi2fit observables, model, max, {nil,probes}, options
  end

  defp chi2fit(observables, {parameters, fun}, max, data, options), do: chi2fit observables, {parameters, fun, &nopenalties/2}, max, data, options
  defp chi2fit(observables, {parameters, fun, penalties}, 0, {{cov,_error}, params}, options) do
    chi = chi2(observables, &(fun.(&1,parameters)), &(penalties.(&1,parameters)), options)
    ranges = params
    |> Enum.flat_map(fn {c,p}-> if(c < chi+1.0, do: [List.to_tuple([c|p])], else: []) end)
    |> Chi2fit.Utilities.unzip
    |> Tuple.to_list
    |> Enum.map(fn x->[Enum.min(x),Enum.max(x)] end)
    |> List.to_tuple

    {chi, cov, parameters, ranges}
  end
  defp chi2fit observables, {parameters, fun, penalties}, 0, {nil,ranges}, options do
    chi2 = chi2(observables, &(fun.(&1,parameters)), &(penalties.(&1,parameters)), options)
    alpha = alpha(observables, {parameters, fun, penalties, options})

    {:ok,cov} = try do
        alpha |> M.inverse
      rescue
        ArithmeticError ->
          throw {:inverse_error, ArithmeticError, chi2, parameters}
      catch
        {:impossible_inverse,error} ->
          throw {:inverse_error, error, chi2, parameters}

        {:failed_to_reach_tolerance,_pars,error} ->
          throw {:failed_to_reach_tolerance, error, chi2, parameters}

      end

    error = cov |> M.diagonal
    chi2fit observables, {parameters, fun, penalties}, 0, {{cov,error},[{chi2,parameters}|ranges]}, options
  end
  defp chi2fit(observables, {parameters, fun, penalties}, max, {preverror,ranges}, options) when max>0 do
    grid? = options[:grid?] || false
    parmax = options[:pariter] || 10

    vecg = gamma(observables, {parameters, fun, penalties, options})
    chi2 = chi2(observables, &(fun.(&1,parameters)), &(penalties.(&1,parameters)),options)
    alpha = alpha(observables, {parameters, fun, penalties,options})

    ranges = [{chi2,parameters}|ranges]

    try do
      {:ok,cov} = alpha |> M.inverse
      error = cov |> M.diagonal

      flags = if grid? do
        List.duplicate(false, length(parameters))
        |> List.replace_at(rem(max,length(parameters)),true) ## grid search
      else
        List.duplicate(true, length(parameters))
      end

      delta = cov |> Enum.map(&(M.dotproduct(&1,vecg)))

      options[:onbegin] && options[:onbegin].(%{step: max, chi: chi2, derivatives: Enum.zip(vecg,M.diagonal(alpha))})

      {params,_chi2,_,ranges} = parameters
      |> vary_params(flags,vecg,M.diagonal(alpha),parmax)
      |> Enum.reduce_while({parameters,chi2,{nil,nil},ranges},
        fn
          (factor,{pars,oldchi,minimum,data}) ->
            dvec = factor |> M.from_diagonal |> Enum.map(&M.dotproduct(&1,delta))
            vec = M.add(parameters,dvec)
            try do
              smoothing? = options[:smoothing] || false
              newchi = chi2smooth observables,vec,{fun,penalties},smoothing?,options
              data = [{newchi,vec}|data]

              new_minimum = case minimum do
                {nil,nil} ->
                  current_par = pars |> filter_param(flags) |> Enum.sum
                  par = vec |> filter_param(flags) |> Enum.sum
                  cond do
                    newchi < oldchi and par <= current_par ->
                      {nil,{pars,oldchi}}

                    newchi < oldchi and par > current_par ->
                      {{pars,oldchi},nil}

                    newchi >= oldchi and par <= current_par ->
                      {{vec,newchi},nil}

                    newchi >= oldchi and par > current_par ->
                      {nil,{vec,newchi}}
                  end

                {left={left_pars,_},nil} ->
                  current_par = pars |> filter_param(flags) |> Enum.sum
                  left_par = left_pars |> filter_param(flags) |> Enum.sum
                  par = vec |> filter_param(flags) |> Enum.sum
                  cond do
                    newchi < oldchi and par < current_par and par > left_par ->
                      {left,{pars,oldchi}}

                    newchi < oldchi and par > current_par ->
                      {{pars,oldchi},nil}

                    newchi > oldchi and par < current_par and par < left_par ->
                      {{vec,newchi},nil}

                    newchi > oldchi and par > current_par ->
                      {left,{vec,newchi}}

                    true ->
                      {left,nil}
                  end

                {nil,right={right_pars,_}} ->
                  current_par = pars |> filter_param(flags) |> Enum.sum
                  right_par = right_pars |> filter_param(flags) |> Enum.sum
                  par = vec |> filter_param(flags) |> Enum.sum
                  cond do
                    newchi < oldchi and par > current_par and par < right_par ->
                      {{pars,oldchi},right}

                    newchi < oldchi and par < current_par ->
                      {nil,{pars,oldchi}}

                    newchi > oldchi and par > current_par and par < right_par ->
                      {nil,{vec,newchi}}

                    newchi > oldchi and par < current_par ->
                      {{vec,newchi},right}

                    true ->
                      {nil,right}
                  end

                local -> local
              end

              case new_minimum do
                {{left_pars,_lchi2},{right_pars,_rchi2}} ->
                  left_par = left_pars |> filter_param(flags) |> Enum.sum
                  right_par = right_pars |> filter_param(flags) |> Enum.sum

                  try do
                    smoothing? = options[:smoothing] || false
                    {_, {left,right},{_vleft,_vright}} = Ma.newton(left_par, right_par, fn x->
                      nvec = parameters |> List.replace_at(rem(max,length(parameters)), {x,1})
                      Ma.der(nvec, fn lp -> chi2smooth(observables,lp,{fun,penalties},smoothing?,options) end,options)
                    end, parmax, options)

                    nvec = pars |> List.replace_at(rem(max,length(parameters)), (left+right)/2)
                    newchi = chi2smooth observables,nvec,{fun,penalties},smoothing?,options
                    data = [{newchi,nvec}|data]

                    if newchi<oldchi do
                      options[:onstep] && options[:onstep].(%{delta: dvec, chi2: newchi, params: vec})
                      {:halt, {nvec,newchi,:done,data}}
                    else
                      {:cont, {pars,oldchi,minimum,data}}
                    end
                  rescue
                    ArgumentError -> {:cont, {pars,oldchi,minimum,data}}
                    ArithmeticError -> {:cont, {pars,oldchi,minimum,data}}
                  end

                _otherwise ->
                  cond do
                    newchi < oldchi ->
                      options[:onstep] && options[:onstep].(%{delta: dvec, chi2: newchi, params: vec})
                      {:cont, {vec,newchi,new_minimum,data}}
                    true ->
                      {:cont, {pars,oldchi,minimum,data}}
                  end
                end
            rescue
              ArithmeticError ->
                Logger.warn "chi2fit: arithmetic error [#{inspect vec}] [#{inspect __STACKTRACE__}]"
                {:cont, {pars,oldchi,minimum,data}}
            end
        end)

      cond do
        Enum.all?(delta, &(&1 == 0)) -> 
          chi2fit observables, {params,fun,penalties}, 0, {{cov,error},ranges}, options

        true ->
          chi2fit observables, {params,fun,penalties}, max-1, {{cov,error},ranges}, options
      end
    rescue
      ArithmeticError ->
        Logger.warn "chi2: arithmetic error"
        IO.puts "#{inspect __STACKTRACE__}"
        chi2fit observables, {parameters,fun,penalties}, 0, {preverror,ranges}, options
    catch
      {:impossible_inverse,error} ->
        Logger.warn "chi2: impossible inverse: #{error}"
        chi2fit observables, {parameters,fun,penalties}, 0, {preverror,ranges}, options
      {:failed_to_reach_tolerance,_pars,error} ->
        Logger.warn "chi2: failed to reach tolerance: #{error}"
        chi2fit observables, {parameters,fun,penalties}, 0, {preverror,ranges}, options
    end

  end

  defp _find_change(enumerable, fun) when is_function(enumerable,2), do: enumerable |> U.subsequences |> Stream.map(fun)
  defp _find_change(enumerable, fun), do: enumerable |> U.subsequences |> Enum.map(fun)

  defp probe_seq(seq, binsize, initial, cdf, options) do
    seq
    |> S.to_bins({binsize,0})
    |> chi2probe(initial, {cdf, &nopenalties/2}, options)
    |> Tuple.to_list
    |> Enum.take(2)
    |> List.to_tuple
  end
  
  @doc """
  Finds the point in the data where the chi-squared has a jump when fitting the model
  """
  @spec find_change(list :: [number()],options :: Keyword.t) :: [{chi :: float,[float]}]
  def find_change(list, options) do
    binsize = options[:bin] || 1
    initial = options[:init]
    model = options[:fitmodel]
    cdf = D.cdf(model)

    list |> _find_change(& probe_seq(&1,binsize,initial,cdf,options))
  end

  @doc """
  Partitions the data list in segments with similar chi-squared values when fitting the model
  """
  @spec find_all(nil | [number],options :: Keyword.t) :: [[float()]]
  def find_all(data,options), do: find_all(data,[],options)
  
  defp find_all(data, acc, options)
  defp find_all(nil, acc, _options), do: Enum.reverse(acc)
  defp find_all(data, _acc, options) do
      threshold = options[:threshold] || 10
      binsize = options[:bin] || 1
      initial = options[:init]
      model = options[:fitmodel]
      tolerance = options[:tolerance] || 10
      cdf = D.cdf(model)

      data
      |> Stream.chunk_while({0,nil,[]},
        fn
          dat, {0,nil,[]} ->
            {chi2, params} = probe_seq([dat],binsize,initial,cdf,options)
            {:cont, {chi2,params,[dat]}}

          dat, {lastchi2, lastparams, saved} ->
            {chi2, params} = probe_seq([dat|saved],binsize,initial,cdf,options)
            if chi2<tolerance*lastchi2 or chi2<threshold do
              {:cont, {chi2, params, [dat|saved]}}
            else
              {newchi2, newparams} = probe_seq([dat],binsize,initial,cdf,options)
              {:cont, {lastchi2, lastparams, Enum.reverse(saved)}, {newchi2, newparams, [dat]}}
            end
        end,
        fn
          {chi,params,saved} -> {:cont,{chi,params,Enum.reverse(saved)},[]}
        end)
      |> (& if is_function(data, 2), do: &1, else: Enum.into(&1, [])).()
  end

end
