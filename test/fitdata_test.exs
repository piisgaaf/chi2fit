defmodule FitTest do

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
  import Chi2fit.Fit
  import Chi2fit.Utilities

  @moduletag distrib_mod: Chi2fit.Distribution
  @moduletag fit: true
  @moduletag p_value: 0.005

  ## See https://www.di-mgt.com.au/chisquare-calculator.html
  @chi2_table  %{
      {0.005,20} => 39.99684,
      {0.005,200} => 255.26415,
      {0.005,2000} => 2166.66430,
      {0.005,20000} => 20518.92189
    }

  doctest Chi2fit.Fit

  setup context do
    distrib = apply context[:distrib_mod], context[:distrib], context[:params]

    sample = 1..context[:size] |> Enum.map(fn _->distrib.() end)
    {cdf,bins,_,_} = get_cdf(sample,1,:wilson)
    {:ok, data: convert_cdf({cdf, bins|>Enum.map(&elem(&1,0))}), init: context[:params_init], cdf: context[:cdf]}
  end

  defp assert_fit expect_params, {data,init,cdf}, options do
    extra_options = options[:extra_options] || []
    fit_options = Keyword.merge([model: :linear],extra_options)

    {chi2,cov,params,_ranges} = chi2fit data, {init,fn x,pars -> apply(options[:distrib_mod],cdf,pars).(x) end}, 100, fit_options
    assert @chi2_table[{options.p_value,options[:size]}] && chi2 < @chi2_table[{options.p_value,50}]

    errors = cov |> Chi2fit.Matrix.diagonal
    Enum.zip([params,errors,expect_params])
    |> Enum.each(fn {par,err,check}->
      ssd = 2.0*:math.sqrt(err)
      assert par-ssd < check and check < par+ssd, "#{par-ssd} < #{check} < #{par+ssd}"
    end)
  end

  ##
  ## Weibull
  ##
  @tag distrib: :weibull
  @tag cdf: :weibullCDF
  @tag size: 20000
  @tag params: [0.4,5.0]
  @tag params_init: [0.2,1.0]
  test "fit - weibull (20000)", context do
    assert_fit context[:params], {context[:data],context[:init],context[:cdf]}, context
  end

  @tag distrib: :weibull
  @tag cdf: :weibullCDF
  @tag size: 2000
  @tag params: [0.4,5.0]
  @tag params_init: [0.2,1.0]
  test "fit - weibull (2000)", context do
    assert_fit context[:params], {context[:data],context[:init],context[:cdf]}, context
  end

  @tag distrib: :weibull
  @tag cdf: :weibullCDF
  @tag size: 200
  @tag params: [0.4,5.0]
  @tag params_init: [0.2,1.0]
  test "fit - weibull (200)", context do
    assert_fit context[:params], {context[:data],context[:init],context[:cdf]}, context
  end

  @tag distrib: :weibull
  @tag cdf: :weibullCDF
  @tag size: 20
  @tag params: [0.4,5.0]
  @tag params_init: [0.2,1.0]
  test "fit - weibull (20)", context do
    assert_fit context[:params], {context[:data],context[:init],context[:cdf]}, context
  end

  ##
  ## Normal
  ##
  @tag distrib: :normal
  @tag cdf: :normalCDF
  @tag size: 20000
  @tag params: [11.2,2.3]
  @tag params_init: [10.0,2.0]
  test "fit - normal (20000)", context do
    assert_fit context[:params], {context[:data],context[:init],context[:cdf]}, context
  end

  @tag distrib: :normal
  @tag cdf: :normalCDF
  @tag size: 2000
  @tag params: [11.2,2.3]
  @tag params_init: [10.0,2.0]
  test "fit - normal (2000)", context do
    assert_fit context[:params], {context[:data],context[:init],context[:cdf]}, context
  end

  @tag distrib: :normal
  @tag cdf: :normalCDF
  @tag size: 200
  @tag params: [11.2,2.3]
  @tag params_init: [10.0,2.0]
  test "fit - normal (200)", context do
    assert_fit context[:params], {context[:data],context[:init],context[:cdf]}, context
  end

  @tag distrib: :normal
  @tag cdf: :normalCDF
  @tag size: 20
  @tag params: [11.2,2.3]
  @tag params_init: [10.0,3.0]
  test "fit - normal (20)", context do
    assert_fit context[:params], {context[:data],context[:init],context[:cdf]}, context
  end

  ##
  ## Exponential
  ##
  @tag distrib: :exponential
  @tag cdf: :exponentialCDF
  @tag size: 20000
  @tag params: [1.7]
  @tag params_init: [1.0]
  test "fit - exponential (20000)", context do
    assert_fit context[:params], {context[:data],context[:init],context[:cdf]}, context
  end

  @tag distrib: :exponential
  @tag cdf: :exponentialCDF
  @tag size: 2000
  @tag params: [1.7]
  @tag params_init: [1.0]
  test "fit - exponential (2000)", context do
    assert_fit context[:params], {context[:data],context[:init],context[:cdf]}, context
  end

  @tag distrib: :exponential
  @tag cdf: :exponentialCDF
  @tag size: 200
  @tag params: [1.7]
  @tag params_init: [1.0]
  test "fit - exponential (200)", context do
    assert_fit context[:params], {context[:data],context[:init],context[:cdf]}, context
  end

  @tag distrib: :exponential
  @tag cdf: :exponentialCDF
  @tag size: 20
  @tag params: [1.7]
  @tag params_init: [1.0]
  test "fit - exponential (20)", context do
    assert_fit context[:params], {context[:data],context[:init],context[:cdf]}, context
  end

  ##
  ## Wald
  ##
  @tag distrib: :wald
  @tag cdf: :waldCDF
  @tag size: 20000
  @tag params: [1.7,6.3]
  @tag params_init: [1.0,1.0]
  test "fit - wald (20000)", context do
    assert_fit context[:params], {context[:data],context[:init],context[:cdf]}, context
  end

  @tag distrib: :wald
  @tag cdf: :waldCDF
  @tag size: 2000
  @tag params: [1.7,6.3]
  @tag params_init: [1.0,1.0]
  test "fit - wald (2000)", context do
    assert_fit context[:params], {context[:data],context[:init],context[:cdf]}, context
  end

  @tag distrib: :wald
  @tag cdf: :waldCDF
  @tag size: 200
  @tag params: [1.7,6.3]
  @tag params_init: [1.0,1.0]
  test "fit - wald (200)", context do
    assert_fit context[:params], {context[:data],context[:init],context[:cdf]}, context
  end

  @tag distrib: :wald
  @tag cdf: :waldCDF
  @tag size: 20
  @tag params: [1.7,6.3]
  @tag params_init: [1.0,6.0]
  test "fit - wald (20)", context do
    assert_fit context[:params], {context[:data],context[:init],context[:cdf]}, context
  end

  ##
  ## Erlang
  ##
  @tag distrib: :erlang
  @tag cdf: :erlangCDF
  @tag size: 20000
  @tag params: [5,1.3]
  @tag params_init: [4,1.0]
  test "fit - erlang (20000)", context do
    assert_fit context[:params], {context[:data],context[:init],context[:cdf]}, context
  end

  @tag distrib: :erlang
  @tag cdf: :erlangCDF
  @tag size: 2000
  @tag params: [5,1.3]
  @tag params_init: [4,1.0]
  test "fit - erlang (2000)", context do
    assert_fit context[:params], {context[:data],context[:init],context[:cdf]}, context
  end

  @tag distrib: :erlang
  @tag cdf: :erlangCDF
  @tag size: 200
  @tag params: [5,1.3]
  @tag params_init: [4,1.0]
  test "fit - erlang (200)", context do
    assert_fit context[:params], {context[:data],context[:init],context[:cdf]}, context
  end

  @tag distrib: :erlang
  @tag cdf: :erlangCDF
  @tag size: 20
  @tag params: [5,1.3]
  @tag params_init: [4,1.0]
  test "fit - erlang (20)", context do
    assert_fit context[:params], {context[:data],context[:init],context[:cdf]}, context
  end

  ##
  ## Various
  ##
  @tag distrib: :normal
  @tag params: [1,2]
  @tag size: 50
  @tag various: :normal50
  @tag :notest
  test "normal (50)",context do
    ## Normal
    distrib = apply context[:distrib_mod], :normal, [33.0,2.4]
    sample = 1..50 |> Enum.map(fn _->distrib.() end)
    {cdf,bins,_,_} = get_cdf(sample,1,:wilson)
    data = convert_cdf({cdf, bins|>Enum.map(&elem(&1,0))})
    {chi2,cov,params,_ranges} = chi2fit data, {[1.0,30.0],fn x,pars -> apply(context[:distrib_mod],:normalCDF,pars).(x) end}, 100, [model: :linear]
    [p1,p2] = params
    [[v1,_],[_,v2]] = cov

    IO.puts "\n\n============================================"
    IO.puts "distrib = normal; mu=33.0, sigma=2.4\n"
    IO.puts "chi2 = #{chi2}"
    IO.puts "dof = 48"
    IO.puts "(mu,   mu_var)    = #{p1},#{v1}"
    IO.puts "(sigma,sigma_var) = #{p2},#{v2}"
  end
  
  @tag distrib: :weibull
  @tag params: [1,2]
  @tag size: 50
  @tag various: :weibull_l50
  @tag :notest
  test "weibull [lean] (50)",context do
    ## Weibull: mode<1.0
    distrib = apply context[:distrib_mod], :weibull, [0.42,33.0]
    sample = 1..50 |> Enum.map(fn _->distrib.() end)
    {cdf,bins,_,_} = get_cdf(sample,1,:wilson)
    data = convert_cdf({cdf, bins|>Enum.map(&elem(&1,0))})
    {chi2,cov,params,_ranges} = chi2fit data, {[1.0,30.0],fn x,pars -> apply(context[:distrib_mod],:weibullCDF,pars).(x) end}, 100, [model: :linear]
    [p1,p2] = params
    [[v1,_],[_,v2]] = cov

    IO.puts "\n\n============================================"
    IO.puts "distrib = weibull; alpha=0.42, beta=33.0\n"
    IO.puts "chi2 = #{chi2}"
    IO.puts "dof = 48"
    IO.puts "(alpha,alpha_var) = #{p1},#{v1}"
    IO.puts "(beta, beta_var)  = #{p2},#{v2}"
  end
  
  @tag distrib: :weibull
  @tag params: [1,2]
  @tag size: 50
  @tag various: :weibull_a50
  @tag :notest
  test "weibull [agile] (50)",context do
    ## Weibull: 1.0<mode<2.0
    distrib = apply context[:distrib_mod], :weibull, [1.42,33.0]
    sample = 1..50 |> Enum.map(fn _->distrib.() end)
    {cdf,bins,_,_} = get_cdf(sample,1,:wilson)
    data = convert_cdf({cdf, bins|>Enum.map(&elem(&1,0))})
    {chi2,cov,params,_ranges} = chi2fit data, {[1.0,30.0],fn x,pars -> apply(context[:distrib_mod],:weibullCDF,pars).(x) end}, 100, [model: :linear]
    [p1,p2] = params
    [[v1,_],[_,v2]] = cov

    IO.puts "\n\n============================================"
    IO.puts "distrib = weibull; alpha=1.42, beta=33.0\n"
    IO.puts "chi2 = #{chi2}"
    IO.puts "dof = 48"
    IO.puts "(alpha,alpha_var) = #{p1},#{v1}"
    IO.puts "(beta, beta_var)  = #{p2},#{v2}"
  end
  
  @tag distrib: :weibull
  @tag params: [1,2]
  @tag size: 50
  @tag various: :weibull_w50
  @tag :notest
  test "weibull [waterfall] (50)",context do
    ## Weibull: mode>2.0
    distrib = apply context[:distrib_mod], :weibull, [2.42,33.0]
    sample = 1..50 |> Enum.map(fn _->distrib.() end)
    {cdf,bins,_,_} = get_cdf(sample,1,:wilson)
    data = convert_cdf({cdf, bins|>Enum.map(&elem(&1,0))})
    {chi2,cov,params,_ranges} = chi2fit data, {[1.0,30.0],fn x,pars -> apply(context[:distrib_mod],:weibullCDF,pars).(x) end}, 100, [model: :linear]
    [p1,p2] = params
    [[v1,_],[_,v2]] = cov

    IO.puts "\n\n============================================"
    IO.puts "distrib = weibull; alpha=2.42, beta=33.0\n"
    IO.puts "chi2 = #{chi2}"
    IO.puts "dof = 48"
    IO.puts "(alpha,alpha_var) = #{p1},#{v1}"
    IO.puts "(beta, beta_var)  = #{p2},#{v2}"
  end

end
