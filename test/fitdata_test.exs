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

  doctest Chi2fit.Fit

  setup context do
    distrib = apply context[:distrib_mod], context[:distrib], context[:params]

    sample = 1..context[:size] |> Enum.map(fn _->distrib.() end)
    {cdf,bins,_,_} = get_cdf(sample,1,:wilson)
    {:ok, data: convert_cdf({cdf, bins|>Enum.map(&elem(&1,0))}), init: context[:params_init], cdf: context[:cdf]}
  end

  defp assert_fit {expect_chi2,expect_params}, {data,init,cdf}, options do
    extra_options = options[:extra_options] || []
    fit_options = Keyword.merge([model: :linear],extra_options)

    {chi2,cov,params,_ranges} = chi2fit data, {init,fn x,pars -> apply(options[:distrib_mod],cdf,pars).(x) end}, 100, fit_options
    assert chi2 < expect_chi2

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
    assert_fit {20.0,context[:params]}, {context[:data],context[:init],context[:cdf]}, context
  end

  @tag distrib: :weibull
  @tag cdf: :weibullCDF
  @tag size: 2000
  @tag params: [0.4,5.0]
  @tag params_init: [0.2,1.0]
  test "fit - weibull (2000)", context do
    assert_fit {20.0,context[:params]}, {context[:data],context[:init],context[:cdf]}, context
  end

  @tag distrib: :weibull
  @tag cdf: :weibullCDF
  @tag size: 200
  @tag params: [0.4,5.0]
  @tag params_init: [0.2,1.0]
  test "fit - weibull (200)", context do
    assert_fit {2.0,context[:params]}, {context[:data],context[:init],context[:cdf]}, context
  end

  @tag distrib: :weibull
  @tag cdf: :weibullCDF
  @tag size: 20
  @tag params: [0.4,5.0]
  @tag params_init: [0.2,1.0]
  test "fit - weibull (20)", context do
    assert_fit {2.0,context[:params]}, {context[:data],context[:init],context[:cdf]}, context
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
    assert_fit {20.0,context[:params]}, {context[:data],context[:init],context[:cdf]}, context
  end

  @tag distrib: :normal
  @tag cdf: :normalCDF
  @tag size: 2000
  @tag params: [11.2,2.3]
  @tag params_init: [10.0,2.0]
  test "fit - normal (2000)", context do
    assert_fit {20.0,context[:params]}, {context[:data],context[:init],context[:cdf]}, context
  end

  @tag distrib: :normal
  @tag cdf: :normalCDF
  @tag size: 200
  @tag params: [11.2,2.3]
  @tag params_init: [10.0,2.0]
  test "fit - normal (200)", context do
    assert_fit {20.0,context[:params]}, {context[:data],context[:init],context[:cdf]}, context
  end

  @tag distrib: :normal
  @tag cdf: :normalCDF
  @tag size: 20
  @tag params: [11.2,2.3]
  @tag params_init: [10.0,3.0]
  test "fit - normal (20)", context do
    assert_fit {2.0,context[:params]}, {context[:data],context[:init],context[:cdf]}, context
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
    assert_fit {20.0,context[:params]}, {context[:data],context[:init],context[:cdf]}, context
  end

  @tag distrib: :exponential
  @tag cdf: :exponentialCDF
  @tag size: 2000
  @tag params: [1.7]
  @tag params_init: [1.0]
  test "fit - exponential (2000)", context do
    assert_fit {20.0,context[:params]}, {context[:data],context[:init],context[:cdf]}, context
  end

  @tag distrib: :exponential
  @tag cdf: :exponentialCDF
  @tag size: 200
  @tag params: [1.7]
  @tag params_init: [1.0]
  test "fit - exponential (200)", context do
    assert_fit {20.0,context[:params]}, {context[:data],context[:init],context[:cdf]}, context
  end

  @tag distrib: :exponential
  @tag cdf: :exponentialCDF
  @tag size: 20
  @tag params: [1.7]
  @tag params_init: [1.0]
  test "fit - exponential (20)", context do
    assert_fit {20.0,context[:params]}, {context[:data],context[:init],context[:cdf]}, context
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
    assert_fit {20.0,context[:params]}, {context[:data],context[:init],context[:cdf]}, context
  end

  @tag distrib: :wald
  @tag cdf: :waldCDF
  @tag size: 2000
  @tag params: [1.7,6.3]
  @tag params_init: [1.0,1.0]
  test "fit - wald (2000)", context do
    assert_fit {20.0,context[:params]}, {context[:data],context[:init],context[:cdf]}, context
  end

  @tag distrib: :wald
  @tag cdf: :waldCDF
  @tag size: 200
  @tag params: [1.7,6.3]
  @tag params_init: [1.0,1.0]
  test "fit - wald (200)", context do
    assert_fit {20.0,context[:params]}, {context[:data],context[:init],context[:cdf]}, context
  end

  @tag distrib: :wald
  @tag cdf: :waldCDF
  @tag size: 20
  @tag params: [1.7,6.3]
  @tag params_init: [1.0,6.0]
  test "fit - wald (20)", context do
    assert_fit {20.0,context[:params]}, {context[:data],context[:init],context[:cdf]}, context
  end

end
