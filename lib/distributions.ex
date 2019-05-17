defmodule Chi2fit.Distribution do

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
  Provides various distributions.
  """

  import Chi2fit.Utilities
  
  @typedoc "A probability distribution function"
  @type distribution() :: ((...) :: term())
  
  @typedoc "Cumulative Distribution function"
  @type cdf() :: ((number) :: number())

  @typedoc "Keyword list containing the CDF function and the number of parameters"
  @type model :: [fun: cdf(), df: pos_integer()]
  
  defmodule UnsupportedDistributionError do
    defexception message: "Unsupported distribution function"
  end

  @doc """
  Returns the model for a name.
  
  The kurtosis is the so-called 'excess kurtosis'.
  
  Supported disributions:
      "wald" - The Wald or Inverse Gauss distribution,
      "weibull" - The Weibull distribution,
      "exponential" - The exponential distribution,
      "poisson" - The Poisson distribution,
      "normal" - The normal or Gaussian distribution,
      "fechet" - The Fréchet distribution,
      "nakagami" - The Nakagami distribution,
      "sep" - The Skewed Exponential Power distribution (Azzalini),
      "erlang" - The Erlang distribution,
      "sep0" - The Skewed Exponential Power distribution (Azzalini) with location parameter set to zero (0).
      
  ## Options
  Available only for the SEP distribution, see 'sepCDF/5'.
  """
  @spec model(name::String.t, options::Keyword.t) :: model()
  def model(name, options \\ []) do
    case name do
      "wald" ->
        d = %Distribution.Wald{}
        [
          fun: Distribution.cdf(d),
          df: Distribution.size(d),
          skewness: Distribution.skewness(d),
          kurtosis: Distribution.kurtosis(d)
        ]
      "weibull" ->
        d = %Distribution.Weibull{}
        [
          fun: Distribution.cdf(d),
          df: Distribution.size(d),
          skewness: Distribution.skewness(d),
          kurtosis: Distribution.kurtosis(d)
        ]
      "exponential" ->
        d = %Distribution.Exponential{}
        [
          fun: Distribution.cdf(d),
          df: Distribution.size(d),
          skewness: Distribution.skewness(d),
          kurtosis: Distribution.kurtosis(d)
        ]
      "frechet" ->
        d = %Distribution.Frechet{}
        [
          fun: Distribution.cdf(d),
          df: Distribution.size(d),
          skewness: Distribution.skewness(d),
          kurtosis: Distribution.kurtosis(d)
        ]
      "nakagami" ->
        d = %Distribution.Nakagami{}
        [
          fun: Distribution.cdf(d),
          df: Distribution.size(d),
          skewness: Distribution.skewness(d),
          kurtosis: Distribution.kurtosis(d)
        ]
      "poisson" ->
        d = %Distribution.Poisson{}
        [
          fun: Distribution.cdf(d),
          df: Distribution.size(d),
          skewness: Distribution.skewness(d),
          kurtosis: Distribution.kurtosis(d)
        ]
      {"poisson", period} when is_number(period) and period>0 ->
        d = %Distribution.Poisson{period: period}
        [
          fun: Distribution.cdf(d),
          df: Distribution.size(d),
          skewness: Distribution.skewness(d),
          kurtosis: Distribution.kurtosis(d)
        ]
      "erlang" ->
        d = %Distribution.Erlang{}
        [
          fun: Distribution.cdf(d),
          df: Distribution.size(d),
          skewness: Distribution.skewness(d),
          kurtosis: Distribution.kurtosis(d)
        ]
      {"erlang", batches} ->
        d = %Distribution.Erlang{batches: batches}
        [
          fun: Distribution.cdf(d),
          df: Distribution.size(d),
          skewness: Distribution.skewness(d),
          kurtosis: Distribution.kurtosis(d)
        ]
      "normal" ->
        d = %Distribution.Normal{}
        [
          fun: Distribution.cdf(d),
          df: Distribution.size(d),
          skewness: Distribution.skewness(d),
          kurtosis: Distribution.kurtosis(d)
        ]
      "sep" ->
        d = %Distribution.SEP{options: options}
        [
          fun: Distribution.cdf(d),
          df: Distribution.size(d),
          skewness: Distribution.skewness(d),
          kurtosis: Distribution.kurtosis(d)
        ]
      "sep0" ->
        d = %Distribution.SEP{offset: 0.0, options: options}
        [
          fun: Distribution.cdf(d),
          df: Distribution.size(d),
          skewness: Distribution.skewness(d),
          kurtosis: Distribution.kurtosis(d)
        ]
      unknown ->
        raise UnsupportedDistributionError, message: "Unsupported cumulative distribution function '#{inspect unknown}'"
    end
  end

  @doc """
  Guesses what distribution is likely to fit the sample data
  """
  @spec guess(sample::[number], n::integer, list::[String.t] | String.t) :: [any]
  def guess(sample,n \\ 100,list \\ ["exponential","poisson","normal","erlang","wald","sep","weibull","frechet","nakagami"])
  def guess(sample,n,list) when is_integer(n) and n>0 and is_list(list) do
    {{skewness,err_s},{kurtosis,err_k}} = sample |> cullen_frey(n) |> cullen_frey_point
    list
    |> Enum.flat_map(
      fn
        distrib ->
          r = sample
          |> guess(n,distrib)
          |> Enum.map(fn {s,k}->((skewness-s)/err_s)*((skewness-s)/err_s) + ((kurtosis-k)/err_k)*((kurtosis-k)/err_k) end)
          |> Enum.min
          [{distrib,r}]
      end)
    |> Enum.sort(fn {_,r1},{_,r2} -> r1<r2 end)
  end
  def guess(_sample,n,distrib) when is_integer(n) and n>0 do
    model = model(distrib)
    params = 1..model[:df]
    1..n
    |> Enum.map(fn _ -> Enum.map(params, fn _ -> 50*:rand.uniform end) end)
    |> Enum.flat_map(fn
      pars ->
        try do
          s = model[:skewness].(pars)
          k = model[:kurtosis].(pars)
          [{s,k}]
        rescue
          _error -> []
        end
    end)
  end

end
