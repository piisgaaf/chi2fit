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
  @spec model(name::String.t, options::Keyword.t) :: any
  def model(name, options \\ []) do
    params = options[:pars] || nil
    
    case name do
      "constant" -> %Distribution.Constant{pars: params}
      "uniform" -> %Distribution.Uniform{pars: params}
      "dice" -> %Distribution.Dice{mode: :regular}
      "dice_gk4" -> %Distribution.Dice{mode: :gk4}
      "wald" -> %Distribution.Wald{pars: params}
      "weibull" -> %Distribution.Weibull{pars: params}
      "exponential" -> %Distribution.Exponential{pars: params}
      "frechet" -> %Distribution.Frechet{pars: params}
      "nakagami" -> %Distribution.Nakagami{pars: params}
      "poisson" -> %Distribution.Poisson{pars: params}
      {"poisson", period} when is_number(period) and period>0 -> %Distribution.Poisson{pars: params,period: period}
      "erlang" -> %Distribution.Erlang{pars: params}
      {"erlang", batches} when is_number(batches) and batches>0 -> %Distribution.Erlang{pars: params,batches: batches}
      "normal" -> %Distribution.Normal{pars: params}
      "sep" -> %Distribution.SEP{pars: params,options: options}
      "sep0" -> %Distribution.SEP{pars: params,offset: 0.0, options: options}
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
    params = 1..Distribution.size(model)
    1..n
    |> Enum.map(fn _ -> Enum.map(params, fn _ -> 50*:rand.uniform end) end)
    |> Enum.flat_map(fn
      pars ->
        try do
          s = Distribution.skewness(model).(pars)
          k = Distribution.kurtosis(model).(pars)
          [{s,k}]
        rescue
          _error -> []
        end
    end)
  end

end
