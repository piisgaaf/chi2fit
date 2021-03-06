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

  defp to_number(string) when is_binary(string) do
    case Integer.parse(string) do
      {val,""} -> val
      _ -> String.to_float(string)
    end
  end
  defp to_numbers(list), do: String.split(list," ") |> Enum.map(& to_number &1)

  @doc ~S"""
  
  ## Examples

      iex> ~M(3 4 5)
      %Distribution.Uniform{pars: [3, 4, 5]}

      iex> ~M(3 4 5)u
      %Distribution.Uniform{pars: [3, 4, 5]}

      iex> ~M()d
      %Distribution.Dice{mode: :regular}

      iex> ~M()dgk
      %Distribution.Dice{mode: :gk4}

      iex> ~M(1.2)p
      %Distribution.Poisson{pars: [1.2], period: 1.0}

      iex> ~M(1.2 5.4)w
      %Distribution.Weibull{pars: [1.2, 5.4]}

      iex> ~M(1.2 5.4)wald
      %Distribution.Wald{pars: [1.2, 5.4]}

  """
  def sigil_M(str, ''), do: %Distribution.Uniform{pars: to_numbers(str)}
  def sigil_M(str, [?u|_]), do: %Distribution.Uniform{pars: to_numbers(str)}
  def sigil_M("", 'coin'), do: %Distribution.Coin{}
  def sigil_M(str, [?c|_]), do: %Distribution.Constant{pars: [to_number(str)]}
  def sigil_M("", 'dgk'), do: %Distribution.Dice{mode: :gk4}
  def sigil_M("", [?d|_]), do: %Distribution.Dice{mode: :regular}

  def sigil_M(str, [?b|_]), do: %Distribution.Bernoulli{pars: [to_number(str)]}

  def sigil_M(str, 'erlangb') do
    [rate, batches] = to_numbers(str)
    %Distribution.Erlang{pars: [rate], batches: batches}
  end
  def sigil_M(str, 'erlang'), do: %Distribution.Erlang{pars: to_numbers(str)}
  def sigil_M(str, [?e|_]), do: %Distribution.Exponential{pars: [to_number(str)]}
  def sigil_M(str, 'pp') do
    [rate, period] = to_numbers(str)
    %Distribution.Poisson{pars: [rate], period: period}
  end
  def sigil_M(str, [?p|_]), do: %Distribution.Poisson{pars: [to_number(str)]}

  def sigil_M(str, 'nakagami'), do: %Distribution.Nakagami{pars: to_numbers(str)}
  def sigil_M(str, [?n|_]), do: %Distribution.Normal{pars: to_numbers(str)}

  def sigil_M(str, 'wald'), do: %Distribution.Wald{pars: to_numbers(str)}
  def sigil_M("", [?w|_]), do: %Distribution.Weibull{}
  def sigil_M(str, [?w|_]), do: %Distribution.Weibull{pars: to_numbers(str)}

  def sigil_M(str, [?f|_]), do: %Distribution.Frechet{pars: to_numbers(str)}
  def sigil_M(str, 'sz'), do: %Distribution.SEP{pars: to_numbers(str), offset: 0.0}
  def sigil_M(str, [?s|_]), do: %Distribution.SEP{pars: to_numbers(str)}

  def sigil_M(_term, modifiers) do
    raise UnsupportedDistributionError, message: "Unsupported modifiers #{inspect modifiers}"
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
