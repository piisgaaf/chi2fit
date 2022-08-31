defmodule Chi2fit.Distribution.Utilities do

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

  import Chi2fit.Statistics
  alias Chi2fit.Distribution, as: D


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
      "frechet" - The Fréchet distribution,
      "nakagami" - The Nakagami distribution,
      "sep" - The Skewed Exponential Power distribution (Azzalini),
      "erlang" - The Erlang distribution,
      "sep0" - The Skewed Exponential Power distribution (Azzalini) with location parameter set to zero (0),
      "tw" - The Tracy-Widom distributions TW1, TW2, and TW4.
      
  ## Options
  Available only for the SEP distribution, see 'sepCDF/5'.
  """
  @spec model(name::String.t, options::Keyword.t) :: any
  def model(name, options \\ []) do
    params = options[:pars] || nil

    case name do
      "constant" -> %D.Constant{pars: params}
      "uniform" -> %D.Uniform{pars: params}
      "dice" -> %D.Dice{mode: :regular}
      "dice_gk4" -> %D.Dice{mode: :gk4}
      "wald" -> %D.Wald{pars: params}
      "weibull" -> %D.Weibull{pars: params}
      "exponential" -> %D.Exponential{pars: params}
      "frechet" -> %D.Frechet{pars: params}
      "nakagami" -> %D.Nakagami{pars: params}
      "poisson" -> %D.Poisson{pars: params}
      {"poisson", period} when is_number(period) and period>0 -> %D.Poisson{pars: params, period: period}
      "erlang" -> %D.Erlang{pars: params}
      {"erlang", batches} when is_number(batches) and batches>0 -> %D.Erlang{pars: params, batches: batches}
      "normal" -> %D.Normal{pars: params}
      "sep" -> %D.SEP{pars: params, options: options}
      "sep0" -> %D.SEP{pars: params, offset: 0.0, options: options}
      "tw" -> %D.TracyWidom{type: 1}
      "tw1" -> %D.TracyWidom{type: 1}
      "tw2" -> %D.TracyWidom{type: 2}
      "tw4" -> %D.TracyWidom{type: 4}
      unknown ->
        raise UnsupportedDistributionError, message: "Unsupported cumulative distribution function '#{inspect unknown}'"
    end
  end

  defp to_number(string) when is_binary(string) do
    case Integer.parse(string) do
      {val, ""} -> val
      _ -> String.to_float(string)
    end
  end
  defp to_numbers(list), do: String.split(list, " ") |> Enum.map(& to_number &1)

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
  def sigil_M(str, ''), do: %D.Uniform{pars: to_numbers(str)}
  def sigil_M(str, [?u|_]), do: %D.Uniform{pars: to_numbers(str)}
  def sigil_M("", 'coin'), do: %D.Coin{}
  def sigil_M(str, [?c|_]), do: %D.Constant{pars: [to_number(str)]}
  def sigil_M("", 'dgk'), do: %D.Dice{mode: :gk4}
  def sigil_M("", [?d|_]), do: %D.Dice{mode: :regular}

  def sigil_M(str, [?b|_]), do: %D.Bernoulli{pars: [to_number(str)]}

  def sigil_M(str, 'erlangb') do
    [rate, batches] = to_numbers(str)
    %D.Erlang{pars: [rate], batches: batches}
  end
  def sigil_M(str, 'erlang'), do: %D.Erlang{pars: to_numbers(str)}
  def sigil_M(str, [?e|_]), do: %D.Exponential{pars: [to_number(str)]}
  def sigil_M(str, 'pp') do
    [rate, period] = to_numbers(str)
    %D.Poisson{pars: [rate], period: period}
  end
  def sigil_M(str, [?p|_]), do: %D.Poisson{pars: [to_number(str)]}

  def sigil_M(str, 'nakagami'), do: %D.Nakagami{pars: to_numbers(str)}
  def sigil_M(str, [?n|_]), do: %D.Normal{pars: to_numbers(str)}

  def sigil_M(str, 'wald'), do: %D.Wald{pars: to_numbers(str)}
  def sigil_M("", [?w|_]), do: %D.Weibull{}
  def sigil_M(str, [?w|_]), do: %D.Weibull{pars: to_numbers(str)}

  def sigil_M(str, [?f|_]), do: %D.Frechet{pars: to_numbers(str)}
  def sigil_M(str, 'sz'), do: %D.SEP{pars: to_numbers(str), offset: 0.0}
  def sigil_M(str, [?s|_]), do: %D.SEP{pars: to_numbers(str)}

  def sigil_M(_str, 'tw'), do: %D.TracyWidom{pars: [], type: 1}
  def sigil_M(_str, 'tww'), do: %D.TracyWidom{pars: [], type: 2}
  def sigil_M(_str, 'twwww'), do: %D.TracyWidom{pars: [], type: 4}
  
  def sigil_M(_term, modifiers) do
    raise UnsupportedDistributionError, message: "Unsupported modifiers #{inspect modifiers}"
  end

  @doc """
  Guesses what distribution is likely to fit the sample data
  """
  @spec guess(sample::[number], n::integer, list::[String.t] | String.t) :: [any]
  def guess(sample,n \\ 100,list \\ ["exponential","poisson","normal","erlang","wald","sep","weibull","frechet","nakagami","tw","tww"])
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
    params = 1..D.size(model)
    1..n
    |> Enum.map(fn _ -> Enum.map(params, fn _ -> 50*:rand.uniform end) end)
    |> Enum.flat_map(fn
      pars ->
        try do
          s = D.skewness(model).(pars)
          k = D.kurtosis(model).(pars)
          [{s,k}]
        rescue
          _error -> []
        end
    end)
  end

end
