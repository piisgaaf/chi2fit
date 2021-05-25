defmodule Chi2fit.MonteCarlo do

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

  import Chi2fit.Statistics

  alias Chi2fit.Matrix, as: M

  @doc """
  Basic Monte Carlo simulation to repeatedly run a simulation multiple times.
  
  ## Options
  
      `:collect_all?` - If true, collects data from each individual simulation run and returns this an the third element of the result tuple
  
  """
  @spec mc(iterations :: pos_integer, fun :: ((pos_integer) -> float), options :: Keyword.t) :: {avg :: float, sd :: float, tries :: [float]} | {avg :: float, sd :: float}
  def mc(iterations,fun,options \\ []) do
    all? = options[:collect_all?] || false

    tries = 1..iterations |> Enum.map(fn _ -> fun.() end)
    avg = moment tries, 1
    sd = :math.sqrt momentc(tries,2,avg)
    if all?, do: {avg,sd,tries}, else: {avg,sd}
  end

  def total_mc(result, fun, mode \\ :use_bounds, iterations \\ 1000) do
    {_, cov, parameters, _} = result

    ranges = case mode do
      :use_bounds ->
        # Pick up the error in the paramater value
        errors = cov |> M.diagonal |> Enum.map(fn x -> x|>abs|>:math.sqrt end)
        Enum.zip(parameters, errors) |> Enum.map(fn {par, err} -> [par - err, par + err] end)

      :use_ranges ->
        {_, _, _, parranges} = result
        parranges |> Tuple.to_list |> tl
    end
    outcomes = ranges
    |> List.foldr([], fn
      pair, [] -> pair
      [left, right], acc -> Enum.flat_map(acc, & [[left|List.wrap(&1)], [right|List.wrap(&1)]])
    end)
    |> Enum.map(& mc(iterations, fun.(&1)))
    |> Enum.map(& elem(&1,0))

    {Enum.min(outcomes), Enum.max(outcomes)}
  end

  @doc """
  Forecasts how many time periods are needed to complete `size` items
  
  Related functions: `forecast_duration/2` and `forecast_items/2`.
  """
  @spec forecast(fun :: (() -> non_neg_integer),size :: pos_integer, tries :: pos_integer, update :: (() -> number)) :: number
  def forecast(fun, size, tries \\ 0,update \\ fn -> 1 end)
  def forecast(fun, size, tries, update) when size>0 do
      forecast(fun, size-fun.(),tries+update.(),update)
  end
  def forecast(_fun,_size,tries,_update), do: tries

  @doc """
  Returns a function for forecasting the duration to complete a number of items.
  
  This function is a wrapper for `forecast/4`.

  ## Arguments
  
      `data` - either a data set to base the forecasting on, or a function that returns (random) numbers
      `size` - the number of items to complete

  """
  @spec forecast_duration(data :: [number] | (()->number), size :: pos_integer) :: (() -> number)
  def forecast_duration(data, size) when is_list(data) do
    fn -> forecast(fn -> Enum.random(data) end, size) end
  end
  def forecast_duration(fun, size) when is_function(fun,0) do
    fn -> forecast(fun, size) end
  end

  @doc """
  Returns a function for forecasting the number of completed items in a number periods.
  
  This function is a wrapper for `forecast/4`.

  ## Arguments
  
      `data` - either a data set to base the forecasting on, or a function that returns (random) numbers
      `periods` - the number of periods to forecast the number of completed items for

  """
  @spec forecast_items(data :: [number] | (()->number), periods :: pos_integer) :: (() -> number)
  def forecast_items(data, periods) when is_list(data) do
    fn -> forecast(fn -> 1 end, periods, 0, fn -> Enum.random(data) end) end
  end
  def forecast_items(fun, periods) when is_function(fun,0) do
    fn -> forecast(fn -> 1 end, periods, 0, fun) end
  end

end
