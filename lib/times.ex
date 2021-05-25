defmodule Chi2fit.Times do

  # Copyright 2015-2021 Pieter Rijken
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

  @hours 24.0
  @default_workday {8.0,18.0}
  @default_epoch ~D[1970-01-01]

  @doc ~s"""
  Adjusts the times to working hours and/or work days.
  
  ## Options
  
      `workhours` - a 2-tuple containing the starting and ending hours of the work day (defaults
          to #{inspect @default_workday})
      `epoch` - the epoch to which all data elements are relative (defaults to #{@default_epoch})
      `saturday` - number of days since the epoch that corresponds to a Saturday (defaults
          to #{13 - Date.day_of_week(@default_epoch)})
      `correct` - whether to correct the times for working hours and weekdays; possible values
          `:worktime`, `:weekday`, `:"weekday+worktime"` (defaults to `false`)

  """
  @spec adjust_times(Enumerable.t, options :: Keyword.t) :: Enumerable.t
  def adjust_times(data, options) do
    {startofday,endofday} = options[:workhours] || @default_workday
    correct = options[:correct] || false
    epoch = options[:epoch] || @default_epoch
    sat = 13 - Date.day_of_week(epoch)
    saturday = options[:saturday] || sat

    data
    |> Stream.map(fn x ->
        case correct do
          :worktime -> map2workhours(x, startofday, endofday)
          :weekday -> map2weekdays(x, saturday)
          :"weekday+worktime" -> x |> map2workhours(startofday, endofday) |> map2weekdays(saturday)
          _ -> x
        end
      end)
    |> Enum.sort(&(&1>&2)) # Sort on new delivery times
  end

  @default_cutoff 0.01

  @doc """
  Returns a list of time differences (assumes an ordered list as input)
  
  ## Options
  
      `cutoff` - time differences below the cutoff are changed to the cutoff value (defaults to `#{@default_cutoff}`)
      `drop?` - whether to drop time differences below the cutoff (defaults to `false`)

  """
  @spec time_diff(data :: Enumrable.t, options :: Keyword.t) :: Enumerable.t
  def time_diff(data, options) do
    cutoff = options[:cutoff] || @default_cutoff
    drop? = options[:drop] || false

    data
    |> Stream.chunk_every(2,1,:discard)
    |> Stream.map(fn [x,y]->x-y end)
    |> Stream.transform(nil,fn x,_acc ->
          {
            cond do
              x < cutoff and drop? -> []
              x < cutoff -> [cutoff]
              true -> [x]
            end,
            nil
          }
        end)
    |> (& if is_function(data, 2), do: &1, else: Enum.into(&1, [])).()
  end

  @doc """
  Maps the time of a day into the working hour period
  
  Scales the resulting part of the day between 0..1.
  
  ## Arguments
  
      `t` - date and time of day as a float; the integer part specifies the day and the fractional part the hour of the day
      `startofday` - start of the work day in hours
      `endofday` - end of the working day in hours
  
  ## Example
  
      iex> map2workhours(43568.1, 8, 18)
      43568.0
  
      iex> map2workhours(43568.5, 8, 18)
      43568.4
  
  """
  @spec map2workhours(t :: number, startofday :: number, endofday :: number) :: number
  def map2workhours(t,startofday,endofday)
    when startofday>0 and startofday<endofday and endofday<@hours do
    frac = t - trunc(t)
    hours = endofday - startofday
    trunc(t) + min(max(0.0,frac - startofday/@hours),hours/@hours) * @hours/hours
  end

  @doc """
  Maps the date to weekdays such that weekends are eliminated; it does so with respect to a given Saturday
  
  ## Example
  
      iex> map2weekdays(43568.123,43566)
      43566.123
  
      iex> map2weekdays(43574.123,43566)
      43571.123
  
  """
  @spec map2weekdays(t :: number, sat :: pos_integer) :: number
  def map2weekdays(t, sat) when is_integer(sat) do
    offset = rem trunc(t)-sat, 7
    weeks = div trunc(t)-sat, 7
    
    part_of_day = t - trunc(t)
    sat + 5*weeks + max(0.0,offset-2.0) + part_of_day
  end

  @doc """
  Returns a `Stream` that generates a stream of dates.
  
  ## Examples
  
      iex> intervals(end: ~D[2019-06-01]) |> Enum.take(4)
      [~D[2019-06-01], ~D[2019-05-16], ~D[2019-05-01], ~D[2019-04-16]]
  
      iex> intervals(end: ~D[2019-06-01], type: :weekly) |> Enum.take(4)
      [~D[2019-06-01], ~D[2019-05-18], ~D[2019-05-04], ~D[2019-04-20]]
  
      iex> intervals(end: ~D[2019-06-01], type: :weekly, weeks: 1) |> Enum.take(4)
      [~D[2019-06-01], ~D[2019-05-25], ~D[2019-05-18], ~D[2019-05-11]]

      iex> intervals(end: ~D[2019-06-01], type: :weekly, weeks: [3,2]) |> Enum.take(4)
      [~D[2019-06-01], ~D[2019-05-11], ~D[2019-04-27], ~D[2019-04-13]]

  """
  @spec intervals(options :: Keyword.t) :: Enumerable.t
  def intervals(options \\ []) do
    type = options[:type] || :half_month
    periods = case options[:weeks] do
      nil -> [2]
      x when is_number(x) -> [x]
      list when is_list(list) -> list
    end
    last = options[:end] || Date.utc_today()

    case type do
      :half_month ->
        recent = case last do
          date = %Date{day: day} when day > 16 -> %Date{date | day: 1, month: date.month+1}
          date = %Date{day: day} when day > 1 -> %Date{date | day: 16}
          date = %Date{day: 1} -> date
        end

        recent |> Stream.iterate(fn
          previous = %Date{day: 16} -> Timex.shift previous, days: -15
          previous = %Date{day: 1} -> Timex.shift previous, days: +15, months: -1
        end)
      
      :weekly ->
        Stream.resource(
          fn -> {last,periods} end,
          fn
            {current,[p]} ->
              next = Timex.shift(current, weeks: -p)
              {[current], {next,[p]}}
            {current,[p|rest]} ->
              next = Timex.shift(current, weeks: -p)
              {[current], {next,rest}}
          end,
          fn _ -> [] end)
    end
  end

  @doc """
  Counts the number of dates (`datelist`) that is between consecutive dates in `intervals` and returns the result as a list of numbers.
  """
  @spec throughput(intervals :: Enumerable.t, datelist :: [NaiveDateTime.t]) :: [number]
  def throughput(intervals, datelist) do
    intervals
    |> Stream.chunk_every(2, 1, :discard)
    |> Stream.transform(datelist, fn
        _, [] ->
          {:halt, []}
        [d1,d2], acc ->
          {left,right} = Enum.split_with(acc, fn d -> Timex.between?(d,d2,d1) end)
          {[{d1,Enum.count(left)}],right}
      end)
    |> Enum.map(fn {_d,count} -> count end)
  end

end
