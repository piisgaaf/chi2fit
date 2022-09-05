defmodule Chi2fit.Utilities do

  # Copyright 2015-2017 Pieter Rijken
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
  Provides various utilities:
  
    * Bootstrapping
    * Creating Cumulative Distribution Functions / Histograms from sample data
    * Autocorrelation coefficients
  
  """

  alias Chi2fit.Distribution, as: D
  alias Chi2fit.Distribution.Utilities, as: U
  alias Chi2fit.Fit, as: F
  alias Chi2fit.Math, as: M
  alias Chi2fit.Matrix, as: Mx
  alias Chi2fit.Statistics, as: S
  
  @typedoc "Average and standard deviationm (error)"
  @type avgsd :: {avg :: float, sd :: float}

  @doc """
  Reads data from a file specified by `filename` and returns a stream with the data parsed as floats.

  It expects a single data point on a separate line and removes entries that:
  
    * are not floats, and
    * smaller than zero (0)

  """
  @spec read_data(filename::String.t) :: Enumerable.t
  def read_data(filename) do
    filename
    |> File.stream!([],:line)
    |> Stream.flat_map(&String.split(&1,"\r",trim: true))
    |> Stream.filter(&is_tuple(Float.parse(&1)))
    |> Stream.map(&elem(Float.parse(&1),0))
    |> Stream.filter(&(&1 >= 0.0))
  end


  @doc """
  Unzips lists of 1-, 2-, 3-, 4-, 5-, 6-, 7-, and 8-tuples.
  """
  @spec unzip(list::[tuple]) :: tuple
  def unzip([]), do: {}
  def unzip([{}|_]), do: {}
  def unzip(list=[{_}|_]), do: {Enum.map(list,fn {x}->x end)}
  def unzip(list=[{_,_}|_]), do: Enum.unzip(list)
  def unzip(list=[{_,_,_}|_]) do
    {
      list |> Enum.map(&elem(&1,0)),
      list |> Enum.map(&elem(&1,1)),
      list |> Enum.map(&elem(&1,2))
    }
  end
  def unzip(list=[{_,_,_,_}|_]) do
    {
      list |> Enum.map(&elem(&1,0)),
      list |> Enum.map(&elem(&1,1)),
      list |> Enum.map(&elem(&1,2)),
      list |> Enum.map(&elem(&1,3))
    }
  end
  def unzip(list=[{_,_,_,_,_}|_]) do
    {
      list |> Enum.map(&elem(&1,0)),
      list |> Enum.map(&elem(&1,1)),
      list |> Enum.map(&elem(&1,2)),
      list |> Enum.map(&elem(&1,3)),
      list |> Enum.map(&elem(&1,4))
    }
  end
  def unzip(list=[{_,_,_,_,_,_}|_]) do
    {
      list |> Enum.map(&elem(&1,0)),
      list |> Enum.map(&elem(&1,1)),
      list |> Enum.map(&elem(&1,2)),
      list |> Enum.map(&elem(&1,3)),
      list |> Enum.map(&elem(&1,4)),
      list |> Enum.map(&elem(&1,5))
    }
  end
  def unzip(list=[{_,_,_,_,_,_,_}|_]) do
    {
      list |> Enum.map(&elem(&1,0)),
      list |> Enum.map(&elem(&1,1)),
      list |> Enum.map(&elem(&1,2)),
      list |> Enum.map(&elem(&1,3)),
      list |> Enum.map(&elem(&1,4)),
      list |> Enum.map(&elem(&1,5)),
      list |> Enum.map(&elem(&1,6))
    }
  end
  def unzip(list=[_|_]) do
    0..tuple_size(hd(list))-1
    |> Enum.reduce({},fn i,tup -> Tuple.append(tup,list |> Enum.map(&elem(&1,i))) end)
  end

  ##
  ## Local functions
  ##

  @doc """
  Outputs and formats the errors that result from a call to `Chi2fit.Fit.chi2/4`
  
  Errors are tuples of length 2 and larger: `{[min1,max1], [min2,max2], ...}`.
  """
  @spec puts_errors(device :: IO.device(), errors :: tuple()) :: none()
  def puts_errors(device \\ :stdio, errors) do
    errors
    |> Tuple.to_list
    |> Enum.with_index
    |> Enum.each(fn
        {[mn,mx],0} -> IO.puts device, "\t\t\tchi2:\t\t#{mn}\t-\t#{mx}"
        {[mn,mx],_} -> IO.puts device, "\t\t\tparameter:\t#{mn}\t-\t#{mx}"
    end)
  end

  @doc """
  Displays results of the function `Chi2fit.Fit.chi2probe/4`
  """
  @spec display(device :: IO.device(), F.chi2probe() | avgsd()) :: none()
  def display(device \\ :stdio, results)
  def display(device,{chi2, parameters, errors, _saved}) do
      IO.puts device,"Initial guess:"
      IO.puts device,"    chi2:\t\t#{chi2}"
      IO.puts device,"    pars:\t\t#{inspect parameters}"
      IO.puts device,"    ranges:\t\t#{inspect errors}\n"
  end
  def display(device,{avg, sd, direction}) when direction in [:-, :+] do
    op = case direction do
      :+ -> &Kernel.+/2
      :- -> &Kernel.-/2
    end
    IO.puts device,"50%    => #{:math.ceil(avg)} units"
    IO.puts device,"84%    => #{:math.ceil(op.(avg,sd))} units"
    IO.puts device,"97.5%  => #{:math.ceil(op.(avg,2*sd))} units"
    IO.puts device,"99.85% => #{:math.ceil(op.(avg,3*sd))} units"
  end
  def display(device, {list, direction}) when is_list(list) and direction in [:-, :+] do
    sorted = case direction do
      :+ -> list |> Enum.sort
      :- -> list |> Enum.sort |> Enum.reverse
    end
    max = length(sorted)
    IO.puts device,"50%    => #{Enum.at(sorted,round(0.5 * max))} units"
    IO.puts device,"84%    => #{Enum.at(sorted,round(0.84 * max))} units"
    IO.puts device,"97.5%  => #{Enum.at(sorted,round(0.975 * max))} units"
    IO.puts device,"99.85% => #{Enum.at(sorted,round(0.9985 * max))} units"
  end

  @doc """
  Pretty prints subsequences.
  """
  @spec display_subsequences(device :: IO.device(), trends :: list(), intervals :: [NaiveDateTime.t]) :: none()
  def display_subsequences(device \\ :stdio, trends, intervals) do
    trends
    |> Stream.transform(1, fn arg={_,_,data}, index -> { [{arg, Enum.at(intervals,index)}], index+length(data)} end)
    |> Stream.each(fn
        {{chimin, [rate], subdata},date} ->
          IO.puts device, "Subsequence ending @#{Timex.format!(date,~S({Mshort}, {D} {YYYY}))}"
          IO.puts device, "----------------------------------"
          IO.puts device, "    chi2@minimum:  #{Float.round(chimin,1)}"
          IO.puts device, "    delivery rate: #{Float.round(rate,1)}"
          IO.puts device, "    subsequence:   #{inspect(subdata, charlists: :as_lists)}"
          IO.puts device, ""
      end)
    |> Stream.run()
  end

  @doc """
  Displays results of the function `Chi2fit.Fit.chi2fit/4`
  """
  @spec display(device :: IO.device(),hdata :: S.ecdf(),model :: U.model(),F.chi2fit(),options :: Keyword.t) :: none()
  def display(device \\ :stdio,hdata,model,{chi2, cov, parameters, errors},options) do
      param_errors = cov |> Mx.diagonal |> Enum.map(fn x->x|>abs|>:math.sqrt end)

      IO.puts device,"Final:"
      IO.puts device,"    chi2:\t\t#{chi2}"
      IO.puts device,"    Degrees of freedom:\t#{length(hdata)-D.size(model)}"
      IO.puts device,"    gradient:\t\t#{inspect M.jacobian(parameters,&F.chi2(hdata,fn x->D.cdf(model).(x,&1) end,fn _->0.0 end,options),options)}"
      IO.puts device,"    parameters:\t\t#{inspect parameters}"
      IO.puts device,"    errors:\t\t#{inspect param_errors}"
      IO.puts device,"    ranges:"
      puts_errors device,errors
  end

  @doc """
  Walks a map structure while applying the function `fun`.
  """
  @spec analyze(map :: %{}, fun :: (([number],Keyword.t) -> Keyword.t), options :: Keyword.t) :: Keyword.t
  def analyze(map = %{}, fun, options) do
      map |> Enum.reduce(%{}, fn {k,v},acc -> Map.put(acc,k,analyze(v,fun,options)) end)
  end
  def analyze(data, fun, options) when is_list(data) do
      cond do
          Keyword.keyword?(data) ->
              Keyword.merge(data, fun.(data,Keyword.put(options,:bin,data[:bin])))
          true ->
              analyze([throughput: data, bin: options[:bin]], fun, options)
      end
  end

  @doc """
  Pretty-prints a nested array-like structure (list or tuple) as a table.
  """
  @spec as_table(rows :: [any], header :: tuple()) :: list()
  def as_table(rows, header) do
    map = 1..tuple_size(header) |> Enum.map(&{&1,0}) |> Enum.into(%{})
    table = [header|rows] |> _to_string()
    map = Enum.reduce(table, map, fn row,acc ->
      row
      |> Enum.with_index(1)
      |> Enum.reduce(acc, fn {str,i},acc2 -> Map.update!(acc2, i, fn v -> max(v,String.length(str)) end) end)
    end)
    table
    |> Enum.with_index()
    |> Enum.each(fn
      {row, 0} ->
        IO.puts row |> Enum.with_index(1) |> Enum.map(fn {str,i} -> String.pad_trailing(str, map[i]) end) |> Enum.join("|")
        IO.puts row |> Enum.with_index(1) |> Enum.map(fn {_,i} -> String.duplicate("-",map[i]) end) |> Enum.join("|")
      {row, _} ->
        IO.puts row |> Enum.with_index(1) |> Enum.map(fn {str,i} -> String.pad_trailing(str, map[i]) end) |> Enum.join("|")
    end)
    
    rows
  end
  
  defp _to_string(list) when is_list(list), do: list |> Enum.map(&_to_string/1)
  defp _to_string(tuple) when is_tuple(tuple), do: tuple |> Tuple.to_list |> _to_string()
  defp _to_string(string) when is_binary(string), do: string
  defp _to_string(float) when is_float(float), do: "#{float}"
  defp _to_string(integer) when is_integer(integer), do: "#{integer}"



  @doc ~S"""

  ## Examples

      iex> subsequences []
      []
  
      iex> subsequences [:a, :b]
      [[:a], [:a, :b]]
  
      iex> Stream.cycle([1,2,3]) |> subsequences |> Enum.take(4)
      [[1], [1, 2], [1, 2, 3], [1, 2, 3, 1]]

  """
  @spec subsequences(Enumerable.t) :: Enumerable.t
  def subsequences(stream) when is_function(stream, 2) do
    stream
    |> Stream.transform([], fn x,acc -> {[Enum.reverse([x|acc])], [x|acc]} end)
  end
  def subsequences(list) do
    {result, _} = list
    |> Enum.reduce({[],[]}, fn  x, {res,acc} -> {[Enum.reverse([x|acc])|res], [x|acc]} end)
    Enum.reverse(result)
  end

end
