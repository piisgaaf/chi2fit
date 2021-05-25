defmodule Chi2fit.CSV do

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

  alias Timex.Format.DateTime.Formatters

  @doc ~S"""
  Reads CSV data, extracts one column, and returns it as a list of `NaiveDateTime`.
  
  ## Examples
  
      iex> csv = ["Done","2019/05/01","2019/06/01"] |> Stream.map(& &1)
      ...> csv_to_list csv, "Done", header?: true
      [~N[2019-06-01 00:00:00], ~N[2019-05-01 00:00:00]]
  
      iex> csv = ["Done","2019/May/01","2019/Jun/01"] |> Stream.map(& &1)
      ...> csv_to_list csv, "Done", header?: true, format: "{YYYY}/{Mshort}/{0D}"
      [~N[2019-06-01 00:00:00], ~N[2019-05-01 00:00:00]]
  
      iex> csv = ["Done","2019/May/01","2019/06/01"] |> Stream.map(& &1)
      ...> csv_to_list csv, "Done", header?: true, format: "{YYYY}/{Mshort}/{0D}"
      [~N[2019-05-01 00:00:00]]
  
      iex> csv = ["Done","2019/May/01","2019/06/01"] |> Stream.map(& &1)
      ...> csv_to_list csv, "Done", header?: true, format: ["{YYYY}/{Mshort}/{0D}","{YYYY}/{0M}/{0D}"]
      [~N[2019-06-01 00:00:00], ~N[2019-05-01 00:00:00]]

      iex> csv = ["Done","2019/May/01","2019/Jun/01"] |> Stream.map(& &1)
      ...> csv_to_list csv, "Done", header?: true, format: ["%Y/%b/%d"], parser: :strftime
      [~N[2019-06-01 00:00:00], ~N[2019-05-01 00:00:00]]

  """
  @spec csv_to_list(csvcata :: Enumerable.t, key :: String.t, options :: Keyword.t) :: [NaiveDateTime.t]
  def csv_to_list(csvdata, key, options \\ []) do
    header? = options[:header?] || false
    format = options[:format] || "{YYYY}/{0M}/{0D}"
    separator = options[:separator] || ?,
    parser = case options[:parser] do
      :strftime -> Formatters.Strftime
      :default -> Formatters.Default
      nil -> Formatters.Default
      _ -> Formatters.Default
    end

    formats = if is_list(format), do: format, else: [format]

    csvdata
    |> CSV.decode!(headers: header?, separator: separator)
    |> Stream.filter(& Map.fetch!(&1, key) != "")
    |> Stream.map(& Map.fetch!(&1, key))
    |> Stream.map(fn t -> Enum.reduce_while(formats,nil,fn f,_acc ->
        case Timex.parse(t, f, parser) do
          {:ok,n} -> {:halt,n}
          {:error,_msg} -> {:cont,nil}
        end
      end)
    end)
    |> Stream.filter(& &1 != nil)
    |> Enum.sort(& NaiveDateTime.compare(&1,&2) === :gt)
  end

end
