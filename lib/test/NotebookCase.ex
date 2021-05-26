defmodule NotebookUnit.Case do

  # Copyright 2019 Pieter Rijken
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

  defmacro __using__(options) do
    quote do
      import unquote(__MODULE__)
      import ExUnit.CaptureIO

      @compile warnings_as_errors: true
      @compile ignore_module_conflict: true

      Module.put_attribute __MODULE__, :nbdir, unquote(options)[:dir] || "."
    end
  end
  
  defmacro nbrun(json, _notebook) do
    quote do
      capture_io :stderr, fn ->
        capture_io fn ->
          try do
            {env, _} = Code.eval_quoted(quote do
              import IEx.Helpers
              __ENV__
            end)

            {:ok, _pid} = Boyle.start_link([])
            result = for {"cells",cells} <- unquote(json) do
              for cell <- cells, cell["cell_type"]=="code" do
                code = cell["source"]
                |> Enum.join("")
                |> String.replace("/app/notebooks/data", "notebooks/data")
                |> String.replace("/app/notebooks/images", "notebooks/images")
                {cell["execution_count"], code}
              end
            end
            |> List.flatten()
            |> Enum.reduce([out: [], binding: [ans: nil, out: %{}], env: env], fn {key,code}, acc ->
              block = { Code.string_to_quoted(code), quote(do: __ENV__) }
              {{result,env}, binding} = Code.eval_quoted(block, acc[:binding], acc[:env])
              new_binding = case result do
                :"do not show this result in output" -> binding
                :"this is an inline image" -> binding
                {:"this is an inline image",_} -> binding
                _ ->
                  binding
                  |> Keyword.update!(:ans,fn _ -> result end)
                  |> Keyword.update!(:out,& Map.put_new(&1,key,result))
              end
              [binding: new_binding, out: [{key,result}|acc[:out]], env: env]
            end)
            send self(), {:execute, result[:out]}
          rescue
            error -> send self(), {:error, error}
          after
            :ok = GenServer.stop Boyle
          end
        end
      end
    end
  end

  @mdheader ~r|#+ .*|
  def md_headers(json) do
    for {"cells",cells} <- json do
      for cell <- cells, cell["cell_type"]=="markdown" do
        cell["source"] |> Enum.filter(& Regex.match?(@mdheader,&1))
      end |> List.flatten()
    end |> List.flatten()
    |> Enum.map(& String.replace_leading(&1,"#",""))
    |> Enum.map(&String.trim/1)
  end

  @toc ~r|\#+ Table [oO]f [cC]ontents ?.*|
  @mdlink ~r|.*\[(?<hdr>.+)\]\(#(?<link>.+)\).*|
  @mdlistitem ~r| *\* +.*|
  def md_toc(json) do
    for {"cells",cells} <- json do
      for cell <- cells, cell["cell_type"]=="markdown", Regex.match?(@toc,hd(cell["source"])) do
        tl(cell["source"])
      end |> List.flatten()
    end |> List.flatten()
    |> Enum.filter(& Regex.match?(@mdlistitem,&1))
    |> Enum.map(& Regex.named_captures(@mdlink,&1))
  end

  defmacro nbtest(notebook) do
    quote do
      @tag notebooks: true
      @tag :notest
      test "Notebook - #{unquote notebook}" do
        json = (@nbdir <> "/" <> unquote notebook) |> File.read!() |> Poison.decode!()

        warns = nbrun json, unquote(notebook)
        assert_received {:execute, _out}
        assert warns == ""
        
        headers = md_headers(json)
        toc = md_toc(json)

        Enum.each(toc, fn %{"hdr" => hdr, "link" => link} ->
          assert hdr in headers, "Table of contents entry '#{hdr}' is missing"
          assert link in Enum.map(headers,& String.replace(&1," ","-")), "Table of contents '#{link}' is missing a corresponding header"
        end)
      end
    end
  end

end
