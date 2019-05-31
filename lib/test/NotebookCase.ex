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

      Module.register_attribute __MODULE__, :out, accumulate: true
      Module.put_attribute __MODULE__, :dir, unquote(options)[:dir] || "."

      @template ~S"""
      defmodule <%= String.replace("Runner_#{name}",~r<[-./]>,"_") %> do
      use NotebookUnit.Case

      <%= for {"cells",cells} <- json do %>
        <%= for cell <- cells, cell["cell_type"]=="code" do %>
      # In[<%= cell["execution_count"] %>]
      runblock <%= cell["execution_count"] %> do
      <%= cell["source"] |> Enum.join("") |> String.replace("/app/notebooks/data", "notebooks/data") %>
      end
        <% end %>
      <% end %>
      
      @out
      end
      """

      require EEx
      EEx.function_from_string :defp, :convert, @template, [:json,:name]
    end
  end
  
  defmacro runblock(key, do: block) do
    quote do
      @out {unquote(key), unquote(block)}
    end
  end

  defmacro nbrun(json, notebook) do
    quote do
      code = unquote(json) |> convert(unquote(notebook))
      capture_io :stderr, fn ->
        capture_io fn ->
          try do
            Code.compiler_options(warnings_as_errors: true)
            {{:module,_mod,_data,out}, _} = Code.eval_string code
            send self(), {:execute, out}
          rescue
            error -> send self(), {:error, error}
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
      test "Notebook - #{unquote notebook}" do
        json = (@dir <> "/" <> unquote notebook) |> File.read!() |> Poison.decode!()

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
