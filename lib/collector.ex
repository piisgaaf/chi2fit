defmodule Chi2fit.Collector do

  # Copyright 2016-2017 Pieter Rijken
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

  use Agent, shutdown: 1_000, restart:

  @moduledoc """
  Implements an agent for collecting data.
  """

  def start_link do
    Agent.start_link(fn -> [] end, name: __MODULE__)
  end

  defp capture(out \\ []) do
    receive do
      {:io_request, from, refid, {:put_chars, :unicode, data}} ->
        send from, {:io_reply, refid, :ok}
        capture [data|out]

      {:stop, from} ->
        send from, {:data, Enum.reverse(out)}
        :ok

      :shutdown ->
        :ok

      msg ->
        raise ArgumentError, "[Collector] Received unexpected message: #{inspect msg}"
    end
  end

  def collect, do: Agent.cast(__MODULE__, fn _ -> capture() end)

  def value do
    send __MODULE__, {:stop, self()}
    receive do
      {:data, data} -> Agent.update(__MODULE__, fn _ -> data end)
    end
    Agent.get(__MODULE__, fn acc -> acc end)
  end

  def stop do
    send __MODULE__, :shutdown
    Agent.stop(__MODULE__)
  end

end
