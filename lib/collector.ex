defmodule Chi2fit.Utilities.Collector do
  use Agent

  def start_link() do
    Agent.start_link(fn -> [] end, name: __MODULE__)
  end

  defp capture(out \\ []) do
    receive do
      {:io_request, from, refid, {:put_chars, :unicode, data}} ->
        send from, {:io_reply, refid, :ok}
        capture [data|out]
          
      {:stop,from} ->
        send from, {:data, Enum.reverse(out)}
        :ok
        
      msg ->
        raise ArgumentError, "[Collector] Received unexpected message: #{inspect msg}"
    end
  end
  
  def collect, do: Agent.cast(__MODULE__, fn _ -> capture() end)

  def value do
    send __MODULE__, {:stop,self()}
    receive do
      {:data, data} -> Agent.update(__MODULE__, fn _ -> data end)
    end
    Agent.get(__MODULE__, fn acc -> acc end)
  end

  def stop, do: Agent.stop(__MODULE__)
end