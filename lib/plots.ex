defmodule Chi2fit.Gnuplotlib do

  # Copyright 2019-2019 Pieter Rijken
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
  Provides various various plots using the [Gnuplot](https://hex.pm/packages/gnuplot) package.
  """

  alias Gnuplot, as: G
  alias Chi2fit.Statistics, as: S
  alias Chi2fit.Utilities, as: U
  
  @imgpath "/app/notebooks/images"
  @terminal "pngcairo"
  @pngoptions ~w(set terminal #{@terminal} transparent enhanced)a

  @doc """
  Captures the output sent by `&Port.open/2` and returns it as a binary
  """
  @timeout 1_000
  @spec capture(out :: binary) :: binary
  def capture(out \\ <<>>) do
    receive do
      {_, {:data, data}} ->
        capture(out <> data)
      {_, :closed} ->
        out
    after
      @timeout ->
        out
    end
  end

  @doc """
  Draws a histogram of the data.
  
  ## Options

      `:bin` - the size of the bins to use,
      `:plottitle` - the title of the plot,
      `:xrange` - the range for the x-values to use in the format '[x1:x2]'
      `:xrange` - the range for the y-values to use in the format '[y1:y2]'
      `:xlabel` - the label to use for the x-axis,
      `:ylabel` - the label to use for the y-axis.

  """
  @spec histogram(data :: [number], options :: Keyword.t) :: none()
  def histogram(data, options \\ []) do
    binsize = options[:bin] || 1
    hist = data |> S.make_histogram(binsize,0) |> Enum.map(&Tuple.to_list/1)

    terminal(options)
      ++ [
            ['width=#{binsize}'],
            ['hist(x,width)=width*floor(x-1)+width/2.0'],
            [:set, :boxwidth, 'width*0.9'],
            [:set, :style, :fill, :solid, 0.5],
            if(options[:plottitle], do: [:set, :title, options[:plottitle]], else: []),
            if(options[:xrange], do: [:set, :xrange, options[:xrange]], else: []),
            if(options[:yrange], do: [:set, :yrange, options[:yrange]], else: []),
            if(options[:xlabel], do: [:set, :xlabel, options[:xlabel]], else: [:set, :xlabel]),
            if(options[:ylabel], do: [:set, :ylabel, options[:ylabel], :rotate, :by, 90], else: [:set, :ylabel]),

            [:plot, "-", :u, '(hist($1,width)):2', :smooth, :freq, :w, :boxes, :lc, 'rgb"green"', :notitle]
          ]
      |> do_output([hist], options)
  end

  @doc """
  Draws a graph of the empirical CDF as steps, the data points with error bars, and the (fitted) function.
  
  ## Options
  
      `:bin` - the size of the bins to use,
      `:plottitle` - the title of the plot,
      `:xrange`- the range for the x-values to use in the format '[x1:x2]'
      `:xrange` - the range for the y-values to use in the format '[y1:y2]'
      `:xlabel` - the label to use for the x-axis,
      `:ylabel` - the label to use for the y-axis,
      `:func` - the data to use for the CDF curve as a list of `[x,y]`,
      `:title` - the title to use for the CDF curve.
      `:bounds` - 2-tuple of functions describing the minimum and maximum error-curves for the CDF

  """
  @type datapoint() :: {x :: number, y :: number, ylow :: number, yhigh :: number}
  @spec ecdf(data :: [datapoint], options :: Keyword.t) :: none()
  def ecdf(data, options) do
    npoints = options[:npoints] || 100
    hist = data |> Enum.map(&Tuple.to_list/1)
    maxx = data |> Enum.map(&elem(&1,0)) |> Enum.max |> Kernel.*(1.2)

    args = [
        [[0,0,0,0]|hist]++[[maxx,1,0,0]],
        hist,
        hist
    ]
    ++ if(options[:func], do: [dofun(npoints,maxx,options[:func])], else: [])
    ++ case options[:bounds] do
         {minrate,maxrate} ->
           [ dofun(npoints,maxx,minrate) |> Enum.zip(dofun(npoints,maxx,maxrate)) |> Enum.map(fn {[x,y1],[x,y2]}->[x,y1,y2] end) ]
         _else ->
           [ ]
        end

    terminal(options)
      ++ [
        [:set, :style, :line, 1,
            :linecolor, :rgb, "#0060ad",
            :linetype, 1, :linewidth, 2,
            :pointtype, 7, :pointsize, 1.5],
        [:set, :style, :line, 2,
            :linecolor, :rgb, "#dd181f",
            :linetype, 1, :linewidth, 2],
        [:set, :style, :line, 3,
            :linecolor, :rgb, "green",
            :linetype, 1, :linewidth, 2],
        ~w(set style fill transparent solid 0.2 noborder)a,
        ~w(set key left top)a,
        if(options[:plottitle], do: [:set, :title, options[:plottitle]], else: []),
        if(options[:xrange], do: [:set, :xrange, options[:xrange]], else: []),
        if(options[:yrange], do: [:set, :yrange, options[:yrange]], else: [:set,:yrange,'[0:1.2]']),
        if(options[:xlabel], do: [:set, :xlabel, options[:xlabel]], else: [:set, :xlabel]),
        if(options[:ylabel], do: [:set, :ylabel, options[:ylabel], :rotate, :by, 90], else: [:set, :ylabel]),

        [:plot, G.list([
                ~w('-' u 1:2 w steps ls 1 notitle)a,
                ~w('' u 1:2 w points ls 1 notitle)a,
                ~w('' u 1:2:3:4 w yerrorbars ls 2 title 'Empirical CDF')a,
                if(options[:func], do: ["", :u, '1:2', :w, :lines, :ls, 3, :title, options[:title]], else: [])
            ] ++ case options[:bounds] do
                  {_,_} -> [
                    ["", :u, '1:2:3', :lc, :rgb, "grey", :w, :filledcurve, :closed, :title, "Error bounds"]
                  ]
                  _else -> []
                end
            )
        ]
      ]
      |> do_output(args, options)
  end

  @doc """
  Draws a graph of the PDF.
  
  ## Options

      `:bin` - the size of the bins to use,
      `:offset` -- the offset of the bin,
      `:plottitle` - the title of the plot,
      `:noerror` -- no error bars,
      `:xrange` - the range for the x-values to use in the format '[x1:x2]'
      `:xrange` - the range for the y-values to use in the format '[y1:y2]'
      `:xlabel` - the label to use for the x-axis,
      `:ylabel` - the label to use for the y-axis,
      `:pdf` - the data to use for the PDF curve as a list of `[x,y]`,
      `:title` - the title to use for the PDF curve.

  """
  @spec pdf(data :: [number], options :: Keyword.t) :: none()
  def pdf(data, options) do
    npoints = options[:npoints] || 100
    bin = options[:bin] || 1
    offset = options[:offset] || 0
    noerrors = options[:noerrors] || false
    maxx = data |> Enum.max |> Kernel.*(1.2)
    
    hist = data
    |> S.make_histogram(bin,offset)
    |> Enum.map(&Tuple.to_list/1)
    |> Enum.map(fn [x,y]->[x*bin,y] end)

    args = [ hist ] ++ if(noerrors, do: [], else: [hist]) ++ [ dofun(npoints,maxx,options[:pdf]) ]

    terminal(options)
      ++ [
        ['count=#{length(data)}'],
        ['width=#{bin}'],
        ['hist(x,width)=width*floor((x-1)/width)+width/2.0'],
        [:set, :boxwidth, 'width*0.9'],
        [:set, :style, :fill, :solid, 0.5],
        if(options[:plottitle], do: [:set, :title, options[:plottitle]], else: []),
        if(options[:xrange], do: [:set, :xrange, options[:xrange]], else: []),
        if(options[:yrange], do: [:set, :yrange, options[:yrange]], else: []),
        if(options[:xlabel], do: [:set, :xlabel, options[:xlabel]], else: [:set,:xlabel]),
        if(options[:ylabel], do: [:set, :ylabel, options[:ylabel], :rotate, :by, 90], else: [:set,:ylabel]),

        [:plot, G.list(
                [~w|'-' u (hist($1,width)):($2/count/#{bin}) smooth freq w boxes lc rgb "green" title "Empirical PDF"|a]
                ++ if(noerrors, do: [], else: [~w|'-' u (hist($1,width)):($2/count/#{bin}):(sqrt($2)/count/#{bin}) w errorbars ls 3 notitle|a])
                ++ [ ["", :u, '1:2', :w, :lines, :ls, 3, :title, options[:title]] ]
            )
        ]
      ]
      |> do_output(args, options)
  end
  
  @doc """
  Plots severals graphs in a multi-plot.
  """
  @spec multi(all :: [{command :: [],data :: []}], options :: Keyword.t) :: none()
  def multi(all, options \\ []) do
    cols = options[:columns] || 2
    rows = options[:rows] || trunc(Float.ceil(length(all)/cols,0))
    {commands,data} = all |> U.unzip

    [
      terminal(options)
      ++ [
           [:set, :multiplot, :layout, '#{rows},#{cols}'] ++ (if options[:title], do: [:title, options[:title], :font, ",14"], else: []),
       ]
    ]
    |> Enum.concat(commands)
    |> Enum.concat
    |> do_output(Enum.concat(data), options)
  end

  @spec surface(data :: [[number()]], options :: Keyword.t) :: none()
  def surface(data, options) do
    dgrid3d = options[:dgrid3d] || ""
    using = options[:parameters] || "1:2"
    
    terminal(options)
      ++ [
        ['set cntrparam levels 50'],
        ['set contour'],
        ['set dgrid3d #{dgrid3d}'],
        ['unset key'],
        ['set cntrlabel start 2 font ",7"'],

        if(options[:plottitle], do: [:set, :title, options[:plottitle]], else: []),
        if(options[:xrange], do: [:set, :xrange, options[:xrange]], else: []),
        if(options[:yrange], do: [:set, :yrange, options[:yrange]], else: []),
        if(options[:zrange], do: [:set, :zrange, options[:zrange]], else: []),
        if(options[:xlabel], do: [:set, :xlabel, options[:xlabel]], else: [:set,:xlabel]),
        if(options[:ylabel], do: [:set, :ylabel, options[:ylabel]], else: [:set,:ylabel]),
        if(options[:zlabel], do: [:set, :zlabel, options[:zlabel], :rotate, :by, 90], else: [:set,:zlabel]),

        ~w(splot '-' u #{using}:3 w lines notitle)a
      ]
      |> do_output([data], options)
  end

  #############################################################################
  ##
  ## Local functions
  ##
  #############################################################################
  
  defp make_terminal(options) do
    List.flatten([
      @pngoptions,
      if(options[:size], do: ~w(size #{options[:size]})a, else: [ ])
    ])
  end

  defp terminal(options) do
    case options[:mode] do
      {:as_file, path} ->
        File.mkdir_p @imgpath
        [ make_terminal(options), ~w(set output '#{@imgpath}/#{path}.png')a ]
      raw when raw==nil or raw==:raw or raw==:raw64 ->
        [ make_terminal(options), ~w(set output)a ]
      _else ->
        [ ]
    end
  end

  defp do_output(commands, datasets, options) do
    case options[:mode] do
      :as_commands ->
        {commands, datasets}
      {:as_file, _} ->
        G.plot(commands, datasets)
      :as_raw ->
        G.plot(commands, datasets)
        capture() |> IO.write
      :as_raw64 ->
        G.plot(commands, datasets)
        capture() |> Base.encode64 |> IO.write
      nil ->
        G.plot(commands, datasets)
        capture() |> Base.encode64 |> IO.write
    end
  end

  defp dofun(npoints,maxx,fun) do
      0..npoints
      |> Enum.map(fn i -> [i*maxx/npoints,fun.(i*maxx/npoints)] end)
  end

end