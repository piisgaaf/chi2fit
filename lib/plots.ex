defmodule Gnuplotlib do
  
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

  alias Chi2fit.Utilities, as: U
  alias Gnuplot, as: G
  
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
    hist = data |> U.make_histogram(binsize,0) |> Enum.map(&Tuple.to_list/1)
    commands = [
      ['width=#{binsize}'],
      ['hist(x,width)=width*floor(x)+width/2.0'],
      [:set, :boxwidth, 'width*0.9'],
      [:set, :style, :fill, :solid, 0.5],
      if(options[:plottitle], do: [:set, :title, options[:plottitle]], else: []),
      if(options[:xrange], do: [:set, :xrange, options[:xrange]], else: []),
      if(options[:yrange], do: [:set, :yrange, options[:yrange]], else: []),
      if(options[:xlabel], do: [:set, :xlabel, options[:xlabel]], else: [:set, :xlabel]),
      if(options[:ylabel], do: [:set, :ylabel, options[:ylabel], :rotate, :by, 90], else: [:set, :ylabel]),

      [:plot, "-", :u, '(hist($1,width)):2', :smooth, :freq, :w, :boxes, :lc, 'rgb"green"', :notitle]
    ]
    if options[:commands], do: {commands,[hist]}, else: G.plot(commands,[hist])
  end

  @doc """
  Draws a graph of the empirical CDF as steps, the data points with error bars, and the (fitted) function.
  
  ## Options
    `:bin` - the size of the bins to use,
    `:plottitle` - the title of the plot,
    `:xrange` - the range for the x-values to use in the format '[x1:x2]'
    `:xrange` - the range for the y-values to use in the format '[y1:y2]'
    `:xlabel` - the label to use for the x-axis,
    `:ylabel` - the label to use for the y-axis,
    `:func` - the data to use for the CDF curve as a list of `[x,y]`,
    `:title` - the title to use for the CDF curve.

  """
  @type datapoint() :: {x :: number, y :: number, ylow :: number, yhigh :: number}
  @spec ecdf(data :: [datapoint], options :: Keyword.t) :: none()
  def ecdf(data, options) do
    npoints = options[:npoints] || 100
    hist = data |> Enum.map(&Tuple.to_list/1)
    maxx = data |> Enum.map(&elem(&1,0)) |> Enum.max |> Kernel.*(1.2)

    commands = [
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
            ])
        ]
    ]

    args = [
        [[0,0,0,0]|hist]++[[maxx,1,0,0]],
        hist,
        hist
    ] ++ if(options[:func], do: [dofun(npoints,maxx,options[:func])], else: [])

    if options[:commands], do: {commands,args}, else: G.plot(commands,args)
  end

  @doc """
  Draws a graph of the PDF.
  
  ## Options
    `:bin` - the size of the bins to use,
    `:plottitle` - the title of the plot,
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
    maxx = data |> Enum.max |> Kernel.*(1.2)
    
    hist = data
    |> U.make_histogram(bin,0)
    |> Enum.map(&Tuple.to_list/1)
    |> Enum.map(fn [x,y]->[(x-0.5)*bin,y] end)

    commands = [
        ['count=#{length(data)}'],
        ['width=#{bin}'],
        ['hist(x,width)=width*floor(x/width)+width/2.0'],
        [:set, :boxwidth, 'width*0.9'],
        [:set, :style, :fill, :solid, 0.5],
        if(options[:plottitle], do: [:set, :title, options[:plottitle]], else: []),
        if(options[:xrange], do: [:set, :xrange, options[:xrange]], else: []),
        if(options[:yrange], do: [:set, :yrange, options[:yrange]], else: []),
        if(options[:xlabel], do: [:set, :xlabel, options[:xlabel]], else: [:set,:xlabel]),
        if(options[:ylabel], do: [:set, :ylabel, options[:ylabel], :rotate, :by, 90], else: [:set,:ylabel]),

        [:plot, G.list([
                ~w|'-' u (hist($1,width)):($2/count/#{bin}) smooth freq w boxes lc rgb "green" title "Empirical PDF"|a,
                ~w|'-' u (hist($1,width)):($2/count/#{bin}):(sqrt($2)/count/#{bin}) w errorbars ls 3 notitle|a,
                ["", :u, '1:2', :w, :lines, :ls, 3, :title, options[:title]]
            ])
        ]
    ]
    
    args = [ hist,hist,dofun(npoints,maxx,options[:pdf]) ]
    
    if options[:commands], do: {commands,args}, else: G.plot(commands,args)
  end
  
  @doc """
  Plots severals graphs in a multi-plot.
  """
  @spec multi(all :: [{command :: [],data :: []}], options :: Keyword.t) :: none()
  def multi(all, options \\ []) do
    cols = options[:columns] || 2
    rows = options[:rows] || trunc(Float.ceil(length(all)/cols,0))
    {commands,data} = all |> U.unzip
    G.plot([ [
            [:set, :terminal, :wxt, :size, '1000,500'],
            [:set, :multiplot, :layout, '#{rows},#{cols}'] ++ (if options[:title], do: [:title, options[:title], :font, ",14"], else: []),
        ] ] |> Enum.concat(commands) |> Enum.concat,
        data |> Enum.concat)
  end

  #############################################################################
  ##
  ## Local functions
  ##
  #############################################################################
  
  defp dofun(npoints,maxx,fun) do
      0..npoints
      |> Enum.map(fn i -> [i*maxx/npoints,fun.(i*maxx/npoints)] end)
  end

end