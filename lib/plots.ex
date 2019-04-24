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
  Provides various various plots using ghe Gnuplot module.  
  """

  alias Chi2fit.Utilities, as: U
  alias Gnuplot, as: G
  
  def histogram(data, options) do
      hist = data |> U.make_histogram(1,0) |> Enum.map(&Tuple.to_list/1)
      commands = [
          ['width=1.'],
          ['hist(x,width)=width*floor(x/width)+width/2.0'],
          [:set, :boxwidth, 'width*0.9'],
          [:set, :style, :fill, :solid, 0.5],
          if(options[:plottitle], do: [:set, :title, options[:plottitle]], else: []),
          if(options[:xrange], do: [:set, :xrange, options[:xrange]], else: []),
          if(options[:yrange], do: [:set, :yrange, options[:yrange]], else: []),
          if(options[:xlabel], do: [:set, :xlabel, options[:xlabel]], else: []),
          if(options[:ylabel], do: [:set, :ylabel, options[:ylabel], :rotate, :by, 90], else: [])
      ]
      G.plot(commands ++
          [[:plot, "-", :u, '(hist($1,width)):2', :smooth, :freq, :w, :boxes, :lc, 'rgb"green"', :notitle]],
          [hist])
  end

  def ecdf(data, options) do
    hist = data|>Enum.map(&Tuple.to_list/1)
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
        [:set, :yrange, '[0:1.2]'],
        if(options[:xlabel], do: [:set, :xlabel, options[:xlabel]], else: []),
        if(options[:ylabel], do: [:set, :ylabel, options[:ylabel], :rotate, :by, 90], else: [])
    ]
    G.plot(commands ++
      [
        [:plot, G.list([
                ~w('-' u 1:2 w steps ls 1 notitle)a,
                ~w('' u 1:2 w points ls 1 notitle)a,
                ~w('' u 1:2:3:4 w yerrorbars ls 2 title 'Empirical CDF')a,
                if(options[:func], do: ["", :u, '1:2', :w, :lines, :ls, 3, :title, options[:title]], else: [])
            ])
        ]
      ],
      [
          [[0,0,0,0]|hist]++[[maxx,1,0,0]],
          hist,
          hist
      ] ++ if(options[:func], do: [options[:func]], else: []))
  end

  def pdf(data, options) do
      bin = options[:bin] || 1
      hist = data
      |> U.make_histogram(bin,0)
      |> Enum.map(&Tuple.to_list/1)
      |> Enum.map(fn [x,y]->[(x-0.5)*bin,y] end)
      commands = [
          ['count=#{length(data)}'],
          ['width=#{options[:bin]||1}.'],
          ['hist(x,width)=width*floor(x/width)+width/2.0'],
          [:set, :boxwidth, 'width*0.9'],
          [:set, :style, :fill, :solid, 0.5],
          if(options[:plottitle], do: [:set, :title, options[:plottitle]], else: []),
          if(options[:yrange], do: [:set, :yrange, options[:yrange]], else: []),
          if(options[:xlabel], do: [:set, :xlabel, options[:xlabel]], else: []),
          if(options[:ylabel], do: [:set, :ylabel, options[:ylabel], :rotate, :by, 90], else: [])
      ]
      G.plot(commands ++
        [
          [:plot, G.list([
                  ~w|'-' u (hist($1,width)):($2/count/#{bin}) smooth freq w boxes lc rgb "green" title "Empirical PDF"|a,
                  ~w|'-' u (hist($1,width)):($2/count/#{bin}):(sqrt($2)/count/#{bin}) w errorbars ls 3 notitle|a,
                  ["", :u, '1:2', :w, :lines, :ls, 3, :title, options[:title]]
              ])
          ]
        ],
        [ hist,hist,options[:pdf] ])
  end
    
end