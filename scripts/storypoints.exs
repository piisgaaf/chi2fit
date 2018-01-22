defmodule StoryPoints do

  def parse_options(args) do
    options = case OptionParser.parse args, strict: [
      help: :boolean,
      sp: :integer,
      min: :float,
      max: :float] do
        {opt,[],[]} -> opt
        {_opt,_,[_rest|_]} -> usage(1)
    end
    if options[:help], do: usage(0)

    options
  end

  def validate_options(context) do
    # Validation of options
    cond do
      !context[:sp] || !context[:min] || !context[:max] ->
        IO.puts :stderr, "Error: please specify 'sp', 'min', and 'max'"
        System.halt(1)

      (context[:min] < 0) or (context[:min] > context[:sp]) ->
        IO.inspect context
        IO.puts :stderr, "Error: 'min' must be smaller than 'sp' and larger than 0"
        System.halt(1)

      context[:max] < context[:sp] ->
        IO.puts :stderr, "Error: 'max' must be larger than 'sp'"
        System.halt(1)

      true ->
        :ok
    end

    context
  end
  
  def usage code do
      IO.puts "Usage: #{__ENV__.file |> String.split("/") |> Enum.reverse |> hd} <options>"
      IO.puts "    --help\t\t\tShows this help"
      System.halt(code)
  end
end
  
##############################################################################################################################
##
## MAIN PROGRAM
##
##############################################################################################################################

import Chi2fit.Utilities
import Chi2fit.Distribution

## Parse and validate options
context = System.argv |> StoryPoints.parse_options |> StoryPoints.validate_options

func = fn k ->
  valmax = weibullCDF(k,context[:sp]/gamma(1+1.0/k)).(context[:max])
  valmin = weibullCDF(k,context[:sp]/gamma(1+1.0/k)).(context[:min])
  valmax-valmin-0.85
end

pdf = fn k,lambda,x ->
  k/lambda*:math.exp(-:math.pow(x/lambda,k))*:math.pow(x/lambda,k-1)
end

{a,b,_} = 0..100 |> Enum.reduce({0.1,10.0,nil}, fn
  _,{mn,mx,:done} -> {mn,mx,:done}
  n,{mn,mx,nil} -> {mn,mx,func.(0.1 + n*(10.0-0.1)/100) > 0}
  n,{mn,mx,sign} ->
    x = 0.1 + n*(10.0-0.1)/100
    val = func.(x)
    case {sign,val>0} do
      {true,true} -> {x,mx,true}
      {true,false} -> {mn,x,:done}
      {false,true} -> {mn,x,:done}
      {false,false} -> {x,mx,false}
    end
end)
{k,_,_} = newton(a,b,func,10, [])
lambda = context[:sp]/gamma(1+1/k)

IO.puts "k      = #{k}"
IO.puts "lambda = #{lambda}"

0..100 |> Enum.each(fn n->x=0.1+n*(5-0.1)/100; IO.puts "#{x},#{pdf.(k,lambda,x)}" end)
