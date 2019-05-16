defmodule Distribution.Normal do

  @moduledoc """
  The normal or Gauss distribution
  """

  defstruct [:pars]
  
  @type t() :: %__MODULE__{
    pars: [number()] | nil
  }

end

defimpl Distribution, for: Distribution.Normal do
  import Distribution.Normal
  alias Distribution.Normal

  defp normal(mean,sigma) when is_number(mean) and is_number(sigma) and sigma>=0 do
    fn () ->
      {w,v1,_} = polar()
      y = v1*:math.sqrt(-2*:math.log(w)/w)
      mean + sigma*y
    end
  end

  defp normalCDF(_mean,sigma) when sigma<0, do: raise ArgumentError
  defp normalCDF(mean,sigma) when is_number(mean) and is_number(sigma) and sigma>=0 do
    fn
      x when (x-mean)/sigma < 4.0 ->
        0.5*:math.erfc(-(x-mean)/sigma/:math.sqrt(2.0)) 
      x ->
        0.5*( 1.0 + :math.erf((x-mean)/sigma/:math.sqrt(2.0)) )
    end
  end

  @spec polar() :: {number(), number(), number()}
  defp polar() do
    v1 = random(-1,1)
    v2 = random(-1,1)
    w = v1*v1 + v2*v2

    cond do
      w >= 1.0 -> polar()
      true -> {w,v1,v2}
    end
  end

  @spec random(min::number(),max::number()) :: number()
  defp random(min,max) when max >= min do
    min + (max-min)*:rand.uniform()
  end
  
  def skewness(%Normal{pars: nil}), do: fn _ -> 0 end
  def kurtosis(%Normal{pars: nil}), do: fn _ -> 0 end
  def size(%Normal{}), do: 2
  
  def cdf(%Normal{pars: nil}), do: fn x, [mu,sigma] -> normalCDF(mu,sigma).(x) end
  def cdf(%Normal{pars: [mu,sigma]}), do: normalCDF(mu,sigma)

  def pdf(%Normal{pars: nil}), do: fn x, [mu,sigma] -> 1/:math.sqrt(2*:math.pi)/sigma*:math.exp(-(x-mu)*(x-mu)/2/sigma/sigma) end

  def random(%Normal{pars: [mu,sigma]}), do: normal(mu,sigma).()
end
