defmodule Distribution.Nakagami do

  @moduledoc """
  The Nakagami distribution.
  """

  defstruct [:pars]
  
  @type t() :: %__MODULE__{
    pars: [number()] | nil
  }

end

defimpl Distribution, for: Distribution.Nakagami do
  import Distribution.Nakagami
  alias Distribution.Nakagami

  import Exboost.Math, only: [tgamma: 1]

  @spec nakagami(scale::number(),shape::number()) :: ((...) -> number)
  defp nakagami(scale,shape) do
    fn ->
      u = :rand.uniform()
      scale*:math.sqrt(Exboost.Math.gamma_p_inv(shape,u)/shape)
    end
  end

  @spec nakagamiCDF(scale :: float,shape :: float) :: (number -> number)
  defp nakagamiCDF(scale,shape) when scale>0 and shape>0 do
    fn
      x ->
        Exboost.Math.tgamma_lower(shape,shape*(x/scale)*(x/scale))
    end
  end
  defp nakagamiCDF(_scale,_shape), do: raise ArithmeticError, "Nakagami is only defined for positive scale and shape"

  def skewness(%Nakagami{pars: nil}) do
    fn [_scale,shape] ->
      g = tgamma(shape)
      g1_2 = tgamma(shape+0.5)
      g1 = tgamma(shape+1.0)
      g3_2 = tgamma(shape+1.5)
      num = 2*g1_2*g1_2*g1_2 + g*g*( g3_2 - 3*shape*g1_2 )
      den = g*g*g*:math.pow(shape*(1.0-shape*g1_2*g1_2/g1/g1),1.5)
      num/den
    end
  end
  def kurtosis(%Nakagami{pars: nil}) do
    fn [_scale,shape] ->
      g = tgamma(shape)
      g1_2 = tgamma(shape+0.5)
      g2 = tgamma(shape+2.0)
      gdouble = tgamma(2*shape)
      num = -6*g1_2*g1_2*g1_2*g1_2 - 3*shape*shape*g*g*g*g + g*g*g*g2 + :math.pow(2,3-4*shape)*(4*shape-1)*:math.pi*gdouble*gdouble
      den = :math.pow(abs(g1_2*g1_2 - shape*g*g),2)
      num/den
    end
  end
  def size(%Nakagami{}), do: 2
  
  def cdf(%Nakagami{pars: nil}), do: fn x, [scale,shape] -> nakagamiCDF(scale,shape).(x) end

  def pdf(%Nakagami{pars: nil}), do: raise Distribution.FunctionNotSupportedError, message: "pdf is not supported for the Nakagami distribution"
  def random(%Nakagami{pars: [scale,shape]}), do: nakagami(scale, shape).()
end
