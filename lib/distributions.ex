defmodule Chi2fit.Distribution do

  # Copyright 2012-2017 Pieter Rijken
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
  Provides various distributions.
  """

  import Chi2fit.Utilities
  
  @type distribution() :: ((...) :: term())
  @type cdf() :: ((number) :: number())

  defmodule UnsupportedDistributionError do
    defexception message: "Unsupported distribution function"
  end

  ###
  ### Standard distributions
  ###

  @doc """
  Uniform distribution.
  """
  @spec uniform(Keyword.t) :: distribution
  def uniform([]), do: uniform(0, 2.0)
  def uniform([avg: average]), do: uniform(0,2*average)
  def uniform(list) when is_list(list), do: fn () -> Enum.random(list) end

  @doc """
  Uniform distribution.
  """
  @spec uniform(min::integer(),max::integer()) :: distribution
  def uniform(min,max) when max>=min, do: fn () -> random(min,max) end

  @doc """
  Constant distribution.
  """
  @spec constant(number | Keyword.t) :: distribution
  def constant([avg: average]), do: fn () -> average end
  def constant(average) when is_number(average), do: fn () -> average end

  @doc """
  The exponential distribution.
  """
  @spec exponential(Keyword.t) :: distribution
  def exponential([avg: average]) do
    fn ->
      u = :rand.uniform()
      -average*:math.log(u)
    end
  end
  def exponential(rate), do: exponential([avg: 1.0/rate])
  def exponentialCDF(rate) when rate > 0.0, do: fn t -> 1.0 - :math.exp(-rate*t) end

  @step 500
  @expstep :math.exp(@step)
  defp iterate(p,r) when p<1 and r>0 and r>@step, do: iterate(p*@expstep,r-@step)
  defp iterate(p,r) when p<1 and r>0, do: iterate(p*:math.exp(r),0)
  defp iterate(p,r), do: {p,r}

  defp _poisson(rate, k \\ 0, p \\ 1.0)
  defp _poisson(rate,k,p) do
    k = k+1
    p = p*:rand.uniform()
    {p,rate} = iterate p,rate
    if p>1, do: _poisson(rate,k,p), else: k-1
  end

  @doc """
  The Poisson distribution.
  
  For the implementation, see https://en.wikipedia.org/wiki/Poisson_distribution, 'Generating Poisson-distributed random variables'
  """
  def poisson(rate), do: fn -> _poisson(rate) end
  def poissonCDF(rate) when rate > 0.0 do
    fn t -> 1.0 - Exboost.Math.gamma_p(Float.floor(t+1.0),rate) end
  end

  @doc """
  The Erlang distribution.
  """
  @spec erlang(k::integer(),lambda::number()) :: distribution
  def erlang(k,lambda) when is_integer(k) and k>0 do
    fn ->
      -1.0/lambda*:math.log(1..k |> Enum.reduce(1.0, fn _,acc -> :rand.uniform()*acc end))
    end
  end
  
  @doc """
  The Erlang cumulative distribution function.
  """
  @spec erlangCDF(k::number(),lambda::number()) :: cdf
  def erlangCDF(k,lambda) when k<0 or lambda<0, do: raise ArithmeticError, "Erlang is only defined for positive shape and mode"
  def erlangCDF(k,lambda) when k>0, do: &Exboost.Math.gamma_p(k,lambda*&1)

  @gamma53 0.902745292950933611297
  @gamma32 0.886226925452758013649

  @doc """
  The Weibull distribution.
  """
  @spec weibull(number, number|Keyword.t) :: distribution
  def weibull(1.0, [avg: average]), do: weibull(1.0, average)
  def weibull(1.5, [avg: average]), do: weibull(1.5, average/@gamma53)
  def weibull(2.0, [avg: average]), do: weibull(2.0, average/@gamma32)
  def weibull(alpha, beta) when is_number(alpha) and is_number(beta) do
    fn ->
      u = :rand.uniform()
      beta*:math.pow(-:math.log(u),1.0/alpha)
    end
  end
  
  @doc """
  The Weibull cumulative distribution function.
  """
  @spec weibullCDF(number,number) :: cdf
  def weibullCDF(k,_) when k<0, do: raise ArithmeticError, "Weibull is only defined for positive shape"
  def weibullCDF(_,lambda) when lambda<0, do: raise ArithmeticError, "Weibull is only defined for positive scale"
  def weibullCDF(k,lambda) when is_number(k) and is_number(lambda) do
    fn
      0 -> 0.0
      0.0 -> 0.0
      x when x<0 -> 0.0
      x ->
        lg = :math.log(x/lambda)*k
        cond do
          lg > 100.0 ->
            0.0
          lg < -18.0 ->
            ## With -18 (x/lambda)^2k < 10^(-16)
            t = :math.pow(x/lambda,k)
            t*(1 - 0.5*t)
          true ->
            1.0 - :math.exp -:math.pow(x/lambda,k)
        end
    end
  end

  @doc """
  The normal or Gauss distribution
  """
  @spec normal(mean::number(),sigma::number()) :: distribution()
  def normal(mean,sigma) when is_number(mean) and is_number(sigma) and sigma>=0 do
    fn () ->
      {w,v1,_} = polar()
      y = v1*:math.sqrt(-2*:math.log(w)/w)
      mean + sigma*y
    end
  end

  @doc """
  The normal or Gauss cumulative distribution
  """
  @spec normalCDF(mean::number(),sigma::number()) :: cdf
  def normalCDF(_mean,sigma) when sigma<0, do: raise ArgumentError
  def normalCDF(mean,sigma) when is_number(mean) and is_number(sigma) and sigma>=0 do
    fn
      x when (x-mean)/sigma < 4.0 ->
        0.5*:math.erfc(-(x-mean)/sigma/:math.sqrt(2.0)) 
      x ->
        0.5*( 1.0 + :math.erf((x-mean)/sigma/:math.sqrt(2.0)) )
    end
  end

  @doc """
  The Bernoulli distribution.
  """
  @spec bernoulli(value :: number) :: distribution
  def bernoulli(value) when is_number(value) do
   fn () ->
       u = :rand.uniform()
       if u <= value, do: 1, else: 0
   end
  end

  @doc """
  Wald or Inverse Gauss distribution.
  """
  @spec wald(mu::number(),lambda::number()) :: distribution
  def wald(mu,lambda) when is_number(mu) and is_number(lambda) do
   fn () ->
       w = normal(0.0,1.0).()
       y = w*w
       x = mu + mu*mu*y/2/lambda - mu/2/lambda*:math.sqrt(4*mu*lambda*y + mu*mu*y*y)

       z = :rand.uniform()
       if z <= mu/(mu+x), do: x, else: mu*mu/x
   end
  end
  def wald([avg: average],lambda), do: wald(average,lambda)

  @doc """
  The Wald (Inverse Gauss) cumulative distribution function.
  """
  @spec waldCDF(number,number) :: cdf
  def waldCDF(mu,_) when mu < 0, do: raise ArithmeticError, "Wald is only defined for positive average"
  def waldCDF(_,lambda) when lambda < 0, do: raise ArithmeticError, "Wald is only defined for positive shape"
  def waldCDF(mu,lambda) do
    fn
      x when x == 0 -> 0.0
      x when x < 0 -> 0.0
      x when x > 0 ->
        phi(:math.sqrt(lambda/x) * (x/mu-1.0)) + :math.exp(2.0*lambda/mu) * phi(-:math.sqrt(lambda/x) * (x/mu+1.0))
    end
  end

  @doc """
  The Skew Exponential Power cumulative distribution (Azzalini).

  ## Options
    `:method` - the integration method to use, :gauss and :romberg types are supported, see below
    `:tolerance` - re-iterate until the tolerance is reached (only for :romberg)
    `:points` - the number of points to use in :gauss method

  ## Integration methods
    `:gauss` - n-point Gauss rule,
    `:gauss2` - n-point Guass rule with tanh transformation,
    `:gauss3` - n-point Gauss rule with linear transformstion,
    `:romberg` - Romberg integration,
    `:romberg2` - Romberg integration with tanh transformation,
    `:romberg3` - Romberg integration with linear transformstion.

  """
  @spec sepCDF(a :: float,b :: float,lambda :: float,alpha :: float, options :: Keyword.t) :: cdf
  def sepCDF(a,b,lambda,alpha,options \\ []) do
    method = options[:method] || :romberg2
    endpoint = if method in [:gauss2,:gauss3,:romberg2,:romberg3], do: :infinity, else: 1000.0
    fn
      x ->
        result2 = integrate(method, sepPDF(a,b,lambda,alpha), 0.0, x, options)
        result3 = integrate(method, sepPDF(a,b,lambda,alpha), 0.0, endpoint, options)
        result2/result3
    end
  end

  @doc """
  Fréhet or inverse Weibull distribution.
  """
  @spec frechet(scale::number(),shape::number()) :: distribution
  def frechet(scale,shape) when is_number(scale) and is_number(shape) do
   fn ->
     u = :rand.uniform()
     scale * :math.pow(-:math.log(u),-1.0/shape)
   end
  end

  @doc """
  The Fréchet distribution, also known inverse Weibull distribution.
  """
  @spec frechetCDF(scale :: float,shape :: float) :: cdf
  def frechetCDF(scale,shape) when scale>0 and shape>0 do
    fn
      x when x==0.0 ->
        0.0
      x ->
        :math.exp(-:math.pow(x/scale,-shape))
    end
  end
  def frechetCDF(_scale,_shape), do: raise ArithmeticError, "Fréchet is only defined for positive scale and shape"

  @doc """
  Nakagami distribution.
  """
  @spec nakagami(scale::number(),shape::number()) :: distribution
  def nakagami(scale,shape) do
    fn ->
      u = :rand.uniform()
      scale*:math.sqrt(Exboost.Math.gamma_p_inv(shape,u)/shape)
    end
  end

  @doc """
  The Nakagami distribution.
  """
  @spec nakagamiCDF(scale :: float,shape :: float) :: cdf
  def nakagamiCDF(scale,shape) when scale>0 and shape>0 do
    fn
      x ->
        Exboost.Math.tgamma_lower(shape,shape*(x/scale)*(x/scale))
    end
  end
  def nakagamiCDF(_scale,_shape), do: raise ArithmeticError, "Nakagami is only defined for positive scale and shape"

  ###
  ### Special distributions
  ###

  @doc """
  Distribution for flipping coins.
  """
  @spec coin() :: distribution
  def coin() do
    fn -> if(:rand.uniform()<0.5, do: :heads, else: :tails) end
  end

  @doc """
  Distribution simulating a dice (1..6)
  """
  @spec dice([] | number) :: distribution
  def dice([]), do: dice(1.0)
  def dice([avg: avg]), do: dice(avg)
  def dice(avg), do: uniform([avg*1,avg*2,avg*3,avg*4,avg*5,avg*6])

  @doc """
  Distribution simulating the dice in the GetKanban V4 simulation game.
  """
  @spec dice_gk4([] | number) :: distribution
  def dice_gk4([]), do: dice_gk4(1.0)
  def dice_gk4([avg: avg]), do: dice_gk4(avg)
  def dice_gk4(avg), do: uniform([avg*3,avg*4,avg*4,avg*5,avg*5,avg*6])

  @doc """
  Returns the model for a name.
  
  The kurtosis is the so-called 'excess kurtosis'.
  
  Supported disributions:
      "wald" - The Wald or Inverse Gauss distribution,
      "weibull" - The Weibull distribution,
      "exponential" - The exponential distribution,
      "poisson" - The Poisson distribution,
      "normal" - The normal or Gaussian distribution,
      "fechet" - The Fréchet distribution,
      "nakagami" - The Nakagami distribution,
      "sep" - The Skewed Exponential Power distribution (Azzalini),
      "erlang" - The Erlang distribution,
      "sep0" - The Skewed Exponential Power distribution (Azzalini) with location parameter set to zero (0).
      
  ## Options
  Available only for the SEP distribution, see 'sepCDF/5'.
  """
  @spec model(name::String.t, options::Keyword.t) :: [fun: cdf, df: pos_integer()]
  def model(name, options \\ []) do
    import Exboost.Math, only: [tgamma: 1]

    case name do
      "wald" -> [
        fun: fn (x,[k,lambda]) -> waldCDF(k,lambda).(x) end,
        df: 2,
        skewness: fn [k,lambda] -> 3*:math.sqrt(k/lambda) end,
        kurtosis: fn [k,lambda] -> 15*k/lambda end
      ]
      "weibull" -> [
        fun: fn (x,[k,lambda]) -> weibullCDF(k,lambda).(x) end,
        df: 2,
        skewness: fn
          [k,lambda] ->
            mu = lambda*tgamma(1+1/k)
            sigma = lambda*:math.sqrt(tgamma(1+2/k) - tgamma(1+1/k)*tgamma(1+1/k))
            tgamma(1+3/k)*:math.pow(lambda/sigma,3) - 3*mu/sigma - :math.pow(mu/sigma,3)
           end,
        kurtosis: fn
          [k,lambda] ->
            mu = lambda*tgamma(1+1/k)
            sigma = lambda*:math.sqrt(tgamma(1+2/k) - tgamma(1+1/k)*tgamma(1+1/k))
            skew = tgamma(1+3/k)*:math.pow(lambda/sigma,3) - 3*mu/sigma - :math.pow(mu/sigma,3)
            tgamma(1+4/k)*:math.pow(lambda/sigma,4) - 4*mu/sigma*skew - 6*:math.pow(mu/sigma,2) - :math.pow(mu/sigma,4) - 3.0
          end
      ]
      "exponential" -> [
        fun: fn (x,[k]) -> exponentialCDF(k).(x) end,
        df: 1,
        skewness: fn _ -> 2 end,
        kurtosis: fn _ -> 6 end
      ]
      "frechet" -> [
        fun: fn (x,[scale,shape]) -> frechetCDF(scale,shape).(x) end,
        df: 2,
        skewness: fn
          [_scale,shape] ->
            g1 = tgamma(1.0-1.0/shape)
            g2 = tgamma(1.0-2.0/shape)
            g3 = tgamma(1.0-3.0/shape)
            (g3 - 3*g2*g1 + 2*g1*g1*g1)/:math.pow(g2 - g1*g1,1.5)
          end,
        kurtosis: fn
          [_scale,shape] ->
            g1 = tgamma(1.0-1.0/shape)
            g2 = tgamma(1.0-2.0/shape)
            g3 = tgamma(1.0-3.0/shape)
            g4 = tgamma(1.0-4.0/shape)
            -6 + (g4 - 4*g3*g1 + 3*g2*g2)/:math.pow(g2 - g1*g1,2.0)
          end
      ]
      "nakagami" -> [
        fun: fn (x,[scale,shape]) -> nakagamiCDF(scale,shape).(x) end,
        df: 2,
        skewness: fn
          [_scale,shape] ->
            g = tgamma(shape)
            g1_2 = tgamma(shape+0.5)
            g1 = tgamma(shape+1.0)
            g3_2 = tgamma(shape+1.5)
            num = 2*g1_2*g1_2*g1_2 + g*g*( g3_2 - 3*shape*g1_2 )
            den = g*g*g*:math.pow(shape*(1.0-shape*g1_2*g1_2/g1/g1),1.5)
            num/den
          end,
        kurtosis: fn
          [_scale,shape] ->
            g = tgamma(shape)
            g1_2 = tgamma(shape+0.5)
            g2 = tgamma(shape+2.0)
            gdouble = tgamma(2*shape)
            num = -6*g1_2*g1_2*g1_2*g1_2 - 3*shape*shape*g*g*g*g + g*g*g*g2 + :math.pow(2,3-4*shape)*(4*shape-1)*:math.pi*gdouble*gdouble
            den = :math.pow(abs(g1_2*g1_2 - shape*g*g),2)
            num/den
          end
      ]
      "poisson" -> [
        fun: fn (x,[lambda]) -> poissonCDF(lambda).(x) end,
        df: 1,
        skewness: fn [lambda] -> 1/:math.sqrt(lambda) end,
        kurtosis: fn [lambda] -> 1/lambda end
      ]
      {"poisson", period} when is_number(period) and period>0 -> [
        fun: fn (x,[lambda]) -> poissonCDF(lambda*period).(x) end,
        df: 1,
        skewness: fn [lambda] -> 1/:math.sqrt(lambda*period) end,
        kurtosis: fn [lambda] -> 1/lambda/period end
      ]
      "erlang" -> [
        fun: fn (x,[k,lambda]) -> erlangCDF(k,lambda).(x) end,
        df: 2,
        skewness: fn [k,_] -> 2/:math.sqrt(k) end,
        kurtosis: fn [k,_] -> 6/k end
      ]
      "normal" -> [
        fun: fn (x,[mu,sigma]) -> normalCDF(mu,sigma).(x) end,
        df: 2,
        skewness: fn _ -> 0 end,
        kurtosis: fn [_mu,_sigma] -> 0 end
      ]
      "sep" -> [
        fun: fn (x,[a,b,lambda,alpha]) -> sepCDF(a,b,lambda,alpha,options).(x) end,
        df: 4,
        skewness: fn
          [_a,_b,lambda,_alpha] ->
            delta = lambda/:math.sqrt(1+lambda*lambda)
            pi = :math.pi()
            0.5*(4-pi)*:math.pow(delta*:math.sqrt(2/pi),3)/:math.pow(1-2*delta*delta/pi,1.5)
          end,
        kurtosis: fn
           [_a,_b,lambda,_alpha] ->
             delta = lambda/:math.sqrt(1+lambda*lambda)
             pi = :math.pi()
             2*(pi-3)*:math.pow(delta*:math.sqrt(2/pi),4)/:math.pow(1-2*delta*delta/pi,2)
          end
      ]
      "sep0" -> [
        fun: fn (x,[b,lambda,alpha]) -> sepCDF(0.0,b,lambda,alpha,options).(x) end,
        df: 3
      ]
      unknown ->
        raise UnsupportedDistributionError, message: "Unsupported cumulative distribution function '#{inspect unknown}'"
    end
  end

  @doc """
  Guesses what distribution is likely to fit the sample data
  """
  @spec guess(sample::[number], n::integer, list::[String.t] | String.t) :: [any]
  def guess(sample,n \\ 100,list \\ ["exponential","poisson","normal","erlang","wald","sep","weibull","frechet","nakagami"])
  def guess(sample,n,list) when is_integer(n) and n>0 and is_list(list) do
    {{skewness,err_s},{kurtosis,err_k}} = sample |> cullen_frey(n) |> cullen_frey_point
    list
    |> Enum.flat_map(
      fn
        distrib ->
          r = sample
          |> guess(n,distrib)
          |> Enum.map(fn {s,k}->((skewness-s)/err_s)*((skewness-s)/err_s) + ((kurtosis-k)/err_k)*((kurtosis-k)/err_k) end)
          |> Enum.min
          [{distrib,r}]
      end)
    |> Enum.sort(fn {_,r1},{_,r2} -> r1<r2 end)
  end
  def guess(_sample,n,distrib) when is_integer(n) and n>0 do
    model = model(distrib)
    params = 1..model[:df]
    1..n
    |> Enum.map(fn _ -> Enum.map(params, fn _ -> 50*:rand.uniform end) end)
    |> Enum.flat_map(fn
      pars ->
        try do
          s = model[:skewness].(pars)
          k = model[:kurtosis].(pars)
          [{s,k}]
        rescue
          _error -> []
        end
    end)
  end

  ##
  ## Local Functions
  ##
  @spec random(min::number(),max::number()) :: number()
  defp random(min,max) when max >= min do
    min + (max-min)*:rand.uniform()
  end
  
  @spec phi(x :: float) :: float
  defp phi(x), do: normalCDF(0.0,1.0).(x)

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

  @spec sepPDF(a::float,b::float,lambda::float,alpha::float) :: cdf
  defp sepPDF(a,b,lambda,alpha) do
    fn x ->
      z = (x-a)/b
      t = :math.pow(abs(z),alpha/2.0)
      w = lambda*:math.sqrt(2.0/alpha)*t

      if z > 0.0 do
        :math.exp(-t*t/alpha) * 0.5 * ( 1.0 + :math.erf(w/:math.sqrt(2.0)) )
      else
        :math.exp(-t*t/alpha) * 0.5 * ( :math.erfc(w/:math.sqrt(2.0)) )
      end
    end
  end

end
