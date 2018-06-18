defmodule Chi2fit.Cli do

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

  @moduledoc """
  Provides a command line interface for fitting data against a known cumulative distribution function.
  
  Tool for fitting particular probability distributions to empirical cumulative distribution functions.
  Distributions supported are Weibull, Wald (Inverse Gauss), Normal, Exponential, Erlang, and Skewed Exponential.

  It uses the Chi-squared Pearson statistic as the likelihood function for fitting. This statistic applies to
  empirical data that is categorial in nature.

  It provides various options for controlling the fitting procedure and assignment of errors. It supports asymmetrical
  errors in fitting the data.

  ## Basic usage: scanning the surface
  
  As described above fitting the parameters is done by minimizing the chi-squared statistic. Usually this is a function of the
  distribution paremeters.
  
  Scanning the surface is a simple way to have an initial guess of the parameters. The following command does a simple scan of
  the chi-squared surface against data:

  ```shell
  $ chi2fit data.csv --ranges '[{0.8,1.2},{0.6,1.2}]' --cdf weibull

  Initial guess:
      chi2:		1399.3190035059733
      pars:		[0.800467783803376, 29.98940654419653]
      errors:		{[0.800467783803376, 0.800467783803376], [29.98940654419653, 29.98940654419653]}
  ```

  and the file `data.csv` is formatted as

  ```
  Lead Time
  26
  0
  105
  69
  3
  36
  ...
  ```

  In this form the command will scan or probe the Chi-squared surface for the parameters within the provided range. It returns the found
  minimum Chi-squared and the parameter values at this minimum. The reported error ranges correspond to a change of Chi-squared of +1.

  Options available:
  
    * `probes` - The number of probes to use for guessing parameter values at initialization
    * `progress` - Shows progress during 'probing' (shows progress every 1000 probes)
      * `c` - Mark progress every 100th probe
      * `x` - Mark progress every 10th probe

  More options are described below and are available using the option `--help`.

  ## Input data options
  
  Several options control how the input data is interpreted. These are:
  
    * `model` - determines how errors are assigned to the data points. Possible values include `simple|asimple|linear`
    * `data` - instead of using the file for data, use this option to pass a list of data points
    * `correction` - Estimate of number of events missed in the right tail of the sample

  An example of specifying data on the command line is:

  ```shell
  $ chi2fit --ranges '[{0.8,1.2},{0.6,1.2}]' --cdf weibull --data '[2,3,4,5,5,4,4,7]'
  ```

  ## Distribution options

  Distributions supported are: Wald, Weibull, Normal, Erlang, Exponential, and SEP (Skewed Exponential: 3 and 4 parameters).

  For the distributions of SEP (4 parameters), and SEP0 (3 parameters) the following options exist:
  
    * `method` -  Supported values are 'gauss|gauss2|gaus3|romberg|romberg2|romberg3'
  
  Romberg integration supports the options:
  
    * `tolerance` - The target precision for Romberg integration
    * `itermax` - The maximum number of iterations to use in Romberg integration
    
  Gauss integration supports the option:

    * `npoints` - The number of points to use in Gauss integration (4, 8, 16, and 32)
  
  ## Fitting options

  AFter probing the surface for an initial guess of the parameters, a fine grained search for the optimum can be done by enabling
  the fit procedure. The algorithm implemented assumes that the initial guess is close enough to the minimum and uses a combination of
  parameter estimation and Monte Carlo methods.
  
  An additional strategy is to use a so-called grid-search by changing only one parameter at a time. It selects the parameters in a
  round robin fashion. Using Romberg iteration and Newton root finding algorithm the parameter value minimizing chi-squared is determined
  while kepping the other parameters constant. Then the other parameters are varied. Especially fitting distributions with 3 or more
  parameters may benefit from this strategy.
  
  Options controlling these are:

    * `fit` - Enables the fine-grained fitting of parameters
    * `iterations` - Number of iterations to use in the optimizing the Likelihood function
    * `grid` - Uses a grid search to fit one parameter at a time in a round robin fashion

  Sometimes the chi-squared surface is not smooth but numerically problematic to get stable. In this case smoothing the surface
  may help. The next option enables this feature:

    * `smoothing` - Smoothing of the likelihood function with a Gaussian kernel

  The fitting procedures uses derivatives (first and second order) to estimate changes in the parameters that will result in
  a better fit. Derivaties are calculated using Romberg differentiation. The accuracy and maximum number of iterations are
  controlled by the options:
  
    * `tolerance` - The target precision for Romberg integration
    * `itermax` - The maximum number of iterations to use in Romberg integration
  
  ## Bootstrapping

  Bootstrapping can be enabled to estimate the errors in the parameters. The supported options are:
  
    * `bootstrap` - Enables bootstrapping. Specifies the number of iterations to perform
    * `sample` - The sample size to use from the empirical distribution

  ## Output options

  These options are useful for printing data for generating charts of the data:

    * `print` - Outputs the empirical input data with errors included
    * `output` - Outputs the fitted distribution function values at the data points
    * `surface` -  Outputs the Chi-squared surface to a file
    * `smoothing` - Smoothing of the likelihood function with a Gaussian kernel

  ## General options

  Options available for scanning, fitting, and bootstrapping:

    * `debug` - Outputs additional data for debugging purposes"
  
  ## References
  
    [1] R.A. Arndt and M.H. MacGregor, Methods in Computational Physics, Vol. 6 (1966) 256-296

    [2] Marius M. Nagels, Baryon-Baryon Scattering in a One-Boson-Exchange Potential Mode, PhD. Thesis, Nijmegen University, 1975

    [3] Richard A. Arndt and Malcolm H. MacGregor, Determination of the Nucleon-Nucleon Elastic-Scattering Matrix. IV. Comparison of Energy-Dependent and Energy-Independent Phase-Shift Analyses, Physical Review Volume 142, Number 3, January 1966
  """

  require Logger

  import Chi2fit.Fit, only: [chi2fit: 4, chi2probe: 4, chi2: 4]
  import Chi2fit.Utilities
  import Chi2fit.Matrix
  import Chi2fit.Distribution

  @datapoints 500
  @maxx 1.1
  @default_iterations 10
  @default_parameter_iterations 10
  @default_probes 100_000
  @default_surface_file "cdf_surface.csv"
  @default_cdf "weibull"
  @default_asymm :linear
  @default_int_method :romberg2
  @default_tolerance 1.0e-6
  @default_npoints 32
  @default_binsize 1
  @default_binoffset 0.5
  @default_error_score :wilson
  @default_bin {@default_binsize,0.5}

  @jac_threshold 0.01
  
  defp penalties(_x,_pars), do: 0.0

  defp probe(data, model, options) do
    penalties = options[:penalties]
    surface = options[:surface]
    surface? = options[:surface?]
    
    {:ok, file} = if surface?, do: File.open(surface, [:write]), else: {:ok,nil}
    options = options |> Keyword.put_new(:surfacefile,file)
    result = chi2probe(data, model[:probe], {model[:fun], penalties}, options)
    if file, do: File.close(file)
    result
  end
  
  defp print_cdf({cdf,[_,maxdur]}, options) do
    0..options[:datapoints]
    |> Stream.map(&(maxdur*options[:maxx]*&1/options[:datapoints]))
    |> Stream.map(fn x -> {x,cdf.(x)} end)
    |> Enum.each(fn ({x,{y,ylow,yhigh}})-> IO.puts("#{x},#{y},#{ylow},#{yhigh}") end)
    System.halt(0)
  end
  
  defp prepare_data(data, options) do
    mcsample = options[:mcsample]
    correction = options[:correction]

    workdata = cond do
      mcsample == :all -> data |> Enum.to_list
      true -> data |> Enum.take_random(mcsample)
    end
    if mcsample != :all && (workdata |> Enum.sum)*(data|>Enum.count)<(data |> Enum.sum)*mcsample/2 do
      IO.puts "WARNING: maximum of sample is smaller than average of complete sample"
      IO.puts "    Sample = #{inspect(workdata)}"
    end

    {cdf,bins,_,_} = get_cdf(workdata,@default_bin,@default_error_score,correction)
    {mindur,_,_} = bins |> hd
    {maxdur,_,_} = bins |> List.last
    if options[:print?], do: print_cdf({cdf,[mindur,maxdur]}, options)

    data = convert_cdf({cdf,bins|>Enum.map(&elem(&1,0))})

    try do
      model = model(options[:name],options) |> Keyword.put(:probe, elem(Code.eval_string(options[:ranges]),0))
      {chi2, parameters,errors} = probe data, model, options
      {data,model, {chi2, parameters,errors}}
    rescue
      _e in Chi2fit.Distribution.UnsupportedDistributionError ->
        IO.puts :stderr, "ERROR: Unsupported distribution '#{options[:name]}'"
        System.halt 1
    end
  end
  
  defp do_output(data, parameters, model, alphainv, options) do
    data |> Enum.sort |> Enum.each(fn
      (x)->
        jac = jacobian parameters, fn (pars)->model[:fun].(x,pars) end, options
        error2 = alphainv |> Enum.map(&(ExAlgebra.Vector.dot(&1,jac))) |> ExAlgebra.Vector.dot(jac)
        try do
          y = model[:fun].(x,parameters)
          error = if abs(error2/y) < 1.0e-6, do: 1.0e-6, else: :math.sqrt(error2)
          IO.puts("#{x},#{y},#{y-error},#{y+error}")
        rescue
          ArithmeticError -> IO.puts "Warning: arithmetic error (probably negative diagonal element (#{error2}) in covariance matrix)"
        end
      end)
  end
  
  defp usage(code) do
    IO.puts "Usage: #{__ENV__.file |> String.split("/") |> Enum.reverse |> hd} <options> <data file>"
    IO.puts "    --help\t\t\t\tShows this help"
    IO.puts ""
    IO.puts "    --ranges \"[{...,...},...]\"\t\tRanges of parameters to search for minimum likelihood"
    IO.puts "    --cdf <cdf>\t\t\t\tThe distribution function (defaults to '#{@default_cdf}') to fit the data."
    IO.puts "    \t\t\t\t\tSupported values are 'wald|weibull|exponential|sep|sep0'"
    IO.puts ""
    IO.puts "    Input data:"
    IO.puts "    --model simple|asimple|linear\tThe model (defaults to '#{@default_asymm}') to use for handling asymmetrical errors in the input data"
    IO.puts "    --data <data>\t\t\tArray of data points to use in fitting"
    IO.puts "    --correction <integer>\t\tEstimate of number of events missed in the right tail of the sample"
    IO.puts ""
    IO.puts "    Fitting data to a CDF:"
    IO.puts "    --guess <number>\t\t\tGuess what distribution fits best. Use <number> of bootstraps."
    IO.puts "    --fit\t\t\t\tTry to fit the parameters"
    IO.puts "    --iterations <number>\t\tNumber of iterations (defaults to '#{@default_iterations}') to use in the optimizing the Likelihood function"
    IO.puts "    --probes <number>\t\t\tThe number of probes (defaults to '#{@default_probes}') to use for guessing parameter values at initialization"
    IO.puts "    --pariter <number>\t\t\tThe maximum number of parameter variations to try (Monte Carlo) and maximum number of iterations in Newton root finding algorithm"
    IO.puts "    --grid\t\t\t\tUses a grid search to fit one parameter at a time in a round robin fashion"
    IO.puts ""
    IO.puts "    Distributions:"
    IO.puts "    --method <method>\t\t\tThe integration method to use (defaults to '#{@default_int_method}'); applies to sep and sep0 only."
    IO.puts "    \t\t\t\t\tSupported values are 'gauss|gauss2|gaus3|romberg|romberg2|romberg3'"
    IO.puts "    --tolerance <tolerance>\t\tThe target precision (defaults to '#{@default_tolerance}') for Romberg integration"
    IO.puts "    --itermax <integer>\t\t\tThe maximum number of iterations to use in Romberg integration"
    IO.puts "    --npoints <points>\t\t\tThe number of points to use in Gauss integration (defaults to '#{@default_npoints}')"
    IO.puts ""
    IO.puts "    Bootstrapping:"
    IO.puts "    --bootstrap <integer>\t\tEnables bootstrapping. Specifies the number of iterations to perform."
    IO.puts "    --sample <size>\t\t\tThe sample size to use from the empirical distribution"
    IO.puts ""
    IO.puts "    Output:"
    IO.puts "    --print\t\t\t\tOutputs the input data"
    IO.puts "    --output\t\t\t\tOutputs the fitted distribution function values at the data points"
    IO.puts "    --surface <file>\t\t\tOutputs the Chi-squared surface to a file (defaults to '#{@default_surface_file}')"
    IO.puts "    --smoothing\t\t\t\tSmoothing of the likelihood function"
    IO.puts ""
    IO.puts "    General options:"
    IO.puts "    --progress\t\t\t\tShows progress during 'probing'"
    IO.puts "    --c\t\t\t\t\tMark progress every 100th probe"
    IO.puts "    --x\t\t\t\t\tMark progress every 10th probe"
    IO.puts "    --debug\t\t\t\tOutputs additional data for debugging purposes"
    System.halt(code)
  end

  defp parse_args args do
    case OptionParser.parse args, strict: [
      help: :boolean,
      debug: :boolean,
      print: :boolean,
      guess: :integer,
      cdf: :string,
      data: :string,
      bootstrap: :integer,
      correction: :integer,
      output: :boolean,
      surface: :string,
      iterations: :integer,
      pariter: :integer,
      model: :string,
      tolerance: :float,
      itermax: :integer,
      npoints: :integer,
      method: :string,
      probes: :integer,
      ranges: :string,
      smoothing: :boolean,
      sample: :integer,
      fit: :boolean,
      grid: :boolean,
      progress: :boolean,
      c: :boolean,
      x: :boolean] do
        {options, [filename], []} -> {options,filename}
        {options, [], []} -> {options,nil}
        _else -> usage(1)
    end
  end

  defp add_defaults(options) do
    options = options
    |> Keyword.put_new(:debug?,     options[:debug] || false)
    |> Keyword.put_new(:print?,     options[:print] || false)
    |> Keyword.put_new(:output?,    options[:output] || false)
    |> Keyword.put_new(:surface?,   options[:surface] || false)
    |> Keyword.put_new(:surface,    @default_surface_file)
    |> Keyword.put_new(:name,       options[:cdf] || @default_cdf)
    |> Keyword.update(:model,       @default_asymm, &String.to_atom/1)
    |> Keyword.update(:method,      @default_int_method, &String.to_atom/1)
    |> Keyword.put_new(:iterations, @default_iterations)
    |> Keyword.put_new(:pariter,    @default_parameter_iterations)
    |> Keyword.put_new(:probes,     @default_probes)
    |> Keyword.put_new(:ranges,     nil)
    |> Keyword.put_new(:correction, options[:correction] || 0)
    |> Keyword.put_new(:smoothing,  false)
    |> Keyword.put_new(:fit?,       options[:fit] || false)
    |> Keyword.put_new(:grid?,      options[:grid] || false)
    |> Keyword.put_new(:progress?,  options[:progress] || false)
    |> Keyword.put_new(:mcsample,   options[:sample] || :all)
    |> Keyword.put_new(:mcbootstrap,options[:bootstrap] || 1)
    |> Keyword.put_new(:mcdata,     options[:data] || false)
    
    options
    |> Keyword.put_new(:mark,       [
      m: fn -> if(options[:progress?], do: IO.write("M")) end,
      c: fn -> if(options[:x] || options[:c], do: IO.write("C")) end,
      x: fn -> if(options[:x], do: IO.write("X")) end,
      *: fn -> if(options[:progress?], do: IO.write("*")) end
      ])
    #
    |> Keyword.put_new(:datapoints, @datapoints)
    |> Keyword.put_new(:maxx,       @maxx)
    |> Keyword.put_new(:penalties,  &penalties/2)
  end

  defp kernel(options) do
    fn sample, wwww ->
      IO.write "#{wwww}/#{options[:mcbootstrap]} Running chi-squared fit: progress:\t"
      {data,model, {_chi2, parameters,_errors}} = prepare_data sample, options
      try do
        IO.write "...fitting..."
        fit = {_,_,pars,_ranges} = chi2fit(data, {parameters, model[:fun], &penalties/2}, options[:iterations], options)
        jac = jacobian(pars,&chi2(data,fn (x)->model[:fun].(x,&1) end,fn (x)->penalties(x,&1) end,options),options)
        |> Enum.map(&(&1*&1))|>Enum.sum|>:math.sqrt
        if jac<@jac_threshold, do: fit, else: {:error, "not in minimum #{jac}"}
      catch
        {:inverse_error, ArithmeticError, chi2, _parameters} ->
          IO.puts "(chi2=#{chi2}; dof=#{length(sample)-model[:df]})"
          {chi2,[],parameters}
      else
        {:error, msg} ->
          IO.puts"..#{msg}...skipping"
          nil
  
        {chi2, alphainv, parameters,_ranges} ->
          IO.puts "(chi2=#{chi2}; dof=#{length(sample)-model[:df]})"
          {chi2, alphainv, parameters}
      end
    end
  end

  defp validate options, filename do
    cond do
      filename == nil && options[:data] == nil ->
        IO.puts :stderr, "ERROR: please specify either a data file or 'data'"
        System.halt 1
      filename != nil && options[:data] != nil ->
        IO.puts :stderr, "ERROR: please specify either a data file or 'data'"
        System.halt 1
      filename && !File.exists?(filename) ->
        IO.puts :stderr, "ERROR: failed to open file '#{filename}'"
        System.halt 1

      options[:guess] != nil && options[:guess] < 100 ->
        IO.puts :stderr, "ERROR: 'guess' needs 100 or more bootstraps"
        System.halt 1

      options[:guess] == nil ->
        model = model(options[:cdf],options)
        ranges = options[:ranges] && elem(Code.eval_string(options[:ranges]),0)

        cond do
          options[:ranges] == nil ->
            IO.puts :stderr, "ERROR: please specify 'ranges' for parameters"
            System.halt 1

          length(ranges) != model[:df] ->
            IO.puts :stderr, "ERROR: 'ranges' must be of length #{model[:df]} for '#{options[:cdf]}'"
            System.halt 1

          true -> options
        end

      true -> options
    end
  end

  ## When called from 'mix run -e ...'
  def main, do: main(System.argv())
  
  ## When called from escript
  def main args do
    {options, filename} = parse_args(args)

    ## Help
    if options[:help], do: usage(0)

    ## Default options
    options = try do
      options |> validate(filename) |> add_defaults
    rescue
      Chi2fit.Distribution.UnsupportedDistributionError ->
        IO.puts :stderr, "ERROR: Unsupported distribution '#{options[:cdf]}'"
        System.halt 1
    end

    ## Read the data
    data = if options[:mcdata], do: elem(Code.eval_string(options[:mcdata]),0), else: read_data(filename)

    ## Guess
    if options[:guess] do
      IO.puts String.pad_trailing(~s(Distribution),20)<>"Score"
      IO.puts String.pad_trailing(~s(------------),20)<>"_____"
      data
      |> guess(options[:guess])
      |> Enum.each(fn {name,score} -> IO.puts String.pad_trailing(name,20)<>"#{score}" end)
      System.halt 0
    end

    cond do
      options[:mcbootstrap]>1 and options[:fit?] ->
        wdata = if options[:mcsample] == :all, do: data, else: data |> Enum.take_random(options[:mcsample])
        boot = bootstrap(options[:mcbootstrap], wdata, kernel(options),options)
        boot = boot |> Enum.filter(&is_tuple/1)

        # Compute average, average sd, sd error, and maximum lag that occured
        model = model(options[:name],options) |> Keyword.put(:probe, elem(Code.eval_string(options[:ranges]),0))
        avgchi2 = (boot |> Stream.map(fn ({chi2,_,_}) -> chi2 end) |> Enum.sum)/length(boot)
        sdchi2 = :math.sqrt((boot |> Stream.map(fn {chi2,_,_}->(chi2-avgchi2)*(chi2-avgchi2) end) |> Enum.sum))/length(boot)

        avgpars = boot |> Stream.map(fn {_,_,pars} -> pars end) |> Stream.map(&List.to_tuple/1) |> Enum.to_list |> :lists.unzip |> Tuple.to_list
        |> Enum.map(&(Enum.sum(&1)/length(boot)))

        sdpars = boot |> Stream.map(fn {_,_,pars} -> pars end) |> Stream.map(&List.to_tuple/1) |> Enum.to_list |> :lists.unzip |> Tuple.to_list
        |> Enum.zip(avgpars) |> Enum.map(fn {parlist,avg} -> :math.sqrt(parlist|>Enum.map(&((&1-avg)*(&1-avg)))|>Enum.sum)/length(parlist) end)

        avgsd = boot |> Stream.map(fn {_,cov,_} -> cov end) |> Stream.filter(&(length(&1)>0)) |> Stream.map(&diagonal/1) |> Stream.map(&(Enum.map(&1,fn x->:math.sqrt(abs(x)) end))) |> Stream.map(&List.to_tuple/1) |> Enum.to_list |> :lists.unzip |> Tuple.to_list
        |> Enum.map(&(Enum.sum(&1)/length(&1)))
      
        IO.puts "Sample:"
        IO.puts "    #{inspect wdata|>Enum.to_list}"
        IO.puts ""
        IO.puts "Final:"
        IO.puts "    chi2:\t\t\t#{avgchi2}"
        IO.puts "    SD (chi2):\t\t\t#{sdchi2}"
        IO.puts "    parameters:\t\t\t#{inspect avgpars}"
        IO.puts "    SD (parameters; sample):\t#{inspect sdpars}"
        IO.puts "    SD (parameters; fit):\t#{inspect avgsd}"
        IO.puts "    Degrees of freedom:\t\t#{length(wdata|>Enum.to_list)-model[:df]}"
        IO.puts "    Total:\t\t\t#{length(boot)}"

        if options[:output?], do: do_output(wdata, avgpars, model, sdpars |> Enum.map(&(&1*&1)) |> from_diagonal, options)

      true ->
        ## TODO: pass the results from probing to the chi2fit function (parameter ranges & chi2 info)
        {data,model, {chi2, parameters,errors}} = prepare_data data, options

        IO.puts "\n\nInitial guess:"
        IO.puts "    chi2:\t\t#{chi2}"
        IO.puts "    pars:\t\t#{inspect parameters}"
        IO.puts "    errors:\t\t#{inspect errors}\n"
  
        if options[:fit?] do
          options = if options[:debug?] do
            options
            |> Keyword.put_new(:onstep, fn %{chi2: newchi, params: vec} ->
              IO.puts "\n[-]\tchi=#{newchi}"
              IO.puts "\tparameters=#{inspect vec}"
            end)
            |> Keyword.put_new(:onbegin, fn %{step: step, chi: chi2, derivatives: ders} ->
              IO.puts "\n[#{step}]\tchi2=#{chi2}"
              IO.puts "\tderivatives(first,second)=#{inspect ders}"
            end)
          else
            options
          end
          {chi2, alphainv, parameters, ranges} = chi2fit(data, {parameters, model[:fun], &penalties/2}, options[:iterations], options)
          IO.puts "Final:"
          IO.puts "    chi2:\t\t#{chi2}"
          IO.puts "    Degrees of freedom:\t#{length(data)-model[:df]}"
          IO.puts "    covariance:\t\t["
          alphainv |> Enum.each(fn row -> IO.puts "    \t\t\t  #{inspect row}" end)
          IO.puts "    \t\t\t]"
          IO.puts "    gradient:\t\t#{inspect jacobian(parameters,&chi2(data,fn (x)->model[:fun].(x,&1) end,fn (x)->penalties(x,&1) end,options),options)}"
          IO.puts "    parameters:\t\t#{inspect parameters}"
          IO.puts "    errors:\t\t#{inspect alphainv |> diagonal |> Enum.map(fn x->x|>abs|>:math.sqrt end)}"
          IO.puts "    ranges:"
          ranges
          |> Tuple.to_list
          |> Enum.with_index
          |> Enum.each(fn
            {[mn,mx],0} -> IO.puts "\t\t\tchi2:\t\t#{mn}\t-\t#{mx}"
            {[mn,mx],_} -> IO.puts "\t\t\tparameter:\t#{mn}\t-\t#{mx}"
          end)

          if options[:output?], do: do_output(Enum.map(data, fn {x,_,_,_}->x end), parameters, model, alphainv, options)
        end
    end
  end
  
end
