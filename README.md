# Chi2fit

Tool for fitting particular probability distributions to empirical cumulative distribution functions.
Distributions supported are Weibull, Wald (Inverse Gauss), Normal, Exponential, Poisson, Erlang, and Skewed Exponential.

It uses the Chi-squared Pearson statistic as the likelihood function for fitting. This statistic applies to
empirical data that is categorial in nature.

It provides various options for controlling the fitting procedure and assignment of errors. It supports asymmetrical
errors in fitting the data.

## Installation

If [available in Hex](https://hex.pm/docs/publish), the package can be installed
by adding `chi2fit` to your list of dependencies in `mix.exs`:

```elixir
def deps do
  [
    {:chi2fit, "~> 0.9"}
  ]
end
```

Documentation can be generated with [ExDoc](https://github.com/elixir-lang/ex_doc)
and published on [HexDocs](https://hexdocs.pm). Once published, the docs can
be found at [https://hexdocs.pm/chi2fit](https://hexdocs.pm/chi2fit).

## Docker & Jupyter Notebooks

`Chi2fit` can be used together with Jupyter Notebooks. The easiest way is to create a docker image and run it.
The docker image is based on [IElixir](https://github.com/pprzetacznik/IElixir).

The image is built using:

```shell
$ docker build -f docker/Dockerfile .
```

Run the image with the command:

```shell
$  docker run -p 8888:8888 --hostname 127.0.0.1 -v /tmp:/app/notebooks <docker image id>
```

In Jupyter use one of the provided example notebooks to learn how Chi2fit is set-up from within a notebook.

## Using the prebuilt docker container

Instead of building the docker image yourself, docker images are available at [https://hub.docker.com/r/pietertje/chi2fit](https://hub.docker.com/r/pietertje/chi2fit). After starting the container the log
shows the url to connect to the jupyter notebook.

## Basic usage: command line

The following command does a simple fit against data:

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

More options are available using the option `--help`.

## Basic usage: jupyter notebook

The repository contains the notebooks:

* `chi2fit.ipynb` - simple template containing the minimal set-up to get started,
* `BacklogForecasting.ipynb` - elaborate example using data to forecast the completion date of a backlog of work items,
* `BacklogForecasting-plots.ipynb` - same as `BacklogForecasting.ipynb` but with plots using `GnuPlot`; see below,
* `BacklogForecasting-non-equilibrium.ipynb` - illustration of binning and changing delivery rate

Plots are supported using the package [:gnuplot](https://hex.pm/packages/gnuplot).
On MacOS execute the following command from the shell to display the GnuPlot window:

```shell
$ socat TCP-LISTEN:6000,reuseaddr,fork UNIX-CLIENT:\"$DISPLAY\"
```

On a Mac using `port` the tool `socat` is installed by the command:

```shell
$ sudo port install socat
```  
## Documentation

For detailed documentation please visit [https://hexdocs.pm/chi2fit](https://hexdocs.pm/chi2fit).
