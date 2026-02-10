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
    {:chi2fit, "~> 2.1"}
  ]
end
```

Documentation can be generated with [ExDoc](https://github.com/elixir-lang/ex_doc)
and published on [HexDocs](https://hexdocs.pm). Once published, the docs can
be found at [https://hexdocs.pm/chi2fit](https://hexdocs.pm/chi2fit).

## Basic usage: command line

The following command does a simple fit against data:

```shell
$ mix escript.build
$ mix escript.install
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

## Documentation

For detailed documentation please visit [https://hexdocs.pm/chi2fit](https://hexdocs.pm/chi2fit).
