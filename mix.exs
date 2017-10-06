defmodule Chi2fit.Mixfile do
  use Mix.Project

  def project do
    [
      app: :chi2fit,
      version: "0.2.0",
      elixir: "~> 1.5",
      start_permanent: Mix.env == :prod,
      deps: deps(),
      ## Hex stuff:
      description: description(),
      package: package(),
      name: "Chi-SquaredFit"
    ]
  end

  # Run "mix help compile.app" to learn about applications.
  def application do
    [
      extra_applications: [:logger]
    ]
  end

  # Run "mix help deps" to learn about dependencies.
  defp deps do
    [
      {:exalgebra, "~> 0.0.5"},
      {:ex_doc, "~> 0.16", only: :dev, runtime: false}
    ]
  end
  
  # Hex Package Manager stuff:
  defp description() do
    """
    Provides functions for fast matrix inversion, creation of empirical CDF from sample data including
    handling of asymmetric errors, and fitting to a funtion using chi-squared. The fitting procedure return
    the full covariance matrix describing the fitted parameters.
    """
  end
  
  defp package() do
    [
      maintainers: ["Pieter Rijken"],
      licenses: ["Apache 2.0"],
      files: ["lib", "mix.exs", "README*", "LICENSE*", "config"],
      links: %{ }
    ]
  end

end
