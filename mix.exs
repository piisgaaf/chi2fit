defmodule Chi2fit.Mixfile do

  # Copyright 2015-2017 Pieter Rijken
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

  use Mix.Project

  def project do
    [
      app: :chi2fit,
      version: "0.6.6",
      elixir: "~> 1.5",
      start_permanent: Mix.env == :prod,
      deps: deps(),
      escript: escript(),
      aliases: aliases(),
      preferred_cli_env: ["test": :test, "test_all": :test, "test_perf": :test],
      ## Hex stuff:
      description: description(),
      package: package(),
      name: "Chi-SquaredFit",
      source_url: "https://github.com/piisgaaf/chi2fit"
    ]
  end

  # Run "mix help compile.app" to learn about applications.
  def application do
    [
      extra_applications: [ :logger ]
    ]
  end

  # Run "mix help deps" to learn about dependencies.
  defp deps do
    [
      {:exalgebra, "~> 0.0.5"},
      {:ex_doc, "~> 0.16", only: :dev, runtime: false}
    ]
  end
  
  defp aliases do
    [
      test_perf: "test --only performance",
      test_all: "test --include performance"
    ]
  end

  defp escript do
    [
      main_module: Chi2fit.Cli
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
      maintainers: [ "Pieter Rijken" ],
      licenses: [ "Apache 2.0" ],
      files: [ "lib", "mix.exs", "README*", "LICENSE*", "config", "chi2fit" ],
      links: %{ "GitHub" => "https://github.com/piisgaaf/chi2fit" }
    ]
  end

end
