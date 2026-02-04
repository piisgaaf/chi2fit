defmodule Chi2fit.MixProject do

  # Copyright 2015-2019 Pieter Rijken
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
      version: "2.0.2",
      elixir: "~> 1.19",
      start_permanent: Mix.env == :prod,
      build_embedded: Mix.env == :prod,
      deps: deps(),
      escript: escript(),
      aliases: aliases(),
      compilers: Mix.compilers ++ if(Mix.env == :docs, do: [:md], else: []), # Add the make compiler
      test_coverage: [tool: ExCoveralls],

      ## Hex stuff:
      description: description(),
      package: package(),
      name: "Chi-SquaredFit",
      source_url: "https://github.com/piisgaaf/chi2fit",

      ## Docs
      docs: docs()
    ]
  end

  # Run "mix help compile.app" to learn about applications.
  def application() do
    [
      extra_applications: [ :logger]
    ]
  end

  # Run "mix help deps" to learn about dependencies.
  defp deps() do
    [
      {:exboost, "~> 0.2"},
      {:graphvix,"~> 1.1"},
      {:csv, "~> 3.2"},
      {:timex, "~> 3.7.13", runtime: false},
      {:stream_data, "~> 1.2"},
      {:gnuplot, "~> 1.22"},
      {:ielixir, "~> 1.0", only: :nb, runtime: false},
      {:tzdata, "~> 1.1", runtime: false},
      {:excoveralls, "~> 0.11.0", only: :test},
      {:ex_doc, "~> 0.40", runtime: false},
      {:credo, "~> 1.7", only: :dev},
      {:mix_test_watch, "~> 1.0", only: :test},
      {:poison, "~> 3.0", only: [:nb,:test]}
    ]
  end

  defp aliases() do
    [
      test_notebook: "test --only notebooks",
      test_perf: "test --only performance",
      test_all: "test --include performance --include notebooks",
      docs: ["docs", &copy_images/1]
    ]
  end

  defp escript() do
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
      licenses: [ "Apache-2.0" ],
      files: [ "lib", "mix.exs", "README*", "LICENSE*", "NOTICE", "config", "chi2fit" ],
      links: %{
        "GitHub" => "https://github.com/piisgaaf/chi2fit",
        "Docker" => "https://hub.docker.com/r/pietertje/chi2fit",
        "Binder" => "https://mybinder.org/v2/gh/piisgaaf/chi2fit/master?filepath=README.ipynb"
      }
    ]
  end

  def cli() do
    [
      preferred_envs: [
        test: :test,
        test_all: :test,
        test_perf: :test,
        test_notebook: :test,
        coveralls: :test,
        docs: :docs,
        "coveralls.html": :test,
        "hex.publish": :docs
      ]
    ]
  end

  defp docs() do
    [
      extras: [
        "README.md",
        "Forecasting-empirical-data.md",
        "Forecasting-fit-to-known-distribution.md",
        "Forecasting-bootstrapping.md",
        "Forecasting-non-equilibrium.md",
        "Example-multi-plot.md",
        "Forecasting-cycle-times.md",
      ] |> Enum.map(& "#{Mix.Project.build_path()}/lib/chi2fit/docs/"<>&1)
    ]
  end

  defp copy_images(_args) do
    "#{Mix.Project.build_path()}/lib/chi2fit/docs/*_files"
    |> Path.wildcard()
    |> Enum.each(& File.cp_r!(&1, "doc/#{Path.basename &1}"))
  end

end

defmodule Mix.Tasks.Compile.Md do
  use Mix.Task.Compiler

  @shortdoc "Compiles with jupyter nbconvert to create markdown from a notebook"

  def run(_) do
    outdir = "#{Mix.Project.build_path()}/lib/chi2fit/docs"
    docs = Chi2fit.Mixfile.project()[:docs][:extras]

    {result, _error_code} = System.cmd("make", docs, stderr_to_stdout: true, env: [{"OUTDIR",outdir}])
    Mix.shell.info result
  end
end
