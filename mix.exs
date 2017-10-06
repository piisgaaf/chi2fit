defmodule Chi2fit.Mixfile do
  use Mix.Project

  def project do
    [
      app: :chi2fit,
      version: "0.2.0",
      elixir: "~> 1.5",
      start_permanent: Mix.env == :prod,
      deps: deps()
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
end
