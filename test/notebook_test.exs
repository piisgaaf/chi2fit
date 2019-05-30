defmodule NotebookTest do

  # Copyright 2019 Pieter Rijken
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

  use ExUnit.Case, async: true
  use NotebookUnit.Case, dir: "notebooks"

  @moduletag timeout: 120000

  nbtest "AgileNXT-Jupiter.ipynb"
  nbtest "Forecasting-empirical-data.ipynb"
  nbtest "Forecasting-fit-to-known-distribution.ipynb"
  nbtest "Forecasting-non-equilibrium.ipynb"
  nbtest "Forecasting-bootstrapping.ipynb"
  nbtest "Forecasting-cycle-times.ipynb"
  nbtest "Forecasting-multiplot.ipynb"

end
