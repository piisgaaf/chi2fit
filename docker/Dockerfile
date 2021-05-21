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

FROM pprzetacznik/ielixir
#FROM pietertje/ielixir

# Set the working directory
WORKDIR /app

USER root

# Environment variables
ENV MIX_ENV=prod
ENV BOOSTLIB=/usr/lib/x86_64-linux-gnu/
ENV ERL_LIBS=/app/chi2fit/_build/${MIX_ENV}/lib

# Install the BOOST++ libraries
RUN set -xe \
  && buildDeps='libstdc++-10-dev \
                libboost1.71-dev \
                libboost-math1.71-dev \
                libboost-random1.71-dev \
                libboost-thread1.71-dev \
                libboost-system1.71-dev \
                libboost1.71-tools-dev \
                libmpfrc++-dev \
                libntl-dev' \
  && apt-get update \
  && apt-get install -y --no-install-recommends $buildDeps
  
# Install the XeTeX distribution with additional fonts
RUN set -xe \
  && buildDeps='texlive-xetex \
                texlive-plain-generic \
                texlive-fonts-recommended' \
  && apt-get update \
  && apt-get install -y --no-install-recommends $buildDeps

# Reinstall anaconda
RUN set -xe \
  && rm -rf /opt/conda \
  && mkdir -p /tmp \
  && wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh \
  && bash /tmp/miniconda.sh -b -p /opt/conda

# Install additional conda packages
RUN set -xe \
  && conda update -n base conda \
  && conda install -c conda-forge rise

# Install additional packages needed for the Chi2fit package
RUN set -xe \
  && buildDeps='gnuplot \
                tini \
                graphviz' \
  && apt-get update \
  && apt-get install -y --no-install-recommends $buildDeps

  # Copy the Chi2fit package
WORKDIR /app/chi2fit
ADD mix.* ./
ADD config config/
ADD lib lib/

# Compile dependencies & Chi2Fit
RUN set -xe \
  && mix local.hex --force \
  && mix deps.get \
  && mix deps.compile

RUN set -xe \
  && mix deps.compile \
  && mix compile

# Add notebooks
WORKDIR /app/notebooks
ADD notebooks/*.ipynb ./
ADD notebooks/data/* ./data/

RUN set -xe \
  && mkdir -p /app/notebooks \
  && mkdir -p /app/notebooks/images \
  && chmod -R o+w /app/notebooks \
  && chown -R jovyan:users /home/jovyan
