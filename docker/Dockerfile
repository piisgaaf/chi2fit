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

#
# To circumvent that the imported base uses arch amd64, but miniforge3
# expects arm64 build using:
# $ docker buildx build --platform=linux/amd64 -f docker/Dockerfile .
#

FROM quay.io/jupyter/base-notebook:aarch64-ubuntu-24.04

# Set the working directory
WORKDIR /app

USER root

# Environment variables
ENV MIX_ENV=prod
ENV ERL_LIBS=/app/chi2fit/_build/${MIX_ENV}/lib

# Install the XeTeX distribution with additional fonts
RUN set -xe \
  && buildDeps='texlive-xetex \
                texlive-plain-generic \
                texlive-fonts-recommended' \
  && apt-get update \
  && apt-get install -y --no-install-recommends $buildDeps

# Install erlang (OTP27) and elixir (1.18):
RUN set -xe \
  && buildDeps='erlang \
                elixir' \
  && apt-get update \
  && apt-get install -y --no-install-recommends $buildDeps

# IElixir:
RUN set -xe \
  && buildDeps='libzmq3-dev \
                libsqlite3-dev \
                unzip \
                make \
                g++ \
                gcc' \
  && apt-get update \
  && apt-get install -y --no-install-recommends $buildDeps

# Source: https://github.com/pprzetacznik/IElixir
ENV CFLAGS="-x c"
RUN set -xe \
  && wget -q https://github.com/pprzetacznik/IElixir/archive/refs/heads/master.zip -O ielixir.zip \
  && unzip ielixir.zip -d ielixir_tmp \
  && mv ielixir_tmp/IElixir-master ielixir \
  && rm -rf ielixir_tmp \
  && cd ielixir \
  && ls -alh \
  && mix local.hex --force \
  && mix local.rebar --force \
  && mix deps.get \
  && mix deps.compile \
  && ./install_script.sh
ADD config/config_ielixir.exs ielixir/config/config.exs
ADD config/prod_ielixir.exs ielixir/config/prod.exs

# Install the BOOST++ libraries
RUN set -xe \
  && buildDeps='libstdc++-11-dev \
                libboost-dev \
                libboost-math-dev \
                libboost-random-dev \
                libboost-thread-dev \
                libboost-system-dev \
                libboost-tools-dev \
                libmpfrc++-dev \
                libntl-dev' \
  && apt-get update \
  && apt-get install -y --no-install-recommends $buildDeps

# Install additional packages needed for the Chi2fit package
RUN set -xe \
  && buildDeps='gnuplot \
                tini \
                graphviz' \
  && apt-get update \
  && apt-get install -y --no-install-recommends $buildDeps

ENV CC=c++
ENV BOOSTLIB=/usr/lib/aarch64-linux-gnu
ENV BOOSTINCLUDE=/usr/include/boost

# Copy the Chi2fit package
WORKDIR /app/chi2fit
ADD mix.* ./
ADD config config/
ADD lib lib/

# Compile dependencies & Chi2Fit
RUN set -xe \
  && mix deps.clean exboost --build \
  && mix clean \
  && mix local.hex --force \
  && mix deps.get \
  && mix deps.compile \
  && mix compile

RUN set -xe \
  && mkdir -p /app/notebooks \
  && mkdir -p /app/notebooks/images \
  && chmod -R o+w /app/notebooks \
  && chown -R jovyan:users /home/jovyan

# Add notebooks
WORKDIR /app/notebooks
ADD notebooks/*.ipynb ./
ADD notebooks/data/* ./data/

RUN set -xe \
  && chown -R jovyan:users /app/*

USER $NB_UID

ENV ELIXIR_ERL_OPTIONS="-smp enable"
CMD ["start-notebook.py"]
