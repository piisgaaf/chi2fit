FROM pprzetacznik/ielixir

USER root

# Set the working directory
WORKDIR /app

# Environment variables
ENV BOOSTLIB=/usr/lib/x86_64-linux-gnu/

# Install the BOOST++ libraries
RUN set -xe \
  && buildDeps='libstdc++-6-dev \
                libboost1.62-dev \
                libboost-math1.62-dev \
                libboost-random1.62-dev \
                libboost-thread1.62-dev \
                libboost-system1.62-dev \
                libboost1.62-tools-dev \
                libmpfrc++-dev \
                libntl-dev' \
  && apt-get update \
  && apt-get install -y --no-install-recommends $buildDeps

WORKDIR /app/chi2fit
ADD mix.* ./
ADD config config/
ADD lib lib/

# Compile dependencies & Chi2Fit
RUN set -xe \
  && mix local.hex --force \
  && mix deps.get \
  && mix deps.compile \
  && mix compile

WORKDIR /app/notebooks
ADD notebooks/*.ipynb ./

RUN set -xe \
  && mkdir -p /app/notebooks \
  && chmod o+w /app/notebooks
