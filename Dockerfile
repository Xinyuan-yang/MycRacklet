FROM debian:stable-slim

MAINTAINER Thibault Roch <thibault.roch@epfl.ch>

# Library dependencies
RUN apt-get -qq update && apt-get install -y \
    gcc cmake \
    libgsl-dev libfftw3-dev libfftw3-openmp-dev \
    python3 python3-dev python3-numpy python3-mpi4py \
  && rm -rf /var/lib/apt/lists/*

# for documentation 
RUN apt-get -qq update && apt-get -qq -y install \
    python3-sphinx python3-sphinx-rtd-theme python3-breathe doxygen graphviz \
    && rm -rf /var/lib/apt/lists/*

# for ci on c4science.ch
RUN apt-get -qq update && apt-get -qq -y install \
    python3-pytest git \
    php-cli php-curl php-xml \
    python3-phabricator python3-click python3-yaml \
    && rm -rf /var/lib/apt/lists/*