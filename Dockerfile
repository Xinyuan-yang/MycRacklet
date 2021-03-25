FROM ubuntu:18.04
MAINTAINER Thibault Roch <thibault.roch@epfl.ch>

ENV LANG C.UTF-8
ENV LC_ALL C.UTF-8

## for apt to be noninteractive
ENV DEBIAN_FRONTEND noninteractive
ENV DEBCONF_NONINTERACTIVE_SEEN true

RUN apt-get -qq update && apt-get -y -qq install \
    clang \
    cmake \
    g++ \
    gfortran \
    curl \
    doxygen \
    git \
    libgsl-dev \
    libfftw3-dev \
    python3-pip \
    python3-breathe \
    python3-dev \
    python3-pytest \
    python3-mpi4py \
    python3-numpy \
    python3-pytest \
    python3-sphinx \
    python3-sphinx-rtd-theme \
    && rm -rf /var/lib/apt/lists/*

# apt-get on one line due to https://docs.docker.com/develop/develop-images/dockerfile_best-practices/#run
