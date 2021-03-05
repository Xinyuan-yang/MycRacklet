FROM debian:stable-slim

MAINTAINER Thibault Roch <thibault.roch@epfl.ch>

# Add contrib and non-free
RUN sed -i 's/main/main contrib non-free/' /etc/apt/sources.list

# Install any needed packages from ubuntu repos
RUN apt-get -qq update && apt-get install -y \
    curl \
    gcc \
    git \
    libboost-dev \
    libpython3-dev \
    libthrust-dev \
    libfftw3-dev \
    libfftw3-mpi-dev \
    python3 \
    python3-dev \
    python3-numpy \
    python3-pip \
    python3-scipy \
    python3-h5py \
    python3-netcdf4 \
    python3-phabricator \
    python3-click \
    python3-yaml \
    python3-matplotlib \
    python3-sphinx \
    python3-breathe \
    python3-sphinx-rtd-theme \
    python3-mpi4py \
    doxygen \
    scons \
  && rm -rf /var/lib/apt/lists/*

RUN pip3 install pytest uvw