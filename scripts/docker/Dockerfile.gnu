FROM ubuntu:18.04

ENV DEBIAN_FRONTEND noninteractive
ENV TZ=Europe/Berlin

# Package dependencies

RUN apt-get -y update && apt-get -y upgrade && \
    apt-get -y install git cmake build-essential autotools-dev autoconf g++ gfortran pkg-config scons \
        mpich libmpich-dev libpthread-stubs0-dev zlib1g-dev libnuma-dev libyaml-cpp-dev libcurl3-dev \
        gdb valgrind rsync && \
    rm -rf /var/lib/apt/lists/*

# Libraries

COPY scripts/XDMF /app/samoa/scripts/XDMF
WORKDIR /app/samoa/scripts/XDMF
RUN ./install_all.sh mpi gnu /app/samoa_xdmf_libs && \
    ./install_all.sh nompi gnu /app/samoa_xdmf_libs && \
    rm -rf /app/samoa/scripts/XDMF/libsrc

# Entrypoint

WORKDIR /app/samoa
CMD scons -j 8 config=my_conf_gnu.py