FROM ubuntu:24.04

ARG HOST_UID=1000
ARG HOST_GID=1000

RUN apt-get update && apt-get install -y \
  gcc-14 \
  g++-14 \
  cmake \
  gdb \
  libboost-dev \
  libtbb-dev \
  git \
  cmake-format \
  clang-format \
  catch2 && \
  rm -rf /var/lib/apt/lists/*


RUN usermod -l user ubuntu && \
    groupmod -n user ubuntu && \
    usermod -d /home/user -m user



RUN if ["$HOST_UID" != "1000"]; then usermod -u $HOST_UID user; fi && \
    if ["$HOST_GID" != "1000"]; then \
        if ! getent group $HOST_GID > /dev/null; then \
            groupmod -g $HOST_GID user; \
        else \
            usermod -g $HOST_GID user; \
        fi \
    fi

WORKDIR /app

RUN chown user:user /app

ENV CC=gcc-14
ENV CXX=g++-14

USER user

CMD /app/scripts/init/init_container.sh
