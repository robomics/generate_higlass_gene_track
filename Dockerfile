# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

FROM ubuntu:22.04 AS downloader

RUN apt-get update \
&&  apt-get install -y git \
&&  git clone 'https://github.com/higlass/clodius.git' --depth 1 --branch v0.19.0 /tmp/clodius \
&&  printf '%s\n' '#!/usr/bin/env python3' > /tmp/exonU.py \
&&  cat /tmp/clodius/scripts/exonU.py >> /tmp/exonU.py \
&&  rm -r /tmp/clodius


FROM ubuntu:22.04 AS final

ARG CONTAINER_VERSION
ARG CONTAINER_TITLE

ARG CLODIUS_VER='0.19.*'
ARG PIP_NO_CACHE_DIR=0

RUN apt-get update \
&&  apt-get install -y gawk \
                       gzip \
                       python3 \
                       python3-pip \
                       zstd \
&&  pip install "clodius==$CLODIUS_VER" \
&&  apt-get remove -y python3-pip \
&&  apt-get autoremove -y \
&&  rm -rf /var/lib/apt/lists/*

COPY --from=downloader /tmp/exonU.py /usr/local/bin/

RUN chmod 755 /usr/local/bin/exonU.py

RUN clodius --help
RUN exonU.py --help

ENTRYPOINT ["/bin/bash"]
WORKDIR /data

LABEL org.opencontainers.image.authors='Roberto Rossini <roberros@uio.no>'
LABEL org.opencontainers.image.url='https://github.com/robomics/generate_higlass_gene_track'
LABEL org.opencontainers.image.documentation='https://github.com/robomics/generate_higlass_gene_track'
LABEL org.opencontainers.image.source='https://github.com/robomics/generate_higlass_gene_track'
LABEL org.opencontainers.image.licenses='MIT'
LABEL org.opencontainers.image.title="${CONTAINER_TITLE:-generate-higlass-gene-track}"
LABEL org.opencontainers.image.version="${CONTAINER_VERSION:-latest}"
