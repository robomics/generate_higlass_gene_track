#!/usr/bin/env bash

# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

set -e
set -x
set -u

nextflow run -c config/example-hg38.config main.nf -resume
