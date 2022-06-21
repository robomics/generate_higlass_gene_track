# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

#!/usr/bin/env bash

set -e
set -x
set -u

nextflow run -c example.config main.nf -resume
