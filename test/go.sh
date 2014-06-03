#!/bin/sh

FFVC=../bin/ffvc_mini

export OMP_NUM_THREADS=2

(mpiexec -n 8 $FFVC --scale=strong --size=128 --division=2x2x2 \
  && ./check.py history_base_0.txt history_base.txt) 2>&1 | tee log
