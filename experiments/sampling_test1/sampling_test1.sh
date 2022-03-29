#!/usr/bin/env bash

Rscript experiments/sampling_test1/sampling_test1.R \
  --method_conf=experiments/sampling_test1/method_conf.json \
  --dim_conf=experiments/sampling_test1/dim_conf.json \
  --sample_path=experiments/sampling_test1/samples \
  --timing_path=experiments/sampling_test1/timings \
  --seed=2022 \
  --n_threads=4