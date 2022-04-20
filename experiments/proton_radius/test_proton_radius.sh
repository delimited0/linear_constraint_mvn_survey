#!/usr/bin/env bash

export R_PROGRESSR_ENABLE=TRUE

Rscript experiments/proton_radius/proton_radius.R \
  --method_conf=experiments/proton_radius/method_conf.json \
  --sample_path=experiments/proton_radius/test_samples \
  --reps=4 \
  --seed=2022 \
  --n_cores=2 \
  --n_blas_threads=1