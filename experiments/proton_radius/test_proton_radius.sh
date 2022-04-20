#!/usr/bin/env bash

export R_PROGRESSR_ENABLE=TRUE

Rscript experiments/proton_radius/proton_radius.R \
  --method_conf=experiments/proton_radius/method_conf.json \
  --sample_path=experiments/proton_radius/test_samples \
  --reps=18 \
  --seed=2022 \
  --n_cores=1 \
  --n_blas_threads=4