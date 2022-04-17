#!/usr/bin/env bash

export R_PROGRESSR_ENABLE=TRUE

Rscript experiments/sample_identity_box/identity_box.R \
  --method_conf=experiments/sample_identity_box/test_method_conf.json \
  --dim_conf=experiments/sample_identity_box/test_dim_conf.json \
  --half_width=3 \
  --sample_path=experiments/sample_identity_box/samples \
  --timing_path=experiments/sample_identity_box/timings \
  --seed=2022 \
  --n_threads=4