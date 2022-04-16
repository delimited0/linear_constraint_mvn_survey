#!/usr/bin/env bash

export R_PROGRESSR_ENABLE=TRUE

Rscript experiments/proton_radius/proton_radius.R \
  --method_conf=experiments/proton_radius/method_conf.json \
  --sample_path=experiments/proton_radius/samples \
  --timing_path=experiments/proton_radius/timings \
  --reps=2 \
  --seed=2022 \
  --n_threads=4