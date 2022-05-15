#!/usr/bin/env bash

Rscript experiments/em_sparse_probit/iq_glasso.R \
  --method_conf=experiments/em_sparse_probit/method_conf.json \
  --output_path=experiments/em_sparse_probit/output \
  --n_reps=18 \
  --n_mc_samples=50 \
  --seed=2022 \
  --n_cores=18 \
  --n_blas_threads=1