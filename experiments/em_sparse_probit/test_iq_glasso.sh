#!/usr/bin/env bash

Rscript experiments/em_sparse_probit/iq_glasso.R \
--method_conf=experiments/em_sparse_probit/method_conf.json \
--output_path=experiments/em_sparse_probit/test_output \
--n_reps=4 \
--n_mc_samples=5 \
--max_iter=10 \
--seed=2022 \
--n_cores=4 \
--n_blas_threads=1