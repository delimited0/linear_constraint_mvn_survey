#!/usr/bin/env bash

Rscript experiments/em_sparse_probit/iq_ep_glasso.R \
--output_path=experiments/em_sparse_probit/test_output \
--n_reps=4 \
--max_iter=20 \
--seed=2022 \
--n_cores=4 \
--n_blas_threads=1 & disown