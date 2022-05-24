#!/usr/bin/env bash

Rscript experiments/em_sparse_probit/iq_ep_glasso.R \
--output_path=experiments/em_sparse_probit/output \
--n_reps=18 \
--max_iter=300 \
--seed=2022 \
--n_cores=18 \
--n_blas_threads=1 & disown