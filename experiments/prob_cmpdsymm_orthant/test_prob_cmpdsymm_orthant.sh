#!/usr/bin/env bash

Rscript experiments/prob_cmpdsymm_orthant/prob_cmpdsymm_orthant.R \
--method_conf=experiments/prob_cmpdsymm_orthant/method_conf.json \
--dim_conf=experiments/prob_cmpdsymm_orthant/test_dim_conf.json \
--corr=.5 \
--result_path=experiments/prob_cmpdsymm_orthant/test_results \
--n_reps=4 \
--seed=2022 \
--n_cores=4 \
--n_blas_threads=1