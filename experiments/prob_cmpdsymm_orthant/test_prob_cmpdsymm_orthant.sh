#!/usr/bin/env bash

Rscript experiments/prob_cmpdsymm_orthant/prob_cmpdsymm_orthant.R \
--method_conf=experiments/prob_cmpdsymm_orthant/method_conf.json \
--dim_conf=experiments/prob_cmpdsymm_orthant/test_dim_conf.json \
--corr=.5 \
--result_path=experiments/prob_cmpdsymm_orthant/test_results \
--seed=2022 \
--n_cores=1 \
--n_blas_threads=4