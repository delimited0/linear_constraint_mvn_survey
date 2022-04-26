#!/usr/bin/env bash

Rscript experiments/prob_exp_covariance/prob_exp_covariance.R \
--method_conf=experiments/prob_exp_covariance/method_conf.json \
--dim_conf=experiments/prob_exp_covariance/test_dim_conf.json \
--dep=.1 \
--result_path=experiments/prob_exp_covariance/test_results \
--seed=2022 \
--n_cores=1 \
--blas_threads=4