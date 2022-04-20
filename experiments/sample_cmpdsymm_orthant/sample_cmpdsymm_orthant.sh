#!/usr/bin/env bash

export R_PROGRESSR_ENABLE=TRUE

Rscript experiments/sample_cmpdsymm_orthant/sample_cmpdsymm_orthant.R \
--method_conf=experiments/sample_cmpdsymm_orthant/method_conf.json \
--dim_conf=experiments/sample_cmpdsymm_orthant/dim_conf.json \
--sample_path=experiments/sample_cmpdsymm_orthant/samples \
--seed=2022 \
--n_threads=18 &