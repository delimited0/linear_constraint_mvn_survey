#!/usr/bin/env bash

Rscript experiments/fernandez_small/fernandez_small.R \
--method_conf=experiments/fernandez_small/lcg_conf.json \
--dim_conf=experiments/fernandez_small/dim_conf.json \
--result_path=experiments/fernandez_small/results \
--seed=2022 \
--n_cores=1 \
--n_threads=18