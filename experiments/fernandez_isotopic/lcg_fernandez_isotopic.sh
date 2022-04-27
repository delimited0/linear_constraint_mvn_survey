#!/usr/bin/env bash

Rscript experiments/fernandez_isotopic/fernandez_isotopic.R \
--method_conf=experiments/fernandez_isotopic/lcg_conf.json \
--dim_conf=experiments/fernandez_isotopic/dim_conf.json \
--result_path=experiments/fernandez_isotopic/results \
--seed=2022 \
--n_cores=1 \
--n_blas_threads=18