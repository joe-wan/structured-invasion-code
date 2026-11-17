#!/usr/bin/env bash
# 00_run_all.sh
# =============
# Author: Joe Wan
# Runs all code from the code repository.

# Make sure we are in the code directory
cd -- "$(dirname -- "${BASH_SOURCE[0]}")"

# Run each script
( cd 01_toy_models/ && Rscript 01_simple_model.R )
( cd 01_toy_models/ && Rscript 02_cyclic_model.R )
( cd 01_toy_models/ && Rscript 03_floquet_theory.R )
( cd 01_toy_models/ && Rscript 04_chaotic_model.R )
( cd 02_case_studies/ && Rscript 01_intraguild_predation.R )
( cd 02_case_studies/ && Rscript 02_plant_soil_feedback.R )
( cd 02_case_studies/ && Rscript 03_plant_pest_parasite.R )
