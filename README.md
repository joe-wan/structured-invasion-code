# Code for: "Invasion analysis for structured populations: analytical principles and numerical workflows"

This repository contains the original analysis code for the following
manuscript:

> Invasion analysis for structured populations: analytical principles and numerical workflows\
> Joe Wan, Hsiang-Chih Lo, Guan-Yu Chen, Shiang-Chi Huang, Jia-Zhen Liu, Chin-Te Tsai, Po-Ju Ke\
> Preprint (bioRxiv). DOI: [10.64898/2026.11.04.763670v1](https://www.biorxiv.org/content/10.64898/2026.11.04.763670)


## File Structure
Each of the two sets of analyses is in its own directory. Run each script with
its containing directory as the working directory.

The contents of this directory are:

- `01_toy_models/` contains the two-stage toy model, cyclic resident model,
  Floquet analysis, and chaotic resident model. These analyses produce main text
  Figures 1, 2, 4, 5, and 6 and Supplemental Figures S8.1–S8.6.
- `02_case_studies/` contains the intraguild predation/cannibalism,
  plant–soil feedback, and plant–pest–parasite analyses, producing main text
  Figures 3b, 3d, and 3f and Supplemental Figure S6.1.
  `00_simulation_igr.R` provides the shared numerical invasion growth rate helper (Algorithm S3.1).
- `00_plot_settings.R` defines shared colors, themes, and figure export settings.
- `00_run_all.sh` runs all scripts in the repository.
- `outputs/` contains the exported figures from both analysis directories.


## Software
All analyses were conducted in [R](https://www.r-project.org/) (version 4.4.0) using the following package versions:

- [cowplot](https://wilkelab.org/cowplot/) (v1.1.3)
- [deSolve](https://cran.r-project.org/package=deSolve) (v1.40)
- [dplyr](https://dplyr.tidyverse.org/) (v1.1.4)
- [expm](https://cran.r-project.org/package=expm) (v1.0-0)
- [forcats](https://forcats.tidyverse.org/) (v1.0.0)
- [ggh4x](https://teunbrand.github.io/ggh4x/) (v0.3.1)
- [ggplot2](https://ggplot2.tidyverse.org/) (v4.0.0)
- [magrittr](https://magrittr.tidyverse.org/) (v2.0.3)
- [purrr](https://purrr.tidyverse.org/) (v1.0.2)
- [scales](https://scales.r-lib.org/) (v1.4.0)
- [stringr](https://stringr.tidyverse.org/) (v1.5.1)
- [tibble](https://tibble.tidyverse.org/) (v3.2.1)
- [tidyr](https://tidyr.tidyverse.org/) (v1.3.0)
- [tidyverse](https://tidyverse.tidyverse.org/) (v2.0.0)


## Author Contributions
JW and HCL edited the source code for public release.
JW and PJK conceptualized all analyses, with input from all authors.
JW designed and implemented the toy model analyses.
HCL implemented the simulation-based IGR algorithm, designed and implemented the intraguild predation model, and standardized the code for case studies.
SCH, JZL, and CTT implemented the plant–soil feedback model.
GYC designed and implemented the plant–pest–parasite model.
