# Licence and provenance of the data in this directory

## Licence

The data in this directory are made available under the
**Creative Commons Attribution 4.0 International licence (CC BY 4.0)**.

Full text: <https://creativecommons.org/licenses/by/4.0/legalcode>
Summary:   <https://creativecommons.org/licenses/by/4.0/>

You are free to share and adapt this material for any purpose, including
commercially, provided you give appropriate credit, link to the licence, and
indicate whether changes were made.

## Contents

| File | What it is |
|---|---|
| `diver_owf_data.RData` | Every input the analysis reads: aerial survey counts, OWF footprints, the prediction grid and its mask, the two reference areas, and the SPDE mesh. The objects are listed in the README. |
| `pred_change_975.tif` | Not survey data. See below. |

## Provenance of the survey data

The counts and the spatial mesh derive from the aerial survey programme and
modelling framework described in:

> Vilela, R., Burger, C., Diederichs, A., Bachl, F.E., Szostek, L., Freund, A.,
> Braasch, A., Bellebaum, J., Beckers, B., Piper, W. and Nehls, G. (2021).
> Use of an INLA latent Gaussian modeling approach to assess bird population
> changes due to the development of offshore wind farms.
> *Frontiers in Marine Science*, 8, 701332.
> <https://doi.org/10.3389/fmars.2021.701332>

The records differ in when they were first released:

- **Spring 2001–2018** were released with Vilela et al. (2021), published open
  access under CC BY 4.0. The licence above follows from that publication.
- **Spring 2019–2021** were not part of that dataset and are released here for
  the first time, under the same licence and with the permission of the survey
  programme.

Counts are **aggregated to the nodes of the spatial mesh**, preserving the
survey technique (HiDef, DAISI, APEM, visual) and the time period. They are not
raw sighting positions. Conventional visual counts have been corrected for
detection probability using distance sampling (`mrds`); see the Methods of the
accompanying manuscript.

## `pred_change_975.tif`

This file is not survey data. It is the change field produced by the earlier
version of the analysis, on the same prediction grid, archived so that the
comparison reported in the response to reviewers can be reproduced. It is a
model output rather than an observation, and is included for verification
only. It is read by `R/diagnose_vs_original_tif.R` and
`R/diagnose_distance_surfaces.R`.

## How to cite

If you use these data, please cite Vilela et al. (2021) above, together with
the manuscript accompanying this repository.
