# Licence and provenance of the data in this directory

## Licence

The datasets in this directory are made available under the
**Creative Commons Attribution 4.0 International licence (CC BY 4.0)**.

Full text: <https://creativecommons.org/licenses/by/4.0/legalcode>
Summary:   <https://creativecommons.org/licenses/by/4.0/>

You are free to share and adapt this material for any purpose, including
commercially, provided you give appropriate credit, link to the licence, and
indicate whether changes were made.

This matches the licence of the original publication in which these data were
first released (Frontiers in Marine Science is fully open access under
CC BY 4.0).

## Provenance

The aerial survey data (`ips_01_21_noFN10.rds`) and the spatial mesh
(`mesh_5km.rds`) derive from the survey programme and modelling framework
described in:

> Vilela, R., Burger, C., Diederichs, A., Bachl, F.E., Szostek, L., Freund, A.,
> Braasch, A., Bellebaum, J., Beckers, B., Piper, W. and Nehls, G. (2021).
> Use of an INLA latent Gaussian modeling approach to assess bird population
> changes due to the development of offshore wind farms.
> *Frontiers in Marine Science*, 8, 701332.
> <https://doi.org/10.3389/fmars.2021.701332>

Counts are **aggregated to the nodes of the spatial mesh**, preserving the
survey method (HiDef, DAISI, APEM, visual) and the time period. They are not
raw sighting positions. Conventional visual counts have been corrected for
detection probability using distance sampling (`mrds`); see the Methods of the
accompanying manuscript.

`farms18_spring_byyear_updt.rds` contains offshore wind farm footprints in the
German Bight with construction and operation dates. `mask_zee_full.*` is the
German EEZ boundary used to clip the prediction grid.

## How to cite

If you use these data, please cite Vilela et al. (2021) above, together with
the manuscript accompanying this repository.
