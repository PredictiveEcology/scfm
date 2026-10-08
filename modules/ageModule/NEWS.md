# ageModule 2.1.0

- Released with the scfm family at 2.1.0.
- Deprecated: the scfm parent module does not load it (projects use scfmDataPrep, scfmIgnition, scfmEscape, scfmSpread and scfmDiagnostics). It will be removed in a later release.
- Stands now age each year. The new ages were computed but never written to `ageMap`, so only burned pixels changed (to 0).
- No longer calls quickPlot (not in `reqdPkgs`, and being retired), which stopped `init` unless another package had attached it; scfmutils is required from `@development`, since its `main` branch is at 0.0.11.

# ageModule 2.0.0

- completed conversion to using `terra` instead of `raster` objects;
- completed conversion to using `sf` instead of `sp` objects;

