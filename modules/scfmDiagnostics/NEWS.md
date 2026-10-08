# scfmDiagnostics 2.1.0

- Released with the scfm family at 2.1.0.
- `terra` and `data.table` are in `reqdPkgs`; the summary used their functions without them and stopped ("could not find function expanse") unless another module had attached them. quickPlot's `clearPlot()` is no longer called.

# scfmDiagnostics 2.0.0

- completed conversion to using `terra` instead of `raster` objects;
- completed conversion to using `sf` instead of `sp` objects;
- updated objects `fireRegimePolys` replaces objects `landscapeAttr`, `scfmDriverPars` and `scfmRegimePars`, which have been removed;