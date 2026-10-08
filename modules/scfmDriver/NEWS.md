# scfmDriver 2.1.0

- Released with the scfm family at 2.1.0.
- Deprecated: scfmDataPrep runs this step (its `eventsToPrepare` parameter), and the scfm parent module does not load this module. It will be removed in a later release.

# scfmDriver 2.0.0

- completed conversion to using `terra` instead of `raster` objects;
- completed conversion to using `sf` instead of `sp` objects;
- removed parameters `bufferLCCYear` and `neighbours`;
- new parameter `dataYear` used to select year for landcover data used for `flammableMap` creation;
- updated objects `fireRegimePolys` replaces objects `landscapeAttr`, `scfmDriverPars` and `scfmRegimePars`, which have been removed;

