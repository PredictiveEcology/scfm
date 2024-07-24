# scfmRegime 2.0.0

- completed conversion to using `terra` instead of `raster` objects;
- completed conversion to using `sf` instead of `sp` objects;
- parameter `fireCause` uses new default `"N"` for naturally occuring ignitions, instead of the previous `"L"` (lightning-caused), due to a recent change in specification in the National Fire Database datasets;
- updated objects `fireRegimePolys` and `fireRegimePolysCalibration` replace objects `landscapeAttr` and `landscapeAttrLarge`, respectively, which have been removed;
- updated object `fireRegimePolys` replaces object `scfmRegimePars`, which has been removed;

