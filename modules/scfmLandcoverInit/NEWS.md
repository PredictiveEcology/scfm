# scfmLandcoverInit 2.0.0

- completed conversion to using `terra` instead of `raster` objects;
- completed conversion to using `sf` instead of `sp` objects;
- new parameter `dataYear` used to select year for landcover data used for `flammableMap` creation;
- new parameter `fireRegimePolysType` used to select data source for `fireRegimePolys` creation;
- parameter `fireCause` uses new default `"N"` for naturally occuring ignitions, instead of the previous `"L"` (lightning-caused), due to a recent change in specification in the National Fire Database datasets;
- updated objects `fireRegimePolys` and `fireRegimePolysCalibration` replace objects `landscapeAttr` and `landscapeAttrLarge`, respectively, which have been removed;
- object `cellsByZone` has been removed, as this attribute is now stored in `fireRegimePolys` and `fireRegimePolysCalibration`;

