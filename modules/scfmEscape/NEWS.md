# scfmEscape 2.1.0

- Released with the scfm family at 2.1.0.

- `spreadState` is `NULL` in a year with no ignitions. It kept the previous year's state, so scfmSpread
  burned the previous year's fires again and counted them again in `burnSummary`.
- Each fire's size cap (`maxBurnCells` of its fire regime polygon) is set in `spreadState` when the fire
  starts, so the escape's pixels count towards it. A cap first given when scfmSpread resumed the fire
  counted each fire as 1 pixel, and was matched to fires in pixel order, not ignition order.

# scfmEscape 2.0.0

- completed conversion to using `terra` instead of `raster` objects;
- completed conversion to using `sf` instead of `sp` objects;
- new parameter `dataYear` used to select year for landcover data used for `flammableMap` creation;
- updated objects `fireRegimePolys` replaces objects `landscapeAttr` and `scfmDriverPars`, which has been removed;
- object `rasterToMatch` was removed;

