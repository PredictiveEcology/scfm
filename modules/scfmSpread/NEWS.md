# scfmSpread (development version)

- `rstCurrentBurn` now starts empty in every burn event. It was rebuilt only in a year with an escaped
  fire, so a year whose ignitions all failed to escape kept the last escape year's pixels, and
  Biomass_regeneration regenerated them again. In such a year each ignition pixel now burns
  (`rstCurrentBurn`, `burnMap`, `burnDT`, `timeSinceFire`), as non-escaped ignitions already did in
  years with escapes; `burnSummary` already counted them.
- Fires are capped at their fire regime polygon's `maxBurnCells`. The caps were looked up in a column
  `polyID` that does not exist (it is `PolyID`), so every cap was `NA` and no fire was capped. scfmEscape
  now sets each fire's cap in `spreadState`, and scfmSpread uses it.
- `timeSinceFire` now ages by 1 in a year with no ignitions; it aged only in years in which something
  burned.

# scfmSpread 2.0.0

- completed conversion to using `terra` instead of `raster` objects;
- completed conversion to using `sf` instead of `sp` objects;
- new parameter `dataYear` used to select year for landcover data used for `flammableMap` creation;
- updated object `fireRegimePolys` replaces objects `landscapeAttr` and `scfmDriverPars`, which have been removed;
