# scfm 2.1.0

This is the first release of scfm as one family under one version number. The repository itself is now the scfm parent module, so `PredictiveEcology/scfm@v2.1.0` loads the five scfm modules in use from this release in one step, and every module in it carries the version 2.1.0. A new scfmDataPrep module prepares the fire regime inputs in one place, the modules use terra and sf throughout, and run names can be added to plot titles and burn maps.

Several fixes change simulated fire. A year whose ignitions all failed to escape no longer repeats the previous year's burned area, fires now respect each fire regime's maximum size, and time since fire ages every year, including years without fire. Fire regime estimates now draw on a longer record of past fires, and the default flammability threshold is lower. ageModule now ages stands, which it never did. ageModule, scfmDriver, scfmLandcoverInit and scfmRegime are deprecated: scfmDataPrep does the data preparation, and the parent loads only the five modules in use. Projects that loaded the scfm modules one by one from `scfm@development/modules/...` can keep doing so.

Each module also has its own NEWS.md in `modules/<module>/`. Changes before 2024-06-03 are not included.
