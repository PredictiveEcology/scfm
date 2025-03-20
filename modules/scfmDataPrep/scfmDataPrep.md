---
title: "scfmDataPrep Manual"
subtitle: "v.0.0.0.9000"
date: "Last updated: 2025-03-17"
output:
  bookdown::html_document2:
    toc: true
    toc_float: true
    theme: sandstone
    number_sections: false
    df_print: paged
    keep_md: yes
editor_options:
  chunk_output_type: console
  bibliography: citations/references_scfmDataPrep.bib
link-citations: true
always_allow_html: true
---

# scfmDataPrep Module

<!-- the following are text references used in captions for LaTeX compatibility -->
(ref:scfmDataPrep) *scfmDataPrep*



[![made-with-Markdown](figures/markdownBadge.png)](https://commonmark.org)

<!-- if knitting to pdf remember to add the pandoc_args: ["--extract-media", "."] option to yml in order to get the badge images -->

#### Authors:

Eliot J B McIntire <eliot.mcintire@nrcan-rncan.gc.ca> [aut], Steve Cumming <stevec@sbf.ulaval.ca> [aut], Ian Eddy <ian.eddy@nrcan-rncan.gc.ca> [aut, cre], Alex M. Chubaty <achubaty@for-cast.ca> [ctb]
<!-- ideally separate authors with new lines, '\n' not working -->

## Module Overview

### Module summary

Generates relevant statistics and estimates parameters for a generic percolation model to simulate fire regime parameters for a given landscape (`studyArea`).

If `scfm` is being parameterized over a larger area (`studyAreaCalibration`), then the following objects must be supplied with identical CRS and resolution, where applicable:

- `studyArea` and `studyAreaCalibration`;
- `rasterToMatch` and `rasterToMatchCalibration`.

The extents should differ between objects and their larger `*Calibration` counterparts.

### Module inputs and parameters

Describe input data required by the module and how to obtain it (e.g., directly from online sources or supplied by other modules)
If `sourceURL` is specified, `downloadData("scfmDataPrep", "..")` may be sufficient.
Table \@ref(tab:moduleInputs-scfmDataPrep) shows the full list of module inputs.

<table class="table" style="color: black; margin-left: auto; margin-right: auto;">
<caption>(\#tab:moduleInputs-scfmDataPrep)(\#tab:moduleInputs-scfmDataPrep)List of (ref:scfmDataPrep) input objects and their description.</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> objectName </th>
   <th style="text-align:left;"> objectClass </th>
   <th style="text-align:left;"> desc </th>
   <th style="text-align:left;"> sourceURL </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> cloudFolderID </td>
   <td style="text-align:left;"> character </td>
   <td style="text-align:left;"> URL for Google-drive-backed cloud cache. Note: turn `cloudCache` on or off with `options('reproducible.useCloud')`. </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> firePoints </td>
   <td style="text-align:left;"> sf </td>
   <td style="text-align:left;"> Historical fire data in point form. Must contain fields 'CAUSE', 'YEAR', and 'SIZE_HA', or pass the parameters to identify those. </td>
   <td style="text-align:left;"> http://cwfis.cfs.nrcan.gc.ca/downloads/nfdb/fire_pnt/current_version/NFDB_point.zip </td>
  </tr>
  <tr>
   <td style="text-align:left;"> fireRegimePolys </td>
   <td style="text-align:left;"> sf </td>
   <td style="text-align:left;"> Areas to calibrate individual fire regime parameters. Defaults to ecozones of Canada. Must have numeric field 'PolyID' or it will be created for individual polygons. </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> fireRegimePolysCalibration </td>
   <td style="text-align:left;"> sf </td>
   <td style="text-align:left;"> if `studyAreaCalibration` is supplied, the corresponding fire regime areas. Requires integer field `PolyID` if supplied. Uses same defaults as `fireRegimePolys`. </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> flammableMap </td>
   <td style="text-align:left;"> SpatRaster </td>
   <td style="text-align:left;"> binary flammability map - defaults to using `LandR::prepInputsLCC` </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> flammableMapCalibration </td>
   <td style="text-align:left;"> SpatRaster </td>
   <td style="text-align:left;"> binary flammability map corresponding to `rasterToMatchCalibration`. It should extent from `studyArea` by &gt;= scfmDriver's `P(sim)$buffDist`. and if unsupplied, will be created using `LandR::prepInputs_NTEMS_LCC_FAO` </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> rasterToMatch </td>
   <td style="text-align:left;"> SpatRaster </td>
   <td style="text-align:left;"> template raster for raster GIS operations. Must be supplied by user. </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> rasterToMatchCalibration </td>
   <td style="text-align:left;"> SpatRaster </td>
   <td style="text-align:left;"> Template raster for `studyAreaCalibration`. Will be created based on `rasterToMatch` if unsupplied. </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> studyArea </td>
   <td style="text-align:left;"> sf </td>
   <td style="text-align:left;"> Polygon to use as the simulation study area (typically buffered). </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> studyAreaCalibration </td>
   <td style="text-align:left;"> sf </td>
   <td style="text-align:left;"> optional larger study area used for parameterization only </td>
   <td style="text-align:left;"> NA </td>
  </tr>
</tbody>
</table>

Provide a summary of user-visible parameters (Table \@ref(tab:moduleParams-scfmDataPrep))


<table class="table" style="color: black; margin-left: auto; margin-right: auto;">
<caption>(\#tab:moduleParams-scfmDataPrep)(\#tab:moduleParams-scfmDataPrep)List of (ref:scfmDataPrep) parameters and their description.</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> paramName </th>
   <th style="text-align:left;"> paramClass </th>
   <th style="text-align:left;"> default </th>
   <th style="text-align:left;"> min </th>
   <th style="text-align:left;"> max </th>
   <th style="text-align:left;"> paramDesc </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> buffDist </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> 20000 </td>
   <td style="text-align:left;"> 1 </td>
   <td style="text-align:left;"> 1e+05 </td>
   <td style="text-align:left;"> Buffer width to mitigate edge effects in fire landscape calibration. If studyAreaCalibration is not supplied, this parameter will also be used to create it via buffering studyArea </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cloudFolderID </td>
   <td style="text-align:left;"> character </td>
   <td style="text-align:left;">  </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> URL for Google-drive-backed cloud cache </td>
  </tr>
  <tr>
   <td style="text-align:left;"> dataYear </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> 2011 </td>
   <td style="text-align:left;"> 1985 </td>
   <td style="text-align:left;"> 2020 </td>
   <td style="text-align:left;"> used to select the year of landcover data used to create flammableMap if the obejct is unsupplied </td>
  </tr>
  <tr>
   <td style="text-align:left;"> empiricalMaxSizeFactor </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> 1.2 </td>
   <td style="text-align:left;"> 1 </td>
   <td style="text-align:left;"> 10 </td>
   <td style="text-align:left;"> scale `xMax` by this if HD estimator fails </td>
  </tr>
  <tr>
   <td style="text-align:left;"> eventsToPrepare </td>
   <td style="text-align:left;"> character </td>
   <td style="text-align:left;"> scfmLand.... </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> which of three possible dataPrep steps to run? scfmLandcoverInit will clean GIS and generate landscape stats regarding flammability in each fire regime poly; scfmRegime will prepare fire regime attributes inc. mean fire size and ignition rate; scfmDriver will estimate the spread probability of flammable pixels </td>
  </tr>
  <tr>
   <td style="text-align:left;"> fireCause </td>
   <td style="text-align:left;"> character </td>
   <td style="text-align:left;"> N, L </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> subset of `c('H', 'H-PB', 'N', 'Re', 'U')` </td>
  </tr>
  <tr>
   <td style="text-align:left;"> fireCauseColumnName </td>
   <td style="text-align:left;"> character </td>
   <td style="text-align:left;"> CAUSE </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Name of the column that has fire cause, consistent with `P(sim)$fireCause`. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> fireEpoch </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> 1971, 2000 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> start of normal period </td>
  </tr>
  <tr>
   <td style="text-align:left;"> fireRegimePolysType </td>
   <td style="text-align:left;"> character </td>
   <td style="text-align:left;"> ECOREGION </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Polygon type to use for scfm `fireRegimePolys`: see `?scfmutils::prepInputsFireRegimePolys` for allowed types. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> fireSizeColumnName </td>
   <td style="text-align:left;"> character </td>
   <td style="text-align:left;"> SIZE_HA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Name of the column that has fire size </td>
  </tr>
  <tr>
   <td style="text-align:left;"> fireYearColumnName </td>
   <td style="text-align:left;"> character </td>
   <td style="text-align:left;"> YEAR </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Name of the column that has fire size </td>
  </tr>
  <tr>
   <td style="text-align:left;"> flammabilityThreshold </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> 0.25 </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> 1 </td>
   <td style="text-align:left;"> Minimum proportion of flammable old pixel needed to define a new pixel as flammable when upscaling the default flammable maps. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> neighbours </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> 8 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Number of immediate cell neighbours </td>
  </tr>
  <tr>
   <td style="text-align:left;"> pJmp </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> 0.23 </td>
   <td style="text-align:left;"> 0.18 </td>
   <td style="text-align:left;"> 0.25 </td>
   <td style="text-align:left;"> default spread prob for degenerate polygons </td>
  </tr>
  <tr>
   <td style="text-align:left;"> pMax </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> 0.253 </td>
   <td style="text-align:left;"> 0.24 </td>
   <td style="text-align:left;"> 0.26 </td>
   <td style="text-align:left;"> maximum spread range for calibration </td>
  </tr>
  <tr>
   <td style="text-align:left;"> pMin </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> 0.185 </td>
   <td style="text-align:left;"> 0.15 </td>
   <td style="text-align:left;"> 0.225 </td>
   <td style="text-align:left;"> minimum spread range for calibration </td>
  </tr>
  <tr>
   <td style="text-align:left;"> scamOptimizer </td>
   <td style="text-align:left;"> character </td>
   <td style="text-align:left;"> bfgs </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> numerical optimization method used in fitting scam model; see `?scam`. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> sliverThreshold </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> 6.25e+08 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> fire regime polygons with area (in m2) less than this number will be merged with their closest non-sliver neighbour using `sf::st_nearest_feature`. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> targetBurnRate </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> 1 </td>
   <td style="text-align:left;"> a named vector giving the proportional annual area burned of each fire regime polygon. These override the default estimate of scfm and are used to estimate a new mean fire size and ignition rate. Names should correspond to `PolyID`. A partial set of polygons is allowed - missing polys are estimated from data. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> targetMaxFireSize </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> a named vector giving the estimated max fire size (in $ha$) of each fire regime polygon. These will override the default estimate of scfm and will be used to estimate a new spread probability. Names should correspond to `PolyID`. A partial set of polygons is allowed - missing polys are estimated from data. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> targetN </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> 4000 </td>
   <td style="text-align:left;"> 1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> target sample size for determining true spread probability </td>
  </tr>
  <tr>
   <td style="text-align:left;"> .plotInitialTime </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Initial time for plotting </td>
  </tr>
  <tr>
   <td style="text-align:left;"> .plotInterval </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Interval between plotting </td>
  </tr>
  <tr>
   <td style="text-align:left;"> .plots </td>
   <td style="text-align:left;"> character </td>
   <td style="text-align:left;"> screen, png </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Used by `Plots` function, which can be optionally used here. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> .saveInitialTime </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Initial time for saving </td>
  </tr>
  <tr>
   <td style="text-align:left;"> .saveInterval </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Interval between save events </td>
  </tr>
  <tr>
   <td style="text-align:left;"> .useCache </td>
   <td style="text-align:left;"> character </td>
   <td style="text-align:left;"> .inputOb.... </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Use caching of events - not recommended as of 10/05/2023 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> .useCloud </td>
   <td style="text-align:left;"> logical </td>
   <td style="text-align:left;"> FALSE </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> should a cloud cache be used for heavy operations </td>
  </tr>
  <tr>
   <td style="text-align:left;"> .useParallelFireRegimePolys </td>
   <td style="text-align:left;"> logical </td>
   <td style="text-align:left;"> FALSE </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> should driver use parallel? Alternatively accepts a numeric argument, i.e., how many cores. </td>
  </tr>
</tbody>
</table>

### Events

Describe what happens for each event type.

### Plotting

Write what is plotted.

### Saving

Write what is saved.

### Module outputs

Description of the module outputs (Table \@ref(tab:moduleOutputs-scfmDataPrep)).

<table class="table" style="color: black; margin-left: auto; margin-right: auto;">
<caption>(\#tab:moduleOutputs-scfmDataPrep)(\#tab:moduleOutputs-scfmDataPrep)List of (ref:scfmDataPrep) outputs and their description.</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> objectName </th>
   <th style="text-align:left;"> objectClass </th>
   <th style="text-align:left;"> desc </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> fireRegimePoints </td>
   <td style="text-align:left;"> sf </td>
   <td style="text-align:left;"> Fire locations. These are filtered according to criteria set in params (i.e. epoch, cause) </td>
  </tr>
  <tr>
   <td style="text-align:left;"> fireRegimePolys </td>
   <td style="text-align:left;"> sf </td>
   <td style="text-align:left;"> `fireRegimePolys` with fire attributes appended. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> fireRegimePolysCalibration </td>
   <td style="text-align:left;"> sf </td>
   <td style="text-align:left;"> `fireRegimePolysCalibration` with attributes appended </td>
  </tr>
  <tr>
   <td style="text-align:left;"> fireRegimeRas </td>
   <td style="text-align:left;"> SpatRaster </td>
   <td style="text-align:left;"> Rasterized version of fireRegimePolys with values representing polygon ID </td>
  </tr>
</tbody>
</table>

### Links to other modules

### Links to other modules

Intended to be run with the `scfm` suite of modules found at <https://github.com/PredictiveEcology/scfm>:

- `ageModule` (optional)
- `scfmDataPrep`
- `scfmEscape`
- `scfmIgnition`
- `scfmSpread`


### Getting help

-   provide a way for people to obtain help (e.g., module repository issues page)
