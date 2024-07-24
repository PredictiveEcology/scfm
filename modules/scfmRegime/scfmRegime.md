---
title: "scfmRegime Manual"
subtitle: "v.2.0.0"
date: "Last updated: 2024-07-03"
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
  bibliography: citations/references_scfmRegime.bib
link-citations: true
always_allow_html: true
---

# scfmRegime Module

<!-- the following are text references used in captions for LaTeX compatibility -->
(ref:scfmRegime) *scfmRegime*



[![made-with-Markdown](figures/markdownBadge.png)](https://commonmark.org)

<!-- if knitting to pdf remember to add the pandoc_args: ["--extract-media", "."] option to yml in order to get the badge images -->

#### Authors:

Steve Cumming <stevec@sbf.ulaval.ca> [aut], Ian Eddy <ian.eddy@nrcan-rncan.gc.ca> [aut], Alex M. Chubaty <achubaty@for-cast.ca> [ctb]
<!-- ideally separate authors with new lines, '\n' not working -->

## Module Overview

### Module summary

Estimates fire regime parameters for a given landscape (`studyArea`).

### Module inputs and parameters

Table \@ref(tab:moduleInputs-scfmRegime) shows the full list of module inputs.

<table class="table" style="color: black; margin-left: auto; margin-right: auto;">
<caption>(\#tab:moduleInputs-scfmRegime)(\#tab:moduleInputs-scfmRegime)List of (ref:scfmRegime) input objects and their description.</caption>
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
   <td style="text-align:left;"> firePoints </td>
   <td style="text-align:left;"> sf </td>
   <td style="text-align:left;"> Historical fire data in point form. Must contain fields 'CAUSE', 'YEAR', and 'SIZE_HA', or pass the parameters to identify those. </td>
   <td style="text-align:left;"> http://cwfis.cfs.nrcan.gc.ca/downloads/nfdb/fire_pnt/current_version/NFDB_point.zip </td>
  </tr>
  <tr>
   <td style="text-align:left;"> fireRegimePolys </td>
   <td style="text-align:left;"> sf </td>
   <td style="text-align:left;"> Areas to calibrate individual fire regime parameters. Defaults to ecoregions. Must have numeric field 'PolyID' or it will be created for individual polygons. Must be a sf object. </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> fireRegimePolysCalibration </td>
   <td style="text-align:left;"> sf </td>
   <td style="text-align:left;"> `sf` polygons object with field 'PolyID' describing unique fire regimes in a larger study area. Not required - but useful if the parameterization region is different from the simulation region. </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> rasterToMatch </td>
   <td style="text-align:left;"> SpatRaster </td>
   <td style="text-align:left;"> template raster for raster GIS operations. Must be supplied by user with same CRS as `studyArea`. </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> rasterToMatchCalibration </td>
   <td style="text-align:left;"> SpatRaster </td>
   <td style="text-align:left;"> large template raster for raster GIS operations. Must be supplied by user with same CRS as `studyAreaCalibration`. </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> studyArea </td>
   <td style="text-align:left;"> sf </td>
   <td style="text-align:left;"> Polygon to use as the simulation study area. Can be a `SpatVector`. </td>
   <td style="text-align:left;"> http://sis.agr.gc.ca/cansis/nsdb/ecostrat/district/ecodistrict_shp.zip </td>
  </tr>
  <tr>
   <td style="text-align:left;"> studyAreaCalibration </td>
   <td style="text-align:left;"> sf </td>
   <td style="text-align:left;"> Polygon to use as the parametrisation study area. Can be a `SpatVector`. Note that `studyAreaCalibration` is only used for parameter estimation, and can be larger than the actual study area used for simulations. </td>
   <td style="text-align:left;"> http://sis.agr.gc.ca/cansis/nsdb/ecostrat/district/ecodistrict_shp.zip </td>
  </tr>
</tbody>
</table>

Table \@ref(tab:moduleParams-scfmRegime) shows the full list of module parameters.

<table class="table" style="color: black; margin-left: auto; margin-right: auto;">
<caption>(\#tab:moduleParams-scfmRegime)(\#tab:moduleParams-scfmRegime)List of (ref:scfmRegime) parameters and their description.</caption>
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
   <td style="text-align:left;"> empiricalMaxSizeFactor </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> 1.2 </td>
   <td style="text-align:left;"> 1 </td>
   <td style="text-align:left;"> 10 </td>
   <td style="text-align:left;"> scale `xMax` by this if HD estimator fails </td>
  </tr>
  <tr>
   <td style="text-align:left;"> fireCause </td>
   <td style="text-align:left;"> character </td>
   <td style="text-align:left;"> N </td>
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
   <td style="text-align:left;"> .useCache </td>
   <td style="text-align:left;"> logical </td>
   <td style="text-align:left;"> FALSE </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Internal. Can be names of events or the whole module name to be cached by SpaDES. </td>
  </tr>
</tbody>
</table>

### Module outputs

Description of the module outputs (Table \@ref(tab:moduleOutputs-scfmRegime)).

<table class="table" style="color: black; margin-left: auto; margin-right: auto;">
<caption>(\#tab:moduleOutputs-scfmRegime)(\#tab:moduleOutputs-scfmRegime)List of (ref:scfmRegime) outputs and their description.</caption>
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
   <td style="text-align:left;"> Fire locations. Points outside `studyArea` are removed </td>
  </tr>
  <tr>
   <td style="text-align:left;"> fireRegimePolys </td>
   <td style="text-align:left;"> sf </td>
   <td style="text-align:left;"> `fireRegimePolys` with fire attributes appended. </td>
  </tr>
</tbody>
</table>

### Links to other modules

Intended to be run with the `scfm` suite of modules found at <https://github.com/PredictiveEcology/scfm>:

- `ageModule` (optional)
- `scfmDriver`
- `scfmEscape`
- `scfmIgnition`
- `scfmLandcoverInit`
- `scfmRegime`
- `scfmSpread`

### Getting help

<https://github.com/PredictiveEcology/scfm/issues>
