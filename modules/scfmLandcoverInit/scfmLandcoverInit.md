---
title: "scfmLandcoverInit Manual"
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
  bibliography: citations/references_scfmLandcoverInit.bib
link-citations: true
always_allow_html: true
---

# scfmLandcoverInit Module

<!-- the following are text references used in captions for LaTeX compatibility -->
(ref:scfmLandcoverInit) *scfmLandcoverInit*



[![made-with-Markdown](figures/markdownBadge.png)](https://commonmark.org)

<!-- if knitting to pdf remember to add the pandoc_args: ["--extract-media", "."] option to yml in order to get the badge images -->

#### Authors:

Eliot J B McIntire <eliot.mcintire@nrcan-rncan.gc.ca> [aut, cre], Steve Cumming <stevec@sbf.ulaval.ca> [aut], Ian Eddy <ian.eddy@nrcan-rncan.gc.ca> [aut], Alex M. Chubaty <achubaty@for-cast.ca> [ctb]
<!-- ideally separate authors with new lines, '\n' not working -->

## Module Overview

### Module summary

Generates some relevant statistics for each fire regime over a `studyArea`.
If scfm is being parameterized over a larger area (`studyAreaCalibration`), then the following objects must be supplied with identical CRS and resolution, where applicable:

- `studyArea` and `studyAreaCalibration`;
- `rasterToMatch` and `rasterToMatchCalibration`.

The extents should differ between objects and their 'large' counterparts.

### Module inputs and parameters

Table \@ref(tab:moduleInputs-scfmLandcoverInit) shows the full list of module inputs.

<table class="table" style="color: black; margin-left: auto; margin-right: auto;">
<caption>(\#tab:moduleInputs-scfmLandcoverInit)(\#tab:moduleInputs-scfmLandcoverInit)List of (ref:scfmLandcoverInit) input objects and their description.</caption>
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
   <td style="text-align:left;"> binary flammability map - defaults to using LandR::prepInputsLCC </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> flammableMapCalibration </td>
   <td style="text-align:left;"> SpatRaster </td>
   <td style="text-align:left;"> binary flammability map - defaults to using `LandR::prepInputsLCC`. This is only necessary if passing `studyAreaCalibration` OR running `scfmDriver`. It should match the extent of `studyAreaCalibration`, and if running `scfmDriver`, it should extend by &gt;= scfmDriver's `P(sim)$buffDist`. </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> rasterToMatch </td>
   <td style="text-align:left;"> SpatRaster </td>
   <td style="text-align:left;"> template raster for raster GIS operations. Must be supplied by user </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> rasterToMatchCalibration </td>
   <td style="text-align:left;"> SpatRaster </td>
   <td style="text-align:left;"> Template raster for raster GIS operations. Only necessary if `studyAreaCalibration` is passed. Must be supplied by user. </td>
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

Table \@ref(tab:moduleParams-scfmLandcoverInit) shows the full list of module parameters.


<table class="table" style="color: black; margin-left: auto; margin-right: auto;">
<caption>(\#tab:moduleParams-scfmLandcoverInit)(\#tab:moduleParams-scfmLandcoverInit)List of (ref:scfmLandcoverInit) parameters and their description.</caption>
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
   <td style="text-align:left;"> dataYear </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> 2011 </td>
   <td style="text-align:left;"> 1985 </td>
   <td style="text-align:left;"> 2020 </td>
   <td style="text-align:left;"> used to select the year of landcover data used to create flammableMap if the obejct is unsupplied </td>
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
   <td style="text-align:left;"> neighbours </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> 8 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Number of immediate cell neighbours </td>
  </tr>
  <tr>
   <td style="text-align:left;"> sliverThreshold </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> 1e+08 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> fire regime polygons with area less than this number will be merged with their closest non-sliver neighbour using `sf::st_nearest_feature`. </td>
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
</tbody>
</table>

### Module outputs

Description of the module outputs (Table \@ref(tab:moduleOutputs-scfmLandcoverInit)).

<table class="table" style="color: black; margin-left: auto; margin-right: auto;">
<caption>(\#tab:moduleOutputs-scfmLandcoverInit)(\#tab:moduleOutputs-scfmLandcoverInit)List of (ref:scfmLandcoverInit) outputs and their description.</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> objectName </th>
   <th style="text-align:left;"> objectClass </th>
   <th style="text-align:left;"> desc </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> fireRegimePolys </td>
   <td style="text-align:left;"> sf </td>
   <td style="text-align:left;"> `fireRegimePolys` with landcover attributes appended </td>
  </tr>
  <tr>
   <td style="text-align:left;"> fireRegimePolysCalibration </td>
   <td style="text-align:left;"> sf </td>
   <td style="text-align:left;"> `fireRegimePolysCalibration` with landcover attributes appended </td>
  </tr>
  <tr>
   <td style="text-align:left;"> fireRegimeRas </td>
   <td style="text-align:left;"> SpatRaster </td>
   <td style="text-align:left;"> Rasterized version of fireRegimePolys with values representing polygon ID </td>
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
