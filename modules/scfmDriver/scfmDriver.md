---
title: "scfmDriver Manual"
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
  bibliography: citations/references_scfmDriver.bib
link-citations: true
always_allow_html: true
---

# scfmDriver Module

<!-- the following are text references used in captions for LaTeX compatibility -->
(ref:scfmDriver) *scfmDriver*



[![made-with-Markdown](figures/markdownBadge.png)](https://commonmark.org)

<!-- if knitting to pdf remember to add the pandoc_args: ["--extract-media", "."] option to yml in order to get the badge images -->

#### Authors:

Steve G Cumming <stevec@sbf.ulaval.ca> [aut, cre], Ian Eddy <ian.eddy@nrcan-rncan.gc.ca> [aut], Alex M Chubaty <achubaty@for-cast.ca> [ctb]
<!-- ideally separate authors with new lines, '\n' not working -->

## Module Overview

### Module summary

Estimates parameters for the generic percolation model.

### Module inputs and parameters

Table \@ref(tab:moduleInputs-scfmDriver) shows the full list of module inputs.

<table class="table" style="color: black; margin-left: auto; margin-right: auto;">
<caption>(\#tab:moduleInputs-scfmDriver)(\#tab:moduleInputs-scfmDriver)List of (ref:scfmDriver) input objects and their description.</caption>
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
   <td style="text-align:left;"> fireRegimePolys </td>
   <td style="text-align:left;"> sf </td>
   <td style="text-align:left;"> Areas to calibrate individual fire regime parameters. Defaults to ecozones of Canada. Must have numeric field 'PolyID' or it will be created for individual polygons. </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> flammableMapCalibration </td>
   <td style="text-align:left;"> SpatRaster </td>
   <td style="text-align:left;"> a flammable map of study area after buffering by `P(sim)$buffDist`. Must be supplied by user if `flammableMap` is also supplied. </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> rasterToMatch </td>
   <td style="text-align:left;"> SpatRaster </td>
   <td style="text-align:left;"> template raster for raster GIS operations. Must be supplied by user. </td>
   <td style="text-align:left;"> NA </td>
  </tr>
</tbody>
</table>

Table \@ref(tab:moduleParams-scfmDriver) shows the full list of module parameters.

<table class="table" style="color: black; margin-left: auto; margin-right: auto;">
<caption>(\#tab:moduleParams-scfmDriver)(\#tab:moduleParams-scfmDriver)List of (ref:scfmDriver) parameters and their description.</caption>
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
   <td style="text-align:left;"> 5000 </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> 1e+05 </td>
   <td style="text-align:left;"> Buffer width for fire landscape calibration </td>
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
   <td style="text-align:left;"> used to select the year of landcover data used to create flammableMapCalibration if the object is unsupplied </td>
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
   <td style="text-align:left;"> 1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> simulation time at which the first plot event should occur </td>
  </tr>
  <tr>
   <td style="text-align:left;"> .plotInterval </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> 1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> simulation time at which the first plot event should occur </td>
  </tr>
  <tr>
   <td style="text-align:left;"> .plots </td>
   <td style="text-align:left;"> character </td>
   <td style="text-align:left;"> screen, png </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Used by Plots function, which can be optionally used here </td>
  </tr>
  <tr>
   <td style="text-align:left;"> .useCache </td>
   <td style="text-align:left;"> logical </td>
   <td style="text-align:left;"> FALSE </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Can be names of events or the whole module name; these will be cached by SpaDES </td>
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

### Module outputs

Description of the module outputs (Table \@ref(tab:moduleOutputs-scfmDriver)).

<table class="table" style="color: black; margin-left: auto; margin-right: auto;">
<caption>(\#tab:moduleOutputs-scfmDriver)(\#tab:moduleOutputs-scfmDriver)List of (ref:scfmDriver) outputs and their description.</caption>
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
   <td style="text-align:left;"> `fireRegimePolys` with driver attributes appended </td>
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
