---
title: "scfmSpread Manual"
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
  bibliography: citations/references_scfmSpread.bib
link-citations: true
always_allow_html: true
---

# scfmSpread Module

<!-- the following are text references used in captions for LaTeX compatibility -->
(ref:scfmSpread) *scfmSpread*



[![made-with-Markdown](figures/markdownBadge.png)](https://commonmark.org)

<!-- if knitting to pdf remember to add the pandoc_args: ["--extract-media", "."] option to yml in order to get the badge images -->

#### Authors:

Steve Cumming <stevec@sbf.ulaval.ca> [aut], Ian Eddy <ian.eddy@nrcan-rncan.gc.ca> [aut], Alex M Chubaty <achubaty@for-cast.ca> [ctb]
<!-- ideally separate authors with new lines, '\n' not working -->

## Module Overview

### Module summary

Simulates wildfire spread on a landscape.

### Module inputs and parameters

Table \@ref(tab:moduleInputs-scfmSpread) shows the full list of module inputs.

<table class="table" style="color: black; margin-left: auto; margin-right: auto;">
<caption>(\#tab:moduleInputs-scfmSpread)(\#tab:moduleInputs-scfmSpread)List of (ref:scfmSpread) input objects and their description.</caption>
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
   <td style="text-align:left;"> `fireRegimePolys` with fire attributes appended. </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> fireRegimeRas </td>
   <td style="text-align:left;"> SpatRaster </td>
   <td style="text-align:left;"> raster with fire regimes from `fireRegimePolys`. </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> flammableMap </td>
   <td style="text-align:left;"> SpatRaster </td>
   <td style="text-align:left;"> binary map of landscape flammability </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> rasterToMatch </td>
   <td style="text-align:left;"> SpatRaster </td>
   <td style="text-align:left;"> template raster for raster GIS operations. Must be supplied by user. </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> spreadState </td>
   <td style="text-align:left;"> data.table </td>
   <td style="text-align:left;"> see `SpaDES.tools::spread2` </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> studyArea </td>
   <td style="text-align:left;"> sf </td>
   <td style="text-align:left;"> Polygon to use as the simulation study area. </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> studyAreaReporting </td>
   <td style="text-align:left;"> sf </td>
   <td style="text-align:left;"> multipolygon (typically smaller/unbuffered than `studyArea`) to use for plotting/reporting. </td>
   <td style="text-align:left;"> NA </td>
  </tr>
</tbody>
</table>

Table \@ref(tab:moduleParams-scfmSpread) shows the full list of module parameters.


<table class="table" style="color: black; margin-left: auto; margin-right: auto;">
<caption>(\#tab:moduleParams-scfmSpread)(\#tab:moduleParams-scfmSpread)List of (ref:scfmSpread) parameters and their description.</caption>
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
   <td style="text-align:left;"> used to select the year of landcover data used to create `flammableMap` if the obejct is unsupplied </td>
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
   <td style="text-align:left;"> pSpread </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> 0.23 </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> 1 </td>
   <td style="text-align:left;"> default spread probability if `fireRegimePolys` is missing attribute `pSpread`. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> returnInterval </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> 1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Time interval between burn events </td>
  </tr>
  <tr>
   <td style="text-align:left;"> startTime </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Simulation time at which to initiate burning </td>
  </tr>
  <tr>
   <td style="text-align:left;"> .plotInterval </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> 10 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> This describes the simulation time at which the first plot event should occur </td>
  </tr>
  <tr>
   <td style="text-align:left;"> .plots </td>
   <td style="text-align:left;"> character </td>
   <td style="text-align:left;"> screen </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Used by Plots function, which can be optionally used here </td>
  </tr>
  <tr>
   <td style="text-align:left;"> .useCache </td>
   <td style="text-align:left;"> character </td>
   <td style="text-align:left;"> .inputOb.... </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Can be names of events or the whole module name; these will be cached by SpaDES </td>
  </tr>
</tbody>
</table>

### Module outputs

Description of the module outputs (Table \@ref(tab:moduleOutputs-scfmSpread)).

<table class="table" style="color: black; margin-left: auto; margin-right: auto;">
<caption>(\#tab:moduleOutputs-scfmSpread)(\#tab:moduleOutputs-scfmSpread)List of (ref:scfmSpread) outputs and their description.</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> objectName </th>
   <th style="text-align:left;"> objectClass </th>
   <th style="text-align:left;"> desc </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> burnDT </td>
   <td style="text-align:left;"> data.table </td>
   <td style="text-align:left;"> data table with pixel IDs of most recent burn </td>
  </tr>
  <tr>
   <td style="text-align:left;"> burnMap </td>
   <td style="text-align:left;"> SpatRaster </td>
   <td style="text-align:left;"> cumulative burn map </td>
  </tr>
  <tr>
   <td style="text-align:left;"> burnSummary </td>
   <td style="text-align:left;"> data.table </td>
   <td style="text-align:left;"> describes details of all burned pixels </td>
  </tr>
  <tr>
   <td style="text-align:left;"> pSpread </td>
   <td style="text-align:left;"> SpatRaster </td>
   <td style="text-align:left;"> spread probability applied to flammability map </td>
  </tr>
  <tr>
   <td style="text-align:left;"> rstCurrentBurn </td>
   <td style="text-align:left;"> SpatRaster </td>
   <td style="text-align:left;"> annual burn map </td>
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
