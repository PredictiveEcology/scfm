---
title: "scfmIgnition Manual"
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
  bibliography: citations/references_scfmIgnition.bib
link-citations: true
always_allow_html: true
---

# scfmIgnition Module

<!-- the following are text references used in captions for LaTeX compatibility -->
(ref:scfmIgnition) *scfmIgnition*



[![made-with-Markdown](figures/markdownBadge.png)](https://commonmark.org)

<!-- if knitting to pdf remember to add the pandoc_args: ["--extract-media", "."] option to yml in order to get the badge images -->

#### Authors:

Steve G Cumming <stevec@sbf.ulaval.ca> [aut, cre]
<!-- ideally separate authors with new lines, '\n' not working -->

## Module Overview

### Module summary

Starts random number of fires based on estimated ignition probabilities for the landscape.

### Module inputs and parameters

Table \@ref(tab:moduleInputs-scfmIgnition) shows the full list of module inputs.

<table class="table" style="color: black; margin-left: auto; margin-right: auto;">
<caption>(\#tab:moduleInputs-scfmIgnition)(\#tab:moduleInputs-scfmIgnition)List of (ref:scfmIgnition) input objects and their description.</caption>
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
   <td style="text-align:left;"> `fireRegimePolys` with ignition rate attribute </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> fireRegimeRas </td>
   <td style="text-align:left;"> SpatRaster </td>
   <td style="text-align:left;"> rasterized version of `fireRegimePolys` </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> flammableMap </td>
   <td style="text-align:left;"> SpatRaster </td>
   <td style="text-align:left;"> map of flammability </td>
   <td style="text-align:left;"> NA </td>
  </tr>
</tbody>
</table>

Table \@ref(tab:moduleParams-scfmIgnition) shows the full list of module parameters.


<table class="table" style="color: black; margin-left: auto; margin-right: auto;">
<caption>(\#tab:moduleParams-scfmIgnition)(\#tab:moduleParams-scfmIgnition)List of (ref:scfmIgnition) parameters and their description.</caption>
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
   <td style="text-align:left;"> pIgnition </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> 0.001 </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> 1 </td>
   <td style="text-align:left;"> default per-cell and time ignition probability if unsupplied. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> startTime </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> simulation time of first ignition </td>
  </tr>
  <tr>
   <td style="text-align:left;"> returnInterval </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> 1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> interval between main events </td>
  </tr>
  <tr>
   <td style="text-align:left;"> .useCache </td>
   <td style="text-align:left;"> logical </td>
   <td style="text-align:left;"> FALSE </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Internal. Can be names of events or the whole module name; these will be cached by SpaDES </td>
  </tr>
</tbody>
</table>

### Module outputs

Description of the module outputs (Table \@ref(tab:moduleOutputs-scfmIgnition)).

<table class="table" style="color: black; margin-left: auto; margin-right: auto;">
<caption>(\#tab:moduleOutputs-scfmIgnition)(\#tab:moduleOutputs-scfmIgnition)List of (ref:scfmIgnition) outputs and their description.</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> objectName </th>
   <th style="text-align:left;"> objectClass </th>
   <th style="text-align:left;"> desc </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> ignitionLoci </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> vector of ignition locations </td>
  </tr>
  <tr>
   <td style="text-align:left;"> pIg </td>
   <td style="text-align:left;"> SpatRaster </td>
   <td style="text-align:left;"> ignition probability raster </td>
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
