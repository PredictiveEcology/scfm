---
title: "ageModule Manual"
subtitle: "v.2.0.0"
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
  bibliography: citations/references_ageModule.bib
link-citations: true
always_allow_html: true
---

# ageModule Module

<!-- the following are text references used in captions for LaTeX compatibility -->
(ref:ageModule) *ageModule*



[![made-with-Markdown](figures/markdownBadge.png)](https://commonmark.org)

<!-- if knitting to pdf remember to add the pandoc_args: ["--extract-media", "."] option to yml in order to get the badge images -->

#### Authors:

Steve G Cumming <stevec@sbf.ulaval.ca> [aut, cre]
<!-- ideally separate authors with new lines, '\n' not working -->

## Module Overview

### Module summary

Creates and maintains a raster called `ageMap`.
This module is optional when running the `scfm` suite of modules.

### Module inputs and parameters

Table \@ref(tab:moduleInputs-ageModule) shows the full list of module inputs.

<table class="table" style="color: black; margin-left: auto; margin-right: auto;">
<caption>(\#tab:moduleInputs-ageModule)(\#tab:moduleInputs-ageModule)List of (ref:ageModule) input objects and their description.</caption>
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
   <td style="text-align:left;"> ageMap </td>
   <td style="text-align:left;"> SpatRaster </td>
   <td style="text-align:left;"> stand age map in study area, default is Canada national stand age map </td>
   <td style="text-align:left;"> http://tree.pfc.forestry.ca/kNN-StructureStandVolume.tar </td>
  </tr>
  <tr>
   <td style="text-align:left;"> studyArea </td>
   <td style="text-align:left;"> sf </td>
   <td style="text-align:left;"> Polygon to use as the simulation study area. </td>
   <td style="text-align:left;"> http://sis.agr.gc.ca/cansis/nsdb/ecostrat/district/ecodistrict_shp.zip </td>
  </tr>
  <tr>
   <td style="text-align:left;"> rasterToMatch </td>
   <td style="text-align:left;"> SpatRaster </td>
   <td style="text-align:left;"> template raster for raster GIS operations. Must be supplied by user. </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> rstCurrentBurn </td>
   <td style="text-align:left;"> SpatRaster </td>
   <td style="text-align:left;"> annual burn map created by `scfmSpread`. </td>
   <td style="text-align:left;"> NA </td>
  </tr>
</tbody>
</table>

Summary of user-visible parameters (Table \@ref(tab:moduleParams-ageModule)).


<table class="table" style="color: black; margin-left: auto; margin-right: auto;">
<caption>(\#tab:moduleParams-ageModule)(\#tab:moduleParams-ageModule)List of (ref:ageModule) parameters and their description.</caption>
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
   <td style="text-align:left;"> initialAge </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> 99 </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> 10000 </td>
   <td style="text-align:left;"> initial age </td>
  </tr>
  <tr>
   <td style="text-align:left;"> maxAge </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> 200 </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> 65535 </td>
   <td style="text-align:left;"> maximum age for plotting </td>
  </tr>
  <tr>
   <td style="text-align:left;"> returnInterval </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> 1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Time interval between aging events </td>
  </tr>
  <tr>
   <td style="text-align:left;"> startTime </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Simulation time at which to initiate aging </td>
  </tr>
  <tr>
   <td style="text-align:left;"> .plotInitialTime </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> This describes the simulation time at which the first plot event should occur </td>
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
   <td style="text-align:left;"> screen, png </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Used by `Plots()`, which can be optionally used here </td>
  </tr>
</tbody>
</table>

### Module outputs

Description of the module outputs (Table \@ref(tab:moduleOutputs-ageModule)).

<table class="table" style="color: black; margin-left: auto; margin-right: auto;">
<caption>(\#tab:moduleOutputs-ageModule)(\#tab:moduleOutputs-ageModule)List of (ref:ageModule) outputs and their description.</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> objectName </th>
   <th style="text-align:left;"> objectClass </th>
   <th style="text-align:left;"> desc </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> ageMap </td>
   <td style="text-align:left;"> SpatRaster </td>
   <td style="text-align:left;"> map of vegetation age </td>
  </tr>
</tbody>
</table>

### Links to other modules

Intended to be run with the `scfm` suite of modules found at <https://github.com/PredictiveEcology/scfm>:

- `ageModule` (optional)
- `scfmDataPrep`
- `scfmEscape`
- `scfmIgnition`
- `scfmSpread`

### Getting help

<https://github.com/PredictiveEcology/scfm/issues>
