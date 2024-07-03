---
title: "scfmEscape Manual"
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
  bibliography: citations/references_scfmEscape.bib
link-citations: true
always_allow_html: true
---

# scfmEscape Module

<!-- the following are text references used in captions for LaTeX compatibility -->
(ref:scfmEscape) *scfmEscape*



[![made-with-Markdown](figures/markdownBadge.png)](https://commonmark.org)

<!-- if knitting to pdf remember to add the pandoc_args: ["--extract-media", "."] option to yml in order to get the badge images -->

#### Authors:

Steven G Cumming <stevec@sbf.ulaval.ca> [aut], Ian MS Eddy <ian.eddy@nrcan-rncan.gc.ca> [aut], Alex M Chubaty <achubaty@for-cast.ca> [ctb]
<!-- ideally separate authors with new lines, '\n' not working -->

## Module Overview

### Module summary

'Escapes' fire(s) from an initial set of loci returned by `scfmIgnition` and prepares the results for use by `scfmSpread`.

### Module inputs and parameters

Table \@ref(tab:moduleInputs-scfmEscape) shows the full list of module inputs.

<table class="table" style="color: black; margin-left: auto; margin-right: auto;">
<caption>(\#tab:moduleInputs-scfmEscape)(\#tab:moduleInputs-scfmEscape)List of (ref:scfmEscape) input objects and their description.</caption>
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
   <td style="text-align:left;"> fire regime polys with ignition rate </td>
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
  <tr>
   <td style="text-align:left;"> ignitionLoci </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> pixel IDs where ignition occurs </td>
   <td style="text-align:left;"> NA </td>
  </tr>
</tbody>
</table>

Table \@ref(tab:moduleParams-scfmEscape) shows the full list of module parameters.

<table class="table" style="color: black; margin-left: auto; margin-right: auto;">
<caption>(\#tab:moduleParams-scfmEscape)(\#tab:moduleParams-scfmEscape)List of (ref:scfmEscape) parameters and their description.</caption>
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
   <td style="text-align:left;"> used to select the year of landcover data used to create `flammableMap` if the obejct is unsupplied. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> neighbours </td>
   <td style="text-align:left;"> integer </td>
   <td style="text-align:left;"> 8 </td>
   <td style="text-align:left;"> 4 </td>
   <td style="text-align:left;"> 8 </td>
   <td style="text-align:left;"> Number of cell immediate neighbours (one of `4L` or `8L`). </td>
  </tr>
  <tr>
   <td style="text-align:left;"> p0 </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> 0.1 </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> 1 </td>
   <td style="text-align:left;"> probability of an ignition spreading to an unburned immediate neighbour </td>
  </tr>
  <tr>
   <td style="text-align:left;"> returnInterval </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> 1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> This specifies the time interval between Escape events </td>
  </tr>
  <tr>
   <td style="text-align:left;"> startTime </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> simulation time of first escape </td>
  </tr>
  <tr>
   <td style="text-align:left;"> .useCache </td>
   <td style="text-align:left;"> character </td>
   <td style="text-align:left;"> .inputOb.... </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Internal. Can be names of events or the whole module name; these will be cached by SpaDES. </td>
  </tr>
</tbody>
</table>

### Module outputs

Description of the module outputs (Table \@ref(tab:moduleOutputs-scfmEscape)).

<table class="table" style="color: black; margin-left: auto; margin-right: auto;">
<caption>(\#tab:moduleOutputs-scfmEscape)(\#tab:moduleOutputs-scfmEscape)List of (ref:scfmEscape) outputs and their description.</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> objectName </th>
   <th style="text-align:left;"> objectClass </th>
   <th style="text-align:left;"> desc </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> spreadState </td>
   <td style="text-align:left;"> data.table </td>
   <td style="text-align:left;"> stores the current fire spread state </td>
  </tr>
  <tr>
   <td style="text-align:left;"> p0 </td>
   <td style="text-align:left;"> SpatRaster </td>
   <td style="text-align:left;"> escape probability raster </td>
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
