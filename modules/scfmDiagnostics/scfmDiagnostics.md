---
title: "scfmDiagnostics Manual"
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
  bibliography: citations/references_scfmDiagnostics.bib
citation-style: citations/ecology-letters.csl
link-citations: true
always_allow_html: true
---

# scfmDiagnostics Module

<!-- the following are text references used in captions for LaTeX compatibility -->
(ref:scfmDiagnostics) *scfmDiagnostics*



[![made-with-Markdown](figures/markdownBadge.png)](http://commonmark.org)

<!-- if knitting to pdf remember to add the pandoc_args: ["--extract-media", "."] option to yml in order to get the badge images -->

#### Authors:

Ian Eddy <ian.eddy@nrcan-rncan.gc.ca> [aut], Steve Cumming <stevec@sbf.ulaval.ca> [aut], Alex M Chubaty <achubaty@for-cast.ca> [aut, cre]
<!-- ideally separate authors with new lines, '\n' not working -->

## Module Overview

### Module summary

Produce diagnostic and summary plots for scfm runs.
Can be run for a single simulation, as part of the main scfm run, or as part of postprocessing.
Inputs objects will be loaded from saved simulation files when in 'multi' mode.
                      
### Module inputs and parameters

Table \@ref(tab:moduleInputs-scfmDiagnostics) shows the full list of module inputs.

<table class="table" style="color: black; margin-left: auto; margin-right: auto;">
<caption>(\#tab:moduleInputs-scfmDiagnostics)(\#tab:moduleInputs-scfmDiagnostics)List of (ref:scfmDiagnostics) input objects and their description.</caption>
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
   <td style="text-align:left;"> burnSummary </td>
   <td style="text-align:left;"> data.table </td>
   <td style="text-align:left;"> describes details of all burned pixels. Required in single mode. </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> burnMap </td>
   <td style="text-align:left;"> SpatRaster </td>
   <td style="text-align:left;"> cumulative burn map from simulation </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> fireRegimePoints </td>
   <td style="text-align:left;"> sf </td>
   <td style="text-align:left;"> Fire locations. Points outside `studyArea` are removed. Required in single mode. </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> fireRegimePolys </td>
   <td style="text-align:left;"> sf </td>
   <td style="text-align:left;"> Areas to calibrate individual fire regime parameters. Defaults to ecozones of Canada. Must have numeric field 'PolyID' or it will be created for individual polygons. Required in single mode. </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> flammableMap </td>
   <td style="text-align:left;"> SpatRaster </td>
   <td style="text-align:left;"> binary flammability map. Required in single mode. </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> studyAreaReporting </td>
   <td style="text-align:left;"> sf </td>
   <td style="text-align:left;"> multipolygon (typically smaller/unbuffered than `studyArea`) to use for plotting/reporting. Required in single mode. </td>
   <td style="text-align:left;"> NA </td>
  </tr>
</tbody>
</table>

Provide a summary of user-visible parameters (Table \@ref(tab:moduleParams-scfmDiagnostics))


<table class="table" style="color: black; margin-left: auto; margin-right: auto;">
<caption>(\#tab:moduleParams-scfmDiagnostics)(\#tab:moduleParams-scfmDiagnostics)List of (ref:scfmDiagnostics) parameters and their description.</caption>
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
   <td style="text-align:left;"> mode </td>
   <td style="text-align:left;"> character </td>
   <td style="text-align:left;"> single </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> use 'single' to run part of an scfm simulation (i.e., along with other scfm modules); use 'multi' to run as part of postprocessing multiple scfm runs. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> reps </td>
   <td style="text-align:left;"> integer </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> number of replicates/runs per study area when running in 'multi' mode. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> simOutPrefix </td>
   <td style="text-align:left;"> character </td>
   <td style="text-align:left;"> mySimOut </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> saved simList file prefix </td>
  </tr>
  <tr>
   <td style="text-align:left;"> simTimes </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> NA, NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Simulation start and end times when running in 'multi' mode. </td>
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
   <td style="text-align:left;"> .studyAreaName </td>
   <td style="text-align:left;"> character </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Human-readable name for the study area used - e.g., a hash of the studyarea obtained using `reproducible::studyAreaName()` </td>
  </tr>
  <tr>
   <td style="text-align:left;"> .useCache </td>
   <td style="text-align:left;"> logical </td>
   <td style="text-align:left;"> FALSE </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Should caching of events or module be used? </td>
  </tr>
</tbody>
</table>

### Events

A single event `diagnosticPlots` produces `ggplot`s (and saves these to disk).

### Module outputs

Description of the module outputs (Table \@ref(tab:moduleOutputs-scfmDiagnostics)).

<table class="table" style="color: black; margin-left: auto; margin-right: auto;">
<caption>(\#tab:moduleOutputs-scfmDiagnostics)(\#tab:moduleOutputs-scfmDiagnostics)List of (ref:scfmDiagnostics) outputs and their description.</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> objectName </th>
   <th style="text-align:left;"> objectClass </th>
   <th style="text-align:left;"> desc </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> scfmSummaryDT </td>
   <td style="text-align:left;"> data.table </td>
   <td style="text-align:left;"> Summary data.table containing diagnostic plot data. Can be used to create customized diagnostic plots; see `?scfmutils::comparePredictions`. </td>
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
