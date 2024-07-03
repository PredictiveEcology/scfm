---
title: "group_scfm Manual"
subtitle: "v."
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
  bibliography: citations/references_group_scfm.bib
citation-style: citations/ecology-letters.csl
link-citations: true
always_allow_html: true
---

# group_scfm Module

<!-- the following are text references used in captions for LaTeX compatibility -->
(ref:group_scfm) *group_scfm*



[![made-with-Markdown](figures/markdownBadge.png)](http://commonmark.org)

<!-- if knitting to pdf remember to add the pandoc_args: ["--extract-media", "."] option to yml in order to get the badge images -->

#### Authors:


<!-- ideally separate authors with new lines, '\n' not working -->

## Module Overview

### Module summary

A module group to run the `scfm` suite of modules found at <https://github.com/PredictiveEcology/scfm>:

- `ageModule` (optional)
- `scfmDriver`
- `scfmEscape`
- `scfmIgnition`
- `scfmLandcoverInit`
- `scfmRegime`
- `scfmSpread`

### Module parameters, inputs and outputs

See documentation for individual modules.

### Getting help

<https://github.com/PredictiveEcology/scfm/issues>
