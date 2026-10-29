# Report and figure design audit (October 2026)

Scope: the HTML reports written by `evaluate` and `compare`, the static
figures in `plots/` and `comparative_plots/`, and the interactive figures in
`interactive_plots/` and `interactive_comparative_plots/`. The audit was done
on the fixture dataset from `tests/conftest.py` (two references, 43 baits,
one 300 bp hole) with a second bait set for `compare`, rendered with headless
Chromium.

## Findings

Findings are grouped by what they cost the reader. D-numbers are referenced
from the CHANGELOG.

### Wrong or misleading chart forms

- D-1 A donut chart showed one proportion (covered against uncovered). A
  single proportion needs a number or one bar, not a ring that asks the eye
  to compare angles.
- D-2 Radar charts carried the quality components and the compare metrics.
  Radar area depends on the arbitrary order of the axes, the axes shared no
  unit, and the compare radar scaled every metric from 0 to its maximum so
  that a 0.025 Gini coefficient drew as a full spoke.
- D-3 The compare "performance heatmap" scaled each column to the 0-1 range
  of the sets, so with two sets every cell was 0 or 1 and the colour said
  only "higher" or "lower"; the raw values were not shown.
- D-4 Effect sizes were plotted for tests with n = 2 references without any
  caveat, and p-values were drawn as bars (a bar implies a quantity that
  starts at zero; p-values are compared with alpha).
- D-5 The depth distribution used fixed thresholds (1x, 5x, 10x, 20x) that
  ignored the user's `--min-coverage` and `--target-coverage`.
- D-6 Count axes showed fractional ticks (0.5 gaps); several histograms
  had degenerate bins when the data had one value.
- D-7 Empty panels were drawn as blank white boxes with no explanation
  (no gap features, no challenging windows, correlations not applicable).
- D-8 Two colorbars overlapped in the reference analysis figure; the
  static coverage heatmap had its own scale per panel.
- D-9 The reports had no coverage profile along the references, the one
  figure a bait designer looks at first, and no per-reference table.

### Colour

- D-10 Three accents competed in one page: a purple gradient header, a blue
  Plotly default colourway and a seaborn palette in the PNGs. Colours did not
  match between the static and the interactive version of the same figure.
- D-11 Each panel picked its own palette (Viridis, Set2, Blues, tab10),
  so the same reference or bait set changed colour between panels.
- D-12 Status colours (red, green) were used for ordinary series, and
  category labels (Excellent, Poor) were conveyed by colour alone.
- D-13 Gradient dashboard cards with white text failed contrast on the
  lighter half of the gradient.

### Text, typography and layout

- D-14 A `\002B` literal was shown instead of a plus sign.
- D-15 Figures had no numbers or captions; tables had no captions; the
  reader could not cite "Figure 3" or know what a value meant without the
  surrounding prose.
- D-16 The compare report mixed two stylesheets and two scripts with the
  evaluate report, with different fonts, widths and section styles.
- D-17 Section numbering, the table of contents and the footer differed
  between the two reports; the footer carried a product credit instead of
  the tool versions and the command that produced the report.
- D-18 Numbers were formatted inconsistently (1234567 bp next to 1.2 Mb,
  0.8375 next to 0.84) and long reference names overflowed table cells.
- D-19 Plot div ids collided with section ids (`gaps`, `quality`,
  `reference`), so Plotly drew the figure into the section element instead
  of the figure box, displacing it above its heading.
- D-20 Legends overlapped the Plotly mode bar and the subplot titles; a
  legend was drawn for single-series figures.

## The design system

All colours, fonts, formatting and layout rules live in `baitUtils/design.py`
and are used by every figure and both reports.

Colour has one job per palette:

- Categorical (identity: references, bait sets, gap size classes), fixed
  order, never cycled: `#2a78d6`, `#eb6834`, `#1baf7a`, `#eda100`,
  `#e87ba4`, `#008300`, `#6250d6`, `#e34948`. A ninth series folds into a
  grey "other". The sequence passes the adjacent-pair colour-vision checks
  on a white surface; the yellow, pink and green entries need a visible
  label beside the mark (they always have one: a legend or a direct label).
- Sequential (magnitude: depth, coverage maps, the comparison matrix): one
  blue ramp, light to dark, white for zero.
- Status (state only): Excellent green, Good blue, Fair amber, Poor red,
  always with the category name beside the swatch.
- Gap fill: one reserved red, used for gaps and uncovered bases and for
  nothing else.

Marks and chrome: 2 px lines, markers of at least 8 px with a surface ring,
bars with a 2 px surface gap, hairline grids in `#e1e0d9`, axis lines in
`#c3c2b7`, secondary ink `#52514e` for labels, primary ink `#0b0b0b` for
values. Text never takes a series colour. Every figure with two or more
series has a legend; single-series figures have none. Panels with nothing to
plot say why.

Report page: a 220 px sticky table of contents, a 960 px column, numbered
sections, numbered figures with captions, numbered tables with captions,
KPI tiles for the headline numbers, a metadata line (inputs, mapper,
thresholds, time), a footer with the versions of baitUtils, Python, numpy,
scipy and Biopython and the command line. The same CSS (`templates/report.css`)
and script (`templates/report.js`) serve both reports; the CSS variables are
emitted from `design.css_variables()` so that the page and the figures share
the tokens.

Plotly figures use `design.plotly_template()` and `design.PLOTLY_CONFIG`
(SVG export, no logo). Matplotlib figures use `design.apply_matplotlib_style`
and `design.style_axes`.

## What changed

### evaluate

Static figures (`plots/`):

- `coverage_profile`: new. Depth along the eight longest references, gaps
  shaded, the minimum depth as a dotted line.
- `coverage_overview`: composition bars (covered, uncovered, N) per
  reference with the overall bar on top, and a cumulative depth curve with
  the minimum and target depth marked. Replaces the donut.
- `depth_distribution`: histogram with integer bins, and median with the
  10th to 90th percentile per reference.
- `gap_analysis`: gap length histogram with round bins (log bins when the
  range spans more than a decade) and the largest gaps.
- `per_reference_coverage`: breadth per reference with the overall breadth
  marked, and breadth against length.
- `coverage_heatmap`: a coverage map (reference by relative position) on
  the sequential scale, instead of a matrix of scaled summary statistics.
- `quality_score`: new. Component scores with their weights and the overall
  score against the category thresholds. Replaces the radar.

Interactive figures (`interactive_plots/`): `coverage_dashboard.html` is
the coverage profile per reference; `detailed_gap_analysis.html`,
`quality_assessment.html` and `reference_analysis.html` keep their names
with the forms above. The empty-panel notes are placed in paper
coordinates, because an axis without traces is not drawn.

Report: Summary tiles, Coverage (profile and cumulative depth, embedded),
Per-reference results (new table: length, GC, N, breadth, breadth over
non-N bases, mean and maximum depth, gaps, hits per strand, difficulty),
Gaps, Quality score, Reference sequence analysis (correlations, breadth by
GC decile, challenging windows), Recommendations, Detailed statistics, Data
files, Methods.

### compare

Static figures (`comparative_plots/`): `coverage_distributions` (breadth,
depth, uniformity and mapping efficiency per set as bars in the set's
colour), `gap_analysis_comparison`, `performance_heatmap` (raw values,
shaded from worst to best within a column, light text on dark cells),
`statistical_comparison` (adjusted p-values as lollipops against alpha,
signed effect sizes, window depth boxplots, paired breadth per reference).

Interactive figures (`interactive_comparative_plots/`):
`comparison_dashboard.html` is the coverage profile per reference with one
line per set; `metric_comparison.html` (replaces `performance_radar.html`)
shows every metric as a fraction of the best set.

Report: the bait sets with their colours, Summary tiles and ranking,
Comparison table (best value per column in bold), Ranking with inline bars,
Coverage profiles, Statistical tests as tables with the significance legend
and a small-n caveat below five references, Gap overlap, Figures and data
files, Recommendations, Pairwise differences, Methods.

### Verification

The fixture reports were rendered with headless Chromium at 1280 px and
inspected page by page for collisions, clipped labels, empty panels and
misplaced figures. `tests/test_reports.py` checks the figure ids, the file
sets and the section text; the rest of the suite covers the report
generators with mocked analyses.
