"""
design.py

One design system for every baitUtils figure and report: colour tokens,
typography, the matplotlib style, the plotly template and the CSS variables
of the HTML reports all come from the dictionaries in this module, so the
static figures, the interactive figures and the report chrome look like one
document.

Colours follow a validated, colour-blind-safe palette: eight categorical
hues in a fixed order (series identity, never cycled past eight), one
sequential blue ramp (magnitude), a blue-red diverging pair and four
status colours that are reserved for the quality categories. Text never
wears a data colour; identity comes from a swatch beside the text.
"""

from __future__ import annotations

from typing import Dict, Iterable, List, Sequence, Tuple

import math
import numpy as np

# ------------------------------------------------------------------ tokens

#: Categorical hues in their validated order: blue, orange, aqua, yellow,
#: magenta, green, violet, red. Used for sets, references and any "which
#: series" encoding. Adjacent pairs clear the colour-vision-deficiency gate;
#: all-pairs forms (scatter) should stay within the first three.
CATEGORICAL: List[str] = [
    "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948",
]

#: Sequential blue ramp, light to dark (steps 100 to 700).
SEQUENTIAL: List[str] = [
    "#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec", "#5598e7", "#3987e5",
    "#2a78d6", "#256abf", "#1c5cab", "#184f95", "#104281", "#0d366b",
]

#: Diverging pair with a neutral midpoint (blue, grey, red).
DIVERGING: Tuple[str, str, str] = ("#2a78d6", "#f0efec", "#e34948")

#: Status colours, reserved for states; always shown with a label.
STATUS: Dict[str, str] = {
    "good": "#0ca30c",
    "warning": "#fab219",
    "serious": "#ec835a",
    "critical": "#d03b3b",
}

#: Quality categories map onto the status scale (Good keeps the accent
#: blue, which reads as "fine" without claiming excellence).
CATEGORY_COLORS: Dict[str, str] = {
    "Excellent": STATUS["good"],
    "Good": "#2a78d6",
    "Fair": STATUS["warning"],
    "Poor": STATUS["critical"],
}

#: Chrome and ink.
CHROME: Dict[str, str] = {
    "surface": "#fcfcfb",       # figure and card surface
    "page": "#f9f9f7",          # page plane behind the cards
    "ink": "#0b0b0b",           # primary text
    "ink_secondary": "#52514e",
    "ink_muted": "#898781",     # axis labels, captions
    "grid": "#e1e0d9",          # hairline gridlines
    "axis": "#c3c2b7",          # baseline and axis rules
    "border": "rgba(11,11,11,0.10)",
    "accent": "#2a78d6",
    "emphasis_gray": "#c3c2b7",  # de-emphasised series
    "gap_fill": "#e34948",      # uncovered regions (status-like, always labelled)
}

FONT_STACK = 'system-ui, -apple-system, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif'
MPL_FONTS = ["Inter", "Helvetica Neue", "Helvetica", "Arial", "DejaVu Sans", "sans-serif"]


# ----------------------------------------------------------------- helpers

def series_color(index: int) -> str:
    """Categorical colour for the series at a position; past eight, the de-emphasis grey."""
    return CATEGORICAL[index] if 0 <= index < len(CATEGORICAL) else CHROME["emphasis_gray"]


def category_color(category: str) -> str:
    """Colour of a quality category (Excellent, Good, Fair, Poor)."""
    return CATEGORY_COLORS.get(str(category), CHROME["ink_muted"])


def label_in_fill(hex_color: str) -> str:
    """Ink for a label drawn inside a filled mark: white, or near-black on light fills."""
    rgb = hex_to_rgb(hex_color)
    luminance = 0.2126 * rgb[0] + 0.7152 * rgb[1] + 0.0722 * rgb[2]
    return CHROME["ink"] if luminance > 0.55 else "#ffffff"


def hex_to_rgb(hex_color: str) -> Tuple[float, float, float]:
    h = hex_color.lstrip("#")
    return tuple(int(h[i:i + 2], 16) / 255.0 for i in (0, 2, 4))


def rgba(hex_color: str, alpha: float) -> str:
    r, g, b = (int(round(c * 255)) for c in hex_to_rgb(hex_color))
    return f"rgba({r},{g},{b},{alpha})"


def fmt_bp(value: float) -> str:
    """Compact base-pair label: 850 bp, 12.5 kb, 3.2 Mb, 1.1 Gb."""
    value = float(value)
    for unit, factor in (("Gb", 1e9), ("Mb", 1e6), ("kb", 1e3)):
        if abs(value) >= factor:
            number = value / factor
            return f"{number:.1f} {unit}" if number < 100 else f"{number:.0f} {unit}"
    return f"{value:.0f} bp"


def fmt_count(value: float) -> str:
    """Thousands separators for counts."""
    return f"{int(round(value)):,}"


def fmt_pct(value: float, digits: int = 1) -> str:
    return f"{value:.{digits}f}%"


def fmt_depth(value: float) -> str:
    """Depth with a multiplication sign, one decimal below 10 and none above."""
    return f"{value:.1f}x" if value < 10 else f"{value:.0f}x"


def truncate(label: str, width: int = 24) -> str:
    """Shorten a long sequence id for an axis without losing its start."""
    label = str(label)
    return label if len(label) <= width else label[:width - 1] + "…"


def bin_coverage(array: Sequence[float], n_bins: int = 600) -> Tuple[np.ndarray, np.ndarray]:
    """
    Mean depth per bin along a reference for a coverage profile. Returns the
    bin centres (bp) and the mean depths; references shorter than n_bins
    are returned base by base.
    """
    values = np.asarray(array, dtype=float)
    length = values.size
    if length == 0:
        return np.array([]), np.array([])
    if length <= n_bins:
        return np.arange(length, dtype=float) + 0.5, values
    edges = np.linspace(0, length, n_bins + 1).astype(int)
    sums = np.add.reduceat(values, edges[:-1])
    widths = np.diff(edges)
    centres = (edges[:-1] + edges[1:]) / 2.0
    return centres, sums / np.maximum(widths, 1)


def select_references(per_reference: Dict[str, Dict], limit: int = 8) -> List[str]:
    """Reference ids ordered by length, longest first, capped at limit."""
    ordered = sorted(per_reference, key=lambda r: per_reference[r].get("length", 0), reverse=True)
    return ordered[:limit]


# ------------------------------------------------------------- matplotlib

def apply_matplotlib_style(dpi: int = 300) -> None:
    """Set the matplotlib rcParams of the design system (idempotent)."""
    import matplotlib as mpl
    from cycler import cycler

    mpl.rcParams.update({
        "figure.facecolor": CHROME["surface"],
        "figure.dpi": 100,
        "savefig.dpi": dpi,
        "savefig.facecolor": CHROME["surface"],
        "savefig.edgecolor": "none",
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.15,
        "axes.facecolor": CHROME["surface"],
        "axes.edgecolor": CHROME["axis"],
        "axes.linewidth": 0.8,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "axes.grid.axis": "y",
        "axes.axisbelow": True,
        "grid.color": CHROME["grid"],
        "grid.linewidth": 0.8,
        "grid.linestyle": "-",
        "axes.titlesize": 11,
        "axes.titleweight": "semibold",
        "axes.titlelocation": "left",
        "axes.titlepad": 10,
        "axes.titlecolor": CHROME["ink"],
        "axes.labelsize": 9.5,
        "axes.labelcolor": CHROME["ink_secondary"],
        "axes.prop_cycle": cycler(color=CATEGORICAL),
        "xtick.labelsize": 8.5,
        "ytick.labelsize": 8.5,
        "xtick.color": CHROME["ink_muted"],
        "ytick.color": CHROME["ink_muted"],
        "xtick.major.size": 3,
        "ytick.major.size": 0,
        "xtick.major.width": 0.8,
        "legend.fontsize": 8.5,
        "legend.frameon": False,
        "legend.labelcolor": CHROME["ink_secondary"],
        "font.family": "sans-serif",
        "font.sans-serif": MPL_FONTS,
        "font.size": 9.5,
        "text.color": CHROME["ink"],
        "figure.titlesize": 12,
        "figure.titleweight": "semibold",
        "lines.linewidth": 2.0,
        "lines.markersize": 5,
        "patch.linewidth": 0,
    })


def style_axes(ax, grid: str = "y", baseline: bool = True) -> None:
    """Recessive axes: no top/right spines, hairline grid on one axis, light rules."""
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(CHROME["axis"])
        ax.spines[side].set_linewidth(0.8)
    ax.grid(False)
    if grid in ("x", "y", "both"):
        ax.grid(True, axis=grid, color=CHROME["grid"], linewidth=0.8, linestyle="-")
    ax.set_axisbelow(True)
    if not baseline:
        ax.spines["left"].set_visible(False)
    ax.tick_params(colors=CHROME["ink_muted"], labelsize=8.5)
    ax.title.set_color(CHROME["ink"])


def empty_panel(ax, message: str) -> None:
    """A panel with nothing to plot says so instead of showing empty axes."""
    ax.set_axis_off()
    ax.text(0.5, 0.5, message, ha="center", va="center", transform=ax.transAxes,
            color=CHROME["ink_muted"], fontsize=9.5)


# ----------------------------------------------------------------- plotly

def plotly_template():
    """The plotly template matching the matplotlib style."""
    import plotly.graph_objects as go

    axis = dict(
        showgrid=True, gridcolor=CHROME["grid"], gridwidth=1, zeroline=False,
        linecolor=CHROME["axis"], linewidth=1, showline=True, ticks="outside",
        tickcolor=CHROME["axis"], ticklen=3, tickfont=dict(color=CHROME["ink_muted"], size=11),
        title=dict(font=dict(color=CHROME["ink_secondary"], size=12)),
    )
    return go.layout.Template(
        layout=go.Layout(
            font=dict(family=FONT_STACK, color=CHROME["ink"], size=12),
            paper_bgcolor=CHROME["surface"],
            plot_bgcolor=CHROME["surface"],
            colorway=CATEGORICAL,
            title=dict(font=dict(size=15, color=CHROME["ink"]), x=0, xanchor="left", y=0.985, yanchor="top"),
            xaxis=axis,
            yaxis=axis,
            legend=dict(bgcolor="rgba(0,0,0,0)", font=dict(color=CHROME["ink_secondary"], size=11),
                        orientation="h", yanchor="bottom", y=1.0, xanchor="right", x=1),
            margin=dict(l=60, r=24, t=64, b=56),
            hoverlabel=dict(bgcolor="#ffffff", bordercolor=CHROME["axis"],
                            font=dict(family=FONT_STACK, color=CHROME["ink"], size=12)),
            colorscale=dict(sequential=[[i / (len(SEQUENTIAL) - 1), c] for i, c in enumerate(SEQUENTIAL)],
                            diverging=[[0, DIVERGING[0]], [0.5, DIVERGING[1]], [1, DIVERGING[2]]]),
        )
    )


def plotly_layout(fig, title: str, height: int, legend_above_subplots: bool = False, **kwargs) -> None:
    """Apply the template, title and height.

    Subplot titles are restyled to the design system. With legend_above_subplots the top margin
    makes room for three bands: the title and mode bar, the legend, and the subplot titles.
    """
    fig.update_layout(template=plotly_template(), title=dict(text=title), height=height, **kwargs)
    for annotation in fig.layout.annotations:
        if annotation.xref == "paper" and annotation.yref == "paper" and annotation.font.size == 16:
            annotation.font.update(size=13, color=CHROME["ink"])
    if legend_above_subplots and "legend" not in kwargs:
        margin = fig.layout.margin
        top = max(margin.t if margin.t is not None else 64, 84)
        bottom = margin.b if margin.b is not None else 56
        plot_height = max(height - top - bottom, 1)
        fig.update_layout(margin=dict(t=top), legend=dict(y=1 + 30 / plot_height))


def subplot_note(fig, text: str, row: int, col: int) -> None:
    """Write a note at the centre of a subplot that has nothing to plot.

    Paper coordinates are used because an axis without traces is not drawn, and an
    annotation bound to it would fall back to the first subplot.
    """
    subplot = fig.get_subplot(row, col)
    x0, x1 = subplot.xaxis.domain
    y0, y1 = subplot.yaxis.domain
    fig.add_annotation(text=text, showarrow=False, xref="paper", yref="paper",
                       x=(x0 + x1) / 2, y=(y0 + y1) / 2, font=dict(color=CHROME["ink_muted"], size=12))


def histogram_bins(values, n_bins: int = 20) -> Dict[str, float]:
    """Plotly xbins with a round bin width, also when every value is the same."""
    values = np.asarray(list(values), dtype=float)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return {"start": 0, "end": 1, "size": 1}
    low, high = float(values.min()), float(values.max())
    span = high - low if high > low else max(abs(high), 1.0)
    raw = span / max(n_bins, 1)
    magnitude = 10 ** math.floor(math.log10(raw)) if raw > 0 else 1
    size = next(m * magnitude for m in (1, 2, 5, 10) if m * magnitude >= raw)
    start = math.floor(low / size) * size
    end = (math.floor(high / size) + 1) * size
    return {"start": start, "end": end, "size": size}


PLOTLY_CONFIG = {
    "displayModeBar": True,
    "modeBarButtonsToRemove": ["lasso2d", "select2d", "autoScale2d"],
    "displaylogo": False,
    "toImageButtonOptions": {"format": "svg", "filename": "baitUtils_figure", "scale": 1},
}


def plotly_sequential_scale() -> List[List]:
    return [[i / (len(SEQUENTIAL) - 1), c] for i, c in enumerate(SEQUENTIAL)]


# -------------------------------------------------------------------- CSS

def css_variables() -> str:
    """CSS custom properties of the design system for the HTML reports."""
    lines = [
        f"  --surface: {CHROME['surface']};",
        f"  --page: {CHROME['page']};",
        f"  --ink: {CHROME['ink']};",
        f"  --ink-2: {CHROME['ink_secondary']};",
        f"  --ink-3: {CHROME['ink_muted']};",
        f"  --grid: {CHROME['grid']};",
        f"  --rule: {CHROME['axis']};",
        f"  --border: {CHROME['border']};",
        f"  --accent: {CHROME['accent']};",
        f"  --accent-soft: {SEQUENTIAL[0]};",
        f"  --font: {FONT_STACK};",
    ]
    for name, color in STATUS.items():
        lines.append(f"  --status-{name}: {color};")
    for name, color in CATEGORY_COLORS.items():
        lines.append(f"  --category-{name.lower()}: {color};")
    for i, color in enumerate(CATEGORICAL, 1):
        lines.append(f"  --series-{i}: {color};")
    return ":root {\n" + "\n".join(lines) + "\n}\n"


def swatch(color: str) -> str:
    """Inline HTML swatch that carries a series or category colour beside text."""
    return (f'<span class="swatch" style="background:{color}" aria-hidden="true"></span>')


def status_badge(category: str) -> str:
    """Category label with its colour dot; the text always carries the meaning."""
    return (f'<span class="badge"><span class="swatch" style="background:{category_color(category)}" '
            f'aria-hidden="true"></span>{category}</span>')


def ordered_unique(values: Iterable[str]) -> List[str]:
    seen: List[str] = []
    for value in values:
        if value not in seen:
            seen.append(value)
    return seen
