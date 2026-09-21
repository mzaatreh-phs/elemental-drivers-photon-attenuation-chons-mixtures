"""Matplotlib helpers that imitate the OriginLab figures in ~/CHON/CHOM_submitted version/figs/ (fig2_mac_lac, fig3_mfp_hvl_tvl, fig4_rpe, fig5_ebf, Z_eff):
boxed axes, outward ticks on the bottom/left only, no grid/title, Arial-like font (Liberation Sans), boxed legend with line+marker,
Origin colour/marker order for MIX 1..10, red zoom rectangle with black connector lines to an inset."""
import os, csv
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, ConnectionPatch
from matplotlib.ticker import FixedLocator, FuncFormatter, NullFormatter, LogLocator, AutoMinorLocator
from decimal import Decimal, ROUND_HALF_UP
from matplotlib.lines import Line2D

FD = os.path.expanduser('~/CHON/CHOM_submitted version/figure_data/')
OUT = os.path.expanduser('~/CHON/CHOM_submitted version/figs_origin_style/')
os.makedirs(OUT, exist_ok=True)

# colours sampled from the legend of the Origin figure fig3_mfp_hvl_tvl.png; marker order as in Origin
COL = ['#515151', '#f14040', '#1a6fdf', '#37ad6b', '#b177de', '#cc9900', '#00cbcc', '#7d4e4e', '#8e8e00', '#fb6501']
MK = ['s', 'o', '^', 'v', 'D', '<', '>', 'h', '*', 'p']
MSF = {'D': 0.8, '*': 1.5, 'p': 1.15, 'h': 1.1}        # marker-size factors so that all marker types look equally big
LAB = [f'MIX {i}' for i in range(1, 11)]
FS, FSL, FSLEG = 14, 17, 10                             # tick labels, axis labels, legend

mpl.rcParams.update({
    'font.family': 'Liberation Sans', 'mathtext.fontset': 'custom', 'mathtext.rm': 'Liberation Sans', 'mathtext.it': 'Liberation Sans:italic',
    'mathtext.bf': 'Liberation Sans:bold', 'mathtext.default': 'regular', 'axes.unicode_minus': False,
    'axes.linewidth': 1.4, 'axes.grid': False, 'axes.labelsize': FSL, 'xtick.labelsize': FS, 'ytick.labelsize': FS,
    'xtick.direction': 'out', 'ytick.direction': 'out', 'xtick.top': False, 'ytick.right': False,
    'xtick.major.size': 6.5, 'ytick.major.size': 6.5, 'xtick.minor.size': 3.5, 'ytick.minor.size': 3.5,
    'xtick.major.width': 1.3, 'ytick.major.width': 1.3, 'xtick.minor.width': 1.1, 'ytick.minor.width': 1.1,
    'xtick.minor.visible': True, 'ytick.minor.visible': True,
    'legend.frameon': True, 'legend.edgecolor': 'black', 'legend.fancybox': False, 'legend.framealpha': 1.0, 'legend.fontsize': FSLEG,
    'legend.handlelength': 1.6, 'legend.handletextpad': 0.5, 'legend.borderpad': 0.25, 'legend.labelspacing': 0.15, 'legend.numpoints': 1,
    'savefig.dpi': 300, 'figure.dpi': 100, 'savefig.facecolor': 'white',
})


# ------------------------------------------------------------------ data
def rd(path):
    """csv in figure_data (relative path) -> (header list, float array; blanks -> nan). Non-numeric columns come back as nan."""
    r = list(csv.reader(open(FD + path)))
    def f(x):
        try: return float(x) if x != '' else np.nan
        except ValueError: return np.nan
    return r[0], np.array([[f(x) for x in row] for row in r[1:]])


def rdtxt(path):
    r = list(csv.reader(open(FD + path))); return r[0], r[1:]


# ------------------------------------------------------------------ axes
def frame(ax):
    for s in ax.spines.values(): s.set_linewidth(1.4)
    ax.tick_params(which='both', top=False, right=False)
    if ax.get_xscale() == 'linear': ax.xaxis.set_minor_locator(AutoMinorLocator(2))     # one minor tick between majors, as in Origin
    if ax.get_yscale() == 'linear': ax.yaxis.set_minor_locator(AutoMinorLocator(2))


def log10_x(ax, xlim, ticks=(0.01, 0.1, 1, 10)):
    ax.set_xscale('log'); ax.set_xlim(*xlim)
    ax.xaxis.set_major_locator(FixedLocator(list(ticks))); ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f'{v:g}'))
    ax.xaxis.set_minor_locator(LogLocator(base=10, subs=np.arange(2, 10), numticks=20)); ax.xaxis.set_minor_formatter(NullFormatter())


def log2_x(ax, xlim, k0, k1, step):
    """Origin 'Log2' axis: major ticks at 2**k0 ... 2**k1 (every `step` octaves), labels with two decimals as Origin prints them."""
    ax.set_xscale('log', base=2); ax.set_xlim(*xlim)
    ax.xaxis.set_major_locator(FixedLocator([2.0 ** k for k in range(k0, k1 + 1, step)]))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: str(Decimal(float(v)).quantize(Decimal('0.01'), rounding=ROUND_HALF_UP))))   # Origin rounds 0.125 -> 0.13
    ax.xaxis.set_minor_locator(FixedLocator([2.0 ** k for k in range(k0 - 2, k1 + 3)])); ax.xaxis.set_minor_formatter(NullFormatter())


def pad_ylim(ax, y, lo_frac=0.04, hi_frac=0.05, zero_floor=True):
    y = np.asarray(y, float); y = y[np.isfinite(y)]; top = y.max(); bot = 0.0 if (zero_floor and y.min() >= 0) else y.min()
    span = top - bot; ax.set_ylim((bot - lo_frac * span) if zero_floor else (bot - lo_frac * span), top + hi_frac * span)


# ------------------------------------------------------------------ series
def plot_mixes(ax, x, Y, ms=5.5, lw=1.0, yerr=None, cols=range(10)):
    """Y: (n_energy, 10). One line + marker per mixture in the Origin colours; optional error bars (small caps as in Origin)."""
    for i in cols:
        m = MSF.get(MK[i], 1.0) * ms
        if yerr is not None:
            ax.errorbar(x, Y[:, i], yerr=yerr[:, i], fmt='none', ecolor=COL[i], elinewidth=0.8, capsize=2.2, capthick=0.8, zorder=2)
        ax.plot(x, Y[:, i], color=COL[i], lw=lw, marker=MK[i], ms=m, mec=COL[i], mfc=COL[i], zorder=3, label=LAB[i])


def mix_handles(ms=4.8):
    return [Line2D([0], [0], color=COL[i], lw=1.0, marker=MK[i], ms=MSF.get(MK[i], 1.0) * ms, mec=COL[i], mfc=COL[i]) for i in range(10)]


def legend(ax, loc='upper right', **kw):
    kw.setdefault('borderaxespad', 0.5)
    return ax.legend(handles=mix_handles(), labels=LAB, loc=loc, **kw)


def fig_legend(fig, loc, anchor, **kw):
    return fig.legend(handles=mix_handles(), labels=LAB, loc=loc, bbox_to_anchor=anchor, **kw)


# ------------------------------------------------------------------ zoom inset (red rectangle + two black connector lines, as in the Origin figures)
def zoom_inset(ax, xlim, ylim, rect, draw, xfix=None, corners=('tl', 'br')):
    axi = ax.inset_axes(rect); draw(axi); axi.set_xlim(*xlim); axi.set_ylim(*ylim); frame(axi)
    axi.tick_params(labelsize=FS)
    if xfix: xfix(axi)
    ax.add_patch(Rectangle((xlim[0], ylim[0]), xlim[1] - xlim[0], ylim[1] - ylim[0], fill=False, ec='#ff6b6b', lw=1.3, zorder=6))
    pts = {'tl': ((xlim[0], ylim[1]), (0, 1)), 'br': ((xlim[1], ylim[0]), (1, 0)), 'tr': ((xlim[1], ylim[1]), (1, 1)), 'bl': ((xlim[0], ylim[0]), (0, 0))}
    for c in corners:
        a, b = pts[c]
        ax.add_artist(ConnectionPatch(xyA=a, coordsA='data', axesA=ax, xyB=b, coordsB='axes fraction', axesB=axi, color='black', lw=0.9, zorder=7, clip_on=False))
    return axi


def save(fig, name):
    fig.savefig(OUT + name, bbox_inches='tight', pad_inches=0.12); plt.close(fig); print('wrote', name)
