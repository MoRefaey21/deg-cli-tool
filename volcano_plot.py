"""Helper function: generate a volcano plot."""

import matplotlib.pyplot as plt
import numpy as np


def volcano_plot(x, y, out_file="volcano_plot.png",
                 log2fc_threshold=np.log2(1.5),
                 xlim=(-1, 1), ylim=(0, 5), alpha=0.05, significant=None):
    """Generate a volcano plot of log2 fold change vs -log10(p-value).

    A gene is coloured only if it meets BOTH the fold change threshold
    AND the significance threshold. Everything else is grey.

    Parameters
    ----------
    x : array-like
        log2 fold change values.
    y : array-like
        -log10 p-values.
    out_file : str
        Path for the saved PNG.
    log2fc_threshold : float
        Threshold for over/under-expression on the log2 fold change axis.
    xlim, ylim : tuple
        Axis limits.
    alpha : float
        Significance threshold; draws the horizontal reference line at
        -log10(alpha) and is used to colour significant genes when
        `significant` is not given.
    significant : array-like of bool, optional
        Per-gene significance call (e.g. FDR q < alpha), aligned with x and y.
        If given it replaces the raw p-value test, and the horizontal
        reference line is drawn at the smallest -log10(p) still called
        significant (omitted if no gene is significant).
    """
    try:
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)

        # Drop non-finite values
        mask = np.isfinite(x) & np.isfinite(y)
        if significant is not None:
            significant = np.asarray(significant, dtype=bool)[mask]
        x = x[mask]
        y = y[mask]

        if significant is None:
            sig_threshold = -np.log10(alpha)
            sig_by_p = y > sig_threshold
        else:
            sig_by_p = significant
            sig_threshold = y[sig_by_p].min() if sig_by_p.any() else None

        over  = (x >  log2fc_threshold) & sig_by_p
        under = (x < -log2fc_threshold) & sig_by_p
        ns    = ~over & ~under

        plt.figure(figsize=(12, 10))
        plt.scatter(x[ns],    y[ns],    color="grey",    s=5, label="Not significant")
        plt.scatter(x[over],  y[over],  color="#D55E00", s=8, label="Over-expressed")
        plt.scatter(x[under], y[under], color="#0072B2", s=8, label="Under-expressed")

        plt.xlabel("log2 Fold Change")
        plt.ylabel("-log10 P-Values")
        plt.title("Volcano Plot: log2 FC vs -log10 P-Values")
        plt.legend(loc="upper right")
        plt.xlim(xlim)
        plt.ylim(ylim)
        if sig_threshold is not None:
            plt.axhline(y=sig_threshold,  color="b", linestyle="--", alpha=0.3)
        plt.axvline(x= log2fc_threshold,  color="b", linestyle="--", alpha=0.3)
        plt.axvline(x=-log2fc_threshold,  color="b", linestyle="--", alpha=0.3)
        plt.tight_layout()
        plt.savefig(out_file, dpi=300)
        plt.close()
        print(f"Volcano plot saved as {out_file}.")
    except Exception as e:
        print(f"Error generating volcano plot: {e}")
