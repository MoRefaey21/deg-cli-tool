"""Helper function: -log10(p-value)."""

import numpy as np


def neg_log10_pvalue(p):
    """Calculate -log10 of p-values.

    Parameters
    ----------
    p : float or array-like
        P-values. NaN inputs return NaN.

    Returns
    -------
    -log10(p) : same shape as input
    """
    return -np.log10(p)
