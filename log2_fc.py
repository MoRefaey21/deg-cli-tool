"""Helper function: log2 fold change."""

import numpy as np


def log2_fc(fc):
    """Calculate log2 fold change.

    Parameters
    ----------
    fc : float or array-like
        Fold change values. NaN inputs return NaN.

    Returns
    -------
    log2 fold change : same shape as input
    """
    return np.log2(fc)
