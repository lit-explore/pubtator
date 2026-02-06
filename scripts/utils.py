"""
Shared utilities for PubTator pipeline scripts
"""

import numpy as np


def numpy_encoder(obj: object) -> int | list:
    """
    JSON encoder to convert numpy int64 elements to generic Python ints.

    Usage:
        json.dumps(data, default=numpy_encoder)
    """
    if isinstance(obj, np.generic):
        return obj.item()
    raise TypeError(f"Object of type {type(obj).__name__} is not JSON serializable")
