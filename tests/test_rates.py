import numpy as np
import foram


def test_k_p1():
    assert np.isclose(np.round(foram.rates.k_p1(25, s=33.77), 3), 0.036)


# test_k_p1()
