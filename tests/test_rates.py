import numpy as np
import foram


def test_k_p1():
    assert np.isclose(np.round(foram.rates.k_p1(25, s=33.77), 3), 0.036)


def test_k_p5():
    assert np.isclose(foram.rates.k_p5(25), 1e10)


def test_k_p6():
    assert np.isclose(foram.rates.k_p6(25), 1.3e-3)


def test_k_m7():
    assert np.isclose(foram.rates.k_m7(25), 1e10)

# test_k_p1()
# test_k_p5()
# test_k_p6()
# test_k_m7()
