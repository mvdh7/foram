import numpy as np
import foram


def test_k_p1():
    assert np.isclose(np.round(foram.rates.k_p1(25, s=33.77), 3), 0.036)


# test_k_p1()

def test_Kw():
    assert np.isclose(np.round(np.log(foram.rates.Kw(25, s=35)), 3), -30.434)
    
def test_k_m4():
    assert np.isclose(np.round(foram.rates.Kw(25, s=35), 3), 3.7*1e-4)