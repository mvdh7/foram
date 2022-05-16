import numpy as np
import foram


def test_k_p1():
    assert np.isclose(np.round(foram.rates.k_p1(25, s=33.77), 3), 0.036)

def test_k_m1():
    assert np.isclose(np.round(foram.rates.k_m1(25, s=33.77), 1), 3.3e4)


def test_K1():
    assert np.isclose(np.round(foram.rates.K1(25, s=35), 4), -13.4847)

foram.rates.K1(25, s=35)    
#test_K1()   
#test_k_m1()

#foram.rates.k_m1(25, s=33.77)