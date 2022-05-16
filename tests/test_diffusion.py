import numpy as np
import foram


def test_Dc_HCO3():  
    assert np.isclose(np.round(np.log(foram.diffusion.Dc(25, a=3, s=0)), 3), 1.18e-9)

