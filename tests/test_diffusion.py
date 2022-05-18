import numpy as np
import foram


# Tests don't work because correct effective radii are not yet available

def test_Dc_HCO3():
    assert np.isclose(
        np.round(foram.diffusion.Dc(25, ion="HCO3", s=0), 11),
        1.18e-9,
        atol=1e-11
    )


def test_Dc_CO3():
    assert np.isclose(
        np.round(foram.diffusion.Dc(25, ion="CO3", s=0), 12),
        0.955e-9,
        atol=1e-12
    )


def test_Dc_H():
    assert np.isclose(
        np.round(foram.diffusion.Dc(25, ion="H", s=0), 2), 
        9.31e-9,
        atol=1e-11
    )


def test_Dc_BOH3():
    assert np.isclose(
        np.round(foram.diffusion.Dc(25, ion="BOH3", s=35), 11),
        1.11e-9,
        atol=1e-11
    )


def test_Dc_BOH4():
    assert np.isclose(
        np.round(foram.diffusion.Dc(25, ion="BOH4", s=35), 11),
        0.97e-9,
        atol=1e-11
    )


# test_Dc_HCO3()
# test_Dc_CO3()
# test_Dc_H()
# test_Dc_BOH3()
# test_Dc_BOH4()
