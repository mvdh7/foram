import numpy as np
import foram


# Tests don't work because correct effective radii are not yet available


def test_Dc_HCO3():
    assert np.isclose(
        np.round(foram.diffusion.Dc(25, ion="HCO3", s=0), 2),
        1.18e-9,
    )


def test_Dc_CO3():
    assert np.isclose(
        np.round(foram.diffusion.Dc(25, ion="CO3", s=0), 3),
        0.955e-9,
    )


def test_Dc_H():
    assert np.isclose(np.round(foram.diffusion.Dc(25, ion="H", s=0), 2), 9.31e-9)


# maybe these are also in 10^-9 but it is not clear from the paper and I don't know enough to say what is sensible
def test_Dc_BOH3():
    assert np.isclose(
        np.round(foram.diffusion.Dc(25, ion="BOH3", s=35), 2),
        1.11,
    )


def test_Dc_BOH4():
    assert np.isclose(
        np.round(foram.diffusion.Dc(25, ion="BOH4", s=35), 2),
        0.97,
    )


# test_Dc_HCO3()
# test_Dc_CO3()
# test_Dc_H()
# test_Dc_BOH3()
# test_Dc_BOH4()
