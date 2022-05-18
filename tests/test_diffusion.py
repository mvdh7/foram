import numpy as np
import foram


# Tests don't work because correct effective radii are not yet available

def test_Dc_HCO3():
    assert np.isclose(
        np.round(foram.diffusion.Dc(25, ion="HCO3", s=0), 11),
        1.18e-9,
        atol=1e-10
    )


def test_Dc_CO3():
    assert np.isclose(
        np.round(foram.diffusion.Dc(25, ion="CO3", s=0), 12),
        0.955e-9,
        atol=1e-12
    )


def test_Dc_H():
    assert np.isclose(
        np.round(foram.diffusion.Dc(25, ion="H", s=0), 11), 
        9.31e-9,
        atol=1e-11
    )
    
def test_Dc_OH():
    assert np.isclose(
        np.round(foram.diffusion.Dc(25, ion="OH", s=0), 11), 
        5.27e-9,
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
        atol=1e-10
    )



print("HCO3:",np.round(foram.diffusion.Dc(25, ion="HCO3", s=0), 11),"   Test:",1.18e-9)
print("CO3:",np.round(foram.diffusion.Dc(25, ion="CO3", s=0), 12),"   Test:",0.955e-9)
print("H:",np.round(foram.diffusion.Dc(25, ion="H", s=0), 11),"   Test:",9.31e-9)
print("OH:",np.round(foram.diffusion.Dc(25, ion="OH", s=0), 11),"   Test:",5.27e-9)
print("BOH4:",np.round(foram.diffusion.Dc(25, ion="BOH4", s=35), 11),"   Test:",0.97e-9)
# test_Dc_HCO3()
# test_Dc_CO3()
# test_Dc_H()
# test_Dc_BOH3()
# test_Dc_BOH4()
