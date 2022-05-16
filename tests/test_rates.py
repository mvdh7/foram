import numpy as np
import foram


def test_kHSO4_FREE_D90a():
    assert np.isclose(np.round(np.log(foram.rates.kHSO4_FREE_D90a(25, s=35)), 2), -2.30)


def test_total_sulfate():
    assert np.isclose(np.round(foram.rates.total_sulfate(s=35), 5), 0.02824)


def test_K1():
    assert np.isclose(
        np.round(
            np.log(foram.rates.K1(25, s=35) * foram.rates.pH_free_to_total(25, s=35)), 4
        ),
        -13.4847,
    )


def test_K2():
    assert np.isclose(
        np.round(
            np.log(foram.rates.K2(25, s=35) * foram.rates.pH_free_to_total(25, s=35)), 4
        ),
        -20.5504,
    )


def test_Kw():
    assert np.isclose(
        np.round(
            np.log(foram.rates.Kw(25, s=35) * foram.rates.pH_free_to_total(25, s=35)), 3
        ),
        -30.434,
    )


def test_k_p1():
    assert np.isclose(np.round(foram.rates.k_p1(25, s=33.77), 3), 0.036)


def test_k_m1():
    assert np.isclose(np.round(foram.rates.k_m1(25, s=33.77), 2), 3.35e4)
    # Had to add extra decimal place - incorrect rounding in paper?


def test_k_p4():
    assert np.isclose(np.round(foram.rates.k_p4(25, s=0), 3), 8500)


def test_k_m4():
    assert np.isclose(np.round(foram.rates.Kw(25, s=35), 3), 3.7 * 1e-4)


def test_k_p5():
    assert np.isclose(foram.rates.k_p5(25), 1e10)


def test_k_m5():
    assert np.isclose(np.round(foram.rates.Kw(25, s=35), 3), 9)


def test_k_p6():
    assert np.isclose(foram.rates.k_p6(25), 1.3e-3)


def test_k_m6():
    assert np.isclose(foram.rates.k_m6(25), 2.8e10)


def test_k_m7():
    assert np.isclose(foram.rates.k_m7(25), 1e10)


# test_kHSO4_FREE_D90a()
# test_total_sulfate()
# test_K1()
# test_K2()
# test_Kw()
# test_k_p1()
test_k_m1()
test_k_p4()
test_k_m4()
# test_k_p5()
test_k_m5()
# test_k_p6()
test_k_m6()
# test_k_m7()
