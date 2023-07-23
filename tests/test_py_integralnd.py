#!/usr/bin/env python
#
# example_simple.py
#
# Mon May 22 16:00:00 2023
# Copyright  2023  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# example_simple.py
# Copyright (C) 2023 Sandro Dias Pinto Vitenti <vitenti@uel.br>
#
# numcosmo is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# numcosmo is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program.  If not, see <http://www.gnu.org/licenses/>.

"""Simple example computing cosmological distances."""

from typing import Tuple, List

import numpy as np
from numpy.testing import assert_almost_equal

from numcosmo_py import Ncm

#
#  Initializing the library objects, this must be called before
#  any other library function.
#
Ncm.cfg_init()


class IntegralND(Ncm.IntegralND):
    """Test class for IntegralND."""

    def __init__(self, w: List[float], **kwargs) -> None:
        """Initialize."""
        super().__init__(**kwargs)
        self.w = np.array(w)

    # pylint: disable-next=arguments-differ
    def do_get_dimensions(self) -> Tuple[int, int]:
        """Get number of dimensions."""
        return 3, 3

    # pylint: disable-next=arguments-differ
    def do_integrand(
        self,
        x_vec: Ncm.Vector,
        dim: int,
        npoints: int,
        _fdim: int,
        fval_vec: Ncm.Vector,
    ) -> None:
        """Integrand function."""

        x = np.asarray(x_vec.dup_array()).reshape((npoints, dim))
        fval = self.w * x * np.sin(self.w * x)

        fval_vec.set_array(fval.flatten())


def test_integral_nd() -> None:
    """Example computing cosmological distances."""
    test_f = IntegralND(w=[1.0, 2.0, 3.0], method=Ncm.IntegralNDMethod.H_V)

    res = Ncm.Vector.new(3)
    err = Ncm.Vector.new(3)

    test_f.set_method(Ncm.IntegralNDMethod.P_V)
    assert test_f.get_method() == Ncm.IntegralNDMethod.P_V
    test_f.eval(
        Ncm.Vector.new_array([0.0, 0.0, 0.0]),
        Ncm.Vector.new_array([1.0, 1.0, 1.0]),
        res,
        err,
    )
    assert_almost_equal(res.dup_array(), np.sin(test_f.w) / test_f.w - np.cos(test_f.w))

    test_f.set_method(Ncm.IntegralNDMethod.H_V)
    assert test_f.get_method() == Ncm.IntegralNDMethod.H_V
    test_f.eval(
        Ncm.Vector.new_array([0.0, 0.0, 0.0]),
        Ncm.Vector.new_array([1.0, 1.0, 1.0]),
        res,
        err,
    )
    assert_almost_equal(res.dup_array(), np.sin(test_f.w) / test_f.w - np.cos(test_f.w))

    test_f.set_method(Ncm.IntegralNDMethod.P)
    assert test_f.get_method() == Ncm.IntegralNDMethod.P
    test_f.eval(
        Ncm.Vector.new_array([0.0, 0.0, 0.0]),
        Ncm.Vector.new_array([1.0, 1.0, 1.0]),
        res,
        err,
    )
    assert_almost_equal(res.dup_array(), np.sin(test_f.w) / test_f.w - np.cos(test_f.w))

    test_f.set_method(Ncm.IntegralNDMethod.H)
    assert test_f.get_method() == Ncm.IntegralNDMethod.H
    test_f.eval(
        Ncm.Vector.new_array([0.0, 0.0, 0.0]),
        Ncm.Vector.new_array([1.0, 1.0, 1.0]),
        res,
        err,
    )
    assert_almost_equal(res.dup_array(), np.sin(test_f.w) / test_f.w - np.cos(test_f.w))


def test_integral_nd_props() -> None:
    """Example computing cosmological distances."""
    test_f = IntegralND(w=[1.0, 2.0, 3.0], method=Ncm.IntegralNDMethod.H_V)

    test_f.set_abstol(1e-6)
    assert_almost_equal(test_f.get_abstol(), 1e-6)

    test_f.set_reltol(1e-6)
    assert_almost_equal(test_f.get_reltol(), 1e-6)

    test_f.set_method(Ncm.IntegralNDMethod.P_V)
    assert test_f.get_method() == Ncm.IntegralNDMethod.P_V

    test_f.set_maxeval(1000)
    assert test_f.get_maxeval() == 1000

    test_f.set_error(Ncm.IntegralNDError.INDIVIDUAL)
    assert test_f.get_error() == Ncm.IntegralNDError.INDIVIDUAL


def test_integral_nd_errors() -> None:
    """Example computing cosmological distances."""
    test_f = IntegralND(w=[1.0, 2.0, 3.0], method=Ncm.IntegralNDMethod.H_V)

    res = Ncm.Vector.new(3)
    err = Ncm.Vector.new(3)

    test_f.set_error(Ncm.IntegralNDError.PAIRWISE)
    assert test_f.get_error() == Ncm.IntegralNDError.PAIRWISE
    test_f.eval(
        Ncm.Vector.new_array([0.0, 0.0, 0.0]),
        Ncm.Vector.new_array([1.0, 1.0, 1.0]),
        res,
        err,
    )
    assert_almost_equal(res.dup_array(), np.sin(test_f.w) / test_f.w - np.cos(test_f.w))

    test_f.set_error(Ncm.IntegralNDError.L2)
    assert test_f.get_error() == Ncm.IntegralNDError.L2
    test_f.eval(
        Ncm.Vector.new_array([0.0, 0.0, 0.0]),
        Ncm.Vector.new_array([1.0, 1.0, 1.0]),
        res,
        err,
    )
    assert_almost_equal(res.dup_array(), np.sin(test_f.w) / test_f.w - np.cos(test_f.w))

    test_f.set_error(Ncm.IntegralNDError.L1)
    assert test_f.get_error() == Ncm.IntegralNDError.L1
    test_f.eval(
        Ncm.Vector.new_array([0.0, 0.0, 0.0]),
        Ncm.Vector.new_array([1.0, 1.0, 1.0]),
        res,
        err,
    )
    assert_almost_equal(res.dup_array(), np.sin(test_f.w) / test_f.w - np.cos(test_f.w))

    test_f.set_error(Ncm.IntegralNDError.LINF)
    assert test_f.get_error() == Ncm.IntegralNDError.LINF
    test_f.eval(
        Ncm.Vector.new_array([0.0, 0.0, 0.0]),
        Ncm.Vector.new_array([1.0, 1.0, 1.0]),
        res,
        err,
    )
    assert_almost_equal(res.dup_array(), np.sin(test_f.w) / test_f.w - np.cos(test_f.w))


def test_integral_nd_serialize() -> None:
    """Example computing cosmological distances."""
    test_f = IntegralND(w=[1.0, 2.0, 3.0], method=Ncm.IntegralNDMethod.H_V)

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.NONE)

    test_f.set_method(Ncm.IntegralNDMethod.H)
    test_f.set_abstol(1e-6)
    test_f.set_reltol(1e-6)
    test_f.set_maxeval(11234)
    test_f.set_error(Ncm.IntegralNDError.LINF)
    test_f_dup = ser.dup_obj(test_f)

    assert isinstance(test_f_dup, IntegralND)

    assert test_f_dup.get_method() == test_f.get_method()
    assert test_f_dup.get_abstol() == test_f.get_abstol()
    assert test_f_dup.get_reltol() == test_f.get_reltol()
    assert test_f_dup.get_maxeval() == test_f.get_maxeval()
    assert test_f_dup.get_error() == test_f.get_error()
