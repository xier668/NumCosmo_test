#
# mcat_calibrate_apes.py
#
# Wed Feb 8 10:00:00 2023
# Copyright  2023  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# mcat_calibrate_apes.py
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

"""NumCosmoPy APES calibration script."""

from pathlib import Path
import math
import numpy as np
import matplotlib.pyplot as plt
import typer

from numcosmo_py import Ncm
from numcosmo_py.plotting.tools import confidence_ellipse
from numcosmo_py.interpolation.stats_dist import (
    create_stats_dist,
    InterpolationMethod,
    InterpolationKernel,
)


Ncm.cfg_init()
app = typer.Typer()


@app.command()
def calibrate_catalog(
    catalog: Path,
    *,
    robust: bool = False,
    interpolation_method: InterpolationMethod = InterpolationMethod.VKDE,
    interpolation_kernel: InterpolationKernel = InterpolationKernel.CAUCHY,
    over_smooth: float = 1.0,
    interpolate: bool = True,
    ntries: int = 100,
):
    """Calibrate the APES sampler using a given catalog."""

    if not catalog.exists() or not catalog.is_file():
        raise typer.BadParameter(f"Catalog {catalog} does not exist.")

    mcat = Ncm.MSetCatalog.new_from_file_ro(catalog.as_posix(), 0)
    mcat_len = mcat.len()
    nwalkers = mcat.nchains()
    nadd_vals = mcat.nadd_vals()
    m2lnL_id = mcat.get_m2lnp_var()  # pylint: disable-msg=invalid-name

    assert mcat_len > nwalkers
    assert over_smooth != 0.0

    last_e = [mcat.peek_row(mcat_len - nwalkers + i) for i in range(nwalkers)]
    ncols = mcat.ncols()
    nvar = ncols - nadd_vals
    params = ["$" + mcat.col_symb(i) + "$" for i in range(nadd_vals, mcat.ncols())]

    sdist = create_stats_dist(
        robust=robust,
        interpolation_method=interpolation_method,
        interpolation_kernel=interpolation_kernel,
        dim=nvar,
        over_smooth=math.fabs(over_smooth),
    )

    if over_smooth < 0.0:
        sdist.set_cv_type(Ncm.StatsDistCV.SPLIT)

    m2lnL = []  # pylint: disable-msg=invalid-name
    for row in last_e:
        m2lnL.append(row.get(m2lnL_id))
        sdist.add_obs(row.get_subvector(nadd_vals, nvar))

    m2lnL_v = Ncm.Vector.new_array(m2lnL)  # pylint: disable-msg=invalid-name
    if interpolate:
        sdist.prepare_interp(m2lnL_v)
    else:
        sdist.prepare()

    ovs = sdist.get_over_smooth()
    print(f"# === Setting over smooth to {ovs}")

    rng = Ncm.RNG.new()
    var_vector = Ncm.Vector.new(nvar)

    try_sample_array = []
    for _ in range(ntries):
        sdist.sample(var_vector, rng)
        try_sample_array.append(var_vector.dup_array())

    try_sample = np.array(try_sample_array)

    weights = np.array(sdist.peek_weights().dup_array())
    weights = weights / np.sum(weights)
    max_w = np.max(weights[np.nonzero(weights)])
    min_w = np.min(weights[np.nonzero(weights)])

    print(f"# === Min weight: {min_w}")
    print(f"# === Max weight: {max_w}")
    print(f"# === Mean weight: {np.mean(weights[np.nonzero(weights)])}")
    print(f"# === Median weight: {np.median(weights[np.nonzero(weights)])}")
    print(f"# === Non-zero weights: {np.count_nonzero(weights)}")

    for a in range(nvar):  # pylint: disable-msg=invalid-name
        for b in range(a + 1, nvar):  # pylint: disable-msg=invalid-name
            indices = np.array([a, b])
            print(indices)

            _, axis = plt.subplots(1, 1, figsize=(16, 8))

            # pylint: disable-next=invalid-name
            for ii in range(0, int(sdist.get_n_kernels())):
                y_i, cov_i, _, w_i = sdist.get_Ki(ii)
                mean = np.array(y_i.dup_array())
                cov = np.array([[cov_i.get(i, j) for j in indices] for i in indices])
                cov = cov * 1.0
                w_i = weights[ii]

                if w_i > 0.0:
                    confidence_ellipse(
                        mean[indices],
                        cov,
                        axis,
                        edgecolor="red",
                        facecolor="red",
                    )

            axis.scatter(try_sample[:, a], try_sample[:, b])
            plt.axis("auto")
            plt.xlabel(params[a])
            plt.ylabel(params[b])
            plt.grid()
            plt.show()


if __name__ == "__main__":
    app()
