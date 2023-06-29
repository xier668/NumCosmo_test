#
# gauss_constraint.py
#
# Wed May 3 10:21:00 2023
# Copyright  2023  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# gauss_constraint.py
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

"""Example of using a Gaussian distribution with positivity constraint to
test the MCMC sampler.
"""

from pathlib import Path
from typing import Optional

import numpy as np

from numcosmo_py import Ncm
from numcosmo_py.sampling.esmcmc import (
    create_esmcmc,
    WalkerTypes,
    InterpolationKernel,
    InterpolationMethod,
)


from numcosmo_py.external.minimax_tilting_sampler import TruncatedMVN


def run_gauss_constraint_mcmc(
    *,
    dim: int = 10,
    ssize: int = 5000000,
    verbose: bool = True,
    fit_first: bool = False,
    robust: bool = False,
    use_apes_interpolation: bool = True,
    use_apes_threads: Optional[bool] = None,
    sampler: WalkerTypes = WalkerTypes.APES,
    interpolation_method: InterpolationMethod = InterpolationMethod.VKDE,
    interpolation_kernel: InterpolationKernel = InterpolationKernel.CAUCHY,
    nwalkers: int = 3000,
    nthreads: int = 4,
    over_smooth: float = 1.1,
    local_fraction: Optional[float] = None,
    init_sampling_scale: float = 1.0,
    start_catalog: Optional[Path] = None,
    tmvn: bool = False,
) -> str:
    """Runs the Funnel MCMC example."""

    mgc = Ncm.ModelMVND.new(dim)
    for i in range(dim):
        if i % 2 == 0:
            mgc.param_set_lower_bound(i, 0.0)
        else:
            mgc.param_set_lower_bound(i, -10.0)
        mgc.param_set_upper_bound(i, 10.0)
        mgc.param_set_scale(i, 0.1)
        mgc.orig_vparam_set(0, i, 0.5)

    mset = Ncm.MSet.empty_new()
    mset.set(mgc)
    mset.param_set_all_ftype(Ncm.ParamType.FREE)
    mset.prepare_fparam_map()

    dgc = Ncm.DataGaussCovMVND.new(dim)
    mean = Ncm.Vector.new(dim)
    mean.set_all(0.0)

    if dim > 50:
        raise ValueError("dim > 50 not supported")

    # Always use the same seed and size to get the same results
    # for all dimensions.
    cov50 = Ncm.Matrix.new(50, 50)
    rng = Ncm.RNG.seeded_new(None, 0)
    cov50.fill_rand_cor(10.0, rng)
    cov = cov50.get_submatrix(0, 0, dim, dim).dup()
    # cov50.log_vals("", "% 22.15g")

    dgc.set_cov_mean(mean, cov)
    m2lnN = dgc.get_log_norma(mset)
    print(f"# Constant normalization {m2lnN}")

    if tmvn:
        Ncm.message_str("# Sampling the likelihood posterior using TruncatedMVN: \n")
        rng = Ncm.RNG.seeded_new(None, 0)
        bounds = np.array(
            [
                [mgc.param_get_lower_bound(i), mgc.param_get_upper_bound(i)]
                for i in range(dim)
            ]
        )
        mset.fparams_set_vector(mean)

        cov_np = np.array(cov.dup_array())
        cov_np.shape = (dim, dim)

        tmvn_sampler = TruncatedMVN(
            mu=np.zeros(dim), cov=cov_np, lb=bounds[:, 0], ub=bounds[:, 1], seed=0
        )

        sv = Ncm.StatsVec.new(dim + 1, Ncm.StatsVecType.COV, True)
        for s in np.transpose(tmvn_sampler.sample(ssize)):
            dgc.peek_mean().set_array(s)
            st = np.concatenate(([dgc.m2lnL_val(mset)], s))
            sv.append(Ncm.Vector.new_array(st), True)

        filename = f"gauss_constraint_{dim}d_tmvn_samples.dat"
        with open(filename, "w", encoding="utf-8") as f:
            f.write(f"#   Number of samples : {sv.nitens}\n")
            f.write(
                f"#   Mean : {' '.join([f'{v: 22.15g}' for v in sv.peek_mean().dup_array()])}\n"
            )
            f.write(
                f"#   Stdev: {' '.join([f'{sv.get_sd(i): 22.15g}' for i in range(dim)])}\n"
            )

            for i in range(sv.nitens):
                f.write(
                    " ".join([f"{v: 22.15g}" for v in sv.peek_row(i).dup_array()])
                    + "\n"
                )

        return filename

    dset = Ncm.Dataset.new()
    dset.append_data(dgc)
    likelihood = Ncm.Likelihood.new(dset)

    start_mcat = None
    if start_catalog is not None:
        start_mcat = Ncm.MSetCatalog.new_from_file_ro(start_catalog.as_posix(), 0)

    esmcmc = create_esmcmc(
        likelihood,
        mset,
        f"gauss_constraint_{dim}d",
        verbose=verbose,
        fit_first=fit_first,
        robust=robust,
        use_apes_interpolation=use_apes_interpolation,
        use_apes_threads=use_apes_threads,
        sampler=sampler,
        interpolation_method=interpolation_method,
        interpolation_kernel=interpolation_kernel,
        nwalkers=nwalkers,
        nthreads=nthreads,
        over_smooth=over_smooth,
        local_fraction=local_fraction,
        init_sampling_scale=init_sampling_scale,
        start_mcat=start_mcat,
    )

    # Running the esmcmc.
    esmcmc.start_run()
    esmcmc.run(ssize // nwalkers)
    esmcmc.end_run()

    fit = esmcmc.peek_fit()
    filename = esmcmc.peek_catalog().peek_filename()
    if verbose:
        esmcmc.mean_covar()
        fit.log_covar()

    return filename
