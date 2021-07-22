/* -*- Mode: C; indent-tabs-mode: nil; c-basic-offset: 2; tab-width: 2 -*-  */

/***************************************************************************
 *            ncm_stats_dist_nd_kde_studentt.c
 *
 *  Wed November 07 17:41:47 2018
 *  Copyright  2018  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * ncm_stats_dist_nd_kde_studentt.c
 * Copyright (C) 2018 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
 *
 * numcosmo is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * numcosmo is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * SECTION:ncm_stats_dist_nd_kde_studentt
 * @title: NcmStatsDistNdKDEStudentt
 * @short_description: An N dimensional probability distributions using Multivariate Student's t KDE
 *
 * This object provides the tools to perform a radial basis interpolation
 * in a multidimensional function, using a multivariate Student's t distribution.
 * The goal is to use a polinomial function generated by the interpolation,
 * such that we are able to sample the original functio from this new interpolation function.
 * For more informations about radial basis interpolation,
 * check #NcmStatsDistNd.
 * A brief description of the Multivariate Studentt function can be found below.
 * For more information, check [[The R Journal Vol. 5/2, December 2013](https://journal.r-project.org/archive/2013/RJ-2013-033/RJ-2013-033.pdf)]
 *
 * We use the Multivariate Student's t function as the radial basis function. The function has the stocastic representation given by
 *
 * \begin{align}
 * \boldsymbol{X}=\boldsymbol{\mu}+\sqrt{W} A \boldsymbol{Z}
 * ,\end{align}
 * where $W= \nu / \chi_{\nu}^{2}$, $\nu$ represents the degree of freedom of a chi-squared distribution $\chi_{\nu}^2$,
 * $\boldsymbol{Z}$ is a p-dimensional random vector, $A$ is a $p \times p$ matrix and $\mu$ is the mean vector.
 *
 * $\boldsymbol{X}$ is fully determined by the covariance matrix $\Sigma = \boldsymbol{A}\boldsymbol{A}^t$ and the mean vector $\mu$.
 * Assuming that the covariance matrix is positive definite, $\boldsymbol{X}$ has the probability density
 *
 * \begin{align}
 * \phi(x) &= \frac{\Gamma[(\nu+p) / 2]}{\Gamma(\nu / 2) \nu^{p / 2} \pi^{p / 2}|\mathbf{\Sigma}|^{1 / 2}}\left[1+\frac{1}{\nu}
 * (\mathbf{x}-\boldsymbol{\mu})^{T} \boldsymbol{\Sigma}^{-1}(\mathbf{x}-\boldsymbol{\mu})\right]^{-(\nu+p) / 2},
 * \end{align}
 * where $p$ is the dimension and $x$ are the points to be evaluated.
 *
 * Once this object is initialized, we may use the methods in the #NcmStasDistNd class to perform the interpolation
 * and to generate a sample from the interpolated polinomal function.
 *
 * The user must provide the following input values: $p$ - ncm_stats_dist_nd_kde_studentt_new(),
 * cv_type - ncm_stats_dist_nd_kde_studentt_new(), $\nu$ - ncm_stats_dist_nd_kde_studentt_new().
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "gsl/gsl_sf_result.h"

#include "math/ncm_stats_dist_nd_kde_studentt.h"
#include "math/ncm_stats_vec.h"
#include "math/ncm_c.h"
#include "math/ncm_lapack.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_blas.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_sf_gamma.h>
#include "levmar/levmar.h"
#endif /* NUMCOSMO_GIR_SCAN */

struct _NcmStatsDistNdKDEStudenttPrivate
{
  gdouble nu;
};

enum
{
  PROP_0,
  PROP_NU,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcmStatsDistNdKDEStudentt, ncm_stats_dist_nd_kde_studentt, NCM_TYPE_STATS_DIST_ND);

static void
ncm_stats_dist_nd_kde_studentt_init (NcmStatsDistNdKDEStudentt *dndt)
{
  NcmStatsDistNdKDEStudenttPrivate * const self = dndt->priv = ncm_stats_dist_nd_kde_studentt_get_instance_private (dndt);
  
  self->nu = 0.0;
}

static void
_ncm_stats_dist_nd_kde_studentt_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmStatsDistNdKDEStudentt *dndt = NCM_STATS_DIST_ND_KDE_STUDENTT (object);
  
  /*NcmStatsDistNdKDEStudenttPrivate * const self = dndt->priv;*/
  
  g_return_if_fail (NCM_IS_STATS_DIST_ND_KDE_STUDENTT (object));
  
  switch (prop_id)
  {
    case PROP_NU:
      ncm_stats_dist_nd_kde_studentt_set_nu (dndt, g_value_get_double (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_stats_dist_nd_kde_studentt_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmStatsDistNdKDEStudentt *dndt = NCM_STATS_DIST_ND_KDE_STUDENTT (object);
  
  /*NcmStatsDistNdKDEStudenttPrivate * const self = dndt->priv;*/
  
  g_return_if_fail (NCM_IS_STATS_DIST_ND_KDE_STUDENTT (object));
  
  switch (prop_id)
  {
    case PROP_NU:
      g_value_set_double (value, ncm_stats_dist_nd_kde_studentt_get_nu (dndt));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_stats_dist_nd_kde_studentt_dispose (GObject *object)
{
  /*NcmStatsDistNdKDEStudentt *dndt               = NCM_STATS_DIST_ND_KDE_STUDENTT (object);*/
  /*NcmStatsDistNdKDEStudenttPrivate * const self = dndt->priv;*/
  
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_stats_dist_nd_kde_studentt_parent_class)->dispose (object);
}

static void
_ncm_stats_dist_nd_kde_studentt_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_stats_dist_nd_kde_studentt_parent_class)->finalize (object);
}

static gdouble _ncm_stats_dist_nd_kde_studentt_get_rot_bandwidth (NcmStatsDistNd *dnd, const guint d, const gdouble n);
static gdouble _ncm_stats_dist_nd_kde_studentt_get_kernel_lnnorm (NcmStatsDistNd *dnd, NcmMatrix *cov_decomp, const guint d, const gdouble n, const NcmVector *href);
static void _ncm_stats_dist_nd_kde_studentt_prepare_IM (NcmStatsDistNd *dnd, GPtrArray *invUsample, const gint d, const gint n, const NcmVector *href, NcmMatrix *IM);

static gdouble _ncm_stats_dist_nd_kde_studentt_eval (NcmStatsDistNd *dnd, NcmVector *weights, NcmVector *invUy, GPtrArray *invUsample, const gint d, const gint n, const NcmVector *href);
static void _ncm_stats_dist_nd_kde_studentt_kernel_sample (NcmStatsDistNd *dnd, NcmMatrix *cov_decomp, const guint d, NcmVector *y, NcmVector *mu, const NcmVector *href, NcmRNG *rng);

static void
ncm_stats_dist_nd_kde_studentt_class_init (NcmStatsDistNdKDEStudenttClass *klass)
{
  GObjectClass *object_class     = G_OBJECT_CLASS (klass);
  NcmStatsDistNdClass *dnd_class = NCM_STATS_DIST_ND_CLASS (klass);
  
  object_class->set_property = &_ncm_stats_dist_nd_kde_studentt_set_property;
  object_class->get_property = &_ncm_stats_dist_nd_kde_studentt_get_property;
  object_class->dispose      = &_ncm_stats_dist_nd_kde_studentt_dispose;
  object_class->finalize     = &_ncm_stats_dist_nd_kde_studentt_finalize;
  
  g_object_class_install_property (object_class,
                                   PROP_NU,
                                   g_param_spec_double ("nu",
                                                        NULL,
                                                        "nu value of the function",
                                                        1.0, G_MAXDOUBLE, 3.0,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  
  dnd_class->get_rot_bandwidth = &_ncm_stats_dist_nd_kde_studentt_get_rot_bandwidth;
  dnd_class->get_kernel_lnnorm = &_ncm_stats_dist_nd_kde_studentt_get_kernel_lnnorm;
  dnd_class->prepare_IM        = &_ncm_stats_dist_nd_kde_studentt_prepare_IM;
  
  dnd_class->eval          = &_ncm_stats_dist_nd_kde_studentt_eval;
  dnd_class->kernel_sample = &_ncm_stats_dist_nd_kde_studentt_kernel_sample;
}

static gdouble
_ncm_stats_dist_nd_kde_studentt_f (NcmStatsDistNdKDEStudenttPrivate * const self, const guint d, gdouble chi2)
{
  return pow (1.0 + chi2 / self->nu, -0.5 * (self->nu + d));
}

static gdouble
_ncm_stats_dist_nd_kde_studentt_get_rot_bandwidth (NcmStatsDistNd *dnd, const guint d, const gdouble n)
{
  NcmStatsDistNdKDEStudentt *dndt               = NCM_STATS_DIST_ND_KDE_STUDENTT (dnd);
  NcmStatsDistNdKDEStudenttPrivate * const self = dndt->priv;
  const gdouble nu                              = (self->nu >= 3.0) ? self->nu : 3.0;
  
  return pow (
    16.0 * gsl_pow_2 (nu - 2) * (1.0 + d + nu) * (3.0 + d + nu) /
    ((2.0 + d) * (d + nu) * (2.0 + d + nu) * (d + 2.0 * nu) * (2.0 + d + 2.0 * nu) * n),
    1.0 / (d + 4.0));
}

static gdouble
_ncm_stats_dist_nd_kde_studentt_get_kernel_lnnorm (NcmStatsDistNd *dnd, NcmMatrix *cov_decomp, const guint d, const gdouble n, const NcmVector *href)
{
  NcmStatsDistNdKDEStudentt *dndt               = NCM_STATS_DIST_ND_KDE_STUDENTT (dnd);
  NcmStatsDistNdKDEStudenttPrivate * const self = dndt->priv;
  
  const gdouble lg_lnnorm   = lgamma (self->nu / 2.0) - lgamma ((self->nu + d) / 2.0);
  const gdouble chol_lnnorm = 0.5 * ncm_matrix_cholesky_lndet (cov_decomp);
  const gdouble nc_lnnorm   = (d / 2.0) * (ncm_c_lnpi () + log (self->nu));
  gdouble href_det          = 1.0;
  gint i;
  
  for (i = 0; i < d; i++)
    href_det *= ncm_vector_fast_get (href, i);
  
  return lg_lnnorm + nc_lnnorm + chol_lnnorm + log (href_det);
}

static void
_ncm_stats_dist_nd_kde_studentt_prepare_IM (NcmStatsDistNd *dnd, GPtrArray *Us, const gint d, const gint n, const NcmVector *href, NcmMatrix *IM)
{
  NcmStatsDistNdKDEStudentt *dndt               = NCM_STATS_DIST_ND_KDE_STUDENTT (dnd);
  NcmStatsDistNdKDEStudenttPrivate * const self = dndt->priv;
  gint i;
  
  for (i = 0; i < n; i++)
  {
    NcmVector *row_i = g_ptr_array_index (Us, i);
    gint j;
    
    ncm_matrix_set (IM, i, i, 1.0);
    
    for (j = i + 1; j < n; j++)
    {
      NcmVector *row_j = g_ptr_array_index (Us, j);
      gdouble m2lnp_ij = 0.0;
      gdouble p_ij;
      gint k;
      
      for (k = 0; k < d; k++)
      {
        m2lnp_ij += gsl_pow_2 ((ncm_vector_fast_get (row_i, k) - ncm_vector_fast_get (row_j, k)) /  ncm_vector_fast_get (href, k));
      }
      
      printf ("result %22.15g", m2lnp_ij);
      p_ij = _ncm_stats_dist_nd_kde_studentt_f (self, d, m2lnp_ij);
      ncm_matrix_set (IM, i, j, p_ij);
      ncm_matrix_set (IM, j, i, p_ij);
      /*printf("valor %22.15g na pos %d  e %d ", p_ij, i, j);*/
    }
  }
}

static gdouble
_ncm_stats_dist_nd_kde_studentt_eval (NcmStatsDistNd *dnd, NcmVector *weights, NcmVector *invUy, GPtrArray *invUsample, const gint d, const gint n, const NcmVector *href)
{
  NcmStatsDistNdKDEStudentt *dndt               = NCM_STATS_DIST_ND_KDE_STUDENTT (dnd);
  NcmStatsDistNdKDEStudenttPrivate * const self = dndt->priv;
  gdouble s                                     = 0.0;
  gdouble c                                     = 0.0;
  gint i;
  
  for (i = 0; i < n; i++)
  {
    NcmVector *row_i = g_ptr_array_index (invUsample, i);
    gdouble e_i, t, chi2_i = 0.0;
    gint k;
    
    for (k = 0; k < d; k++)
    {
      chi2_i += gsl_pow_2 ((ncm_vector_fast_get (row_i, k) - ncm_vector_fast_get (invUy, k)) / ncm_vector_fast_get (href, k));
    }
    
    e_i = ncm_vector_get (weights, i) * _ncm_stats_dist_nd_kde_studentt_f (self, d, chi2_i);
    
    t  = s + e_i;
    c += (s >= e_i) ? ((s - t) + e_i) : ((e_i - t) + s);
    s  = t;
  }
  
  return s;
}

static void
_ncm_stats_dist_nd_kde_studentt_kernel_sample (NcmStatsDistNd *dnd, NcmMatrix *cov_decomp, const guint d, NcmVector *y, NcmVector *mu, const NcmVector *href, NcmRNG *rng)
{
  NcmStatsDistNdKDEStudentt *dndt               = NCM_STATS_DIST_ND_KDE_STUDENTT (dnd);
  NcmStatsDistNdKDEStudenttPrivate * const self = dndt->priv;
  gdouble chi_scale;
  gint i, ret;
  
  for (i = 0; i < d; i++)
  {
    const gdouble u_i = gsl_ran_ugaussian (rng->r);
    
    ncm_vector_set (y, i, u_i * ncm_vector_fast_get (href, i));
  }
  
  /* CblasLower, CblasNoTrans => CblasUpper, CblasTrans */
  ret = gsl_blas_dtrmv (CblasUpper, CblasTrans, CblasNonUnit,
                        ncm_matrix_gsl (cov_decomp), ncm_vector_gsl (y));
  NCM_TEST_GSL_RESULT ("ncm_stats_dist_nd_kde_studentt_sample_mean_scale", ret);
  
  chi_scale = sqrt (self->nu / gsl_ran_chisq (rng->r, self->nu));
  
  ncm_vector_scale (y, chi_scale);
  ncm_vector_add (y, mu);
}

/**
 * ncm_stats_dist_nd_kde_studentt_new:
 * @dim: sample space dimension
 * @cv_type: a #NcmStatsDistNdCV
 * @nu: Student-t parameter $\nu$
 *
 * Creates a new #NcmStatsDistNdKDEStudentt object with sample dimension @dim.
 *
 * Returns: a new #NcmStatsDistNdKDEStudentt.
 */
NcmStatsDistNdKDEStudentt *
ncm_stats_dist_nd_kde_studentt_new (const guint dim, const NcmStatsDistNdCV cv_type, const gdouble nu)
{
  NcmStatsDistNdKDEStudentt *dndt = g_object_new (NCM_TYPE_STATS_DIST_ND_KDE_STUDENTT,
                                                  "dimension", dim,
                                                  "CV-type",   cv_type,
                                                  "nu",        nu,
                                                  NULL);
  
  return dndt;
}

/**
 * ncm_stats_dist_nd_kde_studentt_ref:
 * @dndt: a #NcmStatsDistNdKDEStudentt
 *
 * Increase the reference of @stats_dist_nd_kde_studentt by one.
 *
 * Returns: (transfer full): @stats_dist_nd_kde_studentt.
 */
NcmStatsDistNdKDEStudentt *
ncm_stats_dist_nd_kde_studentt_ref (NcmStatsDistNdKDEStudentt *dndt)
{
  return g_object_ref (dndt);
}

/**
 * ncm_stats_dist_nd_kde_studentt_free:
 * @dndt: a #NcmStatsDistNdKDEStudentt
 *
 * Decrease the reference count of @stats_dist_nd_kde_studentt by one.
 *
 */
void
ncm_stats_dist_nd_kde_studentt_free (NcmStatsDistNdKDEStudentt *dndt)
{
  g_object_unref (dndt);
}

/**
 * ncm_stats_dist_nd_kde_studentt_clear:
 * @dndt: a #NcmStatsDistNdKDEStudentt
 *
 * Decrease the reference count of @stats_dist_nd_kde_studentt by one, and sets the pointer *@stats_dist_nd_kde_studentt to
 * NULL.
 *
 */
void
ncm_stats_dist_nd_kde_studentt_clear (NcmStatsDistNdKDEStudentt **dndt)
{
  g_clear_object (dndt);
}

/**
 * ncm_stats_dist_nd_kde_studentt_set_nu:
 * @dndt: a #NcmStatsDistNdKDEStudentt
 * @nu: the over-smooth factor
 *
 * Sets the over-smooth factor to @nu.
 *
 */
void
ncm_stats_dist_nd_kde_studentt_set_nu (NcmStatsDistNdKDEStudentt *dndt, const gdouble nu)
{
  NcmStatsDistNdKDEStudenttPrivate * const self = dndt->priv;
  
  self->nu = nu;
}

/**
 * ncm_stats_dist_nd_kde_studentt_get_nu:
 * @dndt: a #NcmStatsDistNd
 *
 * Returns: the over-smooth factor.
 */
gdouble
ncm_stats_dist_nd_kde_studentt_get_nu (NcmStatsDistNdKDEStudentt *dndt)
{
  NcmStatsDistNdKDEStudenttPrivate * const self = dndt->priv;
  
  return self->nu;
}

