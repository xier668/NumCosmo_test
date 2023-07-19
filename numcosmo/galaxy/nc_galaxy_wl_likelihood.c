/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 2; tab-width: 2 -*-  */

/***************************************************************************
 *            nc_galaxy_wl_likelihood.c
 *
 *  Mon May 08 16:12:03 2023
 *  Copyright  2023  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_galaxy_wl_likelihood.c
 * Copyright (C) 2023 Caio Lima de Oliveira <caiolimadeoliveira@pm.me>
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
 * SECTION: nc_galaxy_wl_likelihood
 * @title: NcGalaxyWLLikelihood
 * @short_description: Class describing galaxy weak lensing distributions.
 * @stability: Unstable
 *
 *
 * This class describes a galaxy weak lensing distribution.
 * It is composed by three distributions: a shape distribution $P(s)$, a proxy redshift distribution $P(z_p)$, and a position distribution $P(z)P(r)$.
 * The shape distribution is defined by the abstract class #NcGalaxySDShape.
 * The proxy redshift distribution is defined by the abstract class #NcGalaxySDZProxy.
 * The position distribution is defined by the abstract class #NcGalaxySDPosition.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "galaxy/nc_galaxy_wl_likelihood.h"
#include "galaxy/nc_galaxy_sd_shape.h"
#include "galaxy/nc_galaxy_sd_z_proxy.h"
#include "galaxy/nc_galaxy_sd_position.h"
#include "math/ncm_stats_dist1d.h"
#include "math/ncm_stats_dist1d_epdf.h"
#include <math.h>
#include <gsl/gsl_math.h>

struct _NcGalaxyWLLikelihoodPrivate
{
  NcmMatrix *obs;
  NcmMatrix *samples;
  NcGalaxySDZProxy *zp_dist;
  NcGalaxySDPosition *rz_dist;
  NcmStatsDist1dEPDF *kde;
  gdouble cut_fraction;
  gdouble zp_min;
  gdouble zp_max;
  gdouble sigma;
  gdouble scale_cut;
  gint ndata;
  gboolean constructed;
  guint len;
};

enum
{
  PROP_0,
  PROP_OBS,
  PROP_ZP_DIST,
  PROP_RZ_DIST,
  PROP_KDE,
  PROP_ZP_MIN,
  PROP_ZP_MAX,
  PROP_SIGMA,
  PROP_SCALE_CUT,
  PROP_NDATA,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcGalaxyWLLikelihood, nc_galaxy_wl_likelihood, G_TYPE_OBJECT);

static void
nc_galaxy_wl_likelihood_init (NcGalaxyWLLikelihood *gwl)
{
  NcGalaxyWLLikelihoodPrivate * const self = gwl->priv = nc_galaxy_wl_likelihood_get_instance_private (gwl);

  self->obs         = NULL;
  self->zp_dist     = NULL;
  self->rz_dist     = NULL;
  self->kde         = ncm_stats_dist1d_epdf_new_full (20000, NCM_STATS_DIST1D_EPDF_BW_RoT, 0.1, 1.0e-2);
  self->len         = 0;
  self->zp_max      = 0.0;
  self->zp_min      = 0.0;
  self->sigma       = 0.0;
  self->scale_cut   = 0.0;
  self->ndata       = 0;
  self->constructed = FALSE;
  self->samples     = NULL;
}

static void
_nc_galaxy_wl_likelihood_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcGalaxyWLLikelihood *gwl                = NC_GALAXY_WL_LIKELIHOOD (object);
  NcGalaxyWLLikelihoodPrivate * const self = gwl->priv;

  g_return_if_fail (NC_IS_GALAXY_WL_LIKELIHOOD (object));

  switch (prop_id)
  {
    case PROP_OBS:
      nc_galaxy_wl_likelihood_set_obs (gwl, g_value_get_object (value));
      break;
    case PROP_ZP_DIST:
      self->zp_dist = g_value_dup_object (value);
      break;
    case PROP_RZ_DIST:
      self->rz_dist = g_value_dup_object (value);
      break;
    case PROP_ZP_MIN:
      self->zp_min = g_value_get_double (value);

      if (self->constructed)
        g_assert_cmpfloat (self->zp_min, <, self->zp_max);

      break;
    case PROP_ZP_MAX:
      self->zp_max = g_value_get_double (value);

      if (self->constructed)
        g_assert_cmpfloat (self->zp_min, <, self->zp_max);

      break;
    case PROP_SIGMA:
      nc_galaxy_wl_likelihood_set_sigma (gwl, g_value_get_double (value));
      break;
    case PROP_SCALE_CUT:
      nc_galaxy_wl_likelihood_set_scale_cut (gwl, g_value_get_double (value));
      break;
    case PROP_NDATA:
      nc_galaxy_wl_likelihood_set_ndata (gwl, g_value_get_int (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_galaxy_wl_likelihood_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcGalaxyWLLikelihood *gwl                = NC_GALAXY_WL_LIKELIHOOD (object);
  NcGalaxyWLLikelihoodPrivate * const self = gwl->priv;

  g_return_if_fail (NC_IS_GALAXY_WL_LIKELIHOOD (object));

  switch (prop_id)
  {
    case PROP_OBS:
      g_value_set_object (value, nc_galaxy_wl_likelihood_peek_obs (gwl));
      break;
    case PROP_ZP_DIST:
      g_value_set_object (value, self->zp_dist);
      break;
    case PROP_RZ_DIST:
      g_value_set_object (value, self->rz_dist);
      break;
    case PROP_KDE:
      g_value_set_object (value, nc_galaxy_wl_likelihood_peek_kde (gwl));
      break;
    case PROP_ZP_MIN:
      g_value_set_double (value, self->zp_min);
      break;
    case PROP_ZP_MAX:
      g_value_set_double (value, self->zp_max);
      break;
    case PROP_SIGMA:
      g_value_set_double (value, self->sigma);
    case PROP_SCALE_CUT:
      g_value_set_double (value, self->scale_cut);
      break;
    case PROP_NDATA:
      g_value_set_int (value, self->ndata);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_galaxy_wl_likelihood_dispose (GObject *object)
{
  NcGalaxyWLLikelihood *gwl                = NC_GALAXY_WL_LIKELIHOOD (object);
  NcGalaxyWLLikelihoodPrivate * const self = gwl->priv;

  ncm_matrix_clear (&self->obs);
  nc_galaxy_sd_z_proxy_clear (&self->zp_dist);
  nc_galaxy_sd_position_clear (&self->rz_dist);

  ncm_stats_dist1d_epdf_clear (&self->kde);

  G_OBJECT_CLASS (nc_galaxy_wl_likelihood_parent_class)->dispose (object);
}

static void
_nc_galaxy_wl_likelihood_finalize (GObject *object)
{
  G_OBJECT_CLASS (nc_galaxy_wl_likelihood_parent_class)->finalize (object);
}

static void
_nc_galaxy_wl_likelihood_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (nc_galaxy_wl_likelihood_parent_class)->constructed (object);
  {
    NcGalaxyWLLikelihood *gwl                = NC_GALAXY_WL_LIKELIHOOD (object);
    NcGalaxyWLLikelihoodPrivate * const self = gwl->priv;

    g_assert_cmpfloat (self->zp_min, <, self->zp_max);

    self->constructed = TRUE;
  }
}

static void
nc_galaxy_wl_likelihood_class_init (NcGalaxyWLLikelihoodClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &_nc_galaxy_wl_likelihood_set_property;
  object_class->get_property = &_nc_galaxy_wl_likelihood_get_property;
  object_class->dispose      = &_nc_galaxy_wl_likelihood_dispose;
  object_class->finalize     = &_nc_galaxy_wl_likelihood_finalize;
  object_class->constructed  = &_nc_galaxy_wl_likelihood_constructed;

  /**
   * NcGalaxyWLLikelihood:obs:
   *
   * Galaxy weak lensing observables matrix.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_OBS,
                                   g_param_spec_object ("obs",
                                                        NULL,
                                                        "Galaxy weak lensing observables",
                                                        NCM_TYPE_MATRIX,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcGalaxyWLLikelihood:zp-dist:
   *
   * A #NcGalaxySDZProxy object.
   *
   */

  g_object_class_install_property (object_class,
                                   PROP_ZP_DIST,
                                   g_param_spec_object ("zp-dist",
                                                        NULL,
                                                        "Galaxy sample proxy redshift distribution",
                                                        NC_TYPE_GALAXY_SD_Z_PROXY,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcGalaxyWLLikelihood:rz-dist:
   *
   * A #NcGalaxySDZPosition object.
   *
   */

  g_object_class_install_property (object_class,
                                   PROP_RZ_DIST,
                                   g_param_spec_object ("rz-dist",
                                                        NULL,
                                                        "Galaxy sample position distribution",
                                                        NC_TYPE_GALAXY_SD_POSITION,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcGalaxyWLLikelihood:zp-min:
   *
   * Minimum redshift of the weak lensing observables.
   *
   */

  g_object_class_install_property (object_class,
                                   PROP_ZP_MIN,
                                   g_param_spec_double ("zp-min",
                                                        NULL,
                                                        "Minimum redshift of the weak lensing observables",
                                                        0.0, G_MAXDOUBLE, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcGalaxyWLLikelihood:zp-max:
   *
   * Maximum redshift of the weak lensing observables.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_ZP_MAX,
                                   g_param_spec_double ("zp-max",
                                                        NULL,
                                                        "Maximum redshift of the weak lensing observables",
                                                        0.0, G_MAXDOUBLE, 10.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcGalaxyWLLikelihood:sigma:
   *
   * Maximum redshift of the weak lensing observables.
   *
   */

  g_object_class_install_property (object_class,
                                   PROP_SIGMA,
                                   g_param_spec_double ("sigma",
                                                        NULL,
                                                        "Maximum redshift of the weak lensing observables",
                                                        0.0, G_MAXDOUBLE, 10.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcGalaxyWLLikelihood:scale-cut:
   *
   * Scale of probability at which to include a sample in the cut.
   *
   */

  g_object_class_install_property (object_class,
                                   PROP_SCALE_CUT,
                                   g_param_spec_double ("scale-cut",
                                                        NULL,
                                                        "Scale of probability at which to include a sample in the cut",
                                                        0.0, 1.0, 0.00001,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcGalaxyWLLikelihood:ndata:
   *
   * Number of data points to sample for KDE.
   *
   */

  g_object_class_install_property (object_class,
                                   PROP_NDATA,
                                   g_param_spec_int ("ndata",
                                                     NULL,
                                                     "Number of data points to sample for KDE",
                                                     0.0, G_MAXINT, 10000.0,
                                                     G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

/**
 * nc_galaxy_wl_likelihood_new:
 * @s_dist: a #NcGalaxySDShape
 * @zp_dist: a #NcGalaxySDZProxy
 * @rz_dist: a #NcGalaxySDPosition
 *
 * Creates a new galaxy weak lensing object.
 * Requires an instance of #NcGalaxySDShape, #NcGalaxySDZProxy, and #NcGalaxySDPosition.
 *
 * Returns: (transfer full): a new NcGalaxyWLLikelihood.
 */
NcGalaxyWLLikelihood *
nc_galaxy_wl_likelihood_new (NcGalaxySDShape *s_dist, NcGalaxySDZProxy *zp_dist, NcGalaxySDPosition *rz_dist)
{
  NcGalaxyWLLikelihood *gwl = g_object_new (NC_TYPE_GALAXY_WL_LIKELIHOOD,
                                            "zp-dist", zp_dist,
                                            "rz-dist", rz_dist,
                                            NULL);

  return gwl;
}

/**
 * nc_galaxy_wl_likelihood_ref:
 * @gwl: a #NcGalaxyWLLikelihood
 *
 * Increase the reference of @gwl by one.
 *
 * Returns: (transfer full): @gwl.
 */
NcGalaxyWLLikelihood *
nc_galaxy_wl_likelihood_ref (NcGalaxyWLLikelihood *gwl)
{
  return g_object_ref (gwl);
}

/**
 * nc_galaxy_wl_likelihood_free:
 * @gwl: a #NcGalaxyWLLikelihood
 *
 * Decrease the reference count of @gwl by one.
 *
 */
void
nc_galaxy_wl_likelihood_free (NcGalaxyWLLikelihood *gwl)
{
  g_object_unref (gwl);
}

/**
 * nc_galaxy_wl_likelihood_clear:
 * @gwl: a #NcGalaxyWLLikelihood
 *
 * Decrease the reference count of @gwl by one, and sets the pointer *@gwl to
 * NULL.
 *
 */
void
nc_galaxy_wl_likelihood_clear (NcGalaxyWLLikelihood **gwl)
{
  g_clear_object (gwl);
}

/**
 * nc_galaxy_wl_likelihood_set_obs:
 * @gwl: a #NcGalaxyWLLikelihood
 * @obs: a #NcmMatrix
 *
 * Sets the observables matrix @obs.
 */
void
nc_galaxy_wl_likelihood_set_obs (NcGalaxyWLLikelihood *gwl, NcmMatrix *obs)
{
  NcGalaxyWLLikelihoodPrivate * const self = gwl->priv;

  g_assert_cmpuint (ncm_matrix_ncols (obs), ==, 4);
  g_assert_cmpuint (ncm_matrix_nrows (obs), >, 0);

  ncm_matrix_clear (&self->obs);

  self->len = ncm_matrix_nrows (obs);
  self->obs = ncm_matrix_ref (obs);
}

/**
 * nc_galaxy_wl_likelihood_peek_obs:
 * @gwl: a #NcGalaxyWLLikelihood
 *
 * Gets the observables matrix.
 *
 * Returns: (transfer none): the observables matrix.
 */
NcmMatrix *
nc_galaxy_wl_likelihood_peek_obs (NcGalaxyWLLikelihood *gwl)
{
  NcGalaxyWLLikelihoodPrivate * const self = gwl->priv;

  return self->obs;
}

/**
 * nc_galaxy_wl_likelihood_peek_kde:
 * @gwl: a #NcGalaxyWLLikelihood
 *
 * Gets the observables matrix.
 *
 * Returns: (transfer none): the observables matrix.
 */
NcmStatsDist1dEPDF *
nc_galaxy_wl_likelihood_peek_kde (NcGalaxyWLLikelihood *gwl)
{
  NcGalaxyWLLikelihoodPrivate * const self = gwl->priv;

  return self->kde;
}

void
nc_galaxy_wl_likelihood_prepare (NcGalaxyWLLikelihood *gwl, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, const gdouble z_cluster)
{
  NcGalaxyWLLikelihoodPrivate * const self = gwl->priv;
  NcmRNG *rng                              = ncm_rng_new (NULL);
  gdouble in_cut                             = 0;
  gdouble out_cut                            = 0;
  gdouble avg                              = 0.0;
  gdouble sigma                            = 0.0;
  gint i                                   = 0;

  ncm_stats_dist1d_epdf_reset (self->kde);

  while (i < self->ndata)
  {
    gdouble z;
    gdouble zp;

    while (TRUE)
    {
      nc_galaxy_sd_position_gen_z (self->rz_dist, rng, &z);

      if (nc_galaxy_sd_z_proxy_gen (self->zp_dist, rng, z, &zp))
        break;
    }

    {
      if ((zp_i >= self->zp_min) && (self->zp_max >= zp_i))
      {
        i++;
        in_cut++;
      }
      else
      {
        out_cut++;
      }

      ncm_matrix_set (self->samples, i, 0, zp);
      ncm_matrix_set (self->samples, i, 1, z);

      ncm_stats_dist1d_epdf_add_obs (self->kde, zp);
      avg = avg_i;
      sigma = sigma_i;
    }
  }

  self->cut_fraction = in_cut / (in_cut + out_cut);

  // printf ("%f\n", self->cut_fraction);

  ncm_stats_dist1d_prepare (NCM_STATS_DIST1D (self->kde));

  ncm_rng_clear (&rng);
}

/**
 * nc_galaxy_wl_likelihood_eval_m2lnP:
 * @gwl: a #NcGalaxyWLLikelihood
 * @cosmo: a #NcHICosmo
 * @dp: a #NcHaloDensityProfile
 * @smd: a #NcWLSurfaceMassDensity
 * @z_cluster: cluster redshift $z_\mathrm{cl}$
 *
 * Computes the observables probability given the theoretical modeling using
 * integration method.
 *
 *
 * Returns: $-2\ln(P)$.
 */
gdouble
nc_galaxy_wl_likelihood_eval_m2lnP (NcGalaxyWLLikelihood *gwl, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, const gdouble z_cluster)
{
  return 0.0;
}

/**
 * nc_galaxy_wl_likelihood_kde_eval_m2lnP:
 * @gwl: a #NcGalaxyWLLikelihood
 * @cosmo: a #NcHICosmo
 * @dp: a #NcHaloDensityProfile
 * @smd: a #NcWLSurfaceMassDensity
 * @z_cluster: cluster redshift $z_\mathrm{cl}$
 *
 * Computes the observables probability given the theoretical modeling using
 * kernel density estimation method.
 *
 *
 * Returns: $-2\ln(P)$.
 */
gdouble
nc_galaxy_wl_likelihood_kde_eval_m2lnP (NcGalaxyWLLikelihood *gwl, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, const gdouble z_cluster)
{
  NcGalaxyWLLikelihoodPrivate * const self = gwl->priv;
  gdouble res                              = 0.0;
  gint gal_i;

  nc_galaxy_wl_likelihood_prepare (gwl, cosmo, dp, smd, z_cluster);

  for (gal_i = 0; gal_i < self->len; gal_i++)
  {
    const gdouble r_i  = ncm_matrix_get (self->obs, gal_i, 0);
    const gdouble zp_i = ncm_matrix_get (self->obs, gal_i, 1);
    const gdouble et_i = ncm_matrix_get (self->obs, gal_i, 2);
    const gdouble ex_i = ncm_matrix_get (self->obs, gal_i, 3);

    if ((zp_i >= self->zp_min) && (self->zp_max >= zp_i))
    {
      gdouble p = ncm_stats_dist1d_eval_p (NCM_STATS_DIST1D (self->kde), zp_i);

      gint closest = 0;
      gint j;

      for (j = 0; j < self->ndata; j++)
      {
        if (abs (zp_i - ncm_matrix_get (self->samples, j, 0)) < abs (zp_i - ncm_matrix_get (self->samples, closest, 0)))
          closest = j;
      }

      gdouble z_i = ncm_matrix_get (self->samples, closest, 1);

      complex double e_obs    = et_i + I * ex_i;
      complex double shear    = nc_wl_surface_mass_density_reduced_shear (smd, dp, cosmo, r_i, z_i, z_cluster, z_cluster);
      complex double e_source = (e_obs - shear) / (1.0 - conj (shear) * e_obs);

      gdouble et = creal (e_source);
      gdouble ex = cimag (e_source);

      res -= 2 * log (p);
      res += (et * et + ex * ex) / (self->sigma * self->sigma);
      res += 4 * log (sqrt (2 * M_PI) * self->sigma);
    }
  }

  return res;
}

/**
 * nc_galaxy_wl_likelihood_len:
 * @gwll: a #NcGalaxyWL
 *
 * Returns: the number of galaxies in @gwl.
 */
guint
nc_galaxy_wl_likelihood_len (NcGalaxyWLLikelihood *gwll)
{
  NcGalaxyWLLikelihoodPrivate * const self = gwll->priv;

  return self->len;
}

/**
 * nc_galaxy_wl_likelihood_set_cut:
 * @gwl: a #NcGalaxyWL
 * @zp_min: minimum redshift proxy $z_\mathrm{p,min}$
 * @zp_max: maximum redshift proxy $z_\mathrm{p,max}$
 * @r_min: minimum projected radius $r_\mathrm{min}$
 * @r_max: maximum projected radius $r_\mathrm{max}$
 *
 * Sets the cut in the observables.
 *
 */
void
nc_galaxy_wl_likelihood_set_cut (NcGalaxyWLLikelihood *gwl, const gdouble zp_min, const gdouble zp_max)
{
  NcGalaxyWLLikelihoodPrivate * const self = gwl->priv;

  g_assert_cmpfloat (zp_min, <, zp_max);

  self->zp_min = zp_min;
  self->zp_max = zp_max;
}

/**
 * nc_galaxy_wl_likelihood_set_ndata:
 * @gwl: a #NcGalaxyWL
 * @ndata: number of samples to take for KDE
 *
 * Sets the number of samples ndata.
 *
 */
void
nc_galaxy_wl_likelihood_set_ndata (NcGalaxyWLLikelihood *gwl, gint ndata)
{
  NcGalaxyWLLikelihoodPrivate * const self = gwl->priv;

  ncm_matrix_clear (&self->samples);

  self->ndata   = ndata;
  self->samples = ncm_matrix_new (ndata, 2);
}

/**
 * nc_galaxy_wl_likelihood_set_scale_cut:
 * @gwl: a #NcGalaxyWL
 * @scale: scale of the probability at which a point should be included in the cut
 *
 * Sets the scale of the probability at which a point should be included in the cut.
 *
 */
void
nc_galaxy_wl_likelihood_set_scale_cut (NcGalaxyWLLikelihood *gwl, const gdouble scale)
{
  NcGalaxyWLLikelihoodPrivate * const self = gwl->priv;

  self->scale_cut = scale;
}

/**
 * nc_galaxy_wl_likelihood_set_sigma:
 * @gwl: a #NcGalaxyWL
 * @sigma: standard deviation of ellipticity distribution
 *
 * Sets the standard deviation.
 *
 */
void
nc_galaxy_wl_likelihood_set_sigma (NcGalaxyWLLikelihood *gwl, gdouble sigma)
{
  NcGalaxyWLLikelihoodPrivate * const self = gwl->priv;

  self->sigma = sigma;
}

