/***************************************************************************
 *            ncm_fftlog.c
 *
 *  Fri May 18 16:44:23 2012
 *  Copyright  2012  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <sandro@isoftware.com.br>
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
 * SECTION:ncm_fftlog
 * @title: Logarithm Fast Fourier Algorithm
 * @short_description: Object implementing logarithm fast fourier transform
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_fftlog.h"
#include "math/ncm_cfg.h"

#include <math.h>
#include <complex.h>
#ifdef NUMCOSMO_HAVE_FFTW3 
#include <fftw3.h>
#endif /* NUMCOSMO_HAVE_FFTW3 */
#ifndef HAVE_FFTW3_ALLOC
#define fftw_alloc_real(n) (double *) fftw_malloc(sizeof(double) * (n))
#define fftw_alloc_complex(n) (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * (n))
#endif /* HAVE_FFTW3_ALLOC */

enum
{
  PROP_0,
  PROP_R0,
  PROP_K0,
  PROP_L,
  PROP_N,
  PROP_NAME,
};

G_DEFINE_TYPE (NcmFftlog, ncm_fftlog, G_TYPE_OBJECT);

static void
ncm_fftlog_init (NcmFftlog *ncm_fftlog)
{
  ncm_fftlog->lnr0 = 0.0;
  ncm_fftlog->lnk0 = 0.0;
  ncm_fftlog->L    = 0.0;
  ncm_fftlog->N    = 0;
  ncm_fftlog->dr   = 0.0;
  ncm_fftlog->name = NULL;
}

static void
_ncm_fftlog_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmFftlog *fftlog = NCM_FFTLOG (object);
  g_return_if_fail (NCM_IS_FFTLOG (object));

  switch (prop_id)
  {
    case PROP_R0:
      fftlog->lnr0 = log (g_value_get_double (value));
      break;
    case PROP_K0:
      fftlog->lnk0 = log (g_value_get_double (value));
      break;
    case PROP_L:
      fftlog->L = g_value_get_double (value);
      break;
    case PROP_N:
      ncm_fftlog_set_size (fftlog, g_value_get_uint (value));
      break;
    case PROP_NAME:
      ncm_fftlog_set_name (fftlog, g_value_get_string (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_fftlog_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmFftlog *fftlog = NCM_FFTLOG (object);

  g_return_if_fail (NCM_IS_FFTLOG (object));

  switch (prop_id)
  {
    case PROP_R0:
      g_value_set_double (value, exp (fftlog->lnr0));
      break;
    case PROP_K0:
      g_value_set_double (value, exp (fftlog->lnk0));
      break;
    case PROP_L:
      g_value_set_double (value, fftlog->L);
      break;
    case PROP_N:
      g_value_set_uint (value, fftlog->N);
      break;
    case PROP_NAME:
      g_value_set_string (value, fftlog->name);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_fftlog_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (ncm_fftlog_parent_class)->constructed (object);
  {
#ifndef NUMCOSMO_HAVE_FFTW3
    g_error ("Cannot construct FFTLog object fftw3 library is not present, recompile numcosmo with fftw3 support.");
#endif /* NUMCOSMO_HAVE_FFTW3 */
  }
}

static void
ncm_fftlog_finalize (GObject *object)
{
#ifdef NUMCOSMO_HAVE_FFTW3
  NcmFftlog *fftlog = NCM_FFTLOG (object);
  if (fftlog->in != NULL)
    fftw_free (fftlog->in);
  if (fftlog->out != NULL)
    fftw_free (fftlog->out); 
  if (fftlog->p_in2out != NULL)
    fftw_destroy_plan (fftlog->p_in2out);
  if (fftlog->p_out2in != NULL)
    fftw_destroy_plan (fftlog->p_out2in);
#endif /* NUMCOSMO_HAVE_FFTW3 */
  
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_fftlog_parent_class)->finalize (object);
}

static void
ncm_fftlog_class_init (NcmFftlogClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);

  object_class->constructed  = &_ncm_fftlog_constructed;
  object_class->set_property = &_ncm_fftlog_set_property;
  object_class->get_property = &_ncm_fftlog_get_property;
  object_class->finalize     = &ncm_fftlog_finalize;

  g_object_class_install_property (object_class,
                                   PROP_R0,
                                   g_param_spec_double ("r0",
                                                        NULL,
                                                        "Center value for r",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_K0,
                                   g_param_spec_double ("k0",
                                                        NULL,
                                                        "Center value for k",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_L,
                                   g_param_spec_double ("L",
                                                        NULL,
                                                        "Function log-period",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 1.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_N,
                                   g_param_spec_uint ("N",
                                                      NULL,
                                                      "Number of knots",
                                                      0, G_MAXUINT, 10,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_NAME,
                                   g_param_spec_string ("name",
                                                      NULL,
                                                      "FFTW Plan wisdown name",
                                                      "fftlog_default_wisdown",
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

}

/**
 * ncm_fftlog_set_name:
 * @fftlog: FIXME
 * @name: FIXME
 * 
 * FIXME
 * 
 */
void
ncm_fftlog_set_name (NcmFftlog *fftlog, const gchar *name)
{
  g_free (fftlog->name);
  fftlog->name = g_strdup_printf ("%s.fftw3", name);
}

/**
 * ncm_fftlog_peek_name:
 * @fftlog: FIXME
 * 
 * FIXME
 * 
 * Returns: (transfer none): FIXME
 */
gchar *
ncm_fftlog_peek_name (NcmFftlog *fftlog)
{
  return fftlog->name;
}

/**
 * ncm_fftlog_set_size:
 * @fftlog: FIXME
 * @n: FIXME
 * 
 * FIXME
 * 
 */
void
ncm_fftlog_set_size (NcmFftlog *fftlog, guint n)
{
#ifdef NUMCOSMO_HAVE_FFTW3
  fftlog->dr  = fftlog->L / (1.0 * fftlog->N);
  if (n != fftlog->N)
  {
    fftlog->N_2 = fftlog->N % 2 == 1 ? (fftlog->N - 1) / 2 : fftlog->N / 2;
    
    if (fftlog->in != NULL)
      fftw_free (fftlog->in);
    if (fftlog->out != NULL)
      fftw_free (fftlog->out); 
    if (fftlog->p_in2out != NULL)
      fftw_destroy_plan (fftlog->p_in2out);
    if (fftlog->p_out2in != NULL)
      fftw_destroy_plan (fftlog->p_out2in);
      
    fftlog->in  = fftw_alloc_complex (fftlog->N);
    fftlog->out = fftw_alloc_complex (fftlog->N);

    ncm_cfg_load_fftw_wisdom (fftlog->name);
    fftlog->p_in2out = fftw_plan_dft_1d (fftlog->N, fftlog->in, fftlog->out, FFTW_FORWARD, FFTW_PATIENT | FFTW_DESTROY_INPUT);
    fftlog->p_out2in = fftw_plan_dft_1d (fftlog->N, fftlog->out, fftlog->in, FFTW_FORWARD, FFTW_PATIENT | FFTW_DESTROY_INPUT);
    ncm_cfg_save_fftw_wisdom (fftlog->name);
  }
#endif /* NUMCOSMO_HAVE_FFTW3 */
}

/**
 * ncm_fftlog_get_size:
 * @fftlog: FIXME
 * 
 * FIXME
 * 
 * Returns: FIXME
 */
guint
ncm_fftlog_get_size (NcmFftlog *fftlog)
{
  return fftlog->N;
}

/**
 * ncm_fftlog_set_length:
 * @fftlog: FIXME
 * @L: FIXME
 * 
 * FIXME
 * 
 */
void
ncm_fftlog_set_length (NcmFftlog *fftlog, gdouble L)
{
  fftlog->L = L;
  fftlog->dr = fftlog->L / (1.0 * fftlog->N);
}

/**
 * ncm_fftlog_get_length:
 * @fftlog: FIXME
 * 
 * FIXME
 * 
 * Returns: FIXME
 */
gdouble
ncm_fftlog_get_length (NcmFftlog *fftlog)
{
  return fftlog->L;
}

