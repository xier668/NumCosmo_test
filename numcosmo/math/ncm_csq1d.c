/* -*- Mode: C; indent-tabs-mode: nil; c-basic-offset: 2; tab-width: 2 -*-  */
/***************************************************************************
 *            ncm_csq1d.c
 *
 *  Mon September 09 13:56:19 2019
 *  Copyright  2019  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * ncm_csq1d.c
 * Copyright (C) 2019 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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
 * SECTION:ncm_csq1d
 * @title: NcmCSQ1D
 * @short_description: Abstract class for Harmonic Oscillator calculation through complex structure quantization.
 *
 * 
 * The system:
 * \begin{align}
 * q^\prime &= \frac{\Pi_q}{m},
 * \Pi_q^\prime &= m\nu^2q.
 * \end{align}
 * 
 * \begin{equation}
 * \xi = \ln (m\nu)
 * \end{equation}
 * 
 * \begin{equation}
 * F^n = \left(\frac{1}{2\nu}\frac{\partial}{\partial t}\right)^n \xi. 
 * \end{equation}
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include <complex.h>

#include "math/ncm_csq1d.h"
#include "math/ncm_model_ctrl.h"
#include "math/ncm_spline_cubic_notaknot.h"
#include "math/ncm_diff.h"
#include "math/ncm_c.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_roots.h>

#include <nvector/nvector_serial.h>

#include <cvode/cvode.h>
#include <cvode/cvode_ls.h>
#include <arkode/arkode.h>
#include <arkode/arkode_arkstep.h>
#include <arkode/arkode_ls.h>

#include <sundials/sundials_matrix.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>

#include <sundials/sundials_types.h> 
#endif /* NUMCOSMO_GIR_SCAN */

#define PRINT_EVOL TRUE

typedef enum _NcmCSQ1DEvolStop
{
  NCM_CSQ1D_EVOL_STOP_ERROR = -1,
  NCM_CSQ1D_EVOL_STOP_FINISHED = 0,
  NCM_CSQ1D_EVOL_STOP_ADIABATIC_START,
  NCM_CSQ1D_EVOL_STOP_UP_START,
  NCM_CSQ1D_EVOL_STOP_UM_START,
} NcmCSQ1DEvolStop;

typedef enum _NcmCSQ1DEvolState
{
  NCM_CSQ1D_EVOL_STATE_INVALID = 0,
  NCM_CSQ1D_EVOL_STATE_ADIABATIC,
  NCM_CSQ1D_EVOL_STATE_UP,
  NCM_CSQ1D_EVOL_STATE_UM,
} NcmCSQ1DEvolState;

struct _NcmCSQ1DPrivate
{
  gdouble reltol;
  gdouble abstol;
  gdouble k;
  gdouble ti;
  gdouble tf;
  gdouble t;
  NcmCSQ1DEvolState state;
  gdouble adiab_threshold;
  gboolean save_evol;
  NcmModelCtrl *ctrl;
  gpointer cvode;
  gpointer cvode_Up;
  gpointer cvode_Um;
  gpointer arkode;
  gboolean cvode_init;
  gboolean cvode_Up_init;
  gboolean cvode_Um_init;
  gboolean arkode_init;
  N_Vector y;
  N_Vector y_Up;
  N_Vector y_Um;
  SUNMatrix A;
  SUNMatrix A_Up;
  SUNMatrix A_Um;
  SUNLinearSolver LS;
  SUNLinearSolver LS_Up;
  SUNLinearSolver LS_Um;
  NcmSpline *alpha_s;
  NcmSpline *dgamma_s;
  NcmDiff *diff;
};

enum
{
  PROP_0,
  PROP_RELTOL,
  PROP_ABSTOL,
  PROP_K,
  PROP_TI,
  PROP_TF,
  PROP_ADIAB_THRESHOLD,
  PROP_SAVE_EVOL,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcmCSQ1D, ncm_csq1d, G_TYPE_OBJECT);

static void
ncm_csq1d_init (NcmCSQ1D *csq1d)
{
  NcmCSQ1DPrivate * const self = csq1d->priv = ncm_csq1d_get_instance_private (csq1d);

  self->reltol          = 0.0;
  self->abstol          = 0.0;
  self->k               = 0.0;
  self->ti              = 0.0;
  self->tf              = 0.0;
  self->state           = NCM_CSQ1D_EVOL_STATE_INVALID;
  self->adiab_threshold = 0.0;
  self->save_evol       = FALSE;
  self->ctrl            = ncm_model_ctrl_new (NULL);

  self->cvode           = NULL;
  self->cvode_init      = FALSE;
  self->cvode_Up        = NULL;
  self->cvode_Up_init   = FALSE;
  self->cvode_Um        = NULL;
  self->cvode_Um_init   = FALSE;
  self->arkode          = NULL;
  self->arkode_init     = FALSE;

  self->y               = N_VNew_Serial (2);
  self->y_Up            = N_VNew_Serial (2);
  self->y_Um            = N_VNew_Serial (2);
  
  self->A               = SUNDenseMatrix (2, 2);
  NCM_CVODE_CHECK ((gpointer)self->A, "SUNDenseMatrix", 0, );  

  self->A_Up           = SUNDenseMatrix (2, 2);
  NCM_CVODE_CHECK ((gpointer)self->A_Up, "SUNDenseMatrix", 0, );  

  self->A_Um           = SUNDenseMatrix (2, 2);
  NCM_CVODE_CHECK ((gpointer)self->A_Um, "SUNDenseMatrix", 0, );  

  self->LS              = SUNDenseLinearSolver (self->y, self->A);
  NCM_CVODE_CHECK ((gpointer)self->LS, "SUNDenseLinearSolver", 0, );

  self->LS_Up          = SUNDenseLinearSolver (self->y_Up, self->A_Up);
  NCM_CVODE_CHECK ((gpointer)self->LS_Up, "SUNDenseLinearSolver", 0, );

  self->LS_Um          = SUNDenseLinearSolver (self->y_Um, self->A_Um);
  NCM_CVODE_CHECK ((gpointer)self->LS_Um, "SUNDenseLinearSolver", 0, );

  self->alpha_s         = ncm_spline_cubic_notaknot_new ();
  self->dgamma_s        = ncm_spline_cubic_notaknot_new ();

  self->diff            = ncm_diff_new ();
}

static void
_ncm_csq1d_dispose (GObject *object)
{
  NcmCSQ1D *csq1d = NCM_CSQ1D (object);
  NcmCSQ1DPrivate * const self = csq1d->priv;
  
  ncm_model_ctrl_clear (&self->ctrl);
  ncm_spline_clear (&self->alpha_s);
  ncm_spline_clear (&self->dgamma_s);

  ncm_diff_clear (&self->diff);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_csq1d_parent_class)->dispose (object);
}

static void
_ncm_csq1d_finalize (GObject *object)
{
  NcmCSQ1D *csq1d = NCM_CSQ1D (object);
  NcmCSQ1DPrivate * const self = csq1d->priv;

  if (self->cvode != NULL)
  {
    CVodeFree (&self->cvode);
    self->cvode      = NULL;
    self->cvode_init = FALSE;
  }
  if (self->cvode_Up != NULL)
  {
    CVodeFree (&self->cvode_Up);
    self->cvode_Up      = NULL;
    self->cvode_Up_init = FALSE;
  }
  if (self->cvode_Um != NULL)
  {
    CVodeFree (&self->cvode_Um);
    self->cvode_Um      = NULL;
    self->cvode_Um_init = FALSE;
  }
  if (self->arkode != NULL)
  {
    ARKStepFree (&self->arkode);
    self->arkode      = NULL;
    self->arkode_init = FALSE;
  }

  g_clear_pointer (&self->y, N_VDestroy);
  g_clear_pointer (&self->y_Up, N_VDestroy);
  g_clear_pointer (&self->y_Um, N_VDestroy);

  if (self->A != NULL)
    SUNMatDestroy (self->A);

  if (self->A_Up != NULL)
    SUNMatDestroy (self->A_Up);

  if (self->A_Um != NULL)
    SUNMatDestroy (self->A_Um);

  if (self->LS != NULL)
  {
    gint flag = SUNLinSolFree (self->LS);
    NCM_CVODE_CHECK (&flag, "SUNLinSolFree", 1, );
  }

  if (self->LS_Up != NULL)
  {
    gint flag = SUNLinSolFree (self->LS_Up);
    NCM_CVODE_CHECK (&flag, "SUNLinSolFree", 1, );
  }

  if (self->LS_Um != NULL)
  {
    gint flag = SUNLinSolFree (self->LS_Um);
    NCM_CVODE_CHECK (&flag, "SUNLinSolFree", 1, );
  }

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_csq1d_parent_class)->finalize (object);
}

static void
_ncm_csq1d_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmCSQ1D *csq1d = NCM_CSQ1D (object);
  g_return_if_fail (NCM_IS_CSQ1D (object));

  switch (prop_id)
  {
    case PROP_RELTOL:
      ncm_csq1d_set_reltol (csq1d, g_value_get_double (value));
      break;
    case PROP_ABSTOL:
      ncm_csq1d_set_abstol (csq1d, g_value_get_double (value));
      break;
    case PROP_K:
      ncm_csq1d_set_k (csq1d, g_value_get_double (value));
      break;
    case PROP_TI:
      ncm_csq1d_set_ti (csq1d, g_value_get_double (value));
      break;
    case PROP_TF:
      ncm_csq1d_set_tf (csq1d, g_value_get_double (value));
      break;
    case PROP_ADIAB_THRESHOLD:
      ncm_csq1d_set_adiab_threshold (csq1d, g_value_get_double (value));
      break;
    case PROP_SAVE_EVOL:
      ncm_csq1d_set_save_evol (csq1d, g_value_get_boolean (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_csq1d_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmCSQ1D *csq1d = NCM_CSQ1D (object);
  g_return_if_fail (NCM_IS_CSQ1D (object));

  switch (prop_id)
  {
    case PROP_RELTOL:
      g_value_set_double (value, ncm_csq1d_get_reltol (csq1d));
      break;
    case PROP_ABSTOL:
      g_value_set_double (value, ncm_csq1d_get_abstol (csq1d));
      break;
    case PROP_K:
      g_value_set_double (value, ncm_csq1d_get_k (csq1d));
      break;
    case PROP_TI:
      g_value_set_double (value, ncm_csq1d_get_ti (csq1d));
      break;
    case PROP_TF:
      g_value_set_double (value, ncm_csq1d_get_tf (csq1d));
      break;
    case PROP_ADIAB_THRESHOLD:
      g_value_set_double (value, ncm_csq1d_get_adiab_threshold (csq1d));
      break;
    case PROP_SAVE_EVOL:
      g_value_set_boolean (value, ncm_csq1d_get_save_evol (csq1d));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static gdouble _ncm_csq1d_eval_xi (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k);
static gdouble _ncm_csq1d_eval_nu (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k);
static gdouble _ncm_csq1d_eval_F1 (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k);
static gdouble _ncm_csq1d_eval_F2 (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k);
static gdouble _ncm_csq1d_eval_FN (NcmCSQ1D *csq1d, NcmModel *model, const gint n, const gdouble t, const gdouble k);
static gdouble _ncm_csq1d_eval_powspec_factor (NcmCSQ1D *csq1d, NcmModel *model, const gdouble k);

static void
ncm_csq1d_class_init (NcmCSQ1DClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &_ncm_csq1d_set_property;
  object_class->get_property = &_ncm_csq1d_get_property;
  object_class->dispose      = &_ncm_csq1d_dispose;
  object_class->finalize     = &_ncm_csq1d_finalize;

  g_object_class_install_property (object_class,
                                   PROP_RELTOL,
                                   g_param_spec_double ("reltol",
                                                        NULL,
                                                        "Relative tolerance",
                                                        0.0, 1.0, NCM_DEFAULT_PRECISION,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_ABSTOL,
                                   g_param_spec_double ("abstol",
                                                        NULL,
                                                        "Absolute tolerance tolerance",
                                                        0.0, G_MAXDOUBLE, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));  
  g_object_class_install_property (object_class,
                                   PROP_K,
                                   g_param_spec_double ("k", 
                                                        NULL, 
                                                        "Mode k",
                                                        0.0, G_MAXDOUBLE, 1.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_TI,
                                   g_param_spec_double ("ti",
                                                        NULL,
                                                        "The initial time t_i",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_TF,
                                   g_param_spec_double ("tf",
                                                        NULL,
                                                        "The final time t_f",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 1.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_ADIAB_THRESHOLD,
                                   g_param_spec_double ("adiab-threshold",
                                                        NULL,
                                                        "The adiabatic threshold",
                                                        0.0, G_MAXDOUBLE, 1.0e-1,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_SAVE_EVOL,
                                   g_param_spec_boolean ("save-evol",
                                                         NULL,
                                                         "Save the system evolution",
                                                         TRUE,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  klass->eval_xi             = &_ncm_csq1d_eval_xi;
  klass->eval_nu             = &_ncm_csq1d_eval_nu;
  klass->eval_F1             = &_ncm_csq1d_eval_F1;
  klass->eval_F2             = &_ncm_csq1d_eval_F2;
  klass->eval_FN             = &_ncm_csq1d_eval_FN;
  klass->eval_powspec_factor = &_ncm_csq1d_eval_powspec_factor;
}

static gdouble 
_ncm_csq1d_eval_xi (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k)
{
  g_error ("_ncm_csq1d_eval_xi: not implemented."); 
  return 0.0;
}

static gdouble 
_ncm_csq1d_eval_nu (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k)
{
  g_error ("_ncm_csq1d_eval_nu: not implemented."); 
  return 0.0;
}

static gdouble 
_ncm_csq1d_eval_F1 (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k)
{
  g_error ("_ncm_csq1d_eval_F1: not implemented."); 
  return 0.0;
}

static gdouble 
_ncm_csq1d_eval_F2 (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k)
{
  g_error ("_ncm_csq1d_eval_F2: not implemented."); 
  return 0.0;
}

static gdouble 
_ncm_csq1d_eval_FN (NcmCSQ1D *csq1d, NcmModel *model, const gint n, const gdouble t, const gdouble k)
{
  g_error ("_ncm_csq1d_eval_FN: not implemented."); 
  return 0.0;
}

static gdouble 
_ncm_csq1d_eval_powspec_factor (NcmCSQ1D *csq1d, NcmModel *model, const gdouble k)
{
  g_error ("_ncm_csq1d_eval_powspec_factor: not implemented."); 
  return 0.0;
}

/**
 * ncm_csq1d_ref:
 * @csq1d: a #NcmCSQ1D
 *
 * Increases the reference count of @csq1d.
 *
 * Returns: (transfer full): @csq1d.
 */
NcmCSQ1D *
ncm_csq1d_ref (NcmCSQ1D *csq1d)
{
  return g_object_ref (csq1d);
}

/**
 * ncm_csq1d_free:
 * @csq1d: a #NcmCSQ1D
 *
 * Decreases the reference count of @csq1d.
 *
 */
void
ncm_csq1d_free (NcmCSQ1D *csq1d)
{
  g_object_unref (csq1d);
}

/**
 * ncm_csq1d_clear:
 * @csq1d: a #NcmCSQ1D
 *
 * Decreases the reference count of *@csq1d and sets the pointer *@csq1d to NULL.
 *
 */
void
ncm_csq1d_clear (NcmCSQ1D **csq1d)
{
  g_clear_object (csq1d);
}

/**
 * ncm_csq1d_set_reltol:
 * @csq1d: a #NcmCSQ1D
 * @reltol: relative tolerance
 *
 * Sets the relative tolerance to @reltol.
 *
 */
void 
ncm_csq1d_set_reltol (NcmCSQ1D *csq1d, const gdouble reltol)
{
  NcmCSQ1DPrivate * const self = csq1d->priv;
  self->reltol = reltol;
}

/**
 * ncm_csq1d_set_abstol:
 * @csq1d: a #NcmCSQ1D
 * @abstol: absolute tolerance
 *
 * Sets the absolute tolerance to @abstol.
 *
 */
void 
ncm_csq1d_set_abstol (NcmCSQ1D *csq1d, const gdouble abstol)
{
  NcmCSQ1DPrivate * const self = csq1d->priv;
  self->abstol = abstol;
}

/**
 * ncm_csq1d_set_k:
 * @csq1d: a #NcmCSQ1D
 * @k: mode $k$
 *
 * Sets the mode $k$ to @k.
 *
 */
void 
ncm_csq1d_set_k (NcmCSQ1D *csq1d, const gdouble k)
{
  NcmCSQ1DPrivate * const self = csq1d->priv;
  self->k = k;
}

/**
 * ncm_csq1d_set_ti:
 * @csq1d: a #NcmCSQ1D
 * @ti: mode $t_i$
 *
 * Sets the initial time $t_i$ to @ti.
 *
 */
void 
ncm_csq1d_set_ti (NcmCSQ1D *csq1d, const gdouble ti)
{
  NcmCSQ1DPrivate * const self = csq1d->priv;
  self->ti = ti;
}

/**
 * ncm_csq1d_set_tf:
 * @csq1d: a #NcmCSQ1D
 * @tf: mode $t_f$
 *
 * Sets the initial time $t_f$ to @tf.
 *
 */
void 
ncm_csq1d_set_tf (NcmCSQ1D *csq1d, const gdouble tf)
{
  NcmCSQ1DPrivate * const self = csq1d->priv;
  self->tf = tf;
}

/**
 * ncm_csq1d_set_adiab_threshold:
 * @csq1d: a #NcmCSQ1D
 * @adiab_threshold: mode $A_t$
 *
 * Sets the adiabatic threshold $A_t$.
 *
 */
void 
ncm_csq1d_set_adiab_threshold (NcmCSQ1D *csq1d, const gdouble adiab_threshold)
{
  NcmCSQ1DPrivate * const self = csq1d->priv;
  self->adiab_threshold = adiab_threshold;
}

/**
 * ncm_csq1d_set_save_evol:
 * @csq1d: a #NcmCSQ1D
 * @save: whether to save all evolution
 *
 * If true saves all evolution to be evaluted later through FIXME
 *
 */
void 
ncm_csq1d_set_save_evol (NcmCSQ1D *csq1d, gboolean save_evol)
{
  NcmCSQ1DPrivate * const self = csq1d->priv;
  if (self->save_evol != save_evol)
  {
    ncm_model_ctrl_force_update (self->ctrl);    
    self->save_evol = save_evol;
  }
}

/**
 * ncm_csq1d_get_reltol:
 * @csq1d: a #NcmCSQ1D
 *
 * Returns: the relative tolerance.
 */
gdouble 
ncm_csq1d_get_reltol (NcmCSQ1D *csq1d)
{
  NcmCSQ1DPrivate * const self = csq1d->priv;
  return self->reltol;
}

/**
 * ncm_csq1d_get_abstol:
 * @csq1d: a #NcmCSQ1D
 *
 * Returns: the absolute tolerance to @abstol.
 */
gdouble
ncm_csq1d_get_abstol (NcmCSQ1D *csq1d)
{
  NcmCSQ1DPrivate * const self = csq1d->priv;
  return self->abstol;
}

/**
 * ncm_csq1d_get_k:
 * @csq1d: a #NcmCSQ1D
 *
 * Returns: the mode $k$.
 */
gdouble
ncm_csq1d_get_k (NcmCSQ1D *csq1d)
{
  NcmCSQ1DPrivate * const self = csq1d->priv;
  return self->k;
}

/**
 * ncm_csq1d_get_ti:
 * @csq1d: a #NcmCSQ1D
 *
 * Returns: the initial time $t_i$.
 */
gdouble
ncm_csq1d_get_ti (NcmCSQ1D *csq1d)
{
  NcmCSQ1DPrivate * const self = csq1d->priv;
  return self->ti;
}

/**
 * ncm_csq1d_get_tf:
 * @csq1d: a #NcmCSQ1D
 *
 * Returns: the initial time $t_f$.
 */
gdouble
ncm_csq1d_get_tf (NcmCSQ1D *csq1d)
{
  NcmCSQ1DPrivate * const self = csq1d->priv;
  return self->tf;
}

/**
 * ncm_csq1d_get_adiab_threshold:
 * @csq1d: a #NcmCSQ1D
 *
 * Returns: the adiabatic threshold $A_t$.
 */
gdouble
ncm_csq1d_get_adiab_threshold (NcmCSQ1D *csq1d)
{
  NcmCSQ1DPrivate * const self = csq1d->priv;
  return self->adiab_threshold;
}

/**
 * ncm_csq1d_get_save_evol:
 * @csq1d: a #NcmCSQ1D
 *
 * Returns: whether the evolution will be saved.
 */
gboolean
ncm_csq1d_get_save_evol (NcmCSQ1D *csq1d)
{
  NcmCSQ1DPrivate * const self = csq1d->priv;
  return self->save_evol;
}

typedef struct _NcmCSQ1DWS
{
  NcmCSQ1D *csq1d;
  NcmModel *model;
  gdouble reltol;
} NcmCSQ1DWS;

static gint _ncm_csq1d_f (realtype t, N_Vector y, N_Vector ydot, gpointer f_data);
static gint _ncm_csq1d_J (realtype t, N_Vector y, N_Vector fy, SUNMatrix J, gpointer jac_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

static gint _ncm_csq1d_f_Up (realtype t, N_Vector y, N_Vector ydot, gpointer f_data);
static gint _ncm_csq1d_J_Up (realtype t, N_Vector y, N_Vector fy, SUNMatrix J, gpointer jac_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

static gint _ncm_csq1d_f_Um (realtype t, N_Vector y, N_Vector ydot, gpointer f_data);
static gint _ncm_csq1d_J_Um (realtype t, N_Vector y, N_Vector fy, SUNMatrix J, gpointer jac_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

static void
_ncm_csq1d_set_init_cond (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t0)
{
  NcmCSQ1DPrivate * const self = csq1d->priv;

  gdouble alpha, dgamma;

  ncm_csq1d_eval_adiab_at (csq1d, model, t0, &alpha, &dgamma, NULL, NULL);

  NV_Ith_S (self->y, 0) = sinh (alpha); /* chi == sinh (alpha) */
  NV_Ith_S (self->y, 1) = dgamma;

  self->t     = t0;
  self->state = NCM_CSQ1D_EVOL_STATE_ADIABATIC;
}

static void
_ncm_csq1d_prepare_integrator (NcmCSQ1D *csq1d, NcmCSQ1DWS *ws)
{
  NcmCSQ1DPrivate * const self = csq1d->priv;
  gint flag;

  if (!self->cvode_init)
  {
    self->cvode = CVodeCreate (CV_BDF);

    flag = CVodeInit (self->cvode, &_ncm_csq1d_f, self->t, self->y);
    NCM_CVODE_CHECK (&flag, "CVodeInit", 1, );

    self->cvode_init = TRUE;
  }
  else
  {
    flag = CVodeReInit (self->cvode, self->t, self->y);
    NCM_CVODE_CHECK (&flag, "CVodeInit", 1, );
  }

  flag = CVodeSStolerances (self->cvode, self->reltol, self->abstol);
  NCM_CVODE_CHECK (&flag, "CVodeSVtolerances", 1, );

  flag = CVodeSetMaxNumSteps (self->cvode, 0);
  NCM_CVODE_CHECK (&flag, "CVodeSetMaxNumSteps", 1, );

  flag = CVodeSetLinearSolver (self->cvode, self->LS, self->A);
  NCM_CVODE_CHECK (&flag, "CVodeSetLinearSolver", 1, );

  flag = CVodeSetJacFn (self->cvode, &_ncm_csq1d_J);
  NCM_CVODE_CHECK (&flag, "CVodeSetJacFn", 1, );

  flag = CVodeSetInitStep (self->cvode, fabs (self->t) * self->reltol);
  NCM_CVODE_CHECK (&flag, "CVodeSetInitStep", 1, );

  flag = CVodeSetUserData (self->cvode, ws);
  NCM_CVODE_CHECK (&flag, "CVodeSetUserData", 1, );

  flag = CVodeSetStopTime (self->cvode, self->tf);
  NCM_CVODE_CHECK (&flag, "CVodeSetStopTime", 1, );

  flag = CVodeSetMaxStep (self->cvode, (self->tf - self->t) / 100.0);
  NCM_CVODE_CHECK (&flag, "CVodeSetMaxStep", 1, );
}

static void
_ncm_csq1d_prepare_integrator_Up (NcmCSQ1D *csq1d, NcmCSQ1DWS *ws)
{
  NcmCSQ1DPrivate * const self = csq1d->priv;
  gint flag;

  if (!self->cvode_Up_init)
  {
    self->cvode_Up = CVodeCreate (CV_BDF);

    flag = CVodeInit (self->cvode_Up, &_ncm_csq1d_f_Up, self->t, self->y_Up);
    NCM_CVODE_CHECK (&flag, "CVodeInit", 1, );

    self->cvode_Up_init = TRUE;
  }
  else
  {
    flag = CVodeReInit (self->cvode_Up, self->t, self->y_Up);
    NCM_CVODE_CHECK (&flag, "CVodeInit", 1, );
  }

  flag = CVodeSStolerances (self->cvode_Up, self->reltol, self->abstol);
  NCM_CVODE_CHECK (&flag, "CVodeSVtolerances", 1, );

  flag = CVodeSetMaxNumSteps (self->cvode_Up, 100000);
  NCM_CVODE_CHECK (&flag, "CVodeSetMaxNumSteps", 1, );

  flag = CVodeSetLinearSolver (self->cvode_Up, self->LS_Up, self->A_Up);
  NCM_CVODE_CHECK (&flag, "CVodeSetLinearSolver", 1, );

  flag = CVodeSetJacFn (self->cvode_Up, &_ncm_csq1d_J_Up);
  NCM_CVODE_CHECK (&flag, "CVodeSetJacFn", 1, );

  flag = CVodeSetInitStep (self->cvode_Up, fabs (self->t) * self->reltol);
  NCM_CVODE_CHECK (&flag, "CVodeSetInitStep", 1, );  

  flag = CVodeSetUserData (self->cvode_Up, ws);
  NCM_CVODE_CHECK (&flag, "CVodeSetUserData", 1, );

  flag = CVodeSetStopTime (self->cvode_Up, self->tf);
  NCM_CVODE_CHECK (&flag, "CVodeSetStopTime", 1, );

  flag = CVodeSetMaxStep (self->cvode_Up, (self->tf - self->t) / 100.0);
  NCM_CVODE_CHECK (&flag, "CVodeSetMaxStep", 1, );
}

static void
_ncm_csq1d_prepare_integrator_Um (NcmCSQ1D *csq1d, NcmCSQ1DWS *ws)
{
  NcmCSQ1DPrivate * const self = csq1d->priv;
  gint flag;

  if (!self->cvode_Um_init)
  {
    self->cvode_Um = CVodeCreate (CV_BDF);

    flag = CVodeInit (self->cvode_Um, &_ncm_csq1d_f_Um, self->t, self->y_Um);
    NCM_CVODE_CHECK (&flag, "CVodeInit", 1, );

    self->cvode_Um_init = TRUE;
  }
  else
  {
    flag = CVodeReInit (self->cvode_Um, self->t, self->y_Um);
    NCM_CVODE_CHECK (&flag, "CVodeInit", 1, );
  }

  flag = CVodeSStolerances (self->cvode_Um, self->reltol, self->abstol);
  NCM_CVODE_CHECK (&flag, "CVodeSVtolerances", 1, );

  flag = CVodeSetMaxNumSteps (self->cvode_Um, 100000);
  NCM_CVODE_CHECK (&flag, "CVodeSetMaxNumSteps", 1, );

  flag = CVodeSetLinearSolver (self->cvode_Um, self->LS_Um, self->A_Um);
  NCM_CVODE_CHECK (&flag, "CVodeSetLinearSolver", 1, );

  flag = CVodeSetJacFn (self->cvode_Um, &_ncm_csq1d_J_Um);
  NCM_CVODE_CHECK (&flag, "CVodeSetJacFn", 1, );

  flag = CVodeSetInitStep (self->cvode_Um, fabs (self->t) * self->reltol);
  NCM_CVODE_CHECK (&flag, "CVodeSetInitStep", 1, );  

  flag = CVodeSetUserData (self->cvode_Um, ws);
  NCM_CVODE_CHECK (&flag, "CVodeSetUserData", 1, );

  flag = CVodeSetStopTime (self->cvode_Um, self->tf);
  NCM_CVODE_CHECK (&flag, "CVodeSetStopTime", 1, );

  flag = CVodeSetMaxStep (self->cvode_Um, (self->tf - self->t) / 100.0);
  NCM_CVODE_CHECK (&flag, "CVodeSetMaxStep", 1, );
}

static gint
_ncm_csq1d_f (realtype t, N_Vector y, N_Vector ydot, gpointer f_data)
{
  NcmCSQ1DWS *ws = (NcmCSQ1DWS *) f_data;
  NcmCSQ1DPrivate * const self = ws->csq1d->priv;

  const gdouble chi    = NV_Ith_S (y, 0);
  const gdouble dgamma = NV_Ith_S (y, 1);
  const gdouble s1     = sqrt (1.0 + chi * chi);

  const gdouble nu     = ncm_csq1d_eval_nu (ws->csq1d, ws->model, t, self->k);
  const gdouble F1     = ncm_csq1d_eval_F1 (ws->csq1d, ws->model, t, self->k);
  const gdouble twonu  = 2.0 * nu;
  
  NV_Ith_S (ydot, 0) = - twonu * sinh (dgamma) * s1;
  NV_Ith_S (ydot, 1) = + twonu * (-F1 + cosh (dgamma) * chi / s1);

  return 0;
}

static gint
_ncm_csq1d_f_Up (realtype t, N_Vector y, N_Vector ydot, gpointer f_data)
{
  NcmCSQ1DWS *ws = (NcmCSQ1DWS *) f_data;
  NcmCSQ1DPrivate * const self = ws->csq1d->priv;

  const gdouble chi   = NV_Ith_S (y, 0);
  const gdouble Up    = NV_Ith_S (y, 1);

  const gdouble xi      = ncm_csq1d_eval_xi (ws->csq1d, ws->model, t, self->k);
  const gdouble nu      = ncm_csq1d_eval_nu (ws->csq1d, ws->model, t, self->k);
  const gdouble twonu   = 2.0 * nu;
  const gdouble DUp     = Up - xi;
  const gdouble exp_DUp = exp (DUp);

  NV_Ith_S (ydot, 0) = + nu * ((1.0 + chi * chi) / exp_DUp - exp_DUp);
  NV_Ith_S (ydot, 1) = + twonu * chi / exp_DUp;

  return 0;
}

static gint
_ncm_csq1d_f_Um (realtype t, N_Vector y, N_Vector ydot, gpointer f_data)
{
  NcmCSQ1DWS *ws = (NcmCSQ1DWS *) f_data;
  NcmCSQ1DPrivate * const self = ws->csq1d->priv;

  const gdouble chi   = NV_Ith_S (y, 0);
  const gdouble Um    = NV_Ith_S (y, 1);

  const gdouble xi      = ncm_csq1d_eval_xi (ws->csq1d, ws->model, t, self->k);
  const gdouble nu      = ncm_csq1d_eval_nu (ws->csq1d, ws->model, t, self->k);
  const gdouble twonu   = 2.0 * nu;
  const gdouble DUm     = Um + xi;
  const gdouble exp_DUm = exp (DUm);

  NV_Ith_S (ydot, 0) = + nu * (exp_DUm - (1.0 + chi * chi) / exp_DUm);
  NV_Ith_S (ydot, 1) = - twonu * chi / exp_DUm;

  return 0;
}

static gint
_ncm_csq1d_J (realtype t, N_Vector y, N_Vector fy, SUNMatrix J, gpointer jac_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  NcmCSQ1DWS *ws = (NcmCSQ1DWS *) jac_data;
  NcmCSQ1DPrivate * const self = ws->csq1d->priv;

  const gdouble chi    = NV_Ith_S (y, 0);
  const gdouble dgamma = NV_Ith_S (y, 1);
  const gdouble s1     = sqrt (1.0 + chi * chi);
  
  const gdouble nu     = ncm_csq1d_eval_nu (ws->csq1d, ws->model, t, self->k);
  const gdouble twonu  = 2.0 * nu;

  /* - twonu * sinh (dgamma) * s1; */
  SM_ELEMENT_D (J, 0, 0) = - twonu * sinh (dgamma) * chi / s1;
  SM_ELEMENT_D (J, 0, 1) = - twonu * cosh (dgamma) * s1;

  /* + twonu * (-F1 + cosh (dgamma) * chi / s1); */
  SM_ELEMENT_D (J, 1, 0) = + twonu * cosh (dgamma) / gsl_pow_3 (s1);
  SM_ELEMENT_D (J, 1, 1) = + twonu * sinh (dgamma) * chi / s1;
  
  return 0;
}

static gint
_ncm_csq1d_J_Up (realtype t, N_Vector y, N_Vector fy, SUNMatrix J, gpointer jac_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  NcmCSQ1DWS *ws = (NcmCSQ1DWS *) jac_data;
  NcmCSQ1DPrivate * const self = ws->csq1d->priv;

  const gdouble chi     = NV_Ith_S (y, 0);
  const gdouble Up      = NV_Ith_S (y, 1);

  const gdouble xi      = ncm_csq1d_eval_xi (ws->csq1d, ws->model, t, self->k);
  const gdouble nu      = ncm_csq1d_eval_nu (ws->csq1d, ws->model, t, self->k);
  const gdouble twonu   = 2.0 * nu;
  const gdouble DUp     = Up - xi;
  const gdouble exp_DUp = exp (DUp);

  /* nu * ((1.0 + chi * chi) / exp_DUp - exp_DUp); */
  SM_ELEMENT_D (J, 0, 0) = + twonu * chi / exp_DUp;
  SM_ELEMENT_D (J, 0, 1) = - nu * ((1.0 + chi * chi) / exp_DUp + exp_DUp);

  /* + twonu * chi / exp_DUp; */
  SM_ELEMENT_D (J, 1, 0) = + twonu / exp_DUp;
  SM_ELEMENT_D (J, 1, 1) = - twonu * chi / exp_DUp;

  return 0;
}

static gint
_ncm_csq1d_J_Um (realtype t, N_Vector y, N_Vector fy, SUNMatrix J, gpointer jac_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  NcmCSQ1DWS *ws = (NcmCSQ1DWS *) jac_data;
  NcmCSQ1DPrivate * const self = ws->csq1d->priv;

  const gdouble chi     = NV_Ith_S (y, 0);
  const gdouble Um      = NV_Ith_S (y, 1);

  const gdouble xi      = ncm_csq1d_eval_xi (ws->csq1d, ws->model, t, self->k);
  const gdouble nu      = ncm_csq1d_eval_nu (ws->csq1d, ws->model, t, self->k);
  const gdouble twonu   = 2.0 * nu;
  const gdouble DUm     = Um + xi;
  const gdouble exp_DUm = exp (DUm);

  /* nu * (exp_DUm - (1.0 + chi * chi) / exp_DUm); */
  SM_ELEMENT_D (J, 0, 0) = - twonu * chi / exp_DUm;
  SM_ELEMENT_D (J, 0, 1) = + nu * (exp_DUm + (1.0 + chi * chi) / exp_DUm);

  /* - twonu * chi / exp_DUm; */
  SM_ELEMENT_D (J, 1, 0) = - twonu / exp_DUm;
  SM_ELEMENT_D (J, 1, 1) = + twonu * chi / exp_DUm;

  return 0;
}

/**
 * ncm_csq1d_eval_xi: (virtual eval_xi)
 * @csq1d: a #NcmCSQ1D
 * @model: (allow-none): a #NcmModel
 * @t: time $t$
 * @k: mode $k$
 *
 * Returns: $\xi$ 
 */
/**
 * ncm_csq1d_eval_nu: (virtual eval_nu)
 * @csq1d: a #NcmCSQ1D
 * @model: (allow-none): a #NcmModel
 * @t: time $t$
 * @k: mode $k$
 *
 * Returns: $\nu$ 
 */
/**
 * ncm_csq1d_eval_F1: (virtual eval_F1)
 * @csq1d: a #NcmCSQ1D
 * @model: (allow-none): a #NcmModel
 * @t: time $t$
 * @k: mode $k$
 *
 * Returns: $F_1$ 
 */
/**
 * ncm_csq1d_eval_F2: (virtual eval_F2)
 * @csq1d: a #NcmCSQ1D
 * @model: (allow-none): a #NcmModel
 * @t: time $t$
 * @k: mode $k$
 *
 * Returns: $F_2$ 
 */
/**
 * ncm_csq1d_eval_FN: (virtual eval_FN)
 * @csq1d: a #NcmCSQ1D
 * @model: (allow-none): a #NcmModel
 * @n: order $n$
 * @t: time $t$
 * @k: mode $k$
 *
 * Returns: $F_n$ 
 */
/**
 * ncm_csq1d_eval_powspec_factor: (virtual eval_powspec_factor)
 * @csq1d: a #NcmCSQ1D
 * @model: (allow-none): a #NcmModel
 * @k: mode $k$
 *
 * Returns: $F_n$ 
 */

static NcmCSQ1DEvolStop
_ncm_csq1d_evol_adiabatic (NcmCSQ1D *csq1d, NcmModel *model, GArray *asinh_t_a, GArray *alpha_a, GArray *dgamma_a)
{
  NcmCSQ1DPrivate * const self = csq1d->priv;
  NcmCSQ1DEvolStop reason = NCM_CSQ1D_EVOL_STOP_ERROR;
  gdouble last_asinh_t    = asinh (self->t);
  gint flag;

  g_assert (self->state == NCM_CSQ1D_EVOL_STATE_ADIABATIC);

  if (PRINT_EVOL)
    printf ("# ENTER EVOL ADIABATIC: TIME % 22.15g\n", self->t);

  while (TRUE)
  {
    gdouble asinh_t;
    gdouble chi, dgamma;
    gboolean is_finished = FALSE;

    flag = CVode (self->cvode, self->tf, self->y, &self->t, CV_ONE_STEP);
    NCM_CVODE_CHECK (&flag, "CVode[_ncm_csq1d_evol_adiabatic]", 1, NCM_CSQ1D_EVOL_STOP_ERROR);

    asinh_t = asinh (self->t);
    chi     = NV_Ith_S (self->y, 0);
    dgamma  = NV_Ith_S (self->y, 1);

    if (PRINT_EVOL)
    {
      const gdouble xi        = ncm_csq1d_eval_xi (csq1d, model, self->t, self->k);
      const gdouble gamma     = dgamma + xi;
      const gdouble log1pchi2 = log1p (chi * chi);
      const gdouble Up        = + gamma + 0.5 * log1pchi2;
      const gdouble Um        = - gamma + 0.5 * log1pchi2;
      printf ("# E[AD] % 22.15g % 22.15g % 22.15g | % 22.15g % 22.15g\n", self->t, chi, dgamma, Up, Um);
    }

    is_finished = (self->t == self->tf);
    
    if ((fabs ((asinh_t - last_asinh_t) / last_asinh_t) > 1.0e-5) || is_finished)
    {
      const gdouble alpha  = asinh (NV_Ith_S (self->y, 0));

      g_array_append_val (asinh_t_a, asinh_t);
      g_array_append_val (alpha_a,   alpha);
      g_array_append_val (dgamma_a,  dgamma);
      last_asinh_t = asinh_t;
    }
    
    if ((fabs (dgamma) > self->adiab_threshold) || (fabs (chi) > self->adiab_threshold))
    {
      if (dgamma > 0.0)
        reason = NCM_CSQ1D_EVOL_STOP_UP_START;
      else
        reason = NCM_CSQ1D_EVOL_STOP_UM_START;
      break;
    }

    if (is_finished)
    {
      reason = NCM_CSQ1D_EVOL_STOP_FINISHED;
      break;
    }
  }
  return reason;
}

static NcmCSQ1DEvolStop
_ncm_csq1d_evol_Up (NcmCSQ1D *csq1d, NcmModel *model, GArray *asinh_t_a, GArray *alpha_a, GArray *dgamma_a)
{
  NcmCSQ1DPrivate * const self = csq1d->priv;
  NcmCSQ1DEvolStop reason = NCM_CSQ1D_EVOL_STOP_ERROR;
  gdouble last_asinh_t    = asinh (self->t);
  gint flag;

  g_assert (self->state == NCM_CSQ1D_EVOL_STATE_UP);
  if (PRINT_EVOL)
    printf ("# ENTER EVOL UP: TIME % 22.15g\n", self->t);

  while (TRUE)
  {
    gdouble asinh_t;
    gdouble chi, dgamma, gamma, Up, xi;
    gboolean is_finished = FALSE;

    flag = CVode (self->cvode_Up, self->tf, self->y_Up, &self->t, CV_ONE_STEP);
    NCM_CVODE_CHECK (&flag, "CVode[_ncm_csq1d_evol_Up]", 1, NCM_CSQ1D_EVOL_STOP_ERROR);

    asinh_t = asinh (self->t);
    xi      = ncm_csq1d_eval_xi (csq1d, model, self->t, self->k);
    chi     = NV_Ith_S (self->y_Up, 0);
    Up      = NV_Ith_S (self->y_Up, 1);
    gamma   = Up - 0.5 * log1p (chi * chi);
    dgamma  = gamma - xi;

    if (PRINT_EVOL)
    {
      const gdouble log1pchi2 = log1p (chi * chi);
      const gdouble Um        = - Up + log1pchi2;
      printf ("# E[UP] % 22.15g % 22.15g % 22.15g | % 22.15g % 22.15g\n", self->t, chi, dgamma, Up, Um);
    }

    is_finished = (self->t == self->tf);
    
    if ((fabs ((asinh_t - last_asinh_t) / last_asinh_t) > 1.0e-5) || is_finished)
    {
      const gdouble alpha = asinh (chi);
      
      g_array_append_val (asinh_t_a, asinh_t);
      g_array_append_val (alpha_a,   alpha);
      g_array_append_val (dgamma_a,  dgamma);
      last_asinh_t = asinh_t;
    }
    
    if ((fabs (dgamma) < self->adiab_threshold) && (fabs (chi) < self->adiab_threshold))
    {
      reason = NCM_CSQ1D_EVOL_STOP_ADIABATIC_START;
      break;
    }

    if (dgamma < 0.0)
    {
      reason = NCM_CSQ1D_EVOL_STOP_UM_START;
      break;
    }

    if (is_finished)
    {
      reason = NCM_CSQ1D_EVOL_STOP_FINISHED;
      break;
    }
  }
  return reason;
}

static NcmCSQ1DEvolStop
_ncm_csq1d_evol_Um (NcmCSQ1D *csq1d, NcmModel *model, GArray *asinh_t_a, GArray *alpha_a, GArray *dgamma_a)
{
  NcmCSQ1DPrivate * const self = csq1d->priv;
  NcmCSQ1DEvolStop reason = NCM_CSQ1D_EVOL_STOP_ERROR;
  gdouble last_asinh_t    = asinh (self->t);
  gint flag;

  g_assert (self->state == NCM_CSQ1D_EVOL_STATE_UM);
  printf ("# ENTER EVOL UM: TIME % 22.15g\n", self->t);

  while (TRUE)
  {
    gdouble asinh_t;
    gdouble chi, dgamma, gamma, Um, xi;
    gboolean is_finished = FALSE;

    flag = CVode (self->cvode_Um, self->tf, self->y_Um, &self->t, CV_ONE_STEP);
    NCM_CVODE_CHECK (&flag, "CVode[_ncm_csq1d_evol_Um]", 1, NCM_CSQ1D_EVOL_STOP_ERROR);

    asinh_t = asinh (self->t);
    xi      = ncm_csq1d_eval_xi (csq1d, model, self->t, self->k);
    chi     = NV_Ith_S (self->y_Um, 0);
    Um      = NV_Ith_S (self->y_Um, 1);
    gamma   = -Um + 0.5 * log1p (chi * chi);
    dgamma  = gamma - xi;

    if (PRINT_EVOL)
    {
      const gdouble log1pchi2 = log1p (chi * chi);
      const gdouble Up        = - Um + log1pchi2;
      printf ("# E[UM] % 22.15g % 22.15g % 22.15g | % 22.15g % 22.15g\n", self->t, chi, dgamma, Up, Um);
    }

    is_finished = (self->t == self->tf);
    
    if ((fabs ((asinh_t - last_asinh_t) / last_asinh_t) > 1.0e-5) || is_finished)
    {
      const gdouble alpha = asinh (chi);
      
      g_array_append_val (asinh_t_a, asinh_t);
      g_array_append_val (alpha_a,   alpha);
      g_array_append_val (dgamma_a,  dgamma);
      last_asinh_t = asinh_t;
    }
    
    if ((fabs (dgamma) < self->adiab_threshold) && (fabs (chi) < self->adiab_threshold))
    {
      reason = NCM_CSQ1D_EVOL_STOP_ADIABATIC_START;
      break;
    }
    if (dgamma > 0.0)
    {
      reason = NCM_CSQ1D_EVOL_STOP_UP_START;
      break;
    }
    if (is_finished)
    {
      reason = NCM_CSQ1D_EVOL_STOP_FINISHED;
      break;
    }
  }
  return reason;
}


static void
_ncm_csq1d_evol_save (NcmCSQ1D *csq1d, NcmModel *model, NcmCSQ1DWS *ws, GArray *asinh_t_a, GArray *alpha_a, GArray *dgamma_a)
{
  NcmCSQ1DPrivate * const self = csq1d->priv;
  NcmCSQ1DEvolStop stop;

  if (PRINT_EVOL)
    printf ("# ENTERING EVOL SAVE: STATE %d, TIME % 22.15g\n", self->state, self->t);
  
  switch (self->state)
  {
    case NCM_CSQ1D_EVOL_STATE_ADIABATIC:
    {
      _ncm_csq1d_prepare_integrator (csq1d, ws);
      stop = _ncm_csq1d_evol_adiabatic (csq1d, model, asinh_t_a, alpha_a, dgamma_a);
      break;
    }
    case NCM_CSQ1D_EVOL_STATE_UP:
    {
      _ncm_csq1d_prepare_integrator_Up (csq1d, ws);
      stop = _ncm_csq1d_evol_Up (csq1d, model, asinh_t_a, alpha_a, dgamma_a);
      break;
    }
    case NCM_CSQ1D_EVOL_STATE_UM:
    {
      _ncm_csq1d_prepare_integrator_Um (csq1d, ws);
      stop = _ncm_csq1d_evol_Um (csq1d, model, asinh_t_a, alpha_a, dgamma_a);
      break;
    }
    default:
      g_assert_not_reached ();
      break;
  }

  switch (stop)
  {
    case NCM_CSQ1D_EVOL_STOP_UP_START:
    {
      switch (self->state)
      {
        case NCM_CSQ1D_EVOL_STATE_ADIABATIC:
        {
          const gdouble chi    = NV_Ith_S (self->y, 0);
          const gdouble dgamma = NV_Ith_S (self->y, 1);
          const gdouble xi     = ncm_csq1d_eval_xi (csq1d, model, self->t, self->k);
          const gdouble gamma  = xi + dgamma;
          const gdouble Up     = + gamma + 0.5 * log1p (chi * chi);

          NV_Ith_S (self->y_Up, 0) = chi;
          NV_Ith_S (self->y_Up, 1) = Up;
          break;
        }
        case NCM_CSQ1D_EVOL_STATE_UM:
        {
          const gdouble chi = NV_Ith_S (self->y_Um, 0);
          const gdouble Um  = NV_Ith_S (self->y_Um, 1);
          const gdouble Up  = -Um + log1p (chi * chi);

          NV_Ith_S (self->y_Up, 0) = chi;
          NV_Ith_S (self->y_Up, 1) = Up;
          break;
        }
        default:
          g_assert_not_reached ();
          break;
      }
      self->state = NCM_CSQ1D_EVOL_STATE_UP;
      _ncm_csq1d_evol_save (csq1d, model, ws, asinh_t_a, alpha_a, dgamma_a);
      break;
    }
    case NCM_CSQ1D_EVOL_STOP_UM_START:
    {
      switch (self->state)
      {
        case NCM_CSQ1D_EVOL_STATE_ADIABATIC:
        {
          const gdouble chi    = NV_Ith_S (self->y, 0);
          const gdouble dgamma = NV_Ith_S (self->y, 1);
          const gdouble xi     = ncm_csq1d_eval_xi (csq1d, model, self->t, self->k);
          const gdouble gamma  = xi + dgamma;
          const gdouble Um     = - gamma + 0.5 * log1p (chi * chi);

          NV_Ith_S (self->y_Um, 0) = chi;
          NV_Ith_S (self->y_Um, 1) = Um;
          break;
        }
        case NCM_CSQ1D_EVOL_STATE_UP:
        {
          const gdouble chi = NV_Ith_S (self->y_Up, 0);
          const gdouble Up  = NV_Ith_S (self->y_Up, 1);
          const gdouble Um  = -Up + log1p (chi * chi);

          NV_Ith_S (self->y_Um, 0) = chi;
          NV_Ith_S (self->y_Um, 1) = Um;
          break;
        }
        default:
          g_assert_not_reached ();
          break;
      }

      self->state = NCM_CSQ1D_EVOL_STATE_UM;
      _ncm_csq1d_evol_save (csq1d, model, ws, asinh_t_a, alpha_a, dgamma_a);
      break;
    }
    case NCM_CSQ1D_EVOL_STOP_ADIABATIC_START:
    {
      switch (self->state)
      {
        case NCM_CSQ1D_EVOL_STATE_UP:
        {
          const gdouble chi    = NV_Ith_S (self->y_Up, 0);
          const gdouble Up     = NV_Ith_S (self->y_Up, 1);
          const gdouble xi     = ncm_csq1d_eval_xi (csq1d, model, self->t, self->k);
          const gdouble gamma  = Up - 0.5 * log1p (chi * chi);
          const gdouble dgamma = gamma - xi;

          NV_Ith_S (self->y, 0) = chi;
          NV_Ith_S (self->y, 1) = dgamma;
          break;
        }
        case NCM_CSQ1D_EVOL_STATE_UM:
        {
          const gdouble chi    = NV_Ith_S (self->y_Um, 0);
          const gdouble Um     = NV_Ith_S (self->y_Um, 1);
          const gdouble xi     = ncm_csq1d_eval_xi (csq1d, model, self->t, self->k);
          const gdouble gamma  = - Um + 0.5 * log1p (chi * chi);
          const gdouble dgamma = gamma - xi;

          NV_Ith_S (self->y, 0) = chi;
          NV_Ith_S (self->y, 1) = dgamma;
        }
        default:
          g_assert_not_reached ();
          break;
      }

      self->state = NCM_CSQ1D_EVOL_STATE_ADIABATIC;      
      _ncm_csq1d_evol_save (csq1d, model, ws, asinh_t_a, alpha_a, dgamma_a);
      break;
    }
    case NCM_CSQ1D_EVOL_STOP_FINISHED:
      return;
      break;
    default:
      g_error ("_ncm_csq1d_evol_save: error evolving the ode system.");
      break;
  }  
}


/**
 * ncm_csq1d_prepare: (virtual prepare)
 * @csq1d: a #NcmCSQ1D
 * @model: (allow-none): a #NcmModel
 *
 * Prepares the object using @model.
 *
 */
void 
ncm_csq1d_prepare (NcmCSQ1D *csq1d, NcmModel *model)
{
  NcmCSQ1DPrivate * const self = csq1d->priv;
  NcmCSQ1DWS ws = {csq1d, model, 0.0};
  
  if (NCM_CSQ1D_GET_CLASS (csq1d)->prepare != NULL)
  {
    NCM_CSQ1D_GET_CLASS (csq1d)->prepare (csq1d, model);
  }

  g_assert_cmpfloat (self->tf, >, self->ti);

  _ncm_csq1d_set_init_cond (csq1d, model, self->ti);

  if (self->save_evol)
  {
    GArray *asinh_t_a    = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), 1000);
    GArray *alpha_a      = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), 1000);
    GArray *dgamma_a     = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), 1000);

    _ncm_csq1d_evol_save (csq1d, model, &ws, asinh_t_a, alpha_a, dgamma_a);

    ncm_spline_set_array (self->alpha_s,  asinh_t_a, alpha_a,  TRUE);
    ncm_spline_set_array (self->dgamma_s, asinh_t_a, dgamma_a, TRUE);

    g_array_unref (asinh_t_a);
    g_array_unref (alpha_a);
    g_array_unref (dgamma_a);
  }
  else
  {
    g_assert_not_reached ();
  }
}

/**
 * ncm_csq1d_get_time_array:
 * @csq1d: a #NcmCSQ1D
 * @smallest_t: (out) (allow-none): the smallest absolute value of $t$ in the array
 *
 * Returns: (transfer full) (element-type gdouble): the time array of the computed steps.
 */
GArray *
ncm_csq1d_get_time_array (NcmCSQ1D *csq1d, gdouble *smallest_t)
{
  NcmCSQ1DPrivate * const self = csq1d->priv;
  NcmVector *asinh_t_v = ncm_spline_get_xv (self->alpha_s);
  const guint len      = ncm_vector_len (asinh_t_v);
  GArray *t_a          = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), len);
  gdouble s_t          = 1.0e300;
  gint i;
  
  for (i = 0; i < len; i++)
  {
    const gdouble t_i = sinh (ncm_vector_fast_get (asinh_t_v, i));
    const gdouble a_t = fabs (t_i);
    
    g_array_append_val (t_a, t_i);
    s_t = MIN (a_t, s_t);
  }

  if (smallest_t != NULL)
    smallest_t[0] = s_t;

  return t_a;
}

gdouble
_ncm_csq1d_find_adiab_time_limit_f (gdouble t, gpointer params)
{
  NcmCSQ1DWS *ws = (NcmCSQ1DWS *) params;
  /*NcmCSQ1DPrivate * const self = ws->csq1d->priv;*/

  gdouble alpha, dgamma, cmp, alpha_reltol, dgamma_reltol;

  ncm_csq1d_eval_adiab_at (ws->csq1d, ws->model, t, &alpha, &dgamma, &alpha_reltol, &dgamma_reltol);

  cmp = MAX (alpha_reltol, dgamma_reltol);

  /*printf ("NHAC % 22.15g % 22.15g % 22.15g % 22.15g\n", t, fabs (alpha), fabs (dgamma), ncm_csq1d_eval_xi (ws->csq1d, ws->model, t, self->k));*/
  
  return log (cmp / ws->reltol);
}

/**
 * ncm_csq1d_find_adiab_time_limit:
 * @csq1d: a #NcmCSQ1D
 * @model: (allow-none): a #NcmModel
 * @t0: time lower bound $t_0$
 * @t1: time upper bound $t_1$
 * @reltol: relative tolerance
 * @ti: (out): adiabatic time limit $t_i$
 *
 * Computes the time upper limit $t_i \in [t_0, t_1]$ where the adiabatic
 * approximation is satisfied up to @reltol.
 * 
 * Returns: whether the time limit was found.
 */
gboolean 
ncm_csq1d_find_adiab_time_limit (NcmCSQ1D *csq1d, NcmModel *model, gdouble t0, gdouble t1, const gdouble reltol, gdouble *ti)
{
  /*NcmCSQ1DPrivate * const self = csq1d->priv;*/
  gdouble alpha0, dgamma0, alpha_reltol0, dgamma_reltol0;
  gdouble alpha1, dgamma1, alpha_reltol1, dgamma_reltol1;
  
  gboolean adiab0, adiab1;

  g_assert_cmpfloat (reltol, >, 0.0);

  ncm_csq1d_eval_adiab_at (csq1d, model, t0, &alpha0, &dgamma0, &alpha_reltol0, &dgamma_reltol0);
  ncm_csq1d_eval_adiab_at (csq1d, model, t1, &alpha1, &dgamma1, &alpha_reltol1, &dgamma_reltol1);

  adiab0 = ((fabs (alpha_reltol0) < reltol) && (fabs (dgamma_reltol0) < reltol));
  adiab1 = ((fabs (alpha_reltol1) < reltol) && (fabs (dgamma_reltol1) < reltol));
  
  if ((adiab0 && adiab1) || (!adiab0 && !adiab1))
  {
    /*printf ("t0 % 22.15g % 22.15g % 22.15g t1 % 22.15g % 22.15g % 22.15g\n", t0, alpha0, dgamma0, t1, alpha1, dgamma1);*/
    return FALSE;
  }
  else
  {
    NcmCSQ1DWS ws = {csq1d, model, reltol};
    gint iter = 0, max_iter = 1000;
    const gdouble root_reltol = 1.0e-2;
    const gsl_root_fsolver_type *T;
    gsl_root_fsolver *s;
    gsl_function F;
    gint status;

    ti[0] = 0.5 * (t0 + t1);

    F.function = &_ncm_csq1d_find_adiab_time_limit_f;
    F.params   = &ws;

    T = gsl_root_fsolver_brent;
    s = gsl_root_fsolver_alloc (T);
    gsl_root_fsolver_set (s, &F, t0, t1);

    do
    {
      iter++;
      status = gsl_root_fsolver_iterate (s);
      ti[0]  = gsl_root_fsolver_root (s);
      t0     = gsl_root_fsolver_x_lower (s);
      t1     = gsl_root_fsolver_x_upper (s);
      status = gsl_root_test_interval (t0, t1, 0.0, root_reltol);

      /*printf ("%5d [%.7f, %.7f] %.7f %+.7f\n", iter, t0, t1, ti[0], t1 - t0);*/
    }
    while (status == GSL_CONTINUE && iter < max_iter);

    gsl_root_fsolver_free (s);    
  }

  return TRUE;
}

static gdouble 
_ncm_csq1d_F2_func (const gdouble t, gpointer user_data)
{
  NcmCSQ1DWS *ws = (NcmCSQ1DWS *) user_data;
  NcmCSQ1DPrivate * const self = ws->csq1d->priv;
  const gdouble F2 = ncm_csq1d_eval_F2 (ws->csq1d, ws->model, t, self->k);

  return F2;
}

static gdouble 
_ncm_csq1d_lnnu_func (const gdouble t, gpointer user_data)
{
  NcmCSQ1DWS *ws = (NcmCSQ1DWS *) user_data;
  NcmCSQ1DPrivate * const self = ws->csq1d->priv;
  const gdouble nu = ncm_csq1d_eval_nu (ws->csq1d, ws->model, t, self->k);

  return log (nu);
}


/**
 * ncm_csq1d_eval_adiab_at:
 * @csq1d: a #NcmCSQ1D
 * @model: (allow-none): a #NcmModel
 * @t: time $t$
 * @alpha: (out): value of $\alpha(t)$
 * @dgamma: (out): value of $\Delta\gamma(t)$
 * @alpha_reltol: (out) (allow-none): estimated error on $\alpha(t)$
 * @dgamma_reltol: (out) (allow-none): estimated error on $\Delta\gamma(t)$
 *
 * Computes the value of the adiabatic approximation of the variables $\alpha$ and $\Delta\gamma$ at $t$.
 *
 */
void 
ncm_csq1d_eval_adiab_at (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, gdouble *alpha, gdouble *dgamma, gdouble *alpha_reltol, gdouble *dgamma_reltol)
{
  NcmCSQ1DPrivate * const self = csq1d->priv;
  NcmCSQ1DWS ws       = {csq1d, model, 0.0};
  const gdouble F1    = ncm_csq1d_eval_F1 (csq1d, model, t, self->k);
  const gdouble F2    = ncm_csq1d_eval_F2 (csq1d, model, t, self->k);
  const gdouble F1_2  = F1 * F1;
  const gdouble F1_3  = F1_2 * F1;
  const gdouble nu    = ncm_csq1d_eval_nu (csq1d, model, t, self->k);
  const gdouble twonu = 2.0 * nu;
  gdouble err, F3, d2F2, dlnnu, F4, alpha_reltol0, dgamma_reltol0;

  F3             = ncm_diff_rc_d1_1_to_1 (self->diff, t, &_ncm_csq1d_F2_func, &ws, &err) / twonu;
  d2F2           = ncm_diff_rc_d2_1_to_1 (self->diff, t, &_ncm_csq1d_F2_func, &ws, &err);
  dlnnu          = ncm_diff_rc_d1_1_to_1 (self->diff, t, &_ncm_csq1d_lnnu_func, &ws, &err);
  F4             = d2F2 / gsl_pow_2 (twonu) - dlnnu * F3 / twonu;
  alpha_reltol0  = gsl_pow_2 ((F1_3 / 3.0 - F3) / F1);
  dgamma_reltol0 = gsl_pow_2 ((F4 - F1_2 * F2) / F2);
    
  /*printf ("F3(% 22.15g) = % 22.15g +/- % 22.15g\n", t, F3, err / (2.0 * nu));*/
  /*printf ("F1, F2, F3, F4 = % 22.15g % 22.15e % 22.15e % 22.15e % 22.15e % 22.15e % 22.15e\n", t, F1, F2, F3, F4,  alpha_reltol0, dgamma_reltol0);*/
  
  alpha[0]  = + F1 + F1_3 / 3.0 - F3;
  dgamma[0] = - (1.0 + F1_2) * F2 + F4;

  if (alpha_reltol != NULL)
    alpha_reltol[0] = alpha_reltol0;
  if (dgamma_reltol != NULL)
    dgamma_reltol[0] = dgamma_reltol0;
}

/**
 * ncm_csq1d_eval_at:
 * @csq1d: a #NcmCSQ1D
 * @t: time $t$
 * @alpha: (out): value of $\alpha(t)$
 * @dgamma: (out): value of $\Delta\gamma(t)$
 *
 * Computes the value of the variables $\alpha$ and $\Delta\gamma$ at $t$.
 *
 */
void 
ncm_csq1d_eval_at (NcmCSQ1D *csq1d, const gdouble t, gdouble *alpha, gdouble *dgamma)
{
  NcmCSQ1DPrivate * const self = csq1d->priv;
  const gdouble a_t = asinh (t);
  
  alpha[0]  = ncm_spline_eval (self->alpha_s, a_t);
  dgamma[0] = ncm_spline_eval (self->dgamma_s, a_t);
}

/**
 * ncm_csq1d_alpha_dgamma_to_phi_Pphi:
 * @csq1d: a #NcmCSQ1D
 * @model: (allow-none): a #NcmModel
 * @t: time $t$
 * @alpha: value of $\alpha(t)$
 * @dgamma: value of $\Delta\gamma(t)$
 * @phi: (out caller-allocates) (array fixed-size=2): real and imaginary parts of $\phi$, i.e., $[\mathrm{Re}(\phi), \mathrm{Im}(\phi)]$
 * @Pphi: (out caller-allocates) (array fixed-size=2): real and imaginary parts of $\Pi_\phi$, i.e., $[\mathrm{Re}(\Pi_\phi), \mathrm{Im}(\Pi_\phi)]$
 *
 * Computes the value of the variables $\alpha$ and $\Delta\gamma$ at $t$.
 *
 */
void 
ncm_csq1d_alpha_dgamma_to_phi_Pphi (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble alpha, const gdouble dgamma, gdouble *phi, gdouble *Pphi)
{
  NcmCSQ1DPrivate * const self = csq1d->priv;
  const gdouble gamma               = ncm_csq1d_eval_xi (csq1d, model, t, self->k) + dgamma;
  const gdouble exp_gamma_p_alpha_2 = exp (0.5 * (gamma + alpha));
  const gdouble exp_gamma_m_alpha_2 = exp (0.5 * (gamma - alpha));

  /*printf ("=> % 22.15g % 22.15g % 22.15g % 22.15g\n", gamma, alpha, ncm_csq1d_eval_xi (csq1d, model, t, self->k), dgamma);*/
  
  phi[0]  = +0.5 / exp_gamma_m_alpha_2;
  phi[1]  = -0.5 / exp_gamma_p_alpha_2;

  Pphi[0] = -0.5 * exp_gamma_p_alpha_2;
  Pphi[1] = -0.5 * exp_gamma_m_alpha_2;
}

/**
 * ncm_csq1d_get_J_at:
 * @csq1d: a #NcmCSQ1D
 * @model: (allow-none): a #NcmModel
 * @t: time $t$
 * @J11: (out): $J_{11}$
 * @J12: (out): $J_{12}$
 * @J22: (out): $J_{22}$
 * 
 * Computes the complex structure matrix.
 * 
 */
void 
ncm_csq1d_get_J_at (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, gdouble *J11, gdouble *J12, gdouble *J22)
{
  NcmCSQ1DPrivate * const self = csq1d->priv;
  const gdouble a_t = asinh (t);
  
  const gdouble alpha  = ncm_spline_eval (self->alpha_s, a_t);
  const gdouble dgamma = ncm_spline_eval (self->dgamma_s, a_t);
  const gdouble gamma  = ncm_csq1d_eval_xi (csq1d, model, t, self->k) + dgamma;

  J11[0] = cosh (alpha) * exp (-gamma);
  J22[0] = cosh (alpha) * exp (+gamma);
  J12[0] = -sinh (alpha);
}