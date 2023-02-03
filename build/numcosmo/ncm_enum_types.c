
/* This file is generated by glib-mkenums, do not modify it. This code is licensed under the same license as the containing project. Note that it links to GLib, so must comply with the LGPL linking clauses. */


#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

/**
 * SECTION:ncm_enum_types     
 * @title: NcmEnumTypes                           
 * @short_description: Automaticaly generated enum types from NumCosmoMath library.
 *   
 */
/* enumerations from "math/ncm_csq1d.h" */
GType
ncm_cs_q1_devol_state_get_type (void)
{
  static gsize static_g_enum_type_id = 0;

  if (g_once_init_enter (&static_g_enum_type_id))
  {
    static const GEnumValue values[] = {
      { NCM_CSQ1D_EVOL_STATE_INVALID, "NCM_CSQ1D_EVOL_STATE_INVALID", "invalid" },
      { NCM_CSQ1D_EVOL_STATE_ADIABATIC, "NCM_CSQ1D_EVOL_STATE_ADIABATIC", "adiabatic" },
      { NCM_CSQ1D_EVOL_STATE_UP, "NCM_CSQ1D_EVOL_STATE_UP", "up" },
      { NCM_CSQ1D_EVOL_STATE_UM, "NCM_CSQ1D_EVOL_STATE_UM", "um" },
      { 0, NULL, NULL }
    };
    GType g_enum_type_id =
      g_enum_register_static (g_intern_static_string ("NcmCSQ1DEvolState"), values);
    g_once_init_leave (&static_g_enum_type_id, g_enum_type_id);
  }

  return static_g_enum_type_id;
}
/* enumerations from "math/ncm_data_poisson.h" */
GType
ncm_data_poisson_type_get_type (void)
{
  static gsize static_g_enum_type_id = 0;

  if (g_once_init_enter (&static_g_enum_type_id))
  {
    static const GEnumValue values[] = {
      { NCM_DATA_POISSON_INT, "NCM_DATA_POISSON_INT", "int" },
      { 0, NULL, NULL }
    };
    GType g_enum_type_id =
      g_enum_register_static (g_intern_static_string ("NcmDataPoissonType"), values);
    g_once_init_leave (&static_g_enum_type_id, g_enum_type_id);
  }

  return static_g_enum_type_id;
}
/* enumerations from "math/ncm_dataset.h" */
GType
ncm_dataset_bstrap_type_get_type (void)
{
  static gsize static_g_enum_type_id = 0;

  if (g_once_init_enter (&static_g_enum_type_id))
  {
    static const GEnumValue values[] = {
      { NCM_DATASET_BSTRAP_DISABLE, "NCM_DATASET_BSTRAP_DISABLE", "disable" },
      { NCM_DATASET_BSTRAP_PARTIAL, "NCM_DATASET_BSTRAP_PARTIAL", "partial" },
      { NCM_DATASET_BSTRAP_TOTAL, "NCM_DATASET_BSTRAP_TOTAL", "total" },
      { 0, NULL, NULL }
    };
    GType g_enum_type_id =
      g_enum_register_static (g_intern_static_string ("NcmDatasetBStrapType"), values);
    g_once_init_leave (&static_g_enum_type_id, g_enum_type_id);
  }

  return static_g_enum_type_id;
}
/* enumerations from "math/ncm_fit.h" */
GType
ncm_fit_type_get_type (void)
{
  static gsize static_g_enum_type_id = 0;

  if (g_once_init_enter (&static_g_enum_type_id))
  {
    static const GEnumValue values[] = {
      { NCM_FIT_TYPE_GSL_LS, "NCM_FIT_TYPE_GSL_LS", "gsl-ls" },
      { NCM_FIT_TYPE_GSL_MM, "NCM_FIT_TYPE_GSL_MM", "gsl-mm" },
      { NCM_FIT_TYPE_GSL_MMS, "NCM_FIT_TYPE_GSL_MMS", "gsl-mms" },
      { NCM_FIT_TYPE_LEVMAR, "NCM_FIT_TYPE_LEVMAR", "levmar" },
      { NCM_FIT_TYPE_NLOPT, "NCM_FIT_TYPE_NLOPT", "nlopt" },
      { 0, NULL, NULL }
    };
    GType g_enum_type_id =
      g_enum_register_static (g_intern_static_string ("NcmFitType"), values);
    g_once_init_leave (&static_g_enum_type_id, g_enum_type_id);
  }

  return static_g_enum_type_id;
}
GType
ncm_fit_grad_type_get_type (void)
{
  static gsize static_g_enum_type_id = 0;

  if (g_once_init_enter (&static_g_enum_type_id))
  {
    static const GEnumValue values[] = {
      { NCM_FIT_GRAD_ANALYTICAL, "NCM_FIT_GRAD_ANALYTICAL", "analytical" },
      { NCM_FIT_GRAD_NUMDIFF_FORWARD, "NCM_FIT_GRAD_NUMDIFF_FORWARD", "numdiff-forward" },
      { NCM_FIT_GRAD_NUMDIFF_CENTRAL, "NCM_FIT_GRAD_NUMDIFF_CENTRAL", "numdiff-central" },
      { NCM_FIT_GRAD_NUMDIFF_ACCURATE, "NCM_FIT_GRAD_NUMDIFF_ACCURATE", "numdiff-accurate" },
      { 0, NULL, NULL }
    };
    GType g_enum_type_id =
      g_enum_register_static (g_intern_static_string ("NcmFitGradType"), values);
    g_once_init_leave (&static_g_enum_type_id, g_enum_type_id);
  }

  return static_g_enum_type_id;
}
GType
ncm_fit_run_msgs_get_type (void)
{
  static gsize static_g_enum_type_id = 0;

  if (g_once_init_enter (&static_g_enum_type_id))
  {
    static const GEnumValue values[] = {
      { NCM_FIT_RUN_MSGS_NONE, "NCM_FIT_RUN_MSGS_NONE", "none" },
      { NCM_FIT_RUN_MSGS_SIMPLE, "NCM_FIT_RUN_MSGS_SIMPLE", "simple" },
      { NCM_FIT_RUN_MSGS_FULL, "NCM_FIT_RUN_MSGS_FULL", "full" },
      { 0, NULL, NULL }
    };
    GType g_enum_type_id =
      g_enum_register_static (g_intern_static_string ("NcmFitRunMsgs"), values);
    g_once_init_leave (&static_g_enum_type_id, g_enum_type_id);
  }

  return static_g_enum_type_id;
}
/* enumerations from "math/ncm_fit_esmcmc_walker_apes.h" */
GType
ncm_fit_esmcmc_walker_apes_method_get_type (void)
{
  static gsize static_g_enum_type_id = 0;

  if (g_once_init_enter (&static_g_enum_type_id))
  {
    static const GEnumValue values[] = {
      { NCM_FIT_ESMCMC_WALKER_APES_METHOD_KDE, "NCM_FIT_ESMCMC_WALKER_APES_METHOD_KDE", "kde" },
      { NCM_FIT_ESMCMC_WALKER_APES_METHOD_VKDE, "NCM_FIT_ESMCMC_WALKER_APES_METHOD_VKDE", "vkde" },
      { 0, NULL, NULL }
    };
    GType g_enum_type_id =
      g_enum_register_static (g_intern_static_string ("NcmFitESMCMCWalkerAPESMethod"), values);
    g_once_init_leave (&static_g_enum_type_id, g_enum_type_id);
  }

  return static_g_enum_type_id;
}
GType
ncm_fit_esmcmc_walker_apes_ktype_get_type (void)
{
  static gsize static_g_enum_type_id = 0;

  if (g_once_init_enter (&static_g_enum_type_id))
  {
    static const GEnumValue values[] = {
      { NCM_FIT_ESMCMC_WALKER_APES_KTYPE_CAUCHY, "NCM_FIT_ESMCMC_WALKER_APES_KTYPE_CAUCHY", "cauchy" },
      { NCM_FIT_ESMCMC_WALKER_APES_KTYPE_ST3, "NCM_FIT_ESMCMC_WALKER_APES_KTYPE_ST3", "st3" },
      { NCM_FIT_ESMCMC_WALKER_APES_KTYPE_GAUSS, "NCM_FIT_ESMCMC_WALKER_APES_KTYPE_GAUSS", "gauss" },
      { 0, NULL, NULL }
    };
    GType g_enum_type_id =
      g_enum_register_static (g_intern_static_string ("NcmFitESMCMCWalkerAPESKType"), values);
    g_once_init_leave (&static_g_enum_type_id, g_enum_type_id);
  }

  return static_g_enum_type_id;
}
/* enumerations from "math/ncm_fit_gsl_mm.h" */
GType
ncm_fit_gslmm_algos_get_type (void)
{
  static gsize static_g_enum_type_id = 0;

  if (g_once_init_enter (&static_g_enum_type_id))
  {
    static const GEnumValue values[] = {
      { NCM_FIT_GSL_MM_CONJUGATE_FR, "NCM_FIT_GSL_MM_CONJUGATE_FR", "conjugate-fr" },
      { NCM_FIT_GSL_MM_CONJUGATE_PR, "NCM_FIT_GSL_MM_CONJUGATE_PR", "conjugate-pr" },
      { NCM_FIT_GSL_MM_VECTOR_BFGS, "NCM_FIT_GSL_MM_VECTOR_BFGS", "vector-bfgs" },
      { NCM_FIT_GSL_MM_VECTOR_BFGS2, "NCM_FIT_GSL_MM_VECTOR_BFGS2", "vector-bfgs2" },
      { NCM_FIT_GSL_MM_STEEPEST_DESCENT, "NCM_FIT_GSL_MM_STEEPEST_DESCENT", "steepest-descent" },
      { 0, NULL, NULL }
    };
    GType g_enum_type_id =
      g_enum_register_static (g_intern_static_string ("NcmFitGSLMMAlgos"), values);
    g_once_init_leave (&static_g_enum_type_id, g_enum_type_id);
  }

  return static_g_enum_type_id;
}
/* enumerations from "math/ncm_fit_gsl_mms.h" */
GType
ncm_fit_gslmms_algos_get_type (void)
{
  static gsize static_g_enum_type_id = 0;

  if (g_once_init_enter (&static_g_enum_type_id))
  {
    static const GEnumValue values[] = {
      { NCM_FIT_GSL_MMS_NMSIMPLEX2, "NCM_FIT_GSL_MMS_NMSIMPLEX2", "nmsimplex2" },
      { NCM_FIT_GSL_MMS_NMSIMPLEX, "NCM_FIT_GSL_MMS_NMSIMPLEX", "nmsimplex" },
      { NCM_FIT_GSL_MMS_NMSIMPLEX2RAND, "NCM_FIT_GSL_MMS_NMSIMPLEX2RAND", "nmsimplex2rand" },
      { 0, NULL, NULL }
    };
    GType g_enum_type_id =
      g_enum_register_static (g_intern_static_string ("NcmFitGSLMMSAlgos"), values);
    g_once_init_leave (&static_g_enum_type_id, g_enum_type_id);
  }

  return static_g_enum_type_id;
}
/* enumerations from "math/ncm_fit_levmar.h" */
GType
ncm_fit_levmar_algos_get_type (void)
{
  static gsize static_g_enum_type_id = 0;

  if (g_once_init_enter (&static_g_enum_type_id))
  {
    static const GEnumValue values[] = {
      { NCM_FIT_LEVMAR_DER, "NCM_FIT_LEVMAR_DER", "der" },
      { NCM_FIT_LEVMAR_DIF, "NCM_FIT_LEVMAR_DIF", "dif" },
      { NCM_FIT_LEVMAR_BC_DER, "NCM_FIT_LEVMAR_BC_DER", "bc-der" },
      { NCM_FIT_LEVMAR_BC_DIF, "NCM_FIT_LEVMAR_BC_DIF", "bc-dif" },
      { 0, NULL, NULL }
    };
    GType g_enum_type_id =
      g_enum_register_static (g_intern_static_string ("NcmFitLevmarAlgos"), values);
    g_once_init_leave (&static_g_enum_type_id, g_enum_type_id);
  }

  return static_g_enum_type_id;
}
/* enumerations from "math/ncm_fit_mc.h" */
GType
ncm_fit_mc_resample_type_get_type (void)
{
  static gsize static_g_enum_type_id = 0;

  if (g_once_init_enter (&static_g_enum_type_id))
  {
    static const GEnumValue values[] = {
      { NCM_FIT_MC_RESAMPLE_FROM_MODEL, "NCM_FIT_MC_RESAMPLE_FROM_MODEL", "from-model" },
      { NCM_FIT_MC_RESAMPLE_BOOTSTRAP_NOMIX, "NCM_FIT_MC_RESAMPLE_BOOTSTRAP_NOMIX", "bootstrap-nomix" },
      { NCM_FIT_MC_RESAMPLE_BOOTSTRAP_MIX, "NCM_FIT_MC_RESAMPLE_BOOTSTRAP_MIX", "bootstrap-mix" },
      { 0, NULL, NULL }
    };
    GType g_enum_type_id =
      g_enum_register_static (g_intern_static_string ("NcmFitMCResampleType"), values);
    g_once_init_leave (&static_g_enum_type_id, g_enum_type_id);
  }

  return static_g_enum_type_id;
}
/* enumerations from "math/ncm_function_cache.h" */
GType
ncm_function_cache_search_type_get_type (void)
{
  static gsize static_g_enum_type_id = 0;

  if (g_once_init_enter (&static_g_enum_type_id))
  {
    static const GEnumValue values[] = {
      { NCM_FUNCTION_CACHE_SEARCH_BOTH, "NCM_FUNCTION_CACHE_SEARCH_BOTH", "both" },
      { NCM_FUNCTION_CACHE_SEARCH_GT, "NCM_FUNCTION_CACHE_SEARCH_GT", "gt" },
      { NCM_FUNCTION_CACHE_SEARCH_LT, "NCM_FUNCTION_CACHE_SEARCH_LT", "lt" },
      { 0, NULL, NULL }
    };
    GType g_enum_type_id =
      g_enum_register_static (g_intern_static_string ("NcmFunctionCacheSearchType"), values);
    g_once_init_leave (&static_g_enum_type_id, g_enum_type_id);
  }

  return static_g_enum_type_id;
}
/* enumerations from "math/ncm_hoaa.h" */
GType
ncm_hoaa_opt_get_type (void)
{
  static gsize static_g_enum_type_id = 0;

  if (g_once_init_enter (&static_g_enum_type_id))
  {
    static const GEnumValue values[] = {
      { NCM_HOAA_OPT_FULL, "NCM_HOAA_OPT_FULL", "full" },
      { NCM_HOAA_OPT_V_ONLY, "NCM_HOAA_OPT_V_ONLY", "v-only" },
      { NCM_HOAA_OPT_DLNMNU_ONLY, "NCM_HOAA_OPT_DLNMNU_ONLY", "dlnmnu-only" },
      { NCM_HOAA_OPT_INVALID, "NCM_HOAA_OPT_INVALID", "invalid" },
      { 0, NULL, NULL }
    };
    GType g_enum_type_id =
      g_enum_register_static (g_intern_static_string ("NcmHOAAOpt"), values);
    g_once_init_leave (&static_g_enum_type_id, g_enum_type_id);
  }

  return static_g_enum_type_id;
}
GType
ncm_hoaa_sing_type_get_type (void)
{
  static gsize static_g_enum_type_id = 0;

  if (g_once_init_enter (&static_g_enum_type_id))
  {
    static const GEnumValue values[] = {
      { NCM_HOAA_SING_TYPE_ZERO, "NCM_HOAA_SING_TYPE_ZERO", "zero" },
      { NCM_HOAA_SING_TYPE_INF, "NCM_HOAA_SING_TYPE_INF", "inf" },
      { NCM_HOAA_SING_TYPE_INVALID, "NCM_HOAA_SING_TYPE_INVALID", "invalid" },
      { 0, NULL, NULL }
    };
    GType g_enum_type_id =
      g_enum_register_static (g_intern_static_string ("NcmHOAASingType"), values);
    g_once_init_leave (&static_g_enum_type_id, g_enum_type_id);
  }

  return static_g_enum_type_id;
}
GType
ncm_hoaa_var_get_type (void)
{
  static gsize static_g_enum_type_id = 0;

  if (g_once_init_enter (&static_g_enum_type_id))
  {
    static const GEnumValue values[] = {
      { NCM_HOAA_VAR_QBAR, "NCM_HOAA_VAR_QBAR", "qbar" },
      { NCM_HOAA_VAR_PBAR, "NCM_HOAA_VAR_PBAR", "pbar" },
      { NCM_HOAA_VAR_UPSILON, "NCM_HOAA_VAR_UPSILON", "upsilon" },
      { NCM_HOAA_VAR_GAMMA, "NCM_HOAA_VAR_GAMMA", "gamma" },
      { NCM_HOAA_VAR_SYS_SIZE, "NCM_HOAA_VAR_SYS_SIZE", "sys-size" },
      { 0, NULL, NULL }
    };
    GType g_enum_type_id =
      g_enum_register_static (g_intern_static_string ("NcmHOAAVar"), values);
    g_once_init_leave (&static_g_enum_type_id, g_enum_type_id);
  }

  return static_g_enum_type_id;
}
/* enumerations from "math/ncm_lh_ratio1d.h" */
GType
ncm_lh_ratio1d_root_get_type (void)
{
  static gsize static_g_enum_type_id = 0;

  if (g_once_init_enter (&static_g_enum_type_id))
  {
    static const GEnumValue values[] = {
      { NCM_LH_RATIO1D_ROOT_BRACKET, "NCM_LH_RATIO1D_ROOT_BRACKET", "bracket" },
      { NCM_LH_RATIO1D_ROOT_NUMDIFF, "NCM_LH_RATIO1D_ROOT_NUMDIFF", "numdiff" },
      { 0, NULL, NULL }
    };
    GType g_enum_type_id =
      g_enum_register_static (g_intern_static_string ("NcmLHRatio1dRoot"), values);
    g_once_init_leave (&static_g_enum_type_id, g_enum_type_id);
  }

  return static_g_enum_type_id;
}
/* enumerations from "math/ncm_lh_ratio2d.h" */
GType
ncm_lh_ratio2d_root_get_type (void)
{
  static gsize static_g_enum_type_id = 0;

  if (g_once_init_enter (&static_g_enum_type_id))
  {
    static const GEnumValue values[] = {
      { NCM_LH_RATIO2D_ROOT_BRACKET, "NCM_LH_RATIO2D_ROOT_BRACKET", "bracket" },
      { NCM_LH_RATIO2D_ROOT_NUMDIFF, "NCM_LH_RATIO2D_ROOT_NUMDIFF", "numdiff" },
      { 0, NULL, NULL }
    };
    GType g_enum_type_id =
      g_enum_register_static (g_intern_static_string ("NcmLHRatio2dRoot"), values);
    g_once_init_leave (&static_g_enum_type_id, g_enum_type_id);
  }

  return static_g_enum_type_id;
}
/* enumerations from "math/ncm_matrix.h" */
GType
ncm_matrix_internal_get_type (void)
{
  static gsize static_g_enum_type_id = 0;

  if (g_once_init_enter (&static_g_enum_type_id))
  {
    static const GEnumValue values[] = {
      { NCM_MATRIX_SLICE, "NCM_MATRIX_SLICE", "slice" },
      { NCM_MATRIX_GSL_MATRIX, "NCM_MATRIX_GSL_MATRIX", "gsl-matrix" },
      { NCM_MATRIX_MALLOC, "NCM_MATRIX_MALLOC", "malloc" },
      { NCM_MATRIX_GARRAY, "NCM_MATRIX_GARRAY", "garray" },
      { NCM_MATRIX_DERIVED, "NCM_MATRIX_DERIVED", "derived" },
      { 0, NULL, NULL }
    };
    GType g_enum_type_id =
      g_enum_register_static (g_intern_static_string ("NcmMatrixInternal"), values);
    g_once_init_leave (&static_g_enum_type_id, g_enum_type_id);
  }

  return static_g_enum_type_id;
}
/* enumerations from "math/ncm_model_funnel.h" */
GType
ncm_model_funnel_sparams_get_type (void)
{
  static gsize static_g_enum_type_id = 0;

  if (g_once_init_enter (&static_g_enum_type_id))
  {
    static const GEnumValue values[] = {
      { NCM_MODEL_FUNNEL_NU, "NCM_MODEL_FUNNEL_NU", "nu" },
      { 0, NULL, NULL }
    };
    GType g_enum_type_id =
      g_enum_register_static (g_intern_static_string ("NcmModelFunnelSParams"), values);
    g_once_init_leave (&static_g_enum_type_id, g_enum_type_id);
  }

  return static_g_enum_type_id;
}
GType
ncm_model_funnel_vparams_get_type (void)
{
  static gsize static_g_enum_type_id = 0;

  if (g_once_init_enter (&static_g_enum_type_id))
  {
    static const GEnumValue values[] = {
      { NCM_MODEL_FUNNEL_X, "NCM_MODEL_FUNNEL_X", "x" },
      { 0, NULL, NULL }
    };
    GType g_enum_type_id =
      g_enum_register_static (g_intern_static_string ("NcmModelFunnelVParams"), values);
    g_once_init_leave (&static_g_enum_type_id, g_enum_type_id);
  }

  return static_g_enum_type_id;
}
/* enumerations from "math/ncm_model_mvnd.h" */
GType
ncm_model_mvndv_params_get_type (void)
{
  static gsize static_g_enum_type_id = 0;

  if (g_once_init_enter (&static_g_enum_type_id))
  {
    static const GEnumValue values[] = {
      { NCM_MODEL_MVND_MEAN, "NCM_MODEL_MVND_MEAN", "mean" },
      { 0, NULL, NULL }
    };
    GType g_enum_type_id =
      g_enum_register_static (g_intern_static_string ("NcmModelMVNDVParams"), values);
    g_once_init_leave (&static_g_enum_type_id, g_enum_type_id);
  }

  return static_g_enum_type_id;
}
/* enumerations from "math/ncm_model_rosenbrock.h" */
GType
ncm_model_rosenbrock_sparams_get_type (void)
{
  static gsize static_g_enum_type_id = 0;

  if (g_once_init_enter (&static_g_enum_type_id))
  {
    static const GEnumValue values[] = {
      { NCM_MODEL_ROSENBROCK_X1, "NCM_MODEL_ROSENBROCK_X1", "x1" },
      { NCM_MODEL_ROSENBROCK_X2, "NCM_MODEL_ROSENBROCK_X2", "x2" },
      { 0, NULL, NULL }
    };
    GType g_enum_type_id =
      g_enum_register_static (g_intern_static_string ("NcmModelRosenbrockSParams"), values);
    g_once_init_leave (&static_g_enum_type_id, g_enum_type_id);
  }

  return static_g_enum_type_id;
}
/* enumerations from "math/ncm_mpi_job.h" */
GType
ncm_mpi_job_ctrl_msg_get_type (void)
{
  static gsize static_g_enum_type_id = 0;

  if (g_once_init_enter (&static_g_enum_type_id))
  {
    static const GEnumValue values[] = {
      { NCM_MPI_CTRL_SLAVE_INIT, "NCM_MPI_CTRL_SLAVE_INIT", "init" },
      { NCM_MPI_CTRL_SLAVE_FREE, "NCM_MPI_CTRL_SLAVE_FREE", "free" },
      { NCM_MPI_CTRL_SLAVE_KILL, "NCM_MPI_CTRL_SLAVE_KILL", "kill" },
      { NCM_MPI_CTRL_SLAVE_WORK, "NCM_MPI_CTRL_SLAVE_WORK", "work" },
      { 0, NULL, NULL }
    };
    GType g_enum_type_id =
      g_enum_register_static (g_intern_static_string ("NcmMPIJobCtrlMsg"), values);
    g_once_init_leave (&static_g_enum_type_id, g_enum_type_id);
  }

  return static_g_enum_type_id;
}
GType
ncm_mpi_job_ctrl_tag_get_type (void)
{
  static gsize static_g_enum_type_id = 0;

  if (g_once_init_enter (&static_g_enum_type_id))
  {
    static const GEnumValue values[] = {
      { NCM_MPI_CTRL_TAG_CMD, "NCM_MPI_CTRL_TAG_CMD", "cmd" },
      { NCM_MPI_CTRL_TAG_JOB, "NCM_MPI_CTRL_TAG_JOB", "job" },
      { NCM_MPI_CTRL_TAG_WORK_INPUT, "NCM_MPI_CTRL_TAG_WORK_INPUT", "work-input" },
      { NCM_MPI_CTRL_TAG_WORK_RETURN, "NCM_MPI_CTRL_TAG_WORK_RETURN", "work-return" },
      { 0, NULL, NULL }
    };
    GType g_enum_type_id =
      g_enum_register_static (g_intern_static_string ("NcmMPIJobCtrlTag"), values);
    g_once_init_leave (&static_g_enum_type_id, g_enum_type_id);
  }

  return static_g_enum_type_id;
}
/* enumerations from "math/ncm_mset_catalog.h" */
GType
ncm_mset_catalog_sync_get_type (void)
{
  static gsize static_g_enum_type_id = 0;

  if (g_once_init_enter (&static_g_enum_type_id))
  {
    static const GEnumValue values[] = {
      { NCM_MSET_CATALOG_SYNC_DISABLE, "NCM_MSET_CATALOG_SYNC_DISABLE", "disable" },
      { NCM_MSET_CATALOG_SYNC_AUTO, "NCM_MSET_CATALOG_SYNC_AUTO", "auto" },
      { NCM_MSET_CATALOG_SYNC_TIMED, "NCM_MSET_CATALOG_SYNC_TIMED", "timed" },
      { 0, NULL, NULL }
    };
    GType g_enum_type_id =
      g_enum_register_static (g_intern_static_string ("NcmMSetCatalogSync"), values);
    g_once_init_leave (&static_g_enum_type_id, g_enum_type_id);
  }

  return static_g_enum_type_id;
}
GType
ncm_mset_catalog_trim_type_get_type (void)
{
  static gsize static_g_flags_type_id = 0;

  if (g_once_init_enter (&static_g_flags_type_id))
  {
    static const GFlagsValue values[] = {
      { NCM_MSET_CATALOG_TRIM_TYPE_ESS, "NCM_MSET_CATALOG_TRIM_TYPE_ESS", "ess" },
      { NCM_MSET_CATALOG_TRIM_TYPE_HEIDEL, "NCM_MSET_CATALOG_TRIM_TYPE_HEIDEL", "heidel" },
      { NCM_MSET_CATALOG_TRIM_TYPE_CK, "NCM_MSET_CATALOG_TRIM_TYPE_CK", "ck" },
      { NCM_MSET_CATALOG_TRIM_TYPE_ALL, "NCM_MSET_CATALOG_TRIM_TYPE_ALL", "all" },
      { 0, NULL, NULL }
    };
    GType g_flags_type_id =
      g_flags_register_static (g_intern_static_string ("NcmMSetCatalogTrimType"), values);
    g_once_init_leave (&static_g_flags_type_id, g_flags_type_id);
  }

  return static_g_flags_type_id;
}
GType
ncm_mset_catalog_post_norm_method_get_type (void)
{
  static gsize static_g_enum_type_id = 0;

  if (g_once_init_enter (&static_g_enum_type_id))
  {
    static const GEnumValue values[] = {
      { NCM_MSET_CATALOG_POST_LNNORM_METHOD_HYPERBOX, "NCM_MSET_CATALOG_POST_LNNORM_METHOD_HYPERBOX", "hyperbox" },
      { NCM_MSET_CATALOG_POST_LNNORM_METHOD_HYPERBOX_BS, "NCM_MSET_CATALOG_POST_LNNORM_METHOD_HYPERBOX_BS", "hyperbox-bs" },
      { NCM_MSET_CATALOG_POST_LNNORM_METHOD_ELIPSOID, "NCM_MSET_CATALOG_POST_LNNORM_METHOD_ELIPSOID", "elipsoid" },
      { 0, NULL, NULL }
    };
    GType g_enum_type_id =
      g_enum_register_static (g_intern_static_string ("NcmMSetCatalogPostNormMethod"), values);
    g_once_init_leave (&static_g_enum_type_id, g_enum_type_id);
  }

  return static_g_enum_type_id;
}
GType
ncm_mset_catalog_tau_method_get_type (void)
{
  static gsize static_g_enum_type_id = 0;

  if (g_once_init_enter (&static_g_enum_type_id))
  {
    static const GEnumValue values[] = {
      { NCM_MSET_CATALOG_TAU_METHOD_ACOR, "NCM_MSET_CATALOG_TAU_METHOD_ACOR", "acor" },
      { NCM_MSET_CATALOG_TAU_METHOD_AR_MODEL, "NCM_MSET_CATALOG_TAU_METHOD_AR_MODEL", "ar-model" },
      { 0, NULL, NULL }
    };
    GType g_enum_type_id =
      g_enum_register_static (g_intern_static_string ("NcmMSetCatalogTauMethod"), values);
    g_once_init_leave (&static_g_enum_type_id, g_enum_type_id);
  }

  return static_g_enum_type_id;
}
/* enumerations from "math/ncm_mset_trans_kern_cat.h" */
GType
ncm_mset_trans_kern_cat_sampling_get_type (void)
{
  static gsize static_g_enum_type_id = 0;

  if (g_once_init_enter (&static_g_enum_type_id))
  {
    static const GEnumValue values[] = {
      { NCM_MSET_TRANS_KERN_CAT_SAMPLING_CHOOSE, "NCM_MSET_TRANS_KERN_CAT_SAMPLING_CHOOSE", "choose" },
      { NCM_MSET_TRANS_KERN_CAT_SAMPLING_RBF_INTERP, "NCM_MSET_TRANS_KERN_CAT_SAMPLING_RBF_INTERP", "rbf-interp" },
      { NCM_MSET_TRANS_KERN_CAT_SAMPLING_KDE, "NCM_MSET_TRANS_KERN_CAT_SAMPLING_KDE", "kde" },
      { 0, NULL, NULL }
    };
    GType g_enum_type_id =
      g_enum_register_static (g_intern_static_string ("NcmMSetTransKernCatSampling"), values);
    g_once_init_leave (&static_g_enum_type_id, g_enum_type_id);
  }

  return static_g_enum_type_id;
}
/* enumerations from "math/ncm_nnls.h" */
GType
ncm_nnls_umethod_get_type (void)
{
  static gsize static_g_enum_type_id = 0;

  if (g_once_init_enter (&static_g_enum_type_id))
  {
    static const GEnumValue values[] = {
      { NCM_NNLS_UMETHOD_NORMAL, "NCM_NNLS_UMETHOD_NORMAL", "normal" },
      { NCM_NNLS_UMETHOD_NORMAL_LU, "NCM_NNLS_UMETHOD_NORMAL_LU", "normal-lu" },
      { NCM_NNLS_UMETHOD_QR, "NCM_NNLS_UMETHOD_QR", "qr" },
      { NCM_NNLS_UMETHOD_DGELSD, "NCM_NNLS_UMETHOD_DGELSD", "dgelsd" },
      { NCM_NNLS_UMETHOD_GSL, "NCM_NNLS_UMETHOD_GSL", "gsl" },
      { 0, NULL, NULL }
    };
    GType g_enum_type_id =
      g_enum_register_static (g_intern_static_string ("NcmNNLSUMethod"), values);
    g_once_init_leave (&static_g_enum_type_id, g_enum_type_id);
  }

  return static_g_enum_type_id;
}
/* enumerations from "math/ncm_ode_eval.h" */
GType
ncm_ode_eval_return_get_type (void)
{
  static gsize static_g_enum_type_id = 0;

  if (g_once_init_enter (&static_g_enum_type_id))
  {
    static const GEnumValue values[] = {
      { NCM_ODE_EVAL_RETURN_SUCCESS, "NCM_ODE_EVAL_RETURN_SUCCESS", "success" },
      { 0, NULL, NULL }
    };
    GType g_enum_type_id =
      g_enum_register_static (g_intern_static_string ("NcmODEEvalReturn"), values);
    g_once_init_leave (&static_g_enum_type_id, g_enum_type_id);
  }

  return static_g_enum_type_id;
}
/* enumerations from "math/ncm_powspec_filter.h" */
GType
ncm_powspec_filter_type_get_type (void)
{
  static gsize static_g_enum_type_id = 0;

  if (g_once_init_enter (&static_g_enum_type_id))
  {
    static const GEnumValue values[] = {
      { NCM_POWSPEC_FILTER_TYPE_TOPHAT, "NCM_POWSPEC_FILTER_TYPE_TOPHAT", "tophat" },
      { NCM_POWSPEC_FILTER_TYPE_GAUSS, "NCM_POWSPEC_FILTER_TYPE_GAUSS", "gauss" },
      { 0, NULL, NULL }
    };
    GType g_enum_type_id =
      g_enum_register_static (g_intern_static_string ("NcmPowspecFilterType"), values);
    g_once_init_leave (&static_g_enum_type_id, g_enum_type_id);
  }

  return static_g_enum_type_id;
}
/* enumerations from "math/ncm_serialize.h" */
GType
ncm_serialize_opt_get_type (void)
{
  static gsize static_g_flags_type_id = 0;

  if (g_once_init_enter (&static_g_flags_type_id))
  {
    static const GFlagsValue values[] = {
      { NCM_SERIALIZE_OPT_NONE, "NCM_SERIALIZE_OPT_NONE", "none" },
      { NCM_SERIALIZE_OPT_AUTOSAVE_SER, "NCM_SERIALIZE_OPT_AUTOSAVE_SER", "autosave-ser" },
      { NCM_SERIALIZE_OPT_AUTONAME_SER, "NCM_SERIALIZE_OPT_AUTONAME_SER", "autoname-ser" },
      { NCM_SERIALIZE_OPT_CLEAN_DUP, "NCM_SERIALIZE_OPT_CLEAN_DUP", "clean-dup" },
      { 0, NULL, NULL }
    };
    GType g_flags_type_id =
      g_flags_register_static (g_intern_static_string ("NcmSerializeOpt"), values);
    g_once_init_leave (&static_g_flags_type_id, g_flags_type_id);
  }

  return static_g_flags_type_id;
}
/* enumerations from "math/ncm_sparam.h" */
GType
ncm_param_type_get_type (void)
{
  static gsize static_g_enum_type_id = 0;

  if (g_once_init_enter (&static_g_enum_type_id))
  {
    static const GEnumValue values[] = {
      { NCM_PARAM_TYPE_FREE, "NCM_PARAM_TYPE_FREE", "free" },
      { NCM_PARAM_TYPE_FIXED, "NCM_PARAM_TYPE_FIXED", "fixed" },
      { 0, NULL, NULL }
    };
    GType g_enum_type_id =
      g_enum_register_static (g_intern_static_string ("NcmParamType"), values);
    g_once_init_leave (&static_g_enum_type_id, g_enum_type_id);
  }

  return static_g_enum_type_id;
}
/* enumerations from "math/ncm_sphere_map.h" */
GType
ncm_sphere_map_order_get_type (void)
{
  static gsize static_g_enum_type_id = 0;

  if (g_once_init_enter (&static_g_enum_type_id))
  {
    static const GEnumValue values[] = {
      { NCM_SPHERE_MAP_ORDER_NEST, "NCM_SPHERE_MAP_ORDER_NEST", "nest" },
      { NCM_SPHERE_MAP_ORDER_RING, "NCM_SPHERE_MAP_ORDER_RING", "ring" },
      { 0, NULL, NULL }
    };
    GType g_enum_type_id =
      g_enum_register_static (g_intern_static_string ("NcmSphereMapOrder"), values);
    g_once_init_leave (&static_g_enum_type_id, g_enum_type_id);
  }

  return static_g_enum_type_id;
}
GType
ncm_sphere_map_coord_sys_get_type (void)
{
  static gsize static_g_enum_type_id = 0;

  if (g_once_init_enter (&static_g_enum_type_id))
  {
    static const GEnumValue values[] = {
      { NCM_SPHERE_MAP_COORD_SYS_GALACTIC, "NCM_SPHERE_MAP_COORD_SYS_GALACTIC", "galactic" },
      { NCM_SPHERE_MAP_COORD_SYS_ECLIPTIC, "NCM_SPHERE_MAP_COORD_SYS_ECLIPTIC", "ecliptic" },
      { NCM_SPHERE_MAP_COORD_SYS_CELESTIAL, "NCM_SPHERE_MAP_COORD_SYS_CELESTIAL", "celestial" },
      { 0, NULL, NULL }
    };
    GType g_enum_type_id =
      g_enum_register_static (g_intern_static_string ("NcmSphereMapCoordSys"), values);
    g_once_init_leave (&static_g_enum_type_id, g_enum_type_id);
  }

  return static_g_enum_type_id;
}
/* enumerations from "math/ncm_spline_func.h" */
GType
ncm_spline_func_type_get_type (void)
{
  static gsize static_g_enum_type_id = 0;

  if (g_once_init_enter (&static_g_enum_type_id))
  {
    static const GEnumValue values[] = {
      { NCM_SPLINE_FUNCTION_4POINTS, "NCM_SPLINE_FUNCTION_4POINTS", "function-4points" },
      { NCM_SPLINE_FUNCTION_SPLINE, "NCM_SPLINE_FUNCTION_SPLINE", "function-spline" },
      { NCM_SPLINE_FUNCTION_SPLINE_LNKNOT, "NCM_SPLINE_FUNCTION_SPLINE_LNKNOT", "function-spline-lnknot" },
      { NCM_SPLINE_FUNCTION_SPLINE_SINHKNOT, "NCM_SPLINE_FUNCTION_SPLINE_SINHKNOT", "function-spline-sinhknot" },
      { NCM_SPLINE_FUNC_GRID_LINEAR, "NCM_SPLINE_FUNC_GRID_LINEAR", "func-grid-linear" },
      { NCM_SPLINE_FUNC_GRID_LOG, "NCM_SPLINE_FUNC_GRID_LOG", "func-grid-log" },
      { 0, NULL, NULL }
    };
    GType g_enum_type_id =
      g_enum_register_static (g_intern_static_string ("NcmSplineFuncType"), values);
    g_once_init_leave (&static_g_enum_type_id, g_enum_type_id);
  }

  return static_g_enum_type_id;
}
/* enumerations from "math/ncm_spline_func_test.h" */
GType
ncm_spline_func_test_type_get_type (void)
{
  static gsize static_g_enum_type_id = 0;

  if (g_once_init_enter (&static_g_enum_type_id))
  {
    static const GEnumValue values[] = {
      { NCM_SPLINE_FUNC_TEST_TYPE_POLYNOMIAL, "NCM_SPLINE_FUNC_TEST_TYPE_POLYNOMIAL", "polynomial" },
      { NCM_SPLINE_FUNC_TEST_TYPE_COSINE, "NCM_SPLINE_FUNC_TEST_TYPE_COSINE", "cosine" },
      { NCM_SPLINE_FUNC_TEST_TYPE_RBF, "NCM_SPLINE_FUNC_TEST_TYPE_RBF", "rbf" },
      { NCM_SPLINE_FUNC_TEST_TYPE_USER, "NCM_SPLINE_FUNC_TEST_TYPE_USER", "user" },
      { 0, NULL, NULL }
    };
    GType g_enum_type_id =
      g_enum_register_static (g_intern_static_string ("NcmSplineFuncTestType"), values);
    g_once_init_leave (&static_g_enum_type_id, g_enum_type_id);
  }

  return static_g_enum_type_id;
}
GType
ncm_spline_func_test_type_pdf_get_type (void)
{
  static gsize static_g_enum_type_id = 0;

  if (g_once_init_enter (&static_g_enum_type_id))
  {
    static const GEnumValue values[] = {
      { NCM_SPLINE_FUNC_TEST_TYPE_PDF_FLAT, "NCM_SPLINE_FUNC_TEST_TYPE_PDF_FLAT", "flat" },
      { NCM_SPLINE_FUNC_TEST_TYPE_PDF_NORMAL, "NCM_SPLINE_FUNC_TEST_TYPE_PDF_NORMAL", "normal" },
      { 0, NULL, NULL }
    };
    GType g_enum_type_id =
      g_enum_register_static (g_intern_static_string ("NcmSplineFuncTestTypePDF"), values);
    g_once_init_leave (&static_g_enum_type_id, g_enum_type_id);
  }

  return static_g_enum_type_id;
}
/* enumerations from "math/ncm_spline_gsl.h" */
GType
ncm_spline_gsl_type_get_type (void)
{
  static gsize static_g_enum_type_id = 0;

  if (g_once_init_enter (&static_g_enum_type_id))
  {
    static const GEnumValue values[] = {
      { NCM_SPLINE_GSL_LINEAR, "NCM_SPLINE_GSL_LINEAR", "linear" },
      { NCM_SPLINE_GSL_POLYNOMIAL, "NCM_SPLINE_GSL_POLYNOMIAL", "polynomial" },
      { NCM_SPLINE_GSL_CSPLINE, "NCM_SPLINE_GSL_CSPLINE", "cspline" },
      { NCM_SPLINE_GSL_CSPLINE_PERIODIC, "NCM_SPLINE_GSL_CSPLINE_PERIODIC", "cspline-periodic" },
      { NCM_SPLINE_GSL_AKIMA, "NCM_SPLINE_GSL_AKIMA", "akima" },
      { NCM_SPLINE_GSL_AKIMA_PERIODIC, "NCM_SPLINE_GSL_AKIMA_PERIODIC", "akima-periodic" },
      { 0, NULL, NULL }
    };
    GType g_enum_type_id =
      g_enum_register_static (g_intern_static_string ("NcmSplineGslType"), values);
    g_once_init_leave (&static_g_enum_type_id, g_enum_type_id);
  }

  return static_g_enum_type_id;
}
/* enumerations from "math/ncm_spline_rbf.h" */
GType
ncm_spline_rbf_type_get_type (void)
{
  static gsize static_g_enum_type_id = 0;

  if (g_once_init_enter (&static_g_enum_type_id))
  {
    static const GEnumValue values[] = {
      { NCM_SPLINE_RBF_TYPE_POSDEF_GAUSS, "NCM_SPLINE_RBF_TYPE_POSDEF_GAUSS", "posdef-gauss" },
      { NCM_SPLINE_RBF_TYPE_GAUSS, "NCM_SPLINE_RBF_TYPE_GAUSS", "gauss" },
      { 0, NULL, NULL }
    };
    GType g_enum_type_id =
      g_enum_register_static (g_intern_static_string ("NcmSplineRBFType"), values);
    g_once_init_leave (&static_g_enum_type_id, g_enum_type_id);
  }

  return static_g_enum_type_id;
}
/* enumerations from "math/ncm_stats_dist.h" */
GType
ncm_stats_dist_cv_get_type (void)
{
  static gsize static_g_enum_type_id = 0;

  if (g_once_init_enter (&static_g_enum_type_id))
  {
    static const GEnumValue values[] = {
      { NCM_STATS_DIST_CV_NONE, "NCM_STATS_DIST_CV_NONE", "none" },
      { NCM_STATS_DIST_CV_SPLIT, "NCM_STATS_DIST_CV_SPLIT", "split" },
      { 0, NULL, NULL }
    };
    GType g_enum_type_id =
      g_enum_register_static (g_intern_static_string ("NcmStatsDistCV"), values);
    g_once_init_leave (&static_g_enum_type_id, g_enum_type_id);
  }

  return static_g_enum_type_id;
}
/* enumerations from "math/ncm_stats_dist1d_epdf.h" */
GType
ncm_stats_dist1d_epdf_bw_get_type (void)
{
  static gsize static_g_enum_type_id = 0;

  if (g_once_init_enter (&static_g_enum_type_id))
  {
    static const GEnumValue values[] = {
      { NCM_STATS_DIST1D_EPDF_BW_FIXED, "NCM_STATS_DIST1D_EPDF_BW_FIXED", "fixed" },
      { NCM_STATS_DIST1D_EPDF_BW_RoT, "NCM_STATS_DIST1D_EPDF_BW_RoT", "rot" },
      { NCM_STATS_DIST1D_EPDF_BW_AUTO, "NCM_STATS_DIST1D_EPDF_BW_AUTO", "auto" },
      { 0, NULL, NULL }
    };
    GType g_enum_type_id =
      g_enum_register_static (g_intern_static_string ("NcmStatsDist1dEPDFBw"), values);
    g_once_init_leave (&static_g_enum_type_id, g_enum_type_id);
  }

  return static_g_enum_type_id;
}
/* enumerations from "math/ncm_stats_dist_kde.h" */
GType
ncm_stats_dist_kde_cov_type_get_type (void)
{
  static gsize static_g_enum_type_id = 0;

  if (g_once_init_enter (&static_g_enum_type_id))
  {
    static const GEnumValue values[] = {
      { NCM_STATS_DIST_KDE_COV_TYPE_SAMPLE, "NCM_STATS_DIST_KDE_COV_TYPE_SAMPLE", "sample" },
      { NCM_STATS_DIST_KDE_COV_TYPE_FIXED, "NCM_STATS_DIST_KDE_COV_TYPE_FIXED", "fixed" },
      { NCM_STATS_DIST_KDE_COV_TYPE_ROBUST_DIAG, "NCM_STATS_DIST_KDE_COV_TYPE_ROBUST_DIAG", "robust-diag" },
      { NCM_STATS_DIST_KDE_COV_TYPE_ROBUST, "NCM_STATS_DIST_KDE_COV_TYPE_ROBUST", "robust" },
      { 0, NULL, NULL }
    };
    GType g_enum_type_id =
      g_enum_register_static (g_intern_static_string ("NcmStatsDistKDECovType"), values);
    g_once_init_leave (&static_g_enum_type_id, g_enum_type_id);
  }

  return static_g_enum_type_id;
}
/* enumerations from "math/ncm_stats_vec.h" */
GType
ncm_stats_vec_type_get_type (void)
{
  static gsize static_g_enum_type_id = 0;

  if (g_once_init_enter (&static_g_enum_type_id))
  {
    static const GEnumValue values[] = {
      { NCM_STATS_VEC_MEAN, "NCM_STATS_VEC_MEAN", "mean" },
      { NCM_STATS_VEC_VAR, "NCM_STATS_VEC_VAR", "var" },
      { NCM_STATS_VEC_COV, "NCM_STATS_VEC_COV", "cov" },
      { 0, NULL, NULL }
    };
    GType g_enum_type_id =
      g_enum_register_static (g_intern_static_string ("NcmStatsVecType"), values);
    g_once_init_leave (&static_g_enum_type_id, g_enum_type_id);
  }

  return static_g_enum_type_id;
}
GType
ncm_stats_vec_ar_type_get_type (void)
{
  static gsize static_g_enum_type_id = 0;

  if (g_once_init_enter (&static_g_enum_type_id))
  {
    static const GEnumValue values[] = {
      { NCM_STATS_VEC_AR_NONE, "NCM_STATS_VEC_AR_NONE", "none" },
      { NCM_STATS_VEC_AR_FPE, "NCM_STATS_VEC_AR_FPE", "fpe" },
      { NCM_STATS_VEC_AR_AIC, "NCM_STATS_VEC_AR_AIC", "aic" },
      { NCM_STATS_VEC_AR_AICC, "NCM_STATS_VEC_AR_AICC", "aicc" },
      { 0, NULL, NULL }
    };
    GType g_enum_type_id =
      g_enum_register_static (g_intern_static_string ("NcmStatsVecARType"), values);
    g_once_init_leave (&static_g_enum_type_id, g_enum_type_id);
  }

  return static_g_enum_type_id;
}
/* enumerations from "math/ncm_vector.h" */
GType
ncm_vector_internal_get_type (void)
{
  static gsize static_g_enum_type_id = 0;

  if (g_once_init_enter (&static_g_enum_type_id))
  {
    static const GEnumValue values[] = {
      { NCM_VECTOR_SLICE, "NCM_VECTOR_SLICE", "slice" },
      { NCM_VECTOR_GSL_VECTOR, "NCM_VECTOR_GSL_VECTOR", "gsl-vector" },
      { NCM_VECTOR_MALLOC, "NCM_VECTOR_MALLOC", "malloc" },
      { NCM_VECTOR_ARRAY, "NCM_VECTOR_ARRAY", "array" },
      { NCM_VECTOR_DERIVED, "NCM_VECTOR_DERIVED", "derived" },
      { 0, NULL, NULL }
    };
    GType g_enum_type_id =
      g_enum_register_static (g_intern_static_string ("NcmVectorInternal"), values);
    g_once_init_leave (&static_g_enum_type_id, g_enum_type_id);
  }

  return static_g_enum_type_id;
}

/* Generated data ends here */

