/* -*- Mode: C; indent-tabs-mode: nil; c-basic-offset: 2; tab-width: 2 -*-  */

/***************************************************************************
 *            ncm_stats_dist_private.h
 *
 *  Thu July 22 15:12:38 2021
 *  Copyright  2021  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_stats_dist_private.h
 * Copyright (C) 2021 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

#ifndef _NCM_STATS_DIST_PRIVATE_H_
#define _NCM_STATS_DIST_PRIVATE_H_

#include <glib.h>
#include "math/ncm_stats_dist.h"
#include "math/ncm_nnls.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_multimin.h>
#endif /* NUMCOSMO_GIR_SCAN */

G_BEGIN_DECLS

struct _NcmStatsDistPrivate
{
  /*< private >*/
  NcmStatsDistKernel *kernel;
  GPtrArray *sample_array;
  NcmVector *weights;
  NcmVector *wcum;
  gboolean wcum_ready;
  gboolean cv_ready;
  gboolean print_fit;
  gdouble over_smooth;
  NcmStatsDistCV cv_type;
  gboolean use_threads;
  gdouble split_frac;
  gdouble min_m2lnp;
  gdouble max_m2lnp;
  gdouble href;
  gdouble rnorm;
  guint n_obs;
  guint n_kernels;
  guint alloc_n_obs;
  guint alloc_n_kernels;
  gboolean alloc_subs;
  guint d;
  GArray *sampling;
  NcmNNLS *nnls;
  NcmMatrix *IM;
  NcmMatrix *sub_IM;
  NcmVector *sub_x;
  NcmVector *f;
  NcmVector *f1;
  gdouble *levmar_workz;
  guint levmar_n;
  gsl_multimin_fminimizer *fmin;
  GArray *m2lnp_sort;
  NcmRNG *rng;
};

G_END_DECLS

#endif /* _NCM_STATS_DIST_PRIVATE_H_ */

