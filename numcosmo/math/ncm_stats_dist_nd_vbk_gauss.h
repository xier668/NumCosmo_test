/* -*- Mode: C; indent-tabs-mode: nil; c-basic-offset: 2; tab-width: 2 -*-  */
/***************************************************************************
 *            ncm_stats_dist_nd_vbk_gauss.h
 *
 *  Wed November 07 17:41:38 2018
 *  Copyright  2018  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * ncm_stats_dist_nd_vbk_gauss.h
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

#ifndef _NCM_STATS_DIST_ND_VBK_GAUSS_H_
#define _NCM_STATS_DIST_ND_VBK_GAUSS_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_stats_dist_nd_vbk.h>
#include <numcosmo/math/ncm_vector.h>

G_BEGIN_DECLS

#define NCM_TYPE_STATS_DIST_ND_VBK_GAUSS             (ncm_stats_dist_nd_vbk_gauss_get_type ())
#define NCM_STATS_DIST_ND_VBK_GAUSS(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_STATS_DIST_ND_VBK_GAUSS, NcmStatsDistNdVBKGauss))
#define NCM_STATS_DIST_ND_VBK_GAUSS_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_STATS_DIST_ND_VBK_GAUSS, NcmStatsDistNdVBKGaussClass))
#define NCM_IS_STATS_DIST_ND_VBK_GAUSS(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_STATS_DIST_ND_VBK_GAUSS))
#define NCM_IS_STATS_DIST_ND_VBK_GAUSS_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_STATS_DIST_ND_VBK_GAUSS))
#define NCM_STATS_DIST_ND_VBK_GAUSS_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_STATS_DIST_ND_VBK_GAUSS, NcmStatsDistNdVBKGaussClass))

typedef struct _NcmStatsDistNdVBKGaussClass NcmStatsDistNdVBKGaussClass;
typedef struct _NcmStatsDistNdVBKGauss NcmStatsDistNdVBKGauss;
typedef struct _NcmStatsDistNdVBKGaussPrivate NcmStatsDistNdVBKGaussPrivate;


struct _NcmStatsDistNdVBKGaussClass
{
  NcmStatsDistNdVBKClass parent_class;
};

struct _NcmStatsDistNdVBKGauss
{
  NcmStatsDistNdVBK parent_instance;

  NcmStatsDistNdVBKGaussPrivate *priv;
};

GType ncm_stats_dist_nd_vbk_gauss_get_type (void) G_GNUC_CONST;

NcmStatsDistNdVBKGauss *ncm_stats_dist_nd_vbk_gauss_new (const guint dim, const NcmStatsDistNdVBKCV cv_type);

NcmStatsDistNdVBKGauss *ncm_stats_dist_nd_vbk_gauss_ref (NcmStatsDistNdVBKGauss *dndg);
void ncm_stats_dist_nd_vbk_gauss_free (NcmStatsDistNdVBKGauss *dndg);
void ncm_stats_dist_nd_vbk_gauss_clear (NcmStatsDistNdVBKGauss **dndg);

G_END_DECLS

#endif /* _NCM_STATS_DIST_ND_VBK_GAUSS_H_ */
