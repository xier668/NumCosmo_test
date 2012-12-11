/***************************************************************************
 *            nc_data_bao.h
 *
 *  Thu November 22 20:41:23 2012
 *  Copyright  2012  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) 2012 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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

#ifndef _NC_DATA_BAO_H_
#define _NC_DATA_BAO_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/math/ncm_data.h>
#include <numcosmo/nc_distance.h>

G_BEGIN_DECLS

/**
 * NcDataBaoId:
 * @NC_DATA_BAO_A_EISENSTEIN: FIXME
 * @NC_DATA_BAO_DV_EISENSTEIN: FIXME
 * @NC_DATA_BAO_DVDV_PERCIVAL: FIXME
 * @NC_DATA_BAO_RDV_PERCIVAL: FIXME
 *
 * FIXME
 */
typedef enum _NcDataBaoId
{
  NC_DATA_BAO_A_EISENSTEIN = 0,
  NC_DATA_BAO_DV_EISENSTEIN,
  NC_DATA_BAO_DVDV_PERCIVAL,
  NC_DATA_BAO_RDV_PERCIVAL,     /*< private >*/
  NC_DATA_BAO_NSAMPLES,         /*< skip >*/
} NcDataBaoId;

NcmData *nc_data_bao_new (NcDistance *dist, NcDataBaoId id);

G_END_DECLS

#endif /* _NC_DATA_BAO_H_ */