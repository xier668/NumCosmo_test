/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 2; tab-width: 2 -*-  */

/***************************************************************************
 *            nc_galaxy_wl_ellipticity_gauss.h
 *
 *  Mon July 27 11:13:56 2020
 *  Copyright  2020  Sandro Dias Pinto Vitenti & Mariana Penna Lima
 *  <sandro@isoftware.com.br>, <pennalima@gmail.com>
 ****************************************************************************/
/*
 * nc_galaxy_wl_ellipticity_gauss.h
 * Copyright (C) 2020 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
 * Copyright (C) 2020 Mariana Penna Lima <pennalima@gmail.com>
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

#ifndef _NC_GALAXY_WL_ELLIPTICITY_GAUSS_H_
#define _NC_GALAXY_WL_ELLIPTICITY_GAUSS_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_rng.h>
#include <numcosmo/math/ncm_matrix.h>
#include <numcosmo/lss/nc_galaxy_wl_dist.h>

G_BEGIN_DECLS

#define NC_TYPE_GALAXY_WL_ELLIPTICITY_GAUSS             (nc_galaxy_wl_ellipticity_gauss_get_type ())
#define NC_GALAXY_WL_ELLIPTICITY_GAUSS(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_GALAXY_WL_ELLIPTICITY_GAUSS, NcGalaxyWLEllipticityGauss))
#define NC_GALAXY_WL_ELLIPTICITY_GAUSS_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_GALAXY_WL_ELLIPTICITY_GAUSS, NcGalaxyWLEllipticityGaussClass))
#define NC_IS_GALAXY_WL_ELLIPTICITY_GAUSS(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_GALAXY_WL_ELLIPTICITY_GAUSS))
#define NC_IS_GALAXY_WL_ELLIPTICITY_GAUSS_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_GALAXY_WL_ELLIPTICITY_GAUSS))
#define NC_GALAXY_WL_ELLIPTICITY_GAUSS_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_GALAXY_WL_ELLIPTICITY_GAUSS, NcGalaxyWLEllipticityGaussClass))

typedef struct _NcGalaxyWLEllipticityGaussClass NcGalaxyWLEllipticityGaussClass;
typedef struct _NcGalaxyWLEllipticityGauss NcGalaxyWLEllipticityGauss;
typedef struct _NcGalaxyWLEllipticityGaussPrivate NcGalaxyWLEllipticityGaussPrivate;

struct _NcGalaxyWLEllipticityGaussClass
{
  /*< private >*/
  NcGalaxyWLDistClass parent_class;
};

struct _NcGalaxyWLEllipticityGauss
{
  /*< private >*/
  NcGalaxyWLDist parent_instance;
  NcGalaxyWLEllipticityGaussPrivate *priv;
};

/**
 * NcGalaxyWLEllipticityGaussPos:
 * @NC_GALAXY_WL_ELLIPTICITY_GAUSS_POS_ANG: FIXME
 * @NC_GALAXY_WL_ELLIPTICITY_GAUSS_POS_R: FIXME
 *
 * FIXME
 *
 */
typedef enum _NcGalaxyWLEllipticityGaussPos
{
  NC_GALAXY_WL_ELLIPTICITY_GAUSS_POS_ANG,
  NC_GALAXY_WL_ELLIPTICITY_GAUSS_POS_R,
  /* < private > */
  NC_GALAXY_WL_ELLIPTICITY_GAUSS_POS_LEN, /*< skip >*/
} NcGalaxyWLEllipticityGaussPos;

GType nc_galaxy_wl_ellipticity_gauss_get_type (void) G_GNUC_CONST;

NcGalaxyWLEllipticityGauss *nc_galaxy_wl_ellipticity_gauss_new (NcGalaxyWLEllipticityGaussPos pos);
NcGalaxyWLEllipticityGauss *nc_galaxy_wl_ellipticity_gauss_ref (NcGalaxyWLEllipticityGauss *gegauss);

void nc_galaxy_wl_ellipticity_gauss_free (NcGalaxyWLEllipticityGauss *gegauss);
void nc_galaxy_wl_ellipticity_gauss_clear (NcGalaxyWLEllipticityGauss **gegauss);

void nc_galaxy_wl_ellipticity_gauss_set_pos (NcGalaxyWLEllipticityGauss *gegauss, NcGalaxyWLEllipticityGaussPos pos);
NcGalaxyWLEllipticityGaussPos nc_galaxy_wl_ellipticity_gauss_get_pos (NcGalaxyWLEllipticityGauss *gegauss);

void nc_galaxy_wl_ellipticity_gauss_set_obs (NcGalaxyWLEllipticityGauss *gegauss, NcmMatrix *obs);
NcmMatrix *nc_galaxy_wl_ellipticity_gauss_peek_obs (NcGalaxyWLEllipticityGauss *gegauss);

G_END_DECLS

#endif /* _NC_GALAXY_WL_ELLIPTICITY_GAUSS_H_ */

