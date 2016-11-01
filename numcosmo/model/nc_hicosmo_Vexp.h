/***************************************************************************
 *            nc_hicosmo_Vexp.h
 *
 *  Fri October 28 13:27:25 2016
 *  Copyright  2016  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2016 <sandro@isoftware.com.br>
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

#ifndef _NC_HICOSMO_VEXP_H_
#define _NC_HICOSMO_VEXP_H_

#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/nc_hicosmo.h>
#include <nvector/nvector_serial.h>

G_BEGIN_DECLS

#define NC_TYPE_HICOSMO_VEXP             (nc_hicosmo_Vexp_get_type ())
#define NC_HICOSMO_VEXP(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_HICOSMO_VEXP, NcHICosmoVexp))
#define NC_HICOSMO_VEXP_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_HICOSMO_VEXP, NcHICosmoVexpClass))
#define NC_IS_HICOSMO_VEXP(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_HICOSMO_VEXP))
#define NC_IS_HICOSMO_VEXP_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_HICOSMO_VEXP))
#define NC_HICOSMO_VEXP_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_HICOSMO_VEXP, NcHICosmoVexpClass))

typedef struct _NcHICosmoVexpClass NcHICosmoVexpClass;
typedef struct _NcHICosmoVexp NcHICosmoVexp;

struct _NcHICosmoVexpClass
{
  /*< private >*/
  NcHICosmoClass parent_class;
};

/**
 * NcHICosmoVexpParams:
 * @NC_HICOSMO_VEXP_H0:      FIXME
 * @NC_HICOSMO_VEXP_OMEGA_C: FIXME
 * @NC_HICOSMO_VEXP_SIGMA:   FIXME
 * @NC_HICOSMO_VEXP_C:       FIXME
 * @NC_HICOSMO_VEXP_ALPHA0:  FIXME
 *
 * FIXME
 * 
 */
typedef enum _NcHICosmoVexpParams
{
  NC_HICOSMO_VEXP_H0 = 0,
  NC_HICOSMO_VEXP_OMEGA_C,
  NC_HICOSMO_VEXP_SIGMA_PHI,
  NC_HICOSMO_VEXP_D_PHI,
  NC_HICOSMO_VEXP_ALPHA_0,    /*< private >*/
  NC_HICOSMO_VEXP_SPARAM_LEN, /*< skip >*/
} NcHICosmoVexpParams;

struct _NcHICosmoVexp
{
  /*< private >*/
  NcHICosmo parent_instance;
  gpointer cvode_qt;
  gpointer cvode_cl;
  gboolean qt_init;
  gboolean cl_init;
  N_Vector y_qt;
  N_Vector ydot_qt;
  N_Vector y_cl;
};

GType nc_hicosmo_Vexp_get_type (void) G_GNUC_CONST;

NcHICosmoVexp *nc_hicosmo_Vexp_new (void);

#define NC_HICOSMO_VEXP_DEFAULT_H0 (70.0)
#define NC_HICOSMO_VEXP_DEFAULT_OMEGA_C (0.25)
#define NC_HICOSMO_VEXP_DEFAULT_SIGMA_PHI (0.4)
#define NC_HICOSMO_VEXP_DEFAULT_D_PHI (-0.3)
#define NC_HICOSMO_VEXP_DEFAULT_ALPHA_0 (0.1)

G_END_DECLS

#endif /* _NC_HICOSMO_VEXP_H_ */
