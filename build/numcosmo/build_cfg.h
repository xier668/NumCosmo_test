/***************************************************************************
 *            build_cfg.h.in
 *
 *  Wed Apr 21 15:19:05 2010
 *  Copyright  2010  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <vitenti@uel.br>
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
/*
 * File autogenerated by configure script, dont mess with it
 */

#ifndef _NCM_BUILD_CFG_H_
#define _NCM_BUILD_CFG_H_

/**
 * SECTION:build_cfg
 * @title: NumCosmo building configuration
 * @short_description: Macros describing the current building configuration
 *
 * This section shows the macros controling which prerequisite/option was
 * used when compiling NumCosmo.
 * 
 */

#define NUMCOSMO_HAVE_FFTW3 1
#define NUMCOSMO_HAVE_FFTW3F 1
#define NUMCOSMO_HAVE_CFITSIO 1

#define NUMCOSMO_HAVE_NLOPT 1


#define NUMCOSMO_HAVE_INLINE 1
#define NCM_INLINE static inline

#ifndef NUMCOSMO_GIR_SCAN
#include <cblas.h>
#define __MKL_CBLAS_H__
#define __GSL_CBLAS_H__
#define NCM_BLAS_NOT_TYPEDEFED 1
#include <numcosmo/math/ncm_gsl_blas_types.h>

#endif /* NUMCOSMO_GIR_SCAN */

#define NUMCOSMO_MAJOR_VERSION 0
#define NUMCOSMO_MINOR_VERSION 16
#define NUMCOSMO_MICRO_VERSION 0
#define NUMCOSMO_VERSION "0.16.0"

#endif /* _NCM_BUILD_CFG_H_ */
