/* gen1int: compute one-electron integrals using rotational London atomic-orbitals
   Copyright 2009-2012 Bin Gao, and Andreas Thorvaldsen

   gen1int is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   gen1int is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with gen1int. If not, see <http://www.gnu.org/licenses/>.

   This file is the header file of one-electron property integrals.

   2012-05-24, Bin Gao:
   * first version
*/

#if !defined(ONE_PROP_H)
#define ONE_PROP_H

/* QMatrix library for matrix operations */
#include "qmatrix.h"

/* One-electron operator */
typedef struct {
  QInt order_geo_bra;  /* order of partial geometric derivatives on bra center */
  QInt order_geo_ket;  /* order of partial geometric derivatives on ket center */
  QInt order_mag;      /* order of total magnetic derivatives */
  QInt order_mag_bra;  /* order of partial magnetic derivatives on bra center */
  QInt order_mag_ket;  /* order of partial magnetic derivatives on ket center */
  QInt order_ram;      /* order of total derivatives w.r.t. total rotational angular momentum (RAM) */
  QInt order_ram_bra;  /* order of partial derivatives on bra center w.r.t. total RAM */
  QInt order_ram_ket;  /* order of partial derivatives on ket center w.r.t. total RAM */
  GTOType gto_type;    /* type of GTOs */
  QInt add_sr;         /* scalar-relativistic (SR) correction */
  QInt add_so;         /* spin-orbit (SO) correction */
  QInt add_london;     /* transforms the operator by the LAO type gauge-including projector */
} OneProp;

/* functions of one-electron operator */
extern QErrorCode OnePropCreate
extern QErrorCode OnePropSetPartialGeom
extern QErrorCode OnePropSetMag
extern QErrorCode OnePropSetRAM
extern QErrorCode OnePropSetGTO
extern QErrorCode OnePropAddSR
extern QErrorCode OnePropAddSO
extern QErrorCode OnePropGetIntegral
extern QErrorCode OnePropGetFunction
extern QErrorCode OnePropView
extern QErrorCode OnePropDestroy

#endif
