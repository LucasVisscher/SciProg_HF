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

   This file is the header file of Gen1Int used by users.

   2012-05-01, Bin Gao:
   * first version
*/

#if !defined(GEN1INT_H)
#define GEN1INT_H

/* QMatrix library for matrix operations */
#include "qmatrix.h"

/* total geometric derivatives */
#include "geom_total.h"

/* one-electron operator */
#include "one_prop.h"

/* type of Gaussian type orbitals (GTOs)
   (1) NON_LAO: non London atomic orbital
   (2) LONDON: London atomic orbital (LAO)
   (3) ROT_LAO: rotational LAO */
typedef enum {NON_LAO,LONDON,ROT_LAO} GTOType;

/* maximum integer for the non-atomic center */
#define MAX_IDX_NON 0

#endif
