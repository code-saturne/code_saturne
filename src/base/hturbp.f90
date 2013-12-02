!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2013 EDF S.A.
!
! This program is free software; you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free Software
! Foundation; either version 2 of the License, or (at your option) any later
! version.
!
! This program is distributed in the hope that it will be useful, but WITHOUT
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
! FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
! details.
!
! You should have received a copy of the GNU General Public License along with
! this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
! Street, Fifth Floor, Boston, MA 02110-1301, USA.

!-------------------------------------------------------------------------------

!===============================================================================
! Function :
! --------

!> Compute the correction of the exchange coefficient between the fluid and
!> the wall for a turbulent flow.
!>
!> This is function of the dimensionless
!> distance to the wall \f$ y^+ = \dfrac{\centip \centf u_\star}{\nu}\f$.
!>
!> Then the return coefficient reads:
!> \f[
!> h_{tur} = Pr \dfrac{y^+}{T^+}
!> \f]
!>
!> This coefficient is computed thanks to a similarity model between
!> dynamic viscous sub-layer and themal sub-layer.
!>
!> \f$ T^+ \f$ is computed as follows:
!>
!> - For a laminar Prandtl number smaller than 0.1 (such as liquid metals),
!>   the standard model with two sub-layers (Prandtl-Taylor) is used.
!>
!> - For a laminar Prandtl number larger than 0.1 (such as liquids and gaz),
!>   a model with three sub-layers (Arpaci-Larsen) is used.
!>
!> The final exchange coefficient is:
!> \f[
!> h = \dfrac{K}{\centip \centf} h_{tur}
!> \f]
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     prl           laminar Prandtl number
!> \param[in]     prt           turbulent Prandtl number
!> \param[in]     ckarm         Von Karman constant
!> \param[in]     yplus         dimensionless distance to the wall
!> \param[in]     dplus         dimensionless distance for scalable
!>                              wall functions
!> \param[out]    htur          corrected exchange coefficient
!> \param[out]    yplim         value of the limit for \f$ y^+ \f$
!_______________________________________________________________________________

subroutine hturbp &
 ( prl    , prt    , ckarm  , yplus  , dplus, htur , yplim )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use optcal, only: iwallt
use cstnum

!===============================================================================

implicit none

! Arguments

double precision htur
double precision prl,ckarm,prt,yplus, dplus, yplim

! Local variables

double precision tplus
double precision beta2,a2
double precision yp2
double precision prlm1

!===============================================================================

!===============================================================================
! 1. Initializations
!===============================================================================

!===============================================================================

htur = max(yplus-dplus, epzero)/max(yplus, epzero)

prlm1 = 0.1d0

!===============================================================================
! 2. Compute htur for small Prandtl numbers
!===============================================================================

if (prl.le.prlm1) then
  yplim   = prt/(prl*ckarm)
  if (yplus .gt. yplim) then
    tplus = prl*yplim + prt/ckarm * log(yplus/yplim)
    htur = prl*(yplus-dplus)/tplus
  endif

!===============================================================================
! 3. Compute htur for the model with three sub-layers
!===============================================================================

else
  yp2   = ckarm*1000.d0/prt
  yp2   = sqrt(yp2)
  yplim   = (1000.d0/prl)**(1.d0/3.d0)

  a2 = 15.d0*(prl**(2.d0/3.d0))
  beta2 = a2 - 500.d0/ (yp2**2)

  if ((yplus.ge.yplim).and.(yplus.lt.yp2)) then
    tplus = a2 - 500.d0/(yplus*yplus)
    htur = prl*(yplus-dplus)/tplus
  endif

  if ((yplus.ge.yp2)) then
    tplus = beta2 + prt/ckarm*log(yplus/yp2)
    htur = prl*(yplus-dplus)/tplus
  endif

endif

!----
! End
!----

return

end subroutine
