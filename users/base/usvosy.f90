!-------------------------------------------------------------------------------

!VERS


!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2011 EDF S.A., France

!     contact: saturne-support@edf.fr

!     The Code_Saturne Kernel is free software; you can redistribute it
!     and/or modify it under the terms of the GNU General Public License
!     as published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.

!     The Code_Saturne Kernel is distributed in the hope that it will be
!     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
!     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.

!     You should have received a copy of the GNU General Public License
!     along with the Code_Saturne Kernel; if not, write to the
!     Free Software Foundation, Inc.,
!     51 Franklin St, Fifth Floor,
!     Boston, MA  02110-1301  USA

!-------------------------------------------------------------------------------

subroutine usvosy &
!================

 ( nvar   , nscal  , inbcou , ncecpl ,                            &
   iscalt ,                                                       &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   lcecpl , hvol )

!===============================================================================
! Purpose:
! -------

!    User subroutine.

!    Compute a volume exchange coefficient for SYRTHES coupling

!
! Usage
! -----
! The routine is called in cpvosy() for each volume coupling
! therefore it is necessary to test the value of coupling number to separate
! the treatments of the different couplings
!
! Up to now temperature is the only scalar managed for volume couplings.
!

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! inbcou           ! i  ! <-- ! SYRTHES coupling number                        !
! ncecpl           ! i  ! <-- ! number of cells implied for this coupling      !
! iscalt           ! i  ! <-- ! index number of the temperature scalar         !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp              ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (current time step)                           !
! rtpa             ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (preceding time step)                         !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! lcecpl(ncecpl)   ! ri ! <-- ! list of coupled cells                          !
! hvol(ncecpl)     ! ra ! --> ! volume exchange coefficient to compute         !
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use entsor
use optcal
use cstphy
use parall
use period
use mesh

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal
integer          ncecpl
integer          iscalt , inbcou

integer          lcecpl(ncecpl)

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision hvol(ncecpl), tfluid(ncecpl)

! Local variables

character*80     chaine
integer          ivar, iiscvr, iel, iloc, iutile, ivart
integer          ipcrom, ipcvsl, ipcvis, ipccp

double precision cp, mu, lambda, rho, uloc, L, sexcvo
double precision nu, re, pr
double precision hcorr, hvol_cst, lambda_over_cp

!===============================================================================

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================

if(1.eq.1) return

!===============================================================================
! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END


!===============================================================================
! 1. Initialization
!===============================================================================

! --- Index number of the cell properties the propce array
ipcrom = ipproc(irom)
ipcvis = ipproc(iviscl)

if (icp.gt.0) then
   ipccp = ipproc(icp)
else
   ipccp = 0
endif

if (ivisls(iscalt).gt.0) then
   ipcvsl = ipproc(ivisls(iscalt))
else
   ipcvsl = 0
endif

ivart = isca(iscalt)

!===============================================================================
! 2. Example 1 of the computation of a volumic exchange coefficient
!
!    hvol(iel) = cst
!
!===============================================================================

! ----------------------------------------------

! It is quite frequent to forget to remove this example when it is
!  not needed. Therefore the following test is designed to prevent
!  any bad surprise.

iutile = 0

if(iutile.eq.0) return

hvol_cst = 1.0d6

do iloc = 1, ncecpl  ! Loop on coupled cells
   hvol(iloc) = hvol_cst
enddo

!===============================================================================
! 2. Example 2 of the computation of a volumic exchange coefficient
!
!    hvol(iel) =  hsurf(iel) * exchange_surface_by_unit_vol
!
!    with: hsurf = Nusselt * lambda / L
!
!    lambda is the thermal conductivity coefficient
!    L is a characteristic length
!
!    Nusselt is computed by means of the Colburn correlation
!
!    Nu = 0.023 * Re^(0.8) * Pr^(1/3)
!
!    Re is the Reynolds number and Pr is the Prandtl number
!
!===============================================================================

! ----------------------------------------------

! It is quite frequent to forget to remove this example when it is
!  not needed. Therefore the following test is designed to prevent
!  any bad surprise.

iutile = 0

if(iutile.eq.0) return

sexcvo = 36.18d0  ! Surface area where exchanges take place by unit of volume
L = 0.03d0        ! Characteristic length

! No test on the coupling number (inbcou). We assume that the same
! treatment is applied to all volume couplings

do iloc = 1, ncecpl  ! Loop on coupled cells

   iel = lcecpl(iloc)

   ! Get cell properties of the current element

   rho = propce(iel, ipcrom)
   mu = propce(iel, ipcvis)

   if (ipccp.gt.0) then
      cp = propce(iel, ipccp)
   else
      cp = cp0
   endif

   if (ipcvsl.gt.0) then ! lambda/Cp is variable
      lambda_over_cp = propce(iel, ipcvsl)
      lambda =  lambda_over_cp * cp
   else
      lambda_over_cp = visls0(iscalt)
      lambda = lambda_over_cp * cp
   endif

   ! Compute a local molecular Prandtl **(1/3)

   pr = mu / lambda_over_cp

   ! Compute a local Reynolds number

   uloc = sqrt(rtp(iel, iu)**2 + rtp(iel, iv)**2 + rtp(iel, iw)**2)
   re = max(uloc*rho*L/mu, 1.d0) ! To avoid division by zero

   ! Compute Nusselt number thanks to Colburn correlation

   nu = 0.023d0 * re**0.8d0 * pr**(1.d0/3.d0)
   hcorr = nu * lambda / L

   ! Compute hvol
   hvol(iloc) = hcorr * sexcvo

enddo

! ----------------------------------------------

!--------
! Formats
!--------

#if defined(_CS_LANG_FR)

 1000 format(                                                     &
'@                                                            ',/,&
'@ @@ ATTENTION : COUPLAGE VOLUMIQUE SYRTHES AVEC UN SCALAIRE ',/,&
'@      QUI EST DIFFERENT DE LA TEMPERATURE                   ',/,&
'@    =========                                               ',/,&
'@      OPTION NON VALIDE                                     ',/,&
'@                                                            ')

#else

 1000 format(                                                     &
'@                                                            ',/,&
'@ @@ WARNING: SYRTHES VOLUME COUPLING WITH A SCALAR          ',/,&
'@       DIFFERENT FROM TEMPERATURE                           ',/,&
'@    ========                                                ',/,&
'@      OPTION NOT POSSIBLE                                   ',/,&
'@                                                            ')

#endif

!----
! End
!----

return
end subroutine
