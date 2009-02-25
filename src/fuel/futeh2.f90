!-------------------------------------------------------------------------------

!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2008 EDF S.A., France

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

                  subroutine futeh2                               &
!================

 ( ncelet , ncel   , nrtuse ,                                     &
   rtp    , propce , rtuser )

!===============================================================================
! FONCTION :
! --------
! CALCUL DE LA TEMPERATURE DES PARTICULES
!  EN FONCTION DE L'ENTHALPIE DU FOL ET DES CONCENTRATIONS

! Arguments
!__________________.____._____.________________________________________________.
!    nom           !type!mode !                   role                         !
!__________________!____!_____!________________________________________________!
! ncelet           ! e  ! <-- ! nombre d'elements halo compris                 !
! ncel             ! e  ! <-- ! nombre d'elements actifs                       !
! rtp              ! tr ! <-- ! variables de calcul au centre des              !
! (ncelet,*)       !    !     !    cellules (instant courant)                  !
! propce           ! tr ! <-- ! proprietes physiques au centre des             !
! (ncelet,*)       !    !     !    cellules                                    !
! eh0              ! tr ! <-- ! tableau reel de travail                        !
! eh1              ! tr ! <-- ! tableau reel de travail                        !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHAMNUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

implicit none

!==============================================================================
!     DONNEES EN COMMON
!==============================================================================

include "paramx.h"
include "numvar.h"
include "optcal.h"
include "cstphy.h"
include "cstnum.h"
include "entsor.h"
include "ppppar.h"
include "ppthch.h"
include "coincl.h"
include "cpincl.h"
include "fuincl.h"
include "ppincl.h"

!===============================================================================

! Arguments

integer          ncelet, ncel , nrtuse
double precision rtp(ncelet,*), propce(ncelet,*)
double precision rtuser(nrtuse)

! VARIABLES LOCALES

integer          icel , icla
integer          ipcte1 , ipcte2
integer          mode
double precision eh2, xsolid(2)
double precision mkfini,diamgt
double precision masgut,mfgout,mkgout,rhofol

!===============================================================================
! RQ IMPORTANTE : On suppose pour l'instant que H2 = H02 + CP2(T2-TREF)



!===============================================================================
! 1. CALCULS PRELIMINAIRES
!===============================================================================

! --- Initialisation des tableaux

ipcte1 = ipproc(itemp1)

do icla = 1, nclafu

  ipcte2 = ipproc(itemp3(icla))

  do icel = 1, ncel

! --- Initialisation de T2 a T20

    propce(icel,ipcte2) = 373.d0
! 20/09/05 Pour l'instant TinFol en dur, un jour on recuperera dans USFUCL
  enddo
enddo

!===============================================================================
! 2. CALCUL DE LA TEMPERATURE DES PARTICULES
!===============================================================================

mode = 1

do icla = 1, nclafu

  ipcte2 = ipproc(itemp3(icla))

  mkfini = rho0fl*pi/6.d0*dinikf(icla)**3

  do icel = 1, ncel

    rhofol = propce(icel,ipproc(irom3(icla)))
    diamgt = propce(icel,ipproc(idiam3(icla)))
    masgut = rho0fl*pi/6.d0*diamgt**3
    if (diamgt.le.dinikf(icla)) then
      mkgout = masgut
    else
      mkgout = mkfini
    endif
    mfgout = masgut - mkgout
    xsolid(1) = 1.d0-fkc
    xsolid(2) = fkc
    if(masgut.gt.zero) then
      xsolid(1) = mfgout / masgut
      xsolid(2) = mkgout / masgut
    endif
    xsolid(1) = min(1.d0,max(0.d0,xsolid(1)))
    xsolid(2) = min(1.d0,max(0.d0,xsolid(2)))

    if ( rtp(icel,isca(iyfol(icla))) .gt. (3.d3*epsifl) ) then
      eh2 =  rtp(icel,isca(ihlf(icla)))                           &
            /rtp(icel,isca(iyfol(icla)))
      call futhp2                                                 &
!           ============
    ( mode , eh2 , xsolid , propce(icel,ipcte2) )
    else
      propce(icel,ipcte2) = propce(icel,ipcte1)
    endif

  enddo

enddo

!----
! FIN
!----

return
end
