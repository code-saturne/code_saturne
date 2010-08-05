!-------------------------------------------------------------------------------

!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2009 EDF S.A., France

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

subroutine fuphy2 &
!================

 ( ncelet , ncel   ,                                              &
   rtp    , propce )

!===============================================================================
! FONCTION :
! --------

! CALCUL DES PROPRIETES PHYSIQUES DE LA PHASE DISPERSEE
! VALEURS CELLULES
! ----------------

!   FRACTION MASSIQUE DE LIQUIDE
!     ET CLIPPING EVENTUELS
!   DIAMETRE
!   MASSE VOLUMIQUE
!     ET CLIPPING EVENTUELS

! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! ncelet           ! i  ! <-- ! number of extended (real + ghost) cells        !
! ncel             ! i  ! <-- ! number of cells                                !
! rtp              ! tr ! <-- ! variables de calcul au centre des              !
! (ncelet,*)       !    !     !    cellules (instant courant)                  !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use optcal
use cstphy
use entsor
use cstnum
use parall
use ppppar
use ppthch
use coincl
use cpincl
use fuincl
use ppincl

!===============================================================================

implicit none

! Arguments

integer          ncelet , ncel

double precision rtp(ncelet,*) , propce(ncelet,*)

! Local variables

integer          iel
integer          n1     , n2     , ipcdia , icla
double precision xng    ,  d1s3
double precision rhofol , diam3m , diam3x

!===============================================================================

!===============================================================================
! 1. INITIALISATIONS
!===============================================================================

d1s3 = 1.d0/3.d0

!===============================================================================
! 2. CALCUL POUR CHAQUE CLASSE
!    DE LA MASSE VOLUMIQUE DU FOL
!    DE LA FRACTION MASSIQUE DE FOL
!    DU DIAMETRE DU COKE
!===============================================================================

do icla = 1, nclafu

  n1 = 0
  n2 = 0
  diam3m =  1.d0
  diam3x =  0.d0

  ipcdia = ipproc(idiam3(icla))

  do iel = 1, ncel

!        Masse Volumique

    propce(iel,ipproc(irom3(icla))) = rho0fl

    yfol   = rtp(iel,isca(iyfol(icla)))
    xng    = rtp(iel,isca(ing(icla)))
    rhofol = propce(iel,ipproc(irom3(icla)))

! --- Calcul du diametre de la particule

!     Attention !
!     Diam en mm
!     Xng en Ggouttes (en milliards de gouttes/kg de mélange)
!     ça marche mais il faudrait un facteur mis à UN pour s'y retrouver

    if ( yfol .gt. epsifl .and. (xng*yfol) .gt. 0.d0) then

      propce(iel,ipcdia) = ((yfol / rhofol)                       &
                          /(pi/6.d0 * xng) ) ** d1s3

      if ( propce(iel,ipcdia) .gt. dinifl(icla) ) then
        n1 = n1+1
        diam3x = max(diam3x,propce(iel,ipcdia))
        propce(iel,ipcdia) = dinifl(icla)
      endif

      if ( propce(iel,ipcdia) .lt. diniin(icla) ) then
        n2 = n2+1
        diam3m = min(diam3m,propce(iel,ipcdia))
        propce(iel,ipcdia) = diniin(icla)
      endif

    else
      propce(iel,ipcdia) = dinifl(icla)
    endif

  enddo

  if (irangp.ge.0) then

    call parcpt (n1)
    !==========
    call parcpt (n2)
    !==========

    call parmax (diam3x)
    !==========
    call parmin (diam3m)
    !==========
  endif

  if ( n1 .gt. 0 ) then
    write(nfecra,1001) icla, n1, diam3x
  endif
  if ( n2 .gt. 0 ) then
    write(nfecra,1002)  icla, n2, diam3m
  endif

enddo

!----
! FORMATS
!----

 1001 format(/,1X,' CLIPPING EN MAX DU DIAMETRE CLASSE :',I2,           &
       /,10X,' Nombre de points : ',I8,                           &
       /,10X,' Valeur Max       : ',G15.7)
 1002 format(/,1X,' CLIPPING EN MIN DU DIAMETRE CLASSE :',I2,           &
       /,10X,' Nombre de points : ',I8,                           &
       /,10X,' Valeur Min       : ',G15.7)


!----
! FIN
!----

return
end subroutine
