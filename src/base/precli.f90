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

subroutine precli &
!================

 ( nvar   , nscal  ,                                              &
   icodcl ,                                                       &
   propfb ,                                                       &
   rcodcl )

!===============================================================================
! FONCTION :
! --------

!    PREPARATION DU REMPLISSAGE DES CONDITIONS AUX LIMITES

! On dispose des types de faces de bord au pas de temps
!   precedent (sauf au premier pas de temps, ou les tableaux
!   ITYPFB et ITRIFB n'ont pas ete renseignes)

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! icodcl           ! te ! --> ! code de condition limites aux faces            !
!  (nfabor,nvar    !    !     !  de bord                                       !
!                  !    !     ! = 1   -> dirichlet                             !
!                  !    !     ! = 3   -> densite de flux                       !
!                  !    !     ! = 4   -> glissemt et u.n=0 (vitesse)           !
!                  !    !     ! = 5   -> frottemt et u.n=0 (vitesse)           !
!                  !    !     ! = 6   -> rugosite et u.n=0 (vitesse)           !
!                  !    !     ! = 9   -> entree/sortie libre (vitesse          !
!                  !    !     !  entrante eventuelle     bloquee               !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! rcodcl           ! tr ! --> ! valeur des conditions aux limites              !
!  (nfabor,nvar    !    !     !  aux faces de bord                             !
!                  !    !     ! rcodcl(1) = valeur du dirichlet                !
!                  !    !     ! rcodcl(2) = valeur du coef. d'echange          !
!                  !    !     !  ext. (infinie si pas d'echange)               !
!                  !    !     ! rcodcl(3) = valeur de la densite de            !
!                  !    !     !  flux (negatif si gain) w/m2 ou                !
!                  !    !     !  hauteur de rugosite (m) si icodcl=6           !
!                  !    !     ! pour les vitesses (vistl+visct)*gradu          !
!                  !    !     ! pour la pression             dt*gradp          !
!                  !    !     ! pour les scalaires                             !
!                  !    !     !        cp*(viscls+visct/sigmas)*gradt          !
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
use cstnum
use entsor
use pointe
use albase
use ppppar
use ppthch
use ppincl
use mesh

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal

integer          icodcl(nfabor,nvarcl)

double precision propfb(nfabor,*)
double precision rcodcl(nfabor,nvarcl,3)

! Local variables

integer          ifac, ivar, iscal, iut, ivt, iwt

!===============================================================================
!===============================================================================
! 1.  INITIALISATIONS
!===============================================================================


!===============================================================================
! 2.  INITIALISATION DES CONDITIONS LIMITES ET TYPE DE FACES DE BORD
!===============================================================================

!      ICODCL = 0 INDIQUE QUE LA CL N'A PAS ETE RENSEIGNEE
!      ITYPFB = 0 INDIQUE QUE LA CL N'A PAS ETE RENSEIGNEE

!      RINFIN : VALEUR INFINIE

do ifac = 1, nfabor
  itypfb(ifac) = 0
enddo

! Pour toutes les variables, on initialise RCODCL(1)a RINFIN
! Cette valeur sera reinitialisee a zero dans typecl.F

do ivar = 1, nvar
  do ifac = 1, nfabor
    icodcl(ifac,ivar)   = 0
    rcodcl(ifac,ivar,1) = rinfin
    rcodcl(ifac,ivar,2) = rinfin
    rcodcl(ifac,ivar,3) = 0.d0
  enddo
enddo

! Default value for turbulent fluxes
do iscal = 1, nscal
  if (ityturt(iscal).eq.3) then
    iut = nvar + 3*(ifltur(iscal) - 1) + 1
    ivt = nvar + 3*(ifltur(iscal) - 1) + 2
    iwt = nvar + 3*(ifltur(iscal) - 1) + 3
    do ifac = 1, nfabor
      icodcl(ifac,iut)   = 0
      rcodcl(ifac,iut,1) = rinfin
      rcodcl(ifac,iut,2) = rinfin
      rcodcl(ifac,iut,3) = 0.d0
      icodcl(ifac,ivt)   = 0
      rcodcl(ifac,ivt,1) = rinfin
      rcodcl(ifac,ivt,2) = rinfin
      rcodcl(ifac,ivt,3) = 0.d0
      icodcl(ifac,iwt)   = 0
      rcodcl(ifac,iwt,1) = rinfin
      rcodcl(ifac,iwt,2) = rinfin
      rcodcl(ifac,iwt,3) = 0.d0
    enddo
  endif
enddo

! En ALE, on initialise aussi le tableau IALTYB
if (iale.eq.1) then
  do ifac = 1, nfabor
    ialtyb(ifac) = 0
  enddo
endif

! POUR LES PHYSIQUES PARTICULIERES

if (ippmod(iphpar).ge.1) then
  call ppprcl(nvar, izfppp, propfb, rcodcl)
  !==========
endif

!----
! End
!----

return
end subroutine
