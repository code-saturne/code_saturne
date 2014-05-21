!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2014 EDF S.A.
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

subroutine lagesp &
!================

 ( nvar   , nscal  ,                                              &
   nbpmax , nvp    , nvp1   , nvep   , nivep  ,                   &
   ntersl , nvlsta , nvisbr ,                                     &
   itepa  ,                                                       &
   dt     , rtpa   , propce ,                                     &
   ettp   , ettpa  , tepa   , statis , stativ ,                   &
   taup   , tlag   , piil   ,                                     &
   tsuf   , tsup   , bx     , tsfext ,                            &
   vagaus , gradpr , gradvf , brgaus , terbru , romp   , auxl2 ,  &
   vislen  )

!===============================================================================
! Purpose:
! ----------

!   Subroutine of the Lagrangian particle-tracking module :
!   ------------------------------------------------------

!   Integration of particle equations of motion :
!
!   * Standard Model : First order  -> call of subroutine lages1
!                      Second order -> call of subroutine lages2
!
!   * Deposition submodel (Guingo & Minier, 2008) if needed
!
!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! nbpmax           ! e  ! <-- ! nombre max de particulies autorise             !
! nvp              ! e  ! <-- ! nombre de variables particulaires              !
! nvp1             ! e  ! <-- ! nvp sans position, vfluide, vpart              !
! nvep             ! e  ! <-- ! nombre info particulaires (reels)              !
! nivep            ! e  ! <-- ! nombre info particulaires (entiers)            !
! ntersl           ! e  ! <-- ! nbr termes sources de couplage retour          !
! nvlsta           ! e  ! <-- ! nombre de var statistiques lagrangien          !
! nvisbr           ! e  ! <-- ! nombre de statistiques aux frontieres          !
! itepa            ! te ! <-- ! info particulaires (entiers)                   !
! (nbpmax,nivep    !    !     !   (cellule de la particule,...)                !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtpa             ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at previous time step)                       !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! ettp             ! tr ! <-- ! tableaux des variables liees                   !
!  (nbpmax,nvp)    !    !     !   aux particules etape courante                !
! ettpa            ! tr ! <-- ! tableaux des variables liees                   !
!  (nbpmax,nvp)    !    !     !   aux particules etape precedente              !
! tepa             ! tr ! <-- ! info particulaires (reels)                     !
! (nbpmax,nvep)    !    !     !   (poids statistiques,...)                     !
! statis           ! tr ! <-- ! cumul pour les moyennes des                    !
!(ncelet,nvlsta    !    !     !   statistiques volumiques                      !
! stativ           ! tr ! <-- ! cumul pour les variances des                   !
!(ncelet,          !    !     !    statistiques volumiques                     !
!   nvlsta-1)      !    !     !                                                !
! taup(nbpmax)     ! tr ! <-- ! temps caracteristique dynamique                !
! tlag(nbpmax)     ! tr ! <-- ! temps caracteristique fluide                   !
! piil(nbpmax,3    ! tr ! --> ! terme dans l'integration des eds up            !
! tsup(nbpmax,3    ! tr ! --> ! prediction 1er sous-pas pour                   !
!                  !    !     !   la vitesse des particules                    !
! tsuf(nbpmax,3    ! tr ! --> ! prediction 1er sous-pas pour                   !
!                  !    !     !   la vitesse du fluide vu                      !
! bx(nbpmax,3,2    ! tr ! <-- ! caracteristiques de la turbulence              !
! tsfext(nbpmax    ! tr ! <-- ! infos pour le couplage retour                  !
! vagaus           ! tr ! <-- ! variables aleatoires gaussiennes               !
!(nbpmax,nvgaus    !    !     !                                                !
! gradpr(3,ncel)   ! tr ! <-- ! gradient de pression                           !
! gradvf(3,3,ncel) ! tr ! <-- ! gradient de la vitesse du fluide               !
! romp             ! tr ! --- ! masse volumique des particules                 !
! auxl2            ! tr ! --- ! tableau de travail                             !
!    (nbpmax,7)    !    !     !                                                !
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
use cstphy
use cstnum
use optcal
use entsor
use lagpar
use lagran
use ppppar
use ppthch
use mesh

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal
integer          nbpmax , nvp    , nvp1   , nvep  , nivep
integer          ntersl , nvlsta , nvisbr

integer          itepa(nbpmax,nivep)

double precision dt(ncelet) , rtpa(ncelet,nflown:nvar)
double precision propce(ncelet,*)
double precision ettp(nbpmax,nvp) , ettpa(nbpmax,nvp)
double precision tepa(nbpmax,nvep)
double precision statis(ncelet,*),stativ(ncelet,*)
double precision taup(nbpmax) , tlag(nbpmax,3)
double precision piil(nbpmax,3) , bx(nbpmax,3,2)
double precision tsuf(nbpmax,3) , tsup(nbpmax,3)
double precision tsfext(nbpmax)
double precision dlgeo(nfabor,ngeol)
double precision vagaus(nbpmax,*)
double precision gradpr(3,ncelet) , gradvf(3,3,ncelet)
double precision brgaus(nbpmax,*) , terbru(nbpmax)
double precision romp(nbpmax) , auxl2(nbpmax,7)
double precision vislen(nfabor)

! Local variables

integer          ip
integer          iifacl

double precision d3 , aa

double precision, allocatable, dimension(:,:) :: fextla

!===============================================================================

!===============================================================================
! 1.  Initialization
!===============================================================================
! Initialize variables to avoid compiler warnings

iifacl = 0

! Computation of particle density

aa = 6.d0 / pi
do ip = 1,nbpart
  if ( itepa(ip,jisor).gt.0 ) then
    d3 = ettp(ip,jdp) * ettp(ip,jdp) * ettp(ip,jdp)
    romp(ip) = aa * ettp(ip,jmp) / d3
  endif
enddo

!===============================================================================
! 2.  Management of user external force fields
!===============================================================================

! Allocate a temporay array
allocate(fextla(nbpmax,3))

do ip = 1, nbpmax
  fextla(ip,1) = 0.d0
  fextla(ip,2) = 0.d0
  fextla(ip,3) = 0.d0
enddo

call uslafe                                                       &
!==========
 ( nvar   , nscal  ,                                              &
   nbpmax , nvp    , nvp1   , nvep   , nivep  ,                   &
   ntersl , nvlsta , nvisbr ,                                     &
   itepa  ,                                                       &
   dt     ,                                                       &
   ettp   , ettpa  , tepa   , statis , stativ ,                   &
   taup   , tlag   , piil   ,                                     &
   tsuf   , tsup   , bx     , tsfext ,                            &
   vagaus , gradpr , gradvf ,                                     &
   romp   , fextla )


!===============================================================================
! 4.  First order
!===============================================================================

if (nordre.eq.1) then

!=============================================================================
! 4.1 If no deposition sub-model is activated, call of subroutine lages1
!     for every particle
!=============================================================================

  if (idepst.le.0) then

  call lages1                                                     &
  !==========
   ( nbpmax , nvp    , nvep   , nivep  ,                          &
     itepa  ,                                                     &
     rtpa   , propce ,                                            &
     ettp   , ettpa  , tepa   ,                                   &
     taup   , tlag   , piil   ,                                   &
     bx     , vagaus , gradpr , romp   ,                          &
     brgaus , terbru , fextla )


!=============================================================================
! 4.2 Management of the deposition submodel
!=============================================================================

  else

     call lagdep                                                  &
    !==========
   ( nbpmax ,                                                     &
     rtpa   , propce ,                                            &
     taup   , tlag   , piil   ,                                   &
     bx     , vagaus , gradpr , romp   ,                          &
     fextla , vislen)

  endif

!===============================================================================
! 5.  Second order
!===============================================================================

else

  call lages2                                                     &
  !==========
   ( nbpmax , nvp    , nvep   , nivep  ,                          &
     itepa  ,                                                     &
     rtpa   , propce ,                                            &
     ettp   , ettpa  , tepa   ,                                   &
     taup   , tlag   , piil   ,                                   &
     tsuf   , tsup   , bx     , tsfext , vagaus ,                 &
     auxl2  , gradpr ,                                            &
     romp   , brgaus , terbru , fextla )

endif

! Free memory
deallocate(fextla)

!===============================================================================

end subroutine
