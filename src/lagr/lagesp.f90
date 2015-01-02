!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2015 EDF S.A.
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
   ntersl , nvlsta , nvisbr ,                                     &
   dt     , propce ,                                              &
   statis , stativ , taup   , tlag   , piil   ,                   &
   bx     , tsfext ,                                              &
   gradpr , gradvf , terbru , vislen  )

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
! ntersl           ! e  ! <-- ! nbr termes sources de couplage retour          !
! nvlsta           ! e  ! <-- ! nombre de var statistiques lagrangien          !
! nvisbr           ! e  ! <-- ! nombre de statistiques aux frontieres          !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! statis           ! tr ! <-- ! cumul pour les moyennes des                    !
!(ncelet,nvlsta    !    !     !   statistiques volumiques                      !
! stativ           ! tr ! <-- ! cumul pour les variances des                   !
!(ncelet,          !    !     !    statistiques volumiques                     !
!   nvlsta-1)      !    !     !                                                !
! taup(nbpart)     ! tr ! <-- ! temps caracteristique dynamique                !
! tlag(nbpart)     ! tr ! <-- ! temps caracteristique fluide                   !
! piil(nbpart,3)   ! tr ! --> ! terme dans l'integration des eds up            !
! tsup(nbpart,3)   ! tr ! --> ! prediction 1er sous-pas pour                   !
!                  !    !     !   la vitesse des particules                    !
! tsuf(nbpart,3)   ! tr ! --> ! prediction 1er sous-pas pour                   !
!                  !    !     !   la vitesse du fluide vu                      !
! bx(nbpart,3,2)   ! tr ! <-- ! caracteristiques de la turbulence              !
! tsfext(nbpart)   ! tr ! <-- ! infos pour le couplage retour                  !
! gradpr(3,ncel)   ! tr ! <-- ! gradient de pression                           !
! gradvf(3,3,ncel) ! tr ! <-- ! gradient de la vitesse du fluide               !
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
integer          ntersl , nvlsta , nvisbr

double precision dt(ncelet)
double precision propce(ncelet,*)
double precision statis(ncelet,*),stativ(ncelet,*)
double precision taup(nbpart) , tlag(nbpart,3)
double precision piil(nbpart,3) , bx(nbpart,3,2)
double precision tsfext(nbpart)
double precision gradpr(3,ncelet) , gradvf(3,3,ncelet)
double precision terbru(nbpart)
double precision romp(nbpart)
double precision vislen(nfabor)

! Local variables

integer          ip, ivf
integer          iifacl

double precision d3 , aa

double precision, allocatable, dimension(:,:) :: fextla
double precision, allocatable, dimension(:,:) :: tsuf, tsup
double precision, allocatable, dimension(:,:) :: vagaus, brgaus

!===============================================================================

!===============================================================================
! 1.  Initialization
!===============================================================================
! Initialize variables to avoid compiler warnings

iifacl = 0

! Computation of particle density

aa = 6.d0 / pi
do ip = 1,nbpart
  if ( ipepa(jisor,ip).gt.0 ) then
    d3 = eptp(jdp,ip) * eptp(jdp,ip) * eptp(jdp,ip)
    romp(ip) = aa * eptp(jmp,ip) / d3
  endif
enddo

!===============================================================================
! 2.  Management of user external force fields
!===============================================================================

! Allocate temporay arrays
allocate(fextla(nbpart,3))
allocate(vagaus(nbpart,nvgaus))

! Random values

if (idistu.eq.1) then
  do ivf = 1,nvgaus
    call normalen(nbpart, vagaus(:,ivf))
  enddo
else
  do ivf = 1,nvgaus
    do ip = 1,nbpart
      vagaus(ip,ivf) = 0.d0
    enddo
  enddo
endif

! Brownian movement

if (lamvbr.eq.1) then
  allocate(brgaus(nbpart,nbrgau))
  do ivf = 1,nbrgau
    call normalen(nbpart, brgaus(:,ivf))
  enddo
endif

do ip = 1, nbpart
  fextla(ip,1) = 0.d0
  fextla(ip,2) = 0.d0
  fextla(ip,3) = 0.d0
enddo

allocate(tsuf(nbpart,3))
allocate(tsup(nbpart,3))

if (jtsuf(1).gt.0) then
  do ip = 1, nbpart
    tsup(ip,1) = pepa(jtsup(1),ip)
    tsup(ip,2) = pepa(jtsup(2),ip)
    tsup(ip,3) = pepa(jtsup(3),ip)
    tsuf(ip,1) = pepa(jtsuf(1),ip)
    tsuf(ip,2) = pepa(jtsuf(2),ip)
    tsuf(ip,3) = pepa(jtsuf(3),ip)
  enddo
else
  do ip = 1, nbpart
    tsup(ip,1) = 0.d0
    tsup(ip,2) = 0.d0
    tsup(ip,3) = 0.d0
    tsuf(ip,1) = 0.d0
    tsuf(ip,2) = 0.d0
    tsuf(ip,3) = 0.d0
  enddo
endif

call uslafe                                                       &
!==========
 ( nvar   , nscal  ,                                              &
   ntersl , nvlsta , nvisbr ,                                     &
   dt     ,                                                       &
   statis , stativ ,                                              &
   taup   , tlag   , piil   ,                                     &
   tsuf   , tsup   , bx     , tsfext ,                            &
   vagaus , gradpr , gradvf ,                                     &
   romp   , fextla )

deallocate(tsuf)
deallocate(tsup)

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
   ( propce ,                                                     &
     taup   , tlag   , piil   ,                                   &
     bx     , vagaus , gradpr , romp   ,                          &
     brgaus , terbru , fextla )


!=============================================================================
! 4.2 Management of the deposition submodel
!=============================================================================

  else

     call lagdep                                                  &
    !==========
   ( propce ,                                                     &
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
   ( propce ,                                                     &
     taup   , tlag   , piil   ,                                   &
     bx     , tsfext , vagaus ,                                   &
     gradpr ,                                                     &
     romp   , brgaus , terbru , fextla )

endif

! Free memory
if (lamvbr.eq.1) then
  deallocate(brgaus)
endif
deallocate(vagaus)
deallocate(fextla)

!===============================================================================

end subroutine
