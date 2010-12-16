!-------------------------------------------------------------------------------

!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2010 EDF S.A., France

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

subroutine mmtycl &
!================

 ( idbia0 , idbra0 ,                                              &
   nvar   , nscal  , nphas  ,                                     &
   itypfb , icodcl ,                                              &
   ia     ,                                                       &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   rcodcl ,                                                       &
   w1     , w2     , w3     , w4     , w5     , w6     ,          &
   ra     )

!===============================================================================
! FONCTION :
! --------

! TRAITEMENT DES CODES DE CONDITIONS POUR UN MAILLAGE MOBILE
!   LORS D'UN COUPLAGE DE TYPE ROTOR/STATOR

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! nphas            ! i  ! <-- ! number of phases                               !
! itypfb           ! ia ! <-- ! boundary face types                            !
!  (nfabor, nphas) !    !     !                                                !
! icodcl           ! te ! <-- ! code de condition limites aux faces            !
!  (nfabor,nvar    !    !     !  de bord                                       !
!                  !    !     ! = 1   -> dirichlet                             !
!                  !    !     ! = 3   -> densite de flux                       !
! ia(*)            ! ia ! --- ! main integer work array                        !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! rcodcl           ! tr ! <-- ! valeur des conditions aux limites              !
!  (nfabor,nvar    !    !     !  aux faces de bord                             !
!                  !    !     ! rcodcl(1) = valeur du dirichlet                !
!                  !    !     ! rcodcl(2) = valeur du coef. d'echange          !
!                  !    !     !  ext. (infinie si pas d'echange)               !
!                  !    !     ! rcodcl(3) = valeur de la densite de            !
!                  !    !     !  flux (negatif si gain) w/m2                   !
!                  !    !     ! pour les vitesses (vistl+visct)*gradu          !
!                  !    !     ! pour la pression             dt*gradp          !
!                  !    !     ! pour les scalaires                             !
!                  !    !     !        cp*(viscls+visct/sigmas)*gradt          !
! depmob(nnod,3    ! tr ! <-- ! deplacement aux noeuds                         !
! xyzno1(3,nnod    ! tr ! <-- ! coordonnees noeuds maillage initial            !
! w1,2,3,4,5,6     ! ra ! --- ! work arrays                                    !
!  (ncelet)        !    !     !  (computation of pressure gradient)            !
! ra(*)            ! ra ! --- ! main real work array                           !
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
use dimens, only: ndimfb
use numvar
use optcal
use cstnum
use cstphy
use entsor
use pointe
use parall
use mesh

!===============================================================================

implicit none

! Arguments

integer          idbia0 , idbra0
integer          nvar   , nscal  , nphas

integer          itypfb(nfabor,nphas)
integer          icodcl(nfabor,nvar)
integer          ia(*)

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(ndimfb,*)
double precision rcodcl(nfabor,nvar,3)
double precision depmob(nnod,3), xyzno1(3,nnod)
double precision w1(ncelet),w2(ncelet),w3(ncelet)
double precision w4(ncelet),w5(ncelet),w6(ncelet)
double precision ra(*)

! Local variables

integer          idebia, idebra
integer          ifac, iel, iphas, iuiph, iviph, iwiph
integer          ii, inod, icpt
double precision ddepx, ddepy, ddepz
double precision srfbnf, rnx, rny, rnz
double precision rcodcx, rcodcy, rcodcz, rcodsn
double precision vitbox, vitboy, vitboz


!===============================================================================

!===============================================================================
! 1.  INITIALISATIONS
!===============================================================================

idebia = idbia0
idebra = idbra0

!===============================================================================
! 2.  VITESSE DE DEFILEMENT POUR LES PAROIS FLUIDES ET SYMETRIES
!===============================================================================

! Pour les symetries on rajoute toujours la vitesse de maillage, car on
!   ne conserve que la vitesse normale
! Pour les parois, on prend la vitesse de maillage si l'utilisateur n'a
!   pas specifie RCODCL, sinon on laisse RCODCL pour la vitesse tangente
!   et on prend la vitesse de maillage pour la composante normale.
! On se base uniquement sur ITYPFB, a l'utilisateur de gere les choses
!   s'il rentre en CL non standards.

do ifac = 1, nfabor

  iel = ifabor(ifac)

  ! --- En turbomachine on connaît la valeur exacte de la vitesse de maillage


  do iphas = 1, nphas

    iuiph = iu(iphas)
    iviph = iv(iphas)
    iwiph = iw(iphas)

    vitbox = omegay*cdgfbo(3,ifac) - omegaz*cdgfbo(2,ifac)
    vitboy = omegaz*cdgfbo(1,ifac) - omegax*cdgfbo(3,ifac)
    vitboz = omegax*cdgfbo(2,ifac) - omegay*cdgfbo(1,ifac)

    if (itypfb(ifac,iphas).eq.isymet) then
      rcodcl(ifac,iuiph,1) = vitbox
      rcodcl(ifac,iviph,1) = vitboy
      rcodcl(ifac,iwiph,1) = vitboz
    endif

    if (itypfb(ifac,iphas).eq.iparoi) then
      ! Si une des composantes de vitesse de glissement a ete
      !    modifiee par l'utilisateur, on ne fixe que la vitesse
      !    normale
      if (rcodcl(ifac,iuiph,1).gt.rinfin*0.5d0 .and.              &
          rcodcl(ifac,iviph,1).gt.rinfin*0.5d0 .and.              &
          rcodcl(ifac,iwiph,1).gt.rinfin*0.5d0) then
        rcodcl(ifac,iuiph,1) = vitbox
        rcodcl(ifac,iviph,1) = vitboy
        rcodcl(ifac,iwiph,1) = vitboz
      else
      ! On met a 0 les composantes de RCODCL non specifiees
        if (rcodcl(ifac,iuiph,1).gt.rinfin*0.5d0) rcodcl(ifac,iuiph,1) = 0.d0
        if (rcodcl(ifac,iviph,1).gt.rinfin*0.5d0) rcodcl(ifac,iviph,1) = 0.d0
        if (rcodcl(ifac,iwiph,1).gt.rinfin*0.5d0) rcodcl(ifac,iwiph,1) = 0.d0

        srfbnf = ra(isrfbn-1+ifac)
        rnx = surfbo(1,ifac)/srfbnf
        rny = surfbo(2,ifac)/srfbnf
        rnz = surfbo(3,ifac)/srfbnf
        rcodcx = rcodcl(ifac,iuiph,1)
        rcodcy = rcodcl(ifac,iviph,1)
        rcodcz = rcodcl(ifac,iwiph,1)
        rcodsn = (vitbox - rcodcx)*rnx                            &
               + (vitboy - rcodcy)*rny                            &
               + (vitboz - rcodcz)*rnz
        rcodcl(ifac,iuiph,1) = rcodcx + rcodsn*rnx
        rcodcl(ifac,iviph,1) = rcodcy + rcodsn*rny
        rcodcl(ifac,iwiph,1) = rcodcz + rcodsn*rnz
      endif

    endif
  enddo
enddo

!===============================================================================
! FORMATS
!===============================================================================

return
end subroutine
