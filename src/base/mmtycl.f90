!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2019 EDF S.A.
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

subroutine mmtycl &
!================

 ( itypfb , rcodcl )

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
! itypfb           ! ia ! <-- ! boundary face types                            !
! rcodcl           ! tr ! <-- ! valeur des conditions aux limites              !
!                  !    !     !  aux faces de bord                             !
!                  !    !     ! rcodcl(1) = valeur du dirichlet                !
!                  !    !     ! rcodcl(2) = valeur du coef. d'echange          !
!                  !    !     !  ext. (infinie si pas d'echange)               !
!                  !    !     ! rcodcl(3) = valeur de la densite de            !
!                  !    !     !  flux (negatif si gain) w/m2                   !
!                  !    !     ! pour les vitesses (vistl+visct)*gradu          !
!                  !    !     ! pour la pression             dt*gradp          !
!                  !    !     ! pour les scalaires                             !
!                  !    !     !        cp*(viscls+visct/turb_schmidt)*gradt    !
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
use cstnum
use dimens, only: nvar
use rotation
use turbomachinery
use entsor
use parall
use mesh
use field

!===============================================================================

implicit none

! Arguments

integer          itypfb(nfabor)

double precision rcodcl(nfabor,nvar,3)

! Local variables

integer          ifac, iel
integer          idefau

double precision srfbnf, rnx, rny, rnz
double precision rcodcx, rcodcy, rcodcz, rcodsn
double precision vr(3)
double precision visclc, visctc, distbf, hint

double precision, dimension(:), pointer :: viscl, visct

!===============================================================================

!===============================================================================
! Initialization
!===============================================================================

! --- Physical quantities
call field_get_val_s(iviscl, viscl)
call field_get_val_s(ivisct, visct)

!===============================================================================
! VITESSE DE DEFILEMENT POUR LES PAROIS FLUIDES ET SYMETRIES
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

  if (irotce(iel).ne.0) then

    ! --- En turbomachine on connait la valeur exacte de la vitesse de maillage

    call rotation_velocity(irotce(iel), cdgfbo(:,ifac), vr)

    if (itypfb(ifac).eq.isymet) then
      rcodcl(ifac,iu,1) = vr(1)
      rcodcl(ifac,iv,1) = vr(2)
      rcodcl(ifac,iw,1) = vr(3)
    endif

    if (itypfb(ifac).eq.iparoi .or.                        &
        itypfb(ifac).eq.iparug ) then
      ! Si une des composantes de vitesse de glissement a ete
      !    modifiee par l'utilisateur, on ne fixe que la vitesse
      !    normale
      if (rcodcl(ifac,iu,1).gt.rinfin*0.5d0 .and.              &
          rcodcl(ifac,iv,1).gt.rinfin*0.5d0 .and.              &
          rcodcl(ifac,iw,1).gt.rinfin*0.5d0) then
        rcodcl(ifac,iu,1) = vr(1)
        rcodcl(ifac,iv,1) = vr(2)
        rcodcl(ifac,iw,1) = vr(3)
      else
        ! On met a 0 les composantes de RCODCL non specifiees
        if (rcodcl(ifac,iu,1).gt.rinfin*0.5d0) rcodcl(ifac,iu,1) = 0.d0
        if (rcodcl(ifac,iv,1).gt.rinfin*0.5d0) rcodcl(ifac,iv,1) = 0.d0
        if (rcodcl(ifac,iw,1).gt.rinfin*0.5d0) rcodcl(ifac,iw,1) = 0.d0

        srfbnf = surfbn(ifac)
        rnx = surfbo(1,ifac)/srfbnf
        rny = surfbo(2,ifac)/srfbnf
        rnz = surfbo(3,ifac)/srfbnf
        rcodcx = rcodcl(ifac,iu,1)
        rcodcy = rcodcl(ifac,iv,1)
        rcodcz = rcodcl(ifac,iw,1)
        rcodsn = (vr(1) - rcodcx)*rnx                            &
               + (vr(2) - rcodcy)*rny                            &
               + (vr(3) - rcodcz)*rnz
        rcodcl(ifac,iu,1) = rcodcx + rcodsn*rnx
        rcodcl(ifac,iv,1) = rcodcy + rcodsn*rny
        rcodcl(ifac,iw,1) = rcodcz + rcodsn*rnz
      endif

    endif

  endif
enddo

! In case of transient turbomachinery computations, assign the default values
! of coefficients associated to turbulent or rough wall velocity BC, in order
! to update the wall velocity after the geometry update (between prediction and
! correction step).
! Note that the velocity BC update is made only if the user has not specified
! any specific Dirichlet condition for velocity.

if (iturbo.eq.2) then
  do ifac = 1, nfabor

    iel = ifabor(ifac)

    if (rcodcl(ifac,iu,1).gt.rinfin*0.5d0 .and.  &
        rcodcl(ifac,iv,1).gt.rinfin*0.5d0 .and.  &
        rcodcl(ifac,iw,1).gt.rinfin*0.5d0) then
      idefau = 1
    else
      idefau = 0
    endif

    if (idefau.eq.1 .and. irotce(iel).ne.0 .and. &
      (itypfb(ifac).eq.iparoi.or.itypfb(ifac).eq.iparug)) then

      ! --- Physical Properties
      visclc = viscl(iel)
      visctc = visct(iel)

      ! --- Geometrical quantities
      distbf = distb(ifac)

      if (itytur.eq.3) then
        hint =   visclc         /distbf
      else
        hint = ( visclc+visctc )/distbf
      endif

      ! Coefficients associated to laminar wall Dirichlet BC

      coftur(ifac) = 0.d0
      hfltur(ifac) = hint

    else

      ! Large coefficient in others cases (unused)

      coftur(ifac) = rinfin
      hfltur(ifac) = rinfin

    endif

  enddo
endif

!===============================================================================
! FORMATS
!===============================================================================

return
end subroutine
