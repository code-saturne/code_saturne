!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2022 EDF S.A.
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

subroutine stdtcl &
 ( nbzfmx , nozfmx ,                                              &
   icalke , dh     , xintur , itypfb , iznfbr , rcodcl )

!===============================================================================
! FONCTION :
! --------

!    CONDITIONS AUX LIMITES AUTOMATIQUES  EN STANDARD

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nbzfmx           ! e  ! <-- ! nb max de zones de faces de bord               !
! nozfmx           ! e  ! <-- ! numero max de zones de faces de bord           !
! itypfb           ! ia ! <-- ! boundary face types                            !
! iznfbr           ! te ! <-- ! numero de zone de la face de bord              !
! (nfabor)         !    !     !                                                !
! rcodcl           ! tr ! --> ! valeur des conditions aux limites              !
!                  !    !     !  aux faces de bord                             !
!                  !    !     ! rcodcl(1) = valeur du dirichlet                !
!                  !    !     ! rcodcl(2) = valeur du coef. d'echange          !
!                  !    !     !  ext. (infinie si pas d'echange)               !
!                  !    !     ! rcodcl(3) = valeur de la densite de            !
!                  !    !     !  flux (negatif si gain) w/m2 ou                !
!                  !    !     !  hauteur de rugosite (m) si icodcl=6           !
!                  !    !     ! pour les vitesses (vistl+visct)*gradu          !
!                  !    !     ! pour la pression             dt*gradp          !
!                  !    !     ! pour les scalaires                             !
!                  !    !     !        cp*(viscls+visct/turb_schmidt)*gradt    !
!__________________!____!_____!________________________________________________!

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use optcal
use cstphy
use cstnum
use dimens, only: nvar
use ppincl, only: itempk
use entsor
use parall
use mesh
use field
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer          nozfmx
integer          nbzfmx

integer          icalke(nozfmx)
integer          itypfb(nfabor)
integer          iznfbr(nfabor)

double precision dh(nozfmx), xintur(nozfmx)
double precision rcodcl(nfabor,nvar,3)

! Local variables

integer          ifac, izone, ifvu, izonem
integer          nozapm, nzfppp
integer          icke, ii, iel, iok, itk
double precision qisqc, viscla, uref2, rhomoy, dhy, xiturb, brom_loc
double precision, dimension(:), pointer :: brom
double precision, dimension(:), pointer :: viscl
integer, allocatable, dimension(:) :: ilzfbr
double precision, allocatable, dimension(:) :: qcalc

integer          ipass
data             ipass /0/
save             ipass

!===============================================================================
! 1.  INITIALISATIONS
!===============================================================================

allocate(ilzfbr(nbzfmx))
allocate(qcalc(nozfmx))


call field_get_val_s(ibrom, brom)
call field_get_val_s(iviscl, viscl)

! --> On construit une liste des numeros des zones frontieres.
!           (liste locale au processeur, en parallele)
nzfppp = 0
do ifac = 1, nfabor
  ifvu = 0
  do ii = 1, nzfppp
    if (ilzfbr(ii).eq.iznfbr(ifac)) then
      ifvu = 1
    endif
  enddo
  if(ifvu.eq.0) then
    nzfppp = nzfppp + 1
    if(nzfppp.le.nbzfmx) then
      ilzfbr(nzfppp) = iznfbr(ifac)
    else
      write(nfecra,1001) nbzfmx
      write(nfecra,1002)(ilzfbr(ii),ii=1,nbzfmx)
      call csexit(1)
    endif
  endif
enddo

! ---> Plus grand numero de zone

izonem = 0
do ii = 1, nzfppp
  izone = ilzfbr(ii)
  izonem = max(izonem,izone)
enddo
if(irangp.ge.0) then
  call parcmx(izonem)
  !==========
endif
nozapm = izonem

 1001 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : PROBLEME DANS LES CONDITIONS AUX LIMITES    ',/,&
'@    =========                                               ',/,&
'@                ARRET DANS LE SOUS-PROGRAMME STDTCL         ',/,&
'@                                                            ',/,&
'@  The maximum number of boundary zones which can be defined ',/,&
'@    by the user is NBZPPM = ',I10                            ,/,&
'@    It has been exceeded.                                   ',/,&
'@                                                            ',/,&
'@  The calculation will not run.                             ',/,&
'@                                                            ',/,&
'@  Verify the boundary conditions.                           ',/,&
'@                                                            ',/,&
'@  The first NBZPPM boundary zones have the following number:',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 1002 format(i10)

!===============================================================================
! 4.  REMPLISSAGE DU TABLEAU DES CONDITIONS LIMITES
!       ON BOUCLE SUR TOUTES LES FACES D'ENTREE
!                     =========================
!         ON DETERMINE LA FAMILLE ET SES PROPRIETES
!           ON IMPOSE LES CONDITIONS AUX LIMITES
!           POUR LA TURBULENCE
!===============================================================================

do ifac = 1, nfabor

  izone = iznfbr(ifac)

  if (izone .gt. 0) then

    if (     itypfb(ifac).eq.ientre  &
        .or. itypfb(ifac).eq.i_convective_inlet) then

! ----  Traitement automatique de la turbulence

      if ( icalke(izone).ne.0 ) then

!       La turbulence est calculee par defaut si ICALKE different de 0
!          - soit a partir du diametre hydraulique, d'une vitesse
!            de reference adaptes a l'entree courante si ICALKE = 1
!          - soit a partir du diametre hydraulique, d'une vitesse
!            de reference et de l'intensite turvulente
!            adaptes a l'entree courante si ICALKE = 2

        uref2 = rcodcl(ifac,iu,1)**2                       &
              + rcodcl(ifac,iv,1)**2                       &
              + rcodcl(ifac,iw,1)**2
        uref2 = max(uref2,epzero)
        rhomoy = brom(ifac)
        iel    = ifabor(ifac)
        viscla = viscl(iel)
        icke   = icalke(izone)
        dhy    = dh(izone)
        xiturb = xintur(izone)
        if (icke.eq.1) then
          call turbulence_bc_inlet_hyd_diam(ifac,                            &
                                            uref2, dhy, rhomoy, viscla,      &
                                            rcodcl)
        else if (icke.eq.2) then
          call turbulence_bc_inlet_turb_intensity(ifac,                      &
                                                  uref2, xiturb, dhy,        &
                                                  rcodcl)
        endif

      endif

    endif

    if (itypfb(ifac).eq.iparoi) then
        ! condition automatique en paroi pour alpha
        if (iturb.eq.32) then
          rcodcl(ifac,ial,1)  = 0.d0
        endif
    endif

  endif

enddo

! Free memory
deallocate(ilzfbr)
deallocate(qcalc)

!----
! End
!----

return
end subroutine
