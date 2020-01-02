!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2020 EDF S.A.
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

subroutine lwcphy &
!================

 ( mbrom  , izfppp )

!===============================================================================
! FONCTION :
! --------

! ROUTINE PHYSIQUE PARTICULIERE : FLAMME DE PREMELANGE MODELE LWC
! Calcul de RHO adiabatique ou permeatique (transport de H)


! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! mbrom            ! te ! <-- ! indicateur de remplissage de romb              !
! izfppp           ! te ! --> ! numero de zone de la face de bord              !
! (nfabor)         !    !     !  pour le module phys. part.                    !
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
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use mesh
use field

!===============================================================================

implicit none

! Arguments

integer          mbrom
integer          izfppp(nfabor)

! Local variables

integer          igg, iel
integer          izone , ifac
double precision coefg(ngazgm)
double precision nbmol , temsmm
double precision masmg
double precision, dimension(:), pointer :: brom,  crom
double precision, dimension(:), pointer :: bsval
double precision, dimension(:), pointer :: cvar_yfm, cvar_yfp2m
double precision, dimension(:), pointer :: cvar_fm, cvar_fp2m
double precision, dimension(:), pointer :: cvar_coyfp, cpro_ymgg

integer       ipass
data          ipass /0/
save          ipass

!===============================================================================

!===============================================================================
! 0. ON COMPTE LES PASSAGES
!===============================================================================

ipass = ipass + 1

!===============================================================================
! 1. INITIALISATIONS A CONSERVER
!===============================================================================

! --- Initialisation memoire


! ---> Initialisation

do igg = 1, ngazgm
  coefg(igg) = zero
enddo

! ---> Positions des variables, coefficients

call field_get_val_s(icrom, crom)
call field_get_val_s(ibrom, brom)

call field_get_val_s(ivarfl(isca(ifm)), cvar_fm)
call field_get_val_s(ivarfl(isca(ifp2m)), cvar_fp2m)
call field_get_val_s(ivarfl(isca(iyfm)), cvar_yfm)
call field_get_val_s(ivarfl(isca(iyfp2m)), cvar_yfp2m)
if (ippmod(icolwc).ge.2) call field_get_val_s(ivarfl(isca(icoyfp)), cvar_coyfp)

!===============================================================================
! 2. DETERMINATION DES GRANDEURS THERMOCHIMIQUES MOYENNES
!===============================================================================


if ( (ippmod(icolwc).eq.0) .or. (ippmod(icolwc).eq.1) ) then

  call pdflwc                                                     &
  !==========
   ( ncelet        , ncel          ,                              &
     cvar_fm       , cvar_fp2m     ,                              &
     cvar_yfm      , cvar_yfp2m )

endif

if ( (ippmod(icolwc).eq.2) .or. (ippmod(icolwc).eq.3) ) then

  call pdfpp3                                                     &
  !==========
   ( ncelet        , ncel          ,                              &
     cvar_fm       , cvar_fp2m     ,                              &
     cvar_yfm      , cvar_yfp2m    ,                              &
     cvar_coyfp )

endif

if ( (ippmod(icolwc).eq.4).or.(ippmod(icolwc).eq.5) ) then

  call pdfpp4                                                     &
  !==========
   ( ncelet        , ncel          ,                              &
     cvar_fm       , cvar_fp2m     ,                              &
     cvar_yfm      , cvar_yfp2m    ,                              &
     cvar_coyfp )

endif

!===============================================================================
! 3. CALCUL DE RHO ET DES FRACTIONS MASSIQUES DES ESPECES GLOBALES
!    SUR LES BORDS
!===============================================================================

! --> Masse volumique au bord

mbrom = 1

! ---- Masse volumique au bord pour toutes les facettes
!      Les facettes d'entree seront recalculees.

do ifac = 1, nfabor
iel = ifabor(ifac)
  brom(ifac) = crom(iel)
enddo

! ---- Masse volumique au bord pour les facettes d'entree UNIQUEMENT
!      Le test sur IZONE sert pour les reprises de calcul

if ( ipass.gt.1 .or. isuite.eq.1 ) then
  do ifac = 1, nfabor
    izone = izfppp(ifac)
    if(izone.gt.0) then
      if ( ientgb(izone).eq.1 .or. ientgf(izone).eq. 1) then
        coefg(1) = fment(izone)
        coefg(2) = 1.d0-fment(izone)
        coefg(3) = zero
        if ( ientgb(izone).eq.1 ) then
          coefg(1) = max(zero,(fment(izone)-fs(1))/(1.d0-fs(1)))
          coefg(3) = (fment(izone)-coefg(1))/fs(1)
          coefg(2) = 1.d0 - coefg(1) - coefg(3)
        endif
        nbmol = 0.d0
        do igg = 1, ngazg
          nbmol = nbmol + coefg(igg)/wmolg(igg)
        enddo
       masmg = 1.d0/nbmol
       temsmm = tkent(izone)/masmg
       brom(ifac) = p0/(cs_physical_constants_r*temsmm)
      endif
    endif
  enddo
endif

! --> Fractions massiques des especes globales au bord
do igg = 1, ngazg
  call field_get_val_s(ibym(igg), bsval)
  call field_get_val_s(iym(igg), cpro_ymgg)
  do ifac = 1, nfabor
    iel = ifabor(ifac)
    bsval(ifac) = cpro_ymgg(iel)
  enddo
enddo

!----
! FIN
!----

return
end subroutine
