!-------------------------------------------------------------------------------

! This file is part of code_saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2023 EDF S.A.
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

!===============================================================================
! Function:
! ---------

!> \file d3pphy.f90
!>
!> \brief Specific physic subroutine: diffusion flame.
!>
!> Calculation of mean density
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine d3pphy ()

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
use radiat
use mesh
use field
use cs_c_bindings

!===============================================================================

implicit none

! Local variables

integer          if, ih, iel, igg
integer          ifac

double precision coefg(ngazgm), fsir, hhloc, tstoea, tin

integer, allocatable, dimension(:) :: indpdf

double precision, allocatable, dimension(:) :: dirmin, dirmax
double precision, allocatable, dimension(:) :: fdeb, ffin
double precision, allocatable, dimension(:) :: hrec, tpdf
double precision, allocatable, dimension(:) :: w1, w2
double precision, dimension(:), pointer :: bsval
double precision, dimension(:), pointer :: cvar_fm, cvar_fp2m
double precision, dimension(:), pointer :: cpro_ymgg

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

! Allocate temporary arrays
allocate(dirmin(ncelet), dirmax(ncelet))
allocate(fdeb(ncelet), ffin(ncelet))
allocate(hrec(ncelet), tpdf(ncelet))

! Allocate temporary arrays
allocate(w1(ncelet), w2(ncelet))

! --- Initialisation memoire


!     On reserve la memoire pour le tabeau INDPDF : passage ou non
!     par les PDF (on le garde pendant tout le sous-programme

allocate(indpdf(ncelet))

call field_get_val_s(ivarfl(isca(ifm)), cvar_fm)
call field_get_val_s(ivarfl(isca(ifp2m)), cvar_fp2m)

!===============================================================================
! 2. DETERMINATION DES GRANDEURS THERMOCHIMIQUES
!===============================================================================

if ( ipass.le.2 ) then

! Rq : Il faut avoir vu usd3pc.F pour calculer correctement HH et FF

! ---> Calcul de TSTOEA

! ---- Initialisation
  do igg = 1, ngazgm
    coefg(igg) = zero
  enddo

  hstoea = fs(1)*hinfue + (1.d0-fs(1))*hinoxy
  coefg(1) = zero
  coefg(2) = zero
  coefg(3) = 1.d0
  tstoea = cs_gas_combustion_h_to_t(coefg, hstoea)

! ---> Construction d'une table Temperature en fonction de la richesse
!        et de l'enthalpie stoechiometrique
!        de dimension 9X9

  fsir = fs(1)

! ---- Calcul du tableau FF(IF)

  do if = 1, (nmaxf/2+1)
    ff(if) = fsir * dble(2*if-2)/dble(nmaxf-1)
  enddo
  do if = (nmaxf/2+2), nmaxf
    ff(if) = fsir + dble(2*if-nmaxf-1)                            &
                  / dble(nmaxf-1)*(1.d0-fsir)
  enddo

! ----  Remplissage du tableau HH(IH)

  coefg(1) = zero
  coefg(2) = zero
  coefg(3) = 1.d0
  tin = min(tinfue,tinoxy)
  hh(nmaxh) = cs_gas_combustion_t_to_h(coefg, tin)
  hh(1) = hstoea
  do ih = 2, (nmaxh-1)
    hh(ih) = hh(1) + (hh(nmaxh)-hh(1))*                           &
                    dble(ih-1)/dble(nmaxh-1)
  enddo

! ---- Remplissage du tableau TFH(IF,IH)

  do ih = 1, nmaxh
    do if = 1, (nmaxf/2+1)
! ----- Melange pauvre
      coefg(1) = zero
      coefg(2) = (fsir-ff(if))/fsir
      coefg(3) = ff(if)/fsir
      hhloc = hinoxy + dble(2*if-2)/dble(nmaxf-1)                 &
                     * (hh(ih)-hinoxy)
      tfh(if,ih) = cs_gas_combustion_h_to_t(coefg, hhloc)
    enddo
    do if = (nmaxf/2+2), nmaxf
! ----- Melange riche
      coefg(1) = (ff(if)-fsir)/(1.d0-fsir)
      coefg(2) = 0.d0
      coefg(3) = (1.d0-ff(if))/(1.d0-fsir)
      hhloc = ( dble(2*if)*(hinfue-hh(ih))                        &
                + dble(2*nmaxf)*hh(ih)                            &
                - hinfue*dble(nmaxf+1) )                          &
            / dble(nmaxf-1)
      tfh(if,ih) = cs_gas_combustion_h_to_t(coefg, hhloc)
    enddo

  enddo

endif


!===============================================================================
! 3. CALCUL DE PARAMETRES DE LA FONCTION DENSITE DE PROBABILITE
!    POUR LA FRACTION DE MELANGE
!===============================================================================

! --- Definition des bornes min et max de la pdf
!     dans les 2 tableaux de travail W1 et W2

do iel = 1, ncel
  w1(iel) = 0.d0
  w2(iel) = 1.d0
enddo

call pppdfr &
!==========
 ( ncelet, ncel  , indpdf,                                 &
   tpdf  ,                                                 &
   cvar_fm       , cvar_fp2m ,                             &
   w1    , w2    ,                                         &
   dirmin, dirmax, fdeb  , ffin, hrec )



!===============================================================================
! 4. INTEGRATION DE LA FONCTION DENSITE DE PROBABILITE
!    POUR DETERMINER LA TEMPERATURE
!                    LES FRACTIONS MASSIQUES
!                    LA MASSE VOLUMIQUE
!                    qsp RAYONNEMENT
!    Ces variables d'etat sont des champs de type CS_FIELD_PROPERTY
!===============================================================================

call d3pint &
!==========
 ( indpdf ,                                                       &
   dirmin , dirmax , fdeb   , ffin , hrec , tpdf ,                &
   w1 )

! Free memory
deallocate(indpdf)

!===============================================================================
! 4. CALCUL DE DES FRACTIONS MASSIQUES DES ESPECES GLOBALES SUR LES BORDS
!===============================================================================

! --> Fractions massiques des especes globales au bord
do igg = 1, ngazg
  call field_get_val_s(ibym(igg), bsval)
  call field_get_val_s(iym(igg),cpro_ymgg)
  do ifac = 1, nfabor
    iel = ifabor(ifac)
    bsval(ifac) = cpro_ymgg(iel)
  enddo
enddo

! Free memory
deallocate(dirmin, dirmax)
deallocate(fdeb, ffin)
deallocate(hrec, tpdf)
deallocate(w1, w2)

!===============================================================================
! FORMATS
!----


!----
! FIN
!----

return
end subroutine
