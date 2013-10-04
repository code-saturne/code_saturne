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
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[out]    mbrom         indicator of boundary density array filling
!> \param[in]     izfppp        boundary zone index for specific physic
!> \param[in]     rtp           calculated variables at cell centers
!>                               (at current time steps)
!> \param[in,out] propce        physical properties at cell centers
!_______________________________________________________________________________

subroutine d3pphy &
 ( mbrom  , izfppp , rtp    , propce )

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

!===============================================================================

implicit none

! Arguments

integer          mbrom
integer          izfppp(nfabor)

double precision rtp(ncelet,*)
double precision propce(ncelet,*)

! Local variables

integer          if, ih, iel, icg
integer          ifac, mode, izone
integer          ipcycg

double precision coefg(ngazgm), fsir, hhloc, tstoea, tin
double precision temsmm

integer, allocatable, dimension(:) :: indpdf

double precision, allocatable, dimension(:) :: dirmin, dirmax
double precision, allocatable, dimension(:) :: fdeb, ffin
double precision, allocatable, dimension(:) :: hrec, tpdf
double precision, allocatable, dimension(:) :: w1, w2
double precision, dimension(:), pointer :: brom,  crom
double precision, dimension(:), pointer :: bsval
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

!===============================================================================
! 2. DETERMINATION DES GRANDEURS THERMOCHIMIQUES
!===============================================================================

if ( ipass.le.2 ) then

! Rq : Il faut avoir vu usd3pc.F pour calculer correctement HH et FF

! ---> Calcul de TSTOEA

! ---- Initialisation
  do icg = 1, ngazgm
    coefg(icg) = zero
  enddo

  hstoea = fs(1)*hinfue + (1.d0-fs(1))*hinoxy
  coefg(1) = zero
  coefg(2) = zero
  coefg(3) = 1.d0
  mode = 1
  call cothht                                                     &
  !==========
  ( mode   , ngazg , ngazgm  , coefg  ,                           &
    npo    , npot   , th     , ehgazg ,                           &
    hstoea , tstoea )


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
  mode    = -1
  call cothht                                                     &
  !==========
  ( mode      , ngazg , ngazgm  , coefg  ,                        &
    npo       , npot   , th     , ehgazg ,                        &
    hh(nmaxh) , tin    )
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
      mode = 1
      call cothht                                                 &
      !==========
      ( mode   , ngazg , ngazgm  , coefg  ,                       &
        npo    , npot   , th     , ehgazg ,                       &
        hhloc , tfh(if,ih) )
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
      mode = 1
      call cothht                                                 &
      !==========
      ( mode   , ngazg , ngazgm  , coefg  ,                       &
        npo    , npot   , th     , ehgazg ,                       &
        hhloc , tfh(if,ih) )
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
   rtp(1,isca(ifm))      , rtp(1,isca(ifp2m)),             &
   w1    , w2    ,                                         &
   dirmin, dirmax, fdeb  , ffin, hrec )



!===============================================================================
! 4. INTEGRATION DE LA FONCTION DENSITE DE PROBABILITE
!    POUR DETERMINER LA TEMPERATURE
!                    LES FRACTIONS MASSIQUES
!                    LA MASSE VOLUMIQUE
!                    qsp RAYONNEMENT
!    Ces variables d'etat sont dans PROPCE
!===============================================================================

call d3pint &
!==========
 ( indpdf ,                                                       &
   dirmin , dirmax , fdeb   , ffin , hrec , tpdf ,                &
   rtp    , propce , w1 )

! Free memory
deallocate(indpdf)

!===============================================================================
! 4. CALCUL DE RHO ET DES FRACTIONS MASSIQUES DES ESPECES GLOBALES
!    SUR LES BORDS
!===============================================================================

! --> Masse volumique au bord

mbrom = 1
call field_get_val_s(ibrom, brom)
call field_get_val_s(icrom, crom)

! ---- Masse volumique au bord pour toutes les facettes
!      Les facettes d'entree seront recalculees apres
!      a partir des CL (si IPASS > 1). .


! ---- Au premier passage sans suite ou si on n'a pas relu la
!      masse volumique dans le fichier suite, on n'a pas recalcule la
!      masse volumique dans d3pint, pas la peine de la reprojeter aux
!      faces.

if ( ipass.gt.1.or.(isuite.eq.1.and.initro.eq.1)) then

  do ifac = 1, nfabor
    iel = ifabor(ifac)
    brom(ifac) = crom(iel)
  enddo

endif

! ---- Masse volumique au bord pour les facettes d'entree UNIQUEMENT
!      Le test sur IZONE sert pour les reprises de calcul
!      On suppose implicitement que les elements ci-dessus ont ete relus
!      dans le fichier suite (i.e. pas de suite en combustion d'un calcul
!      a froid) -> sera pris en compte eventuellement dans les versions
!      suivantes

if(ipass.gt.1 .or. isuite.eq.1 ) then
  do ifac = 1, nfabor
    izone = izfppp(ifac)
    if(izone.gt.0) then
      if ( ientfu(izone).eq.1 .or. ientox(izone).eq.1 ) then
        temsmm = tinfue/wmolg(1)
        if ( ientox(izone).eq.1 ) temsmm = tinoxy/wmolg(2)
        brom(ifac) = p0/(rr*temsmm)
      endif
    endif
  enddo
endif

! --> Fractions massiques des especes globales au bord
!     Uniquement si rayonnement

if (iirayo.gt.0) then
  do icg = 1, ngazg
    call field_get_val_s(ibym(1), bsval)
    ipcycg = ipproc(iym(icg))
    do ifac = 1, nfabor
      iel = ifabor(ifac)
      bsval(ifac) = propce(iel,ipcycg)
    enddo
  enddo
endif


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
