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

subroutine cppdfr &
!================

 ( ncelet , ncel   ,                                              &
   indpdf , sc     , sd     , sp2m   ,                            &
   dsc    , dsd    , sdeb   , sfin   ,                            &
   hrec   )

!===============================================================================
! FONCTION :
! ----------

! CALCUL DES PARAMETRES DE LA PDF
! PDF RECTANGLE - PICS DE DIRAC CENTREE PPl - AE

! LE RESULTAT EST :
! ---------------
!    CALCUL DES PARAMETRES ASSOCIES AUX FONCTIONS RECTANGLE - DIRAC

!         DSC contient   le Dirac en SC
!         DSD  - - - - - le Dirac en SD
!         SDEB - - - - - l'abcisse de debut du rectangle
!         SFIN - - - - - - - - - - -  fin - - - - - - -
!         HREC - - - - - la hauteur du rectangle
!         INDPDF indique le passage ou non par la pdf

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! ncelet           ! i  ! <-- ! number of extended (real + ghost) cells        !
! ncel             ! i  ! <-- ! number of cells                                !
! sp2m             ! tr ! <-- ! variance de s                                  !
! sc               ! tr ! <-- ! borne min de s                                 !
! sd               ! tr ! <-- ! borne max de s                                 !
! dsc              ! tr !  <- ! dirac en c                                     !
! dsd              ! tr !  <- ! dirac en d                                     !
! sdeb             ! tr !  <- ! abscisse debut rectangle                       !
! sfin             ! tr !  <- ! abscisse fin rectangle                         !
! hrec             ! tr !  <- ! hauteur rectangle                              !
! indpdf           ! tr !  <- ! indicateur passage ou non par pdf              !
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
use parall
use pointe
use ppppar
use ppthch
use coincl
use cpincl
use ppincl

!===============================================================================

implicit none

! Arguments

integer          ncelet , ncel
integer          indpdf(ncelet)

double precision sp2m(ncelet)
double precision sc(ncelet)   , sd(ncelet)
double precision dsc(ncelet)  , dsd(ncelet)
double precision sdeb(ncelet) , sfin(ncelet)
double precision hrec(ncelet)


! Local variables

integer          iel
integer          n1     , n2     , n3     , n4     , n5
#ifdef DEBUG
double precision mom0   , mom1   , mom2
#endif
double precision t1     , t2     , t3     , epsi


!===============================================================================

!===============================================================================
! 1.  CALCULS PRELIMINAIRES
!===============================================================================

epsi = 1.d-06
! Parametre relatif a la variance
t1 = 1.d-05
! Parametre relatif a la moyenne
t2 = 5.d-04
! Parametre supplementaire
t3   = sqrt(3.d0*t1)

do iel = 1, ncel
  if ( indpdf(iel).gt.0 ) then
    if ( (sp2m(iel).lt.t1) .or. (abs(sd(iel)+sc(iel)).lt.t2) )    &
      indpdf(iel) = 0
  endif
enddo


!===============================================================================
! 2.  CALCUL DES PARAMETRES DE LA FONCTION DENSITE DE PROBABILITE

!     DONNEES POUR LES FONCTIONS RECTANGLES
!       SC   contient la borne inferieure de S
!       SD   - - - -  la borne superieure de S
!       SP2M - - - -  la variance de S
!       SM   - - - -  la moyenne de S = 0

!     SC -- DEB ---- 0 ---- FIN ---------- SD

!     PARAMETRES ASSOCIES AUX FONCTIONS RECTANGLES
!       DSC    contient le pic de Dirac en C
!       DSD    - - - -  le pic de Dirac en D
!       SDEB   - - - -  l'abcisse de debut du rectangle
!       SFIN   - - - - - - - - - - - fin   - - - - - -
!       HREC   - - - -  la hauteur du rectangle
!       INDPDF - - - -  le type de PDF : 2,3,11,12,13
!                   (= 0 si pas de passage par les pdf)

!===============================================================================

do iel = 1, ncel

  if (indpdf(iel).gt.0) then

    if ( ( sd(iel)  .ge.(-sc(iel)) .and.                          &
         sp2m(iel).le.(sc(iel)**2)/3.d0  ) .or.                   &
         ( sd(iel)  .lt.(-sc(iel)) .and.                          &
         sp2m(iel).le.(sd(iel)**2)/3.d0  )       ) then

! --> Rectangle seul

      hrec(iel) = sqrt(3.d0*sp2m(iel))
!     ATTENTION HREC n'est pas la hauteur du rectangle
!                    mais un intermediaire de calcul
      dsc(iel) = zero
      dsd(iel) = zero
      sdeb(iel) = - hrec(iel)
      sfin(iel) = + hrec(iel)
    elseif ( sd(iel)  .ge.(-sc(iel)) .and.                        &
           sp2m(iel).le.(-sc(iel)/3.d0*(sc(iel)+2.d0*sd(iel)))    &
           ) then

! --> Rectangle et un Dirac en SC

      sdeb(iel) = sc(iel)
      sfin(iel) = (3.d0*sp2m(iel)+sdeb(iel)**2)/(-2.d0*sdeb(iel))
      dsc(iel)= (sfin(iel)+sdeb(iel))/(sfin(iel)-sdeb(iel))
      dsd(iel) = zero
    elseif ( sd(iel)  .lt.(-sc(iel)) .and.                        &
           sp2m(iel).lt. (-sd(iel)/3.d0*(2.d0*sc(iel)+sd(iel)))   &
           ) then

! --> Rectangle et un Dirac en SD

      sfin(iel) = sd(iel)
      sdeb(iel) = (3.d0*sp2m(iel)+sfin(iel)**2)/(-2.d0*sfin(iel))
      dsd(iel) = (-sfin(iel)-sdeb(iel))/(sfin(iel)-sdeb(iel))
      dsc(iel)= zero
    else

! --> Rectangle et deux Diracs en SC et SD

      sdeb(iel) = sc(iel)
      sfin(iel) = sd(iel)
      dsd(iel) = (3.d0*sp2m(iel)+sdeb(iel)**2                     &
           +2.d0*sdeb(iel)*sfin(iel))                             &
           / ((sfin(iel)-sdeb(iel))**2)
      dsc(iel)= dsd(iel)+(sfin(iel)+sdeb(iel))                    &
           /(sfin(iel)-sdeb(iel))

    endif

    if ( abs(sfin(iel)-sdeb(iel)).gt.epsi ) then
      hrec(iel) = (1.d0-dsc(iel)-dsd(iel))/(sfin(iel)-sdeb(iel))
    else
      sdeb(iel) = min(sd(iel),max(sc(iel),-t3))
      sfin(iel) = min(sd(iel),max(sc(iel),+t3))
      hrec(iel) = (1.d0-dsc(iel)-dsd(iel))/(sfin(iel)-sdeb(iel))
    endif

  endif

enddo


!===============================================================================
! 3. IMPRESSION
!===============================================================================

n1 = 0
n2 = 0
n3 = 0
n4 = 0
n5 = 0

do iel = 1, ncel
  if ( indpdf(iel).ne.0  ) n1 = n1+1
  if ( indpdf(iel).eq.3  ) n2 = n2+1
  if ( indpdf(iel).eq.12 ) n3 = n3+1
  if ( indpdf(iel).eq.13 ) n4 = n4+1
  if ( indpdf(iel).eq.11 ) n5 = n5+1

! ---- Test
#ifdef DEBUG
  if ( indpdf(iel).ne.0  ) then
    mom0 = hrec(iel)*(sfin(iel)-sdeb(iel))                        &
         + dsc(iel)+dsd(iel)
    if ( abs(mom0-1.d0).gt.epsi ) then
      write(NFECRA,*) 'PB CELLULE ',IEL
      write(NFECRA,*) 'MOM0 VC VA', MOM0, 1.D0
    endif
    mom1 = hrec(iel)*(sfin(iel)**2-sdeb(iel)**2)/2.d0             &
         + dsc(iel)*sc(iel)+dsd(iel)*sd(iel)
    if ( mom1.gt.epsi ) then
      write(NFECRA,*) 'PB CELLULE ',IEL
      write(NFECRA,*) 'MOM1 VC VA', MOM1, 0.D0
    endif
    mom2 = hrec(iel)*(sfin(iel)**3-sdeb(iel)**3)/3.d0             &
         + dsc(iel)*sc(iel)**2+dsd(iel)*sd(iel)**2
    if ( abs(mom2-sp2m(iel)).gt.epsi ) then
      write(NFECRA,*) 'PB CELLULE ',IEL
      write(NFECRA,*) 'MOM2 VC VA',  MOM2, SP2M(IEL)
    endif
  endif
#endif
! ---- Fin Test

enddo

if (irangp.ge.0) then
  call parcpt (n1)
  !==========
  call parcpt (n2)
  !==========
  call parcpt (n3)
  !==========
  call parcpt (n4)
  !==========
  call parcpt (n5)
  !==========
endif

write(nfecra,1000) ncel, n1, n2, n3, n4, n5


!----
! FORMATS
!----

 1000 format (/,                                                  &
'CONTROLE DES PARAMETRES DANS CPPDFR.F',/,                  &
'======================================',/,                 &
' Nb de points de calculs                                  = ',   &
   i6,/,                                                    &
' Nb de points turbulents (passage par les PDF)            = ',   &
   i6,/,                                                    &
! ..v.7..1....v    ....2....v....3....v....4....v....5....v....6....v....7.I
' Nb de points turbulents pour lesquels support PDF = I4M  = ',   &
   i6,/,                                                    &
' Nb de points turbulents pour lesquels C app. [I4,L3]     = ',   &
   i6,/,                                                    &
' - - - - - - - - - - - - pour lesquels C app. [I4,L5]     = ',   &
   i6,/,                                                    &
' - - - - - - - - - - - - pour lesquels C app. [L5,I3max]  = ',   &
   i6)


!----
! FIN
!----

return
end subroutine
