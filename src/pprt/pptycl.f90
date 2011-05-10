!-------------------------------------------------------------------------------

!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2009 EDF S.A., France

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

subroutine pptycl &
!================

 ( idbia0 , idbra0 ,                                              &
   nvar   , nscal  ,                                              &
   icodcl , itrifb , itypfb , izfppp ,                            &
   ia     ,                                                       &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  , rcodcl ,                                     &
   w1     , w2     , w3     , w4     , w5     , w6     , coefu  , &
   ra     )

!===============================================================================
! FONCTION :
! --------

!    CONDITIONS AUX LIMITES AUTOMATIQUES

!           AIGUILLAGE SPECIFIQUE AUX PHYSIQUES PARTICULIERES


!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
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
! itrifb           ! ia ! <-- ! indirection for boundary faces ordering        !
! itypfb           ! ia ! <-- ! boundary face types                            !
! izfppp           ! te ! <-- ! numero de zone de la face de bord              !
! (nfabor)         !    !     !  pour le module phys. part.                    !
! ia(*)            ! ia ! --- ! main integer work array                        !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
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
!                  !    !     !        cp*(viscls+visct/sigmas)*gradt          !
! w1,2,3,4,5,6     ! ra ! --- ! work arrays                                    !
!  (ncelet)        !    !     !  (computation of pressure gradient)            !
! coefu            ! ra ! --- ! work array                                     !
!  (nfabor, 3)     !    !     !  (computation of pressure gradient)            !
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

! Arguments

use paramx
use numvar
use optcal
use cstphy
use cstnum
use entsor
use parall
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use atincl
use mesh

!===============================================================================

implicit none

! Arguments

integer          idbia0 , idbra0
integer          nvar   , nscal

integer          icodcl(nfabor,nvar)
integer          itrifb(nfabor), itypfb(nfabor)
integer          izfppp(nfabor)
integer          ia(*)

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision rcodcl(nfabor,nvar,3)
double precision w1(ncelet),w2(ncelet),w3(ncelet)
double precision w4(ncelet),w5(ncelet),w6(ncelet)
double precision coefu(nfabor,ndim)
double precision ra(*)

! Local variables

integer          idebia, idebra
integer          ifac, iok, ifvu, ii, izone, izonem

!===============================================================================
!===============================================================================
! 1.  INITIALISATIONS
!===============================================================================

idebia = idbia0
idebra = idbra0


!===============================================================================
! 2.  LISTE DES ZONES (pour n'importe quel modele)
!===============================================================================

! --> Faces appartiennent toutes a une zone frontiere

iok = 0

do ifac = 1, nfabor
  if(izfppp(ifac).le.0.or.izfppp(ifac).gt.nozppm) then
    iok = iok + 1
    write(nfecra,1000)ifac,nozppm,izfppp(ifac)
  endif
enddo

if(iok.gt.0) then
  call csexit (1)
  !==========
endif

! --> On construit une liste des numeros des zones frontieres.
!           (liste locale au processeur, en parallele)
nzfppp = 0
do ifac = 1, nfabor
  ifvu = 0
  do ii = 1, nzfppp
    if (ilzppp(ii).eq.izfppp(ifac)) then
      ifvu = 1
    endif
  enddo
  if(ifvu.eq.0) then
    nzfppp = nzfppp + 1
    if(nzfppp.le.nbzppm) then
      ilzppp(nzfppp) = izfppp(ifac)
    else
      write(nfecra,1001) nbzppm
      write(nfecra,1002)(ilzppp(ii),ii=1,nbzppm)
      call csexit (1)
      !==========
    endif
  endif
enddo

! ---> Plus grand numero de zone

izonem = 0
do ii = 1, nzfppp
  izone = ilzppp(ii)
  izonem = max(izonem,izone)
enddo
if(irangp.ge.0) then
  call parcmx(izonem)
  !==========
endif
nozapm = izonem

 1000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : PHYSIQUE PARTICULIERE                       ',/,&
'@    =========                                               ',/,&
'@    LES CONDITIONS AUX LIMITES SONT INCOMPLETES OU ERRONEES ',/,&
'@                                                            ',/,&
'@  Le numero de zone associee a la face ',I10   ,' doit etre ',/,&
'@    un entier strictement positif et inferieur ou egal a    ',/,&
'@    NOZPPM = ',I10                                           ,/,&
'@  Ce numero (IZFPPP(IFAC)) vaut ici ',I10                    ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier les conditions aux limites.                      ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 1001 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : PHYSIQUE PARTICULIERE                       ',/,&
'@    =========                                               ',/,&
'@    PROBLEME DANS LES CONDITIONS AUX LIMITES                ',/,&
'@                                                            ',/,&
'@  Le nombre maximal de zones frontieres qui peuvent etre    ',/,&
'@    definies par l''utilisateur est NBZPPM = ',I10           ,/,&
'@    Il a ete depasse.                                       ',/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier les conditions aux limites.                      ',/,&
'@                                                            ',/,&
'@  Les NBZPPM premieres zones frontieres                     ',/,&
'@    portent ici les numeros suivants :                      ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 1002 format(i10)



!===============================================================================
! 5.  REMPLISSAGE DU TABLEAU DES CONDITIONS LIMITES
!       ON BOUCLE SUR TOUTES LES FACES D'ENTREE
!                     =========================
!         ON DETERMINE LA FAMILLE ET SES PROPRIETES
!           ON IMPOSE LES CONDITIONS AUX LIMITES
!           POUR LES SCALAIRES
!    (selon le modele)
!===============================================================================


! ---> Chimie 3 points : USD3PC

if ( ippmod(icod3p).ge.0 ) then

  call d3ptcl                                                     &
  !==========
 ( idebia , idebra ,                                              &
   nvar   , nscal  ,                                              &
   icodcl , itrifb , itypfb , izfppp ,                            &
   ia     ,                                                       &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  , rcodcl ,                                     &
   w1     , w2     , w3     , w4     , w5     , w6     , coefu  , &
   ra     )

! ---> Combustion gaz USEBUC
!      Flamme de premelange modele EBU

elseif ( ippmod(icoebu).ge.0 ) then

  call ebutcl                                                     &
  !==========
 ( idebia , idebra ,                                              &
   nvar   , nscal  ,                                              &
   icodcl , itrifb , itypfb , izfppp ,                            &
   ia     ,                                                       &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  , rcodcl ,                                     &
   w1     , w2     , w3     , w4     , w5     , w6     , coefu  , &
   ra     )

! ---> Combustion gaz USLWCC
!      Flamme de premelange modele LWC

elseif ( ippmod(icolwc).ge.0 ) then

  call lwctcl                                                     &
  !==========
 ( idebia , idebra ,                                              &
   nvar   , nscal  ,                                              &
   icodcl , itrifb , itypfb , izfppp ,                            &
   ia     ,                                                       &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  , rcodcl ,                                     &
   w1     , w2     , w3     , w4     , w5     , w6     , coefu  , &
   ra     )

! ---> Combustion charbon pulverise USCPCL

elseif ( ippmod(icp3pl).ge.0 ) then

  call cpptcl                                                     &
  !==========
 ( idebia , idebra ,                                              &
   nvar   , nscal  ,                                              &
   icodcl , itrifb , itypfb , izfppp ,                            &
   ia     ,                                                       &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  , rcodcl ,                                     &
   w1     , w2     , w3     , w4     , w5     , w6     , coefu  , &
   ra     )

! ---> Combustion charbon pulverise couple Lagrangien USCPLC

elseif ( ippmod(icpl3c).ge.0 ) then

  call cpltcl                                                     &
  !==========
 ( idebia , idebra ,                                              &
   nvar   , nscal  ,                                              &
   icodcl , itrifb , itypfb , izfppp ,                            &
   ia     ,                                                       &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  , rcodcl ,                                     &
   w1     , w2     , w3     , w4     , w5     , w6     , coefu  , &
   ra     )

! ---> Combustion fuel USFUCL

elseif ( ippmod(icfuel).ge.0 ) then

  call fuptcl                                                     &
  !==========
 ( idebia , idebra ,                                              &
   nvar   , nscal  ,                                              &
   icodcl , itrifb , itypfb , izfppp ,                            &
   ia     ,                                                       &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  , rcodcl ,                                     &
   w1     , w2     , w3     , w4     , w5     , w6     , coefu  , &
   ra     )

! ---> Compressible USCFCL

elseif ( ippmod(icompf).ge.0 ) then

  call cfxtcl                                                     &
  !==========
 ( idebia , idebra ,                                              &
   nvar   , nscal  ,                                              &
   icodcl , itrifb , itypfb , izfppp ,                            &
   ia     ,                                                       &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  , rcodcl ,                                     &
   w1     , w2     , w3     , w4     , w5     , w6     , coefu  , &
   ra     )

! ---> Ecoulements atmospheriques

elseif ( ippmod(iatmos).ge.0 ) then

  call attycl                                                     &
  !==========
 ( idebia , idebra ,                                              &
   nvar   , nscal  ,                                              &
   nbmetd , nbmett , nbmetm ,                                     &
   icodcl , itrifb , itypfb , izfppp , iprofm ,                   &
   ia     ,                                                       &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  , rcodcl ,                                     &
   ra(itmmet)      , ra(iztmet)      , ra(izdmet)      ,          &
   ra(ixmet)       , ra(iymet)       , ra(ipmer)       ,          &
   ra(ittmet)      , ra(iqvmet)      , ra(iumet)       ,          &
   ra(ivmet)       , ra(iekmet)      , ra(iepmet)      ,          &
   ra(irmet)       , ra(itpmet)      , ra(iphmet)      ,          &
   w1     , w2     , w3     , w4     , w5     , w6     , coefu  , &
   ra     )

endif
!----
! FORMATS
!----


!----
! FIN
!----

return
end subroutine
