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

subroutine d3pphy &
!================

 ( idbia0 , idbra0 ,                                              &
   nvar   , nscal  , nphas  ,                                     &
   nphmx  ,                                                       &
   ibrom  , izfppp ,                                              &
   ia     ,                                                       &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   dirmin , dirmax , fdeb   , ffin   , hrec   ,                   &
   w1     , w2     , w3     ,                                     &
   ra     )

!===============================================================================
! FONCTION :
! --------

! ROUTINE PHYSIQUE PARTICULIERE : FLAMME DE DIFFUSION
! Calcul de RHO mutualise pour chimie 3 points
!  adiabatique ou permeatique (transport de H)


! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! nphas            ! i  ! <-- ! number of phases                               !
! nphmx            ! e  ! <-- ! nphsmx                                         !
! ibrom            ! te ! <-- ! indicateur de remplissage de romb              !
!        !    !     !                                                !
! izfppp           ! te ! --> ! numero de zone de la face de bord              !
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
! dirmin           ! tr ! --- ! pdf : dirac en fmin                            !
! dirmax           ! tr ! --- ! pdf : dirac en fmax                            !
! fdeb             ! tr ! --- ! pdf : abscisse debut rectangle                 !
! ffin             ! tr ! --- ! pdf : abscisse fin rectangle                   !
! hrec             ! tr ! --- ! pdf : hauteur rectangle                        !
! w1-w3(ncelet)    ! tr ! --- ! tableau de travail                             !
! ra(*)            ! ra ! --- ! main real work array                           !
! rpp              ! tr ! --- ! macro tableau reel pp                          !
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
use radiat
use mesh

!===============================================================================

implicit none

! Arguments

integer          idbia0 , idbra0
integer          iinpdf
integer          nvar   , nscal  , nphas
integer          nphmx

integer          ibrom
integer          izfppp(nfabor)
integer          ia(*)

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision dirmin(ncelet),dirmax(ncelet)
double precision fdeb(ncelet),ffin(ncelet)
double precision hrec(ncelet)
double precision w1(ncelet),w2(ncelet),w3(ncelet)
double precision ra(*)

! Local variables

integer          idebia, idebra, ifinia
integer          if, ih, iel, icg
integer          ifac, mode, izone
integer          iphas, ipbrom, ipcrom, ipbycg, ipcycg
double precision coefg(ngazgm), fsir, hhloc, tstoea, tin
double precision temsmm

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

idebia = idbia0
idebra = idbra0
!     On reserve la memoire pour le tabeau INDPDF : passage ou non
!     par les PDF (on le garde pendant tout le sous-programme

iinpdf = idebia
ifinia = iinpdf + ncelet
call iasize('d3pphy',ifinia)
!==========

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

call pppdfr                                                       &
!==========
 ( ncelet,ncel, ia(iinpdf),                                       &
   rtp(1,isca(ifm)), rtp(1,isca(ifp2m)),                          &
   w1, w2,                                                        &
   dirmin, dirmax, fdeb, ffin, hrec )



!===============================================================================
! 4. INTEGRATION DE LA FONCTION DENSITE DE PROBABILITE
!    POUR DETERMINER LA TEMPERATURE
!                    LES FRACTIONS MASSIQUES
!                    LA MASSE VOLUMIQUE
!                    qsp RAYONNEMENT
!    Ces variables d'etat sont dans PROPCE
!===============================================================================

if ( ippmod(icod3p).eq.1 ) then

  call d3pint                                                     &
  !==========
 ( ncelet,ncel, ia(iinpdf),                                       &
   dirmin,dirmax,fdeb,ffin,hrec,                                  &
   rtp(1,isca(ifm)),rtp(1,isca(ihm)),rtp(1,ipr),               &
   propce,                                                        &
   w1 )

else

  call d3pint                                                     &
  !==========
 ( ncelet,ncel, ia(iinpdf),                                       &
   dirmin,dirmax,fdeb,ffin,hrec,                                  &
   rtp(1,isca(ifm)),w2,rtp(1,ipr),                             &
   propce,                                                        &
   w1 )

endif


!===============================================================================
! 4. CALCUL DE RHO ET DES FRACTIONS MASSIQUES DES ESPECES GLOBALES
!    SUR LES BORDS
!===============================================================================

! --> Masse volumique au bord

iphas = 1
ibrom = 1
ipbrom = ipprob(irom)
ipcrom = ipproc(irom)

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
    propfb(ifac,ipbrom) = propce(iel,ipcrom)
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
        propfb(ifac,ipbrom) = p0/(rr*temsmm)
      endif
    endif
  enddo
endif

! --> Fractions massiques des especes globales au bord
!     Uniquement si rayonnement

if ( iirayo.gt.0 ) then
  do icg = 1, ngazg
    ipbycg = ipprob(iym(icg))
    ipcycg = ipproc(iym(icg))
    do ifac = 1, nfabor
      iel = ifabor(ifac)
      propfb(ifac,ipbycg) = propce(iel,ipcycg)
    enddo
  enddo
endif


!===============================================================================
! FORMATS
!----


!----
! FIN
!----

return
end subroutine
