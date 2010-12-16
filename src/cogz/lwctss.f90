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

subroutine lwctss &
!================

 ( idbia0 , idbra0 ,                                              &
   nvar   , nscal  , nphas  , ncepdp , ncesmp ,                   &
   iscal  ,                                                       &
   icepdc , icetsm , itypsm ,                                     &
   izfppp , ia     ,                            &
   dt     , rtpa   , rtp    , propce , propfa , propfb ,          &
   coefa  , coefb  , ckupdc , smacel ,                            &
   smbrs  , rovsdt ,                                              &
   viscf  , viscb  , xam    ,                                     &
   w1     , w2     , w3     , w4     , w5     ,                   &
   w6     , w7     , w8     , w9     , w10    , w11    ,          &
   ra     )

!===============================================================================
! FONCTION :
! ----------

! ROUTINE PHYSIQUE PARTICULIERE : FLAMME PREMELANGE MODELE LWC
!   ON PRECISE LES TERMES SOURCES POUR UN SCALAIRE PP
!   SUR UN PAS DE TEMPS

! ATTENTION : LE TRAITEMENT DES TERMES SOURCES EST DIFFERENT
! ---------   DE CELUI DE USTSSC.F

! ON RESOUT ROVSDT*D(VAR) = SMBRS

! ROVSDT ET SMBRS CONTIENNENT DEJA D'EVENTUELS TERMES SOURCES
!  UTILISATEUR. IL FAUT DONC LES INCREMENTER ET PAS LES
!  ECRASER

! POUR DES QUESTIONS DE STABILITE, ON NE RAJOUTE DANS ROVSDT
!  QUE DES TERMES POSITIFS. IL N'Y A PAS DE CONTRAINTE POUR
!  SMBRS

! DANS LE CAS D'UN TERME SOURCE EN CEXP + CIMP*VAR ON DOIT
! ECRIRE :
!          SMBRS  = SMBRS  + CEXP + CIMP*VAR
!          ROVSDT = ROVSDT + MAX(-CIMP,ZERO)

! ON FOURNIT ICI ROVSDT ET SMBRS (ILS CONTIENNENT RHO*VOLUME)
!    SMBRS en kg variable/s :
!     ex : pour la vitesse            kg m/s2
!          pour les temperatures      kg degres/s
!          pour les enthalpies        Joules/s
!    ROVSDT en kg /s


!-------------------------------------------------------------------------------
!ARGU                             ARGUMENTS
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! nphas            ! i  ! <-- ! number of phases                               !
! ncepdp           ! i  ! <-- ! number of cells with head loss                 !
! ncesmp           ! i  ! <-- ! number of cells with mass source term          !
! iscal            ! i  ! <-- ! scalar number                                  !
! icepdc(ncelet    ! te ! <-- ! numero des ncepdp cellules avec pdc            !
! icetsm(ncesmp    ! te ! <-- ! numero des cellules a source de masse          !
! itypsm           ! te ! <-- ! type de source de masse pour les               !
! (ncesmp,nvar)    !    !     !  variables (cf. ustsma)                        !
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
! ckupdc           ! tr ! <-- ! tableau de travail pour pdc                    !
!  (ncepdp,6)      !    !     !                                                !
! smacel           ! tr ! <-- ! valeur des variables associee a la             !
! (ncesmp,*   )    !    !     !  source de masse                               !
!                  !    !     !  pour ivar=ipr, smacel=flux de masse           !
! smbrs(ncelet)    ! tr ! --> ! second membre explicite                        !
! rovsdt(ncelet    ! tr ! --> ! partie diagonale implicite                     !
! viscf(nfac)      ! tr ! --- ! tableau de travail    faces internes           !
! viscb(nfabor     ! tr ! --- ! tableau de travail    faces de bord            !
! xam(nfac,2)      ! tr ! --- ! tableau de travail    faces de bord            !
! w1..11(ncelet    ! tr ! --- ! tableau de travail    cellules                 !
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
use numvar
use entsor
use optcal
use cstphy
use cstnum
use parall
use period
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use mesh

!===============================================================================

implicit none

! Arguments

integer          idbia0 , idbra0
integer          nvar   , nscal  , nphas
integer          ncepdp , ncesmp
integer          iscal

integer          icepdc(ncepdp)
integer          icetsm(ncesmp), itypsm(ncesmp,nvar)
integer          izfppp(nfabor)
integer          ia(*)

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision ckupdc(ncepdp,6), smacel(ncesmp,nvar)
double precision smbrs(ncelet), rovsdt(ncelet)
double precision viscf(nfac), viscb(nfabor)
double precision xam(nfac,2)
double precision w1(ncelet), w2(ncelet), w3(ncelet)
double precision w4(ncelet), w5(ncelet), w6(ncelet)
double precision w7(ncelet), w8(ncelet), w9(ncelet)
double precision w10(ncelet), w11(ncelet)
double precision ra(*)

! Local variables

integer          idebia, idebra
integer          ivar, iel, iphas, idirac, ivar0
integer          iphydp
integer          inc , iccocg
integer          ipcvst
integer          ipcrom, ii

integer          iptscl(ndracm), ipfmal(ndracm)
integer          ipfmel(ndracm), iprhol(ndracm)
double precision sum, epsi
double precision tsgrad, tschim, tsdiss

!===============================================================================

!===============================================================================
! 1. INITIALISATION
!===============================================================================

idebia = idbia0
idebra = idbra0
epsi   = 1.0d-10

! --- Numero du scalaire a traiter : ISCAL

! --- Numero de la variable associee au scalaire a traiter ISCAL
ivar = isca(iscal)

! --- Numero de phase associee au scalaire ISCAL
iphas = iphsca(iscal)

! ---
ipcrom = ipproc(irom(iphas))
ipcvst = ipproc(ivisct(iphas))

! --- Numero des grandeurs physiques (voir usclim)
do idirac = 1, ndirac
  iptscl(idirac) = ipproc(itscl(idirac))
  ipfmal(idirac) = ipproc(ifmal(idirac))
  ipfmel(idirac) = ipproc(ifmel(idirac))
  iprhol(idirac) = ipproc(irhol(idirac))
enddo

!===============================================================================
! 2. PRISE EN COMPTE DES TERMES SOURCES
!===============================================================================

if ( ivar.eq.isca(iyfm) ) then

! ---> Terme source pour la fraction massique moyenne de fuel

  do iel = 1, ncel
      sum = zero
      do idirac = 1, ndirac
        sum  = sum + propce(iel,iprhol(idirac))                   &
           *propce(iel,iptscl(idirac))*volume(iel)
      enddo

! terme implicite

      if (rtpa(iel,ivar).gt.epsi) then
        rovsdt(iel) = rovsdt(iel) + max(-sum/rtpa(iel,ivar),zero)
      endif

! terme explicite

       smbrs(iel) =  smbrs(iel) + sum

  enddo

endif

! ---> Terme source pour la variance de la fraction massique moyenne de fuel

if (ivar.eq.isca(iyfp2m)) then

  do iel = 1, ncel
    sum = zero
    do idirac = 1, ndirac
      sum  = sum + (propce(iel,iptscl(idirac))*volume(iel)        &
        *(propce(iel,ipfmal(idirac)) - rtpa(iel,isca(iyfm)))      &
             *propce(iel,iprhol(idirac)))
    enddo
    smbrs(iel) = smbrs(iel) + sum
  enddo

endif

! ---> Terme source pour la covariance

if ( ivar.eq.isca(icoyfp)) then

! --- Calcul du gradient de F
!     =======================

  ii = isca(ifm)
  do iel = 1, ncel
    w10(iel) = rtpa(iel,ii)
  enddo

  ! En periodique et parallele, echange avant calcul du gradient
  if (irangp.ge.0.or.iperio.eq.1) then
    call synsca(w10)
    !==========
  endif

!  IVAR0 = 0 (indique pour la periodicite de rotation que la variable
!     n'est pas la vitesse ni Rij)
  ivar0 = 0
  iphydp = 0
  inc = 1
  iccocg = 1

  call grdcel                                                     &
  !==========
 ( idebia , idebra ,                                              &
   nphas  ,                                                       &
   ivar0  , imrgra , inc    , iccocg , nswrgr(ii) , imligr(ii) ,  &
   iphydp , iwarni(ii) , nfecra ,                                 &
   epsrgr(ii) , climgr(ii) , extrag(ii) ,                         &
   ia     ,                                                       &
   w10    , w10    , w10    ,                                     &
   w10    , coefa(1,iclrtp(ii,icoef))  ,                          &
            coefb(1,iclrtp(ii,icoef))  ,                          &
   w1              , w2              , w3     ,                   &
!        d./dx1          , d./dx2          , d./dx3 ,
   w4     , w5     , w6     ,                                     &
   ra     )

! --- Calcul du gradient de Yfuel
!     ===========================

  ii = isca(iyfm)
  do iel = 1, ncel
    w11(iel) = rtpa(iel,ii)
  enddo

  ! En periodique et parallele, echange avant calcul du gradient
  if (irangp.ge.0.or.iperio.eq.1) then
    call synsca(w11)
    !==========
  endif

!  IVAR0 = 0 (indique pour la periodicite de rotation que la variable
!     n'est pas la vitesse ni Rij)
  ivar0 = 0
  iphydp = 0
  inc = 1
  iccocg = 1

  call grdcel                                                     &
  !==========
 ( idebia , idebra ,                                              &
   nphas  ,                                                       &
   ivar0  , imrgra , inc    , iccocg , nswrgr(ii) , imligr(ii) ,  &
   iphydp , iwarni(ii) , nfecra ,                                 &
   epsrgr(ii) , climgr(ii) , extrag(ii) ,                         &
   ia     ,                                                       &
   w11    , w11    , w11    ,                                     &
   w11    , coefa(1,iclrtp(ii,icoef))  ,                          &
            coefb(1,iclrtp(ii,icoef))  ,                          &
   w7              , w8              , w9     ,                   &
!        d./dx1          , d./dx2          , d./dx3 ,
   w4     , w5     , w6     ,                                     &
   ra     )


! --- Calcul du terme source
!     ======================


! ---> Calcul de K et Epsilon en fonction du modele de turbulence


! ---- TURBULENCE

  if (itytur(iphas).eq.2) then

    do iel = 1, ncel
      w10(iel) = rtpa(iel,ik(iphas))
      w11(iel) = rtpa(iel,iep(iphas))
    enddo

  elseif (itytur(iphas).eq.3) then

    do iel = 1, ncel
      w10(iel) = ( rtpa(iel,ir11(iphas))                          &
                  +rtpa(iel,ir22(iphas))                          &
                  +rtpa(iel,ir33(iphas)) ) / 2.d0
      w11(iel) = rtpa(iel,iep(iphas))
    enddo

  elseif (iturb(iphas).eq.50) then

    do iel = 1, ncel
      w10(iel) = rtpa(iel,ik(iphas))
      w11(iel) = rtpa(iel,iep(iphas))
    enddo

  elseif (iturb(iphas).eq.60) then

    do iel = 1, ncel
      w10(iel) = rtpa(iel,ik(iphas))
      w11(iel) = cmu*rtpa(iel,ik(iphas))*rtpa(iel,iomg(iphas))
    enddo

  endif

  do iel=1,ncel

!  A confirmer :
!   Le terme de dissipation devrait etre implicite
!   Dans le terme de dissipation, il manque une constante Cf
!   Peut-elle etre consideree egale a 1 ?
!   Verifier le signe du terme de production
!-
! terme implicite


    w11(iel) = w11(iel)/(w10(iel)*rvarfl(iscal))                  &
         *volume(iel)*propce(iel,ipcrom)
    rovsdt(iel) = rovsdt(iel) + max(w11(iel),zero)

! terme de gradient

    tsgrad =  (2.0d0                                              &
         * propce(iel,ipcvst)/(sigmas(iscal))                     &
         *(w1(iel)*w7(iel)+w2(iel)*w8(iel)+w3(iel)*w9(iel)))      &
         *volume(iel)


! terme de dissipation

    tsdiss = -w11(iel) * rtpa(iel,ivar)

! terme de chimique

    tschim = zero
    do idirac = 1, ndirac
      tschim =   tschim                                           &
           + (propce(iel,iptscl(idirac))                          &
           *(propce(iel,ipfmel(idirac))-rtpa(iel,isca(ifm)))      &
           *volume(iel))*propce(iel,iprhol(idirac))
    enddo

! --> Somme des termes

    smbrs(iel) = smbrs(iel) + tschim + tsgrad + tsdiss

   enddo

 endif

!----
! FIN
!----

return

end subroutine
