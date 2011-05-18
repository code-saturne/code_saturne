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

subroutine raypun &
!================

 ( idbia0 , idbra0 ,                                              &
   nvar   , nscal  ,                                              &
   itypfb ,                                                       &
   ia     ,                                                       &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   cofrua , cofrub ,                                              &
   flurds , flurdb ,                                              &
   dtr    , viscf  , viscb  ,                                     &
   dam    , xam    ,                                              &
   drtp   , smbrs  , rovsdt ,                                     &
   theta4 , thetaa , sa     ,                                     &
   qx     , qy     , qz     ,                                     &
   qincid , eps    , tparoi ,                                     &
   w1     , w2     , w3     , w4     , w5     ,                   &
   w6     , w7     , w8     , w9     , ckmel  ,                   &
   ra     )

!===============================================================================
! FONCTION :
! ----------

!   SOUS-PROGRAMME DU MODULE DE RAYONNEMENT :
!   -----------------------------------------

!   CALCUL DES FLUX ET DU TERME SOURCE RADIATIFS
!   AVEC L'APPROXIMATION P-1

!-------------------------------------------------------------------------------
!ARGU                             ARGUMENTS
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! itypfb           ! ia ! <-- ! boundary face types                            !
! ia(*)            ! ia ! --- ! main integer work array                        !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
! cofrua,cofrub    ! tr ! --- ! conditions aux limites aux                     !
!(nfabor)          !    !     !    faces de bord pour la luminance             !
! flurds,flurdb    ! tr ! --- ! pseudo flux de masse (faces internes           !
!(nfac)(nfabor)    !    !     !    et faces de bord )                          !
! dtr(ncelet)      ! tr ! --- ! dt*cdtvar                                      !
! viscf(nfac)      ! tr ! --- ! visc*surface/dist aux faces internes           !
! viscb(nfabor     ! tr ! --- ! visc*surface/dist aux faces de bord            !
! dam(ncelet       ! tr ! --- ! tableau de travail pour matrice                !
! xam(nfac,*)      ! tr ! --- ! tableau de travail pour matrice                !
! drtp(ncelet      ! tr ! --- ! tableau de travail pour increment              !
! smbrs(ncelet     ! tr ! --- ! tableau de travail pour sec mem                !
! rovsdt(ncelet    ! tr ! --- ! tableau de travail pour terme instat           !
! theta4(ncelet    ! tr ! --- ! pseudo temperature radiative                   !
! thetaa(ncelet    ! tr ! --- ! pseudo temp rar pdt precedent (nulle)          !
! sa (ncelet)      ! tr ! --> ! part d'absorption du terme source rad          !
! qxqyqz(ncelet    ! tr ! --> ! composante du vecteur densite de flux          !
!                  !    !     ! radiatif explicite                             !
! qincid(nfabor    ! tr ! --> ! densite de flux radiatif aux bords             !
! eps (nfabor)     ! tr ! <-- ! emissivite des facettes de bord                !
! tparoi(nfabor    ! tr ! <-- ! temperature de paroi en kelvin                 !
! w1...9(ncelet    ! tr ! --- ! tableau de travail                             !
! ckmel(ncelet)    ! tr ! <-- ! coeff d'absorption du melange                  !
!                  !    !     !   gaz-particules de charbon                    !
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
use pointe
use ppppar
use ppthch
use cpincl
use ppincl
use radiat
use parall
use period
use mesh

!===============================================================================

implicit none

! Arguments

integer          idbia0 , idbra0
integer          nvar   , nscal

integer          itypfb(nfabor)
integer          ia(*)

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)

double precision cofrua(nfabor), cofrub(nfabor)
double precision flurds(nfac), flurdb(nfabor)

double precision dtr(ncelet)
double precision viscf(nfac), viscb(nfabor)
double precision dam(ncelet), xam(nfac,2)
double precision drtp(ncelet), smbrs(ncelet)
double precision rovsdt(ncelet)

double precision theta4(ncelet), thetaa(ncelet)
double precision sa(ncelet)
double precision qx(ncelet), qy(ncelet), qz(ncelet)
double precision qincid(nfabor), tparoi(nfabor), eps(nfabor)

double precision w1(ncelet), w2(ncelet), w3(ncelet)
double precision w4(ncelet), w5(ncelet), w6(ncelet)
double precision w7(ncelet), w8(ncelet), w9(ncelet)
double precision ckmel(ncelet)
double precision ra(*)


! Local variables

character*80     cnom

integer          idebia, idebra
integer          ifac  , iel
integer          iconv1, idiff1, ndirc1, ireso1
integer          nitmap, nswrsp, nswrgp, iwarnp
integer          imgr1 , imligp, ircflp, ischcp, isstpp, iescap
integer          ncymap, nitmgp
integer          inum
integer          idtva0, ivar0
integer          inc, iccocg
double precision epsrgp, blencp, climgp, epsilp, extrap, epsrsp
double precision aa, aaa, aaaa, relaxp, thetap


!===============================================================================

!===============================================================================
! 0. GESTION MEMOIRE
!===============================================================================

idebia = idbia0
idebra = idbra0

!===============================================================================
! 1. PARAMETRAGE DU SOLVEUR ET INITIALISATION
!===============================================================================

!--> Gradient Conjugue

ireso1 = 0

!--> Parametrage de CODITS

! IVAR0= 0  LA VARIABLE N'EST ICI NI RIJ NI VITESSE
ivar0   = 0
nitmap  = 1000
!     IMRGRA  = 0
nswrsp  = 1
nswrgp  = 100
imligp  = -1
ircflp  = 1
ischcp  = 1
isstpp  = 0
iescap  = 0
imgr1   = 0
ncymap  = 100
nitmgp  = 10
iwarnp  = iimlum
blencp  = zero
epsilp  = 1.d-8
epsrsp  = 1.d-8
epsrgp  = 1.d-5
climgp  = 1.5d0
extrap  = zero
relaxp  = 1.d0

!--> Il y a des dirichlets

ndirc1 = 1

!--> Pas de convection pour le modele P1

iconv1 = 0

!--> Equation de diffusion

idiff1 = 1

!--> Remise a zero des tableaux avant resolution

do iel = 1,ncel
  drtp(iel)   = zero
  theta4(iel) = zero
  thetaa(iel) = zero
enddo

do ifac = 1,nfac
  flurds(ifac) = zero
enddo

do ifac = 1,nfabor
  flurdb(ifac) = zero
enddo

!===============================================================================
! 2. COEFFICIENT DE DIFFUSION AUX FACES
!===============================================================================

do iel = 1,ncel
  ckmel(iel) = 1.d0 / ckmel(iel)
enddo

call viscfa                                                       &
!==========
   ( idebia , idebra ,                                            &
     imvisf ,                                                     &
     ia     ,                                                     &
     ckmel  , viscf  , viscb  ,                                   &
     ra     )

!===============================================================================
! 3.  RESOLUTION
!===============================================================================

!     Parametre pour schemas en temps et stationnaire
thetap = 1.d0
idtva0 = 0

CNOM = ' '
WRITE(CNOM,'(A)') 'Rayon P1'
inum = 1
nomvar(inum) = cnom

call codits                                                       &
!==========
 ( idebia , idebra ,                                              &
   nvar   , nscal  ,                                              &
   idtva0 , ivar0  , iconv1 , idiff1 , ireso1 , ndirc1 , nitmap , &
   imrgra , nswrsp , nswrgp , imligp , ircflp ,                   &
   ischcp , isstpp , iescap ,                                     &
   imgr1  , ncymap , nitmgp , inum   , iwarnp ,                   &
   blencp , epsilp , epsrsp , epsrgp , climgp , extrap ,          &
   relaxp , thetap ,                                              &
   ia     ,                                                       &
   thetaa , thetaa , cofrua , cofrub , cofrua , cofrub ,          &
   flurds , flurdb ,                                              &
   viscf  , viscb  , viscf  , viscb  ,                            &
   rovsdt , smbrs  , theta4 ,                                     &
   dam    , xam    , drtp   ,                                     &
   w1     , w2     , w3     , w4     , w5     ,                   &
   w6     , w7     , w8     , w9     ,                            &
   ra     )

!===============================================================================
! 4. Vecteur densite de flux radiatif
!===============================================================================

!    En periodique et parallele, echange avant calcul du gradient
if (irangp.ge.0.or.iperio.eq.1) then
  call synsca(theta4)
  !==========
endif

!     Calcul de la densite du flux radiatif QX, QY, QZ

inc     = 1
iccocg  = 1
imligp  = -1
iwarnp  = iimlum
epsrgp  = 1.d-8
climgp  = 1.5d0
extrap  = 0.d0
nswrgp  = 100
ivar0   = 0

call grdcel                                                       &
!==========
   ( ivar0  , imrgra , inc    , iccocg , nswrgp , imligp,         &
     iwarnp , nfecra , epsrgp , climgp , extrap ,                 &
     ia     ,                                                     &
     theta4 , cofrua , cofrub ,                                   &
     w1     , w2     , w3     ,                                   &
     w4     , w5     , w6     ,                                   &
     ra     )

aa = - stephn * 4.d0 / 3.d0

do iel = 1,ncel
  aaa = aa * ckmel(iel)
  qx(iel) = w1(iel) * aaa
  qy(iel) = w2(iel) * aaa
  qz(iel) = w3(iel) * aaa
enddo

!===============================================================================
! 5. Terme Source Radiatif d'absorption et densite de flux incident
!===============================================================================

!     Calcul de la part d'absorption du terme Source Radiatif

aa = 4.d0 * stephn
do iel = 1,ncel
  sa(iel) = aa * theta4(iel)
enddo

!     Calcul du flux incident Qincid

do ifac = 1, nfabor

  iel = ifabor(ifac)

  if (itypfb(ifac).eq.iparoi .or.                           &
      itypfb(ifac).eq.iparug ) then

!--> Premiere version plus chere et legerement plus precise

    aaaa = tparoi(ifac)**4

    aaa  = 1.5d0 * distb(ifac) / ckmel(iel)                 &
           * ( 2.d0 /(2.d0-eps(ifac)) -1.d0 )
    aa   = ( aaa * aaaa + theta4(iel) ) / (1.d0 + aaa)

    qincid(ifac) = stephn * (2.d0 * aa - eps(ifac) * aaaa)        &
                       / (2.d0 - eps(ifac))

!--> Deuxieme version plus cheap mais moins precise

!         QINCID(IFAC) = STEPHN *
!    &    (2.D0 * THETA4(IFABOR(IFAC)) - EPS(IFAC) * TPAROI(IFAC)**4)
!    &  / (2.D0 - EPS(IFAC))

  else
    qincid(ifac) = stephn * theta4(iel)                           &
               + ( qx(iel) * surfbo(1,ifac) +                     &
                   qy(iel) * surfbo(2,ifac) +                     &
                   qz(iel) * surfbo(3,ifac) ) /                   &
                   (0.5d0 * surfbn(ifac) )
  endif

enddo

!===============================================================================

!--------
! FORMATS
!--------

!----
! FIN
!----

return

end subroutine
