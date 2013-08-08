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

subroutine raypun &
!================

 ( nvar   , nscal  ,                                              &
   itypfb ,                                                       &
   dt     , rtp    , rtpa   , propce ,                            &
   coefap , coefbp ,                                              &
   cofafp , cofbfp ,                                              &
   flurds , flurdb ,                                              &
   viscf  , viscb  ,                                              &
   smbrs  , rovsdt ,                                              &
   theta4 , thetaa , sa     ,                                     &
   qx     , qy     , qz     ,                                     &
   qincid , eps    , tparoi ,                                     &
   ckmel  )

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
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! itypfb           ! ia ! <-- ! boundary face types                            !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! coefap,coefbp    ! tr ! --- ! conditions aux limites aux                     !
!  cofafp, cofbfp  !    !     !    faces de bord pour la luminance             !
! flurds,flurdb    ! tr ! --- ! pseudo flux de masse (faces internes           !
!(nfac)(nfabor)    !    !     !    et faces de bord )                          !
! viscf(nfac)      ! tr ! --- ! visc*surface/dist aux faces internes           !
! viscb(nfabor     ! tr ! --- ! visc*surface/dist aux faces de bord            !
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
! ckmel(ncelet)    ! tr ! <-- ! coeff d'absorption du melange                  !
!                  !    !     !   gaz-particules de charbon                    !
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

integer          nvar   , nscal

integer          itypfb(nfabor)

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision coefap(nfabor), coefbp(nfabor)
double precision cofafp(nfabor), cofbfp(nfabor)
double precision flurds(nfac), flurdb(nfabor)

double precision viscf(nfac), viscb(nfabor)
double precision smbrs(ncelet)
double precision rovsdt(ncelet)

double precision theta4(ncelet), thetaa(ncelet)
double precision sa(ncelet)
double precision qx(ncelet), qy(ncelet), qz(ncelet)
double precision qincid(nfabor), tparoi(nfabor), eps(nfabor)

double precision ckmel(ncelet)

! Local variables

character*80     cnom

integer          ifac  , iel
integer          iconv1, idiff1, ndirc1, ireso1
integer          nitmap, nswrsp, nswrgp, iwarnp
integer          imgr1 , imligp, ircflp, ischcp, isstpp, iescap
integer          ncymap, nitmgp
integer          inum
integer          idtva0, ivar0
integer          inc, iccocg
integer          imucpp, idftnp, iswdyp
double precision epsrgp, blencp, climgp, epsilp, extrap, epsrsp
double precision aa, aaa, aaaa, relaxp, thetap

double precision rvoid(1)

double precision, allocatable, dimension(:,:) :: grad
double precision, allocatable, dimension(:) :: dpvar

!===============================================================================

!===============================================================================
! 0. GESTION MEMOIRE
!===============================================================================

! Allocate temporary array
allocate(dpvar(ncelet))

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
imucpp  = 0
idftnp  = 1
iswdyp  = 0
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
   ( imvisf ,                                                     &
     ckmel  , viscf  , viscb  )

!===============================================================================
! 3.  RESOLUTION
!===============================================================================

!     Parametre pour schemas en temps et stationnaire
thetap = 1.d0
idtva0 = 0

cnom = ' '
write(cnom,'(A)') 'Rayon P1'
inum = 1
nomvar(inum) = cnom

call codits &
!==========
 ( idtva0 , ivar0  , iconv1 , idiff1 , ireso1 , ndirc1 , nitmap , &
   imrgra , nswrsp , nswrgp , imligp , ircflp ,                   &
   ischcp , isstpp , iescap , imucpp , idftnp , iswdyp ,          &
   imgr1  , ncymap , nitmgp , inum   , iwarnp ,                   &
   blencp , epsilp , epsrsp , epsrgp , climgp , extrap ,          &
   relaxp , thetap ,                                              &
   thetaa , thetaa , coefap , coefbp , cofafp , cofbfp ,          &
   flurds , flurdb ,                                              &
   viscf  , viscb  , rvoid  , viscf  , viscb  , rvoid  ,          &
   rvoid  , rvoid  ,                                              &
   rovsdt , smbrs  , theta4 , dpvar  ,                            &
   rvoid  , rvoid  )

!===============================================================================
! 4. Vecteur densite de flux radiatif
!===============================================================================

! Allocate a temporary array for gradient computation
allocate(grad(ncelet,3))

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

call grdcel &
!==========
   ( ivar0  , imrgra , inc    , iccocg , nswrgp , imligp,         &
     iwarnp , nfecra , epsrgp , climgp , extrap ,                 &
     theta4 , coefap , coefbp ,                                   &
     grad   )

aa = - stephn * 4.d0 / 3.d0

do iel = 1, ncel
  aaa = aa * ckmel(iel)
  qx(iel) = grad(iel,1) * aaa
  qy(iel) = grad(iel,2) * aaa
  qz(iel) = grad(iel,3) * aaa
enddo

! Free memory
deallocate(grad)

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

! Free memory
deallocate(dpvar)

!--------
! FORMATS
!--------

!----
! FIN
!----

return

end subroutine
