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

subroutine raysol &
!================

 ( idbia0 , idbra0 ,                                              &
   nvar   , nscal  , nphas  , iphas  ,                            &
   itypfb ,                                                       &
   ia     ,                                                       &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   cofrua , cofrub ,                                              &
   flurds , flurdb ,                                              &
   dtr    , viscf  , viscb  ,                                     &
   dam    , xam    ,                                              &
   drtp   , smbrs  , rovsdt ,                                     &
   w1     , w2     , w3     , w4     , w5     ,                   &
   w6     , w7     , w8     , w9     , w10    ,                   &
   ru     , rua    ,                                              &
   sa     ,                                                       &
   qx     , qy     , qz     ,                                     &
   qincid , snplus ,                                              &
   ra     )

!===============================================================================
! FONCTION :
! ----------

!   SOUS-PROGRAMME DU MODULE DE RAYONNEMENT :
!   -----------------------------------------

!   CALCUL DES FLUX ET DU TERME SOURCE RADIATIFS

!  1/ DONNEES DES LUMINANCES ENTRANTES AUX LIMITES DU DOMAINE
!        (C.L : REFLEXION ET EMISSION ISOTROPE)

!                               ->  ->           ->
!  2/ CALCUL DE LA LUMINANCE L( X , S ) AU POINT X

!                                    D L
!     PAR RESOLUTION DE L'EQUATION : --- = -TK.L +TS
!                                    D S
!                        ->                o
!     OU ENCORE : DIV (L.S ) = -TK.L + TK.L

!                                  ->   /    ->  ->  ->
!  3/ CALCUL DES DENSITES DE FLUX  Q = /  L( X , S ).S DOMEGA
!                                     /4.PI

!                                       /    ->  ->
!         ET DE L'ABSORPTION       SA= /  L( X , S ).  DOMEGA
!                                     /4.PI

!     PAR INTEGRATION DES LUMINANCES SUR LES ANGLES SOLIDES.

!     N . B : CA SERT A CALCULER LE TAUX D'ECHAUFFEMENT
!     -----
!                                       /    ->  ->  ->  ->
!  4/ CALCUL DU FLUX INCIDENT QINCID = /  L( X , S ).S . N DOMEGA
!                                     /->->
!        ->                          / S.N >0
!        N NORMALE FLUIDE VERS PAROI

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
! iphas            ! i  ! --> ! phase number                                   !
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
! w1...9(ncelet    ! tr ! --- ! tableau de travail                             !
! ru  (ncelet)     ! tr ! --- ! luminance                                      !
! rua (ncelet)     ! tr ! --- ! luminance pdt precedent (nulle)                !
! sa (ncelet)      ! tr ! --> ! part d'absorption du terme source rad          !
! qxqyqz(ncelet    ! tr ! --> ! composante du vecteur densite de flux          !
!                  !    !     ! radiatif explicite                             !
! qincid(nfabor    ! tr ! --> ! densite de flux radiatif aux bords             !
! snplus(nfabor    ! tr ! --- ! integration du demi-espace egale a pi          !
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
use mesh

!===============================================================================

implicit none

! Arguments

integer          idbia0 , idbra0
integer          nvar   , nscal  , nphas  , iphas

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
double precision w1(ncelet), w2(ncelet), w3(ncelet)
double precision w4(ncelet), w5(ncelet), w6(ncelet)
double precision w7(ncelet), w8(ncelet), w9(ncelet)
double precision w10(ncelet)

double precision ru(ncelet), rua(ncelet)
double precision sa(ncelet)
double precision qx(ncelet), qy(ncelet), qz(ncelet)
double precision qincid(nfabor), snplus(nfabor)

double precision ra(*)


! Local variables

character*80     cnom

integer          idebia, idebra
integer          ifac  , iel
integer          iconv1, idiff1, ndirc1, ireso1
integer          nitmap, nswrsp, nswrgp, iwarnp
integer          imgr1 , imligp, ircflp, ischcp, isstpp, iescap
integer          ncymap, nitmgp
integer          idir  , ndirs , kdir  , ipp   , inum
integer          ii, jj, kk, idtva0, ivar0
double precision epsrgp, blencp, climgp, epsilp, extrap, epsrsp
double precision sx, sy, sz, domega
double precision sxt(ndirs8), syt(ndirs8), szt(ndirs8)
double precision aa
double precision relaxp, thetap

!===============================================================================

!===============================================================================
! 0. GESTION MEMOIRE
!===============================================================================

idebia = idbia0
idebra = idbra0

!===============================================================================
! 1. INITIALISATION
!===============================================================================

ivar0   = 0
nitmap  = 1000
!     IMRGRA  = 0
nswrsp  = 2
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

!--> Convection pure

iconv1 = 1

!--> JACOBI

ireso1 = 1

!===============================================================================
! 2. Choix du nombre de directions et calcul des coordonnees SX,SY,SZ
!===============================================================================


ndirs = ndirec/8

domega = pi /(2 *ndirs)

call raydir                                                       &
!==========
     (sxt, syt, szt, ndirs)

!===============================================================================
!                                              / -> ->
! 3. CORRECTION DES C.L. POUR RESPECTER : PI= /  S. N DOMEGA
!                                            /2PI
!===============================================================================

do ifac = 1,nfabor
  snplus(ifac) = zero
enddo

do ii = -1,1,2
  do jj = -1,1,2
    do kk = -1,1,2
      do idir = 1,ndirs

        sx = ii *sxt (idir)
        sy = jj *syt (idir)
        sz = kk *szt (idir)

        do ifac = 1,nfabor
          aa = sx * surfbo(1,ifac)                                &
             + sy * surfbo(2,ifac)                                &
             + sz * surfbo(3,ifac)
          aa = aa / surfbn(ifac)
          snplus(ifac) =snplus(ifac) +0.5d0 *(-aa+abs(aa)) *domega
        enddo

      enddo
    enddo
  enddo
enddo

do ifac = 1,nfabor
  cofrua(ifac) = cofrua(ifac) *(pi /snplus(ifac))
enddo

!===============================================================================
! 4. INITIALISATION POUR INTEGRATION DANS LES BOUCLES SUIVANTES
!===============================================================================

do ifac = 1, nfabor
  qincid(ifac) = zero
  snplus(ifac) = zero
enddo

do iel = 1, ncelet
  sa(iel) = zero
  qx(iel) = zero
  qy(iel) = zero
  qz(iel) = zero
enddo

!--> Stockage du SMBRS dans tableau tampon, il sont recharges
!    a chaque changement de direction

do iel = 1, ncel
  w10(iel) =  smbrs(iel)
enddo

!--> ROVSDT charge une seule fois
do iel = 1, ncel
  rovsdt(iel) = max(rovsdt(iel),zero)
enddo

!--> ON SAUVEGARDE LE NOM DE LA PREMIERE VARIABLE
!    Attention, le passage du nom de variable pour les impressions
!      est fait par common dans NOMVAR (pas tres pratique).
!      Des infos sont egalement remplies en common (residu, nbiter...)
!      Les directions n'etant pas des variables qu'on se donne
!      la possibilite d'imprimer dans le listing, on utilise
!      la position numero 1 qui est une poubelle.

ipp  = 1
NOMVAR(IPP) = 'RayonXXX'

!===============================================================================
! 5. RESOLUTION DE L'EQUATION DES TRANSFERTS RADIATIFS
!===============================================================================

!===============================================================================
! 5.1 DISCRETISATION ANGULAIRE
!===============================================================================

kdir = 0

do ii = -1,1,2
  do jj = -1,1,2
    do kk = -1,1,2
      do idir = 1,ndirs

        sx = ii * sxt(idir)
        sy = jj * syt(idir)
        sz = kk * szt(idir)

        kdir = kdir + 1

        CNOM = ' '
        WRITE(CNOM,'(A5,I3.3)')'Rayon',KDIR
        inum = ipp
        nomvar(inum) = cnom

!===============================================================================
! 5.2 DISCRETISATION SPATIALE
!===============================================================================

!===============================================================================
! 5.1.1 PREPARATION ET PARAMETRAGE DE LA RESOLUTION
!===============================================================================

!--> Terme source explicite

        do iel = 1, ncel
          smbrs(iel) = w10(iel)
        enddo

!--> Terme source implicite (ROVSDT vu plus haut)

!--> Pas de diffusion facette

        idiff1 = 0
        do ifac = 1,nfac
          viscf(ifac) = zero
        enddo
        do ifac = 1,nfabor
          viscb(ifac) = zero
        enddo

        do iel = 1,ncelet
          drtp(iel) = zero
          ru(iel)   = zero
          rua(iel)  = zero
        enddo

        do ifac = 1,nfac
          flurds(ifac) =                                          &
               + sx*surfac(1,ifac)                                &
               + sy*surfac(2,ifac)                                &
               + sz*surfac(3,ifac)
        enddo

        do ifac = 1,nfabor
          flurdb(ifac) =                                          &
               + sx*surfbo(1,ifac)                                &
               + sy*surfbo(2,ifac)                                &
               + sz*surfbo(3,ifac)
        enddo

!===============================================================================
! 5.1.2 RESOLUTION
!===============================================================================

! Dans le cas d'un theta-schema on met theta = 1
! Pas de relaxation en stationnaire non plus

        thetap = 1.0d0
        idtva0 = 0

        call codits                                               &
        !==========
 ( idebia , idebra ,                                              &
   nvar   , nscal  , nphas  ,                                     &
   idtva0 , ivar0  , iconv1 , idiff1 , ireso1 , ndirc1 ,  nitmap ,&
   imrgra , nswrsp , nswrgp , imligp , ircflp ,                   &
   ischcp , isstpp , iescap ,                                     &
   imgr1  , ncymap , nitmgp , inum   , iwarnp ,                   &
   blencp , epsilp , epsrsp , epsrgp , climgp , extrap ,          &
   relaxp , thetap ,                                              &
   ia     ,                                                       &
   rua    , ru     ,                                              &
   cofrua , cofrub , cofrua , cofrub , flurds , flurdb ,          &
   viscf  , viscb  , viscf  , viscb  ,                            &
   rovsdt , smbrs  , ru    ,                                      &
   dam    , xam    , drtp   ,                                     &
   w1     , w2     , w3     , w4     , w5     ,                   &
   w6     , w7     , w8     , w9     ,                            &
   ra     )

!===============================================================================
! 5.2 INTEGRATION DES FLUX ET TERME SOURCE
!===============================================================================

        do iel = 1,ncel
          aa = ru(iel) *domega
          sa(iel) = sa(iel) + aa
          qx(iel) = qx(iel) + aa * sx
          qy(iel) = qy(iel) + aa * sy
          qz(iel) = qz(iel) + aa * sz
        enddo

!===============================================================================
! 5.3 FLUX INCIDENT A LA PAROI
!===============================================================================

        do ifac = 1,nfabor

          aa = sx * surfbo(1,ifac)                                &
             + sy * surfbo(2,ifac)                                &
             + sz * surfbo(3,ifac)
          aa = aa / surfbn(ifac)

          aa = 0.5d0 *(aa+abs(aa)) *domega

          snplus(ifac) = snplus(ifac) + aa

          qincid(ifac) = qincid(ifac) + aa*ru(ifabor(ifac))

        enddo

      enddo
    enddo
  enddo
enddo


!--------
! FORMATS
!--------

!----
! FIN
!----

return

end subroutine
