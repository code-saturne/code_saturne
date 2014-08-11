!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2014 EDF S.A.
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

subroutine raysol &
!================

 ( coefap , coefbp ,                                              &
   cofafp , cofbfp ,                                              &
   flurds , flurdb ,                                              &
   viscf  , viscb  ,                                              &
   smbrs  , rovsdt ,                                              &
   ru     , rua    ,                                              &
   sa     ,                                                       &
   qx     , qy     , qz     ,                                     &
   qincid , snplus )

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
! coefap,coefbp    ! tr ! --- ! conditions aux limites aux                     !
!  cofafp, cofbfp  !    !     !    faces de bord pour la luminance             !
! flurds,flurdb    ! tr ! --- ! pseudo flux de masse (faces internes           !
!(nfac)(nfabor)    !    !     !    et faces de bord )                          !
! viscf(nfac)      ! tr ! --- ! visc*surface/dist aux faces internes           !
! viscb(nfabor     ! tr ! --- ! visc*surface/dist aux faces de bord            !
! smbrs(ncelet     ! tr ! --- ! tableau de travail pour sec mem                !
! rovsdt(ncelet    ! tr ! --- ! tableau de travail pour terme instat           !
! ru  (ncelet)     ! tr ! --- ! luminance                                      !
! rua (ncelet)     ! tr ! --- ! luminance pdt precedent (nulle)                !
! sa (ncelet)      ! tr ! --> ! part d'absorption du terme source rad          !
! qxqyqz(ncelet    ! tr ! --> ! composante du vecteur densite de flux          !
!                  !    !     ! radiatif explicite                             !
! qincid(nfabor    ! tr ! --> ! densite de flux radiatif aux bords             !
! snplus(nfabor    ! tr ! --- ! integration du demi-espace egale a pi          !
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
use mesh

!===============================================================================

implicit none

! Arguments

double precision coefap(nfabor), coefbp(nfabor)
double precision cofafp(nfabor), cofbfp(nfabor)
double precision flurds(nfac), flurdb(nfabor)

double precision viscf(nfac), viscb(nfabor)
double precision smbrs(ncelet)
double precision rovsdt(ncelet)

double precision ru(ncelet), rua(ncelet)
double precision sa(ncelet)
double precision qx(ncelet), qy(ncelet), qz(ncelet)
double precision qincid(nfabor), snplus(nfabor)

! Local variables

character(len=80) :: cnom

integer          ifac  , iel
integer          iconv1, idiff1, ndirc1, ireso1
integer          nitmap, nswrsp, nswrgp, iwarnp
integer          imgr1 , imligp, ircflp, ischcp, isstpp, iescap
integer          ncymap, nitmgp
integer          idir  , kdir  , ipp   , inum
integer          ii, jj, kk, idtva0, ivar0
integer          imucpp, idftnp, iswdyp
integer          icvflb
integer          ivoid(1)

double precision epsrgp, blencp, climgp, epsilp, extrap, epsrsp
double precision sxt, syt, szt, domegat
double precision aa
double precision relaxp, thetap

double precision rvoid(1)

double precision, allocatable, dimension(:) :: rhs0
double precision, allocatable, dimension(:) :: dpvar

!===============================================================================

!===============================================================================
! 0. GESTION MEMOIRE
!===============================================================================

! Allocate a work array
allocate(rhs0(ncelet))
allocate(dpvar(ncelet))

!===============================================================================
! 1. INITIALISATION
!===============================================================================

ivar0   = 0
nitmap  = 1000
nswrsp  = 2
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

!--> Convection pure

iconv1 = 1

!--> JACOBI

ireso1 = 1

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

        sxt = ii *sx (idir)
        syt = jj *sy (idir)
        szt = kk *sz (idir)

        domegat = angsol(idir)

        do ifac = 1,nfabor
          aa = sxt * surfbo(1,ifac)                                &
             + syt * surfbo(2,ifac)                                &
             + szt * surfbo(3,ifac)
          aa = aa / surfbn(ifac)
          snplus(ifac) = snplus(ifac) + 0.5d0 * (-aa+abs(aa)) * domegat
        enddo

      enddo
    enddo
  enddo
enddo

do ifac = 1, nfabor
  coefap(ifac) = coefap(ifac) *(pi /snplus(ifac))
  cofafp(ifac) = cofafp(ifac) *(pi /snplus(ifac))
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
  rhs0(iel) =  smbrs(iel)
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
nomva0 = 'RayonXXX'

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
      do idir = 1, ndirs

        sxt = ii * sx(idir)
        syt = jj * sy(idir)
        szt = kk * sz(idir)
        domegat = angsol(idir)

        kdir = kdir + 1

        cnom = ' '
        write(cnom,'(a5,i3.3)')'Rayon',kdir
        inum = ipp
        nomva0 = cnom

!===============================================================================
! 5.2 DISCRETISATION SPATIALE
!===============================================================================

!===============================================================================
! 5.1.1 PREPARATION ET PARAMETRAGE DE LA RESOLUTION
!===============================================================================

!--> Terme source explicite

        do iel = 1, ncel
          smbrs(iel) = rhs0(iel)
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
          ru(iel)   = zero
          rua(iel)  = zero
        enddo

        do ifac = 1,nfac
          flurds(ifac) =                                          &
               + sxt * surfac(1,ifac)                             &
               + syt * surfac(2,ifac)                             &
               + szt * surfac(3,ifac)
        enddo

        do ifac = 1,nfabor
          flurdb(ifac) =                                          &
               + sxt * surfbo(1,ifac)                             &
               + syt * surfbo(2,ifac)                             &
               + szt * surfbo(3,ifac)
        enddo

!===============================================================================
! 5.1.2 RESOLUTION
!===============================================================================

! Dans le cas d'un theta-schema on met theta = 1
! Pas de relaxation en stationnaire non plus

        thetap = 1.d0
        idtva0 = 0

        ! all boundary convective flux with upwind
        icvflb = 0

        call codits &
        !==========
 ( idtva0 , ivar0  , iconv1 , idiff1 , ireso1 , ndirc1 ,  nitmap ,&
   imrgra , nswrsp , nswrgp , imligp , ircflp ,                   &
   ischcp , isstpp , iescap , imucpp , idftnp , iswdyp ,          &
   imgr1  , ncymap , nitmgp , inum   , iwarnp ,                   &
   blencp , epsilp , epsrsp , epsrgp , climgp , extrap ,          &
   relaxp , thetap ,                                              &
   rua    , ru     ,                                              &
   coefap , coefbp , cofafp , cofbfp , flurds , flurdb ,          &
   viscf  , viscb  , rvoid  , viscf  , viscb  , rvoid  ,          &
   rvoid  , rvoid  ,                                              &
   icvflb , ivoid  ,                                              &
   rovsdt , smbrs  , ru     , dpvar  ,                            &
   rvoid  , rvoid  )

!===============================================================================
! 5.2 INTEGRATION DES FLUX ET TERME SOURCE
!===============================================================================

        do iel = 1,ncel
          aa = ru(iel) * domegat
          sa(iel) = sa(iel) + aa
          qx(iel) = qx(iel) + aa * sxt
          qy(iel) = qy(iel) + aa * syt
          qz(iel) = qz(iel) + aa * szt
        enddo

!===============================================================================
! 5.3 FLUX INCIDENT A LA PAROI
!===============================================================================

        do ifac = 1,nfabor

          aa = sxt * surfbo(1,ifac)                                &
             + syt * surfbo(2,ifac)                                &
             + szt * surfbo(3,ifac)
          aa = aa / surfbn(ifac)

          aa = 0.5d0 *(aa+abs(aa)) *domegat

          snplus(ifac) = snplus(ifac) + aa

          qincid(ifac) = qincid(ifac) + aa*ru(ifabor(ifac))

        enddo

      enddo
    enddo
  enddo
enddo

! Free memory
deallocate(rhs0)
deallocate(dpvar)

!--------
! FORMATS
!--------

!----
! FIN
!----

return

end subroutine
