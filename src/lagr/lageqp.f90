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

subroutine lageqp &
!================

 ( ul     , vl     , wl     , alphal , phi    )

!===============================================================================
! FONCTION :
! ----------

!   SOUS-PROGRAMME DU MODULE LAGRANGIEN :
!   -------------------------------------

!          RESOLUTION D'UNE EQUATION DE POISSON

!            div[ALPHA grad(PHI)] = div(ALPHA <Up>)

!-------------------------------------------------------------------------------
!ARGU                             ARGUMENTS
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! ul,vl,wl(ncelet) ! tr ! <-- ! vitesse lagrangien                             !
! alphal(ncelet)   ! tr ! <-- ! taux de presence                               !
! phi(ncelet)      ! tr ! --> ! terme de correction                            !
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
use lagpar
use lagran
use mesh

!===============================================================================

implicit none

! Arguments

double precision ul(ncelet), vl(ncelet), wl(ncelet)
double precision phi(ncelet), alphal(ncelet)

! Local variables

character(len=80) :: chaine
integer          idtva0, ivar
integer          ifac, iel
integer          nswrgp, imligp, iwarnp , iescap
integer          iconvp, idiffp, ndircp, nitmap
integer          nswrsp, ircflp, ischcp, isstpp
integer          imucpp, idftnp, iswdyp
integer          icvflb
integer          ivoid(1)
double precision epsrgp, climgp, extrap, blencp, epsilp, epsrsp
double precision relaxp, thetap
double precision qimp  , hint, pimp

double precision rvoid(1)

double precision, allocatable, dimension(:) :: viscf, viscb
double precision, allocatable, dimension(:) :: smbrs, rovsdt
double precision, allocatable, dimension(:) :: fmala, fmalb
double precision, allocatable, dimension(:) :: phia
double precision, allocatable, dimension(:) :: w1, w2, w3
double precision, allocatable, dimension(:) :: coefax, coefay, coefaz
double precision, allocatable, dimension(:) :: coefbx, coefby, coefbz
double precision, allocatable, dimension(:) :: coefap, coefbp
double precision, allocatable, dimension(:) :: cofafp, cofbfp
double precision, allocatable, dimension(:) :: dpvar

!===============================================================================

!===============================================================================
! 1. INITIALISATION
!===============================================================================

! Allocate temporary arrays
allocate(viscf(nfac), viscb(nfabor))
allocate(smbrs(ncelet), rovsdt(ncelet))
allocate(fmala(nfac), fmalb(nfabor))
allocate(phia(ncelet))
allocate(dpvar(ncelet))

! Allocate work arrays
allocate(w1(ncelet), w2(ncelet), w3(ncelet))

CHAINE = 'Correction pression'
write(nfecra,1000) chaine(1:19)

!===============================================================================
! 2. TERMES SOURCES
!===============================================================================

! --> Initialisation

do iel = 1, ncel
  smbrs(iel) = 0.d0
enddo
do iel = 1, ncel
  rovsdt(iel) = 0.d0
enddo
do iel = 1, ncel
  phi(iel)  = 0.d0
  phia(iel) = 0.d0
enddo

!     "VITESSE" DE DIFFUSION FACE

  call viscfa                                                     &
  !==========
 ( imvisf ,                                                       &
   alphal ,                                                       &
   viscf  , viscb  )

! CALCUL  de div(Alpha Up) avant correction

do iel = 1, ncel
  w1(iel) = -ul(iel)*alphal(iel)
  w2(iel) = -vl(iel)*alphal(iel)
  w3(iel) = -wl(iel)*alphal(iel)
enddo

! --> Calcul du gradient de W1
!     ========================

! Allocate temporary arrays
allocate(coefax(nfabor), coefay(nfabor), coefaz(nfabor))
allocate(coefbx(nfabor), coefby(nfabor), coefbz(nfabor))

do ifac = 1, nfabor
  iel = ifabor(ifac)

  coefax(ifac) = w1(iel)
  coefbx(ifac) = zero

  coefay(ifac) = w2(iel)
  coefby(ifac) = zero

  coefaz(ifac) = w3(iel)
  coefbz(ifac) = zero

enddo

call diverv                                                       &
!==========
 ( smbrs  , w1     , w2     , w3     ,                            &
   coefax , coefay , coefaz ,                                     &
   coefbx , coefby , coefbz )

! Free memory
deallocate(coefax, coefay, coefaz)
deallocate(coefbx, coefby, coefbz)

! --> Conditions aux limites sur PHI
!     ==============================

! Allocate temporary arrays
allocate(coefap(nfabor), coefbp(nfabor))
allocate(cofafp(nfabor), cofbfp(nfabor))

do ifac = 1, nfabor
  iel = ifabor(ifac)

  hint = alphal(iel)/distb(ifac)

  if (itypfb(ifac).eq.ientre.or.itypfb(ifac).eq.iparoi.or.       &
      itypfb(ifac).eq.iparug.or.itypfb(ifac).eq.isymet) then

    ! Neumann Boundary Conditions
    !----------------------------

    qimp = 0.d0

    call set_neumann_scalar &
         !==================
       ( coefap(ifac), cofafp(ifac),             &
         coefbp(ifac), cofbfp(ifac),             &
         qimp        , hint )

    coefap(ifac) = zero
    coefbp(ifac) = 1.d0

  else if (itypfb(ifac).eq.isolib) then

    ! Dirichlet Boundary Condition
    !-----------------------------

    pimp = phia(iel)

    call set_dirichlet_scalar &
         !====================
       ( coefap(ifac), cofafp(ifac),             &
         coefbp(ifac), cofbfp(ifac),             &
         pimp        , hint        , rinfin )

  else
    write(nfecra,1100) itypfb(ifac)
    call csexit (1)
  endif

enddo

!===============================================================================
! 3. RESOLUTION
!===============================================================================

! Pas de stationnaire
idtva0 = 0
! Pas de terme de convection
iconvp = 0
! Diffusion
idiffp = 1
! Valeur par defaut
ndircp = 1
nitmap = 1000
nswrsp = 2
nswrgp = 10000
imligp = 1
ircflp = 1
ischcp = 1
isstpp = 0
imucpp = 0
idftnp = 1
iswdyp = 0
iwarnp = 10
blencp = 0.d0
epsilp = 1.d-8
epsrsp = 1.d-8
epsrgp = 1.d-5
climgp = 1.5d0
extrap = 0.d0
relaxp = 1.d0
iescap = 0
! all boundary convective flux with upwind
icvflb = 0

nomva0 = 'PoissonL'

!  IVAR = 0 (indique pour la periodicite de rotation que la variable
!     n'est pas la vitesse ni Rij)

ivar = 0

! On annule les flux de masse

do ifac = 1,nfac
  fmala(ifac) = zero
enddo

do ifac = 1,nfabor
  fmalb(ifac) = zero
enddo

! Dans le cas d'un theta-schema on met theta = 1 (ordre 1)

thetap = 1.0d0

call codits &
!==========
 ( idtva0 , ivar   , iconvp , idiffp , ndircp ,                   &
   imrgra , nswrsp , nswrgp , imligp , ircflp ,                   &
   ischcp , isstpp , iescap , imucpp , idftnp , iswdyp ,          &
   iwarnp ,                                                       &
   blencp , epsilp , epsrsp , epsrgp , climgp , extrap ,          &
   relaxp , thetap ,                                              &
   phia   , phia   , coefap , coefbp ,                            &
   cofafp , cofbfp ,                                              &
   fmala  , fmalb  ,                                              &
   viscf  , viscb  , rvoid  , viscf  , viscb  , rvoid  ,          &
   rvoid  , rvoid  ,                                              &
   icvflb , ivoid  ,                                              &
   rovsdt , smbrs  , phi    , dpvar  ,                            &
   rvoid  , rvoid  )


! Free memory
deallocate(viscf, viscb)
deallocate(smbrs, rovsdt)
deallocate(fmala, fmalb)
deallocate(coefap, coefbp)
deallocate(cofafp, cofbfp)
deallocate(phia)
deallocate(w1, w2, w3)
deallocate(dpvar)

!--------
! FORMATS
!--------

 1000 format(/,                                                   &
'   ** RESOLUTION POUR LA VARIABLE ',A19                       ,/,&
'      ---------------------------                            ',/)

 1100 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========                                               ',/,&
'@    ERREUR A LA RESOLUTION DE L''EQUATION DE POISSON :      ',/,&
'@      CONDITIONS AUX LIMITES SUR PHI NON PREVUES (LAGEQP).  ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Contacter l''equipe de developpement.                     ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

!----
! FIN
!----

return

end subroutine
