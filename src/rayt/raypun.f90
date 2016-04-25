!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2016 EDF S.A.
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
 ( itypfb ,                                                       &
   coefap , coefbp ,                                              &
   cofafp , cofbfp ,                                              &
   flurds , flurdb ,                                              &
   viscf  , viscb  ,                                              &
   smbrs  , rovsdt ,                                              &
   theta4 , thetaa , sa     ,                                     &
   q   ,                                                          &
   qincid , eps    , tparoi ,                                     &
   ckmel  , abo    , iband )

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
! itypfb           ! ia ! <-- ! boundary face types                            !
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
! q(3,ncelet)      ! tr ! --> ! vecteur densite de flux radiatif explicite     !
! qincid(nfabor    ! tr ! --> ! densite de flux radiatif aux bords             !
! eps (nfabor)     ! tr ! <-- ! emissivite des facettes de bord                !
! tparoi(nfabor    ! tr ! <-- ! temperature de paroi en kelvin                 !
! ckmel(ncelet)    ! tr ! <-- ! coeff d'absorption du melange                  !
!                  !    !     !   gaz-particules de charbon                    !
! abo              ! ra ! <-- ! Weights of the i-th gray gas at boundaries     !
! iband            ! i  ! <-- ! Number of the i-th grey gas                    !
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
use cs_c_bindings
use field

!===============================================================================

implicit none

! Arguments

integer          iband
integer          itypfb(nfabor)

double precision coefap(nfabor), coefbp(nfabor)
double precision cofafp(nfabor), cofbfp(nfabor)
double precision flurds(nfac), flurdb(nfabor)

double precision viscf(nfac), viscb(nfabor)
double precision smbrs(ncelet)
double precision rovsdt(ncelet)

double precision theta4(ncelet), thetaa(ncelet)
double precision sa(ncelet)
double precision q(3,ncelet)
double precision qincid(nfabor), tparoi(nfabor), eps(nfabor)

double precision ckmel(ncelet)
double precision abo(nfabor,nwsgg)

! Local variables

character(len=80) :: cnom

integer          ifac  , iel
integer          iconv1, idiff1, ndirc1
integer          nswrsp, nswrgp, iwarnp
integer          imligp, ircflp, ischcp, isstpp, iescap
integer          idtva0, ivar0,f_id0
integer          inc, iccocg
integer          imucpp, idftnp, iswdyp, icvflb
integer          ivoid(1)
double precision epsrgp, blencp, climgp, epsilp, extrap, epsrsp
double precision aa, aaa, relaxp, thetap

double precision rvoid(1)

double precision, allocatable, dimension(:) :: dpvar
double precision, dimension(:,:), pointer     :: qinspe

!===============================================================================
if (imoadf.ge.1) then
  ! Pointer to the table
  ! which contains the spectral flux density
  call field_get_val_v(iqinsp,qinspe)
endif
!===============================================================================
! 0. GESTION MEMOIRE
!===============================================================================

! Allocate temporary array
allocate(dpvar(ncelet))

!===============================================================================
! 1. PARAMETRAGE DU SOLVEUR ET INITIALISATION
!===============================================================================

!--> Parametrage de CODITS

! IVAR0= 0  LA VARIABLE N'EST ICI NI RIJ NI VITESSE
ivar0   = 0
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
iwarnp  = iimlum
blencp  = zero
epsilp  = 1.d-8
epsrsp  = 1.d-8
epsrgp  = 1.d-5
climgp  = 1.5d0
extrap  = zero
relaxp  = 1.d0
! all boundary convective flux with upwind
icvflb = 0

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

call viscfa(imvisf, ckmel, viscf, viscb)

!===============================================================================
! 3.  RESOLUTION
!===============================================================================

! Parametre pour schemas en temps et stationnaire
thetap = 1.d0
idtva0 = 0

cnom = ' '
write(cnom,'(a)') 'radiation_p1'
nomva0 = cnom

call codits &
 ( idtva0 , ivar0  , iconv1 , idiff1 , ndirc1 ,                   &
   imrgra , nswrsp , nswrgp , imligp , ircflp ,                   &
   ischcp , isstpp , iescap , imucpp , idftnp , iswdyp ,          &
   iwarnp ,                                                       &
   blencp , epsilp , epsrsp , epsrgp , climgp , extrap ,          &
   relaxp , thetap ,                                              &
   thetaa , thetaa , coefap , coefbp , cofafp , cofbfp ,          &
   flurds , flurdb ,                                              &
   viscf  , viscb  , viscf  , viscb  , rvoid  ,                   &
   rvoid  , rvoid  ,                                              &
   icvflb , ivoid  ,                                              &
   rovsdt , smbrs  , theta4 , dpvar  ,                            &
   rvoid  , rvoid  )

!===============================================================================
! 4. Vecteur densite de flux radiatif
!===============================================================================

!     Calcul de la densite du flux radiatif Q

inc     = 1
iccocg  = 1
imligp  = -1
iwarnp  = iimlum
epsrgp  = 1.d-8
climgp  = 1.5d0
extrap  = 0.d0
nswrgp  = 100
f_id0   = -1

call gradient_s                                                   &
   ( f_id0  , imrgra , inc    , iccocg , nswrgp , imligp,         &
     iwarnp , epsrgp , climgp , extrap ,                          &
     theta4 , coefap , coefbp ,                                   &
     q )

aa = - stephn * 4.d0 / 3.d0

do iel = 1, ncel
  aaa = aa * ckmel(iel)
  q(1,iel) = q(1,iel) * aaa
  q(2,iel) = q(2,iel) * aaa
  q(3,iel) = q(3,iel) * aaa
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
! FIXME
!   aaaa = tparoi(ifac)**4

!   aaa  = 1.5d0 * distb(ifac) / ckmel(iel)                 &
!          * ( 2.d0 /(2.d0-eps(ifac)) -1.d0 )
!   aa   = ( aaa * aaaa + theta4(iel) ) / (1.d0 + aaa)

!   qincid(ifac) = stephn * (2.d0 * aa - eps(ifac) * aaaa)        &
!                      / (2.d0 - eps(ifac))

!--> Deuxieme version plus cheap mais moins precise
!    This expression can be found in the documentation of the P1 model
!    written by S. DAL-SECCO, A. DOUCE, N. MECHITOUA.
    if (imoadf.ge.1) then
      qinspe(iband,ifac) = stephn                                          &
                         * ((2.d0*theta4(iel))+(abo(ifac,iband)*eps(ifac)  &
                            *(tparoi(ifac)**4)))                           &
                         / (2.d0-eps(ifac))
    else
      qincid(ifac) = stephn                                              &
                   * ((2.d0*theta4(iel))+(eps(ifac)*(tparoi(ifac)**4)))  &
                   / (2.d0-eps(ifac))
    endif

  else

    if (imoadf.ge.1) then
      qinspe(iband,ifac) = stephn * theta4(iel)                         &
                         + ( q(1,iel) * surfbo(1,ifac) +                &
                             q(2,iel) * surfbo(2,ifac) +                &
                             q(3,iel) * surfbo(3,ifac) ) /              &
                             (0.5d0 * surfbn(ifac))
    else
      qincid(ifac) = stephn * theta4(iel)                               &
                 + ( q(1,iel) * surfbo(1,ifac) +                        &
                     q(2,iel) * surfbo(2,ifac) +                        &
                     q(3,iel) * surfbo(3,ifac) ) /                      &
                     (0.5d0 * surfbn(ifac) )
    endif
  endif

enddo

!===============================================================================

! Free memory
deallocate(dpvar)

!--------
! Formats
!--------

!----
! End
!----

return

end subroutine
