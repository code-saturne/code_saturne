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

subroutine rijech &
!================

 ( isou   ,                                                       &
   rtpa   ,                                                       &
   produc , smbr   )

!===============================================================================
! FONCTION :
! ----------

! TERMES D'ECHO DE PAROI
!   POUR Rij
! VAR  = R11 R22 R33 R12 R13 R23
! ISOU =  1   2   3   4   5   6

!-------------------------------------------------------------------------------
!ARGU                             ARGUMENTS
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! ivar             ! i  ! <-- ! variable number                                !
! isou             ! e  ! <-- ! numero de passage                              !
! ipp              ! e  ! <-- ! numero de variable pour sorties post           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! produc           ! tr ! <-- ! production                                     !
!  (6,ncelet)      !    !     !                                                !
! smbr(ncelet      ! tr ! <-- ! tableau de travail pour sec mem                !
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
use parall
use period
use mesh
use field
!===============================================================================

implicit none

! Arguments

integer          isou

double precision rtpa(ncelet,*)
double precision produc(6,ncelet)
double precision smbr(ncelet)

! Local variables

integer          ifacpt, iel   , ii    , jj    , kk    , mm
integer          irkm  , irki  , irkj  , iskm  , iski  , iskj
integer          ifac
integer          inc   , iccocg, ivar0

double precision cmu075, distxn, d2s3  , trrij , xk
double precision unssur, vnk   , vnm   , vni   , vnj
double precision deltki, deltkj, deltkm, deltij
double precision aa    , bb    , xnorme

double precision, allocatable, dimension(:,:) :: grad
double precision, allocatable, dimension(:) :: produk, epsk
double precision, allocatable, dimension(:) :: w2, w3, w4, w6
double precision, allocatable, dimension(:) :: coefax, coefbx
double precision, dimension(:), pointer :: crom, cromo

!===============================================================================

!===============================================================================
! 1. INITIALISATION
!===============================================================================

! Allocate temporary arrays
allocate(produk(ncelet), epsk(ncelet))
allocate(w2(ncelet), w3(ncelet), w4(ncelet), w6(ncelet))

! Initialize variables to avoid compiler warnings

ii = 0
jj = 0
irki = 0
irkj = 0
irkm = 0
iski = 0
iskj = 0
iskm = 0

vni = 0.d0
vnj = 0.d0
vnk = 0.d0
vnm = 0.d0

! Memoire

call field_get_val_s(icrom, crom)
if (isto2t.gt.0.and.iroext.gt.0) then
  call field_get_val_prev_s(icrom, cromo)
else
  call field_get_val_s(icrom, cromo)
endif

deltij = 1.0d0
if(isou.gt.3) then
  deltij = 0.0d0
endif

cmu075 = cmu**0.75d0
d2s3   = 2.d0/3.d0

!===============================================================================
! 2. CALCUL AUX CELLULES DES NORMALES EN PAROI CORRESPONDANTES
!===============================================================================

!     On stocke les composantes dans les tableaux de travail W2, W3 et W4

if(abs(icdpar).eq.2) then

!     On connait la face de paroi correspondante

    do iel = 1, ncel
      ifacpt = ifapat(iel)
      unssur = 1.d0/surfbn(ifacpt)
      w2(iel)= surfbo(1,ifacpt)*unssur
      w3(iel)= surfbo(2,ifacpt)*unssur
      w4(iel)= surfbo(3,ifacpt)*unssur
    enddo

elseif(abs(icdpar).eq.1) then

!     La normale est definie comme - gradient de la distance
!       a la paroi

!       La distance a la paroi vaut 0 en paroi
!         par definition et obeit a un flux nul ailleurs

  ! Allocate temporary arrays
  allocate(coefax(nfabor), coefbx(nfabor))

  do ifac = 1, nfabor
    if (itypfb(ifac).eq.iparoi .or. itypfb(ifac).eq.iparug) then
      coefax(ifac) = 0.0d0
      coefbx(ifac) = 0.0d0
    else
      coefax(ifac) = 0.0d0
      coefbx(ifac) = 1.0d0
    endif
  enddo

!       Calcul du gradient

  ! Allcoate a temporary array for the gradient calculation
  allocate(grad(ncelet,3))

  if (irangp.ge.0.or.iperio.eq.1) then
    call synsca(dispar)
    !==========
  endif

  inc    = 1
  iccocg = 1
  ivar0  = 0

  call grdcel                                                     &
  !==========
 ( ivar0  , imrgra , inc    , iccocg , nswrgy , imligy ,          &
   iwarny , nfecra ,                                              &
   epsrgy , climgy , extray ,                                     &
   dispar , coefax , coefbx ,                                     &
   grad   )

  ! Free memory
  deallocate(coefax, coefbx)

!     Normalisation (attention, le gradient peut etre nul, parfois)

  do iel = 1 ,ncel
    xnorme = max(sqrt(grad(iel,1)**2+grad(iel,2)**2+grad(iel,3)**2),epzero)
    w2(iel) = -grad(iel,1)/xnorme
    w3(iel) = -grad(iel,2)/xnorme
    w4(iel) = -grad(iel,3)/xnorme
  enddo

  ! Free memory
  deallocate(grad)

endif

!===============================================================================
! 2. CALCUL DE VARIABLES DE TRAVAIL
!===============================================================================



! ---> Production et k

do iel = 1 , ncel
  produk(iel) = 0.5d0 * (produc(1,iel)  + produc(2,iel)  + produc(3,iel))
  xk          = 0.5d0 * (rtpa(iel,ir11) + rtpa(iel,ir22) + rtpa(iel,ir33))
  epsk(iel)   = rtpa(iel,iep)/xk
enddo



! ---> Indices des tensions

if     ((isou.eq.1).or.(isou.eq.4).or.(isou.eq.5)) then
  ii = 1
elseif ((isou.eq.2).or.(isou.eq.6)) then
  ii = 2
elseif  (isou.eq.3) then
  ii = 3
endif

if     ((isou.eq.3).or.(isou.eq.5).or.(isou.eq.6)) then
  jj = 3
elseif ((isou.eq.2).or.(isou.eq.4)) then
  jj = 2
elseif  (isou.eq.1) then
  jj = 1
endif

! ---> Boucle pour construction des termes sources

do iel = 1, ncel
  w6(iel) = 0.d0
enddo

do kk = 1, 3

! ---> Sommes sur m

  do mm = 1, 3

!   --> Delta km

    if(kk.eq.mm) then
      deltkm = 1.0d0
    else
      deltkm = 0.0d0
    endif

!  --> R km

    if     ((kk*mm).eq.1) then
      irkm = ir11
      iskm = 1
    elseif ((kk*mm).eq.4) then
      irkm = ir22
      iskm = 2
    elseif ((kk*mm).eq.9) then
      irkm = ir33
      iskm = 3
    elseif ((kk*mm).eq.2) then
      irkm = ir12
      iskm = 4
    elseif ((kk*mm).eq.3) then
      irkm = ir13
      iskm = 5
    elseif ((kk*mm).eq.6) then
      irkm = ir23
      iskm = 6
    endif

!  --> Termes en R km et Phi km

    do iel = 1, ncel

      if    (kk.eq.1) then
        vnk    = w2(iel)
      elseif(kk.eq.2) then
        vnk    = w3(iel)
      elseif(kk.eq.3) then
        vnk    = w4(iel)
      endif

      if    (mm.eq.1) then
        vnm    = w2(iel)
      elseif(mm.eq.2) then
        vnm    = w3(iel)
      elseif(mm.eq.3) then
        vnm    = w4(iel)
      endif

      w6(iel) = w6(iel) + vnk*vnm*deltij*(                        &
             crijp1*rtpa(iel,irkm)*epsk(iel)                      &
            -crijp2                                               &
             *crij2*(produc(iskm,iel)-d2s3*produk(iel)*deltkm) )
    enddo

  enddo

! ---> Fin des sommes sur m


!  --> R ki

  if     ((kk*ii).eq.1) then
    irki = ir11
    iski = 1
  elseif ((kk*ii).eq.4) then
    irki = ir22
    iski = 2
  elseif ((kk*ii).eq.9) then
    irki = ir33
    iski = 3
  elseif ((kk*ii).eq.2) then
    irki = ir12
    iski = 4
  elseif ((kk*ii).eq.3) then
    irki = ir13
    iski = 5
  elseif ((kk*ii).eq.6) then
    irki = ir23
    iski = 6
  endif

!  --> R kj

  if     ((kk*jj).eq.1) then
    irkj = ir11
    iskj = 1
  elseif ((kk*jj).eq.4) then
    irkj = ir22
    iskj = 2
  elseif ((kk*jj).eq.9) then
    irkj = ir33
    iskj = 3
  elseif ((kk*jj).eq.2) then
    irkj = ir12
    iskj = 4
  elseif ((kk*jj).eq.3) then
    irkj = ir13
    iskj = 5
  elseif ((kk*jj).eq.6) then
    irkj = ir23
    iskj = 6
  endif

!   --> Delta ki

  if (kk.eq.ii) then
    deltki = 1.d0
  else
    deltki = 0.d0
  endif

!   --> Delta kj

  if (kk.eq.jj) then
    deltkj = 1.d0
  else
    deltkj = 0.d0
  endif

  do iel = 1, ncel

      if    (kk.eq.1) then
        vnk    = w2(iel)
      elseif(kk.eq.2) then
        vnk    = w3(iel)
      elseif(kk.eq.3) then
        vnk    = w4(iel)
      endif

      if    (ii.eq.1) then
        vni    = w2(iel)
      elseif(ii.eq.2) then
        vni    = w3(iel)
      elseif(ii.eq.3) then
        vni    = w4(iel)
      endif

      if    (jj.eq.1) then
        vnj    = w2(iel)
      elseif(jj.eq.2) then
        vnj    = w3(iel)
      elseif(jj.eq.3) then
        vnj    = w4(iel)
      endif

    w6(iel) = w6(iel) + 1.5d0*vnk*(                               &
    -crijp1*(rtpa(iel,irki)*vnj+rtpa(iel,irkj)*vni)*epsk(iel)     &
    +crijp2                                                       &
     *crij2*((produc(iski,iel)-d2s3*produk(iel)*deltki)*vnj       &
            +(produc(iskj,iel)-d2s3*produk(iel)*deltkj)*vni) )

  enddo

enddo


! ---> Distance a la paroi et fonction d'amortissement : W3
!     Pour chaque mode de calcul : meme code, test
!       en dehors de la boucle

if(abs(icdpar).eq.2) then
  do iel = 1 , ncel
    ifacpt = ifapat(iel)
    distxn =                                                      &
          (cdgfbo(1,ifacpt)-xyzcen(1,iel))**2                     &
         +(cdgfbo(2,ifacpt)-xyzcen(2,iel))**2                     &
         +(cdgfbo(3,ifacpt)-xyzcen(3,iel))**2
    distxn = sqrt(distxn)
    trrij  = 0.5d0 * (rtpa(iel,ir11) + rtpa(iel,ir22) + rtpa(iel,ir33))
    aa = 1.d0
    bb = cmu075*trrij**1.5d0/(xkappa*rtpa(iel,iep)*distxn)
    w3(iel) = min(aa, bb)
  enddo
else
  do iel = 1 , ncel
    distxn =  max(dispar(iel),epzero)
    trrij  = 0.5d0 * (rtpa(iel,ir11) + rtpa(iel,ir22) + rtpa(iel,ir33))
    aa = 1.d0
    bb = cmu075*trrij**1.5d0/(xkappa*rtpa(iel,iep)*distxn)
    w3(iel) = min(aa, bb)
  enddo
endif


! ---> Increment du terme source

do iel = 1, ncel
  smbr(iel) = smbr(iel) + cromo(iel)*volume(iel)*w6(iel)*w3(iel)
enddo

! Allocate temporary arrays
deallocate(produk, epsk)
deallocate(w2, w3, w4, w6)


return
end subroutine
