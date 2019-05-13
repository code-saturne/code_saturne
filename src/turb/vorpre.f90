!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2019 EDF S.A.
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

subroutine vorpre

!===============================================================================
! FONCTION :
! --------

!    ROUTINE DE PREPATATION DE LA METHODE DES VORTEX
!    Gestion memoire, connectivites, ...
!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
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
use pointe
use numvar
use optcal
use cstphy
use cstnum
use entsor
use parall
use period
use vorinc
use mesh
use field
!===============================================================================

implicit none

! Arguments

! Local variables

integer          ifac, iel, ii
integer          ient
integer          isurf(nentmx)

double precision xx, yy, zz
double precision xxv, yyv, zzv

double precision, allocatable, dimension(:,:) :: w1x, w1y, w1z, w1v
double precision, dimension(:), pointer :: crom
double precision, dimension(:), pointer :: viscl

!===============================================================================
! 1.  INITIALISATIONS
!===============================================================================

! Allocate temporary arrays
allocate(w1x(icvmax,nnent), w1y(icvmax,nnent), w1z(icvmax,nnent))
allocate(w1v(icvmax,nnent))


nvomax = 0
do ient = 1, nnent
  nvomax = max(nvort(ient),nvomax)
enddo

! NVOMAX = nombre max de vortex (utilise plus tard)

do ient = 1, nnent
  icvor2(ient) = 0
enddo

do ifac = 1, nfabor
  ient = irepvo(ifac)
  if(ient.ne.0) then
    icvor2(ient) = icvor2(ient) + 1
  endif
enddo

! ICVOR2 = compteur du nombre local de faces
!   utilisant des vortex a l'entree IENT

icvmax = 0
if(irangp.ge.0) then
  do ient = 1, nnent
    icvor(ient) = icvor2(ient)
    call parcpt(icvor(ient))
    icvmax = max(icvmax,icvor(ient))
  enddo
else
  do ient = 1, nnent
    icvor(ient) = icvor2(ient)
    icvmax = max(icvmax,icvor(ient))
  enddo
endif

!===============================================================================
! 2. CONSTRUCTION DE LA " GEOMETRIE GOBALE "
!===============================================================================

do ient = 1, nnent
  icvor2(ient) = 0
  xsurfv(ient) = 0.d0
  isurf(ient)  = 0
enddo

! Chaque processeur stocke dans les tableaux 'w1x', ...
! les coordonnees des faces ou il doit ensuite utiliser des vortex

call field_get_val_s(iviscl, viscl)
call field_get_val_s(icrom, crom)
do ifac = 1, nfabor
  ient = irepvo(ifac)
  if(ient.ne.0) then
    iel = ifabor(ifac)
    icvor2(ient) = icvor2(ient) + 1
    w1x(icvor2(ient),ient)= cdgfbo(1,ifac)
    w1y(icvor2(ient),ient)= cdgfbo(2,ifac)
    w1z(icvor2(ient),ient)= cdgfbo(3,ifac)
    w1v(icvor2(ient),ient) = viscl(iel)/crom(iel)
    xsurfv(ient) = xsurfv(ient) + sqrt(surfbo(1,ifac)**2          &
      + surfbo(2,ifac)**2 + surfbo(3,ifac)**2)
!         Vecteur surface d'une face de l'entree
    if (isurf(ient).eq.0) then
      surf(1,ient) = surfbo(1,ifac)
      surf(2,ient) = surfbo(2,ifac)
      surf(3,ient) = surfbo(3,ifac)
      isurf(ient)  = 1
    endif
  endif
enddo

if(irangp.ge.0) then
  do ient = 1, nnent
    call parsom(xsurfv(ient))
  enddo
endif

! -------------
! En parallele
! -------------
if(irangp.ge.0) then
  do ient = 1, nnent
    call paragv &
 ( icvor2(ient) , icvor(ient)    ,  &
   w1x(:,ient)  , xyzv(:,1,ient) )
    call paragv &
 ( icvor2(ient) , icvor(ient)    ,  &
   w1y(:,ient)  , xyzv(:,2,ient) )
    call paragv &
 ( icvor2(ient) , icvor(ient)    ,  &
   w1z(:,ient)  , xyzv(:,3,ient) )
    call paragv &
 ( icvor2(ient) , icvor(ient)  ,  &
   w1v(:,ient)  , visv(:,ient) )
  enddo

!  -> A la fin de cette etape, tous les processeurs connaissent
!     les coordonees des faces d'entree

else
! ----------------------
! Sur 1 seul processeur
! ----------------------
  do ient = 1,nnent
    do ii = 1, icvor(ient)
      xyzv(ii,1,ient) = w1x(ii,ient)
      xyzv(ii,2,ient) = w1y(ii,ient)
      xyzv(ii,3,ient) = w1z(ii,ient)
      visv(ii,ient) = w1v(ii,ient)
    enddo
  enddo
endif

!===============================================================================
! 3. CONSTRUCTION DE LA CONNECTIVITE
!===============================================================================

do ient = 1, nnent
  icvor2(ient) = 0
  do ifac = 1, icvmax
    ifacgl(ifac,ient) = 0
  enddo
enddo

! On cherche ensuite le numero de la ligne du tableau 'xyzv' qui est
! associe a la Ieme face d'entree utilisant des vortex (dans la
! numerotation chronologique que suit ICVOR2).

do ifac = 1, nfabor
  ient = irepvo(ifac)
  if(ient.ne.0) then
    icvor2(ient) = icvor2(ient) + 1
    do ii = 1, icvor(ient)
      xx = cdgfbo(1,ifac)
      yy = cdgfbo(2,ifac)
      zz = cdgfbo(3,ifac)
      xxv = xyzv(ii,1,ient)
      yyv = xyzv(ii,2,ient)
      zzv = xyzv(ii,3,ient)
      if(abs(xxv-xx).lt.epzero.and.abs(yyv-yy).lt.epzero.and.     &
           abs(zzv-zz).lt.epzero) then
        ifacgl(icvor2(ient),ient) = ii
      endif
    enddo
  endif
enddo

! Allocate temporary arrays
deallocate(w1x, w1y, w1z)
deallocate(w1v)

! ---
! FIN
! ---

return
end subroutine
