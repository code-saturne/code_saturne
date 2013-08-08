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

subroutine cspcev &
!================

 ( nptdis , ivar   ,                                              &
   locpts ,                                                       &
   vela   ,                                                       &
   coefav , coefbv ,                                              &
   coopts , rvdis  )

!===============================================================================
! FONCTION :
! --------

! PREPARATION DE L'ENVOI DES VARIABLES DE VITESSE POUR UN COUPLAGE
!   ENTRE DEUX INSTANCES DE CODE_SATURNE VIA LES FACES DE BORD

! L'INFORMATION RECUE SERA TRANSFORMEE EN CONDITION LIMITE DANS
!   LA SUBROUTINE CSC2CL

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
use cplsat
use mesh

!===============================================================================

implicit none

! Arguments

integer          ivar
integer          nptdis

integer          locpts(nptdis)

double precision coopts(3,nptdis), rvdis(3,nptdis)
double precision coefav(3  ,nfabor)
double precision coefbv(3,3,nfabor)
double precision vela(3,ncelet)

! Local variables

integer          ipt    , iel    , isou
integer          inc    , iccocg , iclvar, nswrgp
integer          iwarnp , imligp

double precision epsrgp , climgp , extrap
double precision dx     , dy     , dz

logical ilved

double precision, dimension(:,:,:), allocatable :: gradv

!===============================================================================

! Allocate a temporary array
allocate(gradv(3,3,ncelet))


if (irangp.ge.0.or.iperio.eq.1) then
  call synvin(vela)
  !==========
endif

inc    = 1
iccocg = 1
iclvar = iclrtp(ivar,icoef)
nswrgp = nswrgr(ivar)
imligp = imligr(ivar)
iwarnp = iwarni(ivar)
epsrgp = epsrgr(ivar)
climgp = climgr(ivar)
extrap = extrag(ivar)
ilved  = .true.

call grdvec &
!==========
( ivar   , imrgra , inc    , nswrgp , imligp ,                   &
  iwarnp , nfecra , epsrgp , climgp , extrap ,                   &
  ilved  ,                                                       &
  vela   , coefav , coefbv ,                                     &
  gradv)


! --- Interpolation

do ipt = 1, nptdis

  iel = locpts(ipt)

  dx = coopts(1,ipt) - xyzcen(1,iel)
  dy = coopts(2,ipt) - xyzcen(2,iel)
  dz = coopts(3,ipt) - xyzcen(3,iel)

  do isou = 1, 3
    rvdis(isou,ipt) = vela(isou,iel) + gradv(isou,1,iel)*dx       &
                                     + gradv(isou,2,iel)*dy       &
                                     + gradv(isou,3,iel)*dz
  enddo

enddo

! Free memory
deallocate(gradv)

!--------
! FORMATS
!--------
!----
! FIN
!----

return
end subroutine
