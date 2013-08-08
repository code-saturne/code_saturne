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

subroutine cscpce &
!================

 ( nptdis , ivar   ,                                              &
   locpts ,                                                       &
   rtpa   ,                                                       &
   coefa  , coefb  ,                                              &
   coopts , rvdis  )

!===============================================================================
! FONCTION :
! --------

! PREPARATION DE L'ENVOI DES VARIABLES POUR UN COUPLAGE
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

double precision rtpa(ncelet,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision coopts(3,nptdis), rvdis(nptdis)

! Local variables

integer          ipt    , iel
integer          inc    , iccocg , iclvar, nswrgp
integer          iwarnp , imligp

double precision epsrgp , climgp , extrap
double precision dx     , dy     , dz

double precision, allocatable, dimension(:,:) :: grad

!===============================================================================

! Allocate a temporary array
allocate(grad(ncelet,3))


if (irangp.ge.0.or.iperio.eq.1) then
  call synsca(rtpa(1,ivar))
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

call grdcel                                                       &
!==========
 ( ivar   , imrgra , inc    , iccocg , nswrgp , imligp ,          &
   iwarnp , nfecra ,                                              &
   epsrgp , climgp , extrap ,                                     &
   rtpa(1,ivar)    , coefa(1,iclvar) , coefb(1,iclvar) ,          &
   grad   )


! --- Interpolation

do ipt = 1, nptdis

  iel = locpts(ipt)

  dx = coopts(1,ipt) - xyzcen(1,iel)
  dy = coopts(2,ipt) - xyzcen(2,iel)
  dz = coopts(3,ipt) - xyzcen(3,iel)

  rvdis(ipt) = rtpa(iel,ivar) + grad(iel,1)*dx+grad(iel,2)*dy+grad(iel,3)*dz

enddo

! Free memory
deallocate(grad)

!--------
! FORMATS
!--------
!----
! FIN
!----

return
end subroutine
