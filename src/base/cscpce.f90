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

subroutine cscpce &
!================

 ( idbia0 , idbra0 ,                                              &
   nvar   , nscal  ,                                              &
   nptdis , ityloc ,                                              &
   ivar   ,                                                       &
   locpts ,                                                       &
   ia     ,                                                       &
   dt     , rtpa   , propce , propfa , propfb ,                   &
   coefa  , coefb  ,                                              &
   w1     , w2     , w3     , w4     , w5     , w6     ,          &
   coopts , rvdis  ,                                              &
   ra     )

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

integer          idbia0 , idbra0
integer          nvar   , nscal
integer          ivar
integer          nptdis , ityloc

integer          locpts(nptdis)
integer          ia(*)


double precision dt(ncelet), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision w1(ncelet),w2(ncelet),w3(ncelet)
double precision w4(ncelet),w5(ncelet),w6(ncelet)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision coopts(3,nptdis), rvdis(nptdis)
double precision ra(*)

! Local variables

integer          idebia , idebra , ifinia , ifinra
integer          ipt    , iel
integer          inc    , iccocg , iphydp , iclvar, nswrgp
integer          iwarnp , imligp

double precision epsrgp , climgp , extrap
double precision dx     , dy     , dz

!===============================================================================

idebia = idbia0
idebra = idbra0

if (irangp.ge.0.or.iperio.eq.1) then
  call synsca(rtpa(1,ivar))
  !==========
endif

inc    = 1
iccocg = 1
iphydp = 0
iclvar = iclrtp(ivar,icoef)
nswrgp = nswrgr(ivar)
imligp = imligr(ivar)
iwarnp = iwarni(ivar)
epsrgp = epsrgr(ivar)
climgp = climgr(ivar)
extrap = extrag(ivar)

call grdcel                                                       &
!==========
 ( ivar   , imrgra , inc    , iccocg , nswrgp , imligp , iphydp,  &
   iwarnp , nfecra ,                                              &
   epsrgp , climgp , extrap ,                                     &
   ia     ,                                                       &
   w4     , w4     , w4     ,                                     &
   rtpa(1,ivar)    , coefa(1,iclvar) , coefb(1,iclvar) ,          &
   w1     , w2     , w3     ,                                     &
!        ------   ------   ------
   w4     , w5     , w6     ,                                     &
   ra     )


! --- Interpolation

do ipt = 1, nptdis

  iel = locpts(ipt)

  dx = coopts(1,ipt) - xyzcen(1,iel)
  dy = coopts(2,ipt) - xyzcen(2,iel)
  dz = coopts(3,ipt) - xyzcen(3,iel)

  rvdis(ipt) = rtpa(iel,ivar) + w1(iel)*dx+w2(iel)*dy+w3(iel)*dz

enddo


!--------
! FORMATS
!--------
!----
! FIN
!----

return
end subroutine
