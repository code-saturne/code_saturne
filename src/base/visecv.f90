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

subroutine visecv &
!================

 ( propce ,                            &
   secvif , secvib )

!===============================================================================
! FONCTION :
! ----------

! COMPUTE (K -2/3 MU)

! IN ORDER TO COMPUTE GRAD( (K -2/3 MU) TRACE( GRAD(U)) ) + DIV( MU (GRAD_TRANSPOSE(U)) )

! WITH MU = MU_LAMINAR + MU_TURBULENT
!  AND K = VOLUME VISCOSITY (GENERALLY ZERO)

! GRAD(U) IS A CELL GRADIENT

! REMARKS :
!  - In LES, the tensor <(u-<u>)(u-<u>)> is modeled by mut <S>
!      and not by mut <S> - 2/3 mut Tr(<S>) Id + 2/3 k Id
!      so that no term mut div<u> is needed.
!-------------------------------------------------------------------------------
!ARGU                             ARGUMENTS
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! secvif(nfac)     ! tr ! --- ! lambda*surface at interior faces               !
! secvib(nfabor)   ! tr ! --- ! lambda*surface at boundary faces               !
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
use cstphy
use entsor
use numvar
use optcal
use pointe, only: porosi
use parall
use period
use ppppar
use ppthch
use ppincl
use mesh

!===============================================================================

implicit none

! Arguments

double precision propce(ncelet,*)
double precision secvif(nfac), secvib(nfabor)

! Local variables

integer          iel, ifac, ii, jj
integer          ipcvis, ipcvst
integer          ipcvsv

double precision d2s3m

double precision, allocatable, dimension(:) :: secvis

!===============================================================================

!===============================================================================
! 1.  INITIALIZATION
!===============================================================================

! Allocate temporary arrays
allocate(secvis(ncelet))

ipcvis = ipproc(iviscl)
ipcvst = ipproc(ivisct)

if(ippmod(icompf).ge.0) then
  if(iviscv.gt.0) then
    ipcvsv = ipproc(iviscv)
  else
    ipcvsv = 0
  endif
else
  ipcvsv = -1
endif

! If we extrapolate the source terms, then we take the physical properties at
! time n FIXME check if it is consistent with viscfa.f90
if(isno2t.gt.0) then
  if(iviext.gt.0) then
    ipcvis = ipproc(ivisla)
    ipcvst = ipproc(ivista)
  endif
endif

!===============================================================================
! 2. COMPUTATION OF THE SECOND VISCOSITY: LAMBDA = K -2/3 MU
!===============================================================================

!  Ici pour l'ordre 2 en temps, il faudrait tout prendre en n...

d2s3m = -2.d0/3.d0

! Laminar viscosity
do iel = 1, ncel
  secvis(iel) = d2s3m*propce(iel,ipcvis)
enddo

! Volume viscosity if present
if (ipcvsv.gt.0) then
  do iel = 1, ncel
    secvis(iel) = secvis(iel) + propce(iel,ipcvsv)
  enddo
elseif (ipcvsv.eq.0) then
  do iel = 1, ncel
    secvis(iel) = secvis(iel) + viscv0
  enddo
endif

! Turbulent viscosity (if not in Rij or LES)
if (itytur.ne.3 .and. itytur.ne.4) then
  do iel = 1, ncel
    secvis(iel) = secvis(iel) + d2s3m*propce(iel,ipcvst)
  enddo
endif

! With porosity
if (iporos.ge.1) then
  do iel = 1, ncel
    secvis(iel) = secvis(iel)*porosi(iel)
  enddo
endif


! ---> PARALLELISM AND PERIODICITY TREATMENT

if (irangp.ge.0.or.iperio.eq.1) then
  call synsca(secvis)
  !==========
endif

! --- Interior faces
! TODO we should (re)test the weigthen walue and also add a consistent
! geometrical weigthen if imvisf>0. (see viscfa.f90)

do ifac = 1, nfac
  ii = ifacel(1,ifac)
  jj = ifacel(2,ifac)
  secvif (ifac) = 0.5d0*(secvis(ii)+secvis(jj))
enddo

! --- Boundary faces
! TODO shall we extrapolate this value?

do ifac = 1, nfabor
  ii = ifabor(ifac)
  secvib(ifac) = secvis(ii)
enddo

! --- TODO stresses at the wall?

! Free memory
deallocate(secvis)

return
end subroutine
