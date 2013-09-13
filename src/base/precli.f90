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

!> \file precli.f90
!> \brief Preparation of boudary conditions determination
!> Boundary faces of precedent step are used.
!> Except at first time step, where arrays \ref itypfb and \ref itrifb
!> are undefined.
!>
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Arguments
!------------------------------------------------------------------------------
!   mode          name          role
!------------------------------------------------------------------------------
!> \param[in]     nvar          total number of variables
!> \param[in]     nscal         total number of scalars
!> \param[out]    icodcl        defines the type of boundary condition for the
!>                              variable ivar on the face
!>                               - 1: Dirichlet condition at the face
!>                               - 3: flux condition at the face
!>                               - 4: symmetry condition, for the symmetry
!>                              faces or wall faces without friction.
!>                              Only used for velocity components
!>                               - 5: friction condition, for wall faces
!>                              with friction. This condition can not be applied
!>                              to the pressure.
!>                               - 6: friction condition, for the rough-wall
!>                              faces with friction. This condition can not be
!>                              used with the pressure.
!>                               - 9: free outlet condition for the
!>                              velocity. Only applicable to velocity components
!> \param[in]     propfb        physical properties at boundary face centers
!> \param[out]    rcodcl        gives the numerical values associated with the
!>                              type of boundary condition
!>                              (value of the Dirichlet, of the flux ...).
!>                               - rcodcl(1) = Dirichlet value
!>                               - rcodcl(2) = value of the exchange coefficient
!>                                between the outside and the fluid. Infinite
!>                                value indicates an ideal transfer(default case)
!>                               - rcodcl(3):
!>                                 - if icodcl=6: rugosity height (m)
!>                                 - else value of flux density (w/m2).
!>                                            (negative if gain)
!>                               - for velocity: (vistl+visct)*gradu
!>                               - for pressure:             dt*gradp
!>                               - for scalars:  cp*(viscls+visct/sigmas)*gradt
!______________________________________________________________________________

subroutine precli &
 ( nvar   , nscal  ,                                              &
   icodcl ,                                                       &
   propfb ,                                                       &
   rcodcl )

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use optcal
use cstphy
use cstnum
use entsor
use pointe
use albase
use ppppar
use ppthch
use ppincl
use mesh

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal

integer          icodcl(nfabor,nvarcl)

double precision propfb(nfabor,*)
double precision rcodcl(nfabor,nvarcl,3)

! Local variables

integer          ifac, ivar, iscal, iut, ivt, iwt

!===============================================================================
!===============================================================================
! 1.  INITIALISATIONS
!===============================================================================


!===============================================================================
! 2.  INITIALISATION DES CONDITIONS LIMITES ET TYPE DE FACES DE BORD
!===============================================================================

!      ICODCL = 0 INDIQUE QUE LA CL N'A PAS ETE RENSEIGNEE
!      ITYPFB = 0 INDIQUE QUE LA CL N'A PAS ETE RENSEIGNEE

!      RINFIN : VALEUR INFINIE

do ifac = 1, nfabor
  itypfb(ifac) = 0
enddo

! Pour toutes les variables, on initialise RCODCL(1)a RINFIN
! Cette valeur sera reinitialisee a zero dans typecl.F

do ivar = 1, nvar
  do ifac = 1, nfabor
    icodcl(ifac,ivar)   = 0
    rcodcl(ifac,ivar,1) = rinfin
    rcodcl(ifac,ivar,2) = rinfin
    rcodcl(ifac,ivar,3) = 0.d0
  enddo
enddo

! Default value for turbulent fluxes
do iscal = 1, nscal
  if (ityturt(iscal).eq.3) then
    iut = nvar + 3*(ifltur(iscal) - 1) + 1
    ivt = nvar + 3*(ifltur(iscal) - 1) + 2
    iwt = nvar + 3*(ifltur(iscal) - 1) + 3
    do ifac = 1, nfabor
      icodcl(ifac,iut)   = 0
      rcodcl(ifac,iut,1) = rinfin
      rcodcl(ifac,iut,2) = rinfin
      rcodcl(ifac,iut,3) = 0.d0
      icodcl(ifac,ivt)   = 0
      rcodcl(ifac,ivt,1) = rinfin
      rcodcl(ifac,ivt,2) = rinfin
      rcodcl(ifac,ivt,3) = 0.d0
      icodcl(ifac,iwt)   = 0
      rcodcl(ifac,iwt,1) = rinfin
      rcodcl(ifac,iwt,2) = rinfin
      rcodcl(ifac,iwt,3) = 0.d0
    enddo
  endif
enddo

! En ALE, on initialise aussi le tableau IALTYB
if (iale.eq.1) then
  do ifac = 1, nfabor
    ialtyb(ifac) = 0
  enddo
endif

! POUR LES PHYSIQUES PARTICULIERES

if (ippmod(iphpar).ge.1) then
  call ppprcl(nvar, izfppp, propfb, rcodcl)
  !==========
endif

!----
! End
!----

return
end subroutine
