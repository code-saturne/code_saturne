!-------------------------------------------------------------------------------

!VERS

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

!===============================================================================
! Function:
! ---------

!> \file cs_user_boundary_conditions_ale.f90
!>
!> \brief User subroutine dedicated the use of ALE (Arbitrary Lagrangian
!> Eulerian) Method:
!>  - Fills boundary conditions (ialtyb, icodcl, rcodcl) for mesh velocity.
!>  - This subroutine also enables one to fix displacement on nodes.
!>
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role
!______________________________________________________________________________!
!> \param[in]     itrale        number of iterations for ALE method
!> \param[in]     nvar          total number of variables
!> \param[in]     nscal         total number of scalars
!> \param[out]    icodcl        boundary condition code:
!>                               - 1 Dirichlet
!>                               - 2 Radiative outlet
!>                               - 3 Neumann
!>                               - 4 sliding and
!>                                 \f$ \vect{u} \cdot \vect{n} = 0 \f$
!>                               - 5 smooth wall and
!>                                 \f$ \vect{u} \cdot \vect{n} = 0 \f$
!>                               - 6 rough wall and
!>                                 \f$ \vect{u} \cdot \vect{n} = 0 \f$
!>                               - 9 free inlet/outlet
!>                                 (input mass flux blocked to 0)
!> \param[in,out] itypfb        boundary face types
!> \param[out]    ialtyb        boundary face types for mesh velocity
!> \param[in]     impale        indicator for fixed node displacement
!> \param[in]     dt            time step (per cell)
!> \param[in,out] rcodcl        boundary condition values:
!>                               - rcodcl(1) value of the dirichlet
!>                               - rcodcl(2) value of the exterior exchange
!>                                 coefficient (infinite if no exchange)
!>                               - rcodcl(3) value flux density
!>                                 (negative if gain) in w/m2 or roughness
!>                                 in m if icodcl=6
!>                                 -# for the velocity \f$ (\mu+\mu_T)
!>                                    \gradv \vect{u} \cdot \vect{n}  \f$
!>                                 -# for the pressure \f$ \Delta t
!>                                    \grad P \cdot \vect{n}  \f$
!>                                 -# for a scalar \f$ C_p \left( K +
!>                                     \dfrac{K_T}{\sigma_T} \right)
!>                                     \grad T \cdot \vect{n} \f$
!> \param[in,out] disale        nodes displacement
!> \param[in]     xyzno0        vertex coordinates of initial mesh
!_______________________________________________________________________________

subroutine usalcl &
 ( itrale ,                                                       &
   nvar   , nscal  ,                                              &
   icodcl , itypfb , ialtyb , impale ,                            &
   dt     ,                                                       &
   rcodcl , xyzno0 , disale )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use optcal
use cstphy
use cstnum
use entsor
use parall
use period
use mesh

!===============================================================================

implicit none
!< [arg]
! Arguments

integer          itrale
integer          nvar   , nscal

integer          icodcl(nfabor,nvar)
integer          itypfb(nfabor), ialtyb(nfabor)
integer          impale(nnod)

double precision dt(ncelet)
double precision rcodcl(nfabor,nvar,3)
double precision disale(3,nnod), xyzno0(3,nnod)
!< [arg]
!< [loc_var]
! Local variables

integer          ifac, iel, ii
integer          inod
integer          ilelt, nlelt

double precision delta, deltaa

integer, allocatable, dimension(:) :: lstelt
!< [loc_var]
!===============================================================================


!===============================================================================
! 1.  Initialization
!===============================================================================
!< [allocate_ale]
! Allocate a temporary array for boundary faces selection
allocate(lstelt(nfabor))
!< [allocate_ale]
!===============================================================================
! 2.  Assign boundary conditions to boundary faces here

!     One may use selection criteria to filter boundary case subsets
!       Loop on faces from a subset
!         Set the boundary condition for each face
!===============================================================================


! Calculation of displacement at current time step
!< [calcul]
deltaa = sin(3.141596d0*(ntcabs-1)/50.d0)
delta  = sin(3.141596d0*ntcabs/50.d0)
!< [calcul]

! Example: For boundary faces of color 4 assign a fixed velocity
!< [example_1]

call getfbr('4', nlelt, lstelt)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)
  ! Element adjacent a la face de bord
  iel = ifabor(ifac)

  ialtyb(ifac) = ivimpo
  rcodcl(ifac,iuma,1) = 0.d0
  rcodcl(ifac,ivma,1) = 0.d0
  rcodcl(ifac,iwma,1) = (delta-deltaa)/dt(iel)

enddo

!< [example_1]
! Example: For boundary faces of color 5 assign a fixed displacement on nodes
!< [example_2]

call getfbr('5', nlelt, lstelt)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

  do ii = ipnfbr(ifac), ipnfbr(ifac+1)-1
    inod = nodfbr(ii)
    if (impale(inod).eq.0) then
      disale(1,inod) = 0.d0
      disale(2,inod) = 0.d0
      disale(3,inod) = delta
      impale(inod) = 1
    endif
  enddo

enddo

!< [example_2]
! Example: For boundary faces of color 6 assign a sliding boundary
!< [example_3]

call getfbr('6', nlelt, lstelt)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

  ialtyb(ifac) = igliss

enddo

!< [example_3]
! Example: prescribe elsewhere a fixed boundary
!< [example_4]

call getfbr('not (4 or 5 or 6)', nlelt, lstelt)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

  ialtyb(ifac) = ibfixe

enddo

!< [example_4]
!--------
! Formats
!--------

!----
! End
!----
!< [deallocate_ale]
! Deallocate the temporary array
deallocate(lstelt)

!< [deallocate_ale]
return
end subroutine usalcl
