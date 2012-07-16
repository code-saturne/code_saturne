
! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2012 EDF S.A.
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

! Purpose:
! -------

! User subroutines for input of calculation parameters (Fortran commons).
!   These subroutines are called in all cases.

! If the Code_Saturne GUI is used, this file is not required (but may be
!   used to override parameters entered through the GUI, and to set
!   parameters not accessible through the GUI).

! Several routines are present in the file, each destined to defined
!   specific parameters.

! To modify the default value of parameters which do not appear in the
!   examples provided, code should be placed as follows:
!   - usipsu   for numerical and physical options
!   - usipes   for input-output related options

! As a convention, "specific physics" defers to the following modules only:
!   pulverized coal, gas combustion, electric arcs.

! In addition, specific routines are provided for the definition of some
!   "specific physics" options.
!   These routines are described at the end of this file and will be activated
!   when the corresponding option is selected in the usppmo routine.

!-------------------------------------------------------------------------------


!===============================================================================


subroutine cs_user_scal_drift
!============================

!===============================================================================
! Purpose:
! --------
! Parameters and physical properties for the Nerisson model for aerosol deposition
! based on the Zaichik model for particle transport.

! Note:
! =====
! The "full" Zaichik model is activated with the set of parameters
!               idrift = 2
!               itsttu = 1
!               itstde = 1
!
! In this case, the parameters iiso and ibrow must be prescribed
!
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use entsor
use cstphy
use numvar
use optcal
use mesh

!===============================================================================

implicit none

! Local variables



integer          ii, iel, iscal, ivar
integer          iutile

double precision alpha0, beta0, gamma0
double precision lg

!===============================================================================

!===============================================================================

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================

if (1.eq.1) return

!===============================================================================
! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END

!===============================================================================
! --- itsttu: Activation of the turbophorese modeling
!     ======
!               itsttu = 0: deactivation of the turbophorese modeling
!               itsttu = 1: activation of the turbophorese modeling
!===============================================================================

itsttu = 0

!===============================================================================
! --- itstde: Activation of the particle-trajectory deviation
!     ======
!               itstde = 0: deactivation of the particle-trajectory deviation
!               itstde = 1: activation of the particle-trajectory deviation
!
!CAUTION: UNSTABLE IN THE CURRENT VERSION
!===============================================================================

itstde = 0

!===============================================================================
! --- iiso: Taking into account the extradiagonal terms of the diffusion tensor
!     ===== (1: Yes; 0;: No)
!
!     CAUTION: UNSTABLE IN THE CURRENT VERSION
!===============================================================================

iiso  = 0

!===============================================================================
! --- ibrow: Taking into account the Brownian diffusion
!     ===== (1: Yes; 0;: No)
!===============================================================================

ibrow = 1

!===============================================================================
! --- ithphor: Taking into account the thermophoretic phenomenon
!     ===== (1: Yes; 0;: No)
!
!  CAUTION: TO BE IMPLEMENTED
!===============================================================================

ithphor  = 0

!===============================================================================
! --- ielectro: Taking into account the electrophoretic phenomenon
!     ======== (1: Yes; 0;: No)
!
! CAUTION: TO BE IMPLEMENTED
!===============================================================================

ielectro  = 0

!===============================================================================
! --- idepot: Choice of the wall deposition model
!     ======
!         idepot = 1: Deposition velocity  Vd = k1*Sc^k2 + Vsedim
!         idepot = 2: Lai & Nazaroff model
!         idepot = 3: Sedimentation velocity only
!         idepot = 4: O. Simonin model with y^+ term
!
! CAUTION: TO BE IMPLEMENTED
!===============================================================================

idepot = 4

! No turbulence model: idepot is forced to 3
if (iturb .eq. 0) idepot = 3

! "ndep" parameter for the Simonin model
! ndep = 1:  Wood model
! ndep = 2 : Prandtl model
ndep = 2

!===============================================================================
! Physical properties of the particles to be prescribed by the user
!
!
!==============================================================================

iscal = nscaus + 1
ivar  = isca(iscal)

nomvar(ipprtp(ivar)) = ""

!===============================================================================
! --- idrift: Choice for the model of scalar transport
!     ======
!               idrift = 0: passive tracer
!               idrift = 1: drift-flux model
!               idrift = 2: diffusion-inertia model
!===============================================================================

idrift(iscal) = 2

!-- diamm: Particle diameter (m)
!   =====
!  if polydispersity (ppolyd = 1) prescribe the mean particle diameter

diapart(iscal) = 100.d-6

!-- rhopart: Particle density (kg/m^3)
!   =======

rhopart(iscal) = 5180.d0

!-- sigmas: Turbulent Schmidt number (default = 1.d0):
!   ======

sigmas(iscal) = 1.d0


!--- ppolyd: Particle-size polydispersity (1: Yes; 0: No; default: 0)
!    ======

ppolyd(iscal) = 0

!-- if polydispersity (ppolyd = 1) prescribe the geometric standard deviation
!   denoted sigmag

if (ppolyd(iscal).gt.0) sigmag(iscal) = 1.0d0
!
if (ppolyd(iscal).gt.0) then
   !------------ Example of calculation of particle diameter
   diapart(iscal) = exp(log(diapart(iscal))-(4.d0-8.d0*(iscal-0.5d0)    &
        /nscadr)*log(sigmag(iscal)))
endif

!--- kpart: Particle thermal conductivity W/(m.K)
!    =====
!   (only useful if thermophorese is activated)

if (ithphor.gt.0) kpart(iscal) = 0.0

!--- qpart: Particle charge (C)
!    =====
!   (only useful if electrophorese is activated)

!      Charge des particules (Example of Park et al.)
! q = ce*p
! q = ((1+2*7.0d-8/diamp(iscal))**2.d0 + &
!     (2/(1+2*7.0d-8/diamp(iscal)))*((10-1)/(10+2)))*   &
!     pi*8.84d-12*(diamp(iscal))**2*500000.d0

if (ielectro.gt.0) qpart(iscal) = 0.0

! Electrical field

iutile = 0
if (iutile.eq.1) then

  if (ielectro.gt.0) then

    if (.not.allocated(elfield)) then
      allocate(elfield(ncel,3))
    endif

    do iel =1, ncel
      elfield(iel,1) = 0.d0
      elfield(iel,2) = 0.0d0
      elfield(iel,3) = 500000.0d0
    enddo

  endif

endif ! --- Test on 'iutile'

!----
! Formats
!----

!----
! End
!----

return
end subroutine
