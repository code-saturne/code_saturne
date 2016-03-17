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

!> \file numvar.f90
!> \brief Module for variable numbering

module numvar

  !=============================================================================

  use paramx

  implicit none

  !=============================================================================

  !> \defgroup numvar Module for variable numbering

  !> \addtogroup numvar
  !> \{

  !----------------------------------------------------------------------------
  ! Main variables
  !----------------------------------------------------------------------------

  !> \defgroup main_variables Main variables
  !> \brief Main variable field indices (previously stored in rtp, rtpa).

  !> \addtogroup main_variables
  !> \{

  !> \anchor ipr
  !> pressure
  integer, save :: ipr

  !> \anchor iu
  !> velocity component \f$ u_x \f$
  integer, save :: iu

  !> \anchor iv
  !> velocity component \f$ u_y \f$
  integer, save :: iv

  !> \anchor iw
  !> velocity component \f$ u_z \f$
  integer, save :: iw

  !> \anchor ivoidf
  !> void fraction for cavitation modelling
  integer, save :: ivoidf

  !> \anchor ik
  !> turbulent kinetic energy \f$ k \f$
  integer, save :: ik

  !> \anchor iep
  !> turbulent dissipation \f$ \varepsilon \f$
  integer, save :: iep

  !> \anchor ir11
  !> Reynolds stress component \f$ R_{xx} \f$
  integer, save :: ir11

  !> \anchor ir22
  !> Reynolds stress component \f$ R_{yy} \f$
  integer, save :: ir22

  !> \anchor ir33
  !> Reynolds stress component \f$ R_{zz} \f$
  integer, save :: ir33

  !> \anchor ir12
  !> Reynolds stress component \f$ R_{xy} \f$
  integer, save :: ir12

  !> \anchor ir23
  !> Reynolds stress component \f$ R_{yz} \f$
  integer, save :: ir23

  !> \anchor ir13
  !> Reynolds stress component \f$ R_{zz} \f$
  integer, save :: ir13

  !> \anchor irij
  !> Reynolds stress tenso \f$ R_{ij} \f$
  integer, save :: irij

  !> \anchor iphi
  !> variable \f$ \phi \f$ of the \f$ \phi-f_b \f$ model
  integer, save :: iphi

  !> \anchor ifb
  !> variable \f$ f_b \f$ of the \f$ \phi-f_b \f$ model
  integer, save :: ifb

  !> \anchor ial
  !> variable \f$ \alpha \f$ of the \f$ Bl-v^2-k \f$ model
  integer, save :: ial

  !> \anchor iomg
  !> variable \f$ \omega \f$ of the \f$ k-\omega \f$ SST
  integer, save :: iomg

  !> \anchor inusa
  !> variable \f$ \widetilde{\nu}_T \f$ of the Spalart Allmaras
  integer, save :: inusa

  !> \anchor isca
  !> isca(i) is the index of the scalar i
  integer, save :: isca(nscamx)

  !> \anchor iscapp
  !> iscapp(i) is the index of the specific physics scalar i
  integer, save :: iscapp(nscamx)

  !> \anchor nscaus
  !> number of user scalars
  integer, save :: nscaus

  !> \anchor nscapp
  !> number of specific physics scalars
  integer, save :: nscapp

  !> \anchor nscasp
  !> number of species scalars
  integer, save :: nscasp

  !> \anchor iuma
  !> mesh velocity component \f$ w_x \f$
  integer, save :: iuma

  !> \anchor ivma
  !> mesh velocity component \f$ w_y \f$
  integer, save :: ivma

  !> \anchor iwma
  !> mesh velocity component \f$ w_z \f$
  integer, save :: iwma

  !> \}

  !----------------------------------------------------------------------------
  ! Physical properties
  !----------------------------------------------------------------------------

  !> \defgroup physical_prop Physical properties
  !> \brief Physical properties are stored in propce.
  !> See \ref cs_user_boundary_conditions for some examples.

  !> \addtogroup physical_prop
  !> \{

  !> pointer to cell properties (propce)
  integer, save :: ipproc(npromx)

  !> Density at the current time step
  integer, save :: irom

  !> Density at the previous time step
  integer, save :: iroma

  !> Density at the second previous time step
  integer, save :: iromaa

  !> dynamic molecular viscosity (in kg/(m.s))
  integer, save :: iviscl

  !> dynamic turbulent viscosity
  integer, save :: ivisct

  !> dynamic molecular viscosity (in kg/(m.s)) at the previous time-step
  integer, save :: ivisla

  !> dynamic turbulent viscosity at the previous time-step
  integer, save :: ivista

  !> specific heat \f$ C_p \f$ at the previous time-step
  integer, save :: icpa

  !> error estimator for Navier-Stokes
  integer, save :: iestim(nestmx)

  !> interior and boundary convective mass flux key ids of the variables
  integer, save :: kimasf, kbmasf

  !> constant diffusivity field id key for scalars
  integer, save :: kvisl0

  !> variable diffusivity field id key for scalars
  integer, save :: kivisl

  !> source terms at previous time step for 2nd order
  integer, save :: kstprv

  !> source terms at the current time step (used for limiters)
  integer, save :: kst

  !> convective mass flux of the variables at the previous time-step
  integer, save :: ifluaa(nvarmx)

  !> cell density field ids of the variables
  integer, save :: icrom

  !> boundary density field ids of the variables
  integer, save :: ibrom

  !> cell density at the second previous time step key id of the variables
  integer, save :: icroaa

  !> field ids of the cell porosity
  integer, save :: ipori, iporf

  !> dynamic constant of Smagorinsky
  integer, save :: ismago

  !> field ids of the anisotropic viscosity
  !> \remark turbulent or Darcy module anisotropic diffusion
  integer, save :: ivsten, ivstes

  !> Courant number
  integer, save :: icour

  !> Fourier number
  integer, save :: ifour

  !> Total pressure at cell centers
  !> \f$ P_{tot} = P^\star +\rho \vect{g} \cdot (\vect{x}-\vect{x}_0) \f$
  integer, save :: iprtot

  !> Mesh velocity viscosity for the ALE module
  !> \remark might be orthotropic
  integer, save :: ivisma(3)

  !> pointer for dilatation source terms
  integer, save :: iustdy(nscamx)

  !> pointer for global dilatation source terms
  integer, save :: itsrho

  !> pointer for thermal expansion coefficient
  integer, save :: ibeta

  !> pointer for deduced mass fraction in case of gas mix
  integer, save :: iddgas

  !> pointer for gas mix molar mass
  integer, save :: igmxml

  !> field id of the stresses at boundary  (if post-processed)
  integer, save :: iforbr

  !>  field id of \f$y^+\f$ at boundary (if post-processed)
  integer, save :: iyplbr

  !>  field id of temperature at boundary
  integer, save ::  itempb
  !> \}
  !----------------------------------------------------------------------------
  ! Numerical properties
  !----------------------------------------------------------------------------

  !> \defgroup numerical_prop Numerical properties

  !> \addtogroup numerical_prop
  !> \{

  !> Weighting for gradient calculation on variables
  integer, save :: kwgrec

  !> \}
  !----------------------------------------------------------------------------
  ! Mapping to field structures
  !----------------------------------------------------------------------------

  !> \defgroup field_map Mapping to field structures

  !> \addtogroup field_map
  !> \{

  !> Field id for variable i
  integer, save :: ivarfl(nvarmx)

  !> Field id for property i
  integer, save :: iprpfl(npromx)

  !> Field id for the dttens tensor
  integer, save :: idtten

  !> \}

  !=============================================================================

  !> \}

end module numvar
