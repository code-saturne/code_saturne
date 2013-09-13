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

!> \file cstphy.f90
!> \brief Module for physical constants

module cstphy

  !=============================================================================

  use paramx

  implicit none

  !=============================================================================

  !> \defgroup cstphy Module for physical constants

  !> \addtogroup cstphy
  !> \{

  !> Celsius to Kelvin: positive value of obsolute zero (= +273,15)
  double precision :: tkelvi
  parameter(tkelvi = 273.15d0)
  !> Kelvin to Celsius: negative value of obsolute zero (= -273,15)
  double precision :: tkelvn
  parameter(tkelvn = -273.15d0)

  !> Calories (1 cal = xcal2j J)
  double precision :: xcal2j
  parameter(xcal2j = 4.1855d0)

  !> Stephan Boltzmann (constant of)
  double precision :: stephn
  parameter(stephn = 5.6703d-8)
  !> Perfect gas constant for air (mixture)
  double precision :: rair
  parameter(rair = 287.d0)
  !> Gravity
  double precision, save :: gx, gy, gz
  integer, save :: icorio
  !> Rotation vector
  double precision, save :: omegax, omegay, omegaz
  ! TODO
  double precision, save :: irot(3,3), prot(3,3), qrot(3,3), rrot(3,3)

  !> Constantes physiques du fluide
  !> filling \ref xyzp0 indicator
  integer, save ::          ixyzp0
  !> reference density
  double precision, save :: ro0
  !> reference viscosity
  double precision, save :: viscl0
  !> reference total pressure
  double precision, save :: p0
  !> reference reduced pressure
  double precision, save :: pred0
  !> reference pressure position
  double precision, save :: xyzp0(3)
  !> reference temperature
  double precision, save :: t0
  !> reference specific heat
  double precision, save :: cp0
  !> molar mass of the perfect gas in \f$ kg/mol \f$ (if \ref ieos=1)
  double precision, save :: xmasmr

  !> Uniform thermodynamic pressure for the low-Mach algorithm
  !> Thermodynamic pressure for the current time step
  double precision, save :: pther
  !> Thermodynamic pressure for the previous time step
  double precision, save :: pthera

  !> \defgroup csttur Module for turbulence constants

  !> \addtogroup csttur
  !> \{

  !> constant of Karman (\f$ \kappa = 0.42 \f$)
  double precision, save :: xkappa
  !> constant of log law: \f$ \dfrac{1}{\kappa} \ln(y^+) + cstlog \f$
  !>  (\f$ cstlog = 5.2 \f$)
  double precision, save :: cstlog
  !> yplus limite (\f$ \dfrac{1}{\kappa} \f$ or \f$ 10.88 \f$ if \ref ideuch=2)
  double precision, save :: ypluli
  !> Werner and Wengle coefficient
  double precision, save :: apow
  !> Werner and Wengle coefficient
  double precision, save :: bpow
  !> Werner and Wengle coefficient
  double precision, save :: cpow
  !> Werner and Wengle coefficient
  double precision, save :: dpow
  !> \f$ C_\mu \f$ constant of \f$ k-\varepsilon\f$ model
  double precision, save :: cmu
  !> \f$ C_\mu^\frac{1}{4} \f$
  double precision, save :: cmu025
  !> constant of k-epsilon
  double precision, save :: ce1
  !> constant of k-epsilon
  double precision, save :: ce2
  !> Coefficient of interfacial coefficient in k-eps,
  !> used in Lagrange treatment
  double precision, save :: ce4
  !> constant of k-epsilon
  double precision, save :: sigmak
  !> constant for the Rij-epsilon EBRSM (1.15)
  double precision, save :: sigmae
  !> constant of standard Rij-epsilon (LRR)
  double precision, save :: crij1
  !> constant of standard Rij-epsilon (LRR)
  double precision, save :: crij2
  !> constant of standard Rij-epsilon (LRR)
  double precision, save :: crij3
  !> constant of standard Rij-epsilon (LRR)
  double precision, save :: crijp1
  !> constant of standard Rij-epsilon (LRR)
  double precision, save :: crijp2
  !> specific constant of SSG Rij-epsilon
  double precision, save :: cssge2
  !> specific constant of SSG Rij-epsilon
  double precision, save :: cssgs1
  !> specific constant of SSG Rij-epsilon
  double precision, save :: cssgs2
  !> specific constant of SSG Rij-epsilon
  double precision, save :: cssgr1
  !> specific constant of SSG Rij-epsilon
  double precision, save :: cssgr2
  !> specific constant of SSG Rij-epsilon
  double precision, save :: cssgr3
  !> specific constant of SSG Rij-epsilon
  double precision, save :: cssgr4
  !> specific constant of SSG Rij-epsilon
  double precision, save :: cssgr5
  !> constant of the Rij-epsilon EBRSM
  double precision, save :: cebms1
  !> constant of the Rij-epsilon EBRSM
  double precision, save :: cebms2
  double precision, save :: cebmr1, cebmr2, cebmr3, cebmr4, cebmr5, cebmr6
  !> constant of the Rij-epsilon EBRSM (0.21)
  double precision, save :: csrij
  !> constant of the Rij-epsilon EBRSM
  double precision, save :: cebme2
  !> constant of the Rij-epsilon EBRSM
  double precision, save :: cebmmu
  !> constant of the Rij-epsilon EBRSM
  double precision, save :: xcl
  !> constant in the expression of Ce1' for the Rij-epsilon EBRSM
  double precision, save :: xa1
  !> constant of the Rij-epsilon EBRSM
  double precision, save :: xct
  !> constant of the Rij-epsilon EBRSM
  double precision, save :: xceta

  !> specific constant of v2f "BL-v2k" (or phi-alpha)
  double precision, save :: cpale1
  !> specific constant of v2f "BL-v2k" (or phi-alpha)
  double precision, save :: cpale2
  !> specific constant of v2f "BL-v2k" (or phi-alpha)
  double precision, save :: cpale3
  !> specific constant of v2f "BL-v2k" (or phi-alpha)
  double precision, save :: cpale4
  !> specific constant of v2f "BL-v2k" (or phi-alpha)
  double precision, save :: cpalse
  !> specific constant of v2f "BL-v2k" (or phi-alpha)
  double precision, save :: cpalmu
  !> specific constant of v2f "BL-v2k" (or phi-alpha)
  double precision, save :: cpalc1
  !> specific constant of v2f "BL-v2k" (or phi-alpha)
  double precision, save :: cpalc2
  !> specific constant of v2f "BL-v2k" (or phi-alpha)
  double precision, save :: cpalct
  !> specific constant of v2f "BL-v2k" (or phi-alpha)
  double precision, save :: cpalcl
  !> specific constant of v2f "BL-v2k" (or phi-alpha)
  double precision, save :: cpalet

  !> specific constant of k-omega SST (sigma_k)
  double precision, save :: ckwsk1
  !> specific constant of k-omega SST (sigma_k)
  double precision, save :: ckwsk2
  !> specific constant of k-omega SST (sigma_w)
  double precision, save :: ckwsw1
  !> specific constant of k-omega SST (sigma_w)
  double precision, save :: ckwsw2
  !> specific constant of k-omega SST (beta)
  double precision, save :: ckwbt1
  !> specific constant of k-omega SST (beta)
  double precision, save :: ckwbt2
  !> specific constant of k-omega SST (gamma)
  double precision, save :: ckwgm1
  !> specific constant of k-omega SST (gamma)
  double precision, save :: ckwgm2
  !> specific constant of k-omega SST
  double precision, save :: ckwa1
  !> specific constant of k-omega SST
  double precision, save :: ckwc1
  !> specific constant of Spalart-Allmaras
  double precision, save :: csab1
  !> specific constant of Spalart-Allmaras
  double precision, save :: csab2
  !> specific constant of Spalart-Allmaras
  double precision, save :: csasig
  !> specific constant of Spalart-Allmaras
  double precision, save :: csav1
  !> specific constant of Spalart-Allmaras
  double precision, save :: csaw1
  !> specific constant of Spalart-Allmaras
  double precision, save :: csaw2
  !> specific constant of Spalart-Allmaras
  double precision, save :: csaw3
  !> constant of the Spalart-Shur rotation/curvature correction
  double precision, save :: cssr1
  !> constant of the Spalart-Shur rotation/curvature correction
  double precision, save :: cssr2
  !> constant of the Spalart-Shur rotation/curvature correction
  double precision, save :: cssr3
  !> constants of the Cazalbou rotation/curvature correction
  double precision, save :: ccaze2
  !> constants of the Cazalbou rotation/curvature correction
  double precision, save :: ccazsc
  !> constants of the Cazalbou rotation/curvature correction
  double precision, save :: ccaza
  !> constants of the Cazalbou rotation/curvature correction
  double precision, save :: ccazb
  !> constants of the Cazalbou rotation/curvature correction
  double precision, save :: ccazc
  !> constants of the Cazalbou rotation/curvature correction
  double precision, save :: ccazd
  !> scale of turbulent length
  double precision, save :: almax
  !> reference velocity
  double precision, save :: uref
  !> mixing length for the mixing length model
  double precision, save :: xlomlg
  !> constant used in the definition of LES filtering diameter:
  !> \f$ \delta = \text{xlesfl} . (\text{ales} . volume)^{\text{bles}} \f$
  double precision, save :: xlesfl
  !> constant used in the definition of LES filtering diameter \f$ \delta \f$
  double precision, save :: ales
  !> constant used in the definition of LES filtering diameter \f$ \delta \f$
  double precision, save :: bles
  !> Smagorinsky constant. In theory Smagorinsky constant is 0.18.
  !> For a planar canal plan, 0.065 value is rather taken.
  double precision, save :: csmago
  !> ratio between
  !> explicit and explicit filter width for a dynamic model
  double precision, save :: xlesfd
  !> Maximal desired Smagorinsky constatnt (e.g. 10 times \c csmago)
  double precision, save :: smagmx
  !> Van Driest constant in \f$ (1-\exp^{(-y^+/cdries}) \f$.
  !> (the Van Driest damping is activated  when \c idries=1).
  double precision, save :: cdries
  !> minimal control volume
  double precision, save :: volmin
  !> maximal control volume
  double precision, save :: volmax
  !> total domain volume
  double precision, save :: voltot
  !> constant of the v2f "\f$ \phi \f$-model" (\f$ \overline{f} \f$-bar)
  double precision, save :: cv2fa1
  !> constant of the v2f "\f$ \phi \f$-model" (\f$ \overline{f} \f$-bar)
  double precision, save :: cv2fe2
  !> constant of the v2f "\f$ \phi \f$-model" (\f$ \overline{f} \f$-bar)
  double precision, save :: cv2fmu
  !> constant of the v2f "\f$ \phi \f$-model" (\f$ \overline{f} \f$-bar)
  double precision, save :: cv2fc1
  !> constant of the v2f "\f$ \phi \f$-model" (\f$ \overline{f} \f$-bar)
  double precision, save :: cv2fc2
  !> constant of the v2f "\f$ \phi \f$-model" (\f$ \overline{f} \f$-bar)
  double precision, save :: cv2fct
  !> constant of the v2f "\f$ \phi \f$-model" (\f$ \overline{f} \f$-bar)
  double precision, save :: cv2fcl
  !> constant of the v2f "\f$ \phi \f$-model" (\f$ \overline{f} \f$-bar)
  double precision, save :: cv2fet
  !> constant of the WALE LES method
  double precision, save :: cwale
  !> coefficient of turbulent AFM flow model
  double precision, save :: xiafm
  !> coefficient of turbulent AFM flow model
  double precision, save :: etaafm
  !> coefficient of turbulent DFM flow model
  double precision, save :: c1trit
  !> coefficient of turbulent DFM flow model
  double precision, save :: c2trit
  !> coefficient of turbulent DFM flow model
  double precision, save :: c3trit
  !> coefficient of turbulent DFM flow model
  double precision, save :: c4trit
  !> constant of GGDH and AFM on the thermal scalar
  double precision, save :: cthafm
  !> constant of GGDH and AFM on the thermal scalar
  double precision, save :: cthdfm

  !> \}

  !> \}

  !=============================================================================

end module cstphy
