!-------------------------------------------------------------------------------

!VERS

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

subroutine cfther &
!================

 ( nvar   ,                                                       &
   iccfth , imodif ,                                              &
   rtp    ,                                                       &
   sorti1 , sorti2 , gamagr , bval   , wbfb   )

!===============================================================================
! Purpose:
! -------

!    Define thermodynamic laws (especially for the compressible flow scheme).


! Introduction
! ============

! Avalable thermodynamic laws
! ===========================

!  1. Perfect gas (the molar mass 'xmasml' must be provided)
!  2. Perfect gas with non constant Gamma (example to be adapted)
!  3. Van Der Waals (not yet implemented)


! Implemented calculations
! ========================

! This user subroutine implements the computation of several quantities.
! Each calculation has to be explicitly implemented in the appropriate
! section below (already done for perfect gas).


! Selection of the quantity to return
! ===================================

! When calling the user subroutine, the integer 'iccfth' specifies which
! calculation has to be performed (and which quantity has to be returned).
! The values for 'iccfth' for each case are provided below.
! For some configurations, two systems of references are used for 'iccfth'
! (this is useful to make tests easier to implement in the calling
! subroutines): both systems are explained hereafter for information.

! First system:

!   the variables are referred to using an index i:
!     Variable  P  rho  T   e   h   s  'internal energy - CvT'
!        Index  1   2   3   4   5   6              7

!   iccfth is as follows, depending on which quantity needs to be computed:
!     - compute all variables at cell centers from variable i
!                                              and variable j (i<j):
!               => iccfth = 10*i+j
!     - compute all variables at boundary faces from variable i
!                                                and variable j (i<j):
!               => iccfth = 10*i+j+900

! Second system:

!   the variables are referred to using a different index i:
!     Variable  P  rho  T  e  s
!        Index  2   3   5  7  13

!   iccfth is as follows, depending on which quantity needs to be computed:
!     - compute all variables at cell centers from variable i
!                                              and variable j (i<j):
!               => iccfth = i*j*10000
!     - compute all variables at boundary faces from variable i
!                                                and variable j (i<j):
!               => iccfth = i*j*10000+900

! Other quantities:

!   the variables are referred to using the index of the first system.
!   iccfth is defined as follows:
!     - compute variable i at cell centers (for s and 'internal energy-CvT')
!               => iccfth = i
!                                   \partial(variable i)|
!     - compute partial derivative  --------------------|
!                                   \partial(variable j)|variable k
!               => iccfth = 100*i+10*j+k
!     - compute boundary conditions, resp. symmetry, wall, outlet, inlet, inlet:
!               => iccfth = 90, 91, 93, 94, 95


! Values of iccfth
! ================

! To summarize, the values for iccfth are as follows:

!   Values at the cell centers:

!   -> set calculation options (cst/variable cp)   : iccfth = -1
!   -> set default initialization                  : iccfth =  0
!   -> calculate gamma                             : iccfth =  1
!   -> verification of the density                 : iccfth = -2
!   -> verification of the energy                  : iccfth = -4
!   -> calculation of temperature and energy
!                     from pressure and density    : iccfth =  12 or  60000
!   -> calculation of density and energy
!                     from pressure and temperature: iccfth =  13 or 100000
!   -> calculation of density and temperature
!                     from pressure and energy     : iccfth =  14 or 140000
!   -> calculation of pressure and energy
!                     from density and temperature : iccfth =  23 or 150000
!   -> calculation of pressure and temperature
!                     from density and energy      : iccfth =  24 or 210000
!
!                      2    dP |
!   -> calculation of c  = ----|                   : iccfth = 126
!                          drho|s
!
!                            dP|
!   -> calculation of beta = --|                   : iccfth = 162
!                            ds|rho
!
!                          de|
!   -> calculation of Cv = --|                     : iccfth = 432
!                          dT|rho
!
!   -> calculation of entropie                     : iccfth =   6
!
!
!   Values at the boundary faces
!
!   -> calculation of the boundary conditions:
!     - symmetry                                   : iccfth =  90
!     - wall                                       : iccfth =  91
!     - generalized outlet                         : iccfth =  93
!     - inlet (not yet implemented)                : iccfth =  94
!     - inlet                                      : iccfth =  95
!
!   -> calculation of the variables at the faces for boundary conditions:
!     - temperature and energy
!         from pressure and density                : iccfth = 912 ou  60900
!     - density and energy
!         from pressure and temperature            : iccfth = 913 ou 100900
!     - density and temperature
!         from pressure and energy                 : iccfth = 914 ou 140900
!     - pressure and energy
!         from density and temperature             : iccfth = 923 ou 150900
!     - pressure and temperature
!         from density and energy                  : iccfth = 924 ou 210900


!   Values at the cell centers and at the boundary faces

!   -> calculation of 'internal energy - Cv.T'     : iccfth =   7

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
!    nom           !type!mode !                   role                         !
!__________________!____!_____!________________________________________________!
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! rtp              ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current time step)                        !
! sorti1,2(*)      ! ra ! --> ! output variable (unused if iccfth.lt.0)        !
! gamagr(*)        ! ra ! --> ! equivalent "gamma" constant of the gas         !
!                  !    !     !   (unused if iccfth.lt.0)                      !
!                  !    !     !   (first value only used for perfect gas)      !
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use optcal
use cstphy
use cstnum
use parall
use pointe
use entsor
use ppppar
use ppthch
use ppincl
use mesh
use field

!===============================================================================

implicit none

! Arguments

integer          nvar
integer          iccfth   , imodif

double precision rtp(ncelet,*)
double precision wbfb(nfabor), bval(nfabor,nvar)

double precision sorti1(*), sorti2(*), gamagr(*)

! Local variables

integer          ifac0
integer          ierr
integer          iel    , ifac
integer          itk    , ien    , niter, nitermax
double precision gamagp , xmasml , enint, norm, cosalp
double precision xmach  , xmachi , xmache , dxmach
double precision dir(3)
double precision roi, ro1, pri, ei, uni, un1, y, uns, bc, pinf, ptot, eps, res, rhotot
double precision ci, c1, mi, a, b, sigma1, utxi, utyi, utzi, bMach, old_pstat, pstat

double precision, dimension(:), pointer :: crom, brom

!===============================================================================

!===============================================================================
! 0. Initialization.
!===============================================================================

! Error indicator (stop if non zero)
ierr   = 0

! Rank of the variables in their associated arrays
if (iccfth.ge.0.or.iccfth.le.-2) then
  call field_get_val_s(icrom, crom)
  call field_get_val_s(ibrom, brom)
  itk = isca(itempk)
  ien = isca(ienerg)
endif

! For calculation of values at the cell centers,
!   ifac0 > indicates that the array rtp must be modified
! For calculation of values at the cell faces,
!   ifac0 is the number of the current face
ifac0 = imodif

!===============================================================================
! 1. Cp and Cv variable or not (depends on the equation of state)
!===============================================================================

! Warning: once the thermodynamic law has been chosen,
! =======  the remainder of the user subroutine must be modified

if (iccfth.eq.-1) then
  if (ieos.eq.1) then

! --- Calculation options: constant Cp and Cv (perfect gas)

  ! The value for the isobaric specific heat Cp0 must be provided in
  !   the user subroutine ''usipph''. The value for the isochoric
  !   specific heat Cv0 is calculated in a subsequent section (from Cp0)

    icp = 0
    icv = 0

  endif

  return

endif

!===============================================================================
! 2. Perfect gas
!===============================================================================

if (ieos.eq.1) then

! --- Molar mass of the gas (kg/mol)
    ! The value is put in xmasml

  if (iccfth.ge.0) then
    xmasml = xmasmr
  endif

!===============================================================================
! 2.1. Default laws
!===============================================================================

! --- Calculation of the constant gamagp

  if (iccfth.gt.0) then

    ! Gamagp is supposed to be superior or equal to 1.
    ! It is computed at each call, even if this may seem costly,
    !   to be coherent with the "constant gamma" case for which this
    !   constant is not saved. A ''save'' instruction and a test would
    !   be sufficient to avoid computing gamagp at each call if necessary.

    gamagp = 1.d0 + rr/(xmasml*cp0-rr)

    if (gamagp.lt.1.d0) then
      write(nfecra,1010) gamagp
      call csexit (1)
    endif

    ! Gamma is returned if required

    if (iccfth.eq.1) then
      gamagr(1) = gamagp
    endif

  endif


! --- Default initializations (before uscfxi)

!     T0 is positive (this assumption has been checked in
!       the user programme 'verini')

  if (iccfth.eq.0) then

    cv0 = cp0 - rr/xmasml

    if ( isuite .eq. 0 ) then
      do iel = 1, ncel
        crom(iel) = p0*xmasml/(rr*t0)
        rtp(iel,ien) = cv0*t0
      enddo
    endif


! --- Verification of the density

  elseif (iccfth.eq.-2) then

    ! If the density is lower or equal to zero: clipping, write and stop.
    !   Indeed, if this is the case, the thermodynamic computations will
    !   most probably fail.
    ! This call is done at the end of the density calculation (after
    !   a classical clipping and before parallel communications).

    ierr = 0
    do iel = 1, ncel
      if (rtp(iel,ipr).le.0.d0) then
        rtp(iel,ipr) = epzero
        ierr = ierr + 1
      endif
    enddo
    if (irangp.ge.0) then
      call parcpt (ierr)
    endif
    if (ierr.gt.0) then
      ntmabs = ntcabs
      write(nfecra,8000)ierr, epzero
    endif


! --- Verification of the energy

  elseif (iccfth.eq.-4) then

    ! If the total energy <= zero: clipping, write and stop
    !   Indeed, if this is the case, the thermodynamic computations will
    !   most probably fail.

    ierr = 0
    do iel = 1, ncel
      enint = rtp(iel,ien)                                     &
               - 0.5d0*( rtp(iel,iu)**2                        &
                       + rtp(iel,iv)**2                        &
                       + rtp(iel,iw)**2 )
      if (enint.le.0.d0) then
        rtp(iel,ien) = epzero                                  &
               + 0.5d0*( rtp(iel,iu)**2                        &
                       + rtp(iel,iv)**2                        &
                       + rtp(iel,iw)**2 )
        ierr = ierr + 1
      endif
    enddo
    if (irangp.ge.0) then
      call parcpt (ierr)
    endif
    if (ierr.gt.0) then
      ntmabs = ntcabs
      write(nfecra,8100)ierr, epzero
    endif


! --- Calculation of temperature and energy from pressure and density

  elseif (iccfth.eq.12.or.iccfth.eq.60000) then

    ! Verification of the values of the density
    ierr = 0
    do iel = 1, ncel
      if (crom(iel).le.0.d0) then
        write(nfecra,3010)crom(iel),iel
      endif
    enddo
    ! Stop if a negative value is detected (since the density has been
    ! provided by the user, one potential cause is a wrong user
    ! initialization)
    if (ierr.eq.1) then
      call csexit (1)
    endif

    do iel = 1, ncel
      ! Temperature
      sorti1(iel) = xmasml*rtp(iel,ipr)/(rr*crom(iel))
      ! Total energy
      sorti2(iel) = cv0*sorti1(iel)                        &
           + 0.5d0*( rtp(iel,iu)**2 + rtp(iel,iv)**2 + rtp(iel,iw)**2 )
    enddo

    ! Transfer to the array rtp
    if (imodif.gt.0) then
      do iel = 1, ncel
        rtp(iel,itk) = sorti1(iel)
        rtp(iel,ien) = sorti2(iel)
      enddo
    endif


! --- Calculation of density and energy from pressure and temperature:

  elseif (iccfth.eq.13.or.iccfth.eq.100000) then

    ! Verification of the values of the temperature
    ierr = 0
    do iel = 1, ncel
      if (rtp(iel,itk).le.0.d0) then
        write(nfecra,2010)rtp(iel,itk),iel
      endif
    enddo
    ! Stop if a negative value is detected (since the temperature has been
    ! provided by the user, one potential cause is a wrong user
    ! initialization: a value not provided in Kelvin for example)
    if (ierr.eq.1) then
      call csexit (1)
    endif

    do iel = 1, ncel
      ! Density
      sorti1(iel) = xmasml*rtp(iel,ipr)/(rr*rtp(iel,itk))
      ! Total energy
      sorti2(iel) = cv0*rtp(iel,itk)                    &
           + 0.5d0*( rtp(iel,iu)**2 + rtp(iel,iv)**2 + rtp(iel,iw)**2 )
    enddo

    ! Transfer to the array rtp
    if (imodif.gt.0) then
      do iel = 1, ncel
        crom(iel) = sorti1(iel)
        rtp(iel,ien) = sorti2(iel)
      enddo
    endif


! --- Calculation of density and temperature from pressure and energy

  elseif (iccfth.eq.14.or.iccfth.eq.140000) then

    do iel = 1, ncel
      ! Internal energy (to avoid the need to divide by the temperature
      ! to compute density)
      enint = rtp(iel,ien)                                     &
               - 0.5d0*( rtp(iel,iu)**2                        &
                       + rtp(iel,iv)**2                        &
                       + rtp(iel,iw)**2 )
      ! Density
      sorti1(iel) = rtp(iel,ipr) / ( (gamagp-1.d0) * enint )
      ! Temperature
      sorti2(iel) = xmasml * (gamagp-1.d0) * enint / rr
    enddo

    ! Transfer to the array rtp
    if (imodif.gt.0) then
      do iel = 1, ncel
        crom(iel) = sorti1(iel)
        rtp(iel,itk) = sorti2(iel)
      enddo
    endif


! --- Calculation of pressure and energy from density and temperature

  elseif (iccfth.eq.23.or.iccfth.eq.150000) then

    do iel = 1, ncel
      ! Pressure
      sorti1(iel) = crom(iel)*rtp(iel,itk)*rr/xmasml
      ! Total energy
      sorti2(iel) = cv0*rtp(iel,itk)                    &
           + 0.5d0*( rtp(iel,iu)**2 + rtp(iel,iv)**2 + rtp(iel,iw)**2 )
    enddo

    ! Transfer to the array rtp
    if (imodif.gt.0) then
      do iel = 1, ncel
        rtp(iel,ipr) = sorti1(iel)
        rtp(iel,ien) = sorti2(iel)
      enddo
    endif


! --- Calculation of pressure and temperature from density and energy

  elseif (iccfth.eq.24.or.iccfth.eq.210000) then

    do iel = 1, ncel
      ! Internal energy (to avoid the need to divide by the temperature
      ! to compute density)
      enint = rtp(iel,ien)                                     &
               - 0.5d0*( rtp(iel,iu)**2                        &
                       + rtp(iel,iv)**2                        &
                       + rtp(iel,iw)**2 )
      ! Pressure
      sorti1(iel) = (gamagp-1.d0) * crom(iel) * enint
      ! Temperature
      sorti2(iel) = xmasml * (gamagp-1.d0) * enint / rr
    enddo

    ! Transfer to the array rtp
    if (imodif.gt.0) then
      do iel = 1, ncel
        rtp(iel,ipr) = sorti1(iel)
        rtp(iel,itk) = sorti2(iel)
      enddo
    endif


!                     2                            2         P
! --- Calculation of c from pressure and density: c = gamma*---
!                                                           rho

  elseif (iccfth.eq.126) then

    ! Verification of the values of the density
    !   This test can be discarded to reduce the CPU time (if
    !     density is <= 0, the calculation will simply fail)
    !   It is discarded here with .false.
    if (.false.) then
      ierr = 0
      do iel = 1, ncel
        if (crom(iel).le.0.d0) then
          write(nfecra,4010)crom(iel),iel
        endif
      enddo
      if (ierr.eq.1) then
        call csexit (1)
      endif
    endif

    do iel = 1, ncel
      sorti1(iel) = gamagp * rtp(iel,ipr) / crom(iel)
    enddo


!                                                              gamma
! --- Calculation of beta from pressure and density: beta = rho

  elseif (iccfth.eq.162) then

    ! Verification of the values of the density
    !   This test can be discarded to reduce the CPU time (if
    !     density is <= 0, the calculation will simply fail)
    !   It is discarded here with .false.
    if (.false.) then
      ierr = 0
      do iel = 1, ncel
        if (crom(iel).lt.0.d0) then
          write(nfecra,4020)crom(iel),iel
        endif
      enddo
      if (ierr.eq.1) then
        call csexit (1)
      endif
    endif

    do iel = 1, ncel
      sorti1(iel) = crom(iel)**gamagp
    enddo


! --- Calculation of the isochoric specific heat

    ! It is a constant: nothing to do


!                                                                  P
! --- Calculation of the entropy from pressure and density: s = --------
!                                                                  gamma
!                                                               rho

  elseif (iccfth.eq.6) then

    ! Verification of the values of the density
    !   This test can be discarded to reduce the CPU time (if
    !     density is <= 0, the calculation will simply fail)
    ierr = 0
    do iel = 1, ncel
      if (crom(iel).le.0.d0) then
        write(nfecra,4030)crom(iel),iel
      endif
    enddo
    if (ierr.eq.1) then
      call csexit (1)
    endif

    do iel = 1, ncel
      sorti1(iel) = rtp(iel,ipr) / (crom(iel)**gamagp)
    enddo


! --- Calculation of 'internal energy - Cv.T'

  elseif (iccfth.eq.7) then

    ! It is zero for a perfect gas

    !   At the cell centers
    do iel = 1, ncel
      sorti1(iel) = 0.d0
    enddo

    !   On the boundary faces
    do ifac = 1, nfabor
      sorti2(ifac) = 0.d0
    enddo


! --- Calculation of the boundary conditions on the face ifac = ifac0

!  -- Wall

  elseif (iccfth.eq.91) then

    ifac = ifac0
    iel  = ifabor(ifac)

    ! Calculation of the Mach number at the boundary face, using the
    !   cell center velocity projected on the vector normal to the boundary
    xmach =                                                    &
         ( rtp(iel,iu)*surfbo(1,ifac)                          &
         + rtp(iel,iv)*surfbo(2,ifac)                          &
         + rtp(iel,iw)*surfbo(3,ifac) ) / surfbn(ifac)         &
         / sqrt( gamagp * rtp(iel,ipr) / crom(iel) )

    ! Pressure

    !   A Neumann boundary condition is used. This does not allow to use
    !     the Rusanov scheme, but some stabilization effect is expected.
    !     A test based on the value of coefb at the previous time step
    !     is implemented to avoid oscillating between a rarefaction
    !     situation and a shock configuration from one time step to the
    !     next.

    !   Rarefaction !FIXME with the new cofaf cofbf
    if (xmach.lt.0.d0.and.wbfb(ifac).le.1.d0) then

      if (xmach.gt.2.d0/(1.d0-gamagp)) then
        wbfb(ifac) = (1.d0 + (gamagp-1.d0)/2.d0 * xmach)    &
                     ** (2.d0*gamagp/(gamagp-1.d0))
      else
        ! In case the rarefaction is too strong, a zero Dirichlet value
        !   is used for pressure (the value of wbfb is used here as an
        !   indicator)
        wbfb(ifac) = rinfin
      endif

      !  Shock
    elseif (xmach.gt.0.d0.and.wbfb(ifac).ge.1.d0) then

      wbfb(ifac) = 1.d0 + gamagp*xmach                      &
                   *( (gamagp+1.d0)/4.d0*xmach                           &
                   + sqrt(1.d0 + (gamagp+1.d0)**2/16.d0*xmach**2) )

      !  Oscillation between rarefaction and shock or zero Mach number
    else
      wbfb(ifac) = 1.d0
    endif


!  -- Symmetry
    ! A zero flux condition (homogeneous Neumann condition) is
    !   prescribed by default.
    ! No user input required

  elseif (iccfth.eq.90) then

    ifac = ifac0
    iel  = ifabor(ifac)

!  -- Subsonic inlet with prescribed mass and enthalpy flow rates
    ! The quantities prescribed are rho*u and rho*u*h

    ! The subsonic nature of the inlet is postulated.

    ! This section remains to be implemented: stop for the moment

    ! One may proceed as follows:
    !   Pressure computed with a Newton method
    !   Velocity and density computed from pressure
    !   Total energy computed from enthalpy
    !   (written on paper, to be implemented: contact the user support)

  elseif (iccfth.eq.94) then

    ifac = ifac0
    iel  = ifabor(ifac)

    write(nfecra,7000)

    call csexit (1)
    !==========


!  -- Generalized outlet
  elseif (iccfth.eq.93) then

    ifac = ifac0
    iel  = ifabor(ifac)

    pinf = bval(ifac,ipr)
    pri  = rtp(iel,ipr)
    roi  = crom(iel)

    ci = sqrt(gamagp * pri / roi)
    uni = ( rtp(iel,iu) * surfbo(1,ifac)                             &
          + rtp(iel,iv) * surfbo(2,ifac)                             &
          + rtp(iel,iw) * surfbo(3,ifac) ) / surfbn(ifac)
    ei   = rtp(iel,isca(ienerg)) - 0.5d0 * uni**2

    ! Rarefaction case
    if (pinf.le.pri) then

      ! Computation of the velocity in state 1 using Riemann invariants of the 1-rarefaction
      a = 2 * ci / (gamagp - 1.d0) * (1.d0 - (pinf / pri) ** ((gamagp - 1.d0) / (2.d0 * gamagp)))
      un1 = uni + a

      ! Computation of the density in state 1 using Rieman invariants of the 1-rarefaction
      ro1 = roi * (pinf/pri)**(1.d0 / gamagp)

      ! Subsonic inlet - state 2 should be imposed but too few information is available to compute it
      ! for want of anything better, state 1 is imposed
      if (un1.lt.0.d0) then

        ! Density
        brom(ifac) = ro1
        ! Velocity
        bval(ifac,iu) = rtp(iel,iu) + a * surfbo(1,ifac) / surfbn(ifac)
        bval(ifac,iv) = rtp(iel,iv) + a * surfbo(2,ifac) / surfbn(ifac)
        bval(ifac,iw) = rtp(iel,iw) + a * surfbo(3,ifac) / surfbn(ifac)
        ! Total energy
        bval(ifac,ien) = pinf / ((gamagp - 1.d0) * ro1)                         &
                         + 0.5d0 * (bval(ifac,iu)**2                            &
                                  + bval(ifac,iv)**2                            &
                                  + bval(ifac,iw)**2)

      ! Outlet
      else

        ! Computation of the sound speed in state 1
        c1 = sqrt(gamagp * pinf / ro1)

        ! Subsonic outlet - state 1 is imposed
        if ((un1-c1).lt.0.d0) then

          ! Density
          brom(ifac) = ro1
          ! Velocity
          bval(ifac,iu) = rtp(iel,iu) + a * surfbo(1,ifac) / surfbn(ifac)
          bval(ifac,iv) = rtp(iel,iv) + a * surfbo(2,ifac) / surfbn(ifac)
          bval(ifac,iw) = rtp(iel,iw) + a * surfbo(3,ifac) / surfbn(ifac)
          ! Total energy
          bval(ifac,ien) = pinf / ((gamagp - 1.d0) * ro1)                         &
                           + 0.5d0 * (bval(ifac,iu)**2                            &
                                    + bval(ifac,iv)**2                            &
                                    + bval(ifac,iw)**2)

        ! Sonic outlet
        else if ((uni-ci).lt.0.d0) then

          ! Mach number in the domain
          mi = uni / ci

          b = (gamagp - 1.d0) / (gamagp + 1.d0) * (mi + 2.d0 / (gamagp - 1))

          ! Sonic state pressure
          bval(ifac,ipr) = pri * b ** (2.d0 * gamagp / (gamagp - 1.d0))
          ! Sonic state density
          brom(ifac) = roi * b ** (2.d0 / (gamagp - 1.d0))
          ! Sonic state velocity
          uns = b * ci
          bval(ifac,iu) = uns * surfbo(1,ifac) / surfbn(ifac)
          bval(ifac,iv) = uns * surfbo(2,ifac) / surfbn(ifac)
          bval(ifac,iw) = uns * surfbo(3,ifac) / surfbn(ifac)
          ! Sonic state energy
          bval(ifac,isca(ienerg)) = bval(ifac,ipr) / ((gamagp - 1.d0) * brom(ifac)) + 0.5d0 * uns**2

        ! Supersonic outlet
        else

          ! pb = pri
          bval(ifac,ipr) = pri
          ! ub = uni
          bval(ifac,iu) = rtp(iel,iu)
          bval(ifac,iv) = rtp(iel,iv)
          bval(ifac,iw) = rtp(iel,iw)
          ! rob = roi
          brom(ifac) = roi
          ! eb = ei
          bval(ifac,isca(ienerg)) = rtp(iel,isca(ienerg))

        endif

      endif

    ! Shock case
    else

      ! Computation of the density in state 1 with Rankine-Hugoniot relations
      ro1 = roi * ((gamagp - 1.d0) * pri   + (gamagp + 1.d0) * pinf) /&
                  ((gamagp - 1.d0) * pinf + (gamagp + 1.d0) * pri)

      ! Computation of the velocity in state 1 with Rankine-Hugoniot relations
      ! un1 = un2
      a = sqrt( (pinf - pri) * (1.d0/roi - 1.d0/ro1) )
      un1 = uni - a

      ! Subsonic inlet - state 2 should be imposed but too few information is available to compute it
      ! for want of anything better, state 1 is imposed
      if (un1.le.0d0) then

        ! Density
        brom(ifac) = ro1
        ! Velocity
        bval(ifac,iu) = rtp(iel,iu) - a * surfbo(1,ifac) / surfbn(ifac)
        bval(ifac,iv) = rtp(iel,iv) - a * surfbo(2,ifac) / surfbn(ifac)
        bval(ifac,iw) = rtp(iel,iw) - a * surfbo(3,ifac) / surfbn(ifac)
        ! Total energy
        bval(ifac,ien) = pinf / ((gamagp-1.d0) * brom(ifac))           &
                         + 0.5d0 * (bval(ifac,iu)**2                   &
                                  + bval(ifac,iv)**2                   &
                                  + bval(ifac,iw)**2)

      ! Outlet
      else

        ! Computation of the shock velocity
        sigma1 = (roi * uni - ro1 * un1) / (roi - ro1)

        ! Subsonic outlet - state 1 is imposed
        if (sigma1.le.0.d0) then

          ! Density
          brom(ifac) = ro1
          ! Velocity
          bval(ifac,iu) = rtp(iel,iu) - a * surfbo(1,ifac) / surfbn(ifac)
          bval(ifac,iv) = rtp(iel,iv) - a * surfbo(2,ifac) / surfbn(ifac)
          bval(ifac,iw) = rtp(iel,iw) - a * surfbo(3,ifac) / surfbn(ifac)
          ! Total energy
          bval(ifac,ien) = pinf / ((gamagp-1.d0) * brom(ifac))           &
                           + 0.5d0 * (bval(ifac,iu)**2                   &
                                    + bval(ifac,iv)**2                   &
                                    + bval(ifac,iw)**2)

        ! Supersonic outlet
        else

          ! pb = pri
          bval(ifac,ipr) = pri
          ! unb = uni
          bval(ifac,iu) = rtp(iel,iu)
          bval(ifac,iv) = rtp(iel,iv)
          bval(ifac,iw) = rtp(iel,iw)
          ! rob = roi
          brom(ifac) = roi
          ! eb = ei
          bval(ifac,isca(ienerg)) = rtp(iel,isca(ienerg))

        endif ! test on shock speed sign

      endif ! test on state 1 velocity sign

    endif ! test on pinf-pri sign

  !  -- Boundary condition at prescribed total pressure and enthalpy
  elseif (iccfth.eq.95) then

    niter = 0

    ifac = ifac0
    iel  = ifabor(ifac)

    roi  = crom(iel)
    pri  = rtp(iel,ipr)

    ! Normalize the direction vector given by the user
    norm = sqrt(bval(ifac,iu)**2 + bval(ifac,iv)**2 + bval(ifac,iw)**2)
    if (norm.lt.epzero) then
      write(nfecra,9010)ifac
      call csexit (1)
    endif

    dir(1) = bval(ifac,iu) / norm
    dir(2) = bval(ifac,iv) / norm
    dir(3) = bval(ifac,iw) / norm

    ! Angle between the imposed direction and the inlet normal
    cosalp = ( dir(1)*surfbo(1,ifac) + dir(2)*surfbo(2,ifac)  &
             + dir(3)*surfbo(3,ifac) ) /surfbn(ifac)

    ! If direction vector is outward, warn the user
    if (cosalp.gt.epzero) then
      write(nfecra,9020)ifac
    endif

    ! Computation of the sound speed inside the domain
    ci = sqrt(gamagp * pri / roi)

    uni = ( rtp(iel,iu) * surfbo(1,ifac)                             &
          + rtp(iel,iv) * surfbo(2,ifac)                             &
          + rtp(iel,iw) * surfbo(3,ifac) ) / surfbn(ifac)

    bMach = uni / ci

    utxi = rtp(iel,iu) - uni * surfbo(1,ifac) * surfbn(ifac)
    utyi = rtp(iel,iv) - uni * surfbo(2,ifac) * surfbn(ifac)
    utzi = rtp(iel,iw) - uni * surfbo(3,ifac) * surfbn(ifac)

    ei   = rtp(iel,isca(ienerg)) - 0.5d0 * ( rtp(iel,iu)**2 + rtp(iel,iv)**2 + rtp(iel,iw)**2 )

    ptot = bval(ifac,ipr)
    rhotot = gamagp / (gamagp - 1.d0) * ptot / bval(ifac,isca(ienerg))
    old_pstat = ptot

    nitermax = 100
    eps = epsrsm(ipr)
    res = 1.d0

    do while (niter.le.100.and.res.gt.eps)

      pstat =  ptot * ( 1.d0 + (gamagp - 1.d0) * 0.5d0 * bMach**2 )**(gamagp/(1.d0 - gamagp))
      y = pri / pstat

      ! 1-shock
      if (y.lt.1.d0) then

        ! Computation of the density in state 1 with Rankine-Hugoniot relations
        ro1 = roi * ((gamagp - 1.d0) * pri   + (gamagp + 1.d0) * pstat) / &
             ((gamagp - 1.d0) * pstat + (gamagp + 1.d0) * pri)

        ! Computation of the velocity in state 1 with Rankine-Hugoniot relations
        ! un1 = un2
        un1 = uni - sqrt( (pstat - pri) * (1.d0/roi - 1.d0/ro1) )

        ! Subsonic inlet
        if (un1.le.0.d0) then

          ! unb = u2
          bval(ifac,iu) = un1 / cosalp * dir(1)
          bval(ifac,iv) = un1 / cosalp * dir(2)
          bval(ifac,iw) = un1 / cosalp * dir(3)
          ! rob = ro2
          brom(ifac) = (pstat / ptot)**(1.d0/gamagp) * rhotot
          ! eb = e2
          bval(ifac,isca(ienerg)) = pstat / ((gamagp - 1.d0) * brom(ifac))                    &
                                  + 0.5d0 * ( bval(ifac,iu)**2 + bval(ifac,iv)**2 + bval(ifac,iw)**2 )
        ! Outlet
        else
          ! Computation of the shock velocity
          sigma1 = (roi * uni - ro1 * un1) / (roi - ro1)

          ! subsonic outlet
          if (sigma1.le.0.d0) then

            ! unb = u1
            bval(ifac,iu) = utxi + un1 * surfbo(1,ifac) / surfbn(ifac)
            bval(ifac,iv) = utyi + un1 * surfbo(2,ifac) / surfbn(ifac)
            bval(ifac,iw) = utzi + un1 * surfbo(3,ifac) / surfbn(ifac)
            ! rob = ro1
            brom(ifac) = ro1
            ! eb = e1
            bval(ifac,isca(ienerg)) = ei - 0.5d0 * (pstat + pri) * (1.d0 / ro1 - 1.d0 / roi) &
                                    + 0.5d0 * (un1**2 + utxi**2 + utyi**2 + utzi**2)

          ! supersonic outlet
          else

            ! pb = pri
            pstat = pri
            ! unb = uni
            bval(ifac,iu) = rtp(iel,iu)
            bval(ifac,iv) = rtp(iel,iv)
            bval(ifac,iw) = rtp(iel,iw)
            ! rob = roi
            brom(ifac) = roi
            ! eb = ei
            bval(ifac,isca(ienerg)) = rtp(iel,isca(ienerg))

          endif

        endif

        ! 1-rarefaction
      else

        ! Computation of the velocity in state 1 using Riemann invariants of the 1-rarefaction
        un1 = uni + 2 * ci / (gamagp - 1.d0) * (1.d0 - (pstat / pri) ** ((gamagp - 1.d0) / (2.d0 * gamagp)))

        ! Computation of the density in state 1 using Riemann invariants of the 1-rarefaction
        ro1 = (pstat / pri) ** (1.d0 / gamagp) * roi

        ! Subsonic inlet
        if (un1.le.0.d0) then

          ! unb = u2
          bval(ifac,iu) = un1 / cosalp * dir(1)
          bval(ifac,iv) = un1 / cosalp * dir(2)
          bval(ifac,iw) = un1 / cosalp * dir(3)
          ! rob = ro2
          brom(ifac) = (pstat / ptot)**(1.d0/gamagp) * rhotot
          ! eb = e2
          bval(ifac,isca(ienerg)) = pstat / ((gamagp - 1.d0) * brom(ifac))                    &
                                  + 0.5d0 * ( bval(ifac,iu)**2 + bval(ifac,iv)**2 + bval(ifac,iw)**2 )
        ! Outlet
        else

          ! Computation of the sound speed in state 1
          c1 = sqrt(gamagp * pstat / ro1)

          ! Subsonic outlet
          if ((un1 - c1).lt.0.d0) then

            ! unb = u1
            bval(ifac,iu) = utxi + un1 * surfbo(1,ifac) / surfbn(ifac)
            bval(ifac,iv) = utyi + un1 * surfbo(2,ifac) / surfbn(ifac)
            bval(ifac,iw) = utzi + un1 * surfbo(3,ifac) / surfbn(ifac)
            ! rob = ro1
            brom(ifac) = ro1
            ! eb = e1
            bval(ifac,isca(ienerg)) = pstat / (ro1 * (gamagp - 1.d0))                   &
                 + 0.5d0 * (un1**2 + utxi**2 + utyi**2 + utzi**2)

          ! Supersonic outlet
          else if ((uni - ci).ge.0.d0) then
            !          write(nfecra,*) 'supersonic outlet'
            ! pb = pri
            pstat = pri
            ! ub = uni
            bval(ifac,iu) = rtp(iel,iu)
            bval(ifac,iv) = rtp(iel,iv)
            bval(ifac,iw) = rtp(iel,iw)
            ! rob = roi
            brom(ifac) = roi
            ! eb = ei
            bval(ifac,isca(ienerg)) = rtp(iel,isca(ienerg))

          ! Outlet in sonic state
          else

            ! Mach number in the domain
            mi = uni / ci

            a = (gamagp - 1.d0) / (gamagp + 1.d0) * (mi + 2.d0 / (gamagp - 1))

            ! Sonic state pressure
            pstat = pri * a ** (2.d0 * gamagp / (gamagp - 1.d0))
            ! Sonic state density
            brom(ifac) = roi * a ** (2.d0 / (gamagp - 1.d0))
            ! Sonic state velocity
            uns = a * ci
            bval(ifac,iu) = uns * surfbo(1,ifac) / surfbn(ifac)
            bval(ifac,iv) = uns * surfbo(2,ifac) / surfbn(ifac)
            bval(ifac,iw) = uns * surfbo(3,ifac) / surfbn(ifac)
            ! Sonic state energy
            bval(ifac,isca(ienerg)) = pstat / ((gamagp - 1.d0) * brom(ifac)) + 0.5d0 * uns**2

          endif

        endif

      endif

      bc = sqrt(gamagp * pstat / brom(ifac))
      bMach = ( bval(ifac,iu) * surfbo(1,ifac)                             &
              + bval(ifac,iv) * surfbo(2,ifac)                             &
              + bval(ifac,iw) * surfbo(3,ifac) ) / surfbn(ifac) / bc

      bval(ifac,ipr) = pstat

      ! Pressure residual
      res = abs((pstat - old_pstat) / ptot)

      ! Prepare next iteration
      old_pstat = pstat
      niter = niter + 1

    enddo

    ! Warn the user if fixed point algorithm did not converge
    if (niter.eq.101) then
      write(nfecra,9000)ifac,res
    endif

! --- Calculation of temperature and energy from pressure and density

    ! It is postulated that the pressure and density values are
    !   strictly positive

  elseif (iccfth.eq.912.or.iccfth.eq.60900) then

    ifac = ifac0
    iel  = ifabor(ifac)

    ! Temperature
    bval(ifac, itk) = xmasml*bval(ifac,ipr)/(rr*brom(ifac))

    ! Energie totale
    bval(ifac,ien) = cv0 * bval(ifac,itk)                                       &
                     + 0.5d0 * ( bval(ifac,iu)**2                               &
                               + bval(ifac,iv)**2 + bval(ifac,iw)**2 )


! --- Calculation of density and energy from pressure and temperature

  elseif (iccfth.eq.913.or.iccfth.eq.100900) then

    ifac = ifac0
    iel  = ifabor(ifac)

    ! Density
    brom(ifac) = xmasml * bval(ifac,ipr) / ( rr * bval(ifac,itk) )

    ! Total energy
    bval(ifac,ien) = cv0 * bval(ifac,itk)                                     &
                     + 0.5d0*( bval(ifac,iu)**2                               &
                             + bval(ifac,iv)**2 + bval(ifac,iw)**2 )

! --- Calculation of density and temperature from pressure and total energy

  elseif (iccfth.eq.914.or.iccfth.eq.140900) then

    ifac = ifac0
    iel  = ifabor(ifac)

    ! Density
    brom(ifac) = bval(ifac,ipr)/((gamagp-1.d0)*                                 &
                          ( bval(ifac,ien)                                      &
                          - 0.5d0 * ( bval(ifac,iu)**2                          &
                                    + bval(ifac,iv)**2                          &
                                    + bval(ifac,iw)**2)) )

    ! Temperature
    bval(ifac,itk) = xmasml * bval(ifac,ipr) / ( rr * brom(ifac) )


! --- Calculation of pressure and energy from density and temperature

  elseif (iccfth.eq.923.or.iccfth.eq.150900) then

    ifac = ifac0
    iel  = ifabor(ifac)

    ! Pressure
    bval(ifac,ipr) = brom(ifac)*rr/xmasml * bval(ifac,itk)

    ! Total energy
      bval(ifac,ien) = cv0 * bval(ifac,itk)                                     &
                       + 0.5d0*( bval(ifac,iu)**2                               &
                               + bval(ifac,iv)**2 + bval(ifac,iw)**2 )

! --- Calculation of pressure and temperature from density and energy

  elseif (iccfth.eq.924.or.iccfth.eq.210900) then

    ifac = ifac0
    iel  = ifabor(ifac)

    ! Pressure
    bval(ifac,ipr) = (gamagp-1.d0)*brom(ifac)                                   &
                     *( bval(ifac,ien)                                          &
                      - 0.5d0*( bval(ifac,iu)**2                                &
                              + bval(ifac,iv)**2                                &
                              + bval(ifac,iw)**2 ) )

    ! Temperature
    bval(ifac,itk)= xmasml * bval(ifac,ipr) / ( rr * brom(ifac) )

! --- End of the treatment of the perfect gas
  endif

! --- End of test on the thermodynamic laws
endif


!--------
! Formats
!--------

 1010 format(                                                     &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/,                                                            &
'@ @@ WARNING:    stop in thermodynamics computations',/,         &
'@    =======',/,                                                 &
'@     Error encountered in the user subroutine ''cfther'', ',/,  &
'@       for perfect gas with constant gamma.',/,                 &
'@',/,                                                            &
'@     Gamma = ',e12.4   ,/,                                      &
'@     Gamma must be a real number greater or equal to 1.',/,     &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/)
 2010 format(                                                     &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/,                                                            &
'@ @@ WARNING:    stop in thermodynamics computations',/,         &
'@    =======',/,                                                 &
'@     Error encountered in the user subroutine ''cfther'', ',/,  &
'@       for perfect gas with constant gamma.',/,                 &
'@',/,                                                            &
'@     The computation of density failed.',/,                     &
'@',/,                                                            &
'@     Temperature = ',e12.4   ,' in cell ',i10  ,/,              &
'@     Temperature must be strictly positive.',/,                 &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/)
 3010 format(                                                     &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/,                                                            &
'@ @@ WARNING:    stop in thermodynamics computations',/,         &
'@    =======',/,                                                 &
'@     Error encountered in the user subroutine ''cfther'', ',/,  &
'@       for perfect gas with constant gamma.',/,                 &
'@',/,                                                            &
'@     The computation of temperature failed.',/,                 &
'@',/,                                                            &
'@     Density = ',e12.4   ,' in cell ',i10  ,/,                  &
'@     Density must be strictly positive.',/,                     &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/)
 4010 format(                                                     &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/,                                                            &
'@ @@ WARNING:    stop in thermodynamics computations',/,         &
'@    =======',/,                                                 &
'@     Error encountered in the user subroutine ''cfther'', ',/,  &
'@       for perfect gas with constant gamma.',/,                 &
'@',/,                                                            &
'@     The computation of the squared speed of sound failed.',/,  &
'@',/,                                                            &
'@     Density = ',e12.4   ,' in cell ',i10  ,/,                  &
'@     Density must be strictly positive.',/,                     &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/)
 4020 format(                                                     &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/,                                                            &
'@ @@ WARNING:    stop in thermodynamics computations',/,         &
'@    =======',/,                                                 &
'@     Error encountered in the user subroutine ''cfther'', ',/,  &
'@       for perfect gas with constant gamma.',/,                 &
'@',/,                                                            &
'@     The computation of the variable beta failed.',/,           &
'@',/,                                                            &
'@     Density = ',e12.4   ,' in cell ',i10  ,/,                  &
'@     Density must be strictly positive.',/,                     &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/)
 4030 format(                                                     &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/,                                                            &
'@ @@ WARNING:    stop in thermodynamics computations',/,         &
'@    =======',/,                                                 &
'@     Error encountered in the user subroutine ''cfther'', ',/,  &
'@       for perfect gas with constant gamma.',/,                 &
'@',/,                                                            &
'@     The computation of the entropy failed.',/,                 &
'@',/,                                                            &
'@     Density = ',e12.4   ,' in cell ',i10  ,/,                  &
'@     Density must be strictly positive.',/,                     &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/)
 7000 format (                                                    &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/,                                                            &
'@ @@ WARNING:    stop in thermodynamics computations',/,         &
'@    =======',/,                                                 &
'@     Error encountered in the user subroutine ''cfther'', ',/,  &
'@       for perfect gas with constant gamma.',/,                 &
'@',/,                                                            &
'@     The boundary condition of the type ''prescribed mass',/,   &
'@     and enthalpy flow rates '' is not available in the ',/,    &
'@     current release.',/,                                       &
'@',/,                                                            &
'@     Modify the user subroutine ''cfther''.',/,                 &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/)
 8000 format (                                                    &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/,                                                            &
'@ @@ WARNING:    stop in thermodynamics computations',/,         &
'@    =======',/,                                                 &
'@     Error encountered in the user subroutine ''cfther'', ',/,  &
'@       for perfect gas with constant gamma.',/,                 &
'@',/,                                                            &
'@     Negative values of the density were encountered ',/,       &
'@     in ',i10   ,' cells.',/,                                   &
'@     The density was clipped at ',e12.4  ,/                     &
'@     The run was stopped.',/,                                   &
'@',/,                                                            &
'@     If it is desired to continue the run in spite of this ',/, &
'@     behavior, it is possible to force a standard clipping ',/, &
'@     by setting a minimum value for the density variable in',/, &
'@     the GUI or in the user subroutine ''usipsu'' (set the ',/, &
'@     scamin value associated to the variable ',/,               &
'@     isca(irho).',/,                                     &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/)
 8100 format (                                                    &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/,                                                            &
'@ @@ WARNING:    stop in thermodynamics computations',/,         &
'@    =======',/,                                                 &
'@     Error encountered in the user subroutine ''cfther'', ',/,  &
'@       for perfect gas with constant gamma.',/,                 &
'@',/,                                                            &
'@     Negative values of the internal energy were encountered',/,&
'@     in ',i10   ,' cells.',/,                                   &
'@     The internal energy  was clipped at ',e12.4  ,/            &
'@     The run was stopped.',/,                                   &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/)


! The following formats may be discarded if or when the
! gamma variable option will have been fixed



9000 format(                                                     &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/,                                                            &
'@ @@ WARNING:    in thermodynamics computations',/,              &
'@    =======',/,                                                 &
'@     Error encountered in the user subroutine ''cfther'', ',/,  &
'@       for perfect gas with variable gamma.',/,                 &
'@',/,                                                            &
'@     Fixed point algorithm did not converge when',/,            &
'@     computing the subsonic inlet boundary condition',/,        &
'@     with imposed total pressure and total enthalpy.',/,        &
'@',/,                                                            &
'@     At boundary face ',i10   , /,                              &
'@     boundary Mach number residual = ', e12.4 ,/,               &
'@     maximum number of iterations (100) was reached.' ,/,       &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/)

9010 format(                                                     &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/,                                                            &
'@ @@ WARNING:    in thermodynamics computations',/,              &
'@    =======',/,                                                 &
'@     Error encountered in the user subroutine ''cfther'', ',/,  &
'@       for perfect gas with variable gamma.',/,                 &
'@',/,                                                            &
'@     The computation of the subsonic inlet boundary',/,         &
'@     condition with imposed total pressure and',/,              &
'@     total enthalpy failed.',/,                                 &
'@',/,                                                            &
'@     At boundary face ',i10   , /,                              &
'@     The direction vector given by the user can''t be null.',/, &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/)

9020 format(                                                     &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/,                                                            &
'@ @@ WARNING:    in thermodynamics computations',/,              &
'@    =======',/,                                                 &
'@     In the user subroutine ''cfther'', for perfect gas',/,  &
'@       with variable gamma,',/,                                 &
'@',/,                                                            &
'@       in the computation of the subsonic inlet with',/,        &
'@       imposed total pressure and total enthalpy.',/,           &
'@',/,                                                            &
'@     At boundary face ',i10   , /,                              &
'@     The direction vector points outward the fluid domain.',/,  &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/)

!----
! End
!----

return
end subroutine cfther
