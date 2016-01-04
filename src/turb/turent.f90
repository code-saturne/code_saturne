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

!===============================================================================
! Function:
! --------
!> \file turent.f90
!> \brief Calculation of turbulent inlet conditions for a circular duct flow
!>        with smooth wall.
!>
!> \brief Calculation of \f$ u^\star \f$, \f$ k \f$ and \f$\varepsilon \f$
!>        from a diameter \f$ D_H \f$ and the reference velocity \f$ U_{ref} \f$
!>        for a circular duct flow with smooth wall
!>        (use for inlet boundary conditions).
!>
!> Both \f$ u^\star \f$ and\f$ (k,\varepsilon )\f$ are returned, so that
!> the user may compute other values of \f$ k \f$ and \f$ \varepsilon \f$
!> with the \f$ u^\star \f$.
!>
!> We use the laws coming for Idel'Cik, i.e.
!> the head losses coefficient \f$ \lambda \f$ is defined by:
!> \f[ |\dfrac{\Delta P}{\Delta x}| =
!>                        \dfrac{\lambda}{D_H} \frac{1}{2} \rho U_{ref}^2 \f]
!>
!> then  the relation reads \f$u^\star = U_{ref} \sqrt{\dfrac{\lambda}{8}}\f$.
!> \f$\lambda \f$ depends on the hydraulic Reynolds number
!> \f$ Re = \dfrac{U_{ref} D_H}{ \nu} \f$ and is given by:
!>  - for \f$ Re < 2000 \f$
!>      \f[ \lambda = \dfrac{64}{Re} \f]
!>
!>  - for \f$ Re > 4000 \f$
!>      \f[ \lambda = \dfrac{1}{( 1.8 \log_{10}(Re)-1.64 )^2} \f]
!>
!>  - for \f$ 2000 < Re < 4000 \f$, we complete by a straight line
!>      \f[ \lambda = 0.021377 + 5.3115. 10^{-6} Re \f]
!>
!>  From \f$ u^\star \f$, we can estimate \f$ k \f$ and \f$ \varepsilon\f$
!>  from the well known formulae of developped turbulence
!>
!> \f[ k = \dfrac{u^{\star 2}}{\sqrt{C_\mu}} \f]
!> \f[ \varepsilon = \dfrac{ u^{\star 3}}{(\kappa D_H /10)} \f]
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role
!______________________________________________________________________________!
!> \param[in]     uref2         square of the flow speed of reference
!> \param[in]     dh            hydraulic diameter \f$ D_H \f$
!> \param[in]     xrho          mass density \f$ \rho \f$
!> \param[in]     xmu           dynamic viscosity \f$ \nu \f$
!> \param[in]     cmu           constant \f$ C_\nu \f$
!> \param[in]     xkappa        constant \f$ \kappa \f$
!> \param[out]    ustar2        square of friction speed
!> \param[out]    xk            calculated turbulent intensity \f$ k \f$
!> \param[out]    xeps          calculated turbulent dissipation
!>                              \f$ \varepsilon \f$
!______________________________________________________________________________!


subroutine keendb &
 ( uref2, dh, xrho, xmu , cmu, xkappa, ustar2, xk, xeps )

!===============================================================================
! Module files
!===============================================================================

!===============================================================================

implicit none

! Arguments

double precision uref2, dh, xrho, xmu , ustar2, xk, xeps
double precision cmu, xkappa

! Local variables

double precision re, xlmbda

!===============================================================================

re = sqrt(uref2)*dh*xrho/xmu

if (re.lt.2000) then
  !     in this case we calculate directly \f$u*^2\f$ to avoid an issue with
  !      \f$ xlmbda= \dfrac{64}{Re} \f$ when Re->0

  ustar2 = 8.d0*xmu*sqrt(uref2)/xrho/dh

else if (re.lt.4000) then

  xlmbda = 0.021377d0 + 5.3115d-6*re
  ustar2 = uref2*xlmbda/8.d0

else

  xlmbda = 1/( 1.8d0*log(re)/log(10.d0)-1.64d0)**2
  ustar2 = uref2*xlmbda/8.d0

endif

xk   = ustar2/sqrt(cmu)
xeps = ustar2**1.5d0/(xkappa*dh*0.1d0)

!----
! End
!----

return
end subroutine

!===============================================================================
! Function:
! --------
!> \brief Calculation of \f$ u^\star\f$, \f$ k \f$ and \f$\varepsilon\f$
!>        from a diameter \f$ D_H \f$, a turbulent intensity \f$ I \f$
!>        and the reference velocity \f$ U_{ref} \f$
!>        for a circular duct flow with smooth wall
!>        (use for inlet boundary conditions).
!>
!> \f[ k = 1.5 I {U_{ref}}^2 \f]
!> \f[ \varepsilon = 10 \dfrac{{C_\mu}^{0.75} k^{1.5}}{ \kappa D_H} \f]
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role
!______________________________________________________________________________!
!> \param[in]     uref2         square of the flow velocity of reference
!> \param[in]     xintur        turbulent intensity \f$ I \f$
!> \param[in]     dh            hydraulic diameter \f$ D_H \f$
!> \param[in]     cmu           constant \f$ C_\mu \f$
!> \param[in]     xkappa        constant \f$ \kappa \f$
!> \param[out]    xk            calculated turbulent intensity \f$ k \f$
!> \param[out]    xeps          calculated turbulent disspation
!>                              \f$ \varepsilon \f$
!______________________________________________________________________________!

subroutine keenin &
 ( uref2, xintur, dh, cmu, xkappa, xk, xeps )

!===============================================================================
! Module files
!===============================================================================

!===============================================================================

implicit none

! Arguments

double precision uref2, xintur, dh, cmu, xkappa, xk, xeps

! Local variables

!===============================================================================


xk   = 1.5d0*uref2*xintur**2
xeps = 10.d0*cmu**(0.75d0)*xk**1.5d0/(xkappa*dh)

!----
! End
!----

return
end subroutine
