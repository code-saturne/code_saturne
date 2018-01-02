!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2018 EDF S.A.
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

!> \file cs_f_interfaces.f90
!> Definition of explicit interfaces for Fortran functions

module cs_f_interfaces

  !=============================================================================

  use, intrinsic :: iso_c_binding

  use cs_c_bindings, only: var_cal_opt

  implicit none

  !=============================================================================

  interface

    !---------------------------------------------------------------------------

    subroutine diften &
      (idtvar, ivar, vcopt,                                                    &
      inc, iccocg,                                                             &
      pvar, pvara, coefap, coefbp, cofafp, cofbfp,                             &
      viscf, viscb, viscel,                                                    &
      weighf, weighb,                                                          &
      smbrp)
      use mesh
      use cs_c_bindings, only: var_cal_opt
      integer :: idtvar, ivar
      type(var_cal_opt) :: vcopt
      integer :: inc, iccocg
      double precision, dimension(ncelet) :: pvar, pvara
      double precision, dimension(nfabor) :: coefap, coefbp, cofafp, cofbfp
      double precision, dimension(nfac) :: viscf
      double precision, dimension(nfabor) :: viscb
      double precision, dimension(6,ncelet), target :: viscel
      double precision, dimension(2,nfac) :: weighf
      double precision, dimension(nfabor) :: weighb
      double precision, dimension(ncelet) :: smbrp
    end subroutine diften

    !---------------------------------------------------------------------------

    subroutine itrmav &
     (f_id, init, inc, imrgra, iccocg, nswrgp, imligp, ircflp,                 &
     iphydp, iwgrp, iwarnp, epsrgp, climgp, extrap, frcxt,                     &
     pvar, coefap, coefbp, cofafp, cofbfp, viscf, viscb, viscel,               &
     weighf, weighb, flumas, flumab)
      use mesh
      integer :: f_id, init, inc, imrgra
      integer :: iccocg, nswrgp, imligp, ircflp
      integer :: iwarnp, iphydp
      double precision :: epsrgp , climgp , extrap
      double precision, dimension(ncelet) :: pvar
      double precision, dimension(nfabor) :: coefap, coefbp, cofafp, cofbfp
      double precision, dimension(nfac) :: viscf(nfac)
      double precision, dimension(nfabor) :: viscb(nfabor)
      double precision, dimension(6, ncelet), target :: viscel
      double precision, dimension(2,nfac) :: weighf
      double precision, dimension(nfabor) :: weighb
      double precision, dimension(nfac) :: flumas
      double precision, dimension(nfabor) :: flumab
      double precision, dimension(3, ncelet) :: frcxt
    end subroutine itrmav

    !---------------------------------------------------------------------------

    subroutine itrgrv &
      (f_id, init, inc, imrgra, iccocg, nswrgp, imligp, ircflp,                &
      iphydp, iwarnp,                                                          &
      epsrgp, climgp, extrap, frcxt,                                           &
      pvar, coefap, coefbp, cofafp, cofbfp, viscf, viscb, viscel,              &
      weighf, weighb, diverg)
      use mesh
      integer :: f_id, init, inc, imrgra, iccocg, nswrgp, imligp, ircflp
      integer :: iwarnp , iphydp
      double precision :: epsrgp, climgp, extrap
      double precision, dimension(ncelet) :: pvar
      double precision, dimension(nfabor) :: coefap, coefbp, cofafp, cofbfp
      double precision, dimension(nfac) :: viscf
      double precision, dimension(nfabor) :: viscb
      double precision, dimension(6,ncelet), target :: viscel
      double precision, dimension(2,nfac) :: weighf
      double precision, dimension(nfabor) :: weighb
      double precision, dimension(ncelet) :: diverg
      double precision, dimension(3,ncelet) :: frcxt
    end subroutine itrgrv

    !---------------------------------------------------------------------------

    subroutine matrix &
      (iconvp, idiffp, ndircp, isym, thetap, imucpp, coefbp, cofbfp,           &
      rovsdt, i_massflux, b_massflux, i_visc, b_visc, xcpp, da, xa)
      use mesh
      integer :: iconvp, idiffp, ndircp, isym, imucpp
      double precision :: thetap
      double precision, dimension(ncelet) :: rovsdt, xcpp, da
      double precision, dimension(nfabor) :: coefbp, cofbfp, b_massflux, b_visc
      double precision, dimension(nfac) :: i_massflux, i_visc
      double precision, dimension(2,nfac) :: xa
    end subroutine matrix

    !---------------------------------------------------------------------------

    subroutine matrdt &
      (iconvp, idiffp, isym, coefbp, cofbfp,                                   &
      i_massflux, b_massflux, i_visc, b_visc, da)
      use mesh
      integer :: iconvp, idiffp, isym
      double precision, dimension(ncelet) :: da
      double precision, dimension(nfabor) :: coefbp, cofbfp, b_massflux, b_visc
      double precision, dimension(nfac) :: i_massflux, i_visc
    end subroutine matrdt

    !---------------------------------------------------------------------------

    subroutine post_boundary_thermal_flux &
      (nfbrps, lstfbr, bflux)
      use dimens
      use mesh
      integer, intent(in)                                        :: nfbrps
      integer, dimension(nfbrps), intent(in)                     :: lstfbr
      double precision, dimension(nfbrps), intent(out)           :: bflux
    end subroutine post_boundary_thermal_flux

    !---------------------------------------------------------------------------

    subroutine post_boundary_nusselt &
      (nfbrps, lstfbr, bnussl)
      use dimens
      use mesh
      integer, intent(in)                                        :: nfbrps
      integer, dimension(nfbrps), intent(in)                     :: lstfbr
      double precision, dimension(nfbrps), intent(out)           :: bnussl
    end subroutine post_boundary_nusselt

    !---------------------------------------------------------------------------

    subroutine post_stress &
      (nfbrps, lstfbr, stress)
      use dimens
      use mesh
      integer, intent(in)                                 :: nfbrps
      integer, dimension(nfbrps), intent(in)              :: lstfbr
      double precision, dimension(3, nfbrps), intent(out) :: stress
    end subroutine post_stress

    !---------------------------------------------------------------------------

    subroutine turrij &
      (nvar, nscal, ncepdp, ncesmp, icepdc, icetsm, itypsm,                    &
      dt, tslagr,                                                              &
      ckupdc, smacel)
      use lagran, only: ntersl
      use mesh
      integer :: nvar, nscal, ncepdp, ncesmp
      integer, dimension(ncepdp) :: icepdc
      integer, dimension(ncesmp) :: icetsm
      integer, dimension(ncesmp,nvar), target :: itypsm
      double precision, dimension(ncelet) :: dt
      double precision, dimension(ncelet,ntersl), target :: tslagr
      double precision, dimension(ncepdp,6) :: ckupdc
      double precision, dimension(ncesmp,nvar), target ::  smacel
    end subroutine turrij

    !---------------------------------------------------------------------------

    subroutine vitens &
     (w1, iwarnp, weighf, weighb, viscf, viscb)
      use mesh
      integer :: iwarnp
      double precision, dimension(6, ncelet), target :: w1
      double precision, dimension(2, nfac) :: weighf
      double precision, dimension(nfabor) :: weighb
      double precision, dimension(nfac) :: viscf
      double precision, dimension(nfabor) :: viscb
    end subroutine vitens

    !---------------------------------------------------------------------------

    subroutine vistnv &
      (imvisf, w1, viscf, viscb)
      use mesh
      integer :: imvisf
      double precision, dimension(6,ncelet), target :: w1
      double precision, dimension(3,3,nfac) :: viscf
      double precision, dimension(nfabor) :: viscb
    end subroutine vistnv

    !---------------------------------------------------------------------------

  end interface

  !=============================================================================

end module cs_f_interfaces
