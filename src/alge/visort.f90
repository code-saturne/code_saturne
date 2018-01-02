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

!===============================================================================
! Function :
! --------

!> \file visort.f90
!>
!> \brief Orthotropic diffusion velocity computation.
!>
!> Orthotropic diffusion velocity computation
!> viscf,b = viscosity*surface/distance, in kg/s
!>
!> \f[ viscf,b = (nx^2*visc11\_moy\_face
!>             +  ny^2*visc22\_moy\_face
!>             +  nz^2*visc33\_moy\_face)*surface/distance
!> \f]
!>
!> The viscosity is given by w1, w2, w3
!>
!> \remark There's no need for a reconstruction technique. (To improve if necessary)
!>
!> Please refer to the
!> <a href="../../theory.pdf#visort"><b>visort</b></a> section of the
!> theory guide for more informations.
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments                                                                    !
!______________________________________________________________________________!
!  mode           name          role                                           !
!______________________________________________________________________________!
! \param[in]     imvisf         Face viscosity computation method
!                               = 0 arithmetic
!                               = 1 harmonic
! \param[in]     w1,2,3         Viscosity values
! \param[out]    viscf          visc*surface/dist at internal faces
! \param[out]    viscb          visc*surface/dist at boundary faces
!______________________________________________________________________________!

subroutine visort &
 ( imvisf ,                                                       &
   w1     , w2     , w3     ,                                     &
   viscf  , viscb  )

!===============================================================================
! Module files
!===============================================================================

use paramx
use optcal, only: iporos
use numvar
use parall
use period
use mesh
use field

!===============================================================================

implicit none

! Arguments

integer          imvisf


double precision w1(ncelet), w2(ncelet), w3(ncelet)
double precision viscf(nfac), viscb(nfabor)

! Local variables

integer          ifac, ii, jj, isou
double precision visci(3), viscj(3), surf2(3)
double precision pnd

double precision, dimension(:), pointer :: porosi

!===============================================================================

! ---> Periodicity and parallelism treatment

if (irangp.ge.0.or.iperio.eq.1) then
  call syndia(w1, w2, w3)
endif

! Without porosity
if (iporos.eq.0) then

  ! Arithmetic mean
  if (imvisf.eq.0) then

    do ifac = 1, nfac

      ii = ifacel(1,ifac)
      jj = ifacel(2,ifac)

      visci(1) = w1(ii)
      viscj(1) = w1(jj)
      visci(2) = w2(ii)
      viscj(2) = w2(jj)
      visci(3) = w3(ii)
      viscj(3) = w3(jj)

      do isou = 1, 3
        surf2(isou) = surfac(isou,ifac)**2
      enddo

      viscf(ifac) = 0.5d0*(                                         &
         (visci(1)+viscj(1))*surf2(1)                               &
       + (visci(2)+viscj(2))*surf2(2)                               &
       + (visci(3)+viscj(3))*surf2(3) ) / (surfan(ifac)*dist(ifac))

    enddo

  ! Harmonic mean
  else

    do ifac = 1, nfac

      ii = ifacel(1,ifac)
      jj = ifacel(2,ifac)

      pnd  = pond(ifac)

      visci(1) = w1(ii)
      viscj(1) = w1(jj)
      visci(2) = w2(ii)
      viscj(2) = w2(jj)
      visci(3) = w3(ii)
      viscj(3) = w3(jj)

      do isou = 1, 3
        surf2(isou) = surfac(isou,ifac)**2
      enddo

      viscf(ifac) = &
        ( visci(1)*viscj(1)*surf2(1)/(pnd*visci(1)+(1.d0-pnd)*viscj(1))  &
        + visci(2)*viscj(2)*surf2(2)/(pnd*visci(2)+(1.d0-pnd)*viscj(2))  &
        + visci(3)*viscj(3)*surf2(3)/(pnd*visci(3)+(1.d0-pnd)*viscj(3))  &
        ) / (surfan(ifac)*dist(ifac))

    enddo

  endif

  do ifac = 1, nfabor

    viscb(ifac) = surfbn(ifac)

  enddo

! With porosity
else

  call field_get_val_s(ipori, porosi)

  if (irangp.ge.0.or.iperio.eq.1) then
    call synsca(porosi)
  endif

  ! Arithmetic mean
  if (imvisf.eq.0) then

    do ifac = 1, nfac

      ii = ifacel(1,ifac)
      jj = ifacel(2,ifac)

      visci(1) = w1(ii) * porosi(ii)
      viscj(1) = w1(jj) * porosi(jj)
      visci(2) = w2(ii) * porosi(ii)
      viscj(2) = w2(jj) * porosi(jj)
      visci(3) = w3(ii) * porosi(ii)
      viscj(3) = w3(jj) * porosi(jj)

      do isou = 1, 3
        surf2(isou) = surfac(isou,ifac)**2
      enddo

      viscf(ifac) = 0.5d0*(                                         &
         (visci(1)+viscj(1))*surf2(1)                               &
       + (visci(2)+viscj(2))*surf2(2)                               &
       + (visci(3)+viscj(3))*surf2(3) ) / (surfan(ifac)*dist(ifac))

    enddo

  ! Harmonic mean
  else

    do ifac = 1, nfac

      ii = ifacel(1,ifac)
      jj = ifacel(2,ifac)

      pnd  = pond(ifac)

      visci(1) = w1(ii) * porosi(ii)
      viscj(1) = w1(jj) * porosi(jj)
      visci(2) = w2(ii) * porosi(ii)
      viscj(2) = w2(jj) * porosi(jj)
      visci(3) = w3(ii) * porosi(ii)
      viscj(3) = w3(jj) * porosi(jj)

      do isou = 1, 3
        surf2(isou) = surfac(isou,ifac)**2
      enddo

      viscf(ifac) = &
        ( visci(1)*viscj(1)*surf2(1)/(pnd*visci(1)+(1.d0-pnd)*viscj(1))  &
        + visci(2)*viscj(2)*surf2(2)/(pnd*visci(2)+(1.d0-pnd)*viscj(2))  &
        + visci(3)*viscj(3)*surf2(3)/(pnd*visci(3)+(1.d0-pnd)*viscj(3))  &
        ) / (surfan(ifac)*dist(ifac))

    enddo

  endif

  do ifac = 1, nfabor

    ii = ifabor(ifac)

    viscb(ifac) = surfbn(ifac) * porosi(ii)

  enddo


endif

!----
! End
!----

return

end subroutine
