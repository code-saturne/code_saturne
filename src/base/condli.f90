!-------------------------------------------------------------------------------

! This file is part of code_saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2023 EDF S.A.
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

!> \file condli.f90
!>
!> \brief Translation of the boundary conditions given by
!> cs_user_boundary_conditions
!> in a form that fits to the solver.
!>

!===============================================================================
! Local functions
!===============================================================================

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[out]    coefa         explicit BC coefficient for gradients
!> \param[out]    cofaf         explicit BC coefficient for diffusive flux
!> \param[out]    coefb         implicit BC coefficient for gradients
!> \param[out]    cofbf         implicit BC coefficient for diffusive flux
!> \param[in]     pimp          Dirichlet value to impose
!> \param[in]     hint          Internal exchange coefficient
!> \param[in]     hext          External exchange coefficient (10^30 by default)
!_______________________________________________________________________________

subroutine set_dirichlet_scalar &
 ( coefa , cofaf, coefb , cofbf, pimp  , hint, hext)

!===============================================================================
! Module files
!===============================================================================

use cstnum, only: rinfin

!===============================================================================

implicit none

! Arguments

double precision coefa, cofaf, coefb, cofbf, pimp, hint, hext

! Local variables

double precision heq

!===============================================================================

if (abs(hext).gt.rinfin*0.5d0) then

  ! Gradient BCs
  coefa = pimp
  coefb = 0.d0

  ! Flux BCs
  cofaf = -hint*pimp
  cofbf =  hint

else

  ! Gradient BCs
  coefa = hext*pimp/(hint + hext)
  coefb = hint     /(hint + hext)

  ! Flux BCs
  heq = hint*hext/(hint + hext)
  cofaf = -heq*pimp
  cofbf =  heq

endif

return
end subroutine set_dirichlet_scalar

!===============================================================================

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[out]    coefa         explicit BC coefficient for gradients
!> \param[out]    cofaf         explicit BC coefficient for diffusive flux
!> \param[out]    coefb         implicit BC coefficient for gradients
!> \param[out]    cofbf         implicit BC coefficient for diffusive flux
!> \param[in]     pimpv         Dirichlet value to impose
!> \param[in]     hint          Internal exchange coefficient
!> \param[in]     hextv         External exchange coefficient (10^30 by default)
!_______________________________________________________________________________

subroutine set_dirichlet_vector &
 ( coefa , cofaf, coefb , cofbf, pimpv  , hint , hextv)

!===============================================================================
! Module files
!===============================================================================

use cstnum, only: rinfin

!===============================================================================

implicit none

! Arguments

double precision coefa(3), cofaf(3)
double precision coefb(3,3), cofbf(3,3)
double precision pimpv(3)
double precision hint
double precision hextv(3)

! Local variables

integer          isou  , jsou
double precision heq

!===============================================================================

do isou = 1, 3

  if (abs(hextv(isou)).gt.rinfin*0.5d0) then
    ! Gradient BCs
    coefa(isou) = pimpv(isou)
    do jsou = 1, 3
      coefb(isou,jsou) = 0.d0
    enddo

    ! Flux BCs
    cofaf(isou) = -hint*pimpv(isou)
    do jsou = 1, 3
      if (jsou.eq.isou) then
        cofbf(isou,jsou) = hint
      else
        cofbf(isou,jsou) = 0.d0
      endif
    enddo

  else

    heq = hint*hextv(isou)/(hint + hextv(isou))

    ! Gradient BCs
    coefa(isou) = hextv(isou)*pimpv(isou)/(hint + hextv(isou))
    do jsou = 1, 3
      if (jsou.eq.isou) then
        coefb(isou,jsou) = hint/(hint + hextv(isou))
      else
        coefb(isou,jsou) = 0.d0
      endif
    enddo

    ! Flux BCs
    cofaf(isou) = -heq*pimpv(isou)
    do jsou = 1, 3
      if (jsou.eq.isou) then
        cofbf(isou,jsou) = heq
      else
        cofbf(isou,jsou) = 0.d0
      endif
    enddo

  endif

enddo

return
end subroutine set_dirichlet_vector

!===============================================================================

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[out]    coefa         explicit BC coefficient for gradients
!> \param[out]    cofaf         explicit BC coefficient for diffusive flux
!> \param[out]    coefb         implicit BC coefficient for gradients
!> \param[out]    cofbf         implicit BC coefficient for diffusive flux
!> \param[in]     pimpts        Dirichlet value to impose
!> \param[in]     hint          Internal exchange coefficient
!> \param[in]     hextts        External exchange coefficient (10^30 by default)
!_______________________________________________________________________________

subroutine set_dirichlet_tensor &
 ( coefa , cofaf, coefb , cofbf, pimpts  , hint , hextts)

!===============================================================================
! Module files
!===============================================================================

use cstnum, only: rinfin

!===============================================================================

implicit none

! Arguments

double precision coefa(6), cofaf(6)
double precision coefb(6,6), cofbf(6,6)
double precision pimpts(6)
double precision hint
double precision hextts(6)

! Local variables

integer          isou  , jsou
double precision heq

!===============================================================================

do isou = 1, 6

  if (abs(hextts(isou)).gt.rinfin*0.5d0) then
    ! Gradient BCs
    coefa(isou) = pimpts(isou)
    do jsou = 1, 6
      coefb(isou,jsou) = 0.d0
    enddo

    ! Flux BCs
    cofaf(isou) = -hint*pimpts(isou)
    do jsou = 1, 6
      if (jsou.eq.isou) then
        cofbf(isou,jsou) = hint
      else
        cofbf(isou,jsou) = 0.d0
      endif
    enddo

  else

    heq = hint*hextts(isou)/(hint + hextts(isou))

    ! Gradient BCs
    coefa(isou) = hextts(isou)*pimpts(isou)/(hint + hextts(isou))
    do jsou = 1, 6
      if (jsou.eq.isou) then
        coefb(isou,jsou) = hint/(hint + hextts(isou))
      else
        coefb(isou,jsou) = 0.d0
      endif
    enddo

    ! Flux BCs
    cofaf(isou) = -heq*pimpts(isou)
    do jsou = 1, 6
      if (jsou.eq.isou) then
        cofbf(isou,jsou) = heq
      else
        cofbf(isou,jsou) = 0.d0
      endif
    enddo

  endif

enddo

return
end subroutine set_dirichlet_tensor

!===============================================================================

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[out]    coefa         explicit BC coefficient for gradients
!> \param[out]    cofaf         explicit BC coefficient for diffusive flux
!> \param[out]    coefb         implicit BC coefficient for gradients
!> \param[out]    cofbf         implicit BC coefficient for diffusive flux
!> \param[in]     pimpv         Dirichlet value to impose
!> \param[in]     hint          Internal exchange coefficient
!> \param[in]     hextv         External exchange coefficient (10^30 by default)
!_______________________________________________________________________________

subroutine set_dirichlet_vector_aniso &
 ( coefa , cofaf, coefb , cofbf, pimpv  , hint , hextv)

!===============================================================================
! Module files
!===============================================================================

use cstnum, only: rinfin

!===============================================================================

implicit none

! Arguments

double precision coefa(3), cofaf(3)
double precision coefb(3,3), cofbf(3,3)
double precision pimpv(3)
double precision hint(6)
double precision hextv(3)

! Local variables

integer          isou  , jsou

!===============================================================================

do isou = 1, 3

  if (abs(hextv(isou)).gt.rinfin*0.5d0) then
    ! Gradient BCs
    coefa(isou) = pimpv(isou)
    do jsou = 1, 3
      coefb(isou,jsou) = 0.d0
    enddo

  else

    call csexit(1)

  endif

enddo

! Flux BCs
cofaf(1) = -(hint(1)*pimpv(1) + hint(4)*pimpv(2) + hint(6)*pimpv(3))
cofaf(2) = -(hint(4)*pimpv(1) + hint(2)*pimpv(2) + hint(5)*pimpv(3))
cofaf(3) = -(hint(6)*pimpv(1) + hint(5)*pimpv(2) + hint(3)*pimpv(3))
cofbf(1,1) = hint(1)
cofbf(2,2) = hint(2)
cofbf(3,3) = hint(3)
cofbf(1,2) = hint(4)
cofbf(2,1) = hint(4)
cofbf(2,3) = hint(5)
cofbf(3,2) = hint(5)
cofbf(1,3) = hint(6)
cofbf(3,1) = hint(6)

return
end subroutine set_dirichlet_vector_aniso

!===============================================================================

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[out]    coefa         explicit BC coefficient for gradients
!> \param[out]    cofaf         explicit BC coefficient for diffusive flux
!> \param[out]    coefb         implicit BC coefficient for gradients
!> \param[out]    cofbf         implicit BC coefficient for diffusive flux
!> \param[in]     dimp          Flux value to impose
!> \param[in]     hint          Internal exchange coefficient
!_______________________________________________________________________________

subroutine set_neumann_scalar &
 ( coefa , cofaf, coefb , cofbf, dimp  , hint)

!===============================================================================
! Module files
!===============================================================================

!===============================================================================

implicit none

! Arguments

double precision coefa, cofaf, coefb, cofbf, dimp, hint

! Local variables

!===============================================================================

! Gradient BCs
coefa = -dimp/max(hint, 1.d-300)
coefb = 1.d0

! Flux BCs
cofaf = dimp
cofbf = 0.d0

return
end subroutine set_neumann_scalar

!===============================================================================

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[out]    coefa         explicit BC coefficient for gradients
!> \param[out]    cofaf         explicit BC coefficient for diffusive flux
!> \param[out]    coefb         implicit BC coefficient for gradients
!> \param[out]    cofbf         implicit BC coefficient for diffusive flux
!> \param[in]     qimpv         Flux value to impose
!> \param[in]     hint          Internal exchange coefficient
!_______________________________________________________________________________

subroutine set_neumann_vector &
 ( coefa , cofaf, coefb , cofbf, qimpv  , hint)

!===============================================================================
! Module files
!===============================================================================

!===============================================================================

implicit none

! Arguments

double precision coefa(3), cofaf(3)
double precision coefb(3,3), cofbf(3,3)
double precision qimpv(3)
double precision hint

! Local variables

integer          isou  , jsou

!===============================================================================

do isou = 1, 3

  ! Gradient BCs
  coefa(isou) = -qimpv(isou)/max(hint, 1.d-300)
  do jsou = 1, 3
    if (jsou.eq.isou) then
      coefb(isou,jsou) = 1.d0
    else
      coefb(isou,jsou) = 0.d0
    endif
  enddo

  ! Flux BCs
  cofaf(isou) = qimpv(isou)
  do jsou = 1, 3
    cofbf(isou,jsou) = 0.d0
  enddo

enddo

return
end subroutine set_neumann_vector

!===============================================================================

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[out]    coefa         explicit BC coefficient for gradients
!> \param[out]    cofaf         explicit BC coefficient for diffusive flux
!> \param[out]    coefb         implicit BC coefficient for gradients
!> \param[out]    cofbf         implicit BC coefficient for diffusive flux
!> \param[in]     qimpts         Flux value to impose
!> \param[in]     hint          Internal exchange coefficient
!_______________________________________________________________________________

subroutine set_neumann_tensor &
 ( coefa , cofaf, coefb , cofbf, qimpts  , hint)

!===============================================================================
! Module files
!===============================================================================

!===============================================================================

implicit none

! Arguments

double precision coefa(6), cofaf(6)
double precision coefb(6,6), cofbf(6,6)
double precision qimpts(6)
double precision hint

! Local variables

integer          isou  , jsou

!===============================================================================

do isou = 1, 6

  ! Gradient BCs
  coefa(isou) = -qimpts(isou)/max(hint, 1.d-300)
  do jsou = 1, 6
    if (jsou.eq.isou) then
      coefb(isou,jsou) = 1.d0
    else
      coefb(isou,jsou) = 0.d0
    endif
  enddo

  ! Flux BCs
  cofaf(isou) = qimpts(isou)
  do jsou = 1, 6
    cofbf(isou,jsou) = 0.d0
  enddo

enddo

return
end subroutine set_neumann_tensor

!===============================================================================

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[out]    coefa         explicit BC coefficient for gradients
!> \param[out]    cofaf         explicit BC coefficient for diffusive flux
!> \param[out]    coefb         implicit BC coefficient for gradients
!> \param[out]    cofbf         implicit BC coefficient for diffusive flux
!> \param[in]     qimpv         Flux value to impose
!> \param[in]     hint          Internal exchange coefficient
!_______________________________________________________________________________

subroutine set_neumann_vector_aniso &
 ( coefa , cofaf, coefb , cofbf, qimpv  , hint)

!===============================================================================
! Module files
!===============================================================================

!===============================================================================

implicit none

! Arguments

double precision coefa(3), cofaf(3)
double precision coefb(3,3), cofbf(3,3)
double precision qimpv(3)
double precision hint(6)

! Local variables

integer          isou  , jsou
double precision invh(6), invdet, m(6)

!===============================================================================
m(1) = hint(2)*hint(3) - hint(5)*hint(5)
m(2) = hint(1)*hint(3) - hint(6)*hint(6)
m(3) = hint(1)*hint(2) - hint(4)*hint(4)
m(4) = hint(5)*hint(6) - hint(4)*hint(3)
m(5) = hint(4)*hint(6) - hint(1)*hint(5)
m(6) = hint(4)*hint(5) - hint(2)*hint(6)

invdet = 1.d0/(hint(1)*m(1) + hint(4)*m(4) + hint(6)*m(6))

invh(1) = m(1) * invdet
invh(2) = m(2) * invdet
invh(3) = m(3) * invdet
invh(4) = m(4) * invdet
invh(5) = m(5) * invdet
invh(6) = m(6) * invdet

! Gradient BCs
coefa(1) = -(invh(1)*qimpv(1) + invh(4)*qimpv(2) + invh(6)*qimpv(3))
coefa(2) = -(invh(4)*qimpv(1) + invh(2)*qimpv(2) + invh(5)*qimpv(3))
coefa(3) = -(invh(6)*qimpv(1) + invh(5)*qimpv(2) + invh(3)*qimpv(3))
coefb(1,1) = 1.d0
coefb(2,2) = 1.d0
coefb(3,3) = 1.d0
coefb(1,2) = 0.d0
coefb(2,1) = 0.d0
coefb(2,3) = 0.d0
coefb(3,2) = 0.d0
coefb(1,3) = 0.d0
coefb(3,1) = 0.d0

do isou = 1, 3

  ! Flux BCs
  cofaf(isou) = qimpv(isou)
  do jsou = 1, 3
    cofbf(isou,jsou) = 0.d0
  enddo

enddo

return
end subroutine set_neumann_vector_aniso

!===============================================================================

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[out]    coefa         explicit BC coefficient for gradients
!> \param[out]    cofaf         explicit BC coefficient for diffusive flux
!> \param[out]    coefb         implicit BC coefficient for gradients
!> \param[out]    cofbf         implicit BC coefficient for diffusive flux
!> \param[in]     pimpv         Dirichlet value to impose on the normal
!>                              component
!> \param[in]     qimpv         Flux value to impose on the
!>                              tangential components
!> \param[in]     hint          Internal exchange coefficient
!> \param[in]     normal        normal
!_______________________________________________________________________________

subroutine set_generalized_sym_vector &
 ( coefa , cofaf, coefb , cofbf, pimpv, qimpv, hint, normal)

!===============================================================================
! Module files
!===============================================================================

!===============================================================================

implicit none

! Arguments

double precision coefa(3), cofaf(3)
double precision coefb(3,3), cofbf(3,3)
double precision hint
double precision normal(3)
double precision pimpv(3), qimpv(3)

! Local variables

integer          isou  , jsou

!===============================================================================

do isou = 1, 3

  ! Gradient BCs
  coefa(isou) = - qimpv(isou)/max(hint, 1.d-300)
    ! "[1 -n(x)n] Qimp / hint" is divided into two
  do jsou = 1, 3
    coefa(isou) = coefa(isou) &
      + normal(isou)*normal(jsou)*(pimpv(jsou)+qimpv(jsou)/max(hint, 1.d-300))
    if (jsou.eq.isou) then
      coefb(isou,jsou) = 1.d0 - normal(isou)*normal(jsou)
    else
      coefb(isou,jsou) = - normal(isou)*normal(jsou)
    endif
  enddo

  ! Flux BCs
  cofaf(isou) = qimpv(isou)
    ! "[1 -n(x)n] Qimp" is divided into two
  do jsou = 1, 3
    cofaf(isou) = cofaf(isou) &
      - normal(isou)*normal(jsou)*(hint*pimpv(jsou)+qimpv(jsou))
    cofbf(isou,jsou) = hint*normal(isou)*normal(jsou)
  enddo

enddo

return
end subroutine set_generalized_sym_vector

!===============================================================================

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[out]    coefa         explicit BC coefficient for gradients
!> \param[out]    cofaf         explicit BC coefficient for diffusive flux
!> \param[out]    coefb         implicit BC coefficient for gradients
!> \param[out]    cofbf         implicit BC coefficient for diffusive flux
!> \param[in]     pimpv         Dirichlet value to impose on the normal
!>                              component
!> \param[in]     qimpv         Flux value to impose on the
!>                              tangential components
!> \param[in]     hint          Internal exchange coefficient
!> \param[in]     normal        normal
!_______________________________________________________________________________

subroutine set_generalized_sym_vector_aniso &
 ( coefa , cofaf, coefb , cofbf, pimpv, qimpv, hint, normal)

!===============================================================================
! Module files
!===============================================================================

!===============================================================================

implicit none

! Arguments

double precision coefa(3), cofaf(3)
double precision coefb(3,3), cofbf(3,3)
double precision hint(6)
double precision normal(3)
double precision pimpv(3), qimpv(3)

! Local variables

integer          isou  , jsou
double precision invh(6), invdet, m(6), qshint(3), hintpv(3), hintnm(3)

!===============================================================================

m(1) = hint(2)*hint(3) - hint(5)*hint(5)
m(2) = hint(1)*hint(3) - hint(6)*hint(6)
m(3) = hint(1)*hint(2) - hint(4)*hint(4)
m(4) = hint(5)*hint(6) - hint(4)*hint(3)
m(5) = hint(4)*hint(6) - hint(1)*hint(5)
m(6) = hint(4)*hint(5) - hint(2)*hint(6)

invdet = 1.d0/(hint(1)*m(1) + hint(4)*m(4) + hint(6)*m(6))

invh(1) = m(1) * invdet
invh(2) = m(2) * invdet
invh(3) = m(3) * invdet
invh(4) = m(4) * invdet
invh(5) = m(5) * invdet
invh(6) = m(6) * invdet

qshint(1) = invh(1)*qimpv(1) + invh(4)*qimpv(2) + invh(6)*qimpv(3)
qshint(2) = invh(4)*qimpv(1) + invh(2)*qimpv(2) + invh(5)*qimpv(3)
qshint(3) = invh(6)*qimpv(1) + invh(5)*qimpv(2) + invh(3)*qimpv(3)

hintpv(1) = hint(1)*pimpv(1) + hint(4)*pimpv(2) + hint(6)*pimpv(3)
hintpv(2) = hint(4)*pimpv(1) + hint(2)*pimpv(2) + hint(5)*pimpv(3)
hintpv(3) = hint(6)*pimpv(1) + hint(5)*pimpv(2) + hint(3)*pimpv(3)

hintnm(1) = hint(1)*normal(1) + hint(4)*normal(2) + hint(6)*normal(3)
hintnm(2) = hint(4)*normal(1) + hint(2)*normal(2) + hint(5)*normal(3)
hintnm(3) = hint(6)*normal(1) + hint(5)*normal(2) + hint(3)*normal(3)

do isou = 1, 3

  ! Gradient BCs
  coefa(isou) = - qshint(isou)
    ! "[1 -n(x)n] Qimp / hint" is divided into two
  do jsou = 1, 3
    coefa(isou) = coefa(isou) &
      + normal(isou)*normal(jsou)*(pimpv(jsou)+qshint(jsou))
    if (jsou.eq.isou) then
      coefb(isou,jsou) = 1.d0 - normal(isou)*normal(jsou)
    else
      coefb(isou,jsou) = - normal(isou)*normal(jsou)
    endif
  enddo

  ! Flux BCs
  cofaf(isou) = qimpv(isou)
    ! "[1 -n(x)n] Qimp" is divided into two
  do jsou = 1, 3
    cofaf(isou) = cofaf(isou) &
      - normal(isou)*normal(jsou)*(hintpv(jsou)+qimpv(jsou))
    cofbf(isou,jsou) = hintnm(isou)*normal(jsou)
  enddo

enddo

return
end subroutine set_generalized_sym_vector_aniso

!===============================================================================

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[out]    coefa         explicit BC coefficient for gradients
!> \param[out]    cofaf         explicit BC coefficient for diffusive flux
!> \param[out]    coefb         implicit BC coefficient for gradients
!> \param[out]    cofbf         implicit BC coefficient for diffusive flux
!> \param[in]     pimpv         Dirichlet value to impose on the tangential
!>                              components
!> \param[in]     qimpv         Flux value to impose on the
!>                              normal component
!> \param[in]     hint          Internal exchange coefficient
!> \param[in]     normal        normal
!_______________________________________________________________________________

subroutine set_generalized_dirichlet_vector &
 ( coefa , cofaf, coefb , cofbf, pimpv, qimpv, hint, normal)

!===============================================================================
! Module files
!===============================================================================

!===============================================================================

implicit none

! Arguments

double precision coefa(3), cofaf(3)
double precision coefb(3,3), cofbf(3,3)
double precision hint
double precision normal(3)
double precision pimpv(3), qimpv(3)

! Local variables

integer          isou  , jsou

!===============================================================================

do isou = 1, 3

  ! Gradient BCs
  ! "[1 -n(x)n] Pimp" is divided into two
  coefa(isou) = pimpv(isou)
  do jsou = 1, 3
    coefa(isou) = coefa(isou) &
      - normal(isou)*normal(jsou)*(pimpv(jsou)+qimpv(jsou)/max(hint, 1.d-300))
    coefb(isou,jsou) = normal(isou)*normal(jsou)
  enddo

  ! Flux BCs
  ! "[1 -n(x)n] Pimp" is divided into two
  cofaf(isou) = -hint*pimpv(isou)
  do jsou = 1, 3
    cofaf(isou) = cofaf(isou) &
      + normal(isou)*normal(jsou)*(qimpv(jsou)+pimpv(jsou)*hint)
    if (jsou.eq.isou) then
      cofbf(isou,jsou) = hint*(1.d0-normal(isou)*normal(jsou))
    else
      cofbf(isou,jsou) = -hint*normal(isou)*normal(jsou)
    endif
  enddo

enddo

return
end subroutine set_generalized_dirichlet_vector

!===============================================================================

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[out]    coefa         explicit BC coefficient for gradients
!> \param[out]    cofaf         explicit BC coefficient for diffusive flux
!> \param[out]    coefb         implicit BC coefficient for gradients
!> \param[out]    cofbf         implicit BC coefficient for diffusive flux
!> \param[in]     pimpv         Dirichlet value to impose on the tangential
!>                              components
!> \param[in]     qimpv         Flux value to impose on the
!>                              normal component
!> \param[in]     hint          Internal exchange coefficient
!> \param[in]     normal        normal
!_______________________________________________________________________________

subroutine set_generalized_dirichlet_vector_aniso &
 ( coefa , cofaf, coefb , cofbf, pimpv, qimpv, hint, normal)

!===============================================================================
! Module files
!===============================================================================

!===============================================================================

implicit none

! Arguments

double precision coefa(3), cofaf(3)
double precision coefb(3,3), cofbf(3,3)
double precision hint(6)
double precision normal(3)
double precision pimpv(3), qimpv(3)

! Local variables

integer          isou  , jsou
double precision invh(6), invdet, m(6), qshint(3), hintpv(3), hintnm(3)

!===============================================================================

m(1) = hint(2)*hint(3) - hint(5)*hint(5)
m(2) = hint(1)*hint(3) - hint(6)*hint(6)
m(3) = hint(1)*hint(2) - hint(4)*hint(4)
m(4) = hint(5)*hint(6) - hint(4)*hint(3)
m(5) = hint(4)*hint(6) - hint(1)*hint(5)
m(6) = hint(4)*hint(5) - hint(2)*hint(6)

invdet = 1.d0/(hint(1)*m(1) + hint(4)*m(4) + hint(6)*m(6))

invh(1) = m(1) * invdet
invh(2) = m(2) * invdet
invh(3) = m(3) * invdet
invh(4) = m(4) * invdet
invh(5) = m(5) * invdet
invh(6) = m(6) * invdet

qshint(1) = invh(1)*qimpv(1) + invh(4)*qimpv(2) + invh(6)*qimpv(3)
qshint(2) = invh(4)*qimpv(1) + invh(2)*qimpv(2) + invh(5)*qimpv(3)
qshint(3) = invh(6)*qimpv(1) + invh(5)*qimpv(2) + invh(3)*qimpv(3)

hintpv(1) = hint(1)*pimpv(1) + hint(4)*pimpv(2) + hint(6)*pimpv(3)
hintpv(2) = hint(4)*pimpv(1) + hint(2)*pimpv(2) + hint(5)*pimpv(3)
hintpv(3) = hint(6)*pimpv(1) + hint(5)*pimpv(2) + hint(3)*pimpv(3)

hintnm(1) = hint(1)*normal(1) + hint(4)*normal(2) + hint(6)*normal(3)
hintnm(2) = hint(4)*normal(1) + hint(2)*normal(2) + hint(5)*normal(3)
hintnm(3) = hint(6)*normal(1) + hint(5)*normal(2) + hint(3)*normal(3)

do isou = 1, 3

  ! Gradient BCs
  ! "[1 -n(x)n] Pimp" is divided into two
  coefa(isou) = pimpv(isou)
  do jsou = 1, 3
    coefa(isou) = coefa(isou) &
      - normal(isou)*normal(jsou)*(pimpv(jsou)+qshint(jsou))
    coefb(isou,jsou) = normal(isou)*normal(jsou)
  enddo

  ! Flux BCs
  ! "[1 -n(x)n] Pimp" is divided into two
  cofaf(isou) = -hintpv(isou)
  do jsou = 1, 3
    cofaf(isou) = cofaf(isou) &
      + normal(isou)*normal(jsou)*(qimpv(jsou)+hintpv(jsou))
    if (jsou.eq.isou) then
      cofbf(isou,jsou) = hint(isou)-hintnm(isou)*normal(jsou)
    else
      cofbf(isou,jsou) = -hintnm(isou)*normal(jsou)
    endif
  enddo

enddo

return
end subroutine set_generalized_dirichlet_vector_aniso

!===============================================================================

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[out]    coefa         explicit BC coefficient for gradients
!> \param[out]    cofaf         explicit BC coefficient for diffusive flux
!> \param[out]    coefb         implicit BC coefficient for gradients
!> \param[out]    cofbf         implicit BC coefficient for diffusive flux
!> \param[in]     pimp          Flux value to impose
!> \param[in]     cfl           Local Courant number used to convect
!> \param[in]     hint          Internal exchange coefficient
!_______________________________________________________________________________

subroutine set_convective_outlet_scalar &
 ( coefa , cofaf, coefb , cofbf, pimp  , cfl   , hint)

!===============================================================================
! Module files
!===============================================================================

!===============================================================================

implicit none

! Arguments

double precision coefa, cofaf, coefb, cofbf, pimp, cfl, hint

! Local variables

!===============================================================================

! Gradient BCs
coefb = cfl/(1.d0+cfl)
coefa = (1.d0-coefb)*pimp

! Flux BCs
cofaf = -hint*coefa
cofbf =  hint*(1.d0 - coefb)

return
end subroutine

!===============================================================================

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[out]    coefa         explicit BC coefficient for gradients
!> \param[out]    cofaf         explicit BC coefficient for diffusive flux
!> \param[out]    coefb         implicit BC coefficient for gradients
!> \param[out]    cofbf         implicit BC coefficient for diffusive flux
!> \param[in]     pimpv         Dirichlet value to impose
!> \param[in]     cflv          Local Courant number used to convect
!> \param[in]     hint          Internal exchange coefficient
!_______________________________________________________________________________

subroutine set_convective_outlet_vector &
 ( coefa , cofaf, coefb , cofbf, pimpv  , cflv  , hint)

!===============================================================================
! Module files
!===============================================================================

!===============================================================================

implicit none

! Arguments

double precision coefa(3), cofaf(3)
double precision coefb(3,3), cofbf(3,3)
double precision pimpv(3), cflv(3)
double precision hint

! Local variables

integer          isou  , jsou

!===============================================================================

do isou = 1, 3

  ! Gradient BCs
  do jsou = 1, 3
    if (jsou.eq.isou) then
      coefb(isou,jsou) = cflv(isou)/(1.d0+cflv(isou))
    else
      coefb(isou,jsou) = 0.d0
    endif
  enddo
  coefa(isou) = (1.d0-coefb(isou,isou))*pimpv(isou)

  ! Flux BCs
  cofaf(isou) = -hint*coefa(isou)
  do jsou = 1, 3
    if (jsou.eq.isou) then
      cofbf(isou,jsou) = hint*(1.d0 - coefb(isou,jsou))
    else
      cofbf(isou,jsou) = 0.d0
    endif
  enddo

enddo

return
end subroutine

!===============================================================================

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[out]    coefa         explicit BC coefficient for gradients
!> \param[out]    cofaf         explicit BC coefficient for diffusive flux
!> \param[out]    coefb         implicit BC coefficient for gradients
!> \param[out]    cofbf         implicit BC coefficient for diffusive flux
!> \param[in]     pimpts         Dirichlet value to impose
!> \param[in]     cflts          Local Courant number used to convect
!> \param[in]     hint          Internal exchange coefficient
!_______________________________________________________________________________

subroutine set_convective_outlet_tensor &
 ( coefa , cofaf, coefb , cofbf, pimpts  , cflts  , hint)

!===============================================================================
! Module files
!===============================================================================

!===============================================================================

implicit none

! Arguments

double precision coefa(6), cofaf(6)
double precision coefb(6,6), cofbf(6,6)
double precision pimpts(6), cflts(6)
double precision hint

! Local variables

integer          isou  , jsou

!===============================================================================

do isou = 1, 6

  ! Gradient BCs
  do jsou = 1, 6
    if (jsou.eq.isou) then
      coefb(isou,jsou) = cflts(isou)/(1.d0+cflts(isou))
    else
      coefb(isou,jsou) = 0.d0
    endif
  enddo
  coefa(isou) = (1.d0-coefb(isou,isou))*pimpts(isou)

  ! Flux BCs
  cofaf(isou) = -hint*coefa(isou)
  do jsou = 1, 6
    if (jsou.eq.isou) then
      cofbf(isou,jsou) = hint*(1.d0 - coefb(isou,jsou))
    else
      cofbf(isou,jsou) = 0.d0
    endif
  enddo

enddo

return
end subroutine


!===============================================================================

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[out]    coefa         explicit BC coefficient for gradients
!> \param[out]    cofaf         explicit BC coefficient for diffusive flux
!> \param[out]    coefb         implicit BC coefficient for gradients
!> \param[out]    cofbf         implicit BC coefficient for diffusive flux
!> \param[in]     pimpv         Dirichlet value to impose
!> \param[in]     cflv          Local Courant number used to convect
!> \param[in]     hint          Internal exchange coefficient
!_______________________________________________________________________________

subroutine set_convective_outlet_vector_aniso &
 ( coefa , cofaf, coefb , cofbf, pimpv  , cflv  , hint )

!===============================================================================
! Module files
!===============================================================================

!===============================================================================

implicit none

! Arguments

double precision coefa(3), cofaf(3)
double precision coefb(3,3), cofbf(3,3)
double precision pimpv(3), cflv(3)
double precision hint(6)

! Local variables

integer          isou  , jsou

!===============================================================================

do isou = 1, 3

  ! Gradient BCs
  do jsou = 1, 3
    if (jsou.eq.isou) then
      coefb(isou,jsou) = cflv(isou)/(1.d0+cflv(isou))
    else
      coefb(isou,jsou) = 0.d0
    endif
  enddo
  coefa(isou) = (1.d0-coefb(isou,isou))*pimpv(isou)

enddo

! Flux BCs
cofaf(1) = -(hint(1)*coefa(1) + hint(4)*coefa(2) + hint(6)*coefa(3))
cofaf(2) = -(hint(4)*coefa(1) + hint(2)*coefa(2) + hint(5)*coefa(3))
cofaf(3) = -(hint(6)*coefa(1) + hint(5)*coefa(2) + hint(3)*coefa(3))
cofbf(1,1) = hint(1)*(1.d0 - coefb(1,1))
cofbf(2,2) = hint(2)*(1.d0 - coefb(2,2))
cofbf(3,3) = hint(3)*(1.d0 - coefb(3,3))
cofbf(1,2) = hint(4)*(1.d0 - coefb(1,1))
cofbf(2,1) = hint(4)*(1.d0 - coefb(1,1))
cofbf(2,3) = hint(5)*(1.d0 - coefb(2,2))
cofbf(3,2) = hint(5)*(1.d0 - coefb(2,2))
cofbf(1,3) = hint(6)*(1.d0 - coefb(3,3))
cofbf(3,1) = hint(6)*(1.d0 - coefb(3,3))

return
end subroutine


!===============================================================================

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[out]    coefa         explicit BC coefficient for gradients
!> \param[out]    cofaf         explicit BC coefficient for diffusive flux
!> \param[out]    coefb         implicit BC coefficient for gradients
!> \param[out]    cofbf         implicit BC coefficient for diffusive flux
!> \param[in]     pinf          affine part
!> \param[in]     ratio         linear part
!> \param[in]     hint          internal exchange coefficient
!_______________________________________________________________________________

subroutine set_affine_function_scalar &
 ( coefa , cofaf, coefb, cofbf, pinf , ratio, hint  )

!===============================================================================
! Module files
!===============================================================================

!===============================================================================

implicit none

! Arguments

double precision coefa, cofaf, coefb, cofbf, pinf, ratio, hint

!===============================================================================

! Gradient BCs
coefb = ratio
coefa = pinf

! Flux BCs
cofaf = -hint*coefa
cofbf =  hint*(1.d0 - coefb)

return
end subroutine


!===============================================================================

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[out]    coefa         explicit BC coefficient for gradients
!> \param[out]    cofaf         explicit BC coefficient for diffusive flux
!> \param[out]    coefb         implicit BC coefficient for gradients
!> \param[out]    cofbf         implicit BC coefficient for diffusive flux
!> \param[in]     dimp          Flux value to impose
!> \param[in]     hint          Internal exchange coefficient
!_______________________________________________________________________________

subroutine set_neumann_conv_h_neumann_diff_scalar &
 (coefa , cofaf, coefb, cofbf, dimp, hint)

!===============================================================================
! Module files
!===============================================================================

!===============================================================================

implicit none

! Arguments

double precision coefa, cofaf, coefb, cofbf, dimp, hint

!===============================================================================

! Gradient BCs
call set_neumann_scalar(coefa, cofaf, coefb, cofbf, dimp, hint)

! Flux BCs
cofaf = 0.d0
cofbf = 0.d0

return
end subroutine


!===============================================================================

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[out]    coefa         explicit BC coefficient for gradients
!> \param[out]    cofaf         explicit BC coefficient for diffusive flux
!> \param[out]    coefb         implicit BC coefficient for gradients
!> \param[out]    cofbf         implicit BC coefficient for diffusive flux
!> \param[in]     pinf          affine part
!> \param[in]     ratio         linear part
!> \param[in]     dimp          Flux value to impose
!_______________________________________________________________________________

subroutine set_affine_function_conv_neumann_diff_scalar &
 (coefa, cofaf, coefb, cofbf, pinf, ratio, dimp)

!===============================================================================
! Module files
!===============================================================================

!===============================================================================

implicit none

! Arguments

double precision coefa, cofaf, coefb, cofbf, pinf, ratio, dimp

!===============================================================================

! Gradient BCs
coefb = ratio
coefa = pinf

! Flux BCs
cofaf = dimp
cofbf = 0.d0

return
end subroutine


!===============================================================================

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[out]    coefa         explicit BC coefficient for gradients
!> \param[out]    cofaf         explicit BC coefficient for diffusive flux
!> \param[out]    coefb         implicit BC coefficient for gradients
!> \param[out]    cofbf         implicit BC coefficient for diffusive flux
!> \param[in]     hext          convective flux to be imposed
!> \param[in]     dimp          Flux value to impose
!_______________________________________________________________________________

subroutine set_total_flux &
 ( coefa, cofaf, coefb, cofbf, hext, dimp )

!===============================================================================
! Module files
!===============================================================================

!===============================================================================

implicit none

! Arguments

double precision coefa, cofaf, coefb, cofbf, dimp, hext

! Gradients BCs
coefa = 0.d0
coefb = 1.d0

! Flux BCs
cofaf = dimp
cofbf = hext

return
end subroutine


!===============================================================================

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[out]    coefa         explicit BC coefficient for gradients
!> \param[out]    cofaf         explicit BC coefficient for diffusive flux
!> \param[out]    coefb         implicit BC coefficient for gradients
!> \param[out]    cofbf         implicit BC coefficient for diffusive flux
!> \param[in]     pimp          Dirichlet value to impose
!> \param[in]     dimp          Flux value to impose
!_______________________________________________________________________________

subroutine set_dirichlet_conv_neumann_diff_scalar &
 ( coefa, cofaf, coefb, cofbf, pimp, dimp )

!===============================================================================
! Module files
!===============================================================================

!===============================================================================

implicit none

! Arguments

double precision coefa, cofaf, coefb, cofbf, pimp, dimp

! Gradients BCs
coefa = pimp
coefb = 0.d0

! Flux BCs
cofaf = dimp
cofbf = 0.d0

return
end subroutine

!===============================================================================

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[out]    coefa         explicit BC coefficient for gradients
!> \param[out]    cofaf         explicit BC coefficient for diffusive flux
!> \param[out]    coefb         implicit BC coefficient for gradients
!> \param[out]    cofbf         implicit BC coefficient for diffusive flux
!> \param[in]     pimpv         Dirichlet value to impose
!> \param[in]     qimpv         Flux value to impose
!_______________________________________________________________________________

subroutine set_dirichlet_conv_neumann_diff_vector &
 ( coefa, cofaf, coefb, cofbf, pimpv, qimpv )

!===============================================================================
! Module files
!===============================================================================

!===============================================================================

implicit none

! Arguments

double precision coefa(3), cofaf(3)
double precision coefb(3,3), cofbf(3,3)
double precision pimpv(3), qimpv(3)

! Local variables

integer          isou  , jsou

!===============================================================================

do isou = 1, 3

  ! Gradient BCs
  coefa(isou) = pimpv(isou)
  do jsou = 1, 3
    coefb(isou,jsou) = 0.d0
  enddo

  ! Flux BCs
  cofaf(isou) = qimpv(isou)
  do jsou = 1, 3
    cofbf(isou,jsou) = 0.d0
  enddo

enddo

return
end subroutine

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[out]    coefa         explicit BC coefficient for gradients
!> \param[out]    cofaf         explicit BC coefficient for diffusive flux
!> \param[out]    coefb         implicit BC coefficient for gradients
!> \param[out]    cofbf         implicit BC coefficient for diffusive flux
!> \param[in]     pimpts         Dirichlet value to impose
!> \param[in]     qimpts         Flux value to impose
!_______________________________________________________________________________

subroutine set_dirichlet_conv_neumann_diff_tensor &
 ( coefa, cofaf, coefb, cofbf, pimpts, qimpts )

!===============================================================================
! Module files
!===============================================================================

!===============================================================================

implicit none

! Arguments

double precision coefa(6), cofaf(6)
double precision coefb(6,6), cofbf(6,6)
double precision pimpts(6), qimpts(6)

! Local variables

integer          isou  , jsou

!===============================================================================

do isou = 1, 6

  ! BS test sur hextv ? if (abs(hextv(isou)).gt.rinfin*0.5d0) then

  ! Gradient BCs
  coefa(isou) = pimpts(isou)
  do jsou = 1, 6
    coefb(isou,jsou) = 0.d0
  enddo

  ! Flux BCs
  cofaf(isou) = qimpts(isou)
  do jsou = 1, 6
    cofbf(isou,jsou) = 0.d0
  enddo

enddo

return
end subroutine

