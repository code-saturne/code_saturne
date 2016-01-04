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
! ---------

!> \file bilscv.f90
!>
!> \brief Wrapper to the function which adds the explicit part of the
!> convection/diffusion
!> terms of a transport equation of a vector field \f$ \vect{\varia} \f$.
!>
!> More precisely, the right hand side \f$ \vect{Rhs} \f$ is updated as
!> follows:
!> \f[
!> \vect{Rhs} = \vect{Rhs} - \sum_{\fij \in \Facei{\celli}}      \left(
!>        \dot{m}_\ij \left( \vect{\varia}_\fij - \vect{\varia}_\celli \right)
!>      - \mu_\fij \gradt_\fij \vect{\varia} \cdot \vect{S}_\ij  \right)
!> \f]
!>
!> Remark:
!> if ivisep = 1, then we also take \f$ \mu \transpose{\gradt\vect{\varia}}
!> + \lambda \trace{\gradt\vect{\varia}} \f$, where \f$ \lambda \f$ is
!> the secondary viscosity, i.e. usually \f$ -\frac{2}{3} \mu \f$.
!>
!> Warning:
!> - \f$ \vect{Rhs} \f$ has already been initialized before calling bilscv!
!> - mind the sign minus
!>
!> Options for the diffusive scheme:
!> - idftnp = 1: scalar diffusivity
!> - idftnp = 6: symmetric tensor diffusivity
!>
!> Options for the convective scheme:
!> - blencp = 0: upwind scheme for the advection
!> - blencp = 1: no upwind scheme except in the slope test
!> - ischcp = 0: second order
!> - ischcp = 1: centered
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     idtvar        indicator of the temporal scheme
!> \param[in]     ivar          index of the current variable
!> \param[in]     iconvp        indicator
!>                               - 1 convection,
!>                               - 0 otherwise
!> \param[in]     idiffp        indicator
!>                               - 1 diffusion,
!>                               - 0 otherwise
!> \param[in]     nswrgp        number of reconstruction sweeps for the
!>                               gradients
!> \param[in]     imligp        clipping gradient method
!>                               - < 0 no clipping
!>                               - = 0 by neighboring gradients
!>                               - = 1 by the mean gradient
!> \param[in]     ircflp        indicator
!>                               - 1 flux reconstruction,
!>                               - 0 otherwise
!> \param[in]     ischcp        indicator
!>                               - 1 centered
!>                               - 0 2nd order
!> \param[in]     isstpp        indicator
!>                               - 1 without slope test
!>                               - 0 with slope test
!> \param[in]     inc           indicator
!>                               - 0 when solving an increment
!>                               - 1 otherwise
!> \param[in]     imrgra        indicator
!>                               - 0 iterative gradient
!>                               - 1 least squares gradient
!> \param[in]     ivisep        indicator to take \f$ \divv
!>                               \left(\mu \gradt \transpose{\vect{a}} \right)
!>                               -2/3 \grad\left( \mu \dive \vect{a} \right)\f$
!>                               - 1 take into account,
!>                               - 0 otherwise
!> \param[in]     iwarnp        verbosity
!> \param[in]     idftnp        indicator
!>                               - 1 scalar diffusivity
!>                               - 6 symmetric tensor diffusivity
!> \param[in]     imasac        indicator
!>                               - 1 take mass accumulation into account
!>                               - 0 do not take mass accumulation
!> \param[in]     blencp        fraction of upwinding
!> \param[in]     epsrgp        relative precision for the gradient
!>                               reconstruction
!> \param[in]     climgp        clipping coefficient for the computation of
!>                               the gradient
!> \param[in]     relaxp        coefficient of relaxation
!> \param[in]     thetap        weighting coefficient for the theta-schema,
!>                               - thetap = 0: explicit scheme
!>                               - thetap = 0.5: time-centered
!>                               scheme (mix between Crank-Nicolson and
!>                               Adams-Bashforth)
!>                               - thetap = 1: implicit scheme
!> \param[in]     pvar          solved velocity (current time step)
!> \param[in]     pvara         solved velocity (previous time step)
!> \param[in]     coefav        boundary condition array for the variable
!>                               (explicit part)
!> \param[in]     coefbv        boundary condition array for the variable
!>                               (implicit part)
!> \param[in]     cofafv        boundary condition array for the diffusion
!>                               of the variable (explicit part)
!> \param[in]     cofbfv        boundary condition array for the diffusion
!>                               of the variable (implicit part)
!> \param[in]     flumas        mass flux at interior faces
!> \param[in]     flumab        mass flux at boundary faces
!> \param[in]     viscf         \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
!>                               at interior faces for the r.h.s.
!> \param[in]     viscb         \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
!>                               at boundary faces for the r.h.s.
!> \param[in]     secvif        secondary viscosity at interior faces
!> \param[in]     secvib        secondary viscosity at boundary faces
!> \param[in]     icvflb        global indicator of boundary convection flux
!>                               - 0 upwind scheme at all boundary faces
!>                               - 1 imposed flux at some boundary faces
!> \param[in]     icvfli        boundary face indicator array of convection flux
!>                               - 0 upwind scheme
!>                               - 1 imposed flux
!> \param[in,out] smbr          right hand side \f$ \vect{Rhs} \f$
!_______________________________________________________________________________

subroutine bilscv &
 ( idtvar , ivar   , iconvp , idiffp , nswrgp , imligp , ircflp , &
   ischcp , isstpp , inc    , imrgra , ivisep ,                   &
   iwarnp , idftnp , imasac ,                                     &
   blencp , epsrgp , climgp , relaxp , thetap ,                   &
   pvar   , pvara  ,                                              &
   coefav , coefbv , cofafv , cofbfv ,                            &
   flumas , flumab , viscf  , viscb  , secvif , secvib ,          &
   icvflb , icvfli ,                                              &
   smbr   )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar, only: ivarfl
use pointe
use entsor
use parall
use period
use cplsat
use mesh
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer          idtvar
integer          ivar   , iconvp , idiffp , nswrgp , imligp
integer          ircflp , ischcp , isstpp
integer          inc    , imrgra , ivisep
integer          idftnp , icvflb , imasac
integer          iwarnp
integer          icvfli(nfabor)

double precision blencp , epsrgp , climgp, relaxp , thetap
double precision pvar  (3  ,ncelet)
double precision pvara (3  ,ncelet)
double precision coefav(3  ,nfabor)
double precision cofafv(3  ,nfabor)
double precision coefbv(3,3,nfabor)
double precision cofbfv(3,3,nfabor)
double precision flumas(nfac)  , flumab(nfabor)
double precision viscf (*)  , viscb (nfabor)
double precision secvif(nfac), secvib(nfabor)
double precision smbr(3,ncelet)

! Local variables

integer          idiflc, f_id

type(var_cal_opt) vcopt

!===============================================================================

if (ivar.eq.0) then
  f_id = -1
  vcopt%iwarni = iwarnp
  vcopt%iconv  = iconvp
  vcopt%istat  = -1
  vcopt%idiff  = idiffp
  vcopt%idifft = -1
  vcopt%idften = idftnp
  vcopt%iswdyn = -1
  vcopt%ischcv = ischcp
  vcopt%isstpc = isstpp
  vcopt%nswrgr = nswrgp
  vcopt%nswrsm = -1
  vcopt%imrgra = imrgra
  vcopt%imligr = imligp
  vcopt%ircflu = ircflp
  vcopt%iwgrec = 0
  vcopt%thetav = thetap
  vcopt%blencv = blencp
  vcopt%epsilo = -1.d0
  vcopt%epsrsm = -1.d0
  vcopt%epsrgr = epsrgp
  vcopt%climgr = climgp
  vcopt%extrag = -1.d0
  vcopt%relaxv = relaxp
else
  f_id = ivarfl(ivar)
  call field_get_key_struct_var_cal_opt(f_id, vcopt)
  ! Warning: vcopt is a copy of variable calculation options
  ! and could be locally modified here for explicit time stepping.
  vcopt%thetav = thetap
endif

! Scalar diffusivity
if (idftnp.eq.1) then

  call bilsc4 &
  !==========
   ( idtvar , f_id   , vcopt  ,                                     &
     icvflb , inc    , ifaccp , ivisep , imasac ,                   &
     pvar   , pvara  ,                                              &
     itypfb , icvfli , coefav , coefbv , cofafv , cofbfv ,          &
     flumas , flumab , viscf  , viscb  , secvif ,                   &
     smbr   )

! Symmetric tensor diffusivity
elseif (idftnp.eq.6) then

  ! Nor diffusive part neither secondary viscosity or transpose of gradient
  idiflc = 0
  vcopt%idiff  = idiflc
  ! Convective part
  if (iconvp.eq.1) then

    call bilsc4 &
    !==========
   ( idtvar , f_id   , vcopt  ,                                     &
     icvflb , inc    , ifaccp , idiflc , imasac ,                   &
     pvar   , pvara  ,                                              &
     itypfb , icvfli , coefav , coefbv , cofafv , cofbfv ,          &
     flumas , flumab , viscf  , viscb  , secvif ,                   &
     smbr   )

  endif

  ! Diffusive part (with a 3x3 symmetric diffusivity)
  if (idiffp.eq.1) then

    call diftnv &
    !==========
   ( idtvar , f_id   , vcopt  ,                                     &
     inc    , ifaccp , ivisep ,                                     &
     pvar   , pvara  ,                                              &
     itypfb , coefav , coefbv , cofafv , cofbfv ,                   &
     viscf  , viscb  , secvif ,                                     &
     smbr   )

  endif

endif

!--------
! Formats
!--------

!----
! End
!----

return

end subroutine
