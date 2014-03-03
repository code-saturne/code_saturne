!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2014 EDF S.A.
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

!> \file bilsca.f90
!>
!> \brief Wrapper to the function which adds the explicit part of the
!> convection/diffusion
!> terms of a transport equation of a scalar field \f$ \varia \f$.
!>
!> More precisely, the right hand side \f$ Rhs \f$ is updated as
!> follows:
!> \f[
!> Rhs = Rhs - \sum_{\fij \in \Facei{\celli}}      \left(
!>        \dot{m}_\ij \left( \varia_\fij - \varia_\celli \right)
!>      - \mu_\fij \gradv_\fij \varia \cdot \vect{S}_\ij  \right)
!> \f]
!>
!> Warning:
!> - \f$ Rhs \f$ has already been initialized before calling bilsca!
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
!> - ischcp = 1: centred
!> - imucpp = 0: do not multiply the convective part by \f$ C_p \f$
!> - imucpp = 1: multiply the convective part by \f$ C_p \f$
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
!>                               - 0 sinon
!> \param[in]     idiffp        indicator
!>                               - 1 diffusion,
!>                               - 0 sinon
!> \param[in]     nswrgp        number of reconstruction sweeps for the
!>                               gradients
!> \param[in]     imligp        clipping gradient method
!>                               - < 0 no clipping
!>                               - = 0 thank to neighbooring gradients
!>                               - = 1 thank to the mean gradient
!> \param[in]     ircflp        indicator
!>                               - 1 flux reconstruction,
!>                               - 0 otherwise
!> \param[in]     ischcp        indicator
!>                               - 1 centred
!>                               - 0 2nd order
!> \param[in]     isstpp        indicator
!>                               - 1 without slope test
!>                               - 0 with slope test
!> \param[in]     inc           indicator
!>                               - 0 when solving an increment
!>                               - 1 otherwise
!> \param[in]     imrgra        indicator
!>                               - 0 iterative gradient
!>                               - 1 least square gradient
!> \param[in]     iccocg        indicator
!>                               - 1 re-compute cocg matrix (for iterativ gradients)
!>                               - 0 otherwise
!> \param[in]     iwarnp        verbosity
!> \param[in]     imucpp        indicator
!>                               - 0 do not multiply the convectiv term by Cp
!>                               - 1 do multiply the convectiv term by Cp
!> \param[in]     idftnp        indicator
!>                               - 1 scalar diffusivity
!>                               - 6 symmetric tensor diffusivity
!> \param[in]     blencp        fraction of upwinding
!> \param[in]     epsrgp        relative precision for the gradient
!>                               reconstruction
!> \param[in]     climgp        clipping coeffecient for the computation of
!>                               the gradient
!> \param[in]     extrap        coefficient for extrapolation of the gradient
!> \param[in]     relaxp        coefficient of relaxation
!> \param[in]     thetap        weightening coefficient for the theta-schema,
!>                               - thetap = 0: explicit scheme
!>                               - thetap = 0.5: time-centred
!>                               scheme (mix between Crank-Nicolson and
!>                               Adams-Bashforth)
!>                               - thetap = 1: implicit scheme
!> \param[in]     pvar          solved variable (current time step)
!> \param[in]     pvara         solved variable (previous time step)
!> \param[in]     coefap        boundary condition array for the variable
!>                               (Explicit part)
!> \param[in]     coefbp        boundary condition array for the variable
!>                               (Impplicit part)
!> \param[in]     cofafp        boundary condition array for the diffusion
!>                               of the variable (Explicit part)
!> \param[in]     cofbfp        boundary condition array for the diffusion
!>                               of the variable (Implicit part)
!> \param[in]     flumas        mass flux at interior faces
!> \param[in]     flumab        mass flux at boundary faces
!> \param[in]     viscf         \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
!>                               at interior faces for the r.h.s.
!> \param[in]     viscb         \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
!>                               at border faces for the r.h.s.
!> \param[in]     viscce        symmetric cell tensor \f$ \tens{\mu}_\celli \f$
!> \param[in]     xcpp          array of specific heat (Cp)
!> \param[in]     weighf        internal face weight between cells i j in case
!>                               of tensor diffusion
!> \param[in]     weighb        boundary face weight for cells i in case
!>                               of tensor diffusion
!> \param[in]     icvflb        global indicator of boundary convection flux
!>                               - 0 upwind scheme at all boundary faces
!>                               - 1 imposed flux at some boundary faces
!> \param[in]     icvfli        boundary face indicator array of convection flux
!>                               - 0 upwind scheme
!>                               - 1 imposed flux
!> \param[in,out] smbrp         right hand side \f$ \vect{Rhs} \f$
!_______________________________________________________________________________

subroutine bilsca &
 ( idtvar , ivar   , iconvp , idiffp , nswrgp , imligp , ircflp , &
   ischcp , isstpp , inc    , imrgra , iccocg ,                   &
   iwarnp , imucpp , idftnp ,                                     &
   blencp , epsrgp , climgp , extrap , relaxp , thetap ,          &
   pvar   , pvara  , coefap , coefbp , cofafp , cofbfp ,          &
   flumas , flumab , viscf  , viscb  , viscce , xcpp   ,          &
   weighf , weighb ,                                              &
   icvflb , icvfli ,                                              &
   smbrp  )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use pointe
use entsor
use parall
use numvar, only: ivarfl
use period
use cplsat
use mesh
use cs_f_interfaces
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer          idtvar
integer          ivar   , iconvp , idiffp , nswrgp , imligp
integer          ircflp , ischcp , isstpp
integer          inc    , imrgra , iccocg
integer          iwarnp , imucpp, idftnp, icvflb

integer          icvfli(nfabor)

double precision blencp , epsrgp , climgp, extrap, relaxp , thetap

double precision pvar (ncelet), pvara(ncelet)
double precision coefap(nfabor), coefbp(nfabor)
double precision cofafp(nfabor), cofbfp(nfabor)
double precision flumas(nfac), flumab(nfabor)
double precision viscf (nfac), viscb (nfabor)
double precision weighf(2,nfac), weighb(nfabor)
double precision smbrp(ncelet)
double precision xcpp(ncelet)
double precision viscce(*)

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
  vcopt%idften = -1
  vcopt%iswdyn = -1
  vcopt%ischcv = ischcp
  vcopt%isstpc = isstpp
  vcopt%nswrgr = nswrgp
  vcopt%nswrsm = -1
  vcopt%imrgra = imrgra
  vcopt%imligr = imligp
  vcopt%ircflu = ircflp
  vcopt%thetav = thetap
  vcopt%blencv = blencp
  vcopt%epsilo = -1.d0
  vcopt%epsrsm = -1.d0
  vcopt%epsrgr = epsrgp
  vcopt%climgr = climgp
  vcopt%extrag = extrap
  vcopt%relaxv = relaxp
else
  f_id = ivarfl(ivar)
  call field_get_key_struct_var_cal_opt(f_id, vcopt)
endif

! Scalar diffusivity
if (idftnp.eq.1) then
  if (imucpp.eq.0) then

    call bilsc2 &
    !==========
   ( idtvar , f_id   , vcopt  ,                                     &
     icvflb , inc    , iccocg ,                                     &
     ifaccp ,                                                       &
     pvar   , pvara  ,                                              &
     itypfb , icvfli , coefap , coefbp , cofafp , cofbfp ,          &
     flumas , flumab , viscf  , viscb  ,                            &
     smbrp  )


  ! The convective part is mulitplied by Cp for the Temperature
  else

    call bilsct &
    !==========
   ( idtvar , f_id   , vcopt  ,                                     &
     inc    , iccocg ,                                              &
     ifaccp ,                                                       &
     pvar   , pvara  , itypfb , coefap , coefbp , cofafp , cofbfp , &
     flumas , flumab , viscf  , viscb  , xcpp   ,                   &
     smbrp  )

  endif

! Symmetric tensor diffusivity
elseif (idftnp.eq.6) then

  idiflc = 0
  vcopt%idiff  = idiflc
  ! Convective part
  if (imucpp.eq.0.and.iconvp.eq.1) then

    call bilsc2 &
    !==========
   ( idtvar , f_id   , vcopt  ,                                     &
     icvflb , inc    , iccocg ,                                     &
     ifaccp ,                                                       &
     pvar   , pvara  ,                                              &
     itypfb , icvfli , coefap , coefbp , cofafp , cofbfp ,          &
     flumas , flumab , viscf  , viscb  ,                            &
     smbrp  )


  ! The convective part is mulitplied by Cp for the Temperature
  elseif (imucpp.eq.1.and.iconvp.eq.1) then

    call bilsct &
    !==========
   ( idtvar , f_id   , vcopt  ,                                     &
     inc    , iccocg ,                                              &
     ifaccp ,                                                       &
     pvar   , pvara  , itypfb , coefap , coefbp , cofafp , cofbfp , &
     flumas , flumab , viscf  , viscb  , xcpp   ,                   &
     smbrp  )

  endif

  ! Diffusive part
  if (idiffp.eq.1) then

    call diften &
    !==========
   ( idtvar , f_id   , vcopt  ,                                     &
     inc    , iccocg ,                                              &
     pvar   , pvara  , coefap , coefbp , cofafp , cofbfp ,          &
     viscf  , viscb  , viscce ,                                     &
     weighf , weighb ,                                              &
     smbrp  )

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
