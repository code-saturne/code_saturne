!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2015 EDF S.A.
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

!> \file bilscts.f90
!>
!> \brief Wrapper to the function which adds the explicit part of the
!> convection/diffusion
!> terms of a transport equation of a tensor field \f$ \tens{\varia} \f$.
!>
!> More precisely, the right hand side \f$ \vect{Rhs} \f$ is updated as
!> follows:
!> \f[
!> \tens{Rhs} = \tens{Rhs} - \sum_{\fij \in \Facei{\celli}}      \left(
!>        \dot{m}_\ij \left( \tens{\varia}_\fij - \tens{\varia}_\celli \right)
!>      - \mu_\fij \gradt_\fij \tens{\varia} \cdot \tens{S}_\ij  \right)
!> \f]
!>
!> Warning:
!> - \f$ \tens{Rhs} \f$ has already been initialized before calling bilscts!
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
!> \param[in]     coefa       boundary condition array for the variable
!>                              (Explicit part)
!> \param[in]     coefb       boundary condition array for the variable
!>                              (Impplicit part)
!> \param[in]     cofaf       boundary condition array for the diffusion
!>                              of the variable (Explicit part)
!> \param[in]     cofbf       boundary condition array for the diffusion
!>                              of the variable (Implicit part)
!> \param[in]     flumas        mass flux at interior faces
!> \param[in]     flumab        mass flux at boundary faces
!> \param[in]     viscf         \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
!>                               at interior faces for the r.h.s.
!> \param[in]     viscb         \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
!>                               at boundary faces for the r.h.s.
!> \param[in]     viscce        symmetric cell tensor \f$ \tens{\mu}_\celli \f$
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
!> \param[in,out] smbr          right hand side \f$ \vect{Rhs} \f$
!_______________________________________________________________________________

subroutine bilscts &
 ( idtvar , ivar   , iconvp , idiffp , nswrgp , imligp , ircflp , &
   ischcp , isstpp , inc    , imrgra ,                            &
   iwarnp , idftnp , imasac ,                                     &
   blencp , epsrgp , climgp , relaxp , thetap ,                   &
   pvar   , pvara  ,                                              &
   coefa  , coefb  , cofaf  , cofbf  ,                            &
   flumas , flumab , viscf  , viscb  ,   viscce ,                 &
   weighf , weighb ,                                              &
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
integer          inc    , imrgra
integer          idftnp , icvflb , imasac
integer          iwarnp
integer          icvfli(nfabor)

double precision blencp , epsrgp , climgp, relaxp , thetap
double precision pvar  (6  ,ncelet)
double precision pvara (6  ,ncelet)
double precision coefa(6  ,nfabor)
double precision cofaf(6  ,nfabor)
double precision coefb(6,6,nfabor)
double precision cofbf(6,6,nfabor)
double precision weighf(2,nfac), weighb(nfabor)
double precision flumas(nfac)  , flumab(nfabor)
double precision viscf (*)  , viscb (nfabor)
double precision smbr(6,ncelet)
double precision viscce(6,ncelet)
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
endif

! Scalar diffusivity
if (idftnp.eq.1) then

  call bilsc6 &
  !==========
   ( idtvar , f_id   , vcopt  ,                                     &
     icvflb , inc    , ifaccp , imasac ,                            &
     pvar   , pvara  ,                                              &
     itypfb , coefa , coefb , cofaf , cofbf ,                       &
     flumas , flumab , viscf  , viscb  ,                            &
     smbr   )

! Symmetric tensor diffusivity
elseif (idftnp.eq.6) then

  ! Nor diffusive part neither secondary viscosity or transpose of gradient
  idiflc = 0
  vcopt%idiff  = idiflc
  ! Convective part
  if (iconvp.eq.1) then

    call bilsc6 &
    !==========
   ( idtvar , f_id   , vcopt  ,                                     &
     icvflb , inc    , ifaccp , imasac ,                            &
     pvar   , pvara  ,                                              &
     itypfb , coefa , coefb , cofaf , cofbf ,                       &
     flumas , flumab , viscf  , viscb  ,                            &
     smbr   )
  endif


  ! Diffusive part (with a 6x6 symmetric diffusivity)
  if (idiffp.eq.1) then

    call diftnts &
    !==========
   ( idtvar , f_id   , vcopt  ,                                     &
     inc    ,                                                       &
     pvar   , pvara  , itypfb ,                                     &
     coefa , coefb , cofaf , cofbf ,                                &
     viscf  , viscb  , viscce ,                                     &
     weighf , weighb ,                                              &
     smbr  )

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
