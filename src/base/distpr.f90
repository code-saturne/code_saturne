!-------------------------------------------------------------------------------

! This file is part of code_saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2022 EDF S.A.
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

!> \file distpr.f90
!> \brief Compute distance to wall by solving a 3d diffusion equation.
!> Solve
!>   \f[ -\divs ( \grad \varia ) = 1 \f]
!> with:
!>  - \f$ \varia_|b = 0 \f$  at the wall
!>  - \f$ \grad \varia \cdot \vect{n} = 0 \f$ elsewhere
!> The wall distance is then equal to:
!>  \f[
!>  d \simeq -|\grad \varia |
!>  + \sqrt{ \grad \varia \cdot \grad \varia +2 \varia }
!>  \f]
!>
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Arguments
!------------------------------------------------------------------------------
!   mode          name          role
!------------------------------------------------------------------------------
!> \param[in]     itypfb        boundary face types
!> \param[in]     iterns        iteration number on Navier-Stokes equations
!______________________________________________________________________________

subroutine distpr(itypfb, iterns)

!===============================================================================
! Module files
!===============================================================================

use, intrinsic :: iso_c_binding

use paramx
use numvar
use entsor
use optcal
use cstphy
use cstnum
use ppppar
use parall
use period
use mesh
use field
use field_operator
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer          itypfb(nfabor)
integer          iterns

! Local variables

integer          ndircp, imvisp
integer          iel   , ifac
integer          inc   , iccocg, f_id, f_id_pre
integer          mmprpl, nswrsp
integer          imucpp
integer          icvflb, iescap, ircflp
integer          ivoid(1)
integer          counter

double precision dismax, dismin, hint, pimp, qimp, norm_grad
double precision normp

double precision rvoid(1)

double precision, allocatable, dimension(:) :: viscf, viscb
double precision, allocatable, dimension(:) :: dpvar, smbrp, rovsdt
double precision, allocatable, dimension(:,:) :: grad
double precision, allocatable, dimension(:) :: w1
double precision, allocatable, dimension(:) :: imasfl, bmasfl
double precision, pointer, dimension(:)   :: cvar_var
double precision, pointer, dimension(:)   :: cvara_var
double precision, pointer, dimension(:) :: coefap, coefbp
double precision, pointer, dimension(:) :: cofafp, cofbfp
type(var_cal_opt) :: vcopt
type(var_cal_opt), target :: vcopt_loc
type(var_cal_opt), pointer :: p_k_value
type(c_ptr) :: c_k_value

!===============================================================================

!===============================================================================
! 1. Initialization
!===============================================================================

! Allocate temporary arrays for the species resolution
allocate(viscf(nfac), viscb(nfabor))
allocate(imasfl(nfac), bmasfl(nfabor))

do ifac = 1, nfac
  imasfl(ifac) = 0.d0
enddo

do ifac = 1, nfabor
  bmasfl(ifac) = 0.d0
enddo

allocate(dpvar(ncelet), smbrp(ncelet), rovsdt(ncelet))

! Allocate work arrays
allocate(w1(ncelet))

! Initialize variables to avoid compiler warnings
call field_get_id("wall_distance", f_id)

call field_get_key_struct_var_cal_opt(f_id, vcopt)

call field_get_val_s(f_id, cvar_var)
! Previous value is stored in a specific field beacause
! the solved field is not directly the wall distance
call field_get_id_try("wall_distance_aux_pre", f_id_pre)
if (f_id_pre.ge.0) then
  call field_get_val_s(f_id_pre, cvara_var)
else
  call field_get_val_prev_s(f_id, cvara_var)
endif

!===============================================================================
! 2. Boundary conditions
!===============================================================================

!     Conditions aux limites pour le scalaire resolu T
!       Dirichlet a 0 en paroi
!       Neumann nul ailleurs
!     On test aussi la presence d'un Dirichlet

ndircp = 0

call field_get_coefa_s( f_id, coefap)
call field_get_coefb_s( f_id, coefbp)
call field_get_coefaf_s(f_id, cofafp)
call field_get_coefbf_s(f_id, cofbfp)

do ifac = 1, nfabor
  if (itypfb(ifac).eq.iparoi .or. itypfb(ifac).eq.iparug) then

    ! Dirichlet Boundary Condition
    !-----------------------------

    hint = 1.d0/distb(ifac)
    pimp = 0.d0

    call set_dirichlet_scalar &
         !====================
       ( coefap(ifac), cofafp(ifac),             &
         coefbp(ifac), cofbfp(ifac),             &
         pimp        , hint        , rinfin )


    ndircp = ndircp + 1
  else

    ! Neumann Boundary Conditions
    !----------------------------

    hint = 1.d0/distb(ifac)
    qimp = 0.d0

    call set_neumann_scalar &
         !==================
       ( coefap(ifac), cofafp(ifac),             &
         coefbp(ifac), cofbfp(ifac),             &
         qimp        , hint )

  endif
enddo

if (irangp.ge.0) call parcpt(ndircp)

! If no wall initialization to a big value
if (ndircp.eq.0) then
  do iel = 1, ncel
    cvar_var(iel) = grand
  enddo
  return
endif

!===============================================================================
! 3. Prepare system to solve
!===============================================================================

! -- Diagonal

do iel = 1, ncel
  rovsdt(iel) = 0.d0
enddo

! -- Diffusion at faces

do iel = 1, ncel
  w1(iel) = 1.d0
enddo

imvisp = vcopt%imvisf

call viscfa                                                       &
!==========
 ( imvisp ,                                                       &
   w1     ,                                                       &
   viscf  , viscb  )

!===============================================================================
! 4. Solve system
!===============================================================================

! Distance to wall is initialized to 0 for reconstruction

nswrsp = vcopt%nswrsm
ircflp = vcopt%ircflu
iescap = 0
imucpp = 0
! all boundary convective flux with upwind
icvflb = 0
normp = -1.d0

vcopt_loc = vcopt

vcopt_loc%istat  = -1
vcopt_loc%icoupl = -1
vcopt_loc%ndircl = ndircp
vcopt_loc%idifft = -1
vcopt_loc%iwgrec = 0 ! Warning, may be overwritten if a field
vcopt_loc%blend_st = 0 ! Warning, may be overwritten if a field

p_k_value => vcopt_loc
c_k_value = equation_param_from_vcopt(c_loc(p_k_value))

110 continue

do iel = 1, ncelet
  dpvar(iel)  = 0.d0
enddo

! -- RHS

do iel = 1, ncel
  rovsdt(iel) = 0.d0
  smbrp(iel)  = cell_f_vol(iel)
enddo

call cs_equation_iterative_solve_scalar                           &
 ( idtvar , iterns ,                                              &
   f_id   , c_null_char ,                                         &
   iescap , imucpp , normp  , c_k_value       ,                   &
   cvara_var       , cvara_var       ,                            &
   coefap , coefbp , cofafp , cofbfp ,                            &
   imasfl , bmasfl ,                                              &
   viscf  , viscb  , viscf  , viscb  ,                            &
   rvoid  , rvoid  , rvoid  ,                                     &
   icvflb , ivoid  ,                                              &
   rovsdt , smbrp  , cvar_var        , dpvar  ,                   &
   rvoid  , rvoid  )

! Count clippings
mmprpl = 0
dismin =  grand
do iel = 1, ncel
  if (cvar_var(iel).lt.0.d0) then
    mmprpl = mmprpl + 1
    dismin = min(cvar_var(iel), dismin)
    cvar_var(iel) = epzero*(cell_f_vol(iel))**(1.d0/3.d0)
  endif
enddo

if (irangp.ge.0) then
  call parcpt(mmprpl)
  call parmin(dismin)
endif

! Recompute wall distance without reconstruction (that ensure that it is positive)
if (mmprpl.ge.1) then
  if (nswrsp.gt.0) then
    nswrsp = 0
    ircflp = 0
    ! Reset also in var_cal_opt structure because some basic routines
    ! directly use the field options in vcopt...
    vcopt%nswrsm = nswrsp
    vcopt%ircflu = ircflp
    call field_set_key_struct_var_cal_opt(f_id, vcopt)
    write(nfecra,9000) mmprpl

    ! Reset wall distance
    do iel = 1, ncelet
      cvar_var(iel)  = 0.d0
    enddo

    goto 110
  else

    write(nfecra,9001) dismin

  endif
endif

do iel = 1, ncel
  dpvar(iel) = max(cvar_var(iel), 0.d0)
  ! Save working field for the next time step
  if (f_id_pre.ge.0) then
    cvara_var(iel) = cvar_var(iel)
  endif
enddo

if (f_id_pre.ge.0) then
  call synsca(cvara_var)
endif

!===============================================================================
! 5. Compute distance to wall
!===============================================================================

! Allocate a temporary array for the gradient calculation
allocate(grad(3,ncelet))

! Compute current gradient

inc    = 1
iccocg = 1

call field_gradient_scalar(f_id, 0, 0, inc, iccocg, grad)

counter = 0
do iel = 1, ncel
  norm_grad = grad(1,iel)**2.d0+grad(2,iel)**2.d0+grad(3,iel)**2.d0
  if (norm_grad+2.d0*dpvar(iel).ge.0.d0) then
    cvar_var(iel) = sqrt(norm_grad + 2.d0*dpvar(iel)) - sqrt(norm_grad)
  else
    counter = counter + 1
  endif
enddo

if (irangp.ge.0) then
  call parcpt(counter)
endif

if (counter.gt.0) then
  write(nfecra,8000) counter
endif

! Free memory
deallocate(grad)

!===============================================================================
! 6. Compute bounds and print info
!===============================================================================

dismax = -grand
dismin =  grand

do iel = 1, ncel
  dismin = min(cvar_var(iel),dismin)
  dismax = max(cvar_var(iel),dismax)
enddo

if (irangp.ge.0) then
  call parmin(dismin)
  call parmax(dismax)
endif

write(nfecra,1000) dismin, dismax

! Free memory
deallocate(viscf, viscb)
deallocate(dpvar, smbrp, rovsdt)
deallocate(imasfl, bmasfl)
deallocate(w1)

!===============================================================================
! 7. Formats
!===============================================================================

 1000 format(                                                           &
'                                                             ',/,&
' ** WALL DISTANCE                                            ',/,&
'    -------------                                            ',/,&
'                                                             ',/,&
'  Min distance = ',E14.5    ,' Max distance = ',E14.5         ,/)

 8000   format(                                                         &
'@                                                            ',/,&
'@ @@ WARNING: Wall distance calculation                      ',/,&
'@    ========                                                ',/,&
'@  The associated variable does not converge in ',I10,' cells.'/)

 9000   format(                                                         &
'@'                                                            ,/,&
'@ @@ WARNING: Wall distance calculation'                      ,/,&
'@    ========='                                               ,/,&
'@  The laplacian solution does not respect the maximum'       ,/,&
'@  principle in ', i10,' cells. We recompute the laplacien'   ,/,&
'@  without reconstructions.',/)

 9001   format(                                                         &
'@                                                            ',/,&
'@ @@ WARNING: Wall distance calculation                      ',/,&
'@    =========                                               ',/,&
'@  The laplacian solution does not respect the maximum       ',/,&
'@  principle. (laplacian solution is  negative :', E14.6,')    ',/)

return
end subroutine
