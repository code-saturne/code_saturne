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

!> \file distpr.f90
!> \brief Compute distance to wall by solving a 3d diffusion equation.
!> Solve
!>   \f[ \divs ( \grad \varia ) = -1 \f]
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
!______________________________________________________________________________

subroutine distpr(itypfb)

!===============================================================================
! Module files
!===============================================================================

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


! Local variables

integer          ndircp, iconvp, idiffp
integer          iel   , ifac
integer          inc   , iccocg, f_id
integer          mmprpl, nswrsp
integer          imucpp, idftnp
integer          nswrgp
integer          icvflb, iescap, imligp, ircflp, iswdyp, isstpp, ischcp, iwarnp
integer          ivoid(1)

double precision relaxp, blencp, climgp, epsilp, epsrgp, epsrsp, extrap
double precision dismax, dismin, hint, pimp, qimp, norm_grad, thetap

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
call field_get_val_prev_s(f_id, cvara_var)

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

call viscfa                                                       &
!==========
 ( imvisf ,                                                       &
   w1     ,                                                       &
   viscf  , viscb  )

!===============================================================================
! 4. Solve system
!===============================================================================

! Distance to wall is initialized to 0 for reconstruction

iconvp = vcopt%iconv
idiffp = vcopt%idiff
idftnp = vcopt%idften
nswrsp = vcopt%nswrsm
nswrgp = vcopt%nswrgr
imligp = vcopt%imligr
ircflp = vcopt%ircflu
ischcp = vcopt%ischcv
isstpp = vcopt%isstpc
iescap = 0
imucpp = 0
iswdyp = vcopt%iswdyn
iwarnp = vcopt%iwarni
blencp = vcopt%blencv
epsilp = vcopt%epsilo
epsrsp = vcopt%epsrsm
epsrgp = vcopt%epsrgr
climgp = vcopt%climgr
extrap = vcopt%extrag
relaxp = vcopt%relaxv
thetap = vcopt%thetav
! all boundary convective flux with upwind
icvflb = 0

110 continue

do iel = 1, ncelet
  dpvar(iel)  = 0.d0
enddo

! -- RHS

do iel = 1, ncel
  rovsdt(iel) = 0.d0
  smbrp(iel)  = cell_f_vol(iel)
enddo

call codits &
 ( idtvar , f_id   , iconvp , idiffp , ndircp ,                   &
   imrgra , nswrsp , nswrgp , imligp , ircflp ,                   &
   ischcp , isstpp , iescap , imucpp , idftnp , iswdyp ,          &
   iwarnp ,                                                       &
   blencp , epsilp , epsrsp , epsrgp , climgp , extrap ,          &
   relaxp , thetap ,                                              &
   cvara_var       , cvara_var       ,                            &
   coefap , coefbp , cofafp , cofbfp ,                            &
   imasfl , bmasfl ,                                              &
   viscf  , viscb  , viscf  , viscb  , rvoid  ,                   &
   rvoid  , rvoid  ,                                              &
   icvflb , ivoid  ,                                              &
   rovsdt , smbrp  , cvar_var        , dpvar  ,                   &
   rvoid  , rvoid  )

! Count clippings
mmprpl = 0
dismin =  grand
do iel = 1, ncel
  if (cvar_var(iel).lt.0.d0) then
    mmprpl = mmprpl + 1
    dismin = min(cvar_var(iel),dismin)
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

do iel=1,ncel
  dpvar(iel) = max(cvar_var(iel), 0.d0)
enddo

!===============================================================================
! 5. Compute distance to wall
!===============================================================================

! Allocate a temporary array for the gradient calculation
allocate(grad(3,ncelet))

! Compute current gradient

inc    = 1
iccocg = 1

call field_gradient_scalar(f_id, 0, imrgra, inc, iccocg, grad)

do iel = 1, ncel
  norm_grad = grad(1,iel)**2.d0+grad(2,iel)**2.d0+grad(3,iel)**2.d0
  if (norm_grad+2.d0*dpvar(iel).gt.0.d0) then
    cvar_var(iel) = sqrt(norm_grad + 2.d0*dpvar(iel)) - sqrt(norm_grad)
  else
    write(nfecra,8000)iel, xyzcen(1,iel),xyzcen(2,iel),xyzcen(3,iel)
  endif
enddo

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

#if defined(_CS_LANG_FR)

 1000 format(                                                           &
'                                                             ',/,&
' ** DISTANCE A LA PAROI                                      ',/,&
'    -------------------                                      ',/,&
'                                                             ',/,&
'   Distance min = ',E14.5    ,'  Distance max = ',E14.5       ,/)

 8000   format(                                                         &
'@                                                            ',/,&
'@ @@ ATTENTION : Calcul de la distance a la paroi            ',/,&
'@    =========                                               ',/,&
'@  La variable associee ne converge pas a la cellule ',I10    ,/,&
'@       Coord X      Coord Y      Coord Z                    ',/,&
'@ ',3E13.5                                                    ,/)

 9000   format(                                                         &
'@'                                                            ,/,&
'@ @@ ATTENTION : Calcul de la distance a la paroi'            ,/,&
'@    ========='                                               ,/,&
'@  La solution du laplacien ne respecte pas le principe du'   ,/,&
'@  maximum en ', i10, ' cellules. On recalcule le laplacien'  ,/,&
'@  sans les reconstructions.',/)

 9001   format(                                                         &
'@                                                            ',/,&
'@ @@ ATTENTION : Calcul de la distance a la paroi            ',/,&
'@    =========                                               ',/,&
'@  La solution du laplacien ne respecte pas le principe du   ',/,&
'@  maximum. (lapalcien negatif : ', E14.6,')                 ',/)


#else

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
'@  The associated variable does not converge in cell ',I10    ,/,&
'@       Coord X      Coord Y      Coord Z                    ',/,&
'@ ',3E13.5                                                    ,/)

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

#endif

return
end subroutine
