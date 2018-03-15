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
! Function:
! ---------

!> \file alelav.f90
!>
!> \brief This subroutine performs the solving of a Poisson equation
!> on the mesh velocity for ALE module. It also updates the mesh displacement
!> so that it can be used to update mass fluxes (due to mesh displacement).
!>
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     iterns        Navier-Stokes iteration number
!_______________________________________________________________________________

subroutine alelav &
 ( iterns )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens
use numvar
use entsor
use optcal
use cstnum
use cstphy
use pointe
use albase, only: ialtyb, iortvm, impale, fdiale
use parall
use period
use mesh
use field
use field_operator
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer          iterns

! Local variables

character(len=80) :: chaine
integer          iel   , isou  , jsou  , ifac, inod
integer          iflmas, iflmab
integer          nswrgp, imligp, iwarnp
integer          iconvp, idiffp, ndircp
integer          nswrsp, ircflp, ischcp, isstpp, iescap
integer          ivisep, pot_f_id
integer          iswdyp, idftnp, icvflb
integer          ivoid(1)
integer          inc, iprev

double precision blencp, epsilp, epsrgp, climgp, thetv
double precision epsrsp, prosrf
double precision relaxp
double precision hint, distbf, srfbn2
double precision rinfiv(3), pimpv(3)
double precision rvoid(1)

double precision, allocatable, dimension(:,:) :: dproj
double precision, allocatable, dimension(:,:,:) :: gradm
double precision, allocatable, dimension(:) :: viscf, viscb
double precision, allocatable, dimension(:,:) :: smbr
double precision, allocatable, dimension(:,:,:) :: fimp

double precision, dimension(:), pointer :: imasfl, bmasfl
double precision, dimension(:), pointer :: brom, cpro_vism_s
double precision, dimension(:,:), pointer :: cfaale, claale
double precision, dimension(:,:,:), pointer :: cfbale, clbale
double precision, dimension(:,:), pointer :: mshvel, mshvela, cpro_vism_v
double precision, dimension(:,:), pointer :: disale, disala
double precision, pointer, dimension(:)   :: dt

type(var_cal_opt)  vcopt

!===============================================================================

!===============================================================================
! 1. Initialization
!===============================================================================

! Allocate temporary arrays for the radiative equations resolution
allocate(viscf(nfac), viscb(nfabor))
allocate(smbr(3,ncelet))
allocate(fimp(3,3,ncelet))

rinfiv(1) = rinfin
rinfiv(2) = rinfin
rinfiv(3) = rinfin

call field_get_val_s_by_name('dt', dt)

if (iortvm.eq.0) then
  call field_get_val_s(ivisma, cpro_vism_s)
else
  call field_get_val_v(ivisma, cpro_vism_v)
endif

! The mass flux is necessary to call coditv but not used (ICONV=0)
! Except for the free surface, where it is used as a Boundary condition
call field_get_key_int(ivarfl(iu), kimasf, iflmas)
call field_get_key_int(ivarfl(iu), kbmasf, iflmab)
call field_get_val_s(iflmas, imasfl)
call field_get_val_s(iflmab, bmasfl)

call field_get_val_v(ivarfl(iuma), mshvel)
call field_get_val_prev_v(ivarfl(iuma), mshvela)

call field_get_key_struct_var_cal_opt(ivarfl(iuma), vcopt)

if(vcopt%iwarni.ge.1) then
  write(nfecra,1000)
endif

call field_get_coefa_v(ivarfl(iuma), claale)
call field_get_coefb_v(ivarfl(iuma), clbale)
call field_get_coefaf_v(ivarfl(iuma), cfaale)
call field_get_coefbf_v(ivarfl(iuma), cfbale)

! We compute the boundary condition on the mesh velocity at the free surface
! from the new mass flux.

! Density at the boundary
call field_get_val_s(ibrom, brom)

! The mesh move in the direction of the gravity in case of free-surface
do ifac = 1, nfabor
  if (ialtyb(ifac) .eq. ifresf) then

    iel = ifabor(ifac)
    distbf = distb(ifac)
    srfbn2 = surfbn(ifac)**2
    if (iortvm.eq.0) then
      hint = cpro_vism_s(iel)/distbf
    else !FIXME
      hint = ( cpro_vism_v(1,iel)*surfbo(1,ifac)**2    &
             + cpro_vism_v(2,iel)*surfbo(2,ifac)**2    &
             + cpro_vism_v(3,iel)*surfbo(3,ifac)**2 )  &
           /distbf/srfbn2
    endif

    prosrf = gx*surfbo(1,ifac) + gy*surfbo(2,ifac) + gz*surfbo(3,ifac)
    pimpv(1) = gx*bmasfl(ifac)/(brom(ifac)*prosrf)
    pimpv(2) = gy*bmasfl(ifac)/(brom(ifac)*prosrf)
    pimpv(3) = gz*bmasfl(ifac)/(brom(ifac)*prosrf)

    call set_dirichlet_vector &
         !====================
       ( claale(:,ifac)  , cfaale(:,ifac)  ,             &
         clbale(:,:,ifac), cfbale(:,:,ifac),             &
         pimpv           , hint            , rinfiv )
  endif
enddo

!===============================================================================
! 2. Solving of the mesh velocity equation
!===============================================================================

if (vcopt%iwarni.ge.1) then
  call field_get_name(ivarfl(iuma), chaine)
  write(nfecra,1100) chaine(1:16)
endif

do iel = 1, ncelet
  do isou = 1, 3
    smbr(isou,iel)   = 0.d0
    do jsou = 1, 3
      fimp(isou,jsou,iel) = 0.d0
    enddo
  enddo
enddo

if (iortvm.eq.0) then
  call viscfa &
  !==========
( imvisf ,                                                       &
  cpro_vism_s ,                                                  &
  viscf  , viscb  )
else
  call visort &
  !==========
( imvisf ,                                                       &
  cpro_vism_v(1,:), cpro_vism_v(2,:), cpro_vism_v(3,:),          &
  viscf  , viscb  )
endif

iconvp = vcopt%iconv
idiffp = vcopt%idiff
ndircp = ndircl(iuma)
nswrsp = vcopt%nswrsm
nswrgp = vcopt%nswrgr
imligp = vcopt%imligr
ircflp = vcopt%ircflu
ischcp = vcopt%ischcv
isstpp = vcopt%isstpc
iescap = 0
idftnp = vcopt%idften
iswdyp = vcopt%iswdyn
iwarnp = vcopt%iwarni
blencp = vcopt%blencv
epsilp = vcopt%epsilo
epsrsp = vcopt%epsrsm
epsrgp = vcopt%epsrgr
climgp = vcopt%climgr
relaxp = 1.d0
thetv  = 1.d0
! all boundary convective flux with upwind
icvflb = 0

! we do not take into account the transpose of grad
ivisep = 0

!we do not add gradP to the RHS in coditv
pot_f_id = -1

call coditv &
!==========
 ( idtvar , iterns , ivarfl(iuma)    , iconvp , idiffp , ndircp , &
   imrgra , nswrsp , nswrgp , imligp , ircflp , ivisep ,          &
   ischcp , isstpp , iescap , idftnp , iswdyp ,                   &
   iwarnp ,                                                       &
   blencp , epsilp , epsrsp , epsrgp , climgp ,                   &
   relaxp , thetv  ,                                              &
   mshvela         , mshvela         ,                            &
   claale , clbale , cfaale , cfbale , pot_f_id ,                 &
   imasfl , bmasfl ,                                              &
   viscf  , viscb  , viscf  , viscb  , viscf  , viscb  ,          &
   rvoid  , rvoid  , rvoid  ,                                     &!FIXME do a proper anisotropic version
   icvflb , ivoid  ,                                              &
   fimp   ,                                                       &
   smbr   ,                                                       &
   mshvel ,                                                       &
   rvoid  )

! Free memory
deallocate(viscf, viscb)
deallocate(smbr, fimp)

!===============================================================================
! 3. Update nodes displacement
!===============================================================================

call field_get_val_v(fdiale, disale)
call field_get_val_prev_v(fdiale, disala)

! Allocate a temporary array
allocate(dproj(3,nnod))
allocate(gradm(3,3,ncelet))

inc = 1
iprev = 0

call field_gradient_vector(ivarfl(iuma), iprev, imrgra, inc,      &
                           gradm)

call field_get_coefa_v(ivarfl(iuma), claale)
call field_get_coefb_v(ivarfl(iuma), clbale)

call aledis &
 ( ialtyb ,                                                       &
   mshvel , gradm  ,                                              &
   claale , clbale ,                                              &
   dt     , dproj  )

!FIXME warning if nterup > 1, use itrale ?
! Update mesh displacement only where it is not imposed by the user (ie when impale <> 1)
do inod = 1, nnod
  if (impale(inod).eq.0) then
    do isou = 1, 3
      disale(isou,inod) = disala(isou,inod) + dproj(isou,inod)
    enddo
  endif
enddo

! Free memory
deallocate(dproj)
deallocate(gradm)

!--------
! Formats
!--------

#if defined(_CS_LANG_FR)

 1000    format(/,                                                &
'   ** RESOLUTION DE LA VITESSE DE MAILLAGE'                   ,/,&
'      ------------------------------------'                   ,/)
 1100    format(/,'           RESOLUTION POUR LA VARIABLE ', a16,/)

#else

 1000    format(/,                                                &
'   ** SOLVING MESH VELOCITY'                                  ,/,&
'      ---------------------'                                  ,/)
 1100    format(/,'           SOLVING VARIABLE ', a16          ,/)

#endif

!----
! End
!----

return

end subroutine
