!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2019 EDF S.A.
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
!  Function :
!  --------
!
!> \file cs_coal_physprop.f90
!> \brief Specific physics routine: combustion of pulverized coal
!>        Calculation of \f$ \rho  \f$ of the mixture
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role
!______________________________________________________________________________!
!> \param[in]     mbrom         filling indicator of romb
!> \param[in]     izfppp        area number of the edge face
!>                               for the specific physic module
!______________________________________________________________________________!

subroutine cs_coal_physprop &
 ( mbrom  , izfppp )

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use optcal
use cstphy
use cstnum
use entsor
use parall
use period
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use ppcpfu
use cs_coal_incl
use mesh
use field
use pointe, only:pmapper_double_r1

!===============================================================================

implicit none

! Arguments

integer          mbrom
integer          izfppp(nfabor)

! Local variables

integer          iel, icha, icla
integer          izone, ifac
integer          ioxy , nbclip1,nbclip2
integer          iscdri, keydri, iflid, nfld, keyccl
integer          f_id
integer          iok1,iok2,iok3

character(len=80) :: fname, name

double precision x1sro1, x2sro2, srrom1, uns1pw
double precision x2tot, wmolme, unsro1
double precision ff3min,ff3max,valmin,valmax

integer          ipass
data             ipass /0/
save             ipass

double precision, dimension (:), allocatable :: f1m,f2m,f3m,f4m,f5m
double precision, dimension (:), allocatable :: f6m,f7m,f8m,f9m
double precision, dimension (:), allocatable :: enth1 , fvp2m
double precision, dimension (:), allocatable :: xoxyd,enthox
double precision, dimension(:), pointer :: cpro_taup
double precision, dimension(:), pointer :: cpro_x1
double precision, dimension(:), pointer :: cpro_x2
double precision, dimension(:), pointer :: cpro_rom2, cpro_diam2
double precision, dimension(:), pointer :: brom, crom
double precision, dimension(:,:), pointer :: cvar_vel
double precision, dimension(:,:), pointer :: vg_lim_pi
double precision, dimension(:,:), pointer :: vdc
double precision, dimension(:), pointer :: v_x_pi,v_y_pi,v_z_pi
double precision, dimension(:,:), pointer :: vdp_i
double precision, dimension(:,:), pointer :: vg_pi
double precision, dimension(:), pointer :: cvar_f1m, cvar_f2m
double precision, dimension(:), pointer :: cvar_f4m, cvar_f5m, cvar_f6m
double precision, dimension(:), pointer :: cvar_f7m, cvar_f8m, cvar_f9m
double precision, dimension(:), pointer :: cvar_fvp2m
double precision, dimension(:), pointer :: cvar_hox, cvar_hgas
double precision, dimension(:), pointer :: cpro_rom1
type(pmapper_double_r1), dimension(:) , allocatable :: cpro_x2b, cpro_ro2

!===============================================================================
!
!===============================================================================
! 0. Counting the calls
!===============================================================================

ipass = ipass + 1

!===============================================================================
! 1. Initializations
!===============================================================================

! Massic fraction of gas
call field_get_val_s_by_name("x_c", cpro_x1)

! Number of fields
call field_get_n_fields(nfld)

if ( noxyd .ge. 2 ) then
  call field_get_val_s(ivarfl(isca(if4m)), cvar_f4m)
  if ( noxyd .eq. 3 ) then
    call field_get_val_s(ivarfl(isca(if5m)), cvar_f5m)
  endif
endif
if ( ippmod(iccoal) .ge. 1 ) then
  call field_get_val_s(ivarfl(isca(if6m)), cvar_f6m)
endif
call field_get_val_s(ivarfl(isca(if7m)), cvar_f7m)
if ( ihtco2 .eq. 1 ) then
  call field_get_val_s(ivarfl(isca(if8m)), cvar_f8m)
endif
if ( ihth2o .eq. 1 ) then
  call field_get_val_s(ivarfl(isca(if9m)), cvar_f9m)
endif
call field_get_val_s(ivarfl(isca(ifvp2m)), cvar_fvp2m)
if ( ieqnox .eq. 1 ) then
  call field_get_val_s(ivarfl(isca(ihox)), cvar_hox)
endif
call field_get_val_s(ivarfl(isca(ihgas)), cvar_hgas)

!===============================================================================
! Deallocation dynamic arrays
!----
allocate(f1m(1:ncelet), f2m(1:ncelet), f3m(1:ncelet), STAT=iok1)
allocate(f4m(1:ncelet), f5m(1:ncelet), STAT=iok1)
allocate(f6m(1:ncelet), f7m(1:ncelet), f8m(1:ncelet), f9m(1:ncelet), STAT=iok2)
allocate(enth1(1:ncel), fvp2m(1:ncel), STAT=iok3)
allocate(cpro_x2b(nclacp), cpro_ro2(nclacp))
!----
if ( iok1 > 0 .or. iok2 > 0 .or. iok3 > 0) then
  write(nfecra,*) ' Memory allocation error inside: '
  write(nfecra,*) '     cs_coal_physprop             '
  call csexit(1)
endif
if (ieqnox .eq. 1) then

  allocate(xoxyd(1:ncelet),enthox(1:ncelet),STAT=iok1)

  if (iok1 > 0) then
    write(nfecra,*) ' Memory allocation error inside:         '
    write(nfecra,*) '   cs_coal_physprop for xoxyd and enthox '
  endif
endif

!===============================================================================
! 2. Calculation of the physical properties of the dispersed phase
!                    cell values
!                    -----------
!    Mass fraction of solid
!    Diameter
!    Mass density
!===============================================================================

call cs_coal_physprop2 ( ncelet , ncel )
!=====================

!===============================================================================
! 3. Calculation of the physical properties of the gaseous phase
!                    cell values
!                    -----------
!    Temperature
!    Mass density
!    Concentrations of the gaseous species
!===============================================================================

! --- Calculation of the gas enthalpy  enth1
!        of F1M
!        of F2M
!        of F3M                    in W3=1-F1M-F2M-F4M-F5M-F6M-F7M-F8M-F9M
!        of F4M
!        of F5M
!        of F6M
!        of F7M
!        of F8M
!        of F9M
!        of FVP2M
!
! Initialization of fm and of x2 at 0
f1m( : ) = 0.d0
f2m( : ) = 0.d0
f3m( : ) = 0.d0
f4m( : ) = 0.d0
f5m( : ) = 0.d0
f6m( : ) = 0.d0
f7m( : ) = 0.d0
f8m( : ) = 0.d0
f9m( : ) = 0.d0

do icla = 1, nclacp
  call field_get_val_s(ix2(icla),cpro_x2b(icla)%p)
  call field_get_val_s(irom2(icla),cpro_ro2(icla)%p)
enddo

do iel = 1, ncel
  cpro_x1(iel) = 1.d0
  do icla = 1, nclacp
    cpro_x1(iel) = cpro_x1(iel) - cpro_x2b(icla)%p(iel)
  enddo
enddo

! ---> Handle parallelism and periodicity
!      (periodicity of rotation is not ensured here)
if (irangp.ge.0 .or. iperio.eq.1) then
  call synsca(cpro_x1)
endif

do icha = 1, ncharb
  call field_get_val_s(ivarfl(isca(if1m(icha))), cvar_f1m)
  call field_get_val_s(ivarfl(isca(if2m(icha))), cvar_f2m)
  do iel = 1, ncel
    f1m(iel) =  f1m(iel) + cvar_f1m(iel)
    f2m(iel) =  f2m(iel) + cvar_f2m(iel)
  enddo
enddo

if (ieqnox .eq. 1) then
  do iel = 1, ncel
    xoxyd(iel)= cpro_x1(iel)-f1m(iel)-f2m(iel)
  enddo
endif

ff3min = 1.d+20
ff3max =-1.d+20
nbclip1= 0
nbclip2= 0
valmin = 1.d+20
valmax =-1.d+20
do iel = 1, ncel
  uns1pw = 1.d0/cpro_x1(iel)
  if ( noxyd .ge. 2 ) then
    f4m(iel) = cvar_f4m(iel)
    if ( noxyd .eq. 3 ) then
      f5m(iel) = cvar_f5m(iel)
    endif
  endif

  if ( ippmod(iccoal) .ge. 1 ) then
    f6m(iel) = cvar_f6m(iel)
  endif

  f7m(iel) =  cvar_f7m(iel)

  if ( ihtco2 .eq. 1 ) then
    f8m(iel) = cvar_f8m(iel)
  endif

  if ( ihth2o .eq. 1 ) then
    f9m(iel) = cvar_f9m(iel)
  endif

  fvp2m(iel) = cvar_fvp2m(iel)

  ! Units: [kg scalars / kg gas]
  f1m(iel)  = f1m(iel)    *uns1pw
  f2m(iel)  = f2m(iel)    *uns1pw
  f4m(iel)  = f4m(iel)    *uns1pw
  f5m(iel)  = f5m(iel)    *uns1pw
  f6m(iel)  = f6m(iel)    *uns1pw
  f7m(iel)  = f7m(iel)    *uns1pw
  f8m(iel)  = f8m(iel)    *uns1pw
  f9m(iel)  = f9m(iel)    *uns1pw

  fvp2m(iel)= fvp2m(iel)  *uns1pw

  f3m(iel) = 1.d0                                        &
           -( f1m(iel)+f2m(iel)+f4m(iel)+f5m(iel)        &
             +f6m(iel)+f7m(iel)+f8m(iel)+f9m(iel) )

  ff3max = max(ff3max,f3m(iel))
  ff3min = min(ff3min,f3m(iel))

  if ( ieqnox .eq. 1 ) then
    enthox(iel) = cvar_hox(iel)/xoxyd(iel)
  endif

enddo

if (irangp .ge. 0) then
  call parmin(ff3min)
  call parmax(ff3max)
  call parcpt(nbclip1)
  call parcpt(nbclip2)
  call parmin(valmin)
  call parmax(valmax)
endif

write(nfecra,*) ' Values of F3 min and max: ',ff3min,ff3max
if (nbclip1 .gt. 0) then
  write(nfecra,*) ' Clipping phase gas variance in min:',nbclip1,valmin
endif
if (nbclip2 .gt. 0) then
  write(nfecra,*) ' Clipping phase gas variance in max:',nbclip2,valmax
endif

! ---- Gas Enthalpy h1 (cpro_x1 h1 is transported)
do iel = 1, ncel
  enth1(iel) = cvar_hgas(iel) / cpro_x1(iel)
enddo

call field_get_val_s(irom1, cpro_rom1)

call cs_coal_physprop1 &
!=====================
 ( ncelet , ncel   ,                                      &
   f1m    , f2m    , f3m    , f4m    , f5m    ,           &
   f6m    , f7m    , f8m    , f9m    , fvp2m  ,           &
   enth1  , enthox ,                                      &
   cpro_rom1 )

!===============================================================================
! 4. Calculation of the physical properties of the dispersed phase
!                    cells value
!                    -----------
!    Temperature
!===============================================================================

! --- Transport of H2

call  cs_coal_thfieldconv2 ( ncelet , ncel )
!=========================

!===============================================================================
! 5. Calculation of the physical properties of the mixture
!                    cells value
!                    -----------
!    Mass density
!===============================================================================
! --- Calculation of Rho of the mixture: 1/Rho = X1/Rho1 + Sum(X2/Rho2)
!     We relax when we have a rho n available, ie
!     from the second passage or
!     from the first passage if we are in continuation of the calculation and
!     that we have reread the mass density in the file suite.

call field_get_val_s(icrom, crom)

if (ipass.gt.1.or.(isuite.eq.1.and.initro.eq.1)) then
  srrom1 = srrom
else
  srrom1 = 0.d0
endif

do iel = 1, ncel
  x2sro2 = 0.d0
  do icla = 1, nclacp
    x2sro2 = x2sro2 + cpro_x2b(icla)%p(iel) / cpro_ro2(icla)%p(iel)
  enddo
  x1sro1 = cpro_x1(iel) / cpro_rom1(iel)
  ! ---- Eventual relaxation to give in ppini1.f90
  crom(iel) = srrom1*crom(iel)                  &
            + (1.d0-srrom1)/(x1sro1+x2sro2)
enddo


!===============================================================================
! 6. Calculation of the density of the mixture
!                      face values
!                      -----------
!===============================================================================

mbrom = 1
call field_get_val_s(ibrom, brom)
call field_get_val_s(icrom, crom)

! ---> Mass density on edges for all faces
!      The input faces are recalculated.

do ifac = 1, nfabor
  iel = ifabor(ifac)
  brom(ifac) = crom(iel)
enddo

! ---> Mass density on edge for all ONLY inlet faces
!      The test on izone is used for the calculation

if ( ipass.gt.1 .or. isuite.eq.1 ) then
  do ifac = 1, nfabor

    izone = izfppp(ifac)
    if(izone.gt.0) then
      if ( ientat(izone).eq.1 .or. ientcp(izone).eq.1 ) then
        x2sro2 = 0.d0
        x2tot  = 0.d0
        do icla = 1, nclacp
          x2sro2 = x2sro2 + x20(izone,icla)/rho20(icla)
          x2tot  = x2tot  + x20(izone,icla)
        enddo

        ioxy = inmoxy(izone)
        wmolme =( oxyo2(ioxy)+oxyn2(ioxy)                         &
                 +oxyh2o(ioxy)+oxyco2(ioxy))                      &
               /( wmole(io2) *oxyo2(ioxy)                         &
                 +wmole(in2) *oxyn2(ioxy)                         &
                 +wmole(ih2o)*oxyh2o(ioxy)                        &
                 +wmole(ico2)*oxyco2(ioxy) )

        unsro1 = (wmolme*cs_physical_constants_r*timpat(izone)) / p0
        x1sro1 = (1.d0-x2tot) * unsro1
        brom(ifac) = 1.d0 / (x1sro1+x2sro2)
      endif
    endif

  enddo
endif

!===============================================================================
! 7. Compute the drift velocity if needed
!===============================================================================
if (i_comb_drift.ge.1) then

  ! Get all needed fields
  call field_get_val_v(ivarfl(iu), cvar_vel)

  ! Key id for drift scalar
  call field_get_key_id("drift_scalar_model", keydri)

  ! Key id of the coal scalar class
  call field_get_key_id("scalar_class", keyccl)

  ! 1. Compute the limit velocity
  !------------------------------


  ! Loop over coal particle classes
  ! We only handle here coal class with a drift
  !--------------------------------------------

  do iflid = 0, nfld-1

    ! Index of the scalar class (<0 if the scalar belongs to the gas phase)
    call field_get_key_int(iflid, keyccl, icla)

    call field_get_key_int(iflid, keydri, iscdri)

    ! We only handle here one scalar with a drift per particle class
    if (icla.ge.1.and.btest(iscdri, DRIFT_SCALAR_ADD_DRIFT_FLUX)) then

      call field_get_val_s(irom2(icla), cpro_rom2)
      call field_get_val_s(idiam2(icla), cpro_diam2)
      call field_get_val_s(ix2(icla), cpro_x2)

      ! Position of variables, coefficients
      ! -----------------------------------

      ! Name of the drift scalar
      call field_get_name(iflid, fname)

      ! Index of the corresponding relaxation time (cpro_taup)
      call field_get_id('drift_tau_'//trim(fname), f_id)
      call field_get_val_s(f_id, cpro_taup)

      write(name,'(a,i2.2)')'vg_lim_p_' ,icla
      call field_get_val_v_by_name(name, vg_lim_pi)

      do iel = 1, ncel
        vg_lim_pi(1, iel) = cpro_taup(iel)*gx
        vg_lim_pi(2, iel) = cpro_taup(iel)*gy
        vg_lim_pi(3, iel) = cpro_taup(iel)*gz
      enddo

    endif ! test icla

  enddo ! loop on iflid


  ! 2. Init of the drift velocity of the continuous phase (gas)
  !------------------------------------------------------------
  call field_get_val_v_by_name('vd_c', vdc)
  call field_get_val_s_by_name('x_c', cpro_x1)
  do iel = 1, ncel
    vdc(1, iel) = 0.d0
    vdc(2, iel) = 0.d0
    vdc(3, iel) = 0.d0
  enddo

endif

! 3. Transported particle velocity
!---------------------------------
if (i_comb_drift.eq.1) then

  do icla = 1, nclacp

    write(name,'(a,i2.2)')'v_x_p_' ,icla
    call field_get_val_s_by_name(name, v_x_pi)

    write(name,'(a,i2.2)')'v_y_p_' ,icla
    call field_get_val_s_by_name(name, v_y_pi)

    write(name,'(a,i2.2)')'v_z_p_' ,icla
    call field_get_val_s_by_name(name, v_z_pi)

    write(name,'(a,i2.2)')'vd_p_' ,icla
    call field_get_val_v_by_name(name, vdp_i)

    call field_get_val_s(ix2(icla), cpro_x2)
    do iel = 1,ncel
      ! Vdi = Vpi-Vs
      if (cpro_x2(iel).gt. 1.d-7 ) then
        vdp_i(1, iel) = v_x_pi(iel)-cvar_vel(1,iel)
        vdp_i(2, iel) = v_y_pi(iel)-cvar_vel(2,iel)
        vdp_i(3, iel) = v_z_pi(iel)-cvar_vel(3,iel)
      else
        vdp_i(1, iel) = 0.d0
        vdp_i(2, iel) = 0.d0
        vdp_i(3, iel) = 0.d0
      endif
      vdc(1, iel) = vdc(1, iel) - cpro_x2(iel)*vdp_i(1, iel)
      vdc(2, iel) = vdc(2, iel) - cpro_x2(iel)*vdp_i(2, iel)
      vdc(3, iel) = vdc(3, iel) - cpro_x2(iel)*vdp_i(3, iel)
    enddo

  enddo

  do iel = 1, ncel
    vdc(1, iel) = vdc(1, iel)/cpro_x1(iel)
    vdc(2, iel) = vdc(2, iel)/cpro_x1(iel)
    vdc(3, iel) = vdc(3, iel)/cpro_x1(iel)
  enddo

  do icla = 1 , nclacp

    write(name,'(a,i2.2)')'vd_p_' ,icla
    call field_get_val_v_by_name(name, vdp_i)

    write(name,'(a,i2.2)')'vg_p_' ,icla
    call field_get_val_v_by_name(name, vg_pi)

    do iel = 1, ncel
     vg_pi(1, iel) = vdp_i(1, iel) - vdc(1, iel)
     vg_pi(2, iel) = vdp_i(2, iel) - vdc(2, iel)
     vg_pi(3, iel) = vdp_i(3, iel) - vdc(3, iel)
    enddo

  enddo

! Prescribed drift
!-----------------
elseif (i_comb_drift.gt.1) then

  do icla = 1, nclacp

    write(name,'(a,i2.2)')'vg_lim_p_' ,icla
    call field_get_val_v_by_name(name, vg_lim_pi)

    write(name,'(a,i2.2)')'vg_p_' ,icla
    call field_get_val_v_by_name(name, vg_pi)

    call field_get_val_s(ix2(icla), cpro_x2)

    do iel = 1, ncel

      ! FIXME vg_ is useless!
      vg_pi(1, iel) = vg_lim_pi(1, iel)
      vg_pi(2, iel) = vg_lim_pi(2, iel)
      vg_pi(3, iel) = vg_lim_pi(3, iel)
      vdc(1, iel) = vdc(1, iel) -cpro_x2(iel) * vg_pi(1, iel)
      vdc(2, iel) = vdc(2, iel) -cpro_x2(iel) * vg_pi(2, iel)
      vdc(3, iel) = vdc(3, iel) -cpro_x2(iel) * vg_pi(3, iel)
    enddo
  enddo

  do icla = 1, nclacp
    write(name,'(a,i2.2)')'vg_p_' ,icla
    call field_get_val_v_by_name(name, vg_pi)

    write(name,'(a,i2.2)')'vd_p_' ,icla
    call field_get_val_v_by_name(name, vdp_i)

    do iel = 1, ncel

      vdp_i(1, iel) = vdc(1, iel) + vg_pi(1, iel)
      vdp_i(2, iel) = vdc(2, iel) + vg_pi(2, iel)
      vdp_i(3, iel) = vdc(3, iel) + vg_pi(3, iel)

    enddo

  enddo
endif

!--------
! Formats
!--------

!===============================================================================
! Deallocation dynamic arrays

deallocate(f1m,f2m,f3m,f4m,f5m,STAT=iok1)
deallocate(f6m,f7m,f8m,f9m,    STAT=iok2)
deallocate(enth1,fvp2m,     STAT=iok3)
deallocate(cpro_ro2, cpro_x2b)

if (iok1 > 0 .or. iok2 > 0 .or. iok3 > 0) then
  write(nfecra,*) ' Memory deallocation error inside: '
  write(nfecra,*) '     cs_coal_physprop              '
  call csexit(1)
endif

if (ieqnox .eq. 1) then

  deallocate(xoxyd,enthox)

  if (iok1 > 0) then
    write(nfecra,*) ' Memory deallocation error inside:       '
    write(nfecra,*) '   cs_coal_physprop for xoxyd and enthox '
    call csexit(1)
  endif
endif

!----
! End
!----
return
end subroutine
