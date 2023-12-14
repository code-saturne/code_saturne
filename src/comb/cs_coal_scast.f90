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
! Function:
! ---------
!> \file cs_coal_scast.f90
!> \brief Specific physic routine: pulverized coal flame
!>   Souce terms have to be precised for a scalar PP
!>   on a step of time
!>
!> \warning  the treatement of source terms is different
!> -------   from that of cs_user_source_terms.
!>
!> We solve: \f[ rovsdt D(var) = smbrs \f]
!>
!> rovsdt and smbrs already contain eventual user source terms.
!> So they have to be incremented and not erased.
!>
!> For stability reasons, only positive terms can be added in rovsdt.
!> There is no contraint for smbrs.
!>
!> In the case of a source term in \f$ cexp + cimp var \f$, it has to be written:
!>        - \f$ smbrs  = smbrs  + cexp + cimp var \f$
!>        - \f$ rovsdt = rovsdt + \max(-cimp,0) \f$
!>
!> Here are \f$ rovsdt \f$ and \f$ smbrs \f$ (they contain \f$ \rho volume\f$)
!>    smbrs in kg variable/s :
!>     \c i.e.: - for velocity            \f$ kg . m . s^{-2} \f$
!>              - for temperature         \f$ kg . [degres] . s^{-1} \f$
!>              - for enthalpy            \f$ J . s^{-1} \f$
!>              - rovsdt                  \f$ kg . s^{-1} \f$
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!                ARGUMENTS
!______________________________________________________________________________.
!  mode           name          role
!______________________________________________________________________________!
!> \param[in]     iscal         scalar number
!> \param[in,out] smbrs         explicit second member
!> \param[in,out] rovsdt        implicit diagonal part
!______________________________________________________________________________!

subroutine cs_coal_scast &
 (iscal, smbrs, rovsdt)

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use entsor
use lagran
use optcal
use cstphy
use cstnum
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
use pointe
use field
use radiat
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer          iscal

double precision smbrs(ncelet), rovsdt(ncelet)

! Local variables

character(len=80) :: chaine, fname, name
integer          ivar , iel
integer          numcla , numcha , icla
integer          ifcvsl
integer          mode, ige
integer          icha , ii, jj
integer          itermx,nbpauv,nbrich,nbepau,nberic
integer          iterch,nbpass,nbarre,nbimax
integer          iok1
integer          keyccl, f_id

double precision xnuss
double precision aux
double precision coefe(ngazem)
double precision t1, t2, hh2ov
double precision f1mc(ncharm), f2mc(ncharm)
double precision xhdev1 , xhdev2 , xhco , xho2 , gamdv1 , gamdv2
double precision xhco2  , xhh2   , xhh2o

double precision gamhet , den
double precision xxco,xxo2,xxco2,xxh2o
double precision xkp,xkm,t0p,t0m
double precision anmr,tauchi,tautur
double precision sqh2o , x2
double precision err1mx,err2mx
double precision errch,xden
double precision fn0,fn1,fn2,anmr0,anmr1,anmr2
double precision lnk0p,l10k0e,lnk0m,t0e,xco2eq,xcoeq,xo2eq
double precision xcom,xo2m,xkcequ,xkpequ,visls_0

double precision xw1,char_formation,exp_st,imp_st
double precision xo2,wmel,wmhcn,wmno,wmo2,wmnh3
double precision gmdev1(ncharm),gmdev2(ncharm),gmhet(ncharm)
double precision aux1 , aux2 , aux3
double precision xch,xck,xash,xmx2
double precision tfuelmin,tfuelmax
double precision smbrs1, diam2
double precision, dimension (:), allocatable :: w1,w3,w4,w5
double precision, dimension (:), allocatable :: tfuel
double precision, dimension(:), pointer :: vp_x, vp_y, vp_z
double precision, dimension(:,:), pointer :: vdc
double precision, dimension(:), pointer :: crom
double precision, dimension(:), pointer :: cpro_cp, cpro_viscls
double precision, dimension(:), pointer :: cpro_rovsdt2
double precision, dimension(:), pointer :: cvara_k, cvara_ep
double precision, dimension(:), pointer :: cvara_coke
double precision, dimension(:), pointer :: taup
double precision, dimension(:), pointer :: smbrsh1, rovsdth1
double precision, dimension(:,:), pointer :: vel
double precision, dimension(:,:), pointer ::  vg_lim_pi
double precision, dimension(:), pointer :: cvar_x2h2, cvara_x2h2
double precision, dimension(:), pointer :: cvar_xchcl
double precision, dimension(:), pointer :: cvar_xckcl, cvara_xckcl
double precision, dimension(:), pointer :: cvar_xnpcl
double precision, dimension(:), pointer :: cvar_xwtcl, cvara_xwtcl
double precision, dimension(:), pointer :: cvar_yno, cvara_yno
double precision, dimension(:), pointer :: cvara_yhcn, cvara_ynh3
double precision, dimension(:), pointer :: cvara_var
double precision, dimension(:), pointer :: cpro_temp, cpro_temp2, cpro_rom2
double precision, dimension(:), pointer :: cpro_diam2, cpro_cgd1, cpro_cgd2
double precision, dimension(:), pointer :: cpro_cght, cpro_yox, cpro_yco2
double precision, dimension(:), pointer :: cpro_yco, cpro_yh2o, cpro_rom1
double precision, dimension(:), pointer :: cpro_yn2, cpro_cyf1, cpro_cyf2
double precision, dimension(:), pointer :: cpro_exp1, cpro_exp2, cpro_exp3
double precision, dimension(:), pointer :: cpro_exp4, cpro_exp5, cpro_exprb
double precision, dimension(:), pointer :: cpro_mmel, cpro_cgch, cpro_ghco2
double precision, dimension(:), pointer :: cpro_ghh2o, cpro_x2, cpro_csec
double precision, dimension(:), pointer :: cpro_fhcnr, cpro_fhcnd, cpro_fhcnc
double precision, dimension(:), pointer :: cpro_fnh3d, cpro_fnh3c, cpro_cnorb
double precision, dimension(:), pointer :: cpro_fnoch, cpro_cnohc, cpro_fnohc
double precision, dimension(:), pointer :: cpro_fnonh, cpro_cnonh, cpro_fnoth
type(pmapper_double_r1), dimension(:), allocatable :: cvar_xck, cvara_xck
type(pmapper_double_r1), dimension(:), allocatable :: cvar_xch, cvara_xch
type(pmapper_double_r1), dimension(:), allocatable :: cvar_xnp
type(pmapper_double_r1), dimension(:), allocatable :: cvar_xwt
type(pmapper_double_r1), dimension(:), allocatable :: cpro_x2c, cpro_t2
type(pmapper_double_r1), dimension(:), allocatable :: cpro_gmhet
type(pmapper_double_r1), dimension(:), allocatable :: cpro_ghco2a, cpro_ghh2oa
type(pmapper_double_r1), dimension(:), allocatable :: cpro_gmdv1, cpro_gmdv2

! Auxiliary variables
! -------------------
double precision aux4,aux5,auxhet,auxrb1,auxrb2
double precision core1,core2,core3,para2
double precision ychx
! Heterogeneous combustion
! ------------------------
! (Reaction rate of the heterogeneous combustion of char 1,
! reaction rate of the heterogeneous combustion of char 2)
double precision mckcl1, mckcl2
!
double precision, dimension(:), pointer :: dt

type(var_cal_opt) :: vcopt
logical(kind=c_bool) :: log_active

!===============================================================================
! Interfaces
!===============================================================================

interface

  function cs_coal_thconvers1(xesp, f1mc, f2mc, tp)  result(eh) &
    bind(C, name='cs_coal_thconvers1')
    use, intrinsic :: iso_c_binding
    implicit none
    real(c_double), dimension(*) :: xesp
    real(c_double), dimension(*) :: f1mc, f2mc
    real(c_double), value :: tp
    real(c_double) :: eh
  end function cs_coal_thconvers1

end interface

!===============================================================================
! 1. Initialization
!===============================================================================

! --- Scalar number to treat: iscal
! --- Number of the associated to the scalar to treat variable iscal
ivar = isca(iscal)
call field_get_val_prev_s(ivarfl(isca(iscal)), cvara_var)

! --- Name of the variable associated to the scalar to treat iscal
call field_get_label(ivarfl(ivar), chaine)

! --- Number of the physic quantity
call field_get_val_s(icrom, crom)
if (icp.ge.0) call field_get_val_s(icp, cpro_cp)

call field_get_val_s(itemp, cpro_temp)
call field_get_val_s(irom1, cpro_rom1)

call field_get_val_s(iym1(io2), cpro_yox)
call field_get_val_s(iym1(ico2), cpro_yco2)
call field_get_val_s(iym1(ico), cpro_yco)
call field_get_val_s(iym1(ih2o), cpro_yh2o)
call field_get_val_s(immel, cpro_mmel)

call field_get_key_struct_var_cal_opt(ivarfl(ivar), vcopt)

call field_get_val_v(ivarfl(iu), vel)

! Key id of the coal scalar class
call field_get_key_id("scalar_class", keyccl)

! source term for gas enthalpy due to inter-phase fluxes
call field_get_val_s_by_name('x_h_c_exp_st',smbrsh1)
call field_get_val_s_by_name('x_h_c_imp_st',rovsdth1)

if (ieqnox .eq. 1 .and. ntcabs .gt. 1) then
  call field_get_val_s(ivarfl(isca(iyno)), cvar_yno)
  call field_get_val_prev_s(ivarfl(isca(iyno)), cvara_yno)
  call field_get_val_prev_s(ivarfl(isca(iyhcn)), cvara_yhcn)
  call field_get_val_prev_s(ivarfl(isca(iynh3)), cvara_ynh3)
endif

log_active = cs_log_default_is_active()

!===============================================================================
! Deallocation dynamic arrays
!----
allocate(w1(1:ncel),w3(1:ncel),w4(1:ncel),w5(1:ncel),stat=iok1)
if (iok1 > 0) then
  write(nfecra,*) ' Memory allocation error inside: '
  write(nfecra,*) '     cs_coal_scast               '
  call csexit(1)
endif
allocate(tfuel(1:ncel),stat=iok1)
if (iok1 > 0) then
  write(nfecra,*) ' Memory allocation error inside: '
  write(nfecra,*) '     cs_coal_scast               '
  call csexit(1)
endif

call field_get_val_s_by_name('dt', dt)

!===============================================================================
! 2. Consideration of the source terms for relative variables
!    to the classes of particles
!===============================================================================

! --> Source term for the mass fraction of reactive coal

if (ivar.ge.isca(ixch(1)) .and. ivar.le.isca(ixch(nclacp))) then

  if (vcopt%iwarni.ge.1) then
    write(nfecra,1000) chaine(1:8)
  endif

  ! index of the coal particle class
  call field_get_key_int(ivarfl(ivar), keyccl, numcla)

  call field_get_val_s(igmdch(numcla), cpro_cgch)

  do iel = 1, ncel

    ! ---- Calculation of W1 = - rho.GMDCH > 0

    xw1 = - crom(iel)*cpro_cgch(iel)*cell_f_vol(iel)

    ! ---- Calculation of explicit and implicit parts of source terms

    rovsdt(iel) = rovsdt(iel) + max(xw1,zero)
    smbrs (iel) = smbrs(iel)  - xw1*cvara_var(iel)

  enddo

endif

!===============================================================================
! 2.1 Source term for the mass fraction of coke
!===============================================================================

if (ivar.ge.isca(ixck(1)) .and. ivar.le.isca(ixck(nclacp))) then

  if (vcopt%iwarni.ge.1) then
    write(nfecra,1000) chaine(1:8)
  endif

  ! index of the coal particle class
  call field_get_key_int(ivarfl(ivar), keyccl, numcla)

  call field_get_val_s(ivarfl(isca(ixch(numcla))), cvar_xchcl)
  call field_get_val_prev_s(ivarfl(isca(ixck(numcla))), cvara_xckcl)

  call field_get_val_s(igmdch(numcla), cpro_cgch)
  call field_get_val_s(igmdv1(numcla), cpro_cgd1)
  call field_get_val_s(igmdv2(numcla), cpro_cgd2)
  call field_get_val_s(igmhet(numcla), cpro_cght)

  if (ihtco2 .eq. 1) then
    call field_get_val_s(ighco2(numcla), cpro_ghco2)
  endif
  if (ihth2o .eq. 1) then
    call field_get_val_s(ighh2o(numcla), cpro_ghh2o)
  endif

  do iel = 1, ncel
    exp_st = 0.d0
    ! volatile formation minus coal consuming = char formation
    ! (Coke formation in French)
    ! NB: we take values at current and not previous time step
    !     to be conservative in mass
    char_formation = crom(iel)*cvar_xchcl(iel)*cell_f_vol(iel)      &
                   *(cpro_cgd1(iel)+cpro_cgd2(iel)-cpro_cgch(iel))

    ! Compute the implict part of the Source term
    if (cvara_xckcl(iel) .gt. epsicp) then

      ! Reaction C(s) + O2 ---> 0.5CO
      exp_st = cpro_cght(iel)

      ! Reaction C(s) + CO2 ---> 2CO
      if (ihtco2 .eq. 1) exp_st = exp_st + cpro_ghco2(iel)

      ! Reaction C(s) + H2O ---> CO + H2
      if (ihth2o .eq. 1) exp_st = exp_st + cpro_ghh2o(iel)

      exp_st = -2.d0/3.d0 * crom(iel) * exp_st   &
                          * cell_f_vol(iel) / cvara_xckcl(iel)**(1.d0/3.d0)
    endif

    ! Compute the explicit part of the Source term
    imp_st = -3.d0/2.d0 * exp_st * cvara_xckcl(iel)

    rovsdt(iel) = rovsdt(iel) + max(exp_st, 0.d0)
    smbrs(iel) = smbrs(iel) + char_formation + imp_st

  enddo

endif

!===============================================================================
! 2.2 Source term for the mass fraction of water
!===============================================================================

if (ippmod(iccoal) .eq. 1) then

  if (ivar.ge.isca(ixwt(1)) .and. ivar.le.isca(ixwt(nclacp))) then

    if (vcopt%iwarni.ge.1) then
      write(nfecra,1000) chaine(1:8)
    endif

    ! index of the coal particle class
    call field_get_key_int(ivarfl(ivar), keyccl, numcla)

    numcha = ichcor(numcla)

    call field_get_val_s(igmsec(numcla),cpro_csec)
    call field_get_val_s(ix2(numcla),cpro_x2)

    do iel = 1, ncel

     ! ---- Calculation of explicit and implicit parts of source terms

     if (cvara_var(iel).gt.epsicp .and.                         &
          xwatch(numcha).gt.epsicp) then
       xw1 = crom(iel)*cpro_csec(iel)*cell_f_vol(iel)             &
            *(1.d0/cpro_x2(iel))*(1.d0/xwatch(numcha))

       rovsdt(iel) = rovsdt(iel) + max(xw1,zero)
       smbrs(iel)  = smbrs(iel)  - xw1*cvara_var(iel)
     endif

    enddo

  endif

endif

!===============================================================================
! 2.3 Particle age source term
!===============================================================================

if (i_comb_drift.ge.1) then

  call field_get_name(ivarfl(ivar), fname)

  ! Particle age source term
  if (fname(1:7).eq.'n_p_age') then

    ! index of the coal particle class
    call field_get_key_int(ivarfl(ivar), keyccl, icla)

    call field_get_val_s(ivarfl(isca(inp(icla))), cvar_xnpcl)

    do iel = 1, ncel
      smbrs(iel) = smbrs(iel) + crom(iel) * cell_f_vol(iel)*cvar_xnpcl(iel)
    enddo

  endif

  ! Age of the bulk source term
  if (fname(1:3).eq.'age') then

    do iel = 1, ncel
      smbrs(iel) =  smbrs(iel) + crom(iel) * cell_f_vol(iel)
    enddo

  endif

endif

!===============================================================================
! 2.4 Particle velocity source terms
!===============================================================================

if (i_comb_drift.eq.1) then

  call field_get_name(ivarfl(ivar), fname)

  ! index of the coal particle class
  call field_get_key_int(ivarfl(ivar), keyccl, icla)

  if (icla.ge.1) then

    numcha = ichcor(icla)

    call field_get_val_s(igmhet(icla), cpro_cght)

    if (ihtco2 .eq. 1) then
      call field_get_val_s(ighco2(icla), cpro_ghco2)
    endif
    if (ihth2o .eq. 1) then
      call field_get_val_s(ighh2o(icla), cpro_ghh2o)
    endif

    call field_get_val_prev_s(ivarfl(isca(ixck(icla))) , cvara_coke)
    ! Taup
    write(name,'(a,i2.2)')'n_p_' ,icla
    call field_get_id('drift_tau_'//trim(name), f_id)
    call field_get_val_s(f_id, taup)

    if (fname(1:6).eq.'v_x_p_') then

      write(name,'(a,i2.2)')'v_x_p_' ,icla
      call field_get_val_s_by_name(name, vp_x)

      write(name,'(a,i2.2)')'vg_lim_p_' ,icla
      call field_get_val_v_by_name(name, vg_lim_pi)

      ! Vitesse de deviation de la phase continue
      name='vd_c'
      call field_get_val_v_by_name(name, vdc)

      do iel = 1, ncel
        ! Drying and Devolatilization have no effect on Vp
        ! (mass flux is only exiting the particle)
        ! During Heterogeneous reactions gas molecules, with the velocity Vc,
        ! are absorbed and the product, with the velocity Vp, is released.
        ! For Oxygen oxydation one O2 comes in and two CO go out
        ! For CO2 gasification,one  CO2 comes inxg and two CO get out
        ! For H2O gasification, one H2O comes in and one CO and one H2 get out
        smbrs1 = 0.d0
        if (cvara_coke(iel).gt.epsicp) then
          smbrs1                  = smbrs1 + wmole(io2)/wmolat(iatc)*cpro_cght(iel)
          if (ihtco2.eq.1) smbrs1 = smbrs1 + wmole(ico2)/wmolat(iatc)*cpro_ghco2(iel)
          if (ihth2o.eq.1) smbrs1 = smbrs1 + wmole(ih2o)/wmolat(iatc)*cpro_ghh2o(iel)
          smbrs1                 = smbrs1 * cvara_coke(iel)**(2.d0/3.d0)
        endif

        ! relaxation to drop velocity
        smbrs1 = crom(iel)*cell_f_vol(iel)*(1.d0/taup(iel)+smbrs1)                &
               *(vel(1,iel)+vdc(1,iel)+vg_lim_pi(1, iel)-vp_x(iel))

        smbrs(iel) = smbrs(iel) + smbrs1
        rovsdt(iel) = rovsdt(iel) + crom(iel)*cell_f_vol(iel)/taup(iel)

      enddo !sur iel

    elseif (fname(1:6).eq.'v_y_p_') then

      write(name,'(a,i2.2)')'v_y_p_' ,icla
      call field_get_val_s_by_name(name, vp_y)

      write(name,'(a,i2.2)')'vg_lim_p_' ,icla
      call field_get_val_v_by_name(name, vg_lim_pi)

      ! Vitesse de deviation de la phase continue
      name='vd_c'
      call field_get_val_v_by_name(name, vdc)

      do iel = 1, ncel
        smbrs1 = 0.d0
        if (cvara_coke(iel).gt.epsicp) then
                        smbrs1 = smbrs1 + wmole(io2)/wmolat(iatc)*cpro_cght(iel)
        if(ihtco2.eq.1) smbrs1 = smbrs1 + wmole(ico2)/wmolat(iatc)*cpro_ghco2(iel)
        if(ihth2o.eq.1) smbrs1 = smbrs1 + wmole(ih2o)/wmolat(iatc)*cpro_ghh2o(iel)
        smbrs1 = smbrs1 * cvara_coke(iel)**(2.d0/3.d0)
        endif

        ! relaxation to drop velocity
        smbrs1 = crom(iel)*cell_f_vol(iel)*(1.d0/taup(iel)+smbrs1)                &
               *(vel(2,iel)+vdc(2, iel)+vg_lim_pi(2, iel)-vp_y(iel))
        smbrs(iel) = smbrs(iel) + smbrs1
        rovsdt(iel) = rovsdt(iel) + crom(iel)*cell_f_vol(iel)/taup(iel)

      enddo !sur iel


    elseif (fname(1:6).eq.'v_z_p_') then

      write(name,'(a,i2.2)')'v_z_p_' ,icla
      call field_get_val_s_by_name(name, vp_z)

      write(name,'(a,i2.2)')'vg_lim_p_' ,icla
      call field_get_val_v_by_name(name, vg_lim_pi)

      ! Vitesse de deviation de la phase continue
      name='vd_c'
      call field_get_val_v_by_name(name, vdc)

      do iel = 1, ncel
        smbrs1 = 0.d0
        if (cvara_coke(iel).gt.epsicp) then
                        smbrs1 = smbrs1 + wmole(io2)/wmolat(iatc)*cpro_cght(iel)
        if(ihtco2.eq.1) smbrs1 = smbrs1 + wmole(ico2)/wmolat(iatc)*cpro_ghco2(iel)
        if(ihth2o.eq.1) smbrs1 = smbrs1 + wmole(ih2o)/wmolat(iatc)*cpro_ghh2o(iel)
        smbrs1 = smbrs1 * cvara_coke(iel)**(2.d0/3.d0)
        endif

        ! relaxation to drop velocity
        smbrs1 = crom(iel)*cell_f_vol(iel)*(1.d0/taup(iel)+smbrs1)                &
               *(vel(3,iel)+vdc(3, iel)+vg_lim_pi(3, iel)-vp_z(iel))

        smbrs(iel) = smbrs(iel) + smbrs1
        rovsdt(iel) = rovsdt(iel) + crom(iel)*cell_f_vol(iel)/taup(iel)

      enddo !on iel

    endif !on fname

  endif

endif !on icoal drift

!===============================================================================
! 2.5 Source term for the enthalpies
!===============================================================================

if ((ivar.ge.isca(ih2(1)) .and. ivar.le.isca(ih2(nclacp)))) then

  ! Initialization of the exchange terms for gas enthalpy
  if (ivar.eq.isca(ih2(1))) then
    do iel = 1, ncel
      smbrsh1(iel) = 0.d0
      rovsdth1(iel) = 0.d0
    enddo
  endif

  if (vcopt%iwarni.ge.1) then
    write(nfecra,1000) chaine(1:8)
  endif

  ! index of the coal particle class
  call field_get_key_int(ivarfl(ivar), keyccl, numcla)

  call field_get_val_s(ivarfl(isca(ixch(numcla))), cvar_xchcl)
  call field_get_val_s(ivarfl(isca(ixck(numcla))), cvar_xckcl)
  call field_get_val_prev_s(ivarfl(isca(ixck(numcla))), cvara_xckcl)
  if (ippmod(iccoal) .eq. 1) then
    call field_get_val_s(ivarfl(isca(ixwt(numcla))), cvar_xwtcl)
    call field_get_val_prev_s(ivarfl(isca(ixwt(numcla))), cvara_xwtcl)
  endif
  numcha = ichcor(numcla)

  call field_get_val_s(ix2(numcla),cpro_x2)
  call field_get_val_s(irom2(numcla),cpro_rom2)
  call field_get_val_s(idiam2(numcla),cpro_diam2)
  call field_get_val_s(itemp2(numcla),cpro_temp2)
  call field_get_val_s(igmhet(numcla),cpro_cght)

  if (ihtco2 .eq. 1) then
    call field_get_val_s(ighco2(numcla), cpro_ghco2)
  endif
  if (ihth2o .eq. 1) then
    call field_get_val_s(ighh2o(numcla), cpro_ghh2o)
  endif

  call field_get_val_s(igmdch(numcla), cpro_cgch)
  call field_get_val_s(igmdv1(numcla), cpro_cgd1)
  call field_get_val_s(igmdv2(numcla), cpro_cgd2)

  ! ---- Contribution to the explicit and implicit balance
  !        exchanges by molecular distribution
  !        6 Lambda Nu / diam**2 / Rho2 * Rho * (T1-T2)

  ! ------ Calculation of lambda in W1

  xnuss = 2.d0

  call field_get_key_int (ivarfl(isca(iscalt)), kivisl, ifcvsl)
  if (ifcvsl.ge.0) then
    call field_get_val_s(ifcvsl, cpro_viscls)
  endif

  call field_get_key_double(ivarfl(isca(iscalt)), kvisl0, visls_0)

  do iel = 1, ncel
    if (ifcvsl.ge.0) then
      if (icp.ge.0) then
        w1(iel) = cpro_viscls(iel) * cpro_cp(iel)
      else
        w1(iel) = cpro_viscls(iel) * cp0
      endif
    else
      if (icp.ge.0) then
        w1(iel) = visls_0 * cpro_cp(iel)
      else
        w1(iel) = visls_0 * cp0
      endif
    endif
  enddo

  ! ------ Contribution to the explicit and implicit balance
  !        exchanges by molecular distribution
  !      Remark: We use cpro_x2(iel) because we want X2 at the iteration n
  call field_get_val_s(igmtr(numcla), cpro_rovsdt2)
  do iel = 1, ncel
    ! ------ Calculation of diameter of the particles
    !        d20 = (A0.D0**2+(1-A0)*DCK**2)**0.5
    diam2 =  xashch(numcha)*diam20(numcla)**2 +                &
             (1.d0-xashch(numcha))*cpro_diam2(iel)**2

    aux = 6.d0 * w1(iel) * xnuss / diam2       &
        / cpro_rom2(iel) * crom(iel)       &
        * cell_f_vol(iel)

    smbrs(iel)  = smbrs(iel)-aux*(cpro_temp2(iel)-cpro_temp(iel))*cpro_x2(iel)
    smbrsh1(iel) = smbrsh1(iel)+aux*(cpro_temp2(iel)-cpro_temp(iel))*cpro_x2(iel)

    ! Store the implicite part of the exchange so that we can compute a
    ! conservative exhcange term when computing the gas enthalpy
    cpro_rovsdt2(iel) = aux / cp2ch(numcha)
    rovsdt(iel) = rovsdt(iel) + cpro_rovsdt2(iel)

  enddo

  ! ---- Contribution to the explicit and implicit balances
  !        of exchange term of energy between phases:
  !        gama(dev1) H(mv1,T2)+gama(dev2) H(mv2,T2)

  do iel = 1, ncel

    !        Gama Dev1 et Gama Dev2

    gamdv1 = crom(iel)*cvar_xchcl(iel)                   &
            *cpro_cgd1(iel)

    gamdv2 = crom(iel)*cvar_xchcl(iel)                   &
            *cpro_cgd2(iel)

    !        H(mv1,T2)

    do ige = 1, ngazem
      coefe(ige) = zero
    enddo

    den = a1(numcha)*wmole(ichx1c(numcha))+b1(numcha)*wmole(ico)  &
         +c1(numcha)*wmole(ih2o)          +d1(numcha)*wmole(ih2s) &
         +e1(numcha)*wmole(ihcn)          +f1(numcha)*wmole(inh3)
    coefe(ichx1) = a1(numcha)*wmole(ichx1c(numcha)) /den
    coefe(ico ) = b1(numcha)*wmole(ico)            /den
    coefe(ih2o) = c1(numcha)*wmole(ih2o)           /den
    coefe(ih2s) = d1(numcha)*wmole(ih2s)           /den
    coefe(ihcn) = e1(numcha)*wmole(ihcn)           /den
    coefe(inh3) = f1(numcha)*wmole(inh3)           /den

    t2         = cpro_temp2(iel)
    do icha = 1, ncharm
      f1mc(icha) = zero
      f2mc(icha) = zero
    enddo
    f1mc(numcha) = 1.d0

    xhdev1 = cs_coal_thconvers1(coefe, f1mc, f2mc, t2)

    !        H(mv2,T2)

    do ige = 1, ngazem
      coefe(ige) = zero
    enddo
    den = a2(numcha)*wmole(ichx2c(numcha))+b2(numcha)*wmole(ico)  &
         +c2(numcha)*wmole(ih2o)          +d2(numcha)*wmole(ih2s) &
         +e2(numcha)*wmole(ihcn)          +f2(numcha)*wmole(inh3)
    coefe(ichx2) = a2(numcha)*wmole(ichx2c(numcha)) /den
    coefe(ico ) = b2(numcha)*wmole(ico)            /den
    coefe(ih2o) = c2(numcha)*wmole(ih2o)           /den
    coefe(ih2s) = d2(numcha)*wmole(ih2s)           /den
    coefe(ihcn) = e2(numcha)*wmole(ihcn)           /den
    coefe(inh3) = f2(numcha)*wmole(inh3)           /den

    t2         = cpro_temp2(iel)
    do icha = 1, ncharm
      f1mc(icha) = zero
      f2mc(icha) = zero
    enddo
    f2mc(numcha) = 1.d0

    xhdev2 = cs_coal_thconvers1(coefe, f1mc, f2mc, t2)

    !         Contribution to explicit and implicit balances

     smbrs(iel)   = smbrs(iel)   + (gamdv1*xhdev1+gamdv2*xhdev2)*cell_f_vol(iel)
     smbrsh1(iel) = smbrsh1(iel) - (gamdv1*xhdev1+gamdv2*xhdev2)*cell_f_vol(iel)
  enddo

  ! ------ Heterogeneous combustion: C(s) + 02 ---> 0.5 C0
  !        GamHET * (28/12 H(CO,T2)-16/12 H(O2,T1))

  do iel = 1, ncel

    !        Calculation of HCO(T2)

    do ige = 1, ngazem
      coefe(ige) = zero
    enddo
    coefe(ico) = 1.d0
    do icha = 1, ncharm
      f1mc(icha) = zero
      f2mc(icha) = zero
    enddo

    t2        = cpro_temp2(iel)
    xhco      = cs_coal_thconvers1(coefe, f1mc, f2mc, t2)

    !        Calculation of HO2(T1)

    do ige = 1, ngazem
      coefe(ige) = zero
    enddo
    coefe(io2) = 1.d0
    do icha = 1, ncharm
      f1mc(icha) = zero
      f2mc(icha) = zero
    enddo

    t1        = cpro_temp(iel)
    xho2      = cs_coal_thconvers1(coefe, f1mc, f2mc, t1)

    !         Contribution to explicit and implicit balances

    if (cvara_xckcl(iel) .gt. epsicp) then

      gamhet = crom(iel)*cpro_cght(iel)              &
               * ((cvara_xckcl(iel))**(2.d0/3.d0) +              &
                2.d0/3.d0*(cvar_xckcl(iel)-cvara_xckcl(iel))      &
                 /(cvara_xckcl(iel))**(1.d0/3.d0))

    else
      gamhet = 0.d0
    endif

    gamhet = gamhet *(wmole(ico)*xhco-wmolat(iato)*xho2)/wmolat(iatc) * cell_f_vol(iel)

    smbrs(iel) = smbrs(iel)     + gamhet
    smbrsh1(iel) = smbrsh1(iel) - gamhet
  enddo

  ! ------ Heterogeneous combustion: C(s) + C02 ---> 2 C0
  !        GamHET * (56/12 H(CO,T2)-44/12 H(CO2,T1))

  if (ihtco2 .eq. 1) then
    do iel = 1, ncel

      !        Calculation of HCO(T2)

      do ige = 1, ngazem
        coefe(ige) = zero
      enddo
      coefe(ico) = 1.d0
      do icha = 1, ncharm
        f1mc(icha) = zero
        f2mc(icha) = zero
      enddo

      t2        = cpro_temp2(iel)
      xhco      = cs_coal_thconvers1(coefe, f1mc, f2mc, t2)

      !        Calculation of HCO2(T1)

      do ige = 1, ngazem
        coefe(ige) = zero
      enddo
      coefe(ico2) = 1.d0
      do icha = 1, ncharm
        f1mc(icha) = zero
        f2mc(icha) = zero
      enddo

      t1        = cpro_temp(iel)
      xhco2     = cs_coal_thconvers1(coefe, f1mc, f2mc, t1)

      !         Contribution to explicit and implicit balances

      if (cvara_xckcl(iel) .gt. epsicp) then

        gamhet = crom(iel)*cpro_ghco2(iel)           &
                 * ((cvara_xckcl(iel))**(2.d0/3.d0) +            &
                2.d0/3.d0*(cvar_xckcl(iel)-cvara_xckcl(iel))      &
                 /(cvara_xckcl(iel))**(1.d0/3.d0))

      else
        gamhet = 0.d0
      endif

      gamhet =  gamhet*(2.d0*wmole(ico)*xhco-wmole(ico2)*xhco2)/wmolat(iatc) &
              * cell_f_vol(iel)

      smbrs(iel)   = smbrs(iel)   + gamhet
      smbrsh1(iel) = smbrsh1(iel) - gamhet
    enddo

  endif

  ! ------ Heterogeneous combustion: C(s) + H2O ---> CO + H2
  !        GamHET * (28/12 H(CO,T2)+2/12 H(HY,T2) -18/12 H(H2O,T1))

  if (ihth2o .eq. 1) then
    do iel = 1, ncel

      !        Calculation of HCO(T2)

      do ige = 1, ngazem
        coefe(ige) = zero
      enddo
      coefe(ico) = 1.d0
      do icha = 1, ncharm
        f1mc(icha) = zero
        f2mc(icha) = zero
      enddo

      t2        = cpro_temp2(iel)
      xhco      = cs_coal_thconvers1(coefe, f1mc, f2mc, t2)

      !        Calculation of HH2(T2)

      do ige = 1, ngazem
        coefe(ige) = zero
      enddo
      coefe(ihy) = 1.d0
      do icha = 1, ncharm
        f1mc(icha) = zero
        f2mc(icha) = zero
      enddo

      t2        = cpro_temp2(iel)
      xhh2      = cs_coal_thconvers1(coefe, f1mc, f2mc, t2)

      !        Calculation of HH2O(T1)

      do ige = 1, ngazem
        coefe(ige) = zero
      enddo
      coefe(ih2o) = 1.d0
      do icha = 1, ncharm
        f1mc(icha) = zero
        f2mc(icha) = zero
      enddo

      t1        = cpro_temp(iel)
      xhh2o     = cs_coal_thconvers1(coefe, f1mc, f2mc, t1)

      !         Contribution to explicit and implicit balances

      if (cvara_xckcl(iel) .gt. epsicp) then

        gamhet = crom(iel)*cpro_ghh2o(iel)           &
                 * ((cvara_xckcl(iel))**(2.d0/3.d0) +            &
                2.d0/3.d0*(cvar_xckcl(iel)-cvara_xckcl(iel))      &
                 /(cvara_xckcl(iel))**(1.d0/3.d0))

      else
        gamhet = 0.d0
      endif

      gamhet = gamhet * (wmole(ico)*xhco+wmole(ihy)*xhh2   &
                        -wmole(ih2o)*xhh2o)/wmolat(iatc)   &
                      *cell_f_vol(iel)

      smbrs(iel)   = smbrs(iel)   + gamhet
      smbrsh1(iel) = smbrsh1(iel) - gamhet
    enddo

  endif

  !       --> Source term on H2 (coming from drying)

  if (ippmod(iccoal) .eq. 1) then

    ! ---- Contribution of source term interfacial to explicit and implicit balances

    call field_get_val_s(igmsec(numcla),cpro_csec)
    call field_get_val_s(itemp2(numcla),cpro_temp2)

    do iel = 1, ncel

      !          Calculation of H(H2O) at T2

      do ige = 1, ngazem
        coefe(ige) = zero
      enddo
      coefe(ih2o) = 1.d0
      do icha = 1, ncharm
        f1mc(icha) = zero
        f2mc(icha) = zero
      enddo

      t2 = cpro_temp2(iel)
      if (t2 .gt. 100.d0+tkelvi) then
        t2 = 100.d0+tkelvi
      endif
      mode      = -1

      call cpthp1(mode, hh2ov, coefe, f1mc, f2mc, t2)

      ! Contribution to explicit balance

      if (cvara_xwtcl(iel).gt. epsicp .and.          &
           xwatch(numcha) .gt. epsicp      ) then

        aux = -crom(iel)*cpro_csec(iel)              &
       *(cvar_xwtcl(iel)/cpro_x2(iel))             &
       *(1.d0                    /xwatch(numcha))                 &
             *hh2ov

      else
        aux = 0.d0
      endif

      smbrs(iel)   = smbrs(iel)   + aux*cell_f_vol(iel)
      smbrsh1(iel) = smbrsh1(iel) - aux*cell_f_vol(iel)

    enddo

  endif  ! if : sechage
endif    ! enthalpies des particules

!===============================================================================
! 3. Taking into account source terms for relative variables in the mixture
!===============================================================================

if (ivar .eq. isca(ihgas)) then

  ! source terms from particles (convection and enthalpy drived by mass fluxes)
  do iel = 1, ncel
    smbrs(iel) = smbrs(iel) + smbrsh1(iel)
    rovsdt(iel) = rovsdt(iel) + rovsdth1(iel)
  enddo

  ! Explicit contribution due to implicit source term on particle class enthalpy
  ! TODO adapt it to other classes (fuel..)
  do icla = 1, nclacp
    call field_get_val_s(igmtr(icla), cpro_rovsdt2)
    call field_get_val_s(ivarfl(isca(ih2(icla))), cvar_x2h2)
    call field_get_val_prev_s(ivarfl(isca(ih2(icla))), cvara_x2h2)
    do iel = 1, ncel
      smbrs(iel) = smbrs(iel)  &
                 + cpro_rovsdt2(iel)*(cvar_x2h2(iel)-cvara_x2h2(iel))
    enddo
  enddo

endif

! --> Source term for light volatile materials

if (ivar.ge.isca(if1m(1)) .and. ivar.le.isca(if1m(ncharb))) then

  if (vcopt%iwarni.ge.1) then
    write(nfecra,1000) chaine(1:8)
  endif

! ---- Calculation of GMDEV1 = - Sum (rho.XCH.GMDV1) > 0  --> W1

  numcha = ivar-isca(if1m(1))+1
  do iel = 1, ncel
    w1(iel) = zero
  enddo
  do icla = 1, nclacp
    call field_get_val_s(igmdv1(icla), cpro_cgd1)
    call field_get_val_s(ivarfl(isca(ixch(icla))), cvar_xchcl)
    if (ichcor(icla).eq.numcha) then
      do iel = 1, ncel
        w1(iel) = w1(iel) - crom(iel)*cvar_xchcl(iel)    &
                * cpro_cgd1(iel)
      enddo
    endif
  enddo

  ! Contribution of interfacial source term to explicit and implicit balances

  do iel = 1, ncel
    smbrs(iel)  = smbrs(iel)  + cell_f_vol(iel) * w1(iel)
  enddo

endif


! Source terms for heavy volatil materials

if (ivar.ge.isca(if2m(1)) .and. ivar.le.isca(if2m(ncharb))) then

  if (vcopt%iwarni.ge.1) then
    write(nfecra,1000) chaine(1:8)
  endif

  ! Calculation of GMDEV2 = - Sum (rho.XCH.GMDV2) > 0 --> W1

  numcha = ivar-isca(if2m(1))+1
  do iel = 1, ncel
    w1(iel) = zero
  enddo
  do icla = 1, nclacp
    call field_get_val_s(igmdv2(icla), cpro_cgd2)
    call field_get_val_s(ivarfl(isca(ixch(icla))), cvar_xchcl)
    if (ichcor(icla).eq.numcha) then
      do iel = 1, ncel
        w1(iel) = w1(iel) - crom(iel)*cvar_xchcl(iel)    &
                * cpro_cgd2(iel)
      enddo
    endif
  enddo

  ! Contribution of interfacial source term to explicite balance

  do iel = 1, ncel
    smbrs(iel)  = smbrs(iel)  + cell_f_vol(iel) * w1(iel)
  enddo

endif


! Source term for the tracer 7 (O2) (heterogeneous combustion by C)

if (ivar.eq.isca(if7m)) then

  ! Remark: We take the same source term than for Xck
  !                  to be conservative

  if (vcopt%iwarni.ge.1) then
    write(nfecra,1000) chaine(1:8)
  endif

  do iel = 1, ncel
    w1(iel) = zero
  enddo

  do icla = 1, nclacp
    call field_get_val_s(igmhet(icla),cpro_cght)
    call field_get_val_s(ivarfl(isca(ixck(icla))), cvar_xckcl)
    call field_get_val_prev_s(ivarfl(isca(ixck(icla))), cvara_xckcl)
    do iel = 1, ncel
      if (cvara_xckcl(iel) .gt. epsicp) then
        w1(iel) = w1(iel)                                         &
                 - crom(iel)*cpro_cght(iel)                       &
                 * ((cvara_xckcl(iel))**(2.d0/3.d0) +            &
                    2.d0/3.d0*(cvar_xckcl(iel)-cvara_xckcl(iel))  &
                    /(cvara_xckcl(iel))**(1.d0/3.d0))
      endif
    enddo
  enddo

  ! Contribution of interfacial source term to explicit and implicit balances

  do iel = 1, ncel
    smbrs(iel)  = smbrs(iel)  + cell_f_vol(iel) * w1(iel)
  enddo

endif

! Source term for the tracer 8 (CO2) (heterogeneous combustion by C)

if (ihtco2 .eq. 1) then
  if (ivar.eq.isca(if8m)) then

    ! Remark: We take the same source term than for Xck
    !                  to be conservative

    if (vcopt%iwarni.ge.1) then
      write(nfecra,1000) chaine(1:8)
    endif

    do iel = 1, ncel
      w1(iel) = zero
    enddo

    do icla = 1, nclacp
      call field_get_val_s(ighco2(icla), cpro_ghco2)
      call field_get_val_s(ivarfl(isca(ixck(icla))), cvar_xckcl)
      call field_get_val_prev_s(ivarfl(isca(ixck(icla))), cvara_xckcl)
      do iel = 1, ncel
        if (cvara_xckcl(iel) .gt. epsicp) then
          w1(iel) = w1(iel)                                        &
                   - crom(iel)*cpro_ghco2(iel)         &
                   * ((cvara_xckcl(iel))**(2.d0/3.d0) +           &
                      2.d0/3.d0*(cvar_xckcl(iel)-cvara_xckcl(iel)) &
                      /(cvara_xckcl(iel))**(1.d0/3.d0))
        endif
      enddo
    enddo

    ! ---- Contribution of interfacial source term to explicit and implicit balances

    do iel = 1, ncel
      smbrs(iel)  = smbrs(iel)  + cell_f_vol(iel) * w1(iel)
    enddo

  endif

endif

! --> Source term for the tracer 9 (H2O) (heterogeneous combustion by H2O)

if (ihth2o .eq. 1) then
  if (ivar.eq.isca(if9m)) then

    ! Remark: We take the same source term than for Xck
    !                  to be conservative

    if (vcopt%iwarni.ge.1) then
      write(nfecra,1000) chaine(1:8)
    endif

    do iel = 1, ncel
      w1(iel) = zero
    enddo

    do icla = 1, nclacp
      call field_get_val_s(ighh2o(icla), cpro_ghh2o)
      call field_get_val_s(ivarfl(isca(ixck(icla))), cvar_xckcl)
      call field_get_val_prev_s(ivarfl(isca(ixck(icla))), cvara_xckcl)
      do iel = 1, ncel
        if (cvara_xckcl(iel) .gt. epsicp) then
          w1(iel) = w1(iel)                                        &
                   - crom(iel)*cpro_ghh2o(iel)        &
                   * ((cvara_xckcl(iel))**(2.d0/3.d0) +           &
                      2.d0/3.d0*(cvar_xckcl(iel)-cvara_xckcl(iel)) &
                      /(cvara_xckcl(iel))**(1.d0/3.d0))
        endif
      enddo
    enddo

    ! ---- Contribution of interfacial source term to explicit and implicit balances

    do iel = 1, ncel
      smbrs(iel)  = smbrs(iel)  + cell_f_vol(iel) * w1(iel)
    enddo

  endif

endif

! --> Source term for the fuel variance

if (ivar.eq.isca(ifvp2m)) then

  if (vcopt%iwarni.ge.1) then
    write(nfecra,1000) chaine(1:8)
  endif

  call cs_coal_fp2st(iscal, smbrs, rovsdt)

endif

! --> Source term for the tracer 6 (Water coming from drying)

if (ippmod(iccoal) .eq. 1) then

  if (ivar.eq.isca(if6m)) then


    if (vcopt%iwarni.ge.1) then
      write(nfecra,1000) chaine(1:8)
    endif

    ! ---- Contribution of interfacial source term to explicit and implicit balances

    do iel = 1, ncel
      w1(iel) = zero
    enddo

    do icla = 1, nclacp

      call field_get_val_s(ivarfl(isca(ixwt(icla))), cvar_xwtcl)
      call field_get_val_prev_s(ivarfl(isca(ixwt(icla))), cvara_xwtcl)
      call field_get_val_s(igmsec(icla),cpro_csec)
      call field_get_val_s(ix2(icla),cpro_x2)

      numcha = ichcor(icla)

      do iel = 1, ncel

        if (      cvara_xwtcl(iel).gt. epsicp               &
            .and. xwatch(numcha) .gt. epsicp) then

          w1(iel) =   w1(iel)                               &
                    +   crom(iel)*cpro_csec(iel)            &
                      * (cvar_xwtcl(iel)/cpro_x2(iel))      &
                      * (1.d0 /xwatch(numcha))

        endif

      enddo

    enddo

    do iel = 1, ncel
      smbrs(iel)  = smbrs(iel)  + cell_f_vol(iel) * w1(iel)
    enddo

  endif

endif

! --> Source term for CO2

if (ieqco2 .eq. 1) then

  if (ivar.eq.isca(iyco2)) then

    if (vcopt%iwarni.ge.1) then
      write(nfecra,1000) chaine(1:8)
    endif

    call field_get_val_prev_s(ivarfl(ik), cvara_k)
    call field_get_val_prev_s(ivarfl(iep), cvara_ep)

    allocate(cpro_x2c(nclacp))

    do icla = 1, nclacp
      call field_get_val_s(ix2(icla),cpro_x2c(icla)%p)
    enddo

    ! Contribution of interfacial source term to explicit and implicit balances

    ! Oxydation of CO
    ! ===============

    !  Dryer Glassman : XK0P in (mol/m3)**(-0.75) s-1
    !          XK0P = 1.26D10
    !           XK0P = 1.26D7 * (1.1)**(NTCABS)
    !           IF (XK0P .GT. 1.26D10) XK0P=1.26D10
    !           T0P  = 4807.D0
    !  Howard : XK0P in (mol/m3)**(-0.75) s-1
    !             XK0P = 4.11D9
    !             T0P  = 15090.D0
    !  Westbrook & Dryer

    lnk0p = 23.256d0
    t0p  = 20096.d0

    !  Hawkin et Smith Purdue University Engeneering Bulletin, i
    !  Research series 108 vol 33, n 3n 1949
    !  Kp = 10**(4.6-14833/T)
    !  Equilibrum constant in partial pressure [atm]
    !  XKOE is the decimal log of the pre-exponential constant
    !  TOE is NOT an activation temperature ... there is a lg(e)
    !  to return to Kc and to use concentrations (in mol/m3)
    !  Kc = (1/RT)**variation nb moles * Kp
    !  here Kc = sqrt(0.082*T)*Kp

    l10k0e = 4.6d0
    t0e  = 14833.d0

    ! Dissociation of CO2 (Trinh Minh Chinh)
    ! ===================
    !          XK0M = 5.D8
    !          T0M  = 4807.D0
    !          XK0M = 0.D0
    !  Westbrook & Dryer

    lnk0m = 20.03d0
    t0m  = 20096.d0

    err1mx = 0.d0
    err2mx = 0.d0

    ! Number of iterations
    itermx = 500
    ! Number of convergent points

   nbpauv = 0
   nbepau = 0
   nbrich = 0
   nberic = 0
   nbpass = 0
   nbarre = 0
   nbimax = 0
   ! Precision on the convergence
   errch = 1.d-8

   do iel = 1, ncel

     xxco  = cpro_yco(iel)/wmole(ico)            &
            *cpro_rom1(iel)
     xxo2  = cpro_yox(iel)/wmole(io2)            &
            *cpro_rom1(iel)
     xxco2 = cpro_yco2(iel)/wmole(ico2)          &
            *cpro_rom1(iel)
     xxh2o = cpro_yh2o(iel)/wmole(ih2o)          &
            *cpro_rom1(iel)

     xxco  = max(xxco ,zero)
     xxo2  = max(xxo2 ,zero)
     xxco2 = max(xxco2,zero)
     xxh2o = max(xxh2o,zero)
     sqh2o = sqrt(xxh2o)

     xkp = exp(lnk0p-t0p/cpro_temp(iel))
     xkm = exp(lnk0m-t0m/cpro_temp(iel))

     xkpequ = 10.d0**(l10k0e-t0e/cpro_temp(iel))
     xkcequ = xkpequ                                              &
             /sqrt(8.32d0*cpro_temp(iel)/1.015d5)

     !        initialization by the transported state

     anmr  = xxco2
     xcom  = xxco + xxco2
     xo2m  = xxo2 + 0.5d0*xxco2

     if (cpro_temp(iel) .gt. 1200.d0) then

      !        Search for the equilibrum state
      !        Iterative search without control of convergence
      !         (to keep the parallelisation on meshes)
      !        On the numver of moles of separating reaction
      !         the state before reaction (such as calculated by Cpcym)
      !         of the equilibrum state
      !        anmr has to be the boundary between 0 and Min(XCOM,2.*XO2M)
      !        We look for the solution by dichotomy

       anmr0 = 0.d0
       anmr1 = min(xcom,2.d0*xo2m)
       iterch = 0
       fn2    = 1.d0
       fn0  = -0.5d0                           * anmr0**3         &
            + (     xcom    + xo2m - xkcequ**2) * anmr0**2        &
            - (.5d0*xcom    +2.d0*xo2m)*xcom   * anmr0            &
            +       xcom**2 * xo2m
       fn1  = -0.5d0                           * anmr1**3         &
            + (     xcom    + xo2m - xkcequ**2) * anmr1**2        &
            - (.5d0*xcom    +2.d0*xo2m)*xcom   * anmr1            &
            +       xcom**2 * xo2m

       if (xo2m.gt.1.d-6) then
         do while (iterch.lt.itermx .and. fn2.gt.errch)
           anmr2 = 0.5d0*(anmr0+anmr1)
           fn2  = -0.5d0                            * anmr2**3    &
                + (     xcom    + xo2m - xkcequ**2) * anmr2**2    &
                - (.5d0*xcom    +2.d0*xo2m)*xcom    * anmr2       &
                +       xcom**2 * xo2m
           if(fn0*fn2 .gt. 0.d0) then
             anmr0 = anmr2
             fn0 = fn2
           elseif(fn1*fn2 .gt. 0.d0) then
             anmr1 = anmr2
             fn1 = fn2
           elseif(fn0*fn1 .gt. 0.d0) then
             iterch = itermx
             anmr2 = min(xcom,2.d0*xo2m)
             nbarre = nbarre + 1
           endif
           iterch = iterch + 1
         enddo

         if (iterch .ge. itermx) then
           nberic = nberic + 1
         else
           nbimax = max(nbimax,iterch)
         endif
         err1mx = max(err1mx,fn2)

         xco2eq = anmr2
         xcoeq  = xcom - anmr2
         xo2eq  = xo2m - 0.5d0 * anmr2
       else
         xo2eq  = 0.d0
         xcoeq  = xxco
         xco2eq = 0.d0
       endif

     else

       xco2eq = min(xcom,2.d0*xo2m)
       xo2eq  = xo2m - 0.5d0*xco2eq
       xcoeq  = xcom - xco2eq

     endif

     if (xco2eq.gt.xxco2) then
       ! oxydation
       xden = xkp*sqh2o*(xxo2)**0.25d0
     else
       ! dissociation
       xden = xkm
     endif
     if (xden .ne. 0.d0) then

       tauchi = 1.d0/xden
       tautur = cvara_k(iel)/cvara_ep(iel)

       x2 = 0.d0
       do icla = 1, nclacp
         x2 = x2 + cpro_x2c(icla)%p(iel)
       enddo

       if (ieqco2 .eq. 1) then
       !    We transport CO2

       smbrs(iel)  = smbrs(iel)                                   &
                    +wmole(ico2)/cpro_rom1(iel)        &
         * (xco2eq-xxco2)/(tauchi+tautur)                         &
         * (1.d0-x2)                                              &
         * cell_f_vol(iel) * crom(iel)

       else if (ieqco2 .eq. 2) then
       !    We transport CO

       smbrs(iel)  = smbrs(iel)                                   &
                    +wmole(ico)/cpro_rom1(iel)        &
         * (xco2eq-xxco)/(tauchi+tautur)                         &
         * (1.d0-x2)                                              &
         * cell_f_vol(iel) * crom(iel)
       endif

       w1(iel) = cell_f_vol(iel)*crom(iel)/(tauchi+tautur)
       rovsdt(iel) = rovsdt(iel) +   max(w1(iel),zero)

     else
       rovsdt(iel) = rovsdt(iel) + 0.d0
       smbrs(iel)  = smbrs(iel)  + 0.d0
     endif

   enddo

   if(irangp.ge.0) then
     call parcpt(nberic)
     call parmax(err1mx)
     call parcpt(nbpass)
     call parcpt(nbarre)
     call parcpt(nbarre)
     call parcmx(nbimax)
   endif

   if (log_active .eqv. .true.) then
     write(nfecra,*) ' Max Error = ', err1mx
     write(nfecra,*) ' no Points   ', nberic, nbarre, nbpass
     write(nfecra,*) ' Iter max number ', nbimax
   endif

   !     Source term: heterogeneous combustion by CO2

   if (ihtco2 .eq. 1) then

     ! Arrays of pointers containing the fields values for each class
     ! (loop on cells outside loop on classes)
     allocate(cvara_xck(nclacp))
     allocate(cpro_ghco2a(nclacp))
     do icla = 1,nclacp
       call field_get_val_prev_s(ivarfl(isca(ixck(icla))), cvara_xck(icla)%p)
       call field_get_val_s(ighco2(icla), cpro_ghco2a(icla)%p)
     enddo

     do iel = 1, ncel

       aux = 0.d0
       do icla = 1,nclacp

         aux = aux                                                &
              + crom(iel)*cpro_ghco2a(icla)%p(iel)                &
               *(cvara_xck(icla)%p(iel))**(2.d0/3.d0)*cell_f_vol(iel)

       enddo

       rovsdt(iel) = rovsdt(iel) - aux*(wmole(ico2)/0.012)

     enddo

     deallocate(cvara_xck)
     deallocate(cpro_x2c, cpro_ghco2a)

   endif

  endif

endif

! --> Source term for Enth_Ox
!                       HCN and NO: only from the second iteration

if (ieqnox .eq. 1 .and. ntcabs .gt. 1) then

  ! Terms on Oxydant enthalpy

  if (ivar .eq. isca(ihox)) then

    ! Arrays of pointers containing the fields values for each class
    ! (loop on cells outside loop on classes)
    allocate(cvar_xck(nclacp), cvara_xck(nclacp))
    allocate(cvar_xch(nclacp), cvar_xnp(nclacp))
    allocate(cpro_t2(nclacp),cpro_gmhet(nclacp))

    if (ihtco2 .eq. 1) then
      allocate(cpro_ghco2a(nclacp))
    endif

    if (ihth2o .eq. 1) then
      allocate(cpro_ghh2oa(nclacp))
    endif

    if (ippmod(iccoal) .eq. 1) then
      allocate(cvar_xwt(nclacp))
    endif
    do icla = 1,nclacp
      call field_get_val_s(ivarfl(isca(ixck(icla))), cvar_xck(icla)%p)
      call field_get_val_prev_s(ivarfl(isca(ixck(icla))), cvara_xck(icla)%p)
      call field_get_val_s(ivarfl(isca(ixch(icla))), cvar_xch(icla)%p)
      call field_get_val_s(ivarfl(isca(inp(icla))), cvar_xnp(icla)%p)
      call field_get_val_s(itemp2(icla), cpro_t2(icla)%p)
      call field_get_val_s(igmhet(icla), cpro_gmhet(icla)%p)
      if (ihtco2 .eq. 1) then
        call field_get_val_s(ighco2(icla), cpro_ghco2a(icla)%p)
      endif
      if (ihth2o .eq. 1) then
        call field_get_val_s(ighh2o(icla), cpro_ghh2oa(icla)%p)
      endif
      if (ippmod(iccoal) .eq. 1) then
        call field_get_val_s(ivarfl(isca(ixwt(icla))), cvar_xwt(icla)%p)
      endif
    enddo

    !  Calculation of T2 average on particles

    tfuelmin = 1.d+20
    tfuelmax =-1.d+20
    do iel=1,ncel

      xmx2 = 0.d0
      do icla = 1, nclacp
        xck  = cvar_xck(icla)%p(iel)
        xch  = cvar_xch(icla)%p(iel)
        xash = cvar_xnp(icla)%p(iel)*xmash(icla)
        xmx2   = xmx2 + xch + xck + xash

        !   Taking into account humidity

        if (ippmod(iccoal) .eq. 1) then
          xmx2 = xmx2+cvar_xwt(icla)%p(iel)
        endif
      enddo

      if (xmx2 .gt. 0.d0) then
        tfuel(iel) = 0.d0
        do icla=1,nclacp
          tfuel(iel) = tfuel(iel)                                              &
                      +( cvar_xck(icla)%p(iel)                                 &
                        +cvar_xch(icla)%p(iel)                                 &
                        +cvar_xnp(icla)%p(iel)*xmash(icla))                   &
                        *cpro_t2(icla)%p(iel)

          !  Taking into account humidity

          if (ippmod(iccoal) .eq. 1) then
            tfuel(iel) = tfuel(iel) + (cvar_xwt(icla)%p(iel))              &
                         *cpro_t2(icla)%p(iel)
          endif
        enddo

        tfuel(iel) = tfuel(iel)/xmx2

      else
        tfuel(iel) = cpro_temp(iel)
      endif
      tfuelmin = min(tfuel(iel),tfuelmin)
      tfuelmax = max(tfuel(iel),tfuelmax)
    enddo
    if (irangp .ge. 0) then
      call parmin(tfuelmin)
      call parmax(tfuelmax)
    endif
    if (log_active .eqv. .true.) then
      write(nfecra,*) ' Min max of Tfuel for Hoxy ', tfuelmin, tfuelmax
    endif

    ! Heterogeneous combustion: C + O2 ---> 0.5 CO

    do iel=1,ncel

      !   Calculation of HCO(T2)

      do ige = 1, ngazem
        coefe(ige) = zero
      enddo
      coefe(ico) = 1.d0
      do icha = 1, ncharm
        f1mc(icha) = zero
        f2mc(icha) = zero
      enddo
      t2        = tfuel(iel)
      xhco      = cs_coal_thconvers1(coefe, f1mc, f2mc,  t2)

      !  Calculation of HO2(T1)

      do ige = 1, ngazem
        coefe(ige) = zero
      enddo
      coefe(io2) = 1.d0
      do icha = 1, ncharm
        f1mc(icha) = zero
        f2mc(icha) = zero
      enddo
      t1        = cpro_temp(iel)
      xho2      = cs_coal_thconvers1(coefe, f1mc, f2mc, t1)

      do icla=1,nclacp
        if (cvara_xck(icla)%p(iel) .gt. epsicp) then
          gamhet = crom(iel)*cpro_gmhet(icla)%p(iel)                        &
                   * ((cvara_xck(icla)%p(iel))**(2.d0/3.d0) +              &
                  2.d0/3.d0*(cvar_xck(icla)%p(iel)-cvara_xck(icla)%p(iel))  &
                    /(cvara_xck(icla)%p(iel))**(1.d0/3.d0))
        else
          gamhet = 0.d0
        endif
        smbrs(iel) = smbrs(iel)                                        &
                    -gamhet                                            &
                     *(wmole(ico)*xhco-wmolat(iato)*xho2)/wmolat(iatc)*cell_f_vol(iel)
!
      enddo

    enddo

    !  Heterogeneous combustion: C + CO2 ---> 2 CO

    !  Calculation of HO2(T1)

    if (ihtco2 .eq. 1) then
      do iel=1,ncel

        !  Calculation of HCO(T2)

        do ige = 1, ngazem
          coefe(ige) = zero
        enddo
        coefe(ico) = 1.d0
        do icha = 1, ncharm
          f1mc(icha) = zero
          f2mc(icha) = zero
        enddo

        t2        = tfuel(iel)
        xhco      = cs_coal_thconvers1(coefe, f1mc, f2mc, t2)

        !  Calculation of HCO2(T1)

        do ige = 1, ngazem
          coefe(ige) = zero
        enddo
        coefe(ico2) = 1.d0
        do icha = 1, ncharm
          f1mc(icha) = zero
          f2mc(icha) = zero
        enddo
        t1        = cpro_temp(iel)
        xhco2     = cs_coal_thconvers1(coefe, f1mc, f2mc, t1)

        do icla=1,nclacp
          if (cvara_xck(icla)%p(iel) .gt. epsicp) then
            gamhet = crom(iel)*cpro_ghco2a(icla)%p(iel)              &
                     * ((cvara_xck(icla)%p(iel))**(2.d0/3.d0) +               &
                    2.d0/3.d0*(cvar_xck(icla)%p(iel)-cvara_xck(icla)%p(iel))   &
                      /(cvara_xck(icla)%p(iel))**(1.d0/3.d0))
          else
            gamhet = 0.d0
          endif
          smbrs(iel) = smbrs(iel)                                          &
                      -gamhet                                              &
                       *(2.d0*wmole(ico)*xhco-wmole(ico2)*xhco2)/wmolat(iatc) *cell_f_vol(iel)

        enddo

      enddo

    endif

    !   Heterogeneous combustion: C + H2O ---> CO + H2

    !     Calculation of HO2(T1)

    if (ihth2o .eq. 1) then

      do iel=1,ncel

        !      Calculation of HCO(T2)

        do ige = 1, ngazem
          coefe(ige) = zero
        enddo
        coefe(ico) = 1.d0
        do icha = 1, ncharm
          f1mc(icha) = zero
          f2mc(icha) = zero
        enddo

        t2        = tfuel(iel)
        xhco      = cs_coal_thconvers1(coefe, f1mc, f2mc,  t2)

        !      Calculation of HH2(T2)

        do ige = 1, ngazem
          coefe(ige) = zero
        enddo
        coefe(ihy) = 1.d0
        do icha = 1, ncharm
          f1mc(icha) = zero
          f2mc(icha) = zero
        enddo

        t2        = tfuel(iel)
        xhh2      = cs_coal_thconvers1(coefe, f1mc, f2mc, t2)

        !       Calculation of HH2O(T1)

        do ige = 1, ngazem
          coefe(ige) = zero
        enddo
        coefe(ih2o) = 1.d0
        do icha = 1, ncharm
          f1mc(icha) = zero
          f2mc(icha) = zero
        enddo
        t1        = cpro_temp(iel)
        xhh2o     = cs_coal_thconvers1(coefe, f1mc, f2mc, t1)

        do icla=1,nclacp
          if (cvara_xck(icla)%p(iel) .gt. epsicp) then
            gamhet = crom(iel)*cpro_ghh2oa(icla)%p(iel)                      &
                     * ((cvara_xck(icla)%p(iel))**(2.d0/3.d0) +              &
                    2.d0/3.d0*(cvar_xck(icla)%p(iel)-cvara_xck(icla)%p(iel)) &
                      /(cvara_xck(icla)%p(iel))**(1.d0/3.d0))
          else
            gamhet = 0.d0
          endif
          smbrs(iel) = smbrs(iel)                                          &
                      -gamhet                                              &
                       *(wmole(ico)*xhco+ wmole(ihy)*xhh2                  &
                                         -wmole(ih2o)*xhh2o)               &
                       /wmolat(iatc) *cell_f_vol(iel)

        enddo

      enddo

    endif

    !   Drying

    if (ippmod(iccoal) .eq. 1) then

      do icla=1,nclacp

        numcha = ichcor(icla)

        call field_get_val_s(ivarfl(isca(ixwt(icla))), cvar_xwtcl)
        call field_get_val_prev_s(ivarfl(isca(ixwt(icla))), cvara_xwtcl)
        call field_get_val_s(igmsec(icla),cpro_csec)
        call field_get_val_s(itemp2(icla),cpro_temp2)
        call field_get_val_s(ix2(icla),cpro_x2)

        do iel = 1, ncel

          !  Calculation of H(H2O) at T2

          do ige = 1, ngazem
            coefe(ige) = zero
          enddo
          coefe(ih2o) = 1.d0
          do icha = 1, ncharm
            f1mc(icha) = zero
            f2mc(icha) = zero
          enddo

          t2 = cpro_temp2(iel)
          mode      = -1
          call cpthp1(mode, hh2ov, coefe, f1mc, f2mc, t2)

          !  Contribution to explicit balance

          if (      cvara_xwtcl(iel).gt. epsicp           &
              .and. xwatch(numcha) .gt. epsicp) then

            aux = crom(iel)*cpro_csec(iel)             &
                 *(cvar_xwtcl(iel)/cpro_x2(iel))       &
                 *(1.d0                    /xwatch(numcha))         &
                 *hh2ov

          else
            aux = 0.d0
          endif

          smbrs(iel) = smbrs(iel) - aux*cell_f_vol(iel)

        enddo

      enddo

    endif

    deallocate(cvar_xck, cvara_xck, cvar_xch, cvar_xnp)

    deallocate(cpro_t2, cpro_gmhet)
    if (ihtco2 .eq. 1) then
      deallocate(cpro_ghco2a)
    endif
    if (ihth2o .eq. 1) then
      deallocate(cpro_ghh2oa)
    endif
    if (ippmod(iccoal) .eq. 1) then
      deallocate(cvar_xwt)
    endif

  endif
endif

if (ieqnox .eq. 1 .and. imdnox.eq.0 .and. ntcabs .gt. 1) then

  ! Source terms on Y_HCN and Y_NO

  if (ivar.eq.isca(iyhcn) .or. ivar.eq.isca(iyno)) then

    !  Pointer source terms
    call field_get_val_s(ighcn1,cpro_exp1)
    call field_get_val_s(ighcn2,cpro_exp2)
    call field_get_val_s(ignoth,cpro_exp3)

    !  Mass molar

    wmhcn = wmole(ihcn)
    wmno  = 0.030d0
    wmo2  = wmole(io2 )

    if (ivar.eq.isca(iyhcn)) then

    !  Source term HCN

      if (vcopt%iwarni.ge.1) then
        write(nfecra,1000) chaine(1:8)
      endif

      ! Arrays of pointers containing the fields values for each class
      ! (loop on cells outside loop on classes)
      allocate(cvara_xck(nclacp), cvara_xch(nclacp))
      allocate(cpro_gmhet(nclacp))
      allocate(cpro_gmdv1(nclacp), cpro_gmdv2(nclacp))

      do icla = 1,nclacp
        call field_get_val_prev_s(ivarfl(isca(ixck(icla))), cvara_xck(icla)%p)
        call field_get_val_prev_s(ivarfl(isca(ixch(icla))), cvara_xch(icla)%p)
        call field_get_val_s(igmdv1(icla), cpro_gmdv1(icla)%p)
        call field_get_val_s(igmdv2(icla), cpro_gmdv2(icla)%p)
        call field_get_val_s(igmhet(icla), cpro_gmhet(icla)%p)
      enddo

      do iel=1,ncel
        wmel=cpro_mmel(iel)
        xo2= cpro_yox(iel)*wmel/wmo2

        aux = cell_f_vol(iel)*crom(iel)                      &
             *(cpro_exp2(iel)+cpro_exp1(iel)                &
             *cvara_yno(iel)                             &
             *cpro_mmel(iel)                           &
             /wmno)

        smbrs(iel)  = smbrs(iel)  - aux*cvara_var(iel)
        rovsdt(iel) = rovsdt(iel) + aux

        do icha=1,ncharb
          gmdev1(icha)=0.d0
          gmdev2(icha)=0.d0
          gmhet (icha)=0.d0
        enddo

        do icla=1,nclacp

          icha = ichcor(icla)

          gmdev1(icha) = gmdev1(icha)                              &
               +cpro_gmdv1(icla)%p(iel)                   &
               *crom(iel)                                 &
               *cvara_xch(icla)%p(iel)
          gmdev2(icha) = gmdev2(icha)                              &
               +cpro_gmdv2(icla)%p(iel)                   &
               *crom(iel)                                 &
               *cvara_xch(icla)%p(iel)
          gmhet(icha) = gmhet(icha)                                &
               +cpro_gmhet(icla)%p(iel)                   &
               *crom(iel)                                 &
               *cvara_xck(icla)%p(iel)**(2.d0/3.d0)

        enddo

        do icha=1,ncharb
          !  % of pure nitrogen in the coal

          aux = -cell_f_vol(iel)*fn(icha)*wmhcn/(wmole(in2)/2.d0)        &
                            *(qpr(icha)*(gmdev1(icha)+gmdev2(icha)))
          if(xo2.gt.0.03d0) then
            aux=aux-cell_f_vol(iel)*fn(icha)*wmhcn/(wmole(in2)/2.d0)     &
                               * (1.d0-qpr(icha)*y2ch(icha))         &
                                / (1-y2ch(icha))*gmhet(icha)         &
                                * (1.d0-xashch(icha))
          endif
          smbrs(iel)  = smbrs(iel) + aux
        enddo

      enddo

      deallocate(cvara_xck, cvara_xch)
      deallocate(cpro_gmhet, cpro_gmdv1, cpro_gmdv2)

    endif

    if (ivar.eq.isca(iyno)) then

      !  Source term NO

      call field_get_val_s(iym1(in2),cpro_yn2)

      if (vcopt%iwarni.ge.1) then
        write(nfecra,1000) chaine(1:8)
      endif

      do iel=1,ncel

        wmel=cpro_mmel(iel)

        aux1 = cell_f_vol(iel)*crom(iel)                 &
              *cpro_exp1(iel)*cvara_yhcn(iel)         &
              *cpro_mmel(iel)/wmhcn
        aux2 = cell_f_vol(iel)*crom(iel)                 &
              *cpro_exp2(iel)*cvara_yhcn(iel)         &
              *wmno/wmhcn
        aux3 = cell_f_vol(iel)*crom(iel)**1.5d0                   &
              *cpro_exp3(iel)                                  &
              *cpro_yn2(iel)

        smbrs(iel)  = smbrs(iel) - aux1*cvara_var(iel)            &
                                 + aux2 + aux3
        rovsdt(iel) = rovsdt(iel) + aux1
      enddo

    endif

  endif

endif

if (ieqnox .eq. 1 .and. imdnox.eq.1 .and. ntcabs .gt. 1) then

  ! Source terms on Y_HCN and Y_NO

  if (ivar.eq.isca(iyhcn) .or. ivar.eq.isca(iyno) .or. ivar.eq.isca(iynh3)    &
    ) then

    ! Arrays of pointers containing the fields values for each class
    ! (loop on cells outside loop on classes)
    allocate(cvara_xck(nclacp), cvara_xch(nclacp))
    allocate(cpro_gmhet(nclacp),cpro_t2(nclacp))
    allocate(cpro_gmdv1(nclacp),cpro_gmdv2(nclacp))

    do icla = 1,nclacp
      call field_get_val_prev_s(ivarfl(isca(ixck(icla))), cvara_xck(icla)%p)
      call field_get_val_prev_s(ivarfl(isca(ixch(icla))), cvara_xch(icla)%p)
      call field_get_val_s(itemp2(icla), cpro_t2(icla)%p)
      call field_get_val_s(igmhet(icla), cpro_gmhet(icla)%p)
      call field_get_val_s(igmdv1(icla), cpro_gmdv1(icla)%p)
      call field_get_val_s(igmdv2(icla), cpro_gmdv2(icla)%p)
    enddo

    ! Pointer Source terms NO gas phase
    call field_get_val_s(ighcn1,cpro_exp1)
    call field_get_val_s(ighcn2,cpro_exp2)
    call field_get_val_s(ignoth,cpro_exp3)
    call field_get_val_s(ignh31,cpro_exp4)
    call field_get_val_s(ignh32,cpro_exp5)
    call field_get_val_s(igrb,  cpro_exprb)

    call field_get_val_s(ifnh3d,cpro_cnorb)
    call field_get_val_s(ifnh3c,cpro_fnoch)
    call field_get_val_s(icnohc,cpro_cnohc)
    call field_get_val_s(ifnohc,cpro_fnohc)
    call field_get_val_s(ifnonh,cpro_fnonh)
    call field_get_val_s(ifnoth,cpro_fnoth)
    call field_get_val_s(icnonh,cpro_cnonh)
    call field_get_val_s(ifnh3d,cpro_fnh3d)
    call field_get_val_s(ifnh3c,cpro_fnh3c)

    ! Pointer on CHx1 and CHx2
    call field_get_val_s(iym1(1),cpro_cyf1)
    call field_get_val_s(iym1(2),cpro_cyf2)

    call field_get_val_s(ifhcnr,cpro_fhcnr)
    call field_get_val_s(ifhcnd,cpro_fhcnd)
    call field_get_val_s(ifhcnc,cpro_fhcnc)

    ! Mass molar

    wmhcn = wmole(ihcn)
    wmno  = 0.030d0
    wmnh3 = wmole(inh3)

    if (ivar.eq.isca(iyhcn)) then

      aux = 0.d0

      do iel = 1, ncel
        cpro_fhcnr(iel) = zero
        cpro_fhcnd(iel) = zero
        cpro_fhcnc(iel) = zero
      enddo

      !  Source term HCN

      if (vcopt%iwarni.ge.1) then
        write(nfecra,1000) chaine(1:8)
      endif

      do iel = 1, ncel

        !  Mass molar of the gas mixture
        wmel=cpro_mmel(iel)


        !  Coefficient of reactions HCN + O2 et HCN + NO
        aux = cell_f_vol(iel)*crom(iel)                               &
             *(cpro_exp2(iel)+cpro_exp1(iel)                    &
             *cvara_yno(iel)                                          &
             *wmel/wmno)

        smbrs(iel)  = smbrs(iel)  - aux*cvara_var(iel)
        rovsdt(iel) = rovsdt(iel) + aux

        !  Reburning ?
        !  Chen's model
        if(irb.eq.1) then

          do icha = 1, ncharb

            ychx =   (cpro_cyf1(iel) * wmel/wmchx1)               &
                   + (cpro_cyf2(iel) * wmel/wmchx2)

            aux = cell_f_vol(iel)*wmhcn*cpro_exprb(iel)           &
                * cvara_yno(iel)*wmel/wmno                           &
                * ychx

            smbrs(iel)  = smbrs(iel)  + aux

            cpro_fhcnr(iel) = cpro_fhcnr(iel) + aux

          enddo

        !  Dimitiou's model
        elseif(irb.eq.2) then

          do icha = 1,ncharb

            !  Reburning by CHx1
            if(cpro_cyf1(iel).gt.0.d0) then

              !  Number of point of the temperature discretization
              do ii = 1,7

                !  We look for the interval teno(ii) < Tgaz < teno(ii+1)
                if(cpro_temp(iel).ge.teno(ii).and.cpro_temp(iel).lt.   &
                   teno(ii+1)) then

                  !  JJ indicates the quotient H/C of the fuel (4=CH4;3=CH3,etc.)
                  do jj = 1,4

                    !  We look for the interval jj < chx1(icha) < jj + 1
                    if(chx1(icha).ge.4.d0) then

                    core1 = ka(4,ii) + ((ka(4,ii+1)- ka(4,ii))/(teno(ii+1)-   &
                            teno(ii))) * (cpro_temp(iel) - teno(ii))
                    core2 = kb(4,ii) + ((kb(4,ii+1)- kb(4,ii))/(teno(ii+1)-   &
                            teno(ii))) * (cpro_temp(iel) - teno(ii))
                    para2 = chi2(ii) + ((chi2(ii+1)-chi2(ii))/(teno(ii+1) -   &
                            teno(ii))) * (cpro_temp(iel) - teno(ii))

                    elseif(chx1(icha).le.1.d0) then

                    core1 = ka(1,ii) + ((ka(1,ii+1)- ka(1,ii))/(teno(ii+1)-   &
                            teno(ii))) * (cpro_temp(iel) - teno(ii))
                    core2 = kb(1,ii) + ((kb(1,ii+1)- kb(1,ii))/(teno(ii+1)-   &
                            teno(ii))) * (cpro_temp(iel) - teno(ii))
                    para2 = chi2(ii) + ((chi2(ii+1)-chi2(ii))/(teno(ii+1) -   &
                            teno(ii))) * (cpro_temp(iel) - teno(ii))

                    elseif (chx1(icha).ge.jj.and.chx1(icha).lt.jj+1) then

                    core1 = ka(jj,ii) + ((ka(jj+1,ii+1)- ka(jj,ii))/          &
                            (teno(ii+1)-teno(ii))) * (cpro_temp(iel)      &
                            - teno(ii))
                    core2 = kb(jj,ii) + ((kb(jj+1,ii+1)- kb(jj,ii))/          &
                            (teno(ii+1)-teno(ii))) * (cpro_temp(iel) -    &
                            teno(ii))
                    para2 = chi2(ii) + ((chi2(ii+1)-chi2(ii))/(teno(ii+1)-    &
                            teno(ii))) * (cpro_temp(iel) - teno(ii))

                    endif

                  enddo

                endif

              enddo

              if (chx1(icha).ge.3.d0) then

                aux = (cell_f_vol(iel)*wmhcn)                                  &
                    * ((core1 + core2) * para2)                                &
                    * (cvar_yno(iel)*crom(iel)/wmno)                           &
                    * (cpro_cyf1(iel)*cpro_rom1(iel)/wmchx1)

              else

                aux = (cell_f_vol(iel)*wmhcn)                                  &
                    * (core1 + core2)                                          &
                    * (cvar_yno(iel)*crom(iel)/wmno)                           &
                    * (cpro_cyf1(iel)*cpro_rom1(iel)/wmchx1)

              endif

              smbrs(iel)  = smbrs(iel)  + aux

              cpro_fhcnr(iel) = cpro_fhcnr(iel) + aux

            elseif (cpro_cyf2(iel).gt.0.d0) then

              !  Reburning by CHx2

              !  Number of points of the temperature discretization
              do ii = 1,7

                !  We look for the interval teno(ii) < Tgaz < teno(ii+1)
                if (cpro_temp(iel).ge.teno(ii).and.cpro_temp(iel).lt.   &
                   teno(ii+1)) then

                  !  JJ indicates the quotient H/C of
                  !           the fuel (4=CH4;3=CH3,etc.)
                  do jj = 1,4

                    !  We look for the interval jj < chx1(icha) < jj + 1
                    if(chx2(icha).ge.4.d0) then

                      core1 = ka(4,ii) + ((ka(4,ii+1)- ka(4,ii))/(teno(ii+1)-   &
                              teno(ii))) * (cpro_temp(iel) - teno(ii))
                      core2 = kb(4,ii) + ((kb(4,ii+1)- kb(4,ii))/(teno(ii+1)-   &
                              teno(ii))) * (cpro_temp(iel) - teno(ii))
                      para2 = chi2(ii) + ((chi2(ii+1)-chi2(ii))/(teno(ii+1) -   &
                              teno(ii))) * (cpro_temp(iel) - teno(ii))

                    elseif(chx2(icha).le.1.d0) then

                      core1 = ka(1,ii) + ((ka(1,ii+1)- ka(1,ii))/(teno(ii+1)-   &
                              teno(ii))) * (cpro_temp(iel) - teno(ii))
                      core2 = kb(1,ii) + ((kb(1,ii+1)- kb(1,ii))/(teno(ii+1)-   &
                              teno(ii))) * (cpro_temp(iel) - teno(ii))
                      para2 = chi2(ii) + ((chi2(ii+1)-chi2(ii))/(teno(ii+1) -   &
                              teno(ii))) * (cpro_temp(iel) - teno(ii))

                    elseif (chx2(icha).ge.jj.and.chx2(icha).lt.jj+1) then

                      core1 = ka(jj,ii) + ((ka(jj+1,ii+1)- ka(jj,ii))/          &
                              (teno(ii+1)-teno(ii))) * (cpro_temp(iel) -    &
                              teno(ii))
                      core2 = kb(jj,ii) + ((kb(jj+1,ii+1)- kb(jj,ii))/          &
                              (teno(ii+1)-teno(ii))) * (cpro_temp(iel) -    &
                              teno(ii))
                      para2 = chi2(ii) + ((chi2(ii+1)-chi2(ii))/(teno(ii+1) -   &
                              teno(ii))) * (cpro_temp(iel) - teno(ii))

                    endif

                  enddo

                endif

              enddo

              if(chx2(icha).ge.3.d0) then

                aux = (cell_f_vol(iel)*wmhcn)                                  &
                    * ((core1 + core2) * para2)                                &
                    * (cvar_yno(iel)*crom(iel)/wmno)                           &
                    * (cpro_cyf2(iel)*cpro_rom1(iel)/wmchx2)

              else

                aux = (cell_f_vol(iel)*wmhcn)                                      &
                    * (core1 + core2)                                          &
                    * (cvar_yno(iel)*crom(iel)/wmno)                           &
                    * (cpro_cyf2(iel)*cpro_rom1(iel)/wmchx2)

              endif

              smbrs(iel)  = smbrs(iel)  + aux

              cpro_fhcnr(iel) = cpro_fhcnr(iel) + aux

            endif

          enddo

        endif

        !  Initialization of variables
        do icha = 1, ncharb

          gmdev1(icha)=0.d0
          gmdev2(icha)=0.d0
          gmhet (icha)=0.d0

        enddo

        do icla = 1, nclacp

          icha   = ichcor(icla)

          mckcl1 = (1.d0-y1ch(icha))*a1ch(icha)               &
                   *exp(-e1ch(icha)/( cs_physical_constants_r &
                                     *cpro_t2(icla)%p(iel)))

          mckcl2 = (1.d0-y2ch(icha))*a2ch(icha)               &
                   *exp(-e2ch(icha)/( cs_physical_constants_r &
                                     *cpro_t2(icla)%p(iel)))

          !  Forming rate of the first pyrolisis reaction
          gmdev1(icha) = gmdev1(icha)                                          &
               +cpro_gmdv1(icla)%p(iel)                               &
               *crom(iel)                                             &
               *cvara_xch(icla)%p(iel)

          !  Forming rate of the second pyrolisis reaction
          gmdev2(icha) = gmdev2(icha)                                          &
               +cpro_gmdv2(icla)%p(iel)                               &
               *crom(iel)                                             &
               *cvara_xch(icla)%p(iel)

          if (cvara_xck(icla)%p(iel) .gt. epsicp) then
            !  Reaction rate of the heterogeneous combustion
            gmhet(icha) = gmhet(icha)                                          &
               +cpro_gmhet(icla)%p(iel)                                             &
               *crom(iel)                                             &
               *(cvara_xck(icla)%p(iel)*((1.d0/(mckcl2/mckcl1+1.d0))*yhcnc1(icha)   &
                   +(1.d0/(mckcl1/mckcl2+1.d0))*yhcnc2(icha)))**(2.d0/3.d0)
          endif

        enddo

        !  Modified source term (new model of NOx)

        do icha = 1, ncharb

          !  Release of HCN during devolatilization
          aux = -cell_f_vol(iel)*(gmdev1(icha)*yhcnle(icha)                       &
                +gmdev2(icha)*yhcnlo(icha))

          !   Release of HCN during the heterogeneous combustion according
          !   to the value repnck(icha)

          aux = aux-cell_f_vol(iel)*gmhet(icha)

          smbrs(iel)  = smbrs(iel) + aux

          !  Source terms displaying
          cpro_fhcnd(iel) = cpro_fhcnd(iel)             &
          -cell_f_vol(iel)*(gmdev1(icha)*yhcnle(icha)+gmdev2(icha)*yhcnlo(icha))

          cpro_fhcnc(iel) = cpro_fhcnc(iel)              &
          -cell_f_vol(iel)*gmhet(icha)

        enddo

      enddo

    endif

    if (ivar.eq.isca(iynh3)) then

      aux = 0.d0

      do iel = 1, ncel
        cpro_fnh3d(iel) = zero
      enddo

      !  Source term NH3

      if (vcopt%iwarni.ge.1) then
        write(nfecra,1000) chaine(1:8)
      endif

      do iel = 1, ncel

        !  Mass molar of the gaseous mixture
        wmel = cpro_mmel(iel)

        !  Coefficient of reactions NH3 + O2 and NH3 + NO
        aux  =   cell_f_vol(iel)*crom(iel)                                &
             * (cpro_exp4(iel) + cpro_exp5(iel)                         &
             *   cvara_yno(iel)*wmel/wmno)

        smbrs(iel)  = smbrs(iel)  - aux*cvara_var(iel)
        rovsdt(iel) = rovsdt(iel) + aux

        !  Initialization of variables
        do icha = 1, ncharb

          gmdev1(icha)=0.d0
          gmdev2(icha)=0.d0
          gmhet (icha)=0.d0

        enddo

        do icla = 1, nclacp

          icha = ichcor(icla)

          !  Forming rate of the first pyrolisis reaction
          gmdev1(icha) = gmdev1(icha)                                          &
               +cpro_gmdv1(icla)%p(iel)                               &
               *crom(iel)                                             &
               *cvara_xch(icla)%p(iel)

          !  Forming rate of the second pyrolisis reaction
          gmdev2(icha) = gmdev2(icha)                                          &
               +cpro_gmdv2(icla)%p(iel)                               &
               *crom(iel)                                             &
               *cvara_xch(icla)%p(iel)

        enddo

        do icha = 1, ncharb

          !  Release of NH3 during the devolatization.
          aux = -cell_f_vol(iel)*(gmdev1(icha)*ynh3le(icha)                       &
                +gmdev2(icha)*ynh3lo(icha))

          smbrs(iel)  = smbrs(iel) + aux

          !  Source terms displaying
          cpro_fnh3d(iel) = cpro_fnh3d(iel)             &
          -cell_f_vol(iel)*(gmdev1(icha)*ynh3le(icha)+gmdev2(icha)*ynh3lo(icha))

          cpro_fnh3c(iel) = zero

        enddo

      enddo

    endif

    if (ivar.eq.isca(iyno)) then

      call field_get_val_s(iym1(in2),cpro_yn2)

      do iel = 1, ncel
        cpro_cnorb(iel) = zero
        cpro_fnoch(iel) = zero
      enddo

      !  Source term NO

      if (vcopt%iwarni.ge.1) then
        write(nfecra,1000) chaine(1:8)
      endif

      do iel = 1, ncel

        !  Mass molar of the gaseous mixture
        wmel=cpro_mmel(iel)

        !  Coefficient of reaction HCN + NO
        aux1 = cell_f_vol(iel)*crom(iel)                              &
              *cpro_exp1(iel)*cvara_yhcn(iel)                      &
              *wmel/wmhcn

        cpro_cnohc(iel) = aux1*cvara_var(iel)

        !  Coefficient of reaction HCN + O2
        aux2 = cell_f_vol(iel)*crom(iel)                              &
              *cpro_exp2(iel)*cvara_yhcn(iel)                      &
              *wmno/wmhcn

        cpro_fnohc(iel) = aux2

        !  Coefficient of thermal NO
        aux3 = cell_f_vol(iel)*crom(iel)**1.5d0                           &
              *cpro_exp3(iel)                                               &
!FIXME       Pourquoi la fraction massique d'azote n'a ete pas transforme dans une
!       fraction molaire ?
              *cpro_yn2(iel)

        cpro_fnoth(iel) = aux3

        !  Coefficient of reaction NH3 + O2 --> NO + ...
        aux4 = cell_f_vol(iel)*crom(iel)                               &
              *cpro_exp4(iel)*cvara_ynh3(iel)                          &
              *wmno/wmnh3

        cpro_fnonh(iel) = aux4

        !  Coefficient of reaction NH3 + NO --> N2 + ...
        aux5 = cell_f_vol(iel)*crom(iel)                               &
              *cpro_exp5(iel)*cvara_ynh3(iel)                          &
              *wmel/wmnh3

        cpro_cnonh(iel) = aux5*cvara_var(iel)

        !  Reburning ?
        !  Chen's model
        if (irb.eq.1) then

          do icha = 1,ncharb

             ychx = (cpro_cyf1(iel) * wmel/wmchx1)                       &
                  + (cpro_cyf2(iel) * wmel/wmchx2)

             aux = cell_f_vol(iel)*wmhcn*cpro_exprb(iel)                   &
                 * cvara_yno(iel) * wmel/wmno  * ychx

             smbrs(iel)  = smbrs(iel)  - aux

             cpro_cnorb(iel) = cpro_cnorb(iel) + aux

          enddo

        !  Dimitiou's model
        else if (irb.eq.2) then

          do icha = 1, ncharb

            !  Reburning by CHx1
            if (cpro_cyf1(iel).gt.0.d0) then

              !  Number of point of the temperature discretization
              do ii = 1,7

                !  We look for the interval teno(ii) < Tgaz < teno(ii+1)
                if(cpro_temp(iel).ge.teno(ii).and.cpro_temp(iel).lt.   &
                   teno(ii+1)) then

                  !  JJ indicates the quotient H/C of the fuel (4=CH4;3=CH3,etc.)
                  do jj = 1,4

                    !  We look for the interval jj < chx1(icha) < jj + 1
                    if(chx1(icha).ge.4.d0) then

                      core1 = ka(4,ii) + ((ka(4,ii+1)- ka(4,ii))/(teno(ii+1)-   &
                              teno(ii))) * (cpro_temp(iel) - teno(ii))
                      core2 = kb(4,ii) + ((kb(4,ii+1)- kb(4,ii))/(teno(ii+1)-   &
                              teno(ii))) * (cpro_temp(iel) - teno(ii))
                      core3 = kc(4,ii) + ((kc(4,ii+1)- kc(4,ii))/(teno(ii+1)-   &
                              teno(ii))) * (cpro_temp(iel) - teno(ii))
                      para2 = chi2(ii) + ((chi2(ii+1)-chi2(ii))/(teno(ii+1) -   &
                              teno(ii))) * (cpro_temp(iel) - teno(ii))

                    elseif(chx1(icha).le.1.d0) then

                      core1 = ka(1,ii) + ((ka(1,ii+1)- ka(1,ii))/           &
                              (teno(ii+1)-teno(ii))) * (cpro_temp(iel)  -   &
                              teno(ii))
                      core2 = kb(1,ii) + ((kb(1,ii+1)- kb(1,ii))/(teno(ii+1)-   &
                              teno(ii))) * (cpro_temp(iel) - teno(ii))
                      core3 = kc(1,ii) + ((kc(1,ii+1)- kc(1,ii))/           &
                              (teno(ii+1)-teno(ii))) * (cpro_temp(iel) -    &
                              teno(ii))
                      para2 = chi2(ii) + ((chi2(ii+1)-chi2(ii)) /           &
                              (teno(ii+1) - teno(ii))) * (cpro_temp(iel) -  &
                              teno(ii))

                    elseif (chx1(icha).ge.jj.and.chx1(icha).lt.jj+1) then

                      core1 = ka(jj,ii) + ((ka(jj+1,ii+1)- ka(jj,ii))/      &
                              (teno(ii+1)-teno(ii))) * (cpro_temp(iel) -    &
                              teno(ii))
                      core2 = kb(jj,ii) + ((kb(jj+1,ii+1)- kb(jj,ii))/      &
                              (teno(ii+1)-teno(ii))) * (cpro_temp(iel) -    &
                              teno(ii))
                      core3 = kc(jj,ii) + ((kc(jj+1,ii+1)- kc(jj,ii))/      &
                              (teno(ii+1)-teno(ii))) * (cpro_temp(iel) -    &
                              teno(ii))
                      para2 = chi2(ii) + ((chi2(ii+1)-chi2(ii))/            &
                              (teno(ii+1) - teno(ii))) * (cpro_temp(iel) -  &
                              teno(ii))

                    endif

                  enddo

                endif

              enddo

              if (chx1(icha).ge.3.d0) then

                auxrb1 = (cell_f_vol(iel)*wmno)                                &
                       * ((core1 + core2 + core3) * para2)                     &
                       * (cpro_cyf1(iel)*cpro_rom1(iel)/wmchx1)         &
                       * (cvar_yno(iel) * crom(iel)/wmno)

              else

                auxrb1 = (cell_f_vol(iel)*wmno)                                &
                       * (core1 + core2 + core3)                               &
                       * (cpro_cyf1(iel)*cpro_rom1(iel)/wmchx1)         &
                       * (cvar_yno(iel) * crom(iel)/wmno)

              endif

              smbrs(iel)  = smbrs(iel)  - auxrb1

              cpro_cnorb(iel) = cpro_cnorb(iel) + auxrb1

            !  Reburning by CHx2
            elseif(cpro_cyf2(iel).gt.0.d0) then

              !  Number of point of the temperature discretization
              do ii = 1,7

                !  We look for the interval teno(ii) < Tgaz < teno(ii+1)
                if(cpro_temp(iel).ge.teno(ii).and.cpro_temp(iel).lt.   &
                   teno(ii+1)) then

                  !  JJ incates the quotient H/C of the fuel (4=CH4;3=CH3,etc.)
                  do jj = 1,4

                    !  We look for the interval jj < chx1(icha) < jj + 1
                    if(chx2(icha).ge.4.d0) then

                      core1 = ka(4,ii) + ((ka(4,ii+1)- ka(4,ii))/               &
                              (teno(ii+1)-teno(ii))) * (cpro_temp(iel) -    &
                              teno(ii))
                      core2 = kb(4,ii) + ((kb(4,ii+1)- kb(4,ii))/               &
                              (teno(ii+1)-teno(ii))) * (cpro_temp(iel) -    &
                              teno(ii))
                      core3 = kc(4,ii) + ((kc(4,ii+1)- kc(4,ii))/               &
                              (teno(ii+1)-teno(ii))) * (cpro_temp(iel) -    &
                              teno(ii))
                      para2 = chi2(ii) + ((chi2(ii+1)-chi2(ii))/                &
                              (teno(ii+1) - teno(ii))) * (cpro_temp(iel) -  &
                              teno(ii))

                    elseif(chx2(icha).le.1.d0) then

                      core1 = ka(1,ii) + ((ka(1,ii+1)- ka(1,ii))/(teno(ii+1)-   &
                              teno(ii))) * (cpro_temp(iel) -                &
                              teno(ii))
                      core2 = kb(1,ii) + ((kb(1,ii+1)- kb(1,ii))/               &
                              (teno(ii+1)-teno(ii))) * (cpro_temp(iel) -    &
                              teno(ii))
                      core3 = kc(1,ii) + ((kc(1,ii+1)- kc(1,ii))/               &
                              (teno(ii+1)-teno(ii))) * (cpro_temp(iel) -    &
                              teno(ii))
                      para2 = chi2(ii) + ((chi2(ii+1)-chi2(ii)) /               &
                              (teno(ii+1) - teno(ii))) * (cpro_temp(iel) -  &
                              teno(ii))
                    elseif (chx2(icha).ge.jj.and.chx2(icha).lt.jj+1) then

                      core1 = ka(jj,ii) + ((ka(jj+1,ii+1)- ka(jj,ii))/          &
                              (teno(ii+1)-teno(ii))) * (cpro_temp(iel) -    &
                              teno(ii))
                      core2 = kb(jj,ii) + ((kb(jj+1,ii+1)- kb(jj,ii))/          &
                              (teno(ii+1)-teno(ii))) * (cpro_temp(iel) -    &
                              teno(ii))
                      core3 = kc(jj,ii) + ((kc(jj+1,ii+1)- kc(jj,ii))/          &
                              (teno(ii+1)-teno(ii))) * (cpro_temp(iel) -    &
                              teno(ii))
                      para2 = chi2(ii) + ((chi2(ii+1)-chi2(ii))/                &
                              (teno(ii+1) - teno(ii))) * (cpro_temp(iel) -  &
                              teno(ii))

                    endif

                  enddo

                endif

              enddo

              if(chx2(icha).ge.3.d0) then

                auxrb2 = (cell_f_vol(iel)*wmno)                                &
                       * ((core1 + core2 + core3) * para2)                     &
                       * (cpro_cyf2(iel)*cpro_rom1(iel)/wmchx2)         &
                       * (cvar_yno(iel) * crom(iel)/wmno)

              else

                auxrb2 = (cell_f_vol(iel)*wmno)                                &
                       * (core1 + core2 + core3)                               &
                       * (cpro_cyf2(iel)*cpro_rom1(iel)/wmchx2)         &
                       * (cvar_yno(iel) * crom(iel)/wmno)

              endif

              smbrs(iel)  = smbrs(iel)  - auxrb2

              cpro_cnorb(iel) = cpro_cnorb(iel) + auxrb2

            endif

          enddo

        endif


        !  Initialization
        do icha=1,ncharb

          gmhet (icha)=0.d0

        enddo

        do icla=1,nclacp

          icha   = ichcor(icla)

          mckcl1 = (1.d0-y1ch(icha))*a1ch(icha)               &
                   *exp(-e1ch(icha)/( cs_physical_constants_r &
                                     *cpro_t2(icla)%p(iel)))

          mckcl2 = (1.d0-y2ch(icha))*a2ch(icha)               &
                   *exp(-e2ch(icha)/( cs_physical_constants_r &
                                     *cpro_t2(icla)%p(iel)))

          !  Reaction rate of the heterogeneous combustion
          if (cvara_xck(icla)%p(iel) .gt. epsicp) then
          gmhet(icha) = gmhet(icha)                                            &
               +cpro_gmhet(icla)%p(iel)                               &
               *crom(iel)                                             &
               *(cvara_xck(icla)%p(iel)*((1.d0/(mckcl2/mckcl1+1.d0))*ynoch1(icha)   &
                   +(1.d0/(mckcl1/mckcl2+1.d0))*ynoch2(icha)))**(2.d0/3.d0)
          endif

        enddo

        !  Coefficient of released NO during the heterogeneous combustion

        do icha=1,ncharb

          auxhet = -cell_f_vol(iel)*gmhet(icha)
          cpro_fnoch(iel) = cpro_fnoch(iel) + auxhet


          smbrs(iel)  = smbrs(iel) - aux1*cvara_var(iel)                       &
                                   - aux5*cvara_var(iel)                       &
                                   + aux2 + aux3 + aux4 + auxhet

          rovsdt(iel) = rovsdt(iel) + aux1 + aux5

        enddo

      enddo

    endif

    deallocate(cvara_xck, cvara_xch)
    deallocate(cpro_gmhet,cpro_t2)
    deallocate(cpro_gmdv1,cpro_gmdv2)

  endif

endif

!--------
! Formats
!--------

 1000 format(' Specific physic source term for the variable ', a8, /)

!----
! End
!----

! Deallocate work arrays
deallocate(w1, w3, w4, w5, tfuel, stat=iok1)

if (iok1 > 0) then
  write(nfecra,*) ' cs_coal_scast: memory deallocation error'
  call csexit(1)
endif

end subroutine
