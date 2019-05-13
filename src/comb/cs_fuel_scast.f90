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
! Function:
! ---------
!> \file cs_fuel_scast.f90
!> \brief Specific physic routine: fuel oil flame.
!>   We indicate the source terms for a scalar PP
!>   on a step time
!>
!>
!> \warning   The treatment of source terms is different from
!>            the treatment in ustssc.f
!>
!> we solve \f$ rovsdt D(var) = smbrs \f$
!>
!> \f$ rovsdt \f$ and \f$ smbrs \f$ already contain eventual user source term.
!>  So they have to be incremented and be erased
!>
!> For stability reasons, we only add in rovsdt positive terms.
!>  There is no stress for smbrs.
!>
!> In the case of a source term in \f$ cexp + cimp \varia \f$ we must
!> write:
!>        \f[ smbrs  = smbrs  + cexp + cimp\cdot \varia\f]
!>        \f[ rovsdt = rovsdt + Max(-cimp,0)\f]
!>
!> We provide here rovsdt and smbrs (they contain rho*volume)
!>    smbrs in \f$kg\cdot [variable] \cdot s^{-1}\f$:
!>     ex : for velocity               \f$kg\cdot m \cdot s^{-2}\f$
!>          for temperatures           \f$kg \cdot [degres] \cdot s^{-1}\f$
!>          for enthalpies             \f$J \cdot s^{-1} \f$
!>    rovsdt in \f$kg \cdot s^{-1}\f$
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!            ARGUMENTS
!______________________________________________________________________________.
!  mode           name          role
!______________________________________________________________________________!
!> \param[in]     iscal         scalar number
!> \param[in,out] smbrs         second explicit member
!> \param[in,out] rovsdt        implicit diagonal part
!______________________________________________________________________________!

subroutine cs_fuel_scast &
 ( iscal  ,                                                       &
   smbrs  , rovsdt )
!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use entsor
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
use cs_fuel_incl
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
integer          ivar , iel, icla , numcla
integer          imode  , iesp
integer          itermx,nbpauv,nbrich,nbepau,nberic
integer          nbarre,nbimax,nbpass
integer          iterch
integer          keyccl, f_id

double precision aux, rhovst
double precision rom
double precision hfov
double precision ho2,hco,xesp(ngazem),t2mt1
double precision gmech,gmvap,gmhet
double precision xxco,xxo2,xxco2,xxh2o
double precision xkp,xkm,t0p,t0m
double precision aux1 , aux2 , aux3 , w1
double precision anmr,tauchi,tautur
double precision sqh2o , x2 , wmhcn , wmno ,wmo2
double precision err1mx,err2mx
double precision errch,fn,qpr
double precision auxmax,auxmin
double precision ymoy
double precision fn0,fn1,fn2,anmr0,anmr1,anmr2
double precision lnk0p,l10k0e,lnk0m,t0e,xco2eq,xcoeq,xo2eq
double precision xcom,xo2m,xkcequ,xkpequ,xden
double precision smbrs1
double precision, dimension(:), pointer :: vp_x, vp_y, vp_z
double precision, dimension(:,:), pointer :: vdc
double precision, dimension(:), pointer :: crom
double precision, dimension(:), pointer :: cvara_k, cvara_ep
double precision, dimension(:), pointer :: cvar_yfolcl, cvara_yfolcl
double precision, dimension(:), pointer :: cvara_yno, cvara_yhcn
double precision, dimension(:), pointer :: cvara_var
type(pmapper_double_r1), dimension(:), allocatable :: cvara_yfol
double precision, dimension(:), pointer :: taup
double precision, dimension(:), pointer :: smbrsh1, rovsdth1
double precision, dimension(:,:), pointer :: vel
double precision, dimension(:,:), pointer ::  vg_lim_pi
double precision, dimension(:), pointer :: cpro_temp1, cpro_temp2, cpro_rom2
double precision, dimension(:), pointer :: cpro_diam2, cpro_cgev, cpro_cght
double precision, dimension(:), pointer :: cpro_chgl, cpro_yox, cpro_yco2
double precision, dimension(:), pointer :: cpro_yco, cpro_yh2o, cpro_rom1
double precision, dimension(:), pointer :: cpro_exp1, cpro_exp2, cpro_exp3
double precision, dimension(:), pointer :: cpro_mmel, cpro_yn2

type(var_cal_opt) :: vcopt

!===============================================================================
! 1. Initialization
!===============================================================================

! --- Number of the scalar to treat: iscal

! --- Number of the variable associated to the scalar to treat iscal
ivar = isca(iscal)
call field_get_val_prev_s(ivarfl(isca(iscal)), cvara_var)

! --- Name of the variable associated to scalar to treat iscal
call field_get_label(ivarfl(ivar), chaine)

! --- Number of the physic bulks (Cf cs_user_boundary_conditions)
call field_get_val_s(icrom, crom)

! --- Gas phase temperature

call field_get_val_s(itemp1, cpro_temp1)
call field_get_val_s(irom1, cpro_rom1)

call field_get_val_prev_s(iym1(io2), cpro_yox)
call field_get_val_prev_s(iym1(ico2), cpro_yco2)
call field_get_val_prev_s(iym1(ico), cpro_yco)
call field_get_val_prev_s(iym1(ih2o), cpro_yh2o)

call field_get_val_v(ivarfl(iu), vel)

! Key id of the coal scalar class
call field_get_key_id("scalar_class", keyccl)

! source term for gas enthalpy due to inter-phase fluxes
call field_get_val_s_by_name('x_h_c_exp_st',smbrsh1)
call field_get_val_s_by_name('x_h_c_imp_st',rovsdth1)

call field_get_key_struct_var_cal_opt(ivarfl(ivar), vcopt)

!===============================================================================
! 2. Taking into account the source terms for the relative particles
!    in the particles classes
!===============================================================================

!===============================================================================
! 2.1 Source term for the enthalpies
!===============================================================================

if ( ivar .ge. isca(ih2(1)) .and. ivar .le. isca(ih2(nclafu)) ) then

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

  ! index of the droplet class
  call field_get_key_int(ivarfl(ivar), keyccl, numcla)

  call field_get_val_prev_s(ivarfl(isca(iyfol(numcla))), cvara_yfolcl)
  call field_get_val_s(irom2(numcla),cpro_rom2)
  call field_get_val_s(idiam2(numcla),cpro_diam2)
  call field_get_val_s(itemp2(numcla),cpro_temp2)
  call field_get_val_s(igmeva(numcla),cpro_cgev)
  call field_get_val_s(igmhtf(numcla),cpro_cght)
  call field_get_val_s(ih1hlf(numcla),cpro_chgl)

  !       The variable is the liquid enthalpy for the mixture mass
  !       The interfacial flux contribute to the variation of the liquid
  !       enthalpy
  !       The vapor takes away its enthalpy
  !       flux = cpro_cgev(iel)
  !       massic enthalpy reconstructed from ehgaze(ifov )
  !       at the drop temperature
  !       The heterogeneous oxidation contains an input flux of O2
  !       an output flux of CO
  !       The net flux is the carbon flux
  !       fluxIN  = 16/12 * cpro_cght(iel)
  !       fluxOUT = 28/12 * cpro_cght(iel)
  !       Input enthalpy reconstructed from ehgaze(IO2 )
  !       at the surrounding gas temperature
  !       Output enthalpy reconstructed from ehgaze(ico )
  !       at the grain temperature

  imode = -1
  do iel = 1, ncel

    if (cvara_yfolcl(iel) .gt. epsifl) then

      rom = crom(iel)

      do iesp = 1, ngazem
        xesp(iesp) = zero
      enddo

      xesp(ifov) = 1.d0
      call cs_fuel_htconvers1(imode,hfov,xesp,cpro_temp2(iel))

      xesp(ifov) = zero
      xesp(io2)  = 1.d0
      call cs_fuel_htconvers1(imode,ho2 ,xesp,cpro_temp1(iel))

      xesp(io2)  = zero
      xesp(ico)  = 1.d0
      call cs_fuel_htconvers1(imode,hco,xesp,cpro_temp2(iel))

      t2mt1 = cpro_temp2(iel)-cpro_temp1(iel)

      gmech = -cpro_chgl(iel)*t2mt1
      gmvap = cpro_cgev(iel)*hfov*t2mt1
      gmhet = 16.d0/12.d0*cpro_cght(iel)*ho2                    &
             -28.d0/12.d0*cpro_cght(iel)*hco

      smbrs(iel) = smbrs(iel) + (gmech+gmvap+gmhet)*rom*volume(iel)
      ! FIXME better time stepping?
      smbrsh1(iel) = smbrsh1(iel) -(gmech+gmvap+gmhet)*rom*volume(iel)
      rhovst = ( cpro_chgl(iel)                                 &
                -cpro_cgev(iel)*hfov )/cp2fol                   &
              *rom*volume(iel)
      rovsdt(iel) = rovsdt(iel) +  max(zero,rhovst)

    endif

  enddo

! --> Source terme for the liquid mass

elseif ( ivar .ge. isca(iyfol(1))     .and.                       &
         ivar .le. isca(iyfol(nclafu))        ) then

  if (vcopt%iwarni.ge.1) then
    write(nfecra,1000) chaine(1:8)
  endif

  ! index of the droplet class
  call field_get_key_int(ivarfl(ivar), keyccl, numcla)

  call field_get_val_s(irom2(numcla),cpro_rom2)
  call field_get_val_s(idiam2(numcla),cpro_diam2)
  call field_get_val_s(itemp2(numcla),cpro_temp2)
  call field_get_val_s(igmeva(numcla),cpro_cgev)
  call field_get_val_s(igmhtf(numcla),cpro_cght)
  call field_get_val_s(ih1hlf(numcla),cpro_chgl)

  do iel = 1, ncel

    t2mt1 =  cpro_temp2(iel)-cpro_temp1(iel)
    gmvap = -cpro_cgev(iel)*t2mt1
    gmhet = -cpro_cght(iel)

    smbrs(iel) = smbrs(iel)                                       &
         - crom(iel)*volume(iel)*(gmvap+gmhet)
    if ( cvara_var(iel).gt.epsifl ) then
      rhovst = crom(iel)*volume(iel)*(gmvap + gmhet)              &
              / cvara_var(iel)
    else
      rhovst = 0.d0
    endif
    rovsdt(iel) = rovsdt(iel) + max(zero,rhovst)

  enddo

endif

!===============================================================================
! 2.2 Particle velocity source terms
!===============================================================================

if (i_comb_drift.eq.1) then

  call field_get_name(ivarfl(ivar), fname)

  ! index of the particle class
  call field_get_key_int(ivarfl(ivar), keyccl, icla)

  if (icla.ge.1) then

    ! Taup
    write(name,'(a,i2.2)')'nd_fuel_' ,icla!FIXME change name?
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
        ! TODO PP
        smbrs1 = 0.d0

        ! relaxation to drop velocity
        smbrs1 = crom(iel)*volume(iel)*(1.d0/taup(iel)+smbrs1)                &
               *(vel(1,iel)+vdc(1,iel)+vg_lim_pi(1, iel)-vp_x(iel))

        smbrs(iel) = smbrs(iel) + smbrs1
        rovsdt(iel) = rovsdt(iel) + crom(iel)*volume(iel)/taup(iel)

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
        ! TODO PP
        smbrs1 = 0.d0
        ! relaxation to drop velocity
        smbrs1 = crom(iel)*volume(iel)*(1.d0/taup(iel)+smbrs1)                &
               *(vel(2,iel)+vdc(2, iel)+vg_lim_pi(2, iel)-vp_y(iel))
        smbrs(iel) = smbrs(iel) + smbrs1
        rovsdt(iel) = rovsdt(iel) + crom(iel)*volume(iel)/taup(iel)

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
        ! TODO PP
        smbrs1 = 0.d0

        ! relaxation to drop velocity
        smbrs1 = crom(iel)*volume(iel)*(1.d0/taup(iel)+smbrs1)                &
               *(vel(3,iel)+vdc(3, iel)+vg_lim_pi(3, iel)-vp_z(iel))

        smbrs(iel) = smbrs(iel) + smbrs1
        rovsdt(iel) = rovsdt(iel) + crom(iel)*volume(iel)/taup(iel)

      enddo !on iel

    endif !on fname

  endif

endif !on icoal drift

!===============================================================================
! 3. Taking into account source terms for relative variables in the mixture
!===============================================================================

if (ivar .eq. isca(ihgas)) then

  ! source terms from particles (convection and enthalpy drived by mass fluxes)
  do iel = 1, ncel
    smbrs(iel) = smbrs(iel) + smbrsh1(iel)
    rovsdt(iel) = rovsdt(iel) + rovsdth1(iel)
  enddo

  ! terme source de rayonnement pour l'enthalpie du gaz
  ! serait mieux dans ... raysca
  if (iirayo.ge.1) then
    !TODO PP
  endif
endif


if ( ivar .eq. isca(ifvap) ) then

  if (vcopt%iwarni.ge.1) then
    write(nfecra,1000) chaine(1:8)
  endif

  do icla = 1, nclafu

    call field_get_val_s(ivarfl(isca(iyfol(icla))), cvar_yfolcl)
    call field_get_val_prev_s(ivarfl(isca(iyfol(icla))), cvara_yfolcl)
    call field_get_val_s(itemp2(icla),cpro_temp2)
    call field_get_val_s(igmeva(icla),cpro_cgev)

    do iel = 1, ncel

      t2mt1 = cpro_temp2(iel)-cpro_temp1(iel)
      if ( cvara_yfolcl(iel) .gt. epsifl ) then
        gmvap = -cpro_cgev(iel)*t2mt1*cvar_yfolcl(iel)        &
                / cvara_yfolcl(iel)
      else
        gmvap = -cpro_cgev(iel)*t2mt1
      endif

      smbrs(iel) = smbrs(iel)                                     &
                 + gmvap*crom(iel)*volume(iel)
    enddo

  enddo

! --> Source term for the C tracer ex heterogeneous reaction

elseif ( ivar .eq. isca(if7m) ) then

  if (vcopt%iwarni.ge.1) then
    write(nfecra,1000) chaine(1:8)
  endif

  do icla = 1, nclafu

    call field_get_val_s(ivarfl(isca(iyfol(icla))), cvar_yfolcl)
    call field_get_val_prev_s(ivarfl(isca(iyfol(icla))), cvara_yfolcl)
    call field_get_val_s(igmhtf(icla),cpro_cght)

    do iel = 1, ncel
      if (cvara_yfolcl(iel) .gt. epsifl) then
        smbrs(iel) = smbrs(iel)                                        &
             -crom(iel)*cpro_cght(iel)*volume(iel)                 &
                                *cvar_yfolcl(iel)                      &
                                /cvara_yfolcl(iel)
      else
        smbrs(iel) = smbrs(iel)                                        &
                    -crom(iel)*cpro_cght(iel)*volume(iel)
      endif

    enddo

  enddo

endif

! --> Source term for the variance of the tracer 4 (Air)

if ( ivar.eq.isca(ifvp2m) ) then

  if (vcopt%iwarni.ge.1) then
    write(nfecra,1000) chaine(1:8)
  endif

  ! ---- Calculation of the source the explicit and implicit source terms
  !      relative to interfacial exchanges between phases

  call cs_fuel_fp2st &
 !==================
 ( iscal  ,                                                        &
   smbrs  , rovsdt )

endif


! --> Source term for CO2

if ( ieqco2 .ge. 1 ) then

  if ( ivar.eq.isca(iyco2) ) then

    if (vcopt%iwarni.ge.1) then
      write(nfecra,1000) chaine(1:8)
    endif

    call field_get_val_prev_s(ivarfl(ik), cvara_k)
    call field_get_val_prev_s(ivarfl(iep), cvara_ep)

    ! Arrays of pointers containing the fields values for each class
    ! (loop on cells outside loop on classes)
    allocate(cvara_yfol(nclafu))
    do icla = 1, nclafu
      call field_get_val_prev_s(ivarfl(isca(iyfol(icla))), cvara_yfol(icla)%p)
    enddo

    ! ---- Contribution of the interfacial source term to the explicit and implicit balances

    ! Oxidation of CO
    ! ===============

    !  Dryer Glassman : XK0P in (mol/m3)**(-0.75) s-1
    !          XK0P = 1.26D10
    !          XK0P = 1.26D7 * (1.1)**(NTCABS)
    !          IF ( XK0P .GT. 1.26D10 ) XK0P=1.26D10
    !          T0P  = 4807.D0
    !  Howard : XK0P en [(moles/m3)**(-0.75) s-1]
    !          XK0P = 4.11D9
    !          T0P  = 15090.D0
    !  Westbrook & Dryer

    lnk0p = 23.256d0
    t0p  = 20096.d0

    !  Hawkin et Smith Purdue University Engeneering Bulletin, i
    !  Research series 108 vol 33, n 3n 1949
    !  Kp = 10**(4.6-14833/T)
    !  Equilibrum constant in partial pressure [atm           !]
    !  XKOE is the decimal log of the pre-exponential constant
    !  TOE  is NOT an activation temperature  ... it remains a lg(e)
    !  to go back in Kc and to use concentrations [moles/m3]
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
   ! Precision for convergence
   errch = 1.d-8

   do iel = 1, ncel

     xxco  = cpro_yco(iel)/wmole(ico)           &
            *cpro_rom1(iel)
     xxo2  = cpro_yox(iel)/wmole(io2)           &
            *cpro_rom1(iel)
     xxco2 = cpro_yco2(iel)/wmole(ico2)         &
            *cpro_rom1(iel)
     xxh2o = cpro_yh2o(iel)/wmole(ih2o)         &
            *cpro_rom1(iel)

     xxco  = max(xxco ,zero)
     xxo2  = max(xxo2 ,zero)
     xxco2 = max(xxco2,zero)
     xxh2o = max(xxh2o,zero)
     sqh2o = sqrt(xxh2o)

     xkp = exp(lnk0p-t0p/cpro_temp1(iel))
     xkm = exp(lnk0m-t0m/cpro_temp1(iel))

     xkpequ = 10.d0**(l10k0e-t0e/cpro_temp1(iel))
     xkcequ = xkpequ                                              &
             /sqrt(8.32d0*cpro_temp1(iel)/1.015d5)

     !        initialization per transported state

     anmr  = xxco2
     xcom  = xxco + xxco2
     xo2m  = xxo2 + 0.5d0*xxco2

     if ( cpro_temp1(iel) .gt. 1200.d0 ) then

     !           Search for the equilibrum state
     !           Iterative search with convergence control
     !            (to preserve parallelism on the meshes)
     !            on the number of reaction mols which separate
     !            the state before reaction (as calculated by Cpcym)
     !            of the equilibrum state
     !           anmr must be confined between 0 and Min(xcom,2.*xo2m)
     !           We look for the solution by dichotomy

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

       if ( xo2m.gt.1.d-6) then
         do while ( iterch.lt.itermx .and. fn2.gt.errch )
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

         if ( iterch .ge. itermx) then
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

     if ( xco2eq.gt.xxco2 ) then
       !           oxidation
       xden = xkp*sqh2o*(xxo2)**0.25d0
     else
       !           dissociation
       xden = xkm
     endif
     if ( xden .ne. 0.d0 ) then

       tauchi = 1.d0/xden
       tautur = cvara_k(iel)/cvara_ep(iel)

       x2 = 0.d0
       do icla = 1, nclafu
         x2 = x2 + cvara_yfol(icla)%p(iel)
       enddo

       !    We transport CO2

       smbrs(iel)  = smbrs(iel)                                   &
                    +wmole(ico2)/cpro_rom1(iel)        &
         * (xco2eq-xxco2)/(tauchi+tautur)                         &
         * (1.d0-x2)                                              &
         * volume(iel) * crom(iel)

       w1 = volume(iel)*crom(iel)/(tauchi+tautur)
       rovsdt(iel) = rovsdt(iel) +   max(w1,zero)

     else
       rovsdt(iel) = rovsdt(iel) + 0.d0
       smbrs(iel)  = smbrs(iel)  + 0.d0
     endif

   enddo

   deallocate(cvara_yfol)

   if(irangp.ge.0) then
     call parcpt(nberic)
     call parmax(err1mx)
     call parcpt(nbpass)
     call parcpt(nbarre)
     call parcpt(nbarre)
     call parcmx(nbimax)
   endif

   write(nfecra,*) ' Max Error = ', err1mx
   write(nfecra,*) ' no Points   ', nberic, nbarre, nbpass
   write(nfecra,*) ' Iter max number ', nbimax


  endif

endif


! --> Source term for HCN and NO: only from the second
!                                   iteration

if ( ieqnox .eq. 1 .and. ntcabs .gt. 1) then

  if ( ivar.eq.isca(iyhcn) .or. ivar.eq.isca(iyno) ) then

    call field_get_val_s(ighcn1,cpro_exp1)
    call field_get_val_s(ighcn2,cpro_exp2)
    call field_get_val_s(ignoth,cpro_exp3)
    call field_get_val_s(immel, cpro_mmel)
    call field_get_val_s(iym1(in2),cpro_yn2)

    ! QPR= %N released during the evaporation/average volatile materials
    !          rate

    qpr = 1.3d0

    ! YMOY = % output vapor

    ymoy = 0.7d0

    ! Azote in the fuel oil

    fn = 0.015

    ! Molar mass
    wmhcn = wmole(ihcn)
    wmno  = 0.030d0
    wmo2  = wmole(io2)

    if ( ivar.eq.isca(iyhcn) ) then

      !        Source term HCN

      if (vcopt%iwarni.ge.1) then
        write(nfecra,1000) chaine(1:8)
      endif

      call field_get_val_prev_s(ivarfl(isca(iyno)), cvara_yno)

      auxmin = 1.d+20
      auxmax =-1.d+20

      do iel=1,ncel

        xxo2 = cpro_yox(iel)                       &
              *cpro_mmel(iel)/wmo2

        aux = volume(iel)*crom(iel)                                &
             *( cpro_exp2(iel)                                  &
               +cpro_exp1(iel)*cvara_yno(iel)                   &
                                 *cpro_mmel(iel)/wmno )

        smbrs(iel)  = smbrs(iel)  - aux*cvara_var(iel)
        rovsdt(iel) = rovsdt(iel) + aux

        gmvap = 0.d0
        gmhet = 0.d0

        do icla=1,nclafu

          call field_get_val_s(itemp2(icla),cpro_temp2)
          call field_get_val_s(igmeva(icla),cpro_cgev)
          call field_get_val_s(igmhtf(icla),cpro_cght)

          gmvap = gmvap                                           &
                 + crom(iel)*cpro_cgev(iel)                   &
                  *(cpro_temp2(iel)-cpro_temp1(iel))

          gmhet = gmhet                                           &
                 +crom(iel)*cpro_cght(iel)

        enddo
        if ( xxo2 .gt. 0.03d0 ) then
          aux = -volume(iel)*fn*wmhcn/(wmole(in2)/2.d0)           &
                *( qpr*gmvap+(1.d0-qpr*ymoy)/(1.d0-ymoy)*gmhet )
        else
          aux = -volume(iel)*fn*wmhcn/(wmole(in2)/2.d0)           &
                            *(qpr*gmvap)
        endif
        smbrs(iel)  = smbrs(iel) + aux

      enddo

    endif

    if ( ivar.eq.isca(iyno) ) then

      !        Source term NO

      if (vcopt%iwarni.ge.1) then
        write(nfecra,1000) chaine(1:8)
      endif

      call field_get_val_prev_s(ivarfl(isca(iyhcn)), cvara_yhcn)

      do iel=1,ncel

        aux1 = volume(iel)*crom(iel)                     &
              *cpro_exp1(iel)*cvara_yhcn(iel)            &
              *cpro_mmel(iel)/wmhcn
        aux2 = volume(iel)*crom(iel)                     &
              *cpro_exp2(iel)*cvara_yhcn(iel)            &
              *wmno/wmhcn
        aux3 = volume(iel)*crom(iel)**1.5d0              &
              *cpro_exp3(iel)                            &
              *cpro_yn2(iel)

        smbrs(iel)  = smbrs(iel) - aux1*cvara_var(iel)   &
                               + aux2 + aux3
        rovsdt(iel) = rovsdt(iel) + aux1
      enddo

    endif

  endif

endif

!--------
! Formats
!--------

 1000 format(' Specific physic source term for the variable '  &
       ,a8,/)

!----
! End
!----

return

end subroutine
