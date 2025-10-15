!-------------------------------------------------------------------------------

! This file is part of code_saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2025 EDF S.A.
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
!> \file rayso.f90
!>
!-------------------------------------------------------------------------------
!> \brief Compute reflexion and transmission
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role
!______________________________________________________________________________!
!> \param[in]   pioc        Albedo of simple diffusion for cloud (water)
!> \param[in]   piaero      Albedo of simple diffusion for aerosol
!> \param[in]   gasym       Asymmetry factor for clouds
!> \param[in]   gaero       Asymmetry factor for aerosols
!> \param[in]   tauc        Optical depth for clouds
!> \param[in]   taua        Optical depth for aersols
!> \param[out]  ref         Reflexion
!> \param[out]  tra         Transmission
!> \param[in]   epsc        clipping threshold
!> \param[in]   dqqv        Optical depth for Water vapor
!> \param[in]   mui         mui (Gauss quadrature)
!> \param[in]   muzero_cor  mu0 corrected
!_______________________________________________________________________________

subroutine reftra  &
    (pioc, piaero, gasym, gaero, tauc, taua, &
    ref, tra, epsc,dqqv, mui, muzero_cor)

  !===========================================================================

  implicit none

  ! Arguments

  double precision, intent(in) :: pioc, piaero,  gasym, gaero
  double precision, intent(in) :: tauc, taua,  dqqv, epsc
  double precision, intent(inout) :: ref, tra
  double precision, intent(in) :: mui, muzero_cor

  ! Local
  double precision ::gas, fas, kt, gama1, gama2, tln
  double precision :: drt, extlnp, extlnm
  double precision :: pic, tau
  double precision :: gama3,gama4,a1,a2,ktmu,expmuo,exnmuo,rnum,tnum
  logical          :: glob = .true. ! solve global and direct
  !===========================================================================

  tau = tauc +taua + dqqv
  ! For 0 optical depth
  if (tau .lt. epsc) then
    ref = 0.d0
    tra = 1.d0
  else

    ! We solve Global and Direct
    if (glob) then
      ! Pure diffusion atmosphere (pioc=1)
      if (pioc.ge.(1.d0-epsc)) then !TODO check .and. (taua .le. epsc))
        gama1=(1.d0-gasym)/(2.d0*mui)
        ref = gama1*tau/(1.d0+gama1*tau)
        tra = 1.d0/(1.d0+gama1*tau)

      else
        pic =(pioc*tauc+piaero*taua)/tau
        ! Pure absorbing atmosphere (pioc=0)
        if (pic .lt. epsc) then
          gama1=1.d0/mui
          ref = 0.d0
          tra = dexp(-gama1*tau)
        else

          gas=(pioc*tauc*gasym+piaero*taua*gaero) / (pic*tau)
          ! Only forward diffusion
          if(gas.ge.(1.d0-epsc)) then
            gama1=(1.d0-pic)/mui
            ref=0.d0
            tra=dexp(-gama1*tau)
          else
            ! Joseph correction
            fas=gas*gas
            tau=(1.d0-pic*fas)*tau
            pic=pic*(1.d0-fas)/(1.d0-pic*fas)
            gas=(gas-fas)/(1.d0-fas)
            gama1=(1.d0-pic*(1.d0+gas)/2.d0)/mui
            gama2=pic*(1.d0-gas)/(2.d0*mui)
            kt=dsqrt(gama1*gama1-gama2*gama2)
            tln=kt*tau
            extlnp=dexp(tln)
            extlnm=dexp(-tln)
            drt=(kt+gama1)*extlnp+(kt-gama1)*extlnm
            ref=gama2*(extlnp-extlnm)/drt
            tra=2.d0*kt/drt
          endif
        endif

      endif

    ! We solve Diffuse and Direct
    else

      ! Pure diffusion atmosphere (pioc=1)
      if (pioc.ge.(1.d0-epsc)) then
        gama1=0.5d0 * (1.d0-gasym)/mui
        gama3= 0.5d0 * (1.d0-gasym)
        exnmuo=dexp(-tau/muzero_cor)
        ref=(gama1*tau+(gama3-gama1*muzero_cor)*(1.d0-exnmuo))/(1.d0+gama1*tau)
        tra = 1.d0-ref

      else
        pic =(pioc*tauc+piaero*taua)/tau
        ! Pure absorbing atmosphere (pioc=0)
        if (pic .lt. epsc) then
          ref = 0.d0
          tra = dexp(-tau/muzero_cor)
        else

          gas=(pioc*tauc*gasym+piaero*taua*gaero) / (pic*tau)
          ! Only forward diffusion
          if(gas.ge.(1.d0-epsc)) then
            gama1=(1.d0-pic)/mui
            gama2=0.d0
            gama3=0.5d0 * (1.d0-gas)
            gama4=0.5d0 * (1.d0+gas)
            expmuo=dexp(tau/muzero_cor*(1.d0-gama1*muzero_cor))
            exnmuo=dexp(-tau/muzero_cor*(1.d0+gama1*muzero_cor))

            ref=gama3*pic*(1.d0-exnmuo)/(1.d0+gama1*muzero_cor)

            tra= dexp(-tau/muzero_cor)*(1.d0-pic*gama4*expmuo &
              /(1.d0-gama1*muzero_cor))
          else
            ! Joseph correction
            fas=gas*gas
            tau=(1.d0-pic*fas)*tau
            pic=pic*(1.d0-fas)/(1.d0-pic*fas)
            gas=(gas-fas)/(1.d0-fas)

            gama1=(1.d0-pic*(1.d0+gas)/2.d0)/mui
            gama2=pic*(1.d0-gas)/(2.d0*mui)
            gama3=0.5d0*(1.d0 - 3.d0* gas*muzero_cor * mui)
            gama4=0.5d0*(1.d0 + 3.d0* gas*muzero_cor * mui)
            a1=gama1*gama4+gama2*gama3
            a2=gama1*gama3+gama2*gama4

            kt=dsqrt(gama1*gama1-gama2*gama2)
            ktmu=kt*muzero_cor
            tln=kt*tau
            extlnp=dexp(tln)
            extlnm=dexp(-tln)
            expmuo=dexp(tau/muzero_cor)
            exnmuo=dexp(-tau/muzero_cor)
            drt=(1.d0-ktmu*ktmu)*((kt+gama1)*extlnp+(kt-gama1)*extlnm)

            rnum=(1.d0-ktmu)*(a2+kt*gama3)*extlnp &
               - (1.d0+ktmu)*(a2-kt*gama3)*extlnm &
               - 2.d0*kt*(gama3-a2*muzero_cor) * exnmuo
            ref=pic*rnum/drt

            tnum=(1.d0+ktmu)*(a1+kt*gama4)*extlnp &
                -(1.d0-ktmu)*(a1-kt*gama4)*extlnm &
                -2.d0*kt*(gama4+a1*muzero_cor)*expmuo
            tra=exnmuo*(1.d0-pic*tnum/drt)
          endif
        endif
      endif

    endif
  endif

end subroutine reftra

!-------------------------------------------------------------------------------
!> \brief Compute solar fluxes for both clear and cloudy atmosphere following
!> Lacis and Hansen (1974). The multiple diffusion is taken into account by an
!> addition method and overlapping between water vapor and liquid water with k
!> distribution method.
!> Some improvements from original version concerns:
!> - introduction of cloud fraction with hazardous recovering
!> - introduction of aerosol diffusion in the same way as for cloud droplets
!>   but with specific optical properties for aerosols.
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role
!______________________________________________________________________________!
!> \param[in]   ivertc      index of vertical profile
!> \param[in]   k1          index for ground level
!> \param[in]   kmray       vertical levels number for radiation
!> \param[in]   heuray      Universal time (Hour)
!> \param[in]   imer1       sea index
!> \param[in]   albe        albedo
!> \param[in]   qqv         optical depth for water vapor (0,zqq)
!> \param[in]   qqqv        idem for intermediate levels
!> \param[in]   qqvinf      idem qqv but for altitude above 11000m
!> \param[in]   zqq         vertical levels (interfaces)
!> \param[in]   zray        vertical levels (volumes)
!> \param[in]   qvray       specific humidity for water vapor
!> \param[in]   qlray       specific humidity for liquid water
!> \param[in]   fneray      cloud fraction
!> \param[in]   romray      air density
!> \param[in]   preray      pressure
!> \param[in]   temray      temperature
!> \param[out]  fos         global downward solar flux at the ground
!> \param[out]  rayst       flux divergence of solar radiation
!_______________________________________________________________________________

subroutine rayso  &
 (ivertc, k1, kmray, heuray, imer1, albe,        &
  qqv, qqqv, qqvinf, zqq,                        &
  zray, qvray, qlray, fneray,                    &
  romray, preray, temray, fos, rayst, ncray)

!===============================================================================
! Module files
!===============================================================================

use optcal
use cstphy
use parall
use ppincl
use cs_c_bindings
use mesh
use field
use atincl, only: kmx, nbmett, sigc, squant, xlat, xlon, sold, &
                  solu, piaero_o3,piaero_h2o, &
                  black_carbon_frac,zaero, gaero_o3, gaero_h2o, &
                  aod_h2o_tot, aod_o3_tot, cp_a, cp_v, rad_atmo_model
use cstnum, only: epzero, pi

!===============================================================================

implicit none

procedure() :: reftra

! Arguments

integer ivertc, k1, kmray, imer1
double precision albe, heuray, fos
double precision qqv(kmx+1), qqqv(kmx+1), qqvinf, zqq(kmx+1)
double precision qlray(kmx), fneray(kmx), zray(kmx)
double precision qvray(kmx), preray(kmx)
double precision rayst(kmx), romray(kmx)
double precision temray(kmx)
double precision ncray(kmx)

! Local variables

integer i,k,n,l,iaer,iaero_top,iaero_topr
integer itop,ibase
integer          ifac, iz1, iz2, f_id, c_id, iel
integer          k0
double precision za, muzero, muzero_cor,fo,m,mbar,mbarh2o,rabar,rabar2,rbar
double precision rabarc,rbarc, refx, trax, refx0, trax0
double precision qqvtot,y,ystar
double precision yp1,ystarp1
double precision zqm1,zq,x,xstar
double precision zqp1,xp1,xstarp1
double precision rrbar,rrbar2s
double precision tauctot,wh2ol,rm,req,deltaz
double precision pioc,zbas,dud1
double precision gasym,drt,tln,drtt1,dtrb1
double precision kn(8),pkn(8),dqqv
double precision cphum,qureel
double precision taua(kmx+1),tauatot
double precision rabara,rbara
double precision gama1,gama2,kt,gas,fas
double precision omega, var, zent
double precision cpvcpa
double precision dzx, dy
! For postprocessing
logical          is_active

integer, dimension(:), pointer :: itypfb
type(c_ptr) :: c_itypfb

double precision, allocatable:: fabsh2o(:),fabso3(:),tauc(:)
double precision, allocatable:: tau(:,:),pic(:,:),ref(:,:)
double precision, allocatable:: reft(:,:),trat(:,:)
double precision, allocatable:: refb(:,:),upwf(:,:)
double precision, allocatable:: fabso3c(:,:),tra(:,:)
double precision, allocatable:: dow(:,:)
double precision, allocatable:: dowf(:,:),atln(:,:),absn(:,:)
double precision, allocatable :: ckup_sir_f(:), ckdown_sir_r(:), ckdown_sir_f(:)
double precision, allocatable :: ckup_suv_f(:), ckdown_suv_r(:), ckdown_suv_f(:)
double precision, allocatable :: w0_suv(:), w0_sir(:), g_apc_suv(:), g_apc_sir(:)
double precision, allocatable:: fnebmax(:),fneba(:)

double precision, allocatable, dimension(:,:) :: dowd, trad, trard
double precision, allocatable, dimension(:) :: dffsh2o, dffso3
double precision, allocatable, dimension(:) :: ddfsh2o, ddfso3
double precision, allocatable, dimension(:) :: dddfsh2o, dddfso3
double precision, allocatable, dimension(:) :: drfs, ddfs

double precision, allocatable, dimension(:) ::  dfsh2o, ufsh2o
double precision, allocatable, dimension(:) ::  dfso3, ufso3
double precision, allocatable, dimension(:,:) ::  ufso3c
double precision, allocatable, dimension(:) ::  dfs, ufs
double precision, dimension(:,:), pointer :: bpro_rad_inc
double precision, dimension(:,:), pointer :: cpro_ck_up
double precision, dimension(:,:), pointer :: cpro_ck_down
double precision, dimension(:,:), pointer :: cpro_w0
double precision, dimension(:,:), pointer :: cpro_gapc
! For computing albedo PIC, PIOC
double precision epsc
double precision pioco3,pioch2o,gasymo3,gasymh2o
double precision pic_o3(kmx+1),pic_h2o(kmx+1),gco3(kmx+1),gch2o(kmx+1)
double precision pioco3_1, pioco3_2,pioch2o_1, pioch2o_2
double precision rayst_h2o(kmx),rayst_o3(kmx)
double precision tabara, tabarc
! 5 minor gas (NO, NO2, CO2, O2, CO)
double precision Tmg,amg(5),bmg(5),cmg(5),dmg(5),umg(5)
! For black carbon, 12 bands
double precision piocv(12), copioc(12), omega0(12), beta1(12)
double precision beta2(12), beta3(12),beta4(12),copioc20(12)
double precision nu0, dm0, dm, coeff_E_o3(12), coeff_E_h2o(12)
double precision pioco3C, pioch2oC
double precision tauao3(kmx+1) , tauah2o(kmx+1)
double precision corp
double precision mui,tauapc,ckapcd,ckapcf,picapc,gapc
double precision ck_aero_h2of,ck_aero_h2od,ck_aero_o3f,ck_aero_o3d
! data for pkn and kn distribution
data kn/4.d-5,0.002d0,0.035d0,0.377d0,1.95d0,9.40d0,44.6d0,190.d0/
data pkn/0.647d0,0.0698d0,0.1443d0,0.0584d0,0.0335d0,0.0225d0,0.0158d0,0.0087d0/

! Data from Chuang 2002 for calculation of cloud SSA  taking into account black carbon
data beta1/0.2382803,0.2400113,0.2471480,0.2489583,0.2542476,0.2588392,0.2659081,0.2700860,0.2783093,0.2814346,0.2822860,0.1797007/
data beta2/0.2940957,0.2936845,0.2880274,0.2871209,0.2824498,0.2775943,0.2698008,0.265296,0.2564840,0.2535739,0.2487382,0.1464709/
data beta3/61.83657,58.25082,52.79042,50.06907,45.75322,42.43440,37.03823,32.32349,25.99426,20.05043,12.76966,3.843661/
data beta4/574.2111,565.0809,519.0711,494.0088,448.3519,409.9063,348.9051,297.9909,233.7397,175.4385,112.8208,39.24047/
data omega0/5189239.d-11,2261712.d-11,1264190.d-11,9446845.d-11,6090293.d-12,3794524.d-12,1735499.d-12,1136807.d-12,2261422.d-12,&
     1858815.d-11,5551822.d-9,2325124.d-7/
data coeff_E_o3/0.0,0.0,0.0,0.0,0.029193795335,0.045219606416,0.16880411522,0.186215099078,0.5705673839,0.0,0.0,0.0/
data coeff_E_h2o/0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.27068650951091927, 0.6844296138881737, 0.044883876600907306/

! Data for calculation of Tmg - Transmission function for minor gases
data amg/0.0721d0,0.0062d0,0.0326d0,0.0192d0,0.0003d0/
data bmg/377.890d0,243.670d0,107.413d0,166.095d0,476.934d0/
data cmg/0.5855d0,0.4246d0,0.5501d0,0.4221d0,0.4892d0/
data dmg/3.1709d0,1.7222d0,0.9093d0,0.7186d0,0.1261d0/
data umg/390.0d0,0.075d0,0.28d0,1.6d0,209500.d0/

!========================================================================

allocate(fabsh2o(kmx+1),fabso3(kmx+1),tauc(kmx+1))
! Note ref(0) = albedo
allocate(tau(kmx+1,8),pic(kmx+1,8),ref(0:kmx+1,8))
allocate(reft(kmx+1,8),trat(kmx+1,8))
! Note refb(0) = albedo
allocate(refb(0:kmx+1,8),upwf(kmx+1,8))
! Note tra(0) for the ground
allocate(fabso3c(kmx+1,2),tra(0:kmx+1,8))
allocate(dow(kmx+1,8))
allocate(dowf(kmx+1,8),atln(kmx+1,8),absn(kmx+1,8))
allocate(fnebmax(kmx+1),fneba(kmx+1))

allocate(ufso3c(kmx+1,2))
allocate(dowd(kmx+1,8), trad(kmx+1,8), trard(0:kmx+1,8))

allocate(dfsh2o(kmx+1), ufsh2o(kmx+1))
allocate(dfso3(kmx+1), ufso3(kmx+1))
allocate(dfs(kmx+1), ufs(kmx+1))
allocate(ckup_sir_f(kmx), ckdown_sir_r(kmx), ckdown_sir_f(kmx))
allocate(ckup_suv_f(kmx), ckdown_suv_r(kmx), ckdown_suv_f(kmx))
allocate(w0_sir(kmx), w0_suv(kmx), g_apc_sir(kmx), g_apc_suv(kmx))

allocate(dffsh2o(kmx+1), dffso3(kmx+1))
allocate(ddfsh2o(kmx+1), ddfso3(kmx+1))
allocate(dddfsh2o(kmx+1), dddfso3(kmx+1))
allocate(drfs(kmx+1), ddfs(kmx+1))

! 1 - local initializations
! ===========================

cpvcpa = cp_v / cp_a

iaer = 1 ! has aerosols, always, remove
ibase = 0
epsc=1.d-8

do k = 1,kmray
  w0_sir(k) = 0.d0
  w0_suv(k) = 0.d0
  g_apc_sir(k) = 0.d0
  g_apc_suv(k) = 0.d0
  ! TODO useful ?
  if (qlray(k).lt.epsc) then
    qlray(k) = 0.d0
    fneray(k)=0.d0
  endif
enddo

do k = 1, kmx+1
  fabsh2o(k) = 0.d0
  fabso3(k) = 0.d0
  tauc(k) = 0.d0
  taua(k) = 0.d0
  tauao3(k)=0.d0
  tauah2o(k)=0.d0
  fnebmax(k) = 0.d0
  gco3(k)=0.d0
  gch2o(k)=0.d0
  pic_o3(k)=0.d0
  pic_h2o(k)=0.d0
  dfs(k) = 0.d0
  ddfs(k) = 0.d0
  drfs(k) = 0.d0
  dffsh2o(k) =0.d0
  dffso3(k) = 0.d0
  ddfsh2o(k) = 0.d0
  ddfso3(k) = 0.d0
  dddfsh2o(k) = 0.d0
  dddfso3(k) = 0.d0
  dfsh2o(k) = 0.d0
  ufsh2o(k) = 0.d0
  dfso3(k) = 0.d0
  ufso3(k) = 0.d0

  if(iaer.ge.1) then
    fneba(k) = 1.d0
  endif
  do n = 1, 2
    fabso3c(k,n) = 0.d0
  enddo
  do n = 1, 8
    tau(k,n) = 0.d0
    pic(k,n) = 0.d0
    atln(k,n) = 0.d0
    absn(k,n) = 0.d0
    ! Note for slice l->l+1
    ! reflexion stored at level l
    ! transmission stored at level l
    ref(k,n) = 0.d0
    reft(k,n) = 0.d0
    refb(k,n) = 0.d0
    trat(k,n) = 1.d0
    tra(k,n) = 1.d0
    upwf(k,n) = 0.d0
    dow(k,n) = 0.d0
    dowf(k,n) = 0.d0
  enddo
enddo

do k = 1, kmx+1
  do n = 1, 8
    trad(k,n) = 1.d0
    dowd(k,n) = 1.d0
  enddo
enddo
refx=0.d0
trax=0.d0
refx0=0.d0
trax0=0.d0
!initialisation variables for multiple diffusion
drt=0.d0
gas=0.d0
fas=0.d0
gama1=0.d0
gama2=0.d0
kt=0.d0
tln=0.d0
! Leighton 1980
! id of the top of the aerosol layer
iaero_top=0
! remaining aerosol quantities above this level
iaero_topr=0
k0 = k1-1

!initialisation Chuang calculations for  black carbon in droplets
dm0=20.d0 !micrometres
! Volume fraction of Black Carbon (module)
nu0=1.d-8
!Initialisation of data used for the calculation of SSA using Chuang 2002
pioco3C=0.d0
pioch2oC=0.d0
do k = 1, 12
  piocv(k)=0.d0
  copioc(k)=0.d0
  copioc20(k)=0.d0
enddo

!  2 - calculation for muzero and solar constant fo
!  ===================================================
!        muzero = cosin of zenithal angle
!  (corrected to take curvature of the Earth)
!        fo = solar constant in watt/m2

!
!  careful : 0. h < heuray < 24. h
!  ---------

qureel = float(squant)
call raysze(xlat, xlon, qureel, heuray, imer1, albe, za, muzero_cor, omega, fo)
! if muzero is negative, it is night and solar radiation is not
! computed
muzero = dcos(za)

if (muzero.gt.epzero) then

  ! Optical air mass
  ! Corrected value for very low angle
  ! cf. Kasten, F., Young, A.T., 1989. Revised optical air mass tables and approximation formula.
  m = 1.d0/muzero_cor

  ! m coefficient for O3 (5/3 or 1.9 for LH74)

  ! as we suppressed m and mbar in dzx and dzy we can eventually keep 1.9 for O3 to tested
  mbar = 1.9d0
  ! For H2O mbar is 5/3
  mbarh2o = 5.d0 / 3.d0
  ! for LH74 quadrature method in two-stream
  mui=1.d0/dsqrt(3.d0)

  !  3 -  albedos for O3 and Rayleigh diffusion

  ! Note LH74 equation (18)
  rabar = 0.219d0/(1.d0 + 0.816d0*muzero)
  rabar2 = 0.144d0
  ! Note LH74: eq (15)
  rbar = rabar + (1.d0 - rabar)*(1.d0 - rabar2)*albe/(1.d0 - rabar2*albe)
  rrbar = 0.28d0/(1.d0 + 6.43d0*muzero)
  rrbar2s = 0.0685d0

  !  4 - addition of one level for solar radiation

  qqvtot = qqvinf + qqv(kmray)
  qqv(kmray+1) = qqvtot - qqvinf

  ! Transmission for minor gases
  Tmg = 1.d0
  do i = 1, 5
    Tmg = Tmg* (1.d0 - (amg(i) * m * umg(i)) / &
      ((1.d0+bmg(i)*m*umg(i))**cmg(i)+dmg(i)*m*umg(i)))
  enddo

  ! introduction of absorption by minor gases
  fo=fo*Tmg

  ! 5 - Solar radiation calculation for cloudy sky
  ! In order to take into account cloud fraction, multiple diffusion is achieved
  ! for both cloudy (index 1) and clear (index 2) sky

  !  5.1 cloud level determination (top for the top of the higher cloud,
  !  base for the bottom of the lower cloud)

  itop = -1
  do i = kmray, k1, -1
    if(qlray(i).gt.epsc) then
      !FIXME to be coherent with 3D
      ! à supprimer si l'on veut garder la nébulosité fractionnaire(devrait avoir peu deffet sur le cas parisfog)
      ! FIX fneray(i) = 1.d0
      if (itop.eq.-1) then
        itop = i+1
      endif
      ibase = i
    endif
  enddo

  ! if itop = -1, there is no cloud but, nevertheless, it is possible to execute
  ! the adding method for the water vapor (SIR band) only

  if (itop.eq.-1) then
    itop = k1
    ibase = k1
  endif


  ! 5.2 calculation for optical parameters of clouds and aerosols
  ! (single scattering albedo, optical depth, radius, asymmetry factor)

  fnebmax(kmray+1) = 0.d0
  tauctot = 0.d0
  tauatot = 0.d0
  do i = kmray, k1, -1
    if((i.ge.ibase).and.(i.lt.itop)) then
      ! liquid water density in g/m3 in the layers
      wh2ol = 1.d3*(romray(i)*qlray(i))
      ! mean droplet radius in µm
      rm = 30.d0*wh2ol + 2.d0
      !  the max of the mean radius is fixed at 10 µm in the considered domain
      !  and at 2 µm above
      if(i.le.nbmett) then
        rm = min(10.d0,rm)
      else
        rm = min(2.d0, rm)
      endif

      ! Efficient radius
      if (ncray(i).gt.epsc.and.qlray(i).gt.epsc) then
        ! Simplification:
        ! Note Old formula
        ! req = 1.d6*( (3.d0*romray(i)*qlray(i)) /                 &
        !              (4.*pi*1000.*ncray(i)*1.d6))**(1./3.)       &
        !      *dexp(sigc**2)
        req = 1.d3*( (3.d0*romray(i)*qlray(i)) /                 &
          (4.d0*pi*ncray(i)))**(1.d0/3.d0)   &
          *dexp(sigc**2)
      else
        req = 1.5d0 * rm
      endif

      ! Clippling: Climatological limits for effective radius
      if (req .gt. 20.d0) req=20.d0
      if (req .lt. 1.d0) req=1.d0

      ! Cloud optical thickness
      deltaz = zqq(i+1) - zqq(i)

      ! req has to be in µm
      tauc(i) = 1.5d0 * wh2ol * deltaz / req
      tauctot = tauctot + tauc(i)
    else
      tauc(i) = 0.d0
      req = 0.d0
    endif

    ! Calculation of aerosol optical depth AOD
    fneba(i) = 0.d0
    if((iaer.eq.1).and.(zray(i).le.zaero)) then
      iaero_top=max(i+1,iaero_top)
      fneba(i) = 1.d0
      deltaz = zqq(i+1) - zqq(i)

      ! Distribution of AOD on the vertical homogeneous between 0 and zaero
      ! Note, we used a formula based on concentration before v6.2
      tauao3(i) = aod_o3_tot*deltaz/zqq(iaero_top)
      tauah2o(i) = aod_h2o_tot*deltaz/zqq(iaero_top)
    endif
    ! Estimation of the law of the cloud fraction

    ! Calculation of SSA and Asymmetry factor for clouds using- Nielsen 2014
    ! Note only the first two bands are taken, the third only is
    ! approximately 0
    ! Usefull for gas assymmetry, for albedo, either Chuang if BC, or Nielsen
    pioco3_1 = 1.d0 - 33.d-9*req !FIXME MF : ce n'est pas une bonne idée en terme de precision numérique...
    pioco3_2 = 1.d0 - 1.d-7*req

    pioch2o_1 = 0.99999d0 -149.d-7*req
    pioch2o_2 = 0.9985d0 -92.d-5*req

    fnebmax(i) = max(fnebmax(i+1),fneray(i))
    if(black_carbon_frac .gt. epsc) then
      ! Calculation of SSA for clouds taking into account black carbon fraction - Chuang 2002
      dm = req*4.d0/3.d0 !mean diameter
      do k = 5, 12 !5 to 12 because we there is no energy in the
        !4 first spectral band defined by Chuang 2002
        copioc20(k)=omega0(k) &
          + beta1(k)*(1.d0-dexp(-beta3(k)*(black_carbon_frac-nu0))) &
          + beta2(k)*(1.d0-dexp(-beta4(k)*(black_carbon_frac-nu0)))
        copioc(k) = (copioc20(k)*dm/dm0) &
          / (1.d0 + 1.8d0*copioc20(k)*(dm/dm0 - 1.d0))
        piocv(k) = 1.d0 - copioc(k)
        pioco3C = pioco3C+coeff_E_o3(k)*piocv(k)
        pioch2oC = pioch2oC+coeff_E_h2o(k)*piocv(k)
      enddo

      pic_h2o(i)=pioch2oC
      pic_o3(i)=pioco3C
    else

      ! Calculation of SSA and Asymmetry factor for clouds using- Nielsen 2014
      ! Note only the first two bands are taken, the third only is
      ! approximately 0
      pioco3=pioco3_1*0.24d0+pioco3_2*0.76d0

      pioch2o=0.60d0*pioch2o_1+0.40d0*pioch2o_2

      pic_o3(i)=pioco3
      pic_h2o(i)=pioch2o
    endif

    ! Gas asymmetry: Nielsen
    gasymo3 =( 0.868d0 + 14.d-5*req  &
      - 61.d-4*dexp(-0.25*req))*pioco3_1*0.24d0 &
      + ( 0.868d0 + 25.d-5*req  &
      - 63.d-4*dexp(-0.25*req))*pioco3_2*0.76d0

    gasymh2o = ( 0.867d0 + 31.d-5*req  &
      - 78.d-4*dexp(-0.195d0*req))*0.60d0*pioch2o_1 &
      + ( 0.864d0 + 54.d-5*req &
      - 0.133d0*dexp(-0.194d0*req))*0.40d0*pioch2o_2

    gco3(i)=gasymo3
    gch2o(i)=gasymh2o

  enddo

  tauc(kmray+1) = 0.d0

  ! 5.3 O3 absorption in presence of clouds

  ! Calculation of the different albedos for O3 (LH 74)

  ! Asymmetry factor and SSA for liquid water
  gasym=gco3(itop)

  pioc=pic_o3(itop)
  !Calculation for cloudy layers
  call reftra  &
    (pioc, 0.d0, gasym, 0.d0, tauctot, 0.d0, &
    refx, trax, epsc, 0.d0, mui, muzero_cor)

  rabarc=refx
  tabarc=trax
  ! Note: LH74 equation (15)
  rbarc = rabarc+tabarc*tabarc*albe/(1.d0 - rabarc*albe)
  !Calculation for aerosol layers
  call reftra  &
    (0.d0, piaero_o3, 0.d0, gaero_o3, 0.d0, aod_o3_tot, &
    refx, trax, epsc,0.d0, mui, muzero_cor)
  rabara=refx
  tabara=trax

  rbara = rabara+tabara*tabara*albe/(1.d0 - rabara*albe)

  ! in the case there is an aerosol layer above the cloud layer

  if((iaer.eq.1).and.(iaero_top.gt.itop)) then
    itop=iaero_top
    rbar=rbara
    rrbar2s=rabara
  endif

  ! Calculation above the top of the cloud or the aerosol layer

  do l = itop, kmray   !calculation have to start at the first level
    ! ( itop=1  in the case without aerosols nor clouds)
    zq = zqq(l)
    zqp1 = zqq(l+1)

    ! The ozone amount traversed by the direct solar beam in reaching the
    ! l^th layer
    xp1 = m*rayuoz(zqp1)
    x = m*rayuoz(zq)

    ! Calculation of heat and radiation fluxes during cloudy sky
    zbas = zqq(itop)
    ! Path traversed by the diffuse radiation illuminating the l^th layer
    ! from below
    xstarp1 = m*rayuoz(zbas) + mbar*(rayuoz(zbas) - rayuoz(zqp1))
    xstar   = m*rayuoz(zbas) + mbar*(rayuoz(zbas) - rayuoz(zq))

    ! Heat
    ! Note upwf(l) - dow(l) -(upwf(l+1) - dow(l+1)
    dud1=1.d0 / (1.d0-rrbar2s*albe)
    fabso3c(l,1) = muzero*fo*(raysoz(x) - raysoz(xp1)*dud1  &
                            + rbarc*(raysoz(xstarp1) - raysoz(xstar)))

    ! Direct downward radiation
    ddfso3(l) = muzero*fo*(0.647d0-rrbar-raysoz(x))
    ! Diffuse Downward radiation (factor is for Rayleight)
    dddfso3(l) = ddfso3(l) * (dud1 - 1.d0)
    ! Global Downward radiation
    dfso3(l) = ddfso3(l) * dud1
    ! Upward (diffuse) radiation
    ! Note: rbarc correspond to reflexion of cloud + ground
    ! Note 1, correspond to cloudy sky
    ufso3c(l,1) = muzero*fo*(0.647d0-rrbar-raysoz(xstar))*rbarc

    ! Calculation of heat and radiation fluxes during  Clear sky
    zbas = zqq(k1)
    xstarp1 = m*rayuoz(zbas) + mbar*(rayuoz(zbas) - rayuoz(zqp1))
    xstar = m*rayuoz(zbas) + mbar*(rayuoz(zbas) - rayuoz(zq))
    ! heat LH74 (14), ,2) correspond to clear sky
    fabso3c(l,2) = muzero*fo*((raysoz(x) - raysoz(xp1))*dud1 &
                             + albe*dud1*(raysoz(xstarp1) - raysoz(xstar)))
    ! Upward (diffuse) radiation
    ufso3c(l,2) = muzero*fo*(0.647d0-rrbar-raysoz(xstar))*albe*dud1

    ! Summation depending on cloud fraction
    fabso3(l) = fnebmax(k1)*fabso3c(l,1) + (1.d0-fnebmax(k1))*fabso3c(l,2)

    ufso3(l) = fnebmax(k1)*ufso3c(l,1)+ (1.d0-fnebmax(k1))*ufso3c(l,2)
    ! Calculation of absorption coefficient ckup and ckdown
    ! useful for 3D simulation
    dzx = drayuoz(zq)

    ckdown_suv_r(l)=dzxaoz(x,dzx)/(0.647d0-rrbar-raysoz(x))
    ckdown_suv_f(l)=dzxaoz(x,dzx)/(0.647d0-rrbar-raysoz(x))
    ckup_suv_f(l)=dzxaoz(xstar,dzx)/(0.647d0-rrbar-raysoz(xstar))

  enddo

  ! Calculation under the top of the cloud or the aerosol layer, the adding
  ! Method with multiple diffusion is used
  n = 1 !no O3 overlapping

  ! Top boundary conditions
  tra(kmray+1,n)= 1.d0
  ! value from 44km to 11km
  trard(kmray+1,n)= 1.d0
  trad(kmray+1,n)= 1.d0
  ref(kmray+1,n)= 0.d0
  trat(kmray+1,n)=1.d0
  reft(kmray+1,n)=ref(kmray+1,n)

  ! Bottom boundary conditions
  tra(k0,n) = 0.d0
  trard(k0,n) = 0.d0
  ref(k0,n) = albe

  do l = k1, kmray
    tau(l,n)=tauc(l) + tauao3(l)
    gasym=gco3(l)
    pioc=pic_o3(l)
    ! In the cloud layers
    call reftra  &
      (pioc, piaero_o3, gasym, gaero_o3, tauc(l), tauao3(l), &
      refx, trax, epsc, 0.d0, mui, muzero_cor)

    ref(l,n)=fneray(l)*refx
    tra(l,n)=fneray(l)*trax

    ! In the aerosol layers
    call reftra  &
      ( 0.d0, piaero_o3, 0.d0, gaero_o3, 0.d0, tauao3(l), &
      refx, trax, epsc, 0.d0, mui, muzero_cor)

    ref(l,n)=ref(l,n) + (1.d0 - fneray(l))*refx
    tra(l,n)=tra(l,n) + (1.d0 - fneray(l))*trax

    trard(l,n) = (fneray(l)* dexp(-m*tauc(l)) + (1.d0 - fneray(l)))  &
                * dexp(-m*tauao3(l))

  end do
  ! Downward addition of layers
  do l = kmray, k1, -1
    ! Equations 33 of LH74
    drtt1=1.d0/(1.d0-reft(l+1,n)*ref(l,n))
    ! Note:
    ! R(top->l) = R(top->l+1)
    !           + T(top->l+1) R(l) T*(top->l+1) / (1 - R*(top->l+1)- R(l))
    ! Equations 34 of LH74
    ! Note R*(top->l) = R(top->l)
    !      T*(top->l+1) = T(top->l+1)
    reft(l,n)=reft(l+1,n) + trat(l+1,n)*ref(l,n)*trat(l+1,n) * drtt1

    trat(l,n)=trat(l+1,n)*tra(l,n) * drtt1

    ! Trad transmission for direct radiation
    trad(l,n) = trad(l+1,n)*trard(l,n)
  end do

  ! Upward addition of layers
  refb(k0,n)=ref(k0,n)

  do l = k1, kmray+1
    dtrb1=1.d0-ref(l,n)*refb(l-1,n)
    refb(l,n)=ref(l,n)+tra(l,n)*refb(l-1,n)*tra(l,n)/dtrb1

    ! Calculation of upward and downward fluxes and absorption
    dud1=1.d0-reft(l,n)*refb(l-1,n)
    ! Direct
    dowd(l,n)= trad(l,n)
    ! Diffuse
    dowf(l,n)=trat(l,n)/dud1-trad(l,n)
    ! Global == dowf(l,n) + dowd(l,n)
    dow(l,n) = trat(l,n)/dud1
    upwf(l,n)=refb(l-1,n)*dow(l,n)
    ! Absorption from top to level l
    atln(l,n)=1.d0-reft(k1,n)+upwf(l,n)-dow(l,n)
  enddo

  ! If there is a cloud
  if (itop .gt. k1) then
    do l = k1, itop - 1
      ! addition of ozone absorption for heating in the layers when adding method is used
      zq = zqq(l)
      zqp1 = zqq(l+1)
      zbas = zqq(itop)
      ! Note LH74 compute absn as:
      ! absn(l,n) = atln(l,n) - atln(l+1,n)
      ! but can be simplified as
      ! upwf(l,n)-dow(l,n) - upwf(l+1,n) + dow(l+1,n)
      absn(l,n) = dow(l,n) * (refb(l-1,n) - 1.d0) &
                - dow(l+1,n) * (refb(l,n) - 1.d0)

      x    = m*rayuoz(zbas) + mbar*(rayuoz(zq) - rayuoz(zbas))
      xp1  = m*rayuoz(zbas) + mbar*(rayuoz(zqp1) - rayuoz(zbas))
      xstar= m*rayuoz(zbas) + mbar*(rayuoz(zqq(k1)) - rayuoz(zbas)) &
                            + mbar*(rayuoz(zqq(k1)) - rayuoz(zq))
      xstarp1= m*rayuoz(zbas) + mbar*(rayuoz(zqq(k1)) - rayuoz(zbas)) &
                              + mbar*(rayuoz(zqq(k1)) - rayuoz(zqp1))

      !taking into account ozone absorption in the layers with clouds or aerosols
      fabso3(l) = muzero*fo*(raysoz(x) - raysoz(xp1) +                        &
                            (0.647d0- raysoz(m*rayuoz(zbas)))*absn(l,n)  &
                            +    rbar *(raysoz(xstarp1) - raysoz(xstar)))
      ! fluxes calculation taking into account ozone absorption
      ! Direct
      ddfso3(l)=muzero*fo*(0.647d0-rrbar-raysoz(x))*dowd(l,n)
      ! Diffuse:
      ! 2 contributions: one which transforms direct into diffuse (dowf)
      ! the other which is the diffuse coming from the top.
      dddfso3(l)=muzero*fo*(0.647d0-rrbar-raysoz(x)) * ( &
        dowf(l,n) &
       +dow(l,n)*albe*rrbar2s/(1.d0-rrbar2s*albe))
      ! Global: note this verifies:
      ! dfso3(l)= ddfso3(l)+dddfso3(l)
      dfso3(l)=muzero*fo*(0.647d0-rrbar-raysoz(x))*dow(l,n)/(1.d0-rrbar2s*albe)

      ufso3(l)=muzero*fo*(0.647d0-rrbar-raysoz(xstar))*upwf(l,n)

      dzx = drayuoz(zq)
      ! calculation of absorption coefficient ckup and ckdown useful for 3D simulation
      ckdown_suv_r(l)= dzxaoz(x,dzx)/(0.647d0-rrbar-raysoz(x))
      ! the optical depths for gases have to be changed in order to take
      ! into account the transformation of direct radiation in diffuse
      ! radiation under the top of the cloud
      ckdown_suv_f(l)= dzxaoz(x,dzx)/(0.647d0-rrbar-raysoz(x))
      ckup_suv_f(l)=dzxaoz(xstar,dzx)/(0.647d0-rrbar-raysoz(xstar))

    enddo

    ! Calculation of upward flux above cloud or aerosol layers taking into account
    ! the upward flux transmitted by cloud or aerosol layers
    ! if there is no cloud and no aerosol (itop=k1) this term have not to be add
    do l = itop, kmray+1
      zq = zqq(l)
      zbas = zqq(k1)
      xstar = m*rayuoz(zbas) + mbar*(rayuoz(zbas) - rayuoz(zq))
      ufso3(l)=muzero*fo*(0.647d0-rrbar-raysoz(xstar))*upwf(itop,n)
    enddo
  endif

  ! 6.4 Absorption by water vapor and liquid water (H20 band, SIR)

  ! In that case we have to solve multiple diffusion. This is achieved by means
  ! of the adding method following Lacis et Hansen, 1974
  ! calculation of reflexivity and transmissivity for each vertical layer

  do n = 1, 8

    ! Layer 44km -> 11km
    dqqv = kn(n)*(qqvtot - qqv(kmray))/10.d0

    ! Tod and ground BCs
    tau(kmray+1,n) = dqqv
    tra(kmray+1,n) = dexp(-m*tau(kmray+1,n))

    ! For direct radiation
    trard(kmray+1,n) = dexp(-m*tau(kmray+1,n))
    trad(kmray+1,n) = trard(kmray+1,n)

    ref(kmray+1,n) = 0.d0
    tra(k0,n) = 0.d0
    ref(k0,n) = albe

    trat(kmray+1,n) = tra(kmray+1,n)
    reft(kmray+1,n) = ref(kmray+1,n)

    ! For direct radiation
    trad(kmray+1,n) = trard(kmray+1,n)

    do l = k1, kmray
      gasym=gch2o(l)
      pioc=pic_h2o(l)
      ! Note /10 is for conversion purpose.
      dqqv = kn(n)*(qqv(l+1) - qqv(l))/10.d0 !Note qqv(1) = 0

      ! In the cloud  layers
      tau(l,n) = tauc(l) + dqqv + tauah2o(l)

      if (qlray(l).ge.epsc) then
        call reftra &
          (pioc, piaero_h2o, gasym, gaero_h2o, tauc(l) , tauah2o(l), &
          refx, trax, epsc, dqqv, mui, muzero_cor)

        ref(l,n)=fneray(l)*refx
        tra(l,n)=fneray(l)*trax + (1.d0 - fneray(l))*dexp(-mbarh2o * dqqv)

        if (iaer.eq.1) then
          call reftra &
            (0.d0, piaero_h2o, 0.d0, gaero_h2o, 0.d0 , tauah2o(l), &
            refx0, trax0, epsc, dqqv, mui, muzero_cor)

          ref(l,n) = fneray(l)*refx + (1.d0 - fneray(l))*refx0
          tra(l,n) = fneray(l)*trax + (1.d0 - fneray(l))*trax0
        endif

        ! trard transmissivity for direct radiation
        trard(l,n) = (fneray(l)*dexp(-m*tauc(l))+(1.d0-fneray(l))) &
                     *dexp(-m*(dqqv+tauah2o(l)))

      else

        ! in the clear sky layers
        ref(l,n) = 0.d0
        tra(l,n) = dexp(-mbarh2o * tau(l,n))

        trard(l,n)=dexp(-m*(dqqv+tauah2o(l)))

        if(l.ge.itop) tra(l,n) = dexp(-m*tau(l,n))
        if (iaer.eq.1) then
          call reftra  &
            (0.d0, piaero_h2o, 0.d0, gaero_h2o, 0.d0, tauah2o(l), &
            refx, trax, epsc,dqqv, mui, muzero_cor)

          ref(l,n)=fneba(l)*refx
          tra(l,n)=fneba(l)*trax+(1.d0-fneba(l))*dexp(-mbarh2o * dqqv)

        endif

      endif

    enddo

    ! downward addition of layers
    do l = kmray, k1, -1
      drtt1 = 1.d0 / (1.d0 - reft(l+1,n)*ref(l,n))
      ! Note R(l->top) = R*(l->top)
      ! R*(l->top) = R(l) + T*(l)*T(l) * R*(l+1>top)/(1 - R*(l+1->top).R(l))
      reft(l,n) = reft(l+1,n) &
        + trat(l+1,n)*ref(l,n)*trat(l+1,n)*drtt1
      ! T(l->top) = T(l+1->top).T(l) /(1 - R*(l+1->top).R(l))
      ! Note also it is equal to
      ! T*(l->top) = T*(l+1->top)*T(l) /(1 - R*(l+1->top).R(l))
      trat(l,n) = trat(l+1,n)*tra(l,n)*drtt1

      ! trad for direct radiation
      trad(l,n)=trad(l+1,n)*trard(l,n)

    enddo

    ! upward layer addition
    refb(k0,n) = ref(k0,n)

    ! Note LH74 equation (33)
    ! T*(l) = T(l)
    do l = k1, kmray
      dtrb1 = 1.d0/(1.d0 - refb(l-1,n)*ref(l,n))
      refb(l,n) = ref(l,n) + tra(l,n)*refb(l-1,n)*tra(l,n)*dtrb1
    enddo

    ! calculation of downward and upward fluxes
    do l = kmray+1, k1, -1
      ! downward fluxes for direct radiation
      dowd(l,n) = trad(l,n)

      dud1 = 1.d0 - reft(l,n)*refb(l-1,n)
      ! Diffuse
      dowf(l,n) = trat(l,n)/dud1-trad(l,n)

      ! Global
      dow(l,n) = trat(l,n)/dud1
      ! Diffuse up = global up
      upwf(l,n) = refb(l-1,n)*dow(l,n)

      ! calculation of absorption from 16km to level l
      ! Note can be factorized
      atln(l,n) =  pkn(n)*(1.d0 - reft(k1,n) + upwf(l,n) - dow(l,n))
    enddo

    ! absorption in individual layers
    do l = kmray, k1, -1
      ! Note LH74 compute absn as:
      ! absn(l,n) = atln(l,n) - atln(l+1,n)
      ! but can be simplified as
      ! upwf(l,n)-dow(l,n) - upwf(l+1,n) + dow(l+1,n)
      absn(l,n) = pkn(n) * ( dow(l,n) * (refb(l-1,n) - 1.d0)                   &
                           - dow(l+1,n) * (refb(l,n) - 1.d0))
    enddo
  enddo

  ! summation over frequencies and estimation of absorption integrated
  !  on the whole spectrum
  do l = kmray, k1, -1
    fabsh2o(l) = 0.d0
    do n = 1, 8
      fabsh2o(l) = fabsh2o(l)+absn(l,n)
    enddo
    fabsh2o(l) = fabsh2o(l)*fo*muzero
  enddo

  ! In the case with no clouds and aerosol in order to have exactly
  ! the same expressions in 1D and 3D for the heating rate
  if (itop.eq.k1) then
    do l = k1, kmray
      y = m*(qqvtot - qqv(l)) ! Note qqv(1) = 0
      yp1 = m*(qqvtot - qqv(l+1))
      ystar = m*qqvtot + mbarh2o * qqv(l)
      ystarp1 = m*qqvtot + mbarh2o * qqv(l+1)
      fabsh2o(l) = muzero*fo*(raysve(y) - raysve(yp1) + albe*(raysve(ystarp1) &
                                                             -raysve(ystar)))
    enddo
  endif
  ! 5.5 heating in the layers
  do i = k1, kmray
    deltaz = zqq(i+1) - zqq(i)

    cphum = cp0*(1.d0 + (cpvcpa - 1.d0)*qvray(i))
    rayst(i) = (fabsh2o(i) + fabso3(i))/deltaz/romray(i)/cphum

    rayst_h2o(i) = fabsh2o(i)/deltaz/romray(i)/cphum
    rayst_o3(i) =  fabso3(i) /deltaz/romray(i)/cphum
  enddo

  ! 5.6 calculation of solar fluxes
  !  for global radiation, fd for direct radiation for the water vapor band
  ! H20: SIR band

  do i = k1, kmray
    dfsh2o(i) = 0.d0
    ufsh2o(i) = 0.d0
    ddfsh2o(i) = 0.d0
    dddfsh2o(i) = 0.d0

    do n = 2,8
      ufsh2o(i) = ufsh2o(i) + pkn(n)*upwf(i,n)
      ddfsh2o(i) = ddfsh2o(i) + pkn(n)*dowd(i,n)
      dddfsh2o(i) = dddfsh2o(i) + pkn(n)*dowf(i,n)
    enddo

    ufsh2o(i) = fo*muzero*ufsh2o(i)
    ddfsh2o(i) = fo*muzero*ddfsh2o(i)
    dddfsh2o(i) = fo*muzero*dddfsh2o(i)
    ! Global
    dfsh2o(i)=ddfsh2o(i)+dddfsh2o(i)

    ! the optical depth for gases have to be changed to take into account
    ! the transformation of direct radiation in diffuse radiation under
    ! the top of the cloud
    ystar = m*(qqvtot-qqv(itop)) + mbarh2o * (qqv(itop)+qqv(i))
    if (i.ge.itop)then
      y = m*(qqvtot - qqv(i))
    else
      y = m*(qqvtot - qqv(itop)) + mbarh2o * (qqv(itop) - qqv(i))
    endif


    ! to test in the case with no clouds and aerosol in order to have exactly
    ! the same expressions in 1D and 3D for the fluxes
    if (itop.eq.k1) then
      ! Diffuse
      dddfsh2o(i) = 0.d0
      ! Global
      dfsh2o(i) = fo*muzero*(0.353d0-raysve(y))
      ufsh2o(i) = fo*muzero*(0.353d0-raysve(ystar))*albe
      ! Direct
      ddfsh2o(i) = fo*muzero*(0.353d0-raysve(y))
    endif

    ! (p_i/p_k) * (p_k/p0) = p_i/p0
    !We keep 1013.15 for standard atmosphere ref:LH74
    corp = preray(i) / 101315.d0
    dy = romray(i)*(qvray(i)*corp*sqrt(tkelvi/(temray(i) + tkelvi)))
    ! Calculation of absorption coefficient ckup and ckdown useful
    ! for 3D calculations
    ! Note: gas contribution
    ckdown_sir_r(i) = dzyama(y,dy)/(0.353d0-raysve(y))
    ckup_sir_f(i) = dzyama(ystar,dy)/(0.353d0-raysve(ystar))
    ckdown_sir_f(i) = dzyama(y,dy)/(0.353d0-raysve(y))
  enddo

  ! 6. Calculation of solar fluxes For the whole spectrum
  do k = k1, kmray
    ! Global (down and up) fluxes
    dfs(k) = dfsh2o(k) + dfso3(k)
    ufs(k) = ufsh2o(k) + ufso3(k)

    ! direct radiation and diffuse mod (sum of vapor water band and O3 band)
    ! Direct
    drfs(k) = ddfsh2o(k)+ddfso3(k)
    ! Diffuse
    ddfs(k) = dddfsh2o(k)+dddfso3(k)

    solu(k,ivertc) = ufs(k)
    sold(k,ivertc) = dfs(k)
  enddo

  ! Note: Multiplication by transmission function for minor gases
  ! Tmg is now taking into account by fo=fo*Tmg

  ! solar heating of the ground surface by the downward global flux
  fos=dfs(k1)*(1.d0-albe)


  ! Calculation of absorption coefficient ckup and ckdown useful
  ! for 3D simulation
  ! For water vapor without clouds and aerosols downward is only direct,
  ! upward is only diffuse

  ! SIR band
  do k = k1, kmray
    tauapc=tauah2o(k)+tauc(k)
    deltaz = zqq(k+1) - zqq(k)
    if (tauapc.lt.epsc)then
      ckapcd=0.d0
      ckapcf=0.d0
      ck_aero_h2of=0.d0
      ck_aero_h2od=0.d0
      w0_sir(k) = 0.d0
      g_apc_sir(k) = 0.d0
    else
      picapc=(pic_h2o(k)*tauc(k)+piaero_h2o*tauah2o(k))/tauapc

      ! if we take into account asymmetry factor for forward diffuse radiation
      ! Note apc means aerosols+clouds
      gapc=(pic_h2o(k)*tauc(k)*gch2o(k)+piaero_h2o*tauah2o(k)*gaero_h2o)/(tauapc*picapc)

      ! Save values without Joseph correction for 3D
      g_apc_sir(k) = gapc
      w0_sir(k) = picapc
      ! absorption and forward diffusion
      ckapcf = tauapc/deltaz

      ! direct do not take Joseph correction into account
      ckapcd=tauapc/(deltaz*muzero_cor)
      ck_aero_h2of = tauah2o(k)/deltaz
      ck_aero_h2od=tauah2o(k)/(deltaz*muzero_cor)
    endif

    ckup_sir_f(k) = ckup_sir_f(k) + ck_aero_h2of*(1.d0-fneray(k)) &
                  + ckapcf * fneray(k)
    ckdown_sir_r(k) = ckdown_sir_r(k) + ck_aero_h2od*(1.d0-fneray(k))  &
                    + ckapcd*fneray(k)
    ckdown_sir_f(k) = ckdown_sir_f(k) + ck_aero_h2of*(1.d0-fneray(k))  &
                    + ckapcf * fneray(k)
  enddo

  ! SUV band
  do k = k1, kmray
    tauapc=tauao3(k)+tauc(k)
    deltaz = zqq(k+1) - zqq(k)
    if(tauapc.lt.epsc)then
      ckapcd=0.d0
      ckapcf=0.d0
      ck_aero_o3f=0.d0
      ck_aero_o3d=0.d0
      w0_suv(k) = 0.d0
      g_apc_suv(k) = 0.d0
    else
      picapc=(pic_o3(k)*tauc(k) + piaero_o3*tauao3(k))/tauapc
      ! if we take into account asymmetry factor for forward diffuse radiation
      gapc=(pic_o3(k)*tauc(k)*gco3(k)+piaero_o3*tauao3(k)*gaero_o3)/(tauapc*picapc)
      ! direct do not take Joseph correction into account
      ckapcd=tauapc/(deltaz*muzero_cor)
      w0_suv(k) = picapc
      g_apc_suv(k) = gapc
      ! absorption and forward diffusion
      ckapcf=tauapc/deltaz
      ck_aero_o3f=tauao3(k)/deltaz
      ck_aero_o3d=tauao3(k)/(deltaz*muzero_cor)
    endif

    ckup_suv_f(k) = ckup_suv_f(k) + ck_aero_o3f*(1.d0-fneray(k))&
                  + ckapcf * fneray(k)
    ckdown_suv_r(k) = ckdown_suv_r(k) + ck_aero_o3d*(1.d0-fneray(k)) &
                    + ckapcd * fneray(k)

    ckdown_suv_f(k) = ckdown_suv_f(k) + ck_aero_o3f*(1.d0-fneray(k))&
                    + ckapcf * fneray(k)
  enddo

  !In addition a source term has to be added in 3D for diffuse radiation

  ! if mui=1/sqrt(3) quadrature method as LH74 the source term is added for
  ! both downward and upward radiation
  ! where  g1= 31/2(1-wo) , g3=(1-31/2µo)/2, g4=(1+31/2µo)/2 and t the total
  ! optical depth (gas + aerosol +
  ! cloud) t = tg + ta + tc.
  ! if mui=muzero_cor delta method the source term is added only in the downward
  !diffuse radiation
  ! In that condition the two-stream equations can be written:
  ! where  g1=(1-wo)/µo, g4=1 and t the total optical depth (gas + aerosol + cloud) t = tg + ta + tc.

! if muzero < 0, it is night
else

  muzero = 0.d0
  do k = k1, kmray
    rayst(k) = 0.d0

    rayst_h2o(k) = 0.d0
    rayst_o3(k) = 0.d0

    solu(k,ivertc) = 0.d0
    sold(k,ivertc) = 0.d0

    ckup_sir_f(k) = 0.d0
    ckdown_sir_r(k) = 0.d0
    ckdown_sir_f(k) = 0.d0
    ckup_suv_f(k)=0.d0
    ckdown_suv_r(k)=0.d0
    ckdown_suv_f(k)=0.d0
    w0_sir(k) = 0.d0
    w0_suv(k) = 0.d0
    g_apc_sir(k) = 0.d0
    g_apc_suv(k) = 0.d0

  enddo
  fos = 0.d0
endif

! Compute Boundary conditions for the 3D (Director diFfuse) Solar radiance
! at the top of the CFD domain
! and the absorption coefficients
call field_get_id_try("spectral_rad_incident_flux", f_id)

is_active = cs_rad_time_is_active()

if (f_id.ge.0.and.is_active) then

  call cs_f_boundary_conditions_get_pointers(c_itypfb)
  call c_f_pointer(c_itypfb, itypfb, [nfabor])

  call field_get_val_v(f_id, bpro_rad_inc)

  call field_get_val_v_by_name("rad_absorption_coeff_up", cpro_ck_up)
  call field_get_val_v_by_name("rad_absorption_coeff_down", cpro_ck_down)
  call field_get_val_v_by_name("simple_diffusion_albedo", cpro_w0)
  call field_get_val_v_by_name("asymmetry_factor", cpro_gapc)

  c_id = 0
  ! Direct Solar (denoted by _r) (for Solar IR band absorbed by H20)
  if (iand(rad_atmo_model, 1).eq.1) then

    c_id = c_id + 1

    ! Store the incident radiation of the 1D model
    do ifac = 1, nfabor
      if(itypfb(ifac).eq.iparug.or.itypfb(ifac).eq.iparoi) then
        bpro_rad_inc(c_id, ifac) = 0.d0
      else

        ! Interpolate at zent
        zent = cdgfbo(3, ifac)

        call intprz &
          (kmray, zqq,                                               &
          ddfsh2o, zent, iz1, iz2, var )

        ! TODO do not multiply and divide by cos(zenital) = muzero
        if (muzero.gt.epzero) then
          bpro_rad_inc(c_id, ifac) = var / muzero_cor
        else
          bpro_rad_inc(c_id, ifac) = 0.d0
        endif
      endif
    enddo

    ! Store the (downward) absorption coefficient of the 1D model
    do iel = 1, ncel

      ! Interpolate at zent
      zent = xyzcen(3, iel)

      call intprz &
        (kmray, zray,                                               &
        ckdown_sir_r, zent, iz1, iz2, var )

      cpro_ck_down(c_id, iel) = var
    enddo

  endif

  ! Direct Solar (denoted by _r) (for visible UV (SUV) band absorbed by O3)
  if (iand(rad_atmo_model, 2).eq.2) then

    c_id = c_id + 1

    ! Store the incident radiation of the 1D model
    do ifac = 1, nfabor
      if(itypfb(ifac).eq.iparug.or.itypfb(ifac).eq.iparoi) then
        bpro_rad_inc(c_id, ifac) = 0.d0
      else

        ! Interpolate at zent
        zent = cdgfbo(3, ifac)

        call intprz &
          (kmray, zqq,                                               &
          ddfso3, zent, iz1, iz2, var )

        ! TODO do not multiply and divide by cos(zenital) = muzero
        if (muzero.gt.epzero) then
          bpro_rad_inc(c_id, ifac) = var / muzero_cor
        else
          bpro_rad_inc(c_id, ifac) = 0.d0
        endif
      endif
    enddo

    ! Store the (downward) absortion coefficient of the 1D model
    do iel = 1, ncel

      ! Interpolate at zent
      zent = xyzcen(3, iel)

      call intprz &
        (kmray, zray,                                               &
        ckdown_suv_r, zent, iz1, iz2, var )

      cpro_ck_down(c_id, iel) = var
    enddo

    ! Direct Solar (denoted by _r): if O3 band not activated it is added
    ! to total solar band (ie the one of H20)
  else if (iand(rad_atmo_model, 1).eq.1) then

    ! Store the incident radiation of the 1D model
    do ifac = 1, nfabor
      if(itypfb(ifac).eq.iparug.or.itypfb(ifac).eq.iparoi) then
        bpro_rad_inc(c_id, ifac) = 0.d0
      else

        ! Interpolate at zent
        zent = cdgfbo(3, ifac)

        call intprz &
          (kmray, zqq,                                               &
          ddfso3, zent, iz1, iz2, var )

        ! TODO do not multiply and divide by cos(zenital) = muzero
        if (muzero.gt.epzero) then
          bpro_rad_inc(c_id, ifac) = bpro_rad_inc(c_id, ifac) + var / muzero_cor
        else
          bpro_rad_inc(c_id, ifac) = 0.d0
        endif
      endif
    enddo

    ! Store the (downward) absortion coefficient of the 1D model
    do iel = 1, ncel

      ! Interpolate at zent
      zent = xyzcen(3, iel)

      call intprz &
        (kmray, zray,                                               &
        ckdown_suv_r, zent, iz1, iz2, var )

      cpro_ck_down(c_id, iel) = cpro_ck_down(c_id, iel) + var
    enddo


  endif

  ! Diffuse solar radiation incident up and down (SIR band)
  if (iand(rad_atmo_model, 4).eq.4) then

    c_id = c_id + 1
    do ifac = 1, nfabor

      if(itypfb(ifac).eq.iparug.or.itypfb(ifac).eq.iparoi) then
        bpro_rad_inc(c_id, ifac) = 0.d0
      else
        ! Interpolate at zent
        zent = cdgfbo(3, ifac)

        call intprz &
          (kmray, zqq,                                               &
          dddfsh2o, zent, iz1, iz2, var )

        bpro_rad_inc(c_id, ifac) = var
      endif

    enddo

    ! Store the (downward and upward) absorption coefficient of the 1D model
    do iel = 1, ncel

      ! Interpolate at zent
      zent = xyzcen(3, iel)

      ! Note: it is zero
      call intprz &
        (kmray, zray,                                               &
        ckdown_sir_f, zent, iz1, iz2, var )

      cpro_ck_down(c_id, iel) = var

      call intprz &
        (kmray, zray,                                               &
        ckup_sir_f, zent, iz1, iz2, var )

      cpro_ck_up(c_id, iel) = var

      ! Simple diffusion albedo w0

      call intprz &
        (kmray, zray,                                               &
        w0_sir, zent, iz1, iz2, var )

      cpro_w0(c_id, iel) = var

      ! Asymmetry factor

      call intprz &
        (kmray, zray,                                               &
        g_apc_sir, zent, iz1, iz2, var )

      cpro_gapc(c_id, iel) = var

    enddo

  endif

  ! Diffuse solar radiation incident up and down (SUV - O3 band)
  if (iand(rad_atmo_model, 8).eq.8) then

    c_id = c_id + 1
    do ifac = 1, nfabor

      if(itypfb(ifac).eq.iparug.or.itypfb(ifac).eq.iparoi) then
        bpro_rad_inc(c_id, ifac) = 0.d0
      else
        ! Interpolate at zent
        zent = cdgfbo(3, ifac)

        call intprz &
          (kmray, zqq,                                               &
          dddfso3, zent, iz1, iz2, var )

        bpro_rad_inc(c_id, ifac) = var
      endif

    enddo

    ! Store the (downward and upward) absortion coefficient of the 1D model
    do iel = 1, ncel

      ! Interpolate at zent
      zent = xyzcen(3, iel)

      call intprz &
        (kmray, zray,                                               &
        ckdown_suv_f, zent, iz1, iz2, var )

      cpro_ck_down(c_id, iel) = var

      call intprz &
        (kmray, zray,                                               &
        ckup_suv_f, zent, iz1, iz2, var )

      cpro_ck_up(c_id, iel) = var

      ! Simple diffusion albedo w0

      call intprz &
        (kmray, zray,                                               &
        w0_suv, zent, iz1, iz2, var )

      cpro_w0(c_id, iel) = var

      ! Asymmetry factor

      call intprz &
        (kmray, zray,                                               &
        g_apc_suv, zent, iz1, iz2, var )

      cpro_gapc(c_id, iel) = var
    enddo

    ! If Diffuse solar O3 band in 3D is not activated -> add in the total
    ! diffuse solar band
  elseif (iand(rad_atmo_model, 4).eq.4) then

    do ifac = 1, nfabor

      if(itypfb(ifac).eq.iparug.or.itypfb(ifac).eq.iparoi) then
        bpro_rad_inc(c_id, ifac) = 0.d0
      else
        ! Interpolate at zent
        zent = cdgfbo(3, ifac)

        call intprz &
          (kmray, zqq,                                               &
          dddfso3, zent, iz1, iz2, var )

        bpro_rad_inc(c_id, ifac) = bpro_rad_inc(c_id, ifac) + var
      endif

    enddo

    ! Store the (downward and upward) absortion coefficient of the 1D model
    do iel = 1, ncel

      ! Interpolate at zent
      zent = xyzcen(3, iel)

      call intprz &
        (kmray, zray,                                               &
        ckdown_suv_f, zent, iz1, iz2, var )

      cpro_ck_down(c_id, iel) = cpro_ck_down(c_id, iel) + var

      call intprz &
        (kmray, zray,                                               &
        ckup_suv_f, zent, iz1, iz2, var )

      cpro_ck_up(c_id, iel) = cpro_ck_up(c_id, iel) + var

      ! Simple diffusion albedo w0

      call intprz &
        (kmray, zray,                                               &
        w0_suv, zent, iz1, iz2, var )

      cpro_w0(c_id, iel) = cpro_w0(c_id, iel) + var

      ! Asymmetry factor

      call intprz &
        (kmray, zray,                                               &
        g_apc_suv, zent, iz1, iz2, var )

      cpro_gapc(c_id, iel) = cpro_gapc(c_id, iel) + var
    enddo

  endif

endif

deallocate(fabsh2o,fabso3,tauc)
deallocate(tau,pic,ref)
deallocate(reft)
deallocate(refb,upwf)
deallocate(fabso3c,tra)
deallocate(dow)
deallocate(dowf,atln,absn)
deallocate(fnebmax,fneba)

deallocate(ufso3c)
deallocate(dowd, trad, trard)

deallocate(dfsh2o, ufsh2o, dfso3, ufso3, dfs, ufs)
deallocate(dffsh2o, dffso3)
deallocate(ddfsh2o, ddfso3)
deallocate(dddfsh2o, dddfso3)
deallocate(drfs, ddfs)
deallocate(ckup_sir_f, ckdown_sir_r, ckdown_sir_f)
deallocate(ckup_suv_f, ckdown_suv_r, ckdown_suv_f)
deallocate(g_apc_suv, g_apc_sir, w0_suv, w0_sir)

return

!===============================================================================

contains

  !-----------------------------------------------------------------------------

  !> \brief Computes ozone amount above a given altitude
  ! See LH74 equation (17)

  !> \param[in]   zh          altitude

  function rayuoz(zh)

    !===========================================================================

    implicit none

    ! Arguments

    double precision, intent(in) :: zh  ! absolute altitude

    ! Local

    double precision ::  rayuoz
    double precision ::  a, b, c

    !===========================================================================

    a = 0.4d0
    b = 20000.d0
    c = 5000.d0

    rayuoz = a*(1.d0 + dexp(-b/c))/(1.d0 + dexp((zh-b)/c))

  end function rayuoz

  !-----------------------------------------------------------------------------

  !> \brief Absorption function of the solar radiation by water vapor

  !> \param[in]       y       optical depth for water vapor

  function raysve(y)

    !===========================================================================

    implicit none

    ! Arguments

    double precision, intent(in) :: y        ! specific humidity

    ! Local

    double precision :: raysve

    !===========================================================================

    raysve = 0.29d0*y/((1.d0 + 14.15d0*y)**0.635d0 + 0.5925d0*y)

  end function raysve

  !-----------------------------------------------------------------------------

  !> \brief Absorption derivative-function of the solar radiation by water vapor

  !> \param[in]       y       optical depth for water vapor
  !> \param[in]       dy      TODO?

  function dzyama(y, dy)

    implicit none

    ! Arguments

    double precision, intent(in):: y, dy

    double precision:: num,denum
    double precision:: dzyama

    num = 14.15d0*0.635d0*dy*(1.0d0 + 14.15d0*y)**(0.635d0-1.0d0) + 0.5925d0*dy
    denum = (1.0d0 + 14.15d0*y)**0.635d0 + 0.5925d0*y
    dzyama= 0.29d0*dy/denum - 0.29d0*y*num/(denum**2.0d0)

  end function dzyama

  !-----------------------------------------------------------------------------

  !> \brief Absorption function of the solar radiation by ozone

  !> \param[in]       x       optical depth for ozone

  function raysoz(x)

    !===========================================================================

    implicit none

    ! Arguments

    double precision, intent(in) :: x

    ! Local

    double precision :: raysoz, a, b, c, d, e, f, g, h
    double precision :: ao3vis, ao3uv

    !===========================================================================

    a = 0.02118d0
    b = 0.042d0
    c = 0.000323d0
    d = 1.082d0
    e = 138.6d0
    f = 0.805d0
    g = 0.0658d0
    h = 103.6d0
    ao3vis = a*x/(1.d0 + (b + c*x)*x)
    ao3uv =  d*x/(1.d0 + e*x)**f                     &
           + g*x/(1.d0 + (h*x)**3)
    raysoz = ao3vis + ao3uv

  end function raysoz

  !-----------------------------------------------------------------------------

  !> \brief Absorption derivative-function of the solar radiation by ozone

  !> \param[in]       x
  !> \param[in]       dx

  function dzxaoz(x, dx)

    implicit none

    ! Arguments

    double precision, intent(in):: x, dx

    double precision:: num1,denum1,num2,denum2,num3,denum3
    double precision:: dzxaoz

    num1=0.02118d0*(1.d0-0.000323d0*x**2.d0)
    denum1=(1.d0+0.042d0*x+0.000323d0*x**2.d0)**2.d0
    num2=0.0658d0*(1.d0-2.d0*(103.6d0*x)**3.d0)
    denum2=(1.d0+(103.6d0*x)**3.d0)**2.d0
    num3=1.082d0*((1.d0+138.6d0*x)**(0.805d0)-x*(0.805d0*138.6d0 &
         *(1.d0+138.6d0*x)**(0.805d0-1.d0)))
    denum3=(1.d0+138.6d0*x)**(2.d0*0.805d0)

    dzxaoz= dx*(num1/denum1+num2/denum2+num3/denum3)

  end function dzxaoz

function drayuoz(zh)

  implicit none
  double precision, intent(in) :: zh  ! absolute altitude
  double precision ::  drayuoz
  double precision ::  a, b, c

  a = 0.4d0
  b = 20000.d0
  c = 5000.d0

  drayuoz = a / c *(1.d0 + dexp(-b/c))              &
             * dexp((zh-b)/c)                        &
             / ((1.d0 + dexp((zh-b)/c))**(2.0d0))

end function drayuoz

end subroutine rayso
