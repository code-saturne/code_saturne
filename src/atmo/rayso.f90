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
!_______________________________________________________________________________

subroutine reftra  &
    (pioc, piaero, gasym, gaero, tauc, taua, &
    ref, tra, epsc,dqqv, mui)

  !===========================================================================

  implicit none

  ! Arguments

  double precision, intent(in) :: pioc, piaero,  gasym, gaero
  double precision, intent(in) :: tauc, taua,  dqqv, epsc
  double precision, intent(inout) :: ref, tra

  ! Local
  double precision ::gas, fas, kt, gama1, gama2, tln
  double precision :: drt, extlnp, extlnm
  double precision :: pic, tau, mui

  !===========================================================================

  tau = tauc +taua + dqqv
  ! For 0 optical depth
  if (tau .lt. epsc) then
    ref = 0.d0
    tra = 1.d0
  else

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

        gas=(pioc*tauc*gasym+piaero*taua*gaero)&
          /(pic*tau)
        ! Only forward diffusion
        if(gas.ge.(1.d0-epsc)) then
          gama1=(1.d0-pic)/mui
          ref=0.d0
          tra=dexp(-gama1*tau)
        else
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
!> \param[in]   qqv         optical depth for water vapor (0,z)
!> \param[in]   qqqv        idem for intermediate levels
!> \param[in]   qqvinf      idem qqv but for altitude above 11000m
!> \param[in]   zqq         vertical levels
!> \param[in]   zray        altitude (physical mesh)
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

use paramx
use numvar
use optcal
use cstphy
use entsor
use parall
use period
use ppppar
use ppthch
use ppincl
use cs_c_bindings
use mesh
use field
use atincl, only: kmx, nbmett, sigc, squant, xlat, xlon, sold, &
                  solu, piaero_o3,piaero_h2o, &
                  black_carbon_frac,zaero, gaero_o3, gaero_h2o, &
                  aod_h2o_tot, aod_o3_tot
use ctincl, only: cp_a, cp_v
use cstnum, only: epzero, pi
use pointe, only: itypfb
use radiat

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

integer i,k,n,l,k1p1,iaer,iaero_top,iaero_topr
integer itop,ibase,itopp1,itopp2,ibasem1
integer          ifac, iz1, iz2, f_id, c_id, iel
double precision za, muzero, muzero_cor,fo,m,mbar,rabar,rabar2,rbar
double precision rabarc,rbarc, refx, trax, refx0, trax0
double precision qqvtot,y,ystar
double precision yp1,ystarp1
double precision zqm1,zq,xm1,x,xstar,xstarm1
double precision rrbar,rrbar2s
double precision tauctot,wh2ol,rm,req,deltaz
double precision pioc,zbas,dud1
double precision gasym,drt,tln,drtt1,dtrb1
double precision kn(8),pkn(8),dqqv
double precision cphum,qureel
double precision taua(kmx+1),tauatot
double precision rabara,rbara
double precision tauca(kmx+1,8)
double precision gama1,gama2,kt,gas,fas
double precision omega, var, zent
double precision cpvcpa
double precision dzx, dy
! For postprocessing
double precision soil_direct_flux , soil_global_flux
double precision soil_direct_flux_h2o,  soil_global_flux_h2o
double precision soil_direct_flux_o3,  soil_global_flux_o3
logical          is_active

double precision, allocatable:: fabsh2o(:),fabso3(:),tauc(:)
double precision, allocatable:: tau(:,:),pic(:,:),ref(:,:)
double precision, allocatable:: reft(:,:),trat(:,:),refts(:,:)
double precision, allocatable:: refs(:,:),tras(:,:),trats(:,:)
double precision, allocatable:: refb(:,:),trab(:,:),upw(:,:)
double precision, allocatable:: refbs(:,:),fabso3c(:,:),tra(:,:)
double precision, allocatable:: dow(:,:),atln(:,:),absn(:,:)
double precision, allocatable :: ckup_sir_f(:), ckdown_sir_r(:), ckdown_sir_f(:)
double precision, allocatable :: ckup_suv_f(:), ckdown_suv_r(:), ckdown_suv_f(:)
double precision, allocatable :: w0_suv(:), w0_sir(:), g_apc_suv(:), g_apc_sir(:)
double precision, allocatable:: fnebmax(:),fneba(:)

double precision, allocatable, dimension(:,:) :: dowd, trad, trard, ddfso3c
double precision, allocatable, dimension(:) :: dffsh2o, dffso3
double precision, allocatable, dimension(:) :: ddfsh2o, ddfso3
double precision, allocatable, dimension(:) :: drfs, dffs, ddfs, dddsh2o, dddso3

double precision, allocatable, dimension(:) ::  dfsh2o, ufsh2o
double precision, allocatable, dimension(:) ::  dfso3, ufso3
double precision, allocatable, dimension(:,:) ::  dfso3c, ufso3c
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
double precision corp, rov
double precision mui,tauapc,ckapcd,ckapcf,picapc,gapc
double precision ck_aero_h2of,ck_aero_h2od,ck_aero_o3f,ck_aero_o3d
! data for pkn and kn distribution
data kn/4.d-5,0.002,0.035,0.377,1.95,9.40,44.6,190./
data pkn/0.6470,0.0698,0.1443,0.0584,0.0335,0.0225,0.0158,0.0087/

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
allocate(tau(kmx+1,8),pic(kmx+1,8),ref(kmx+1,8))
allocate(reft(kmx+1,8),trat(kmx+1,8),refts(kmx+1,8))
allocate(refs(kmx+1,8),tras(kmx+1,8),trats(kmx+1,8))
allocate(refb(kmx+1,8),trab(kmx+1,8),upw(kmx+1,8))
allocate(refbs(kmx+1,8),fabso3c(kmx+1,2),tra(kmx+1,8))
allocate(dow(kmx+1,8),atln(kmx+1,8),absn(kmx+1,8))
allocate(fnebmax(kmx+1),fneba(kmx+1))

allocate(dfso3c(kmx+1,2), ufso3c(kmx+1,2), ddfso3c(kmx+1,2))
allocate(dowd(kmx+1,8), trad(kmx+1,8), trard(kmx+1,8))

allocate(dfsh2o(kmx+1), ufsh2o(kmx+1))
allocate(dfso3(kmx+1), ufso3(kmx+1))
allocate(dfs(kmx+1), ufs(kmx+1))
allocate(ckup_sir_f(kmx), ckdown_sir_r(kmx), ckdown_sir_f(kmx))
allocate(ckup_suv_f(kmx), ckdown_suv_r(kmx), ckdown_suv_f(kmx))
allocate(w0_sir(kmx), w0_suv(kmx), g_apc_sir(kmx), g_apc_suv(kmx))

allocate(dffsh2o(kmx+1), dffso3(kmx+1))
allocate(ddfsh2o(kmx+1), ddfso3(kmx+1))
allocate(drfs(kmx+1), dffs(kmx+1), ddfs(kmx+1), dddsh2o(kmx+1), dddso3(kmx+1))

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
  ! TODO usefull ?
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
  dddsh2o(k) = 0.d0
  dddso3(k) = 0.d0
  ddfs(k) = 0.d0
  dffs(k) = 0.d0
  drfs(k) = 0.d0
  dffsh2o(k) =0.d0
  dffso3(k) = 0.d0
  ddfsh2o(k) = 0.d0
  ddfso3(k) = 0.d0
  dfsh2o(k) = 0.d0
  ufsh2o(k) = 0.d0
  dfso3(k) = 0.d0
  ufso3(k) = 0.d0

  if(iaer.ge.1) then
    fneba(k) = 1.d0
  endif
  do l = 1, 2
    fabso3c(k,l) = 0.d0
  enddo
  do l = 1, 8
    tau(k,l) = 0.d0
    tauca(k,l)=0.d0
    pic(k,l) = 0.d0
    atln(k,l) = 0.d0
    absn(k,l) = 0.d0
    ref(k,l) = 0.d0
    reft(k,l) = 0.d0
    refts(k,l) = 0.d0
    refs(k,l) = 0.d0
    refb(k,l) = 0.d0
    refbs(k,l) = 0.d0
    trat(k,l) = 0.d0
    tras(k,l) = 0.d0
    trab(k,l) = 0.d0
    trats(k,l) = 0.d0
    tra(k,l) = 0.d0
    upw(k,l) = 0.d0
    dow(k,l) = 0.d0
  enddo
enddo

do k = 1, kmx+1
  do l = 1, 8
    trad(k,l) = 0.d0
    dowd(k,l) = 0.d0
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
k1p1 = k1+1

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
  ! TODO test 5/3 but not coherent with LH74...

  ! as we suppressed m and mbar in dzx and dzy we can eventually keep 1.9 for O3 to tested
  mbar = 1.9d0
  ! for LH74 quadrature method in two-stream
  mui=1.d0/dsqrt(3.d0)

  !  3 -  albedos for O3 and Rayleigh diffusion

  rabar = 0.219d0/(1.d0 + 0.816d0*muzero)
  rabar2 = 0.144d0
  rbar = rabar + (1.d0 - rabar)*(1.d0 - rabar2)*albe/(1.d0 - rabar2*albe)
  rrbar = 0.28d0/(1.d0 + 6.43d0*muzero)
  rrbar2s = 0.0685d0

  !  4 - addition of one level for solar radiation

  zqq(kmray+1) = 16000.d0
  qqvtot = qqvinf + qqqv(kmray)
  qqv(kmray+1) = qqvtot - qqvinf/4.d0

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

  itop = 0

  do i = kmray, k1p1, -1
    if(qlray(i).gt.epsc) then
      if (itop.eq.0) then
        itop = i
      endif
      ibase = i
    endif
  enddo

  ! if itop = 0, there is no cloud but, nevertheless, it is possible to execute
  ! the adding method for the water vapor (SIR band) only

  if (itop.eq.0) then
    itop = k1
    ibase = k1
  endif
  itopp1 = itop +1
  itopp2 = itop +2
  ibasem1 = ibase -1

  ! 5.2 calculation for optical parameters of clouds and aerosols
  ! (single scattering albedo, optical depth, radius, asymmetry factor)

  fnebmax(kmray+1) = 0.d0
  tauctot = 0.d0
  tauatot = 0.d0
  do i = kmray, k1p1, -1
    if((i.gt.ibasem1).and.(i.lt.itopp1)) then
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

      deltaz = zqq(i+1) - zqq(i)
      ! Cloud optical thickness
      if (i.eq.k1p1) deltaz = zqq(i+1) - zray(k1)
      ! req has to be in µm
      tauc(i) = 1.5d0 * wh2ol * deltaz / req
      tauctot = tauctot + tauc(i)
    else
      tauc(i) = 0.d0
      req = 0.d0
    endif

    ! Calculation of aerosol optical depth AOD
    fneba(i) = 0.d0
    if((iaer.eq.1).and.(zqq(i).le.zaero)) then
      iaero_top=max(i,iaero_top)
      fneba(i) = 1.d0
      deltaz = zqq(i+1) - zqq(i)
      if(i.eq.k1p1) deltaz=zqq(i+1) - zray(k1)
      ! Distribution of AOD on the vertical
      ! Note, we used a formula based on concentration before v6.2
      tauao3(i) = aod_o3_tot*deltaz/zqq(iaero_top+1)
      tauah2o(i) = aod_h2o_tot*deltaz/zqq(iaero_top+1)
    endif
    ! Estimation of the law of the cloud fraction

    ! Calculation of SSA and Asymmetry factor for clouds using- Nielsen 2014
    ! Note only the first two bands are taken, the third only is
    ! approximately 0
    ! Usefull for gas assymmetry, for albedo, either Chuang if BC, or Nielsen
    pioco3_1 = 1.d0 - 33.d-9*req
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

    ! Gas assymmmetry: Nielsen
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

  fnebmax(k1) = fnebmax(k1p1)

  tauc(kmray+1) = 0.d0

  ! 5.3 O3 absorption in presence of clouds

  ! Calculation of the different albedos for O3 (LH 74)

  ! Asymmetry factor and SSA for liquid water
  gasym=gco3(itop)

  pioc=pic_o3(itop)
  !Calculation for cloudy layers
  call reftra  &
    (pioc, 0.d0, gasym, 0.d0, tauctot, 0.d0, &
    refx, trax, epsc, 0.d0, mui)

  rabarc=refx
  tabarc=trax
  rbarc = rabarc+tabarc*tabarc*albe/(1.d0 - rabarc*albe)
  !Calculation for aerosol layers
  call reftra  &
    (0.d0, piaero_o3, 0.d0, gaero_o3, 0.d0, aod_o3_tot, &
    refx, trax, epsc,0.d0, mui)
  rabara=refx
  tabara=trax

  rbara = rabara+tabara*tabara*albe/(1.d0 - rabara*albe)

  ! in the case there is an aerosol layer above the cloud layer

  if((iaer.eq.1).and.(iaero_top.gt.itop)) then
    itop=iaero_top
    itopp1=itop+1
    rbar=rbara
    rrbar2s=rabara
  endif

  ! Calculation above the top of the cloud or the aerosol layer

  do i = itop, kmray   !calculation have to start at the first level
    ! ( itop=1  in the case without aerosols and clouds)
    zqm1 = zqq(i)

    zq = zqq(i+1)
    ! For the ground, we have to get fabso3c(1,1)=0. and ddfso3c(1,1)= incomming
    ! flux. It is the cas with zq=zqq(k1)(=zray(k1))
    if (i.eq.k1) zq = zray(k1)
    xm1 = m*rayuoz(zqm1)
    x = m*rayuoz(zq)

    ! Calculation of heat and radiation fluxes during cloudy sky
    zbas = zray(itop)
    xstarm1 = m*rayuoz(zbas) + mbar*(rayuoz(zbas) - rayuoz(zqm1))
    xstar = m*rayuoz(zbas) + mbar*(rayuoz(zbas) - rayuoz(zq))

    ! Heat
    fabso3c(i,1) = muzero*fo*(raysoz(xm1) - raysoz(x))  &
      + rbarc*(raysoz(xstar) - raysoz(xstarm1))

    ! Direct downward radiation
    ddfso3c(i,1) = muzero*fo*(0.647d0-rrbar-raysoz(x))
    ! Downward radiation
    dfso3c(i,1) = muzero*fo*(0.647d0-rrbar-raysoz(x))
    ! Upward (diffuse) radiation
    ufso3c(i,1) = muzero*fo*(0.647d0-rrbar-raysoz(xstar))*rabarc

    ! Calculation ofheat and radiation fluxes during  Clear sky
    zbas = zray(k1)
    xstarm1 = m*rayuoz(zbas) + mbar*(rayuoz(zbas) - rayuoz(zqm1))
    xstar = m*rayuoz(zbas) + mbar*(rayuoz(zbas) - rayuoz(zq))
    !heat
    fabso3c(i,2) = muzero*fo*(raysoz(xm1) - raysoz(x)) &
      +rbar*(raysoz(xstar) - raysoz(xstarm1))

    dfso3c(i,2) = muzero*fo*(0.647d0-rrbar-raysoz(x)) &
      /(1.d0-rrbar2s*albe)

    ddfso3c(i,2) = muzero*fo*(0.647-rrbar-raysoz(x))
    ddfso3c(i,2) = min(ddfso3c(i,2),dfso3c(i,2))
    ufso3c(i,2) = muzero*fo*(0.647d0-rrbar-raysoz(xstar))*albe         &
      / (1.d0-rrbar2s*albe)
    ! Summation depending on cloud fraction
    fabso3(i) = fnebmax(k1p1)*fabso3c(i,1) + (1.d0-fnebmax(k1p1))*fabso3c(i,2)

    dfso3(i) = fnebmax(k1p1)*dfso3c(i,1)+ (1.d0-fnebmax(k1p1))*dfso3c(i,2)
    ddfso3(i) = fnebmax(k1p1)*ddfso3c(i,1)+ (1.d0-fnebmax(k1p1))*ddfso3c(i,2)
    ufso3(i) = fnebmax(k1p1)*ufso3c(i,1)+ (1.d0-fnebmax(k1p1))*ufso3c(i,2)
    ! Calculation of absorption coefficient ckup and ckdown
    ! useful for 3D simulation
    dzx = drayuoz(zq)

    ckdown_suv_r(i)=dzxaoz(x,dzx)/(0.647d0-rrbar-raysoz(x))
    ckdown_suv_f(i)= dzxaoz(x,dzx)/(0.647d0-rrbar-raysoz(x))
    ckup_suv_f(i)=dzxaoz(xstar,dzx)/(0.647d0-rrbar-raysoz(xstar))

  enddo

  ! Calculation under the top of the cloud or the aerosol layer, the adding
  ! Method with multiple diffusion is used
  n = 1 !no O3 overlapping
  do l=k1p1,itopp1
    tau(l,n)=tauc(l) + tauao3(l)
    tauca(l,n)= tauao3(l)
    gasym=gco3(l)
    pioc=pic_o3(l)
    ! In the cloud layers
    if (tauc(l).gt.epsc) then
      call reftra  &
        (pioc, piaero_o3, gasym, gaero_o3, tauc(l), tauao3(l), &
        refx, trax, epsc, 0.d0, mui)

      ref(l,n)=fneray(l)*refx
      tra(l,n)=fneray(l)*trax

    endif
    ! In the aerosol layers
    call reftra  &
      ( 0.d0, piaero_o3, 0.d0, gaero_o3, 0.d0, tauao3(l), &
      refx, trax, epsc, 0.d0, mui)

    ref(l,n)=ref(l,n)+(1.d0-fneray(l))*refx
    tra(l,n)=tra(l,n) + (1.d0 -fneray(l))*trax

    refs(l,n)=ref(l,n)
    tras(l,n)=tra(l,n)

    trard(l,n)=fneray(l)* dexp(-m*(tauc(l)+tauao3(l))) &
      + (1.d0 - fneray(l))*dexp(-m*tauao3(l))

  end do
  ! Top boundary conditions
  tau(itop+1,n)= 0.d0
  tra(itop+1,n)= 1.d0
  trard(itop+1,n)= 1.d0
  trad(itop+1,n)= trard(itop+1,n)
  ref(itop+1,n)=0.d0
  tras(itop+1,n)=tra(itop+1,n)
  refs(itop+1,n)=ref(itop+1,n)
  trat(itop+1,n)=tra(itop+1,n)
  reft(itop+1,n)=ref(itop+1,n)
  refts(itop+1,n)=refs(itop+1,n)
  trats(itop+1,n)=tras(itop+1,n)
  ! Bottom boundary conditions
  tra(k1,n) = 0.d0
  trard(k1,n) = 0.d0
  ref(k1,n) = albe
  tras(k1,n) = 0.d0
  refs(k1,n) =0.0

  ! Downward addition of layers
  do l=itop,k1,-1
    ! Equations 34 of LH74
    drtt1=1.d0-refts(l+1,n)*ref(l,n)
    reft(l,n)=reft(l+1,n)+trat(l+1,n)*ref(l,n) &
      *trats(l+1,n)/drtt1

    trat(l,n)=trat(l+1,n)*tra(l,n)/drtt1

    ! Trad transmission for direct radiation
    trad(l,n) = trad(l+1,n)*trard(l,n)
    if(l.gt.k1) then
      refts(l,n)=refs(l,n)+tras(l,n)*refts(l+1,n)*tra(l,n)/drtt1
      trats(l,n)=trats(l+1,n)*tras(l,n)/drtt1
    end if
  end do

  ! Upward addition of layers
  refb(k1,n)=ref(k1,n)
  refbs(k1,n)=refs(k1,n)
  do l=k1p1,itop
    dtrb1=1.d0-refb(l-1,n)*refs(l,n)
    refb(l,n)=ref(l,n)+tra(l,n)*refb(l-1,n)*tras(l,n)/dtrb1
  end do

  ! Calculation of upward and downward fluxes and absorption
  do l=itopp1,k1p1,-1
    dud1=1.d0-refts(l,n)*refb(l-1,n)
    if(dud1.gt.1.d-30) then

      upw(l,n)=trat(l,n)*refb(l-1,n)/dud1
      dow(l,n)=trat(l,n)/dud1
    else
      upw(l,n)= trat(l,n)*refb(l-1,n)
      dow(l,n)= trat(l,n)
    endif
    dowd(l,n)= trad(l,n)
    atln(l,n)=((1.d0-reft(k1,n))+upw(l,n)-dow(l,n))
  enddo
  do l = itop, k1p1, -1
    absn(l,n) = atln(l,n) - atln(l+1,n)
    ! Fux divergence
    !flux coming from the top
    fabso3(l)=fo*muzero*(0.647d0- raysoz(m*rayuoz(zray(itop))))*absn(l,n)
  enddo
  ! If there is a cloud
  if (itop .gt. k1) then
    do i = k1, itop - 1
      ! addition of ozone absorption for heating in the layers when adding method is used
      zqm1 = zqq(i)
      zq = zqq(i+1)
      if(i.eq.k1) zq = zqm1
      xm1 = m*rayuoz(zqm1)
      x = m*rayuoz(zq)
      zbas = zray(k1)
      xstarm1 = m*rayuoz(zbas) + mbar*(rayuoz(zbas) - rayuoz(zqm1))
      xstar = m*rayuoz(zbas) + mbar*(rayuoz(zbas) - rayuoz(zq))
      !taking into account ozone absorption in the layers with clouds or aerosols
      fabso3(i) =fabso3(i)+ muzero*fo*(raysoz(xm1) - raysoz(x) + rbar                  &
        *(raysoz(xstar) - raysoz(xstarm1)))
      ! fluxes calculation taking into account ozone absorption
      dfso3(i)=muzero*fo*(0.647d0-rrbar-raysoz(x))*dow(i+1,n)
      ddfso3(i)=muzero*fo*(0.647d0-rrbar-raysoz(x))*dowd(i+1,n)
      ufso3(i)=muzero*fo*(0.647d0-rrbar-raysoz(xstar))*upw(i+1,n)

      dzx = drayuoz(zq)
      ! calculation of absorption coefficient ckup and ckdown useful for 3D simulation
      ckdown_suv_r(i)= dzxaoz(x,dzx)/(0.647d0-rrbar-raysoz(x))
      ckdown_suv_f(i)= dzxaoz(x,dzx)/(0.647d0-rrbar-raysoz(x))
      ckup_suv_f(i)=dzxaoz(xstar,dzx)/(0.647d0-rrbar-raysoz(xstar))

    enddo

    ! Calculation of upward flux above cloud or aerosol layers taking into account
    ! the upward flux transmitted by cloud or aerosol layers
    ! if there is no cloud and no aerosol (itop=k1) this term have not to be add
    do i = itop, kmray
      zq = zqq(i)
      zbas = zray(itop-1)
      xstar = mbar*(rayuoz(zbas) - rayuoz(zq))
      ufso3(i) = ufso3(itop-1)*(1.d0 -raysoz(xstar))
    enddo
  endif

  ! 6.4 Absorption by water vapor and liquid water (H20 band, SIR)

  ! In that case we have to solve multiple diffusion. This is achieved by means
  ! of the adding method following Lacis et Hansen, 1974
  ! calculation of reflexivity and transmissivity for each vertical layer

  do n = 1, 8
    do l = k1p1, kmray
      gasym=gch2o(l)
      pioc=pic_h2o(l)
      dqqv = kn(n)*(qqv(l+1) - qqv(l))/10.d0
      if (l.eq.k1p1) dqqv = kn(n)*qqv(l+1)/10.d0

      ! In the cloud  layers
      tau(l,n) = tauc(l) + dqqv + tauah2o(l)

      if(qlray(l).ge.epsc) then
        call reftra &
          (pioc, piaero_h2o, gasym, gaero_h2o, tauc(l) , tauah2o(l), &
          refx, trax, epsc, dqqv, mui)

        ref(l,n)=fneray(l)*refx
        tra(l,n)=fneray(l)*trax + (1.d0 - fneray(l))*dexp(-5.d0*dqqv/3.d0)

        if (iaer.eq.1) then
          call reftra &
            (0.d0, piaero_h2o, 0.d0, gaero_h2o, 0.d0 , tauah2o(l), &
            refx0, trax0, epsc, dqqv, mui)

          ref(l,n) = fneray(l)*refx + (1.d0 - fneray(l))*refx0
          tra(l,n) = fneray(l)*trax + (1.d0 - fneray(l))*trax0
        endif

        refs(l,n) = ref(l,n)
        tras(l,n) = tra(l,n)

        ! trard transmissivity for direct radiation
        trard(l,n) = fneray(l)*dexp(-m*(dqqv+tauc(l)+tauah2o(l))) &
          +(1.d0-fneray(l))*dexp(-m*(dqqv+tauah2o(l)))

      else

        ! in the clear sky layers
        ref(l,n) = 0.d0
        tra(l,n) = dexp(-5.d0*tau(l,n)/3.d0)
        refs(l,n) = ref(l,n)
        tras(l,n) = tra(l,n)

        trard(l,n)=dexp(-m*(dqqv+tauah2o(l)))

        if(l.ge.itopp1) tra(l,n) = dexp(-m*tau(l,n))
        if (iaer.eq.1) then
          call reftra  &
            (0.d0, piaero_h2o, 0.d0, gaero_h2o, 0.d0, tauah2o(l), &
            refx, trax, epsc,dqqv, mui)

          ref(l,n)=fneba(l)*refx
          tra(l,n)=fneba(l)*trax+(1.d0-fneba(l))*dexp(-5.d0*dqqv/3.d0)

        endif

      endif

    enddo

    tau(kmray+1,n) = 0.d0
    tra(kmray+1,n) = 1.d0

    ! For direct radiation
    trard(kmray+1,n) = 1.d0
    trad(kmray+1,n) = trard(kmray+1,n)

    ref(kmray+1,n) = 0.d0
    tras(kmray+1,n) = tra(kmray+1,n)
    refs(kmray+1,n) = ref(kmray+1,n)
    tra(k1,n) = 0.d0
    ref(k1,n) = albe
    tras(k1,n) = 0.d0
    refs(k1,n) = 0.d0

    trat(kmray+1,n) = tra(kmray+1,n)
    reft(kmray+1,n) = ref(kmray+1,n)
    refts(kmray+1,n) = refs(kmray+1,n)
    trats(kmray+1,n) = tras(kmray+1,n)
    fneray(k1) = 0.d0

    ! For direct radiation
    trad(kmray+1,n) = trard(kmray+1,n)
    ! downward addition of layers
    do l = kmray, k1, -1
      drtt1 = 1.d0 - refts(l+1,n)*ref(l,n)
      reft(l,n) = reft(l+1,n) &
        + trat(l+1,n)*ref(l,n)*trats(l+1,n)/drtt1
      trat(l,n) = trat(l+1,n)*tra(l,n)/drtt1

      ! trad for direct radiation
      trad(l,n)=trad(l+1,n)*trard(l,n)

      if(l.gt.k1) then
        refts(l,n) = refs(l,n) &
          + tras(l,n)*refts(l+1,n)*tra(l,n)/drtt1
        trats(l,n) = trats(l+1,n)*tras(l,n)/drtt1
      endif
    enddo

    ! upward layer addition
    refb(k1,n) = ref(k1,n)
    refbs(k1,n) = refs(k1,n)

    do l = k1p1, kmray
      dtrb1 = 1.d0 - refb(l-1,n)*refs(l,n)
      refb(l,n) = ref(l,n) + tra(l,n)*refb(l-1,n)*tras(l,n)/dtrb1
    enddo

    ! calculation of downward and upward fluxes
    do l = kmray+1, k1p1, -1
      dud1 = 1.d0 - refts(l,n)*refb(l-1,n)
      upw(l,n) = trat(l,n)*refb(l-1,n)/dud1
      dow(l,n) = trat(l,n)/dud1

      ! downward fluxes for direct radiation
      dowd(l,n) = trad(l,n)

      !calculation of absorption
      atln(l,n) =  pkn(n)*((1.d0 - reft(k1,n)) + upw(l,n) &
        - dow(l,n))
    enddo

    ! absorption in individual layers
    do l = kmray, k1p1, -1
      absn(l,n) = atln(l,n) - atln(l+1,n)
    enddo
  enddo

  ! summation over frequencies and estimation of absorption integrated
  !  on the whole spectrum
  do l = kmray, k1p1, -1
    fabsh2o(l) = 0.d0
    do n = 1, 8
      fabsh2o(l) = fabsh2o(l)+absn(l,n)
    enddo
    fabsh2o(l) = fabsh2o(l)*fo*muzero
  enddo

  ! In the case with no clouds and aerosol in order to have exactly
  ! the same expressions in 1D and 3D for the heating rate
  if (itop.eq.k1) then
    do i = k1p1, kmray
      y = m*(qqvtot - qqv(i))
      if(i.eq.k1p1) y = m*qqvtot
      yp1 = m*(qqvtot - qqv(i+1))
      ystar = m*qqvtot + 5.d0/3.d0*qqv(i)
      if(i.eq.k1p1) ystar = m*qqvtot
      ystarp1 = m*qqvtot + 5.d0/3.d0*qqv(i+1)
      fabsh2o(i) = muzero*fo*(raysve(y) - raysve(yp1) + albe*(raysve(ystarp1) &
                                                             -raysve(ystar)))
    enddo
  endif
  ! 5.5 heating in the layers
  rayst(k1) = 0.d0
  rayst_h2o(k1) = 0.d0
  rayst_o3(k1) = 0.d0
  do i = k1p1, kmray
    deltaz = zqq(i+1) - zqq(i)
    if(i.eq.k1p1) deltaz = zqq(i+1) - zray(k1)
    cphum = cp0*(1.d0 + (cpvcpa - 1.d0)*qvray(i))
    rayst(i) = (fabsh2o(i) + fabso3(i))/deltaz/romray(i)/cphum

    rayst_h2o(i) = (fabsh2o(i))/deltaz/romray(i)/cphum
    rayst_o3(i) = ( fabso3(i))/deltaz/romray(i)/cphum
  enddo

  ! 5.6 calculation of solar fluxes
  !  for global radiation, fd for direct radiation for the water vapor band

  do i = k1, kmray
    dfsh2o(i) = 0.d0
    ufsh2o(i) = 0.d0
    ddfsh2o(i) = 0.d0

    do n = 2,8
      dfsh2o(i) = dfsh2o(i) + pkn(n)*dow(i+1,n)
      ufsh2o(i) = ufsh2o(i) + pkn(n)*upw(i+1,n)
      ddfsh2o(i) = ddfsh2o(i) + pkn(n)*dowd(i+1,n)
    enddo

    dfsh2o(i) = fo*muzero*dfsh2o(i)
    ufsh2o(i) = fo*muzero*ufsh2o(i)
    ddfsh2o(i) = fo*muzero*ddfsh2o(i)

    y = m*(qqvtot - qqv(i))
    if (i.eq.k1) y = m*qqvtot
    ystar = m*qqvtot + 5.d0/3.d0*qqv(i)
    if(i.eq.k1) ystar = m*qqvtot

    ! to test in the case with no clouds and aerosol in order to have exactly
    ! the same expressions in 1D and 3D for the fluxes
    if (itop.eq.k1) then
      dfsh2o(i) = fo*muzero*(0.353d0-raysve(y))
      ufsh2o(i) = fo*muzero*(0.353d0-raysve(ystar))*albe
      ddfsh2o(i) = fo*muzero*(0.353d0-raysve(y))
    endif

    ! (p_i/p_k) * (p_k/p0) = p_i/p0
    !We keep 1013.15 for standard atmosphere ref:LH74
    corp = preray(i) / 101315.d0
    rov = romray(i)*(qvray(i)*corp*sqrt(tkelvi/(temray(i) + tkelvi)))
    dy = rov
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

    ! direct radiation mod (sum of vapor water band and O3 band)
    drfs(k) = ddfsh2o(k)+ddfso3(k)
    ! diffuse radiation (estmated by difference between global and direct)
    dddsh2o(k) = dfsh2o(k)-ddfsh2o(k)
    dddso3(k) = dfso3(k)-ddfso3(k)
    ddfs(k) = dddsh2o(k)+dddso3(k)

    solu(k,ivertc) = ufs(k)
    sold(k,ivertc) = dfs(k)
  enddo

  ! Note: Mutiplication by transmission function for minor gases
  ! Tmg is now taking into account by fo=fo*Tmg

  ! solar heating of the ground surface by the downward global flux
  fos=dfs(k1)*(1.d0-albe)

  soil_direct_flux=drfs(k1)
  soil_global_flux=dfs(k1)
  soil_direct_flux_h2o=ddfsh2o(k1)
  soil_global_flux_h2o=dfsh2o(k1)
  soil_direct_flux_o3=ddfso3(k1)
  soil_global_flux_o3=dfso3(k1)

  ! Calculation of absorption coefficient ckup and ckdown useful
  ! for 3D simulation
  ! For water vapor without clouds and aerosols downward is only direct,
  ! upward is only diffuse

  ! SIR band
  do k = k1p1, kmray
    tauapc=tauah2o(k)+tauc(k)
    deltaz = zqq(k+1) - zqq(k)
    if(k.eq.k1p1) deltaz=zqq(k+1) - zray(k1)
    if(tauapc.lt.epsc)then
      ckapcd=0.d0
      ckapcf=0.d0
      ck_aero_h2of=0.d0
      ck_aero_h2od=0.d0
      w0_sir(k) = 0.d0
      g_apc_sir(k) = 0.d0
    else
      picapc=(pic_h2o(k)*tauc(k)+piaero_h2o*tauah2o(k))/tauapc
      w0_sir(k) = picapc
      ! if we take into account asymmetry factor for forward diffuse radiation
      ! Note apc means aerosols+clouds
      gapc=(pic_h2o(k)*tauc(k)*gch2o(k)+piaero_h2o*tauah2o(k)*gaero_h2o)/tauapc*picapc
      g_apc_sir(k) = gapc
      ! absorption and forward diffusion
      ckapcf=(1.d0-picapc*(1.d0+gapc)/2.d0)*tauapc/(deltaz*mui)
      ckapcd=tauapc/(deltaz*muzero_cor)
      ck_aero_h2of=(1.d0-piaero_h2o)*tauah2o(k)/(deltaz*mui)
      ck_aero_h2od=tauah2o(k)/(deltaz*muzero_cor)
    endif

    ckup_sir_f(k) = (ckup_sir_f(k) + ck_aero_h2of)*(1.d0-fneray(k)) + &
      (ckup_sir_f(k)+ckapcf)*(fneray(k))
    ckdown_sir_r(k) = (ckdown_sir_r(k) + ck_aero_h2od)*(1.d0-fneray(k)) + &
      (ckdown_sir_r(k)+ckapcd)*(fneray(k))
    ckdown_sir_f(k) = (ckdown_sir_f(k) + ck_aero_h2of)*(1.d0-fneray(k)) + &
      (ckdown_sir_f(k)+ckapcf)*(fneray(k))
  enddo

  ! SUV band
  do k = k1p1, kmray
    tauapc=tauao3(k)+tauc(k)
    deltaz = zqq(k+1) - zqq(k)
    if(k.eq.k1p1) deltaz=zqq(k+1) - zray(k1)
    if(tauapc.lt.epsc)then
      ckapcd=0.d0
      ckapcf=0.d0
      ck_aero_o3f=0.d0
      ck_aero_o3d=0.d0
      w0_suv(k) = 0.d0
      g_apc_suv(k) = 0.d0
    else
      picapc=(pic_o3(k)*tauc(k)+piaero_o3*tauao3(k))/tauapc
      w0_suv(k) = picapc
      ! if we take into account asymmetry factor for forward diffuse radiation
      gapc=(pic_o3(k)*tauc(k)*gco3(k)+piaero_o3*tauao3(k)*gaero_o3)/tauapc*picapc
      g_apc_suv(k) = gapc
      ! absorption and forward diffusion
      ckapcf=(1.d0-picapc*(1.d0+gapc)/2.d0)*tauapc/(deltaz*mui)
      ckapcd=tauapc/(deltaz*muzero_cor)
      ck_aero_o3f=(1.d0-piaero_o3)*tauao3(k)/(deltaz*mui)
      ck_aero_o3d=tauao3(k)/(deltaz*muzero_cor)
    endif

    ckup_suv_f(k) = (ckup_suv_f(k)+ck_aero_o3f)*(1.d0-fneray(k))+&
      (ckup_suv_f(k)+ckapcf)*(fneray(k))
    ckdown_suv_r(k) = (ckdown_suv_r(k) +  ck_aero_o3d)*(1.d0-fneray(k))+&
      (ckdown_suv_r(k)+ckapcd)*(fneray(k))

    ckdown_suv_f(k) = (ckdown_suv_f(k) +  ck_aero_o3f)*(1.d0-fneray(k))+&
      (ckdown_suv_f(k)+ckapcf)*(fneray(k))
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
  soil_direct_flux = 0.d0
  soil_global_flux = 0.d0
  soil_direct_flux_h2o=0.0
  soil_global_flux_h2o=0.0
  soil_direct_flux_o3=0.0
  soil_global_flux_o3=0.0
endif

! Compute Boundary conditions for the 3D (Director diFfuse) Solar radiance
! at the top of the CFD domain
! and the absorption coefficients
call field_get_id_try("spectral_rad_incident_flux", f_id)

is_active = cs_rad_time_is_active()

if (f_id.ge.0.and.is_active) then
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
          dddsh2o, zent, iz1, iz2, var )

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
          dddso3, zent, iz1, iz2, var )

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
          dddso3, zent, iz1, iz2, var )

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
deallocate(reft,refts)
deallocate(refs,tras,trats)
deallocate(refb,trab,upw)
deallocate(refbs,fabso3c,tra)
deallocate(dow,atln,absn)
deallocate(fnebmax,fneba)

deallocate(dfso3c, ufso3c, ddfso3c)
deallocate(dowd, trad, trard)

deallocate(dfsh2o, ufsh2o, dfso3, ufso3, dfs, ufs)
deallocate(dffsh2o, dffso3)
deallocate(ddfsh2o, ddfso3)
deallocate(drfs, dffs, ddfs, dddsh2o, dddso3)
deallocate(ckup_sir_f, ckdown_sir_r, ckdown_sir_f)
deallocate(ckup_suv_f, ckdown_suv_r, ckdown_suv_f)
deallocate(g_apc_suv, g_apc_sir, w0_suv, w0_sir)

return

!===============================================================================

contains

  !-----------------------------------------------------------------------------

  !> \brief Computes ozone concentration for a given altitude

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

  !> \brief Aborption function of the solar radiation by water vapor

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

  !> \brief Aborption derivative-function of the solar radiation by water vapor

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

  !> \brief Aborption function of the solar radiation by ozone

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
