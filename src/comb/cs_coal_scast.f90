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
!> \file cs_coal_scast.f90
!> \brief Specific physic routine: pulverized coal flame
!>   Souce terms have to be precised for a scalar PP
!>   on a step of time
!>
!> \warning  the treatement of source terms is different
!> -------   from that of ustssc.f
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
!> \param[in]     rtp, rtpa     calculated variables at cell centers
!>                               (at current and previous time steps)
!> \param[in,out] propce        physic properties at cell centers
!> \param[in,out] smbrs         explicit second member
!> \param[in,out] rovsdt        implicit diagonal part
!______________________________________________________________________________!

subroutine cs_coal_scast &
 ( iscal  ,                                                       &
   rtpa   , rtp    , propce ,                                     &
   smbrs  , rovsdt )

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use entsor
use optcal
use dimens, only: nvar
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
use field

!===============================================================================

implicit none

! Arguments

integer          iscal

double precision rtp(ncelet,nflown:nvar), rtpa(ncelet,nflown:nvar)
double precision propce(ncelet,*)
double precision smbrs(ncelet), rovsdt(ncelet)

! Local variables

character*80     string
character*80     fname
character*80     name
integer          ivar , iel
integer          numcla , numcha , icla
integer          ipcgch , ipcgd1 , ipcgd2 , ipcght , ipcsec
integer          ipghco2 , ipghh2o
integer          ixchcl , ixckcl
integer          ipcro2 , ipcte1 , ipcte2 , ipcvsl
integer          ipcdia
integer          mode, ige
integer          ipcx2c , icha , ii, jj
integer          itermx,nbpauv,nbrich,nbepau,nberic,ipghc2
integer          iterch,nbpass,nbarre,nbimax
integer          iexp1 , iexp2 , iexp3
integer          iok1
integer          keyccl

double precision xnuss
double precision aux, rhovst
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
double precision xcom,xo2m,xkcequ,xkpequ

double precision xw1,xw2,xw3,xw4
double precision xo2,wmel,wmhcn,wmno,wmo2,wmnh3
double precision gmdev1(ncharm),gmdev2(ncharm),gmhet(ncharm)
double precision aux1 , aux2 , aux3
double precision xch,xck,xash,xmx2
double precision tfuelmin,tfuelmax
double precision auxdev,auxht3,auxco2,auxh2o,auxwat

double precision, dimension (:), allocatable :: w1,w2,w3,w4,w5
double precision, dimension (:), allocatable :: tfuel
double precision, dimension(:), pointer :: gamvlei, gamvloi, xchcpi, gaheto2i, xckcpi
double precision, dimension(:), pointer :: gaseci, frmcpi, agei, gahetco2i
double precision, dimension(:), pointer :: gaheth2oi, xwtcpi, xacpip
double precision, dimension(:), pointer ::  crom
double precision, dimension(:), pointer :: cpro_cp
double precision, dimension(:), pointer :: cvara_k, cvara_ep

!Lacal variables
!===============
!
! Pointers of variables of state
! ------------------------------
! (Temperature of a particle of iclas,
!  Pointer on CHx1,
!  Pointer on CHx2,
!  Mass density of the gaseous phase,
!  Mass density of the mixture)
integer ipctem, ipcyf1, ipcyf2, idgaz
!
! Kinetic laminar constant
! ------------------------
! (NH3 + O2,
!  NH3 + NO,
!  reburning)
integer iexp4,iexp5,iexprb
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
!===============================================================================
!
!===============================================================================
! 1. Initialization
!===============================================================================
! --- Scalar number to treat: iscal
! --- Number of the associated to the scalar to treat variable iscal
ivar = isca(iscal)

! --- Name of the variable associated to the scalar to treat iscal
call field_get_label(ivarfl(ivar), chaine)

! --- Number of the physic quantity
call field_get_val_s(icrom, crom)
if (icp.gt.0) call field_get_val_s(iprpfl(icp), cpro_cp)
ipcte1 = ipproc(itemp1)

call field_get_val_prev_s(ivarfl(ik), cvara_k)
call field_get_val_prev_s(ivarfl(iep), cvara_ep)

! Key id of the coal scalar class
call field_get_key_id("scalar_class", keyccl)

!===============================================================================
! Deallocation dynamic arrays
!----
allocate(w1(1:ncel),w2(1:ncel),w3(1:ncel),w4(1:ncel),w5(1:ncel),stat=iok1)
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
!===============================================================================

!===============================================================================
! 2. Consideration of the source terms for relative variables
!    to the classes of particles
!===============================================================================

! --> Source term for the mass fraction of reactive coal

if ( ivar.ge.isca(ixch(1)) .and. ivar.le.isca(ixch(nclacp)) ) then

  if (iwarni(ivar).ge.1) then
    write(nfecra,1000) chaine(1:8)
  endif
  numcla = ivar-isca(ixch(1))+1
  ipcgch = ipproc(igmdch(numcla))

  do iel = 1, ncel

    ! ---- Calculation of W1 = - rho.GMDCH > 0

    xw1 = - crom(iel)*propce(iel,ipcgch)*volume(iel)

    ! ---- Calculation of explicit and implicit parts of source terms

    rovsdt(iel) = rovsdt(iel) + max(xw1,zero)
    smbrs (iel) = smbrs(iel)  - xw1*rtpa(iel,ivar)

  enddo

endif


! --> Source term for the mass fraction of coke
if ( ivar.ge.isca(ixck(1)) .and. ivar.le.isca(ixck(nclacp)) ) then

  if (iwarni(ivar).ge.1) then
    write(nfecra,1000) chaine(1:8)
  endif

  numcla = ivar-isca(ixck(1))+1
  ipcgch = ipproc(igmdch(numcla))
  ipcgd1 = ipproc(igmdv1(numcla))
  ipcgd2 = ipproc(igmdv2(numcla))
  ixchcl = isca(ixch(numcla))
  ixckcl = isca(ixck(numcla))
  ipcght = ipproc(igmhet(numcla))
  if ( ihtco2 .eq. 1 ) then
    ipghco2 = ipproc(ighco2(numcla))
  endif
  if ( ihth2o .eq. 1 ) then
    ipghh2o = ipproc(ighh2o(numcla))
  endif

  do iel = 1, ncel

    ! ---- Calculation of W1 = - rho.Xch.GMDCH.Volume > 0

    xw1=-crom(iel)*rtp(iel,ixchcl)*propce(iel,ipcgch)     &
                                           *volume(iel)

    ! AE : We consider rtp(iel,ixchcl) and not rtpa(iel,ixchcl) to
    !      preserve the mass.

    ! ---- Calculation of W2 = rho.Xch.(GMDV1+GMDV2)Volume < 0

    xw2 = crom(iel)*rtp(iel,ixchcl)                  &
         *(propce(iel,ipcgd1)+propce(iel,ipcgd2))*volume(iel)

    if ( rtpa(iel,ixckcl) .gt. epsicp ) then

      ! Reaction C(s) + O2 ---> 0.5CO
      ! =============================

      ! ---- Calculation of the implicit part > 0 of source terms
      !      relative to gmhet

      xw3 = -2.d0/3.d0*crom(iel)*propce(iel,ipcght) &
           /(rtpa(iel,ixckcl))**(1.d0/3.d0)*volume(iel)

      ! ---- Calculation of the explicit part < 0 of source terms
      !      relative to gmhet

      xw4 = crom(iel)*propce(iel,ipcght)             &
                * (rtpa(iel,ixckcl))**(2.d0/3.d0)*volume(iel)

    else
      xw3 = 0.d0
      xw4 = 0.d0
    endif

    ! ---- Calculation of the explicit and implicit parts of source terms

    rovsdt(iel) = rovsdt(iel) + max(xw3,zero)
    smbrs(iel)  = smbrs(iel)  + xw1 + xw2 + xw4

  enddo

  if ( ihtco2 .eq. 1 ) then

    do iel = 1, ncel

      if ( rtpa(iel,ixckcl) .gt. epsicp ) then

        ! Reaction C(s) + CO2 ---> 2CO
        ! =============================

        ! ---- Calculation of the implicit part > 0 of source terms
        !      relative to gmhet

        xw3 = -2.d0/3.d0*crom(iel)*propce(iel,ipghco2)     &
              /(rtpa(iel,ixckcl))**(1.d0/3.d0)*volume(iel)

        ! ---- Calculation of the explicit part < 0 of source terms
        !      relative to gmhet

        xw4 = crom(iel)*propce(iel,ipghco2)                 &
             *(rtpa(iel,ixckcl))**(2.d0/3.d0)*volume(iel)

      else
        xw3 = 0.d0
        xw4 = 0.d0
      endif

      ! ---- Calculation of explicit and implicit parts of source terms

      rovsdt(iel) = rovsdt(iel) + max(xw3,zero)
      smbrs(iel)  = smbrs(iel)  + xw4

    enddo

  endif

  if ( ihth2o .eq. 1 ) then

    do iel = 1, ncel

      if ( rtpa(iel,ixckcl) .gt. epsicp ) then

        ! Reaction C(s) + CO2 ---> 2CO
        ! =============================

        ! ---- Calculation of the explicit part > 0 of source terms
        !      relative to gmhet

        xw3 = -2.d0/3.d0*crom(iel)*propce(iel,ipghh2o)     &
              /(rtpa(iel,ixckcl))**(1.d0/3.d0)*volume(iel)

        ! ---- Calculation of the explicit part < 0 of source terms
        !      relative to gmhet

        xw4 = crom(iel)*propce(iel,ipghh2o)                &
             *(rtpa(iel,ixckcl))**(2.d0/3.d0)*volume(iel)

      else
        xw3 = 0.d0
        xw4 = 0.d0
      endif

      ! ---- Calculation of the explicit and implicit part of source terms

      rovsdt(iel) = rovsdt(iel) + max(xw3,zero)
      smbrs(iel)  = smbrs(iel)  + xw4

    enddo

  endif

endif


! --> Source term for the mass fraction of water

if ( ippmod(iccoal) .eq. 1 ) then

  if ( ivar.ge.isca(ixwt(1)) .and.                                &
       ivar.le.isca(ixwt(nclacp)) ) then

    if (iwarni(ivar).ge.1) then
      write(nfecra,1000) chaine(1:8)
    endif
    numcla = ivar-isca(ixwt(1))+1
    numcha = ichcor(numcla)

    ipcsec = ipproc(igmsec(numcla))
    ipcx2c = ipproc(ix2(numcla))

    do iel = 1, ncel

     ! ---- Calculation of explicit and implicit parts of source term

     if ( rtpa(iel,ivar).gt. epsicp .and.                         &
          xwatch(numcha).gt. epsicp       ) then
       xw1 = crom(iel)*propce(iel,ipcsec)*volume(iel)    &
            *(1.d0/propce(iel,ipcx2c))*(1.d0/xwatch(numcha))

       rovsdt(iel) = rovsdt(iel) + max(xw1,zero)
       smbrs(iel)  = smbrs(iel)  - xw1*rtpa(iel,ivar)
     endif

    enddo

  endif

endif

if (i_coal_drift.eq.1) then

  call field_get_name(ivarfl(ivar), fname)

  ! Particle age source term
  if (fname(1:10).eq.'x_age_coal') then

    ! Index of the coal particle class
    call field_get_key_int(ivarfl(ivar), keyccl, icla)

    ! Array values at previous time step
    call field_get_val_prev_s_by_name(fname, xacpip)

    ! Light volatile's source term
    write(name,'(a,i2.2)') 'm_transfer_v1_coal_', icla
    call field_get_val_s_by_name(name, gamvlei)

    ! Heavy volatile's source term
    write(name,'(a,i2.2)') 'm_transfer_v2_coal_', icla
    call field_get_val_s_by_name(name, gamvloi)

    ! Coal mass fraction of icla
    write(name,'(a,i2.2)') 'x_coal_', icla
    call field_get_val_s_by_name(name, xchcpi)

    ! Heterogeneous combustion temporary scale
    write(name,'(a,i2.2)') 'het_ts_o2_coal_', icla
    call field_get_val_s_by_name(name, gaheto2i)

    ! Char mass fraction of icla
    write(name,'(a,i2.2)') 'w_ck_coal_', icla
    call field_get_val_s_by_name(name, xckcpi)

    ! Coal indicator of class icla
    numcha = ichcor(icla)

    ! Dryer temporary scale
    if (ippmod(iccoal) .eq. 1) then
      write(name,'(a,i2.2)') 'dry_ts_coal_', icla
      call field_get_val_s_by_name(name, gaseci)
    endif

    ! Mass fraction of solid phase
    write(name,'(a,i2.2)') 'w_solid_coal_', icla
    call field_get_val_s_by_name(name, frmcpi)

    ! Particle age per cell
    write(name,'(a,i2.2)') 'age_coal_', icla
    call field_get_val_s_by_name(name, agei)

    ! Gazefiation temporary scale by CO2
    if (ihtco2 .eq. 1) then
      write(name,'(a,i2.2)') 'het_ts_co2_coal_', icla
      call field_get_val_s_by_name(name, gahetco2i)
    endif

    ! Gazefiation temporary scale by H2O
    if (ihth2o .eq. 1) then
      write(name,'(a,i2.2)') 'het_ts_h2o_coal_', icla
      call field_get_val_s_by_name(name, gaheth2oi)
    endif

    if (ippmod(iccoal).eq.1) then
      write(name,'(a,i2.2)') 'xwt_coal_', icla
      call field_get_val_s_by_name(name, xwtcpi)
    endif

    do iel = 1, ncel
      ! Mass flux: Devolatilization
      auxdev =  -(gamvlei(iel)+gamvloi(iel))*xchcpi(iel)
      ! Coal consumpsion by heterogeneous combustion
      if (xckcpi(iel) .gt. epsicp) then
        auxht3 = -gaheto2i(iel) * (xckcpi(iel))**(2.d0/3.d0)
      else
        auxht3 = 0.d0
      endif
      ! Mass flux: Gazefication by CO2
      if (ihtco2 .eq. 1) then
        if (xckcpi(iel) .gt. epsicp ) then
          auxco2 = -gahetco2i(iel)* (xckcpi(iel))**(2.d0/3.d0)
        else
          auxco2 = 0.d0
        endif
      else
        auxco2 = 0.d0
      endif
      ! Mass flux: Gazefication by H2O
      if (ihth2o .eq. 1) then
        if (xckcpi(iel) .gt. epsicp) then
          auxh2o = -gaheth2oi(iel)*(xckcpi(iel))**(2.d0/3.d0)
        else
          auxh2o = 0.d0
        endif
      else
        auxh2o = 0.d0
      endif
      !  Mass flux: Drying
      if (ippmod(iccoal) .eq. 1) then
        if (xwtcpi(iel).gt.epsicp .and.xwatch(numcha).gt.epsicp) then
          auxwat = gaseci(iel)/(frmcpi(iel)*xwatch(numcha))
        else
          auxwat = 0.d0
        endif
      else
        auxwat = 0.d0
      endif

      if (frmcpi(iel).gt.epsicp) then
         smbrs(iel) =  smbrs(iel) + (crom(iel) * volume(iel) *          &
                       (frmcpi(iel)                                    -          &
                       ((auxdev+auxht3+auxco2+auxh2o+auxwat)*xacpip(iel))) )
         rovsdt(iel) = rovsdt(iel) + max(((auxdev+auxht3+auxco2+auxh2o+auxwat)/   &
                                         frmcpi(iel)), zero)
      endif

    enddo
  endif

  ! Gas age source term
  if (fname(1:9).eq.'x_age_gas') then

    ! Loop over particle classes
    do icla = 1, nclacp
      ! Light volatile's source term
      write(name,'(a,i2.2)') 'm_transfer_v1_coal_', icla
      call field_get_val_s_by_name(name, gamvlei)

      ! Heavy volatile's source term
      write(name,'(a,i2.2)') 'm_transfer_v2_coal_', icla
      call field_get_val_s_by_name(name, gamvloi)

      ! Coal mass fraction of icla
      write(name,'(a,i2.2)') 'x_coal_', icla
      call field_get_val_s_by_name(name, xchcpi)

      ! Heterogeneous combustion temporary scale
      write(name,'(a,i2.2)') 'het_ts_o2_coal_', icla
      call field_get_val_s_by_name(name, gaheto2i)

      ! Coal mass fraction of icla
      write(name,'(a,i2.2)') 'xck_cp_', icla
      call field_get_val_s_by_name(name, xckcpi)

      ! Coal indicator of class icla
      numcha = ichcor(icla)

      ! Drying temporary scale
      if (ippmod(iccoal) .eq. 1) then
        write(name,'(a,i2.2)') 'dry_ts_coal_', icla
        call field_get_val_s_by_name(name, gaseci)
      endif

      ! Mass fraction of the solid phase
      write(name,'(a,i2.2)') 'w_solid_coal_', icla
      call field_get_val_s_by_name(name, frmcpi)

      ! Particle age by cell
      write(name,'(a,i2.2)') 'age_coal_', icla
      call field_get_val_s_by_name(name, agei)

      ! Gazefication temporary scale by CO2
      if (ihtco2 .eq. 1) then
        write(name,'(a,i2.2)') 'het_ts_co2_coal_', icla
        call field_get_val_s_by_name(name, gahetco2i)
      endif

      ! Gazefication temporary scale by H2O
      if (ihth2o .eq. 1) then
        write(name,'(a,i2.2)') 'het_ts_h2o_coal_', icla
        call field_get_val_s_by_name(name, gaheth2oi)
      endif

      if (ippmod(iccoal).eq.1) then
        write(name,'(a6,i2.2)') 'xwt_coal_', icla
        call field_get_val_s_by_name(name, xwtcpi)
      endif

      do iel = 1, ncel
        ! Mass flux: Devolatilization
        auxdev = -(gamvlei(iel)+gamvloi(iel))*xchcpi(iel)
        ! Coal consumption by heterogeneous combustion
        if (xckcpi(iel) .gt. epsicp) then
          auxht3 = -gaheto2i(iel) * (xckcpi(iel))**(2.d0/3.d0)
        else
          auxht3 = 0.d0
        endif
        ! Mass flux: Gazefication by CO2
        if (ihtco2 .eq. 1) then
          if (xckcpi(iel) .gt. epsicp ) then
            auxco2 = -gahetco2i(iel)* (xckcpi(iel))**(2.d0/3.d0)
          else
            auxco2 = 0.d0
          endif
        else
          auxco2 = 0.d0
        endif
        ! Mass flux: Gazefication by H2O
        if (ihth2o .eq. 1) then
          if (xckcpi(iel) .gt. epsicp) then
            auxh2o = -gaheth2oi(iel)*(xckcpi(iel))**(2.d0/3.d0)
          else
            auxh2o = 0.d0
          endif
        else
          auxh2o = 0.d0
        endif
        !  Mass flux: Drying
        if (ippmod(iccoal) .eq. 1) then
          if (xwtcpi(iel).gt.epsicp .and.xwatch(numcha).gt.epsicp) then
            auxwat = gaseci(iel)/(frmcpi(iel)*xwatch(numcha))
          else
            auxwat = 0.d0
          endif
        else
          auxwat = 0.d0
        endif

        smbrs(iel) = smbrs(iel) + crom(iel)*volume(iel)               &
                                *( ( auxdev+auxht3+auxco2+auxh2o+auxwat)       &
                                   * agei(iel)                                 &
                                 - frmcpi(iel)                                 &
                                 )
      enddo

    enddo

    ! Finalization
    do iel = 1, ncel
      ! The formula is (1- Sum X2)
      smbrs(iel) =  smbrs(iel) + crom(iel) * volume(iel)
    enddo

  endif

endif

! --> Source term for the enthalpy of the solid

if ( ivar.ge.isca(ih2(1)) .and. ivar.le.isca(ih2(nclacp)) ) then

  if (iwarni(ivar).ge.1) then
    write(nfecra,1000) chaine(1:8)
  endif

  numcla = ivar-isca(ih2(1))+1
  numcha = ichcor(numcla)
  ixchcl = isca(ixch(numcla))
  ixckcl = isca(ixck(numcla))
  ipcx2c = ipproc(ix2(numcla))
  ipcro2 = ipproc(irom2(numcla ))
  ipcdia = ipproc(idiam2(numcla))
  ipcte2 = ipproc(itemp2(numcla))
  ipcght = ipproc(igmhet(numcla))
  if ( ihtco2 .eq. 1 ) then
    ipghco2 = ipproc(ighco2(numcla))
  endif
  if ( ihth2o .eq. 1 ) then
    ipghh2o = ipproc(ighh2o(numcla))
  endif
  ipcgd1 = ipproc(igmdv1(numcla))
  ipcgd2 = ipproc(igmdv2(numcla))
  ipcgch = ipproc(igmdch(numcla))
  ixchcl = isca(ixch(numcla))
  ixckcl = isca(ixck(numcla))

  ! ---- Contribution to the explicit and implicit balance
  !        exchanges by molecular distribution
  !        6 Lambda Nu / diam**2 / Rho2 * Rho * (T1-T2)

  ! ------ Calculation of lambda in W1

  xnuss = 2.d0
  do iel = 1, ncel
    if ( ivisls(iscalt).gt.0 ) then
      ipcvsl = ipproc(ivisls(iscalt))
      if (icp.gt.0) then
        w1(iel) = propce(iel,ipcvsl) * cpro_cp(iel)
      else
        w1(iel) = propce(iel,ipcvsl) * cp0
      endif
    else
      if (icp.gt.0) then
        w1(iel) = visls0(iscalt) * cpro_cp(iel)
      else
        w1(iel) = visls0(iscalt) * cp0
      endif
    endif
  enddo

  ! ------ Calculation of diameter of the particles in W2
  !        d20 = (A0.D0**2+(1-A0)*DCK**2)**0.5

  do iel = 1, ncel
    w2(iel) = ( xashch(numcha)*diam20(numcla)**2 +                &
                (1.d0-xashch(numcha))*propce(iel,ipcdia)**2       &
              )**0.5
  enddo

  ! ------ Contribution to the explicit and implicit balance
  !        exchanges by molecular distribution
  !      Remark : We use propce(iel,IPCX2C) because we want X2 at the iteration n

  do iel = 1, ncel
    if ( propce(iel,ipcx2c) .gt. epsicp ) then
      aux         = 6.d0 * w1(iel) * xnuss / w2(iel)**2           &
                  / propce(iel,ipcro2) * crom(iel)       &
                  * propce(iel,ipcx2c) * volume(iel)
      rhovst      = aux / cp2ch(numcha) /propce(iel,ipcx2c)

      smbrs(iel)  = smbrs(iel)-aux*(propce(iel,ipcte2)-propce(iel,ipcte1))
      rovsdt(iel) = rovsdt(iel)+ max(zero,rhovst)
    endif

  enddo


  ! ---- Contribution to the explicit and implicit balances
  !        of exchange term of energy between phases:
  !        gama(dev1) H(mv1,T2)+gama(dev2) H(mv2,T2)

  do iel = 1, ncel

    !        Gama Dev1 et Gama Dev2

    gamdv1 = crom(iel)*rtp(iel,ixchcl)                   &
            *propce(iel,ipcgd1)

    gamdv2 = crom(iel)*rtp(iel,ixchcl)                   &
            *propce(iel,ipcgd2)

    !        H(mv1,T2)

    do ige = 1, ngazem
      coefe(ige) = zero
    enddo

    den = a1(numcha)*wmole(ichx1c(numcha))+b1(numcha)*wmole(ico)  &
         +c1(numcha)*wmole(ih2o)          +d1(numcha)*wmole(ih2s) &
         +e1(numcha)*wmole(ihcn)          +f1(numcha)*wmole(inh3)
    coefe(ichx1) = a1(numcha)*wmole(ichx1c(numcha)) /den
    coefe(ico  ) = b1(numcha)*wmole(ico)            /den
    coefe(ih2o ) = c1(numcha)*wmole(ih2o)           /den
    coefe(ih2s ) = d1(numcha)*wmole(ih2s)           /den
    coefe(ihcn ) = e1(numcha)*wmole(ihcn)           /den
    coefe(inh3 ) = f1(numcha)*wmole(inh3)           /den

    t2         = propce(iel,ipcte2)
    do icha = 1, ncharm
      f1mc(icha) = zero
      f2mc(icha) = zero
    enddo
    f1mc(numcha) = 1.d0

    mode      = -1
    call cs_coal_htconvers1 &
    !======================
    ( mode , xhdev1 , coefe , f1mc , f2mc , t2 )

    !        H(mv2,T2)

    do ige = 1, ngazem
      coefe(ige) = zero
    enddo
    den = a2(numcha)*wmole(ichx2c(numcha))+b2(numcha)*wmole(ico)  &
         +c2(numcha)*wmole(ih2o)          +d2(numcha)*wmole(ih2s) &
         +e2(numcha)*wmole(ihcn)          +f2(numcha)*wmole(inh3)
    coefe(ichx2) = a2(numcha)*wmole(ichx2c(numcha)) /den
    coefe(ico  ) = b2(numcha)*wmole(ico)            /den
    coefe(ih2o ) = c2(numcha)*wmole(ih2o)           /den
    coefe(ih2s ) = d2(numcha)*wmole(ih2s)           /den
    coefe(ihcn ) = e2(numcha)*wmole(ihcn)           /den
    coefe(inh3 ) = f2(numcha)*wmole(inh3)           /den

    t2         = propce(iel,ipcte2)
    do icha = 1, ncharm
      f1mc(icha) = zero
      f2mc(icha) = zero
    enddo
    f2mc(numcha) = 1.d0

    mode      = -1
    call cs_coal_htconvers1 &
    ( mode , xhdev2 , coefe , f1mc , f2mc , t2 )

    !         Contribution to explicit and implicit balances

    smbrs(iel) = smbrs(iel)+(gamdv1*xhdev1+gamdv2*xhdev2)*volume(iel)

  enddo

  ! ------ Heterogeneous combustion: C(s) + 02 ---> 0.5 C0
  !        GamHET * (28/12 H(CO,T2)-16/12 H(O2,T1) )

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

    t2        = propce(iel,ipcte2)
    mode      = -1
    call cs_coal_htconvers1 &
    ( mode , xhco , coefe , f1mc , f2mc , t2 )

    !        Calculation of HO2(T1)

    do ige = 1, ngazem
      coefe(ige) = zero
    enddo
    coefe(io2) = 1.d0
    do icha = 1, ncharm
      f1mc(icha) = zero
      f2mc(icha) = zero
    enddo

    t1        = propce(iel,ipcte1)
    mode      = -1
    call cs_coal_htconvers1 &
    ( mode , xho2 , coefe , f1mc , f2mc , t1 )

    !         Contribution to explicit and implicit balances

    if ( rtpa(iel,ixckcl) .gt. epsicp ) then

      gamhet = crom(iel)*propce(iel,ipcght)              &
               * ( (rtpa(iel,ixckcl))**(2.d0/3.d0) +              &
                2.d0/3.d0*(rtp(iel,ixckcl)-rtpa(iel,ixckcl))      &
                 /(rtpa(iel,ixckcl))**(1.d0/3.d0) )

    else
      gamhet = 0.d0
    endif

   smbrs(iel) = smbrs(iel)                                        &
                 +gamhet                                          &
                  *(28.d0/12.d0*xhco-16.d0/12.d0*xho2)            &
                  *volume(iel)

  enddo

  ! ------ Heterogeneous combustion: C(s) + C02 ---> 2 C0
  !        GamHET * (56/12 H(CO,T2)-44/12 H(CO2,T1) )

  if ( ihtco2 .eq. 1 ) then
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

      t2        = propce(iel,ipcte2)
      mode      = -1
      call cs_coal_htconvers1 &
      ( mode , xhco , coefe , f1mc , f2mc , t2  )

      !        Calculation of HCO2(T1)

      do ige = 1, ngazem
        coefe(ige) = zero
      enddo
      coefe(ico2) = 1.d0
      do icha = 1, ncharm
        f1mc(icha) = zero
        f2mc(icha) = zero
      enddo

      t1        = propce(iel,ipcte1)
      mode      = -1
      call cs_coal_htconvers1 &
      !======================
      ( mode , xhco2 , coefe , f1mc , f2mc , t1    )

      !         Contribution to explicit and implicit balances

      if ( rtpa(iel,ixckcl) .gt. epsicp ) then

        gamhet = crom(iel)*propce(iel,ipghco2)           &
                 * ( (rtpa(iel,ixckcl))**(2.d0/3.d0) +            &
                2.d0/3.d0*(rtp(iel,ixckcl)-rtpa(iel,ixckcl))      &
                 /(rtpa(iel,ixckcl))**(1.d0/3.d0) )

      else
        gamhet = 0.d0
      endif

     smbrs(iel) = smbrs(iel)                                      &
                   +gamhet                                        &
                    *(56.d0/12.d0*xhco-44.d0/12.d0*xhco2)         &
                    *volume(iel)

    enddo

  endif
  ! ------ Heterogeneous combustion: C(s) + H2O ---> CO + H2
  !        GamHET * (28/12 H(CO,T2)+2/12 H(HY,T2) -18/12 H(H2O,T1) )

  if ( ihth2o .eq. 1 ) then
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

      t2        = propce(iel,ipcte2)
      mode      = -1
      call cs_coal_htconvers1 &
      ( mode , xhco , coefe , f1mc , f2mc , t2 )

      !        Calculation of HH2(T2)

      do ige = 1, ngazem
        coefe(ige) = zero
      enddo
      coefe(ihy) = 1.d0
      do icha = 1, ncharm
        f1mc(icha) = zero
        f2mc(icha) = zero
      enddo

      t2        = propce(iel,ipcte2)
      mode      = -1
      call cs_coal_htconvers1 &
      ( mode , xhh2 , coefe , f1mc , f2mc , t2    )

      !        Calculation of HH2O(T1)

      do ige = 1, ngazem
        coefe(ige) = zero
      enddo
      coefe(ih2o) = 1.d0
      do icha = 1, ncharm
        f1mc(icha) = zero
        f2mc(icha) = zero
      enddo

      t1        = propce(iel,ipcte1)
      mode      = -1
      call cs_coal_htconvers1 &
      ( mode , xhh2o , coefe , f1mc , f2mc , t1    )

      !         Contribution to explicit and implicit balances

      if ( rtpa(iel,ixckcl) .gt. epsicp ) then

        gamhet = crom(iel)*propce(iel,ipghh2o)           &
                 * ( (rtpa(iel,ixckcl))**(2.d0/3.d0) +            &
                2.d0/3.d0*(rtp(iel,ixckcl)-rtpa(iel,ixckcl))      &
                 /(rtpa(iel,ixckcl))**(1.d0/3.d0) )

      else
        gamhet = 0.d0
      endif

     smbrs(iel) = smbrs(iel)                                    &
                 +gamhet                                        &
                  *(28.d0/12.d0*xhco+ 2.d0/12.d0*xhh2           &
                                    -18.d0/12.d0*xhh2o )        &
                    *volume(iel)

    enddo

  endif

  !       --> Source term on H2 (coming from drying)

  if ( ippmod(iccoal) .eq. 1 ) then

    ! ---- Contribution of source term interfacial to explicit and implicit balances


    ipcsec = ipproc(igmsec(numcla))
    ipcte2 = ipproc(itemp2(numcla))

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

      t2 = propce(iel,ipcte2)
      if ( t2 .gt. 100.d0+tkelvi ) then
        t2 = 100.d0+tkelvi
      endif
      mode      = -1
      call cpthp1 &
      !==========
      ( mode , hh2ov , coefe , f1mc , f2mc ,  t2  )

!         Contribution to explicit balance

      if ( rtpa(iel,isca(ixwt(numcla))).gt. epsicp .and.          &
           xwatch(numcha) .gt. epsicp       ) then

        aux = -crom(iel)*propce(iel,ipcsec)              &
       *(rtp(iel,isca(ixwt(numcla)))/propce(iel,ipcx2c))          &
       *(1.d0                    /xwatch(numcha))                 &
             *hh2ov

      else
        aux = 0.d0
      endif

      smbrs(iel) = smbrs(iel) + aux*volume(iel)

    enddo

  endif
endif

!===============================================================================
! 3. Taking into account source terms for relative variables in the mixture
!===============================================================================

! --> Source term for light volatile materials

if ( ivar.ge.isca(if1m(1)) .and. ivar.le.isca(if1m(ncharb)) ) then

  if (iwarni(ivar).ge.1) then
    write(nfecra,1000) chaine(1:8)
  endif

! ---- Calculation of GMDEV1 = - Sum (rho.XCH.GMDV1) > 0  --> W1

  numcha = ivar-isca(if1m(1))+1
  do iel = 1, ncel
    w1(iel) = zero
  enddo
  do icla = 1, nclacp
    ipcgd1 = ipproc(igmdv1(icla))
    ixchcl = isca(ixch(icla))
    if ( ichcor(icla).eq.numcha ) then
      do iel = 1, ncel
        w1(iel) = w1(iel) - crom(iel)*rtp(iel,ixchcl)    &
                * propce(iel,ipcgd1)
      enddo
    endif
  enddo

! ---- Contribution of interfacial source term to explicit and implicit balances

  do iel = 1, ncel
    smbrs(iel)  = smbrs(iel)  + volume(iel) * w1(iel)
  enddo

endif


! --> Source terms for heavy volatil materials

if ( ivar.ge.isca(if2m(1)) .and. ivar.le.isca(if2m(ncharb)) ) then

  if (iwarni(ivar).ge.1) then
    write(nfecra,1000) chaine(1:8)
  endif

! ---- Calculation of GMDEV2 = - Sum (rho.XCH.GMDV2) > 0 --> W1

  numcha = ivar-isca(if2m(1))+1
  do iel = 1, ncel
    w1(iel) = zero
  enddo
  do icla = 1, nclacp
    ipcgd2 = ipproc(igmdv2(icla))
    ixchcl = isca(ixch(icla))
    if ( ichcor(icla).eq.numcha ) then
      do iel = 1, ncel
        w1(iel) = w1(iel) - crom(iel)*rtp(iel,ixchcl)    &
                * propce(iel,ipcgd2)
      enddo
    endif
  enddo

! ---- Contribution of interfacial source term to explicite balance

  do iel = 1, ncel
    smbrs(iel)  = smbrs(iel)  + volume(iel) * w1(iel)
  enddo

endif


! --> Source term for the tracer 7 (O2) (heterogeneous combustion by C)

if ( ivar.eq.isca(if7m) ) then

  ! Remark: We take the same source term than for Xck
  !                  to be conservative

  if (iwarni(ivar).ge.1) then
    write(nfecra,1000) chaine(1:8)
  endif

  do iel = 1, ncel
    w1(iel) = zero
  enddo

  do icla = 1, nclacp
    ipcght = ipproc(igmhet(icla))
    ixckcl = isca(ixck(icla))
    do iel = 1, ncel
      if ( rtpa(iel,ixckcl) .gt. epsicp ) then
        w1(iel) = w1(iel)                                         &
                 - crom(iel)*propce(iel,ipcght)          &
                 * ( (rtpa(iel,ixckcl))**(2.d0/3.d0) +            &
                    2.d0/3.d0*(rtp(iel,ixckcl)-rtpa(iel,ixckcl))  &
                    /(rtpa(iel,ixckcl))**(1.d0/3.d0) )
      endif
    enddo
  enddo

  ! ---- Contribution of interfacial source term to explicit and implicit balances

  do iel = 1, ncel
    smbrs(iel)  = smbrs(iel)  + volume(iel) * w1(iel)
  enddo

endif



! --> Source term for the tracer 8 (CO2) (heterogeneous combustion by C)

if ( ihtco2 .eq. 1 ) then
  if ( ivar.eq.isca(if8m) ) then

    ! Remark: We take the same source term than for Xck
    !                  to be conservative

    if (iwarni(ivar).ge.1) then
      write(nfecra,1000) chaine(1:8)
    endif

    do iel = 1, ncel
      w1(iel) = zero
    enddo

    do icla = 1, nclacp
      ipcght = ipproc(ighco2(icla))
      ixckcl = isca(ixck(icla))
      do iel = 1, ncel
        if ( rtpa(iel,ixckcl) .gt. epsicp ) then
          w1(iel) = w1(iel)                                        &
                   - crom(iel)*propce(iel,ipcght)         &
                   * ( (rtpa(iel,ixckcl))**(2.d0/3.d0) +           &
                      2.d0/3.d0*(rtp(iel,ixckcl)-rtpa(iel,ixckcl)) &
                      /(rtpa(iel,ixckcl))**(1.d0/3.d0) )
        endif
      enddo
    enddo

    ! ---- Contribution of interfacial source term to explicit and implicit balances

    do iel = 1, ncel
      smbrs(iel)  = smbrs(iel)  + volume(iel) * w1(iel)
    enddo

  endif

endif

! --> Source term for the tracer 9 (H2O) (heterogeneous combustion by H2O)

if ( ihth2o .eq. 1 ) then
  if ( ivar.eq.isca(if9m) ) then

    ! Remark: We take the same source term than for Xck
    !                  to be conservative


    if (iwarni(ivar).ge.1) then
      write(nfecra,1000) chaine(1:8)
    endif

    do iel = 1, ncel
      w1(iel) = zero
    enddo

    do icla = 1, nclacp
      ipcght = ipproc(ighh2o(icla))
      ixckcl = isca(ixck(icla))
      do iel = 1, ncel
        if ( rtpa(iel,ixckcl) .gt. epsicp ) then
          w1(iel) = w1(iel)                                        &
                   - crom(iel)*propce(iel,ipcght)        &
                   * ( (rtpa(iel,ixckcl))**(2.d0/3.d0) +           &
                      2.d0/3.d0*(rtp(iel,ixckcl)-rtpa(iel,ixckcl)) &
                      /(rtpa(iel,ixckcl))**(1.d0/3.d0) )
        endif
      enddo
    enddo

    ! ---- Contribution of interfacial source term to explicit and implicit balances

    do iel = 1, ncel
      smbrs(iel)  = smbrs(iel)  + volume(iel) * w1(iel)
    enddo

  endif

endif


! --> Source term for the fuel variance

if ( ivar.eq.isca(ifvp2m) ) then

  if (iwarni(ivar).ge.1) then
    write(nfecra,1000) chaine(1:8)
  endif

  call cs_coal_fp2st                                               &
 !==================
 ( iscal  ,                                                        &
   rtpa   , rtp    , propce ,                                      &
   smbrs  , rovsdt )

endif

! --> Source term for the tracer 6 (Water coming from drying)

if ( ippmod(iccoal) .eq. 1 ) then

  if ( ivar.eq.isca(if6m) ) then


    if (iwarni(ivar).ge.1) then
      write(nfecra,1000) chaine(1:8)
    endif

    ! ---- Contribution of interfacial source term to explicit and implicit balances

    do iel = 1, ncel
      w1(iel) = zero
    enddo

    do icla = 1, nclacp

      ipcsec = ipproc(igmsec(icla))
      ipcx2c = ipproc(ix2(icla))
      numcha = ichcor(icla)

      do iel = 1, ncel

        if (  rtpa(iel,isca(ixwt(icla))).gt. epsicp               &
            .and.                                                 &
              xwatch(numcha) .gt. epsicp       ) then

          w1(iel) = w1(iel)                                       &
      + crom(iel)*propce(iel,ipcsec)                     &
       *(rtp(iel,isca(ixwt(icla)))/propce(iel,ipcx2c))            &
       *(1.d0                  /xwatch(numcha))

        endif

      enddo

    enddo

    do iel = 1, ncel
      smbrs(iel)  = smbrs(iel)  + volume(iel) * w1(iel)
    enddo

  endif

endif

! --> Source term for CO2

if ( ieqco2 .eq. 1 ) then

  if ( ivar.eq.isca(iyco2) ) then

    if (iwarni(ivar).ge.1) then
      write(nfecra,1000) chaine(1:8)
    endif

    ! ---- Contribution of interfacial source term to explicit and implicit balances

    ! Oxydation of CO
    ! ===============

    !  Dryer Glassman : XK0P in (mol/m3)**(-0.75) s-1
    !          XK0P = 1.26D10
    !           XK0P = 1.26D7 * (1.1)**(NTCABS)
    !           IF ( XK0P .GT. 1.26D10 ) XK0P=1.26D10
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

     xxco  = propce(iel,ipproc(iym1(ico  )))/wmole(ico)           &
            *propce(iel,ipproc(irom1))
     xxo2  = propce(iel,ipproc(iym1(io2  )))/wmole(io2)           &
            *propce(iel,ipproc(irom1))
     xxco2 = propce(iel,ipproc(iym1(ico2 )))/wmole(ico2)          &
            *propce(iel,ipproc(irom1))
     xxh2o = propce(iel,ipproc(iym1(ih2o )))/wmole(ih2o)          &
            *propce(iel,ipproc(irom1))

     xxco  = max(xxco ,zero)
     xxo2  = max(xxo2 ,zero)
     xxco2 = max(xxco2,zero)
     xxh2o = max(xxh2o,zero)
     sqh2o = sqrt(xxh2o)

     xkp = exp(lnk0p-t0p/propce(iel,ipproc(itemp1)))
     xkm = exp(lnk0m-t0m/propce(iel,ipproc(itemp1)))

     xkpequ = 10.d0**(l10k0e-t0e/propce(iel,ipproc(itemp1)))
     xkcequ = xkpequ                                              &
             /sqrt(8.32d0*propce(iel,ipproc(itemp1))/1.015d5)

     !        initialization by the transported state

     anmr  = xxco2
     xcom  = xxco + xxco2
     xo2m  = xxo2 + 0.5d0*xxco2

     if ( propce(iel,ipproc(itemp1)) .gt. 1200.d0 ) then

      !        Search for the equilibrum state
      !        Iterative search without control of convergence
      !         (to keep the parallelisation on meshes)
      !        On the numver of moles of separating reaction
      !         the state before reaction (such as calculated by Cpcym)
      !         of the equilibrum state
      !        anmr has to be the border between 0 and Min(XCOM,2.*XO2M)
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
           oxydation
       xden = xkp*sqh2o*(xxo2)**0.25d0
     else
           dissociation
       xden = xkm
     endif
     if ( xden .ne. 0.d0 ) then

       tauchi = 1.d0/xden
       tautur = cvara_k(iel)/cvara_ep(iel)

       x2 = 0.d0
       do icla = 1, nclacp
         x2 = x2 + propce(iel,ipproc(ix2(icla)))
       enddo

       !    We transport CO2

       smbrs(iel)  = smbrs(iel)                                   &
                    +wmole(ico2)/propce(iel,ipproc(irom1))        &
         * (xco2eq-xxco2)/(tauchi+tautur)                         &
         * (1.d0-x2)                                              &
         * volume(iel) * crom(iel)

       w1(iel) = volume(iel)*crom(iel)/(tauchi+tautur)
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

   write(nfecra,*) ' Max Error = ',ERR1MX
   write(nfecra,*) ' no Points   ',NBERIC,NBARRE,NBPASS
   write(nfecra,*) ' Iter max number ',NBIMAX

   !     Source term: heterogeneous combustion by CO2

   if ( ihtco2 .eq. 1) then

     do iel = 1, ncel

       aux = 0.d0
       do icla = 1,nclacp

         ixckcl = isca(ixck(icla))
         ipghc2 = ipproc(ighco2(icla))

         aux = aux                                                &
              + crom(iel)*propce(iel,ipghc2)             &
               *(rtpa(iel,ixckcl))**(2.d0/3.d0)*volume(iel)

       enddo

       rovsdt(iel) = rovsdt(iel) - aux*(wmole(ico2)/0.012)

     enddo

   endif

  endif

endif

! --> Source term for Enth_Ox
!                       HCN and NO: only from the second iteration

if ( ieqnox .eq. 1 .and. ntcabs .gt. 1) then

  ! Terms on Oxydant enthalpy

  if ( ivar .eq. isca(ihox) ) then

    !  Calculation of T2 average on particles

    tfuelmin = 1.d+20
    tfuelmax =-1.d+20
    do iel=1,ncel

      xmx2 = 0.d0
      do icla = 1, nclacp
        xck  = rtp(iel,isca(ixck(icla)))
        xch  = rtp(iel,isca(ixch(icla)))
        xash = rtp(iel,isca(inp (icla)))*xmash(icla)
        xmx2   = xmx2 + xch + xck + xash

        !   Taking into account humidity

        if ( ippmod(iccoal) .eq. 1 ) then
          xmx2 = xmx2+rtp(iel,isca(ixwt(icla)))
        endif
      enddo

      if ( xmx2 .gt. 0.d0 ) then
        tfuel(iel) = 0.d0
        do icla=1,nclacp
          ipcte2=ipproc(itemp2(icla))
          tfuel(iel) = tfuel(iel)                                              &
                      +( rtp(iel,isca(ixck(icla)))                             &
                        +rtp(iel,isca(ixch(icla)))                             &
                        +rtp(iel,isca(inp (icla)))*xmash(icla) )               &
                        *propce(iel,ipcte2)

          !  Taking into account humidity

          if ( ippmod(iccoal) .eq. 1 ) then
            tfuel(iel) = tfuel(iel) + (rtp(iel,isca(ixwt(icla))))              &
                         *propce(iel,ipcte2)
          endif
        enddo

        tfuel(iel) = tfuel(iel)/xmx2

      else
        tfuel(iel) = propce(iel,ipcte1)
      endif
      tfuelmin = min(tfuel(iel),tfuelmin)
      tfuelmax = max(tfuel(iel),tfuelmax)
    enddo
    if (irangp .ge. 0 ) then
      call parmin(tfuelmin)
      call parmax(tfuelmax)
    endif
    write(nfecra,*) ' Min max de Tfuel pour Hoxy ',tfuelmin,tfuelmax

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
      mode      = -1
      call cs_coal_htconvers1 &
    ( mode  , xhco    , coefe  , f1mc   , f2mc   ,  t2    )

      !  Calculation of HO2(T1)

      do ige = 1, ngazem
        coefe(ige) = zero
      enddo
      coefe(io2) = 1.d0
      do icha = 1, ncharm
        f1mc(icha) = zero
        f2mc(icha) = zero
      enddo
      t1        = propce(iel,ipcte1)
      mode      = -1
      call cs_coal_htconvers1 &
      ( mode  , xho2    , coefe  , f1mc   , f2mc   , t1    )

      do icla=1,nclacp
        ixckcl = isca(ixck(icla))
        ipcght = ipproc(igmhet(icla))
        if ( rtpa(iel,ixckcl) .gt. epsicp ) then
          gamhet = crom(iel)*propce(iel,ipcght)              &
                   * ( (rtpa(iel,ixckcl))**(2.d0/3.d0) +              &
                  2.d0/3.d0*(rtp(iel,ixckcl)-rtpa(iel,ixckcl))        &
                    /(rtpa(iel,ixckcl))**(1.d0/3.d0) )
        else
          gamhet = 0.d0
        endif
        smbrs(iel) = smbrs(iel)                                        &
                    -gamhet                                            &
                     *(28.d0/12.d0*xhco-16.d0/12.d0*xho2)*volume(iel)

      enddo

    enddo

    !  Heterogeneous combustion: C + CO2 ---> 2 CO

    !  Calculation of HO2(T1)

    if ( ihtco2 .eq. 1 ) then
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
        mode      = -1
        call cs_coal_htconvers1 &
        ( mode  , xhco    , coefe  , f1mc   , f2mc   , t2    )


        !  Calculation of HCO2(T1)

        do ige = 1, ngazem
          coefe(ige) = zero
        enddo
        coefe(ico2) = 1.d0
        do icha = 1, ncharm
          f1mc(icha) = zero
          f2mc(icha) = zero
        enddo
        t1        = propce(iel,ipcte1)
        mode      = -1
        call cs_coal_htconvers1 &
        ( mode  , xhco2    , coefe  , f1mc   , f2mc   , t1    )

        do icla=1,nclacp
          ixckcl  = isca(ixck(icla))
          ipghco2 = ipproc(ighco2(icla))
          if ( rtpa(iel,ixckcl) .gt. epsicp ) then
            gamhet = crom(iel)*propce(iel,ipghco2)              &
                     * ( (rtpa(iel,ixckcl))**(2.d0/3.d0) +               &
                    2.d0/3.d0*(rtp(iel,ixckcl)-rtpa(iel,ixckcl))         &
                      /(rtpa(iel,ixckcl))**(1.d0/3.d0) )
          else
            gamhet = 0.d0
          endif
          smbrs(iel) = smbrs(iel)                                          &
                      -gamhet                                              &
                       *(56.d0/12.d0*xhco-44.d0/12.d0*xhco2) *volume(iel)

        enddo

      enddo

    endif

    !   Heterogeneous combustion: C + H2O ---> CO + H2

    !     Calculation of HO2(T1)

    if ( ihth2o .eq. 1 ) then

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
        mode      = -1
        call cs_coal_htconvers1 &
        ( mode  , xhco    , coefe  , f1mc   , f2mc   ,  t2    )

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
        mode      = -1
        call cs_coal_htconvers1 &
        ( mode  , xhh2    , coefe  , f1mc   , f2mc   , t2    )

        !       Calculation of HH2O(T1)

        do ige = 1, ngazem
          coefe(ige) = zero
        enddo
        coefe(ih2o) = 1.d0
        do icha = 1, ncharm
          f1mc(icha) = zero
          f2mc(icha) = zero
        enddo
        t1        = propce(iel,ipcte1)
        mode      = -1
        call cs_coal_htconvers1 &
      ( mode  , xhh2o    , coefe  , f1mc   , f2mc   , t1    )

        do icla=1,nclacp
          ixckcl  = isca(ixck(icla))
          ipghh2o = ipproc(ighh2o(icla))
          if ( rtpa(iel,ixckcl) .gt. epsicp ) then
            gamhet = crom(iel)*propce(iel,ipghh2o)           &
                     * ( (rtpa(iel,ixckcl))**(2.d0/3.d0) +            &
                    2.d0/3.d0*(rtp(iel,ixckcl)-rtpa(iel,ixckcl))      &
                      /(rtpa(iel,ixckcl))**(1.d0/3.d0) )
          else
            gamhet = 0.d0
          endif
          smbrs(iel) = smbrs(iel)                                           &
                      -gamhet                                               &
                       *(28.d0/12.d0*xhco+ 2.d0/12.d0*xhh2                  &
                                         -18.d0/12.d0*xhh2o ) *volume(iel)

        enddo

      enddo

    endif

    !   Drying

    if ( ippmod(iccoal) .eq. 1 ) then

      do icla=1,nclacp

        numcha = ichcor(icla)

        ipcsec = ipproc(igmsec(icla))
        ipcte2 = ipproc(itemp2(icla))
        ipcx2c = ipproc(ix2(icla))

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

          t2 = propce(iel,ipcte2)
          mode      = -1
          call cpthp1 &
          ( mode  , hh2ov    , coefe  , f1mc   , f2mc   ,t2    )

          !  Contribution to explicit balance

          if ( rtpa(iel,isca(ixwt(icla))).gt. epsicp .and.          &
               xwatch(numcha) .gt. epsicp       ) then

            aux = crom(iel)*propce(iel,ipcsec)             &
                 *(rtp(iel,isca(ixwt(icla)))/propce(iel,ipcx2c))    &
                 *(1.d0                    /xwatch(numcha))         &
                 *hh2ov

          else
            aux = 0.d0
          endif

          smbrs(iel) = smbrs(iel) - aux*volume(iel)

        enddo

      enddo

    endif

  endif
endif

if ( ieqnox .eq. 1 .and. imdnox.eq.0 .and. ntcabs .gt. 1) then

  ! Source terms on Y_HCN and Y_NO

  if ( ivar.eq.isca(iyhcn) .or. ivar.eq.isca(iyno) ) then

    !  Pointer source terms

    iexp1  = ipproc(ighcn1)
    iexp2  = ipproc(ighcn2)
    iexp3  = ipproc(ignoth)


    !  Mass molar

    wmhcn = wmole(ihcn)
    wmno  = 0.030d0
    wmo2  = wmole(io2  )

    if ( ivar.eq.isca(iyhcn) ) then

    !  Source term HCN

      if (iwarni(ivar).ge.1) then
        write(nfecra,1000) chaine(1:8)
      endif

      do iel=1,ncel
        wmel=propce(iel,ipproc(immel))
        xo2= propce(iel,ipproc(iym1(io2)))*wmel/wmo2

        aux = volume(iel)*crom(iel)                      &
             *(propce(iel,iexp2)+propce(iel,iexp1)                &
             *rtpa(iel,isca(iyno))                                &
             *propce(iel,ipproc(immel))                           &
             /wmno)

        smbrs(iel)  = smbrs(iel)  - aux*rtpa(iel,ivar)
        rovsdt(iel) = rovsdt(iel) + aux

        do icha=1,ncharb
          gmdev1(icha)=0.d0
          gmdev2(icha)=0.d0
          gmhet (icha)=0.d0
        enddo

        do icla=1,nclacp

          icha = ichcor(icla)

          gmdev1(icha) = gmdev1(icha)                              &
               +propce(iel,ipproc(igmdv1(icla)))                   &
               *crom(iel)                                 &
               *rtpa(iel,isca(ixch(icla)))
          gmdev2(icha) = gmdev2(icha)                              &
               +propce(iel,ipproc(igmdv2(icla)))                   &
               *crom(iel)                                 &
               *rtpa(iel,isca(ixch(icla)))
          gmhet(icha) = gmhet(icha)                                &
               +propce(iel,ipproc(igmhet(icla)))                   &
               *crom(iel)                                 &
               *rtpa(iel,isca(ixck(icla)))**(2.d0/3.d0)

        enddo

        do icha=1,ncharb
          !  % of pure nitrogen in the coal

          aux = -volume(iel)*fn(icha)*wmhcn/(wmole(in2)/2.d0)        &
                            *(qpr(icha)*(gmdev1(icha)+gmdev2(icha)))
          if(xo2.gt.0.03d0) then
            aux=aux-volume(iel)*fn(icha)*wmhcn/(wmole(in2)/2.d0)     &
                               * (1.d0-qpr(icha)*y2ch(icha))         &
                                / (1-y2ch(icha))*gmhet(icha)         &
                                * (1.d0-xashch(icha))
          endif
          smbrs(iel)  = smbrs(iel) + aux
        enddo

      enddo

    endif

    if ( ivar.eq.isca(iyno) ) then

      !  Source term NO

      if (iwarni(ivar).ge.1) then
        write(nfecra,1000) chaine(1:8)
      endif

      do iel=1,ncel

        wmel=propce(iel,ipproc(immel))

        aux1 = volume(iel)*crom(iel)                     &
              *propce(iel,iexp1)*rtpa(iel,isca(iyhcn))            &
              *propce(iel,ipproc(immel))/wmhcn
        aux2 = volume(iel)*crom(iel)                     &
              *propce(iel,iexp2)*rtpa(iel,isca(iyhcn))            &
              *wmno/wmhcn
        aux3 = volume(iel)*crom(iel)**1.5d0              &
              *propce(iel,iexp3)                                  &
              *propce(iel,ipproc(iym1(in2)))

        smbrs(iel)  = smbrs(iel) - aux1*rtpa(iel,ivar)            &
                                 + aux2 + aux3
        rovsdt(iel) = rovsdt(iel) + aux1
      enddo

    endif

  endif

endif

if ( ieqnox .eq. 1 .and. imdnox.eq.1 .and. ntcabs .gt. 1) then

  ! Source terms on Y_HCN and Y_NO

  if ( ivar.eq.isca(iyhcn) .or. ivar.eq.isca(iyno) .or. ivar.eq.isca(iynh3)    &
     ) then

    ! Pointer Source terms NO gas phase
    iexp1  = ipproc(ighcn1)
    iexp2  = ipproc(ighcn2)
    iexp3  = ipproc(ignoth)
    iexp4  = ipproc(ignh31)
    iexp5  = ipproc(ignh32)
    iexprb = ipproc(igrb)
    ! Pointer on CHx1 and CHx2
    ipcyf1 = ipproc(iym1(1))
    ipcyf2 = ipproc(iym1(2))
    idgaz  = ipproc(irom1)

    ! Mass molar

    wmhcn = wmole(ihcn)
    wmno  = 0.030d0
    wmnh3 = wmole(inh3)

    if ( ivar.eq.isca(iyhcn) ) then

    aux = 0.d0

      do iel = 1,ncel

         propce(iel,ipproc(ifhcnr))  = zero
         propce(iel,ipproc(ifhcnd))  = zero
         propce(iel,ipproc(ifhcnc)) = zero

      enddo

      !  Source term HCN

      if (iwarni(ivar).ge.1) then
        write(nfecra,1000) chaine(1:8)
      endif

      do iel=1,ncel

        !  Mass molar of the gas mixture
        wmel=propce(iel,ipproc(immel))


        !  Coefficient of reactions HCN + O2 et HCN + NO
        aux = volume(iel)*crom(iel)                                   &
             *(propce(iel,iexp2)+propce(iel,iexp1)                             &
             *rtpa(iel,isca(iyno))                                             &
             *wmel/wmno)

        smbrs(iel)  = smbrs(iel)  - aux*rtpa(iel,ivar)
        rovsdt(iel) = rovsdt(iel) + aux

        !  Reburning ?
        !  Chen's model
        if(irb.eq.1) then

          do icha = 1,ncharb

             ychx = ( propce(iel,ipcyf1) * wmel/wmchx1 )                       &
                  + ( propce(iel,ipcyf2) * wmel/wmchx2 )

             aux = volume(iel)*wmhcn*propce(iel,iexprb)                        &
                 * rtpa(iel,isca(iyno))*wmel/wmno                              &
                 * ychx

             smbrs(iel)  = smbrs(iel)  + aux

             propce(iel,ipproc(ifhcnr)) = propce(iel,ipproc(ifhcnr)) + aux

          enddo

        !  Dimitiou's model
        elseif(irb.eq.2) then

          do icha = 1,ncharb

            !  Reburning by CHx1
            if(propce(iel,ipcyf1).gt.0.d0) then

              !  Number of point of the temperature discretization
              do ii = 1,7

                !  We look for the interval teno(ii) < Tgaz < teno(ii+1)
                if(propce(iel,ipcte1).ge.teno(ii).and.propce(iel,ipcte1).lt.   &
                   teno(ii+1)) then

                  !  JJ indicates the quotient H/C of the fuel (4=CH4;3=CH3,etc.)
                  do jj = 1,4

                    !  We look for the interval jj < chx1(icha) < jj + 1
                    if(chx1(icha).ge.4.d0) then

                    core1 = ka(4,ii) + ( (ka(4,ii+1)- ka(4,ii))/(teno(ii+1)-   &
                            teno(ii)) ) * (propce(iel,ipcte1) - teno(ii))
                    core2 = kb(4,ii) + ( (kb(4,ii+1)- kb(4,ii))/(teno(ii+1)-   &
                            teno(ii)) ) * (propce(iel,ipcte1) - teno(ii))
                    para2 = chi2(ii) + ( (chi2(ii+1)-chi2(ii))/(teno(ii+1) -   &
                            teno(ii)) ) * (propce(iel,ipcte1) - teno(ii))

                    elseif(chx1(icha).le.1.d0) then

                    core1 = ka(1,ii) + ( (ka(1,ii+1)- ka(1,ii))/(teno(ii+1)-   &
                            teno(ii)) ) * (propce(iel,ipcte1) - teno(ii))
                    core2 = kb(1,ii) + ( (kb(1,ii+1)- kb(1,ii))/(teno(ii+1)-   &
                            teno(ii)) ) * (propce(iel,ipcte1) - teno(ii))
                    para2 = chi2(ii) + ( (chi2(ii+1)-chi2(ii))/(teno(ii+1) -   &
                            teno(ii)) ) * (propce(iel,ipcte1) - teno(ii))

                    elseif (chx1(icha).ge.jj.and.chx1(icha).lt.jj+1) then

                    core1 = ka(jj,ii) + ( (ka(jj+1,ii+1)- ka(jj,ii))/          &
                            (teno(ii+1)-teno(ii)) ) * (propce(iel,ipcte1)      &
                            - teno(ii))
                    core2 = kb(jj,ii) + ( (kb(jj+1,ii+1)- kb(jj,ii))/          &
                            (teno(ii+1)-teno(ii)) ) * (propce(iel,ipcte1) -    &
                            teno(ii))
                    para2 = chi2(ii) + ( (chi2(ii+1)-chi2(ii))/(teno(ii+1)-    &
                            teno(ii)) ) * (propce(iel,ipcte1) - teno(ii))

                    endif

                  enddo

                endif

              enddo

              if(chx1(icha).ge.3.d0) then

              aux = ( volume(iel)*wmhcn )                                      &
                  * ( (core1 + core2) * para2 )                                &
                  * ( rtp(iel,isca(iyno))*crom(iel)/wmno )                     &
                  * ( propce(iel,ipcyf1)*propce(iel,idgaz)/wmchx1 )

              else

              aux = ( volume(iel)*wmhcn )                                      &
                  * ( core1 + core2 )                                          &
                  * ( rtp(iel,isca(iyno))*crom(iel)/wmno )                     &
                  * ( propce(iel,ipcyf1)*propce(iel,idgaz)/wmchx1 )

              endif

              smbrs(iel)  = smbrs(iel)  + aux

              propce(iel,ipproc(ifhcnr)) = propce(iel,ipproc(ifhcnr)) + aux

            elseif(propce(iel,ipcyf2).gt.0.d0) then

              !  Reburning by CHx2

              !  Number of points of the temperature discretization
              do ii = 1,7

                !  We look for the interval teno(ii) < Tgaz < teno(ii+1)
                if(propce(iel,ipcte1).ge.teno(ii).and.propce(iel,ipcte1).lt.   &
                   teno(ii+1)) then

                  !  JJ indicates the quotient H/C of
                  !           the fuel (4=CH4;3=CH3,etc.)
                  do jj = 1,4

                    !  We look for the interval jj < chx1(icha) < jj + 1
                    if(chx2(icha).ge.4.d0) then

                    core1 = ka(4,ii) + ( (ka(4,ii+1)- ka(4,ii))/(teno(ii+1)-   &
                            teno(ii)) ) * (propce(iel,ipcte1) - teno(ii))
                    core2 = kb(4,ii) + ( (kb(4,ii+1)- kb(4,ii))/(teno(ii+1)-   &
                            teno(ii)) ) * (propce(iel,ipcte1) - teno(ii))
                    para2 = chi2(ii) + ( (chi2(ii+1)-chi2(ii))/(teno(ii+1) -   &
                            teno(ii)) ) * (propce(iel,ipcte1) - teno(ii))

                    elseif(chx2(icha).le.1.d0) then

                    core1 = ka(1,ii) + ( (ka(1,ii+1)- ka(1,ii))/(teno(ii+1)-   &
                            teno(ii)) ) * (propce(iel,ipcte1) - teno(ii))
                    core2 = kb(1,ii) + ( (kb(1,ii+1)- kb(1,ii))/(teno(ii+1)-   &
                            teno(ii)) ) * (propce(iel,ipcte1) - teno(ii))
                    para2 = chi2(ii) + ( (chi2(ii+1)-chi2(ii))/(teno(ii+1) -   &
                            teno(ii)) ) * (propce(iel,ipcte1) - teno(ii))

                    elseif (chx2(icha).ge.jj.and.chx2(icha).lt.jj+1) then

                    core1 = ka(jj,ii) + ( (ka(jj+1,ii+1)- ka(jj,ii))/          &
                            (teno(ii+1)-teno(ii)) ) * (propce(iel,ipcte1) -    &
                            teno(ii))
                    core2 = kb(jj,ii) + ( (kb(jj+1,ii+1)- kb(jj,ii))/          &
                            (teno(ii+1)-teno(ii)) ) * (propce(iel,ipcte1) -    &
                            teno(ii))
                    para2 = chi2(ii) + ( (chi2(ii+1)-chi2(ii))/(teno(ii+1) -   &
                            teno(ii)) ) * (propce(iel,ipcte1) - teno(ii))

                    endif

                  enddo

                endif

              enddo

              if(chx2(icha).ge.3.d0) then

              aux = ( volume(iel)*wmhcn )                                      &
                  * ( (core1 + core2) * para2 )                                &
                  * ( rtp(iel,isca(iyno))*crom(iel)/wmno )                     &
                  * ( propce(iel,ipcyf2)*propce(iel,idgaz)/wmchx2 )

              else

              aux = ( volume(iel)*wmhcn )                                      &
                  * ( core1 + core2 )                                          &
                  * ( rtp(iel,isca(iyno))*crom(iel)/wmno )                     &
                  * ( propce(iel,ipcyf2)*propce(iel,idgaz)/wmchx2 )

              endif

              smbrs(iel)  = smbrs(iel)  + aux

              propce(iel,ipproc(ifhcnr)) = propce(iel,ipproc(ifhcnr)) + aux

            endif

          enddo

        endif

        !  Initialization of variables
        do icha=1,ncharb

          gmdev1(icha)=0.d0
          gmdev2(icha)=0.d0
          gmhet (icha)=0.d0

        enddo

        do icla=1,nclacp

          icha   = ichcor(icla)
          ipctem = ipproc(itemp2(icla))
          ixckcl = isca(ixck(icla))
          ipcght = ipproc(igmhet(icla))

          mckcl1 = (1.d0-y1ch(icha))*a1ch(icha)                                &
                   *exp(-e1ch(icha)/(rr*propce(iel,ipctem)))

          mckcl2 = (1.d0-y2ch(icha))*a2ch(icha)                                &
                   *exp(-e2ch(icha)/(rr*propce(iel,ipctem)))

          !  Forming rate of the first pyrolisis reaction
          gmdev1(icha) = gmdev1(icha)                                          &
               +propce(iel,ipproc(igmdv1(icla)))                               &
               *crom(iel)                                             &
               *rtpa(iel,isca(ixch(icla)))

          !  Forming rate of the second pyrolisis reaction
          gmdev2(icha) = gmdev2(icha)                                          &
               +propce(iel,ipproc(igmdv2(icla)))                               &
               *crom(iel)                                             &
               *rtpa(iel,isca(ixch(icla)))

          if ( rtpa(iel,ixckcl) .gt. epsicp ) then
            !  Reaction rate of the heterogeneous combustion
            gmhet(icha) = gmhet(icha)                                          &
               +propce(iel,ipcght)                                             &
               *crom(iel)                                             &
               *( rtpa(iel,ixckcl)*((1.d0/(mckcl2/mckcl1+1.d0))*yhcnc1(icha)   &
                   +(1.d0/(mckcl1/mckcl2+1.d0))*yhcnc2(icha)) )**(2.d0/3.d0)
          endif

        enddo

        !  Modified source term (new model of NOx)

        do icha=1,ncharb

           !  Release of HCN during devolatilization
           aux = -volume(iel)*(gmdev1(icha)*yhcnle(icha)                       &
                 +gmdev2(icha)*yhcnlo(icha))

           !   Release of HCN during the heterogeneous combustion according
           !   to the value repnck(icha)

           aux = aux-volume(iel)*gmhet(icha)

           smbrs(iel)  = smbrs(iel) + aux

           !  Source terms displaying
           propce(iel,ipproc(ifhcnd)) = propce(iel,ipproc(ifhcnd))             &
           -volume(iel)*(gmdev1(icha)*yhcnle(icha)+gmdev2(icha)*yhcnlo(icha))

           propce(iel,ipproc(ifhcnc))= propce(iel,ipproc(ifhcnc))              &
           -volume(iel)*gmhet(icha)

        enddo

      enddo

    endif


    if ( ivar.eq.isca(iynh3)) then

       aux = 0.d0

      do iel = 1, ncel

         propce(iel,ipproc(ifnh3d)) = zero

      enddo

      !  Source term NH3

      if (iwarni(ivar).ge.1) then
        write(nfecra,1000) chaine(1:8)
      endif

      do iel=1,ncel

        !  Mass molar of the gaseous mixture
        wmel = propce(iel,ipproc(immel))

        !  Coefficient of reactions NH3 + O2 and NH3 + NO
        aux  =   volume(iel)*crom(iel)                                &
             * ( propce(iel,iexp4) + propce(iel,iexp5)                         &
             *   rtpa(iel,isca(iyno))*wmel/wmno )

        smbrs(iel)  = smbrs(iel)  - aux*rtpa(iel,ivar)
        rovsdt(iel) = rovsdt(iel) + aux

        !  Initialization of variables
        do icha=1,ncharb

          gmdev1(icha)=0.d0
          gmdev2(icha)=0.d0
          gmhet (icha)=0.d0

        enddo

        do icla=1,nclacp

          icha = ichcor(icla)

          !  Forming rate of the first pyrolisis reaction
          gmdev1(icha) = gmdev1(icha)                                          &
               +propce(iel,ipproc(igmdv1(icla)))                               &
               *crom(iel)                                             &
               *rtpa(iel,isca(ixch(icla)))

          !  Forming rate of the second pyrolisis reaction
          gmdev2(icha) = gmdev2(icha)                                          &
               +propce(iel,ipproc(igmdv2(icla)))                               &
               *crom(iel)                                             &
               *rtpa(iel,isca(ixch(icla)))

        enddo

        do icha=1,ncharb

           !  Release of NH3 during the devolatization.
           aux = -volume(iel)*(gmdev1(icha)*ynh3le(icha)                       &
                 +gmdev2(icha)*ynh3lo(icha))

           smbrs(iel)  = smbrs(iel) + aux

           !  Source terms displaying
           propce(iel,ipproc(ifnh3d)) = propce(iel,ipproc(ifnh3d))             &
           -volume(iel)*(gmdev1(icha)*ynh3le(icha)+gmdev2(icha)*ynh3lo(icha))

           propce(iel,ipproc(ifnh3c)) = zero

        enddo



      enddo

    endif

    if ( ivar.eq.isca(iyno) ) then

      do iel = 1 ,ncel

        propce(iel,ipproc(icnorb))  = zero
        propce(iel,ipproc(ifnoch)) = zero

      enddo

      !  Source term NO

      if (iwarni(ivar).ge.1) then
        write(nfecra,1000) chaine(1:8)
      endif

      do iel=1,ncel

        !  Mass molar of the gaseous mixture
        wmel=propce(iel,ipproc(immel))

        !  Coefficient of reaction HCN + NO
        aux1 = volume(iel)*crom(iel)                                  &
              *propce(iel,iexp1)*rtpa(iel,isca(iyhcn))                         &
              *wmel/wmhcn

        propce(iel,ipproc(icnohc)) = aux1*rtpa(iel,ivar)

        !  Coefficient of reaction HCN + O2
        aux2 = volume(iel)*crom(iel)                                  &
              *propce(iel,iexp2)*rtpa(iel,isca(iyhcn))                         &
              *wmno/wmhcn

        propce(iel,ipproc(ifnohc)) = aux2

        !  Coefficient of thermal NO
        aux3 = volume(iel)*crom(iel)**1.5d0                           &
              *propce(iel,iexp3)                                               &
!       Pourquoi la fraction massique d'azote n'a ete pas transforme dans une
!       fraction molaire ?
              *propce(iel,ipproc(iym1(in2)))

        propce(iel,ipproc(ifnoth)) = aux3

        !  Coefficient of reaction NH3 + O2 --> NO + ...
        aux4 = volume(iel)*crom(iel)                                  &
              *propce(iel,iexp4)*rtpa(iel,isca(iynh3))                         &
              *wmno/wmnh3

        propce(iel,ipproc(ifnonh)) = aux4

        !  Coefficient of reaction NH3 + NO --> N2 + ...
        aux5 = volume(iel)*crom(iel)                                  &
              *propce(iel,iexp5)*rtpa(iel,isca(iynh3))                         &
              *wmel/wmnh3

        propce(iel,ipproc(icnonh)) = aux5*rtpa(iel,ivar)

        !  Reburning ?
        !  Chen's model
        if(irb.eq.1) then

          do icha = 1,ncharb

             ychx = ( propce(iel,ipcyf1) * wmel/wmchx1 )                       &
                  + ( propce(iel,ipcyf2) * wmel/wmchx2 )

             aux = volume(iel)*wmhcn*propce(iel,iexprb)                        &
                 * rtpa(iel,isca(iyno)) * wmel/wmno  * ychx

             smbrs(iel)  = smbrs(iel)  - aux

             propce(iel,ipproc(icnorb)) = propce(iel,ipproc(icnorb)) + aux

          enddo

        !  Dimitiou's model
        elseif(irb.eq.2) then

          do icha = 1,ncharb

            !  Reburning by CHx1
            if (propce(iel,ipcyf1).gt.0.d0) then

              !  Number of point of the temperature discretization
              do ii = 1,7

                !  We look for the interval teno(ii) < Tgaz < teno(ii+1)
                if(propce(iel,ipcte1).ge.teno(ii).and.propce(iel,ipcte1).lt.   &
                   teno(ii+1)) then

                  !  JJ indicates the quotient H/C of the fuel (4=CH4;3=CH3,etc.)
                  do jj = 1,4

                    !  We look for the interval jj < chx1(icha) < jj + 1
                    if(chx1(icha).ge.4.d0) then

                    core1 = ka(4,ii) + ( (ka(4,ii+1)- ka(4,ii))/(teno(ii+1)-   &
                            teno(ii)) ) * (propce(iel,ipcte1) - teno(ii))
                    core2 = kb(4,ii) + ( (kb(4,ii+1)- kb(4,ii))/(teno(ii+1)-   &
                            teno(ii)) ) * (propce(iel,ipcte1) - teno(ii))
                    core3 = kc(4,ii) + ( (kc(4,ii+1)- kc(4,ii))/(teno(ii+1)-   &
                            teno(ii)) ) * (propce(iel,ipcte1) - teno(ii))
                    para2 = chi2(ii) + ( (chi2(ii+1)-chi2(ii))/(teno(ii+1) -   &
                            teno(ii)) ) * (propce(iel,ipcte1) - teno(ii))

                    elseif(chx1(icha).le.1.d0) then

                    core1 = ka(1,ii) + ( (ka(1,ii+1)- ka(1,ii))/               &
                            (teno(ii+1)-teno(ii)) ) * (propce(iel,ipcte1)  -   &
                            teno(ii))
                    core2 = kb(1,ii) + ( (kb(1,ii+1)- kb(1,ii))/(teno(ii+1)-   &
                            teno(ii)) ) * (propce(iel,ipcte1) - teno(ii))
                    core3 = kc(1,ii) + ( (kc(1,ii+1)- kc(1,ii))/               &
                            (teno(ii+1)-teno(ii)) ) * (propce(iel,ipcte1) -    &
                            teno(ii))
                    para2 = chi2(ii) + ( (chi2(ii+1)-chi2(ii)) /               &
                            (teno(ii+1) - teno(ii)) ) * (propce(iel,ipcte1) -  &
                            teno(ii))

                    elseif (chx1(icha).ge.jj.and.chx1(icha).lt.jj+1) then

                    core1 = ka(jj,ii) + ( (ka(jj+1,ii+1)- ka(jj,ii))/          &
                            (teno(ii+1)-teno(ii)) ) * (propce(iel,ipcte1) -    &
                            teno(ii))
                    core2 = kb(jj,ii) + ( (kb(jj+1,ii+1)- kb(jj,ii))/          &
                            (teno(ii+1)-teno(ii)) ) * (propce(iel,ipcte1) -    &
                            teno(ii))
                    core3 = kc(jj,ii) + ( (kc(jj+1,ii+1)- kc(jj,ii))/          &
                            (teno(ii+1)-teno(ii)) ) * (propce(iel,ipcte1) -    &
                            teno(ii))
                    para2 = chi2(ii) + ( (chi2(ii+1)-chi2(ii))/                &
                            (teno(ii+1) - teno(ii)) ) * (propce(iel,ipcte1) -  &
                            teno(ii))

                    endif

                  enddo

                endif

              enddo

              if(chx1(icha).ge.3.d0) then

              auxrb1 = ( volume(iel)*wmno )                                    &
                     * ( (core1 + core2 + core3) * para2 )                     &
                     * ( propce(iel,ipcyf1)*propce(iel,idgaz)/wmchx1 )         &
                     * ( rtp(iel,isca(iyno)) * crom(iel)/wmno )

              else

              auxrb1 = ( volume(iel)*wmno )                                    &
                     * ( core1 + core2 + core3 )                               &
                     * ( propce(iel,ipcyf1)*propce(iel,idgaz)/wmchx1 )         &
                     * ( rtp(iel,isca(iyno)) * crom(iel)/wmno )

              endif

              smbrs(iel)  = smbrs(iel)  - auxrb1

              propce(iel,ipproc(icnorb)) = propce(iel,ipproc(icnorb)) + auxrb1

            !  Reburning by CHx2
            elseif(propce(iel,ipcyf2).gt.0.d0) then

              !  Number of point of the temperature discretization
              do ii = 1,7

                !  We look for the interval teno(ii) < Tgaz < teno(ii+1)
                if(propce(iel,ipcte1).ge.teno(ii).and.propce(iel,ipcte1).lt.   &
                   teno(ii+1)) then

                  !  JJ incates the quotient H/C of the fuel (4=CH4;3=CH3,etc.)
                  do jj = 1,4

                    !  We look for the interval jj < chx1(icha) < jj + 1
                    if(chx2(icha).ge.4.d0) then

                    core1 = ka(4,ii) + ( (ka(4,ii+1)- ka(4,ii))/               &
                            (teno(ii+1)-teno(ii)) ) * (propce(iel,ipcte1) -    &
                            teno(ii))
                    core2 = kb(4,ii) + ( (kb(4,ii+1)- kb(4,ii))/               &
                            (teno(ii+1)-teno(ii)) ) * (propce(iel,ipcte1) -    &
                            teno(ii))
                    core3 = kc(4,ii) + ( (kc(4,ii+1)- kc(4,ii))/               &
                            (teno(ii+1)-teno(ii)) ) * (propce(iel,ipcte1) -    &
                            teno(ii))
                    para2 = chi2(ii) + ( (chi2(ii+1)-chi2(ii))/                &
                            (teno(ii+1) - teno(ii)) ) * (propce(iel,ipcte1) -  &
                            teno(ii))

                    elseif(chx2(icha).le.1.d0) then

                    core1 = ka(1,ii) + ( (ka(1,ii+1)- ka(1,ii))/(teno(ii+1)-   &
                            teno(ii)) ) * (propce(iel,ipcte1) -                &
                            teno(ii))
                    core2 = kb(1,ii) + ( (kb(1,ii+1)- kb(1,ii))/               &
                            (teno(ii+1)-teno(ii)) ) * (propce(iel,ipcte1) -    &
                            teno(ii))
                    core3 = kc(1,ii) + ( (kc(1,ii+1)- kc(1,ii))/               &
                            (teno(ii+1)-teno(ii)) ) * (propce(iel,ipcte1) -    &
                            teno(ii))
                    para2 = chi2(ii) + ( (chi2(ii+1)-chi2(ii)) /               &
                            (teno(ii+1) - teno(ii)) ) * (propce(iel,ipcte1) -  &
                            teno(ii))
                    elseif (chx2(icha).ge.jj.and.chx2(icha).lt.jj+1) then

                    core1 = ka(jj,ii) + ( (ka(jj+1,ii+1)- ka(jj,ii))/          &
                            (teno(ii+1)-teno(ii)) ) * (propce(iel,ipcte1) -    &
                            teno(ii))
                    core2 = kb(jj,ii) + ( (kb(jj+1,ii+1)- kb(jj,ii))/          &
                            (teno(ii+1)-teno(ii)) ) * (propce(iel,ipcte1) -    &
                            teno(ii))
                    core3 = kc(jj,ii) + ( (kc(jj+1,ii+1)- kc(jj,ii))/          &
                            (teno(ii+1)-teno(ii)) ) * (propce(iel,ipcte1) -    &
                            teno(ii))
                    para2 = chi2(ii) + ( (chi2(ii+1)-chi2(ii))/                &
                            (teno(ii+1) - teno(ii)) ) * (propce(iel,ipcte1) -  &
                            teno(ii))

                    endif

                  enddo

                endif

              enddo

              if(chx2(icha).ge.3.d0) then

              auxrb2 = ( volume(iel)*wmno )                                    &
                     * ( (core1 + core2 + core3) * para2 )                     &
                     * ( propce(iel,ipcyf2)*propce(iel,idgaz)/wmchx2 )         &
                     * ( rtp(iel,isca(iyno)) * crom(iel)/wmno )

              else

              auxrb2 = ( volume(iel)*wmno )                                    &
                     * ( core1 + core2 + core3 )                               &
                     * ( propce(iel,ipcyf2)*propce(iel,idgaz)/wmchx2 )         &
                     * ( rtp(iel,isca(iyno)) * crom(iel)/wmno )

              endif

              smbrs(iel)  = smbrs(iel)  - auxrb2

              propce(iel,ipproc(icnorb)) = propce(iel,ipproc(icnorb)) + auxrb2

            endif

          enddo

        endif


        !  Initialization
        do icha=1,ncharb

          gmhet (icha)=0.d0

        enddo

        do icla=1,nclacp

          icha   = ichcor(icla)
          ipctem = ipproc(itemp2(icla))
          ixckcl = isca(ixck(icla))

          mckcl1 = (1.d0-y1ch(icha))*a1ch(icha)                                &
                   *exp(-e1ch(icha)/(rr*propce(iel,ipctem)))

          mckcl2 = (1.d0-y2ch(icha))*a2ch(icha)                                &
                   *exp(-e2ch(icha)/(rr*propce(iel,ipctem)))

          !  Reaction rate of the heterogeneous combustion
          if ( rtpa(iel,ixckcl) .gt. epsicp ) then
          gmhet(icha) = gmhet(icha)                                            &
               +propce(iel,ipproc(igmhet(icla)))                               &
               *crom(iel)                                             &
               *( rtpa(iel,ixckcl)*((1.d0/(mckcl2/mckcl1+1.d0))*ynoch1(icha)   &
                   +(1.d0/(mckcl1/mckcl2+1.d0))*ynoch2(icha)) )**(2.d0/3.d0)
          endif

        enddo

        !  Coefficient of released NO during the heterogeneous combustion

        do icha=1,ncharb

          auxhet = -volume(iel)*gmhet(icha)
          propce(iel,ipproc(ifnoch)) = propce(iel,ipproc(ifnoch)) + auxhet


          smbrs(iel)  = smbrs(iel) - aux1*rtpa(iel,ivar)                       &
                                   - aux5*rtpa(iel,ivar)                       &
                                   + aux2 + aux3 + aux4 + auxhet

          rovsdt(iel) = rovsdt(iel) + aux1 + aux5

        enddo

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

!===============================================================================
! Deallocation dynamic arrays
!----
deallocate(w1,w2,w3,w4,w5,STAT=iok1)
!----
if ( iok1 > 0 ) THEN
  write(nfecra,*) ' Memory deallocation error inside: '
  write(nfecra,*) '      cs_coal_scast                '
  call csexit(1)
endif
deallocate(tfuel,STAT=iok1)
if ( iok1 > 0 ) THEN
  write(nfecra,*) ' Memory deallocation error inside: '
  write(nfecra,*) '      cs_coal_scast                '
  call csexit(1)
endif
!===============================================================================

return

end subroutine
