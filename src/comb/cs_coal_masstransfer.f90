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
! --------
!> \file cs_coal_masstransfert.f90
!>
!> \brief Calculation of the terms of mass transfert
!>        between the continous phase and the dispersed phase
!>

!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!   mode           name          role
!______________________________________________________________________________!
!> \param[in]      ncelet        number of extended (real + ghost) cells
!> \param[in]      ncel          number of cells
!> \param[in]      rtpa          calculation variables at the centre of cells
!>                                 (previous instant)
!> \param[in,out]  propce        physic properties at cell centers
!> \param[in]      volume        cell volumes
!______________________________________________________________________________!

subroutine cs_coal_masstransfer &
 ( ncelet , ncel , rtpa , propce , volume )

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use optcal
use dimens, only: nvar
use cstphy
use cstnum
use entsor
use parall
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use ppcpfu
use cs_coal_incl
use field
!===============================================================================

implicit none

! Arguments

integer          ncelet , ncel

double precision rtpa(ncelet,nflown:nvar), propce(ncelet,*)
double precision volume(ncelet)

! Local variables

integer          iel    , icha   , icla
integer          ipcte1 , ipctem , ipcro2 , ipcdia
integer          ipcgd1 , ipcgd2 , ipcgch , ipcght , ipcyox
integer          ipcsec
integer          ipcvsl , iromf , ipcte2
integer          npoin1,npoin2,npoin3,npoin4,npoin63,npoint
integer          npyv, modntl

double precision xx2     , xch    , xck    , xash   , xnp , xuash
double precision pparo2 , xdfchi , xdfext , xdftot0 , xdftot1
double precision devto1(ncharm)  , devto2(ncharm) , coxck ,  den
double precision diacka
double precision dp , lv, yvs, yv , tebl , shrd , xmeau, xmgaz
double precision xnuss, tlimit , tmini
double precision tmin, tmax, yvmin, yvmax, yymax

integer          ipyco2,ipyh2o
double precision pprco2,pprh2o

integer          iok1
double precision, dimension (:), allocatable :: x2, x2srho2, rho1, w1
double precision, dimension(:), pointer ::  crom
double precision, dimension(:), pointer :: cpro_cp

!===============================================================================
! 1. Initialization and preliminary computations
!===============================================================================

!===============================================================================
! --- Memory allocation
allocate(x2(1:ncel),x2srho2(1:ncel),rho1(1:ncel),w1(1:ncel),STAT=iok1)
if ( iok1 > 0 ) then
  write(nfecra,*) ' Memory allocation error inside : '
  write(nfecra,*) '     cs_coal_masstransfer         '
  call csexit(1)
endif
!===============================================================================

! --- Initialization of mass transfert terms

do icla = 1, nclacp
  ipcgd1 = ipproc(igmdv1(icla))
  ipcgd2 = ipproc(igmdv2(icla))
  ipcgch = ipproc(igmdch(icla))
  ipcght = ipproc(igmhet(icla))
  do iel = 1, ncel
    propce(iel,ipcgd1) = zero
    propce(iel,ipcgd2) = zero
    propce(iel,ipcgch) = zero
    propce(iel,ipcght) = zero
  enddo
enddo

! --- Calculation of mass density of the gas mixture

call field_get_val_s(icrom, crom)

! ---- Calculation of X2=somme(X2i) and X2SRO2 = sum(X2i/Rho2i)

x2     ( : ) = zero
x2srho2( : ) = zero
do icla = 1, nclacp
  ipcro2 = ipproc(irom2(icla))
  do iel = 1, ncel
    xck  = rtpa(iel,isca(ixck(icla)))
    xch  = rtpa(iel,isca(ixch(icla)))
    xash = rtpa(iel,isca(inp (icla)))*xmash(icla)
    xx2   = xch + xck + xash

   ! ---- Taking into account humidity

    if ( ippmod(iccoal) .eq. 1 ) then
      xx2 = xx2+rtpa(iel,isca(ixwt(icla)))
    endif

    x2    (iel) = x2(iel)     + xx2
    x2srho2(iel) = x2srho2(iel) + xx2 / propce(iel,ipcro2)
  enddo
enddo

! ---- Rho 1

do iel = 1, ncel
  rho1(iel) = (1.d0-x2(iel)) / (1.d0/crom(iel)-x2srho2(iel))
enddo



!===============================================================================
! 2. Mass transfert by devolatilization
!===============================================================================

ipcte1 = ipproc(itemp1)

do icla = 1, nclacp

  ipcgd1 = ipproc(igmdv1(icla))
  ipcgd2 = ipproc(igmdv2(icla))
  ipcgch = ipproc(igmdch(icla))
  ipctem = ipproc(itemp2(icla))

  do iel = 1, ncel

  ! --- Mass transfert due to degagement of light mass density (s-1) < 0

    propce(iel,ipcgd1) = -y1ch(ichcor(icla))*a1ch(ichcor(icla))   &
      * exp(-e1ch(ichcor(icla))/(rr*propce(iel,ipctem)))

  ! --- Mass transfert due to degagement of heavy mass density (s-1) < 0

    propce(iel,ipcgd2) = -y2ch(ichcor(icla))*a2ch(ichcor(icla))   &
      * exp(-e2ch(ichcor(icla))/(rr*propce(iel,ipctem)))

  ! --- Rate of disappearance of reactive coal (s-1) < 0

    propce(iel,ipcgch) = - a1ch(ichcor(icla))                     &
      * exp(-e1ch(ichcor(icla))/(rr*propce(iel,ipctem)))          &
                         - a2ch(ichcor(icla))                     &
      * exp(-e2ch(ichcor(icla))/(rr*propce(iel,ipctem)))

  enddo

enddo


!===============================================================================
! 3. Calculation of average RHO_COKE for each coal
!    It is assumed for the calculation of the mass density of the coke
!     that the devolatization takes place at constant volume
!===============================================================================

! --- Initialization

do icha = 1, ncharb
  devto1(icha) = zero
  devto2(icha) = zero
  rhock(icha)  = rho0ch(icha)
enddo

! --- Calculation of the integral of GMDEV1 and GMDEV2 for each coal

call field_get_val_s(icrom, crom)
do icla = 1, nclacp
  ipcgd1 = ipproc(igmdv1(icla))
  ipcgd2 = ipproc(igmdv2(icla))
  do iel = 1, ncel
    xch = rtpa(iel,isca(ixch(icla)))
    devto1(ichcor(icla)) = devto1(ichcor(icla)) -                 &
      ( propce(iel,ipcgd1)*xch*crom(iel)                 &
        *volume(iel) )
    devto2(ichcor(icla)) = devto2(ichcor(icla)) -                 &
      ( propce(iel,ipcgd2)*xch*crom(iel)                 &
        *volume(iel) )
  enddo
  if(irangp.ge.0) then
    call parsom(devto1(ichcor(icla)))
    call parsom(devto2(ichcor(icla)))
  endif
enddo

! --- Calculation of average mass density of coke

do icha = 1, ncharb
  den = y2ch(icha)*devto1(icha)+y1ch(icha)*devto2(icha)
  if ( den.gt.epsicp ) then
    rhock(icha) = rho0ch(icha) * ( 1.d0 -                         &
     ( y1ch(icha)*y2ch(icha)*(devto1(icha)+devto2(icha))/den ) )
  endif
enddo

!===============================================================================
! 4. Mass transfert by heterogeneous combustion with O2
!===============================================================================

ipcyox = ipproc(iym1(io2))
ipcte1 = ipproc(itemp1)

do icla = 1, nclacp

  ipcght = ipproc(igmhet(icla))
  ipcdia = ipproc(idiam2(icla))
  icha   = ichcor(icla)
  ipctem = ipproc(itemp2(icla))

  do iel = 1, ncel

    xnp   = rtpa(iel,isca(inp(icla)))
    xuash = xnp*(1.d0-xashch(icha))*xmp0(icla)

    ! --- Calculation of the partial pressure of oxygen [atm]
    !                                                 ---
    !       PO2 = RHO1*RR*T*YO2/MO2

    pparo2 = rho1(iel)*rr*propce(iel,ipcte1)*propce(iel,ipcyox)/wmole(io2)
    pparo2 = pparo2 / prefth

    ! --- Coefficient of chemical kinetics of formation of CO
    !       in [kg.m-2.s-1.atm(-n)]

    xdfchi = ahetch(ichcor(icla))                                 &
      * exp(-ehetch(ichcor(icla))*4185.d0                         &
             / (rr * propce(iel,ipctem)) )

    ! --- Diffusion coefficient  in kg/m2/s/[atm] : XDFEXT
    !     Global coefficient for n=0.5 in kg/m2/s : XDFTOT0
    !     Global coefficient for n=1   in Kg/m2/s : XDFTOT1

    diacka = propce(iel,ipcdia)/diam20(icla)
    if ( diacka .gt. epsicp ) then
      xdfext = 2.53d-7*((propce(iel,ipctem))**0.75d0)             &
              / propce(iel,ipcdia)*2.d0
      xdftot1 = pparo2 / ( 1.d0/xdfchi + 1.d0/xdfext )
      xdftot0 = -(xdfchi**2)/(2.d0*xdfext)+(pparo2*xdfchi**2      &
                + (xdfchi**4)/(4.d0*xdfext**2))**0.5d0
    else
      xdftot1 = xdfchi*pparo2
      xdftot0 = xdfchi*pparo2**0.5d0
    endif

    ! Remark AE: For now, we limit ourselves to this test
    !         The introduction of a blowing correlation will come later

    ! --- Calculation of coxck such as: Se = coxck * Xck**(2/3)

    coxck = zero
    if ( xuash.gt.epsicp ) then
      coxck = pi*(diam20(icla)**2)*                               &
            (rho20(icla)/(rhock(icha)*xuash))**(2.d0/3.d0)
    endif

    ! --- Calculation of propce(iel,ipcght) = - coxck*Xdftoto*PPARO2*Xnp < 0
    ! --- or  propce(iel,ipcght) = - coxck*XDFTOT1*PPARO2*Xnp < 0

    if (iochet(icha).eq.1) then
      propce(iel,ipcght) = - xdftot1*coxck*xnp
    else
      propce(iel,ipcght) = - xdftot0*coxck*xnp
    endif

  enddo

enddo

!===============================================================================
! 5. Mass transfert by heterogeneous combustion with CO2
!===============================================================================

if ( ihtco2 .eq. 1) then

  ipyco2 = ipproc(iym1(ico2))
  ipcte1 = ipproc(itemp1)

  do icla = 1, nclacp

    ipcght = ipproc(ighco2(icla))
    ipcdia = ipproc(idiam2(icla))
    icha   = ichcor(icla)
    ipctem = ipproc(itemp2(icla))

    do iel = 1, ncel

      xnp   = rtpa(iel,isca(inp(icla)))
      xuash = xnp*(1.d0-xashch(icha))*xmp0(icla)

      ! --- Calculation of partial pressure of CO2 [atm]
      !                                                 ---
      !       PCO2 = RHO1*RR*T*YCO2/MCO2

      pprco2 = rho1(iel)*rr*propce(iel,ipcte1)                    &
                           *propce(iel,ipyco2)/wmole(ico2)
      pprco2 = pprco2 / prefth

      ! --- Coefficient of chemical kinetics of formation of CO
      !       in kg.m-2.s-1.atm(-n)

      xdfchi = ahetc2(ichcor(icla))                               &
        * exp(-ehetc2(ichcor(icla))*4185.d0                       &
               / (rr * propce(iel,ipctem)) )

      ! --- Diffusion coefficient in kg/m2/s/[atm] : XDFEXT
      !     Global coefficient for n=0.5 in kg/m2/s : XDFTOT0
      !     Glabal coefficient for n=1   in [g/m2/s : XDFTOT1

      diacka = propce(iel,ipcdia)/diam20(icla)
      if ( diacka .gt. epsicp ) then
        xdfext = 2.53d-7*((propce(iel,ipctem))**0.75d0)           &
                / propce(iel,ipcdia)*2.d0
        xdftot1 = pprco2 / ( 1.d0/xdfchi + 1.d0/xdfext )
        xdftot0 = -(xdfchi**2)/(2.d0*xdfext)+(pprco2*xdfchi**2    &
                  +(xdfchi**4)/(4.d0*xdfext**2))**0.5d0
      else
        xdftot1 = xdfchi*pprco2
        xdftot0 = xdfchi*pprco2**0.5d0
      endif

      ! Remark AE: For now, we limit ourselves to this test
      !         The introduction of a blowing correlation will come later

      ! --- Calculation of coxck such as: Se = coxck * Xck**(2/3)

      coxck = zero
      if ( xuash.gt.epsicp ) then
        coxck = pi*(diam20(icla)**2)*                             &
              (rho20(icla)/(rhock(icha)*xuash))**(2.d0/3.d0)
      endif

      ! --- Calculation of propce(iel,ipcght) = - coxck*XDFTOT0*PPRCO2*Xnp < 0
      ! --- or  propce(iel,ipcght) = - coxck*XDFTOT1*PPRCO2*Xnp < 0

      if (ioetc2(icha).eq.1) then
        propce(iel,ipcght) = - xdftot1*coxck*xnp
      else
        propce(iel,ipcght) = - xdftot0*coxck*xnp
      endif

    enddo

  enddo

endif

!===============================================================================
! 6. Mass transfert by heterogeneous combustion with H2O
!===============================================================================

if ( ihth2o .eq. 1) then

  ipyh2o = ipproc(iym1(ih2o))
  ipcte1 = ipproc(itemp1)

  do icla = 1, nclacp

    ipcght = ipproc(ighh2o(icla))
    ipcdia = ipproc(idiam2(icla))
    icha   = ichcor(icla)
    ipctem = ipproc(itemp2(icla))

    do iel = 1, ncel

      xnp   = rtpa(iel,isca(inp(icla)))
      xuash = xnp*(1.d0-xashch(icha))*xmp0(icla)

      ! --- Calculation of partial pressure of H2O [atm]
      !                                                 ---
      !       PH2O = RHO1*RR*T*YH2O/MH2O

      pprh2o = rho1(iel)*rr*propce(iel,ipcte1)*propce(iel,ipyh2o)/wmole(ih2o)
      pprh2o = pprh2o/ prefth

      ! --- Coefficient of chemical kinetics of formation of CO
      !       in kg.m-2.s-1.atm(-n)

      xdfchi = ahetwt(ichcor(icla))                               &
        * exp(-ehetwt(ichcor(icla))*4185.d0                       &
               / (rr * propce(iel,ipctem)) )

      ! --- Diffusion coefficient in kg/m2/s/[atm]: XDFEXT
      !     Global coefficient for n=0.5 in kg/m2/s: XDFTOT0
      !     Global coefficient for n=1 in kg/m2/s: XDFTOT1

      diacka = propce(iel,ipcdia)/diam20(icla)
      if ( diacka .gt. epsicp ) then
        xdfext = 2.53d-7*((propce(iel,ipctem))**0.75d0)           &
                / propce(iel,ipcdia)*2.d0
        xdftot1 = pprh2o / ( 1.d0/xdfchi + 1.d0/xdfext )
        xdftot0 = -(xdfchi**2)/(2.d0*xdfext)+(pprh2o*xdfchi**2    &
                  +(xdfchi**4)/(4.d0*xdfext**2))**0.5d0
      else
        xdftot1 = xdfchi*pprh2o
        xdftot0 = xdfchi*pprh2o**0.5d0
      endif

      ! Reamrk AE: For now, we limit ourselves to this test
      !         The introduction of a blowing correlation will come later

      ! --- Calculation of coxck such as: Se = coxck * Xck**(2/3)

      coxck = zero
      if ( xuash.gt.epsicp ) then
        coxck = pi*(diam20(icla)**2)*                             &
              (rho20(icla)/(rhock(icha)*xuash))**(2.d0/3.d0)
      endif

      ! --- Calculation of propce(iel,ipcght) = - coxck*XDFTOT0*PPRH2O*Xnp < 0
      ! --- or propce(iel,ipcght) = - coxck*XDFTOT1*PPRH2O*Xnp < 0

      if (ioetwt(icha).eq.1) then
        propce(iel,ipcght) = - xdftot1*coxck*xnp
      else
        propce(iel,ipcght) = - xdftot0*coxck*xnp
      endif

    enddo

  enddo

endif

!===============================================================================
! 7. Mass transfert during the dryer phase
!===============================================================================

if ( ippmod(iccoal) .ge. 1 ) then

  !      Latente heat in J/kg
  lv     = 2.263d+6
  tebl   = 100.d0+tkelvi

  xnuss = 2.d0
  shrd  = 2.d0
  xmeau = 0.018d0

  tlimit = 302.24d0
  tmini   = tlimit*(1.d0-tlimit/(lv*xmeau))

  if (icp.gt.0) call field_get_val_s(iprpfl(icp), cpro_cp)

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

  if(ntlist.gt.0) then
    modntl = mod(ntcabs,ntlist)
  elseif(ntlist.eq.-1.and.ntcabs.eq.ntmabs) then
    modntl = 0
  else
    modntl = 1
  endif

  do icla = 1, nclacp

    npoin1 = 0
    npoin2 = 0
    npoin3 = 0
    npoin4 = 0
    npoin63 = 0

    icha   = ichcor(icla)

    iromf  = ipproc(irom1)
    ipcte1 = ipproc(itemp1)
    ipcte2 = ipproc(itemp2(icla))
    ipcsec = ipproc(igmsec(icla))
    ipcro2 = ipproc(irom2(icla))

    ! -------- Calculation of the diameter of particles in W2
    !          d20 = (A0.D0**2+(1-A0)*DCK**2)**0.5


    tmax = -1.d+20
    tmin =  1.d+20
    yvmin = 1.d+20
    yvmax =-1.d+20

    npyv  = 0
    yymax = 0.d0

    do iel = 1, ncel

      propce(iel,ipcsec) = 0.d0

      if ( rtpa(iel,isca(ixwt(icla))) .gt. epsicp  ) then

        npoin1 = npoin1 + 1

        xnp = rtpa(iel,isca(inp(icla)))

        !             Calculation of the diameter of the particules in W2
        !             d20 = (A0.D0**2+(1-A0)*DCK**2)**0.5

        dp  = ( xashch(icha)*diam20(icla)**2 +                    &
                 (1.d0-xashch(icha))*propce(iel,ipcdia)**2        &
                  )**0.5d0

        npoin2 = npoin2 + 1

        if ( propce(iel,ipcte2) .gt. tlimit )then
          xmgaz = propce(iel,ipproc(immel))
          yvs   = xmeau/xmgaz                                   &
                 *exp( lv*xmeau                                 &
                      *(1.d0/tebl-1.d0/propce(iel,ipcte2))      &
                      /rr )

        else
          xmgaz = propce(iel,ipproc(immel))
          yvs   = xmeau/xmgaz                                   &
                 *exp( lv*xmeau                                 &
                      *(1.d0/tebl-1.d0/tlimit)                  &
                      /rr )                                     &
                 *(lv*xmeau*(propce(iel,ipcte2)-tmini))         &
                 /(tlimit*tlimit)
        endif

        yv = yvs

        if ( yv .lt. 0.d0 ) then
          write(nfecra,*) yv,yvs,propce(iel,ipproc(iym1(ih2o)))
          write(nfecra,*) propce(iel,ipcte2),tmini
          call csexit(1)
        endif
        if ( yv .ge. 1.d0) then
          yymax = max(yymax,yv)
          npyv  = npyv +1
          yv    = 0.99d0
        endif

        yvmin=min(yvmin,yv)
        yvmax=max(yvmax,yv)

        if ( yv .gt. epsicp .and. yv .lt. 1.d0) then

          npoin3 = npoin3 + 1

          propce(iel,ipcsec) = pi*dp*propce(iel,iromf)                    &
                              *diftl0*shrd*xnp                            &
                   *log((1.D0-propce(iel,ipproc(iym1(ih2o))))/(1.d0-yv))

          if ( propce(iel,ipcsec) .lt. 0.d0 ) then
            propce(iel,ipcsec) = 0.d0
            npoin63 = npoin63 + 1
          endif
        else
          npoin4 = npoin4 + 1
          propce(iel,ipcsec) = 0.d0
        endif

      else

        propce(iel,ipcsec) = 0.d0

      endif

      tmax = max(propce(iel,ipcsec),tmax)
      tmin = min(propce(iel,ipcsec),tmin)

    enddo

    npoint = ncel

    if ( irangp .ge. 0 ) then

      call parmin(tmin)
      call parmax(tmax)
      call parmin(yvmin)
      call parmax(yvmax)
      call parcpt(npyv)
      call parmax(yymax)

      call parcpt(npoint)
      call parcpt(npoin1)
      call parcpt(npoin2)
      call parcpt(npoin3)
      call parcpt(npoin4)
      call parcpt(npoin63)

    endif


    if (modntl.eq.0) then
      WRITE(NFECRA,*) ' For the class = ',ICLA,npoint
      WRITE(NFECRA,*) '   Number of pts Xwat >0   =',NPOIN1
      WRITE(NFECRA,*) '               T < TEBL  =',NPOIN2
      WRITE(NFECRA,*) '               YV     >0 =',NPOIN3
      WRITE(NFECRA,*) '               YV     <0 =',NPOIN4
      WRITE(NFECRA,*) '   Annul G   T < TEBL=   ', npoin63
      WRITE(NFECRA,*) '   Min Max ST = ',TMIN,TMAX
      WRITE(NFECRA,*) '   Min Max YV = ',YVMIN,YVMAX

      WRITE(NFECRA,*) ' Clipping of YV at Max ',NPYV,YYMAX
    endif

  enddo

endif

!===============================================================================
! Deallocation dynamic arrays
!----
deallocate(x2,x2srho2,rho1,w1,STAT=iok1)
if ( iok1 > 0 ) then
  write(nfecra,*) ' Deallocation memory error in : '
  write(nfecra,*) '     cs_coal_masstransfer    '
  call csexit(1)
endif
!===============================================================================

!--------
! Formats
!--------
!----
! End
!----

return
end subroutine
