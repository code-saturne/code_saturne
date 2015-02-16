!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2015 EDF S.A.
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
! Function :
! --------
!> \file cs_coal_physprop1.f90
!> \brief Calculation of the physic propeties in gaseous phase

!> Cell values
!> -----------
!> Temperature, mass density and average concentrations
!> (use of a pdf rectangle-dirac)
!> ==> fast chemistry model in 3 points
!> Extension to 3 combustibles for pulverized coal
!>
!> Heterogeneous reactions
!>  - Pyrolysis
!>    Elementary composition of the mol of volatile materials
!>    The reactive coal is written C(1)H(ALPHA)O(BETA)
!>      -(k1)-> ALPHA/4 CH4  + BETA CO + (1-ALPHA/4-BETA)    Coke
!>    Reactive coal
!>      -(k2)-> ALPHA/Y CXHY + BETA CO + (1-ALPHA/RYSX-BETA) Coke
!>      With RYSX = Y/X
!>  - Heterogeneous combustion
!>    Coke + 1/2 (O2 + XSI N2) -> CO + XSI/2 N2
!>  - Gas-phase reaction
!> (4/(4-RYSX)) CH4 + (O2 + XSI N2)   -(1)->  4/X/(4-RYSX)*CXHY + 2 H2O
!>                                            + XSI N2
!> CXHY + X/4*(2+RYSX) (O2 + XSI N2)  -(2)->  X CO + Y/2 H2O
!>                                           + X/4*(2+RYSX)*XSI N2
!>            CO + 1/2 (O2 + XSI N2)  -(3)->  CO2 + XSI/2 N2
!> Variable choice
!> - F1 is the mass fractions of volatile materials: CH4  + CO
!> - F2 is the mass fractions of volatile materials: CXHY + CO
!> - F3 is the mass fraction of carbon coming from the homogeneous combustion
!> Let Y be the mass fractions and Z be the concentrations (moles/kg)
!> of index f before reaction, b final
!>
!> Joint PDF degenerated into a 1D PDF of type RECTANGLE - DIRAC
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role
!______________________________________________________________________________!
!> \param[in]     ncelet        number of extended (real + ghost) cells
!> \param[in]     ncel          number of cells
!> \param[in]     f1m           average of tracer 1 mvl [CHx1+CO]
!> \param[in]     f2m           average of tracer 2 mvl [CHx2+CO]
!> \param[in]     f3m           average of tracer 3 (CO heterogeneous comb.)
!> \param[in]     f4m           average of tracer 4 (oxyd2 mass fraction)
!> \param[in]     f5m           average of tracer 5 (oxyd3 mass fraction)
!> \param[in]     f6m           average of tracer 6 (water mass fraction)
!> \param[in]     f7m           average of tracer 7 (coal oxidyzed by O2
!>                              fraction)
!> \param[in]     f8m           average of tracer 8 (coal gasified by CO2
!>                              fraction)
!> \param[in]     f9m           average of tracer 9 (coal gasified by H2O
!>                              fraction)
!> \param[in]     fvp2m         f1f2 variance
!> \param[in]     enth          enthalpy in \f$ j . kg^{-1} \f$  either of
!>                                         the gas or of the mixture
!> \param[in]     enthox        oxydant enthalpy
!> \param[in,out] propce        physical properties at cell centers
!> \param[out]    rom1          gas density
!______________________________________________________________________________!

subroutine cs_coal_physprop1 &
 ( ncelet , ncel   ,                                      &
   f1m    , f2m    , f3m    , f4m    , f5m    ,           &
   f6m    , f7m    , f8m    , f9m    , fvp2m  ,           &
   enth   , enthox ,                                      &
   propce , rom1   )

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
use pointe
use field

!===============================================================================

implicit none

! Arguments
integer          ncelet , ncel

double precision f1m(ncelet), f2m(ncelet) , f3m(ncelet)
double precision f4m(ncelet), f5m(ncelet) , f6m(ncelet)
double precision f7m(ncelet), f8m(ncelet) , f9m(ncelet)
double precision fvp2m(ncelet)
double precision enth(ncelet),enthox(ncelet)
double precision propce(ncelet,*)
double precision rom1(ncelet)

! Local variables

integer          iel , ii , icha ,ice , icla
integer          ipcyf1,ipcyf2,ipcyf3,ipcyf4,ipcyf5,ipcyf6,ipcyf7
integer          ipcyox,ipcyp1,ipcyp2,ipcyp3,ipcyin
integer          ipcte1
integer          iok1 , iok2 , iok3 , iok4 , iok5

double precision xch    , xck    , xash , xwat
double precision zchx10 , zchx20
double precision den1   , den2 , f1mc , f2mc
double precision wmolme
double precision somch , somck , chxc , chxh , chxo , ckxc , ckxh , ckxo

integer          , dimension ( : )     , allocatable :: intpdf
double precision , dimension ( : )     , allocatable :: fmini,fmaxi,ffuel
double precision , dimension ( : )     , allocatable :: dfuel,doxyd,pdfm1,pdfm2,hrec
double precision , dimension ( : )     , allocatable :: cx1m,cx2m,wmchx1,wmchx2
double precision , dimension ( : , : ) , allocatable :: af1    , af2
double precision , dimension ( : )     , allocatable :: fs3no  , fs4no
double precision , dimension ( : , : ) , allocatable :: yfs4no
double precision, allocatable, dimension(:) :: tpdf
double precision, dimension(:), pointer :: cpro_x1
double precision, dimension(:), pointer :: cvar_xchcl, cvar_xckcl, cvar_xnpcl
double precision, dimension(:), pointer :: cvar_xwtcl
type(pmapper_double_r1), dimension(:), allocatable :: cvar_f1m, cvar_f2m

integer          ipass
data ipass / 0 /

!===============================================================================
! 0. Memory allocation
!===============================================================================

! Massic fraction of gas
call field_get_val_s_by_name("x_c", cpro_x1)

allocate(intpdf(1:ncel)                                 ,stat=iok1)
allocate(fmini(1:ncel)      ,fmaxi(1:ncel),ffuel(1:ncel),stat=iok2)
allocate(dfuel(1:ncel)      ,doxyd(1:ncel),pdfm1(1:ncel),stat=iok3)
allocate(pdfm2(1:ncel)      ,hrec(1:ncel)               ,stat=iok3)
allocate(cx1m(1:ncel) ,cx2m(1:ncel) ,stat=iok4)
allocate(wmchx1(1:ncel)     ,wmchx2(1:ncel)             ,stat=iok4)
allocate(af1(1:ncel,1:ngazg),af2(1:ncel,1:ngazg)        ,stat=iok5)
!----
if ( iok1 > 0 .or. iok2 > 0 .or. iok3 > 0 .or. iok4 > 0 .or. iok5 > 0 ) then
  write(nfecra,*) ' Memory allocation error inside : '
  write(nfecra,*) '     cs_coal_physprop1            '
  call csexit(1)
endif

if ( ieqnox .eq. 1 ) then
  allocate(fs3no(1:ncel) , fs4no(1:ncel),stat=iok1)
  if ( iok1 > 0 ) then
    write(nfecra,*) ' Memory allocation error inside: '
    write(nfecra,*) '     cs_coal_physprop1           '
    write(nfecra,*) ' for fs3no and fs4no arrays      '
    call csexit(1)
  endif
  allocate(yfs4no(1:ncel,1:ngazg),stat=iok2)
  if ( iok2 > 0 ) then
    write(nfecra,*) ' Memory allocation error inside: '
    write(nfecra,*) '     cs_coal_physprop1           '
    write(nfecra,*) ' for fs4 concentrations array    '
    call csexit(1)
  endif
endif
!
!===============================================================================
! 1. Initialization
!===============================================================================

! Arrays of pointers containing the fields values for each class
! (loop on cells outside loop on classes)
allocate(cvar_f1m(ncharb), cvar_f2m(ncharb))
do icha = 1, ncharb
  call field_get_val_s(ivarfl(isca(if1m(icha))), cvar_f1m(icha)%p)
  call field_get_val_s(ivarfl(isca(if2m(icha))), cvar_f2m(icha)%p)
enddo

! pointer
ipcyf1 = ipproc(iym1(ichx1))
ipcyf2 = ipproc(iym1(ichx2))
ipcyf3 = ipproc(iym1(ico  ))
ipcyf4 = ipproc(iym1(ih2s ))
ipcyf5 = ipproc(iym1(ihy  ))
ipcyf6 = ipproc(iym1(ihcn ))
ipcyf7 = ipproc(iym1(inh3 ))
ipcyox = ipproc(iym1(io2  ))
ipcyp1 = ipproc(iym1(ico2 ))
ipcyp2 = ipproc(iym1(ih2o ))
ipcyp3 = ipproc(iym1(iso2 ))
ipcyin = ipproc(iym1(in2  ))

ipass = ipass + 1

!===============================================================================
! 2. Determining the type of pdf
!===============================================================================

do iel = 1, ncel
  !  min et max boundary of the pdf
  fmini(iel) = 0
  fmaxi(iel) = 1.d0
  ! Sum F1+F2
  ffuel(iel)=f1m(iel)+f2m(iel)
enddo

allocate(tpdf(ncelet))

call pppdfr &
!==========
 ( ncelet , ncel  , intpdf ,                                      &
   tpdf   ,                                                       &
   ffuel  , fvp2m ,                                               &
   fmini  , fmaxi ,                                               &
   doxyd  , dfuel , pdfm1 , pdfm2 , hrec )

! Free memory
deallocate(tpdf)

!===============================================================================
! 2.Calculation of average concentrations
!===============================================================================

! Calculation of F8MC = sum(F8M(icha))
! Calculation of F9MC = sum(F9M(icha))
do iel = 1, ncel

  do ii=1,ngazg
    af1(iel,ii) = 0.d0
    af2(iel,ii) = 0.d0
  enddo

enddo

do iel = 1, ncel

  zchx10 = zero
  zchx20 = zero
  cx1m(iel) = zero
  cx2m(iel) = zero

  do icha = 1, ncharb

    f1mc = cvar_f1m(icha)%p(iel) / cpro_x1(iel)
    f2mc = cvar_f2m(icha)%p(iel) / cpro_x1(iel)

    den1 = 1.d0                                                    &
          / ( a1(icha)*wmole(ichx1c(icha))+b1(icha)*wmole(ico)     &
             +c1(icha)*wmole(ih2o)        +d1(icha)*wmole(ih2s)    &
             +e1(icha)*wmole(ihcn)        +f1(icha)*wmole(inh3) )
    cx1m(iel) = cx1m(iel)                                          &
               +den1*f1mc*a1(icha)*wmole(ichx1c(icha))
    af1(iel,ichx1)= af1(iel,ichx1)                                 &
                   +den1*f1mc*a1(icha)
    af1(iel,ico)  = af1(iel,ico)                                   &
                   +den1*f1mc*b1(icha)
    af1(iel,ih2o) = af1(iel,ih2o)                                  &
                   +den1*f1mc*c1(icha)
    af1(iel,ih2s) = af1(iel,ih2s)                                  &
                   +den1*f1mc*d1(icha)
    af1(iel,ihcn) = af1(iel,ihcn)                                  &
                   +den1*f1mc*e1(icha)
    af1(iel,inh3) = af1(iel,inh3)                                  &
                   +den1*f1mc*f1(icha)

    den2 = 1.d0                                                    &
          / ( a2(icha)*wmole(ichx2c(icha))+b2(icha)*wmole(ico)     &
             +c2(icha)*wmole(ih2o)        +d2(icha)*wmole(ih2s)    &
             +e2(icha)*wmole(ihcn)        +f2(icha)*wmole(inh3) )
    cx2m(iel) = cx2m(iel)                                          &
               +den2*f2mc*a2(icha)*wmole(ichx2c(icha))
    af2(iel,ichx2)=  af2(iel,ichx2)                                &
                    +den2*f2mc*a2(icha)
    af2(iel,ico)  = af2(iel,ico)                                   &
                   +den2*f2mc*b2(icha)
    af2(iel,ih2o) = af2(iel,ih2o)                                  &
                   +den2*f2mc*c2(icha)
    af2(iel,ih2s) = af2(iel,ih2s)                                  &
                   +den2*f2mc*d2(icha)
    af2(iel,ihcn) = af2(iel,ihcn)                                  &
                   +den2*f2mc*e2(icha)
    af2(iel,inh3) = af2(iel,inh3)                                  &
                   +den2*f2mc*f2(icha)

  enddo
  if ( af1(iel,ichx1).gt.epzero ) then
    cx1m(iel) = ( cx1m(iel)/af1(iel,ichx1)-wmolat(iatc) ) / wmolat(iath)
  else
    cx1m(iel) = 4.D0
  endif
  if ( af2(iel,ichx2).gt.epzero )  then
    cx2m(iel) = ( cx2m(iel)/af2(iel,ichx2)-wmolat(iatc) ) / wmolat(iath)
  else
    cx2m(iel) = 2.D0
  endif
  if ( f1m(iel) .gt. zero ) then
    do ii=1,ngazg
      af1(iel,ii)= af1(iel,ii)/f1m(iel)
    enddo
  else
    do ii=1,ngazg
      af1(iel,ii)= 0.d0
    enddo
  endif
  if ( f2m(iel) .gt. zero ) then
    do ii=1,ngazg
      af2(iel,ii)= af2(iel,ii)/f2m(iel)
    enddo
  else
    do ii=1,ngazg
      af2(iel,ii)= 0.d0
    enddo
  endif

  wmchx1(iel) = wmolat(iatc) + cx1m(iel)*wmolat(iath)
  wmchx2(iel) = wmolat(iatc) + cx2m(iel)*wmolat(iath)

enddo

deallocate(cvar_f1m, cvar_f2m)

call cs_gascomb &
!==============
 ( ncelet , ncel   , ichx1 , ichx2 ,                              &
   intpdf ,                                                       &
   f1m    , f2m , f3m , f4m , f5m , f6m , f7m , f8m , f9m ,       &
   pdfm1  , pdfm2  , doxyd    , dfuel  , hrec ,                   &
   af1    , af2    , cx1m     , cx2m   , wmchx1   , wmchx2 ,      &
   propce(1,ipcyf1) , propce(1,ipcyf2) , propce(1,ipcyf3) ,       &
   propce(1,ipcyf4) , propce(1,ipcyf5) , propce(1,ipcyf6) ,       &
   propce(1,ipcyf7) , propce(1,ipcyox) , propce(1,ipcyp1) ,       &
   propce(1,ipcyp2) , propce(1,ipcyp3) , propce(1,ipcyin) ,       &
   fs3no , fs4no , yfs4no )

! --> Eventual clipping of mass fractions

do iel = 1, ncel
  do ice = 1, ngazg

    if ( propce(iel,ipproc(iym1(ice))) .lt. zero )  then
       propce(iel,ipproc(iym1(ice))) = zero
    endif

    if ( propce(iel,ipproc(iym1(ice))) .gt. 1.d0 )  then
       propce(iel,ipproc(iym1(ice))) = 1.d0
    endif

    if ( abs(propce(iel,ipproc(iym1(ice)))) .lt. epsicp )  then
       propce(iel,ipproc(iym1(ice))) = zero
    endif
  enddo
enddo

!===============================================================================
! 4. Calculation of temperature and density
!===============================================================================

ipcte1 = ipproc(itemp1)

call cs_coal_thfieldconv1 &
!========================
 ( ncelet , ncel   ,                                              &
   enth   ,                                                       &
   propce(1,ipcyf1), propce(1,ipcyf2), propce(1,ipcyf3),          &
   propce(1,ipcyf4), propce(1,ipcyf5), propce(1,ipcyf6),          &
   propce(1,ipcyf7), propce(1,ipcyox), propce(1,ipcyp1),          &
   propce(1,ipcyp2), propce(1,ipcyp3), propce(1,ipcyin),          &
   propce(1,ipcte1) )

do iel = 1, ncel
  wmolme = propce(iel,ipcyf1)/wmchx1(iel)                         &
         + propce(iel,ipcyf2)/wmchx2(iel)                         &
         + propce(iel,ipcyf3)/wmole(ico )                         &
         + propce(iel,ipcyf4)/wmole(ih2s)                         &
         + propce(iel,ipcyf5)/wmole(ihy )                         &
         + propce(iel,ipcyf6)/wmole(ihcn)                         &
         + propce(iel,ipcyf7)/wmole(inh3)                         &
         + propce(iel,ipcyox)/wmole(io2 )                         &
         + propce(iel,ipcyp1)/wmole(ico2)                         &
         + propce(iel,ipcyp2)/wmole(ih2o)                         &
         + propce(iel,ipcyp3)/wmole(iso2)                         &
         + propce(iel,ipcyin)/wmole(in2 )

  ! storage of the mixture molar mass

  propce(iel,ipproc(immel)) = 1.d0 / wmolme

  ! ---- We do not include the mecanical pressure IPR
  !      but P0

  rom1(iel) = p0 / (wmolme*rr*propce(iel,ipcte1))
enddo

!===============================================================================
! 5. Nox's model
!===============================================================================


! Nox's model: Not used at the first relative iteration

if ( ieqnox .eq. 1 .and. ipass .gt. 1 ) then

  call cs_coal_noxst &

 ( ncelet , ncel   ,                                              &
   intpdf ,                                                       &
   pdfm1  , pdfm2  , doxyd  , dfuel  , hrec ,                     &
   f3m    , f4m    , f5m    , f6m    , f7m  , f8m , f9m ,         &
   fs3no  , fs4no  , yfs4no , enthox ,                            &
   propce  )

else if ( ieqnox .eq. 1 ) then

  do iel = 1, ncel
    propce(iel,ipproc(ighcn1)) = 0.d0
    propce(iel,ipproc(ighcn2)) = 0.d0
    propce(iel,ipproc(ignoth)) = 0.d0
    propce(iel,ipproc(ignh31)) = 0.d0
    propce(iel,ipproc(ignh32)) = 0.d0
    propce(iel,ipproc(ifhcnr)) = 0.d0
    propce(iel,ipproc(ifhcnd)) = 0.d0
    propce(iel,ipproc(ifhcnc)) = 0.d0
    propce(iel,ipproc(ifnh3d)) = 0.d0
    propce(iel,ipproc(ifnh3c)) = 0.d0
    propce(iel,ipproc(ifnohc)) = 0.d0
    propce(iel,ipproc(ifnonh)) = 0.d0
    propce(iel,ipproc(ifnoch)) = 0.d0
    propce(iel,ipproc(ifnoth)) = 0.d0
    propce(iel,ipproc(icnohc)) = 0.d0
    propce(iel,ipproc(icnonh)) = 0.d0
    propce(iel,ipproc(ifhcnr)) = 0.d0
    propce(iel,ipproc(icnorb)) = 0.d0
    propce(iel,ipproc(igrb)) = 0.d0
  enddo

endif

!===============================================================================
! 6. Calculation of balances in C, O et H
!===============================================================================

do iel=1,ncel
  propce(iel,ipproc(ibcarbone )) = cpro_x1(iel)                        &
            *( propce(iel,ipcyf1)*wmolat(iatc)/wmchx1(iel)             &
              +propce(iel,ipcyf2)*wmolat(iatc)/wmchx2(iel)             &
              +propce(iel,ipcyf3)*wmolat(iatc)/wmole(ico )             &
              +propce(iel,ipcyf6)*wmolat(iatc)/wmole(ihcn)             &
              +propce(iel,ipcyp1)*wmolat(iatc)/wmole(ico2) )

  propce(iel,ipproc(iboxygen  )) = cpro_x1(iel)                        &
            *( propce(iel,ipcyf3)*     wmolat(iato)/wmole(ico )        &
              +propce(iel,ipcyox)*2.d0*wmolat(iato)/wmole(io2 )        &
              +propce(iel,ipcyp1)*2.d0*wmolat(iato)/wmole(ico2)        &
              +propce(iel,ipcyp2)*     wmolat(iato)/wmole(ih2o)        &
              +propce(iel,ipcyp3)*2.d0*wmolat(iato)/wmole(iso2) )

  propce(iel,ipproc(ibhydrogen)) = cpro_x1(iel)                        &
            *( propce(iel,ipcyf1)*cx1m(iel)*wmolat(iath)/wmchx1(iel)   &
              +propce(iel,ipcyf2)*cx2m(iel)*wmolat(iath)/wmchx2(iel)   &
              +propce(iel,ipcyf4)*2.d0     *wmolat(iath)/wmole(ih2s)   &
              +propce(iel,ipcyf5)*2.d0     *wmolat(iath)/wmole(ihy )   &
              +propce(iel,ipcyf6)*          wmolat(iath)/wmole(ihcn)   &
              +propce(iel,ipcyf7)*3.d0     *wmolat(iath)/wmole(inh3)   &
              +propce(iel,ipcyp2)*2.d0     *wmolat(iath)/wmole(ih2o) )

enddo

do icla = 1, nclacp
  icha = ichcor(icla)
  call field_get_val_s(ivarfl(isca(ixch(icla))), cvar_xchcl)
  call field_get_val_s(ivarfl(isca(ixck(icla))), cvar_xckcl)
  call field_get_val_s(ivarfl(isca(inp(icla))), cvar_xnpcl)
  if ( ippmod(iccoal) .eq. 1 ) then
    call field_get_val_s(ivarfl(isca(ixwt(icla))), cvar_xwtcl)
  endif
  do iel = 1, ncel
    xch  = cvar_xchcl(iel)
    xck  = cvar_xckcl(iel)
    xash = cvar_xnpcl(iel)*xmash(icla)
    if ( ippmod(iccoal) .eq. 1 ) then
      xwat = cvar_xwtcl(iel)
    else
      xwat = 0.d0
    endif

    somch = cch(icha)+hch(icha)+och(icha)+sch(icha)+nch(icha)
    chxc  = cch(icha)/somch
    chxh  = hch(icha)/somch
    chxo  = och(icha)/somch

    somck = cck(icha)+hck(icha)+ock(icha)+sck(icha)+nck(icha)
    if ( somck .gt. 0.d0 ) then
      ckxc  = cck(icha)/somck
      ckxh  = hck(icha)/somck
      ckxo  = ock(icha)/somck
    else
      ckxc  = 1.d0
      ckxh  = 0.d0
      ckxo  = 0.d0
    endif

    propce(iel,ipproc(ibcarbone )) = propce(iel,ipproc(ibcarbone ))  &
            +xch*chxc                                                &
            +xck*ckxc

    propce(iel,ipproc(iboxygen )) = propce(iel,ipproc(iboxygen ))    &
            +xch*chxo                                                &
            +xck*ckxo                                                &
            +xwat            *wmolat(iato)/wmole(ih2o)

    propce(iel,ipproc(ibhydrogen)) = propce(iel,ipproc(ibhydrogen))  &
            +xch*chxh                                                &
            +xck*ckxh                                                &
            +xwat*2.d0       *wmolat(iath)/wmole(ih2o)

  enddo
enddo

!--------
! Formats
!--------

! Deallocation dynamic arrays
deallocate(intpdf,stat=iok1)
deallocate(fmini,fmaxi,ffuel,stat=iok2)
deallocate(dfuel,doxyd,pdfm1,pdfm2,hrec,stat=iok3)
deallocate(cx1m,cx2m,wmchx1,wmchx2,stat=iok4)
deallocate(af1, af2,stat=iok5)

if ( iok1 > 0 .or. iok2 > 0 .or. iok3 > 0 .or. iok4 > 0 .or. iok5 > 0 ) then
  write(nfecra,*) ' Memory deallocation error inside: '
  write(nfecra,*) '     cs_coal_physprop1             '
  call csexit(1)
endif
if ( ieqnox .eq. 1 ) then

  deallocate(fs3no,fs4no,stat=iok1)

  if ( iok1 > 0 ) then
    write(nfecra,*) ' Memory deallocation error inside: '
    write(nfecra,*) '     cs_coal_physprop1             '
    write(nfecra,*) ' for fs3no and fs4no arrays        '
    call csexit(1)
  endif
  deallocate(yfs4no,stat=iok2)
  if ( iok2 > 0 ) then
    write(nfecra,*) ' Memory deallocation error inside: '
    write(nfecra,*) '     cs_coal_physprop1             '
    write(nfecra,*) ' For fs4 concentrations array      '
    call csexit(1)
  endif
endif

!----
! End
!----

return
end subroutine
