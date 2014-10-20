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
!> \file cs_fuel_physprop1.f90
!> \brief Calculation of physic properties of the gaseous phase
!>
!> Cell values
!> -----------
!> Temperature, mass density and average concentrations
!> (use of a PDF RECTANGLE-DIRAC)
!> ==> Fast chemistry model in 3 points
!> Extension for 3 fuels for pulverized coal
!>
!> Heterogeneous reactions
!>   - Pyrolysis
!>     Elementary composition of the mol of volatile materials
!>     The reactive coal is written C(1)H(ALPHA)O(BETA)
!>       -(k1)-> ALPHA/4 CH4  + BETA CO + (1-ALPHA/4-BETA)    Coke
!>     Reactive coal
!>       -(k2)-> ALPHA/Y CXHY + BETA CO + (1-ALPHA/RYSX-BETA) Coke
!>       With RYSX = Y/X
!>   - Heterogeneous combustion
!>     Coke + 1/2 (O2 + XSI N2) -> CO + XSI/2 N2
!>   - Reactions in gaaseous phase
!> (4/(4-RYSX)) CH4 + (O2 + XSI N2)   -(1)->  4/X/(4-RYSX)*CXHY + 2 H2O
!>                                           + XSI N2
!> CXHY + X/4*(2+RYSX) (O2 + XSI N2)  -(2)->  X CO + Y/2 H2O
!>                                           + X/4*(2+RYSX)*XSI N2
!>           CO + 1/2 (O2 + XSI N2)  -(3)->  CO2 + XSI/2 N2
!> Variables choice
!> - f1 is the mass fractions of volatile materials: CH4  + CO
!> - f2 is the mass fractions of volatile materials: CXHY + CO
!> - f3 is the mass fraction of carbone coming from the heterogeneous
!>    combustion
!>
!> Let Y  be the mass fractions and Z be the concentrations [moles/kg]
!> index f before reaction, b final
!>
!> Joint PDF degenerated into a PDF 1D of type RECTANGLE - DIRAC
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role
!______________________________________________________________________________!
!> \param[in]     ncelet        number of extended (real + ghost) cells
!> \param[in]     ncel          number of cells
!> \param[in]     nitbcp        size of the macro table \f$C_p\f$ integers
!> \param[in]     nrtbcp        size of the macro table \f$C_p\f$ reals
!> \param[in]     nitbmc        size of the macro table mc integers
!> \param[in]     nrtbmc        size of the macro table mc reals
!> \param[in]     nitbwo        size of the macro table work integers
!> \param[in]     nrtbwo        size of the macro table work reals
!> \param[in]     pa            absolute pressure in [pascals]
!> \param[in]     f1m           average of tracer 1 mvl [chx1+CO]
!> \param[in]     f2m           average of tracer 2 mvl [chx2+CO]
!> \param[in]     f3m           average of tracer 3 (CO heterogeneous comb.)
!> \param[in]     f4m           average of tracer 4 (air)
!> \param[in]     f4m           average of tracer 5 (H2O)
!> \param[in]     f4p2m         variance of tracer 4 (air)
!> \param[in]     enth          enthalpy in \f$j . kg^{-1}\f$  either gaz
!>                              or mixture
!> \param[in,out] propce        physical properties at cell centers
!______________________________________________________________________________!

subroutine cs_fuel_physprop1 &
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
use dimens, only: nvar
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
use cs_fuel_incl
use field

!===============================================================================

implicit none

! Arguments
integer          ncelet , ncel
!
double precision f1m(ncelet), f2m(ncelet) , f3m(ncelet)
double precision f4m(ncelet), f5m(ncelet) , f6m(ncelet)
double precision f7m(ncelet), f8m(ncelet) , f9m(ncelet)
double precision fvp2m(ncelet)
double precision enth(ncelet),enthox(ncelet)
double precision propce(ncelet,*)
double precision rom1(ncelet)
! Local variables
integer          iel , ii ,ice , icla
integer          ipcyf1,ipcyf2,ipcyf3,ipcyf4,ipcyf5,ipcyf6,ipcyf7
integer          ipcyox,ipcyp1,ipcyp2,ipcyp3,ipcyin
integer          ipcte1

double precision wmolme
double precision somch , cfolc , cfolh , cfolo

integer          iok1 , iok2 , iok3 , iok4 , iok5
integer          , dimension ( : )     , allocatable :: intpdf
double precision , dimension ( : )     , allocatable :: fmini,fmaxi,ffuel
double precision , dimension ( : )     , allocatable :: dfuel,doxyd,pdfm1,pdfm2,hrec
double precision , dimension ( : )     , allocatable :: cx1m,cx2m,wmf1,wmf2
double precision , dimension ( : , : ) , allocatable :: af1    , af2
double precision , dimension ( : )     , allocatable :: fs3no  , fs4no
double precision , dimension ( : , : ) , allocatable :: yfs4no
double precision, allocatable, dimension(:) :: tpdf
double precision, dimension(:), pointer :: cpro_x1
double precision, dimension(:), pointer :: cvar_yfolcl

integer          ipass
data ipass / 0 /

!===============================================================================
! 0. Memory allocation
!===============================================================================

! Massic fraction of gas
call field_get_val_s_by_name("x_c", cpro_x1)

allocate(intpdf(1:ncel)                                       ,stat=iok1)
allocate(fmini(1:ncel)      ,fmaxi(1:ncel)      ,ffuel(1:ncel),stat=iok2)
allocate(dfuel(1:ncel)      ,doxyd(1:ncel)      ,pdfm1(1:ncel),stat=iok3)
allocate(pdfm2(1:ncel)      ,hrec(1:ncel)                     ,stat=iok3)
allocate(cx1m(1:ncel)       ,cx2m(1:ncel) ,stat=iok4)
allocate(wmf1(1:ncel)       ,wmf2(1:ncel)                     ,stat=iok4)
allocate(af1(1:ncel,1:ngazg),af2(1:ncel,1:ngazg)              ,stat=iok5)
!----
if ( iok1 > 0 .or. iok2 > 0 .or. iok3 > 0 .or. iok4 > 0 .or. iok5 > 0 ) then
  write(nfecra,*) ' Memory allocation error inside: '
  write(nfecra,*) '     cs_fuel_physprop1           '
  call csexit(1)
endif

if ( ieqnox .eq. 1 ) then

  allocate(fs3no(1:ncel) , fs4no(1:ncel),stat=iok1)

  if ( iok1 > 0 ) then
    write(nfecra,*) ' Memory allocation error inside: '
    write(nfecra,*) '     cs_fuel_physprop1           '
    write(nfecra,*) ' for fs3no and fs4no arrays      '
    call csexit(1)
  endif

  allocate(yfs4no(1:ncel,1:ngazg),stat=iok2)

  if ( iok2 > 0 ) then
    write(nfecra,*) ' Memory allocation error inside: '
    write(nfecra,*) '     cs_fuel_physprop1           '
    write(nfecra,*) ' for fs4 concentrations array    '
    call csexit(1)
  endif
endif
!
!===============================================================================
! 1. Initialization
!===============================================================================

! pointer
ipcyf1 = ipproc(iym1(ifo0))
ipcyf2 = ipproc(iym1(ifov))
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
  !  Sum of F1+F2 (for the fuel, f1=0)
  ffuel(iel)=f1m(iel)+f2m(iel)
enddo

allocate(tpdf(ncelet))

call pppdfr &
!==========
 ( ncelet ,ncel   , intpdf ,                                      &
   tpdf   ,                                                       &
   ffuel  , fvp2m ,                                               &
   fmini  , fmaxi ,                                               &
   doxyd  , dfuel , pdfm1 , pdfm2 , hrec )

! Free memory
deallocate(tpdf)

!===============================================================================
! 2.Calculation of average concentrations
!===============================================================================
! - The molar masses of gas fuel are calculated in
!   cs_fuel_readata.f90.

do iel = 1, ncel
  wmf1(iel) = wmole(ifo0)
  wmf2(iel) = wmole(ifov)
enddo

do iel=1,ncel

  do ii=1,ngazg
    af1(iel,ii) = zero
    af2(iel,ii) = zero
  enddo

  af2(iel,ifov) = chfov/wmole(ifov)
  af2(iel,ico ) = cofov/wmole(ico )
  af2(iel,ih2s) = hsfov/wmole(ih2s)

  cx1m(iel) = 0.d0
  cx2m(iel) = nhcfov

enddo

call cs_gascomb &
!==============
 ( ncelet , ncel   , ifo0 , ifov ,                                &
   intpdf ,                                                       &
   f1m    , f2m , f3m , f4m , f5m , f6m , f7m , f8m , f9m ,       &
   pdfm1  , pdfm2  , doxyd    , dfuel  , hrec ,                   &
   af1    , af2    , cx1m     , cx2m   , wmf1   , wmf2 ,          &
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
! 4. Calculation of temperature and mass density
!===============================================================================

ipcte1 = ipproc(itemp1)

call cs_fuel_thfieldconv1 &
!========================
 ( ncelet , ncel   ,                                              &
   enth   ,                                                       &
   propce(1,ipcyf1), propce(1,ipcyf2), propce(1,ipcyf3),          &
   propce(1,ipcyf4), propce(1,ipcyf5), propce(1,ipcyf6),          &
   propce(1,ipcyf7), propce(1,ipcyox), propce(1,ipcyp1),          &
   propce(1,ipcyp2), propce(1,ipcyp3), propce(1,ipcyin),          &
   propce(1,ipcte1) )

do iel = 1, ncel
  wmolme = propce(iel,ipcyf1)/wmole(ifo0)                         &
         + propce(iel,ipcyf2)/wmole(ifov)                         &
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

  call cs_fuel_noxst &
 !==================
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
  enddo

endif

!===============================================================================
! 6. Calculation of balances in C, O et H
!===============================================================================

do iel=1,ncel
  propce(iel,ipproc(ibcarbone )) = cpro_x1(iel)                        &
            *( propce(iel,ipcyf1)*wmolat(iatc)/wmole(ifo0)             &
              +propce(iel,ipcyf2)*wmolat(iatc)/wmole(ifov)             &
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
            *( propce(iel,ipcyf1)*nhcfov*wmolat(iath)/wmole(ifo0)      &
              +propce(iel,ipcyf2)*nhcfov*wmolat(iath)/wmole(ifov)      &
              +propce(iel,ipcyf4)*2.d0     *wmolat(iath)/wmole(ih2s)   &
              +propce(iel,ipcyf5)*2.d0     *wmolat(iath)/wmole(ihy )   &
              +propce(iel,ipcyf6)*          wmolat(iath)/wmole(ihcn)   &
              +propce(iel,ipcyf7)*3.d0     *wmolat(iath)/wmole(inh3)   &
              +propce(iel,ipcyp2)*2.d0     *wmolat(iath)/wmole(ih2o) )

enddo

do icla = 1, nclafu
  call field_get_val_s(ivarfl(isca(iyfol(icla))), cvar_yfolcl)
  do iel = 1, ncel

    somch = cfol + hfol + ofol
    cfolc  = cfol/somch
    cfolh  = hfol/somch
    cfolo  = ofol/somch

    propce(iel,ipproc(ibcarbone )) = propce(iel,ipproc(ibcarbone ))    &
                                    +cvar_yfolcl(iel)*cfolc

    propce(iel,ipproc(iboxygen )) = propce(iel,ipproc(iboxygen ))    &
                                   +cvar_yfolcl(iel)*cfolo

    propce(iel,ipproc(ibhydrogen)) = propce(iel,ipproc(ibhydrogen))  &
                                    +cvar_yfolcl(iel)*cfolh

  enddo
enddo

!----
! Formats
!----

! Deallocation dynamic arrays
deallocate(intpdf,                      stat=iok1)
deallocate(fmini,fmaxi,ffuel,           stat=iok2)
deallocate(dfuel,doxyd,pdfm1,pdfm2,hrec,stat=iok3)
deallocate(cx1m,cx2m,wmf1,wmf2,      stat=iok4)
deallocate(af1, af2,                    stat=iok5)

if ( iok1 > 0 .or. iok2 > 0 .or. iok3 > 0 .or. iok4 > 0 .or. iok5 > 0 ) then
  write(nfecra,*) ' Memory deallocation error inside: '
  write(nfecra,*) '     cs_fuel_physprop1             '
  call csexit(1)
endif
if ( ieqnox .eq. 1 ) then

  deallocate(fs3no,fs4no,stat=iok1)

  if ( iok1 > 0 ) then
    write(nfecra,*) ' Memory deallocation error inside: '
    write(nfecra,*) '     cs_fuel_physprop1             '
    write(nfecra,*) ' for fs3no and fs4no arrays        '
    call csexit(1)
  endif

  deallocate(yfs4no,stat=iok2)

  if ( iok2 > 0 ) then
    write(nfecra,*) ' Memory deallocation error inside: '
    write(nfecra,*) '     cs_fuel_physprop1             '
    write(nfecra,*) ' for fs4 concentrations array      '
    call csexit(1)
  endif
endif

!----
! End
!----

return
end subroutine
