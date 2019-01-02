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
! --------
!> \file cs_fuel_physprop1.f90
!> \brief Calculation of physic properties of the gaseous phase
!>
!> Cell values
!> -----------
!> Temperature, mass density and average concentrations
!> (use of a PDF rectangle-dirac)
!> ==> fast chemistry model in 3 points
!> Extension for 3 fuels for pulverized coal
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
!> - F3 is the mass fraction of carbon coming from the heterogeneous
!>    combustion
!>
!> Let Y  be the mass fractions and Z be the concentrations [moles/kg]
!> index f before reaction, b final
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
!> \param[in]     f4m           average of tracer 4 (air)
!> \param[in]     f5m           average of tracer 5 (H2O)
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
!> \param[out]    rom1          gas density
!______________________________________________________________________________!

subroutine cs_fuel_physprop1 &
 ( ncelet , ncel   ,                                      &
   f1m    , f2m    , f3m    , f4m    , f5m    ,           &
   f6m    , f7m    , f8m    , f9m    , fvp2m  ,           &
   enth   , enthox ,                                      &
   rom1   )

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
use cs_fuel_incl
use field
use cs_c_bindings
use pointe

!===============================================================================

implicit none

! Arguments
integer          ncelet , ncel

double precision f1m(ncelet), f2m(ncelet) , f3m(ncelet)
double precision f4m(ncelet), f5m(ncelet) , f6m(ncelet)
double precision f7m(ncelet), f8m(ncelet) , f9m(ncelet)
double precision fvp2m(ncelet)
double precision enth(ncelet),enthox(ncelet)
double precision rom1(ncelet)

! Local variables

integer          iel , ii ,ice , icla
integer          iok1 , iok2 , iok3 , iok4 , iok5

double precision wmolme
double precision somch , cfolc , cfolh , cfolo

integer          , dimension ( : )     , allocatable :: intpdf
double precision , dimension ( : )     , allocatable :: fmini,fmaxi,ffuel
double precision , dimension ( : )     , allocatable :: dfuel,doxyd,pdfm1,pdfm2,hrec
double precision , dimension ( : )     , allocatable :: cx1m,cx2m,wmf1,wmf2
double precision , dimension ( : , : ) , allocatable :: af1    , af2
double precision , dimension ( : )     , allocatable :: fs3no  , fs4no
double precision , dimension ( : , : ) , allocatable :: yfs4no
double precision, allocatable, dimension(:) :: tpdf
double precision, dimension(:), pointer :: cpro_x1, cpro_mmel, cpro_temp1
double precision, dimension(:), pointer :: cvar_yfolcl
double precision, dimension(:), pointer :: cpro_cyf1, cpro_cyf2, cpro_cyf3
double precision, dimension(:), pointer :: cpro_cyf4, cpro_cyf5, cpro_cyf6
double precision, dimension(:), pointer :: cpro_cyf7, cpro_cyox, cpro_cyp1
double precision, dimension(:), pointer :: cpro_cyp2, cpro_cyp3, cpro_cyin
double precision, dimension(:), pointer :: cpro_ghcn1, cpro_ghcn2, cpro_gnoth
double precision, dimension(:), pointer :: cpro_bcarbone, cpro_boxygen, cpro_bhydrogen
type(pmapper_double_r1), dimension(:), allocatable :: cpro_ym1

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
allocate(wmf1(1:ncel)       ,wmf2(1:ncel)               ,stat=iok4)
allocate(af1(1:ncel,1:ngazg),af2(1:ncel,1:ngazg)        ,stat=iok5)
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
call field_get_val_s(iym1(ifo0 ),cpro_cyf1)
call field_get_val_s(iym1(ifov ),cpro_cyf2)
call field_get_val_s(iym1(ico  ),cpro_cyf3)
call field_get_val_s(iym1(ih2s ),cpro_cyf4)
call field_get_val_s(iym1(ihy  ),cpro_cyf5)
call field_get_val_s(iym1(ihcn ),cpro_cyf6)
call field_get_val_s(iym1(inh3 ),cpro_cyf7)
call field_get_val_s(iym1(io2  ),cpro_cyox)
call field_get_val_s(iym1(ico2 ),cpro_cyp1)
call field_get_val_s(iym1(ih2o ),cpro_cyp2)
call field_get_val_s(iym1(iso2 ),cpro_cyp3)
call field_get_val_s(iym1(in2  ),cpro_cyin)

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
! - The molar masses of gas fuel are calculated in
!   cs_fuel_readata.f90.

do iel = 1, ncel
  wmf1(iel) = wmole(ifo0)
  wmf2(iel) = wmole(ifov)
enddo

do iel = 1, ncel

  do ii = 1, ngazg
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
   cpro_cyf1 , cpro_cyf2 , cpro_cyf3 ,                            &
   cpro_cyf4 , cpro_cyf5 , cpro_cyf6 ,                            &
   cpro_cyf7 , cpro_cyox , cpro_cyp1 ,                            &
   cpro_cyp2 , cpro_cyp3 , cpro_cyin ,                            &
   fs3no , fs4no , yfs4no )

! --> Eventual clipping of mass fractions

allocate(cpro_ym1(ngazg))
do ice = 1, ngazg
  call field_get_val_s(iym1(ice), cpro_ym1(ice)%p)
enddo

do iel = 1, ncel
  do ice = 1, ngazg

    if ( cpro_ym1(ice)%p(iel) .lt. zero )  then
       cpro_ym1(ice)%p(iel) = zero
    endif

    if ( cpro_ym1(ice)%p(iel) .gt. 1.d0 )  then
       cpro_ym1(ice)%p(iel) = 1.d0
    endif

    if ( abs(cpro_ym1(ice)%p(iel)) .lt. epsicp )  then
       cpro_ym1(ice)%p(iel) = zero
    endif
  enddo
enddo

!===============================================================================
! 4. Calculation of temperature and density
!===============================================================================

call field_get_val_s(itemp1,cpro_temp1)
call field_get_val_s(immel,cpro_mmel)

call cs_fuel_thfieldconv1(MESH_LOCATION_CELLS, enth, cpro_temp1)

do iel = 1, ncel
  wmolme = cpro_cyf1(iel)/wmole(ifo0)                         &
         + cpro_cyf2(iel)/wmole(ifov)                         &
         + cpro_cyf3(iel)/wmole(ico )                         &
         + cpro_cyf4(iel)/wmole(ih2s)                         &
         + cpro_cyf5(iel)/wmole(ihy )                         &
         + cpro_cyf6(iel)/wmole(ihcn)                         &
         + cpro_cyf7(iel)/wmole(inh3)                         &
         + cpro_cyox(iel)/wmole(io2 )                         &
         + cpro_cyp1(iel)/wmole(ico2)                         &
         + cpro_cyp2(iel)/wmole(ih2o)                         &
         + cpro_cyp3(iel)/wmole(iso2)                         &
         + cpro_cyin(iel)/wmole(in2 )

  ! storage of the mixture molar mass

  cpro_mmel(iel) = 1.d0 / wmolme

  ! ---- We do not include the mecanical pressure IPR
  !      but P0

  rom1(iel) = p0 / (wmolme*cs_physical_constants_r*cpro_temp1(iel))
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
   fs3no  , fs4no  , yfs4no , enthox )

else if ( ieqnox .eq. 1 ) then

  call field_get_val_s(ighcn1,cpro_ghcn1)
  call field_get_val_s(ighcn2,cpro_ghcn2)
  call field_get_val_s(ignoth,cpro_gnoth)

  do iel = 1, ncel
    cpro_ghcn1(iel) = 0.d0
    cpro_ghcn2(iel) = 0.d0
    cpro_gnoth(iel) = 0.d0
  enddo

endif

!===============================================================================
! 6. Calculation of balances in C, O et H
!===============================================================================

call field_get_val_s(ibcarbone,cpro_bcarbone)
call field_get_val_s(iboxygen,cpro_boxygen)
call field_get_val_s(ibhydrogen,cpro_bhydrogen)

do iel=1,ncel
  cpro_bcarbone(iel) = cpro_x1(iel)                                &
            *( cpro_cyf1(iel)*wmolat(iatc)/wmole(ifo0)             &
              +cpro_cyf2(iel)*wmolat(iatc)/wmole(ifov)             &
              +cpro_cyf3(iel)*wmolat(iatc)/wmole(ico )             &
              +cpro_cyf6(iel)*wmolat(iatc)/wmole(ihcn)             &
              +cpro_cyp1(iel)*wmolat(iatc)/wmole(ico2) )

  cpro_boxygen(iel) = cpro_x1(iel)                                 &
            *( cpro_cyf3(iel)*     wmolat(iato)/wmole(ico )        &
              +cpro_cyox(iel)*2.d0*wmolat(iato)/wmole(io2 )        &
              +cpro_cyp1(iel)*2.d0*wmolat(iato)/wmole(ico2)        &
              +cpro_cyp2(iel)*     wmolat(iato)/wmole(ih2o)        &
              +cpro_cyp3(iel)*2.d0*wmolat(iato)/wmole(iso2) )

  cpro_bhydrogen(iel) = cpro_x1(iel)                               &
            *( cpro_cyf1(iel)*nhcfov*wmolat(iath)/wmole(ifo0)      &
              +cpro_cyf2(iel)*nhcfov*wmolat(iath)/wmole(ifov)      &
              +cpro_cyf4(iel)*2.d0     *wmolat(iath)/wmole(ih2s)   &
              +cpro_cyf5(iel)*2.d0     *wmolat(iath)/wmole(ihy )   &
              +cpro_cyf6(iel)*          wmolat(iath)/wmole(ihcn)   &
              +cpro_cyf7(iel)*3.d0     *wmolat(iath)/wmole(inh3)   &
              +cpro_cyp2(iel)*2.d0     *wmolat(iath)/wmole(ih2o) )

enddo

do icla = 1, nclafu
  call field_get_val_s(ivarfl(isca(iyfol(icla))), cvar_yfolcl)
  do iel = 1, ncel

    somch = cfol + hfol + ofol
    cfolc  = cfol/somch
    cfolh  = hfol/somch
    cfolo  = ofol/somch

    cpro_bcarbone(iel) = cpro_bcarbone(iel)    &
                                    +cvar_yfolcl(iel)*cfolc

    cpro_boxygen(iel) = cpro_boxygen(iel)      &
                                   +cvar_yfolcl(iel)*cfolo

    cpro_bhydrogen(iel) = cpro_bhydrogen(iel)  &
                                    +cvar_yfolcl(iel)*cfolh

  enddo
enddo

!--------
! Formats
!--------

! Deallocation dynamic arrays
deallocate(intpdf,stat=iok1)
deallocate(fmini,fmaxi,ffuel,stat=iok2)
deallocate(dfuel,doxyd,pdfm1,pdfm2,hrec,stat=iok3)
deallocate(cx1m,cx2m,wmf1,wmf2,stat=iok4)
deallocate(af1, af2,stat=iok5)

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
