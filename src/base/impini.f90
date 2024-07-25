!-------------------------------------------------------------------------------

! This file is part of code_saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2024 EDF S.A.
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

subroutine impini () &
  bind(C, name='cs_f_impini')
!================


!===============================================================================
! Purpose:
!  ---------

! Print computation parameters after user changes in cs_user_parameters.f90

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
!__________________!____!_____!________________________________________________!

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use cstnum
use dimens
use numvar
use optcal
use cstphy
use entsor
use albase
use parall
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use radiat
use lagran
use mesh
use field
use cavitation
use cs_c_bindings
use vof

!===============================================================================

implicit none

! Arguments


! Local variables

character        chaine*80
integer          ii
integer          kscmin, kscmax, keyvar
integer          f_id, n_fields
integer          igg, ige
integer          kturt, turb_flux_model

character(len=3), dimension(3) :: nomext3
character(len=4), dimension(3) :: nomext63

!===============================================================================
! 1. Introduction
!===============================================================================

nomext3 = (/'[X]', '[Y]', '[Z]'/)
nomext63 = (/'[11]', '[22]', '[33]'/)

call field_get_n_fields(n_fields)

! Key ids for clipping
call field_get_key_id("min_scalar_clipping", kscmin)
call field_get_key_id("max_scalar_clipping", kscmax)

call field_get_key_id("variable_id", keyvar)

write(nfecra,1000)

if (ippmod(icod3p).ne.-1) then
  write(nfecra,1010)
  write(nfecra,1020) ippmod(icod3p)
  write(nfecra,1070) namgas, pcigas
  write(nfecra,1080) trim(nomcog(igfuel(1))), &
  -nreact(igoxy(1)), trim(nomcog(igoxy(1))), trim(nomcog(igprod(1)))
  write(nfecra,'(a20,10(1x,a14))') "Mass composition", "Fuel", "Oxydizer", "Products"
  write(nfecra,'(a20,10(1x,a14))') "----------------", "----", "--------", "--------"
  do ige = 1, ngaze
    write(nfecra,'(a15,10(1x,f14.5))') trim(nomcoe(ige)), (coefeg(ige,igg), igg=1, ngazg)
  enddo
  write(nfecra,1000)
  write(nfecra,'(a20,10(1x,a14))') "Molar composition", "Fuel", "Oxydizer", "Products"
  write(nfecra,'(a20,10(1x,a14))') "-----------------", "----", "--------", "--------"
  do ige = 1, ngaze
    write(nfecra,'(a15,10(1x,f14.5))') trim(nomcoe(ige)), (compog(ige,igg), igg=1, ngazg)
  enddo
else if (ippmod(islfm).ne.-1) then
  write(nfecra,1010)
  write(nfecra,1040) ippmod(islfm)
else if (ippmod(icoebu).ne.-1) then
  write(nfecra,1010)
  write(nfecra,1030) ippmod(icoebu), cebu
endif


 1000 format(                                                     &
                                                                /,&
' ===========================================================', /,&
                                                                /,&
'               CALCULATION PARAMETERS SUMMARY',                /,&
'               ==============================',                /,&
                                                                /,&
' -----------------------------------------------------------', /)

 9900 format(                                                     &
                                                                /,&
' -----------------------------------------------------------', /)
 1010 format(                                                     &
                                                                /,&
' ** SPECIFIC PHYSICS:',                                        /,&
'    ----------------',                                         /)
 1020 format(                                                     &
' --- Diffusion Flame: 3 Point Chemistry',                      /,&
'       OPTION = ',4x,i10                                       /)
 1030 format(                                                     &
' --- Premixed Flame: EBU Model',                               /,&
'       OPTION = ',4x,i10,                                      /,&
'       CEBU   = ',e14.5                                        /)
 1040 format(                                                     &
' --- Diffusion Flame: Steady laminar flamelet model',          /,&
'       OPTION = ',4x,i10,                                      /)
 1070 format(                                                     &
' --- Combustible characteristics',                             /,&
'       Combustible : ',4x,a,                                   /,&
'       PCI = ',4x,e14.5,  ' J/kg',                             /)
 1080 format(                                                     &
" --- Chemical reaction: ",                                     /,&
"       ", a," + ",f6.3," (",a,") --> ",a,                      /)

!===============================================================================
! 2. DEFINITION GENERALE DU CAS
!===============================================================================

! --- Dimensions

write(nfecra,1500)
write(nfecra,1520) nvar,nscal,nscaus,nscapp

write(nfecra,9900)

 1500 format(                                                     &
                                                                /,&
' ** DIMENSIONS',                                               /,&
'    ----------',                                               /)
 1520 format(                                                     &
' --- Physics',                                                 /,&
'       NVAR   = ',4x,i10,    ' (Nb variables                )',/,&
'       NSCAL  = ',4x,i10,    ' (Nb scalars                  )',/,&
'       NSCAUS = ',4x,i10,    ' (Nb user scalars             )',/,&
'       NSCAPP = ',4x,i10,    ' (Nb specific physics scalars )',/)

!===============================================================================
! 3. MODELISATION PHYSIQUE
!===============================================================================

! --- Homogeneous VoF Model

write(nfecra,2101)
write(nfecra,2111) ivofmt

if (ivofmt.gt.0) then

  write(nfecra,2121) rho1, mu1
  write(nfecra,2131) rho2, mu2
  write(nfecra,2141) sigmaS
  write(nfecra,2151) idrift, kdrift, cdrift

  ! --- cavitation model

  write(nfecra,2100)

  if (iand(ivofmt,VOF_MERKLE_MASS_TRANSFER).ne.0) then

    write(nfecra,2120)
    write(nfecra,2130)

    write(nfecra,2140) presat, linf, uinf
    if (itytur.eq.2 .or. itytur.eq.5 .or. iturb.eq.60 .or. iturb.eq.70) then
      write(nfecra,2150) icvevm, mcav
    endif

  endif

endif

write(nfecra,9900)


 2100 format(                                                     &
                                                                /,&
' ** CAVITATION MODEL'                        ,                 /,&
'    ----------------------------------------',                 /)
 2120 format(                                                     &
'  -- Liquid phase: fluid 1',                                   /)
 2130 format(                                                     &
'  -- Gas phase: fluid 2',                                      /)
 2140 format(                                                     &
'  -- Vaporization/condensation model (Merkle)',                /,&
'       PRESAT = ', e14.5,    ' (Saturation pressure         )',/,&
'       LINF   = ', e14.5,    ' (Reference length scale      )',/,&
'       UINF   = ', e14.5,    ' (Reference velocity          )',/)
 2150 format(                                                     &
'  -- Eddy-viscosity correction (Reboud correction)',           /,&
'       ICVEVM = ',4x,i10,    ' (Activated (1) or not (0)    )',/,&
'       MCAV   = ', e14.5,    ' (mcav constant               )',/)

 2101 format(                                                     &
                                                                /,&
' ** HOMOGENEOUS MIXTURE MODEL VoF'           ,                 /,&
'    ----------------------------------------',                 /)
 2111 format(                                                     &
'       IVOFMT = ',4x,i10,    ' ( 0  : disabled              )',/,&
'                               ( > 0: enabled               )',/)
 2121 format(                                                     &
'  -- Fluid 1:',                                                /,&
'       RHO1   = ', e14.5,    ' (Reference density           )',/,&
'       MU1    = ', e14.5,    ' (Ref. molecular dyn. visc.   )',/)
 2131 format(                                                     &
'  -- Fluid 2:',                                                /,&
'       RHO2   = ', e14.5,    ' (Reference density           )',/,&
'       MU2    = ', e14.5,    ' (Ref. molecular dyn. visc.   )',/)
 2141 format(                                                     &
'  -- Surface tension:',                                        /,&
'       SIGMA  = ', e14.5,    ' (Surface tension             )',/)
 2151 format(                                                     &
'  -- Drift velocity:',                                         /,&
'       IDRIFT = ',4x,i10,    ' (0: disabled; > 0: enabled   )',/,&
'       KDRIFT = ', e14.5,    ' (Diffusion effect coeff.     )',/,&
'       CDRIFT = ', e14.5,    ' (Drift flux coeff.           )',/)

!===============================================================================
! 4. DISCRETISATION DES EQUATIONS
!===============================================================================

! --- Marche en temps

if (idtvar.ge.0) then

  ! Coefficient de relaxation de la masse volumique

  if (ippmod(icod3p).ge.0 .or. ippmod(islfm).ge.0 .or. &
      ippmod(icoebu).ge.0 .or. ippmod(icolwc).ge.0 .or. &
      ippmod(iccoal).ge.0) then
    write(nfecra,3000)
    write(nfecra,3050) srrom
    write(nfecra,9900)

  endif

endif

 3000 format(                                                     &
                                                                /,&
' ** TIME STEPPING',                                            /,&
'    -------------',                                            /)
 3050 format(                                                     &
'--- Relaxation coefficient',                                   /,&
'    RHO(n+1)=SRROM*RHO(n)+(1-SRROM)*RHO(n+1)',                 /,&
'       SRROM  = ',e14.5,                                       /)

return
end subroutine
