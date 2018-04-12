!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2018 EDF S.A.
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

subroutine impini
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

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
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

character        name*300, chaine*80
integer          ii    , iiesca, iest
integer          iwar  , kval
integer          kscmin, kscmax, keyvar
integer          f_id, n_fields
integer          igg, ige
double precision scmaxp, scminp
double precision turb_schmidt

character(len=3), dimension(3) :: nomext3
character(len=4), dimension(3) :: nomext63

type(var_cal_opt) :: vcopt

!===============================================================================

interface

  subroutine syr_coupling_log_setup()  &
      bind(C, name='cs_syr_coupling_log_setup')
    use, intrinsic :: iso_c_binding
    implicit none
  end subroutine syr_coupling_log_setup

end interface

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
  write(nfecra,1060) indjon
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
else if (ippmod(icoebu).ne.-1) then
  write(nfecra,1010)
  write(nfecra,1030) ippmod(icoebu), cebu
  write(nfecra,1060) indjon
else if (ippmod(icfuel).ne.-1) then
  write(nfecra,1010)
  write(nfecra,1050) ippmod(icfuel)
endif



#if defined(_CS_LANG_FR)

 1000 format(                                                     &
                                                                /,&
' ===========================================================', /,&
                                                                /,&
'               RESUME DES PARAMETRES DE CALCUL',               /,&
'               ===============================',               /,&
                                                                /,&
' -----------------------------------------------------------', /)

 9900 format(                                                     &
                                                                /,&
' -----------------------------------------------------------', /)
 1010 format(                                                     &
                                                                /,&
' ** PHYSIQUE PARTICULIERE :',                                  /,&
'    ---------------------',                                    /)
 1020 format(                                                     &
' --- Flamme de diffusion : Chimie 3 points',                   /,&
'       OPTION = ',4x,i10                                       /)
 1030 format(                                                     &
' --- Flamme premelangee : Modele EBU',                         /,&
'       OPTION = ',4x,i10,                                      /,&
'       CEBU   = ',e14.5                                        /)
 1050 format(                                                     &
' --- Fuel              : Modele Combustible moyen local',       /&
'       OPTION = ',4x,i10                                       /)
 1060 format(                                                     &
' --- Janaf ou non (dans ce cas tabulation utilisateur)',       /,&
'       INDJON = ',4x,i10,    ' (1: Janaf, 0: utilisateur)',    /)
 1070 format(                                                     &
' --- Caracteristiques du combustible',                         /,&
'       Combustible : ',4x,a,                                   /,&
'       PCI = ',4x,e14.5,  ' J/kg',                             /)
 1080 format(                                                     &
" --- Reaction chimique : ",                                    /,&
"       ", a," + ",f6.3," (",a,") --> ",a,                      /)

#else

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
 1050 format(                                                     &
' --- Fuel:            Local Mean Combustible Model',           /,&
'       OPTION = ',4x,i10                                       /)
 1060 format(                                                     &
' --- Janaf or not (user tabulation required in this case)',    /,&
'       INDJON = ',4x,i10,    ' (1: Janaf, 0: user)',           /)
 1070 format(                                                     &
' --- Combustible characteristics',                             /,&
'       Combustible : ',4x,a,                                   /,&
'       PCI = ',4x,e14.5,  ' J/kg',                             /)
 1080 format(                                                     &
" --- Chemical reaction: ",                                     /,&
"       ", a," + ",f6.3," (",a,") --> ",a,                      /)

#endif


!===============================================================================
! 2. DEFINITION GENERALE DU CAS
!===============================================================================

! --- Dimensions

write(nfecra,1500)
write(nfecra,1520) nvar,nscal,nscaus,nscapp

write(nfecra,9900)


#if defined(_CS_LANG_FR)

 1500 format(                                                     &
                                                                /,&
' ** DIMENSIONS',                                               /,&
'    ----------',                                               /)
 1520 format(                                                     &
' --- Physique',                                                /,&
'       NVAR   = ',4x,i10,    ' (Nb de variables             )',/,&
'       NSCAL  = ',4x,i10,    ' (Nb de scalaires             )',/,&
'       NSCAUS = ',4x,i10,    ' (Nb de scalaires utilisateur )',/,&
'       NSCAPP = ',4x,i10,    ' (Nb de scalaires phys. part. )',/)

#else

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

#endif

!===============================================================================
! 3. MODELISATION PHYSIQUE
!===============================================================================

! --- VOF Model

write(nfecra,2101)
write(nfecra,2111) ivofmt

if (ivofmt.ge.0) then

  write(nfecra,2121) rho1, mu1
  write(nfecra,2131) rho2,  mu2

endif

! --- Modele diphasique homogene de cavitation

write(nfecra,2100)
write(nfecra,2110) icavit

if (icavit.ge.0) then

  write(nfecra,2120)
  write(nfecra,2130)

  if (icavit.eq.1) then
    write(nfecra,2140) presat, linf, uinf
  endif
  if (itytur.eq.2 .or. itytur.eq.5 .or. iturb.eq.60 .or. iturb.eq.70) then
    write(nfecra,2150) icvevm, mcav
  endif

endif

write(nfecra,9900)


#if defined(_CS_LANG_FR)

 2100 format(                                                     &
                                                                /,&
' ** MODELE DIPHASIQUE HOMOGENE DE CAVITATION',                 /,&
'    ----------------------------------------',                 /)
 2110 format(                                                     &
'       ICAVIT = ',4x,i10,    ' (-1: ecoulement monophasique )',/,&
'                               ( 0: sans transfert de masse )',/,&
'                               ( 1: modele de Merkle        )',/)
 2120 format(                                                     &
'  -- Phase liquide : fluide 1',                                /)
 2130 format(                                                     &
'  -- Phase gaz : fluide 2',                                     /)
 2140 format(                                                     &
'  -- Modele de vaporisation/condensation (Merkle)',            /,&
'       PRESAT = ', e14.5,    ' (Pression de saturation      )',/,&
'       LINF   = ', e14.5,    ' (Longueur de reference       )',/,&
'       UINF   = ', e14.5,    ' (Vitesse de reference        )',/)
 2150 format(                                                     &
'  -- Correction de viscosite turbulente de Reboud',            /,&
'       ICVEVM = ',4x,i10,    ' (Active (1) ou non     (0)   )',/,&
'       MCAV   = ', e14.5,    ' (Constante mcav              )',/)

 2101 format(                                                     &
                                                                /,&
' ** METHODE VOF                             ',                 /,&
'    ----------------------------------------',                 /)
 2111 format(                                                     &
'       IVOFMT = ',4x,i10,    ' ( -1: inactif                 )',/,&
'                               (  0: actif                   )',/)
 2121 format(                                                     &
'  -- Fluide 1 :',                                              /,&
'       RHO1     = ', e14.5,  ' (Masse volumique     de ref.)', /,&
'       MU1      = ', e14.5,  ' (Visc. molec. dynam. de ref.)', /)
 2131 format(                                                     &
'  -- Fluide 2 :',                                              /,&
'       RHO2     = ', e14.5,  ' (Masse volumique     de ref.)', /,&
'       MU2      = ', e14.5,  ' (Visc. molec. dynam. de ref.)', /)

#else

 2100 format(                                                     &
                                                                /,&
' ** HOMOGENEOUS MIXTURE MODEL FOR CAVITATION',                 /,&
'    ----------------------------------------',                 /)
 2110 format(                                                     &
'       ICAVIT = ',4x,i10,    ' (-1: single phase flow       )',/,&
'                               ( 0: no vap./cond. model     )',/,&
'                               ( 1: Merkle''s model        )' ,/)
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
' ** VOF METHOD                              ',                 /,&
'    ----------------------------------------',                 /)
 2111 format(                                                     &
'       IVOFMT = ',4x,i10,    ' ( -1: inactive                )',/,&
'                               (  0: active                  )',/)
 2121 format(                                                     &
'  -- Fluid 1:',                                                /,&
'       RHO1   = ', e14.5,  ' (Reference density          )',   /,&
'       MU1    = ', e14.5,  ' (Ref. molecular dyn. visc.  )',   /)
 2131 format(                                                     &
'  -- Fluid 2:',                                                /,&
'       RHO2   = ', e14.5,  ' (Reference density          )',   /,&
'       MU2    = ', e14.5,  ' (Ref. molecular dyn. visc.  )',   /)

#endif

! --- Compressible

if (ippmod(icompf).ge.0) then
  write(nfecra,2700)
  write(nfecra,2710) icv, iviscv, viscv0, icfgrp

  write(nfecra,9900)

endif

#if defined(_CS_LANG_FR)

 2700 format(                                                     &
                                                                /,&
' ** COMPRESSIBLE : donnees complementaires',                   /,&
'    ------------',                                             /)
 2710 format(                                                     &
' --- Phase continue :',                                        /,&
'       ICV    = ',4x,i10,    ' (0 : Cv cst ; 1 : variable   )',/,&
'       IVISCV = ',4x,i10,    ' (0 : kappa cst ; 1 : variable', /,&
'                                kappa : viscosite en volume',  /,&
'                                en kg/(m s)                 )',/,&
'       VISCV0 = ',e14.5,     ' (Valeur de kappa si cst      )',/,&
'       ICFGRP = ',4x,i10,    ' (1 : C.L. pression avec effet', /,&
'                                hydrostatique dominant      )',/)

#else

 2700 format(                                                     &
                                                                /,&
' ** COMPRESSIBLE: additional data',                            /,&
'    ------------',                                             /)
 2710 format(                                                     &
' --- Continuous phase :',                                      /,&
'       ICV    = ',4x,i10,    ' (0: Cv cst; 1: variable      )',/,&
'       IVISCV = ',4x,i10,    ' (0: kappa cst; 1: variable',    /,&
'                                kappa: volume viscosity',      /,&
'                                in kg/(m.s)                 )',/,&
'       VISCV0 = ',e14.5,     ' (kappa value if constant     )',/,&
'       ICFGRP = ',4x,i10,    ' (1: pressure BC with dominant', /,&
'                                hydrostatic effect          )',/)

#endif

!===============================================================================
! 4. DISCRETISATION DES EQUATIONS
!===============================================================================

! --- Marche en temps

write(nfecra,3000)

! Instationnaire
if (idtvar.ge.0) then

!   - Coef multiplicatif du pas de temps

  write(nfecra,3040)
  do f_id = 0, n_fields-1
    call field_get_key_int(f_id, keyvar, ii)
    if (ii.lt.0) cycle
    call field_get_label(f_id, chaine)
    call field_get_key_struct_var_cal_opt(f_id, vcopt)
    write(nfecra,3041) chaine(1:16),vcopt%istat,cdtvar(ii)
  enddo
  write(nfecra,3042)


!   - Coefficient de relaxation de la masse volumique

  if (ippmod(iphpar).ge.2) write(nfecra,3050) srrom

!   - Ordre du schema en temps

  write(nfecra,3060)
  write(nfecra,3061) ischtp
  write(nfecra,3062)

endif

write(nfecra,9900)


#if defined(_CS_LANG_FR)

 3000 format(                                                     &
                                                                /,&
' ** MARCHE EN TEMPS',                                          /,&
'    ---------------',                                          /)
 3040 format(                                                     &
' --- Proprietes par variable',                                 /,&
                                                                /,&
'------------------------------------',                         /,&
' Variable          ISTAT      CDTVAR',                         /,&
'------------------------------------'                            )
 3041 format(                                                     &
 1x,    a16,    i7,      e12.4                                    )
 3042 format(                                                     &
'----------------------------',                                 /,&
                                                                /,&
'       ISTAT  =  0 ou  1       (1 pour instationnaire       )',/,&
'       CDTVAR >  0             (coeff mult. du pas de temps )',/)

 3050 format(                                                     &
'--- Coefficient de relaxation',                                /,&
'    RHO(n+1)=SRROM*RHO(n)+(1-SRROM)*RHO(n+1)',                 /,&
'       SRROM  = ',e14.5,                                       /)

 3060 format(                                                     &
' --- Ordre du schema en temps de base'                          )
 3061 format(                                                     &
'       ISCHTP = ',4x,i10,    ' (1 : ordre 1 ; 2 : ordre 2   )'  )
 3062 format(                                                     &
'                                                             '  )

#else

 3000 format(                                                     &
                                                                /,&
' ** TIME STEPPING',                                            /,&
'    -------------',                                            /)
 3040 format(                                                     &
' --- Per-variable properties',                                 /,&
                                                                /,&
'------------------------------------',                         /,&
' Variable          ISTAT      CDTVAR',                         /,&
'------------------------------------'                            )
 3041 format(                                                     &
 1x,    a16,    i7,      e12.4                                    )
 3042 format(                                                     &
'----------------------------',                                 /,&
                                                                /,&
'       ISTAT  =  0 ou  1       (1 for unsteady              )',/,&
'       CDTVAR >  0             (time step multiplier        )',/)

 3050 format(                                                     &
'--- Relaxation coefficient',                                   /,&
'    RHO(n+1)=SRROM*RHO(n)+(1-SRROM)*RHO(n+1)',                 /,&
'       SRROM  = ',e14.5,                                       /)

 3060 format(                                                     &
' --- Order of base time stepping scheme'                        )
 3061 format(                                                     &
'       ISCHTP = ',4x,i10,    ' (1: order 1; 2: order 2      )'  )
 3062 format(                                                     &
'                                                             '  )

#endif

! --- Stokes
write(nfecra,4114)istmpf,thetfl,     &
     iroext,thetro,                  &
     iviext,thetvi,                  &
     icpext,thetcp,                  &
     thetsn,thetst,epsup

write(nfecra,9900)


#if defined(_CS_LANG_FR)

 4114 format(                                                     &
                                                                /,&
' ** STOKES',                                                   /,&
'    ------',                                                   /,&
                                                                /,&
'  -- Phase continue :',                                        /,&
                                                                /,&
'       ISTMPF = ',4x,i10,    ' (schema en temps pour le flux', /,&
'                ',14x,       ' (0 : explicite (THETFL = 0   )',/,&
'                ',14x,       ' (1 : schema std (Saturne 1.0 )',/,&
'                ',14x,       ' (2 : ordre 2   (THETFL = 0.5 )',/,&
'       THETFL = ', e14.5,    ' (theta pour flux de masse    )',/,&
'       IROEXT = ',4x,i10,    ' (extrap. masse volumique',      /,&
'                ',14x,       ' (0 : explicite',                /,&
'                ',14x,       ' (1 : n+thetro avec thetro=1/2', /,&
'                ',14x,       ' (2 : n+thetro avec thetro=1',   /,&
'       THETRO = ', e14.5,    ' (theta pour masse volumique',   /,&
'                               ((1+theta)nouveau-theta ancien',/,&
'       IVIEXT = ',4x,i10,    ' (extrap. viscosite totale',     /,&
'                ',14x,       ' (0 : explicite',                /,&
'                ',14x,       ' (1 : n+thetvi avec thetro=1/2', /,&
'                ',14x,       ' (2 : n+thetvi avec thetro=1',   /,&
'       THETVI = ', e14.5,    ' (theta pour viscosite totale',  /,&
'                               ((1+theta)nouveau-theta ancien',/,&
'       ICPEXT = ',4x,i10,    ' (extrap. chaleur specifique',   /,&
'                ',14x,       ' (0 : explicite',                /,&
'                ',14x,       ' (1 : n+thetcp avec thetro=1/2', /,&
'                ',14x,       ' (2 : n+thetcp avec thetro=1',   /,&
'       THETCP = ', e14.5,    ' (theta schema chaleur spec',    /,&
'                               ((1+theta)nouveau-theta ancien',/,&
'       THETSN = ', e14.5,    ' (theta schema T.S. Nav-Stokes)',/,&
'                               ((1+theta)nouveau-theta ancien',/,&
'       THETST = ', e14.5,    ' (theta schema T.S. Turbulence)',/,&
'                               ((1+theta)nouveau-theta ancien',/,&
'       EPSUP  = ', e14.5,    ' (Test d''arret du couplage',    /,&
'                ',14x,       '  vitesse/pression            )',/)


#else

 4114 format(                                                     &
                                                                /,&
' ** STOKES',                                                   /,&
'    ------',                                                   /,&
                                                                /,&
'  -- Phase continue :',                                        /,&
                                                                /,&
'       ISTMPF = ',4x,i10,    ' (time scheme for flow',         /,&
'                ',14x,       ' (0: explicit (THETFL = 0     )',/,&
'                ',14x,       ' (1: std scheme (Saturne 1.0  )',/,&
'                ',14x,       ' (2: 2nd-order (THETFL = 0.5  )',/,&
'       THETFL = ', e14.5,    ' (theta for mass flow         )',/,&
'       IROEXT = ',4x,i10,    ' (density extrapolation',        /,&
'                ',14x,       ' (0: explicit',                  /,&
'                ',14x,       ' (1: n+thetro with thetro=1/2',  /,&
'                ',14x,       ' (2: n+thetro with thetro=1',    /,&
'       THETRO = ', e14.5,    ' (theta for density',            /,&
'                               ((1+theta).new-theta.old',      /,&
'       IVIEXT = ',4x,i10,    ' (total viscosity extrapolation',/,&
'                ',14x,       ' (0: explicit',                  /,&
'                ',14x,       ' (1: n+thetvi with thetro=1/2',  /,&
'                ',14x,       ' (2: n+thetvi with thetro=1',    /,&
'       THETVI = ', e14.5,    ' (theta for total viscosity',    /,&
'                               ((1+theta).new-theta.old',      /,&
'       ICPEXT = ',4x,i10,    ' (specific heat extrapolation',  /,&
'                ',14x,       ' (0: explicit',                  /,&
'                ',14x,       ' (1: n+thetcp with thetro=1/2',  /,&
'                ',14x,       ' (2: n+thetcp with thetro=1',    /,&
'       THETCP = ', e14.5,    ' (specific heat theta-scheme',   /,&
'                               ((1+theta).new-theta.old',      /,&
'       THETSN = ', e14.5,    ' (Nav-Stokes S.T. theta scheme)',/,&
'                               ((1+theta).new-theta.old',      /,&
'       THETST = ', e14.5,    ' (Turbulence S.T. theta-scheme)',/,&
'                               ((1+theta).new-theta.old',      /,&
'       EPSUP  = ', e14.5,    ' (Velocity/pressure coupling',   /,&
'                ',14x,       '  stop test                   )',/)

#endif

! --- Estimateurs d'erreurs pour Navier-Stokes

iiesca = 0
do iest = 1, nestmx
  if(iescal(iest).gt.0) then
    iiesca = 1
  endif
enddo

if(iiesca.gt.0) then
  write(nfecra,4820)
  write(nfecra,4821)
  do iest = 1, nestmx
    write(nfecra,4822)iest, iescal(iest)
  enddo
  write(nfecra,4823)
  write(nfecra,4824)iespre,iesder,iescor,iestot
  write(nfecra,9900)
endif

! --- Calcul de la distance a la paroi

if(ineedy.eq.1) then

  write(nfecra,4950) icdpar
  if(abs(icdpar).eq.1) then
    write(nfecra,4951) ntcmxy, coumxy, epscvy, yplmxy
  endif
  write(nfecra,9900)

endif


#if defined(_CS_LANG_FR)

 4820 format(                                                     &
                                                                /,&
' ** ESTIMATEURS D''ERREUR POUR NAVIER-STOKES',                 /,&
'    ----------------------------------------',                 /)
 4821 format(                                                     &
'----------------------------------------',                     /,&
' Estimateur      IESCAL (mode de calcul)',                     /,&
'----------------------------------------'                       )
 4822 format(                                                     &
 1x,     i10,2x,    i10                                          )
 4823 format(                                                     &
'----------------------------------------'                       )
 4824 format(                                                     &
                                                                /,&
' Estimateurs possibles :',                                     /,&
' ',i2,' =IESPRE : prediction',                                 /,&
'            L''estimateur est base sur la grandeur',           /,&
'            I = rho_n (u*-u_n)/dt + rho_n u_n grad u*',        /,&
'              - rho_n div (mu+mu_t)_n grad u* + grad P_n',     /,&
'              - reste du smb(u_n, P_n, autres variables_n)',   /,&
' ',i2,' =IESDER : derive',                                     /,&
'            L''estimateur est base sur la grandeur',           /,&
'            I = div (flux de masse corrige apres etape',       /,&
'                                               de pression)',  /,&
'            Idealement nul quand l''equation de Poisson est',  /,&
'              resolue exactement',                             /,&
' ',i2,' =IESCOR : correction',                                 /,&
'            L''estimateur est base sur la grandeur',           /,&
'            I = div (rho_n u_(n+1))',                          /,&
'            Idealement nul quand l''equation de Poisson est',  /,&
'              resolue exactement et que le passage des flux',  /,&
'              de masse aux faces vers les vitesses au centre', /,&
'              se fait dans un espace de fonctions',            /,&
'              a divergence nulle',                             /,&
' ',i2,' =IESTOT : total',                                      /,&
'            L''estimateur est base sur la grandeur',           /,&
'            I = rho_n (u_(n+1)-u_n)/dt',                       /,&
'                                 + rho_n u_(n+1) grad u_(n+1)',/,&
'              - rho_n div (mu+mu_t)_n grad u_(n+1)',           /,&
'                                               + gradP_(n+1)', /,&
'              - reste du smb(u_(n+1), P_(n+1),',               /,&
'                                          autres variables_n)',/,&
'             Le flux du terme convectif est calcule a partir', /,&
'               de u_(n+1) pris au centre des cellules (et',    /,&
'               non pas a partir du flux de masse aux faces',   /,&
'               actualise)',                                    /,&
                                                                /,&
' On evalue l''estimateur selon les valeurs de IESCAL :',       /,&
'   IESCAL = 0 : l''estimateur n''est pas calcule',             /,&
'   IESCAL = 1 : l''estimateur    est     calcule,',            /,&
'               sans contribution du volume  (on prend abs(I))',/,&
'   IESCAL = 2 : l''estimateur    est     calcule,',            /,&
'               avec contribution du volume ("norme L2")',      /,&
'               soit abs(I)*SQRT(Volume_cellule),',             /,&
'               sauf pour IESCOR : on calcule',                 /,&
'                 abs(I)*Volume_cellule pour mesurer',          /,&
'                 un ecart en kg/s',                            /)

 4950 format(                                                     &
                                                                /,&
' ** CALCUL DE LA DISTANCE A LA PAROI',                         /,&
'    --------------------------------',                         /,&
                                                                /,&
'       ICDPAR = ',4x,i10,    ' ( 1: std et relu      si suite',/,&
'                               (-1: std et recalcule si suite',/,&
'                               ( 2: old et relu      si suite',/,&
'                               (-2: old et recalcule si suite',/)
4951  format(                                                     &
                                                                /,&
'       NTCMXY = ',4x,i10,    ' (Nb iter pour convection stat.',/,&
                                                                /,&
'       COUMXY = ',e14.5,     ' (Courant max pour convection )',/,&
'       EPSCVY = ',e14.5,     ' (Precision pour convect. stat.',/,&
'       YPLMXY = ',e14.5,     ' (y+ max avec influence amort.)',/)

#else

 4820 format(                                                     &
                                                                /,&
' ** ERROR ESTIMATORS FOR NAVIER-STOKES',                       /,&
'    ----------------------------------',                       /)
 4821 format(                                                     &
'------------------------------------------',                   /,&
' Estimateur      IESCAL (calculation mode)',                   /,&
'------------------------------------------'                     )
 4822 format(                                                     &
 1x,     i10,2x,    i10                                          )
 4823 format(                                                     &
'----------------------------------------'                       )
 4824 format(                                                     &
                                                                /,&
' Possible estimators:',                                        /,&
' ',i2,' =IESPRE: prediction',                                  /,&
'            The estimatore is based on the quantity',          /,&
'            I = rho_n (u*-u_n)/dt + rho_n u_n grad u*',        /,&
'              - rho_n div (mu+mu_t)_n grad u* + grad P_n',     /,&
'              - remainder of rhs(u_n, P_n, other variables_n)',/,&
' ',i2,' =IESDER: drift',                                       /,&
'            The estimator is based on quantity',               /,&
'            I = div (mass flow corrected after pressure step)',/,&
'            Ideally zero when Poisson''s equation is',         /,&
'              resolved exactly',                               /,&
' ',i2,' =IESCOR: correction',                                  /,&
'            The estimator is based on quantity',               /,&
'            I = div (rho_n u_(n+1))',                          /,&
'            Ideally zero when Poisson''s equation is',         /,&
'              resolved exactly and the passage from mass flow',/,&
'              at faces to velocity at cell centers is done',   /,&
'              in a function space with zero divergence',       /,&
' ',i2,' =IESTOT: total',                                       /,&
'            Estimator is based on the quantity',               /,&
'            I = rho_n (u_(n+1)-u_n)/dt',                       /,&
'                                 + rho_n u_(n+1) grad u_(n+1)',/,&
'              - rho_n div (mu+mu_t)_n grad u_(n+1)',           /,&
'                                               + gradP_(n+1)', /,&
'              - rmainder of rhs(u_(n+1), P_(n+1),',            /,&
'                                          other variables_n)',/, &
'             The convective term flow is calculated from',     /,&
'               u_(n+1) taken at the cell centers (and not',    /,&
'               from the updated mass flow at faces)',          /,&
                                                                /,&
' We evaluate the estimator based on values of IESCAL:',        /,&
'   IESCAL = 0: the estimator is not calculated',               /,&
'   IESCAL = 1: the estimator is calculated, with no',          /,&
'               contribution from the volume (we take abs(I))', /,&
'   IESCAL = 2: the estimator is calculated,',                  /,&
'               with contribution from the volume ("L2 norm")', /,&
'               that is abs(I)*SQRT(Cell_volume),',             /,&
'               except for IESCOR: we calculate',               /,&
'                 abs(I)*Cell_volume to measure',               /,&
'                 a difference in kg/s',                        /)

 4950 format(                                                     &
                                                                /,&
' ** WALL DISTANCE COMPUTATION',                                /,&
'    -------------------------',                                /,&
                                                                /,&
'       ICDPAR = ',4x,i10,    ' ( 1: std, reread if restart',   /,&
'                               (-1: std, recomputed if restrt',/,&
'                               ( 2: old, reread if restart',   /,&
'                               (-2: old, recomputed if restrt',/)
4951  format(                                                     &
                                                                /,&
'       NTCMXY = ',4x,i10,    ' (Nb iter for steady convect. )',/,&
                                                                /,&
'       COUMXY = ',e14.5,     ' (Max CFL for convection      )',/,&
'       EPSCVY = ',e14.5,     ' (Precision for steady conv.  )',/,&
'       YPLMXY = ',e14.5,     ' (y+ max w. damping influence )',/)

#endif


!===============================================================================
! 5. SOLVEURS
!===============================================================================

! --- Solveurs iteratifs de base

write(nfecra,5010)
do f_id = 0, n_fields-1
  call field_get_key_int(f_id, keyvar, ii)
  if (ii.lt.0) cycle
  call field_get_label(f_id, chaine)
  call field_get_key_struct_var_cal_opt(f_id, vcopt)
  write(nfecra,5020) chaine(1:16), vcopt%epsilo, idircl(ii)
enddo
write(nfecra,5030)

write(nfecra,9900)


#if defined(_CS_LANG_FR)

 5010 format(                                                     &
                                                                /,&
' ** SOLVEURS ITERATIFS DE BASE',                               /,&
'    --------------------------',                               /,&
                                                                /,&
'------------------------------------',                         /,&
' Variable              EPSILO IDIRCL',                         /,&
'------------------------------------'                           )
 5020 format(                                                     &
 1x,    a16,    e12.4,    i7                                      )
 5030 format(                                                     &
'------------------------------------',                         /,&
                                                                /,&
'       EPSILO =                (precision de la resolution)',  /,&
'       IDIRCL = 0 ou 1         (decalage de la diagonale si',  /,&
'                                ISTAT=0 et pas de Dirichlet)', /)

#else

 5010 format(                                                     &
                                                                /,&
' ** BASE ITERATIVE SOLVERS',                                   /,&
'    ----------------------',                                   /,&
                                                                /,&
'------------------------------------',                         /,&
' Variable              EPSILO IDIRCL',                         /,&
'------------------------------------'                          )
 5020 format(                                                     &
 1x,    a16,    e12.4,    i7                                      )
 5030 format(                                                     &
'------------------------------------',                          /,&
                                                                 /,&
'       EPSILO =                (resolution precision)',         /,&
'       IDIRCL = 0 ou 1         (shift diagonal if',             /,&
'                                ISTAT=0 and no Dirichlet)',     /)

#endif

!===============================================================================
! 6. SCALAIRES
!===============================================================================

! --- Scalaires

if (nscal.ge.1) then
  write(nfecra,6000)
  write(nfecra,6010)itbrrb
  write(nfecra,6011)
  do ii = 1, nscal
    f_id = ivarfl(isca(ii))
    call field_get_label(f_id, chaine)
    call field_get_key_double(f_id, ksigmas, turb_schmidt)
    write(nfecra,6021) chaine(1:16),ii,iscacp(ii),      &
                       iturt(ii),visls0(ii),turb_schmidt
  enddo
  write(nfecra,6031)
  write(nfecra,6012)
  do ii = 1, nscal
    f_id = ivarfl(isca(ii))
    call field_get_label(f_id, chaine)
    write(nfecra,6022) chaine(1:16),ii, rvarfl(ii)
  enddo
  write(nfecra,6032)
  write(nfecra,6013)
  do ii = 1, nscal
    ! Get the min clipping
    f_id = ivarfl(isca(ii))
    call field_get_key_double(f_id, kscmin, scminp)
    call field_get_key_double(f_id, kscmax, scmaxp)
    call field_get_label(f_id, chaine)
    write(nfecra,6023) chaine(1:16),ii,iclvfl(ii),      &
                       scminp,scmaxp
  enddo
  write(nfecra,6033)
  write(nfecra,6030)
  write(nfecra,6040)
  do ii = 1, nscal
    write(nfecra,6041) ii,thetss(ii),ivsext(ii),thetvs(ii)
  enddo
  write(nfecra,6042)

  write(nfecra,9900)

endif


#if defined(_CS_LANG_FR)

 6000 format( &
                                                                       /,&
' ** SCALAIRES',                                                       /,&
'    ---------',/)
 6010 format(                                                            &
'       ITBRRB = ',4x,i10,    ' (Reconstruction T ou H au brd)',/)
 6011 format(                                                            &
'-------------------------------------------------------------',       /,&
' Variable         Numero ISCACP  ITURT     VISLS0      SIGMAS',       /,&
'-------------------------------------------------------------'  )
 6021 format( &
 1x,    a16,    i7,    i7,    i7,      e12.4,      e12.4  )
 6031 format( &
'---------------------------------------------------------------------',/)
 6012 format( &
'------------------------------------',                         /,&
' Variable         Number      RVARFL',                         /,&
'------------------------------------' )
 6022 format( &
 1x,    a16,           i7,    e12.4 )
 6032 format(                                                     &
'-------------------------------------------',                   /)
 6013 format(                                                     &
'-------------------------------------------------------',      /,&
' Variable         Numero ICLVFL      SCAMIN      SCAMAX',      /,&
'-------------------------------------------------------'         )
 6023 format(                                                     &
 1x,    a16,    i7,    i7,      e12.4,      e12.4         )
 6033 format(                                                     &
'-------------------------------------------------------',       /)
 6030 format(                                                     &
'-------------------------------------------------------------',/,&
                                                                /,&
'       Le numero indique pour chaque scalaire le rang',        /,&
'         dans la liste de tous les scalaires. Les scalaires',  /,&
'         utilisateurs sont places en tete, de 1 a NSCAUS. Les',/,&
'         scalaires physique particuliere sont a la fin, de',   /,&
'         NSCAUS+1 a NSCAPP+NSCAUS=NSCAL.',                     /,&
                                                                /,&
'       ISCACP = 0 ou 1         (Utilisation de Cp ou non    )',/,&
'       VISLS0 = >0             (Viscosite de reference      )',/,&
'       SIGMAS = >0             (Schmidt                     )',/,&
'       RVARFL = >0             (Rf, cf dissipation variance )',/,&
'       ICLVFL = 0, 1 ou 2      (Mode de clipping variance   )',/,&
'       SCAMIN =                (Valeur min autorisee        )',/,&
'       SCAMAX =                (Valeur max autorisee        )',/,&
'        Pour les variances, SCAMIN est ignore et SCAMAX n est',/,&
'          pris en compte que si ICLVFL = 2',                   /)
 6040 format(                                                     &
'------------------------------------------------------',       /,&
'   Scalaire      THETSS    IVSEXT      THETVS',                /,&
'------------------------------------------------------'         )
 6041 format(                                                     &
 1x,     i10,      e12.4,      i10,      e12.4                   )
 6042 format(                                                     &
'------------------------------------------------------',       /,&
                                                                /,&
'       THETSS =                (theta pour termes sources   )',/,&
'                               ((1+theta)nouveau-theta ancien',/,&
'       IVSEXT =                (extrap. viscosite totale    )',/,&
'                               (0 : explicite               )',/,&
'                               (1 : n+thetvs avec thetvs=1/2', /,&
'                               (2 : n+thetvs avec thetvs=1  )',/,&
'       THETVS =                (theta pour diffusiv. scalaire',/,&
'                               ((1+theta)nouveau-theta ancien',/)

#else

 6000 format(                                                     &
                                                                /,&
' ** SCALARS',                                                  /,&
'    -------',                                                  /)
 6010 format(                                                     &
'       ITBRRB = ',4x,i10,    ' (T or H reconstruction at bdy)',/)
 6011 format(                                                     &
'-------------------------------------------------------------',/,&
' Variable         Number ISCACP  ITURT     VISLS0      SIGMAS',/,&
'-------------------------------------------------------------'  )
 6021 format( &
 1x,    a16,    i7,    i7,    i7,     e12.4,      e12.4  )
 6031 format( &
'-------------------------------------------------------------',/)
 6012 format( &
'------------------------------------',                         /,&
' Variable         Number      RVARFL',                         /,&
'------------------------------------' )
 6022 format( &
 1x,    a16,           i7,    e12.4 )
 6032 format(                                                     &
'-------------------------------------------',                   /)
 6013 format(                                                     &
'-------------------------------------------------------',      /,&
' Variable         Number ICLVFL      SCAMIN      SCAMAX',      /,&
'-------------------------------------------------------'        )
 6023 format(                                                     &
 1x,    a16,    i7,    i7,      e12.4,      e12.4         )
 6033 format(                                                     &
'-------------------------------------------------------',       /)
 6030 format(                                                     &
'-------------------------------------------------------------',/,&
                                                                /,&
'       For each scalar, the number indicates it''s rank',      /,&
'         in the list of all scalars. User scalars are placed', /,&
'         first, from 1 to NSCAUS. Specific physics scalars',   /,&
'         are placed at the end, from',                         /,&
'         NSCAUS+1 to NSCAPP+NSCAUS=NSCAL.',                    /,&
                                                                /,&
'       ISCACP = 0 or 1     2   (use Cp or not               )',/,&
'       VISLS0 = >0             (Reference viscosity         )',/,&
'       SIGMAS = >0             (Schmidt                     )',/,&
'       RVARFL = >0             (Rf, cf variance dissipation )',/,&
'       ICLVFL = 0, 1 or 2      (Variance clipping mode      )',/,&
'       SCAMIN =                (Min authorized value        )',/,&
'       SCAMAX =                (Max authorized value        )',/,&
'        For variances, SCAMIN is ignored and SCAMAX is used',  /,&
'          only if ICLVFL = 2',                                 /)
 6040 format(                                                     &
'------------------------------------------------------',       /,&
'   Scalar        THETSS    IVSEXT      THETVS',                /,&
'------------------------------------------------------'         )
 6041 format(                                                     &
 1x,     i10,      e12.4,      i10,      e12.4                   )
 6042 format(                                                     &
'------------------------------------------------------',       /,&
                                                                /,&
'       THETSS =                (theta for source terms      )',/,&
'                               ((1+theta).new-theta.old     )',/,&
'       IVSEXT =                (extrap. total viscosity     )',/,&
'                               (0: explicit                 )',/,&
'                               (1: n+thetvs with thetvs=1/2 )',/,&
'                               (2: n+thetvs with thetvs=1   )',/,&
'       THETVS =                (theta for scalar diffusivity', /,&
'                               ((1+theta).new-theta.old     )',/)

#endif

!===============================================================================
! 7. GESTION DU CALCUL
!===============================================================================

! --- Gestion du calcul

write(nfecra,7000)

!   - Suite de calcul

write(nfecra,7010) isuite, ileaux, iecaux

!   - Duree du calcul

write(nfecra,7110) inpdt0,ntmabs

!   - Marge en temps CPU

write(nfecra,7210) tmarus

write(nfecra,9900)


#if defined(_CS_LANG_FR)

 7000 format(                                                     &
                                                                /,&
' ** GESTION DU CALCUL',                                        /,&
'    -----------------',                                        /)
 7010 format(                                                     &
' --- Suite de calcul',                                         /,&
'       ISUITE = ',4x,i10,    ' (1 : suite de calcul         )',/,&
'       ILEAUX = ',4x,i10,    ' (1 : lecture  de restart/auxiliary)',/,&
'       IECAUX = ',4x,i10,    ' (1 : ecriture de checkpoint/auxiliary)',/,&
                                                                /)
 7110 format(                                                     &
' --- Duree du calcul',                                         /,&
'     La numerotation des pas de temps et la mesure du temps',  /,&
'       physique simule sont des valeurs absolues',             /,&
'       et non pas des valeurs relatives au calcul en cours.',  /,&
                                                                /,&
'       INPDT0 = ',4x,i10,    ' (1 : calcul a zero pas de tps)',/,&
'       NTMABS = ',4x,i10,    ' (Pas de tps final demande    )',/)
 7210 format(                                                     &
' --- Marge en temps CPU',                                      /,&
'       TMARUS = ', e14.5,    ' (Marge CPU avant arret       )',/)
#else

 7000 format(                                                     &
                                                                /,&
' ** CALCULATION MANAGEMENT',                                   /,&
'    ----------------------',                                   /)
 7010 format(                                                     &
' --- Restarted calculation',                                   /,&
'       ISUITE = ',4x,i10,    ' (1: restarted calculation    )',/,&
'       ILEAUX = ',4x,i10,    ' (1: read  restart/auxiliary  )',/,&
'       IECAUX = ',4x,i10,    ' (1: write checkpoint/auxiliary)',/,&
                                                                /)
 7110 format(                                                     &
' --- Calculation time',                                        /,&
'     The numbering of time steps and the measure of simulated',/,&
'       physical time are absolute values, and not values',     /,&
'       relative to the current calculation.',                  /,&
                                                                /,&
'       INPDT0 = ',4x,i10,    ' (1: 0 time step calcuation   )',/,&
'       NTMABS = ',4x,i10,    ' (Final time step required    )',/)
 7210 format(                                                     &
' --- CPU time margin',                                         /,&
'       TMARUS = ', e14.5,    ' (CPU time margin before stop )',/)

#endif

!===============================================================================
! 8. ENTREES SORTIES
!===============================================================================

write(nfecra,7500)

!   - Fichier suite

write(nfecra,7510) ntsuit

!   - Fichiers historiques
write(nfecra,7530) nthist,frhist
write(nfecra,7532)

!   - Fichiers listing

write(nfecra,7540) ntlist
do f_id = 0, n_fields-1
  call field_get_key_int(f_id, keyvar, ii)
  if (ii.ge.1) then
    call field_get_key_struct_var_cal_opt(f_id, vcopt)
    iwar = vcopt%iwarni
  else
    iwar = -999
  endif
  call field_get_key_int(f_id, keylog, kval)
  if (kval.eq.1) then
    call field_get_label(f_id, name)
    write(nfecra,7531) name(1:16), iwar
  endif
enddo
write(nfecra,7532)

!   - Post-traitement automatique (bord)

write(nfecra,7550)   'ipstfo', ipstdv(ipstfo),                     &
                     'ipstyp', ipstdv(ipstyp),                     &
                     'ipsttp', ipstdv(ipsttp),                     &
                     'ipstft', ipstdv(ipstft),                     &
                     'ipstnu', ipstdv(ipstnu)

write(nfecra,9900)


#if defined(_CS_LANG_FR)

 7500 format(                                                     &
                                                                /,&
' ** ENTREES SORTIES',                                          /,&
'    ---------------',                                          /)
 7510 format(                                                     &
' --- Fichier suite',                                           /,&
'       NTSUIT = ',4x,i10,    ' (Periode de sauvegarde)',       /)
 7530 format(                                                     &
' --- Fichiers historiques',                                    /,&
'       NTHIST = ',4x,i10,    ' (Periode de sortie    )',       /,&
'       FRHIST = ',4x,e11.5,  ' (Periode de sortie (s))')
 7531 format(1X,          A16,6X,         i10                )
 7532 format(                                                     &
'         --           --                --',                   /)
 7540 format(                                                     &
' --- Fichiers listing',                                        /,&
'       NTLIST = ',4x,i10,    ' (Periode de sortie    )',       /,&
                                                                /,&
'       Numero Nom                 Niveau d''impression IWARNI',/,&
'                                      (-999 : non applicable)',/)
 7550 format(                                                     &
' --- Variables supplementaires en post-traitement (ipstdv)',   /,&
'       ',a6,' = ',4x,i10,    ' (Force exercee par',            /,&
'       ',6x,'   ',4x,10x,    '   le fluide sur le bord)',      /,&
'       ',a6,' = ',4x,i10,    ' (y+ au bord)',                  /,&
'       ',a6,' = ',4x,i10,    ' (T+ au bord)',                  /,&
'       ',a6,' = ',4x,i10,    ' (Flux thermique au bord)',      /,&
'       ',a6,' = ',4x,i10,    ' (Flux thermique',               /,&
'       ',6x,'   ',4x,10x,    '  sans dimension au bord)',      /)

#else

 7500 format(                                                     &
                                                                /,&
' ** INPUT-OUTPUT',                                             /,&
'    ------------',                                             /)
 7510 format(                                                     &
' --- Restart file',                                            /,&
'       NTSUIT = ',4x,i10,    ' (Checkpoint frequency )',       /)
 7530 format(                                                     &
' --- Probe history files',                                     /,&
'       NTHIST = ',4x,i10,    ' (Output frequency     )',       /,&
'       FRHIST = ',4x,e11.5,  ' (Output frequency (s) )')
 7531 format(1X,          A16,6X,         i10                )
 7532 format(                                                     &
'         --           --                --',                   /)
 7540 format(                                                     &
' --- Log files',                                               /,&
'       NTLIST = ',4x,i10,    ' (Output frequency     )',       /,&
                                                                /,&
'       Number Name                IWARNI verbosity level',     /,&
'                                      (-999: not applicable)', /)
 7550 format(                                                     &
' --- Additional post-processing variables (ipstdv)',           /,&
'       ',a6,' = ',4x,i10,    ' (Force exerted by the',         /,&
'       ',6x,'   ',4x,10x,    '       fluid on the boundary)',  /,&
'       ',a6,' = ',4x,i10,    ' (y+ at boundary)',              /,&
'       ',a6,' = ',4x,i10,    ' (T+ at boundary)',              /,&
'       ',a6,' = ',4x,i10,    ' (Thermal flux   at boundary)',  /,&
'       ',a6,' = ',4x,i10,    ' (Dimensionless thermal',        /,&
'       ',6x,'   ',4x,10x,    '            flux at boundary)',  /)

#endif


!===============================================================================
! 9. COUPLAGES
!===============================================================================

call syr_coupling_log_setup

!===============================================================================
! 10. METHODE ALE
!===============================================================================
! --- Activation de la methode ALE

write(nfecra,8210)
write(nfecra,8220) iale, nalinf, iflxmw

write(nfecra,9900)


#if defined(_CS_LANG_FR)

 8210 format(                                                     &
                                                                /,&
' ** METHODE ALE (MAILLAGE MOBILE)',                            /,&
'    -----------',                                              /)
 8220 format(                                                     &
'       IALE   = ',4x,i10,    ' (1 : activee                 )',/ &
'       NALINF = ',4x,i10,    ' (Iterations d''initialisation', / &
'                                                   du fluide)',/ &
'       IFLXMW = ',4x,i10,    ' (Calcul du flux de masse ALE',  / &
'                                0 : déplacement des sommets',  / &
'                                1 : vitesse ALE)',/)

#else

 8210 format(                                                     &
                                                                /,&
' ** ALE METHOD (MOVING MESH)',                                 /,&
'    -----------',                                              /)
 8220 format(                                                     &
'       IALE   = ',4x,i10,    ' (1: activated                )',/ &
'       NALINF = ',4x,i10,    ' (Fluid initialization',         / &
'                                                  iterations)',/ &
'       IFLXMW = ',4x,i10,    ' (ALE mass flux computation',    / &
'                                0: thanks to vertices',        / &
'                                1: thanks to mesh velocity)',/)

#endif

!===============================================================================
! 11. FIN
!===============================================================================

return
end subroutine
