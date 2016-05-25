!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2016 EDF S.A.
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

!===============================================================================

implicit none

! Arguments


! Local variables

character        name*300, chaine*80
integer          iokss , iokcaz
integer          ii    , iiesca, iest
integer          ipp   , iwar  , kval
integer          nbccou, nbsucp, nbvocp, issurf, isvol
integer          kscmin, kscmax, keypp, keyvar
integer          c_id, f_id, f_dim, n_fields
integer          igg, ige
double precision scmaxp, scminp

character(len=3), dimension(3) :: nomext3
character(len=4), dimension(3) :: nomext63

!===============================================================================


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
call field_get_key_id("post_id", keypp)

write(nfecra,1000)

if (ippmod(icod3p).ne.-1) then
  write(nfecra,1010)
  write(nfecra,1020) ippmod(icod3p)
  write(nfecra,1060) indjon
  write(nfecra,1070) namgas, pcigas
  write(nfecra,1080) trim(nomcog(igfuel(1))), &
  -nreact(igoxy(1)), trim(nomcog(igoxy(1))), trim(nomcog(igprod(1)))
  write(nfecra,'(a15,10(1x,a14))') "Composition", "Fuel", "Oxydizer", "Products"
  write(nfecra,'(a15,10(1x,a14))') "-----------", "----", "--------", "--------"
  do ige = 1, ngaze
    write(nfecra,'(a15,10(1x,f14.5))') trim(nomcoe(ige)), (coefeg(ige,igg), igg=1, ngazg)
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
write(nfecra,1520) nvar,nscal,nscaus,nscapp,nproce

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
'       NSCAPP = ',4x,i10,    ' (Nb de scalaires phys. part. )',/,&
'       NPROCE = ',4x,i10,    ' (Nb de proprietes (cellules) )',/)

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
'       NSCAPP = ',4x,i10,    ' (Nb specific physics scalars )',/,&
'       NPROCE = ',4x,i10,    ' (Nb cell properties          )',/)

#endif

!===============================================================================
! 3. MODELISATION PHYSIQUE
!===============================================================================

! --- Proprietes physiques

write(nfecra,2000)
write(nfecra,2010) gx,gy,gz
write(nfecra,2011) icorio

write(nfecra,2020) ro0, viscl0, cp0, icp, p0, pred0, t0,  &
                   irovar,ivivar, (xyzp0(ii),ii=1,3)

if (ippmod(iphpar).ge.1) write(nfecra,2030) diftl0


write(nfecra,9900)


#if defined(_CS_LANG_FR)

 2000 format(                                                     &
                                                                /,&
' ** PROPRIETES PHYSIQUES',                                     /,&
'    --------------------',                                     /)
 2010 format(                                                     &
'       GX     = ', e14.5,    ' (Composante x de la gravite  )',/,&
'       GY     = ', e14.5,    ' (Composante y de la gravite  )',/,&
'       GZ     = ', e14.5,    ' (Composante z de la gravite  )',/)
 2011 format(                                                     &
'       ICORIO = ', i10,      ' (Termes source de Coriolis   )',/)
 2020 format(                                                     &
'  -- Phase continue :',                                        /,&
                                                                /,&
'       RO0    = ', e14.5,    ' (Masse volumique     de ref. )',/,&
'       VISCL0 = ', e14.5,    ' (Visc. molec. dynam. de ref. )',/,&
'       CP0    = ', e14.5,    ' (Chal. Spec.     de reference)',/,&
'       ICP    = ',4x,i10,    ' (> 0 : CP variable   (usphyv))',/,&
'       P0     = ', e14.5,    ' (Pression totale de reference)',/,&
'       PRED0  = ', e14.5,    ' (Press. reduite  de reference)',/,&
'       T0     = ', e14.5,    ' (Temperature     de reference)',/,&
                                                                /,&
'       IROVAR = ',4x,i10,    ' (Masse vol.  cst (0) ou non(1)',/,&
'       IVIVAR = ',4x,i10,    ' (Visc molec. cst (0) ou non(1)',/,&
/,                                                          &
'       Point de reference initial pour la pression',           /,&
'       XYZP0  = ', e14.5, e14.5, e14.5                          )
 2030 format(                                                     &
'       DIFTL0 = ', e14.5,    ' (Diff. dynam.    de reference)',/)

#else

 2000 format(                                                     &
                                                                /,&
' ** PHYSICAL PROPERTIES',                                      /,&
'    -------------------',                                      /)
 2010 format(                                                     &
'       GX     = ', e14.5,    ' (Gravity x component         )',/,&
'       GY     = ', e14.5,    ' (Gravity y component         )',/,&
'       GZ     = ', e14.5,    ' (Gravity z component         )',/)
 2011 format(                                                     &
'       ICORIO = ', i10,      ' (Coriolis source terms       )',/)
 2020 format(                                                     &
'  -- Continuous phase:',                                       /,&
                                                                /,&
'       RO0    = ', e14.5,    ' (Reference density           )',/,&
'       VISCL0 = ', e14.5,    ' (Ref. molecular dyn. visc.   )',/,&
'       CP0    = ', e14.5,    ' (Ref. specific heat          )',/,&
'       ICP    = ',4x,i10,    ' (> 0: variable CP (usphyv)   )',/,&
'       P0     = ', e14.5,    ' (Ref. total pressure         )',/,&
'       PRED0  = ', e14.5,    ' (Ref. reduced pressure       )',/,&
'       T0     = ', e14.5,    ' (Ref. temperature            )',/,&
                                                                /,&
'       IROVAR = ',4x,i10,    ' (Density constant(0) or not(1)',/,&
'       IVIVAR = ',4x,i10,    ' (Molec. visc cst.(0) or not(1)',/,&
/,                                                                &
'       Initial reference point for pressure',                  /,&
'       XYZP0  = ', e14.5, e14.5, e14.5                          )
 2030 format(                                                     &
'       DIFTL0 = ', e14.5,    ' (Ref. dynamic diffusivity    )',/)

#endif

! --- Modele diphasique homogene de cavitation

write(nfecra,2100)
write(nfecra,2110) icavit

if (icavit.ge.0) then

  write(nfecra,2120) rol, mul
  write(nfecra,2130) rov, muv
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
'                               ( 1: modele de Merkle        )' ,/)
 2120 format(                                                     &
'  -- Phase lique :',                                           /,&
'       ROL    = ', e14.5,    ' (Masse volumique     de ref. )',/,&
'       MUL    = ', e14.5,    ' (Visc. molec. dynam. de ref. )',/)
 2130 format(                                                     &
'  -- Phase gazeuse :',                                         /,&
'       ROV    = ', e14.5,    ' (Masse volumique     de ref. )',/,&
'       MUV    = ', e14.5,    ' (Visc. molec. dynam. de ref. )',/)
 2140 format(                                                     &
'  -- Modele de vaporisation/condensation (Merkle)',            /,&
'       PRESAT = ', e14.5,    ' (Pression de saturation      )',/,&
'       LINF   = ', e14.5,    ' (Longueur de reference       )',/,&
'       UINF   = ', e14.5,    ' (Vitesse de reference        )',/)
 2150 format(                                                     &
'  -- Correction de viscosite turbulente de Reboud',            /,&
'       ICVEVM = ',4x,i10,    ' (Active (1) ou non     (0)   )',/,&
'       MCAV   = ', e14.5,    ' (Constante mcav              )',/)

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
'  -- Liquid phase:',                                           /,&
'       ROL    = ', e14.5,    ' (Reference density           )',/,&
'       MUL    = ', e14.5,    ' (Ref. molecular dyn. visc.   )',/)
 2130 format(                                                     &
'  -- Gas phase:',                                              /,&
'       ROV    = ', e14.5,    ' (Reference density           )',/,&
'       MUV    = ', e14.5,    ' (Ref. molecular dyn. visc.   )',/)
 2140 format(                                                     &
'  -- Vaporization/condensation model (Merkle)',                /,&
'       PRESAT = ', e14.5,    ' (Saturation pressure         )',/,&
'       LINF   = ', e14.5,    ' (Reference length scale      )',/,&
'       UINF   = ', e14.5,    ' (Reference velocity          )',/)
 2150 format(                                                     &
'  -- Eddy-viscosity correction (Reboud correction)',           /,&
'       ICVEVM = ',4x,i10,    ' (Activated (1) or not (0)    )',/,&
'       MCAV   = ', e14.5,    ' (mcav constant               )',/)

#endif

! --- Thermique

write(nfecra,2500) itherm, itpscl, iscalt

#if defined(_CS_LANG_FR)

 2500 format(                                                     &
                                                                /,&
' ** MODELE THERMIQUE',                                         /,&
'    ----------------',                                         /,&
' --- Phase continue :',                                        /,&
                                                                /,&
'   - Communs',                                                 /,&
'       ITHERM = ',4x,i10,    ' (0: pas de thermique         )',/,&
'                               (1: temperature              )',/,&
'                               (2: enthalpie                )',/,&
'                               (3: energie totale           )',/,&
'       ITPSCL = ',4x,i10,    ' (0: pas de temperature       )',/,&
'                               (1: temperature en Kelvin    )',/,&
'                               (2: temperature en Celsius   )',/,&
'       ISCALT = ',4x,i10,    ' (Numero du scalaire thermique)',/)

#else

 2500 format(                                                     &
                                                                /,&
' ** THERMAL MODEL',                                            /,&
'    -------------',                                            /,&
' --- Continuous phase:',                                       /,&
                                                                /,&
'   - Commons',                                                 /,&
'       ITHERM = ',4x,i10,    ' (0: no thermal model         )',/,&
'                               (1: temperature              )',/,&
'                               (2: enthalpy                 )',/,&
'                               (3: total energy             )',/,&
'       ITPSCL = ',4x,i10,    ' (0: none                     )',/,&
'                               (1: temperature in Kelvin    )',/,&
'                               (2: temperature in Celsius   )',/,&
'       ISCALT = ',4x,i10,    ' (Thermal scalar number       )',/)

#endif

! --- Turbulence

write(nfecra,2510)

!   - Modeles

write(nfecra,2515) iturb, iwallf, iwallt, ypluli, igrhok
if(iturb.eq.10) then
  write(nfecra,2516)                                            &
       xlomlg
elseif(iturb.eq.20) then
  write(nfecra,2517)                                            &
       almax, uref,                                             &
       iclkep,ikecou,igrake
  if (ikecou.eq.0 .and. idtvar.ge.0) then
    write(nfecra,2527) relaxv(ik),relaxv(iep)
  else
    write(nfecra,2550)
  endif
elseif(iturb.eq.21) then
  write(nfecra,2518) almax, uref, iclkep,ikecou,igrake
  if (ikecou.eq.0.and. idtvar.ge.0) then
    write(nfecra,2527) relaxv(ik),relaxv(iep)
  else
    write(nfecra,2550)
  endif
elseif(iturb.eq.30) then
  write(nfecra,2519)                                            &
       almax, uref, irijco,                                     &
       irijnu,irijrb,irijec,                                    &
       idifre,igrari,iclsyr,iclptr
elseif(iturb.eq.31) then
  write(nfecra,2520) almax, uref,irijco, irijnu,irijrb, igrari,iclsyr,iclptr
elseif(iturb.eq.32) then
  write(nfecra,2525)                                            &
    almax, uref,reinit_turb, irijco,                            &
    irijnu,irijrb,                                              &
    igrari,iclsyr,iclptr
elseif(itytur.eq.4) then
  write(nfecra,2521)                                            &
       csmago,cwale, xlesfl,ales,                               &
       bles,idries,cdries, xlesfd,                              &
       smagmx, ivrtex
elseif(iturb.eq.50) then
  write(nfecra,2522) almax, uref, iclkep,ikecou,igrake
  if (ikecou.eq.0 .and. idtvar.ge.0) then
    write(nfecra,2527) relaxv(ik),relaxv(iep)
  else
    write(nfecra,2550)
  endif
elseif(iturb.eq.51) then
  write(nfecra,2524) almax, uref, iclkep,ikecou,igrake
  if (ikecou.eq.0 .and. idtvar.ge.0) then
    write(nfecra,2527) relaxv(ik),relaxv(iep)
  else
    write(nfecra,2529)
  endif
elseif(iturb.eq.60) then
  write(nfecra,2523) almax, uref, ikecou,igrake
  if (ikecou.eq.0 .and. idtvar.ge.0) then
    write(nfecra,2528) relaxv(ik),relaxv(iomg)
  else
    write(nfecra,2550)
  endif
elseif(iturb.eq.70) then
  write(nfecra,2529) almax,  uref,  relaxv(inusa)
endif
if (itytur.eq.2.or.itytur.eq.5.or.iturb.eq.60.or.iturb.eq.70) then
  write(nfecra,2540) irccor
endif

!   - Constantes

write(nfecra,2530)xkappa,cstlog,apow,bpow

if (iturb.eq.20) then
  write(nfecra,2531)ce1,ce2,sigmak,sigmae,cmu
endif
if (iturb.eq.21) then
  write(nfecra,2532)ce1,ce2,sigmak,sigmae,cmu
endif
if (iturb.eq.30) then
  write(nfecra,2533)ce1,ce2,crij1,crij2,crij3,sigmae,csrij,       &
                    crijp1,crijp2,cmu
endif
if (iturb.eq.31) then
  write(nfecra,2534)cssgs1,cssgs2,cssgr1,cssgr2,cssgr3,cssgr4,    &
       cssgr5,csrij,crij3,ce1,cssge2,sigmae,cmu
endif
if (iturb.eq.32) then
  write(nfecra,2539)cebms1,cebmr1,cebmr2,cebmr3,cebmr4,cebmr5,    &
                    csrij,cebmr6,cebme2,ce1,sigmae,xa1,sigmak,    &
                    xceta,xct
endif
if (iturb.eq.50) then
  write(nfecra,2535) cv2fa1,cv2fe2,sigmak,sigmae,cv2fmu,cv2fct,   &
       cv2fcl,cv2fet,cv2fc1,cv2fc2
endif
if (iturb.eq.51) then
  write(nfecra,2538) cpale1,cpale2,cpale3,cpale4,sigmak,cpalse,cpalmu,cpalct, &
       cpalcl,cpalet,cpalc1,cpalc2
endif
if (iturb.eq.60) then
  write(nfecra,2536) ckwsk1,ckwsk2,ckwsw1,ckwsw2,ckwbt1,ckwbt2,   &
       ckwgm1,ckwgm2,ckwa1,ckwc1,cmu
endif
if (iturb.eq.70) then
  write(nfecra,2537) csab1,csab2,csasig,csav1,csaw1,csaw2,csaw3
endif

iokss = 0
iokcaz = 0
if (irccor.eq.1) then
  if (itycor.eq.1) then
    iokcaz = 1
  elseif (itycor.eq.2) then
    iokss = 1
  endif
endif
if (iokcaz.gt.0) then
  write(nfecra,2541) ccaze2,ccazsc,ccaza,ccazb,ccazc,ccazd
endif
if (iokss.gt.0) then
  write(nfecra,2542) cssr1,cssr2,cssr3
endif

write(nfecra,9900)

#if defined(_CS_LANG_FR)

 2510 format(                                                     &
                                                                /,&
' ** TURBULENCE',                                               /,&
'    ----------',                                               /)
 2515 format(                                                     &
' --- Phase continue :',                                        /,&
                                                                /,&
'   - Communs',                                                 /,&
'       ITURB  = ',4x,i10,    ' (Modele de turbulence        )',/,&
'       IWALLF = ',4x,i10,    ' (loi de paroi                )',/,&
'                               (0: non activee              )',/,&
'                               (1: une echelle loi puissance', /,&
'                                   (interdite en k-epsilon) )',/,&
'                               (2: une echelle loi log      )',/,&
'                               (3: deux echelles loi log    )',/,&
'                               (4: loi de paroi scalable    )',/,&
'                               (5: deux echelles V. Driest  )',/,&
'                               (6: deux echelles lisse/rug. )',/,&
'       IWALLT = ',4x,i10,    ' (correlation coeff. echange  )',/,&
'                               (0: non activee              )',/,&
'                               (1: activee                  )',/,&
'       YPLULI = ', e14.5,    ' (Y plus limite               )',/,&
'                               (1: loi log une echelle      )',/,&
'       IGRHOK = ',4x,i10,    ' (1: Grad (rho k ) calcule    )',/)
 2516 format(                                                     &
'   - Longueur de melange (ITURB = 10)',                        /,&
'       XLOMLG = ', e14.5,    ' (Longueur caracteristique    )',/)
 2517 format(                                                     &
'   - k-epsilon           (ITURB = 20)',                        /,&
'       ALMAX  = ', e14.5,    ' (Longueur caracteristique    )',/,&
'       UREF   = ', e14.5,    ' (Vitesse  caracteristique    )',/,&
'       ICLKEP = ',4x,i10,    ' (Mode de clipping k-epsilon  )',/,&
'       IKECOU = ',4x,i10,    ' (Mode de couplage k-epsilon  )',/,&
'       IGRAKE = ',4x,i10,    ' (Prise en compte de gravite  )')
 2518 format(                                                     &
'   - k-epsilon production lineaire (ITURB = 21)',              /,&
'       ALMAX  = ', e14.5,    ' (Longueur caracteristique    )',/,&
'       UREF   = ', e14.5,    ' (Vitesse  caracteristique    )',/,&
'       ICLKEP = ',4x,i10,    ' (Mode de clipping k-epsilon  )',/,&
'       IKECOU = ',4x,i10,    ' (Mode de couplage k-epsilon  )',/,&
'       IGRAKE = ',4x,i10,    ' (Prise en compte de gravite  )')
 2519 format(                                                     &
'   - Rij-epsilon         (ITURB = 30)',                        /,&
'       ALMAX  = ', e14.5,    ' (Longueur caracteristique    )',/,&
'       UREF   = ', e14.5,    ' (Vitesse  caracteristique    )',/,&
'       IRIJCO = ',4x,i10,    ' (Resolution couplee          )',/,&
'       IRIJNU = ',4x,i10,    ' (Stabilisation matricielle   )',/,&
'       IRIJRB = ',4x,i10,    ' (Reconstruction aux bords    )',/,&
'       IRIJEC = ',4x,i10,    ' (Termes d echo de paroi      )',/,&
'       IDIFRE = ',4x,i10,    ' (Traitmnt du tenseur de diff.)',/,&
'       IGRARI = ',4x,i10,    ' (Prise en compte de gravite  )',/,&
'       ICLSYR = ',4x,i10,    ' (Implicitation en symetrie   )',/,&
'       ICLPTR = ',4x,i10,    ' (Implicitation en paroi      )',/)
 2520 format(                                                     &
'   - Rij-epsilon SSG     (ITURB = 31)',                        /,&
'       ALMAX  = ', e14.5,    ' (Longueur caracteristique    )',/,&
'       UREF   = ', e14.5,    ' (Vitesse  caracteristique    )',/,&
'       IRIJCO = ',4x,i10,    ' (Resolution couplee          )',/,&
'       IRIJNU = ',4x,i10,    ' (Stabilisation matricielle   )',/,&
'       IRIJRB = ',4x,i10,    ' (Reconstruction aux bords    )',/,&
'       IGRARI = ',4x,i10,    ' (Prise en compte de gravite  )',/,&
'       ICLSYR = ',4x,i10,    ' (Implicitation en symetrie   )',/,&
'       ICLPTR = ',4x,i10,    ' (Implicitation en paroi      )',/)
 2521 format(                                                     &
'   - LES                 (ITURB = 40, 41, 42)',                /,&
'                               (Modele de sous-maille       )',/,&
'                               (40 Modele de Smagorinsky    )',/,&
'                               (41 Modele dynamique         )',/,&
'                               (42 Modele WALE              )',/,&
'       CSMAGO = ', e14.5,    ' (Constante de Smagorinski    )',/,&
'       CWALE  = ', e14.5,    ' (Constante du modele WALE    )',/,&
'       XLESFL = ', e14.5,    ' (La largeur du filtre en une )',/,&
'       ALES   = ', e14.5,    ' (cellule s''ecrit            )',/,&
'       BLES   = ', e14.5,    ' (XLESFL*(ALES*VOLUME)**(BLES))',/,&
'       IDRIES = ',4x,i10,    ' (=1 Amortissement Van Driest )',/,&
'       CDRIES = ', e14.5,    ' (Constante de Van Driest     )',/,&
'       XLESFD = ', e14.5,    ' (Rapport entre le filtre     )',/,&
'                               (explicite et le filtre LES  )',/,&
'                               (valeur conseillee 1.5       )',/,&
'       SMAGMX = ', e14.5,    ' (Smagorinsky max dans le cas )',/,&
'                               (du modele dynamique         )',/,&
'       IVRTEX = ',4x,i10,    ' (Utilisation de la methode   )',/,&
'                               (des vortex                  )')
 2522 format(                                                     &
'   - v2f phi-model       (ITURB = 50)',                        /,&
'       ALMAX  = ', e14.5,    ' (Longueur caracteristique    )',/,&
'       UREF   = ', e14.5,    ' (Vitesse  caracteristique    )',/,&
'       ICLKEP = ',4x,i10,    ' (Mode de clipping k-epsilon  )',/,&
'       IKECOU = ',4x,i10,    ' (Mode de couplage k-epsilon  )',/,&
'       IGRAKE = ',4x,i10,    ' (Prise en compte de gravite  )')
 2524 format(                                                     &
'   - v2f BL-v2/k         (ITURB = 51)',                        /,&
'       ALMAX  = ', e14.5,    ' (Longueur caracteristique    )',/,&
'       UREF   = ', e14.5,    ' (Vitesse  caracteristique    )',/,&
'       ICLKEP = ',4x,i10,    ' (Mode de clipping k-epsilon  )',/,&
'       IKECOU = ',4x,i10,    ' (Mode de couplage k-epsilon  )',/,&
'       IGRAKE = ',4x,i10,    ' (Prise en compte de gravite  )')
 2523 format(                                                     &
'   - k-omega SST         (ITURB = 60)',                        /,&
'       ALMAX  = ', e14.5,    ' (Longueur caracteristique    )',/,&
'       UREF   = ', e14.5,    ' (Vitesse  caracteristique    )',/,&
'       IKECOU = ',4x,i10,    ' (Mode de couplage k-omega    )',/,&
'       IGRAKE = ',4x,i10,    ' (Prise en compte de gravite  )')
 2525 format(                                                     &
'   - Rij-epsilon EBRSM     (ITURB = 32)',                      /,&
'       ALMAX  = ', e14.5,    ' (Longueur caracteristique    )',/,&
'       UREF   = ', e14.5,    ' (Vitesse  caracteristique    )',/,&
'       REINIT_                 (Reinitialisation de la',       /,&
'        TURB  = ',4x,i10,    '  turbulence                  )',/,&
'       IRIJCO = ',4x,i10,    ' (Resolution couplee          )',/,&
'       IRIJNU = ',4x,i10,    ' (Stabilisation matricielle   )',/,&
'       IRIJRB = ',4x,i10,    ' (Reconstruction aux bords    )',/,&
'       IGRARI = ',4x,i10,    ' (Prise en compte de gravite  )',/,&
'       ICLSYR = ',4x,i10,    ' (Implicitation en symetrie   )',/,&
'       ICLPTR = ',4x,i10,    ' (Implicitation en paroi      )',/)
 2527 format(                                                     &
'       RELAXV = ', e14.5,    ' pour k       (Relaxation)',     /,&
'       RELAXV = ', e14.5,    ' pour epsilon (Relaxation)',     /)
 2528 format(                                                     &
'       RELAXV = ', e14.5,    ' pour k     (Relaxation)',       /,&
'       RELAXV = ', e14.5,    ' pour omega (Relaxation)',       /)
 2529 format(                                                     &
'   - Spalart-Allmares    (ITURB = 70)',                        /,&
'       ALMAX  = ', e14.5,    ' (Longueur caracteristique    )',/,&
'       UREF   = ', e14.5,    ' (Vitesse  caracteristique    )',/,&
'       RELAXV = ', e14.5,    ' pour nu (Relaxation)',          /)

 2530 format(                                                     &
' --- Constantes',                                              /,&
                                                                /,&
'   - Communs',                                                 /,&
'       XKAPPA = ', e14.5,    ' (Constante de Von Karman     )',/,&
'       CSTLOG = ', e14.5,    ' (U+=Log(y+)/kappa +CSTLOG    )',/,&
'       APOW   = ', e14.5,    ' (U+=APOW (y+)**BPOW (W&W law))',/,&
'       BPOW   = ', e14.5,    ' (U+=APOW (y+)**BPOW (W&W law))',/)
 2531 format(                                                     &
'   - k-epsilon           (ITURB = 20)',                        /,&
'       Ce1    = ', e14.5,    ' (Cepsilon 1 : coef de Prod.  )',/,&
'       CE2    = ', e14.5,    ' (Cepsilon 2 : coef de Diss.  )',/,&
'       SIGMAK = ', e14.5,    ' (Prandtl relatif a k         )',/,&
'       SIGMAE = ', e14.5,    ' (Prandtl relatif a epsilon   )',/,&
'       CMU    = ', e14.5,    ' (Constante Cmu               )',/)
 2532 format(                                                     &
'   - k-epsilon production lineaire (ITURB = 21)',              /,&
'       Ce1    = ', e14.5,    ' (Cepsilon 1 : coef de Prod.  )',/,&
'       CE2    = ', e14.5,    ' (Cepsilon 2 : coef de Diss.  )',/,&
'       SIGMAK = ', e14.5,    ' (Prandtl relatif a k         )',/,&
'       SIGMAE = ', e14.5,    ' (Prandtl relatif a epsilon   )',/,&
'       CMU    = ', e14.5,    ' (Constante Cmu               )',/)
 2533 format(                                                     &
'   - Rij-epsilon std     (ITURB = 30)',                        /,&
'       Ce1    = ', e14.5,    ' (Cepsilon 1 : coef de Prod.  )',/,&
'       CE2    = ', e14.5,    ' (Cepsilon 2 : coef de Diss.  )',/,&
'       CRIJ1  = ', e14.5,    ' (Coef terme lent             )',/,&
'       CRIJ2  = ', e14.5,    ' (Coef terme rapide           )',/,&
'       CRIJ3  = ', e14.5,    ' (Coef terme de gravite       )',/,&
'       SIGMAE = ', e14.5,    ' (Coef sigma_eps              )',/,&
'       CSRIJ  = ', e14.5,    ' (Coef diffusion Rij          )',/,&
'       CRIJP1 = ', e14.5,    ' (Coef lent pour echo de paroi)',/,&
'       CRIJP2 = ', e14.5,    ' (Coef rapide    echo de paroi)',/,&
'       CMU    = ', e14.5,    ' (Constante Cmu               )',/)
 2534 format(                                                     &
'   - Rij-epsilon SSG     (ITURB = 31)',                        /,&
'       CSSGS1 = ', e14.5,    ' (Coef Cs1                    )',/,&
'       CSSGS2 = ', e14.5,    ' (Coef Cs2                    )',/,&
'       CSSGR1 = ', e14.5,    ' (Coef Cr1                    )',/,&
'       CSSGR2 = ', e14.5,    ' (Coef Cr2                    )',/,&
'       CSSGR3 = ', e14.5,    ' (Coef Cr3                    )',/,&
'       CSSGR4 = ', e14.5,    ' (Coef Cr4                    )',/,&
'       CSSGR5 = ', e14.5,    ' (Coef Cr5                    )',/,&
'       CSRIJ  = ', e14.5,    ' (Coef Cs diffusion de Rij    )',/,&
'       CRIJ3  = ', e14.5,    ' (Coef terme de gravite       )',/,&
'       Ce1    = ', e14.5,    ' (Coef Ceps1                  )',/,&
'       CSSGE2 = ', e14.5,    ' (Coef Ceps2                  )',/,&
'       SIGMAE = ', e14.5,    ' (Coef sigma_eps              )',/,&
'       CMU    = ', e14.5,    ' (Constante Cmu               )',/)
 2535 format(                                                     &
'   - v2f phi-model       (ITURB = 50)',                        /,&
'       CV2FA1 = ', e14.5,    ' (a1 pour calculer Cepsilon1  )',/,&
'       CV2FE2 = ', e14.5,    ' (Cepsilon 2 : coef de Diss.  )',/,&
'       SIGMAK = ', e14.5,    ' (Prandtl relatif a k         )',/,&
'       SIGMAE = ', e14.5,    ' (Prandtl relatif a epsilon   )',/,&
'       CV2FMU = ', e14.5,    ' (Constante Cmu               )',/,&
'       CV2FCT = ', e14.5,    ' (Constante CT                )',/,&
'       CV2FCL = ', e14.5,    ' (Constante CL                )',/,&
'       CV2FET = ', e14.5,    ' (Constante C_eta             )',/,&
'       CV2FC1 = ', e14.5,    ' (Constante C1                )',/,&
'       CV2FC2 = ', e14.5,    ' (Constante C2                )',/)
 2536 format(                                                     &
'   - k-omega SST         (ITURB = 60)',                        /,&
'       CKWSK1 = ', e14.5,    ' (Constante sigma_k1          )',/,&
'       CKWSK2 = ', e14.5,    ' (Constante sigma_k2          )',/,&
'       CKWSW1 = ', e14.5,    ' (Constante sigma_omega1      )',/,&
'       CKWSW2 = ', e14.5,    ' (Constante sigma_omega2      )',/,&
'       CKWBT1 = ', e14.5,    ' (Constante beta1             )',/,&
'       CKWBT2 = ', e14.5,    ' (Constante beta2             )',/,&
'       CKWGM1 = ', e14.5,    ' (Constante gamma1            )',/,&
'       CKWGM2 = ', e14.5,    ' (Constante gamma2            )',/,&
'       CKWA1  = ', e14.5,    ' (Cste a1 pour calculer mu_t  )',/,&
'       CKWC1  = ', e14.5,    ' (Cste c1 pour limiteur prod  )',/,&
'       CMU    = ', e14.5,    ' (Cste Cmu (ou Beta*) pour    )',/,&
'                                    conversion omega/epsilon)',/)
 2537 format( &
'   - Spalart-Allmaras    (ITURB = 70)',                        /,&
'       CSAB1  = ', e14.5,    ' (Constante b1                )',/,&
'       CSAB2  = ', e14.5,    ' (Constante b2                )',/,&
'       CSASIG = ', e14.5,    ' (Constante sigma             )',/,&
'       CSAV1  = ', e14.5,    ' (Constante v1                )',/,&
'       CSAW1  = ', e14.5,    ' (Constante w1                )',/,&
'       CSAW2  = ', e14.5,    ' (Constante w2                )',/,&
'       CSAW3  = ', e14.5,    ' (Constante w3                )',/)
 2538 format( &
'   - v2f BL-v2/k         (ITURB = 51)',                        /,&
'       CPALe1 = ', e14.5,    ' (Cepsilon 1 : coef de Prod.  )',/,&
'       CPALE2 = ', e14.5,    ' (Cepsilon 2 : coef de Diss.  )',/,&
'       CPALE3 = ', e14.5,    ' (Cepsilon 3 : coef terme E   )',/,&
'       CPALE4 = ', e14.5,    ' (Cepsilon 4 : coef Diss. mod.)',/,&
'       SIGMAK = ', e14.5,    ' (Prandtl relatif a k         )',/,&
'       CPALSE = ', e14.5,    ' (Prandtl relatif a epsilon   )',/,&
'       CPALMU = ', e14.5,    ' (Constante Cmu               )',/,&
'       CPALCT = ', e14.5,    ' (Constante CT                )',/,&
'       CPALCL = ', e14.5,    ' (Constante CL                )',/,&
'       CPALET = ', e14.5,    ' (Constante C_eta             )',/,&
'       CPALC1 = ', e14.5,    ' (Constante C1                )',/,&
'       CPALC2 = ', e14.5,    ' (Constante C2                )',/)
 2539 format( &
'   - Rij-epsilon EBRSM     (ITURB = 32)',                      /,&
'       CEBMS1 = ', e14.5,    ' (Coef Cs1                    )',/,&
'       CEBMR1 = ', e14.5,    ' (Coef Cr1                    )',/,&
'       CEBMR2 = ', e14.5,    ' (Coef Cr2                    )',/,&
'       CEBMR3 = ', e14.5,    ' (Coef Cr3                    )',/,&
'       CEBMR4 = ', e14.5,    ' (Coef Cr4                    )',/,&
'       CEBMR5 = ', e14.5,    ' (Coef Cr5                    )',/,&
'       CSRIJ  = ', e14.5,    ' (Coef Cs diffusion de Rij    )',/,&
'       CEBMR6 = ', e14.5,    ' (Coef terme de gravite       )',/,&
'       CEBME2 = ', e14.5,    ' (Coef Ceps2                  )',/,&
'       Ce1    = ', e14.5,    ' (Coef Ceps1                  )',/,&
'       SIGMAE = ', e14.5,    ' (Coef sigma_eps              )',/,&
'       XA1    = ', e14.5,    ' (Coef A1                     )',/,&
'       SIGMAK = ', e14.5,    ' (Coef sigma_k                )',/,&
'       XCETA  = ', e14.5,    ' (Coef Ceta                   )',/,&
'       XCT    = ', e14.5,    ' (Coef CT                     )',/)

 2540 format( &
'   - Correction rotation/courbure'                            ,/,&
'       IRCCOR = ',4X,I10,    ' (0: desactivee               )',/,&
'                               (1: activee                  )',/)
 2541 format( &
'   - Correction rotation/courbure (Cazalbou)'                 ,/,&
'       CCAZE2 = ', E14.5,    ' (Coef Ce2^0                  )',/,&
'       CCAZSC = ', E14.5,    ' (Coef Csc                    )',/,&
'       CCAZA  = ', E14.5,    ' (Coef a                      )',/,&
'       CCAZB  = ', E14.5,    ' (Coef b                      )',/,&
'       CCAZC  = ', E14.5,    ' (Coef c                      )',/,&
'       CCAZD  = ', E14.5,    ' (Coef d                      )',/)
 2542 format( &
'   - Correction rotation/courbure (Spalart-Shur)'             ,/,&
'       CSSR1  = ', E14.5,    ' (Coef c_r1                   )',/,&
'       CSSR2  = ', E14.5,    ' (Coef c_r2                   )',/,&
'       CSSR3  = ', E14.5,    ' (Coef c_r3                   )',/)

 2550 format(/)

#else

 2510 format(                                                     &
                                                                /,&
' ** TURBULENCE',                                               /,&
'    ----------',                                               /)
 2515 format(                                                     &
' --- Continuous phase:',                                       /,&
                                                                /,&
'   - Commons',                                                 /,&
'       ITURB  = ',4x,i10,    ' (Turbulence model            )',/,&
'       IWALLF = ',4x,i10,    ' (wall function               )',/,&
'                               (0: disabled                 )',/,&
'                               (1: one scale power law',       /,&
'                                   (forbidden for k-epsilon))',/,&
'                               (2: one scale log law        )',/,&
'                               (3: two scales log law       )',/,&
'                               (4: scalable wall function   )',/,&
'                               (5: two scales V. Driest     )',/,&
'                               (6: two scales smooth/rough  )',/,&
'       IWALLT = ',4x,i10,    ' (Exch. coeff. correlation    )',/,&
'                               (0: not activated            )',/,&
'                               (1: activated                )',/,&
'       YPLULI = ', e14.5,    ' (Limit Y+                    )',/,&
'       IGRHOK = ',4x,i10,    ' (1: computed Grad(rho k)     )',/)
 2516 format(                                                     &
'   - Mixing length       (ITURB = 10)',                        /,&
'       XLOMLG = ', e14.5,    ' (Characteristic length       )',/)
 2517 format(                                                     &
'   - k-epsilon           (ITURB = 20)',                        /,&
'       ALMAX  = ', e14.5,    ' (Characteristic length       )',/,&
'       UREF   = ', e14.5,    ' (Characteristic velocity     )',/,&
'       ICLKEP = ',4x,i10,    ' (k-epsilon clipping model    )',/,&
'       IKECOU = ',4x,i10,    ' (k-epsilon coupling mode     )',/,&
'       IGRAKE = ',4x,i10,    ' (Account for gravity         )')
 2518 format(                                                     &
'   - Linear production k-epsilon (ITURB = 21)',                /,&
'       ALMAX  = ', e14.5,    ' (Characteristic length       )',/,&
'       UREF   = ', e14.5,    ' (Characteristic velocity     )',/,&
'       ICLKEP = ',4x,i10,    ' (k-epsilon clipping model    )',/,&
'       IKECOU = ',4x,i10,    ' (k-epsilon coupling mode     )',/,&
'       IGRAKE = ',4x,i10,    ' (Account for gravity         )')
 2519 format(                                                     &
'   - Rij-epsilon         (ITURB = 30)',                        /,&
'       ALMAX  = ', e14.5,    ' (Characteristic length       )',/,&
'       UREF   = ', e14.5,    ' (Characteristic velocity     )',/,&
'       IRIJCO = ',4x,i10,    ' (Coupled resolution          )',/,&
'       IRIJNU = ',4x,i10,    ' (Matrix stabilization        )',/,&
'       IRIJRB = ',4x,i10,    ' (Reconstruct at boundaries   )',/,&
'       IRIJEC = ',4x,i10,    ' (Wall echo terms             )',/,&
'       IDIFRE = ',4x,i10,    ' (Handle diffusion tensor     )',/,&
'       IGRARI = ',4x,i10,    ' (Prise en compte de gravite  )',/,&
'       ICLSYR = ',4x,i10,    ' (Symmetry implicitation      )',/,&
'       ICLPTR = ',4x,i10,    ' (Wall implicitation          )',/)
 2520 format(                                                     &
'   - SSG Rij-epsilon     (ITURB = 31)',                        /,&
'       ALMAX  = ', e14.5,    ' (Characteristic length       )',/,&
'       UREF   = ', e14.5,    ' (Characteristic velocity     )',/,&
'       IRIJCO = ',4x,i10,    ' (Coupled resolution          )',/,&
'       IRIJNU = ',4x,i10,    ' (Matrix stabilization        )',/,&
'       IRIJRB = ',4x,i10,    ' (Reconstruct at boundaries   )',/,&
'       IGRARI = ',4x,i10,    ' (Account for gravity         )',/,&
'       ICLSYR = ',4x,i10,    ' (Symmetry implicitation      )',/,&
'       ICLPTR = ',4x,i10,    ' (Wall implicitation          )',/)
 2521 format(                                                     &
'   - LES                 (ITURB = 40, 41, 42)',                /,&
'                               (Sub-grid scale model        )',/,&
'                               (40 Smagorinsky model        )',/,&
'                               (41 Dynamic model            )',/,&
'                               (42 WALE model               )',/,&
'       CSMAGO = ', e14.5,    ' (Smagorinsky constant        )',/,&
'       CWALE  = ', e14.5,    ' (WALE model constant         )',/,&
'       XLESFL = ', e14.5,    ' (Filter with in a cell is    )',/,&
'       ALES   = ', e14.5,    ' (written as                  )',/,&
'       BLES   = ', e14.5,    ' (XLESFL*(ALES*VOLUME)**(BLES))',/,&
'       IDRIES = ',4x,i10,    ' (=1 Van Driest damping       )',/,&
'       CDRIES = ', e14.5,    ' (Van Driest constant         )',/,&
'       XLESFD = ', e14.5,    ' (Ratio between the explicit  )',/,&
'                               (filter and LES filter       )',/,&
'                               (recommended value: 1.5      )',/,&
'       SMAGMX = ', e14.5,    ' (Max Smagonsky in the        )',/,&
'                               (dynamic model case          )',/,&
'       IVRTEX = ',4x,i10,    ' (Use the vortex method       )')
 2522 format(                                                     &
'   - v2f phi-model       (ITURB = 50)',                        /,&
'       ALMAX  = ', e14.5,    ' (Characteristic length       )',/,&
'       UREF   = ', e14.5,    ' (Characteristic velocity     )',/,&
'       ICLKEP = ',4x,i10,    ' (k-epsilon clipping model    )',/,&
'       IKECOU = ',4x,i10,    ' (k-epsilon coupling mode     )',/,&
'       IGRAKE = ',4x,i10,    ' (Account for gravity         )')
 2523 format(                                                     &
'   - k-omega SST         (ITURB = 60)',                        /,&
'       ALMAX  = ', e14.5,    ' (Characteristic length       )',/,&
'       UREF   = ', e14.5,    ' (Characteristic velocity     )',/,&
'       IKECOU = ',4x,i10,    ' (k-epsilon coupling mode     )',/,&
'       IGRAKE = ',4x,i10,    ' (Account for gravity         )')
 2524 format(                                                     &
'   - v2f BL-v2/k         (ITURB = 51)',                        /,&
'       ALMAX  = ', e14.5,    ' (Characteristic length       )',/,&
'       UREF   = ', e14.5,    ' (Characteristic velocity     )',/,&
'       ICLKEP = ',4x,i10,    ' (k-epsilon clipping model    )',/,&
'       IKECOU = ',4x,i10,    ' (k-epsilon coupling mode     )',/,&
'       IGRAKE = ',4x,i10,    ' (Account for gravity         )')

 2525 format(                                                     &
'   - Rij-epsilon EBRSM     (ITURB = 32)',                      /,&
'       ALMAX  = ', e14.5,    ' (Characteristic length       )',/,&
'       UREF   = ', e14.5,    ' (Characteristic velocity     )',/,&
'       REINIT_                 (Reinitialization of the',      /,&
'        TURB  = ',4x,i10,    '  turbulence                  )',/,&
'       IRIJCO = ',4x,i10,    ' (Coupled resolution          )',/,&
'       IRIJNU = ',4x,i10,    ' (Matrix stabilization        )',/,&
'       IRIJRB = ',4x,i10,    ' (Reconstruct at boundaries   )',/,&
'       IGRARI = ',4x,i10,    ' (Account for gravity         )',/,&
'       ICLSYR = ',4x,i10,    ' (Symmetry implicitation      )',/,&
'       ICLPTR = ',4x,i10,    ' (Wall implicitation          )',/)
 2527 format(                                                     &
'       RELAXV = ', e14.5,    ' for k        (Relaxation)',     /,&
'       RELAXV = ', e14.5,    ' for epsilon  (Relaxation)',     /)
 2528 format(                                                     &
'       RELAXV = ', e14.5,    ' for k      (Relaxation)',       /,&
'       RELAXV = ', e14.5,    ' for omega  (Relaxation)',       /)
 2529 format(                                                     &
'   - Spalart-Allmaras    (ITURB = 70)',                        /,&
'       ALMAX  = ', e14.5,    ' (Characteristic length       )',/,&
'       UREF   = ', e14.5,    ' (Characteristic velocity     )',/,&
'       RELAXV = ', e14.5,    ' for nu (Relaxation)',           /)

 2530 format(                                                     &
' --- Constants',                                               /,&
                                                                /,&
'   - Commons',                                                 /,&
'       XKAPPA = ', e14.5,    ' (Von Karman constant         )',/,&
'       CSTLOG = ', e14.5,    ' (U+=Log(y+)/kappa +CSTLOG    )',/,&
'       APOW   = ', e14.5,    ' (U+=APOW (y+)**BPOW (W&W law))',/,&
'       BPOW   = ', e14.5,    ' (U+=APOW (y+)**BPOW (W&W law))',/)
 2531 format(                                                     &
'   - k-epsilon           (ITURB = 20)',                        /,&
'       Ce1    = ', e14.5,    ' (Cepsilon 1: production coef.)',/,&
'       CE2    = ', e14.5,    ' (Cepsilon 2: dissipat.  coef.)',/,&
'       SIGMAK = ', e14.5,    ' (Prandtl relative to k       )',/,&
'       SIGMAE = ', e14.5,    ' (Prandtl relative to epsilon )',/,&
'       CMU    = ', e14.5,    ' (Cmu constant                )',/)
 2532 format(                                                     &
'   - Linear production k-epsilon (ITURB = 21)',                /,&
'       Ce1    = ', e14.5,    ' (Cepsilon 1: production coef.)',/,&
'       CE2    = ', e14.5,    ' (Cepsilon 2: dissipat.  coef.)',/,&
'       SIGMAK = ', e14.5,    ' (Prandtl relative to k       )',/,&
'       SIGMAE = ', e14.5,    ' (Prandtl relative to epsilon )',/,&
'       CMU    = ', e14.5,    ' (Cmu constant                )',/)
 2533 format(                                                     &
'   - Rij-epsilon         (ITURB = 30)',                        /,&
'       Ce1    = ', e14.5,    ' (Cepsilon 1: production coef.)',/,&
'       CE2    = ', e14.5,    ' (Cepsilon 2: dissipat.  coef.)',/,&
'       CRIJ1  = ', e14.5,    ' (Slow term coefficient       )',/,&
'       CRIJ2  = ', e14.5,    ' (Fast term coefficient       )',/,&
'       CRIJ3  = ', e14.5,    ' (Gravity term coefficient    )',/,&
'       SIGMAE = ', e14.5,    ' (sigma_eps coeff.            )',/,&
'       CSRIJ  = ', e14.5,    ' (Rij diffusion coeff.        )',/,&
'       CRIJP1 = ', e14.5,    ' (Slow coeff. for wall echo   )',/,&
'       CRIJP2 = ', e14.5,    ' (Fast coeff. for wall echo   )',/,&
'       CMU    = ', e14.5,    ' (Cmu constant                )',/)
 2534 format(                                                     &
'   - SSG Rij-epsilon     (ITURB = 31)',                        /,&
'       CSSGS1 = ', e14.5,    ' (Cs1 coeff.                  )',/,&
'       CSSGS2 = ', e14.5,    ' (Cs2 coeff.                  )',/,&
'       CSSGR1 = ', e14.5,    ' (Cr1 coeff.                  )',/,&
'       CSSGR2 = ', e14.5,    ' (Cr2 coeff.                  )',/,&
'       CSSGR3 = ', e14.5,    ' (Cr3 coeff.                  )',/,&
'       CSSGR4 = ', e14.5,    ' (Cr4 coeff.                  )',/,&
'       CSSGR5 = ', e14.5,    ' (Cr5 coeff.                  )',/,&
'       CSRIJ  = ', e14.5,    ' (Rij Cs diffusion coeff.     )',/,&
'       CRIJ3  = ', e14.5,    ' (Gravity term coeff.         )',/,&
'       Ce1    = ', e14.5,    ' (Ceps1 coeff.                )',/,&
'       CSSGE2 = ', e14.5,    ' (Ceps2 coeff.                )',/,&
'       SIGMAE = ', e14.5,    ' (sigma_eps coeff.            )',/,&
'       CMU    = ', e14.5,    ' (Cmu constant                )',/)
 2535 format(                                                     &
'   - v2f phi-model       (ITURB = 50)',                        /,&
'       CV2FA1 = ', e14.5,    ' (a1 to calculate Cepsilon1   )',/,&
'       CV2FE2 = ', e14.5,    ' (Cepsilon 2: dissip. coeff.  )'/, &
'       SIGMAK = ', e14.5,    ' (Prandtl relative to k       )',/,&
'       SIGMAE = ', e14.5,    ' (Prandtl relative to epsilon )',/,&
'       CV2FMU = ', e14.5,    ' (Cmu constant                )',/,&
'       CV2FCT = ', e14.5,    ' (CT constant                 )',/,&
'       CV2FCL = ', e14.5,    ' (CL constant                 )',/,&
'       CV2FET = ', e14.5,    ' (C_eta constant              )',/,&
'       CV2FC1 = ', e14.5,    ' (C1 constant                 )',/,&
'       CV2FC2 = ', e14.5,    ' (C2 constant                 )',/)
 2536 format(                                                     &
'   - k-omega SST         (ITURB = 60)',                        /,&
'       CKWSK1 = ', e14.5,    ' (sigma_k1 constant           )',/,&
'       CKWSK2 = ', e14.5,    ' (sigma_k2 constant           )',/,&
'       CKWSW1 = ', e14.5,    ' (sigma_omega1 constant       )',/,&
'       CKWSW2 = ', e14.5,    ' (sigma_omega2 constant       )',/,&
'       CKWBT1 = ', e14.5,    ' (beta1 constant              )',/,&
'       CKWBT2 = ', e14.5,    ' (beta2 constant              )',/,&
'       CKWGM1 = ', e14.5,    ' (gamma1 constant             )',/,&
'       CKWGM2 = ', e14.5,    ' (gamma2 constant             )',/,&
'       CKWA1  = ', e14.5,    ' (a1 constant to compute mu_t )',/,&
'       CKWC1  = ', e14.5,    ' (c1 const. for prod. limiter )',/,&
'       CMU    = ', e14.5,    ' (Cmu (or Beta*) constant for )',/,&
'                                    omega/epsilon conversion)',/)
 2537 format(                                                     &
'   - Spalart-Allmaras    (ITURB = 70)',                        /,&
'       CSAB1  = ', e14.5,    ' (b1 constant                 )',/,&
'       CSAB2  = ', e14.5,    ' (b2 constant                 )',/,&
'       CSASIG = ', e14.5,    ' (sigma constant              )',/,&
'       CSAV1  = ', e14.5,    ' (v1 constant                 )',/,&
'       CSAW1  = ', e14.5,    ' (w1 constant                 )',/,&
'       CSAW2  = ', e14.5,    ' (w2 constant                 )',/,&
'       CSAW3  = ', e14.5,    ' (w3 constant                 )',/)
 2538 format( &
'   - v2f BL-v2/k         (ITURB = 51)',                        /,&
'       CPALe1 = ', e14.5,    ' (Cepsilon 1 : Prod. coeff.   )',/,&
'       CPALE2 = ', e14.5,    ' (Cepsilon 2 : Diss. coeff.   )',/,&
'       CPALE3 = ', e14.5,    ' (Cepsilon 3 : E term coeff.  )',/,&
'       CPALE4 = ', e14.5,    ' (Cepsilon 4 : Mod Diss. coef.)',/,&
'       SIGMAK = ', e14.5,    ' (Prandtl relative to k       )',/,&
'       CPALSE = ', e14.5,    ' (Prandtl relative to epsilon )',/,&
'       CPALMU = ', e14.5,    ' (Cmu constant               )',/,&
'       CPALCT = ', e14.5,    ' (CT constant                )',/,&
'       CPALCL = ', e14.5,    ' (CL constant                )',/,&
'       CPALET = ', e14.5,    ' (C_eta constant             )',/,&
'       CPALC1 = ', e14.5,    ' (C1 constant                )',/,&
'       CPALC2 = ', e14.5,    ' (C2 constant                )',/)
 2539  format( &
'   - EBRSM Rij-epsilon     (ITURB = 32)',                      /,&
'       CEBMS1 = ', e14.5,    ' (Cs1 coeff.                  )',/,&
'       CEBMR1 = ', e14.5,    ' (Cr1 coeff.                  )',/,&
'       CEBMR2 = ', e14.5,    ' (Cr2 coeff.                  )',/,&
'       CEBMR3 = ', e14.5,    ' (Cr3 coeff.                  )',/,&
'       CEBMR4 = ', e14.5,    ' (Cr4 coeff.                  )',/,&
'       CEBMR5 = ', e14.5,    ' (Cr5 coeff.                  )',/,&
'       CSRIJ  = ', e14.5,    ' (Rij Cs diffusion coeff.     )',/,&
'       CEBMR6 = ', e14.5,    ' (Gravity term coeff.         )',/,&
'       CEBME2 = ', e14.5,    ' (Coef Ceps2                  )',/,&
'       Ce1    = ', e14.5,    ' (Coef Ceps1                  )',/,&
'       SIGMAE = ', e14.5,    ' (Coef sigma_eps              )',/,&
'       XA1    = ', e14.5,    ' (Coef A1                     )',/,&
'       SIGMAK = ', e14.5,    ' (Coef sigma_k                )',/,&
'       XCETA  = ', e14.5,    ' (Coef Ceta                   )',/,&
'       XCT    = ', e14.5,    ' (Coef CT                     )',/)

 2540 format( &
'   - Rotation/curvature correction'                           ,/,&
'       IRCCOR = ',4X,I10,    ' (0: desactivated             )',/,&
'                               (1: activated                )',/)
 2541 format( &
'   - Rotation/curvature correction (Cazalbou)'                ,/,&
'       CCAZE2 = ', E14.5,    ' (Coef Ce2^0                  )',/,&
'       CCAZSC = ', E14.5,    ' (Coef Csc                    )',/,&
'       CCAZA  = ', E14.5,    ' (Coef a                      )',/,&
'       CCAZB  = ', E14.5,    ' (Coef b                      )',/,&
'       CCAZC  = ', E14.5,    ' (Coef c                      )',/,&
'       CCAZD  = ', E14.5,    ' (Coef d                      )',/)
 2542 format( &
'   - Rotation/curvature correction (Spalart-Shur)'            ,/,&
'       CSSR1  = ', E14.5,    ' (Coef c_r1                   )',/,&
'       CSSR2  = ', E14.5,    ' (Coef c_r2                   )',/,&
'       CSSR3  = ', E14.5,    ' (Coef c_r3                   )',/)

 2550 format(/)

#endif

! --- Viscosite secondaire

write(nfecra,2610)
write(nfecra,2620) ivisse

write(nfecra,9900)


#if defined(_CS_LANG_FR)

 2610 format(                                                     &
                                                                /,&
' ** VISCOSITE SECONDAIRE',                                     /,&
'    --------------------',                                     /)
 2620 format(                                                     &
' --- Phase continue :',                                        /,&
'       IVISSE = ',4x,i10,    ' (1 : pris en compte          )',/)

#else

 2610 format(                                                     &
                                                                /,&
' ** SECONDARY VISCOSITY',                                      /,&
'    -------------------',                                      /)
 2620 format(                                                     &
' --- Continuous phase:',                                       /,&
'       IVISSE = ',4x,i10,    ' (1: accounted for            )',/)

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

!     Stationnaire
if (idtvar.lt.0) then

!   - Parametres du pas de temps

  write(nfecra,3010) idtvar, relxst

!   - Champ de vitesse fige

  write(nfecra,3030) iccvfg

!   - Coefficient de relaxation

  write(nfecra,3011)

  do f_id = 0, n_fields-1
    call field_get_key_int(f_id, keyvar, ii)
    if (ii.lt.0) cycle
    call field_get_label(f_id, chaine)
    write(nfecra,3012) chaine(1:16),relaxv(ii)
  enddo

  write(nfecra,3013)

!     Instationnaire
else

!   - Parametres du pas de temps

  write(nfecra,3020) idtvar,iptlro,coumax,foumax,                 &
       varrdt,dtmin,dtmax,dtref

!   - Champ de vitesse fige

  write(nfecra,3030) iccvfg

!   - Coef multiplicatif du pas de temps

  write(nfecra,3040)
  do f_id = 0, n_fields-1
    call field_get_key_int(f_id, keyvar, ii)
    if (ii.lt.0) cycle
    call field_get_label(f_id, chaine)
    write(nfecra,3041) chaine(1:16),istat(ii),cdtvar(ii)
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
 3010 format(                                                     &
'    ALGORITHME STATIONNAIRE',                                  /,&
                                                                /,&
' --- Parametres globaux',                                      /,&
                                                                /,&
'       IDTVAR = ',4x,i10,    ' (-1: algorithme stationnaire )',/,&
'       RELXST = ', e14.5,    ' (Coef relaxation de reference)',/,&
                                                                /)
 3011 format(                                                     &
' --- Coefficient de relaxation par variable',                  /,&
                                                                /,&
'-----------------------------',                                /,&
' Variable          RELAXV',                                    /,&
'-----------------------------'                                   )
 3012 format(                                                     &
 1x,    a16,      e12.4                                           )
 3013 format(                                                     &
'----------------------------',                                 /,&
                                                                /,&
'       RELAXV =  [0.,1.]       (coeff de relaxation         )', /)
 3020 format(                                                     &
'    ALGORITHME INSTATIONNAIRE',                                /,&
                                                                /,&
' --- Parametres du pas de temps',                              /,&
                                                                /,&
'       IDTVAR = ',4x,i10,    ' (0 cst;1,2 var(tps,tps-espace)',/,&
'       IPTLRO = ',4x,i10,    ' (1 : clipping de DT lie a rho)',/,&
'       COUMAX = ', e14.5,    ' (Courant maximum cible       )',/,&
'       FOUMAX = ', e14.5,    ' (Fourier maximum cible       )',/,&
'       VARRDT = ', e14.5,    ' (En DT var, accroissement max)',/,&
'       DTMIN  = ', e14.5,    ' (Pas de temps min            )',/,&
'       DTMAX  = ', e14.5,    ' (Pas de temps max            )',/,&
'       DTREF  = ', e14.5,    ' (Pas de temps de reference   )',/,&
                                                                /,&
'       En pas de temps non constant (IDTVAR = 1 ou 2),',       /,&
'         lorsque la valeur de COUMAX ou FOUMAX est negative',  /,&
'         ou nulle, la limitation du pas de temps associee (au',/,&
'         nombre de Courant et de Fourier, respectivement)',    /,&
'         n ''est pas prise en compte.',                        /)
 3030 format(                                                     &
' --- Champ de vitesse fige',                                   /,&
                                                                /,&
'       ICCVFG = ',4x,i10,    ' (1 : champ de vitesse fige   )',/)
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
 3010 format(                                                     &
'    STEADY ALGORITHM',                                         /,&
                                                                /,&
' --- Global parameters',                                       /,&
                                                                /,&
'       IDTVAR = ',4x,i10,    ' (-1: steady algorithm        )',/,&
'       RELXST = ', e14.5,    ' (Reference relaxation coeff. )',/,&
                                                                /)
 3011 format(                                                     &
' --- Per variable relaxation coefficient',                     /,&
                                                                /,&
'-----------------------------',                                /,&
' Variable          RELAXV',                                    /,&
'-----------------------------'                                   )
 3012 format(                                                     &
 1x,    a16,      e12.4                                           )
 3013 format(                                                     &
'----------------------------',                                 /,&
                                                                /,&
'       RELAXV =  [0.,1.]       (relaxation coefficient      )',/)
 3020 format(                                                     &
'    UNSTEADY ALGORITHM',                                       /,&
                                                                /,&
' --- Time step parameters',                                    /,&
                                                                /,&
'       IDTVAR = ',4x,i10,    ' (0 cst; 1,2 var (t, t-space  )',/,&
'       IPTLRO = ',4x,i10,    ' (1: rho-related DT clipping  )',/,&
'       COUMAX = ', e14.5,    ' (Maximum target CFL          )',/,&
'       FOUMAX = ', e14.5,    ' (Maximum target Fourier      )',/,&
'       VARRDT = ', e14.5,    ' (For var. DT, max. increase  )',/,&
'       DTMIN  = ', e14.5,    ' (Minimum time step           )',/,&
'       DTMAX  = ', e14.5,    ' (Maximum time step           )',/,&
'       DTREF  = ', e14.5,    ' (Reference time step         )',/,&
                                                                /,&
'       With a non-constant time step (IDTVAR = 1 or 2),',      /,&
'         when the value of COUMAX or FOUMAX is negative',      /,&
'         or zero, the associated time step limitation (for',   /,&
'         CFL and Fourier respectively) is ignored.',           /)
 3030 format(                                                     &
' --- Frozen velocity field',                                   /,&
                                                                /,&
'       ICCVFG = ',4x,i10,    ' (1: frozen velocity field    )',/)
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

! --- Convection Diffusion

write(nfecra,4000)

write(nfecra,4010)

do f_id = 0, n_fields-1
  call field_get_key_int(f_id, keyvar, ii)
  if (ii.lt.0) cycle
  call field_get_label(f_id, chaine)
  write(nfecra,4020) chaine(1:16),                                &
                     iconv(ii),idiff(ii),idifft(ii),              &
                     ischcv(ii),isstpc(ii),                       &
                     blencv(ii),thetav(ii)
enddo
write(nfecra,4030)

write(nfecra,9900)

! --- Stokes

write(nfecra,4110) idilat,iporos,iphydr,icalhy,iprco,ipucou,nterup
write(nfecra,4111) irevmc
if (idtvar.ge.0) then
  write(nfecra,4112) relaxv(ipr),arak
else
  write(nfecra,4113) arak*relaxv(iu)
endif
write(nfecra,4114)istmpf,thetfl,     &
     iroext,thetro,                  &
     iviext,thetvi,                  &
     icpext,thetcp,                  &
     thetsn,thetst,epsup

write(nfecra,9900)


#if defined(_CS_LANG_FR)

 4000 format(                                                     &
                                                                /,&
' ** CONVECTION - DIFFUSION',                                   /,&
'    ----------------------',                                   /)
 4010 format(                                                             &
'---------------------------------------------------------------------',/,&
' Variable          ICONV  IDIFF IDIFFT ISCHCV ISSTPC   BLENCV  THETAV',/,&
'---------------------------------------------------------------------'  )
 4020 format(                                                     &
 1x,    a16,    i7,    i7,    i7,    i7,    i7,    e9.2,    e9.2  )
 4030 format(                                                     &
'-------------------------------------------------------------',/,&
                                                                /,&
'       ICONV  =  0 ou  1       (1 pour convection branchee  )',/,&
'       IDIFF  =  0 ou  1       (1 pour diff. tot branchee   )',/,&
'       IDIFFT =  0 ou  1       (1 pour diff. turb. branchee )',/,&
'       ISCHCV =  0 ou  1       (SOLU ou CD                  )',/,&
'       ISSTPC =  0,1,2,3       (0 : test de pente'            ,/,&
'                                1 : sans test de pente'       ,/,&
'                                2 : min/max limier'           ,/,&
'                                3 : Roe Sweby limiter       )',/,&
'       BLENCV =  [0.;1.]       (1-proportion d upwind       )',/,&
'       THETAV =  [0.;1.]       (0.5 Crank-Nicolson/AB       )',/,&
'                               (theta pour les termes de    )',/,&
'                               (convection diffusion utilise)',/,&
'                               ((1-theta)ancien+theta nouveau',/)

 4110 format(                                                     &
                                                                /,&
' ** STOKES',                                                   /,&
'    ------',                                                   /,&
'       IDILAT = ',4x,i10,  ' (1 : sans prise en compte du',    /,&
'                ',14x,     '      terme instationnaire dans',  /,&
'                ',14x,     '      l''equation de continuite',  /,&
'                ',14x,     '  2 : avec prise en compte du',    /,&
'                ',14x,     '      terme instationnaire dans',  /,&
'                ',14x,     '      l''equation de continuite',  /,&
'       IPOROS = ',4x,i10,  ' (0 : sans modelisation poreuse',  /,&
'                ',14x,     '  1 : avec modelisation poreuse)', /,&
'       IPHYDR = ',4x,i10,  ' (1 : prise en compte explicite',  /,&
'                ',14x,     '      de l''equilibre entre grad', /,&
'                ',14x,     '      de pression et termes',      /,&
'                ',14x,     '      sources de gravite et de',   /,&
'                ',14x,     '      pertes de charge',           /,&
'                ',14x,     '  2 : prise en compte explicite',  /,&
'                ',14x,     '      du desequilibre entre grad', /,&
'                ',14x,     '      de pression et termes',      /,&
'                ',14x,     '      sources de gravite        )',/,&
'       ICALHY = ',4x,i10,  ' (1 : calcul de la pression',      /,&
'                ',14x,     '      hydrostatique pour les',     /,&
'                ',14x,     '      conditions de Dirichlet en', /,&
'                ',14x,     '      sortie sur la pression    )',/,&
'       IPRCO  = ',4x,i10,  ' (1 : avec pression-continuite  )',/,&
'       IPUCOU = ',4x,i10,  ' (1 : avec couplage U-P renforce)',/,&
'       NTERUP = ',4x,i10,  ' (n : avec n sweep sur navsto',    /,&
'                ',14x,     '      pour couplage vites/pressio',/)
 4111 format(                                                     &
'  -- Phase continue :',                                        /,&
                                                                /,&
'       IREVMC = ',4x,i10,    ' (Mode de reconstruction vites)',/)
 4112 format(                                                     &
'       RELAXV = ', e14.5,    ' pour la pression (relaxation)', /,&
'       ARAK   = ', e14.5,    ' (Facteur d Arakawa           )',/)
 4113 format(                                                     &
'       ARAK   = ', e14.5,    ' (Facteur d Arakawa           )',/)
 4114 format(                                                     &
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

 4000 format(                                                     &
                                                                /,&
' ** CONVECTION - DIFFUSION',                                   /,&
'    ----------------------',                                   /)
 4010 format(                                                             &
'---------------------------------------------------------------------',/,&
' Variable          ICONV  IDIFF IDIFFT ISCHCV ISSTPC   BLENCV  THETAV',/,&
'---------------------------------------------------------------------'  )
 4020 format(                                                     &
 1x,    a16,    i7,    i7,    i7,    i7,    i7,    e9.2,    e9.2  )
 4030 format(                                                     &
'-------------------------------------------------------------',/,&
                                                                /,&
'       ICONV  =  0 ou  1       (1 for convection active     )',/,&
'       IDIFF  =  0 ou  1       (1 for total diffusion active)',/,&
'       IDIFFT =  0 ou  1       (1 for turbulent diff. active)',/,&
'       ISCHCV =  0 ou  1       (SOLU or CD                  )',/,&
'       ISSTPC =  0,1,2,3       (0: slope test'                ,/,&
'                                1: no slope test'             ,/,&
'                                2: min/max limiter'           ,/,&
'                                3: Roe Sweby limiter        )',/,&
'       BLENCV =  [0.;1.]       (1-upwind proportion         )',/,&
'       THETAV =  [0.;1.]       (0.5 Crank-Nicolson/AB       )',/,&
'                               (theta for convection-       )',/,&
'                               (diffusion terms uses        )',/,&
'                               ((1-theta).old+theta.new     )',/)

 4110 format(                                                     &
                                                                /,&
' ** STOKES',                                                   /,&
'    ------',                                                   /,&
'       IDILAT = ',4x,i10,  ' (1 : without unsteady term',      /,&
'                ',14x,     '      in the continuity equation', /,&
'                ',14x,     '  2 : with unsteady term in ',     /,&
'                ',14x,     '      the continuity equation)',   /,&
'       IPOROS = ',4x,i10,  ' (0 : without porous media',       /,&
'                ',14x,     '  1 : with porous media)',         /,&
'       IPHYDR = ',4x,i10,  ' (1: account for explicit',        /,&
'                ',14x,     '     balance between pressure',    /,&
'                ',14x,     '     gradient, gravity source',    /,&
'                ',14x,     '     terms, and head losses     )',/,&
'       ICALHY = ',4x,i10,  ' (1: compute hydrastatic',        /, &
'                ',14x,     '     pressure for Dirichlet',      /,&
'                ',14x,     '     conditions for pressure',     /,&
'                ',14x,     '     on outlet                  )',/,&
'       IPRCO  = ',4x,i10,  ' (1: pressure-continuity        )',/,&
'       IPUCOU = ',4x,i10,  ' (1: reinforced U-P coupling    )',/,&
'       NTERUP = ',4x,i10,  ' (n: n sweeps on navsto for',      /,&
'                ',14x,     '     velocity/pressure coupling )',/)
 4111 format(                                                     &
'  -- Continuous phase:',                                       /,&
                                                                /,&
'       IREVMC = ',4x,i10,    ' (Velocity reconstruction mode)',/)
 4112 format(                                                     &
'       RELAXV = ', e14.5,    ' for pressure (relaxation)',     /,&
'       ARAK   = ', e14.5,    ' (Arakawa factor              )',/)
 4113 format(                                                     &
'       ARAK   = ', e14.5,    ' (Arakawa factor              )',/)
 4114 format(                                                     &
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

! --- Calcul des gradients

write(nfecra,4500)

write(nfecra,4510) imrgra, anomax

do f_id = 0, n_fields-1
  call field_get_key_int(f_id, keyvar, ii)
  if (ii.lt.0) cycle
  call field_get_label(f_id, chaine)
  write(nfecra,4520) chaine(1:16),                                 &
    nswrgr(ii),nswrsm(ii),epsrgr(ii),epsrsm(ii),extrag(ii)
enddo
write(nfecra,4511)
do f_id = 0, n_fields-1
  call field_get_key_int(f_id, keyvar, ii)
  if (ii.lt.0) cycle
  call field_get_label(f_id, chaine)
  write(nfecra,4521) chaine(1:16),                                 &
    ircflu(ii),imligr(ii),climgr(ii)
enddo
write(nfecra,4530)

write(nfecra,9900)

! --- Interpolation face des viscosites

write(nfecra,4810) imvisf

write(nfecra,9900)

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
    write(nfecra,4951)                                            &
        nitmay, nswrsy, nswrgy, imligy, ircfly, ischcy,           &
        isstpy, iwarny, ntcmxy,                                   &
        blency, epsily, epsrsy, epsrgy, climgy, extray, coumxy,   &
        epscvy, yplmxy
  endif
  write(nfecra,9900)

endif


#if defined(_CS_LANG_FR)

 4500 format(                                                     &
                                                                /,&
' ** CALCUL DES GRADIENTS',                                     /,&
'    --------------------',                                     /)
 4510 format(                                                     &
'       IMRGRA = ',4x,i10,    ' (Mode de reconstruction      )',/,&
'       ANOMAX = ',e14.5,     ' (Angle de non ortho. limite  )',/,&
'                               (pour moindres carres etendu )',/,&
                                                                /,&
'-------------------------------------------------------------------',  /,&
' Variable         NSWRGR NSWRSM      EPSRGR      EPSRSM      EXTRAG',  /,&
'-------------------------------------------------------------------'    )
 4520 format(                                                     &
 1x,    a16,    i7,    i7,      e12.4,      e12.4,      e12.4     )
 4511 format(                                                     &
'-----------------------------------------------------------',  /,&
                                                                /,&
'-------------------------------------------',                  /,&
' Variable         IRCFLU IMLIGR      CLIMGR',                  /,&
'-------------------------------------------'                     )
 4521 format(                                                     &
 1x,    a16,    i7,    i7,      e12.4                             )
 4530 format(                                                     &
'-----------------------------------',                          /,&
                                                                /,&
'       NSWRGR =                (nb sweep reconstruction grad)',/,&
'       NSWRSM =                (nb sweep reconstruction smb )',/,&
'       EPSRGR =                (precision reconstruction gra)',/,&
'       EPSRSM =                (precision reconstruction smb)',/,&
'       EXTRAG =  [0.;1.]       (extrapolation des gradients )',/,&
'       IRCFLU =  0 ou  1       (reconstruction des flux     )',/,&
'       IMLIGR =  < 0, 0 ou 1   (methode de limit. des grad  )',/,&
'       CLIMGR =  > 1 ou 1      (coef de limitation des grad )',/)

 4810 format(                                                     &
                                                                /,&
' ** INTERPOLATION FACE',                                       /,&
'    ------------------',                                       /,&
'       IMVISF = ',4x,i10,    ' (0 arithmetique              )',/)

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
'       NITMAY = ',4x,i10,    ' (Nb iter pour resolution iter.',/,&
'       NSWRSY = ',4x,i10,    ' (Nb iter pour reconstr. smb. )',/,&
'       NSWRGY = ',4x,i10,    ' (Nb iter pour reconstr. grd. )',/,&
'       IMLIGY = ',4x,i10,    ' (Methode de limitation grd.  )',/,&
'       IRCFLY = ',4x,i10,    ' (Reconst. flux conv. diff.   )',/,&
'       ISCHCY = ',4x,i10,    ' (Schema convectif            )',/,&
'       ISSTPY = ',4x,i10,    ' (Utilisation test de pente   )',/,&
'       IWARNY = ',4x,i10,    ' (Niveau d''impression        )',/,&
'       NTCMXY = ',4x,i10,    ' (Nb iter pour convection stat.',/,&
                                                                /,&
'       BLENCY = ',e14.5,     ' (Prop. ordre 2 schema convect.',/,&
'       EPSILY = ',e14.5,     ' (Precision solveur iteratif  )',/,&
'       EPSRSY = ',e14.5,     ' (Precision reconst. smb.     )',/,&
'       EPSRGY = ',e14.5,     ' (Precision reconst. grd.     )',/,&
'       CLIMGY = ',e14.5,     ' (Coeff. pour limitation grd. )',/,&
'       EXTRAY = ',e14.5,     ' (Coeff. pour extrapolation grd',/,&
'       COUMXY = ',e14.5,     ' (Courant max pour convection )',/,&
'       EPSCVY = ',e14.5,     ' (Precision pour convect. stat.',/,&
'       YPLMXY = ',e14.5,     ' (y+ max avec influence amort.)',/)

#else

 4500 format(                                                     &
                                                                /,&
' ** GRADIENTS CALCULATION',                                    /,&
'    ---------------------',                                    /)
 4510 format(                                                     &
'       IMRGRA = ',4x,i10,    ' (Reconstruction mode         )',/,&
'       ANOMAX = ',e14.5,     ' (Non-ortho angle: limit for  )',/,&
'                               (least squares ext. neighbors)',/,&
                                                                /,&
'-------------------------------------------------------------------',  /,&
' Variable         NSWRGR NSWRSM      EPSRGR      EPSRSM      EXTRAG',  /,&
'-------------------------------------------------------------------'    )
 4520 format(                                                     &
 1x,    a16,    i7,    i7,      e12.4,      e12.4,      e12.4     )
 4511 format(                                                     &
'-----------------------------------------------------------',  /,&
                                                                /,&
'-------------------------------------------',                  /,&
' Variable         IRCFLU IMLIGR      CLIMGR',                  /,&
'-------------------------------------------'                     )
 4521 format(                                                     &
 1x,    a16,    i7,    i7,      e12.4                             )
 4530 format(                                                     &
'-----------------------------------',                          /,&
                                                                /,&
'       NSWRGR =                (nb sweep gradient reconstr. )',/,&
'       NSWRSM =                (nb sweep rhs reconstrcution )',/,&
'       EPSRGR =                (grad. reconstruction prec.  )',/,&
'       EPSRSM =                (rhs   reconstruction prec.  )',/,&
'       EXTRAG =  [0.;1.]       (gradients extrapolation     )',/,&
'       IRCFLU =  0 ou  1       (flow reconstruction         )',/,&
'       IMLIGR =  < 0, 0 ou 1   (gradient limitation method  )',/,&
'       CLIMGR =  > 1 ou 1      (gradient limitation coeff.  )',/)

 4810 format(                                                     &
                                                                /,&
' ** FACE INTERPOLATION',                                       /,&
'    ------------------',                                       /,&
'       IMVISF = ',4x,i10,    ' (0 arithmetic                )',/)

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
'       NITMAY = ',4x,i10,    ' (Nb iter for iter resolution )',/,&
'       NSWRSY = ',4x,i10,    ' (Nb iter for rhs reconstr.   )',/,&
'       NSWRGY = ',4x,i10,    ' (Nb iter for grad. reconstr. )',/,&
'       IMLIGY = ',4x,i10,    ' (Gradient limitation method  )',/,&
'       IRCFLY = ',4x,i10,    ' (Conv. Diff. flow reconstr.  )',/,&
'       ISCHCY = ',4x,i10,    ' (Convective scheme           )',/,&
'       ISSTPY = ',4x,i10,    ' (Slope tet use               )',/,&
'       IWARNY = ',4x,i10,    ' (Verbosity level             )',/,&
'       NTCMXY = ',4x,i10,    ' (Nb iter for steady convect. )',/,&
                                                                /,&
'       BLENCY = ',e14.5,     ' (2nd order conv. scheme prop.)',/,&
'       EPSILY = ',e14.5,     ' (Iterative solver precision  )',/,&
'       EPSRSY = ',e14.5,     ' (rhs reconstruction precision)',/,&
'       EPSRGY = ',e14.5,     ' (Gradient reconstr. precision)',/,&
'       CLIMGY = ',e14.5,     ' (Coeff. for grad. limitation )',/,&
'       EXTRAY = ',e14.5,     ' (Coeff. for grad. extrapolat.)',/,&
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
  write(nfecra,5020) chaine(1:16), epsilo(ii), idircl(ii)
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
'       IDIRCL = 0 ou 1         (decalage de la diagonale si',  /,&
'                                ISTAT=0 et pas de Dirichlet )',/)

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
'       IRESOL =            -1  (automatic solver choice     )',/,&
'                IPOL*1000 + 0  (p conjuguate gradient       )',/,&
'                            1  (Jacobi                      )',/,&
'                IPOL*1000 + 2  (bicgstab                    )',/,&
'                  avec IPOL    (preconditioning degree      )',/,&
'       NITMAX =                (max number of iterations    )',/,&
'       EPSILO =                (resolution precision        )',/,&
'       IDIRCL = 0 ou 1         (shift diagonal if   ',         /,&
'                                ISTAT=0 and no Dirichlet    )',/)

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
    write(nfecra,6021) chaine(1:16),ii,iscacp(ii),      &
                       iturt(ii),visls0(ii),sigmas(ii)
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
write(nfecra,7530) nthist,frhist,ncapt
do f_id = 0, n_fields-1
  call field_get_key_int(f_id, keypp, ipp)
  if (ipp.le.1) cycle
  call field_get_dim (f_id, f_dim)
  do c_id = 1, min(f_dim, 3)
    ii = ipp + c_id - 1
    if (ihisvr(ii,1).ne.0) then
      call field_get_label(f_id, name)
      if (f_dim .eq. 3) then
        name = trim(name) // nomext3(c_id)
      else if (f_dim .eq. 6) then
        name = trim(name) // nomext63(c_id)
      endif
      write(nfecra,7531) ii,name  (1:16),ihisvr(ii,1)
    endif
  enddo
enddo
write(nfecra,7532)

!   - Fichiers listing

write(nfecra,7540) ntlist
do f_id = 0, n_fields-1
  call field_get_key_int(f_id, keypp, ipp)
  if (ipp.lt.1) cycle
  call field_get_key_int(f_id, keyvar, ii)
  if (ii.ge.1) then
    iwar = iwarni(ii)
  else
    iwar = -999
  endif
  call field_get_key_int(f_id, keylog, kval)
  if (kval.eq.1) then
    call field_get_label(f_id, name)
    write(nfecra,7531) ipp,name  (1:16),iwar
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
'       FRHIST = ',4x,e11.5,  ' (Periode de sortie (s))',       /,&
'       NCAPT  = ',4x,i10,    ' (Nombre de capteurs   )',       /,&
                                                                /,&
'       Numero Nom                   Nb. sondes (-1 : toutes)'   )
 7531 format(i10,1X,          A16,6X,         i10                )
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
'       FRHIST = ',4x,e11.5,  ' (Output frequency (s) )',       /,&
'       NCAPT  = ',4x,i10,    ' (Number of probes     )',       /,&
                                                                /,&
'       Number Name                  Nb. probes (-1: all)'       )
 7531 format(i10,1X,          A16,6X,         i10                )
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


! --- Couplage SYRTHES

!     RECUPERATION DU NOMBRE DE CAS DE COUPLAGE

call nbcsyr (nbccou)
!==========

if (nbccou .ge. 1) then

  write(nfecra,8000)
  write(nfecra,8010) nbccou

  nbsucp = 0
  nbvocp = 0

  do ii = 1, nbccou

     ! Add a new surface coupling if detected
     issurf = 0
     call tsursy(ii, issurf)
     nbsucp = nbsucp + issurf

     ! Add a new volume coupling if detected
     isvol = 0
     call tvolsy(ii, isvol)
     nbvocp = nbvocp + isvol

  enddo

  write(nfecra,8020) nbsucp, nbvocp
  write(nfecra,8030)
  do ii = 1, nscal
    f_id = ivarfl(isca(ii))
    call field_get_label(f_id, chaine)
    write(nfecra,8031) chaine(1:16),ii,icpsyr(ii)
  enddo
  write(nfecra,8032)

  write(nfecra,9900)

endif


#if defined(_CS_LANG_FR)

 8000 format(                                                     &
                                                                /,&
' ** COUPLAGE SYRTHES',                                         /,&
'    ----------------',                                         /)
 8010 format(                                                     &
'       NBCCOU = ',4x,i10,    ' (Nombre de couplages         )',/)
 8020 format(                                                     &
'       dont', 8x,i10, ' couplage(s) surfacique(s)',/,            &
'       dont', 8x,i10, ' couplage(s) volumique(s)',/)
 8030 format(                                                     &
                                                                /,&
'  -- Scalaires couples',                                       /,&
'-------------------------------',                              /,&
' Scalaire         Numero ICPSYR',                              /,&
'-------------------------------'                                )
 8031 format(                                                     &
 1x,    a16,    i7,    i7                                        )
 8032 format(                                                     &
'-----------------------',                                      /,&
                                                                /,&
'       ICPSYR = 0 ou 1         (1 : scalaire couple SYRTHES )',/)

#else

 8000 format(                                                     &
                                                                /,&
' ** SYRTHES COUPLING',                                         /,&
'    ----------------',                                         /)
 8010 format(                                                     &
'       NBCCOU = ',4x,i10,    ' (Number of couplings         )',/)
 8020 format(                                                     &
'       with', 8x,i10, ' surface coupling(s)',/,                  &
'       with', 8x,i10, ' volume coupling(s)',/)
 8030 format(                                                     &
                                                                /,&
'  -- Coupled scalars',                                         /,&
'-------------------------------',                              /,&
' Scalar           Number ICPSYR',                              /,&
'-------------------------------'                                )
 8031 format(                                                     &
 1x,    a16,    i7,    i7                                        )
 8032 format(                                                     &
'-----------------------',                                      /,&
                                                                /,&
'       ICPSYR = 0 or 1         (1: scalar coupled to SYRTHES)',/)

#endif

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
