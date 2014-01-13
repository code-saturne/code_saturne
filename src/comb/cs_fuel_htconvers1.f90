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

subroutine cs_fuel_htconvers1 &
!============================
 ( mode  , eh     , xesp   , tp  )

!===============================================================================
! FONCTION :
! --------
! CALCUL DE LA TEMPERATURE DU GAZ
!  EN FONCTION DE L'ENTHALPIE DU GAZ ET DES CONCENTRATIONS
!  SI MODE = 1
! CALCUL DE L'ENTHALPIE DU GAZ
!  EN FONCTION DE LA TEMPERATURE DU GAZ ET DES CONCENTRATIONS
!  SI MODE = -1

! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! eh               ! tr ! <-- ! enthalpie du gaz                               !
!                  !    !     ! (j/kg de melange gazeux)                       !
! xesp             ! tr ! <-- ! fraction massique des especes                  !
! tp               ! tr ! --> ! temperature du gaz (kelvin)                    !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHAMNUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

!==============================================================================
! Module files
!==============================================================================

use paramx
use numvar
use optcal
use cstphy
use cstnum
use entsor
use ppppar
use ppthch
use coincl
use cpincl
use cs_fuel_incl
use ppincl
use ppcpfu

!===============================================================================

implicit none

! Arguments

integer          mode
double precision eh,tp
double precision xesp(ngazem)

! Local variables

integer          ii

double precision eh0 , eh1

!===============================================================================
! 1. CALCUL DE LA TEMPERATURE A PARTIR DE l'ENTHALPIE
!===============================================================================

if ( mode .eq. 1 ) then

  ii = npo


! --- Clipping eventuel de TP a TH(NPO) si EH > EH1



  eh1 =  xesp(ifov )*ehgaze(ifov,ii)                               &
       + xesp(ico  )*ehgaze(ico ,ii)                               &
       + xesp(ih2s )*ehgaze(ih2s,ii)                               &
       + xesp(ihy  )*ehgaze(ihy ,ii)                               &
       + xesp(ihcn )*ehgaze(ihcn,ii)                               &
       + xesp(io2  )*ehgaze(io2 ,ii)                               &
       + xesp(ico2 )*ehgaze(ico2,ii)                               &
       + xesp(ih2o )*ehgaze(ih2o,ii)                               &
       + xesp(iso2 )*ehgaze(iso2,ii)                               &
       + xesp(inh3 )*ehgaze(inh3,ii)                               &
       + xesp(in2  )*ehgaze(in2 ,ii)

  if ( eh .ge. eh1 ) then
    tp = th(ii)
    goto 501
  endif

  ii = 1

! --- Clipping eventuel de TP a TH(1) si EH < EH0

  eh0 =  xesp(ifov )*ehgaze(ifov,ii)                               &
       + xesp(ico  )*ehgaze(ico ,ii)                               &
       + xesp(ih2s )*ehgaze(ih2s,ii)                               &
       + xesp(ihy  )*ehgaze(ihy ,ii)                               &
       + xesp(ihcn )*ehgaze(ihcn,ii)                               &
       + xesp(io2  )*ehgaze(io2 ,ii)                               &
       + xesp(ico2 )*ehgaze(ico2,ii)                               &
       + xesp(ih2o )*ehgaze(ih2o,ii)                               &
       + xesp(iso2 )*ehgaze(iso2,ii)                               &
       + xesp(inh3 )*ehgaze(inh3,ii)                               &
       + xesp(in2  )*ehgaze(in2 ,ii)

  if ( eh .le. eh0 ) then
    tp= th(ii)
    goto 501
  endif


 500    continue
  ii = ii + 1

! --- Clipping eventuel de TP a TH(II-1) si EH < EH0

  eh0 =  xesp(ifov )*ehgaze(ifov,ii-1)                             &
       + xesp(ico  )*ehgaze(ico ,ii-1)                             &
       + xesp(ih2s )*ehgaze(ih2s,ii-1)                             &
       + xesp(ihy  )*ehgaze(ihy ,ii-1)                             &
       + xesp(ihcn )*ehgaze(ihcn,ii-1)                             &
       + xesp(io2  )*ehgaze(io2 ,ii-1)                             &
       + xesp(ico2 )*ehgaze(ico2,ii-1)                             &
       + xesp(ih2o )*ehgaze(ih2o,ii-1)                             &
       + xesp(iso2 )*ehgaze(iso2,ii-1)                             &
       + xesp(inh3 )*ehgaze(inh3,ii-1)                             &
       + xesp(in2  )*ehgaze(in2 ,ii-1)

  eh1 =  xesp(ifov )*ehgaze(ifov,ii)                               &
       + xesp(ico  )*ehgaze(ico ,ii)                               &
       + xesp(ih2s )*ehgaze(ih2s,ii)                               &
       + xesp(ihy  )*ehgaze(ihy ,ii)                               &
       + xesp(ihcn )*ehgaze(ihcn,ii)                               &
       + xesp(io2  )*ehgaze(io2 ,ii)                               &
       + xesp(ico2 )*ehgaze(ico2,ii)                               &
       + xesp(ih2o )*ehgaze(ih2o,ii)                               &
       + xesp(iso2 )*ehgaze(iso2,ii)                               &
       + xesp(inh3 )*ehgaze(inh3,ii)                               &
       + xesp(in2  )*ehgaze(in2 ,ii)

  if ( eh .ge. eh0  .and. eh .le. eh1  ) then
    tp = th(ii-1) + (eh-eh0) *                                     &
                 (th(ii)-th(ii-1))/(eh1-eh0)
    goto 501
  endif
  goto 500
 501    continue

!===============================================================================
! 1. CALCUL DE L'ENTHALPIE A PARTIR DE LA TEMPERATURE
!===============================================================================

else if ( mode .eq. -1 ) then

  ii = npo

! --- Calcul en Max

  eh =  xesp(ifov )*ehgaze(ifov,ii)                               &
      + xesp(ico  )*ehgaze(ico ,ii)                               &
      + xesp(ih2s )*ehgaze(ih2s,ii)                               &
      + xesp(ihy  )*ehgaze(ihy ,ii)                               &
      + xesp(ihcn )*ehgaze(ihcn,ii)                               &
      + xesp(io2  )*ehgaze(io2 ,ii)                               &
      + xesp(ico2 )*ehgaze(ico2,ii)                               &
      + xesp(ih2o )*ehgaze(ih2o,ii)                               &
      + xesp(iso2 )*ehgaze(iso2,ii)                               &
      + xesp(inh3 )*ehgaze(inh3,ii)                               &
      + xesp(in2  )*ehgaze(in2 ,ii)

  if (tp.gt.th(ii)) goto 601

! Clipping en Min

  ii = 1


! --- Clipping eventuel de TP a TH(1) si EH < EH0

  eh =  xesp(ifov )*ehgaze(ifov,ii)                               &
      + xesp(ico  )*ehgaze(ico ,ii)                               &
      + xesp(ih2s )*ehgaze(ih2s,ii)                               &
      + xesp(ihy  )*ehgaze(ihy ,ii)                               &
      + xesp(ihcn )*ehgaze(ihcn,ii)                               &
      + xesp(io2  )*ehgaze(io2 ,ii)                               &
      + xesp(ico2 )*ehgaze(ico2,ii)                               &
      + xesp(ih2o )*ehgaze(ih2o,ii)                               &
      + xesp(iso2 )*ehgaze(iso2,ii)                               &
      + xesp(inh3 )*ehgaze(inh3,ii)                               &
      + xesp(in2  )*ehgaze(in2 ,ii)

  if (tp.lt.th(ii)) goto 601

! Interpolation dans la table

  ii = 1
 600    continue

  ii = ii + 1
  if ( tp .le. th(ii) ) then

! --- Calcul de l'enthalpie de l'espece gazeuse TH(II-1)

    eh0 = xesp(ifov )*ehgaze(ifov,ii-1)                           &
        + xesp(ico  )*ehgaze(ico ,ii-1)                           &
        + xesp(ih2s )*ehgaze(ih2s,ii-1)                           &
        + xesp(ihy  )*ehgaze(ihy ,ii-1)                           &
        + xesp(ihcn )*ehgaze(ihcn,ii-1)                           &
        + xesp(io2  )*ehgaze(io2 ,ii-1)                           &
        + xesp(ico2 )*ehgaze(ico2,ii-1)                           &
        + xesp(ih2o )*ehgaze(ih2o,ii-1)                           &
        + xesp(iso2 )*ehgaze(iso2,ii-1)                           &
        + xesp(inh3 )*ehgaze(inh3,ii-1)                           &
        + xesp(in2  )*ehgaze(in2 ,ii-1)

! --- Calcul de l'enthalpie de l'espece gazeuse TH(I)

    eh1 =  xesp(ifov )*ehgaze(ifov,ii)                            &
        + xesp(ico  )*ehgaze(ico ,ii)                             &
        + xesp(ih2s )*ehgaze(ih2s,ii)                             &
        + xesp(ihy  )*ehgaze(ihy ,ii)                             &
        + xesp(ihcn )*ehgaze(ihcn,ii)                             &
        + xesp(io2  )*ehgaze(io2 ,ii)                             &
        + xesp(ico2 )*ehgaze(ico2,ii)                             &
        + xesp(ih2o )*ehgaze(ih2o,ii)                             &
        + xesp(iso2 )*ehgaze(iso2,ii)                             &
        + xesp(inh3 )*ehgaze(inh3,ii)                             &
        + xesp(in2  )*ehgaze(in2 ,ii)


    eh =  eh0                                                     &
         +(eh1-eh0)*(tp-th(ii-1))/(th(ii)-th(ii-1))
    goto 601
  endif
  goto 600

 601    continue

else
  write(nfecra,1000) mode
  call csexit(1)
endif

!--------
! Formats
!--------

 1000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ERREUR DANS CS_FUEL_HTCONVERS1         ',/,&
'@    =========                                               ',/,&
'@    VALEUR INCORRECTE DE L''ARGUMENT MODE                   ',/,&
'@    CE DOIT ETRE UN ENTIER EGAL A 1 OU -1                   ',/,&
'@    IL VAUT ICI ',I10                                        ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

!----
! End
!----

return
end subroutine

