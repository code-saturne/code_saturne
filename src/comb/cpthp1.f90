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

subroutine cpthp1 &
!================

 ( mode  , eh     , xesp   , f1mc   , f2mc   ,                    &
   tp     )

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
! f1mc             ! tr ! <-- ! f1 moyen                                       !
! f2mc             ! tr ! <-- ! f2 moyen                                       !
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
use ppincl

!===============================================================================

implicit none

! Arguments

integer          mode
double precision eh,tp
double precision xesp(ngazem)
double precision f1mc(*),f2mc(*) ! size ncharm or ncharm2

! Local variables

integer          i , icha

double precision ychx10 , ychx20 , ehchx1 , ehchx2
double precision den1   , den2
double precision eh0 , eh1

!===============================================================================
!===============================================================================
! 1. CALCUL DE LA TEMPERATURE A PARTIR DE l'ENTHALPIE
!===============================================================================

if ( mode .eq. 1 ) then

  i = npo

! --- Calcul de l'enthalpie de l'espece gazeuse CHx1m
!                                            et CHx2m a TH(NPO)
  ehchx1 = zero
  ehchx2 = zero
  ychx10 = zero
  ychx20 = zero
  do icha = 1, ncharb
    den1   = 1.d0                                                 &
         / ( a1(icha)*wmole(ichx1c(icha))+b1(icha)*wmole(ico))
    ychx10 = ychx10 + den1 *                                      &
         ( f1mc(icha)*a1(icha)*wmole(ichx1c(icha)) )
    ehchx1 = ehchx1 + den1 *                                      &
         ( ehgaze(ichx1c(icha),i)*                                &
         f1mc(icha)*a1(icha)*wmole(ichx1c(icha)) )
    den2   = 1.d0                                                 &
         / ( a2(icha)*wmole(ichx2c(icha))+b2(icha)*wmole(ico))
    ychx20 = ychx20 + den2 *                                      &
         ( f2mc(icha)*a2(icha)*wmole(ichx2c(icha)) )
    ehchx2 = ehchx2 + den2 *                                      &
         ( ehgaze(ichx2c(icha),i)*                                &
         f2mc(icha)*a2(icha)*wmole(ichx2c(icha)) )
  enddo
  if ( ychx10.gt.epzero ) then
    ehchx1 = ehchx1 / ychx10
  else
    ehchx1 = ehgaze(ichx1,i)
  endif
  if ( ychx20.gt.epzero ) then
    ehchx2 = ehchx2 / ychx20
  else
    ehchx2 = ehgaze(ichx2,i)
  endif

! --- Clipping eventuel de TP a TH(NPO) si EH > EH1

  eh1 = xesp(ichx1)*ehchx1                                        &
      + xesp(ichx2)*ehchx2                                        &
      + xesp(ico  )*ehgaze(ico ,i)                                &
      + xesp(io2  )*ehgaze(io2 ,i)                                &
      + xesp(ico2 )*ehgaze(ico2,i)                                &
      + xesp(ih2o )*ehgaze(ih2o,i)                                &
      + xesp(in2  )*ehgaze(in2 ,i)

  if ( eh .ge. eh1 ) then
    tp = th(i)
    goto 501
  endif

  i = 1

! --- Calcul de l'enthalpie de l'espece gazeuse CHx1m
!                                            et CHx2m a TH(1)
  ehchx1 = zero
  ehchx2 = zero
  ychx10 = zero
  ychx20 = zero
  do icha = 1, ncharb
    den1   = 1.d0                                                 &
         / ( a1(icha)*wmole(ichx1c(icha))+b1(icha)*wmole(ico))
    ychx10 = ychx10 + den1 *                                      &
     ( f1mc(icha)*a1(icha)*wmole(ichx1c(icha)) )
    ehchx1 = ehchx1 + den1 *                                      &
     ( ehgaze(ichx1c(icha),i)*                                    &
       f1mc(icha)*a1(icha)*wmole(ichx1c(icha)) )
    den2   = 1.d0                                                 &
         / ( a2(icha)*wmole(ichx2c(icha))+b2(icha)*wmole(ico))
    ychx20 = ychx20 + den2 *                                      &
     ( f2mc(icha)*a2(icha)*wmole(ichx2c(icha)) )
    ehchx2 = ehchx2 + den2 *                                      &
     ( ehgaze(ichx2c(icha),i)*                                    &
       f2mc(icha)*a2(icha)*wmole(ichx2c(icha)) )
  enddo
  if ( ychx10.gt.epzero ) then
    ehchx1 = ehchx1 / ychx10
  else
    ehchx1 = ehgaze(ichx1,i)
  endif
  if ( ychx20.gt.epzero ) then
    ehchx2 = ehchx2 / ychx20
  else
    ehchx2 = ehgaze(ichx2,i)
  endif

! --- Clipping eventuel de TP a TH(1) si EH < EH0

  eh0 = xesp(ichx1)*ehchx1                                        &
      + xesp(ichx2)*ehchx2                                        &
      + xesp(ico  )*ehgaze(ico ,i)                                &
      + xesp(io2  )*ehgaze(io2 ,i)                                &
      + xesp(ico2 )*ehgaze(ico2,i)                                &
      + xesp(ih2o )*ehgaze(ih2o,i)                                &
      + xesp(in2  )*ehgaze(in2 ,i)

  if ( eh .le. eh0 ) then
    tp= th(i)
    goto 501
  endif


 500    continue
  i = i + 1

! --- Calcul de l'enthalpie de l'espece gazeuse CHx1m
!                                            et CHx2m pour TH(I-1)
  ehchx1 = zero
  ehchx2 = zero
  ychx10 = zero
  ychx20 = zero
  do icha = 1, ncharb
    den1   = 1.d0                                                 &
           / ( a1(icha)*wmole(ichx1c(icha))+b1(icha)*wmole(ico))
    ychx10 = ychx10 + den1 *                                      &
      ( f1mc(icha)*a1(icha)*wmole(ichx1c(icha)) )
    ehchx1 = ehchx1 + den1 *                                      &
      ( ehgaze(ichx1c(icha),i-1)*                                 &
        f1mc(icha)*a1(icha)*wmole(ichx1c(icha)) )
    den2   = 1.d0                                                 &
           / ( a2(icha)*wmole(ichx2c(icha))+b2(icha)*wmole(ico))
    ychx20 = ychx20 + den2 *                                      &
       ( f2mc(icha)*a2(icha)*wmole(ichx2c(icha)) )
    ehchx2 = ehchx2 + den2 *                                      &
       ( ehgaze(ichx2c(icha),i-1)*                                &
         f2mc(icha)*a2(icha)*wmole(ichx2c(icha)) )
  enddo
  if ( ychx10.gt.epzero ) then
    ehchx1 = ehchx1 / ychx10
  else
    ehchx1 = ehgaze(ichx1,i-1)
  endif
  if ( ychx20.gt.epzero ) then
    ehchx2 = ehchx2 / ychx20
  else
    ehchx2 = ehgaze(ichx2,i-1)
  endif
  eh0 = xesp(ichx1)*ehchx1                                        &
      + xesp(ichx2)*ehchx2                                        &
      + xesp(ico  )*ehgaze(ico ,i-1)                              &
      + xesp(io2  )*ehgaze(io2 ,i-1)                              &
      + xesp(ico2 )*ehgaze(ico2,i-1)                              &
      + xesp(ih2o )*ehgaze(ih2o,i-1)                              &
      + xesp(in2  )*ehgaze(in2 ,i-1)

! --- Calcul de l'enthalpie de l'espece gazeuse CHx1m
!                                            et CHx2m pour TH(I)
  ehchx1 = zero
  ehchx2 = zero
  ychx10 = zero
  ychx20 = zero
  do icha = 1, ncharb
    den1   = 1.d0                                                 &
           / ( a1(icha)*wmole(ichx1c(icha))+b1(icha)*wmole(ico))
    ychx10 = ychx10 + den1 *                                      &
      ( f1mc(icha)*a1(icha)*wmole(ichx1c(icha)) )
    ehchx1 = ehchx1 + den1 *                                      &
      ( ehgaze(ichx1c(icha),i)*                                   &
        f1mc(icha)*a1(icha)*wmole(ichx1c(icha)) )
    den2   = 1.d0                                                 &
           / ( a2(icha)*wmole(ichx2c(icha))+b2(icha)*wmole(ico))
    ychx20 = ychx20 + den2 *                                      &
       ( f2mc(icha)*a2(icha)*wmole(ichx2c(icha)) )
    ehchx2 = ehchx2 + den2 *                                      &
       ( ehgaze(ichx2c(icha),i)*                                  &
         f2mc(icha)*a2(icha)*wmole(ichx2c(icha)) )
  enddo
  if ( ychx10.gt.epzero ) then
    ehchx1 = ehchx1 / ychx10
  else
    ehchx1 = ehgaze(ichx1,i)
  endif
  if ( ychx20.gt.epzero ) then
    ehchx2 = ehchx2 / ychx20
  else
    ehchx2 = ehgaze(ichx2,i)
  endif

  eh1 = xesp(ichx1)*ehchx1                                        &
      + xesp(ichx2)*ehchx2                                        &
      + xesp(ico  )*ehgaze(ico ,i)                                &
      + xesp(io2  )*ehgaze(io2 ,i)                                &
      + xesp(ico2 )*ehgaze(ico2,i)                                &
      + xesp(ih2o )*ehgaze(ih2o,i)                                &
      + xesp(in2  )*ehgaze(in2 ,i)

  if ( eh .ge. eh0  .and. eh .le. eh1  ) then
    tp = th(i-1) + (eh-eh0) *                                     &
                 (th(i)-th(i-1))/(eh1-eh0)
    goto 501
  endif
  goto 500
 501    continue

!===============================================================================
! 1. CALCUL DE L'ENTHALPIE A PARTIR DE LA TEMPERATURE
!===============================================================================

else if ( mode .eq. -1 ) then

  i = npo

! --- Calcul en Max

  if ( tp .ge. th(i) ) then
    ehchx1 = zero
    ehchx2 = zero
    ychx10 = zero
    ychx20 = zero
    do icha = 1, ncharb
      den1   = 1.d0                                               &
           / ( a1(icha)*wmole(ichx1c(icha))+b1(icha)*wmole(ico))
      ychx10 = ychx10 + den1 *                                    &
           ( f1mc(icha)*a1(icha)*wmole(ichx1c(icha)) )
      ehchx1 = ehchx1 + den1 *                                    &
           ( ehgaze(ichx1c(icha),i)*                              &
           f1mc(icha)*a1(icha)*wmole(ichx1c(icha)) )
      den2   = 1.d0                                               &
           / ( a2(icha)*wmole(ichx2c(icha))+b2(icha)*wmole(ico))
      ychx20 = ychx20 + den2 *                                    &
           ( f2mc(icha)*a2(icha)*wmole(ichx2c(icha)) )
      ehchx2 = ehchx2 + den2 *                                    &
           ( ehgaze(ichx2c(icha),i)*                              &
           f2mc(icha)*a2(icha)*wmole(ichx2c(icha)) )
    enddo
    if ( ychx10.gt.epzero ) then
      ehchx1 = ehchx1 / ychx10
    else
      ehchx1 = ehgaze(ichx1,i)
    endif
    if ( ychx20.gt.epzero ) then
      ehchx2 = ehchx2 / ychx20
    else
      ehchx2 = ehgaze(ichx2,i)
    endif
    eh = xesp(ichx1)*ehchx1                                       &
       + xesp(ichx2)*ehchx2                                       &
       + xesp(ico  )*ehgaze(ico ,i)                               &
       + xesp(io2  )*ehgaze(io2 ,i)                               &
       + xesp(ico2 )*ehgaze(ico2,i)                               &
       + xesp(ih2o )*ehgaze(ih2o,i)                               &
       + xesp(in2  )*ehgaze(in2 ,i)
    goto 601
  endif

! Clipping en Min

  i = 1

  if ( tp .le. th(i) ) then
    ehchx1 = zero
    ehchx2 = zero
    ychx10 = zero
    ychx20 = zero
    do icha = 1, ncharb
      den1   = 1.d0                                               &
           / ( a1(icha)*wmole(ichx1c(icha))+b1(icha)*wmole(ico))
      ychx10 = ychx10 + den1 *                                    &
           ( f1mc(icha)*a1(icha)*wmole(ichx1c(icha)) )
      ehchx1 = ehchx1 + den1 *                                    &
           ( ehgaze(ichx1c(icha),i)*                              &
           f1mc(icha)*a1(icha)*wmole(ichx1c(icha)) )
      den2   = 1.d0                                               &
           / ( a2(icha)*wmole(ichx2c(icha))+b2(icha)*wmole(ico))
      ychx20 = ychx20 + den2 *                                    &
           ( f2mc(icha)*a2(icha)*wmole(ichx2c(icha)) )
      ehchx2 = ehchx2 + den2 *                                    &
           ( ehgaze(ichx2c(icha),i)*                              &
           f2mc(icha)*a2(icha)*wmole(ichx2c(icha)) )
    enddo
    if ( ychx10.gt.epzero ) then
      ehchx1 = ehchx1 / ychx10
    else
      ehchx1 = ehgaze(ichx1,i)
    endif
    if ( ychx20.gt.epzero ) then
      ehchx2 = ehchx2 / ychx20
    else
      ehchx2 = ehgaze(ichx2,i)
    endif

! --- Clipping eventuel de TP a TH(1) si EH < EH0

    eh = xesp(ichx1)*ehchx1                                       &
       + xesp(ichx2)*ehchx2                                       &
       + xesp(ico  )*ehgaze(ico ,i)                               &
       + xesp(io2  )*ehgaze(io2 ,i)                               &
       + xesp(ico2 )*ehgaze(ico2,i)                               &
       + xesp(ih2o )*ehgaze(ih2o,i)                               &
       + xesp(in2  )*ehgaze(in2 ,i)
    goto 601
  endif

! Interpolation dans la table

  i = 1
 600    continue

  i = i + 1
  if ( tp .le. th(i) ) then

! --- Calcul de l'enthalpie de l'espece gazeuse TH(I-1)

    ehchx1 = zero
    ehchx2 = zero
    ychx10 = zero
    ychx20 = zero
    do icha = 1, ncharb
      den1   = 1.d0                                               &
             / ( a1(icha)*wmole(ichx1c(icha))+b1(icha)*wmole(ico))
      ychx10 = ychx10 + den1 *                                    &
        ( f1mc(icha)*a1(icha)*wmole(ichx1c(icha)) )
      ehchx1 = ehchx1 + den1 *                                    &
        ( ehgaze(ichx1c(icha),i-1)*                               &
          f1mc(icha)*a1(icha)*wmole(ichx1c(icha)) )
      den2   = 1.d0                                               &
             / ( a2(icha)*wmole(ichx2c(icha))+b2(icha)*wmole(ico))
      ychx20 = ychx20 + den2 *                                    &
         ( f2mc(icha)*a2(icha)*wmole(ichx2c(icha)) )
      ehchx2 = ehchx2 + den2 *                                    &
         ( ehgaze(ichx2c(icha),i-1)*                              &
           f2mc(icha)*a2(icha)*wmole(ichx2c(icha)) )
    enddo
    if ( ychx10.gt.epzero ) then
      ehchx1 = ehchx1 / ychx10
    else
      ehchx1 = ehgaze(ichx1,i-1)
    endif
    if ( ychx20.gt.epzero ) then
      ehchx2 = ehchx2 / ychx20
    else
      ehchx2 = ehgaze(ichx2,i-1)
    endif

    eh0 = xesp(ichx1)*ehchx1                                      &
        + xesp(ichx2)*ehchx2                                      &
        + xesp(ico  )*ehgaze(ico ,i-1)                            &
        + xesp(io2  )*ehgaze(io2 ,i-1)                            &
        + xesp(ico2 )*ehgaze(ico2,i-1)                            &
        + xesp(ih2o )*ehgaze(ih2o,i-1)                            &
        + xesp(in2  )*ehgaze(in2 ,i-1)

! --- Calcul de l'enthalpie de l'espece gazeuse TH(I)

    ehchx1 = zero
    ehchx2 = zero
    ychx10 = zero
    ychx20 = zero
    do icha = 1, ncharb
      den1   = 1.d0                                               &
             / ( a1(icha)*wmole(ichx1c(icha))+b1(icha)*wmole(ico))
      ychx10 = ychx10 + den1 *                                    &
        ( f1mc(icha)*a1(icha)*wmole(ichx1c(icha)) )
      ehchx1 = ehchx1 + den1 *                                    &
        ( ehgaze(ichx1c(icha),i)*                                 &
          f1mc(icha)*a1(icha)*wmole(ichx1c(icha)) )
      den2   = 1.d0                                               &
             / ( a2(icha)*wmole(ichx2c(icha))+b2(icha)*wmole(ico))
      ychx20 = ychx20 + den2 *                                    &
         ( f2mc(icha)*a2(icha)*wmole(ichx2c(icha)) )
      ehchx2 = ehchx2 + den2 *                                    &
         ( ehgaze(ichx2c(icha),i)*                                &
           f2mc(icha)*a2(icha)*wmole(ichx2c(icha)) )
    enddo
    if ( ychx10.gt.epzero ) then
      ehchx1 = ehchx1 / ychx10
    else
      ehchx1 = ehgaze(ichx1,i)
    endif
    if ( ychx20.gt.epzero ) then
      ehchx2 = ehchx2 / ychx20
    else
      ehchx2 = ehgaze(ichx2,i)
    endif

    eh1 = xesp(ichx1)*ehchx1                                      &
        + xesp(ichx2)*ehchx2                                      &
        + xesp(ico  )*ehgaze(ico ,i)                              &
        + xesp(io2  )*ehgaze(io2 ,i)                              &
        + xesp(ico2 )*ehgaze(ico2,i)                              &
        + xesp(ih2o )*ehgaze(ih2o,i)                              &
        + xesp(in2  )*ehgaze(in2 ,i)

    eh =  eh0                                                     &
           + (eh1-eh0)*(tp-th(i-1))/(th(i)-th(i-1))
    goto 601
  endif
  goto 600

 601    continue

else
  write(nfecra,1000) mode
  call csexit(1)
endif

!--------
! FORMATS
!--------

 1000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ERREUR DANS CPTHP1                          ',/,&
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
! FIN
!----

return
end subroutine
