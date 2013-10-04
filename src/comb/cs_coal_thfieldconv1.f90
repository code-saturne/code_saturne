!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2013 EDF S.A.
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

subroutine cs_coal_thfieldconv1 &
!==============================
 ( ncelet , ncel   ,                                              &
   eh     , x2     , rtp    ,                                     &
   fuel1  , fuel2  , fuel3  , fuel4 , fuel5 , fuel6 , fuel7 ,     &
   oxyd   , prod1  , prod2  , prod3 , xiner ,                     &
   tp     )
!===============================================================================
! FONCTION :
! --------
! CALCUL DE LA TEMPERATURE DU GAZ
!  EN FONCTION DE L'ENTHALPIE DU GAZ ET DES CONCENTRATIONS
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! ncelet           ! i  ! <-- ! number of extended (real + ghost) cells        !
! ncel             ! i  ! <-- ! number of cells                                !
! ntbmci           ! e  ! <-- ! taille du macro tableau mc entiers             !
! ntbmcr           ! e  ! <-- ! taille du macro tableau mc reels               !
! eh               ! tr ! <-- ! enthalpie du gaz                               !
!                  !    !     ! (j/kg de melange gazeux)                       !
! fuel1            ! tr ! <-- ! fraction massique chx1                         !
! fuel2            ! tr ! <-- ! fraction massique chx2                         !
! fuel3            ! tr ! <-- ! fraction massique co                           !
! oxyd             ! tr ! <-- ! fraction massique o2                           !
! prod1            ! tr ! <-- ! fraction massique co2                          !
! prod2            ! tr ! <-- ! fraction massique h2o                          !
! xiner            ! tr ! <-- ! fraction massique n2                           !
! tp               ! tr ! --> ! temperature du gaz (kelvin)                    !
! tbmci            ! tr ! <-- ! macro tableau entier mc travail                !
! tbmcr            ! tr ! <-- ! macro tableau reel   mc travail                !
! eh0              ! tr ! <-- ! tableau reel de travail                        !
! eh1              ! tr ! <-- ! tableau reel de travail                        !
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
use ppcpfu

!===============================================================================

implicit none

! Arguments

integer          ncelet , ncel

double precision eh(ncelet)   , x2(ncelet)   , rtp(ncelet,*)
double precision fuel1(ncelet), fuel2(ncelet), fuel3(ncelet)
double precision fuel4(ncelet), fuel5(ncelet), fuel6(ncelet) , fuel7(ncelet)
double precision oxyd(ncelet), xiner(ncelet)
double precision prod1(ncelet),prod2(ncelet),prod3(ncelet)
double precision tp(ncelet)

! Local variables

integer          i, icel, icha

double precision ychx10 , ychx20 , ehchx1 , ehchx2
double precision den1   , den2 , eh0 , eh1

integer          iok
double precision , dimension ( : , : )     , allocatable :: f1mc,f2mc

!===============================================================================
!
!===============================================================================
! Deallocation dynamic arrays
!----
allocate(f1mc(1:ncel,1:ncharb),f2mc(1:ncel,1:ncharb),stat=iok)
!----
if ( iok > 0  ) then
  write(nfecra,*) ' Memory allocation error inside: '
  write(nfecra,*) '    cs_coal_thfieldconv1         '
  call csexit(1)
endif
!===============================================================================
!
do icel = 1, ncel
  do icha = 1, ncharb
    f1mc(icel,icha) = rtp(icel,isca(if1m(icha))) /(1.d0-x2(icel))
    f2mc(icel,icha) = rtp(icel,isca(if2m(icha))) /(1.d0-x2(icel))
  enddo
enddo
!
i = npo-1
do icel = 1, ncel

! --- Calcul de l'enthalpie de l'espece gazeuse CHx1m
!                                            et CHx2m a TH(NPO)
  ehchx1 = zero
  ehchx2 = zero
  ychx10 = zero
  ychx20 = zero
  do icha = 1, ncharb
!
    den1   = 1.d0                                                  &
         / ( a1(icha)*wmole(ichx1c(icha))                          &
            +b1(icha)*wmole(ico)                                   &
            +c1(icha)*wmole(ih2o)                                  &
            +d1(icha)*wmole(ih2s)                                  &
            +e1(icha)*wmole(ihcn)                                  &
            +f1(icha)*wmole(inh3) )
    ychx10 = ychx10                                                &
            +den1*(f1mc(icel,icha)*a1(icha)*wmole(ichx1c(icha)) )
    ehchx1 = ehchx1                                                &
            +den1*( ehgaze(ichx1c(icha),i+1)                       &
                   *f1mc(icel,icha)*a1(icha)*wmole(ichx1c(icha)) )
    den2   = 1.d0                                                  &
         / ( a2(icha)*wmole(ichx2c(icha))                          &
            +b2(icha)*wmole(ico)                                   &
            +c2(icha)*wmole(ih2o)                                  &
            +d2(icha)*wmole(ih2s)                                  &
            +e2(icha)*wmole(ihcn)                                  &
            +f2(icha)*wmole(inh3) )
    ychx20 = ychx20                                                &
            +den2 *(f2mc(icel,icha)*a2(icha)*wmole(ichx2c(icha)) )
    ehchx2 = ehchx2 + den2 *                                       &
         ( ehgaze(ichx2c(icha),i+1)                                &
          *f2mc(icel,icha)*a2(icha)*wmole(ichx2c(icha)) )
  enddo
  if ( ychx10.gt.epzero ) then
    ehchx1 = ehchx1 / ychx10
  else
    ehchx1 = ehgaze(ichx1,i+1)
  endif
  if ( ychx20.gt.epzero ) then
    ehchx2 = ehchx2 / ychx20
  else
    ehchx2 = ehgaze(ichx2,i+1)
  endif

! --- Clipping eventuel de TP a TH(NPO) si EH > EH1

  eh1 = fuel1(icel)*ehchx1                                  &
       +fuel2(icel)*ehchx2                                  &
       +fuel3(icel)*ehgaze(ico ,i+1)                        &
       +fuel4(icel)*ehgaze(ih2s,i+1)                        &
       +fuel5(icel)*ehgaze(ihy ,i+1)                        &
       +fuel6(icel)*ehgaze(ihcn,i+1)                        &
       +fuel7(icel)*ehgaze(inh3,i+1)                        &
       +oxyd(icel) *ehgaze(io2 ,i+1)                        &
       +prod1(icel)*ehgaze(ico2,i+1)                        &
       +prod2(icel)*ehgaze(ih2o,i+1)                        &
       +prod3(icel)*ehgaze(iso2,i+1)                        &
       +xiner(icel)*ehgaze(in2 ,i+1)

  if ( eh(icel).ge.eh1 ) tp(icel)= th(i+1)
enddo

i = 1
do icel = 1, ncel

! --- Calcul de l'enthalpie de l'espece gazeuse CHx1m
!                                            et CHx2m a TH(1)
  ehchx1 = zero
  ehchx2 = zero
  ychx10 = zero
  ychx20 = zero
  do icha = 1, ncharb
    den1   = 1.d0                                                    &
         / ( a1(icha)*wmole(ichx1c(icha))                            &
            +b1(icha)*wmole(ico)                                     &
            +c1(icha)*wmole(ih2o)                                    &
            +d1(icha)*wmole(ih2s)                                    &
            +e1(icha)*wmole(ihcn)                                    &
            +f1(icha)*wmole(inh3) )
    ychx10 = ychx10                                                  &
            +den1*(f1mc(icel,icha)*a1(icha)*wmole(ichx1c(icha)) )
    ehchx1 = ehchx1                                                  &
            +den1*( ehgaze(ichx1c(icha),i)                           &
                   *f1mc(icel,icha)*a1(icha)*wmole(ichx1c(icha)) )
    den2   = 1.d0                                                     &
         / ( a2(icha)*wmole(ichx2c(icha))                             &
            +b2(icha)*wmole(ico)                                      &
            +c2(icha)*wmole(ih2o)                                     &
            +d2(icha)*wmole(ih2s)                                     &
            +e2(icha)*wmole(ihcn)                                     &
            +f2(icha)*wmole(inh3) )
    ychx20 = ychx20                                                   &
            +den2*(f2mc(icel,icha)*a2(icha)*wmole(ichx2c(icha)) )
    ehchx2 = ehchx2                                                   &
            +den2*( ehgaze(ichx2c(icha),i)                            &
                   *f2mc(icel,icha)*a2(icha)*wmole(ichx2c(icha)) )
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

  eh0  = fuel1(icel)*ehchx1                                  &
        +fuel2(icel)*ehchx2                                  &
        +fuel3(icel)*ehgaze(ico ,i)                          &
        +fuel4(icel)*ehgaze(ih2s,i)                          &
        +fuel5(icel)*ehgaze(ihy ,i)                          &
        +fuel6(icel)*ehgaze(ihcn,i)                          &
        +fuel7(icel)*ehgaze(inh3,i)                          &
        +oxyd(icel) *ehgaze(io2 ,i)                          &
        +prod1(icel)*ehgaze(ico2,i)                          &
        +prod2(icel)*ehgaze(ih2o,i)                          &
        +prod3(icel)*ehgaze(iso2,i)                          &
        +xiner(icel)*ehgaze(in2 ,i)

  if ( eh(icel) .le. eh0 ) tp(icel)= th(i)

enddo


do i = 1, npo-1
  do icel = 1, ncel

! --- Calcul de l'enthalpie de l'espece gazeuse CHx1m
!                                            et CHx2m pour TH(I)
    ehchx1 = zero
    ehchx2 = zero
    ychx10 = zero
    ychx20 = zero
    do icha = 1, ncharb
      den1   = 1.d0                                                 &
             / ( a1(icha)*wmole(ichx1c(icha))                       &
                +b1(icha)*wmole(ico)                                &
                +c1(icha)*wmole(ih2o)                               &
                +d1(icha)*wmole(ih2s)                               &
                +e1(icha)*wmole(ihcn)                               &
                +f1(icha)*wmole(inh3) )
      ychx10 = ychx10                                               &
              +den1*f1mc(icel,icha)*a1(icha)*wmole(ichx1c(icha))
      ehchx1 = ehchx1                                               &
              +den1*( ehgaze(ichx1c(icha),i)                        &
                     *f1mc(icel,icha)*a1(icha)*wmole(ichx1c(icha)) )
      den2   = 1.d0                                                 &
             / ( a2(icha)*wmole(ichx2c(icha))                       &
                +b2(icha)*wmole(ico)                                &
                +c2(icha)*wmole(ih2o)                               &
                +d2(icha)*wmole(ih2s)                               &
                +e2(icha)*wmole(ihcn)                               &
                +f2(icha)*wmole(inh3) )
      ychx20 = ychx20                                               &
              +den2*(f2mc(icel,icha)*a2(icha)*wmole(ichx2c(icha)) )
      ehchx2 = ehchx2                                               &
              +den2*( ehgaze(ichx2c(icha),i)                        &
                     *f2mc(icel,icha)*a2(icha)*wmole(ichx2c(icha)) )
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
    eh0 = fuel1(icel)*ehchx1                                  &
         +fuel2(icel)*ehchx2                                  &
         +fuel3(icel)*ehgaze(ico ,i)                          &
         +fuel4(icel)*ehgaze(ih2s,i)                          &
         +fuel5(icel)*ehgaze(ihy ,i)                          &
         +fuel6(icel)*ehgaze(ihcn,i)                          &
         +fuel7(icel)*ehgaze(inh3,i)                          &
         +oxyd(icel) *ehgaze(io2 ,i)                          &
         +prod1(icel)*ehgaze(ico2,i)                          &
         +prod2(icel)*ehgaze(ih2o,i)                          &
         +prod3(icel)*ehgaze(iso2,i)                          &
         +xiner(icel)*ehgaze(in2 ,i)

! --- Calcul de l'enthalpie de l'espece gazeuse CHx1m
!                                            et CHx2m pour TH(I+1)
    ehchx1 = zero
    ehchx2 = zero
    ychx10 = zero
    ychx20 = zero
    do icha = 1, ncharb
      den1   = 1.d0                                                   &
           / ( a1(icha)*wmole(ichx1c(icha))                           &
              +b1(icha)*wmole(ico)                                    &
              +c1(icha)*wmole(ih2o)                                   &
              +d1(icha)*wmole(ih2s)                                   &
              +e1(icha)*wmole(ihcn)                                   &
              +f1(icha)*wmole(inh3) )
      ychx10 = ychx10                                                 &
              +den1*f1mc(icel,icha)*a1(icha)*wmole(ichx1c(icha))
      ehchx1 = ehchx1                                                 &
              +den1*( ehgaze(ichx1c(icha),i+1)                        &
                     *f1mc(icel,icha)*a1(icha)*wmole(ichx1c(icha)) )
      den2   = 1.d0                                                   &
             / ( a2(icha)*wmole(ichx2c(icha))                         &
                +b2(icha)*wmole(ico)                                  &
                +c2(icha)*wmole(ih2o)                                 &
                +d2(icha)*wmole(ih2s)                                 &
                +e2(icha)*wmole(ihcn)                                 &
                +f2(icha)*wmole(inh3) )
      ychx20 = ychx20                                                 &
              +den2*f2mc(icel,icha)*a2(icha)*wmole(ichx2c(icha))
      ehchx2 = ehchx2                                                 &
              +den2*( ehgaze(ichx2c(icha),i+1)                        &
                     *f2mc(icel,icha)*a2(icha)*wmole(ichx2c(icha)) )
    enddo
    if ( ychx10.gt.epzero ) then
      ehchx1 = ehchx1 / ychx10
    else
      ehchx1 = ehgaze(ichx1,i+1)
    endif
    if ( ychx20.gt.epzero ) then
      ehchx2 = ehchx2 / ychx20
    else
      ehchx2 = ehgaze(ichx2,i+1)
    endif

    eh1 = fuel1(icel)*ehchx1                                    &
         +fuel2(icel)*ehchx2                                    &
         +fuel3(icel)*ehgaze(ico ,i+1)                          &
         +fuel4(icel)*ehgaze(ih2s,i+1)                          &
         +fuel5(icel)*ehgaze(ihy ,i+1)                          &
         +fuel6(icel)*ehgaze(ihcn,i+1)                          &
         +fuel7(icel)*ehgaze(inh3,i+1)                          &
         +oxyd(icel) *ehgaze(io2 ,i+1)                          &
         +prod1(icel)*ehgaze(ico2,i+1)                          &
         +prod2(icel)*ehgaze(ih2o,i+1)                          &
         +prod3(icel)*ehgaze(iso2,i+1)                          &
         +xiner(icel)*ehgaze(in2 ,i+1)

    if ( eh(icel).ge.eh0 .and. eh(icel).le.eh1 ) then
      tp(icel)= th(i) + (eh(icel)-eh0) *                        &
                        (th(i+1)-th(i))/(eh1-eh0)
    endif
  enddo
enddo

!----
! End
!----

return
end subroutine
