!-------------------------------------------------------------------------------

!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2009 EDF S.A., France

!     contact: saturne-support@edf.fr

!     The Code_Saturne Kernel is free software; you can redistribute it
!     and/or modify it under the terms of the GNU General Public License
!     as published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.

!     The Code_Saturne Kernel is distributed in the hope that it will be
!     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
!     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.

!     You should have received a copy of the GNU General Public License
!     along with the Code_Saturne Kernel; if not, write to the
!     Free Software Foundation, Inc.,
!     51 Franklin St, Fifth Floor,
!     Boston, MA  02110-1301  USA

!-------------------------------------------------------------------------------

subroutine cptehm &
!================

 ( ncelet , ncel   , ntbmci , ntbmcr ,                            &
   eh     , rtp    ,                                              &
   fuel1  , fuel2  , fuel3  , oxyd   , prod1  , prod2  ,          &
   xiner   ,                                                      &
   tp     ,                                                       &
   tbmci  , tbmcr  ,                                              &
   eh0    , eh1    , x1     )

!===============================================================================
! FONCTION :
! --------
! CALCUL DE LA TEMPERATURE DU MELANGE
!  EN FONCTION DE L'ENTHALPIE DU MELANGE ET DES CONCENTRATIONS

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
! rtp              ! tr ! <-- ! variables de calcul au centre des              !
! (ncelet,*)       !    !     !    cellules (instant courant)                  !
! fuel1            ! tr ! <-- ! fraction massique chx1                         !
! fuel2            ! tr ! <-- ! fraction massique chx2                         !
! fuel3            ! tr ! <-- ! fraction massique co                           !
! oxyd             ! tr ! <-- ! fraction massique o2                           !
! prod1            ! tr ! <-- ! fraction massique co2                          !
! prod2            ! tr ! <-- ! fraction massique h2o                          !
! xiner            ! tr ! <-- ! fraction massique n2                           !
! tp               ! tr ! --> ! temperature du melange (kelvin)                !
! tbmci            ! tr ! <-- ! macro tableau entier mc travail                !
! tbmcr            ! tr ! <-- ! macro tableau reel   mc travail                !
! eh0              ! tr ! <-- ! tableau reel de travail                        !
! eh1              ! tr ! <-- ! tableau reel de travail                        !
! x1               ! tr ! <-- ! tableau reel de travail                        !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHAMNUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

implicit none

!==============================================================================
! Common blocks
!==============================================================================

include "paramx.f90"
include "numvar.f90"
include "optcal.f90"
include "cstphy.f90"
include "cstnum.f90"
include "entsor.f90"
include "ppppar.f90"
include "ppthch.f90"
include "coincl.f90"
include "cpincl.f90"
include "ppincl.f90"

!===============================================================================

! Arguments

integer          ncelet , ncel
integer          ntbmci , ntbmcr
integer          tbmci(ncelet,ntbmci)

double precision eh(ncelet) , rtp(ncelet,*)
double precision fuel1(ncelet), fuel2(ncelet), fuel3(ncelet)
double precision oxyd(ncelet), xiner(ncelet)
double precision prod1(ncelet),prod2(ncelet)
double precision tp(ncelet)
double precision tbmcr(ncelet,ntbmcr)
double precision eh0(ncelet), eh1(ncelet), x1(ncelet)

! Local variables

integer          i, icel, icha , icla

double precision ychx10 , ychx20 , ehchx1 , ehchx2
double precision den1   , den2
double precision xch , xck , xash , x2 , x2h2

!===============================================================================

!===============================================================================
! 1. CALCULS PRELIMINAIRES
!===============================================================================


! RQ IMPORTANTE : Pour la definition des pointeurs pour TBMC
!                 VOIR PPINI1.F

! --- Calcul de X1

do icel = 1, ncel
  x1(icel) = 1.d0
enddo

do icla = 1, nclacp
  do icel = 1, ncel
    xck  = rtp(icel,isca(ixck(icla)))
    xch  = rtp(icel,isca(ixch(icla)))
    xash = rtp(icel,isca(inp(icla)))*xmash(icla)
    x2  = xch + xck + xash
    x1(icel) = x1(icel) - x2
  enddo
enddo


!===============================================================================
! 2. CALCUL DE LA TEMPERATURE DU MELANGE GAZEUX
!===============================================================================

i = npo-1

do 100 icel = 1, ncel

! --- Calcul de l'enthalpie de l'espece gazeuse CHx1m
!                                            et CHx2m a TH(NPO)
  ehchx1 = zero
  ehchx2 = zero
  ychx10 = zero
  ychx20 = zero
  do icha = 1, ncharb
    den1   = 1.d0                                                 &
           / ( a1(icha)*wmole(ichx1c(icha))                       &
              +b1(icha)*wmole(ico)                                &
              +c1(icha)*wmole(ih2o) )
    ychx10 = ychx10 + den1 *                                      &
       ( tbmcr(icel,if1mc(icha))*a1(icha)*wmole(ichx1c(icha)) )
    ehchx1 = ehchx1 + den1 *                                      &
       ( ehgaze(ichx1c(icha),i+1)*                                &
         tbmcr(icel,if1mc(icha))*a1(icha)*wmole(ichx1c(icha)) )
    den2   = 1.d0                                                 &
           / ( a2(icha)*wmole(ichx2c(icha))                       &
              +b2(icha)*wmole(ico)                                &
              +c2(icha)*wmole(ih2o) )
    ychx20 = ychx20 + den2 *                                      &
       ( tbmcr(icel,if2mc(icha))*a2(icha)*wmole(ichx2c(icha)) )
    ehchx2 = ehchx2 + den2 *                                      &
       ( ehgaze(ichx2c(icha),i+1)*                                &
         tbmcr(icel,if2mc(icha))*a2(icha)*wmole(ichx2c(icha)) )
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

! --- Calcul de Somme (X2.H2) a TH(NPO)

  x2h2 = zero
  do icla = 1, nclacp
    xck  = rtp(icel,isca(ixck(icla)))
    xch  = rtp(icel,isca(ixch(icla)))
    xash = rtp(icel,isca(inp(icla)))*xmash(icla)
    x2h2 = x2h2                                                   &
         + xch  * ehsoli(ich(ichcor(icla) ),i+1)                  &
         + xck  * ehsoli(ick(ichcor(icla) ),i+1)                  &
         + xash * ehsoli(iash(ichcor(icla)),i+1)
  enddo

! --- Clipping eventuel de TP a TH(NPO) si EH > EH1

  eh1(icel) =  (  fuel1(icel)*ehchx1                              &
                + fuel2(icel)*ehchx2                              &
                + fuel3(icel)*ehgaze(ico ,i+1)                    &
                + oxyd(icel) *ehgaze(io2 ,i+1)                    &
                + prod1(icel)*ehgaze(ico2,i+1)                    &
                + prod2(icel)*ehgaze(ih2o,i+1)                    &
                + xiner(icel)*ehgaze(in2 ,i+1) ) * x1(icel)       &
             + x2h2

  if ( eh(icel).ge.eh1(icel) ) tp(icel)= th(i+1)
 100  continue



i = 1
do 200 icel = 1, ncel

! --- Calcul de l'enthalpie de l'espece gazeuse CHx1m
!                                            et CHx2m a TH(1)
  ehchx1 = zero
  ehchx2 = zero
  ychx10 = zero
  ychx20 = zero
  do icha = 1, ncharb
    den1   = 1.d0                                                 &
           / ( a1(icha)*wmole(ichx1c(icha))                       &
              +b1(icha)*wmole(ico)                                &
              +c1(icha)*wmole(ih2o) )
    ychx10 = ychx10 + den1 *                                      &
       ( tbmcr(icel,if1mc(icha))*a1(icha)*wmole(ichx1c(icha)) )
    ehchx1 = ehchx1 + den1 *                                      &
       ( ehgaze(ichx1c(icha),i)*                                  &
         tbmcr(icel,if1mc(icha))*a1(icha)*wmole(ichx1c(icha)) )
    den2   = 1.d0                                                 &
           / ( a2(icha)*wmole(ichx2c(icha))                       &
              +b2(icha)*wmole(ico)                                &
              +c2(icha)*wmole(ih2o) )
    ychx20 = ychx20 + den2 *                                      &
       ( tbmcr(icel,if2mc(icha))*a2(icha)*wmole(ichx2c(icha)) )
    ehchx2 = ehchx2 + den2 *                                      &
       ( ehgaze(ichx2c(icha),i)*                                  &
         tbmcr(icel,if2mc(icha))*a2(icha)*wmole(ichx2c(icha)) )
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

! --- Calcul de Somme (X2.H2) a TH(1)

  x2h2 = zero
  do icla = 1, nclacp
    xck  = rtp(icel,isca(ixck(icla)))
    xch  = rtp(icel,isca(ixch(icla)))
    xash = rtp(icel,isca(inp(icla)))*xmash(icla)
    x2h2 = x2h2                                                   &
         + xch  * ehsoli(ich(ichcor(icla) ),i)                    &
         + xck  * ehsoli(ick(ichcor(icla) ),i)                    &
         + xash * ehsoli(iash(ichcor(icla)),i)
  enddo

! --- Clipping eventuel de TP a TH(1) si EH < EH0

  eh0(icel) = (  fuel1(icel)*ehchx1                               &
               + fuel2(icel)*ehchx2                               &
               + fuel3(icel)*ehgaze(ico ,i)                       &
               + oxyd(icel) *ehgaze(io2 ,i)                       &
               + prod1(icel)*ehgaze(ico2,i)                       &
               + prod2(icel)*ehgaze(ih2o,i)                       &
               + xiner(icel)*ehgaze(in2 ,i) ) * x1(icel)          &
             + x2h2

  if ( eh(icel).le.eh0(icel) ) tp(icel)= th(i)
 200  continue


do i = 1, npo-1
  do icel = 1, ncel

! --- Calcul de l'enthalpie de l'espece gazeuse CHx1m
!                                            et CHx2m pour TH(I)
    ehchx1 = zero
    ehchx2 = zero
    ychx10 = zero
    ychx20 = zero
    do icha = 1, ncharb
      den1   = 1.d0                                               &
             / ( a1(icha)*wmole(ichx1c(icha))                     &
                +b1(icha)*wmole(ico)                              &
                +c1(icha)*wmole(ih2o) )
      ychx10 = ychx10 + den1 *                                    &
        ( tbmcr(icel,if1mc(icha))*a1(icha)*wmole(ichx1c(icha)) )
      ehchx1 = ehchx1 + den1 *                                    &
        ( ehgaze(ichx1c(icha),i)*                                 &
          tbmcr(icel,if1mc(icha))*a1(icha)*wmole(ichx1c(icha)) )
      den2   = 1.d0                                               &
             / ( a2(icha)*wmole(ichx2c(icha))                     &
                +b2(icha)*wmole(ico)                              &
                +c2(icha)*wmole(ih2o) )
      ychx20 = ychx20 + den2 *                                    &
         ( tbmcr(icel,if2mc(icha))*a2(icha)*wmole(ichx2c(icha)) )
      ehchx2 = ehchx2 + den2 *                                    &
         ( ehgaze(ichx2c(icha),i)*                                &
           tbmcr(icel,if2mc(icha))*a2(icha)*wmole(ichx2c(icha)) )
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

! --- Calcul de Somme (X2.H2) pour TH(I)

    x2h2 = zero
    do icla = 1, nclacp
      xck  = rtp(icel,isca(ixck(icla)))
      xch  = rtp(icel,isca(ixch(icla)))
      xash = rtp(icel,isca(inp(icla)))*xmash(icla)
      x2h2 = x2h2                                                 &
           + xch  * ehsoli(ich(ichcor(icla) ),i  )                &
           + xck  * ehsoli(ick(ichcor(icla) ),i  )                &
           + xash * ehsoli(iash(ichcor(icla)),i  )
    enddo

    eh0(icel) = (  fuel1(icel)*ehchx1                             &
                 + fuel2(icel)*ehchx2                             &
                 + fuel3(icel)*ehgaze(ico ,i  )                   &
                 + oxyd(icel) *ehgaze(io2 ,i  )                   &
                 + prod1(icel)*ehgaze(ico2,i  )                   &
                 + prod2(icel)*ehgaze(ih2o,i  )                   &
                 + xiner(icel)*ehgaze(in2 ,i  ) ) * x1(icel)      &
               + x2h2

! --- Calcul de l'enthalpie de l'espece gazeuse CHx1m
!                                            et CHx2m pour TH(I+1)
    ehchx1 = zero
    ehchx2 = zero
    ychx10 = zero
    ychx20 = zero
    do icha = 1, ncharb
      den1   = 1.d0                                               &
             / ( a1(icha)*wmole(ichx1c(icha))                     &
                +b1(icha)*wmole(ico)                              &
                +c1(icha)*wmole(ih2o) )
      ychx10 = ychx10 + den1 *                                    &
        ( tbmcr(icel,if1mc(icha))*a1(icha)*wmole(ichx1c(icha)) )
      ehchx1 = ehchx1 + den1 *                                    &
        ( ehgaze(ichx1c(icha),i+1)*                               &
          tbmcr(icel,if1mc(icha))*a1(icha)*wmole(ichx1c(icha)) )
      den2   = 1.d0                                               &
             / ( a2(icha)*wmole(ichx2c(icha))                     &
                +b2(icha)*wmole(ico)                              &
                +c2(icha)*wmole(ih2o) )
      ychx20 = ychx20 + den2 *                                    &
         ( tbmcr(icel,if2mc(icha))*a2(icha)*wmole(ichx2c(icha)) )
      ehchx2 = ehchx2 + den2 *                                    &
         ( ehgaze(ichx2c(icha),i+1)*                              &
           tbmcr(icel,if2mc(icha))*a2(icha)*wmole(ichx2c(icha)) )
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

! --- Calcul de Somme (X2.H2) pour TH(I+1)

    x2h2 = zero
    do icla = 1, nclacp
      xck  = rtp(icel,isca(ixck(icla)))
      xch  = rtp(icel,isca(ixch(icla)))
      xash = rtp(icel,isca(inp(icla)))*xmash(icla)
      x2h2 = x2h2                                                 &
           + xch  * ehsoli(ich(ichcor(icla) ),i+1)                &
           + xck  * ehsoli(ick(ichcor(icla) ),i+1)                &
           + xash * ehsoli(iash(ichcor(icla)),i+1)
    enddo

    eh1(icel) = (  fuel1(icel)*ehchx1                             &
                 + fuel2(icel)*ehchx2                             &
                 + fuel3(icel)*ehgaze(ico ,i+1)                   &
                 + oxyd(icel) *ehgaze(io2 ,i+1)                   &
                 + prod1(icel)*ehgaze(ico2,i+1)                   &
                 + prod2(icel)*ehgaze(ih2o,i+1)                   &
                 + xiner(icel)*ehgaze(in2 ,i+1) ) * x1(icel)      &
               + x2h2

    if ( eh(icel).ge.eh0(icel) .and. eh(icel).le.eh1(icel) ) then
      tp(icel)= th(i) + (eh(icel)-eh0(icel)) *                    &
                        (th(i+1)-th(i))/(eh1(icel)-eh0(icel))
    endif

  enddo
enddo


!----
! FIN
!----

return
end subroutine
