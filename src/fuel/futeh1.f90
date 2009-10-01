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

subroutine futeh1 &
!================

 ( ncelet , ncel   ,                                              &
   eh     ,                                                       &
   fuel1  , fuel3  , oxyd   , prod1  , prod2  ,                   &
   xiner  , xh2s   , xso2   ,                                     &
   tp     ,                                                       &
   eh0    , eh1    )

!===============================================================================
! FONCTION :
! --------
! CALCUL DE LA TEMPERATURE DU GAZ
!  EN FONCTION DE L'ENTHALPIE DU GAZ ET DES CONCENTRATIONS

! Arguments
!__________________.____._____.________________________________________________.
!    nom           !type!mode !                   role                         !
!__________________!____!_____!________________________________________________!
! ncelet           ! e  ! <-- ! nombre d'elements halo compris                 !
! ncel             ! e  ! <-- ! nombre d'elements actifs                       !
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
! eh0              ! tr ! <-- ! tableau reel de travail                        !
! eh1              ! tr ! <-- ! tableau reel de travail                        !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHAMNUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

implicit none

!==============================================================================
!     DONNEES EN COMMON
!==============================================================================

include "paramx.h"
include "numvar.h"
include "optcal.h"
include "cstphy.h"
include "cstnum.h"
include "entsor.h"
include "ppppar.h"
include "ppthch.h"
include "coincl.h"
include "cpincl.h"
include "fuincl.h"
include "ppincl.h"

!===============================================================================

! Arguments

integer          ncelet , ncel

double precision eh(ncelet)
double precision fuel1(ncelet), fuel3(ncelet)
double precision oxyd(ncelet), xiner(ncelet)
double precision prod1(ncelet),prod2(ncelet)
double precision xh2s(ncelet) , xso2(ncelet)
double precision tp(ncelet)
double precision eh0(ncelet), eh1(ncelet)

! VARIABLES LOCALES

integer          ii, icel

!===============================================================================

ii = npo-1
do icel = 1, ncel

! --- Clipping eventuel de TP a TH(NPO) si EH > EH1

  eh1(icel) = fuel1(icel)*ehgaze(ifov,ii+1)                       &
       + fuel3(icel)*ehgaze(ico ,ii+1)                            &
       + oxyd(icel) *ehgaze(io2 ,ii+1)                            &
       + prod1(icel)*ehgaze(ico2,ii+1)                            &
       + prod2(icel)*ehgaze(ih2o,ii+1)                            &
       + xiner(icel)*ehgaze(in2 ,ii+1)                            &
       + xh2s (icel)*ehgaze(ih2s,ii+1)                            &
       + xso2 (icel)*ehgaze(iso2,ii+1)

  if ( eh(icel).ge.eh1(icel) ) tp(icel)= th(ii+1)
enddo

ii = 1
do icel = 1, ncel

! --- Clipping eventuel de TP a TH(1) si EH < EH0

  eh0(icel) = fuel1(icel)*ehgaze(ifov,ii)                         &
            + fuel3(icel)*ehgaze(ico ,ii)                         &
            + oxyd(icel) *ehgaze(io2 ,ii)                         &
            + prod1(icel)*ehgaze(ico2,ii)                         &
            + prod2(icel)*ehgaze(ih2o,ii)                         &
            + xiner(icel)*ehgaze(in2 ,ii)                         &
            + xh2s (icel)*ehgaze(ih2s,ii)                         &
            + xso2 (icel)*ehgaze(iso2,ii)

  if ( eh(icel).le.eh0(icel) ) tp(icel)= th(1)
enddo

do ii = 1, npo-1
  do icel = 1, ncel

    eh0(icel) = fuel1(icel)*ehgaze(ifov,ii)                       &
               +fuel3(icel)*ehgaze(ico ,ii)                       &
               +oxyd(icel) *ehgaze(io2 ,ii)                       &
               +prod1(icel)*ehgaze(ico2,ii)                       &
               +prod2(icel)*ehgaze(ih2o,ii)                       &
               +xiner(icel)*ehgaze(in2 ,ii)                       &
               +xh2s (icel)*ehgaze(ih2s,ii)                       &
               +xso2 (icel)*ehgaze(iso2,ii)

    eh1(icel) = fuel1(icel)*ehgaze(ifov,ii+1)                     &
               +fuel3(icel)*ehgaze(ico ,ii+1)                     &
               +oxyd(icel) *ehgaze(io2 ,ii+1)                     &
               +prod1(icel)*ehgaze(ico2,ii+1)                     &
               +prod2(icel)*ehgaze(ih2o,ii+1)                     &
               +xiner(icel)*ehgaze(in2 ,ii+1)                     &
               +xh2s (icel)*ehgaze(ih2s,ii+1)                     &
               +xso2 (icel)*ehgaze(iso2,ii+1)

    if ( eh(icel).ge.eh0(icel) .and. eh(icel).le.eh1(icel) ) then
      tp(icel)= th(ii) + (eh(icel)-eh0(icel))                     &
                       *(th(ii+1)-th(ii))/(eh1(icel)-eh0(icel))
    endif
  enddo
enddo

!----
! FIN
!----

return
end subroutine
