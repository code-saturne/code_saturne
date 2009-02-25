!-------------------------------------------------------------------------------

!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2008 EDF S.A., France

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

                  subroutine cpteh2                               &
!================

 ( ncelet , ncel   ,                                              &
   rtp    , propce ,                                              &
   eh0    , eh1     )

!===============================================================================
! FONCTION :
! --------
! CALCUL DE LA TEMPERATURE DES PARTICULES
!  EN FONCTION DE L'ENTHALPIE DU SOLIDE ET DES CONCENTRATIONS

! Arguments
!__________________.____._____.________________________________________________.
!    nom           !type!mode !                   role                         !
!__________________!____!_____!________________________________________________!
! ncelet           ! e  ! <-- ! nombre d'elements halo compris                 !
! ncel             ! e  ! <-- ! nombre d'elements actifs                       !
! rtp              ! tr ! <-- ! variables de calcul au centre des              !
! (ncelet,*)       !    !     !    cellules (instant courant)                  !
! propce           ! tr ! <-- ! proprietes physiques au centre des             !
! (ncelet,*)       !    !     !    cellules                                    !
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
include "ppincl.h"

!===============================================================================

! Arguments

integer          ncelet, ncel
double precision rtp(ncelet,*), propce(ncelet,*)
double precision eh0(ncelet), eh1(ncelet)

! VARIABLES LOCALES

integer          i      , icla   , icha   , icel
integer          ipcte1 , ipcte2
integer          ihflt2
double precision h2     , x2     , xch    , xck
double precision xash   , xnp    , xtes   , xwat

!===============================================================================
! RQ IMPORTANTE : On suppose pour l'instant que H2 = H02 + CP2(T2-TREF)

ihflt2 = 1

!===============================================================================
! 1. CALCULS PRELIMINAIRES
!===============================================================================

! --- Initialisation des tableaux

do icel = 1, ncel
  eh0(icel) = zero
  eh1(icel) = zero
enddo

! --- Initialisation de T2 a T1

ipcte1 = ipproc(itemp1)
do icla = 1, nclacp
  ipcte2 = ipproc(itemp2(icla))
  do icel = 1, ncel
    propce(icel,ipcte2) = propce(icel,ipcte1)
  enddo
enddo

!===============================================================================
! 2. CALCUL DE LA TEMPERATURE DES PARTICULES
!===============================================================================

if ( ihflt2.eq.0 ) then

! --> H2 fonction lineaire de T2

  do icla = 1, nclacp
    ipcte2 = ipproc(itemp2(icla))
    icha = ichcor(icla)
    do icel = 1, ncel
      propce(icel,ipcte2) =                                       &
            (rtp(icel,isca(ih2(icla)))-h02ch(icha))               &
            / cp2ch(icha) + trefth
    enddo
  enddo

else

! --> H2 tabule

  do icla = 1, nclacp

    ipcte2 = ipproc(itemp2(icla))

    i = npoc-1
    do icel = 1, ncel
      xch  = rtp(icel,isca(ixch(icla)))
      xnp  = rtp(icel,isca(inp(icla)))
      xck  = rtp(icel,isca(ixck(icla)))
      xash = xmash(icla)*xnp
      if ( ippmod(icp3pl) .eq. 1 ) then
        xwat = rtp(icel,isca(ixwt(icla)))
      else
        xwat = 0.d0
      endif

      x2   = xch + xck + xash + xwat

      if ( x2 .gt. epsicp*100.d0 ) then

        h2   = rtp(icel,isca(ih2(icla)))/x2

        xtes = xmp0(icla)*xnp

        if ( xtes.gt.epsicp .and. x2.gt.epsicp*100.d0 ) then
          eh1(icel) = xch /x2 * ehsoli(ich(ichcor(icla) ),i+1)    &
                    + xck /x2 * ehsoli(ick(ichcor(icla) ),i+1)    &
                    + xash/x2 * ehsoli(iash(ichcor(icla)),i+1)    &
                    + xwat/x2 * ehsoli(iwat(ichcor(icla)),i+1)
          if ( h2.ge.eh1(icel) ) then
            propce(icel,ipcte2) = thc(i+1)

          endif
        endif

      endif

    enddo

    i = 1
    do icel = 1, ncel
      xch  = rtp(icel,isca(ixch(icla)))
      xnp  = rtp(icel,isca(inp(icla)))
      xck  = rtp(icel,isca(ixck(icla)))
      xash = xmash(icla)*xnp
      if ( ippmod(icp3pl) .eq. 1 ) then
        xwat = rtp(icel,isca(ixwt(icla)))
      else
        xwat = 0.d0
      endif

      x2   = xch + xck + xash + xwat

      if ( x2 .gt. epsicp*100.d0 ) then

        h2   = rtp(icel,isca(ih2(icla)))/x2

        xtes = xmp0(icla)*xnp

        if ( xtes.gt.epsicp .and. x2.gt.epsicp*100.d0 ) then
          eh0(icel) = xch /x2 * ehsoli(ich(ichcor(icla) ),i)      &
                    + xck /x2 * ehsoli(ick(ichcor(icla) ),i)      &
                    + xash/x2 * ehsoli(iash(ichcor(icla)),i)      &
                    + xwat/x2 * ehsoli(iwat(ichcor(icla)),i)
          if ( h2.le.eh0(icel) ) then
            propce(icel,ipcte2) = thc(i)
          endif
        endif

      endif

    enddo

    do i = 1, npoc-1
      do icel = 1, ncel
        xch  = rtp(icel,isca(ixch(icla)))
        xnp  = rtp(icel,isca(inp(icla)))
        xck  = rtp(icel,isca(ixck(icla)))
        xash = xmash(icla)*xnp
        if ( ippmod(icp3pl) .eq. 1 ) then
          xwat = rtp(icel,isca(ixwt(icla)))
        else
          xwat = 0.d0
        endif

        x2   = xch + xck + xash + xwat

        if ( x2 .gt. epsicp*100.d0 ) then

          h2   = rtp(icel,isca(ih2(icla)))/x2

          xtes = xmp0(icla)*xnp

          if ( xtes.gt.epsicp .and. x2.gt.epsicp*100.d0 ) then
            eh0(icel) = xch /x2 * ehsoli(ich(ichcor(icla) ),i  )  &
                      + xck /x2 * ehsoli(ick(ichcor(icla) ),i  )  &
                      + xash/x2 * ehsoli(iash(ichcor(icla)),i  )  &
                      + xwat/x2 * ehsoli(iwat(ichcor(icla)),i  )

            eh1(icel) = xch /x2 * ehsoli(ich(ichcor(icla) ),i+1)  &
                      + xck /x2 * ehsoli(ick(ichcor(icla) ),i+1)  &
                      + xash/x2 * ehsoli(iash(ichcor(icla)),i+1)  &
                      + xwat/x2 * ehsoli(iwat(ichcor(icla)),i+1)

            if ( h2.ge.eh0(icel) .and. h2.le.eh1(icel) ) then
              propce(icel,ipcte2) = thc(i) + (h2-eh0(icel)) *     &
                    (thc(i+1)-thc(i))/(eh1(icel)-eh0(icel))
            endif
          endif

        endif

      enddo
    enddo

  enddo

endif

!----
! FIN
!----

return
end
