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

subroutine causta &
!================

 ( iface , imp    , xkappa , cstlog , ypluli ,                    &
   apow  , bpow   , dpow   ,                                      &
   uu    , dp     , xnu    , uet      )

!===============================================================================
!     FONCTION  :
!     --------

!  CALCUL DE LA VITESSE DE FROTTEMENT UET PAR UNE METHODE DE NEWTON
!   EN SUPPOSANT QUE UU VITESSE TANGENTIELLE REPOND A LA LOI
!    UU/UET = 1/XKAPPA LOG (YPLUS) + CSTLOG
!    YPLUS = UET*DP/XNU
!-------------------------------------------------------------------------------
!ARGU                         ARGUMENTS
!__________________.____._____.________________________________________________.
!    nom           !type!mode !                   role                         !
!__________________!____!_____!________________________________________________!
! iface            ! e  ! <-- ! numero d'entite traitee                        !
! imp              ! e  ! <-- ! niveau d'impression                            !
!                  !    !     ! > 0 en cas d'erreur fatale                     !
!                  !    !     ! > 1 non respect d'un critere                   !
!                  !    !     ! > 2 mini pour bonne analyse                    !
!                  !    !     ! > 3 intervention inhabituelle du code          !
!                  !    !     ! > 4 impression pour info                       !
!                  !    !     ! > 5 idem + complements                         !
! xkappa           ! r  ! <-- ! cst de karman                                  !
! cstlog           ! r  ! <-- ! cst de la loi log                              !
! ypluli           ! r  ! <-- ! yplus limite                                   !
! *pow             ! r  ! <-- ! coefficients de la loi de werner               !
! uu               ! r  ! <-- ! vitesse tangentielle a la paroi                !
! dp               ! r  ! <-- ! dist a la paroi du pt ou est pris uu           !
! xnu              ! r  ! <-- ! viscosite cinematique moleculaire              !
! uet              ! r  ! --> ! vitesse de frottement                          !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

implicit none

!==============================================================================
!     DONNEES EN COMMON
!==============================================================================

include "paramx.h"
include "entsor.h"

!===============================================================================

! Arguments

integer          iface, imp
double precision xkappa, cstlog, ypluli, apow, bpow, dpow
double precision uu, dp, xnu, uet

! VARIABLES LOCALES

integer          ipos, nitm, nit
double precision eps, rr, yp, uetwer, uetmin, ueta, ysurnu

!===============================================================================

!     -------------------
!     1 - INITIALISATIONS
!     -------------------

eps  = 0.001d0
nitm = 100
nit  = 0

ipos = 0
if ( imp .gt. 4 ) then
   write ( nfecra,1000 ) iface
   ipos = 1
   write ( nfecra,1100 ) uu, dp
endif

!     --------------------------
!   2 - CALCUL DU REYNOLDS LOCAL
!     --------------------------

ysurnu = dp / xnu
rr = uu * ysurnu

!     ------------------------------------------
!     3 - CALCUL DE LA VITESSE DE FROTTEMENT UET
!     ------------------------------------------

!       SOUS-COUCHE VISQUEUSE

if ( rr .le. ypluli**2 ) then
  uet = sqrt( uu / ysurnu )
else

!       RESOLUTION ITERATIVE PAR UNE METHODE DE NEWTON

!     On initialise par Werner ou la borne inf de uet assurant la cv
!     de la methode
  uetwer = ( abs(uu)/apow/(ysurnu**bpow) )**dpow
  uetmin = exp(-cstlog*xkappa)/ysurnu
  uet = max(uetwer,uetmin)

 300    continue
  nit = nit+1

  ueta = uet
  uet = (xkappa*uu + uet)/(log(ysurnu*uet)+xkappa*cstlog+1.d0)

  if ( abs ( uet-ueta ) .le. eps * ueta ) then
    if ( imp .gt. 5 ) write ( nfecra,3000 ) nit, eps
  else if ( nit .ge. nitm ) then
    if ( imp .gt. 1 ) then
      if ( ipos .le. 0 ) write ( nfecra,1000 ) iface
      ipos = 1
      write ( nfecra,3100 ) nitm, eps
    endif
  else
    go to 300
  endif

endif

!     ------------------------------------------------------
!     4 -  POSITION DU POINT INTERIEUR DANS LA COUCHE LIMITE
!     ------------------------------------------------------

if ( imp .gt. 4 ) then
  yp = ysurnu*uet
  if ( yp .le. ypluli ) then
    write ( nfecra,4000 ) yp
  else
    write ( nfecra,4100 ) yp
  endif
endif

!--------
! FORMATS
!--------

#if defined(_CS_LANG_FR)

 1000 FORMAT ( 5X,'APPEL DU SOUS-PROGRAMME CAUSTA A LA FACE IFAC=',I9)
 1100 FORMAT ( 10X,'CONDITIONS D''ENTREE : UU=',E10.3,2X,'DP=',E10.3 )
 3000 FORMAT ( 10X,'CONVERGENCE DU CALCUL DE UET APRES NIT=',I4,        &
 1X,'ITERATIONS, PRECISION DEMANDEE EPS=',E10.3 )
 3100 FORMAT ( 10X,'NOMBRE MAXIMAL D''ITERATIONS ATTEINT POUR LE',      &
 1X,'CALCUL DE UET : NITM=',I4,', PRECISION',                     &
 1X,'DEMANDEE EPS=',E10.3 )
 4000 FORMAT ( 10X,'LE PREMIER POINT INTERIEUR EST DANS LA SOUS-COUCHE',&
 1X,'VISQUEUSE, Y+=',E10.3 )
 4100 FORMAT ( 10X,'LE PREMIER POINT INTERIEUR EST A LA DISTANCE DE',   &
 1X,'LA PAROI, Y+=',E10.3 )

#else

 1000 FORMAT ( 5X,'CAUSTA SUBROUTINE CALLED FOR FACE IFAC=',I9)
 1100 FORMAT ( 10X,'INPUT CONDITIONS: UU=',E10.3,2X,'DP=',E10.3 )
 3000 FORMAT ( 10X,'UET COMPUTATION CONVERGENCE AFTER NIT=',I4,         &
 1X,'ITERATIONS, DESIRED PRECISION EPS=',E10.3 )
 3100 FORMAT ( 10X,'MAXIMUM NUMBER OF ITERATIONS REACHED FOR THE',      &
 1X,'COMPUTATION OF UET: NITM=',I4,', DESIRED PRECISION',         &
 1X,'EPS=',E10.3 )
 4000 FORMAT ( 10X,'THE FIRST POINT IS IN THE VISCOUS SUBLAYER',        &
 1X,' Y+=',E10.3 )
 4100 FORMAT ( 10X,'THE FIRST POINT IS AT A WALL-DISTANCE OF',          &
 1X,' Y+=',E10.3 )

#endif

!----
! FIN
!----

return

end subroutine
