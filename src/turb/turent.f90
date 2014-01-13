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

subroutine keendb &
!================

 ( uref2, dh, xrho, xmu , cmu, xkappa, ustar2, xk, xeps )

!===============================================================================
! FONCTION :
! --------

!    CALCUL DE U*, K ET EPSILON A PARTIR D'UN DIAMETRE ET D'UNE
!      VITESSE DEBITANTE POUR DES ECOULEMENTS EN CONDUITE
!      CIRCULAIRE A PAROI LISSE
!    -> UTILISE POUR LES CONDITIONS AUX LIMITES D'ENTREE

!    ON SORT A LA FOIS U* ET (XK,XEPS) POUR PERMETTRE
!      A L'UTILISATEUR DE CALCULER DES XK ET XEPS DIFFERENTS
!      AVEC LE U*, S'IL LE SOUHAITE


!    ON UTILISE DES LOIS TIREES DE IDEL'CIK
!    LE COEFFICIENT DE PERTE DE CHARGE XLMBDA EST DEFINI PAR
!     |dP/dx| = XLMBDA/DH * 1/2*XRHO*UREF**2

!     PUIS U*=UREF*SQRT(XLMBDA/8)

!    POUR Re < 2000
!      XLMBDA = 64/Re

!    POUR Re > 4000
!      XLMBDA = 1/( 1.8*LOG10(Re)-1.64 )**2

!    POUR 2000 < Re < 4000, ON COMPLETE PAR UNE DROITE
!      XLMBDA = 0.021377 + 5.3115D-6*Re


!    A PARTIR DE U*, ON ESTIME XK ET XEPS A PARTIR DE FORMULES
!      CLASSIQUES DE TURBULENCE DEVELOPPEE

!    XK   = USTAR2/SQRT(CMU)
!    XEPS = USTAR2**1.5D0/(XKAPPA*DH*0.1D0)
!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! uref2            ! r  ! <-- ! carre de la vitesse debitante de               !
!                  !    !     !                            reference           !
! dh               ! r  ! <-- ! diametre hydraulique                           !
! xrho             ! r  ! <-- ! masse volumique                                !
! xmu              ! r  ! <-- ! viscosite dynamique                            !
! cmu              ! r  ! <-- ! constante cmu                                  !
! xkappa           ! r  ! <-- ! constante kappa                                !
! ustar2           ! r  ! --> ! carre de la vitesse de frottement              !
! xk               ! r  ! --> ! intensite turbulente calculee                  !
! xeps             ! r  ! --> ! dissipation turbulente calculee                !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

!===============================================================================

implicit none

! Arguments

double precision uref2, dh, xrho, xmu , ustar2, xk, xeps
double precision cmu, xkappa

! Local variables

double precision re, xlmbda

!===============================================================================

re = sqrt(uref2)*dh*xrho/xmu

if (re.lt.2000) then
!     dans ce cas on calcule directement u*^2 pour eviter un probleme
!      sur xlmbda=64/Re quand Re->0

  ustar2 = 8.d0*xmu*sqrt(uref2)/xrho/dh

else if (re.lt.4000) then

  xlmbda = 0.021377d0 + 5.3115d-6*re
  ustar2 = uref2*xlmbda/8.d0

else

  xlmbda = 1/( 1.8d0*log(re)/log(10.d0)-1.64d0)**2
  ustar2 = uref2*xlmbda/8.d0

endif

xk   = ustar2/sqrt(cmu)
xeps = ustar2**1.5d0/(xkappa*dh*0.1d0)

!----
! FIN
!----

return
end subroutine
subroutine keenin &
!================

 ( uref2, xintur, dh, cmu, xkappa, xk, xeps )

!===============================================================================
! FONCTION :
! --------

!    CALCUL DE U*, K ET EPSILON A PARTIR D'UN DIAMETRE ET D'UNE
!      INTENSITE TURBULENTE
!    -> UTILISE POUR LES CONDITIONS AUX LIMITES D'ENTREE


!    ON ESTIME XK ET XEPS A PARTIR DE FORMULES
!      CLASSIQUES DE TURBULENCE DEVELOPPEE

!    XK   = USTAR2/SQRT(CMU)
!    XEPS = USTAR2**1.5D0/(XKAPPA*DH*0.1D0)
!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! uref2            ! r  ! <-- ! carre de la vitesse debitante de               !
!                  !    !     !                            reference           !
! xintur           ! r  ! <-- ! intensite turbulente                           !
! dh               ! r  ! <-- ! diametre hydraulique                           !
! cmu              ! r  ! <-- ! constante cmu                                  !
! xkappa           ! r  ! <-- ! constante kappa                                !
! xk               ! r  ! --> ! intensite turbulente calculee                  !
! xeps             ! r  ! --> ! dissipation turbulente calculee                !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

!===============================================================================

implicit none

! Arguments

double precision uref2, xintur, dh, cmu, xkappa, xk, xeps

! Local variables

!===============================================================================


xk   = 1.5d0*uref2*xintur**2
xeps = 10.d0*cmu**(0.75d0)*xk**1.5d0/(xkappa*dh)

!----
! End
!----

return
end subroutine
