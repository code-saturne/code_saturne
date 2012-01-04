!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2012 EDF S.A.
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

subroutine hturbp &
!================

 ( prl    , prt    , ckarm  , yplus  , htur   )

!===============================================================================

! FONCTION :
! --------
! 1) CALCUL DU COEFFICIENT CORRECTEUR DU COEFFICIENT D'ECHANGE
!   ENTRE LE FLUIDE ET LA PAROI POUR UN ECOULEMENT TURBULENT
!   EN FONCTION DE LA DISTANCE ADIMENSIONELLE YPLUS = USTAR*DP/RNU
!   HTUR = PR*YPLUS/TPLUS


! CE COEFFICIENT EST CALCULE A L'AIDE D'UN MODELE DE SIMILITUDE
! ENTRE COUCHE LIMITE DYNAMIQUE ET COUCHE LIMITE THERMIQUE

! LE  TPLUS EST CALCULE :

! - POUR UN NOMBRE DE PRANDTL  << 0.1 (METAUX LIQUIDES) :
!   PAR LE MODELE STANDARD A DEUX COUCHES (PRANDTL-TAYLOR)

! - POUR UN NOMBRE DE PRANDTL  >> 0.1  (LIQUIDES et GAZ):
!   PAR UN MODELE A TROIS COUCHES (ARPACI-LARSEN)

! -->> LE COEFFICIENT D'ECHANGE FINAL : H = (K/dp)*htur

!-------------------------------------------------------------------------------
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! yplus            ! r  ! <-- ! distance a la paroi adimensionnelle            !
! ckarm            ! r  ! <-- ! constante de karman                            !
! prt              ! r  ! <-- ! nombre de prandtl turbulent                    !
! prl              ! r  ! <-- ! nombre de prandtl moleculaire                  !
! htur             ! r  ! --> ! coefficient correcteur d'echange(adim          !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!-------------------------------------------------------------------------------
!===============================================================================

implicit none

!===============================================================================

! Arguments

double precision htur
double precision prl,ckarm,prt,yplus

! Local variables

double precision tplus
double precision beta2,a2
double precision yp0,yp1,yp2
double precision prlm1

!============================================================================

!     1)INITIALISATIONS
!     -----------------

htur = 1.d0

prlm1 = 0.1d0

yp0   = prt/(prl*ckarm)
yp2   = ckarm*1000.d0/prt
yp2   = sqrt(yp2)
yp1   = (1000.d0/prl)**(1.d0/3.d0)

!     ====================================================
!     2) CALCUL DU COEFFICIENT CORRECTEUR
!        POUR LES NOMBRES DE PRANDTL TRES PETITS
!     ====================================================
if( prl .le. prlm1) then
  if(yplus .gt. yp0) then
    tplus = prl*yp0 + prt/ckarm * log(yplus/yp0)
    htur = prl*yplus/tplus
  endif

endif


!     ====================================================
!     3) CALCUL DU COEFFICIENT CORRECTEUR
!        POUR UN MODELE A TROIS COUCHES
!     ====================================================
if( prl .gt. prlm1) then

 a2 = 15.d0*(prl**(2.d0/3.d0))
 beta2 = a2 - 500.d0/ (yp2**2)

 if( (yplus .ge. yp1) .and. (yplus.lt.yp2) )then
    tplus = a2 - 500.d0/(yplus*yplus)
    htur = prl*yplus/tplus
 endif

 if( (yplus .ge. yp2) )then
    tplus = beta2 + prt/ckarm*log(yplus/yp2)
    htur = prl*yplus/tplus
 endif

endif

!----
! FIN
!----

return

end subroutine
