!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2020 EDF S.A.
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

subroutine vorlgv &
 ( ncevor , ient   , dtref  ,                                     &
   yzc    , xu     , xv     , xw     )

!===============================================================================
!  FONCTION  :
!  ---------

! GENRATION DES FLUCUTUATIONS DE VITESSE DANS LA DIRECTION
! PRINCIPALE A PARTIR DE L'EQUATION DE LANGEVIN

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! ncevor           ! e  ! <-- ! nombre de face a l'entree ou est               !
!                  !    !     ! utilise la methode                             !
! ient             ! e  ! <-- ! numero de l'entree                             !
! dtref            ! r  ! <-- ! pas de temps                                   !
! yzc              ! tr ! <-- ! coordonnees des faces d'entree dans            !
!   (icvmax ,2)    !    !     ! le referentiel local                           !
! xu(icvmax)       ! tr ! <-- ! composante de vitesse principale               !
! xv(icvmax)       ! tr ! <-- ! composantes de vitesse transverses             !
! xw(icvmax)       ! tr ! <-- !                                                !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use entsor
use vorinc

use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer          ncevor, ient

double precision dtref
double precision yzc(icvmax ,2)
double precision xu(icvmax ), xv(icvmax ), xw(icvmax )

!     VARIABLES LOCALES

integer          ii, iii, iun
double precision dw1(1)

double precision cst1, cst2
parameter       (cst1 = 1.8d0)
parameter       (cst2 = 0.6d0)

double precision sinth, costh
double precision norme, ufluc, ek_vor, ee_vor, u_vor, du_vor
double precision phidat, vfluc, yy, zz

!===============================================================================
! 1. CALCUL DES FLUCTUATIONS DE VITESSE (SELON LA DIRECTION X)
!===============================================================================

iun = 1

do ii = 1, ncevor
  yy = yzc(ii,1)
  zz = yzc(ii,2)

  iii = 0
  u_vor = phidat(nfecra,icas(ient),ndat(ient),yy,zz,              &
          ydat(1,ient),zdat(1,ient),udat(1,ient),iii)

  if(icas(ient).eq.2) then
    du_vor =  phidat(nfecra,icas(ient),ndat(ient),yy,zz,          &
              ydat(1,ient),zdat(1,ient),dudat(1,ient),iii)
    ek_vor =  phidat(nfecra,icas(ient),ndat(ient),yy,zz,          &
              ydat(1,ient),zdat(1,ient),kdat(1,ient),iii)
    ee_vor =  phidat(nfecra,icas(ient),ndat(ient),yy,zz,          &
              ydat(1,ient),zdat(1,ient),epsdat(1,ient),iii)

    ufluc = xu(ii) - u_vor

    norme = sqrt(yzc(ii,1)**2+yzc(ii,2)**2)
    costh = yzc(ii,1) / norme
    sinth = yzc(ii,2) / norme
    vfluc = - costh*xv(ii) - sinth*xw(ii)

! le signe - vient du fait que l'on veut
! la fluctuation dans la direction normale
! a la paroi

    call cs_random_normal(1, dw1)
    ufluc = ( ufluc -                                             &
            (1.d0-2.d0/3.d0*cst2)*du_vor*vfluc*dtref +            &
             2.d0* sqrt(2.d0/3.d0*(cst1-1.d0)*ee_vor*dtref)       &
            *dw1(1))
    ufluc = ufluc/(1.d0+(0.5d0*cst1*dtref*ee_vor                  &
                 /ek_vor))
    xu(ii) = u_vor + ufluc
  else
    xu(ii) = u_vor
  endif
enddo

! ---
! FIN
! ---

return

end subroutine
