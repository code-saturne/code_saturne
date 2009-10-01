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

subroutine vorvit &
!================

 ( ncevor , nvor   , ient   , ivorce , visco  ,                   &
   yzcel  , xu     , xv     , xw     ,                            &
   yzvor  , signv  , sigma  , gamma  , temps )

!===============================================================================
!  FONCTION  :
!  ---------

! METHODE DES VORTEX POUR LES CONDITIONS AUX LIMITES D'ENTREE EN
!   L.E.S. : CALCUL DES LA VITESSE DES VORTEX

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
!    nom           !type!mode !                   role                         !
!__________________!____!_____!________________________________________________!
! ncevor           ! e  ! <-- ! nombre de face a l'entree ou est               !
!                  !    !     ! utilise la methode                             !
! nvor             ! e  !     ! nombre de vortex a l'entree ou est             !
!                  !    !     ! utilisee la methode                            !
! ient             ! e  ! <-- ! numero de l'entree ou est utilisee             !
!                  !    !     ! la methode                                     !
! ivorce           ! te ! <-- ! numero du vortex le plus proche d'une          !
!     (nvomax)     !    !     ! face donnee                                    !
! visco            ! tr ! <-- ! viscosite cinematique sur les faces            !
!(icvmax,nnent)    !    !     ! d'entree                                       !
! yzcel            ! tr ! <-- ! coordonnees des faces d'entree dans            !
!   (icvmax ,2)    !    !     ! le referentiel local                           !
! xu(icvmax)       ! tr ! <-- ! composante de vitesse principale               !
! xv(icvmax)       ! tr ! <-- ! composantes de vitesse transverses             !
! xw(icvmax)       ! tr ! <-- !                                                !
! yzvor            ! tr ! <-- ! coordonnees du centre des vortex               !
!   (nvomax,2)     !    !     !                                                !
! signv(nvomax)    ! tr ! <-- ! sens de rotation des vortex                    !
! sigma            ! tr ! <-- ! taille des vortex                              !
!(nvomax,nnent)    !    !     !                                                !
! gamma            ! tr ! <-- ! intensite (circulation) des vortex             !
!(nvomax,2,nnen    !    !     ! dans les deux directions du plan               !
! temps            ! tr ! <-- ! temps ecoule depuis la creation                !
!     (nvomax)     !    !     ! du vortex                                      !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

implicit none

!===============================================================================
!     DONNEES EN COMMON
!===============================================================================

include "paramx.h"
include "cstnum.h"
include "cstphy.h"
include "entsor.h"
include "vortex.h"

!===============================================================================

! Arguments

integer          ncevor , nvor   , ient
integer          ivorce(nvomax)

double precision yzcel(icvmax,2) , visco(icvmax)
double precision xu(icvmax )     , xv(icvmax )     , xw(icvmax )
double precision yzvor(nvomax,2) , signv(nvomax)   , sigma(nvomax)
double precision gamma(nvomax,2) , temps(nvomax)

! VARIABLES LOCALES


integer          ii, jj, iii

double precision yy, zz, xvisc
double precision norme, vv, ww, theta, dv, dw
double precision yparoi, zparoi, yperio, zperio, ysym, zsym
double precision phidat
double precision alpha, ek_vor, ee_vor, v_vor, w_vor
double precision lt, lk, deuxpi

!===============================================================================
! 1. CALCUL DE GAMMA
!===============================================================================

alpha = 4.d0 * sqrt(pi * xsurfv(ient)/                            &
       (3.d0 * nvor*(2.d0*log(3.d0)-3.d0*log(2.d0))))

do ii = 1, nvor
  yy = yzvor(ii,1)
  zz = yzvor(ii,2)
  iii = 0
  ek_vor =  phidat(nfecra,icas(ient),ndat(ient),yy,zz,            &
       ydat(1,ient),zdat(1,ient),kdat(1,ient),iii)
  gamma(ii,1) = alpha*sqrt(ek_vor)*signv(ii)
  gamma(ii,2) = gamma(ii,1)
enddo

!===============================================================================
! 2. CALCUL DE SIGMA
!===============================================================================

if(isgmvo(ient).eq.1) then
  do ii = 1, nvor
    sigma(ii) = xsgmvo(ient)
  enddo
elseif (isgmvo(ient).eq.2) then
  do ii = 1, nvor
    yy = yzvor(ii,1)
    zz = yzvor(ii,2)
    iii = 0
    ek_vor =  phidat(nfecra,icas(ient),ndat(ient),yy,zz,          &
         ydat(1,ient),zdat(1,ient),kdat(1,ient),iii)
    ee_vor =  phidat(nfecra,icas(ient),ndat(ient),yy,zz,          &
         ydat(1,ient),zdat(1,ient),epsdat(1,ient),iii)
    sigma(ii) = (cmu**(0.75d0))*(ek_vor**(1.5d0))/ee_vor
  enddo
elseif (isgmvo(ient).eq.3) then
  do ii = 1, nvor
    yy = yzvor(ii,1)
    zz = yzvor(ii,2)
    iii = 0
    ek_vor =  phidat(nfecra,icas(ient),ndat(ient),yy,zz,          &
         ydat(1,ient),zdat(1,ient),kdat(1,ient),iii)
    ee_vor =  phidat(nfecra,icas(ient),ndat(ient),yy,zz,          &
         ydat(1,ient),zdat(1,ient),epsdat(1,ient),iii)
    xvisc  =  visco(ivorce(ii))
    lt = sqrt(5.d0*xvisc*ek_vor/ee_vor)
    lk = 200.d0*(xvisc**3/ee_vor)**(0.25d0)
    sigma(ii) = max(lt,lk)
  enddo
endif

!===============================================================================
! 3. CALCUL DU CHAMP DE VITESSE INDUIT (AU CENTRE DES CELLULES)
!===============================================================================

deuxpi = 2.d0*pi

do ii = 1, ncevor
  vv = 0.d0
  ww = 0.d0
  do jj = 1, nvor
    yy = yzvor(jj,1)-yzcel(ii,1)
    zz = yzvor(jj,2)-yzcel(ii,2)
    norme = yy**2+zz**2
    alpha = sigma(jj)
    if(norme.gt.epzero.and.norme.lt.4.d0*alpha) then
      theta = -norme/(2.d0*alpha**2)
      theta = exp(theta)
      dv = zz/norme*(1.d0-theta)*gamma(jj,1)*theta/deuxpi
      dw = yy/norme*(1.d0-theta)*gamma(jj,2)*theta/deuxpi
      vv = vv - dv
      ww = ww + dw
    endif

!===============================================================================
! 4. TRAITEMENT DES CONDITIONS AUX LIMITES
!===============================================================================
! suite des boucles en DO
! -----------------------
! Conduite rectangulaire
! -----------------------
    if(icas(ient).eq.1) then

! Periodicites

      if(iclvor(1,ient).eq.3) then
        if(yzvor(jj,1).gt.0.d0) then
          yperio = yzvor(jj,1) - lly(ient)
          zperio = yzvor(jj,2)
        else
          yperio = yzvor(jj,1) + lly(ient)
          zperio = yzvor(jj,2)
        endif
        yy = yperio - yzcel(ii,1)
        zz = zperio - yzcel(ii,2)
        norme = yy**2+zz**2
        alpha = sigma(jj)
        if(norme.gt.epzero.and.norme.lt.4.d0*alpha) then
          theta = -norme/(2.d0*alpha**2)
          theta = exp(theta)
          dv = zz/norme*(1.d0-theta)*gamma(jj,1)*theta/deuxpi
          dw = yy/norme*(1.d0-theta)*gamma(jj,2)*theta/deuxpi
          vv = vv - dv
          ww = ww + dw
        endif
      endif

      if(iclvor(2,ient).eq.3) then
        if(yzvor(jj,2).gt.0.d0) then
          yperio = yzvor(jj,1)
          zperio = yzvor(jj,2) - llz(ient)
        else
          yperio = yzvor(jj,1)
          zperio = yzvor(jj,2) + llz(ient)
        endif
        yy = yperio - yzcel(ii,1)
        zz = zperio - yzcel(ii,2)
        norme = yy**2+zz**2
        alpha = sigma(jj)
        if(norme.gt.epzero.and.norme.lt.4.d0*alpha) then
          theta = -norme/(2.d0*alpha**2)
          theta = exp(theta)
          dv = zz/norme*(1.d0-theta)*gamma(jj,1)*theta/deuxpi
          dw = yy/norme*(1.d0-theta)*gamma(jj,2)*theta/deuxpi
          vv = vv - dv
          ww = ww + dw
        endif
      endif

! Parois

      if(iclvor(1,ient).eq.1) then
        yparoi = lly(ient)/2.d0
        zparoi = yzcel(ii,2)
        yparoi = 2.d0*yparoi - yzvor(jj,1)
        zparoi = 2.d0*zparoi - yzvor(jj,2)
        yy = yparoi - yzcel(ii,1)
        zz = zparoi - yzcel(ii,2)
        norme = yy**2+zz**2
        alpha = sigma(jj)
        if(norme.gt.epzero.and.norme.lt.4.d0*alpha) then
          theta = -norme/(2.d0*alpha**2)
          theta = exp(theta)
          dv = zz/norme*(1.d0-theta)*gamma(jj,1)*theta/deuxpi
          dw = yy/norme*(1.d0-theta)*gamma(jj,2)*theta/deuxpi
          vv = vv - dv
          ww = ww + dw
        endif
      endif

      if(iclvor(3,ient).eq.1) then
        yparoi = - lly(ient)/2.d0
        zparoi = yzcel(ii,2)
        yparoi = 2.d0*yparoi - yzvor(jj,1)
        zparoi = 2.d0*zparoi - yzvor(jj,2)
        yy = yparoi - yzcel(ii,1)
        zz = zparoi - yzcel(ii,2)
        norme = yy**2+zz**2
        alpha = sigma(jj)
        if(norme.gt.epzero.and.norme.lt.4.d0*alpha) then
          theta = -norme/(2.d0*alpha**2)
          theta = exp(theta)
          dv = zz/norme*(1.d0-theta)*gamma(jj,1)*theta/deuxpi
          dw = yy/norme*(1.d0-theta)*gamma(jj,2)*theta/deuxpi
          vv = vv - dv
          ww = ww + dw
        endif
      endif

      if(iclvor(2,ient).eq.1) then
        yparoi = yzcel(ii,1)
        zparoi = llz(ient)/2.d0
        yparoi = 2.d0*yparoi - yzvor(jj,1)
        zparoi = 2.d0*zparoi - yzvor(jj,2)
        yy = yparoi - yzcel(ii,1)
        zz = zparoi - yzcel(ii,2)
        norme = yy**2+zz**2
        alpha = sigma(jj)
        if(norme.gt.epzero.and.norme.lt.4.d0*alpha) then
          theta = -norme/(2.d0*alpha**2)
          theta = exp(theta)
          dv = zz/norme*(1.d0-theta)*gamma(jj,1)*theta/deuxpi
          dw = yy/norme*(1.d0-theta)*gamma(jj,2)*theta/deuxpi
          vv = vv - dv
          ww = ww + dw
        endif
      endif

      if(iclvor(4,ient).eq.1) then
        yparoi = yzcel(ii,1)
        zparoi = -llz(ient)/2.d0
        yparoi = 2.d0*yparoi - yzvor(jj,1)
        zparoi = 2.d0*zparoi - yzvor(jj,2)
        yy = yparoi - yzcel(ii,1)
        zz = zparoi - yzcel(ii,2)
        norme = yy**2+zz**2
        alpha = sigma(jj)
        if(norme.gt.epzero.and.norme.lt.4.d0*alpha) then
          theta = -norme/(2.d0*alpha**2)
          theta = exp(theta)
          dv = zz/norme*(1.d0-theta)*gamma(jj,1)*theta/deuxpi
          dw = yy/norme*(1.d0-theta)*gamma(jj,2)*theta/deuxpi
          vv = vv - dv
          ww = ww + dw
        endif
      endif

! Symetries

      if(iclvor(1,ient).eq.2) then
        ysym = lly(ient)/2.d0
        zsym = yzcel(ii,2)
        ysym = 2.d0*ysym - yzvor(jj,1)
        zsym = 2.d0*zsym - yzvor(jj,2)
        yy = ysym - yzcel(ii,1)
        zz = zsym - yzcel(ii,2)
        norme = yy**2+zz**2
        alpha = sigma(jj)
        if(norme.gt.epzero.and.norme.lt.4.d0*alpha) then
          theta = -norme/(2.d0*alpha**2)
          theta = exp(theta)
          dv = zz/norme*(1.d0-theta)*gamma(jj,1)*theta/deuxpi
          vv = vv - dv
        endif
      endif

      if(iclvor(3,ient).eq.2) then
        ysym = - lly(ient)/2.d0
        zsym = yzcel(ii,2)
        ysym = 2.d0*ysym - yzvor(jj,1)
        zsym = 2.d0*zsym - yzvor(jj,2)
        yy = ysym - yzcel(ii,1)
        zz = zsym - yzcel(ii,2)
        norme = yy**2+zz**2
        alpha = sigma(jj)
        if(norme.gt.epzero.and.norme.lt.4.d0*alpha) then
          theta = -norme/(2.d0*alpha**2)
          theta = exp(theta)
          dv = zz/norme*(1.d0-theta)*gamma(jj,1)*theta/deuxpi
          vv = vv - dv
        endif
      endif

      if(iclvor(2,ient).eq.2) then
        ysym = yzcel(ii,1)
        zsym = llz(ient)/2.d0
        ysym = 2.d0*ysym - yzvor(jj,1)
        zsym = 2.d0*zsym - yzvor(jj,2)
        yy = ysym - yzcel(ii,1)
        zz = zsym - yzcel(ii,2)
        norme = yy**2+zz**2
        alpha = sigma(jj)
        if(norme.gt.epzero.and.norme.lt.4.d0*alpha) then
          theta = -norme/(2.d0*alpha**2)
          theta = exp(theta)
          dw = yy/norme*(1.d0-theta)*gamma(jj,2)*theta/deuxpi
          ww = ww + dw
        endif
      endif

      if(iclvor(4,ient).eq.2) then
        ysym = yzcel(ii,1)
        zsym = -llz(ient)/2.d0
        ysym = 2.d0*ysym - yzvor(jj,1)
        zsym = 2.d0*zsym - yzvor(jj,2)
        yy = ysym - yzcel(ii,1)
        zz = zsym - yzcel(ii,2)
        norme = yy**2+zz**2
        alpha = sigma(jj)
        if(norme.gt.epzero.and.norme.lt.4.d0*alpha) then
          theta = -norme/(2.d0*alpha**2)
          theta = exp(theta)
          dw = yy/norme*(1.d0-theta)*gamma(jj,2)*theta/deuxpi
          ww = ww + dw
        endif
      endif

! --------------------
! Conduite circulaire
! --------------------
    elseif(icas(ient).eq.2) then

      yparoi = yzcel(ii,1)*((lld(ient)/2.0d0)                     &
           /sqrt(yzcel(ii,1)**2 + yzcel(ii,2)**2))
      zparoi = yzcel(ii,2)*((lld(ient)/2.0d0)                     &
           /sqrt(yzcel(ii,1)**2 + yzcel(ii,2)**2))
      yparoi = 2.d0*yparoi - yzvor(jj,1)
      zparoi = 2.d0*zparoi - yzvor(jj,2)

      yy = yparoi - yzcel(ii,1)
      zz = zparoi - yzcel(ii,2)
      norme = yy**2+zz**2
      alpha = sigma(jj)

      if(norme.gt.epzero.and.norme.lt.4.d0*alpha) then
        theta = -norme/(2.d0*alpha**2)
        theta = exp(theta)
        dv = zz/norme*(1.d0-theta)*gamma(jj,1)*theta/deuxpi
        dw = yy/norme*(1.d0-theta)*gamma(jj,2)*theta/deuxpi
        vv = vv - dv
        ww = ww + dw
      endif

    endif

  enddo
  xv(ii) = vv
  xw(ii) = ww
enddo

!===============================================================================
! 5. AJOUT DE LA VITESSE MOYENNE POUR V ET W
!===============================================================================
if(icas(ient).eq.1.or.icas(ient).eq.2.or.icas(ient).eq.3) then
  do ii = 1, ncevor
    yy = yzcel(ii,1)
    zz = yzcel(ii,2)
    iii = 0
    v_vor = phidat(nfecra,icas(ient),ndat(ient),yy,zz,            &
         ydat(1,ient),zdat(1,ient),vdat(1,ient),iii)
    w_vor = phidat(nfecra,icas(ient),ndat(ient),yy,zz,            &
         ydat(1,ient),zdat(1,ient),wdat(1,ient),iii)
    xv(ii) = xv(ii) + v_vor
    xw(ii) = xw(ii) + w_vor
  enddo
endif

!----
! FIN
!----

return
end subroutine
