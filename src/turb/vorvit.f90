!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2018 EDF S.A.
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

subroutine vorvit &
 ( ncevor , nvor   , ient   , ivocel , visco  ,                   &
   yzc    , xv     , xw     ,                                     &
   yzv    , xsignv , xsigma , xgamma )

!===============================================================================
!  FONCTION  :
!  ---------

! METHODE DES VORTEX POUR LES CONDITIONS AUX LIMITES D'ENTREE EN
!   L.E.S. : CALCUL DES LA VITESSE DES VORTEX

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! ncevor           ! e  ! <-- ! nombre de face a l'entree ou est               !
!                  !    !     ! utilise la methode                             !
! nvor             ! e  !     ! nombre de vortex a l'entree ou est             !
!                  !    !     ! utilisee la methode                            !
! ient             ! e  ! <-- ! numero de l'entree ou est utilisee             !
!                  !    !     ! la methode                                     !
! ivocel           ! te ! <-- ! numero du vortex le plus proche d'une          !
!     (nvomax)     !    !     ! face donnee                                    !
! visco            ! tr ! <-- ! viscosite cinematique sur les faces            !
!(icvmax,nnent)    !    !     ! d'entree                                       !
! yzc              ! tr ! <-- ! coordonnees des faces d'entree dans            !
!   (icvmax ,2)    !    !     ! le referentiel local                           !
! xv(icvmax)       ! tr ! <-- ! composantes de vitesse transverses             !
! xw(icvmax)       ! tr ! <-- !                                                !
! yzv              ! tr ! <-- ! coordonnees du centre des vortex               !
!   (nvomax,2)     !    !     !                                                !
! xsignv(nvomax)   ! tr ! <-- ! sens de rotation des vortex                    !
! xsigma           ! tr ! <-- ! taille des vortex                              !
!(nvomax,nnent)    !    !     !                                                !
! xgamma           ! tr ! <-- ! intensite (circulation) des vortex             !
!(nvomax,2,nnen    !    !     ! dans les deux directions du plan               !
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
use cstnum
use cstphy
use entsor
use vorinc

!===============================================================================

implicit none

! Arguments

integer          ncevor , nvor   , ient
integer          ivocel(nvomax)

double precision yzc(icvmax,2) , visco(icvmax)
double precision xv(icvmax )   , xw(icvmax )
double precision yzv(nvomax,2) , xsignv(nvomax)  , xsigma(nvomax)
double precision xgamma(nvomax,2)

! Local variables


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
  yy = yzv(ii,1)
  zz = yzv(ii,2)
  iii = 0
  ek_vor =  phidat(nfecra,icas(ient),ndat(ient),yy,zz,            &
       ydat(1,ient),zdat(1,ient),kdat(1,ient),iii)
  xgamma(ii,1) = alpha*sqrt(ek_vor)*xsignv(ii)
  xgamma(ii,2) = xgamma(ii,1)
enddo

!===============================================================================
! 2. CALCUL DE SIGMA
!===============================================================================

if(isgmvo(ient).eq.1) then
  do ii = 1, nvor
    xsigma(ii) = xsgmvo(ient)
  enddo
elseif (isgmvo(ient).eq.2) then
  do ii = 1, nvor
    yy = yzv(ii,1)
    zz = yzv(ii,2)
    iii = 0
    ek_vor =  phidat(nfecra,icas(ient),ndat(ient),yy,zz,          &
         ydat(1,ient),zdat(1,ient),kdat(1,ient),iii)
    ee_vor =  phidat(nfecra,icas(ient),ndat(ient),yy,zz,          &
         ydat(1,ient),zdat(1,ient),epsdat(1,ient),iii)
    xsigma(ii) = (cmu**(0.75d0))*(ek_vor**(1.5d0))/ee_vor
  enddo
elseif (isgmvo(ient).eq.3) then
  do ii = 1, nvor
    yy = yzv(ii,1)
    zz = yzv(ii,2)
    iii = 0
    ek_vor =  phidat(nfecra,icas(ient),ndat(ient),yy,zz,          &
         ydat(1,ient),zdat(1,ient),kdat(1,ient),iii)
    ee_vor =  phidat(nfecra,icas(ient),ndat(ient),yy,zz,          &
         ydat(1,ient),zdat(1,ient),epsdat(1,ient),iii)
    xvisc  =  visco(ivocel(ii))
    lt = sqrt(5.d0*xvisc*ek_vor/ee_vor)
    lk = 200.d0*(xvisc**3/ee_vor)**(0.25d0)
    xsigma(ii) = max(lt,lk)
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
    yy = yzv(jj,1)-yzc(ii,1)
    zz = yzv(jj,2)-yzc(ii,2)
    norme = yy**2+zz**2
    alpha = xsigma(jj)
    if(norme.gt.epzero.and.norme.lt.4.d0*alpha) then
      theta = -norme/(2.d0*alpha**2)
      theta = exp(theta)
      dv = zz/norme*(1.d0-theta)*xgamma(jj,1)*theta/deuxpi
      dw = yy/norme*(1.d0-theta)*xgamma(jj,2)*theta/deuxpi
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
        if(yzv(jj,1).gt.0.d0) then
          yperio = yzv(jj,1) - lly(ient)
          zperio = yzv(jj,2)
        else
          yperio = yzv(jj,1) + lly(ient)
          zperio = yzv(jj,2)
        endif
        yy = yperio - yzc(ii,1)
        zz = zperio - yzc(ii,2)
        norme = yy**2+zz**2
        alpha = xsigma(jj)
        if(norme.gt.epzero.and.norme.lt.4.d0*alpha) then
          theta = -norme/(2.d0*alpha**2)
          theta = exp(theta)
          dv = zz/norme*(1.d0-theta)*xgamma(jj,1)*theta/deuxpi
          dw = yy/norme*(1.d0-theta)*xgamma(jj,2)*theta/deuxpi
          vv = vv - dv
          ww = ww + dw
        endif
      endif

      if(iclvor(2,ient).eq.3) then
        if(yzv(jj,2).gt.0.d0) then
          yperio = yzv(jj,1)
          zperio = yzv(jj,2) - llz(ient)
        else
          yperio = yzv(jj,1)
          zperio = yzv(jj,2) + llz(ient)
        endif
        yy = yperio - yzc(ii,1)
        zz = zperio - yzc(ii,2)
        norme = yy**2+zz**2
        alpha = xsigma(jj)
        if(norme.gt.epzero.and.norme.lt.4.d0*alpha) then
          theta = -norme/(2.d0*alpha**2)
          theta = exp(theta)
          dv = zz/norme*(1.d0-theta)*xgamma(jj,1)*theta/deuxpi
          dw = yy/norme*(1.d0-theta)*xgamma(jj,2)*theta/deuxpi
          vv = vv - dv
          ww = ww + dw
        endif
      endif

! Parois

      if(iclvor(1,ient).eq.1) then
        yparoi = lly(ient)/2.d0
        zparoi = yzc(ii,2)
        yparoi = 2.d0*yparoi - yzv(jj,1)
        zparoi = 2.d0*zparoi - yzv(jj,2)
        yy = yparoi - yzc(ii,1)
        zz = zparoi - yzc(ii,2)
        norme = yy**2+zz**2
        alpha = xsigma(jj)
        if(norme.gt.epzero.and.norme.lt.4.d0*alpha) then
          theta = -norme/(2.d0*alpha**2)
          theta = exp(theta)
          dv = zz/norme*(1.d0-theta)*xgamma(jj,1)*theta/deuxpi
          dw = yy/norme*(1.d0-theta)*xgamma(jj,2)*theta/deuxpi
          vv = vv - dv
          ww = ww + dw
        endif
      endif

      if(iclvor(3,ient).eq.1) then
        yparoi = - lly(ient)/2.d0
        zparoi = yzc(ii,2)
        yparoi = 2.d0*yparoi - yzv(jj,1)
        zparoi = 2.d0*zparoi - yzv(jj,2)
        yy = yparoi - yzc(ii,1)
        zz = zparoi - yzc(ii,2)
        norme = yy**2+zz**2
        alpha = xsigma(jj)
        if(norme.gt.epzero.and.norme.lt.4.d0*alpha) then
          theta = -norme/(2.d0*alpha**2)
          theta = exp(theta)
          dv = zz/norme*(1.d0-theta)*xgamma(jj,1)*theta/deuxpi
          dw = yy/norme*(1.d0-theta)*xgamma(jj,2)*theta/deuxpi
          vv = vv - dv
          ww = ww + dw
        endif
      endif

      if(iclvor(2,ient).eq.1) then
        yparoi = yzc(ii,1)
        zparoi = llz(ient)/2.d0
        yparoi = 2.d0*yparoi - yzv(jj,1)
        zparoi = 2.d0*zparoi - yzv(jj,2)
        yy = yparoi - yzc(ii,1)
        zz = zparoi - yzc(ii,2)
        norme = yy**2+zz**2
        alpha = xsigma(jj)
        if(norme.gt.epzero.and.norme.lt.4.d0*alpha) then
          theta = -norme/(2.d0*alpha**2)
          theta = exp(theta)
          dv = zz/norme*(1.d0-theta)*xgamma(jj,1)*theta/deuxpi
          dw = yy/norme*(1.d0-theta)*xgamma(jj,2)*theta/deuxpi
          vv = vv - dv
          ww = ww + dw
        endif
      endif

      if(iclvor(4,ient).eq.1) then
        yparoi = yzc(ii,1)
        zparoi = -llz(ient)/2.d0
        yparoi = 2.d0*yparoi - yzv(jj,1)
        zparoi = 2.d0*zparoi - yzv(jj,2)
        yy = yparoi - yzc(ii,1)
        zz = zparoi - yzc(ii,2)
        norme = yy**2+zz**2
        alpha = xsigma(jj)
        if(norme.gt.epzero.and.norme.lt.4.d0*alpha) then
          theta = -norme/(2.d0*alpha**2)
          theta = exp(theta)
          dv = zz/norme*(1.d0-theta)*xgamma(jj,1)*theta/deuxpi
          dw = yy/norme*(1.d0-theta)*xgamma(jj,2)*theta/deuxpi
          vv = vv - dv
          ww = ww + dw
        endif
      endif

! Symetries

      if(iclvor(1,ient).eq.2) then
        ysym = lly(ient)/2.d0
        zsym = yzc(ii,2)
        ysym = 2.d0*ysym - yzv(jj,1)
        zsym = 2.d0*zsym - yzv(jj,2)
        yy = ysym - yzc(ii,1)
        zz = zsym - yzc(ii,2)
        norme = yy**2+zz**2
        alpha = xsigma(jj)
        if(norme.gt.epzero.and.norme.lt.4.d0*alpha) then
          theta = -norme/(2.d0*alpha**2)
          theta = exp(theta)
          dv = zz/norme*(1.d0-theta)*xgamma(jj,1)*theta/deuxpi
          vv = vv - dv
        endif
      endif

      if(iclvor(3,ient).eq.2) then
        ysym = - lly(ient)/2.d0
        zsym = yzc(ii,2)
        ysym = 2.d0*ysym - yzv(jj,1)
        zsym = 2.d0*zsym - yzv(jj,2)
        yy = ysym - yzc(ii,1)
        zz = zsym - yzc(ii,2)
        norme = yy**2+zz**2
        alpha = xsigma(jj)
        if(norme.gt.epzero.and.norme.lt.4.d0*alpha) then
          theta = -norme/(2.d0*alpha**2)
          theta = exp(theta)
          dv = zz/norme*(1.d0-theta)*xgamma(jj,1)*theta/deuxpi
          vv = vv - dv
        endif
      endif

      if(iclvor(2,ient).eq.2) then
        ysym = yzc(ii,1)
        zsym = llz(ient)/2.d0
        ysym = 2.d0*ysym - yzv(jj,1)
        zsym = 2.d0*zsym - yzv(jj,2)
        yy = ysym - yzc(ii,1)
        zz = zsym - yzc(ii,2)
        norme = yy**2+zz**2
        alpha = xsigma(jj)
        if(norme.gt.epzero.and.norme.lt.4.d0*alpha) then
          theta = -norme/(2.d0*alpha**2)
          theta = exp(theta)
          dw = yy/norme*(1.d0-theta)*xgamma(jj,2)*theta/deuxpi
          ww = ww + dw
        endif
      endif

      if(iclvor(4,ient).eq.2) then
        ysym = yzc(ii,1)
        zsym = -llz(ient)/2.d0
        ysym = 2.d0*ysym - yzv(jj,1)
        zsym = 2.d0*zsym - yzv(jj,2)
        yy = ysym - yzc(ii,1)
        zz = zsym - yzc(ii,2)
        norme = yy**2+zz**2
        alpha = xsigma(jj)
        if(norme.gt.epzero.and.norme.lt.4.d0*alpha) then
          theta = -norme/(2.d0*alpha**2)
          theta = exp(theta)
          dw = yy/norme*(1.d0-theta)*xgamma(jj,2)*theta/deuxpi
          ww = ww + dw
        endif
      endif

! --------------------
! Conduite circulaire
! --------------------
    elseif(icas(ient).eq.2) then

      yparoi = yzc(ii,1)*((lld(ient)/2.0d0)                     &
           /sqrt(yzc(ii,1)**2 + yzc(ii,2)**2))
      zparoi = yzc(ii,2)*((lld(ient)/2.0d0)                     &
           /sqrt(yzc(ii,1)**2 + yzc(ii,2)**2))
      yparoi = 2.d0*yparoi - yzv(jj,1)
      zparoi = 2.d0*zparoi - yzv(jj,2)

      yy = yparoi - yzc(ii,1)
      zz = zparoi - yzc(ii,2)
      norme = yy**2+zz**2
      alpha = xsigma(jj)

      if(norme.gt.epzero.and.norme.lt.4.d0*alpha) then
        theta = -norme/(2.d0*alpha**2)
        theta = exp(theta)
        dv = zz/norme*(1.d0-theta)*xgamma(jj,1)*theta/deuxpi
        dw = yy/norme*(1.d0-theta)*xgamma(jj,2)*theta/deuxpi
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
    yy = yzc(ii,1)
    zz = yzc(ii,2)
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
