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

subroutine vordep &
 ( ncevor , nvor   , ient   , dtref  ,                            &
   ivocel , yzc    , xv     , xw     ,                            &
   yzv    , yzva   , xsignv , xtps   , xtpsli )

!===============================================================================
!  FONCTION  :
!  ----------

! METHODE DES VORTEX POUR LES ENTREES EN L.E.S. :
!   - DEPLACEMENT DES VORTEX
!   - TRAITEMENT DES VORTEX MORTS ET SORTANTS

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! ncevor           ! e  ! <-- ! nombre de face a l'entree ou est               !
!                  !    !     ! utilise la methode                             !
! nvor             ! e  ! <-- ! nombre de vortex a l'entree                    !
! ient             ! e  ! <-- ! numero de l'entree                             !
! dtref            ! r  ! <-- ! pas de temps                                   !
! ivocel           ! te ! <-- ! numero du vortex le plus proche d'une          !
!     (nvomax)     !    !     ! face donnee                                    !
! yzc              ! tr ! <-- ! coordonnees des faces d'entree dans            !
!   (icvmax ,2)    !    !     ! le referentiel local                           !
! xv(icvmax)       ! tr ! <-- ! composantes de vitesse transverses             !
! xw(icvmax)       ! tr ! <-- !                                                !
! yzv              ! tr ! --> ! nouvelles coordonnees du centre                !
!   (nvomax,2)     !    !     ! des vortex                                     !
! yzva             ! tr ! <-- ! anciennes coordonnees du centre                !
!   (nvomax,2)     !    !     ! des vortex                                     !
! xsignv(nvomax)   ! tr ! <-- ! sens de rotation des vortex                    !
! xtps             ! tr ! <-- ! temps ecoule depuis la creation                !
!     (nvomax)     !    !     ! du vortex                                      !
! xtpsli           ! tr ! <-- ! duree de vie du vortex                         !
!     (nvomax)     !    !     !                                                !
!__________________.____._____.________________________________________________.

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
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer          ncevor , nvor   , ient
integer          ivocel(nvomax)

double precision dtref
double precision yzc(icvmax ,2)
double precision xv(icvmax)      , xw(icvmax)
double precision yzv(nvomax,2) , yzva(nvomax,2)
double precision xsignv(nvomax)
double precision xtps(nvomax)   , xtpsli(nvomax)

! Local variables

integer          ii, jj, kk, iii

double precision sens, dy, dz
double precision xxv, xxw
double precision drand(1), u_vor, ek_vor, ee_vor
double precision dd, yy, zz
double precision phidat

!===============================================================================
! 1. DEPLACEMENT DES VORTEX
!===============================================================================

do ii = 1, nvor
  xtps(ii) = xtps(ii) + dtref
enddo

! - Deplacement aleatoire

if(idepvo(ient).eq.1)then

  do ii = 1, nvor
    yzva(ii,1) = yzv(ii,1)
    yzva(ii,2) = yzv(ii,2)
  enddo
  do ii = 1, nvor
    sens = 1.d0
    call cs_random_uniform(1, drand)
    if (drand(1).lt.0.5d0) sens = -1.d0
    call cs_random_uniform(1, drand)
    dy = drand(1) * ud(ient) * dtref
    yzv(ii,1) = yzv(ii,1) + sens * dy
    sens = 1.d0
    call cs_random_uniform(1, drand)
    if (drand(1).lt.0.5d0) sens = -1.d0
    call cs_random_uniform(1, drand)
    dz = drand(1) * ud(ient) * dtref
    yzv(ii,2) = yzv(ii,2) + sens * dz
  enddo

! - Convection des vortex

elseif(idepvo(ient).eq.2) then

  do ii = 1, nvor
    yzva(ii,1) = yzv(ii,1)
    yzva(ii,2) = yzv(ii,2)
  enddo

  do ii = 1, nvor
    kk = ivocel(ii)
    xxv = xv(kk)
    xxw = xw(kk)
    yzv(ii,1) = yzv(ii,1) + dtref * xxv
    yzv(ii,2) = yzv(ii,2) + dtref * xxw
  enddo

endif

!===============================================================================
! 2. GESTION DES VORTEX SORTANT DU DOMAINE
!===============================================================================

if(icas(ient).eq.1) then

  if(iclvor(1,ient).eq.1.or.iclvor(1,ient).eq.2) then
    do ii = 1, nvor
      if(yzv(ii,1).gt.(lly(ient)/2.d0)) then
        yzv(ii,1) = yzva(ii,1)
      endif
    enddo
  elseif(iclvor(1,ient).eq.3) then
    do ii = 1, nvor
      if(yzv(ii,1).gt.(lly(ient)/2.d0).and.                     &
        yzv(ii,1).lt.(3.d0*lly(ient)/2.d0)) then
        yzv(ii,1) = yzv(ii,1) - lly(ient)
      elseif(yzv(ii,1).lt.-(lly(ient)/2.d0).and.                &
        yzv(ii,1).gt.-(3.d0*lly(ient)/2.d0)) then
        yzv(ii,1) = yzv(ii,1) + lly(ient)
      elseif(yzv(ii,1).gt.(3.d0*lly(ient)/2.d0).or.             &
        yzv(ii,1).lt.-(3.d0*lly(ient)/2.d0)) then
        yzv(ii,1) = yzva(ii,1)
      endif
    enddo
  endif

  if(iclvor(2,ient).eq.1.or.iclvor(2,ient).eq.2) then
    do ii = 1, nvor
      if(yzv(ii,2).gt.(llz(ient)/2.d0)) then
        yzv(ii,2) = yzva(ii,2)
      endif
    enddo
  elseif(iclvor(2,ient).eq.3) then
    do ii = 1, nvor
      if(yzv(ii,2).gt.(llz(ient)/2.d0).and.                     &
         yzv(ii,2).lt.(3.d0*llz(ient)/2.d0)) then
        yzv(ii,2) = yzv(ii,2) - llz(ient)
      elseif(yzv(ii,2).lt.-(llz(ient)/2.d0).and.                &
        yzv(ii,2).gt.-(3.d0*llz(ient)/2.d0)) then
        yzv(ii,2) = yzv(ii,2) + llz(ient)
      elseif(yzv(ii,2).gt.(3.d0*llz(ient)/2.d0).or.             &
        yzv(ii,2).lt.-(3.d0*llz(ient)/2.d0)) then
        yzv(ii,2) = yzva(ii,2)
      endif
    enddo
  endif

  if(iclvor(3,ient).eq.1.or.iclvor(3,ient).eq.2) then
    do ii = 1, nvor
      if(yzv(ii,1).lt.-(lly(ient)/2.d0)) then
        yzv(ii,1) = yzva(ii,1)
      endif
    enddo
  endif

  if(iclvor(4,ient).eq.1.or.iclvor(4,ient).eq.2) then
    do ii = 1, nvor
      if(yzv(ii,2).lt.-(llz(ient)/2.d0)) then
        yzv(ii,2) = yzva(ii,2)
      endif
    enddo
  endif

elseif(icas(ient).eq.2) then
  do ii = 1, nvor
    if((yzv(ii,1)**2+yzv(ii,2)**2).gt.                        &
      (lld(ient)/2.0d0)**2)then
      yzv(ii,1) = yzva(ii,1)
      yzv(ii,2) = yzva(ii,2)
    endif
  enddo

elseif(icas(ient).eq.3.or.icas(ient).eq.4) then
  do ii = 1, nvor
    if(yzv(ii,1).lt.ymin(ient).or.                              &
      yzv(ii,1).gt.ymax(ient)) then
      yzv(ii,1) = yzva(ii,1)
    endif
    if(yzv(ii,2).lt.zmin(ient).or.                              &
      yzv(ii,2).gt.zmax(ient)) then
      yzv(ii,2) = yzva(ii,2)
    endif
  enddo
endif

!===============================================================================
! 3. REGENERATION DES VORTEX MORTS
!===============================================================================

do ii = 1, nvor
  if(xtps(ii).gt.xtpsli(ii))then
    xtps(ii) = 0.d0

! - Position

    if(icas(ient).eq.1) then
      call cs_random_uniform(1, drand)
      yzv(ii,1) = lly(ient) * drand(1) - lly(ient) / 2.d0
      call cs_random_uniform(1, drand)
      yzv(ii,2) = llz(ient) * drand(1) - llz(ient) / 2.d0
    elseif(icas(ient).eq.2) then
 15   continue
      call cs_random_uniform(1, drand)
      yzv(ii,1) = lld(ient) * drand(1) - lld(ient) / 2.0d0
      call cs_random_uniform(1, drand)
      yzv(ii,2) = lld(ient) * drand(1) - lld(ient) / 2.0d0
      if((yzv(ii,1)**2+yzv(ii,2)**2).gt.                      &
         (lld(ient)/2.0d0)**2) then
        goto 15
      endif
    elseif(icas(ient).eq.3.or.icas(ient).eq.4) then
      call cs_random_uniform(1, drand)
      yzv(ii,1) = ymin(ient) + lly(ient) * drand(1)
      call cs_random_uniform(1, drand)
      yzv(ii,2) = zmin(ient) + llz(ient) * drand(1)
    endif

! - Duree de vie

    if(itlivo(ient).eq.1) then
      xtpsli(ii) = tlimvo(ient)
    elseif(itlivo(ient).eq.2) then
      yy = yzv(ii,1)
      zz = yzv(ii,2)
      iii = 0
      u_vor  = phidat(nfecra,icas(ient),ndat(ient),yy,zz,         &
               ydat(1,ient),zdat(1,ient),udat(1,ient),iii)
      ek_vor = phidat(nfecra,icas(ient),ndat(ient),yy,zz,         &
               ydat(1,ient),zdat(1,ient),kdat(1,ient),iii)
      ee_vor = phidat(nfecra,icas(ient),ndat(ient),yy,zz,         &
               ydat(1,ient),zdat(1,ient),epsdat(1,ient),iii)
      xtpsli(ii) = 5.d0*cmu*ek_vor**(3.d0/2.d0)/ee_vor
      xtpsli(ii) = xtpsli(ii)/u_vor
    endif

! - Sens de rotation

    call cs_random_uniform(1, drand)
    if(drand(1).gt.0.5d0) xsignv(ii) = -1.0d0*xsignv(ii)
  endif
enddo

!===============================================================================
! 4. RECHERCHE DE LA FACE LA PLUS PROCHE DE CHAQUE VORTEX
!===============================================================================

do ii = 1, nvor
   kk = 0
   dd = grand
   do jj = 1, ncevor
     if(((yzc(jj,1)-yzv(ii,1))**2+                            &
         (yzc(jj,2)-yzv(ii,2))**2).lt.dd)then
       dd = (yzc(jj,1)-yzv(ii,1))**2                          &
           +(yzc(jj,2)-yzv(ii,2))**2
       kk = jj
     endif
   enddo
   ivocel(ii) = kk
enddo

! ---
! FIN
! ---

return

end subroutine
