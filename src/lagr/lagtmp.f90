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

subroutine lagtmp                                                              &
!================

 ( nbpmax , nvp    , nvp1  , nvep  , nivep ,                                   &
   npt    ,                                                                    &
   itepa  , propce , ettp  , ettpa , tepa  , tempct ,                          &
   rayon  , mlayer , phith , temp  , tsvar , volume_couche )

!===============================================================================
! FONCTION :
! ----------

!   SOUS-PROGRAMME DU MODULE LAGRANGIEN :
!   -------------------------------------

!     EVOLUTION DE LA TEMPERATURE D'UNE PARTICULE DE CHARBON OU DE BIOMASSE

!     1. CALCUL DES FLUX THERMIQUES (RAYONNEMENT ET CONDUCTION)

!     2. RESOLUTION DE L'EQUATION DE DIFFUSION DE LA CHALEUR AU SEIN D'UNE
!        PARTICULE DE CHARBON OU DE BIOMASSE AVEC PRISE EN COMPTE DES FLUX
!        VOLUMIQUES (REACTION CHIMIQUE)

!        rho cp dT/dt = div( lambda grad(T) ) + phi

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nbpmax           ! e  ! <-- ! nombre max de particulies autorise             !
! nvp              ! e  ! <-- ! nombre de variables particulaires              !
! nvp1             ! e  ! <-- ! nvp sans position, vfluide, vpart              !
! nvep             ! e  ! <-- ! nombre info particulaires (reels)              !
! nivep            ! e  ! <-- ! nombre info particulaires (entiers)            !
! npt              ! e  ! <-- ! numero de la particule a traiter               !
! itepa            ! te ! <-- ! info particulaires (entiers)                   !
! (nbpmax,nivep    !    !     !   (cellule de la particule,...)                !
! propce(ncelet, *)! tr ! <-- ! physical properties at cell centers            !
! ettp             ! tr ! <-- ! tableaux des variables liees                   !
!  (nbpmax,nvp)    !    !     !   aux particules etape courante                !
! ettpa            ! tr ! <-- ! tableaux des variables liees                   !
!  (nbpmax,nvp)    !    !     !   aux particules etape precedente              !
! tepa             ! tr ! <-- ! info particulaires (reels)                     !
! (nbpmax,nvep)    !    !     !   (poids statistiques,...)                     !
! tempct           ! tr ! <-- ! temps caracteristique thermique                !
!  (nbpmax,2)      !    !     !                                                !
! rayon            ! tr ! <-- ! rayons frontieres des differentes couches      !
!  (nlayer)        !    !     !   (en m) (1 par couche)                        !
! mlayer           ! tr ! <-- ! masse des differentes couches (en kg)          !
!  (nlayer)        !    !     !   (1 par couche)                               !
! phith            ! tr ! <-- ! termes sources thermiques (en W) issu des      !
!  (nlayer)        !    !     !   reaction chimique (1 par couche)             !
! temp             ! tr ! --> ! temperature de la particule apres evolution    !
!  (nlayer)        !    !     !                                                !
! tsvar            ! tr ! <-- ! prediction 1er sous-pas pour la                !
! (nbpmax,nvp1)    !    !     !   variable ivar, utilise pour la               !
!                  !    !     !   correction au 2eme sous-pas                  !
! volume_couche    ! r  ! <-- ! volume occuppe par une couche (en m^3)         !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use numvar
use cstnum
use entsor
use lagpar
use lagran
use ppppar
use ppthch
use cpincl
use radiat
use mesh

!===============================================================================

implicit none

! Arguments
integer          nbpmax , nvp, nvp1, nvep, nivep
integer          itepa(nbpmax,nivep)
integer          npt
double precision propce(ncelet,*)
double precision ettp(nbpmax,nvp) , ettpa(nbpmax,nvp)
double precision tepa(nbpmax,nvep)
double precision tempct(nbpmax,2)
double precision rayon(nlayer), mlayer(nlayer)
double precision phith(nlayer), temp(nlayer)
double precision tsvar(nbpmax,nvp1), volume_couche

! Local variables
integer          ilayer, iel, icha
double precision delray(nlayer-1), rayond(nlayer)
double precision tpscara, coefh, phirayo, temprayo, aux1, aux2
double precision dd2, diamp2
double precision lambda, rho(nlayer), f
double precision a(2:nlayer), b(nlayer), c(1:nlayer-1), d(nlayer)
double precision w1(nlayer-1), w2(nlayer)


!===============================================================================
! INITIALISATION DES VARIABLES
!===============================================================================

iel  = itepa(npt,jisor)
icha = itepa(npt,jinch)

!===============================================================================
! A. RESOLUTION MONOCOUCHE (SI NLAYER > 1)
!===============================================================================
if (nlayer.gt.1) then

  !=============================================================================
  ! A.1. INITIALISATION DES VARIABLES MULTICOUCHES
  !=============================================================================
  ! Calcul des rayons demi et des delta rayon
  do ilayer=1,nlayer
    if (ilayer.eq.nlayer) then
      rayond(ilayer) = (rayon(ilayer-1)+rayon(ilayer))/2.d0
    elseif (ilayer.eq.1) then
      rayond(ilayer) = rayon(ilayer)/2.d0
      delray(ilayer) = rayon(ilayer+1)/2.d0
    else
      rayond(ilayer) = (rayon(ilayer-1)+rayon(ilayer))/2.d0
      delray(ilayer) = (rayon(ilayer+1)-rayon(ilayer-1))/2.d0
    endif
  enddo

  ! Calcul de la masse volumique des couches
  do ilayer=1,nlayer
    rho(ilayer) = mlayer(ilayer)/volume_couche
    if (rho(ilayer).le.0.0d0) then
      write(nfecra,1000) npt,ilayer,rho(ilayer)
      call csexit (1)
    endif
  enddo


  !=============================================================================
  ! A.2. CALCUL DES TRANSFERTS CONDUCTIFS ET RADIATIFS
  !=============================================================================
  ! Calcul des termes sources conducto-convectif et radiatifs

  ! Conduction au sein de la particule
  lambda = thcdch(icha)

  ! Conducto-conductif
  dd2 = ettp(npt,jdp)**2
  diamp2 = xashch(icha)*tepa(npt,jrd0p)*tepa(npt,jrd0p)                        &
           +(1.d0-xashch(icha))*tepa(npt,jrdck)*tepa(npt,jrdck)
  tpscara = tempct(npt,1)*diamp2/dd2
  coefh = ettpa(npt,jmp)*ettpa(npt,jcp)/(tpscara*pi*diamp2)

  ! Temperature equivalente de rayonnement
  temprayo = (propce(iel,ipproc(ilumin))/(4.0d0*stephn))**0.25


  !=============================================================================
  ! A.3. CONSTRUCTION DU SYSTEME
  !=============================================================================
  ! Le phith donné par lagich est en W

  do ilayer=1,nlayer
    if (ilayer.eq.1) then
      b(ilayer) = 1.d0 + 4.0d0*(lambda*dtp)/(rho(ilayer)*ettpa(npt,jcp))       &
                         * ( 1.d0  + 1.d0/(rayon(ilayer+1)*rayon(ilayer))      &
                   + 2.d0/(rayon(ilayer+1)*(rayon(ilayer) + rayon(ilayer+1))) )

      c(ilayer) = - 4.0d0*(lambda*dtp)/(rho(ilayer)*ettpa(npt,jcp))            &
                         * ( 1.d0 + 1.d0/(rayon(ilayer+1)*rayon(ilayer))       &
                   + 2.d0/(rayon(ilayer+1)*(rayon(ilayer) + rayon(ilayer+1))) )
      d(ilayer) = ettp(npt,jhp(ilayer)) +                                      &
                  (phith(ilayer)*dtp)/(mlayer(ilayer)*ettpa(npt,jcp))

    elseif (ilayer.eq.nlayer) then
      f = stephn*( temprayo**2 + (ettp(npt,jhp(nlayer)))**2 )                  &
                *( temprayo    +  ettp(npt,jhp(nlayer))     )
      a(ilayer) = - (lambda*dtp)/(rho(ilayer)*ettpa(npt,jcp)*delray(ilayer-1)) &
                    * ( 1.0d0/delray(ilayer-1)-1.0d0/rayond(ilayer) )
      !b(ilayer) = 1.d0 + (lambda*dtp)/(rho(ilayer)                             &
      !                   * ettpa(npt,jcp)*delray(ilayer-1))                    &
      !                   * ( 1.0d0/delray(ilayer-1)-1.0d0/rayond(ilayer) )
      !d(ilayer) = ettp(npt,jhp(ilayer))+dtp/(mlayer(ilayer)*ettpa(npt,jcp))*   &
      !            (phith(ilayer)+(phicc+phirayo)*volume_couche                 &
      !                           *(1.0d0/delray(ilayer-1)+1.0d0/rayond(ilayer)))
      b(ilayer) = 1.d0 + (lambda*dtp)/(rho(ilayer)                             &
                         * ettpa(npt,jcp)*delray(ilayer-1))                    &
                         * ( 1.0d0/delray(ilayer-1)-1.0d0/rayond(ilayer) )     &
                       + ( dtp*(coefh+f)/(rho(ilayer)*ettpa(npt,jcp))          &
                          * ( 1.0d0/delray(ilayer-1)+1.0d0/rayond(ilayer) ) )
      d(ilayer) = ettp(npt,jhp(ilayer))+dtp/(mlayer(ilayer)*ettpa(npt,jcp))*   &
                  (phith(ilayer) + ( coefh*(ettp(npt,jtf)+tkelvi)+f*temprayo)  &
                                    * volume_couche                            &
                              *(1.0d0/delray(ilayer-1)+1.0d0/rayond(ilayer)) )
    else
      f = (lambda*dtp)/(rho(ilayer)*ettpa(npt,jcp)                             &
                        *delray(ilayer-1)*delray(ilayer))
      a(ilayer) = - f* ( 2.0d0*delray(ilayer)/(delray(ilayer-1)+delray(ilayer))&
                        - (delray(ilayer)/rayond(ilayer)) )
      b(ilayer) = 1.d0 + f* ( 2.0d0 - ((delray(ilayer)-delray(ilayer-1))       &
                                        /rayond(ilayer)) )
      c(ilayer) = -f*( 2.0d0*delray(ilayer-1)/(delray(ilayer-1)+delray(ilayer))&
                        + (delray(ilayer-1)/rayond(ilayer)) )
      d(ilayer) = ettp(npt,jhp(ilayer)) +                                      &
                  (phith(ilayer)*dtp)/(mlayer(ilayer)*ettpa(npt,jcp))
    endif
  enddo


  !=============================================================================
  ! A.4. RESOLUTION DU SYSTEME
  !=============================================================================
  ! On cherche a appliquer l'algorithme de Thomas
  ! a_i T_i-1 + b_i T_i + c_i T_i+1 = d_i
  ! T_i = w2_i - w1_i T_i+1

  do ilayer=1,nlayer
    if (ilayer.eq.1) then
      ! La relation entre T_1 et T_2 est connue
      w1(ilayer) = c(ilayer)/b(ilayer)
      w2(ilayer) = d(ilayer)/b(ilayer)
    elseif (ilayer.eq.nlayer) then
      w2(ilayer) = (d(ilayer)-w2(ilayer-1)*a(ilayer)) /       &
                   (b(ilayer)-w1(ilayer-1)*a(ilayer))
    else
      ! On calcule w1_i et w2_i par reccurrence
      w1(ilayer) = c(ilayer)/(b(ilayer)-w1(ilayer-1)*a(ilayer))
      w2(ilayer) = (d(ilayer)-w2(ilayer-1)*a(ilayer)) /       &
                   (b(ilayer)-w1(ilayer-1)*a(ilayer))
    endif
  enddo

  do ilayer=nlayer,1,-1
    if (ilayer.eq.nlayer) then
      ! La valeur de la temperature en nlayer est connue
      temp(ilayer) = w2(ilayer)
    else
      ! On calcule la temperature en ilayer par reccurrence
      temp(ilayer) = w2(ilayer)-w1(ilayer)*temp(ilayer+1)
    endif
  enddo

!===============================================================================
! B. RESOLUTION MONOCOUCHE (SI NLAYER = 1)
!===============================================================================
else if (nlayer.eq.1) then

  dd2 = ettp(npt,jdp)**2
  diamp2 = xashch(icha)*tepa(npt,jrd0p)*tepa(npt,jrd0p)                        &
           +(1.d0-xashch(icha))*tepa(npt,jrdck)*tepa(npt,jrdck)
  tpscara = tempct(npt,1)*diamp2/dd2

  !     Rayonnement
  phirayo = ( propce(iel,ipproc(ilumin))/4.d0                                  &
              - stephn*((ettp(npt,jhp(nlayer)))**4) )

  aux1 = ettp(npt,jtf)+tkelvi + tpscara*(phirayo*pi*diamp2+phith(nlayer))      &
                                /(ettp(npt,jmp)*ettp(npt,jcp))
  aux2 = exp(-dtp/tpscara)

  if (nor.eq.1) then
    tsvar(npt,jhp(nlayer)) = 0.5d0*ettpa(npt,jhp(nlayer))*aux2                 &
                             + (-aux2+(1.0d0-aux2)*tpscara/dtp)*aux1
    temp(nlayer) = ettpa(npt,jhp(nlayer))*aux2 +                               &
                   (1.0d0-aux2)*aux1

  else if (nor.eq.2) then
    temp(nlayer) = tsvar(npt,jhp(nlayer)) + 0.5d0*ettpa(npt,jhp(nlayer))*aux2  &
                   + (1.0d0-(1.0d0-aux2)*tpscara/dtp)*aux1

  endif

endif

!--------
! FORMATS
!--------

1000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========                                               ',/,&
'@    LAGTMP: UNE COUCHE DE MASSE VOLUMIQUE NULLE A ETE       ',/,&
'@      DETECTEE                                              ',/,&
'@                                                            ',/,&
'@       PARTICULE    = ', I10                                 ,/,&
'@       COUCHE       = ',I10                                  ,/,&
'@       RHO DE LA COUCHE = ',E14.5                            ,/,&
'@                                                            ',/,&
'@  Le calcul de la conduction thermique au sein de la        ',/,&
'@   particule est impossible car le coefficient de diffusion ',/,&
'@   thermique est infini                                     ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier la valeur de xashch dans la subroutine USLAG2    ',/,&
'@  ( xashch > 0 pour éviter ce problème                      ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)


!----
! End
!----

end subroutine
