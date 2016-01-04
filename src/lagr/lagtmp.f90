!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2016 EDF S.A.
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

 ( npt    ,                                                                    &
   propce , tempct ,                                                           &
   rayon  , mlayer , phith , temp  , volume_couche )

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
! npt              ! e  ! <-- ! numero de la particule a traiter               !
! propce(ncelet, *)! tr ! <-- ! physical properties at cell centers            !
! tempct           ! tr ! <-- ! temps caracteristique thermique                !
!  (nbpart,2)      !    !     !                                                !
! rayon            ! tr ! <-- ! rayons frontieres des differentes couches      !
!  (nlayer)        !    !     !   (en m) (1 par couche)                        !
! mlayer           ! tr ! <-- ! masse des differentes couches (en kg)          !
!  (nlayer)        !    !     !   (1 par couche)                               !
! phith            ! tr ! <-- ! termes sources thermiques (en W) issu des      !
!  (nlayer)        !    !     !   reaction chimique (1 par couche)             !
! temp             ! tr ! --> ! temperature de la particule apres evolution    !
!  (nlayer)        !    !     !                                                !
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
integer          npt
double precision propce(ncelet,*)
double precision tempct(nbpart,2)
double precision rayon(nlayer), mlayer(nlayer)
double precision phith(nlayer), temp(nlayer)
double precision volume_couche

! Local variables
integer          ilayer, iel, icha, jtshp0
double precision delray(nlayer-1), rayond(nlayer)
double precision tpscara, coefh, phirayo, temprayo, aux1, aux2
double precision dd2, diamp2
double precision lambda, rho(nlayer), f
double precision a(2:nlayer), b(nlayer), c(1:nlayer-1), d(nlayer)
double precision w1(nlayer-1), w2(nlayer)


!===============================================================================
! INITIALISATION DES VARIABLES
!===============================================================================

iel  = ipepa(jisor,npt)
icha = ipepa(jinch,npt)

if (associated(ptsvar)) then
  jtshp0 = jhp(1) - 1
else
  jtshp0 = -1
endif

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
  dd2 = eptp(jdp,npt)**2
  diamp2 = xashch(icha)*pepa(jrd0p,npt)*pepa(jrd0p,npt)                        &
           +(1.d0-xashch(icha))*pepa(jrdck,npt)*pepa(jrdck,npt)
  tpscara = tempct(npt,1)*diamp2/dd2
  coefh = eptpa(jmp,npt)*eptpa(jcp,npt)/(tpscara*pi*diamp2)

  ! Temperature equivalente de rayonnement
  temprayo = (propce(iel,ipproc(ilumin))/(4.0d0*stephn))**0.25


  !=============================================================================
  ! A.3. CONSTRUCTION DU SYSTEME
  !=============================================================================
  ! Le phith donné par lagich est en W

  do ilayer=1,nlayer
    if (ilayer.eq.1) then
      b(ilayer) = 1.d0 + 4.0d0*(lambda*dtp)/(rho(ilayer)*eptpa(jcp,npt))       &
                         * ( 1.d0  + 1.d0/(rayon(ilayer+1)*rayon(ilayer))      &
                   + 2.d0/(rayon(ilayer+1)*(rayon(ilayer) + rayon(ilayer+1))) )

      c(ilayer) = - 4.0d0*(lambda*dtp)/(rho(ilayer)*eptpa(jcp,npt))            &
                         * ( 1.d0 + 1.d0/(rayon(ilayer+1)*rayon(ilayer))       &
                   + 2.d0/(rayon(ilayer+1)*(rayon(ilayer) + rayon(ilayer+1))) )
      d(ilayer) = eptp(jhp(ilayer),npt) +                                      &
                  (phith(ilayer)*dtp)/(mlayer(ilayer)*eptpa(jcp,npt))

    elseif (ilayer.eq.nlayer) then
      f = stephn*( temprayo**2 + (eptp(jhp(nlayer),npt))**2 )                  &
                *( temprayo    +  eptp(jhp(nlayer),npt)     )
      a(ilayer) = - (lambda*dtp)/(rho(ilayer)*eptpa(jcp,npt)*delray(ilayer-1)) &
                    * ( 1.0d0/delray(ilayer-1)-1.0d0/rayond(ilayer) )
      !b(ilayer) = 1.d0 + (lambda*dtp)/(rho(ilayer)                             &
      !                   * eptpa(jcp,npt)*delray(ilayer-1))                    &
      !                   * ( 1.0d0/delray(ilayer-1)-1.0d0/rayond(ilayer) )
      !d(ilayer) = eptp(jhp(ilayer),npt)+dtp/(mlayer(ilayer)*eptpa(jcp,npt))*   &
      !            (phith(ilayer)+(phicc+phirayo)*volume_couche                 &
      !                           *(1.0d0/delray(ilayer-1)+1.0d0/rayond(ilayer)))
      b(ilayer) = 1.d0 + (lambda*dtp)/(rho(ilayer)                             &
                         * eptpa(jcp,npt)*delray(ilayer-1))                    &
                         * ( 1.0d0/delray(ilayer-1)-1.0d0/rayond(ilayer) )     &
                       + ( dtp*(coefh+f)/(rho(ilayer)*eptpa(jcp,npt))          &
                          * ( 1.0d0/delray(ilayer-1)+1.0d0/rayond(ilayer) ) )
      d(ilayer) = eptp(jhp(ilayer),npt)+dtp/(mlayer(ilayer)*eptpa(jcp,npt))*   &
                  (phith(ilayer) + ( coefh*(eptp(jtf,npt)+tkelvi)+f*temprayo)  &
                                    * volume_couche                            &
                              *(1.0d0/delray(ilayer-1)+1.0d0/rayond(ilayer)) )
    else
      f = (lambda*dtp)/(rho(ilayer)*eptpa(jcp,npt)                             &
                        *delray(ilayer-1)*delray(ilayer))
      a(ilayer) = - f* ( 2.0d0*delray(ilayer)/(delray(ilayer-1)+delray(ilayer))&
                        - (delray(ilayer)/rayond(ilayer)) )
      b(ilayer) = 1.d0 + f* ( 2.0d0 - ((delray(ilayer)-delray(ilayer-1))       &
                                        /rayond(ilayer)) )
      c(ilayer) = -f*( 2.0d0*delray(ilayer-1)/(delray(ilayer-1)+delray(ilayer))&
                        + (delray(ilayer-1)/rayond(ilayer)) )
      d(ilayer) = eptp(jhp(ilayer),npt) +                                      &
                  (phith(ilayer)*dtp)/(mlayer(ilayer)*eptpa(jcp,npt))
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

  dd2 = eptp(jdp,npt)**2
  diamp2 = xashch(icha)*pepa(jrd0p,npt)*pepa(jrd0p,npt)                        &
           +(1.d0-xashch(icha))*pepa(jrdck,npt)*pepa(jrdck,npt)
  tpscara = tempct(npt,1)*diamp2/dd2

  !     Rayonnement
  phirayo = ( propce(iel,ipproc(ilumin))/4.d0                                  &
              - stephn*((eptp(jhp(nlayer),npt))**4) )

  aux1 = eptp(jtf,npt)+tkelvi + tpscara*(phirayo*pi*diamp2+phith(nlayer))      &
                                /(eptp(jmp,npt)*eptp(jcp,npt))
  aux2 = exp(-dtp/tpscara)

  if (nor.eq.1) then
    if (jtshp0.ge.0) then
      ptsvar(jtshp0+ilayer,npt) =    0.5d0*eptpa(jhp(nlayer),npt)*aux2         &
                                  + (-aux2+(1.0d0-aux2)*tpscara/dtp)*aux1
    endif
    temp(nlayer) = eptpa(jhp(nlayer),npt)*aux2 + (1.0d0-aux2)*aux1

  else if (nor.eq.2) then
    temp(nlayer) =   ptsvar(jtshp0+ilayer,npt)                                 &
                   + 0.5d0*eptpa(jhp(nlayer),npt)*aux2                         &
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
