!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2019 EDF S.A.
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

!!===============================================================================
! Function:
! ---------

!> \file vortex.f90
!>
!> \brief L.E.S Vortex method data management.
!>
!> Please refer to the <a href="../../theory.pdf#vortex"><b>vortex</b></a>
!> section of the theory guide for more informations.
!-------------------------------------------------------------------------------

subroutine vortex

!===============================================================================
! Module files
!===============================================================================

use paramx
use entsor
use optcal
use vorinc

!===============================================================================

implicit none

! Arguments

! Local variables

character        ficsui*32

integer          ii, ient
integer          ipass
data             ipass /0/
save             ipass

!===============================================================================
! 1. INITIALISATION
!===============================================================================


! L'EQUATION DE LANGEVIN RESTE A TRAVAILLER POUR UN CAS 3D QUELCONQUE
! OU A MODIFIER    . VU LE PEU D'IMPORTANCE QU'ELLE A SUR CE QUI SE PASSE
! EN AVAL DE L'ENTREE (L'IMPORTANT ETANT D'IMPOSER V' ET W'), ON NE VA
! PAS PLUS LOIN (ON ANNULE CES CONTRIBUTION POUR LE MOMENT POUR LES
! CAS 3 ET 4)

ipass = ipass + 1

do ient = 1, nnent

  if (ipass.eq.1)then

    call vorini &
    !==========
 ( icvor(ient)     , nvort(ient)     ,                            &
   ient   , ivorce(1,ient)  ,                                     &
   xyzv(1,1,ient)  , yzcel(1,1,ient) ,                            &
   uvort(1,ient)   ,                                              &
   yzvor(1,1,ient) , signv(1,ient)   , temps(1,ient)   ,          &
   tpslim(1,ient)  )

  endif

!===============================================================================
! 2. DEPLACEMENT DU VORTEX
!===============================================================================

  call vordep &
  !==========
 ( icvor(ient)     , nvort(ient)     , ient   , dtref  ,          &
   ivorce(1,ient)  , yzcel(1,1,ient) ,                            &
   vvort(1,ient)   , wvort(1,ient)   ,                            &
   yzvor(1,1,ient) , yzvora(1,1,ient), signv(1,ient)   ,          &
   temps(1,ient)   , tpslim(1,ient)  )

!===============================================================================
! 3. CALCUL DE LA VITESSE
!===============================================================================

  call vorvit                                                     &
  !==========
 ( icvor(ient)     , nvort(ient)     , ient   ,                   &
   ivorce(1,ient)  , visv(1,ient)    ,                            &
   yzcel(1,1,ient) , vvort(1,ient)   , wvort(1,ient)   ,          &
   yzvor(1,1,ient) , signv(1,ient)   ,                            &
   sigma(1,ient)   , gamma(1,1,ient) )

!===============================================================================
! 4. CALCUL DES FLUCTUATIONS DANS LE SENS DE L'ECOULEMENT
!===============================================================================

  call vorlgv                                                     &
  !==========
 ( icvor(ient)     , ient   , dtref  ,                            &
   yzcel(1,1,ient) ,                                              &
   uvort(1,ient)   , vvort(1,ient)   , wvort(1,ient)   )

enddo

!===============================================================================
! 5. ECRITURE DU FICHIER SUITE
!===============================================================================

! on ecrit a tous les pas de temps pour eviter
! les mauvaises surprises en cas de fin prematuree.
! Il ne faut pas mettre cette partie dans la boucle sur IENT
! car on accede deja a l'unite IMPMVO(=IMPDVO) dans VORINI.
! Seul le premier processeur ecrit (test avant l'appel à VORTEX)

ficsui = 'checkpoint/vortex'
open(unit=impvvo,file=ficsui)
rewind(impvvo)
do ient = 1, nnent
  write(impvvo,100) ient
  write(impvvo,100) nvort(ient)
  do ii = 1, nvort(ient)
    write(impvvo,200) yzvor(ii,1,ient),yzvor(ii,2,ient),          &
         temps(ii,ient), tpslim(ii,ient), signv(ii,ient)
  enddo
enddo
close(impvvo)

!===============================================================================
! 6. FIN
!===============================================================================

! FORMATS

 100  format(i10)
 200  format(5e13.5)

return

end subroutine
