!-------------------------------------------------------------------------------

!VERS

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

subroutine usthht &
!================

 ( mode   , enthal , temper  )

!===============================================================================
! FONCTION :
! --------

!    ROUTINE UTILISATEUR
!    LOI ENTHALPIE   -> TEMPERATURE (MODE =  1)
!    LOI TEMPERATURE -> ENTHALPIE   (MODE = -1)

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! mode             ! e  ! <-- !  -1 : t -> h  ;   1 : h -> t                   !
! enthal           ! r  ! <-- ! enthalpie                                      !
! temper           ! r  ! <-- ! temperature                                    !
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use entsor
use parall
use period

!===============================================================================

implicit none

! Arguments

integer          mode

double precision enthal, temper

! Local variables

integer         it

!     Pour le second exemple ci-dessous,
!        donnee de NTAB > 1 valeurs (fictives) de H (=HT) tabulees
!        en fonction de NTAB valeurs de T (=TH) (attention a l'unite K ou C)

integer          ntab
parameter       (ntab=5)
double precision ht(ntab), th(ntab)
data             ht /100000.d0,200000.d0,300000.d0,               &
                               400000.d0,500000.d0 /
data             th /   100.d0,   200.d0,   300.d0,               &
                                  400.d0,   500.d0 /

!===============================================================================


!     ATTENTION, CE SOUS-PROGRAMME EST APPELE DANS DES BOUCLES :
!     =========                               ================

!       EVITER LES IMPRESSIONS
!       ======

!       EVITER LES OPERATIONS FAISANT APPEL AU PARALLELISME
!       ======

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START

write(nfecra,9000)
call csexit (1)
!==========

!----
! FORMATS
!----

 9000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET LORS DU CALCUL DE LA TEMPERATURE      ',/,&
'@    =========                                               ',/,&
'@    LES TABLES ENTHALPIE TEMPERATURE NE SONT PAS DISPONIBLES',/,&
'@                                                            ',/,&
'@  Le sous-programme utilisateur usthht doit etre complete.  ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Le couplage avec SYRTHES necessite la donne d''une        ',/,&
'@    temperature de paroi.                                   ',/,&
'@  Le scalaire choisi pour le couplage SYRTHES est ici une   ',/,&
'@    enthalpie.                                              ',/,&
'@  La loi donnant la temperature en fonction de l''enthalpie ',/,&
'@    doit etre fournie par l''utilisateur dans le            ',/,&
'@    sous-programme usthht.                                  ',/ &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END

!===============================================================================
! EXEMPLES
!===============================================================================

! Premier exemple, correspondant a H=CpT avec Cp = 4000
! =====================================================

! --- Mode H -> T
if (mode .eq.  1) then
  temper = enthal / 4000.d0

! --- Mode T -> H
else
  enthal = temper * 4000.d0

endif

return












! Second exemple, correspondant a une interpolation simple
! === a partir d'une tabulation H=f(T) entree en DATA ====
!     ================================================

! --- Mode H -> T
if (mode .eq.  1) then

!       Initialisation par defaut
  temper = 0.d0

!       Si H plus petit que la plus petite enthalpie tabulee
!         on limite arbitrairement T a la plus petite temperature
  if(enthal.le.ht(1)) then
    temper = th(1)

!       Si H plus grand que la plus grande enthalpie tabulee
!         on limite arbitrairement T a la plus grande temperature
  elseif(enthal.ge.ht(ntab)) then
    temper = th(ntab)

!       Sinon, on interpole lineairement
  else
    do it = 2, ntab
      if(enthal.le.ht(it)) then
        temper = th(it-1)                                         &
          +(enthal-ht(it-1))*(th(it)-th(it-1))/(ht(it)-ht(it-1))
      endif
    enddo
  endif

! --- Mode T -> H
else

!       Initialisation par defaut
  enthal = 0.d0

!       Si T plus petit que la plus petite temperature tabulee
!         on limite arbitrairement H a la plus petite enthalpie
  if(temper.le.th(1)) then
    enthal = ht(1)

!       Si T plus grand que la plus grande temperature tabulee
!         on limite arbitrairement H a la plus grande enthalpie
  elseif(temper.ge.th(ntab)) then
    enthal = ht(ntab)

!       Sinon, on interpole lineairement
  else
    do it = 2, ntab
      if(temper.le.th(it)) then
        enthal = ht(it-1)                                         &
          +(temper-th(it-1))*(ht(it)-ht(it-1))/(th(it)-th(it-1))
      endif
    enddo
  endif
endif

return

!----
! FIN
!----

end subroutine usthht
