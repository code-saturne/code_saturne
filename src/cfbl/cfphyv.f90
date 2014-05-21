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

subroutine cfphyv &
!================

 ( nvar   , nscal  ,                                              &
   dt     , propce )

!===============================================================================
! FONCTION :
! --------

! ROUTINE PHYSIQUE PARTICULIERE : COMPRESSIBLE SANS CHOC

! Calcul des proprietes physiques variables


! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
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
use numvar
use optcal
use cstphy
use cstnum
use entsor
use ppppar
use ppthch
use ppincl
use mesh

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal

double precision dt(ncelet)
double precision propce(ncelet,*)

! Local variables

integer          iel

integer          ipass
data             ipass /0/
save             ipass

!===============================================================================
!===============================================================================
! 1. INITIALISATIONS A CONSERVER
!===============================================================================

! --- Initialisation memoire



!===============================================================================
! 2. ON DONNE LA MAIN A L'UTILISATEUR
!===============================================================================

call uscfpv                                                       &
!==========
 ( nvar   , nscal  ,                                              &
   dt     )

!===============================================================================
! 3. MISE A JOUR DE LAMBDA/CV
!===============================================================================

! On a vérifié auparavant que CV0 était non nul.
! Si CV variable est nul, c'est une erreur utilisateur. On fait
!     un test à tous les passages (pas optimal), sachant que pour
!     le moment, on est en gaz parfait avec CV constant : si quelqu'un
!     essaye du CV variable, ce serait dommage que cela lui explose à la
!     figure pour de mauvaises raisons.
! Si IVISLS(IENERG).EQ.0, on a forcement IVISLS(ITEMPK).EQ.0
!     et ICV.EQ.0, par construction de IVISLS(IENERG) dans
!     le sous-programme cfvarp

if(ivisls(ienerg).gt.0) then

  if(ivisls(itempk).gt.0) then

    do iel = 1, ncel
      propce(iel,ipproc(ivisls(ienerg))) =               &
           propce(iel,ipproc(ivisls(itempk)))
    enddo

  else
    do iel = 1, ncel
      propce(iel,ipproc(ivisls(ienerg))) =               &
           visls0(itempk)
    enddo

  endif

  if(icv.gt.0) then

    do iel = 1, ncel
      if(propce(iel,ipproc(icv)).le.0.d0) then
        write(nfecra,2000)iel,propce(iel,ipproc(icv))
        call csexit (1)
        !==========
      endif
    enddo

    do iel = 1, ncel
      propce(iel,ipproc(ivisls(ienerg))) =               &
           propce(iel,ipproc(ivisls(ienerg)))            &
           / propce(iel,ipproc(icv))
    enddo

  else

    do iel = 1, ncel
      propce(iel,ipproc(ivisls(ienerg))) =               &
           propce(iel,ipproc(ivisls(ienerg)))            &
           / cv0
    enddo

  endif

else

  visls0(ienerg) = visls0(itempk)/cv0

endif

!--------
! FORMATS
!--------

 2000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION (MODULE COMPRESSIBLE)  ',/,&
'@    =========                                               ',/,&
'@                                                            ',/,&
'@  La capacité calorifique à volume constant présente (au    ',/,&
'@    moins) une valeur négative ou nulle :                   ',/,&
'@    cellule ',I10,   '  Cv = ',E18.9                         ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier uscfpv.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)


!----
! FIN
!----

return
end subroutine
