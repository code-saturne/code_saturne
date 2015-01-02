!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2015 EDF S.A.
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

subroutine lagini &
!================

 ( lndnod )

!===============================================================================
! FONCTION :
! ----------

!   SOUS-PROGRAMME DU MODULE LAGRANGIEN :
!   -------------------------------------

!   Remplissage de LNDNOD pour l'allocation de memoire dans MEMLA1

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! lndnod           ! e  ! --> ! dim. connect. cellules->faces                  !
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
use period
use lagpar
use lagran
use mesh

!===============================================================================

implicit none

! Arguments

integer          lndnod

! Local variables

integer          iel , ifac , ip

integer, allocatable, dimension(:) :: nbrfac

!===============================================================================
!===============================================================================
! 0. GESTION DE LA MEMOIRE
!===============================================================================

! Allocate a temporary array
allocate(nbrfac(ncelet))


!===============================================================================
! 1. Calcul de la dimension du tableau de connectivites cellules->faces
!===============================================================================

!-->Initialisation

do iel = 1,ncelet
  nbrfac(iel) = 0
enddo

!-->Pour chaque cellule on compte le nombre de faces internes qui l'entourent

do ifac = 1,nfac
  do ip = 1,2
    iel = ifacel(ip,ifac)
    nbrfac(iel) = nbrfac(iel) + 1
  enddo
enddo

!-->Pour chaque cellule on compte le nombre de faces de bord qui l'entourent

do ifac = 1,nfabor
  iel = ifabor(ifac)
  nbrfac(iel) = nbrfac(iel) + 1
enddo

!-->verif : chaque cellule ne peut avoir moins de 4 faces en regard
!   Rem : d'apres la routine LEDGEO on a toujours NDIM = 3

ip = 0
do iel = 1,ncel
  if (nbrfac(iel).lt.4) then
    ip = ip + 1
  endif
enddo
if (ip.gt.0) then
  write(nfecra,9000) ip
  call csexit (1)
endif


!-->Calcul de la dimension du tableau de connectivite

lndnod = 0
do iel = 1,ncelet
  lndnod = lndnod + nbrfac(iel)
enddo

! Free memory
deallocate(nbrfac)

!===============================================================================
! 2. Ouverture du fichier listing specifique lagrangien
!===============================================================================

open ( unit=implal, file=ficlal,                                  &
     STATUS='UNKNOWN', FORM='FORMATTED',                          &
     ACCESS='SEQUENTIAL')
rewind(implal)

!--------
! FORMATS
!--------

 9000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========                                               ',/,&
'@  Il y a ',I10,' cellules qui ont moins de 4 faces.         ',/,&
'@   Erreur rencontree dans LAGINI (module Lagrangien).       ',/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier le maillage.                                     ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)


!----
! FIN
!----

end subroutine
