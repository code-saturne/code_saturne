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

subroutine ledgeo &
!================

     ( ndim   , ncelet , ncel   , nfac   , nfabor ,               &
       nprfml , nfml   , nsom   , lndfac , lndfbr )

!===============================================================================
! FONCTION :
! ---------

! LECTURE DES DIMENSIONS DES TABLEAUX ENTITES GEOMETRIQUES

!-------------------------------------------------------------------------------
!ARGU                             ARGUMENTS
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! ndim             ! e  ! --> ! dimension de l'espace (=3)                     !
! ncelet           ! e  ! --> ! nombre d'elements halo compris                 !
! ncel             ! e  ! --> ! nombre d'elements actifs                       !
! nfac             ! e  ! --> ! nombre de faces internes                       !
! nfabor           ! e  ! --> ! nombre de faces de bord                        !
! nprfml           ! e  ! --> ! nombre de propietes des familles               !
!                  !    !     ! de faces de bord                               !
! nfml             ! e  ! --> ! nombre de familles de faces de bord            !
! nsom             ! e  ! --> ! nombre de sommets du maillage                  !
! lndfac           ! e  ! --> ! longueur du tableau nodfac (optionnel          !
! lndfbr           ! e  ! --> ! longueur du tableau nodfbr (optionnel          !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

implicit none

!===============================================================================
! Common blocks
!===============================================================================

include "paramx.h"
include "entsor.h"

!===============================================================================

! Arguments

integer          ndim   , ncelet , ncel   , nfac   , nfabor
integer          nprfml , nfml   , nsom   , lndfac , lndfbr

! Local variables

integer          iel , ifac , ifml , nn , n1 , ip
integer          ntetra , npyram , nprism , nhexae

!===============================================================================

ndim = 0
ncel = 0
nfac = 0
nfabor = 0
nprfml = 0
nfml   = 0
nsom   = 0

lndfac = 0
lndfbr = 0

!===============================================================================
! 1.  Ouverture de FICGEO et lecture des entetes principales
!===============================================================================

!---> Ouverture

OPEN (FILE=FICGEO,UNIT=IMPGEO,FORM='formatted',ERR=9000)

!---> Ecriture avec rubrique ... avant

ndim = 3

read (impgeo,   *,err=9002,end=9003)
read (impgeo,   *,err=9002,end=9003)
read (impgeo,1100,err=9002,end=9003) ncel, nfac, nfabor, nsom

read (impgeo,   *,err=9002,end=9003)
read (impgeo,   *,err=9002,end=9003)
read (impgeo,1100,err=9002,end=9003) ntetra,npyram,nprism,nhexae

read (impgeo,   *,err=9002,end=9003)
read (impgeo,   *,err=9002,end=9003)
read (impgeo,1100,err=9002,end=9003) nprfml, nfml

!===============================================================================
! 2.  Positionnement dans FICGEO pour la lecture des tableaux optionnels
!     de connectivites faces sommets (IPNFAC, IPNFBR, NODFAC, NODFBR)
!===============================================================================

!-->  Positionnement dans le fichier avant lecture

!     Connectivite faces internes - cellules
read(impgeo,*)
read(impgeo,*)
do ifac = 1,nfac
  read(impgeo,*)
enddo

!     Connectivite faces de bord - cellule
read(impgeo,*)
read(impgeo,*)
do ifac = 1, nfabor
  read(impgeo,*)
enddo

!     Coordonnees du centre des cellules
read(impgeo,*)
read(impgeo,*)
do iel = 1,ncel
  read(impgeo,*)
enddo

!     Surfaces des faces internes
read(impgeo,*)
read(impgeo,*)
do ifac = 1,nfac
  read(impgeo,*)
enddo

!     Surfaces des faces de bord
read(impgeo,*)
read(impgeo,*)
do ifac = 1,nfabor
  read(impgeo,*)
enddo

!     Coordonnees du centre des faces
read(impgeo,*)
read(impgeo,*)
do ifac = 1,nfac
  read(impgeo,*)
enddo

!     Coordonnees du centre des faces de bord
read(impgeo,*)
read(impgeo,*)
do ifac = 1,nfabor
  read(impgeo,*)
enddo

!     Familles des faces de bord
read(impgeo,*)
read(impgeo,*)
do ifac = 1,nfabor
  read(impgeo,*)
enddo

!     Proprietes des familles
read(impgeo,*)
read(impgeo,*)
do ifml = 1,nfml
  read(impgeo,*)
enddo

!     Coordonnees des noeuds
read(impgeo,*)
read(impgeo,*)
do n1 = 1,nsom
  read(impgeo,*)
enddo

!     Connectivite cellules points
read(impgeo,*)
read(impgeo,*)
do n1 = 1,ncel
  read(impgeo,*)
enddo

!===============================================================================
! 3.  Calcul des dimensions LNDFAC et LNDFBR si possible
!===============================================================================

!     Remarque : si les tableaux ne sont pas disponibles,
!                on aura une erreur de fin de fichier.
!                Dans ce cas, on sort en fermant normalement
!                le fichier, en mettant LNDFAC et LNDFBR à 0.

!-->  Dimension du tableau de connectivite faces de bord -> points

read(impgeo,*,err=9000,end=100)
read(impgeo,*,err=9000,end=100)
lndfbr = 0
do n1 = 1,nfabor
  read(impgeo,1200,err=9000,end=100) nn , ip
  lndfbr = lndfbr + ip
enddo

!-->  Dimension du tableau de connectivite faces internes -> points

read(impgeo,*,err=9000,end=100)
read(impgeo,*,err=9000,end=100)
lndfac = 0
do n1 = 1,nfac
  read(impgeo,1200,err=9000,end=100) nn , ip
  lndfac = lndfac + ip
enddo

!-->  Si l'on a pu lire les tableaux de connectivite faces -> points,
!     on conserve les dimensions LNDFBR et LNDFAC calculees ;
!     sinon, on les remet a zero

goto 200

 100  lndfbr = 0
lndfac = 0

 200  continue

!===============================================================================
! 3.  Fermeture du fichier et mise a jour des structures C
!===============================================================================

close(impgeo)

!---> MISE A JOUR STRUCTURES C

ncelet = ncel

call dimgeo                                                       &
!==========
     (ndim  , ncelet, ncel  , nfac  , nfabor , nsom  ,            &
      lndfac, lndfbr, nfml  , nprfml,                             &
      ntetra, npyram, nprism, nhexae)

!---> FORMATS

 1100 format(20i10)
 1200 format(i10,i10)

return

! ERREURS

 9000 continue
write(nfecra,8000)ficgeo
call csexit (1)
 9002 continue
write(nfecra,8002)                                                &
            ficgeo,ndim,ncel,nfac,nfabor,nprfml,nfml  ,nsom
call csexit (1)
 9003 continue
write(nfecra,8003)                                                &
            ficgeo,ndim,ncel,nfac,nfabor,nprfml,nfml  ,nsom
call csexit (1)

! FORMATS

#if defined(_CS_LANG_FR)

 8000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''OUVERTURE DU FICHIER GEOMETRIE   ',/,&
'@    =========                                               ',/,&
'@      ERREUR DANS ledgeo POUR LE FICHIER ',A6                ,/,&
'@                                                            ',/,&
'@      Le calcul ne peut etre execute.                       ',/,&
'@                                                            ',/,&
'@      Verifier le fichier geometrie (existence, droits      ',/,&
'@        d''acces, recopie correcte dans le repertoire de    ',/,&
'@        travail).                                           ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8002 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER GEOMETRIE     ',/,&
'@    =========                                               ',/,&
'@      ERREUR DANS ledgeo POUR LE FICHIER ',A6                ,/,&
'@                                                            ',/,&
'@      L''etat des dimensions lues est le suivant :          ',/,&
'@        NDIM      NCEL      NFAC    NFABOR                  ',/,&
'@  ',4I10                                                     ,/,&
'@        NPRFML    NFML      NSOM                            ',/,&
'@  ',3I10                                                     ,/,&
'@                                                            ',/,&
'@      Le calcul ne peut etre execute.                       ',/,&
'@                                                            ',/,&
'@      Verifier le fichier geometrie (existence, droits      ',/,&
'@        d''acces, format...).                               ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8003 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER GEOMETRIE     ',/,&
'@    =========                                               ',/,&
'@      FIN PREMATUREE DANS ledgeo POUR LE FICHIER ',A6        ,/,&
'@                                                            ',/,&
'@      L''etat des dimensions lues est le suivant :          ',/,&
'@        NDIM      NCEL      NFAC    NFABOR                  ',/,&
'@  ',4I10                                                     ,/,&
'@        NPRFML    NFML      NSOM                            ',/,&
'@  ',3I10                                                     ,/,&
'@                                                            ',/,&
'@      Le calcul ne peut etre execute.                       ',/,&
'@                                                            ',/,&
'@      Verifier le fichier geometrie (existence, droits      ',/,&
'@        d''acces, format...).                               ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

#else

 8000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: ABORT WHILE OPENING THE GEOMETRY FILE          ',/,&
'@    ========                                                ',/,&
'@      ERROR IN ledgeo FOR FILE ',A6                          ,/,&
'@                                                            ',/,&
'@      The calculation will not be run.                      ',/,&
'@                                                            ',/,&
'@      Verify the geometry file (existence, access           ',/,&
'@        permission, correct copy in the working directory)  ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8002 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: ABORT WHILE READING THE GEOMETRY FILE          ',/,&
'@    ========                                                ',/,&
'@      ERROR IN ledgeo FOR FILE ',A6                          ,/,&
'@                                                            ',/,&
'@      The state of the read dimensions is:                  ',/,&
'@        NDIM      NCEL      NFAC    NFABOR                  ',/,&
'@  ',4I10                                                     ,/,&
'@        NPRFML    NFML      NSOM                            ',/,&
'@  ',3I10                                                     ,/,&
'@                                                            ',/,&
'@      The calculation will not be run.                      ',/,&
'@                                                            ',/,&
'@      Verify the geometry file (existence, access           ',/,&
'@        permission, format...)                              ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8003 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: ABORT WHILE READING THE GEOMETRY FILE          ',/,&
'@    ========                                                ',/,&
'@      UNEXPECTED END OF FILE in ledgeo FOR FILE ',A6         ,/,&
'@                                                            ',/,&
'@      The state of the read dimensions is:                  ',/,&
'@        NDIM      NCEL      NFAC    NFABOR                  ',/,&
'@  ',4I10                                                     ,/,&
'@        NPRFML    NFML      NSOM                            ',/,&
'@  ',3I10                                                     ,/,&
'@                                                            ',/,&
'@      The calculation will not be run.                      ',/,&
'@                                                            ',/,&
'@      Verify the geometry file (existence, access           ',/,&
'@        permission, format...)                              ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

#endif

end subroutine
