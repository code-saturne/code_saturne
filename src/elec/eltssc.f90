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

subroutine eltssc &
!================

 ( iscal  ,                                                       &
   propce ,                                                       &
   smbrs  )

!===============================================================================
! FONCTION :
! ----------

! ROUTINE PHYSIQUE PARTICULIERE POUR LE MODULE ELECTRIQUE
! CALCUL DES TERMES SOURCE ELECTRIQUES POUR L'ENTHALPIE
! CALCUL DES TERMES SOURCES POUR LE POTENTIEL VECTEUR

! ATTENTION : LE TRAITEMENT DES TERMES SOURCES EST DIFFERENT
! ---------   DE CELUI DE USTSSC

!-------------------------------------------------------------------------------
!ARGU                             ARGUMENTS
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! iscal            ! i  ! <-- ! scalar number                                  !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! smbrs(ncelet)    ! tr ! --> ! second membre explicite                        !
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
use entsor
use optcal
use cstphy
use cstnum
use parall
use period
use ppppar
use ppthch
use ppincl
use elincl
use mesh
use field

!===============================================================================

implicit none

! Arguments

integer          iscal

double precision propce(ncelet,*)
double precision smbrs(ncelet)

! Local variables

character*80     chaine
integer          ivar , iel
integer          ipcdc1, ipcdc2, ipcdc3
integer          ipcefj

double precision valmin,valmax

double precision, allocatable, dimension(:) :: w1

!===============================================================================
!===============================================================================
! 1. INITIALISATION
!===============================================================================

! Allocate a temporary array
allocate(w1(ncelet))


! --- Numero du scalaire a traiter : ISCAL

! --- Numero de la variable associee au scalaire a traiter ISCAL
ivar = isca(iscal)

! --- Nom de la variable associee au scalaire a traiter ISCAL
call field_get_label(ivarfl(ivar), chaine)

!===============================================================================
! 2. PRISE EN COMPTE DES TERMES SOURCES ET VARIABLES STANDARD ET
!   COMMUNES A TOUTES LES VERSIONS ELECTRIQUES : EFFET JOULE, E, J
!         CAS IELJOU >= 1 ou IELARC >= 1 ou IELION >= 1
!===============================================================================

!   2.1 Terme source pour l'enthalpie :
!  ----------------------------------

if ( ivar.eq.isca(iscalt) ) then

  if (iwarni(ivar).ge.1) then
    write(nfecra,1000) chaine(1:8)
  endif

  ipcefj = ipproc(iefjou)

!     Pour eviter de prendre en compte des TS trop violents
!     sur les premiers pas de temps (gradient sur le potentiel)
!     on ne prend en compte les termes sources qu'a partir du
!     troisieme pas de temps. Sinon on laisse ROVSDT et SMBRS
!     inchanges (pas de TS si ce n'est ceux de ustssc)

  if ( ntcabs .ge. 3 ) then
    do iel = 1, ncel
      w1(iel) =  propce(iel,ipcefj)*volume(iel)
    enddo

!        Arc electrique : rajout du TS radiatif si l'utilisateur
!        l'a donne dans le fichier dp_elec (purement explicite)

    if( ippmod(ielarc).ge. 1 ) then
      if ( ixkabe .eq. 2 ) then
        do iel = 1, ncel
          w1(iel) = w1(iel)                                       &
                       -propce(iel,ipproc(idrad))*volume(iel)
        enddo
      endif
    endif

    do iel = 1, ncel
       smbrs(iel) = smbrs(iel) + w1(iel)
    enddo

    if (iwarni(ivar).ge.2) then
      valmin = w1(1)
      valmax = w1(1)
      do iel = 1, ncel
        valmin = min(valmin,w1(iel))
        valmax = max(valmax,w1(iel))
      enddo
      if (irangp.ge.0) then
        call parmax(valmax)
        call parmin(valmin)
      endif
      write(nfecra,2000) valmin,valmax
    endif
  endif
endif

 2000  format(                                                          &
 ' Termes Sources sur H  min= ',E14.5,', max= ',E14.5)

!   2.1  Terme source pour les composantes potentiel vecteur : ARC ELECTRIQUE
!  ---------------------------------------------------------

if( ippmod(ielarc).ge. 2 ) then


! -->   Terme source pour les composantes potentiel vecteur
!============================================================

  if ( ivar .eq. isca(ipotva(1)) ) then

    if (iwarni(ivar).ge.1) then
      write(nfecra,1000) chaine(1:8)
    endif

    ipcdc1 = ipproc(idjr(1))
    do iel = 1, ncel
      smbrs(iel) = smbrs(iel) +                                   &
            permvi*propce(iel,ipcdc1)*volume(iel)
    enddo

  else if ( ivar .eq. isca(ipotva(2)) ) then

    if (iwarni(ivar).ge.1) then
      write(nfecra,1000) chaine(1:8)
    endif

    ipcdc2 = ipproc(idjr(2))
    do iel = 1, ncel
      smbrs(iel) = smbrs(iel) +                                   &
            permvi*propce(iel,ipcdc2)*volume(iel)
    enddo

  else if ( ivar .eq. isca(ipotva(3)) ) then

    if (iwarni(ivar).ge.1) then
      write(nfecra,1000) chaine(1:8)
    endif

    ipcdc3 = ipproc(idjr(3))
    do iel = 1, ncel
      smbrs(iel) = smbrs(iel) +                                   &
            permvi*propce(iel,ipcdc3)*volume(iel)
    enddo
  endif

endif

! Free memory
deallocate(w1)

!--------
! FORMATS
!--------

 1000 format(                                                           &
'  Calcul des termes sources pour la variable : ',A8             )



!----
! FIN
!----

return
end subroutine
