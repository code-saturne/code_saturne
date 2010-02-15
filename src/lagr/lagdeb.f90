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

subroutine lagdeb &
!================

 ( idbia0 , idbra0 ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndnod , lndfac , lndfbr , ncelbr ,                   &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   icocel , itycel ,                                              &
   idevel , ituser , ia     ,                                     &
   rdevel , rtuser , ra     )

!===============================================================================
! FONCTION :
! ----------

!   SOUS-PROGRAMME DU MODULE LAGRANGIEN :
!   -------------------------------------

! Ce sous-programme est lu lors du premier passage dans le module
! Lagrangien (IPLAR=1).

! Construction du stockage compact (tableaux "positions"/"valeurs")
! des connectivites cellules -> faces
! (pour tous formats de maillage).
! On part des connectivites faces->cellule IFACEL et IFABOR.
! Les faces de bord sont stockees avec un numero negatif.


! On rappelle que la lecture des connectivites face->sommets
! pour un maillage .slc est faite dans LAGGET, appelee par LETGEO.



! ----------------------------------------------------------------
! Schema d'association entre les tableaux ITYCEL "positions"
! et ICOCEL "valeurs"
! ----------------------------------------------------------------


!         Tableau des valeurs ICOCEL (de dimension LNDNOD)


! .---.------..------.------.------..------.------..------.------.
! |   | ...  ||IFACn | ...  |      || IFACm| ...  || ...  |      |
! `---'------'`------'------'------'`------'------'`------'------'
!  1           iVal         jVal-1   jVal                  LNDNOD
!                |                     |                      |
!                |                     |                      |
!                `----------.      .---'         .------------'
!                           |      |             |
!                           |      |             |
!        .-----.-------.------.------.-------.------.--------.
!        |  1  |  ...  | iVal | jVal |  ...  |LNDNOD|LNDNOD+1|
!        `-----'-------'------'------'-------'------'--------'
!           1            iPos  iPos+1         nPos-1 nPos=NCELET+1


!         Tableau des positions ITYCEL (de dimension NCELET+1)


!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
! ndim             ! i  ! <-- ! spatial dimension                              !
! ncelet           ! i  ! <-- ! number of extended (real + ghost) cells        !
! ncel             ! i  ! <-- ! number of cells                                !
! nfac             ! i  ! <-- ! number of interior faces                       !
! nfabor           ! i  ! <-- ! number of boundary faces                       !
! nfml             ! i  ! <-- ! number of families (group classes)             !
! nprfml           ! i  ! <-- ! number of properties per family (group class)  !
! nnod             ! i  ! <-- ! number of vertices                             !
! lndnod           ! e  ! <-- ! dim. connectivite cellules->faces              !
! lndfac           ! e  ! <-- ! longueur du tableau nodfac                     !
! lndfbr           ! e  ! <-- ! longueur du tableau nodfbr                     !
! ncelbr           ! i  ! <-- ! number of cells with faces on boundary         !
! nideve, nrdeve   ! i  ! <-- ! sizes of idevel and rdevel arrays              !
! nituse, nrtuse   ! i  ! <-- ! sizes of ituser and rtuser arrays              !
! ifacel(2, nfac)  ! ia ! <-- ! interior faces -> cells connectivity           !
! ifabor(nfabor)   ! ia ! <-- ! boundary faces -> cells connectivity           !
! ifmfbr(nfabor)   ! ia ! <-- ! boundary face family numbers                   !
! ifmcel(ncelet)   ! ia ! <-- ! cell family numbers                            !
! iprfml           ! te ! <-- ! proprietes d'une famille                       !
!  (nfml,nprfml    !    !     !                                                !
! ipnfac           ! te ! <-- ! position du premier noeud de chaque            !
!   (lndfac)       !    !     !  face interne dans nodfac                      !
! nodfac           ! te ! <-- ! connectivite faces internes/noeuds             !
!   (nfac+1)       !    !     !                                                !
! ipnfbr           ! te ! <-- ! position du premier noeud de chaque            !
!   (lndfbr)       !    !     !  face de bord dans nodfbr                      !
! nodfbr           ! te ! <-- ! connectivite faces de bord/noeuds              !
!   (nfabor+1)     !    !     !                                                !
! icocel           ! te ! --> ! connectivite cellules -> faces                 !
! (lndnod)         !    !     !    face de bord si numero negatif              !
! itycel           ! te ! --> ! connectivite cellules -> faces                 !
! (ncelet+1)       !    !     !    pointeur du tableau icocel                  !
! idevel(nideve)   ! ia ! <-> ! integer work array for temporary development   !
! ituser(nituse)   ! ia ! <-> ! user-reserved integer work array               !
! ia(*)            ! ia ! --- ! main integer work array                        !
! rdevel(nrdeve)   ! ra ! <-> ! real work array for temporary development      !
! rtuser(nrtuse)   ! ra ! <-> ! user-reserved real work array                  !
! ra(*)            ! ra ! --- ! main real work array                           !
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
include "numvar.h"
include "optcal.h"
include "entsor.h"
include "pointe.h"
include "period.h"
include "parall.h"
include "lagpar.h"
include "lagran.h"

!===============================================================================

! Arguments

integer          idbia0 , idbra0
integer          ndim   , ncelet , ncel   , nfac   , nfabor
integer          nfml   , nprfml
integer          nnod   , lndnod , lndfac , lndfbr , ncelbr
integer          nideve , nrdeve , nituse , nrtuse
integer          ifacel(2,nfac) , ifabor(nfabor)
integer          ifmfbr(nfabor) , ifmcel(ncelet)
integer          iprfml(nfml,nprfml)
integer          ipnfac(nfac+1) , nodfac(lndfac)
integer          ipnfbr(nfabor+1) , nodfbr(lndfbr)
integer          icocel(lndnod) , itycel(ncelet+1)
integer          idevel(nideve) , ituser(nituse)
integer          ia(*)

double precision rdevel(nrdeve) , rtuser(nrtuse) , ra(*)


! Local variables

integer          idebia, idebra

integer          iel , ifac , ip , n1
integer          nbfac1 , nbfac2

!===============================================================================

!===============================================================================
! 0.  GESTION MEMOIRE
!===============================================================================

idebia = idbia0
idebra = idbra0

!===============================================================================
! 1. Connectivite cellules -> faces
!    Les faces de bord sont stockees avec des valeurs negatives
!===============================================================================

!-->Initialisation

do n1 = 1,ncelet+1
  itycel(n1) = 0
enddo
do n1 = 1,lndnod
  icocel(n1) = 0
enddo

!-->Premier passage : on compte le nombre de faces par cellules

do ifac = 1,nfac
  do ip = 1,2
    iel = ifacel(ip,ifac)
    itycel(iel) = itycel(iel) + 1
  enddo
enddo

do ifac = 1,nfabor
  iel = ifabor(ifac)
  itycel(iel) = itycel(iel) + 1
enddo

!-->Verification

ip = 0
do iel = 1,ncelet
  ip = ip + itycel(iel)
enddo
if (ip.ne.lndnod) then
  write(nfecra,5000) lndnod , ip
  call csexit (1)
endif

!-->Deuxieme passage : on calcule les positions du tableau ICOCEL

nbfac2 = itycel(1)
itycel(1) = 1
do iel = 2, ncelet+1
  nbfac1 = nbfac2
  nbfac2 = itycel(iel)
  itycel(iel) = itycel(iel-1) + nbfac1
enddo

!-->Verification

if (itycel(ncelet+1).ne.lndnod+1) then
  write(nfecra,1000) itycel(ncelet+1) , lndnod+1
  call csexit (1)
endif

!-->Troisieme passage : on remplit le tableau de connectivites ICOCEL
!   Rem : les faces de bord sont stockees avec des valeurs negatives

do ifac = 1,nfac
  do ip = 1,2
    iel = ifacel(ip,ifac)
    n1 = itycel(iel)
    do while (icocel(n1).ne.0)
      n1 = n1 + 1
      if (n1.eq.itycel(iel+1)) then
        write(nfecra,2000)
        call csexit (1)
      endif
    enddo
    icocel(n1) = ifac
  enddo
enddo

do ifac = 1,nfabor
  iel = ifabor(ifac)
  n1 = itycel(iel)
  do while (icocel(n1).ne.0)
    n1 = n1 + 1
    if (n1.eq.itycel(iel+1)) then
      write(nfecra,3000)
      call csexit (1)
    endif
  enddo
  icocel(n1) = -ifac
enddo

!-->Archive : ancienne version du 3e passage NE PAS SUPPRIMER SVP
!   Attention ne marche pas dans le cas de la TURBINE,
!   peut-etre a cause du fait que deux cellules ont en commun
!   plusieurs faces

!     N1 = 0
!     DO IEL = 1,NCEL
!       DO IFAC = 1,NFAC
!         DO IP = 1,2
!           IF (IFACEL(IP,IFAC).EQ.IEL) THEN
!             N1 = N1 + 1
!             ICOCEL(N1) = ICOCEL(N1) - IFAC
!           ENDIF
!         ENDDO
!       ENDDO
!       DO IFAC = 1,NFABOR
!         IF (IFABOR(IFAC).EQ.IEL) THEN
!           N1 = N1 + 1
!           ICOCEL(N1) = ICOCEL(N1)  - IFAC
!         ENDIF
!       ENDDO
!     ENDDO


!-->Verification

do n1 = 1,lndnod
  if (icocel(n1).eq.0) then
    write(nfecra,4000)
    call csexit (1)
  endif
enddo

!===============================================================================

!--------
! FORMATS
!--------

 1000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========                                               ',/,&
'@    Erreur dans la construction de la connectivite          ',/,&
'@    cellules -> faces (LAGDEB)                              ',/,&
'@                                                            ',/,&
'@    Mauvais remplissage de ITYCEL : le dernier              ',/,&
'@    element de ITYCEL ne vaut pas LNDNOD+1                  ',/,&
'@                                                            ',/,&
'@    ITYCEL(NCELET+1) = ',I10                                 ,/,&
'@    LNDNOD+1         = ',I10                                 ,/,&
'@                                                            ',/,&
'@    Verifier le remplissage de ITYCEL dans LAGDEB           ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 2000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========                                               ',/,&
'@    Erreur dans la construction de la connectivite          ',/,&
'@    cellules -> faces (LAGDEB)                              ',/,&
'@                                                            ',/,&
'@    Mauvais remplissage de ITYCEL avec les faces internes   ',/,&
'@                                                            ',/,&
'@    Verifier le remplissage de ITYCEL dans LAGDEB           ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 3000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========                                               ',/,&
'@    Erreur dans la construction de la connectivite          ',/,&
'@    cellules -> faces (LAGDEB)                              ',/,&
'@                                                            ',/,&
'@    Mauvais remplissage de ITYCEL avec les faces de bord    ',/,&
'@                                                            ',/,&
'@    Verifier le remplissage de ITYCEL dans LAGDEB           ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 4000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========                                               ',/,&
'@    Erreur dans la construction de la connectivite          ',/,&
'@    cellules -> faces (LAGDEB)                              ',/,&
'@                                                            ',/,&
'@    Mauvais remplissage de ICOCEL                           ',/,&
'@    (ICOCEL contient au moins une valeur nulle)             ',/,&
'@                                                            ',/,&
'@    Verifier le remplissage de ICOCEL dans LAGDEB           ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 5000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========                                               ',/,&
'@    Erreur dans la construction de la connectivite          ',/,&
'@    cellules -> faces (LAGDEB)                              ',/,&
'@                                                            ',/,&
'@    La taille du tableau ICOCEL n''est plus identique a     ',/,&
'@    celle calculee dans LAGINI                              ',/,&
'@                                                            ',/,&
'@     - LAGINI : LNDNOD = ',I10                               ,/,&
'@     - LAGDEB : LNDNOD = ',I10                               ,/,&
'@                                                            ',/,&
'@    Verifier le calcul de LNDNOD dans LAGINI et LAGDEB      ',/,&
'@    ou l''integrite des tableaux IFACEL et IFABOR           ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

!----
! FIN
!----

return
end subroutine
