!-------------------------------------------------------------------------------

!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2008 EDF S.A., France

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

subroutine invers &
!================

 ( cnom   , idbia0 , idbra0 ,                                     &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nideve , nrdeve , nituse , nrtuse ,                            &
   isym   , ipol   , ireslp , nitmap , imgrp  ,                   &
   ncymxp , nitmfp ,                                              &
   iwarnp , nfecra , niterf , icycle , iinvpe ,                   &
   epsilp , rnorm  , residu ,                                     &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dam    , xam    , smbrp  , vx     ,                            &
   w1     , w2     , w3     , w4     , w5     , w6     ,          &
   rdevel , rtuser , ra     )

!===============================================================================
! FONCTION :
! ----------

! APPEL AUX ROUTINE D'INVERSION DE SYSTEMES LINEAIRE
!  MULTIGRILLE + GRADCO OU JACOBI OU BI-CGSTAB
!  GRADCO
!  JACOBI
!  BI-CGSTAB
!  ON SUPPOSE VX INITIALISE EN ENTREE

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
!    nom           !type!mode !                   role                         !
!__________________!____!_____!________________________________________________!
! cnom             ! a  ! <-- ! nom de la variable                             !
! idbia0           ! e  ! <-- ! numero de la 1ere case libre dans ia           !
! idbra0           ! e  ! <-- ! numero de la 1ere case libre dans ra           !
! ndim             ! e  ! <-- ! dimension de l'espace                          !
! ncelet           ! e  ! <-- ! nombre d'elements halo compris                 !
! ncel             ! e  ! <-- ! nombre d'elements actifs                       !
! nfac             ! e  ! <-- ! nombre de faces internes                       !
! nfabor           ! e  ! <-- ! nombre de faces de bord                        !
! nfml             ! e  ! <-- ! nombre de familles d entites                   !
! nprfml           ! e  ! <-- ! nombre de proprietese des familles             !
! nnod             ! e  ! <-- ! nombre de sommets                              !
! lndfac           ! e  ! <-- ! longueur du tableau nodfac (optionnel          !
! lndfbr           ! e  ! <-- ! longueur du tableau nodfbr (optionnel          !
! ncelbr           ! e  ! <-- ! nombre d'elements ayant au moins une           !
!                  !    !     ! face de bord                                   !
! nideve nrdeve    ! e  ! <-- ! longueur de idevel rdevel                      !
! nituse nrtuse    ! e  ! <-- ! longueur de ituser rtuser                      !
! isym             ! e  ! <-- ! indicateur = 1 matrice sym                     !
!                  !    !     !            = 2 matrice non sym                 !
! ipol             ! e  ! <-- ! degre du polynome pour precond                 !
!                  !    !     !         (0 -> diagonal)                        !
! ireslp           ! e  ! <-- ! indicateur = 0 gradco                          !
!                  !    !     !            = 1 jacobi                          !
!                  !    !     !            = 2 cgstab                          !
! nitmap           ! e  ! <-- ! nombre max d'iter pour resol iterativ          !
! imgrp            ! e  ! <-- ! indicateur = 0 pas de multigrille              !
! ncymxp           ! e  ! <-- ! nombre de cycles max pour multigrille          !
! nitmfp           ! e  ! <-- ! nombre d iter sur maillage fin                 !
! iwarnp           ! e  ! <-- ! niveau d'impression                            !
! nfecra           ! e  ! <-- ! unite du fichier sortie std                    !
! niterf           ! e  ! --> ! nombre d'iterations effectuees                 !
!                  !    !     !  (non multigrille)                             !
! icycle           ! e  ! --> ! nombre de cycles mgm effectues                 !
! iinvpe           ! e  ! <-- ! indicateur pour annuler les increment          !
!                  !    !     ! en periodicite de rotation (=2) ou             !
!                  !    !     ! pour les echanger normalement de               !
!                  !    !     ! maniere scalaire (=1)                          !
! epsilp           ! r  ! <-- ! precision pour resol iter                      !
! rnorm            ! r  ! <-- ! normalisation du residu                        !
! residu           ! r  ! --> ! residu final non norme                         !
! ifacel           ! te ! <-- ! elements voisins d'une face interne            !
! (2, nfac)        !    !     !                                                !
! ifabor           ! te ! <-- ! element  voisin  d'une face de bord            !
! (nfabor)         !    !     !                                                !
! ifmfbr           ! te ! <-- ! numero de famille d'une face de bord           !
! (nfabor)         !    !     !                                                !
! ifmcel           ! te ! <-- ! numero de famille d'une cellule                !
! (ncelet)         !    !     !                                                !
! iprfml           ! te ! <-- ! proprietes d'une famille                       !
! nfml  ,nprfml    !    !     !                                                !
! ipnfac           ! te ! <-- ! position du premier noeud de chaque            !
!   (lndfac)       !    !     !  face interne dans nodfac (optionnel)          !
! nodfac           ! te ! <-- ! connectivite faces internes/noeuds             !
!   (nfac+1)       !    !     !  (optionnel)                                   !
! ipnfbr           ! te ! <-- ! position du premier noeud de chaque            !
!   (lndfbr)       !    !     !  face de bord dans nodfbr (optionnel)          !
! nodfbr           ! te ! <-- ! connectivite faces de bord/noeuds              !
!   (nfabor+1)     !    !     !  (optionnel)                                   !
! idevel(nideve    ! te ! <-- ! tab entier complementaire developemt           !
! ituser(nituse    ! te ! <-- ! tab entier complementaire utilisateur          !
! ia(*)            ! tr ! --- ! macro tableau entier                           !
! xyzcen           ! tr ! <-- ! point associes aux volumes de control          !
! (ndim,ncelet     !    !     !                                                !
! surfac           ! tr ! <-- ! vecteur surface des faces internes             !
! (ndim,nfac)      !    !     !                                                !
! surfbo           ! tr ! <-- ! vecteur surface des faces de bord              !
! (ndim,nfabor)    !    !     !                                                !
! cdgfac           ! tr ! <-- ! centre de gravite des faces internes           !
! (ndim,nfac)      !    !     !                                                !
! cdgfbo           ! tr ! <-- ! centre de gravite des faces de bord            !
! (ndim,nfabor)    !    !     !                                                !
! xyznod           ! tr ! <-- ! coordonnes des noeuds (optionnel)              !
! (ndim,nnod)      !    !     !                                                !
! volume           ! tr ! <-- ! volume d'un des ncelet elements                !
! (ncelet          !    !     !                                                !
! dam(ncelet       ! tr ! <-- ! diagonale (maillage fin si mgm)                !
! xam(nfac,isym    ! tr ! <-- ! extradiagonale (maillage fin si mgm)           !
! smbrp(ncelet     ! tr ! <-- ! second membre (maillage fin si mgm)            !
! vx   (ncelet     ! tr ! <-- ! solution du systeme                            !
! w1,2,3,4,5,6     ! tr ! --- ! auxiliaires de travail                         !
!      (ncelet     !    !     !                                                !
! rdevel(nrdeve    ! tr ! <-- ! tab reel complementaire developemt             !
! rtuser(nrtuse    ! tr ! <-- ! tab reel complementaire utilisateur            !
! ra(*)            ! tr ! --- ! macro tableau reel                             !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

implicit none

!===============================================================================
!     DONNEES EN COMMON
!===============================================================================

include "paramx.h"
include "optcal.h"

!===============================================================================

! Arguments

character*8      cnom
integer          idbia0 , idbra0
integer          ndim   , ncelet , ncel   , nfac   , nfabor
integer          nfml   , nprfml
integer          nnod   , lndfac , lndfbr , ncelbr
integer          nideve , nrdeve , nituse , nrtuse
integer          isym   , ipol   , ireslp , nitmap
integer          imgrp  , ncymxp , nitmfp
integer          iwarnp , nfecra
integer          niterf , icycle , iinvpe
double precision epsilp , rnorm  , residu

integer          ifacel(2,nfac) , ifabor(nfabor)
integer          ifmfbr(nfabor) , ifmcel(ncelet)
integer          iprfml(nfml,nprfml)
integer          ipnfac(nfac+1), nodfac(lndfac)
integer          ipnfbr(nfabor+1), nodfbr(lndfbr)
integer          idevel(nideve), ituser(nituse)
integer          ia(*)

double precision xyzcen(ndim,ncelet)
double precision surfac(ndim,nfac), surfbo(ndim,nfabor)
double precision cdgfac(ndim,nfac), cdgfbo(ndim,nfabor)
double precision xyznod(ndim,nnod), volume(ncelet)
double precision dam(ncelet), xam(nfac ,2)
double precision smbrp(ncelet)
double precision vx(ncelet)
double precision w1(ncelet),w2(ncelet),w3(ncelet),w4(ncelet)
double precision w5(ncelet),w6(ncelet)
double precision rdevel(nrdeve), rtuser(nrtuse), ra(*)


! VARIABLES LOCALES

integer          idebia, idebra
integer          lnom
integer          iresds, iresas, nitmds, nitmas

!===============================================================================

! INITIALISATIONS

lnom = len(cnom)

icycle = 0
niterf = 0

idebia = idbia0
idebra = idbra0

! RESOLUTION

if( imgrp.eq.1 ) then

  iresds = ireslp
  iresas = ireslp

  nitmds = nitmfp
  nitmas = nitmfp

  call resmgr                                                     &
  !==========
 ( cnom   , lnom   , ncelet , ncel   , nfac   ,                   &
   isym   , iresds , iresas , ireslp , ipol   ,                   &
   ncymxp , nitmds , nitmas , nitmap , iinvpe ,                   &
   iwarnp , icycle , niterf , epsilp , rnorm  , residu ,          &
   ifacel , smbrp  , vx     )

elseif(imgrp.eq.0) then

  if (ireslp.ge.0 .and. ireslp.le. 2) then

    call reslin                                                   &
    !==========
 ( cnom   , lnom   , ncelet , ncel   , nfac   ,                   &
   isym   , ireslp , ipol   , nitmap , iinvpe ,                   &
   iwarnp , niterf , epsilp , rnorm  , residu ,                   &
!                 ------                     ------
   ifacel , dam    , xam    , smbrp  , vx     )
!                                   -----

  else
    write(nfecra,1000) cnom, ireslp
    call csexit (1)
  endif

endif


#if defined(_CS_LANG_FR)

 1000 format('INVERS APPELE POUR ',A8,' AVEC IRESOL = ', I10)

#else

 1000 format('INVERS CALLED FOR ',A8,' WITH IRESOL = ', I10)

#endif

!----
! FIN
!----

return

end
