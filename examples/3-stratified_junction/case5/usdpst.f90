!-------------------------------------------------------------------------------

!                      Code_Saturne version 2.0.0-beta1
!                      --------------------------

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

subroutine  usdpst &
!=================

 ( idbia0 , idbra0 ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   lstcel , lstfac , lstfbr ,                                     &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   rdevel , rtuser , ra     )

!===============================================================================
! FONCTION :
! --------

! ROUTINE UTILISATEUR POUR LOCALISER DES CELLULES, DES FACES
! INTERNES ET/OU DES FACES DE BORD DEFINISSANT UN MAILLAGE DE
! POST-TRAITEMENT.

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
!    nom           !type!mode !                   role                         !
!__________________!____!_____!________________________________________________!
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
! nideve nrdeve    ! e  ! <-- ! longueur de idevel rdevel                      !
! nituse nrtuse    ! e  ! <-- ! longueur de ituser rtuser                      !
! ifacel           ! te ! <-- ! elements voisins d'une face interne            !
! (2, nfac)        !    !     !                                                !
! ifabor           ! te ! <-- ! element  voisin  d'une face de bord            !
! (nfabor)         !    !     !                                                !
! ifmfbr           ! te ! <-- ! numero de famille d'une face de bord           !
! (nfabor)         !    !     !                                                !
! ifmcel           ! te ! <-- ! numero de famille d'une cellule                !
! (ncelet)         !    !     !                                                !
! iprfml           ! te ! <-- ! proprietes d'une famille                       !
! (nfml,nprfml)    !    !     !                                                !
! ipnfac           ! te ! <-- ! position du premier noeud de chaque            !
!   (nfac+1)       !    !     !  face interne dans nodfac (optionnel)          !
! nodfac           ! te ! <-- ! connectivite faces internes/noeuds             !
!   (lndfac)       !    !     !  (optionnel)                                   !
! ipnfbr           ! te ! <-- ! position du premier noeud de chaque            !
!  (nfabor+1)      !    !     !  face de bord dans nodfbr (optionnel)          !
! nodfbr           ! te ! <-- ! connectivite faces de bord/noeuds              !
!   (lndfbr  )     !    !     !  (optionnel)                                   !
! lstcel           ! te ! --- ! tableau de travail (liste des                  !
! (ncelet)         !    !     !  cellules d'un maillage de sortie)             !
! lstfac           ! te ! --- ! tableau de travail (liste des faces            !
! (nfac)           !    !     !  internes d'un maillage de sortie)             !
! lstfbr           ! te ! --- ! tableau de travail (liste des faces            !
! (nfabor)         !    !     !  de bord d'un maillage de sortie)              !
! ia(*)            ! te ! --- ! macro tableau entier                           !
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
! ra(*)            ! tr ! --- ! macro tableau reel                             !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

implicit none

!===============================================================================

!===============================================================================
!     DONNEES EN COMMON
!===============================================================================

include "paramx.h"
include "optcal.h"
include "entsor.h"
include "parall.h"
include "period.h"

!===============================================================================

! Arguments

integer          idbia0 , idbra0
integer          ndim   , ncelet , ncel   , nfac   , nfabor
integer          nfml   , nprfml
integer          nnod   , lndfac , lndfbr
integer          nideve , nrdeve , nituse , nrtuse

integer          ifacel(2,nfac) , ifabor(nfabor)
integer          ifmfbr(nfabor) , ifmcel(ncelet)
integer          iprfml(nfml,nprfml)
integer          ipnfac(nfac+1), nodfac(lndfac)
integer          ipnfbr(nfabor+1), nodfbr(lndfbr)
integer          lstcel(ncelet), lstfac(nfac), lstfbr(nfabor)
integer          idevel(nideve), ituser(nituse)
integer          ia(*)

double precision xyzcen(ndim,ncelet)
double precision surfac(ndim,nfac), surfbo(ndim,nfabor)
double precision cdgfac(ndim,nfac), cdgfbo(ndim,nfabor)
double precision xyznod(ndim,nnod), volume(ncelet)
double precision rdevel(nrdeve), rtuser(nrtuse)
double precision ra(*)

! VARIABLES LOCALES

integer          indmod, icas, nbcas, ipart, nbpart, ipref
integer          ntchrl

integer          nlcel, nlfac , nlfbr
integer          iel, ifac  , ii
integer          idebia, idebra
integer          icoul , icoul1, icoul2, iel1  , iel2
character*32     nomcas, nomfmt, nommai
character*96     nomrep, optfmt

double precision xfac  , yfac  , zfac


!===============================================================================



nbcas  = 0
nbpart = 0

! Entiers "pointeurs" sur la premiere case libre de IA et RA

idebia = idbia0
idebra = idbra0

!===============================================================================
!     CREATION DES GESTIONNAIRES D'ECRITURE POUR LE POST TRAITEMENT
!         (UN PAR CAS ET PAR FORMAT, A RENSEIGNER PAR L'UTILISATEUR)
!===============================================================================

!     NOMBRE DE GESTIONNAIRES (case au sens EnSight, etude au sens MED,
!                              ou racine d'une arborescence CGNS)

nbcas = 1

do icas = 1, nbcas

!       INITIALISATIONS DIVERSES

  do ii = 1, len(nomcas)
    NOMCAS (II:II) = ' '
  enddo
  do ii = 1, len(nomrep)
    NOMREP (II:II) = ' '
  enddo
  do ii = 1, len(nomfmt)
    NOMFMT (II:II) = ' '
  enddo
  do ii = 1, len(optfmt)
    OPTFMT (II:II) = ' '
  enddo

!       DEFINITION UTILISATEUR :

!       NOMCAS et NOMREP indiquent respectivement le prefixe du nom
!       des fichiers et le repertoire correspondant.
!       Si NOMREP est de la forme xxxx.ensight ou xxxx.med, le lanceur le
!       rapatriera automatiquement sous le nom XXXX.ENSIGHT.$DATE ou
!       XXXX.MED.$DATE dans le repertoire RESU. Si NOMREP est d'une autre
!       forme, il faudra gerer son rapatriement a la main.

!       NOMFMT permet de choisir le format de sortie
!       ("EnSight Gold", "MED_fichier", ou "CGNS").

!       OPTFMT permet de fournir des options specifiques au format de
!       sortie (separees par des virgules) ;
!         Pour EnSight : "text" ou "binary" (defaut),
!         Pour EnSight, MED, ou CGNS :
!                        "discard_polygons" pour supprimer les polygones,
!                        "discard_polyhedra" pour supprimer les polyedres.
!         Pour EnSight  ou MED :
!                        "divide_polygons" pour découper les polygones,
!                        "divide_polyhedra" pour découper les polyedres.

!       INDMOD indique si les maillages ecrits seront :
!         0 : fixes,
!         1 : deformables a topologie constante,
!         2 : modifiables (pourront etre completement redefinis en
!             cours de calcul via le sous-programme USMPST).
!        10 : comme INDMOD = 0, avec champ de déplacement
!        11 : comme INDMOD = 1, avec champ de déplacement
!        12 : comme INDMOD = 2, avec champ de déplacement

!       NTCHRL donne la frequence de sortie par defaut associee,
!       (la sortie a un pas de temps donne pouvant etre forcee ou
!       empechee via le sous-programme utilisateur USNPST).

  if (icas .eq. 1) then

    NOMCAS = 'chr'
    NOMREP = 'TINF21.ensight'
    NOMFMT = 'EnSight Gold'
    OPTFMT = 'binary, discard_polygons'
    indmod = 2
    ntchrl = 5

  endif

!       DEFINITION EFFECTIVE

  call pstcwr (icas  , nomcas, nomrep, nomfmt, optfmt, &
  !==========
               indmod, ntchrl)

enddo

!===============================================================================
!     NOMBRE DE MAILLAGES EXTRAITS POUR POST TRAITEMENT
!         A RENSEIGNER PAR L'UTILISATEUR
!===============================================================================

!   NBPART est le nombre de "parts" qui seront generees
!   (au sens EnSight ; les équivalents MED et CGNS sont le maillage
!    et la base respectivement)

!   Une "part" peut etre tout volume ou surface que l'on definira par
!   l'identification des cellules ou faces du maillage


nbpart = 2

!===============================================================================
!     DEBUT DE LA BOUCLE SUR LES PARTS DEFINIES PAR L'UTILISATEUR
!===============================================================================

do ipart = 1, nbpart


!===============================================================================
!       INITIALISATIONS DIVERSES
!         PAS D'INTERVENTION UTILISATEUR REQUISE
!===============================================================================

  nlcel = 0
  nlfac = 0
  nlfbr = 0
  do iel = 1, ncelet
    lstcel(iel) = 0
  enddo
  do ifac = 1, nfac
    lstfac(ifac) = 0
  enddo
  do ifac = 1, nfabor
    lstfbr(ifac) = 0
  enddo

  do ii = 1, len(nommai)
    NOMMAI(II:II) = ' '
  enddo

!===============================================================================
!       REPERAGE DES CELLULES OU FACES INCLUSES DANS LE MAILLAGE
!         A RENSEIGNER PAR L'UTILISATEUR
!===============================================================================

!       Ce sous programme est appele avant la definition des
!        conditions aux limites


!       POUR LA 1ere COUPE (PART 1) : coupe exemple

!         Exemple : on selectionne
!                   les faces internes sur le plan median

  if (ipart .eq. 1) then

    NOMMAI = 'Coupe'

!         Pour les faces internes

    do ifac = 1, nfac

!           Determination si la face appartient a la coupe

      if (abs(cdgfac(2,ifac)).lt.1.d-4) then
        nlfac = nlfac+1
        lstfac(nlfac)= ifac
      endif

    enddo




!       POUR LA 2eme COUPE (PART 2) : cellules a T<21 degres

!         Exemple : on initialise a NCELET, le choix
!                   des cellules sera fait dans usmpst.F

  else if (ipart .eq. 2) then

    NOMMAI = 'celTinf21'

!         Si NLCEL = NCELET, pas la peine de definir LSTCEL

    nlcel = ncelet

  endif

!===============================================================================
!       CREATION DES STRUCTURES CONSERVANT LES DONNEES DES PARTS
!         PAS D'INTERVENTION UTILISATEUR REQUISE
!===============================================================================

  call pstcma (ipart, nommai, nlcel, nlfac, nlfbr, lstcel, lstfac, lstfbr)
  !==========

!===============================================================================
!       IDENTIFICATION DU MAILLAGE EXTRAIT ET GESTION DE SORTIE
!         A RENSEIGNER PAR L'UTILISATEUR
!===============================================================================

  if (ipart .eq. 1) then

!         Le maillage 1 est associe a la sortie standard (ICAS=-1)
    icas = -1
    call pstass(ipart, icas)
    !==========

  else if (ipart .eq. 2) then

!         Le maillage 2 est associe au cas cree ici (ICAS=1)
    icas = 1
    call pstass(ipart, icas)
    !==========

  endif

!===============================================================================
!     FIN   DE LA BOUCLE SUR LES PARTS DEFINIES PAR L'UTILISATEUR
!===============================================================================

enddo

return

!===============================================================================
!     FORMATS
!===============================================================================

end

