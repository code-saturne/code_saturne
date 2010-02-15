!-------------------------------------------------------------------------------

!VERS


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

subroutine usdpst &
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
! lndfac           ! i  ! <-- ! size of nodfac indexed array                   !
! lndfbr           ! i  ! <-- ! size of nodfbr indexed array                   !
! nideve, nrdeve   ! i  ! <-- ! sizes of idevel and rdevel arrays              !
! nituse, nrtuse   ! i  ! <-- ! sizes of ituser and rtuser arrays              !
! ifacel(2, nfac)  ! ia ! <-- ! interior faces -> cells connectivity           !
! ifabor(nfabor)   ! ia ! <-- ! boundary faces -> cells connectivity           !
! ifmfbr(nfabor)   ! ia ! <-- ! boundary face family numbers                   !
! ifmcel(ncelet)   ! ia ! <-- ! cell family numbers                            !
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
! xyzcen           ! ra ! <-- ! cell centers                                   !
!  (ndim, ncelet)  !    !     !                                                !
! surfac           ! ra ! <-- ! interior faces surface vectors                 !
!  (ndim, nfac)    !    !     !                                                !
! surfbo           ! ra ! <-- ! boundary faces surface vectors                 !
!  (ndim, nfabor)  !    !     !                                                !
! cdgfac           ! ra ! <-- ! interior faces centers of gravity              !
!  (ndim, nfac)    !    !     !                                                !
! cdgfbo           ! ra ! <-- ! boundary faces centers of gravity              !
!  (ndim, nfabor)  !    !     !                                                !
! xyznod           ! ra ! <-- ! vertex coordinates (optional)                  !
!  (ndim, nnod)    !    !     !                                                !
! volume(ncelet)   ! ra ! <-- ! cell volumes                                   !
! ra(*)            ! ra ! --- ! main real work array                           !
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================

implicit none

!===============================================================================

!===============================================================================
! Common blocks
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

! Local variables

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


! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================

if(1.eq.1) return

!===============================================================================
! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END

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

nbcas = 4

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
    NOMREP = 'EnSight'
    NOMFMT = 'EnSight Gold'
    OPTFMT = 'binary, discard_polygons'
    indmod = 0
    ntchrl = 4

  else if (icas .eq. 2) then

    NOMCAS = 'chr'
    NOMREP = 'EnSight_text'
    NOMFMT = 'ensight'
    OPTFMT = 'text, discard_polyhedra'
    indmod = 1
    ntchrl = ntchr

  else if (icas .eq. 3) then

    NOMCAS = 'modif'
    NOMREP = 'EnSight'
    NOMFMT = 'ensight'
    OPTFMT = 'discard_polyhedra'
    indmod = 2
    ntchrl = ntchr
    ntchrl = 2

  else if (icas .eq. 4) then

    NOMCAS = 'CHR'
    NOMREP = ' '
    NOMFMT = 'MED'
    OPTFMT = ' '
    indmod = 1
    ntchrl = ntchr

  endif

!       DEFINITION EFFECTIVE

  call pstcwr (icas  , nomcas, nomrep, nomfmt, optfmt, indmod, ntchrl)
  !==========

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

!   Exemple : 4 "parts", correspondant respectivement à une coupe
!             mixte faces de bord / faces internes, puis à une coupe
!             contenant uniquement des faces internes, puis à deux extraits
!             evolutifs du maillage global. On ajoutera une 5eme "part",
!             alias de la 2eme.

nbpart = 4

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
!                   les faces internes situees entre une cellule de
!                   couleur 2 et une cellule de couleur 3
!                   et les faces de bord  de couleur 4

  if (ipart .eq. 1) then

    NOMMAI = 'Coupe 1'

!         Pour les faces internes

    do ifac = 1, nfac

!           Elements voisins
      iel1 = ifacel(1,ifac)
      iel2 = ifacel(2,ifac)

!           Couleur des elements voisins
      icoul1 = iprfml(ifmcel(iel1),1)
      icoul2 = iprfml(ifmcel(iel2),1)

!           Determination si la face appartient a la coupe

      if ((icoul1.eq.2.and.icoul2.eq.3).or.(icoul1.eq.3.and.icoul2.eq.2)) then
        nlfac = nlfac+1
        lstfac(nlfac)= ifac
      endif

    enddo

!         Pour les faces de bord

    do ifac = 1, nfabor

!           Couleur de la face
      icoul = iprfml(ifmfbr(ifac),1)

!           Determination si la face appartient a la coupe
      if (icoul.eq.4) then
        nlfbr = nlfbr+1
        lstfbr(nlfbr)= ifac
      endif

    enddo



!       POUR LA 2eme COUPE (PART 2) : coupe exemple

!         Exemple : on selectionne
!                   les faces internes situees a y = 0.5

  else if (ipart .eq. 2) then

    NOMMAI = 'Coupe 2'

!         Pour les faces internes

    do ifac = 1, nfac

!           Determination si la face appartient a la coupe

!           Centre de gravite de la face
      xfac = cdgfac(1,ifac)
      yfac = cdgfac(2,ifac)
      zfac = cdgfac(3,ifac)

      if ((yfac .gt. 0.4999).and.(yfac .lt. 0.5001)) then
        nlfac = nlfac+1
        lstfac(nlfac)= ifac
      endif

    enddo

!       POUR LA PART NUMERO 3 (EXTRAIT DU DOMAINE FLUIDE)

!         Par defaut : on selectionne toutes les cellules, on modifiera
!                      la selection dans USMPST

  else if (ipart .eq. 3) then

    NOMMAI = 'Volume v > 0.5'

    nlcel = ncel
    nlfac = 0
    nlfbr = 0

!       POUR LA PART NUMERO 4 (EXTRAIT DU DOMAINE FLUIDE)

!         Par defaut : on selectionne toutes les faces de bord,
!                      on modifiera la selection dans USMPST

  else if (ipart .eq. 4) then

    NOMMAI = 'Surface "iso" v'

    nlcel = 0
    nlfac = 0
    nlfbr = nfabor

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

  if ((ipart .eq. 1) .or. (ipart .eq. 2)) then

!         Les maillages extraits numero 1 et 2 sont associes aux cas 1 et 2
    icas = 1
    call pstass(ipart, icas)
    !==========
    icas = 2
    call pstass(ipart, icas)
    !==========

  else if ((ipart .eq. 3) .or. (ipart .eq. 4)) then

!         Les maillages extraits numero 3 et 4 sont associes au cas 3
    icas = 3
    call pstass(ipart, icas)
    !==========

  endif

!===============================================================================
!     FIN   DE LA BOUCLE SUR LES PARTS DEFINIES PAR L'UTILISATEUR
!===============================================================================

enddo

!===============================================================================
!     TRAITEMENT DES ALIAS EVENTUELS ; UN ALIAS EST SURTOUT UTILE
!     LORSQU'ON LUI AFFECTE UN 'WRITER' DIFFERENT DU MAILLAGE AUQUEL
!     IL FAIT REFERENCE, AVEC LEQUEL ON PEUT SORTIR DES VARIABLES
!     DIFFERENTES, SANS DUPLICATION EN MEMOIRE DU MAILLAGE POST
!===============================================================================

!     PART NUMERO 7 : ALIAS DE LA PART 2 (POUR SORTIES DE VARIABLES
!                       ASSOCIEES INDEPENDANTE)

ipart = 5
ipref = 2

call pstalm(ipart, ipref)
!==========

!     Pour le maillage numero 5 : writer 3
icas = 3
call pstass(ipart, icas)
!==========


return

!===============================================================================
!     FORMATS
!===============================================================================

end subroutine

