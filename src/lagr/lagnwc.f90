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

subroutine lagnwc &
!================

 ( idbia0 , idbra0 ,                                              &
   lndnod ,                                                       &
   nbpmax , nvp    , nvp1   , nvep   , nivep  ,                   &
   npt    , nptnew , new    ,                                     &
   itycel , icocel ,                                              &
   ifrlag , isorti , iworkp ,                                     &
   ia     ,                                                       &
   surfbn , ettp   ,                                              &
   ra     )


!===============================================================================
! FONCTION :
! ----------

!   SOUS-PROGRAMME DU MODULE LAGRANGIEN :
!   -------------------------------------

!       INJECTION EN CONTINUE DES PARTICULES

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
!                  !    !     !                                                !
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
! lndnod           ! e  ! <-- ! dim. connectivite cellules->faces              !
! nbpmax           ! e  ! <-- ! nombre max de particulies autorise             !
! nvp              ! e  ! <-- ! nombre de variables particulaires              !
! nvp1             ! e  ! <-- ! nvp sans position, vfluide, vpart              !
! nvep             ! e  ! <-- ! nombre info particulaires (reels)              !
! nivep            ! e  ! <-- ! nombre info particulaires (entiers)            !
! npt              ! e  ! --> ! nombre courant de particules                   !
! nptnew           ! e  ! <-- ! nombre total de nouvelles particules           !
!                  !    !     ! pour toutes les zones d'injection              !
! new              ! e  ! <-- ! nombre de nouvelles part a injecter            !
!                  !    !     ! pour la zone d'injection courante              !
! icocel           ! te ! <-- ! connectivite cellules -> faces                 !
!   (lndnod)       !    !     !    face de bord si numero negatif              !
! itycel           ! te ! <-- ! connectivite cellules -> faces                 !
!   (ncelet+1)     !    !     !    pointeur du tableau icocel                  !
! isorti           ! te ! <-- ! pour chaque particule :                        !
!   (nbpmax)       !    !     !    * numero de sa cellule                      !
!                  !    !     !    * 0 si sortie du domaine                    !
! iworkp(npbmax    ! te ! <-- ! numero de la face d'injection                  !
! ia(*)            ! ia ! --- ! main integer work array                        !
! surfbn(nfabor    ! tr ! <-- ! surface des faces de bord                      !
! ettp             ! tr ! <-- ! tableaux des variables liees                   !
!  (nbpmax,nvp)    !    !     !   aux particules etape courante                !
! ra(*)            ! ra ! --- ! main real work array                           !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

!==============================================================================
! Module files
!==============================================================================

use paramx
use numvar
use optcal
use entsor
use period
use lagpar
use lagran
use mesh

!==============================================================================

implicit none

! Arguments

integer          idbia0 , idbra0
integer          lndnod
integer          nbpmax , nvp    , nvp1   , nvep  , nivep
integer          npt    , nptnew , new

integer          icocel(lndnod) , itycel(ncelet+1)
integer          isorti(nbpmax)
integer          ifrlag(nfabor) , iworkp(nbpmax)
integer          ia(*)

double precision surfbn(nfabor)
double precision ettp(nbpmax,nvp)
double precision ra(*)

! Local variables

integer          idebia, idebra, ifinia

integer          np , iel , n1  , ifac , kfac , nbp
integer          ii , jj  , in  , isort
integer          indian , ifaold , ifanew
integer          itypfo , iconfo(100)
integer          idehor , ierrie , icecpt
integer          icelcr, ipercr ,itepas, iper

double precision rd(1)
double precision xf, yf, zf
double precision up, vp, wp, uf, vf, wf
double precision pta(3), ptb(3), vect(3), vectn(3)

!===============================================================================
! 0.  GESTION MEMOIRE
!===============================================================================

idebia = idbia0
idebra = idbra0

!===============================================================================
! 1. INITIALISATION
!===============================================================================

! Initialize variables to avoid compiler warnings

icelcr = 0
ipercr = 0

! Traitement de la periodicite

if (iperio.eq.1) then

  icelcr = idebia
  ipercr = icelcr + ncelet-ncel
  ifinia = ipercr + ncelet-ncel
  call iasize('lagnwc', ifinia)
  !==========

  do iel = 1,ncelet-ncel
    ia(icelcr+iel-1) = 0
    ia(ipercr+iel-1) = -1
  enddo
  call perloc(ia(icelcr), ia(ipercr))
  !==========

endif

!===============================================================================
! 2. Injection des particules
!===============================================================================

!     BOUCLE SUR LES NOUVELLES PARTICULES :

do np = 1,new

!   RAZ du compteur de cellules traversees

  icecpt = 0

!-->IDEHOR = 1 si exterieur du domaine de calcul

  idehor = 0

!      Incrementation du pointeur sur les particules

  npt = npt + 1

!      On sauvegarde la maille de depart

  isort = isorti(npt)

!      Tirage aleatoire

  n1 = 1
  call zufall(n1,rd)

!      Nouvelle position : on le fait dans la
!      direction de la vitesse mais possibilite
!      de le faire par rapport a la normale

  xf = ettp(npt,jxp) + rd(1) *dtp *ettp(npt,jup)
  yf = ettp(npt,jyp) + rd(1) *dtp *ettp(npt,jvp)
  zf = ettp(npt,jzp) + rd(1) *dtp *ettp(npt,jwp)

  up = ettp(npt,jup)
  vp = ettp(npt,jvp)
  wp = ettp(npt,jwp)

  uf = ettp(npt,juf)
  vf = ettp(npt,jvf)
  wf = ettp(npt,jwf)

  if (  xf.eq. ettp(npt,jxp)                                      &
  .and. yf.eq. ettp(npt,jyp)                                      &
  .and. zf.eq. ettp(npt,jzp) ) goto 300

! Est-ce que le point est dans le domaine (Clone de LAGCEL)

!      Numero de face d 'injection dans l'element

  ifanew = iworkp(npt)

 100    continue

    iel    = isort
    ifaold = ifanew
    indian = 0
    icecpt = icecpt + 1

!         ---> Elimination des particules qui posent problemes
!              la particule reste au niveau de la face d'entree
!              (boucles infinies)

    if (icecpt.gt.30) then
      idehor = 1
      goto 200
    endif

!      --> balayage des KFAC faces entourant la cellule IEL
!           Elles sont stockees entre ITYCEL(IEL) et ITYCEL(IEL+1)-1
!          (donc KFAC ne peut pas valoir ITYCEL(IEL+1)...)

    kfac = itycel(iel)-1

    do while (indian.eq.0)

      kfac = kfac + 1

!           --> Erreur : la particule reste au niveau de la face
!                        d'entree

      if (kfac.eq.itycel(iel+1)) then
        idehor = 1
        goto 200
      endif

      ifac = icocel(kfac)

!-->boucle sur les faces internes

!-->si la face interne a deja ete traitee dans la cellule precedente
!   son numero est dans IFAOLD et on ne la retraite pas une seconde fois

      if (ifac.gt.0 .and. ifac.ne.ifaold) then

        in = 0
        do nbp = ipnfac(ifac),ipnfac(ifac+1)-1
          in = in + 1
          iconfo(in) = nodfac(nbp)
        enddo
        itypfo = ipnfac(ifac+1) - ipnfac(ifac) + 1
        iconfo(itypfo) = iconfo(1)

        call ouestu                                               &
        !==========
   (    nfecra , ndim   , nnod ,                                  &
        ierrie ,                                                  &
        ettp(npt,jxp)  , ettp(npt,jyp)  , ettp(npt,jzp)  ,        &
        xf                , yf                , zf               ,&
        cdgfac(1,ifac)    , cdgfac(2,ifac)    , cdgfac(3,ifac)   ,&
        xyzcen(1,iel)     , xyzcen(2,iel)     , xyzcen(3,iel)    ,&
        itypfo , iconfo , xyznod ,                                &
        indian )

        if (ierrie.eq.1) then
          idehor = 1
          goto 200

!-->si la particule passe dans la cellule voisine

        else if (indian.eq.1) then

!-->si la particule NPT est dans la cellule II alors le voisin ne
!   peut etre que la cellule JJ, et vice versa, et inversement, et
!   ainsi de suite.

          ifanew = ifac

          ii = ifacel(1,ifac)
          jj = ifacel(2,ifac)

          if (iel.eq.ii) then
            isort = jj
          else if (iel.eq.jj) then
            isort = ii
          endif

! Traitement de la periodicite
! Meme commentaire que dans lagcel.F

          if (isort.gt.ncel) then

            itepas = isort
            isort = ia(icelcr+itepas-ncel-1)

!                On recupere les informations sur la periodicite

            iper  = ia(ipercr+itepas-ncel-1)

!                 MODIFICATION DE LA POSITION

            pta(1)= xf
            pta(2)= yf
            pta(3)= zf

            call lagper(iper,pta, ptb)
            !==========

            xf = ptb(1)
            yf = ptb(2)
            zf = ptb(3)

!                 MODIFICATION DE LA VITESSE

            vect(1) = up
            vect(2) = vp
            vect(3) = wp

            call lagvec(iper, vect, vectn)
            !==========

            up = vectn(1)
            vp = vectn(2)
            wp = vectn(3)

!                 MODIFICATION DE LA VITESSE FLUIDE VUE

            vect(1) = uf
            vect(2) = vf
            vect(3) = wf

            call lagvec(iper, vect, vectn)
            !==========

            uf = vectn(1)
            vf = vectn(2)
            wf = vectn(3)

            ifanew = 0

          endif

!--> Retour pour balayage des face de la cellule suivante

          goto 100

        endif

!--> Balayage des faces de bord (reperees par leur valeur negative
!    dans ICOCEL)

!    resultat : INDIAN =  0 le rayon PQ ne sort pas de la cellule par
!    ~~~~~~~~               cette face
!               INDIAN = -1 meme cellule
!               INDIAN =  1 interaction avec la frontiere

      else if (ifac.lt.0 .and. ifac.ne.ifaold) then

        ifac = -ifac

        in = 0
        do nbp = ipnfbr(ifac),ipnfbr(ifac+1)-1
          in = in + 1
          iconfo(in) = nodfbr(nbp)
        enddo
        itypfo = ipnfbr(ifac+1) - ipnfbr(ifac) + 1
        iconfo(itypfo) = iconfo(1)

        call ouestu                                               &
        !==========
   (    nfecra , ndim   , nnod ,                                  &
        ierrie ,                                                  &
        ettp(npt,jxp)  , ettp(npt,jyp)  , ettp(npt,jzp)  ,        &
        xf                , yf                , zf               ,&
        cdgfbo(1,ifac)    , cdgfbo(2,ifac)    , cdgfbo(3,ifac)   ,&
        xyzcen(1,iel)     , xyzcen(2,iel)     , xyzcen(3,iel)    ,&
        itypfo , iconfo , xyznod ,                                &
        indian )

!-->si la trajectoire de la particule traverse la face de bord
!    alors on laisse la particule au niveau de la face d'entree

        if (ierrie.eq.1 .or. indian.eq.1) then
          idehor = 1
          goto 200
        endif

      endif

! fin de DO WHILE
    enddo

!-->Fin de la boucle principale sur les particules

 200      continue

! Si le point est dans le domaine, alors on injecte la
! particule du point sinon on laisse la particule au niveau
! de l'entree

  if (idehor.eq.0) then
    ettp(npt,jxp) = xf
    ettp(npt,jyp) = yf
    ettp(npt,jzp) = zf
    isorti(npt)   = isort
  endif

 300      continue

enddo

!==============================================================================

!----
! FIN
!----

return

end subroutine
