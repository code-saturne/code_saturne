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

subroutine recvmc &
!================

 ( idbia0 , idbra0 ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   rom    , flumas , flumab ,                                     &
   ux     , uy     , uz     ,                                     &
   bx     , by     , bz     , cocg   ,                            &
   rdevel , rtuser , ra     )

!===============================================================================
! FONCTION :
! ----------

! RECONSTRUCTION DE LA VITESSE A PARTIR DU FLUX DE MASSE
!     PAR MOINDRES CARRES (VITESSE CONSTANTE PAR ELEMENT)

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
! ncelbr           ! e  ! <-- ! nombre d'elements ayant au moins une           !
!                  !    !     ! face de bord                                   !
! nvar             ! e  ! <-- ! nombre total de variables                      !
! nscal            ! e  ! <-- ! nombre total de scalaires                      !
! nphas            ! e  ! <-- ! nombre de phases                               !
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
! rom(ncelet       ! tr ! <-- ! masse volumique aux cellules                   !
! flumas(nfac)     ! tr ! <-- ! flux de masse aux faces internes               !
! flumab(nfabor    ! tr ! <-- ! flux de masse aux faces de bord                !
! ux   uy          ! tr ! --> ! vitesse reconstruite                           !
! uz   (ncelet     ! tr !     !                                                !
! bx,y,z(ncelet    ! tr ! --- ! tableau de travail                             !
! cocg             ! tr ! --- ! tableau de travail                             !
!   (ncelet,3,3    !    !     !                                                !
! rdevel(nrdeve    ! tr ! <-- ! tab reel complementaire developemt             !
! rtuser(nrtuse    ! tr ! <-- ! tab reel complementaire utilisateur            !
! ra(*)            ! tr ! --- ! macro tableau reel                             !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!-------------------------------------------------------------------------------
!===============================================================================

implicit none

!===============================================================================
!     DONNEES EN COMMON
!===============================================================================

include "vector.h"

!===============================================================================

! Arguments

integer          idbia0 , idbra0
integer          ndim   , ncelet , ncel   , nfac   , nfabor
integer          nfml   , nprfml
integer          nnod   , lndfac , lndfbr , ncelbr
integer          nvar   , nscal  , nphas
integer          nideve , nrdeve , nituse , nrtuse

integer          ifacel(2,nfac) , ifabor(nfabor)
integer          ifmfbr(nfabor) , ifmcel(ncelet)
integer          iprfml(nfml,nprfml)
integer          ipnfac(nfac+1), nodfac(lndfac)
integer          ipnfbr(nfabor+1), nodfbr(lndfbr)
integer          idevel(nideve), ituser(nituse), ia(*)

double precision xyzcen(ndim,ncelet)
double precision surfac(ndim,nfac), surfbo(ndim,nfabor)
double precision cdgfac(ndim,nfac), cdgfbo(ndim,nfabor)
double precision xyznod(ndim,nnod), volume(ncelet)
double precision rom(ncelet)
double precision flumas(nfac), flumab(nfabor)
double precision ux  (ncelet), uy  (ncelet), uz  (ncelet)
double precision bx(ncelet),   by(ncelet),   bz(ncelet)
double precision cocg(ncelet,3,3)
double precision rdevel(nrdeve), rtuser(nrtuse), ra(*)

! VARIABLES LOCALES

integer          lbloc
parameter       (lbloc = 1024)

integer          idebia, idebra, ii, jj, iel, ifac
integer          ibloc, nbloc, irel, idim1, idim2
double precision aa(lbloc,3,3)
double precision a11, a22, a33, a12, a13, a23, unsdet
double precision cocg11, cocg12, cocg13, cocg21, cocg22, cocg23
double precision cocg31, cocg32, cocg33
double precision smbx, smby, smbz, unsrho
double precision vecfac, pfacx, pfacy, pfacz

!===============================================================================

idebia = idbia0
idebra = idbra0

!===============================================================================
! 1. CALCUL DE LA MATRICE
!===============================================================================

!   INITIALISATION

do ii = 1, 3
  do jj = 1, 3
    do iel = 1, ncelet
      cocg(iel,ii,jj) = 0.d0
    enddo
  enddo
enddo

!   ASSEMBLAGE A PARTIR DES FACETTES FLUIDES

do idim1 = 1, 3
  do idim2 = idim1, 3

    if (ivecti.eq.1) then

!CDIR NODEP
      do ifac = 1, nfac
        ii = ifacel(1,ifac)
        jj = ifacel(2,ifac)
        vecfac = surfac(idim1,ifac)*surfac(idim2,ifac)
        cocg(ii,idim1,idim2) = cocg(ii,idim1,idim2) + vecfac
        cocg(jj,idim1,idim2) = cocg(jj,idim1,idim2) + vecfac
      enddo

    else

! VECTORISATION NON FORCEE
      do ifac = 1, nfac
        ii = ifacel(1,ifac)
        jj = ifacel(2,ifac)
        vecfac = surfac(idim1,ifac)*surfac(idim2,ifac)
        cocg(ii,idim1,idim2) = cocg(ii,idim1,idim2) + vecfac
        cocg(jj,idim1,idim2) = cocg(jj,idim1,idim2) + vecfac
      enddo

    endif

    if (ivectb.eq.1) then

!CDIR NODEP
      do ifac = 1, nfabor
        ii = ifabor(ifac)
        cocg(ii,idim1,idim2) = cocg(ii,idim1,idim2)               &
                         + surfbo(idim1,ifac)*surfbo(idim2,ifac)
      enddo

    else

! VECTORISATION NON FORCEE
      do ifac = 1, nfabor
        ii = ifabor(ifac)
        cocg(ii,idim1,idim2) = cocg(ii,idim1,idim2)               &
                         + surfbo(idim1,ifac)*surfbo(idim2,ifac)
      enddo

    endif

  enddo
enddo


!   SYMETRISATION

do iel = 1, ncel
  cocg(iel,2,1) = cocg(iel,1,2)
  cocg(iel,3,1) = cocg(iel,1,3)
  cocg(iel,3,2) = cocg(iel,2,3)
enddo

!===============================================================================
! 2. INVERSION DE LA MATRICE
!===============================================================================


nbloc = ncel/lbloc
if (nbloc.gt.0) then
  do ibloc = 1, nbloc
    do ii = 1, lbloc
      iel = (ibloc-1)*lbloc+ii

      cocg11 = cocg(iel,1,1)
      cocg12 = cocg(iel,1,2)
      cocg13 = cocg(iel,1,3)
      cocg21 = cocg(iel,2,1)
      cocg22 = cocg(iel,2,2)
      cocg23 = cocg(iel,2,3)
      cocg31 = cocg(iel,3,1)
      cocg32 = cocg(iel,3,2)
      cocg33 = cocg(iel,3,3)

      a11=cocg22*cocg33-cocg32*cocg23
      a12=cocg32*cocg13-cocg12*cocg33
      a13=cocg12*cocg23-cocg22*cocg13
      a22=cocg11*cocg33-cocg31*cocg13
      a23=cocg21*cocg13-cocg11*cocg23
      a33=cocg11*cocg22-cocg21*cocg12

      unsdet = 1.d0/(cocg11*a11+cocg21*a12+cocg31*a13)

      aa(ii,1,1) = a11 *unsdet
      aa(ii,1,2) = a12 *unsdet
      aa(ii,1,3) = a13 *unsdet
      aa(ii,2,2) = a22 *unsdet
      aa(ii,2,3) = a23 *unsdet
      aa(ii,3,3) = a33 *unsdet

    enddo

    do ii = 1, lbloc
      iel = (ibloc-1)*lbloc+ii
      cocg(iel,1,1) = aa(ii,1,1)
      cocg(iel,1,2) = aa(ii,1,2)
      cocg(iel,1,3) = aa(ii,1,3)
      cocg(iel,2,2) = aa(ii,2,2)
      cocg(iel,2,3) = aa(ii,2,3)
      cocg(iel,3,3) = aa(ii,3,3)
    enddo

  enddo

endif

irel = mod(ncel,lbloc)
if (irel.gt.0) then
  ibloc = nbloc + 1
  do ii = 1, irel
    iel = (ibloc-1)*lbloc+ii

    cocg11 = cocg(iel,1,1)
    cocg12 = cocg(iel,1,2)
    cocg13 = cocg(iel,1,3)
    cocg21 = cocg(iel,2,1)
    cocg22 = cocg(iel,2,2)
    cocg23 = cocg(iel,2,3)
    cocg31 = cocg(iel,3,1)
    cocg32 = cocg(iel,3,2)
    cocg33 = cocg(iel,3,3)

    a11=cocg22*cocg33-cocg32*cocg23
    a12=cocg32*cocg13-cocg12*cocg33
    a13=cocg12*cocg23-cocg22*cocg13
    a22=cocg11*cocg33-cocg31*cocg13
    a23=cocg21*cocg13-cocg11*cocg23
    a33=cocg11*cocg22-cocg21*cocg12

    unsdet = 1.d0/(cocg11*a11+cocg21*a12+cocg31*a13)

    aa(ii,1,1) = a11 *unsdet
    aa(ii,1,2) = a12 *unsdet
    aa(ii,1,3) = a13 *unsdet
    aa(ii,2,2) = a22 *unsdet
    aa(ii,2,3) = a23 *unsdet
    aa(ii,3,3) = a33 *unsdet

  enddo

  do ii = 1, irel
    iel = (ibloc-1)*lbloc+ii
    cocg(iel,1,1) = aa(ii,1,1)
    cocg(iel,1,2) = aa(ii,1,2)
    cocg(iel,1,3) = aa(ii,1,3)
    cocg(iel,2,2) = aa(ii,2,2)
    cocg(iel,2,3) = aa(ii,2,3)
    cocg(iel,3,3) = aa(ii,3,3)
  enddo
endif


!         MATRICE SYMETRIQUE

do iel = 1, ncel
  cocg(iel,2,1) = cocg(iel,1,2)
  cocg(iel,3,1) = cocg(iel,1,3)
  cocg(iel,3,2) = cocg(iel,2,3)
enddo


!===============================================================================
! 3. CALCUL DU SECOND MEMBRE
!===============================================================================

do iel = 1, ncelet
  bx(iel) = 0.d0
  by(iel) = 0.d0
  bz(iel) = 0.d0
enddo


!     ASSEMBLAGE A PARTIR DES FACETTES FLUIDES

if (ivecti.eq.1) then

!CDIR NODEP
  do ifac = 1,nfac
    ii = ifacel(1,ifac)
    jj = ifacel(2,ifac)
    pfacx = flumas(ifac)*surfac(1,ifac)
    pfacy = flumas(ifac)*surfac(2,ifac)
    pfacz = flumas(ifac)*surfac(3,ifac)
    bx(ii) = bx(ii) + pfacx
    by(ii) = by(ii) + pfacy
    bz(ii) = bz(ii) + pfacz
    bx(jj) = bx(jj) + pfacx
    by(jj) = by(jj) + pfacy
    bz(jj) = bz(jj) + pfacz
  enddo

else

! VECTORISATION NON FORCEE
  do ifac = 1,nfac
    ii = ifacel(1,ifac)
    jj = ifacel(2,ifac)
    pfacx = flumas(ifac)*surfac(1,ifac)
    pfacy = flumas(ifac)*surfac(2,ifac)
    pfacz = flumas(ifac)*surfac(3,ifac)
    bx(ii) = bx(ii) + pfacx
    by(ii) = by(ii) + pfacy
    bz(ii) = bz(ii) + pfacz
    bx(jj) = bx(jj) + pfacx
    by(jj) = by(jj) + pfacy
    bz(jj) = bz(jj) + pfacz
  enddo

endif


!     ASSEMBLAGE A PARTIR DES FACETTES DE BORD

if (ivectb.eq.1) then

!CDIR NODEP
  do ifac = 1,nfabor
    ii = ifabor(ifac)
    bx(ii) = bx(ii) + flumab(ifac)*surfbo(1,ifac)
    by(ii) = by(ii) + flumab(ifac)*surfbo(2,ifac)
    bz(ii) = bz(ii) + flumab(ifac)*surfbo(3,ifac)
  enddo

else

! VECTORISATION NON FORCEE
  do ifac = 1,nfabor
    ii = ifabor(ifac)
    bx(ii) = bx(ii) + flumab(ifac)*surfbo(1,ifac)
    by(ii) = by(ii) + flumab(ifac)*surfbo(2,ifac)
    bz(ii) = bz(ii) + flumab(ifac)*surfbo(3,ifac)
  enddo

endif

!===============================================================================
! 4. RESOLUTION
!===============================================================================


do iel = 1, ncel
  unsrho = 1.d0/rom(iel)
  smbx = bx(iel)
  smby = by(iel)
  smbz = bz(iel)
  ux  (iel) = (cocg(iel,1,1)*smbx+cocg(iel,1,2)*smby              &
              +cocg(iel,1,3)*smbz)*unsrho
  uy  (iel) = (cocg(iel,2,1)*smbx+cocg(iel,2,2)*smby              &
              +cocg(iel,2,3)*smbz)*unsrho
  uz  (iel) = (cocg(iel,3,1)*smbx+cocg(iel,3,2)*smby              &
              +cocg(iel,3,3)*smbz)*unsrho
enddo

!----
! FIN
!----

return

end subroutine
