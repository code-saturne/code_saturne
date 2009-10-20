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

subroutine ppcabs &
!================

 ( idbia0 , idbra0 ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , iphas  ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml , itypfb ,          &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   w1     , w2     , w3     ,                                     &
   rdevel , rtuser ,                                              &
   ra     )

!===============================================================================
! FONCTION :
! ----------


!   SOUS-PROGRAMME PHYSIQUES PARTICULIERES

!  DONNE LA VALEUR DU COEFFICIENT D'ABSORPTION POUR
!    LE MELANGE GAZEUX ET LES PARTICULES POUR LE CP.

!  NOTER QUE IPHAS VAUT 1.

!  DANS LE CAS DU MODELE P-1 ON VERIFIE QUE LA LONGUEUR OPTIQUE
!    DU MILIEU EST AU MINIMUM DE L'ORDRE DE L'UNITE


!-------------------------------------------------------------------------------
!ARGU                             ARGUMENTS
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
! iphas            ! e  ! <-- ! numero de la phase courante                    !
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
! itypfb(nfabor    ! te ! <-- ! type des faces de bord                         !
!  nphas      )    !    !     !                                                !
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
! dt(ncelet)       ! tr ! <-- ! pas de temps                                   !
! rtp, rtpa        ! tr ! <-- ! variables de calcul au centre des              !
! (ncelet,*)       !    !     !    cellules (instant courant ou prec)          !
! propce           ! tr ! <-- ! proprietes physiques au centre des             !
! (ncelet,*)       !    !     !    cellules                                    !
! propfa           ! tr ! <-- ! proprietes physiques au centre des             !
!  (nfac,*)        !    !     !    faces internes                              !
! propfb           ! tr ! <-- ! proprietes physiques au centre des             !
!  (nfabor,*)      !    !     !    faces de bord                               !
! coefa, coefb     ! tr ! <-- ! conditions aux limites aux                     !
!  (nfabor,*)      !    !     !    faces de bord                               !
! w1...3(ncelet    ! tr ! --- ! tableau de travail                             !
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

include "paramx.h"
include "numvar.h"
include "entsor.h"
include "optcal.h"
include "cstphy.h"
include "cstnum.h"
include "pointe.h"
include "ppppar.h"
include "ppthch.h"
include "coincl.h"
include "cpincl.h"
include "fuincl.h"
include "ppincl.h"
include "parall.h"
include "radiat.h"

!===============================================================================

! Arguments

integer          idbia0 , idbra0
integer          ndim   , ncelet , ncel   , nfac   , nfabor
integer          nfml   , nprfml
integer          nnod   , lndfac , lndfbr , ncelbr
integer          nvar   , nscal  , iphas
integer          nideve , nrdeve , nituse , nrtuse

integer          ifacel(2,nfac) , ifabor(nfabor)
integer          ifmfbr(nfabor) , ifmcel(ncelet)
integer          iprfml(nfml,nprfml) , itypfb(nfabor)
integer          ipnfac(nfac+1), nodfac(lndfac)
integer          ipnfbr(nfabor+1), nodfbr(lndfbr)
integer          idevel(nideve), ituser(nituse), ia(*)

double precision xyzcen(ndim,ncelet)
double precision surfac(ndim,nfac), surfbo(ndim,nfabor)
double precision cdgfac(ndim,nfac), cdgfbo(ndim,nfabor)
double precision xyznod(ndim,nnod), volume(ncelet)
double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision w1(ncelet), w2(ncelet), w3(ncelet)

double precision rdevel(nrdeve), rtuser(nrtuse), ra(*)


! VARIABLES LOCALES

integer          idebia, idebra, iel, ifac, icla, ipck, icha, iok
double precision xm, d2, vv, sf, xlc, xkmin, pp

!===============================================================================

!===============================================================================
! 0 - GESTION MEMOIRE
!===============================================================================

idebia = idbia0
idebra = idbra0

!===============================================================================
!  1 - COEFFICIENT D'ABSORPTION DU MELANGE GAZEUX (m-1)
!===============================================================================

if ( ippmod(icod3p).ge.0 .or. ippmod(icoebu).ge.0 ) then

! ----> Combustion gaz : Flamme de diffusion
!                        Flamme de premelange (Modele EBU)

  if (imodak.eq.1) then

    do iel = 1, ncel
      xm = 1.d0/ (  propce(iel,ipproc(iym(1)))/wmolg(1)                 &
                  + propce(iel,ipproc(iym(2)))/wmolg(2)                 &
                  + propce(iel,ipproc(iym(3)))/wmolg(3) )
      w1(iel) = propce(iel,ipproc(iym(3)))*xm/wmolg(3)*xco2
      w2(iel) = propce(iel,ipproc(iym(3)))*xm/wmolg(3)*xh2o
      w3(iel) = 0.d0
    enddo
    call raydak(ncel,ncelet,                                      &
    !==========
      propce(1,ipproc(icak(1))),w1,w2,w3,propce(1,ipproc(itemp)))

    write(NFECRA,*) ' a verifier '
    write(NFECRA,*) ' a finir   : raydak '
    write(NFECRA,*) ' Le codage est a terminer par le groupe I81'
    write(NFECRA,*) '                     13-10-03 22:38:03      '
    call csexit(1)

  else
    do iel = 1, ncel
      propce(iel,ipproc(icak(1))) = propce(iel,ipproc(ickabs))
    enddo
  endif

else if ( ippmod(icp3pl).ge.0 ) then

! ---->  Charbon

  if (imodak.eq.1) then

    do iel = 1,ncel
! concentration volumique en CO2
      w1(iel) = propce(iel,ipproc(immel))/wmole(ico2)             &
               *propce(iel,ipproc(iym1(ico2)))
! concentration volumique en H20
      w2(iel) = propce(iel,ipproc(immel))/wmole(ih2o)             &
               *propce(iel,ipproc(iym1(ih2o)))
! fraction volumique de suies
      w3(iel) = 0.d0

    enddo

    call raydak(ncel,ncelet,                                      &
    !==========
     propce(1,ipproc(icak(1))),w1,w2,w3,propce(1,ipproc(itemp1)))

  else
    do iel = 1, ncel
      propce(iel,ipproc(icak(1))) = ckabs1
    enddo
  endif

else if ( ippmod(icfuel).ge.0 ) then

! ---->  Fuel

  if (imodak.eq.1) then

    do iel = 1,ncel
! concentration volumique en CO2
      w1(iel) = propce(iel,ipproc(immel))/wmole(ico2)             &
               *propce(iel,ipproc(iym1(ico2)))
! concentration volumique en H20
      w2(iel) = propce(iel,ipproc(immel))/wmole(ih2o)             &
               *propce(iel,ipproc(iym1(ih2o)))
! fraction volumique de suies
      w3(iel) = 0.d0

    enddo

    call raydak(ncel,ncelet,                                      &
    !==========
     propce(1,ipproc(icak(1))),w1,w2,w3,propce(1,ipproc(itemp1)))

  else
    do iel = 1, ncel
      propce(iel,ipproc(icak(1))) = ckabs1
    enddo
  endif

endif


!===============================================================================
!  2 - COEFFICIENT D'ABSORPTION DES PARTICULES PAR CLASSE K2/X2 (m-1)
!===============================================================================

! ---->  Charbon

if ( ippmod(icp3pl).ge.0 ) then

  do icla = 1, nclacp

    ipck = 1 + icla
    icha = ichcor(icla)

    do iel = 1, ncel

! ---> Calcul du diametre des particules

      d2 = ( xashch(icha)*diam20(icla)**2 +                       &
           ( 1.d0-xashch(icha))                                   &
             *propce(iel,ipproc(idiam2(icla)))**2 )**0.5d0

! ---> Calcul du coeficient d'absorption des particules K2/X2
!         3./2. ROM/(ROM2*D2)

      propce(iel,ipproc(icak(ipck))) =                            &
                         1.5d0*propce(iel,ipproc(irom(iphas)))    &
                       / ( propce(iel,ipproc(irom2(icla)))*d2)

    enddo

  enddo

endif

! ---->  Fuel

if ( ippmod(icfuel).ge.0 ) then

  do icla = 1, nclafu

    ipck = 1 + icla

    do iel = 1, ncel

! ---> Calcul du coeficient d'absorption des particules K2/X2
!         3./2. ROM/(ROM2*D2)

      propce(iel,ipproc(icak(ipck))) =                            &
                         1.5d0*propce(iel,ipproc(irom(iphas)))    &
                 / ( propce(iel,ipproc(irom3(icla)))              &
                    *propce(iel,ipproc(idiam3(icla))) )

    enddo

  enddo

endif

!===============================================================================
!  3 - COEFFICIENT D'ABSORPTION GAZ (Arc Electrique)
!===============================================================================


if ( ippmod(ielarc).ge.1 ) then

  do iel = 1, ncel

! ---> Directement donne par le fichier dp_elec

      propce(iel,ipproc(icak(1))) = propce(iel,ipproc(idrad))

  enddo

endif

!===============================================================================
!  4 - CLIPPING DU COEFFICIENT D'ABSORPTION DANS LA CAS DE L'APPROX P-1
!===============================================================================


!--> MODELE P-1 : Controle standard des valeurs du coefficient
!                 d'absorption. Ce coefficient doit assurer une
!                 longueur optique au minimum de l'ordre de l'unite.


  if (iirayo.eq.2) then

!         Coefficient d'absorption du melange gaz-particules de charbon

     do iel = 1, ncel
       w3(iel) =  propce(iel,ipproc(icak(1)))
     enddo

     if ( ippmod(icp3pl).ge.0 ) then
       do icla = 1,nclacp
         ipck = 1+icla
         do iel = 1,ncel
           w3(iel) = w3(iel)                                      &
                    + ( propce(iel,ipproc(ix2(icla)))             &
                      * propce(iel,ipproc(icak(ipck))) )
         enddo
       enddo
     elseif ( ippmod(icfuel).ge.0 ) then
       do icla = 1,nclafu
         ipck = 1+icla
         do iel = 1,ncel
           w3(iel) = w3(iel)                                      &
                    + ( propce(iel,ipproc(iyfol(icla)))           &
                      * propce(iel,ipproc(icak(ipck))) )
         enddo
       enddo
     endif

!         Calcul de la longueur caractéristique XLC du domaine de calcul

    sf = 0.d0
    vv = 0.d0

    do ifac = 1,nfabor
       sf = sf + sqrt(surfbo(1,ifac)**2 + surfbo(2,ifac)**2 + surfbo(3,ifac)**2)
    enddo
    if (irangp.ge.0) then
      call parsom(sf)
      !==========
    endif

    do iel = 1,ncel
       vv = vv + volume(iel)
    enddo
    if (irangp.ge.0) then
      call parsom(vv)
      !==========
    endif

    xlc = 3.6d0 * vv / sf

!         Clipping de CK a XKMIN

    xkmin = 1.d0 / xlc

    iok = 0
    do iel = 1,ncel
      if (w3(iel).lt.xkmin) then
        iok = iok +1
      endif
    enddo
    if (irangp.ge.0) then
      call parcpt(iok)
    endif

!     Arret en fin de pas de temps si epaisseur optique trop grande
    pp = xnp1mx/100.0d0
    if (dble(iok).gt.pp*dble(ncelgb)) then
       write(nfecra,1000) xkmin, dble(iok)/dble(ncelgb)*100.d0, xnp1mx
       istpp1 = 1
    endif

  endif

! -------
! FORMAT
! -------

 1000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : RAYONNEMENT APPROXIMATION P-1 (PPCABS)      ',/,&
'@    =========                                               ',/,&
'@                                                            ',/,&
'@    LA LONGUEUR OPTIQUE DU MILIEU SEMI-TRANSPARENT          ',/,&
'@      DOIT AU MOINS ETRE DE L''ORDRE DE L''UNITE POUR ETRE  ',/,&
'@      DANS LE DOMAINE D''APPLICATION DE L''APPROXIMATION P-1',/,&
'@    CELA NE SEMBLE PAS ETRE LE CAS ICI.                     ',/,&
'@                                                            ',/,&
'@    LE COEFFICIENT D''ABSORPTION MINIMUM POUR ASSURER CETTE ',/,&
'@      LONGUEUR OPTIQUE EST XKMIN = ',E10.4                   ,/,&
'@    CETTE VALEUR N''EST PAS ATTEINTE POUR ', E10.4,'%       ',/,&
'@      DES CELLULES DU MAILLAGE.                             ',/,&
'@    LE POURCENTAGE DE CELLULES DU MAILLAGE POUR LESQUELLES  ',/,&
'@      ON ADMET QUE CETTE CONDITION SOIT VIOLEE EST IMPOSE   ',/,&
'@      PAR DEFAUT OU DANS USINI1 A XNP1MX = ', E10.4,'%      ',/,&
'@                                                            ',/,&
'@    Le calcul est interrompu.                               ',/,&
'@                                                            ',/,&
'@    Verifier les valeurs du coefficient d''absorption CK    ',/,&
'@      dans PPCABS, USRAY3 ou Fichier thermochimie.          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

return

end subroutine
