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

subroutine ebuphy &
!================

 ( idbia0 , idbra0 ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   nideve , nrdeve , nituse , nrtuse , nphmx  ,                   &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr , ibrom  , izfppp ,          &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   yfuegf , yoxygf , yprogf ,                                     &
   yfuegb , yoxygb , yprogb ,                                     &
   temp   , masmel ,                                              &
   rdevel , rtuser , ra     )

!===============================================================================
! FONCTION :
! --------

! ROUTINE PHYSIQUE PARTICULIERE : FLAMME DE PREMELAMGE MODELE EBU
! Calcul de RHO adiabatique ou permeatique (transport de H)


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
! nphmx            ! e  ! <-- ! nphsmx                                         !
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
! ibrom            ! te ! <-- ! indicateur de remplissage de romb              !
!   (nphmx   )     !    !     !                                                !
! izfppp           ! te ! <-- ! numero de zone de la face de bord              !
! (nfabor)         !    !     !  pour le module phys. part.                    !
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
! w1...8(ncelet    ! tr ! --- ! tableau de travail                             !
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
include "numvar.h"
include "optcal.h"
include "cstphy.h"
include "cstnum.h"
include "entsor.h"
include "ppppar.h"
include "ppthch.h"
include "coincl.h"
include "cpincl.h"
include "ppincl.h"
include "radiat.h"

!===============================================================================

! Arguments

integer          idbia0 , idbra0
integer          ndim   , ncelet , ncel   , nfac   , nfabor
integer          nfml   , nprfml
integer          nnod   , lndfac , lndfbr , ncelbr
integer          nvar   , nscal  , nphas
integer          nideve , nrdeve , nituse , nrtuse , nphmx

integer          ifacel(2,nfac) , ifabor(nfabor)
integer          ifmfbr(nfabor) , ifmcel(ncelet)
integer          iprfml(nfml,nprfml)
integer          ipnfac(nfac+1), nodfac(lndfac)
integer          ipnfbr(nfabor+1), nodfbr(lndfbr), ibrom(nphmx)
integer          izfppp(nfabor)
integer          idevel(nideve), ituser(nituse), ia(*)

double precision xyzcen(ndim,ncelet)
double precision surfac(ndim,nfac), surfbo(ndim,nfabor)
double precision cdgfac(ndim,nfac), cdgfbo(ndim,nfabor)
double precision xyznod(ndim,nnod), volume(ncelet)
double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision yfuegf(ncelet),yoxygf(ncelet),yprogf(ncelet)
double precision yfuegb(ncelet),yoxygb(ncelet),yprogb(ncelet)
double precision temp(ncelet),masmel(ncelet)
double precision rdevel(nrdeve), rtuser(nrtuse), ra(*)

! VARIABLES LOCALES

integer          idebia, idebra
integer          ipctem, ipcfue, ipcoxy, ipcpro
integer          igg, iphas, iel, ipcrom
integer          ipckab, ipt4, ipt3
integer          ifac, izone
integer          ipbrom, ipbycg, ipcycg, mode
double precision coefg(ngazgm), ygfm, ygbm, epsi
double precision nbmol , temsmm , fmel , ckabgf, ckabgb
double precision masmgb, hgb, tgb, masmgf, masmg

integer       ipass
data          ipass /0/
save          ipass

!===============================================================================

!===============================================================================
! 0. ON COMPTE LES PASSAGES
!===============================================================================

ipass = ipass + 1

!===============================================================================
! 1. INITIALISATIONS A CONSERVER
!===============================================================================

! --- Initialisation memoire

idebia = idbia0
idebra = idbra0

! ---> Initialisation

do igg = 1, ngazgm
  coefg(igg) = zero
enddo

! ---> Positions des variables, coefficients

iphas = 1
ipcrom = ipproc(irom(iphas))
ipbrom = ipprob(irom(iphas))
ipctem = ipproc(itemp)
ipcfue = ipproc(iym(1))
ipcoxy = ipproc(iym(2))
ipcpro = ipproc(iym(3))
if ( iirayo.gt.0 ) then
  ipckab = ipproc(ickabs)
  ipt4 = ipproc(it4m)
  ipt3 = ipproc(it3m)
endif


!===============================================================================
! 2. DETERMINATION DES GRANDEURS THERMOCHIMIQUES
!===============================================================================

! ---> Grandeurs GAZ FRAIS

! ----- Fournies par l'utilisateur
!       FMEL         --> Taux de melange
!                        constant pour les options 0 et 1
!                        variable sinon
!       TGF          --> Temperature gaz frais en K identique
!                        pour premelange frais et dilution
! ----- Deduites
!       YFUEGF(    .)    --> Fraction massique fuel gaz frais
!       YOXYGF(    .)    --> Fraction massique oxydant gaz frais
!       YPROGF(    .)    --> Fraction massique produits gaz frais
!       HGF          --> Enthalpie massique gaz frais identique
!                        pour premelange frais et dilution
!       MASMGF       --> Masse molaire gaz frais
!       CKABGF       --> Coefficient d'absorption

! ---> Grandeurs GAZ BRULES

! ----- Deduites
!       TGB          --> Temperature gaz brules en K
!       YFUEGB(    .)    --> Fraction massique fuel gaz brules
!       YOXYGB(    .)    --> Fraction massique oxydant gaz brules
!       YPROGB(    .)    --> Fraction massique produits gaz brules
!       MASMGB       --> Masse molaire gaz brules
!       CKABGB       --> Coefficient d'absorption

! ---> Grandeurs MELANGE

!       MASMEL           --> Masse molaire du melange
!       PROPCE(    .,IPCTEM) --> Temperature du melange
!       PROPCE(    .,IPCROM) --> Masse volumique du melange
!       PROPCE(,.F,O,P ) --> Fractions massiques en F, O, P
!       PROPCE(    .,IPCKAB) --> Coefficient d'absorption
!       PROPCE(    .,IPT4  ) --> terme T^4
!       PROPCE(    .,IPT3  ) --> terme T^3


! ---> Fractions massiques des gaz frais et brules en F, O, P

do iel = 1, ncel

  if ( ippmod(icoebu).eq.0 .or. ippmod(icoebu).eq.1 ) then
    fmel = frmel
  else
    fmel = rtp(iel,isca(ifm))
  endif

  yfuegf(iel) = fmel
  yoxygf(iel) = 1.d0-fmel
  yprogf(iel) = 0.d0

  yfuegb(iel) = max(zero,(fmel-fs(1))/(1.d0-fs(1)))
  yprogb(iel) = (fmel-yfuegb(iel))/fs(1)
  yoxygb(iel) = 1.d0 - yfuegb(iel) - yprogb(iel)

enddo

epsi = 1.d-06

do iel = 1, ncel

! ---> Coefficients d'absorption des gaz frais et brules

  if ( iirayo.gt.0 ) then
     ckabgf = yfuegf(iel)*ckabsg(1) + yoxygf(iel)*ckabsg(2)       &
            + yprogf(iel)*ckabsg(3)
     ckabgb = yfuegb(iel)*ckabsg(1) + yoxygb(iel)*ckabsg(2)       &
            + yprogb(iel)*ckabsg(3)
  endif

! ---> Masse molaire des gaz frais

  coefg(1) = yfuegf(iel)
  coefg(2) = yoxygf(iel)
  coefg(3) = yprogf(iel)
  nbmol = 0.d0
  do igg = 1, ngazg
    nbmol = nbmol + coefg(igg)/wmolg(igg)
  enddo
  masmgf = 1.d0/nbmol

! ---> Calcul de l'enthalpie des gaz frais

  mode    = -1
  call cothht                                                     &
  !==========
  ( mode   , ngazg , ngazgm  , coefg  ,                           &
    npo    , npot   , th     , ehgazg ,                           &
    hgf    , tgf    )

! ---> Masse molaire des gaz brules

  coefg(1) = yfuegb(iel)
  coefg(2) = yoxygb(iel)
  coefg(3) = yprogb(iel)
  nbmol = 0.d0
  do igg = 1, ngazg
    nbmol = nbmol + coefg(igg)/wmolg(igg)
  enddo
  masmgb = 1.d0/nbmol

  ygfm = rtp(iel,isca(iygfm))
  ygbm = 1.d0 - ygfm

! ---> Masse molaire du melange

  masmel(iel) = 1.d0 / ( ygfm/masmgf + ygbm/masmgb )

! ---> Calcul Temperature des gaz brules

  if ( ippmod(icoebu).eq.0 .or. ippmod(icoebu).eq.2 ) then
! ---- EBU Standard et modifie en conditions adiabatiques (sans H)
    hgb = hgf
  else if ( ippmod(icoebu).eq.1 .or. ippmod(icoebu).eq.3 ) then
! ---- EBU Standard et modifie en conditions permeatiques (avec H)
    hgb = hgf
    if ( ygbm.gt.epsi ) then
      hgb = ( rtp(iel,isca(ihm))-hgf*ygfm ) / ygbm
    endif
  endif

  mode = 1
  call cothht                                                     &
  !==========
  ( mode   , ngazg , ngazgm  , coefg  ,                           &
    npo    , npot   , th     , ehgazg ,                           &
    hgb    , tgb    )

  if ( ippmod(icoebu).eq.0 .or. ippmod(icoebu).eq.2 ) then
! ---- EBU Standard et modifie en conditions adiabatiques (sans H)
    tgbad = tgb
  endif

! ---> Temperature du melange
!      Rq PPl : Il serait plus judicieux de ponderer par les CP (GF et GB)
  propce(iel,ipctem) = ygfm*tgf + ygbm*tgb

! ---> Temperature / Masse molaire

  temsmm = ygfm*tgf/masmgf + ygbm*tgb/masmgb

! ---> Masse volumique du melange

  if (ipass.gt.1.or.(isuite.eq.1.and.initro(iphas).eq.1)) then
    propce(iel,ipcrom) = srrom*propce(iel,ipcrom)                 &
                       + (1.d0-srrom)*                            &
                         ( p0(iphas)/(rr*temsmm) )
  endif

! ---> Fractions massiques des especes globales

  propce(iel,ipcfue) = yfuegf(iel)*ygfm + yfuegb(iel)*ygbm
  propce(iel,ipcoxy) = yoxygf(iel)*ygfm + yoxygb(iel)*ygbm
  propce(iel,ipcpro) = yprogf(iel)*ygfm + yprogb(iel)*ygbm

! ---> Grandeurs relatives au rayonnement

  if ( iirayo.gt.0 ) then
    propce(iel,ipckab) = ygfm*ckabgf + ygbm*ckabgb
    propce(iel,ipt4)   = ygfm*tgf**4 + ygbm*tgb**4
    propce(iel,ipt3)   = ygfm*tgf**3 + ygbm*tgb**3
  endif

enddo


!===============================================================================
! 3. CALCUL DE RHO ET DES FRACTIONS MASSIQUES DES ESPECES GLOBALES
!    SUR LES BORDS
!===============================================================================

! --> Masse volumique au bord

iphas = 1
ibrom(iphas) = 1

! ---- Masse volumique au bord pour toutes les facettes
!      Les facettes d'entree seront recalculees apres
!      a partir des CL (si IPASS > 2).

! ---- Au premier passage sans suite ou si on n'a pas relu la
!      masse volumique dans le fichier suite, on n'a pas recalcule la
!      masse volumique ci-dessus, pas la peine de la reprojeter aux
!      faces.

if (ipass.gt.1.or.(isuite.eq.1.and.initro(iphas).eq.1)) then

  do ifac = 1, nfabor
    iel = ifabor(ifac)
    propfb(ifac,ipbrom) = propce(iel,ipcrom)
  enddo

endif


! ---- Masse volumique au bord pour les facettes d'entree UNIQUEMENT
!      Le test sur IZONE sert pour les reprises de calcul
!      On suppose implicitement que les elements ci-dessus ont ete relus
!      dans le fichier suite (i.e. pas de suite en combustion d'un calcul
!      a froid) -> sera pris en compte eventuellement dans les versions
!      suivantes

if ( ipass.gt.1 .or. isuite.eq.1 ) then
  do ifac = 1, nfabor
    izone = izfppp(ifac)
    if(izone.gt.0) then
      if ( ientgb(izone).eq.1 .or. ientgf(izone).eq.1 ) then
        coefg(1) = fment(izone)
        coefg(2) = 1.d0-fment(izone)
        coefg(3) = zero
        if ( ientgb(izone).eq.1 ) then
          coefg(1) = max(zero,(fment(izone)-fs(1))/(1.d0-fs(1)))
          coefg(3) = (fment(izone)-coefg(1))/fs(1)
          coefg(2) = 1.d0 - coefg(1) - coefg(3)
        endif
        nbmol = 0.d0
        do igg = 1, ngazg
          nbmol = nbmol + coefg(igg)/wmolg(igg)
        enddo
       masmg = 1.d0/nbmol
       temsmm = tkent(izone)/masmg
       propfb(ifac,ipbrom) = p0(iphas)/(rr*temsmm)
      endif
    endif
  enddo
endif


! --> Fractions massiques des especes globales au bord
!     Uniquement si transport de H

if ( iirayo.gt.0 ) then
  do igg = 1, ngazg
    ipbycg = ipprob(iym(igg))
    ipcycg = ipproc(iym(igg))
    do ifac = 1, nfabor
      iel = ifabor(ifac)
      propfb(ifac,ipbycg) = propce(iel,ipcycg)
    enddo
  enddo
endif

!===============================================================================
! FORMATS
!----


!----
! FIN
!----

return
end
