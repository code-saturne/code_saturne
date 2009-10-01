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

subroutine cplphy &
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
   w1     , w2     , w3     , w4     ,                            &
   w5     , w6     , w7     , w8     ,                            &
   rdevel , rtuser , ra     )

!===============================================================================
! FONCTION :
! --------

!   SOUS-PROGRAMME DU MODULE LAGRANGIEN COUPLE CHARBON PULVERISE :
!   --------------------------------------------------------------

!    ROUTINE UTILISATEUR POUR PHYSIQUE PARTICULIERE

!      COMBUSTION EULERIENNE DE CHARBON PULVERISE ET
!      TRANSPORT LAGRANGIEN DES PARTICULES DE CHARBON

! Calcul de RHO de la phase gazeuse


! ATTENTION :
! =========


! Il est INTERDIT de modifier la viscosite turbulente VISCT ici
!        ========
!  (une routine specifique est dediee a cela : usvist)


!  Il FAUT AVOIR PRECISE ICP(IPHAS) = 1
!     ==================
!    dans usini1 si on souhaite imposer une chaleur specifique
!    CP variable pour la phase IPHAS (sinon: ecrasement memoire).


!  Il FAUT AVOIR PRECISE IVISLS(Numero de scalaire) = 1
!     ==================
!     dans usini1 si on souhaite une diffusivite VISCLS variable
!     pour le scalaire considere (sinon: ecrasement memoire).




! Remarques :
! ---------

! Cette routine est appelee au debut de chaque pas de temps

!    Ainsi, AU PREMIER PAS DE TEMPS (calcul non suite), les seules
!    grandeurs initialisees avant appel sont celles donnees
!      - dans usini1 :
!             . la masse volumique (initialisee a RO0(IPHAS))
!             . la viscosite       (initialisee a VISCL0(IPHAS))
!      - dans usppiv :
!             . les variables de calcul  (initialisees a 0 par defaut
!             ou a la valeur donnee dans usiniv)

! On peut donner ici les lois de variation aux cellules
!     - de la masse volumique                      ROM    kg/m3
!         (et eventuellememt aux faces de bord     ROMB   kg/m3)
!     - de la viscosite moleculaire                VISCL  kg/(m s)
!     - de la chaleur specifique associee          CP     J/(kg degres)
!     - des "diffusivites" associees aux scalaires VISCLS kg/(m s)


! On dispose des types de faces de bord au pas de temps
!   precedent (sauf au premier pas de temps, ou les tableaux
!   ITYPFB et ITRIFB n'ont pas ete renseignes)


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
! lndfac           ! e  ! <-- ! longueur du tableau nodfac                     !
! lndfbr           ! e  ! <-- ! longueur du tableau nodfbr                     !
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
!   (lndfac)       !    !     !  face interne dans nodfac                      !
! nodfac           ! te ! <-- ! connectivite faces internes/noeuds             !
!   (nfac+1)       !    !     !                                                !
! ipnfbr           ! te ! <-- ! position du premier noeud de chaque            !
!   (lndfbr)       !    !     !  face de bord dans nodfbr                      !
! nodfbr           ! te ! <-- ! connectivite faces de bord/noeuds              !
!   (nfabor+1)     !    !     !                                                !
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
! xyznod           ! tr ! <-- ! coordonnes des noeuds                          !
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
include "parall.h"
include "ppppar.h"
include "ppthch.h"
include "coincl.h"
include "cpincl.h"
include "ppincl.h"
include "ppcpfu.h"

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
double precision w1(ncelet),w2(ncelet),w3(ncelet),w4(ncelet)
double precision w5(ncelet),w6(ncelet),w7(ncelet),w8(ncelet)
double precision rdevel(nrdeve), rtuser(nrtuse), ra(*)

! VARIABLES LOCALES

integer          idebia, idebra
integer          ntbcpi, icpwi, ntbcpr, icpwr
integer          ntbmci, imcwi, ntbmcr, imcwr
integer          ntbwoi, iwori, ntbwor, iworr
integer          ifinia, ifinra
integer          iel, icha, iphas, ipcrom
integer          izone, ifac
integer          ipbrom
double precision srrom1
double precision wmolme

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

! --- Initialisation des tableaux de travail

do iel = 1, ncel
  w1(iel) = zero
  w2(iel) = zero
  w3(iel) = zero
  w4(iel) = zero
  w5(iel) = zero
  w6(iel) = zero
  w7(iel) = zero
  w8(iel) = zero
enddo

!===============================================================================
! 2. CALCUL DES PROPRIETES PHYSIQUES DE LA PHASE GAZEUSE
!                    VALEURS CELLULES
!                    ----------------
!    TEMPERATURE
!    MASSE VOLUMIQUE
!    CONCENTRATIONS DES ESPECES GAZEUSES
!===============================================================================

! --- Calcul de l'enthalpie du gaz     dans W8
!            de F1M                    dans W2
!            de F2M                    dans W3
!            de F3M                    dans W4
!            de F4M                    dans W5
!            de F4P2M                  dans W7


! ---- W2 = F1M = SOMME(F1M(ICHA))
!      W3 = F2M = SOMME(F2M(ICHA))
!      W4 = F3M
!      W5 = F4M = 1. - F1M - F2M - F3M
!      W7 = F4P2M

do icha = 1, ncharb
  do iel = 1, ncel
    w2(iel) =  w2(iel) + rtp(iel,isca(if1m(icha)))
    w3(iel) =  w3(iel) + rtp(iel,isca(if2m(icha)))
  enddo
enddo

do iel = 1,ncel
  w4(iel) =  rtp(iel,isca(if3m))
  w5(iel) = 1.d0 - w2(iel) - w3(iel) -w4(iel)
  w7(iel) =  rtp(iel,isca(if4p2m))
  w8(iel) = rtp(iel,isca(ihm))
enddo

! --- Gestion memoire
!     Autres tableaux

! ------ Macro tableau d'entiers TBCPI : NTBCPI
!        Macro tableau de reels  TBCPR : NTBCPR
!        Macro tableau d'entiers TBMCI : NTBMCI
!        Macro tableau de reels  TBMCR : NTBMCR
!        Macro tableau d'entiers TBWOI : NTBWOI
!        Macro tableau de reels  TBWOR : NTBWOR

ntbcpi = 1
ntbcpr = 9
ntbmci = 0
ntbmcr = 2 + 2*ncharb + 4

!  Ce sont en fait X1M, X2M,
!                  F1M(ICHA) et F2M(ICHA) pour chaque charbon
!                  ACHX1F1, ACHX2F2, ACOF1, ACOF2
ntbwoi = 1
ntbwor = 4

call memcp1                                                       &
!==========
 ( idebia , idebra ,                                              &
   nvar   , ncelet , ncel   , nfac   , nfabor ,                   &
   ntbcpi , icpwi  ,                                              &
   ntbcpr , icpwr  ,                                              &
   ntbmci , imcwi  ,                                              &
   ntbmcr , imcwr  ,                                              &
   ntbwoi , iwori  ,                                              &
   ntbwor , iworr  ,                                              &
   ifinia , ifinra )

call cplph1                                                       &
!==========
 ( ifinia , ifinra ,                                              &
   ncelet , ncel   ,                                              &
   ntbcpi , ntbcpr , ntbmci , ntbmcr , ntbwoi , ntbwor ,          &
   w2     , w3     , w4     , w5     , w6     , w7     ,          &
!         F1M      F2M      F3M      F4M      F3P2M    F4P2M
   w8     ,                                                       &
!         ENTH
   rtp    , propce  , w1 ,                                        &
!                          ----
!                 ATTENTION W1 contient RHO1
   ia(icpwi) , ra(icpwr) ,                                        &
   ia(imcwi) , ra(imcwr) ,                                        &
   ia(iwori) , ra(iworr)  )


!===============================================================================
! 3. Relaxation de la masse volumique de la phase gazeuse
!===============================================================================

! --- Calcul de Rho avec relaxation

iphas = 1
ipcrom = ipproc(irom(iphas))

if (ipass.gt.1.or.(isuite.eq.1.and.initro(iphas).eq.1)) then
  srrom1 = srrom
else
  srrom1 = 1.d0
endif


do iel = 1, ncel
! ---- Sous relaxation eventuelle a donner dans ppini1.F
  propce(iel,ipcrom) = srrom1*propce(iel,ipcrom)                  &
                     + (1.d0-srrom1)*w1(iel)
enddo


!===============================================================================
! 4. CALCUL DE RHO DE LA PHASE GAZEUSE

!                             VALEURS FACES
!                             -------------
!===============================================================================

iphas = 1
ibrom(iphas) = 1
ipbrom = ipprob(irom(iphas))
ipcrom = ipproc(irom(iphas))

! ---> Masse volumique au bord pour toutes les faces
!      Les faces d'entree seront recalculees.

do ifac = 1, nfabor
  iel = ifabor(ifac)
  propfb(ifac,ipbrom) = propce(iel,ipcrom)
enddo

! ---> Masse volumique au bord pour les faces d'entree UNIQUEMENT
!     Le test sur IZONE sert pour les reprises de calcul

if ( ipass.gt.1 .or. isuite.eq.1 ) then
  do ifac = 1, nfabor

    izone = izfppp(ifac)
    if(izone.gt.0) then
      if ( ientat(izone).eq.1 ) then
        wmolme = (1.d0+xsi) / (wmole(io2)+xsi*wmole(in2))
        propfb(ifac,ipbrom) = p0(iphas)                           &
                             /(wmolme*rr*timpat(izone))
      endif
    endif

  enddo
endif

!----
! FIN
!----

return
end subroutine
