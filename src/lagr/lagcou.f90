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

subroutine lagcou &
!================

 ( idbia0 , idbra0 ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   ,          &
   nprfml , nvar   , nscal  , nphas  ,                            &
   nbpmax , nvp    , nvp1   , nvep   , nivep  ,                   &
   ntersl , nvlsta , nvisbr ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &
   itepa  , indep  , ibord  ,                                     &
   idevel , ituser , ia     ,                                     &
   volume , rtp    , propce ,                                     &
   ettp   , ettpa  , tepa   , taup   , tempct , tsfext , tslagr , &
   cpgd1  , cpgd2  , cpght  ,                                     &
   tslag  , volp   , volm   ,                                     &
   auxl1  , auxl2  , auxl3  ,                                     &
   rdevel , rtuser , ra     )

!===============================================================================
! FONCTION :
! ----------

!      SOUS-PROGRAMME DU MODULE LAGRANGIEN :
!      -----------------------------------

!     CALCUL DES TERMES SOURCES DU COUPLAGE RETOUR

!     Remarque : les termes sources sont calcules pour
!                la cellule de depart de la particule
!                lors de l'iteration courante. Attention, meme
!                si la particule est sortante du domaine de
!                calcul (peu importe la maniere) on doit calculer
!                un terme source qui correspond a ce qu'echange le
!                fluide porteur et la particule au debut du pas de
!                temps. Si NORDRE = 2 et que la particule est en
!                interaction avec la frontiere, alors les termes
!                source sont calcules comme si NORDRE=1
!                (on oublie le pre-remplissage de TSFEXT dans
!ONFC                 LAGES2).


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
! nvar             ! e  ! <-- ! nombre total de variables                      !
! nscal            ! e  ! <-- ! nombre total de scalaires                      !
! nphas            ! e  ! <-- ! nombre de phases                               !
! nbpmax           ! e  ! <-- ! nombre max de particulies autorise             !
! nvp              ! e  ! <-- ! nombre de variables particulaires              !
! nvp1             ! e  ! <-- ! nvp sans position, vfluide, vpart              !
! nvep             ! e  ! <-- ! nombre info particulaires (reels)              !
! nivep            ! e  ! <-- ! nombre info particulaires (entiers)            !
! ntersl           ! e  ! <-- ! nbr termes sources de couplage retour          !
! nvlsta           ! e  ! <-- ! nombre de var statistiques lagrangien          !
! nvisbr           ! e  ! <-- ! nombre de statistiques aux frontieres          !
! nideve nrdeve    ! e  ! <-- ! longueur de idevel rdevel                      !
! nituse nrtuse    ! e  ! <-- ! longueur de ituser rtuser                      !
! itepa            ! te ! <-- ! info particulaires (entiers)                   !
! (nbpmax,nivep    !    !     !   (cellule de la particule,...)                !
! indep            ! te ! <-- ! pour chaque particule :                        !
!  (nbpmax)        !    !     !    numero de la cellule de depart              !
! ibord            ! te ! <-- ! contient le numero de la                       !
!   (nbpmax)       !    !     !   face d'interaction part/frontiere            !
! idevel(nideve    ! te ! <-- ! tab entier complementaire developemt           !
! ituser(nituse    ! te ! <-- ! tab entier complementaire utilisateur          !
! ia(*)            ! tr ! --- ! macro tableau entier                           !
! volume(ncelet    ! tr ! <-- ! volume d'un des ncelet elements                !
! rtp              ! tr ! <-- ! variables de calcul au centre des              !
! (ncelet,*)       !    !     !    cellules                                    !
! propce           ! tr ! <-- ! proprietes physiques au centre des             !
! (ncelet,*)       !    !     !    cellules                                    !
! ettp             ! tr ! <-- ! tableaux des variables liees                   !
!  (nbpmax,nvp)    !    !     !   aux particules etape courante                !
! ettpa            ! tr ! <-- ! tableaux des variables liees                   !
!  (nbpmax,nvp)    !    !     !   aux particules etape precedente              !
! tepa             ! tr ! <-- ! info particulaires (reels)                     !
! (nbpmax,nvep)    !    !     !   (poids statistiques,...)                     !
! taup(nbpmax)     ! tr ! <-- ! temps caracteristique dynamique                !
! tsfext(nbpmax    ! tr ! <-- ! forces externes                                !
! tempct           ! tr ! <-- ! temps caracteristique thermique                !
!  (nbpmax,2)      !    !     !                                                !
! tslagr(nbpmax    ! tr ! --> ! termes sources de couplage retour              !
!     ntersl)      !    !     !                                                !
! cpgd1,cpgd2,     ! tr ! <-- ! termes de devolatilisation 1 et 2 et           !
!  cpght(nbpmax    !    !     !   de combusion heterogene (charbon             !
!                  !    !     !   avec couplage retour thermique)              !
! tslag(nbpmax,    ! tr ! --- ! tableau de travail                             !
!     ntersl)      !    !     !                                                !
! volp(ncelet)     ! tr ! --- ! fraction volumique des particules              !
! volm(ncelet)     ! tr ! --- ! fraction massique des particules               !
! auxl1(nbpmax)    ! tr ! --- ! tableau de travail                             !
! auxl2(nbpmax)    ! tr ! --- ! tableau de travail                             !
! auxl3(nbpmax)    ! tr ! --- ! tableau de travail                             !
! rdevel(nrdeve    ! tr ! <-- ! tab reel complementaire developemt             !
! rtuser(nrtuse    ! tr ! <-- ! tab reel complementaire utilisateur            !
! ra(*)            ! tr ! --- ! macro tableau reel                             !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail

!===============================================================================

!===============================================================================

implicit none

!===============================================================================
!     DONNEES EN COMMON
!===============================================================================

include "paramx.h"
include "cstnum.h"
include "numvar.h"
include "optcal.h"
include "entsor.h"
include "cstphy.h"
include "pointe.h"
include "period.h"
include "parall.h"
include "lagpar.h"
include "lagran.h"
include "ppppar.h"
include "ppthch.h"
include "ppincl.h"
include "radiat.h"

!===============================================================================

! Arguments

integer          idbia0 , idbra0
integer          ndim   , ncelet , ncel   , nfac   , nfabor
integer          nfml   , nprfml
integer          nvar   , nscal  , nphas
integer          nbpmax , nvp    , nvp1   , nvep  , nivep
integer          ntersl , nvlsta , nvisbr
integer          nideve , nrdeve , nituse , nrtuse
integer          itepa(nbpmax,nivep), indep(nbpmax), ibord(nbpmax)
integer          idevel(nideve), ituser(nituse)
integer          ia(*)

double precision volume(ncelet) , propce(ncelet,*) , rtp(ncelet,*)
double precision ettp(nbpmax,nvp) , ettpa(nbpmax,nvp)
double precision tepa(nbpmax,nvep)
double precision taup(nbpmax) , tempct(nbpmax,2)
double precision tsfext(nbpmax)
double precision cpgd1(nbpmax) , cpgd2(nbpmax) , cpght(nbpmax)
double precision tslag(ncelet,ntersl)
double precision volp(ncelet) , volm(ncelet)
double precision tslagr(ncelet,ntersl)
double precision auxl1(nbpmax) , auxl2(nbpmax) , auxl3(nbpmax)
double precision rdevel(nrdeve) , rtuser(nrtuse)
double precision ra(*)

! VARIABLES LOCALES

integer          npt , iel , ivar , icha , iphas
double precision tvmax , tauv , taum , aux1
double precision uuf , vvf , wwf , mf

!===============================================================================

!===============================================================================
! 1. INITIALISATION
!===============================================================================

iphas = ilphas

tvmax = 0.8d0

!   Nombre de passage pour les termes sources en stationnaire

if (isttio.eq.1 .and. iplas.ge.nstits)  npts = npts + 1

ntxerr = 0
vmax   = 0.d0
tmamax = 0.d0

do iel=1,ncel
  volp(iel) = 0.d0
  volm(iel) = 0.d0
enddo

do ivar = 1,ntersl
  do iel = 1,ncel
    tslag(iel,ivar) = 0.d0
  enddo
enddo

!===============================================================================
! 2. CALCULS PRELIMINAIRES
!===============================================================================

! Finalisation des forces externes (Si la particule a interagit avec
! une frontiere du domaine de calcul, on degenere a l'ordre 1).


do npt = 1,nbpart
  aux1 = dtp/taup(npt)
  if (nordre.eq.1 .or. ibord(npt).gt.0) then
    tsfext(npt)= (1.d0-exp(-aux1)) *ettp(npt,jmp) *taup(npt)
  else
    tsfext(npt) = tsfext(npt)                                     &
                + (1.d0- (1.d0-exp(-aux1)) /aux1 ) * taup(npt)    &
                * ettp(npt,jmp)
  endif
enddo

do npt = 1,nbpart
  auxl1(npt) = tepa(npt,jrpoi)*                                   &
        ( ettp(npt,jmp)  * ettp(npt,jup)                          &
         -ettpa(npt,jmp) * ettpa(npt,jup)                         &
         -gx*tsfext(npt)  ) / dtp
  auxl2(npt) = tepa(npt,jrpoi)*                                   &
        ( ettp(npt,jmp)  * ettp(npt,jvp)                          &
         -ettpa(npt,jmp)* ettpa(npt,jvp)                          &
         -gy*tsfext(npt) ) / dtp
  auxl3(npt) = tepa(npt,jrpoi)*                                   &
        ( ettp(npt,jmp)  * ettp(npt,jwp)                          &
         -ettpa(npt,jmp)* ettpa(npt,jwp)                          &
         -gz*tsfext(npt) ) / dtp
enddo

!===============================================================================
! 3. TERMES SOURCES DE QUANTITE DE MOUVEMENT
!===============================================================================

if (ltsdyn.eq.1) then

  do npt = 1,nbpart

    iel = indep(npt)

! Volume et masse des particules dans la maille

    volp(iel) = volp(iel)                                         &
              + tepa(npt,jrpoi)*pi*(ettpa(npt,jdp)**3)/6.d0
    volm(iel) = volm(iel)                                         &
              + tepa(npt,jrpoi)*ettpa(npt,jmp)

! TS de QM

    tslag(iel,itsvx) = tslag(iel,itsvx) - auxl1(npt)
    tslag(iel,itsvy) = tslag(iel,itsvy) - auxl2(npt)
    tslag(iel,itsvz) = tslag(iel,itsvz) - auxl3(npt)
    tslag(iel,itsli) = tslag(iel,itsli)                           &
                     - 2.d0*tepa(npt,jrpoi)*ettp(npt,jmp)         &
                     / taup(npt)

  enddo

!===============================================================================
! 4. TERMES SOURCES SUR LA TURBULENCE
!===============================================================================

  if (itytur(iphas).eq.2 .or. iturb(iphas).eq.50                  &
       .or. iturb(iphas).eq.60 ) then
! En v2f (ITURB=50) les TS lagrangiens influent uniquement sur k et eps
! (difficile d'ecrire quoi que ce soit sur v2, qui perd son sens de
!  "composante de Rij")

    do npt = 1,nbpart

      iel = indep(npt)

      uuf = 0.5d0 * ( ettpa(npt,juf) + ettp(npt,juf) )
      vvf = 0.5d0 * ( ettpa(npt,jvf) + ettp(npt,jvf) )
      wwf = 0.5d0 * ( ettpa(npt,jwf) + ettp(npt,jwf) )

      tslag(iel,itske) = tslag(iel,itske)                         &
                       - uuf * auxl1(npt)                         &
                       - vvf * auxl2(npt)                         &
                       - wwf * auxl3(npt)

    enddo

    do iel = 1,ncel

      tslag(iel,itske) = tslag(iel,itske)                         &
                       - ettp(npt,juf) * tslag(iel,itsvx)         &
                       - ettp(npt,jvf) * tslag(iel,itsvy)         &
                       - ettp(npt,jwf) * tslag(iel,itsvz)

    enddo

  else if (itytur(iphas).eq.3) then

    do npt = 1,nbpart

      iel = indep(npt)

      uuf = 0.5d0 * ( ettpa(npt,juf) + ettp(npt,juf) )
      vvf = 0.5d0 * ( ettpa(npt,jvf) + ettp(npt,jvf) )
      wwf = 0.5d0 * ( ettpa(npt,jwf) + ettp(npt,jwf) )

      tslag(iel,itsr11) = tslag(iel,itsr11)                       &
                        - 2.d0 * uuf * auxl1(npt)

      tslag(iel,itsr12) = tslag(iel,itsr12)                       &
                        - uuf * auxl2(npt)                        &
                        - vvf * auxl1(npt)

      tslag(iel,itsr13) = tslag(iel,itsr13)                       &
                        - uuf * auxl3(npt)                        &
                        - wwf * auxl1(npt)

      tslag(iel,itsr22) = tslag(iel,itsr22)                       &
                        - 2.d0 * vvf * auxl2(npt)

      tslag(iel,itsr23) = tslag(iel,itsr23)                       &
                        - vvf * auxl3(npt)                        &
                        - wwf * auxl2(npt)

      tslag(iel,itsr33) = tslag(iel,itsr33)                       &
                        - 2.d0 * wwf * auxl3(npt)

    enddo

    do iel = 1,ncel

      tslag(iel,itsr11) = tslag(iel,itsr11)                       &
                 - 2.d0 * rtp(iel,iu(iphas)) * tslag(iel,itsvx)

      tslag(iel,itsr12) = tslag(iel,itsr12)                       &
                        - rtp(iel,iu(iphas)) * tslag(iel,itsvy)   &
                        - rtp(iel,iv(iphas)) * tslag(iel,itsvx)

      tslag(iel,itsr13) = tslag(iel,itsr13)                       &
                        - rtp(iel,iu(iphas)) * tslag(iel,itsvz)   &
                        - rtp(iel,iw(iphas)) * tslag(iel,itsvx)

      tslag(iel,itsr22) = tslag(iel,itsr22)                       &
                 - 2.d0 * rtp(iel,iv(iphas)) * tslag(iel,itsvy)

      tslag(iel,itsr23) = tslag(iel,itsr23)                       &
                        - rtp(iel,iv(iphas)) * tslag(iel,itsvz)   &
                        - rtp(iel,iw(iphas)) * tslag(iel,itsvy)

      tslag(iel,itsr33) = tslag(iel,itsr33)                       &
                 - 2.d0 * rtp(iel,iw(iphas)) * tslag(iel,itsvz)

    enddo

  endif

endif

!===============================================================================
! 5. TERME SOURCE MASSIQUES
!===============================================================================

if ( ltsmas.eq.1 .and. (impvar.eq.1 .or. idpvar.eq.1) ) then

  do npt = 1,nbpart

! Dans saturne TSmasse > 0 ===> Apport de masse sur le fluide

    iel = indep(npt)

    tslag(iel,itsmas) = tslag(iel,itsmas) - tepa(npt,jrpoi)       &
     * ( ettp(npt,jmp) - ettpa(npt,jmp) ) /dtp

  enddo

endif

!===============================================================================
! 6. TERMES SOURCES THERMIQUE
!===============================================================================

if (ltsthe.eq.1) then

  if (iphyla.eq.1 .and. itpvar.eq.1) then

    do npt = 1,nbpart

      iel = indep(npt)

      tslag(iel,itste) = tslag(iel,itste)                         &
     -( ettp(npt,jmp)  *ettp(npt,jtp) *ettp(npt,jcp)              &
        -ettpa(npt,jmp) *ettpa(npt,jtp)                           &
         *ettpa(npt,jcp) ) / dtp * tepa(npt,jrpoi)

      tslag(iel,itsti) = tslag(iel,itsti)                         &
                       + tempct(npt,2) * tepa(npt,jrpoi)

    enddo

    if (iirayo.gt.0) then

      do npt = 1,nbpart

        iel = indep(npt)

        aux1 = pi *ettp(npt,jdp) *ettp(npt,jdp) *tepa(npt,jreps)  &
                *(propce(iel,ipproc(ilumin))                      &
                -4.d0 *stephn *ettp(npt,jtp)**4 )

        tslag(iel,itste) =tslag(iel,itste)+aux1*tepa(npt,jrpoi)

      enddo

    endif

  else if (iphyla.eq.2) then

    do npt = 1,nbpart

      iel = indep(npt)
      icha = itepa(npt,jinch)

      tslag(iel,itste) = tslag(iel,itste)                         &
                -( ettp(npt,jmp)  *ettp(npt,jhp)                  &
                    *ettp(npt,jcp)                                &
                  -ettpa(npt,jmp)*ettpa(npt,jhp)                  &
                    *ettpa(npt,jcp) )                             &
                 /dtp*tepa(npt,jrpoi)

      tslag(iel,itsti) = tslag(iel,itsti)                         &
                       + tempct(npt,2) * tepa(npt,jrpoi)

      tslag(iel,itsmv1(icha)) = tslag(iel,itsmv1(icha))           &
                              + tepa(npt,jrpoi) * cpgd1(npt)

      tslag(iel,itsmv2(icha)) = tslag(iel,itsmv2(icha))           &
                              + tepa(npt,jrpoi) * cpgd2(npt)

      tslag(iel,itsco)  = tslag(iel,itsco)                        &
                        + tepa(npt,jrpoi) * cpght(npt)

      tslag(iel,itsfp4) = 0.d0

    enddo

  endif

endif

!===============================================================================
! 7. Verif que le taux volumique maximal TVMAX admissible de particules
!    ne soit pas depasse dans quelques cellules.
!===============================================================================

do iel = 1,ncel

  mf   = volume(iel) * propce(iel,ipproc(irom(iphas)))
  tauv = volp(iel) / volume(iel)
  taum = volm(iel) / mf

  if (tauv.gt.tvmax) then

    ntxerr = ntxerr + 1

    do ivar =1,ntersl
      tslagr(iel,ivar) = 0.d0
    enddo

  endif

  vmax   = max(tauv,vmax)
  tmamax = max(tmamax,taum)

enddo

!===============================================================================
! 8. MOYENNE TEMPORELLE DES TERMES SOURCES
!===============================================================================

if (isttio.eq.1 .and. npts.gt.0) then

  do ivar = 1,ntersl
    do iel = 1,ncel
      tslagr(iel,ivar) =                                          &
   ( tslag(iel,ivar) + (npts-1.d0)*tslagr(iel,ivar) ) / dble(npts)
    enddo
  enddo

else

  do ivar = 1,ntersl
    do iel = 1,ncel
      tslagr(iel,ivar) = tslag(iel,ivar)
    enddo
  enddo

endif

!===============================================================================

!----
! FIN
!----

end subroutine
