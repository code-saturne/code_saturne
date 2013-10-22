!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2013 EDF S.A.
!
! This program is free software; you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free Software
! Foundation; either version 2 of the License, or (at your option) any later
! version.
!
! This program is distributed in the hope that it will be useful, but WITHOUT
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
! FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
! details.
!
! You should have received a copy of the GNU General Public License along with
! this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
! Street, Fifth Floor, Boston, MA 02110-1301, USA.

!-------------------------------------------------------------------------------

subroutine lagcou &
!================

 ( nbpmax ,                                                       &
   ntersl ,                                                       &
   indep  , ibord  ,                                              &
   rtp    , propce ,                                              &
   taup   , tempct , tsfext ,                                     &
   cpgd1  , cpgd2  , cpght  ,                                     &
   tslag  , volp   , volm   ,                                     &
   auxl1  , auxl2  , auxl3  )

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
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nbpmax           ! e  ! <-- ! nombre max de particulies autorise             !
! ntersl           ! e  ! <-- ! nbr termes sources de couplage retour          !
! indep            ! te ! <-- ! pour chaque particule :                        !
!  (nbpmax)        !    !     !    numero de la cellule de depart              !
! ibord            ! te ! <-- ! contient le numero de la                       !
!   (nbpmax)       !    !     !   face d'interaction part/frontiere            !
! rtp              ! tr ! <-- ! variables de calcul au centre des              !
! (ncelet,*)       !    !     !    cellules                                    !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! taup(nbpmax)     ! tr ! <-- ! temps caracteristique dynamique                !
! tsfext(nbpmax    ! tr ! <-- ! forces externes                                !
! tempct           ! tr ! <-- ! temps caracteristique thermique                !
!  (nbpmax,2)      !    !     !                                                !
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
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use cstnum
use numvar
use optcal
use entsor
use cstphy
use pointe
use parall
use period
use lagpar
use lagran
use ppppar
use ppthch
use ppincl
use radiat
use mesh
use field

!===============================================================================

implicit none

! Arguments

integer          nbpmax
integer          ntersl
integer          indep(nbpmax), ibord(nbpmax)

double precision propce(ncelet,*) , rtp(ncelet,*)
double precision taup(nbpmax) , tempct(nbpmax,2)
double precision tsfext(nbpmax)
double precision cpgd1(nbpmax) , cpgd2(nbpmax) , cpght(nbpmax)
double precision tslag(ncelet,ntersl)
double precision volp(ncelet) , volm(ncelet)
double precision auxl1(nbpmax) , auxl2(nbpmax) , auxl3(nbpmax)

! Local variables

integer          npt , iel , ivar , icha
double precision tvmax , tauv , taum , aux1
double precision uuf , vvf , wwf , mf

double precision, dimension(:), pointer ::  crom

!===============================================================================

!===============================================================================
! 1. INITIALISATION
!===============================================================================

tvmax = 0.8d0

call field_get_val_s(icrom, crom)

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

  if (itytur.eq.2 .or. iturb.eq.50                  &
       .or. iturb.eq.60 ) then
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

      tslag(iel,itske) = tslag(iel,itske)                       &
                       - rtp(iel,iu) * tslag(iel,itsvx)         &
                       - rtp(iel,iv) * tslag(iel,itsvy)         &
                       - rtp(iel,iw) * tslag(iel,itsvz)

    enddo

  else if (itytur.eq.3) then

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
                 - 2.d0 * rtp(iel,iu) * tslag(iel,itsvx)

      tslag(iel,itsr12) = tslag(iel,itsr12)                       &
                        - rtp(iel,iu) * tslag(iel,itsvy)   &
                        - rtp(iel,iv) * tslag(iel,itsvx)

      tslag(iel,itsr13) = tslag(iel,itsr13)                       &
                        - rtp(iel,iu) * tslag(iel,itsvz)   &
                        - rtp(iel,iw) * tslag(iel,itsvx)

      tslag(iel,itsr22) = tslag(iel,itsr22)                       &
                 - 2.d0 * rtp(iel,iv) * tslag(iel,itsvy)

      tslag(iel,itsr23) = tslag(iel,itsr23)                       &
                        - rtp(iel,iv) * tslag(iel,itsvz)   &
                        - rtp(iel,iw) * tslag(iel,itsvy)

      tslag(iel,itsr33) = tslag(iel,itsr33)                       &
                 - 2.d0 * rtp(iel,iw) * tslag(iel,itsvz)

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

    if (nlayer.gt.1) then

        ! Couplage thermique non-fonctionnel en multi-layer

      write(1001,*)
      call csexit (1)

    else

      do npt = 1,nbpart

        iel = indep(npt)
        icha = itepa(npt,jinch)

        tslag(iel,itste) = tslag(iel,itste)                         &
             -( ettp(npt,jmp)  *ettp(npt,jhp(1))                    &
             *ettp(npt,jcp)                                         &
             -ettpa(npt,jmp)*ettpa(npt,jhp(1))                      &
             *ettpa(npt,jcp) )                                      &
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
endif

!===============================================================================
! 7. Verif que le taux volumique maximal TVMAX admissible de particules
!    ne soit pas depasse dans quelques cellules.
!===============================================================================

do iel = 1,ncel

  mf   = volume(iel) * crom(iel)
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
