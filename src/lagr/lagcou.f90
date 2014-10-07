!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2014 EDF S.A.
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

 ( ntersl ,                                                       &
   propce ,                                                       &
   taup   , tempct , tsfext ,                                     &
   cpgd1  , cpgd2  , cpght  ,                                     &
   volp   , volm   )

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
! ntersl           ! e  ! <-- ! nbr termes sources de couplage retour          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! taup(nbpart)     ! tr ! <-- ! temps caracteristique dynamique                !
! tsfext(nbpart)   ! tr ! <-- ! forces externes                                !
! tempct           ! tr ! <-- ! temps caracteristique thermique                !
!  (nbpart,2)      !    !     !                                                !
! cpgd1,cpgd2,     ! tr ! <-- ! termes de devolatilisation 1 et 2 et           !
!  cpght(nbpart)   !    !     !   de combusion heterogene (charbon             !
!                  !    !     !   avec couplage retour thermique)              !
! volp(ncelet)     ! tr ! --- ! fraction volumique des particules              !
! volm(ncelet)     ! tr ! --- ! fraction massique des particules               !
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

integer          ntersl

double precision propce(ncelet,*)
double precision taup(nbpart) , tempct(nbpart,2)
double precision tsfext(nbpart)
double precision cpgd1(nbpart) , cpgd2(nbpart) , cpght(nbpart)
double precision volp(ncelet) , volm(ncelet)

! Local variables

integer          npt , iel , ivar , icha
double precision tvmax , tauv , taum , aux1
double precision uuf , vvf , wwf , mf

double precision, dimension(:), pointer ::  crom
double precision, dimension(:,:), pointer :: vel
double precision, allocatable, dimension(:,:) :: tslag
double precision, allocatable, dimension(:) :: auxl1, auxl2, auxl3

!===============================================================================

! Map field arrays
call field_get_val_v(ivarfl(iu), vel)

!===============================================================================
! 1. INITIALISATION
!===============================================================================

tvmax = 0.8d0

call field_get_val_s(icrom, crom)

allocate(tslag(ncelet,ntersl))
allocate(auxl1(nbpart) , auxl2(nbpart) , auxl3(nbpart))

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
  if (nordre.eq.1 .or. ipepa(jord1,npt).gt.0) then
    tsfext(npt)= (1.d0-exp(-aux1)) *eptp(jmp,npt) *taup(npt)
  else
    tsfext(npt) = tsfext(npt)                                     &
                + (1.d0- (1.d0-exp(-aux1)) /aux1 ) * taup(npt)    &
                * eptp(jmp,npt)
  endif
enddo

do npt = 1,nbpart
  auxl1(npt) = pepa(jrpoi,npt)*                                   &
        ( eptp(jmp,npt)  * eptp(jup,npt)                          &
         -eptpa(jmp,npt) * eptpa(jup,npt)                         &
         -gx*tsfext(npt)  ) / dtp
  auxl2(npt) = pepa(jrpoi,npt)*                                   &
        ( eptp(jmp,npt)  * eptp(jvp,npt)                          &
         -eptpa(jmp,npt)* eptpa(jvp,npt)                          &
         -gy*tsfext(npt) ) / dtp
  auxl3(npt) = pepa(jrpoi,npt)*                                   &
        ( eptp(jmp,npt)  * eptp(jwp,npt)                          &
         -eptpa(jmp,npt)* eptpa(jwp,npt)                          &
         -gz*tsfext(npt) ) / dtp
enddo

!===============================================================================
! 3. TERMES SOURCES DE QUANTITE DE MOUVEMENT
!===============================================================================

if (ltsdyn.eq.1) then

  do npt = 1,nbpart

    iel = ipepa(jisora,npt)

! Volume et masse des particules dans la maille

    volp(iel) = volp(iel)                                         &
              + pepa(jrpoi,npt)*pi*(eptpa(jdp,npt)**3)/6.d0
    volm(iel) = volm(iel)                                         &
              + pepa(jrpoi,npt)*eptpa(jmp,npt)

! TS de QM

    tslag(iel,itsvx) = tslag(iel,itsvx) - auxl1(npt)
    tslag(iel,itsvy) = tslag(iel,itsvy) - auxl2(npt)
    tslag(iel,itsvz) = tslag(iel,itsvz) - auxl3(npt)
    tslag(iel,itsli) = tslag(iel,itsli)                           &
                     - 2.d0*pepa(jrpoi,npt)*eptp(jmp,npt)         &
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

      iel = ipepa(jisora,npt)

      uuf = 0.5d0 * ( eptpa(juf,npt) + eptp(juf,npt) )
      vvf = 0.5d0 * ( eptpa(jvf,npt) + eptp(jvf,npt) )
      wwf = 0.5d0 * ( eptpa(jwf,npt) + eptp(jwf,npt) )

      tslag(iel,itske) = tslag(iel,itske)                         &
                       - uuf * auxl1(npt)                         &
                       - vvf * auxl2(npt)                         &
                       - wwf * auxl3(npt)

    enddo

    do iel = 1,ncel

      tslag(iel,itske) = tslag(iel,itske)                       &
                       - vel(1,iel) * tslag(iel,itsvx)         &
                       - vel(2,iel) * tslag(iel,itsvy)         &
                       - vel(3,iel) * tslag(iel,itsvz)

    enddo

  else if (itytur.eq.3) then

    do npt = 1,nbpart

      iel = ipepa(jisora,npt)

      uuf = 0.5d0 * ( eptpa(juf,npt) + eptp(juf,npt) )
      vvf = 0.5d0 * ( eptpa(jvf,npt) + eptp(jvf,npt) )
      wwf = 0.5d0 * ( eptpa(jwf,npt) + eptp(jwf,npt) )

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
                 - 2.d0 * vel(1,iel) * tslag(iel,itsvx)

      tslag(iel,itsr12) = tslag(iel,itsr12)                       &
                        - vel(1,iel) * tslag(iel,itsvy)   &
                        - vel(2,iel) * tslag(iel,itsvx)

      tslag(iel,itsr13) = tslag(iel,itsr13)                       &
                        - vel(1,iel) * tslag(iel,itsvz)   &
                        - vel(3,iel) * tslag(iel,itsvx)

      tslag(iel,itsr22) = tslag(iel,itsr22)                       &
                 - 2.d0 * vel(2,iel) * tslag(iel,itsvy)

      tslag(iel,itsr23) = tslag(iel,itsr23)                       &
                        - vel(2,iel) * tslag(iel,itsvz)   &
                        - vel(3,iel) * tslag(iel,itsvy)

      tslag(iel,itsr33) = tslag(iel,itsr33)                       &
                 - 2.d0 * vel(3,iel) * tslag(iel,itsvz)

    enddo

  endif

endif

!===============================================================================
! 5. TERME SOURCE MASSIQUES
!===============================================================================

if ( ltsmas.eq.1 .and. (impvar.eq.1 .or. idpvar.eq.1) ) then

  do npt = 1,nbpart

! Dans saturne TSmasse > 0 ===> Apport de masse sur le fluide

    iel = ipepa(jisora,npt)

    tslag(iel,itsmas) = tslag(iel,itsmas) - pepa(jrpoi,npt)       &
     * ( eptp(jmp,npt) - eptpa(jmp,npt) ) /dtp

  enddo

endif

!===============================================================================
! 6. TERMES SOURCES THERMIQUE
!===============================================================================

if (ltsthe.eq.1) then

  if (iphyla.eq.1 .and. itpvar.eq.1) then

    do npt = 1,nbpart

      iel = ipepa(jisora,npt)

      tslag(iel,itste) = tslag(iel,itste)                         &
     -( eptp(jmp,npt)  *eptp(jtp,npt) *eptp(jcp,npt)              &
        -eptpa(jmp,npt) *eptpa(jtp,npt)                           &
         *eptpa(jcp,npt) ) / dtp * pepa(jrpoi,npt)

      tslag(iel,itsti) = tslag(iel,itsti)                         &
                       + tempct(npt,2) * pepa(jrpoi,npt)

    enddo

    if (iirayo.gt.0) then

      do npt = 1,nbpart

        iel = ipepa(jisora,npt)

        aux1 = pi *eptp(jdp,npt) *eptp(jdp,npt) *pepa(jreps,npt)  &
                *(propce(iel,ipproc(ilumin))                      &
                -4.d0 *stephn *eptp(jtp,npt)**4 )

        tslag(iel,itste) =tslag(iel,itste)+aux1*pepa(jrpoi,npt)

      enddo

    endif

  else if (iphyla.eq.2) then

    if (nlayer.gt.1) then

        ! Couplage thermique non-fonctionnel en multi-layer

      write(1001,*)
      call csexit (1)

    else

      do npt = 1,nbpart

        iel = ipepa(jisora,npt)
        icha = ipepa(jinch,npt)

        tslag(iel,itste) = tslag(iel,itste)                         &
             -( eptp(jmp,npt)  *eptp(jhp(1),npt)                    &
             *eptp(jcp,npt)                                         &
             -eptpa(jmp,npt)*eptpa(jhp(1),npt)                      &
             *eptpa(jcp,npt) )                                      &
             /dtp*pepa(jrpoi,npt)

        tslag(iel,itsti) = tslag(iel,itsti)                         &
             + tempct(npt,2) * pepa(jrpoi,npt)

        tslag(iel,itsmv1(icha)) = tslag(iel,itsmv1(icha))           &
             + pepa(jrpoi,npt) * cpgd1(npt)

        tslag(iel,itsmv2(icha)) = tslag(iel,itsmv2(icha))           &
             + pepa(jrpoi,npt) * cpgd2(npt)

        tslag(iel,itsco)  = tslag(iel,itsco)                        &
             + pepa(jrpoi,npt) * cpght(npt)

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

deallocate(auxl1, auxl2, auxl3)
deallocate(tslag)

!----
! FIN
!----

end subroutine lagcou
