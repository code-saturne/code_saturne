!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2015 EDF S.A.
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

subroutine lagsec                                                              &
!================

 ( npt   ,                                                                     &
   propce , tempct ,                                                           &
   rayon , mlayer , mwater , mwat_max , volume_couche  , sherw , fwat   )


!===============================================================================
! FONCTION :
! ----------

!   SOUS-PROGRAMME DU MODULE LAGRANGIEN :
!   -------------------------------------

!     CALCUL DU FLUX D'EVAPORATION D'UNE PARTICULE PAR UN MODELE D'EQUILIBRE
!       DES PRESSIONS
!         - CALCUL DE LA PRESSION PARTIELLE SATURANTE A LA TEMPERATURE DE LA
!             PARTICULE
!         - CALCUL DU DEBIT DE VAPEUR QUITTANT LA PARTICULE
!
!     LIMITATION EVENTUELLE DU FLUX (LA PARTICULE A UN COMPORTEMENT FIXE AU
!       COURS DU TEMPS: SOIT ELLE VAPORISE, SOIT ELLE CONDENSE)

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! npt              ! e  ! <-- ! numero de la particule a traiter               !
! propce(ncelet, *)! tr ! <-- ! physical properties at cell centers            !
! tempct           ! tr ! <-- ! temps caracteristique thermique                !
!  (nbpart,2)      !    !     !                                                !
! rayon            ! tr ! <-- ! rayons frontieres des differentes couches      !
!  (nlayer)        !    !     !   (en m) (1 par couche)                        !
! mlayer           ! tr ! <-- ! masse des differentes couches (en kg)          !
!  (nlayer)        !    !     !   (1 par couche)                               !
! mwater           ! tr ! --> ! masse d'eau dans la particule (en kg) pour     !
!  (nlayer)        !    !     !   chaque couche                                !
! mwat_max         ! r  ! <-- ! masse maximum d'eau presente sur une couche    !
! volume_couche    ! r  ! <-- ! volume occuppe par une couche (en m^3)         !
! sherw            ! r  ! <-- ! nombre de Sherwood de la particule             !
! fwat             ! tr ! --> ! flux de sechage (en kg/s) pour la chaque       !
!  (nlayer)        !    !     !   couche                                       !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use numvar
use cstnum
use entsor
use lagpar
use lagran
use ppppar
use ppthch
use ppincl
use cpincl
use mesh

!===============================================================================

implicit none

! Arguments
integer          npt

double precision propce(ncelet,*)
double precision tempct(nbpart,2)
double precision rayon(nlayer), mlayer(nlayer), mwater(nlayer)
double precision mwat_max, volume_couche, sherw, fwat(nlayer)

! Local variables
logical          limitateur
integer          ilayer, ilayer_wat, iel, jtshp0
double precision tpk, aux1, aux2, aux3, fwattot , fwat_restant
double precision tsat, fwatsat(nlayer), phith(nlayer), temp(nlayer)
double precision tssauv(nlayer)

double precision precis, lv, tebl, tlimit, tmini
precis = 1.d-15
lv = 2.263d+6
tebl = 100.d0 + tkelvi
tlimit = 302.24d0
tmini = tlimit*(1.d0-tlimit*cs_physical_constants_r/(lv*wmole(ih2o)))

if (associated(ptsvar)) then
  jtshp0 = jhp(1) - 1
else
  jtshp0 = -1
endif

!===============================================================================
! 1. CALCUL DU FLUX DE VAPEUR POUR LA COUCHE ILAYER
!===============================================================================

! --- Initialisation
fwattot = 0.d0
do ilayer=1,nlayer
  fwat(ilayer) = 0.d0
  fwatsat(ilayer) = 0.d0
enddo
iel = ipepa(jisor,npt)

! --- Reperage de la couche
ilayer_wat = 1
do ilayer=1,nlayer
  if (mwater(ilayer).gt.0.0d0) then
    ilayer_wat = ilayer
  endif
enddo

tpk = eptp(jhp(ilayer_wat),npt)

! --- Calcul de la fraction massique d'eau saturante
if (tpk.ge.tmini) then
  if (tpk.ge.tlimit) then
    aux1 = wmole(ih2o) / propce(iel,ipproc(immel))
    aux2 = aux1 * exp(  lv * wmole(ih2o) * (1.0d0/tebl - 1.0d0/tpk) &
                      / cs_physical_constants_r )
  else
    ! On linearise la fraction massique d'eau saturante entre tmini et Tlimit
    ! En Tlimit, la fraction massique d'eau saturante est nulle
    aux1 = wmole(ih2o) / propce(iel,ipproc(immel))
    aux2 =  aux1 * exp(  lv * wmole(ih2o) * (1.0d0/tebl - 1.0d0/tlimit) &
                       / cs_physical_constants_r )                      &
          * (lv*wmole(ih2o) / (cs_physical_constants_r*tlimit**2))      &
          * (tpk - tmini)
  endif
  ! --- Calcul du terme source d eau diffusee
  aux3 = max(1.0d0 - aux2, precis)
  fwattot = pi*eptpa(jdp,npt)*diftl0*sherw*                                    &
               log((1.0d0-propce(iel,ipproc(iym1(ih2o))))/aux3)
else
  ! Le flux est nul
  fwattot = 0.0d0
endif


!===============================================================================
! 2. REPARTITION DE CE FLUX SUR LES CELLULES VOISINES
!===============================================================================
! On repartit ce flux de sechage fwattot sur les differentes cellules voisines

fwat_restant = fwattot

if (fwattot.ge.0.0d0) then
  ! La particule seche vers le coeur
  do ilayer=ilayer_wat,1,-1
    ! On ne peut pas secher plus que l'eau presente
    fwat(ilayer) = min(mwater(ilayer)/dtp,fwat_restant)
    ! Mise a jour du flux restant a evaporer
    fwat_restant = max(0.0d0,fwat_restant-fwat(ilayer))
  enddo
else
  ! La particule se condense vers l'exterieur
  do ilayer=ilayer_wat,nlayer
    if (ilayer.eq.nlayer) then
      ! En nlayer, on condense tout le flux restant
      fwat(ilayer)=fwat_restant
    else
      ! On ne peut pas condenser plus que la qte d'eau sur une couche
      fwat(ilayer) = max(-(mwat_max-mwater(ilayer))/dtp,fwat_restant)
    endif
    ! Flux restant a condenser
    fwat_restant = min(0.0d0,fwat_restant-fwat(ilayer))
  enddo
endif


!===============================================================================
! 3. CALCUL DES FLUX DE SATURATION
!===============================================================================
! Limitation du flux par rapport à la température de saturation
! On limite le flux de sechage pour que, à la fin d'un pas de temps,
! l'enthalpie de la particule soit suffisament élevée pour que sa pression
! saturante en eau soit superieure a la pression partielle d'eau dans l'air
! qui l'entoure

! Calcul de tsat, temperature saturante à la fraction partielle de l'air
if (propce(iel,ipproc(iym1(ih2o))) .gt. precis) then
  aux1 = wmole(ih2o) / propce(iel,ipproc(immel))
  tsat = 1 / (  1/tebl                                    &
              - cs_physical_constants_r                   &
               *log(propce(iel,ipproc(iym1(ih2o)))/aux1)  &
               /(lv * wmole(ih2o)) )
  if (tsat .lt. tlimit) then
    tsat = tmini + propce(iel,ipproc(iym1(ih2o))) / (aux1                  &
                   *exp( lv*wmole(ih2o)*(1.0d0/tebl-1.0d0/tlimit)          &
                        /cs_physical_constants_r)                          &
                   *(lv*wmole(ih2o)/(cs_physical_constants_r*tlimit**2)) )
  endif
else
  tsat = tmini
endif

! On calcule la temperature a la fin du pas de temps sans chimie
do ilayer=1,nlayer
  ! On debranche tous les termes sources thermiques volumiques pour ce calcul
  phith(ilayer) = 0.0d0
enddo

! On sauvegarde le tableau de correction pour le 2e ordre
if (jtshp0.ge.0) then
  do ilayer=1,nlayer
    tssauv(ilayer) = ptsvar(jtshp0+ilayer,npt)
  enddo
endif

call lagtmp                                                                    &
!==========
( npt    ,                                                                     &
  propce , tempct ,                                                            &
  rayon  , mlayer , phith , temp  , volume_couche )

! On remet le tableau de correction pour le 2e ordre
if (jtshp0.ge.0) then
  do ilayer=1,nlayer
    ptsvar(jtshp0+ilayer,npt) = tssauv(ilayer)
  enddo
endif

! On calcule le flux d'evaporation/condensation tel que T_i=Tsat
do ilayer=1,nlayer
  fwatsat(ilayer) = mlayer(ilayer)*eptpa(jcp,npt)                              &
                    *(temp(ilayer)-tsat)/(lv*dtp)
enddo

!===============================================================================
! 4. LIMITATION EVENTUELLE DE VAPEUR
!===============================================================================
! On compare les resultats de 2. et 3. et on arbitre

limitateur = .false.
if (fwattot.ge.0.0d0) then

  ! La particule seche vers le coeur
  do ilayer=nlayer,1,-1
    if (limitateur.eqv..false.) then

      ! On vérifie que la couche n'a pas un comportement opposé
      if (fwatsat(ilayer).lt.0.0d0) then
        ! On bloque toutes les couches suivantes
        limitateur = .true.
      endif

      ! On limite le flux
      if (fwat(ilayer).gt.fwatsat(ilayer)) then
        fwat(ilayer) = max(0.0d0 , fwatsat(ilayer))
      endif

    else
      ! Le limitateur bloque les couches plus internes
      fwat(ilayer) = 0.0d0
    endif
  enddo

else if (fwattot.lt.0.0d0) then

  ! On vérifie que les couches ext n'ont pas un comportement opposé
  do ilayer=nlayer,ilayer_wat
    if (fwatsat(ilayer).gt.0.0d0) then
      ! On bloque toutes les couches suivantes
      limitateur = .true.
    endif
  enddo

  ! La particule se condense vers l'exterieur
  do ilayer=ilayer_wat,nlayer
    if (limitateur.eqv..false.) then

      ! On limite le flux
      if (fwatsat(ilayer).gt.fwat(ilayer)) then
        fwat(ilayer) = min(0.0d0 , fwatsat(ilayer))
      endif

    else
      ! Le limitateur bloque
      fwat(ilayer) = 0.0d0
    endif
  enddo

endif


!----
! End
!----

end subroutine
