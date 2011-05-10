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

subroutine lagout &
!================

 ( idbia0 , idbra0 ,                                              &
   lndnod ,                                                       &
   nvar   , nscal  , nphas  ,                                     &
   nbpmax , nvp    , nvp1   , nvep   , nivep  ,                   &
   ntersl , nvlsta , nvisbr ,                                     &
   icocel , itycel , itepa  ,                                     &
   ia     ,                                                       &
   dt     , rtpa   , rtp    , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   ettp   , tepa   , parbor , statis , stativ , tslagr ,          &
   ra     )

!===============================================================================
! FONCTION :
! ----------

!   SOUS-PROGRAMME DU MODULE LAGRANGIEN :
!   -------------------------------------

! 1. Ecriture du fichier suite 'lagava' :
!     * variables sur les particules (ETTP)
!     * informations sur les particules (ITEPA, TEPA)

! 2. Ecriture du fichier suite statistiques et termes sources
!     'lasava' :
!     * statistiques volumiques (STATIS)
!     * statistiques aux frontieres (PARBOR)
!     * termes sources de couplage retour (TSLAGR)

! 3. Finalisation des sorties graphiques

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
! lndnod           ! e  ! <-- ! dim. connectivite cellules->faces              !
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! nphas            ! i  ! <-- ! number of phases                               !
! nbpmax           ! e  ! <-- ! nombre max de particulies autorise             !
! nvp              ! e  ! <-- ! nombre de variables particulaires              !
! nvp1             ! e  ! <-- ! nvp sans position, vfluide, vpart              !
! nvep             ! e  ! <-- ! nombre info particulaires (reels)              !
! nivep            ! e  ! <-- ! nombre info particulaires (entiers)            !
! ntersl           ! e  ! <-- ! nbr termes sources de couplage retour          !
! nvlsta           ! e  ! <-- ! nombre de var statistiques lagrangien          !
! nvisbr           ! e  ! <-- ! nombre de statistiques aux frontieres          !
! icocel           ! te ! <-- ! connectivite cellules -> faces                 !
! (lndnod)         !    !     !    face de bord si numero negatif              !
! itycel           ! te ! <-- ! connectivite cellules -> faces                 !
! (ncelet+1)       !    !     !    pointeur du tableau icocel                  !
! itepa            ! te ! <-- ! info particulaires (entiers)                   !
! (nbpmax,nivep    !    !     !   (cellule de la particule,...)                !
! ia(*)            ! ia ! --- ! main integer work array                        !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
! ettp             ! tr ! <-- ! tableaux des variables liees                   !
!  (nbpmax,nvp)    !    !     !   aux particules etape courante                !
! tepa             ! tr ! <-- ! info particulaires (reels)                     !
! (nbpmax,nvep)    !    !     !   (poids statistiques,...)                     !
! statis           ! tr ! <-- ! cumul pour les moyennes des                    !
!(ncelet,nvlsta    !    !     !   statistiques volumiques                      !
! stativ           ! tr ! <-- ! cumul pour les variances des                   !
!(ncelet,          !    !     !    statistiques volumiques                     !
!   nvlsta-1)      !    !     !                                                !
! tslagr           ! tr ! <-- ! terme de couplage retour du                    !
!(ncelet,ntersl    !    !     !   lagrangien sur la phase porteuse             !
! parbor           ! tr ! <-- ! infos sur interaction des particules           !
!(nfabor,nvisbr    !    !     !   aux faces de bord                            !
! ra(*)            ! ra ! --- ! main real work array                           !
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
use cpincl
use radiat
use mesh

!===============================================================================

implicit none

! Arguments

integer          idbia0 , idbra0
integer          lndnod
integer          nvar   , nscal  , nphas
integer          nbpmax , nvp    , nvp1   , nvep  , nivep
integer          ntersl , nvlsta , nvisbr

integer          icocel(lndnod) , itycel(ncelet+1)
integer          itepa(nbpmax,nivep)
integer          ia(*)

double precision dt(ncelet) , rtp(ncelet,*) , rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*) , propfb(nfabor,*)
double precision coefa(nfabor,*) , coefb(nfabor,*)
double precision ettp(nbpmax,nvp) ,  tepa(nbpmax,nvep)
double precision statis(ncelet,nvlsta)
double precision stativ(ncelet,nvlsta-1)
double precision tslagr(ncelet,ntersl)
double precision parbor(nfabor,nvisbr)
double precision ra(*)

! Local variables


character        rubriq*64 , car4*4
character        nomnvl(nvplmx)*60 , nomtsl(nvplmx)*60
character        nomite(nvplmx)*64 , nomrte(nvplmx)*64
character        ficsui*32
integer          idebia , idebra
integer          ifinia , ifinra
integer          ierror , irtyp  , itysup , nbval
integer          ivers  , ilecec
integer          nfin   , iforce , icha   , ii
integer          itrav1
integer          ipas   , jj
integer          impavl , impvls

!===============================================================================
!===============================================================================
! 0.  GESTION MEMOIRE
!===============================================================================

idebia = idbia0
idebra = idbra0

!===============================================================================
! 1. ECRITURE DU FICHIER SUITE : VARIABLES LIEES AUX PARTICULES
!===============================================================================

! ---> Ouverture (et on saute si erreur)
!     ILECEC = 2 : ecriture

write(nfecra,6010)

ilecec = 2
ficsui = 'lagrangian'
call opnsui(ficsui, len(ficsui), ilecec, impavl, ierror)
!==========
if (ierror.ne.0) then
  write(nfecra,9010) ficsui
  goto 9998
endif

write(nfecra,6011)


! Entete et Infos sur le calcul ou on saute si erreur
!     On inclut une rubrique destinee a distinguer ce fichier
!       d'un autre fichier suite
!     Pour le moment, IVERS n'est pas utilise

itysup = 0
nbval  = 1

ivers  = 111
RUBRIQ = 'version_fichier_suite_Lagrangien_variables'
irtyp  = 1
call ecrsui(impavl,rubriq,len(rubriq),itysup,nbval,irtyp,ivers,   &
            ierror)
if (ierror.ne.0) then
  write(nfecra,9010)
  goto 9998
endif

! Temps (par securite)

RUBRIQ = 'nombre_iterations_Lagrangiennes'
irtyp  = 1
call ecrsui(impavl,rubriq,len(rubriq),itysup,nbval,irtyp,iplas,   &
            ierror)
if(ierror.ne.0) then
  write(nfecra,9020)                                              &
  'nombre_iterations_Lagrangiennes                             ', &
  'IPLAS', IPLAS
endif

RUBRIQ = 'temps_physique_Lagrangien'
irtyp  = 2
call ecrsui(impavl,rubriq,len(rubriq),itysup,nbval,irtyp,ttclag,  &
            ierror)
if(ierror.ne.0) then
  write(nfecra,9021)                                              &
  'temps_physique_Lagrangien                                   ', &
  'TTCLAG', TTCLAG
endif

! Infos sur le suivi du calcul

irtyp  = 1

RUBRIQ = 'nombre_courant_particules'
call ecrsui(impavl,rubriq,len(rubriq),itysup,nbval,irtyp,nbpart,  &
            ierror)
if(ierror.ne.0) then
  write(nfecra,9030)                                              &
  'nombre_courant_particules                                   ', &
  'NBPART', NBPART
  goto 9998
endif

RUBRIQ = 'nombre_total_particules'
call ecrsui(impavl,rubriq,len(rubriq),itysup,nbval,irtyp,nbptot,  &
            ierror)
if(ierror.ne.0) then
  write(nfecra,9020)                                              &
  'nombre_total_particules                                     ', &
  'NBPTOT', NBPTOT
endif

RUBRIQ = 'nombre_particules_perdues'
call ecrsui(impavl,rubriq,len(rubriq),itysup,nbval,irtyp,nbpert,  &
            ierror)
if(ierror.ne.0) then
  write(nfecra,9020)                                              &
  'nombre_particules_perdues                                   ', &
  'NBPERT', NBPERT
endif

RUBRIQ = 'indicateur_physique_particules'
call ecrsui(impavl,rubriq,len(rubriq),itysup,nbval,irtyp,iphyla,  &
            ierror)
if(ierror.ne.0) then
  write(nfecra,9030)                                              &
  'indicateur_physique_particules                              ', &
  'IPHYLA', IPHYLA
  goto 9998
endif

RUBRIQ = 'indicateur_temperature_particules'
call ecrsui(impavl,rubriq,len(rubriq),itysup,nbval,irtyp,itpvar,  &
            ierror)
if(ierror.ne.0) then
  write(nfecra,9030)                                              &
  'indicateur_temperature_particules                           ', &
  'ITPVAR', ITPVAR
  goto 9998
endif

RUBRIQ = 'indicateur_diametre_particules'
call ecrsui(impavl,rubriq,len(rubriq),itysup,nbval,irtyp,idpvar,  &
            ierror)
if(ierror.ne.0) then
  write(nfecra,9020)                                              &
  'indicateur_diametre_particules                              ', &
  'IDPVAR', IDPVAR
endif

RUBRIQ = 'indicateur_masse_particules'
call ecrsui(impavl,rubriq,len(rubriq),itysup,nbval,irtyp,impvar,  &
            ierror)
if(ierror.ne.0) then
  write(nfecra,9020)                                              &
  'indicateur_masse_particules                                 ', &
  'IMPVAR', IMPVAR
endif

RUBRIQ = 'nombre_variables_utilisateur'
call ecrsui(impavl,rubriq,len(rubriq),itysup,nbval,irtyp,nvls,    &
            ierror)
if(ierror.ne.0) then
  write(nfecra,9020)                                              &
  'nombre_variables_utilisateur                                ', &
  'NVLS', NVLS
endif

write(nfecra,6012)

! Variables particulaires

NOMNVL(JXP) = 'variable_positionX_particule'
NOMNVL(JYP) = 'variable_positionY_particule'
NOMNVL(JZP) = 'variable_positionZ_particule'
NOMNVL(JUP) = 'variable_vitesseU_particule'
NOMNVL(JVP) = 'variable_vitesseV_particule'
NOMNVL(JWP) = 'variable_vitesseW_particule'
NOMNVL(JUF) = 'variable_vitesseU_fluide_vu'
NOMNVL(JVF) = 'variable_vitesseV_fluide_vu'
NOMNVL(JWF) = 'variable_vitesseW_fluide_vu'
NOMNVL(JMP) = 'variable_masse_particule'
NOMNVL(JDP) = 'variable_diametre_particule'
if (iphyla.eq.1 .and. itpvar.eq.1) then
  NOMNVL(JTP) = 'variable_temperature_particule'
  NOMNVL(JTF) = 'variable_temperature_fluide_vu'
  NOMNVL(JCP) = 'variable_chaleur_specifique_particule'
elseif (iphyla.eq.2) then
  NOMNVL(JHP) = 'variable_temperature_particule'
  NOMNVL(JTF) = 'variable_temperature_fluide_vu'
  NOMNVL(JMCH) = 'variable_masse_charbon_reactif'
  NOMNVL(JMCK) = 'variable_masse_coke'
  NOMNVL(JCP) = 'variable_chaleur_specifique_particule'
endif
if (nvls.gt.0) then
  do ii = 1,nvls
    WRITE(CAR4,'(I4.4)') II
    NOMNVL(JVLS(II)) = 'variable_supplementaire_'//CAR4
  enddo
endif

itysup = 0
nbval  = nbpart
irtyp  = 2

do ii = jmp,jwf
  rubriq = nomnvl(ii)
  call ecrsui(impavl,rubriq,len(rubriq),itysup,nbval,irtyp,       &
              ettp(1,ii),ierror)
  if(ierror.ne.0) then
!         advienne que pourra sur le format
    write(nfecra,9100) rubriq
    goto 9998
  endif
enddo

do ii = 1,jmp-1
  rubriq = nomnvl(ii)
  call ecrsui(impavl,rubriq,len(rubriq),itysup,nbval,irtyp,       &
              ettp(1,ii),ierror)
  if(ierror.ne.0) then
!         advienne que pourra sur le format
    write(nfecra,9101) rubriq
  endif
enddo

! Caracteristiques et infos particulaires (ENTIERS)

NOMITE(JISOR) = 'numero_cellule_particules'
if (iphyla.eq.2) then
  NOMITE(JINCH) = 'numero_charbon'
endif

itysup = 0
nbval  = nbpart
irtyp  = 1

rubriq = nomite(jisor)
call ecrsui(impavl,rubriq,len(rubriq),itysup,nbval,irtyp,         &
            itepa(1,jisor),ierror)
if(ierror.ne.0) then
!         advienne que pourra sur le format
  write(nfecra,9100) rubriq
  goto 9998
endif

! groupe statistique particules

if (nbclst .gt. 0 ) then
  NOMITE(JCLST) = 'numero_groupe_statistiques'

  itysup = 0
  nbval  = nbpart
  irtyp  = 1

  rubriq = nomite(jclst)
  call ecrsui(impavl,rubriq,len(rubriq),itysup,nbval,irtyp,       &
              itepa(1,jclst),ierror)
  if(ierror.ne.0) then
!           advienne que pourra sur le format
    write(nfecra,9100) rubriq
    goto 9998
  endif
endif

! Numero du charbon des particules

if (iphyla.eq.2) then
  rubriq = nomite(jinch)
  call ecrsui(impavl,rubriq,len(rubriq),itysup,nbval,irtyp,       &
              itepa(1,jinch),ierror)
  if(ierror.ne.0) then
!         advienne que pourra sur le format
    write(nfecra,9101) rubriq
  endif
endif

! Caracteristiques et infos particulaires (REELS)

NOMRTE(JRTSP) = 'temps_sejour_particules'
NOMRTE(JRPOI) = 'poids_statistiques_particules'
if (iphyla.eq.1 .and. itpvar.eq.1 .and.iirayo.gt.0) then
  NOMRTE(JREPS) = 'emissivite_particules'
endif
if (iphyla.eq.2) then
  NOMRTE(JRDCK) = 'diametre_coeur_retrecissant_charbon'
  NOMRTE(JRD0P) = 'diametre_initial_charbon'
  NOMRTE(JRR0P) = 'masse_volumique_initial_charbon'
endif

itysup = 0
nbval  = nbpart
irtyp  = 2

do ii = 1, nvep
  rubriq = nomrte(ii)
  call ecrsui(impavl,rubriq,len(rubriq),itysup,nbval,irtyp,       &
              tepa(1,ii),ierror)
  if(ierror.ne.0) then
!         advienne que pourra sur le format
    write(nfecra,9101) rubriq
  endif
enddo

write(nfecra,6013)

! ---> Fermeture du fichier suite
call clssui(impavl,ierror)

if (ierror.ne.0) then
  write(nfecra,9140) ficavl
endif

! ---> En cas d'erreur, on continue quand meme
 9998 continue

write(nfecra,6014)

!===============================================================================
! 2. ECRITURE DU FICHIER SUITE STATISTIQUES ET TERMES SOURCES
!    DE COUPLAGE RETOUR
!===============================================================================


if ( (istala.eq.1 .and. iplas.ge.idstnt) .or.                     &
      iilagr.eq.2                        .or.                     &
     (iensi3.eq.1 .and. nvisbr.gt.0)          ) then


! ---> Ouverture (et on saute si erreur)
!     ILECEC = 2 : ecriture

  write(nfecra,7010)

  ilecec = 2
  ficsui = 'lagrangian_stats'
  call opnsui(ficsui, len(ficsui), ilecec, impvls, ierror)
  !==========
  if (ierror.ne.0) then
    write(nfecra,9510) ficsui
    goto 9999
  endif

  write(nfecra,7011)


! Entete et Infos sur le calcul ou on saute si erreur
!     On inclut une rubrique destinee a distinguer ce fichier
!       d'un autre fichier suite
!     Pour le moment, IVERS n'est pas utilise

  itysup = 0
  nbval  = 1

  ivers  = 111
  RUBRIQ = 'version_fichier_suite_Lagrangien_statistiques'
  irtyp  = 1
  call ecrsui(impvls,rubriq,len(rubriq),itysup,nbval,irtyp,ivers, &
              ierror)

  if(ierror.ne.0) then
    write(nfecra,9510)
    goto 9999
  endif

! ---> On ecrit ISTTIO c'est utile dans tous les cas

  RUBRIQ = 'indicateur_ecoulement_stationnaire'
  irtyp  = 1
  call ecrsui(impvls,rubriq,len(rubriq),itysup,nbval,irtyp,       &
              isttio, ierror)

  if(ierror.ne.0) then
    write(nfecra,9510)
    goto 9999
  endif

! --> En premier, on ecrit les statistiques volumiques

  if (istala.eq.1 .and. iplas.ge.idstnt) then

    RUBRIQ = 'iteration_debut_statistiques'
    irtyp  = 1
    call ecrsui(impvls,rubriq,len(rubriq),itysup,nbval,irtyp,     &
                idstnt,ierror)
    if(ierror.ne.0) then
      write(nfecra,9520)                                          &
  'iteration_debut_statistiques                                ', &
      'IDSTNT', IDSTNT
    endif

    RUBRIQ = 'iteration_debut_statistiques_stationnaires'
    irtyp  = 1
    call ecrsui(impvls,rubriq,len(rubriq),itysup,nbval,irtyp,     &
                nstist,ierror)
    if(ierror.ne.0)  then
      write(nfecra,9520)                                          &
  'iteration_debut_statistiques_stationnaires                  ', &
      'NSTIST', NSTIST
    endif

    RUBRIQ = 'nombre_iterations_statistiques_stationnaires'
    irtyp  = 1
    call ecrsui(impvls,rubriq,len(rubriq),itysup,nbval,irtyp,     &
                npst,ierror)
    if(ierror.ne.0) then
      write(nfecra,9520)                                          &
  'nombre_iterations_statistiques_stationnaires                ', &
      'NPST', NPST
    endif

    RUBRIQ = 'temps_statistiques_stationnaires'
    irtyp  = 2
    call ecrsui(impvls,rubriq,len(rubriq),itysup,nbval,irtyp,     &
                tstat,ierror)
    if(ierror.ne.0) then
      write(nfecra,9521)                                          &
  'temps_statistiques_stationnaires                            ', &
      'TSTAT', TSTAT
    endif

    RUBRIQ = 'classe_statistique_particules'
    irtyp  = 1
    call ecrsui(impvls,rubriq,len(rubriq),itysup,nbval,irtyp,     &
                nbclst,ierror)
    if(ierror.ne.0) then
      write(nfecra,9520)                                          &
  'classes_statistiques                                        ', &
      'NBCLST', NBCLST
    endif

    RUBRIQ = 'nombre_statistiques_utilisateur'
    irtyp  = 1
    call ecrsui(impvls,rubriq,len(rubriq),itysup,nbval,irtyp,     &
                nvlsts,ierror)
    if(ierror.ne.0) then
      write(nfecra,9520)                                          &
  'nombre_statistiques_utilisateur                             ', &
     'NVLSTS', NVLSTS
    endif

!  Statistiques volumiques

    itysup = 1
    irtyp  = 2
    nbval  = 1

    do ipas  = 0,nbclst
      do jj = 1,nvlsta

        ii = ipas*nvlsta +jj
        if (ipas.gt.0) then
          WRITE(CAR4,'(I4.4)') IPAS
          RUBRIQ = 'moy_stat_vol_groupe_'//CAR4//'_'//NOMLAG(II)
        else
          RUBRIQ = 'moy_stat_vol_'//NOMLAG(II)
        endif
        call ecrsui(impvls,rubriq,len(rubriq),itysup,nbval,irtyp, &
                  statis(1,ii),ierror)

        if(ierror.ne.0) then
!         advienne que pourra sur le format
          write(nfecra,9550) rubriq
        endif
      enddo

      do jj = 1,nvlsta-1

        ii = ipas*nvlsta +jj
        if (ipas.gt.0) then
          WRITE(CAR4,'(I4.4)') IPAS
          RUBRIQ = 'var_stat_vol_groupe_'//CAR4//'_'//NOMLAV(II)
        else
          RUBRIQ = 'var_stat_vol_'//NOMLAV(II)
        endif
        call ecrsui(impvls,rubriq,len(rubriq),itysup,nbval,irtyp, &
                  stativ(1,ii),ierror)

        if(ierror.ne.0) then
!         advienne que pourra sur le format
          write(nfecra,9550) rubriq
        endif

      enddo

    enddo

  endif

! --> En second, c'est le tour des statistiques aux frontieres

  if (iensi3.eq.1 .and. nvisbr.gt.0) then

    itysup = 0
    nbval  = 1

    RUBRIQ = 'iteration_debut_stats_frontieres_stationnaires'
    irtyp  = 1
    call ecrsui(impvls,rubriq,len(rubriq),itysup,nbval,irtyp,     &
                nstbor,ierror)
    if(ierror.ne.0)  then
      write(nfecra,9520)                                          &
  'iteration_debut_stats_frontieres_stationnaires              ', &
      'NSTBOR', NSTBOR
    endif

    RUBRIQ = 'nombre_iterations_stats_frontieres'
    irtyp  = 1
    call ecrsui(impvls,rubriq,len(rubriq),itysup,nbval,irtyp,     &
                npstft,ierror)
    if(ierror.ne.0) then
      write(nfecra,9520)                                          &
  'nombre_iterations_stats_frontieres                          ', &
      'NPSTFT', NPSTFT
    endif

    RUBRIQ = 'nombre_iterations_stats_frontieres_stationnaires'
    irtyp  = 1
    call ecrsui(impvls,rubriq,len(rubriq),itysup,nbval,irtyp,     &
                npstf,ierror)
    if(ierror.ne.0) then
      write(nfecra,9520)                                          &
  'nombre_iterations_stats_frontieres_stationnaires            ', &
      'NPSTF', NPSTF
    endif

    RUBRIQ = 'temps_stats_frontieres_stationnaires'
    irtyp  = 2
    call ecrsui(impvls,rubriq,len(rubriq),itysup,nbval,irtyp,     &
                tstatp,ierror)
    if(ierror.ne.0) then
      write(nfecra,9521)                                          &
  'temps_stats_frontieres_stationnaires                        ', &
      'TSTATP', TSTATP
    endif

    RUBRIQ = 'nombre_stats_frontieres_utilisateur'
    irtyp  = 1
    call ecrsui(impvls,rubriq,len(rubriq),itysup,nbval,irtyp,     &
                nusbor,ierror)
    if(ierror.ne.0) then
      write(nfecra,9521)                                          &
  'nombre_stats_frontieres_utilisateur                         ', &
      'NUSBOR', NUSBOR
    endif

!  Statistiques aux frontieres

    itysup = 3
    nbval  = 1
    irtyp  = 2

    do ii = 1,nvisbr
      RUBRIQ = 'stat_bord_'//NOMBRD(II)
      call ecrsui(impvls,rubriq,len(rubriq),itysup,nbval,irtyp,   &
                  parbor(1,ii),ierror)
      if(ierror.ne.0) then
!         advienne que pourra sur le format
        write(nfecra,9550) rubriq
      endif
    enddo

  endif

! --> Enfin, en cas de couplage retour, on ecrit les termes sources

  if (iilagr.eq.2) then

    itysup = 0
    nbval  = 1

    RUBRIQ = 'iteration_debut_termes_sources_stationnaires'
    irtyp  = 1
    call ecrsui(impvls,rubriq,len(rubriq),itysup,nbval,irtyp,     &
                nstits,ierror)
    if(ierror.ne.0) then
      write(nfecra,9520)                                          &
  'iteration_debut_termes_sources_stationnaires                ', &
      'NSTITS', NSTITS
    endif

    RUBRIQ = 'nombre_iterations_termes_sources_stationnaires'
    irtyp  = 1
    call ecrsui(impvls,rubriq,len(rubriq),itysup,nbval,irtyp,     &
                npts,ierror)
    if(ierror.ne.0) then
      write(nfecra,9520)                                          &
  'nombre_iterations_termes_sources_stationnaires              ', &
      'NPTS', NPTS
    endif

    RUBRIQ = 'modele_turbulence_termes_sources'
    irtyp  = 1
    call ecrsui(impvls,rubriq,len(rubriq),itysup,nbval,irtyp,     &
                iturb,ierror)
    if(ierror.ne.0) then
      write(nfecra,9520)                                          &
  'modele_turbulence_termes_sources                            ', &
      'ITURB',ITURB
    endif

!       On donne des labels au different TS pour les noms de rubriques
!       On donne le meme label au keps, au v2f et au k-omega (meme variable k)

    if (ltsdyn.eq.1) then
      NOMTSL(ITSVX) = 'terme_source_vitesseX'
      NOMTSL(ITSVY) = 'terme_source_vitesseY'
      NOMTSL(ITSVZ) = 'terme_source_vitesseZ'
      NOMTSL(ITSLI) = 'terme_source_vitesse_implicite'
      if (itytur.eq.2 .or. iturb.eq.50              &
           .or. iturb.eq.60) then
        NOMTSL(ITSKE) = 'terme_source_turbulence_keps'
      else if (itytur.eq.3) then
        NOMTSL(ITSR11) = 'terme_source_turbulence_R11'
        NOMTSL(ITSR12) = 'terme_source_turbulence_R12'
        NOMTSL(ITSR13) = 'terme_source_turbulence_R13'
        NOMTSL(ITSR22) = 'terme_source_turbulence_R22'
        NOMTSL(ITSR23) = 'terme_source_turbulence_R23'
        NOMTSL(ITSR33) = 'terme_source_turbulence_R33'
      endif
    endif
    if (ltsmas.eq.1) then
      NOMTSL(ITSMAS) = 'terme_source_masse'
    endif
    if (ltsthe.eq.1) then
      if (iphyla.eq.1 .and. itpvar.eq.1) then
        NOMTSL(ITSTE) = 'terme_source_thermique_explicite'
        NOMTSL(ITSTI) = 'terme_source_thermique_implicite'
      else if (iphyla.eq.2) then
        NOMTSL(ITSTE) = 'terme_source_thermique_explicite'
        NOMTSL(ITSTI) = 'terme_source_thermique_implicite'
        do icha = 1,ncharb
          WRITE(CAR4,'(I4.4)') ICHA
          NOMTSL(ITSMV1(ICHA)) = 'terme_source_legeres_F1_'//CAR4
          NOMTSL(ITSMV2(ICHA)) = 'terme_source_lourdes_F2_'//CAR4
        enddo
        NOMTSL(ITSCO) = 'terme_source_F3'
        NOMTSL(ITSFP4) = 'terme_source_variance_traceur_air'
      endif
    endif

!  Termes source de couplage retour

    itysup = 1
    nbval  = 1
    irtyp  = 2

    do ii = 1,ntersl
      rubriq = nomtsl(ii)
      call ecrsui(impvls,rubriq,len(rubriq),itysup,nbval,irtyp,   &
                  tslagr(1,ii),ierror)
      if(ierror.ne.0) then
!         advienne que pourra sur le format
        write(nfecra,9550) rubriq
      endif
    enddo

!  Dans le cas specifique de la combustion de grains de charbon
!  avec un couplage retour sur une combustion gaz en phase porteuse

!      --> A verifier l'utilite de cette sauvegarde pour une suite...

    if (ippmod(icpl3c).ge.0) then
      do ii = 1, nsalpp
        icha = nsalto-nsalpp+ii
        itysup = 1
        nbval  = 1
        irtyp  = 2
        WRITE(CAR4,'(I4.4)') II
        RUBRIQ = 'scalaires_physiques_pariculieres_charbon'//CAR4
        call ecrsui(impvls,rubriq,len(rubriq),itysup,nbval,irtyp, &
                    propce(1,ipproc(icha)),ierror)
        if(ierror.ne.0) write(nfecra,9550) rubriq
      enddo

    endif

  endif

  write(nfecra,7013)

! ---> Fermeture du fichier suite
  call clssui(impvls,ierror)

  if(ierror.ne.0) then
    write(nfecra,9700) ficvls
  endif

! ---> En cas d'erreur, on continue quand meme
 9999   continue

  write(nfecra,7014)

endif

!===============================================================================
! 3. Visualisations
!===============================================================================

if (ntcabs.lt.ntmabs) return

nfin = 1

!-->Stockage des trajectoires au format Ensight Gold

if (iensi1.eq.1) then

  iforce = 0
  call enslag                                                     &
  !==========
   ( idebia , idebra ,                                            &
     nbpmax , nvp    , nvp1   , nvep   , nivep  ,                 &
     nfin   , iforce ,                                            &
     itepa  ,                                                     &
     ettp   , tepa   , ra )

endif

!-->Stockage des deplacements au format Ensight Gold

if (iensi2.eq.1) then

  ifinia = idebia
  itrav1 = idebra
  ifinra = itrav1 + 3*nbpmax
  call rasize ('lagout',ifinra)
  !==========

  call enswaf                                                     &
  !==========
   ( nbpmax , nvp    , nvp1   , nvep   , nivep  ,                 &
     nfin   ,                                                     &
     itepa  ,                                                     &
     ettp   , tepa , ra(itrav1)   )

endif

!===============================================================================
! 4. FIN
!===============================================================================

return

!===============================================================================

!--------
! FORMATS
!--------

 6010 FORMAT (/, 3X,'** INFORMATIONS SUR LE CALCUL LAGRANGIEN       ',/,&
           3X,'   -------------------------------------       ',/,&
           3X,' Ecriture d''un fichier suite                  ',/,&
           3X,'   sur les variables liees aux particules      ',/)


 6011 FORMAT (   3X,'   Debut de l''ecriture                        '  )
 6012 FORMAT (   3X,'   Fin de l''ecriture des infos sur le calcul  '  )
 6013 FORMAT (   3X,'   Fin de l''ecriture des infos particulaires  '  )
 6014 FORMAT (   3X,' Fin de l''ecriture du fichier suite           ',/,&
           3X,'   sur les variables liees aux particules      ',/)

 7010 FORMAT (/, 3X,'** INFORMATIONS SUR LE CALCUL LAGRANGIEN       ',/,&
           3X,'   -------------------------------------       ',/,&
           3X,' Ecriture d''un fichier suite                  ',/,&
           3X,'   sur les statistiques volumiques et aux      ',/,&
           3X,'   fontieres, ainsi que les termes sources     ',/,&
           3X,'   de couplage retour                          ',/)


 7011 FORMAT (   3X,'   Debut de l''ecriture des stats et TS        '  )
 7013 FORMAT (   3X,'   Fin de l''ecriture des statistiques et TS   '  )
 7014 FORMAT (   3X,' Fin de l''ecriture du fichier suite           ',/,&
           3X,'   sur les statistiques et TS couplage retour  ',/)

 9010 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION: A L''ECRITURE DU FICHIER SUITE LAGRANGIEN    ',/,&
'@    =========    SUR LES CARACTERISTIQUES DES PARTICULES    ',/,&
'@                                                            ',/,&
'@    ERREUR A L''OUVERTURE EN ECRITURE DU FICHIER SUITE      ',/,&
'@      (',A13,')                                             ',/,&
'@                                                            ',/,&
'@  Le calcul se termine mais ne fournira pas de fichier      ',/,&
'@    suite sur les caracteristiques des particules.          ',/,&
'@                                                            ',/,&
'@  Verifier que le repertoire de travail est accessible en   ',/,&
'@    ecriture et que le fichier suite peut y etre cree.      ',/,&
'@  Voir le sous-programme LAGOUT.                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 9020 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION: A L''ECRITURE DU FICHIER SUITE LAGRANGIEN    ',/,&
'@    =========    SUR LES CARACTERISTIQUES DES PARTICULES    ',/,&
'@                                                            ',/,&
'@      ERREUR A L''ECRITURE DE LA RUBRIQUE                   ',/,&
'@ ',A60                                                       ,/,&
'@                                                            ',/,&
'@      LA VALEUR DU MOT CLE CONCERNE VAUT :                  ',/,&
'@        ',A10   ,'  = ',I10                                  ,/,&
'@                                                            ',/,&
'@    Le calcul continue...                                   ',/,&
'@                                                            ',/,&
'@    Contacter l''equipe de developpement.                   ',/,&
'@    Voir le sous-programme LAGOUT.                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 9021 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION: A L''ECRITURE DU FICHIER SUITE LAGRANGIEN    ',/,&
'@    =========    SUR LES CARACTERISTIQUES DES PARTICULES    ',/,&
'@                                                            ',/,&
'@      ERREUR A L''ECRITURE DE LA RUBRIQUE                   ',/,&
'@ ',A60                                                       ,/,&
'@                                                            ',/,&
'@      LA VALEUR DU MOT CLE CONCERNE VAUT :                  ',/,&
'@        ',A10   ,'  = ',E14.5                                ,/,&
'@                                                            ',/,&
'@    Le calcul continue...                                   ',/,&
'@                                                            ',/,&
'@    Contacter l''equipe de developpement.                   ',/,&
'@    Voir le sous-programme LAGOUT.                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 9030 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION: A L''ECRITURE DU FICHIER SUITE LAGRANGIEN    ',/,&
'@    =========    SUR LES CARACTERISTIQUES DES PARTICULES    ',/,&
'@                                                            ',/,&
'@      ERREUR A L''ECRITURE DE LA RUBRIQUE                   ',/,&
'@ ',A60                                                       ,/,&
'@                                                            ',/,&
'@      LA VALEUR DU MOT CLE CONCERNE VAUT :                  ',/,&
'@        ',A10   ,'  = ',I10                                  ,/,&
'@                                                            ',/,&
'@    Le calcul continue mais ne fournira pas de fichier      ',/,&
'@      suite sur les caracteristiques des particules.        ',/,&
'@                                                            ',/,&
'@    Contacter l''equipe de developpement.                   ',/,&
'@    Voir le sous-programme LAGOUT.                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 9100 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION: A L''ECRITURE DU FICHIER SUITE LAGRANGIEN    ',/,&
'@    =========    SUR LES CARACTERISTIQUES DES PARTICULES    ',/,&
'@                                                            ',/,&
'@      ERREUR A L''ECRITURE DE LA RUBRIQUE                   ',/,&
'@ ',A60                                                       ,/,&
'@                                                            ',/,&
'@    Le calcul continue mais ne fournira pas de fichier      ',/,&
'@      suite sur les caracteristiques des particules.        ',/,&
'@                                                            ',/,&
'@    Contacter l''equipe de developpement.                   ',/,&
'@    Voir le sous-programme LAGOUT.                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 9101 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION: A L''ECRITURE DU FICHIER SUITE LAGRANGIEN    ',/,&
'@    =========    SUR LES CARACTERISTIQUES DES PARTICULES    ',/,&
'@                                                            ',/,&
'@      ERREUR A L''ECRITURE DE LA RUBRIQUE                   ',/,&
'@ ',A60                                                       ,/,&
'@                                                            ',/,&
'@    Le calcul continue...                                   ',/,&
'@                                                            ',/,&
'@    Contacter l''equipe de developpement.                   ',/,&
'@    Voir le sous-programme LAGOUT.                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 9140 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION: A LA FERMETURE DU FICHIER SUITE LAGRANGIEN   ',/,&
'@    =========    SUR LES CARACTERISTIQUES DES PARTICULES    ',/,&
'@                                                            ',/,&
'@    Probleme sur le fichier de nom : ',A13                   ,/,&
'@                                                            ',/,&
'@    Le calcul continue...                                   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 9510 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION: A L''ECRITURE DU FICHIER SUITE LAGRANGIEN    ',/,&
'@    =========    SUR LES STATISTIQUES ET LES TERMES SOURCES ',/,&
'@                 DE COUPLAGE RETOUR                         ',/,&
'@                                                            ',/,&
'@    ERREUR A L''OUVERTURE EN ECRITURE DU FICHIER SUITE      ',/,&
'@      (',A13,')                                             ',/,&
'@                                                            ',/,&
'@  Le calcul continue mais ne fournira pas de fichier        ',/,&
'@    suite sur les statistiques volumiques et aux frontieres ',/,&
'@    ainsi que sur les termes sources de couplage retour     ',/,&
'@                                                            ',/,&
'@  Verifier que le repertoire de travail est accessible en   ',/,&
'@    ecriture et que le fichier suite peut y etre cree.      ',/,&
'@  Voir le sous-programme LAGOUT.                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 9520 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION: A L''ECRITURE DU FICHIER SUITE LAGRANGIEN    ',/,&
'@    =========    SUR LES STATISTIQUES ET LES TERMES SOURCES ',/,&
'@                 DE COUPLAGE RETOUR                         ',/,&
'@                                                            ',/,&
'@      ERREUR A L''ECRITURE DE LA RUBRIQUE                   ',/,&
'@ ',A60                                                       ,/,&
'@                                                            ',/,&
'@      LA VALEUR DU MOT CLE CONCERNE VAUT :                  ',/,&
'@        ',A10   ,'  = ',I10                                  ,/,&
'@                                                            ',/,&
'@    Le calcul continue...                                   ',/,&
'@                                                            ',/,&
'@    Contacter l''equipe de developpement.                   ',/,&
'@    Voir le sous-programme LAGOUT.                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 9521 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION: A L''ECRITURE DU FICHIER SUITE LAGRANGIEN    ',/,&
'@    =========    SUR LES STATISTIQUES ET LES TERMES SOURCES ',/,&
'@                 DE COUPLAGE RETOUR                         ',/,&
'@                                                            ',/,&
'@      ERREUR A L''ECRITURE DE LA RUBRIQUE                   ',/,&
'@ ',A60                                                       ,/,&
'@                                                            ',/,&
'@      LA VALEUR DU MOT CLE CONCERNE VAUT :                  ',/,&
'@        ',A10   ,'  = ',E14.5                                ,/,&
'@                                                            ',/,&
'@    Le calcul continue...                                   ',/,&
'@                                                            ',/,&
'@    Contacter l''equipe de developpement.                   ',/,&
'@    Voir le sous-programme LAGOUT.                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 9550 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION: A L''ECRITURE DU FICHIER SUITE LAGRANGIEN    ',/,&
'@    =========    SUR LES STATISTIQUES ET LES TERMES SOURCES ',/,&
'@                 DE COUPLAGE RETOUR                         ',/,&
'@                                                            ',/,&
'@      ERREUR A L''ECRITURE DE LA RUBRIQUE                   ',/,&
'@ ',A60                                                       ,/,&
'@                                                            ',/,&
'@    Le calcul continue...                                   ',/,&
'@                                                            ',/,&
'@    Contacter l''equipe de developpement.                   ',/,&
'@    Voir le sous-programme LAGOUT.                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 9700 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION: A LA FERMETURE DU FICHIER SUITE LAGRANGIEN   ',/,&
'@    =========    SUR LES STATISTIQUES ET LES TERMES SOURCES ',/,&
'@                 DE COUPLAGE RETOUR                         ',/,&
'@                                                            ',/,&
'@    Probleme sur le fichier de nom : ',A13                   ,/,&
'@                                                            ',/,&
'@    Le calcul continue...                                   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

!----
! FIN
!----

end subroutine
