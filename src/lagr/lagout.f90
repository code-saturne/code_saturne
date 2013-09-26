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

subroutine lagout &
!================

 ( lndnod ,                                                       &
   nvar   , nscal  ,                                              &
   nbpmax , nvp    , nvp1   , nvep   , nivep  ,                   &
   ntersl , nvlsta , nvisbr ,                                     &
   icocel , itycel , itepa  ,                                     &
   dt     , rtpa   , rtp    , propce , propfb ,                   &
   coefa  , coefb  ,                                              &
   ettp   , tepa   , parbor , statis , stativ , tslagr )

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
! lndnod           ! i  ! <-- ! dim. connectivite cellules->faces              !
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! nbpmax           ! i  ! <-- ! nombre max de particulies autorise             !
! nvp              ! i  ! <-- ! nombre de variables particulaires              !
! nvp1             ! i  ! <-- ! nvp sans position, vfluide, vpart              !
! nvep             ! i  ! <-- ! nombre info particulaires (reels)              !
! nivep            ! i  ! <-- ! nombre info particulaires (entiers)            !
! ntersl           ! i  ! <-- ! nbr termes sources de couplage retour          !
! nvlsta           ! i  ! <-- ! nombre de var statistiques lagrangien          !
! nvisbr           ! i  ! <-- ! nombre de statistiques aux frontieres          !
! icocel           ! ia ! <-- ! connectivite cellules -> faces                 !
! (lndnod)         !    !     !    face de bord si numero negatif              !
! itycel           ! ia ! <-- ! connectivite cellules -> faces                 !
! (ncelet+1)       !    !     !    pointeur du tableau icocel                  !
! itepa            ! ia ! <-- ! info particulaires (entiers)                   !
! (nbpmax,nivep    !    !     !   (cellule de la particule,...)                !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
! ettp             ! ra ! <-- ! tableaux des variables liees                   !
!  (nbpmax,nvp)    !    !     !   aux particules etape courante                !
! tepa             ! ra ! <-- ! info particulaires (reels)                     !
! (nbpmax,nvep)    !    !     !   (poids statistiques,...)                     !
! statis           ! ra ! <-- ! cumul pour les moyennes des                    !
!(ncelet,nvlsta    !    !     !   statistiques volumiques                      !
! stativ           ! ra ! <-- ! cumul pour les variances des                   !
!(ncelet,          !    !     !    statistiques volumiques                     !
!   nvlsta-1)      !    !     !                                                !
! tslagr           ! ra ! <-- ! terme de couplage retour du                    !
!(ncelet,ntersl    !    !     !   lagrangien sur la phase porteuse             !
! parbor           ! ra ! <-- ! infos sur interaction des particules           !
!(nfabor,nvisbr    !    !     !   aux faces de bord                            !
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
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

integer          lndnod
integer          nvar   , nscal
integer          nbpmax , nvp    , nvp1   , nvep  , nivep
integer          ntersl , nvlsta , nvisbr

integer          icocel(lndnod) , itycel(ncelet+1)
integer          itepa(nbpmax,nivep)

double precision dt(ncelet) , rtp(ncelet,*) , rtpa(ncelet,*)
double precision propce(ncelet,*), propfb(nfabor,*)
double precision coefa(nfabor,*) , coefb(nfabor,*)
double precision ettp(nbpmax,nvp) ,  tepa(nbpmax,nvep)
double precision statis(ncelet,nvlsta)
double precision stativ(ncelet,nvlsta-1)
double precision tslagr(ncelet,ntersl)
double precision parbor(nfabor,nvisbr)

! Local variables


character        rubriq*64 , car4*4
character        nomnvl(nvplmx)*60 , nomtsl(nvplmx)*60
character        nomite(nvplmx)*64 , nomrte(nvplmx)*64
character        ficsui*32
integer          ierror , nbval, itysup , irtyp  , irfsup, idbase
integer          ivers  , ilecec
integer          icha   , ii , ilayer
integer          itrav1
integer          ipas   , jj
integer          impavl , impvls, inmcoo, ipasup

integer, allocatable, dimension(:) :: icepar
double precision, allocatable, dimension(:,:) :: coopar

!===============================================================================

!===============================================================================
! Output restart file: variables related to particles
!===============================================================================

! Open restart file

write(nfecra,6010)

ilecec = 2
ficsui = 'lagrangian'
call opnsui(ficsui, len(ficsui), ilecec, impavl, ierror)
!==========
if (ierror.ne.0) then
  goto 9998
endif

write(nfecra,6011)

! Entete et Infos sur le calcul ou on saute si erreur
!     On inclut une rubrique destinee a distinguer ce fichier
!       d'un autre fichier suite
!     Pour le moment, IVERS n'est pas utilise

itysup = 0
nbval  = 1

ivers  = 112
rubriq = 'version_fichier_suite_Lagrangien_variables'
irtyp  = 1
call ecrsui(impavl,rubriq,len(rubriq),itysup,nbval,irtyp,ivers)

! Temps (par securite)

rubriq = 'nombre_iterations_Lagrangiennes'
irtyp  = 1
call ecrsui(impavl,rubriq,len(rubriq),itysup,nbval,irtyp,iplas)

rubriq = 'temps_physique_Lagrangien'
irtyp  = 2
call ecrsui(impavl,rubriq,len(rubriq),itysup,nbval,irtyp,ttclag)

! Infos sur le suivi du calcul

irtyp  = 1
inmcoo = 0

allocate(icepar(nbpart))
allocate(coopar(3,nbpart))

do ii = 1, nbpart
  icepar(ii) = abs(itepa(ii,jisor))
  coopar(1,ii) = ettp(ii,jxp)
  coopar(2,ii) = ettp(ii,jyp)
  coopar(3,ii) = ettp(ii,jzp)
enddo

rubriq = 'particles'
call ecpsui(impavl,rubriq,len(rubriq),inmcoo,nbpart,icepar,coopar,ipasup)

deallocate(coopar)

rubriq = 'nombre_total_particules'
call ecrsui(impavl,rubriq,len(rubriq),itysup,nbval,irtyp,nbptot)

rubriq = 'nombre_particules_perdues'
call ecrsui(impavl,rubriq,len(rubriq),itysup,nbval,irtyp,nbpert)

rubriq = 'indicateur_physique_particules'
call ecrsui(impavl,rubriq,len(rubriq),itysup,nbval,irtyp,iphyla)

rubriq = 'indicateur_temperature_particules'
call ecrsui(impavl,rubriq,len(rubriq),itysup,nbval,irtyp,itpvar)

rubriq = 'indicateur_diametre_particules'
call ecrsui(impavl,rubriq,len(rubriq),itysup,nbval,irtyp,idpvar)

rubriq = 'indicateur_masse_particules'
call ecrsui(impavl,rubriq,len(rubriq),itysup,nbval,irtyp,impvar)

rubriq = 'nombre_variables_utilisateur'
call ecrsui(impavl,rubriq,len(rubriq),itysup,nbval,irtyp,nvls)

write(nfecra,6012)

! Particle flags (currently: stuck or not)

do ii = 1, nbpart
  if (itepa(ii,jisor) .lt. 0) then
    icepar(ii) = 1
  else
    icepar(ii) = 0
  endif
enddo

itysup = ipasup
nbval  = 1
irtyp  = 1

rubriq = 'particle_status_flag'
call ecrsui(impavl, rubriq, len(rubriq), itysup, nbval, irtyp, icepar)

deallocate(icepar)

! Variables particulaires

nomnvl(jup) = 'variable_vitesseU_particule'
nomnvl(jvp) = 'variable_vitesseV_particule'
nomnvl(jwp) = 'variable_vitesseW_particule'
nomnvl(juf) = 'variable_vitesseU_fluide_vu'
nomnvl(jvf) = 'variable_vitesseV_fluide_vu'
nomnvl(jwf) = 'variable_vitesseW_fluide_vu'
nomnvl(jmp) = 'variable_masse_particule'
nomnvl(jdp) = 'variable_diametre_particule'
if (iphyla.eq.1 .and. itpvar.eq.1) then
  nomnvl(jtp) = 'variable_temperature_particule'
  nomnvl(jtf) = 'variable_temperature_fluide_vu'
  nomnvl(jcp) = 'variable_chaleur_specifique_particule'
elseif (iphyla.eq.2) then
  do ilayer = 1, nlayer
    write(nomnvl(jhp(ilayer)),'(A38,I4.4)') 'variable_temperature_particule_couche_',ilayer
  enddo
  nomnvl(jtf) = 'variable_temperature_fluide_vu'
  nomnvl(jmwat) = 'variable_masse_humidite'
  do ilayer = 1, nlayer
    write(nomnvl(jmch(ilayer)),'(A38,I4.4)') 'variable_masse_charbon_reactif_couche_',ilayer
  enddo
  do ilayer = 1, nlayer
    write(nomnvl(jmck(ilayer)),'(A27,I4.4)') 'variable_masse_coke_couche_',ilayer
  enddo
  nomnvl(jcp) = 'variable_chaleur_specifique_particule'
endif
if (nvls.gt.0) then
  do ii = 1,nvls
    write(car4,'(i4.4)') ii
    nomnvl(jvls(ii)) = 'variable_supplementaire_'//car4
  enddo
endif

itysup = ipasup
nbval  = 1
irtyp  = 2

do ii = jmp,jwf
  if (ii .lt. jxp .or. ii.gt.jzp) then
    rubriq = nomnvl(ii)
    call ecrsui(impavl,rubriq,len(rubriq),itysup,nbval,irtyp,ettp(1,ii))
  endif
enddo

do ii = 1,jmp-1
  rubriq = nomnvl(ii)
  call ecrsui(impavl,rubriq,len(rubriq),itysup,nbval,irtyp,ettp(1,ii))
enddo

! Caracteristiques et infos particulaires (ENTIERS)

nomite(jisor) = 'indicateur_'

if (iphyla.eq.2) then
  nomite(jinch) = 'numero_charbon'
endif

! Deposition submodel
if (idepst.eq.1) then
  nomite(jimark) = 'indicateur_de_saut'
  nomite(jdiel) = 'diel_particules'
  nomite(jdfac) = 'dfac_particules'
  nomite(jdifel) = 'difel_particules'
  nomite(jtraj) = 'traj_particules'
  nomite(jptdet) = 'ptdet_particules'
  nomite(jinjst) = 'indic_stat'
  nomite(jdepo) = 'part_depo'
endif

if (ireent.eq.1) then
   nomite(jnbasg) = 'nb_ls_aspe'
   nomite(jnbasp) = 'nb_sms_aspe'
endif

itysup = ipasup
nbval  = 1
irtyp  = 1

do ii = 1, nivep
  if (ii .ne. jisor) then
    rubriq = nomite(ii)
    if (ii.eq.jdfac) then
      idbase = 1
      irfsup = 3
      call ecisui(impavl, rubriq, len(rubriq), itysup, irfsup, idbase,  &
                  itepa(1,ii))
    else
      call ecrsui(impavl, rubriq, len(rubriq), itysup, nbval, irtyp,  &
                  itepa(1,ii))
    endif
  endif
enddo

! groupe statistique particules

if (nbclst .gt. 0 ) then
  nomite(jclst) = 'numero_groupe_statistiques'

  itysup = ipasup
  nbval  = 1
  irtyp  = 1

  rubriq = nomite(jclst)
  call ecrsui(impavl,rubriq,len(rubriq),itysup,nbval,irtyp,itepa(1,jclst))
endif

! Numero du charbon des particules

if (iphyla.eq.2) then
  rubriq = nomite(jinch)
  call ecrsui(impavl,rubriq,len(rubriq),itysup,nbval,irtyp,itepa(1,jinch))
endif

! Caracteristiques et infos particulaires (REELS)

nomrte(jrval) = 'random_value'
nomrte(jrtsp) = 'temps_sejour_particules'
nomrte(jrpoi) = 'poids_statistiques_particules'
if (iphyla.eq.1 .and. itpvar.eq.1 .and.iirayo.gt.0) then
  nomrte(jreps) = 'emissivite_particules'
endif
if (iphyla.eq.2) then
  nomrte(jrdck) = 'diametre_coeur_retrecissant_charbon'
  nomrte(jrd0p) = 'diametre_initial_charbon'
  do ilayer = 1, nlayer
    write(nomrte(jrhock(ilayer)),'(A28,I4.4)') 'masse_volumique_coke_couche_',ilayer
  enddo
endif

! Deposition submodel
if (idepst.eq.1) then
  nomrte(jryplu) = 'yplus_particules'
  nomrte(jrinpf) = 'dx_particules'
endif

if (ireent.eq.1) then
   nomrte(jfadh) = 'force_adhesion'
   nomrte(jmfadh) = 'moment_adhesion'
   nomrte(jndisp) = 'disp_norm'
endif

itysup = ipasup
nbval  = 1
irtyp  = 2

do ii = 1, nvep
  rubriq = nomrte(ii)
  call ecrsui(impavl,rubriq,len(rubriq),itysup,nbval,irtyp,tepa(1,ii))
enddo

write(nfecra,6013)

! ---> Fermeture du fichier suite
call clssui(impavl,ierror)

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
  rubriq = 'version_fichier_suite_Lagrangien_statistiques'
  irtyp  = 1
  call ecrsui(impvls,rubriq,len(rubriq),itysup,nbval,irtyp,ivers)

! ---> On ecrit ISTTIO c'est utile dans tous les cas

  rubriq = 'indicateur_ecoulement_stationnaire'
  irtyp  = 1
  call ecrsui(impvls,rubriq,len(rubriq),itysup,nbval,irtyp,isttio)

! --> En premier, on ecrit les statistiques volumiques

  if (istala.eq.1 .and. iplas.ge.idstnt) then

    rubriq = 'iteration_debut_statistiques'
    irtyp  = 1
    call ecrsui(impvls,rubriq,len(rubriq),itysup,nbval,irtyp,idstnt)

    rubriq = 'iteration_debut_statistiques_stationnaires'
    irtyp  = 1
    call ecrsui(impvls,rubriq,len(rubriq),itysup,nbval,irtyp,nstist)

    rubriq = 'nombre_iterations_statistiques_stationnaires'
    irtyp  = 1
    call ecrsui(impvls,rubriq,len(rubriq),itysup,nbval,irtyp,npst)

    rubriq = 'temps_statistiques_stationnaires'
    irtyp  = 2
    call ecrsui(impvls,rubriq,len(rubriq),itysup,nbval,irtyp,tstat)

    rubriq = 'classe_statistique_particules'
    irtyp  = 1
    call ecrsui(impvls,rubriq,len(rubriq),itysup,nbval,irtyp,nbclst)

    rubriq = 'nombre_statistiques_utilisateur'
    irtyp  = 1
    call ecrsui(impvls,rubriq,len(rubriq),itysup,nbval,irtyp,nvlsts)

    !  Statistiques volumiques

    itysup = 1
    irtyp  = 2
    nbval  = 1

    do ipas  = 0,nbclst
      do jj = 1,nvlsta

        ii = ipas*nvlsta +jj
        if (ipas.gt.0) then
          write(car4,'(i4.4)') ipas
          rubriq = 'moy_stat_vol_groupe_'//car4//'_'//nomlag(ii)
        else
          rubriq = 'moy_stat_vol_'//nomlag(ii)
        endif
        call ecrsui(impvls,rubriq,len(rubriq),itysup,nbval,irtyp,statis(1,ii))

      enddo

      do jj = 1,nvlsta-1

        ii = ipas*nvlsta +jj
        if (ipas.gt.0) then
          write(car4,'(i4.4)') ipas
          rubriq = 'var_stat_vol_groupe_'//car4//'_'//nomlav(ii)
        else
          rubriq = 'var_stat_vol_'//nomlav(ii)
        endif
        call ecrsui(impvls,rubriq,len(rubriq),itysup,nbval,irtyp, &
                    stativ(1,ii))

      enddo

    enddo

  endif

! --> En second, c'est le tour des statistiques aux frontieres

  if (iensi3.eq.1 .and. nvisbr.gt.0) then

    itysup = 0
    nbval  = 1

    rubriq = 'iteration_debut_stats_frontieres_stationnaires'
    irtyp  = 1
    call ecrsui(impvls,rubriq,len(rubriq),itysup,nbval,irtyp,nstbor)

    rubriq = 'nombre_iterations_stats_frontieres'
    irtyp  = 1
    call ecrsui(impvls,rubriq,len(rubriq),itysup,nbval,irtyp,npstft)

    rubriq = 'nombre_iterations_stats_frontieres_stationnaires'
    irtyp  = 1
    call ecrsui(impvls,rubriq,len(rubriq),itysup,nbval,irtyp,npstf)

    rubriq = 'temps_stats_frontieres_stationnaires'
    irtyp  = 2
    call ecrsui(impvls,rubriq,len(rubriq),itysup,nbval,irtyp,tstatp)

    rubriq = 'nombre_stats_frontieres_utilisateur'
    irtyp  = 1
    call ecrsui(impvls,rubriq,len(rubriq),itysup,nbval,irtyp,nusbor)

    !  Statistiques aux frontieres

    itysup = 3
    nbval  = 1
    irtyp  = 2

    do ii = 1,nvisbr
      rubriq = 'stat_bord_'//nombrd(II)
      call ecrsui(impvls,rubriq,len(rubriq),itysup,nbval,irtyp,parbor(1,ii))
    enddo

  endif

  ! --> Enfin, en cas de couplage retour, on ecrit les termes sources

  if (iilagr.eq.2) then

    itysup = 0
    nbval  = 1

    rubriq = 'iteration_debut_termes_sources_stationnaires'
    irtyp  = 1
    call ecrsui(impvls,rubriq,len(rubriq),itysup,nbval,irtyp,nstits)

    rubriq = 'nombre_iterations_termes_sources_stationnaires'
    irtyp  = 1
    call ecrsui(impvls,rubriq,len(rubriq),itysup,nbval,irtyp,npts)

    rubriq = 'modele_turbulence_termes_sources'
    irtyp  = 1
    call ecrsui(impvls,rubriq,len(rubriq),itysup,nbval,irtyp,iturb)

    ! On donne des labels au different TS pour les noms de rubriques
    ! On donne le meme label au keps, au v2f et au k-omega (meme variable k)

    if (ltsdyn.eq.1) then
      nomtsl(itsvx) = 'terme_source_vitesseX'
      nomtsl(itsvy) = 'terme_source_vitesseY'
      nomtsl(itsvz) = 'terme_source_vitesseZ'
      nomtsl(itsli) = 'terme_source_vitesse_implicite'
      if (itytur.eq.2 .or. iturb.eq.50 .or. iturb.eq.60) then
        nomtsl(itske) = 'terme_source_turbulence_keps'
      else if (itytur.eq.3) then
        nomtsl(itsr11) = 'terme_source_turbulence_R11'
        nomtsl(itsr12) = 'terme_source_turbulence_R12'
        nomtsl(itsr13) = 'terme_source_turbulence_R13'
        nomtsl(itsr22) = 'terme_source_turbulence_R22'
        nomtsl(itsr23) = 'terme_source_turbulence_R23'
        nomtsl(itsr33) = 'terme_source_turbulence_R33'
      endif
    endif
    if (ltsmas.eq.1) then
      nomtsl(itsmas) = 'terme_source_masse'
    endif
    if (ltsthe.eq.1) then
      if (iphyla.eq.1 .and. itpvar.eq.1) then
        nomtsl(itste) = 'terme_source_thermique_explicite'
        nomtsl(itsti) = 'terme_source_thermique_implicite'
      else if (iphyla.eq.2) then
        nomtsl(itste) = 'terme_source_thermique_explicite'
        nomtsl(itsti) = 'terme_source_thermique_implicite'
        do icha = 1,ncharb
          write(car4,'(i4.4)') icha
          nomtsl(itsmv1(icha)) = 'terme_source_legeres_F1_'//car4
          nomtsl(itsmv2(icha)) = 'terme_source_lourdes_F2_'//car4
        enddo
        nomtsl(itsco) = 'terme_source_F3'
        nomtsl(itsfp4) = 'terme_source_variance_traceur_air'
      endif
    endif

    ! Termes source de couplage retour

    itysup = 1
    nbval  = 1
    irtyp  = 2

    do ii = 1,ntersl
      rubriq = nomtsl(ii)
      call ecrsui(impvls,rubriq,len(rubriq),itysup,nbval,irtyp,tslagr(1,ii))
    enddo

    ! Dans le cas specifique de la combustion de grains de charbon
    ! avec un couplage retour sur une combustion gaz en phase porteuse

    ! --> A verifier l'utilite de cette sauvegarde pour une suite...

    if (ippmod(icpl3c).ge.0) then
      do ii = 1, nsalpp
        icha = nsalto-nsalpp+ii
        itysup = 1
        nbval  = 1
        irtyp  = 2
        write(car4,'(i4.4)') II
        rubriq = 'scalaires_physiques_pariculieres_charbon'//car4
        call ecrsui(impvls,rubriq,len(rubriq),itysup,nbval,irtyp, &
                    propce(1,ipproc(icha)))
      enddo

    endif

  endif

  write(nfecra,7013)

  ! close restart file
  call clssui(impvls,ierror)

  ! In case of error not leading to an abort in lower-level layers, continue
9999 continue

  write(nfecra,7014)

endif

!===============================================================================
! End
!===============================================================================

return

!===============================================================================

!--------
! Formats
!--------

#if defined(_CS_LANG_FR)

 6010 format(3x,'** Ecriture du fichier suite lagrangien',/,      &
             3x,'   ------------------------------------',/)

 6011 format(3x,'   Debut de l''ecriture')
 6012 format(3x,'   Fin de l''ecriture des infos sur le calcul')
 6013 format(3x,'   Fin de l''ecriture des infos particulaires')
 6014 format(3x,' Fin de l''ecriture du fichier suite',         /,&
             3x,'   sur les variables liees aux particules',/)

 7010 format(/, 3x,'** INFORMATIONS SUR LE CALCUL LAGRANGIEN',  /,&
                3x,'   -------------------------------------',  /,&
                3x,' Ecriture d''un fichier suite',             /,&
                3x,'   sur les statistiques volumiques et aux', /,&
                3x,'   fontieres, ainsi que les termes sources',/,&
                3x,'   de couplage retour',/)


 7011 format(3x,'   Debut de l''ecriture des stats et TS')
 7013 format(3x,'   Fin de l''ecriture des statistiques et TS')
 7014 format(3x,' Fin de l''ecriture du fichier suite',         /,&
             3x,'   sur les statistiques et TS couplage retour',/)

#else

 6010 format(3x,'** Writing the Lagrangian restart file',/,       &
             3x,'   -----------------------------------',/)

 6011 format(3x,'   Start writing')
 6012 format(3x,'   End writing info on calculation')
 6013 format(3x,'   End writing of specific info')
 6014 format(3x,' End writing of restart file',                 /,&
             3x,'   on particle-based variables',/)

 7010 format(/, 3x,'** INFORMATION ON LAGRANGIAN CALCULATION',  /,&
                3x,'   -------------------------------------',  /,&
                3x,' Writing a restart file',                   /,&
                3x,'   for volume and boundary statistics',     /,&
                3x,'   as well as for return coupling',         /,&
                3x,'   source terms',/)

 7011 format(3x,'   Start writing statistics and ST')
 7013 format(3x,'   End writign statistics and ST')
 7014 format(3x,' End writing of restart file',                 /,&
             3x,'   on statistics and return coupling ST',/)

#endif

!----
! End
!----

end subroutine
