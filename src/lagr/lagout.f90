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

subroutine lagout &
!================

 ( nvep   , nivep  ,                                              &
   ntersl , nvlsta , nvisbr ,                                     &
   propce )

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
! nvep             ! i  ! <-- ! nombre info particulaires (reels)              !
! nivep            ! i  ! <-- ! nombre info particulaires (entiers)            !
! ntersl           ! i  ! <-- ! nbr termes sources de couplage retour          !
! nvlsta           ! i  ! <-- ! nombre de var statistiques lagrangien          !
! nvisbr           ! i  ! <-- ! nombre de statistiques aux frontieres          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use, intrinsic :: iso_c_binding

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
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer          nvep  , nivep
integer          ntersl , nvlsta , nvisbr

double precision propce(ncelet,*)

! Local variables

character        rubriq*64 , car4*4
character        nomnvl(nvplmx)*60 , nomtsl(nvplmx)*60
character        nomite(nvplmx)*64 , nomrte(nvplmx)*64
character        ficsui*32
integer          nbval, itysup , irfsup, idbase
integer          ivers  , ilecec, n_sec
integer          icha   , ii , ilayer
integer          ipas   , jj
integer          inmcoo
integer          ival(1)
double precision rval(1)

type(c_ptr) :: rp

integer, allocatable, dimension(:) :: icepar
double precision, allocatable, dimension(:,:) :: coopar

!===============================================================================

!===============================================================================
! Output restart file: variables related to particles
!===============================================================================

! Open restart file

write(nfecra,6010)

ficsui = 'lagrangian'
call restart_create(ficsui, '', 1, rp)

write(nfecra,6011)

! Entete et Infos sur le calcul ou on saute si erreur
!     On inclut une rubrique destinee a distinguer ce fichier
!       d'un autre fichier suite
!     Pour le moment, IVERS n'est pas utilise

itysup = 0
nbval  = 1

ivers  = 32000
rubriq = 'version_fichier_suite_Lagrangien_variables'
ival(1) = ivers
call restart_write_section_int_t(rp,rubriq,itysup,nbval,ival)

! Temps (par securite)

rubriq = 'nombre_iterations_Lagrangiennes'
ival(1) = iplas
call restart_write_section_int_t(rp,rubriq,itysup,nbval,ival)

rubriq = 'temps_physique_Lagrangien'
rval(1) = ttclag
call restart_write_section_real_t(rp,rubriq,itysup,nbval,rval)

! Infos sur le suivi du calcul

ival(1) = nbptot
rubriq = 'nombre_total_particules'
call restart_write_section_int_t(rp,rubriq,itysup,nbval,ival)

ival(1) = nbpert
rubriq = 'nombre_particules_perdues'
call restart_write_section_int_t(rp,rubriq,itysup,nbval,ival)

ival(1) = iphyla
rubriq = 'indicateur_physique_particules'
call restart_write_section_int_t(rp,rubriq,itysup,nbval,ival)

ival(1) = itpvar
rubriq = 'indicateur_temperature_particules'
call restart_write_section_int_t(rp,rubriq,itysup,nbval,ival)

ival(1) = idpvar
rubriq = 'indicateur_diametre_particules'
call restart_write_section_int_t(rp,rubriq,itysup,nbval,ival)

ival(1) = impvar
rubriq = 'indicateur_masse_particules'
call restart_write_section_int_t(rp,rubriq,itysup,nbval,ival)

ival(1) = nvls
rubriq = 'nombre_variables_utilisateur'
call restart_write_section_int_t(rp,rubriq,itysup,nbval,ival)

write(nfecra,6012)

n_sec = lagr_restart_write_particle_data(rp)

write(nfecra,6013)

! ---> Fermeture du fichier suite
call restart_destroy(rp)

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
  call restart_create(ficsui, '', 1, rp)

  write(nfecra,7011)

! Entete et Infos sur le calcul ou on saute si erreur
!     On inclut une rubrique destinee a distinguer ce fichier
!       d'un autre fichier suite
!     Pour le moment, IVERS n'est pas utilise

  itysup = 0
  nbval  = 1

  ivers  = 111
  rubriq = 'version_fichier_suite_Lagrangien_statistiques'
  ival(1) = ivers
  call restart_write_section_int_t(rp,rubriq,itysup,nbval,ival)

! ---> On ecrit ISTTIO c'est utile dans tous les cas

  rubriq = 'indicateur_ecoulement_stationnaire'
  ival(1) = isttio
  call restart_write_section_int_t(rp,rubriq,itysup,nbval,ival)

! --> En premier, on ecrit les statistiques volumiques

  if (istala.eq.1 .and. iplas.ge.idstnt) then

    rubriq = 'iteration_debut_statistiques'
    ival(1) = idstnt
    call restart_write_section_int_t(rp,rubriq,itysup,nbval,ival)

    rubriq = 'iteration_debut_statistiques_stationnaires'
    ival(1) = nstist
    call restart_write_section_int_t(rp,rubriq,itysup,nbval,ival)

    rubriq = 'nombre_iterations_statistiques_stationnaires'
    ival(1) = npst
    call restart_write_section_int_t(rp,rubriq,itysup,nbval,ival)

    rubriq = 'temps_statistiques_stationnaires'
    rval(1) = tstat
    call restart_write_section_real_t(rp,rubriq,itysup,nbval,rval)

    rubriq = 'classe_statistique_particules'
    ival(1) = nbclst
    call restart_write_section_int_t(rp,rubriq,itysup,nbval,ival)

    rubriq = 'nombre_statistiques_utilisateur'
    ival(1) = nvlsts
    call restart_write_section_int_t(rp,rubriq,itysup,nbval,ival)

    !  Statistiques volumiques

    itysup = 1
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
        call restart_write_section_real_t(rp,rubriq,itysup,nbval,statis(:,ii))

      enddo

      do jj = 1,nvlsta-1

        ii = ipas*nvlsta +jj
        if (ipas.gt.0) then
          write(car4,'(i4.4)') ipas
          rubriq = 'var_stat_vol_groupe_'//car4//'_'//nomlav(ii)
        else
          rubriq = 'var_stat_vol_'//nomlav(ii)
        endif
        call restart_write_section_real_t(rp,rubriq,itysup,nbval,stativ(:,ii))

      enddo

    enddo

  endif

! --> En second, c'est le tour des statistiques aux frontieres

  if (iensi3.eq.1 .and. nvisbr.gt.0) then

    itysup = 0
    nbval  = 1

    rubriq = 'iteration_debut_stats_frontieres_stationnaires'
    ival(1) = nstbor
    call restart_write_section_int_t(rp,rubriq,itysup,nbval,ival)

    rubriq = 'nombre_iterations_stats_frontieres'
    ival(1) = npstft
    call restart_write_section_int_t(rp,rubriq,itysup,nbval,ival)

    rubriq = 'nombre_iterations_stats_frontieres_stationnaires'
    ival(1) = npstf
    call restart_write_section_int_t(rp,rubriq,itysup,nbval,ival)

    rubriq = 'temps_stats_frontieres_stationnaires'
    rval(1) = tstatp
    call restart_write_section_real_t(rp,rubriq,itysup,nbval,rval)

    rubriq = 'nombre_stats_frontieres_utilisateur'
    ival(1) = nusbor
    call restart_write_section_int_t(rp,rubriq,itysup,nbval,ival)

    !  Statistiques aux frontieres

    itysup = 3
    nbval  = 1

    do ii = 1,nvisbr
      rubriq = 'stat_bord_'//nombrd(II)
      call restart_write_section_real_t(rp,rubriq,itysup,nbval,parbor(:,ii))
    enddo

  endif

  ! --> Enfin, en cas de couplage retour, on ecrit les termes sources

  if (iilagr.eq.2) then

    itysup = 0
    nbval  = 1

    rubriq = 'iteration_debut_termes_sources_stationnaires'
    ival(1) = nstits
    call restart_write_section_int_t(rp,rubriq,itysup,nbval,ival)

    rubriq = 'nombre_iterations_termes_sources_stationnaires'
    ival(1) = npts
    call restart_write_section_int_t(rp,rubriq,itysup,nbval,ival)

    rubriq = 'modele_turbulence_termes_sources'
    ival(1) = iturb
    call restart_write_section_int_t(rp,rubriq,itysup,nbval,ival)

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

    do ii = 1,ntersl
      rubriq = nomtsl(ii)
      call restart_write_section_real_t(rp,rubriq,itysup,nbval,tslagr(:,ii))
    enddo

    ! Dans le cas specifique de la combustion de grains de charbon
    ! avec un couplage retour sur une combustion gaz en phase porteuse

    ! --> A verifier l'utilite de cette sauvegarde pour une suite...

    if (ippmod(icpl3c).ge.0) then
      do ii = 1, nsalpp
        icha = nsalto-nsalpp+ii
        itysup = 1
        nbval  = 1
        write(car4,'(i4.4)') II
        rubriq = 'scalaires_physiques_pariculieres_charbon'//car4
        call restart_write_section_real_t(rp,rubriq,itysup,nbval,  &
                                          propce(:,ipproc(icha)))
      enddo

    endif

  endif

  write(nfecra,7013)

  ! close restart file
  call restart_destroy(rp)

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
