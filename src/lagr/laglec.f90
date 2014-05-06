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

subroutine laglec &
!================

 ( ncelet , ncel   , nfabor ,                                     &
   nbpmax , nvp    , nvep   , nivep  ,                            &
   ntersl , nvlsta , nvisbr ,                                     &
   itepa  ,                                                       &
   rtpa   , propce ,                                              &
   ettp   , tepa   , statis , stativ , parbor , tslagr )

!===============================================================================
! FONCTION :
! ----------

!   SOUS-PROGRAMME DU MODULE LAGRANGIEN :
!   -------------------------------------

!    Lecture des fichiers suite Lagrangien "lagamo" et "lasamo"
!    contenant les informations sur les particule, les statistiques
!    volumiques et aux frontieres, ainsi que les termes sources
!    de couplage retour.

!    Tous les tableaux sont initialise a zero avant d'être remplis
!    dans le cas d'une suite (sinon ils restent a zero).
!    On realise donc ici l'initialisation des tableaux ouverts
!    dans MEMLA1, ce qui termine l'etape d'initialisation debutee
!    dans LAGOPT.

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! ncelet           ! i  ! <-- ! number of extended (real + ghost) cells        !
! ncel             ! i  ! <-- ! number of cells                                !
! nfabor           ! i  ! <-- ! number of boundary faces                       !
! nbpmax           ! e  ! <-- ! nombre max de particulies autorise             !
! nvp              ! e  ! <-- ! nombre de variables particulaires              !
! nvep             ! e  ! <-- ! nombre info particulaires (reels)              !
! nivep            ! e  ! <-- ! nombre info particulaires (entiers)            !
! ntersl           ! e  ! <-- ! nbr termes sources de couplage retour          !
! nvlsta           ! e  ! <-- ! nombre de var statistiques lagrangien          !
! nvisbr           ! e  ! <-- ! nombre de statistiques aux frontieres          !
! itepa            ! te ! <-- ! info particulaires (entiers)                   !
! (nbpmax,nivep    !    !     !   (cellule de la particule,...)                !
! rtpa             ! tr ! <-- ! variables de calcul au centre des              !
! (ncelet,*)       !    !     !    cellules instant precedent                  !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! ettp             ! tr ! <-- ! tableaux des variables liees                   !
!  (nbpmax,nvp)    !    !     !   aux particules etape courante                !
! tepa             ! tr ! <-- ! info particulaires (reels)                     !
! (nbpmax,nvep)    !    !     !   (poids statistiques,...)                     !
! statis           ! tr ! <-- ! cumul pour les moyennes des                    !
!(ncelet,nvlsta    !    !     !   statistiques volumiques                      !
! stativ           ! tr ! <-- ! cumul pour les variances des                   !
!(ncelet,          !    !     !    statistiques volumiques                     !
!   nvlsta-1)      !    !     !                                                !
! parbor           ! tr ! <-- ! infos sur interaction des particules           !
!(nfabor,nvisbr    !    !     !   aux faces de bord                            !
! tslagr           ! tr ! <-- ! terme de couplage retour du                    !
!(ncelet,ntersl    !    !     !   lagrangien sur la phase porteuse             !
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
use cstphy
use numvar
use optcal
use dimens, only: nvar
use entsor
use parall
use period
use lagpar
use lagran
use ppppar
use ppthch
use ppincl
use cpincl
use radiat

!===============================================================================

implicit none

! Arguments

integer          ncelet , ncel   , nfabor
integer          nbpmax , nvp    , nvep  , nivep
integer          ntersl , nvlsta , nvisbr
integer          itepa(nbpmax,nivep)

double precision rtpa(ncelet,nflown:nvar) , propce(ncelet,*)
double precision ettp(nbpmax,nvp) , tepa(nbpmax,nvep)
double precision statis(ncelet,nvlsta)
double precision stativ(ncelet,nvlsta-1)
double precision tslagr(ncelet,ntersl)
double precision parbor(nfabor,nvisbr)

! Local variables

character        rubriq*64 , car4*4, car8*8, kar8*8
character        nomnvl(nvplmx)*60 , nomtsl(nvplmx)*60
character        nomite(nvplmx)*64 , nomrte(nvplmx)*64
character        ficsui*32
integer          ncelok , nfaiok , nfabok , nsomok
integer          ierror , irtyp  , itysup , ipasup, nbval
integer          irfsup , idbase
integer          ilecec , nberro , ivers
integer          mvls   , ivar   , ip     , icha  , ilayer
integer          ifac   , iel    , ii     , iok
integer          jphyla , jtpvar , jdpvar , jmpvar
integer          jsttio , jdstnt , mstist , mvlsts
integer          mstbor , musbor , mstits , jturb, jtytur
integer          mode   , ipas   , ivl    , nclsto
integer          impaml , impmls

integer, allocatable, dimension(:) :: iflpar
double precision, allocatable, dimension(:,:) :: coopar

!===============================================================================

!===============================================================================
! 1. Initialisations par defaut
!===============================================================================

!---> Il faut faire dans cette routine les initialisations des
!     tableaux lagrangiens ouverts dans la routine MEMLA1
!     (sauf ITYCEL et ICOCEL qui sont initialises dans LAGDEB),

do ivar = 1,nvp
  do ip = 1,nbpmax
    ettp(ip,ivar) = 0.d0
  enddo
enddo

do ivar = 1,nivep
  do ip = 1,nbpmax
    itepa(ip,ivar) = 0
  enddo
enddo

do ivar = 1,nvep
  do ip = 1,nbpmax
    tepa(ip,ivar) = 0.d0
  enddo
enddo

if (istala.eq.1) then
  do ipas  = 0,nbclst
    do ivl = 1,nvlsta
      ivar = ipas*nvlsta +ivl
      do iel = 1,ncel
        statis(iel,ivar) = 0.d0
      enddo
    enddo
  enddo
  do ipas  = 0,nbclst
    do ivl = 1,nvlsta-1
      ivar = ipas*(nvlsta-1) +ivl
      do iel = 1,ncel
        stativ(iel,ivar) = 0.d0
      enddo
    enddo
  enddo
endif

if (iilagr.eq.2) then
  do ivar = 1,ntersl
    do iel = 1,ncel
      tslagr(iel,ivar) = 0.d0
    enddo
  enddo
endif

if (iensi3.eq.1 .and. nvisbr.gt.0) then
  do ivar = 1,nvisbr
    do ifac = 1,nfabor
      parbor(ifac,ivar) = 0.d0
    enddo
  enddo
endif

if (isuila.eq.0) return

!===============================================================================
! 2. LECTURE DU FICHIER SUITE : VARIABLES LIEES AUX PARTICULES
!===============================================================================

!  ---> Ouverture

write(nfecra,6000)

!     (ILECEC=1:lecture)
ilecec = 1
ficsui = 'lagrangian'
call opnsui(ficsui,len(ficsui),ilecec,impaml,ierror)
!==========
if(ierror.ne.0) then
  write(nfecra,9010) ficsui
  call csexit (1)
endif

write(nfecra,6010)

!  ---> Type de fichier suite
!        Pourrait porter le numero de version si besoin.

itysup = 0
nbval  = 1
irtyp  = 1
rubriq = 'version_fichier_suite_Lagrangien_variables'
call lecsui(impaml,rubriq,len(rubriq),itysup,nbval,irtyp,         &
            ivers,ierror)

if(ierror.ne.0) then
  write(nfecra,9020) ficsui
  call csexit (1)
endif

!  ---> Tests

iok = 0

!     Dimensions des supports

call tstsui(impaml,ncelok,nfaiok,nfabok,nsomok)
!==========
if(ncelok.eq.0) then
  write(nfecra,9030) ficsui
  iok = iok + 1
endif

if(nfaiok.eq.0) write(nfecra,9031) ficsui,'internes','internes'

if(nfabok.eq.0) write(nfecra,9031) ficsui,'de bord ','de bord '

!     Nombre de particules dans le domaine du calcul

itysup = 0
nbval  = 1
irtyp  = 1

rubriq = 'particles'
call lipsui(impaml,rubriq,len(rubriq),nbpart,ipasup)

if(ierror.ne.0) then
  write(nfecra,9040) ficsui,                                      &
  'nombre_courant_particules                                   '
  iok = iok + 1
endif
if(nbpart.gt.nbpmax) then
  write(nfecra,9050) ficsui, nbpart, nbpmax
  iok = iok + 1
endif

!     Physique associee aux particules

rubriq = 'indicateur_physique_particules'
irtyp  = 1
call lecsui(impaml,rubriq,len(rubriq),itysup,nbval,irtyp,         &
            jphyla,ierror)
if(ierror.ne.0) then
  write(nfecra,9040) ficsui, rubriq
  iok = iok + 1
endif

rubriq = 'indicateur_temperature_particules'
irtyp  = 1
call lecsui(impaml,rubriq,len(rubriq),itysup,nbval,irtyp,         &
            jtpvar,ierror)
if(ierror.ne.0) then
  write(nfecra,9040) ficsui,                                      &
  'indicateur_temperature_particules                           '
  iok = iok + 1
endif

!     Arret
if(iok.ne.0) then
  call csexit (1)
endif

rubriq = 'indicateur_diametre_particules'
irtyp  = 1
call lecsui(impaml,rubriq,len(rubriq),itysup,nbval,irtyp,         &
            jdpvar,ierror)
if(ierror.ne.0) then
  write(nfecra,9062) ficsui,                                      &
  'indicateur_diametre_particules                              '
  jdpvar = idpvar
endif

rubriq = 'indicateur_masse_particules'
irtyp  = 1
call lecsui(impaml,rubriq,len(rubriq),itysup,nbval,irtyp,         &
            jmpvar,ierror)
if(ierror.ne.0) then
  write(nfecra,9062) ficsui,                                      &
  'indicateur_masse_particules                                 '
  jmpvar = impvar
endif

! ---> On previent si des parametres sont differents

if ( jphyla.ne.iphyla .or.                                        &
     jtpvar.ne.itpvar .or.                                        &
     jdpvar.ne.idpvar .or.                                        &
     jmpvar.ne.impvar      ) then
  write(nfecra,9070) ficsui,                                      &
                     jphyla, jtpvar, jdpvar, jmpvar,              &
                     iphyla, itpvar, idpvar, impvar
endif

! ---> Verification de la compatibilite si changement de thermique

if (jphyla.ne.0 .and. iphyla.eq.0) then
  write(nfecra,9071) ficsui
endif

if (itpvar.eq.1 .and. jtpvar.eq.0) then
  write(nfecra,9072) ficsui, tpart, cppart
endif

if (iphyla.eq.2 .and. jphyla.ne.2) then
  write(nfecra,9073) ficsui
  call csexit (1)
endif

if ( (jphyla.eq.2 .and. iphyla.eq.1) .or.                         &
     (jphyla.eq.1 .and. iphyla.eq.2)      ) then
  write(nfecra,9074) ficsui
  call csexit (1)
endif

! ---> Infos suivi du calcul

rubriq = 'nombre_iterations_Lagrangiennes'
irtyp  = 1
call lecsui(impaml,rubriq,len(rubriq),itysup,nbval,irtyp,         &
            iplas,ierror)
if(ierror.ne.0) then
  write(nfecra,9060) ficsui,                                      &
  'nombre_iterations_Lagrangiennes                             ', &
  'IPLAS',iplas
endif

if(istala.eq.1 .and. isuist.eq.0 .and. iplas.ge.idstnt) then
  write(nfecra,9065) ficsui, isuist, iplas +1, idstnt
  call csexit (1)
endif

if(iensi3.eq.1 .and. isuist.eq.0 .and. iplas.ge.nstbor) then
  write(nfecra,9066) ficsui, isuist, iplas +1, nstbor
  call csexit (1)
endif

rubriq = 'temps_physique_Lagrangien'
irtyp  = 2
call lecsui(impaml,rubriq,len(rubriq),itysup,nbval,irtyp,         &
            ttclag,ierror)
if(ierror.ne.0) then
  write(nfecra,9061) ficsui,                                      &
  'temps_physique_Lagrangien                                   ', &
  'TTCLAG',ttclag
endif

rubriq = 'nombre_total_particules'
irtyp  = 1
call lecsui(impaml,rubriq,len(rubriq),itysup,nbval,irtyp,         &
            nbptot,ierror)
if(ierror.ne.0) then
  write(nfecra,9060) ficsui,                                      &
  'nombre_total_particules                                     ', &
  'NBPTOT',nbptot
endif

rubriq = 'nombre_particules_perdues'
irtyp  = 1
call lecsui(impaml,rubriq,len(rubriq),itysup,nbval,irtyp,         &
            nbpert,ierror)
if(ierror.ne.0) then
  write(nfecra,9060) ficsui,                                      &
  'nombre_particules_perdues                                   ', &
  'NBPERT',nbpert
endif

rubriq = 'nombre_variables_utilisateur'
irtyp  = 1
call lecsui(impaml,rubriq,len(rubriq),itysup,nbval,irtyp,         &
            mvls,ierror)
if(ierror.ne.0) then
  mvls = 0
  if (nvls.gt.0) then
    write(nfecra,9062) ficsui,                                    &
  'nombre_variables_utilisateur                                '
  endif
endif

if (nvls.lt.mvls) then
  write(nfecra,9080) ficsui, mvls, nvls, nvls, nvls
  mvls = nvls
elseif (nvls.gt.mvls ) then
  write(nfecra,9080) ficsui, mvls, nvls, nvls, nvls
endif

! --> Caracteristiques et infos particulaires (ENTIERS)

nberro = 0

nomite(jord1) = 'order_1'

if (nbclst.gt.0) then
  nomite(jclst) = 'numero_groupe_statistiques'
endif
if (iphyla.eq.2) then
  nomite(jinch) = 'numero_charbon'
endif
! deposition submodel
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

allocate(coopar(3,nbpart))

call lepsui(impaml,itepa(1,jisor),coopar,itysup,ierror)

nberro = nberro+ierror

if (ierror.eq.0) then
  do ii = 1, nbpart
    ettp(ii,jxp) = coopar(1,ii)
    ettp(ii,jyp) = coopar(2,ii)
    ettp(ii,jzp) = coopar(3,ii)
  enddo
endif

deallocate(coopar)

allocate(iflpar(nbpart))

irtyp  = 1

rubriq = 'particle_status_flag'
call lecsui(impaml,rubriq,len(rubriq),itysup,nbval,irtyp,       &
            iflpar,ierror)
nberro = nberro+ierror

if (ierror.eq.0) then
  do ii = 1, nbpart
    if (iand(iflpar(ii), 1) .eq. 1) itepa(ii,jisor) = -itepa(ii,jisor)
  enddo
endif

deallocate(iflpar)

irtyp  = 1

do ivar = 1, nivep
  if (ivar.ne.jisor) then
    rubriq = nomite(ivar)
    if (ivar.eq.jdfac) then
      idbase = 1
      irfsup = 3
      call leisui(impaml,rubriq,len(rubriq),itysup,irfsup,idbase,   &
                  itepa(1,ivar),ierror)
    else
      call lecsui(impaml,rubriq,len(rubriq),itysup,nbval,irtyp,     &
                  itepa(1,ivar),ierror)
    endif
    nberro = nberro+ierror
  endif
enddo

! --> Caracteristiques et infos particulaires (REELS)

nomrte(jrval) = 'random_value'
nomrte(jrtsp) = 'temps_sejour_particules'
nomrte(jrpoi) = 'poids_statistiques_particules'
if (iphyla.eq.1 .and. itpvar.eq.1 .and. iirayo.gt.0) then
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

irtyp  = 2

do ivar = 1, nvep
  rubriq = nomrte(ivar)
  call lecsui(impaml,rubriq,len(rubriq),itysup,nbval,irtyp,       &
              tepa(1,ivar),ierror)
  if (ivar.eq.jrval .and. ierror.ne.0) then
    do ii = 1, nbpart
      call random_number(tepa(ii,jrval))
    enddo
  endif
  nberro = nberro+ierror
enddo

dnbpar = 0.d0
do ip  = 1, nbpart
  dnbpar = dnbpar + tepa(ip,jrpoi)
enddo

! --> Variables particulaires

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
if (mvls.gt.0) then
  do ip = 1,mvls
    write(car4,'(i4.4)') ip
    nomnvl(jvls(ip)) = 'variable_supplementaire_'//car4
  enddo
endif

irtyp  = 2

do ivar = jmp, jwf
  if (ivar .lt. jxp .or. ivar.gt.jzp) then
    rubriq = nomnvl(ivar)
    call lecsui(impaml,rubriq,len(rubriq),itysup,nbval,irtyp,       &
                ettp(1,ivar),ierror)
    nberro = nberro+ierror
  endif
enddo

if (iphyla.eq.1 .and. itpvar.eq.1) then

  if (itpvar.eq.1 .and. jtpvar.eq.0) then

    do ip = 1,nbpart
      iel = itepa(ip,jisor)
      ettp(ip,jtp) = tpart
      ettp(ip,jcp) = cppart
      if ( ippmod(iccoal).ge.0 .or.                               &
           ippmod(icpl3c).ge.0 .or.                               &
           ippmod(icfuel).ge.0      ) then
        ettp(ip,jtf) = propce(iel,ipproc(itemp1)) -tkelvi
      else if ( ippmod(icod3p).ge.0 .or.                          &
                ippmod(icoebu).ge.0 .or.                          &
                ippmod(ielarc).ge.0 .or.                          &
                ippmod(ieljou).ge.0      ) then
        ettp(ip,jtf) = propce(iel,ipproc(itemp)) -tkelvi
      else if (itherm.eq.1 .and. itpscl.eq.2) then
        ettp(ip,jtf) = rtpa(iel,isca(iscalt)) -tkelvi
      else if (itherm.eq.1 .and. itpscl.eq.1) then
        ettp(ip,jtf) = rtpa(iel,isca(iscalt))
      else if (itherm.eq.2) then
        mode = 1
        call usthht(mode, rtpa(iel,isca(iscalt)), ettp(ip,jtf))
      endif
    enddo

  else if (itpvar.eq.1 .and. jtpvar.eq.1) then
    do ivar = jtp,jcp
      rubriq = nomnvl(ivar)
      call lecsui(impaml,rubriq,len(rubriq),itysup,nbval,irtyp,   &
                  ettp(1,ivar),ierror)
      nberro = nberro+ierror
    enddo
  endif

else if (iphyla.eq.2) then
  do ivar = jhp(1),jcp
    rubriq = nomnvl(ivar)
    call lecsui(impaml,rubriq,len(rubriq),itysup,nbval,irtyp,     &
                ettp(1,ivar),ierror)
    nberro = nberro+ierror
  enddo
endif

if (mvls.gt.0) then
  do ivar = 1,mvls
    rubriq = nomnvl(jvls(ivar))
    call lecsui(impaml,rubriq,len(rubriq),itysup,nbval,irtyp,     &
                ettp(1,jvls(ivar)),ierror)
    nberro = nberro+ierror
  enddo
endif

!  ---> Si pb : arret


if(nberro.ne.0) then
  write(nfecra,9041) ficsui
  call csexit (1)
endif

write(nfecra,6011)

!  ---> Fermeture du fichier suite

call clssui(impaml,ierror)

if(ierror.ne.0) then
  write(nfecra,9090) ficsui
endif

write(nfecra,6099)


!===============================================================================
! 3. LECTURE DU FICHIER SUITE STATISTIQUES ET TERMES SOURCES
!    DE COUPLAGE RETOUR
!===============================================================================

if (isuist.eq.1) then

!  ---> Ouverture

  write(nfecra,7000)

  ! (ILECEC=1:lecture)
  ilecec = 1
  ficsui = 'lagrangian_stats'
  call opnsui(ficsui,len(ficsui),ilecec,impmls,ierror)
  !==========
  if(ierror.ne.0) then
    write(nfecra,9010) ficsui
    call csexit (1)
  endif

  write(nfecra,7010)

!  ---> Type de fichier suite
!        Pourrait porter le numero de version si besoin.
!        On ne se sert pas de IVERS pour le moment

  itysup = 0
  nbval  = 1
  irtyp  = 1

  rubriq = 'version_fichier_suite_Lagrangien_statistiques'
  call lecsui(impmls,rubriq,len(rubriq),itysup,nbval,irtyp,       &
              ivers,ierror)
  if(ierror.ne.0) then
    write(nfecra,9020) ficsui
    call csexit (1)
  endif

  rubriq = 'indicateur_ecoulement_stationnaire'
  call lecsui(impmls,rubriq,len(rubriq),itysup,nbval,irtyp,       &
              jsttio,ierror)
  if (ierror.ne.0) then
    write(nfecra,9040) ficsui,                                    &
  'indicateur_ecoulement_stationnaire                          '
    call csexit (1)
  endif

!     Dimensions des supports

  call tstsui(impmls,ncelok,nfaiok,nfabok,nsomok)
  !==========
  if(ncelok.eq.0) then
    write(nfecra,9030) ficsui
    call csexit (1)
  endif

  if(nfaiok.eq.0) write(nfecra,9031) ficsui,'internes','internes'

  if(nfabok.eq.0) write(nfecra,9031) ficsui,'de bord ','de bord '

! --> Est-on cense lire une suite de stats volumiques ?

  if (istala.eq.1 .and. iplas.ge.idstnt) then

    nberro = 0
    itysup = 0
    nbval  = 1

    rubriq = 'iteration_debut_statistiques'
    irtyp  = 1
    call lecsui(impmls,rubriq,len(rubriq),itysup,nbval,irtyp,     &
                jdstnt,ierror)
    nberro = nberro+ierror

    rubriq = 'iteration_debut_statistiques_stationnaires'
    irtyp  = 1
    call lecsui(impmls,rubriq,len(rubriq),itysup,nbval,irtyp,     &
                mstist,ierror)
    nberro = nberro+ierror

!  ---> S'il y a des erreurs, on suppose que c'est parce que le fichier
!         suite ne contient pas d'infos sur les stats volumiques.
!         Dans ce cas, si on est en instationnaire on se dit que c'est
!         pas grave, on saute l'etape et on continue. Par contre si on
!         est dans une configuration de calcul de stats volumiques
!         en stationnaire on stoppe.

    if(nberro.ne.0) then
      if ( isttio.eq.0 .or.                                       &
          (isttio.eq.1 .and. iplas.lt.nstist) ) then
        write(nfecra,9110) ficsui, isttio, idstnt, nstist, iplas+1
        goto 9991
      else
        write(nfecra,9120) ficsui, isttio, idstnt, nstist, iplas+1
        call csexit (1)
      endif
    endif

! --> A partir d'ici on considere que le fichier suite contient
!       des stats volumiques

    if ( jsttio.ne.isttio .or.                                    &
         jdstnt.ne.idstnt .or.                                    &
         mstist.ne.nstist     ) then
      write (nfecra,9130) ficsui,                                 &
                          jsttio, jdstnt, mstist,                 &
                          isttio, idstnt, nstist
    endif

!  --> Lecture de l'avancement du calcul stats volumiques

    rubriq = 'nombre_iterations_statistiques_stationnaires'
    irtyp  = 1
    call lecsui(impmls,rubriq,len(rubriq),itysup,nbval,irtyp,     &
                npst,ierror)
    if(ierror.ne.0) then
      write(nfecra,9060) ficsui,                                  &
  'nombre_iterations_statistiques_stationnaires                ', &
      'NPST',npst
    endif

    rubriq = 'temps_statistiques_stationnaires'
    irtyp  = 2
    call lecsui(impmls,rubriq,len(rubriq),itysup,nbval,irtyp,     &
                tstat,ierror)
    if(ierror.ne.0) then
      write(nfecra,9061) ficsui,                                  &
  'temps_statistiques_stationnaires                            ', &
      'TSTAT',tstat
    endif

    rubriq = 'classe_statistique_particules'
    irtyp  = 1
    call lecsui(impmls,rubriq,len(rubriq),itysup,nbval,irtyp,     &
                nclsto,ierror)
    if(ierror.ne.0) then
      write(nfecra,9061) ficsui,                                  &
  'classes_statistiques                                        ', &
      'NBCLST',nclsto
    endif

!  --> Verif de coherence de l'avancement du calcul avec les
!       indicateurs de calcul de la suite actuelle :

!    1) Amont Instationnaire -> Actuel Instationnaire : OK
!         (NPST = 0)
!                            -> Actuel Stationnaire : Exit sauf debut

!    2) Amont Stationnaire   -> Actuel Instationnaire : OK
!         (NPST > 0)                            (pertes Stats amont)
!                            -> Actuel Stationnaire : OK si IDSTNT et
!                                 NSTIST n'ont pas change, sinon Exit)

    if (npst.eq.0 .and. (isttio.eq.1 .and. nstist.le.iplas)) then
      write(nfecra,9140) ficsui, iplas+1, nstist
      call csexit (1)
    endif

    if ( npst.gt.0 .and.                                          &
        ( (isttio.eq.1 .and. iplas.le.nstist) .or.                &
           isttio.eq.0)                              ) then
      write(nfecra,9141) ficsui
    endif

    if (npst.gt.0 .and. (isttio.eq.1 .and. iplas.ge.nstist)) then
     if (  jdstnt.ne.idstnt .or.                                  &
           mstist.ne.nstist      ) then
        write(nfecra,9142) ficsui
        call csexit (1)
      endif
    endif

    if ( nbclst .ne. nclsto ) then
      write(nfecra,9143) ficsui
      call csexit (1)
    endif

! --> Stats supplementaires utilisateurs

    rubriq = 'nombre_statistiques_utilisateur'
    irtyp  = 1
    call lecsui(impmls,rubriq,len(rubriq),itysup,nbval,irtyp,     &
                mvlsts,ierror)

    if (nvlsts.lt.mvlsts) then
      write(nfecra,9150) ficsui, mvlsts, nvlsts, nvlsts, nvlsts
    endif

!  --> Lecture des Statistiques volumiques. Pas de traitement d'erreurs,
!        on suppose qu'elles sont dues a un changement de physique.

    itysup = 1
    irtyp  = 2
    nbval  = 1

    do ipas  = 0,nbclst
      do ivl = 1,nvlsta
        ivar  = ipas*nvlsta +ivl
        if (ipas.gt.0) then
          write(car4,'(i4.4)') ipas
          rubriq = 'moy_stat_vol_groupe_'//car4//'_'//nomlag(ivar)
        else
          rubriq = 'moy_stat_vol_'//nomlag(ivar)
        endif
        call lecsui(impmls,rubriq,len(rubriq),itysup,nbval,irtyp, &
                    statis(1,ivar),ierror)
      enddo

      do ivl = 1,nvlsta-1
        ivar  = ipas*nvlsta +ivl
        if (ipas.gt.0) then
          write(car4,'(i4.4)') ipas
          rubriq = 'var_stat_vol_groupe_'//car4//'_'//nomlav(ivar)
        else
          rubriq = 'var_stat_vol_'//nomlav(ivar)
        endif
        call lecsui(impmls,rubriq,len(rubriq),itysup,nbval,irtyp, &
                    stativ(1,ivar),ierror)
      enddo
    enddo

  endif


 9991   continue

  if (iensi3.eq.1 .and. nvisbr.gt.0 .and. nfabok.ne.0) then

    itysup = 0
    nbval  = 1

    rubriq = 'iteration_debut_stats_frontieres_stationnaires'
    irtyp  = 1
    call lecsui(impmls,rubriq,len(rubriq),itysup,nbval,irtyp,     &
                mstbor,ierror)

!  ---> S'il y a une erreur, on suppose que c'est parce que le fichier
!         suite ne contient pas d'infos sur les stats aux frontieres.
!         Dans ce cas, si on est en instationnaire on se dit que c'est
!         pas grave, on saute l'etape et on continue. Par contre si on
!         est dans une configuration de calcul de stats aux frontieres
!         en stationnaire on stoppe.

    if(ierror.ne.0) then
      if ( isttio.eq.0 .or.                                       &
          (isttio.eq.1 .and. iplas.lt.nstbor) ) then
        write(nfecra,9210) ficsui, isttio, nstbor, iplas+1
        goto 9992
      else
        write(nfecra,9220) ficsui, isttio, nstbor, iplas+1
        call csexit (1)
      endif
    endif

! --> A partir d'ici on considere que le fichier suite contient
!       des stats volumiques

    if ( jsttio.ne.isttio .or.                                    &
         mstbor.ne.nstbor     ) then
      write (nfecra,9230) ficsui,                                 &
                          jsttio, mstbor,                         &
                          isttio, nstbor
    endif

!  --> Lecture de l'avancement du calcul stats aux frontieres

    rubriq = 'nombre_iterations_stats_frontieres'
    irtyp  = 1
    call lecsui(impmls,rubriq,len(rubriq),itysup,nbval,irtyp,     &
                npstft,ierror)
    if(ierror.ne.0) then
      write(nfecra,9060) ficsui,                                  &
  'nombre_iterations_stats_frontieres                          ', &
      'NPSTFT',npstft
    endif

    rubriq = 'nombre_iterations_stats_frontieres_stationnaires'
    irtyp  = 1
    call lecsui(impmls,rubriq,len(rubriq),itysup,nbval,irtyp,     &
                npstf,ierror)
    if(ierror.ne.0) then
      write(nfecra,9060) ficsui,                                  &
  'nombre_iterations_stats_frontieres_stationnaires            ', &
      'NPSTF',npstf
    endif

    rubriq = 'temps_stats_frontieres_stationnaires'
    irtyp  = 2
    call lecsui(impmls,rubriq,len(rubriq),itysup,nbval,irtyp,     &
                tstatp,ierror)
    if(ierror.ne.0) then
      write(nfecra,9060) ficsui,                                  &
  'temps_stats_frontieres_stationnaires                        ', &
      'TSTATP',tstatp
    endif

!  --> Verif de coherence de l'avancement du calcul avec les
!       indicateurs de calcul de la suite actuelle :


    if (npstf.eq.0 .and. (isttio.eq.1 .and. nstbor.le.iplas)) then
      write(nfecra,9240) ficsui, iplas+1, nstbor
      call csexit (1)
    endif

    if ( npstf.gt.0 .and.                                         &
        ( (isttio.eq.1 .and. iplas.le.nstbor) .or.                &
           isttio.eq.0)                             ) then
      write(nfecra,9241) ficsui
    endif

    if (npstf.gt.0 .and. (isttio.eq.1 .and. iplas.ge.nstbor)) then
     if (mstbor.ne.nstbor) then
        write(nfecra,9242) ficsui
        call csexit (1)
      endif
    endif

! --> Stats supplementaires utilisateurs

    rubriq = 'nombre_stats_frontieres_utilisateur'
    irtyp  = 1
    call lecsui(impmls,rubriq,len(rubriq),itysup,nbval,irtyp,     &
                musbor,ierror)

    if (nusbor.lt.musbor) then
      write(nfecra,9250) ficsui, musbor, nusbor, nusbor, nusbor
    endif

!  --> Lecture des stats aux frontieres. Pas de traitement d'erreurs,
!        on suppose qu'elles sont dues a un changement de physique.

    itysup = 3
    irtyp  = 2
    nbval  = 1

    do ivar = 1,nvisbr
      rubriq = 'stat_bord_'//nombrd(ivar)
      call lecsui(impmls,rubriq,len(rubriq),itysup,nbval,irtyp,   &
                  parbor(1,ivar),ierror)
    enddo

  endif

 9992   continue

  if (iilagr.eq.2) then

    itysup = 0
    nbval  = 1

    rubriq = 'iteration_debut_termes_sources_stationnaires'
    irtyp  = 1
    call lecsui(impmls,rubriq,len(rubriq),itysup,nbval,irtyp,     &
                mstits,ierror)
    if(ierror.ne.0) then
      write(nfecra,9020) ficsui,                                  &
  'iteration_debut_termes_sources_stationnaires                ', &
      'NSTITS',mstits
    endif

!  ---> S'il y a une erreur, on suppose que c'est parce que le fichier
!         suite ne contient pas d'infos sur les TS de couplage retour.
!         Dans ce cas, si on est en instationnaire on se dit que c'est
!         pas grave, on saute l'etape et on continue. Par contre si on
!         est dans une configuration de calcul de stats aux frontieres
!         en stationnaire on stoppe.

    if(ierror.ne.0) then
      if ( isttio.eq.0 .or.                                       &
          (isttio.eq.1 .and. iplas.lt.nstits) ) then
        write(nfecra,9310) ficsui, isttio, nstits, iplas+1
        goto 9993
      else
        write(nfecra,9320) ficsui, isttio, nstits, iplas+1
        call csexit (1)
      endif
    endif

! --> A partir d'ici on considere que le fichier suite contient
!       des stats volumiques

    rubriq = 'modele_turbulence_termes_sources'
    irtyp  = 1
    call lecsui(impmls,rubriq,len(rubriq),itysup,nbval,irtyp,     &
                jturb,ierror)

    jtytur = jturb/10

    if ( jsttio.ne.isttio .or.                                    &
         mstits.ne.nstits     ) then
      if (jtytur.eq.2) car8 = 'k-eps'
      if (jtytur.eq.3) car8 = 'Rij-eps'
      if (jturb.eq.50) car8 = 'v2f'
      if (jturb.eq.60) car8 = 'k-omega'
      if (itytur.eq.2) kar8 = 'k-eps'
      if (itytur.eq.3) kar8 = 'Rij-eps'
      if (iturb.eq.50) kar8 = 'v2f'
      if (iturb.eq.60) kar8 = 'k-omega'
      write (nfecra,9330) ficsui,                                 &
                          jsttio, mstits, car8,                   &
                          isttio, nstits, kar8
    endif


!  --> Lecture de l'avancement du couplage retour

    rubriq = 'nombre_iterations_termes_sources_stationnaires'
    irtyp  = 1
    call lecsui(impmls,rubriq,len(rubriq),itysup,nbval,irtyp,     &
                npts,ierror)
    if(ierror.ne.0) then
      write(nfecra,9060) ficsui,                                  &
  'nombre_iterations_termes_sources_stationnaires              ', &
      'NPTS',npts
    endif

!  --> Verif de coherence de l'avancement du calcul avec les
!       indicateurs de calcul de la suite actuelle :

    if (npts.eq.0 .and. (isttio.eq.1 .and. nstits.le.iplas)) then
      write(nfecra,9340) ficsui, iplas+1, nstits
      call csexit (1)
    endif

    if ( npts.gt.0 .and.                                          &
        ( (isttio.eq.1 .and. iplas.le.nstits) .or.                &
           isttio.eq.0)                             ) then
      write(nfecra,9341) ficsui
    endif

    if (npts.gt.0 .and. (isttio.eq.1 .and. iplas.ge.nstits)) then
     if (mstits.ne.nstits) then
        write(nfecra,9342) ficsui
        call csexit (1)
      endif
    endif

!       On donne des labels au different TS pour les noms de rubriques
!       On donne le meme label au keps, au v2f et au k-omega (meme variable k)

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

!  Termes source de couplage retour

    itysup = 1
    nbval  = 1
    irtyp  = 2

    do ivar = 1,ntersl
      rubriq = nomtsl(ivar)
      call lecsui(impmls,rubriq,len(rubriq),itysup,nbval,irtyp,   &
                  tslagr(1,ivar),ierror)
    enddo


!  Dans le cas specifique de la combustion de grains de charbon
!  avec un couplage retour sur une combustion gaz en phase porteuse

!      --> A verifier l'utilite de cette lecture pour une suite...

    if (ippmod(icpl3c).ge.0) then
      do ivar = 1, nsalpp
        icha = nsalto-nsalpp+ivar
        itysup = 1
        nbval  = 1
        irtyp  = 2
        write(car4,'(i4.4)') ivar
        rubriq = 'scalaires_physiques_pariculieres_charbon'//car4
        call lecsui(impmls,rubriq,len(rubriq),itysup,nbval,irtyp, &
                    propce(1,ipproc(icha)),ierror)
      enddo
    endif

  endif

 9993   continue

  write(nfecra,6011)

!  ---> Fermeture du fichier suite

  call clssui(impmls,ierror)

  if(ierror.ne.0) write(nfecra,9090) ficsui, ficsui

! ---> En cas d'erreur, on continue quand meme

  write(nfecra,7099)

endif

write(nfecra,2000)

!===============================================================================

!--------
! FORMATS
!--------


 2000 format(                                                     &
'                                                             ',/,&
'-------------------------------------------------------------',/)

 6000 FORMAT (/, 3X,'** INFORMATIONS SUR LE CALCUL LAGRANGIEN',/,  &
           3X,'   -------------------------------------      ',/,  &
           3X,' Lecture d''un fichier suite                  ',/,  &
           3X,'   sur les variables liees aux particules     '  )
 6010 FORMAT (   3X,'   Debut de la lecture                  '  )
 6011 FORMAT (   3X,'   Fin   de la lecture                  '  )
 6099 FORMAT (   3X,' Fin de la lecture du fichier suite     ',/,  &
           3X,'   sur les variables liees aux particules    ',/)

 7000 FORMAT (/, 3X,'** INFORMATIONS SUR LE CALCUL LAGRANGIEN',/,  &
           3X,'   -------------------------------------      ',/,  &
           3X,' Lecture d''un fichier suite                  ',/,  &
           3X,'   sur les statistiques et TS couplage retour'  )
 7010 FORMAT (   3X,'   Debut de la lecture                  '  )
 7099 FORMAT (   3X,' Fin de la lecture du fichier suite     ',/,  &
           3X,'   sur les statistiques et TS couplage retour'  )

 9010 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''OUVERTURE D''UN FICHIER SUITE    ',/,&
'@    =========     LAGRANGIEN ',A13                           ,/,&
'@                                                            ',/,&
'@    Le calcul ne peut pas etre execute.                     ',/,&
'@                                                            ',/,&
'@    Verifier l''existence de ce fichier suite dans le       ',/,&
'@        sous-repertoire ''restart'' du repertoire de travail.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 9020 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,&
'@    =========     LAGRANGIEN ',A13                           ,/,&
'@      TYPE DE FICHIER INCORRECT                             ',/,&
'@                                                            ',/,&
'@    Ce fichier ne semble pas etre un fichier                ',/,&
'@      suite Lagrangien.                                     ',/,&
'@                                                            ',/,&
'@    Le calcul ne peut etre execute.                         ',/,&
'@                                                            ',/,&
'@    Verifier que le fichier suite utilise correspond bien   ',/,&
'@        a un fichier suite Lagrangien.                      ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 9030 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,&
'@    =========     LAGRANGIEN ',A13                           ,/,&
'@      DONNEES AMONT ET ACTUELLES INCOHERENTES               ',/,&
'@                                                            ',/,&
'@    Le nombre de cellules a ete modifie                     ',/,&
'@                                                            ',/,&
'@    Le calcul ne peut etre execute.                         ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 9031 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : A LA LECTURE DU FICHIER SUITE               ',/,&
'@    =========     LAGRANGIEN ',A13                           ,/,&
'@      DONNEES AMONT ET ACTUELLES DIFFERENTES                ',/,&
'@                                                            ',/,&
'@    Le nombre de faces ',A8  ,' a ete modifie.              ',/,&
'@                                                            ',/,&
'@    Le calcul peut etre execute mais les donnees            ',/,&
'@      sur les faces ',A8  ,' ne seront pas relues           ',/,&
'@      dans le fichier suite.                                ',/,&
'@    Elles seront initialisees par des valeurs par defaut.   ',/,&
'@                                                            ',/,&
'@    Le calcul se poursuit...                                ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 9040 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,&
'@    =========     LAGRANGIEN ',A13                           ,/,&
'@                                                            ',/,&
'@      ERREUR A LA LECTURE DE LA RUBRIQUE                    ',/,&
'@ ',A60                                                       ,/,&
'@                                                            ',/,&
'@    Le calcul ne peut pas etre execute.                     ',/,&
'@                                                            ',/,&
'@    Verifier que ce fichier suite                           ',/,&
'@        correspond bien a un fichier suite Lagrangien,      ',/,&
'@        et qu''il n''a pas ete endommage.                   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 9041 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,&
'@    =========     LAGRANGIEN ',A13                           ,/,&
'@                                                            ',/,&
'@      ERREUR A LA LECTURE DES VARIABLES PARTICULAIRES       ',/,&
'@                                                            ',/,&
'@    Le calcul ne peut pas etre execute.                     ',/,&
'@                                                            ',/,&
'@    Verifier que ce fichier suite                           ',/,&
'@        correspond bien a un fichier suite Lagrangien,      ',/,&
'@        et qu''il n''a pas ete endommage.                   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 9050 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,&
'@    =========     LAGRANGIEN ',A13                           ,/,&
'@      DONNEES AMONT ET ACTUELLES INCOHERENTES               ',/,&
'@                                                            ',/,&
'@    Le nombre de particules contenu dans le fichier suite   ',/,&
'@      est superieur au nombre de particule maximal autorise ',/,&
'@      NBPART = ',I10                                         ,/,&
'@      NBPMAX = ',I10                                         ,/,&
'@                                                            ',/,&
'@    Le calcul ne peut etre execute.                         ',/,&
'@                                                            ',/,&
'@    Augmenter la valeur de NBPMAX dans USLAG1.              ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 9060 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : A LA LECTURE DU FICHIER SUITE               ',/,&
'@    =========     LAGRANGIEN ',A13                           ,/,&
'@                                                            ',/,&
'@      ERREUR A LA LECTURE DE LA RUBRIQUE                    ',/,&
'@ ',A60                                                       ,/,&
'@                                                            ',/,&
'@    Le mot cle est initialise avec sa valeur par defaut     ',/,&
'@      ou celle donnee dans le sous-programme USLAG1 :       ',/,&
'@        ',A10   ,'  = ',I10                                  ,/,&
'@                                                            ',/,&
'@    Le calcul se poursuit...                                ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 9061 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : A LA LECTURE DU FICHIER SUITE               ',/,&
'@    =========     LAGRANGIEN ',A13                           ,/,&
'@                                                            ',/,&
'@      ERREUR A LA LECTURE DE LA RUBRIQUE                    ',/,&
'@ ',A60                                                       ,/,&
'@                                                            ',/,&
'@    Le mot cle est initialise avec sa valeur par defaut     ',/,&
'@      ou celle donnee dans le sous-programme USLAG1 :       ',/,&
'@        ',A10   ,'  = ',E14.5                                ,/,&
'@                                                            ',/,&
'@    Le calcul se poursuit...                                ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 9062 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : A LA LECTURE DU FICHIER SUITE               ',/,&
'@    =========     LAGRANGIEN ',A13                           ,/,&
'@                                                            ',/,&
'@      ERREUR A LA LECTURE DE LA RUBRIQUE                    ',/,&
'@ ',A60                                                       ,/,&
'@                                                            ',/,&
'@    Le calcul se poursuit...                                ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 9065 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : A LA LECTURE DU FICHIER SUITE               ',/,&
'@    =========     LAGRANGIEN ',A13                           ,/,&
'@      DONNEES AMONT ET ACTUELLES INCOHERENTES               ',/,&
'@                                                            ',/,&
'@    L''INDICATEUR DE CALCUL DES STATISTIQUES VOLUMIQUES     ',/,&
'@       A UNE VALEUR NON PERMISE (LAGLEC).                   ',/,&
'@                                                            ',/,&
'@    LORSQU''IL N''Y A PAS DE SUITE DE CALCUL SUR LES        ',/,&
'@    STATISTIQUES VOLUMIQUES (ISUIST = ',  I3, '),'           ,/,&
'@    IDSTNT DOIT ETRE UN ENTIER SUPERIEUR AU NUMERO          ',/,&
'@       DE L''ITERATION LAGRANGIENNE DE REDEMARRAGE ', I10    ,/,&
'@       IL VAUT ICI IDSTNT = ', I10                           ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier la valeur de IDSTNT dans la subroutine USLAG1.   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 9066 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : A LA LECTURE DU FICHIER SUITE               ',/,&
'@    =========     LAGRANGIEN ',A13                           ,/,&
'@      DONNEES AMONT ET ACTUELLES INCOHERENTES               ',/,&
'@                                                            ',/,&
'@    L''INDICATEUR DE CALCUL STATIONNAIRES DES STATISTIQUES  ',/,&
'@       AUX FRONTIERES A UNE VALEUR NON PERMISE (LAGLEC).    ',/,&
'@                                                            ',/,&
'@    LORSQU''IL N''Y A PAS DE SUITE DE CALCUL SUR LES        ',/,&
'@    STATISTIQUES AUX FRONTIERES (ISUIST = ',  I3, '),'       ,/,&
'@    NSTBOR DOIT ETRE UN ENTIER SUPERIEUR AU NUMERO          ',/,&
'@       DE L''ITERATION LAGRANGIENNE DE REDEMARRAGE ', I10    ,/,&
'@       IL VAUT ICI NSTBOR = ', I10                           ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier la valeur de NSTBOR dans la subroutine USLAG1.   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 9070 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : A LA LECTURE DU FICHIER SUITE               ',/,&
'@    =========     LAGRANGIEN ',A13                           ,/,&
'@      DONNEES AMONT ET ACTUELLES DIFFERENTES                ',/,&
'@                                                            ',/,&
'@    Les indicateurs concernant la physique associee         ',/,&
'@      aux particules sont modifies :                        ',/,&
'@                                                            ',/,&
'@              IPHYLA    ITPVAR    IDPVAR    IMPVAR          ',/,&
'@  AMONT : ',4I10                                             ,/,&
'@  ACTUEL: ',4I10                                             ,/,&
'@                                                            ',/,&
'@    Le calcul se poursuit...                                ',/,&
'@                                                            ',/,&
'@    Il est conseille de verifier ces indicateurs dans       ',/,&
'@      le sous-programme USLAG1.                             ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 9071 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : A LA LECTURE DU FICHIER SUITE               ',/,&
'@    =========     LAGRANGIEN ',A13                           ,/,&
'@      DONNEES AMONT ET ACTUELLES DIFFERENTES                ',/,&
'@                                                            ',/,&
'@    Aucune selection de physique associee aux particules    ',/,&
'@      n''est active. Les donnees amont sont perdues.        ',/,&
'@                                                            ',/,&
'@    Le calcul se poursuit...                                ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 9072 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : A LA LECTURE DU FICHIER SUITE               ',/,&
'@    =========     LAGRANGIEN ',A13                           ,/,&
'@      DONNEES AMONT ET ACTUELLES DIFFERENTES                ',/,&
'@                                                            ',/,&
'@    Une equation sur la temperature des particules est      ',/,&
'@      enclenchee en cours de calcul.                        ',/,&
'@    Initialisation par defaut :                             ',/,&
'@       Temperature TPART = ', E14.5                          ,/,&
'@       Chaleur massique CPPART = ', E14.5                    ,/,&
'@                                                            ',/,&
'@    Le calcul se poursuit...                                ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 9073 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,&
'@    =========     LAGRANGIEN ',A13                           ,/,&
'@      DONNEES AMONT ET ACTUELLES INCOHERENTES               ',/,&
'@                                                            ',/,&
'@    L''indicateur d''un calcul Lagrangien de grains         ',/,&
'@      de charbon est enclenche (IPHYLA = 2).                ',/,&
'@    Ce fichier suite ne correspond pas                      ',/,&
'@      a un calcul Lagrangien de grains de charbon.          ',/,&
'@                                                            ',/,&
'@    Le calcul ne peut etre execute.                         ',/,&
'@                                                            ',/,&
'@    Verifier IPHYLA dans le sous-programme USLAG1.          ',/,&
'@    Verifier le fichier suite utilise.                      ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 9074 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,&
'@    =========     LAGRANGIEN ',A13                           ,/,&
'@      DONNEES AMONT ET ACTUELLES INCOHERENTES               ',/,&
'@                                                            ',/,&
'@    Ce fichier suite correspond                             ',/,&
'@      a un calcul Lagrangien de grains de charbon.          ',/,&
'@    L''indicateur de physique actuel associee aux particules',/,&
'@      a une valeur non permise dans le cadre d''une suite   ',/,&
'@      d''un calcul Lagrangien de grains de charbon.         ',/,&
'@                                                            ',/,&
'@    Le calcul ne peut etre execute.                         ',/,&
'@                                                            ',/,&
'@    Verifier IPHYLA dans le sous-programme USLAG1.          ',/,&
'@    Verifier le fichier suite utilise.                      ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 9080 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : A LA LECTURE DU FICHIER SUITE               ',/,&
'@    =========     LAGRANGIEN ',A13                           ,/,&
'@      DONNEES AMONT ET ACTUELLES DIFFERENTES                ',/,&
'@                                                            ',/,&
'@    L''indicateur du  nombre de variables supplementaires   ',/,&
'@      utilisateur est modifie, ou n''a pas pu etre relu.    ',/,&
'@                                                            ',/,&
'@              NVLS                                          ',/,&
'@    AMONT : ',I10   ,'      ACTUEL : ',I10                   ,/,&
'@                                                            ',/,&
'@    Si ACTUEL > AMONT, on initialise les ',I10   ,' 1eres   ',/,&
'@      variables supplementaires actuelles avec celles       ',/,&
'@      du fichier suite, les autres sont initialisees a zero.',/,&
'@                                                            ',/,&
'@    Si ACTUEL < AMONT, on initialise les ',I10   ,' 1eres   ',/,&
'@      variables supplementaires actuelles avec les 1eres    ',/,&
'@      du fichier suite, le reste des variables du fichier   ',/,&
'@      suite sont perdues.                                   ',/,&
'@                                                            ',/,&
'@    Le calcul se poursuit...                                ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 9090 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : A LA FERMETURE DU FICHIER SUITE             ',/,&
'@    =========     LAGRANGIEN ',A13                           ,/,&
'@                                                            ',/,&
'@    Le calcul continue...                                   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 9110 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : A LA LECTURE DU FICHIER SUITE               ',/,&
'@    =========     LAGRANGIEN ',A13                           ,/,&
'@                                                            ',/,&
'@      LES STATISTIQUES VOLUMIQUES DU CALCUL AMONT           ',/,&
'@        NE PEUVENT PAS ETRE RELUES OU SONT ABSENTES         ',/,&
'@        DU FICHIER SUITE                                    ',/,&
'@                                                            ',/,&
'@    Les indicateurs concernant le calcul                    ',/,&
'@      instationnaire/stationnaire des                       ',/,&
'@      statistiques volumiques sont positionnes en mode      ',/,&
'@      instationnaire ou en debut de calcul stationnaire :   ',/,&
'@                                                            ',/,&
'@          ISTTIO    IDSTNT    NSTIST    Iter de redemarrage ',/,&
'@      ',4I10                                                 ,/,&
'@                                                            ',/,&
'@    Elles seront initialisees par des valeurs par defaut.   ',/,&
'@                                                            ',/,&
'@    Le calcul se poursuit...                                ',/,&
'@                                                            ',/,&
'@    Il est conseille de verifier ces indicateurs dans       ',/,&
'@      le sous-programme USLAG1.                             ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 9120 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,&
'@    =========     LAGRANGIEN ',A13                           ,/,&
'@                                                            ',/,&
'@      LES STATISTIQUES VOLUMIQUES DU CALCUL AMONT           ',/,&
'@        NE PEUVENT PAS ETRE RELUES OU SONT ABSENTES         ',/,&
'@        DU FICHIER SUITE                                    ',/,&
'@                                                            ',/,&
'@    Les indicateurs concernant le calcul                    ',/,&
'@      instationnaire/stationnaire des                       ',/,&
'@      statistiques volumiques sont positionnes              ',/,&
'@      en mode stationnaire :                                ',/,&
'@                                                            ',/,&
'@          ISTTIO    IDSTNT    NSTIST    Iter de redemarrage ',/,&
'@      ',4I10                                                 ,/,&
'@                                                            ',/,&
'@    Le calcul ne peut pas etre execute.                     ',/,&
'@                                                            ',/,&
'@    Verifier ces indicateurs dans le sous-programme USLAG1. ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 9130 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : A LA LECTURE DU FICHIER SUITE               ',/,&
'@    =========     LAGRANGIEN ',A13                           ,/,&
'@      DONNEES AMONT ET ACTUELLES DIFFERENTES                ',/,&
'@                                                            ',/,&
'@    Les indicateurs concernant le calcul                    ',/,&
'@      instationnaire/stationnaire des                       ',/,&
'@      statistiques volumiques sont modifies :               ',/,&
'@                                                            ',/,&
'@              ISTTIO    IDSTNT    NSTIST                    ',/,&
'@  AMONT : ',3I10                                             ,/,&
'@  ACTUEL: ',3I10                                             ,/,&
'@                                                            ',/,&
'@    Le calcul se poursuit...                                ',/,&
'@                                                            ',/,&
'@    Il est conseille de verifier ces indicateurs dans       ',/,&
'@      le sous-programme USLAG1.                             ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 9140 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,&
'@    =========     LAGRANGIEN ',A13                           ,/,&
'@      DONNEES AMONT ET ACTUELLES INCOHERENTES               ',/,&
'@                                                            ',/,&
'@    LE CALCUL ACTUEL DES STATISTIQUES VOLUMIQUES            ',/,&
'@      EST EN MODE STATIONNAIRE, ALORS QUE LE FICHIER        ',/,&
'@      SUITE CONTIENT DES STATISTIQUES INSTATIONNAIRES.      ',/,&
'@                                                            ',/,&
'@    NSTIST devrait etre un entier superieur ou egal         ',/,&
'@      a l''iteration Lagrangienne absolue de redemarrage    ',/,&
'@      du calcul (iteration : ',I10   ,')                    ',/,&
'@                                                            ',/,&
'@      Il vaut ici NSTIST = ',I10                             ,/,&
'@                                                            ',/,&
'@    Le calcul ne peut pas etre execute.                     ',/,&
'@                                                            ',/,&
'@    Verifier ces indicateurs dans le sous-programme USLAG1. ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 9141 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : A LA LECTURE DU FICHIER SUITE               ',/,&
'@    =========     LAGRANGIEN ',A13                           ,/,&
'@      DONNEES AMONT ET ACTUELLES DIFFERENTES                ',/,&
'@                                                            ',/,&
'@    LE CALCUL ACTUEL DES STATISTIQUES VOLUMIQUES            ',/,&
'@      EST EN MODE INSTATIONNAIRE, ALORS QUE LE FICHIER      ',/,&
'@      SUITE CONTIENT DES STATISTIQUES STATIONNAIRES.        ',/,&
'@                                                            ',/,&
'@    Les statistiques volumiques stationnaires amont         ',/,&
'@      seront remises a zero.                                ',/,&
'@                                                            ',/,&
'@    Le calcul se poursuit...                                ',/,&
'@                                                            ',/,&
'@    Il est conseille de verifier ces indicateurs dans       ',/,&
'@      le sous-programme USLAG1.                             ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 9142 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,&
'@    =========     LAGRANGIEN ',A13                           ,/,&
'@      DONNEES AMONT ET ACTUELLES INCOHERENTES               ',/,&
'@                                                            ',/,&
'@    LE CALCUL SE POURSUIT AVEC UN CALCUL DE                 ',/,&
'@      STATISTIQUES VOLUMIQUES EN MODE STATIONNAIRE          ',/,&
'@      MAIS LES INDICATEURS DE CONTROLES DES STATISTIQUES    ',/,&
'@      ON ETE MODIFIEES.                                     ',/,&
'@                                                            ',/,&
'@    Pour eviter les incoherences dans le calcul             ',/,&
'@      IDSTNT et NSTIST ne devraient pas etre modifies entre ',/,&
'@      deux calculs de statistiques volumiques stationnaires.',/,&
'@                                                            ',/,&
'@    Le calcul ne peut pas etre execute.                     ',/,&
'@                                                            ',/,&
'@    Verifier ces indicateurs dans le sous-programme USLAG1. ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 9143 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,&
'@    =========     LAGRANGIEN ',A13                           ,/,&
'@      DONNEES AMONT ET ACTUELLES INCOHERENTES               ',/,&
'@                                                            ',/,&
'@    LE CALCUL SE POURSUIT AVEC UN CALCUL DE                 ',/,&
'@      STATISTIQUES VOLUMIQUES EN MODE STATIONNAIRE          ',/,&
'@      MAIS LES INDICATEURS DE CONTROLES DES STATISTIQUES    ',/,&
'@      ON ETE MODIFIEES.                                     ',/,&
'@                                                            ',/,&
'@    Pour eviter les incoherences dans le calcul             ',/,&
'@      NBCLST ne devrait pas etre modifies entre             ',/,&
'@      deux calculs de statistiques volumiques.              ',/,&
'@                                                            ',/,&
'@    Le calcul ne peut pas etre execute.                     ',/,&
'@                                                            ',/,&
'@    Verifier ces indicateurs dans le sous-programme USLAG1. ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 9150 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : A LA LECTURE DU FICHIER SUITE               ',/,&
'@    =========     LAGRANGIEN ',A13                           ,/,&
'@      DONNEES AMONT ET ACTUELLES DIFFERENTES                ',/,&
'@                                                            ',/,&
'@    L''indicateur du  nombre de statistiques volumiques     ',/,&
'@      supplementaires utilisateur est modifie,              ',/,&
'@      ou n''a pas pu etre relu.                             ',/,&
'@                                                            ',/,&
'@              NVLSTS                                        ',/,&
'@    AMONT : ',I10   ,'      ACTUEL : ',I10                   ,/,&
'@                                                            ',/,&
'@    Si ACTUEL > AMONT, on initialise les ',I10   ,' 1eres   ',/,&
'@      statistiques supplementaires actuelles avec celles    ',/,&
'@      du fichier suite, les autres sont initialisees a zero.',/,&
'@                                                            ',/,&
'@    Si ACTUEL < AMONT, on initialise les ',I10   ,' 1eres   ',/,&
'@      statistiques supplementaires actuelles avec les 1eres ',/,&
'@      du fichier suite, le reste des statistiques du fichier',/,&
'@      suite sont perdues.                                   ',/,&
'@                                                            ',/,&
'@    Le calcul se poursuit...                                ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 9210 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : A LA LECTURE DU FICHIER SUITE               ',/,&
'@    =========     LAGRANGIEN ',A13                           ,/,&
'@                                                            ',/,&
'@      LES STATISTIQUES AUX FRONTIERES DU CALCUL AMONT       ',/,&
'@        NE PEUVENT PAS ETRE RELUES OU SONT ABSENTES         ',/,&
'@        DU FICHIER SUITE                                    ',/,&
'@                                                            ',/,&
'@    Les indicateurs concernant le calcul                    ',/,&
'@      instationnaire/stationnaire des                       ',/,&
'@      statistiques aux frontieres sont positionnes en mode  ',/,&
'@      instationnaire ou en debut de calcul stationnaire :   ',/,&
'@                                                            ',/,&
'@          ISTTIO    NSTBOR    Iter de redemarrage           ',/,&
'@      ',3I10                                                 ,/,&
'@                                                            ',/,&
'@    Elles seront initialisees par des valeurs par defaut.   ',/,&
'@                                                            ',/,&
'@    Le calcul se poursuit...                                ',/,&
'@                                                            ',/,&
'@    Il est conseille de verifier ces indicateurs dans       ',/,&
'@      le sous-programme USLAG1.                             ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 9220 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,&
'@    =========     LAGRANGIEN ',A13                           ,/,&
'@                                                            ',/,&
'@      LES STATISTIQUES AUX FRONTIERES DU CALCUL AMONT       ',/,&
'@        NE PEUVENT PAS ETRE RELUES OU SONT ABSENTES         ',/,&
'@        DU FICHIER SUITE                                    ',/,&
'@                                                            ',/,&
'@    Les indicateurs concernant le calcul                    ',/,&
'@      instationnaire/stationnaire des                       ',/,&
'@      statistiques aux frontieres sont positionnes          ',/,&
'@      en mode stationnaire :                                ',/,&
'@                                                            ',/,&
'@          ISTTIO    NSTBOR    Iter de redemarrage           ',/,&
'@      ',3I10                                                 ,/,&
'@                                                            ',/,&
'@    Le calcul ne peut pas etre execute.                     ',/,&
'@                                                            ',/,&
'@    Verifier ces indicateurs dans le sous-programme USLAG1. ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 9230 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : A LA LECTURE DU FICHIER SUITE               ',/,&
'@    =========     LAGRANGIEN ',A13                           ,/,&
'@      DONNEES AMONT ET ACTUELLES DIFFERENTES                ',/,&
'@                                                            ',/,&
'@    Les indicateurs concernant le calcul                    ',/,&
'@      instationnaire/stationnaire des                       ',/,&
'@      statistiques aux frontieres sont modifies :           ',/,&
'@                                                            ',/,&
'@                        ISTTIO    NSTBOR                    ',/,&
'@            AMONT : ',2I10                                   ,/,&
'@            ACTUEL: ',2I10                                   ,/,&
'@                                                            ',/,&
'@    Le calcul se poursuit...                                ',/,&
'@                                                            ',/,&
'@    Il est conseille de verifier ces indicateurs dans       ',/,&
'@      le sous-programme USLAG1.                             ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 9240 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,&
'@    =========     LAGRANGIEN ',A13                           ,/,&
'@      DONNEES AMONT ET ACTUELLES INCOHERENTES               ',/,&
'@                                                            ',/,&
'@    LE CALCUL ACTUEL DES STATISTIQUES AUX FRONTIERES        ',/,&
'@      EST EN MODE STATIONNAIRE, ALORS QUE LE FICHIER        ',/,&
'@      SUITE CONTIENT DES STATISTIQUES INSTATIONNAIRES.      ',/,&
'@                                                            ',/,&
'@    NSTBOR devrait etre un entier superieur ou egal         ',/,&
'@      a l''iteration Lagrangienne absolue de redemarrage    ',/,&
'@      du calcul (iteration : ',I10   ,')                    ',/,&
'@                                                            ',/,&
'@      Il vaut ici NSTBOR = ',I10                             ,/,&
'@                                                            ',/,&
'@    Le calcul ne peut pas etre execute.                     ',/,&
'@                                                            ',/,&
'@    Verifier ces indicateurs dans le sous-programme USLAG1. ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 9241 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : A LA LECTURE DU FICHIER SUITE               ',/,&
'@    =========     LAGRANGIEN ',A13                           ,/,&
'@      DONNEES AMONT ET ACTUELLES DIFFERENTES                ',/,&
'@                                                            ',/,&
'@    LE CALCUL ACTUEL DES STATISTIQUES AUX FRONTIERES        ',/,&
'@      EST EN MODE INSTATIONNAIRE, ALORS QUE LE FICHIER      ',/,&
'@      SUITE CONTIENT DES STATISTIQUES STATIONNAIRES.        ',/,&
'@                                                            ',/,&
'@    Les statistiques aux frontieres stationnaires amont     ',/,&
'@      seront remises a zero.                                ',/,&
'@                                                            ',/,&
'@    Le calcul se poursuit...                                ',/,&
'@                                                            ',/,&
'@    Il est conseille de verifier ces indicateurs dans       ',/,&
'@      le sous-programme USLAG1.                             ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 9242 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,&
'@    =========     LAGRANGIEN ',A13                           ,/,&
'@      DONNEES AMONT ET ACTUELLES INCOHERENTES               ',/,&
'@                                                            ',/,&
'@    LE CALCUL SE POURSUIT AVEC UN CALCUL DE                 ',/,&
'@      STATISTIQUES AUX FRONTIERES EN MODE STATIONNAIRE      ',/,&
'@      MAIS LES INDICATEURS DE CONTROLES DES STATISTIQUES    ',/,&
'@      ON ETE MODIFIEES.                                     ',/,&
'@                                                            ',/,&
'@    Pour eviter les incoherences dans le calcul             ',/,&
'@      NSTBOR ne devrait pas etre modifie entre deux calculs ',/,&
'@      de statistiques aux frontieres stationnaires.         ',/,&
'@                                                            ',/,&
'@    Le calcul ne peut pas etre execute.                     ',/,&
'@                                                            ',/,&
'@    Verifier ces indicateurs dans le sous-programme USLAG1. ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 9250 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : A LA LECTURE DU FICHIER SUITE               ',/,&
'@    =========     LAGRANGIEN ',A13                           ,/,&
'@      DONNEES AMONT ET ACTUELLES DIFFERENTES                ',/,&
'@                                                            ',/,&
'@    L''indicateur du  nombre de statistiques aux frontieres ',/,&
'@      supplementaires utilisateur est modifie,              ',/,&
'@      ou n''a pas pu etre relu.                             ',/,&
'@                                                            ',/,&
'@              NUSBOR                                        ',/,&
'@    AMONT : ',I10   ,'      ACTUEL : ',I10                   ,/,&
'@                                                            ',/,&
'@    Si ACTUEL > AMONT, on initialise les ',I10   ,' 1eres   ',/,&
'@      statistiques supplementaires actuelles avec celles    ',/,&
'@      du fichier suite, les autres sont initialisees a zero.',/,&
'@                                                            ',/,&
'@    Si ACTUEL < AMONT, on initialise les ',I10   ,' 1eres   ',/,&
'@      statistiques supplementaires actuelles avec les 1eres ',/,&
'@      du fichier suite, le reste des statistiques du fichier',/,&
'@      suite sont perdues.                                   ',/,&
'@                                                            ',/,&
'@    Le calcul se poursuit...                                ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 9310 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : A LA LECTURE DU FICHIER SUITE               ',/,&
'@    =========     LAGRANGIEN ',A13                           ,/,&
'@                                                            ',/,&
'@      LES TERMES SOURCES DE COUPLAGE RETOUR DU CALCUL AMONT ',/,&
'@        NE PEUVENT PAS ETRE RELUES OU SONT ABSENTS          ',/,&
'@        DU FICHIER SUITE                                    ',/,&
'@                                                            ',/,&
'@    Les indicateurs concernant le calcul                    ',/,&
'@      instationnaire/stationnaire des termes sources        ',/,&
'@      de couplage retour sont positionnes en mode           ',/,&
'@      instationnaire ou en debut de calcul stationnaire :   ',/,&
'@                                                            ',/,&
'@          ISTTIO    NSTITS    Iter de redemarrage           ',/,&
'@      ',3I10                                                 ,/,&
'@                                                            ',/,&
'@    Ils seront initialisees par des valeurs par defaut.     ',/,&
'@                                                            ',/,&
'@    Le calcul se poursuit...                                ',/,&
'@                                                            ',/,&
'@    Il est conseille de verifier ces indicateurs dans       ',/,&
'@      le sous-programme USLAG1.                             ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 9320 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,&
'@    =========     LAGRANGIEN ',A13                           ,/,&
'@                                                            ',/,&
'@      LES TERMES SOURCES DE COUPLAGE RETOUR DU CALCUL AMONT ',/,&
'@        NE PEUVENT PAS ETRE RELUES OU SONT ABSENTS          ',/,&
'@        DU FICHIER SUITE                                    ',/,&
'@                                                            ',/,&
'@    Les indicateurs concernant le calcul                    ',/,&
'@      instationnaire/stationnaire des termes sources        ',/,&
'@      de couplage retour sont positionnes                   ',/,&
'@      en mode stationnaire :                                ',/,&
'@                                                            ',/,&
'@          ISTTIO    NSTITS    Iter de redemarrage           ',/,&
'@      ',3I10                                                 ,/,&
'@                                                            ',/,&
'@    Le calcul ne peut pas etre execute.                     ',/,&
'@                                                            ',/,&
'@    Verifier ces indicateurs dans le sous-programme USLAG1. ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 9330 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : A LA LECTURE DU FICHIER SUITE               ',/,&
'@    =========     LAGRANGIEN ',A13                           ,/,&
'@      DONNEES AMONT ET ACTUELLES DIFFERENTES                ',/,&
'@                                                            ',/,&
'@    Les indicateurs concernant le calcul                    ',/,&
'@      instationnaire/stationnaire des termes sources        ',/,&
'@      de couplage retour sont modifies :                    ',/,&
'@                                                            ',/,&
'@                   ISTTIO    NSTITS    Turbulence           ',/,&
'@       AMONT : ',2I10,A13                                    ,/,&
'@       ACTUEL: ',2I10,A13                                    ,/,&
'@                                                            ',/,&
'@    Le calcul se poursuit...                                ',/,&
'@                                                            ',/,&
'@    Il est conseille de verifier ces indicateurs dans       ',/,&
'@      le sous-programme USLAG1 et USINI1                    ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 9340 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,&
'@    =========     LAGRANGIEN ',A13                           ,/,&
'@      DONNEES AMONT ET ACTUELLES INCOHERENTES               ',/,&
'@                                                            ',/,&
'@    LE CALCUL ACTUEL DES TERMES SOURCES DE COUPLAGE RETOUR  ',/,&
'@      EST EN MODE STATIONNAIRE, ALORS QUE LE FICHIER        ',/,&
'@      SUITE CONTIENT DES TERMES SOURCES INSTATIONNAIRES.    ',/,&
'@                                                            ',/,&
'@    NSTITS devrait etre un entier superieur ou egal         ',/,&
'@      a l''iteration Lagrangienne absolue de redemarrage    ',/,&
'@      du calcul (iteration : ',I10   ,')                    ',/,&
'@                                                            ',/,&
'@      Il vaut ici NSTITS = ',I10                             ,/,&
'@                                                            ',/,&
'@    Le calcul ne peut pas etre execute.                     ',/,&
'@                                                            ',/,&
'@    Verifier ces indicateurs dans le sous-programme USLAG1. ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 9341 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : A LA LECTURE DU FICHIER SUITE               ',/,&
'@    =========     LAGRANGIEN ',A13                           ,/,&
'@      DONNEES AMONT ET ACTUELLES DIFFERENTES                ',/,&
'@                                                            ',/,&
'@    LE CALCUL ACTUEL DES TERMES SOURCES DE COUPLAGE RETOUR  ',/,&
'@      EST EN MODE INSTATIONNAIRE, ALORS QUE LE FICHIER      ',/,&
'@      SUITE CONTIENT DES TERMES SOURCES STATIONNAIRES.      ',/,&
'@                                                            ',/,&
'@    Les termes sources de couplage retour stationnaires     ',/,&
'@      amont seront remises a zero.                          ',/,&
'@                                                            ',/,&
'@    Le calcul se poursuit...                                ',/,&
'@                                                            ',/,&
'@    Il est conseille de verifier ces indicateurs dans       ',/,&
'@      le sous-programme USLAG1.                             ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 9342 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,&
'@    =========     LAGRANGIEN ',A13                           ,/,&
'@      DONNEES AMONT ET ACTUELLES INCOHERENTES               ',/,&
'@                                                            ',/,&
'@    LE CALCUL SE POURSUIT AVEC UN CALCUL DES TERMES         ',/,&
'@      SOURCES DE COUPLAGE RETOUR EN MODE STATIONNAIRE       ',/,&
'@      MAIS LES INDICATEURS DE CONTROLES DES TERMES SOURCES  ',/,&
'@      ON ETE MODIFIEES.                                     ',/,&
'@                                                            ',/,&
'@    Pour eviter les incoherences dans le calcul             ',/,&
'@      NSTITS ne devrait pas etre modifie entre deux calculs ',/,&
'@      de termes sources de couplage retour stationnaires.   ',/,&
'@                                                            ',/,&
'@    Le calcul ne peut pas etre execute.                     ',/,&
'@                                                            ',/,&
'@    Verifier ces indicateurs dans le sous-programme USLAG1. ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

!----
! FIN
!----

end subroutine
