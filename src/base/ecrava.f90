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

subroutine ecrava &
!================

 ( ndim   , ncelet , ncel   , nfabor ,                            &
   nvar   , nscal  ,                                              &
   xyzcen , cdgfbo ,                                              &
   dt     , rtp    , propce ,                                     &
   frcxt  , prhyd  )

!===============================================================================

! FONCTION :
! ----------
! ECRITURE D'UN FICHIER SUITE

! PAS DE STOP (RETURN SI ON ECHOUE)

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! ncelet           ! i  ! <-- ! number of extended (real + ghost) cells        !
! ncel             ! i  ! <-- ! number of cells                                !
! nfabor           ! i  ! <-- ! number of boundary faces                       !
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp              ! tr ! <-- ! variables de calcul au centre des              !
! (ncelet,*)       !    !     !    cellules (instant courant)                  !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! frcxt(3,ncelet)  ! tr ! <-- ! force exterieure generant la pression          !
!                  !    !     !  hydrostatique                                 !
! prhyd(ncelet)    ! tr ! <-- ! hydrostatic pressure predicted                 !
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
use cstphy
use entsor
use pointe
use optcal
use albase
use alstru
use alaste
use ppppar
use ppthch
use ppincl
use coincl
use cpincl
use cs_fuel_incl
use elincl
use ppcpfu
use cplsat
use field
use atincl, only: init_at_chem
use atchem, only: ichemistry
use siream, only: iaerosol
use mesh, only: isympa
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer          ndim   , ncelet , ncel   , nfabor
integer          nvar   , nscal

double precision xyzcen(ndim,ncelet)
double precision cdgfbo(ndim,nfabor)
double precision dt(ncelet), rtp(ncelet,nflown:nvar)
double precision propce(ncelet,*)
double precision frcxt(3,ncelet), prhyd(ncelet)

! Local variables

integer          nbmom2
parameter       (nbmom2=nbmomx*2)

character*80     fname
character        rubriq*64,car2*2,car4*4,car54*54
character        cindfp*2,cindfs*4,cindff*4,cindfm*4
character        cindfc*2,cindfl*4
character        cphase*2 , cscal(nscamx)*4
character        cflu  (nvarmx)*4
character        nomflu(nvarmx)*18, nomrtp(nvarmx)*20
character        nomcli(nvarmx)*18
character        cstruc(nstrmx)*2, cindst*2
character        ficsui*32
logical          lprev
integer          nphas
integer          ivar  , iscal , f_id
integer          idecal, iclapc, icha  , icla
integer          ii    , ivers
integer          iptsna, iptsta, iptsca
integer          ierror, irtyp , itysup, nbval
integer          ipcefj, ipcla1, ipcla2, ipcla3
integer          nfmtsc, nfmtfl, nfmtmo, nfmtch, nfmtcl
integer          nfmtst
integer          nbflu , ilecec, iecr
integer          ngbstr(2)
integer          ifac, iel, istr, isou
integer          impava, impavx, nfld, iflmas, iflmab
double precision tmpstr(27)

integer, allocatable, dimension(:) :: mflnum
double precision, allocatable, dimension(:) :: w1
double precision, allocatable, dimension(:,:) :: l_velocity

double precision, dimension(:,:), pointer :: xut
double precision, dimension(:), pointer :: sval
double precision, dimension(:,:), pointer :: vel
double precision, dimension(:,:), pointer :: val_vp

!===============================================================================
!     A noter :
!        Lorsque qu'il est necessaire d'utiliser un ordre implicite
!        de rangement des variables, on a choisi :
!          P,
!          (U, V, W, turbulence, scalaires)_phase1,
!          (...)_phase2,
!          (...)_...
!          scalaires

!          avec turbulence = k, epsilon
!                       ou   R11, R22, R33, R12, R13, R23, epsilon
!                       ou   k, epsilon, phi, f_barre
!                       ou   k, omega

!        Ceci est par exemple utilise pour relier les flux de masse aux
!        variables


!===============================================================================
! 1. INITIALISATION
!===============================================================================


!===============================================================================
! 1. VERIFICATIONS DE BASE ET CODAGE DES CHAINES DE CARACTERES
!===============================================================================

!  --->  On code en chaine le numero des phases et scalaires
!        ----------------------------------------------------

!     Nombre de scalaires, de flux, de moments et de charbons
!       max pour les formats choisis
nfmtsc = 9999
nfmtfl = 9999
nfmtmo = 9999
nfmtch = 99
nfmtcl = 9999

!     Indefini (on met qqch de different de lecamo (pour generer une
!       erreur a la lecture)
CINDFP = 'XX'
CINDFS = 'XXXX'
CINDFF = 'XXXX'
CINDFM = 'XXXX'
CINDFC = 'XX'
CINDFL = 'XXXX'


!     Codage en chaine de caracteres du numero de la phase
write(cphase,'(I2.2)') 1

!     Codage en chaine de caracteres du numero du scalaire
do iscal = 1, min(nscal ,nfmtsc)
  write(cscal(iscal),'(I4.4)') iscal
enddo
do iscal = min(nscal ,nfmtsc)+1,nscal
  cscal(iscal) = cindfs
enddo

!     Codage en chaine de caracteres du numero du flux de masse
do ivar = 1, min(nvar  ,nfmtfl)
  write(cflu(ivar),'(I4.4)') ivar
enddo
do ivar = min(nvar  ,nfmtfl)+1,nvar
  cflu(ivar) = cindff
enddo

!     Verifications pour les formats et les numeros
!       de scalaire en chaine.
!     Avertissement (certaines infos passent a la trappe)
if (nscamx.gt.nfmtsc) then
  write(nfecra,7001)nfmtsc,nscamx
endif
if (nvarmx.gt.nfmtfl) then
  write(nfecra,7002)nfmtfl,nvarmx
endif
if (nbmomx.gt.nfmtmo) then
  write(nfecra,7003)nfmtmo,nbmomx
endif
if (ncharm.gt.nfmtch) then
  write(nfecra,7004)nfmtch,ncharm
endif
if (ncpcmx.gt.nfmtcl) then
  write(nfecra,7005)nfmtcl,ncpcmx
endif


!===============================================================================
! 2. OUVERTURE FICHIER SUITE DE BASE
!===============================================================================
! ILECEC = 2 : ecriture

write(nfecra,1000)

ilecec = 2

ficsui = 'main'
call opnsui(ficsui, len(ficsui), ilecec, impava, ierror)
!==========
if (ierror.ne.0) then
  write(nfecra,8000) ficsui
  return
endif


!===============================================================================
! 3. ECRITURE FICHIER SUITE DE BASE
!===============================================================================

write(nfecra,1100)


! 3.0 VERSION  : Rubrique "fichier suite ppal"
!=============   Pourrait porter le numero de version si besoin.
!                On ne se sert pas de IVERS (=1.2.0) pour le moment.

ivers  = 120
itysup = 0
nbval  = 1
irtyp  = 1
rubriq = 'version_fichier_suite_principal'
call ecrsui(impava,rubriq,len(rubriq),itysup,nbval,irtyp,ivers)


! 3.1 DIMENSIONS : les dimensions geometriques sont ecrites
!===============   automatiquement lors de l'ouverture du fichier
!                  on ecrit ici les nombres de variables

! The checkpoint file is now single-phase by default
nphas = 1

itysup = 0
nbval  = 1
irtyp  = 1

rubriq = 'nombre_variables'
call ecrsui(impava,rubriq,len(rubriq),itysup,nbval,irtyp,nvar)

rubriq = 'nombre_scalaires'
call ecrsui(impava,rubriq,len(rubriq),itysup,nbval,irtyp,nscal)

rubriq = 'nombre_scalaires_us'
call ecrsui(impava,rubriq,len(rubriq),itysup,nbval,irtyp,nscaus)

rubriq = 'nombre_scalaires_pp'
call ecrsui(impava,rubriq,len(rubriq),itysup,nbval,irtyp,nscapp)

rubriq = 'nombre_phases'
call ecrsui(impava,rubriq,len(rubriq),itysup,nbval,irtyp,nphas)

#if defined(_CS_LANG_FR)
car54 =' Fin de l''ecriture des dimensions                    '
#else
car54 =' End writing the dimensions                           '
#endif
write(nfecra,1110) car54

! 3.2 OPTIONS (Celles servant a donner le nombre de tableaux a lire)
!============================================================================
! Remarque : ces variables ne sont pas toutes utilies pour la relecture
!            en revanche, elles peuvent completer les infos au sein
!            du fichier suite

!  ---> Nombre de pas de temps, instant precedent
rubriq = 'nbre_pas_de_temps'
itysup = 0
nbval  = 1
irtyp  = 1
call ecrsui(impava,rubriq,len(rubriq),itysup,nbval,irtyp,ntcabs)

rubriq = 'instant_precedent'
itysup = 0
nbval  = 1
irtyp  = 2
call ecrsui(impava,rubriq,len(rubriq),itysup,nbval,irtyp,ttcabs)

!  ---> Modeles de turbulence
rubriq = 'modele_turbulence_phase'//cphase
itysup = 0
nbval  = 1
irtyp  = 1
call ecrsui(impava,rubriq,len(rubriq),itysup,nbval,irtyp,iturb)

!  ---> Methode ALE
rubriq = 'methode_ALE'
itysup = 0
nbval  = 1
irtyp  = 1
call ecrsui(impava,rubriq,len(rubriq),itysup,nbval,irtyp,iale)

!  ---> Cavitation
rubriq = 'cavitation'
itysup = 0
nbval  = 1
irtyp  = 1
call ecrsui(impava,rubriq,len(rubriq),itysup,nbval,irtyp,icavit)

rubriq = 'instant_mobile_precedent'
itysup = 0
nbval  = 1
irtyp  = 2
call ecrsui(impava,rubriq,len(rubriq),itysup,nbval,irtyp,ttcmob)


#if defined(_CS_LANG_FR)
car54 =' Fin de l''ecriture des options                       '
#else
car54 =' End writing the options                              '
#endif
write(nfecra,1110) car54

! 3.3 VARIABLES "PRINCIPALES"
!====================================

nomrtp(ipr)='pression_ce_phase'//cphase
nomrtp(iu)='vitesse_u_ce_phase'//cphase
nomrtp(iv)='vitesse_v_ce_phase'//cphase
nomrtp(iw)='vitesse_w_ce_phase'//cphase
if (itytur == 2) then
  nomrtp(ik)='k_ce_phase'//cphase
  nomrtp(iep)='eps_ce_phase'//cphase
elseif (itytur == 3) then
  nomrtp(ir11)='R11_ce_phase'//cphase
  nomrtp(ir22)='R22_ce_phase'//cphase
  nomrtp(ir33)='R33_ce_phase'//cphase
  nomrtp(ir12)='R12_ce_phase'//cphase
  nomrtp(ir13)='R13_ce_phase'//cphase
  nomrtp(ir23)='R23_ce_phase'//cphase
  nomrtp(iep)='eps_ce_phase'//cphase
  if (iturb.eq.32) then
    nomrtp(ial)='alp_ce_phase'//cphase
  endif
elseif (itytur == 5) then
  nomrtp(ik)='k_ce_phase'//cphase
  nomrtp(iep)='eps_ce_phase'//cphase
  nomrtp(iphi)='phi_ce_phase'//cphase
  if (iturb.eq.50) then
    nomrtp(ifb)='fb_ce_phase'//cphase
  elseif (iturb.eq.51) then
    nomrtp(ial)='al_ce_phase'//cphase
  endif
elseif (iturb == 60) then
  nomrtp(ik)='k_ce_phase'//cphase
  nomrtp(iomg)='omega_ce_phase'//cphase
elseif (iturb.eq.70) then
  nomrtp(inusa)='nusa_ce_phase'//cphase
endif
if (nscal.gt.0) then
  do iscal = 1, nscal
    !  ---> Turbulent flux model
    rubriq = 'turbulent_flux_model'//cscal(iscal)
    itysup = 0
    nbval  = 1
    irtyp  = 1
    call ecrsui(impava,rubriq,len(rubriq),itysup,nbval,irtyp,iturt(iscal))
    nomrtp(isca(iscal))='scalaire_ce_'//cscal(iscal)
  enddo
endif
if (iale.eq.1) then
  nomrtp(iuma)='vit_maillage_u_ce'
  nomrtp(ivma)='vit_maillage_v_ce'
  nomrtp(iwma)='vit_maillage_w_ce'
endif
if (icavit.ge.0) then
  nomrtp(ivoidf)='taux_vide_ce'
endif

!     Dans le cas ou il y a plusieurs phases,
!       on ne veut ecrire la pression qu'une seule fois
!       mais ca tombe bien, car on ne la verra qu'une seule fois
!       dans la liste des variables

! TODO: improve ecrsui for field structure
call field_get_val_v(ivarfl(iu), vel)
allocate(l_velocity(ncelet, 3))
do iel = 1, ncelet
  do isou = 1, 3
    l_velocity(iel,isou) = vel(isou,iel)
  enddo
enddo

do ivar = 1, 3

  itysup = 1
  nbval  = 1
  irtyp  = 2
  rubriq = nomrtp(ivar)
  call ecrsui(impava,rubriq,len(rubriq),itysup,nbval,irtyp,l_velocity(1,ivar))

enddo
deallocate(l_velocity)

do ivar = nflown, nvar

  itysup = 1
  nbval  = 1
  irtyp  = 2
  rubriq = nomrtp(ivar)
  call ecrsui(impava,rubriq,len(rubriq),itysup,nbval,irtyp,rtp(1,ivar))

enddo


do iscal = 1, nscal
  if (ityturt(iscal).eq.2 .or. ityturt(iscal).eq.3) then
    !  ---> Turbulent flux model
    ivar = isca(iscal)
    call field_get_name(ivarfl(ivar), fname)
    ! Index of the corresponding turbulent flux
    call field_get_id(trim(fname)//'_turbulent_flux', f_id)
    call field_get_val_v(f_id, xut)
    rubriq = trim(fname)//'_turbulent_flux_ce'
    itysup = 1
    nbval  = 3
    irtyp  = 2
    call ecrsui(impava,rubriq,len(rubriq),itysup,nbval,irtyp,xut)
  endif
enddo

if (ichemistry.gt.0.or.iaerosol.gt.0) then
  rubriq = 'atmospheric_chem'
  itysup = 0
  nbval  = 1
  irtyp  = 1
  call ecrsui(impava,rubriq,len(rubriq),itysup,nbval,irtyp,init_at_chem)
endif

#if defined(_CS_LANG_FR)
car54 =' Fin de l''ecriture des variables principales         '
#else
car54 =' End writing the main variables                       '
#endif
write(nfecra,1110)car54


!===============================================================================
! 4. FERMETURE FICHIER SUITE DE BASE
!===============================================================================

!     Fermeture du fichier suite principal
call clssui(impava,ierror)

if (ierror.ne.0) then
  write(nfecra,8010) ficsui
endif

write(nfecra,1200)

!===============================================================================
! 5. ECRITURE FICHIER SUITE AUXILIAIRE
!===============================================================================

!     Si  l'ecriture du fichier suite auxiliaire est demandee
if (iecaux.eq.1) then

! 5.0. OUVERTURE FICHIER SUITE AUXILIAIRE
!================================================

  write(nfecra,2000)

  ilecec = 2
  ficsui = 'auxiliary'
  call opnsui(ficsui, len(ficsui), ilecec, impavx, ierror)
  !==========
  if (ierror.ne.0) then
    write(nfecra,8001) ficsui
    return
  endif

  write(nfecra,1100)

! 5.1 VERSION  : Rubrique "fichier suite auxiliaire"
!=============   Pourrait porter le numero de version si besoin.
!                On ne se sert pas de IVERS (=1.2.0) pour le moment.

  ivers  = 120
  itysup = 0
  nbval  = 1
  irtyp  = 1
  rubriq = 'version_fichier_suite_auxiliaire'
  call ecrsui(impavx,rubriq,len(rubriq),itysup,nbval,irtyp,ivers)

! 5.2 DIMENSIONS : les dimensions geometriques sont ecrites
!===============   automatiquement lors de l'ouverture du fichier


!  ---> Nombre de phases
!       On les reecrit ici car on en aura besoin a la relecture
  rubriq = 'nombre_phases'
  itysup = 0
  nbval  = 1
  irtyp  = 1
  call ecrsui(impavx,rubriq,len(rubriq),itysup,nbval,irtyp,nphas)

!  ---> Nombre de pas de temps, instant precedent
!       On les reecrit ici car on en aura besoin a la relecture
  rubriq = 'nbre_pas_de_temps'
  itysup = 0
  nbval  = 1
  irtyp  = 1
  call ecrsui(impavx,rubriq,len(rubriq),itysup,nbval,irtyp,ntcabs)

!  ---> Indicateur de pas de temps variable
  rubriq = 'indic_dt_variable'
  itysup = 0
  nbval  = 1
  irtyp  = 1
  call ecrsui(impavx,rubriq,len(rubriq),itysup,nbval,irtyp,idtvar)

!  ---> Modeles de turbulence
!       On les reecrit ici car on en aura besoin a la relecture
  rubriq = 'modele_turbulence_phase'//cphase
  itysup = 0
  nbval  = 1
  irtyp  = 1
  call ecrsui(impavx,rubriq,len(rubriq),itysup,nbval,irtyp,iturb)

  rubriq = 'methode_ALE'
  itysup = 0
  nbval  = 1
  irtyp  = 1
  call ecrsui(impavx,rubriq,len(rubriq),itysup,nbval,irtyp,iale)

  !  ---> Cavitation
  rubriq = 'cavitation'
  itysup = 0
  nbval  = 1
  irtyp  = 1
  call ecrsui(impavx,rubriq,len(rubriq),itysup,nbval,irtyp,icavit)


#if defined(_CS_LANG_FR)
  car54 =' Fin de l''ecriture des dimensions et des options     '
#else
  car54 =' End writing the dimensions and options               '
#endif
  write(nfecra,1110)car54

! 5.3 ECRITURE DES VARIABLES
!===================================

! --->  Proprietes physiques

  !     Point de reference pour la pression totale
  !     On n'ecrit que si XYZP0 a ete specifie par l'utilisateur ou
  !       calcule a partir de faces de sorties ou de Dirichlet
  if (ixyzp0.eq.1) then
    rubriq = 'ref_presstot'//cphase
    itysup = 0
    nbval  = 3
    irtyp  = 2
    call ecrsui(impavx,rubriq,len(rubriq),itysup,nbval,irtyp,xyzp0(1))
  endif

  ! The physical variables herebelow are required for the low-Mach algorithm

  if (idilat.eq.3) then

    !the reference density updated with the low-Mach algorithm
    rubriq = 'ro0'//cphase
    itysup = 0
    nbval  = 1
    irtyp  = 2
    call ecrsui(impavx,rubriq,len(rubriq),itysup,nbval,irtyp,ro0)

    ! the thermodynamic pressure for the previous time step
    rubriq = 'pther'//cphase
    itysup = 0
    nbval  = 1
    irtyp  = 2
    call ecrsui(impavx,rubriq,len(rubriq),itysup,nbval,irtyp,pther)
  endif

  !     Masse volumique si elle est variable uniquement
  !     La masse volumique est egalement ecrite pour le modele de cavitation
  if (irovar.eq.1.or.icavit.ge.0) then
    !          Masse volumique - cellules
    rubriq = 'rho_ce_phase'//cphase
    itysup = 1
    nbval  = 1
    irtyp  = 2
    call field_get_val_s(icrom, sval)
    call ecrsui(impavx,rubriq,len(rubriq),itysup,nbval,irtyp,sval)

    ! Masse volumique du pdt precedent - cellules (uniquement pour cavitation)
    if (icavit.ge.0) then
      rubriq = 'rho_old_ce_phase'//cphase
      itysup = 1
      nbval  = 1
      irtyp  = 2
      call field_get_val_prev_s(icrom, sval)
      call ecrsui(impavx,rubriq,len(rubriq),itysup,nbval,irtyp,sval)
    endif

    !          Masse volumique - faces de bord
    rubriq = 'rho_fb_phase'//cphase
    itysup = 3
    nbval  = 1
    irtyp  = 2
    call field_get_val_s(ibrom, sval)
    call ecrsui(impavx,rubriq,len(rubriq),itysup,nbval,irtyp,sval)
  endif

  !     On n'ecrit les proprietes physiques que si on les extrapole ou
  !     pour le modele de cavitation
  !       On pourrait les ecrire a tous les coups en prevision d'une
  !       suite avec extrapolation, mais
  !          - c'est rare
  !          - si on demarre un calcul a l'ordre deux a partir d'un calcul
  !            a l'ordre 1, on peut estimer que les premiers pas de temps
  !            sont a jeter de toute facon.
  !       Une exception : on ecrit egalement Cp en effet joule pour
  !         pouvoir calculer la temperature H/Cp en debut de calcul

  if (iviext.gt.0.or.icavit.ge.0) then
    !         Viscosite moleculaire - cellules (si variable ou cavitation)
    if (ivivar.eq.1.or.icavit.ge.0) then
      rubriq = 'viscl_ce_phase'//cphase
      itysup = 1
      nbval  = 1
      irtyp  = 2
      call ecrsui(impavx,rubriq,len(rubriq),itysup,nbval,irtyp, &
                  propce(1,ipproc(iviscl)))
    endif

    if (iviext.gt.0) then
      !         Viscosite turbulente ou de sous-maille - cellules
      rubriq = 'visct_ce_phase'//cphase
      itysup = 1
      nbval  = 1
      irtyp  = 2
      call ecrsui(impavx,rubriq,len(rubriq),itysup,nbval,irtyp,   &
                  propce(1,ipproc(ivisct)))
    endif
  endif

  if ((icpext.gt.0.and.icp.gt.0).or.              &
       (ippmod(ieljou).ge.1.and.icp.gt.0))  then
    !         Chaleur massique - cellules
    rubriq = 'cp_ce_phase'//cphase
    itysup = 1
    nbval  = 1
    irtyp  = 2
    call ecrsui(impavx,rubriq,len(rubriq),itysup,nbval,irtyp,   &
                propce(1,ipproc(icp)))
  endif

!     Si on a des scalaires, on ecrit leur model de flux et
!     leur diffusivite (on ne l'ecrit pas pour les variances)
  if (nscal.gt.0) then
    do iscal = 1, nscal
      if (ivsext(iscal).gt.0.and.ivisls(iscal).gt.0.and.           &
         (iscavr(iscal).le.0.or.iscavr(iscal).gt.nscal) ) then
        ! Cell diffusivity
        rubriq = 'visls_ce_scalaire'//cscal(iscal)
        itysup = 1
        nbval  = 1
        irtyp  = 2
        call ecrsui(impavx,rubriq,len(rubriq),itysup,nbval,irtyp, &
                    propce(1,ipproc(ivisls(iscal))))
      endif
    enddo
  endif

#if defined(_CS_LANG_FR)
  car54 =' Fin de l''ecriture des proprietes physiques          '
#else
  car54 =' End writing the physical properties                  '
#endif
  write(nfecra,1110)car54

! ---> Pas de temps

  if (idtvar.eq.2) then
    rubriq = 'dt_variable_espace_ce'
    itysup = 1
    nbval  = 1
    irtyp  = 2
    call ecrsui(impavx,rubriq,len(rubriq),itysup,nbval,irtyp,dt)
  elseif (idtvar.eq.1) then
    rubriq = 'dt_variable_temps'
    itysup = 0
    nbval  = 1
    irtyp  = 2
    call ecrsui(impavx,rubriq,len(rubriq),itysup,nbval,irtyp,dt)
  endif

#if defined(_CS_LANG_FR)
  car54 =' Fin de l''ecriture du pas de temps                   '
#else
  car54 =' End writing the time step                            '
#endif
  write(nfecra,1110)car54

! ---> Flux de masse

!     Pour garder la memoire de la correspondance entre les variables
!     et les flux de masse, on memorise le nom de chaque variable
!     (nomflu(i)= nom de la ieme variable)
!     Ensuite, pour chaque variable, on ecrit son nom et le numero
!     local du flux de masse correspondant (en pratique 1 ou 2)

!       Initialisation des tableaux de travail

  call field_get_n_fields(nfld)

  allocate(mflnum(nfld))

  nbflu = 0

  do ii = 1, nfld
    mflnum(ii) = 0
  enddo

  nomflu(ipr)='fm_p_phase'//cphase
  nomflu(iu)='fm_u_phase'//cphase
  nomflu(iv)='fm_v_phase'//cphase
  nomflu(iw)='fm_w_phase'//cphase
  if (itytur.eq.2) then
    nomflu(ik)='fm_k_phase'//cphase
    nomflu(iep)='fm_eps_phase'//cphase
  elseif (itytur.eq.3) then
    nomflu(ir11)='fm_R11_phase'//cphase
    nomflu(ir22)='fm_R22_phase'//cphase
    nomflu(ir33)='fm_R33_phase'//cphase
    nomflu(ir12)='fm_R12_phase'//cphase
    nomflu(ir13)='fm_R13_phase'//cphase
    nomflu(ir23)='fm_R23_phase'//cphase
    nomflu(iep)='fm_eps_phase'//cphase
    if (iturb.eq.32) then
      nomflu(ial)='fm_alp_phase'//cphase
    endif
  elseif (itytur.eq.5) then
    nomflu(ik)='fm_k_phase'//cphase
    nomflu(iep)='fm_eps_phase'//cphase
    nomflu(iphi)='fm_phi_phase'//cphase
    ! On n'utilise pas le flux de masse pour fb/al en fait mais on le laisse
    ! ici, car ca ne change rien (le flux n'est ecrit qu'une seule fois)
    if (iturb.eq.50) then
      nomflu(ifb)='fm_fb_phase'//cphase
    elseif (iturb.eq.51) then
      nomflu(ial)='fm_al_phase'//cphase
    endif
  elseif (iturb.eq.60) then
    nomflu(ik)='fm_k_phase'//cphase
    nomflu(iomg)='fm_omega_phase'//cphase
  elseif (iturb.eq.70) then
    nomflu(inusa)='fm_nusa_phase'//cphase
  endif
  if (nscal.gt.0) then
    do iscal = 1, nscal
      nomflu(isca(iscal))='fm_scalaire'//cscal(iscal)
    enddo
  endif
  if (iale.eq.1) then
    nomflu(iuma)='fm_vit_maill_u'
    nomflu(ivma)='fm_vit_maill_v'
    nomflu(iwma)='fm_vit_maill_w'
  endif
  if (icavit.ge.0) then
    nomflu(ivoidf)='fm_taux_vide'
  endif

  ! For variables

  do ivar = 1, nvar

    f_id = ivarfl(ivar)

    ! If the variable is not associated with a mass flux, do nothing

    call field_get_key_int(f_id, kimasf, iflmas) ! interior mass flux

    if (iflmas.ge.0) then

      if (mflnum(iflmas).eq.0) then

        ! Flux has not been written yet

        call field_get_val_s(iflmas, sval)

        nbflu=nbflu+1
        mflnum(iflmas) = nbflu

        ! Write mass flux at interior faces
        rubriq = 'flux_masse_fi_'//cflu(nbflu)
        itysup = 2
        nbval  = 1
        irtyp  = 2
        call ecrsui(impavx, rubriq, len(rubriq), itysup, nbval, irtyp, sval)

        call field_get_key_int(f_id, kbmasf, iflmab) ! boundary mass flux
        call field_get_val_s(iflmab, sval)

        ! Write mass flux at boundary faces
        rubriq = 'flux_masse_fb_'//cflu(nbflu)
        itysup = 3
        nbval  = 1
        irtyp  = 2
        call ecrsui(impavx, rubriq, len(rubriq), itysup, nbval, irtyp, sval)

      endif

      ! Whether flux has been written or not, associate variable with flux
      rubriq = nomflu(ivar)
      itysup = 0
      nbval  = 1
      irtyp  = 1
      call ecrsui(impavx, rubriq, len(rubriq), itysup, nbval, irtyp,   &
                  mflnum(iflmas))

    endif

  enddo

  ! Do the same for mass fluxes at previous time

  nbflu = 0

  do ii = 1, nfld
    mflnum(ii) = 0
  enddo

  nomflu(ipr)='fm_a_p_phase'//cphase
  nomflu(iu)='fm_a_u_phase'//cphase
  nomflu(iv)='fm_a_v_phase'//cphase
  nomflu(iw)='fm_a_w_phase'//cphase
  if (itytur.eq.2) then
    nomflu(ik)='fm_a_k_phase'//cphase
    nomflu(iep)='fm_a_eps_phase'//cphase
  elseif (itytur.eq.3) then
    nomflu(ir11)='fm_a_R11_phase'//cphase
    nomflu(ir22)='fm_a_R22_phase'//cphase
    nomflu(ir33)='fm_a_R33_phase'//cphase
    nomflu(ir12)='fm_a_R12_phase'//cphase
    nomflu(ir13)='fm_a_R13_phase'//cphase
    nomflu(ir23)='fm_a_R23_phase'//cphase
    nomflu(iep)='fm_a_eps_phase'//cphase
    if (iturb.eq.32) then
      nomflu(ial)='fm_a_alp_phase'//cphase
    endif
  elseif (itytur.eq.5) then
    nomflu(ik)='fm_a_k_phase'//cphase
    nomflu(iep)='fm_a_eps_phase'//cphase
    nomflu(iphi)='fm_a_phi_phase'//cphase
    if (iturb.eq.50) then
      nomflu(ifb)='fm_a_fb_phase'//cphase
    elseif (iturb.eq.51) then
      nomflu(ial)='fm_a_al_phase'//cphase
    endif
  elseif (iturb.eq.60) then
    nomflu(ik)='fm_a_k_phase'//cphase
    nomflu(iomg)='fm_a_omega_phase'//cphase
  elseif (iturb.eq.70) then
    nomflu(inusa)='fm_a_nusa_phase'//cphase
  endif
  if (nscal.gt.0) then
    do iscal = 1, nscal
      nomflu(isca(iscal))='fm_a_scalaire'//cscal(iscal)
    enddo
  endif
  if (iale.eq.1) then
    nomflu(iuma)='fm_a_vit_maill_u'
    nomflu(ivma)='fm_a_vit_maill_v'
    nomflu(iwma)='fm_a_vit_maill_w'
  endif
  if (icavit.ge.0) then
    nomflu(ivoidf)='fm_a_taux_vide'
  endif

  do ivar = 1, nvar

    f_id = ivarfl(ivar)

    ! If the variable is not associated with a mass flux, do nothing

    call field_get_key_int(f_id, kimasf, iflmas) ! interior mass flux

    if (mflnum(iflmas).eq.0) then

      ! Flux has not been written yet

      if (mflnum(iflmas).eq.0) then

        call field_have_previous(iflmas, lprev)

        if (.not. lprev) cycle ! skip to next loop variable

        ! Flux has not been written yet

        call field_get_val_prev_s(iflmas, sval)

        nbflu=nbflu+1
        mflnum(iflmas) = nbflu

        ! Write mass flux at interior faces
        rubriq = 'flux_masse_a_fi_'//cflu(nbflu)
        itysup = 2
        nbval  = 1
        irtyp  = 2
        call ecrsui(impavx, rubriq, len(rubriq), itysup, nbval, irtyp, sval)

        call field_get_key_int(f_id, kbmasf, iflmab) ! boundary mass flux
        call field_get_val_prev_s(iflmab, sval)

        ! Write mass flux at boundary faces
        rubriq = 'flux_masse_a_fb_'//cflu(nbflu)
        itysup = 3
        nbval  = 1
        irtyp  = 2
        call ecrsui(impavx, rubriq, len(rubriq), itysup, nbval, irtyp, sval)

      endif

      ! Whether flux has been written or not, associate variable with flux
      rubriq = nomflu(ivar)
      itysup = 0
      nbval  = 1
      irtyp  = 1
      call ecrsui(impavx, rubriq, len(rubriq), itysup, nbval, irtyp,   &
                  mflnum(iflmas))

    endif

  enddo

  deallocate(mflnum)

#if defined(_CS_LANG_FR)
  car54 =' Fin de l''ecriture des flux de masse                 '
#else
  car54 =' End writing the mass fluxes                          '
#endif
  write(nfecra,1110)car54

! ---> Conditions aux limites

  nomcli(IPR)='_p_phase'//cphase
  nomcli(IU)='_u_phase'//cphase
  nomcli(IV)='_v_phase'//cphase
  nomcli(IW)='_w_phase'//cphase
  if (itytur.eq.2) then
    nomcli(IK)='_k_phase'//cphase
    nomcli(IEP)='_eps_phase'//cphase
  elseif (itytur.eq.3) then
    nomcli(IR11)='_R11_phase'//cphase
    nomcli(IR22)='_R22_phase'//cphase
    nomcli(IR33)='_R33_phase'//cphase
    nomcli(IR12)='_R12_phase'//cphase
    nomcli(IR13)='_R13_phase'//cphase
    nomcli(IR23)='_R23_phase'//cphase
    nomcli(IEP)='_eps_phase'//cphase
    if (iturb.eq.32) then
      nomcli(ial)='_alp_phase'//cphase
    endif
  elseif (itytur.eq.5) then
    nomcli(IK)='_k_phase'//cphase
    nomcli(IEP)='_eps_phase'//cphase
    nomcli(IPHI)='_phi_phase'//cphase
    if (iturb.eq.50) then
      NOMCLI(IFB)='_fb_phase'//CPHASE
    elseif (iturb.eq.51) then
      nomcli(ial)='_al_phase'//cphase
    endif
  elseif (iturb.eq.60) then
    nomcli(ik)='_k_phase'//cphase
    nomcli(iomg)='_omega_phase'//cphase
  elseif (iturb.eq.70) then
    nomcli(inusa)='_nusa_phase'//cphase
  endif
  if (nscal.gt.0) then
    do iscal = 1, nscal
      nomcli(isca(iscal))='_scalaire'//cscal(iscal)
    enddo
  endif
  if (iale.eq.1) then
    nomcli(iuma)='_vit_maillage_u'
    nomcli(ivma)='_vit_maillage_v'
    nomcli(iwma)='_vit_maillage_w'
  endif
  if (icavit.ge.0) then
    nomcli(ivoidf)='_taux_vide'
  endif

  call restart_write_bc_coeffs(impavx)

  ! Symmetry flag (used for least-squares gradients,
  ! with extrapolation at boundary).

  rubriq = 'isympa_fb_phase01'
  itysup = 3
  nbval  = 1
  irtyp  = 1
  call ecrsui(impavx,rubriq,len(rubriq),itysup,nbval,irtyp,isympa)

#if defined(_CS_LANG_FR)
  car54 =' Fin de l''ecriture des conditions aux limites        '
#else
  car54 =' End writing the boundary conditions                  '
#endif
  write(nfecra,1110)car54

  ! Backward differential scheme in time
  if (ibdtso.gt.1) then

    ! Warning: must be adapted if ivar.ne.iu
    !          must be adapted if ibdtso.gt.2
    ivar = iu
    f_id = ivarfl(ivar)
    call field_get_val_prev_v(f_id, val_vp)
    rubriq = 'velocity_prev'
    itysup = 1 ! cells location
    nbval = 3 ! interleaved velocity
    irtyp  = 2 ! double precision
    call ecrsui(impavx,rubriq,len(rubriq),itysup,nbval,irtyp,   &
                val_vp)

  endif

! ---> Termes sources
!      Lorsqu'ils sont extrapoles (pour les versions elec, voir plus bas)

  iecr = 0

! ---> Termes sources Navier-Stokes (plus taux de vide si cavitation)

  !     Si les termes sont a l'ordre 2
  if (isno2t.gt.0) then

    iecr = 1

    iptsna = ipproc(itsnsa)
    itysup = 1
    nbval  = 1
    irtyp  = 2

    rubriq = 'tsource_ns_ce_x_phase'//cphase
    call ecrsui(impavx,rubriq,len(rubriq),itysup,nbval,irtyp,   &
                propce(1,iptsna))

    rubriq = 'tsource_ns_ce_y_phase'//cphase
    call ecrsui(impavx,rubriq,len(rubriq),itysup,nbval,irtyp,   &
                propce(1,iptsna+1))

    rubriq = 'tsource_ns_ce_z_phase'//cphase
    call ecrsui(impavx,rubriq,len(rubriq),itysup,nbval,irtyp,   &
                propce(1,iptsna+2))

    if (icavit.ge.0) then
      rubriq = 'tsource_vf_ce_phase'//cphase
      call ecrsui(impavx,rubriq,len(rubriq),itysup,nbval,irtyp,   &
                propce(1,iptsna+3))
    endif

  endif

! ---> Termes sources turbulence

  !        Si les termes sont a l'ordre 2
  if (isto2t.gt.0) then

    iecr = 1

    iptsta = ipproc(itstua)
    itysup = 1
    nbval  = 1
    irtyp  = 2

    !          En k-eps
    if (itytur.eq.2) then

      rubriq = 'tsource_tu_ce_k_phase'//cphase
      call ecrsui(impavx,rubriq,len(rubriq),itysup,nbval,irtyp, &
                  propce(1,iptsta))

      rubriq = 'tsource_tu_ce_eps_phase'//cphase
      call ecrsui(impavx,rubriq,len(rubriq),itysup,nbval,irtyp, &
                  propce(1,iptsta+1))

      !          En Rij
    elseif (itytur.eq.3) then

      rubriq = 'tsource_tu_ce_R11_phase'//cphase
      call ecrsui(impavx,rubriq,len(rubriq),itysup,nbval,irtyp, &
                  propce(1,iptsta))
      rubriq = 'tsource_tu_ce_R22_phase'//cphase
      call ecrsui(impavx,rubriq,len(rubriq),itysup,nbval,irtyp, &
                  propce(1,iptsta+1))
      rubriq = 'tsource_tu_ce_R33_phase'//cphase
      call ecrsui(impavx,rubriq,len(rubriq),itysup,nbval,irtyp, &
                  propce(1,iptsta+2))
      rubriq = 'tsource_tu_ce_R12_phase'//cphase
      call ecrsui(impavx,rubriq,len(rubriq),itysup,nbval,irtyp, &
                  propce(1,iptsta+3))
      rubriq = 'tsource_tu_ce_R13_phase'//cphase
      call ecrsui(impavx,rubriq,len(rubriq),itysup,nbval,irtyp, &
                  propce(1,iptsta+4))
      rubriq = 'tsource_tu_ce_R23_phase'//cphase
      call ecrsui(impavx,rubriq,len(rubriq),itysup,nbval,irtyp, &
                  propce(1,iptsta+5))

      rubriq = 'tsource_tu_ce_eps_phase'//cphase
      call ecrsui(impavx,rubriq,len(rubriq),itysup,nbval,irtyp, &
                  propce(1,iptsta+6))

      if (iturb.eq.32) then
        rubriq = 'tsource_tu_ce_alp_phase'//cphase
        call ecrsui(impavx,rubriq,len(rubriq),itysup,nbval,irtyp, &
                    propce(1,iptsta+7))
      endif

      !          En v2f
    elseif (itytur.eq.5) then

      rubriq = 'tsource_tu_ce_k_phase'//cphase
      call ecrsui(impavx,rubriq,len(rubriq),itysup,nbval,irtyp, &
                  propce(1,iptsta))

      rubriq = 'tsource_tu_ce_eps_phase'//cphase
      call ecrsui(impavx,rubriq,len(rubriq),itysup,nbval,irtyp, &
           propce(1,iptsta+1))

      rubriq = 'tsource_tu_ce_phi_phase'//cphase
      call ecrsui(impavx,rubriq,len(rubriq),itysup,nbval,irtyp, &
                  propce(1,iptsta+2))

      if (iturb.eq.50) then
        rubriq = 'tsource_tu_ce_fb_phase'//cphase
        call ecrsui(impavx,rubriq,len(rubriq),itysup,nbval,irtyp, &
                    propce(1,iptsta+3))
      elseif (iturb.eq.51) then
        rubriq = 'tsource_tu_ce_al_phase'//cphase
        call ecrsui(impavx,rubriq,len(rubriq),itysup,nbval,irtyp, &
                    propce(1,iptsta+3))
      endif

      !          En k-omega
    elseif (iturb.eq.60) then

      rubriq = 'tsource_tu_ce_k_phase'//cphase
      call ecrsui(impavx,rubriq,len(rubriq),itysup,nbval,irtyp, &
                  propce(1,iptsta))

      rubriq = 'tsource_tu_ce_omega_phase'//cphase
      call ecrsui(impavx,rubriq,len(rubriq),itysup,nbval,irtyp, &
                  propce(1,iptsta+1))
      !          En Spalart Allmaras
    elseif (iturb.eq.70) then

      rubriq = 'tsource_tu_ce_nusa_phase'//cphase
      call ecrsui(impavx,rubriq,len(rubriq),itysup,nbval,irtyp, &
                  propce(1,iptsta))

    endif
  endif

! ---> Termes sources scalaires

!     Boucle sur les scalaires
  do iscal = 1, nscal
!     Si le terme est a l'ordre 2
    if (isso2t(iscal).gt.0) then

      iecr = 1

      iptsca = ipproc(itssca(iscal))
      itysup = 1
      nbval  = 1
      irtyp  = 2
      rubriq = 'tsource_sc_ce_scalaire'//cscal(iscal)
      call ecrsui(impavx,rubriq,len(rubriq),itysup,nbval,irtyp,   &
                  propce(1,iptsca))

    endif
  enddo


  if (iecr.ne.0) then
#if defined(_CS_LANG_FR)
    car54=' Fin de l''ecriture des termes sources                '
#else
    car54=' End writing the source terms                         '
#endif
    write(nfecra,1110)car54
  endif


! ---> Moyennes (cumuls)

  call time_moment_restart_write(impavx)

! ---> Distance a la paroi
!      On pourra ecrire ici la distance a la paroi

  iecr = 0

  if (ineedy.eq.1) then
!     Ancien mode de calcul. On ecrit aussi la distance a la paroi,
!       au cas ou on fait une suite en ICDPAR=1.
    if (abs(icdpar).eq.2) then
      iecr   = 1
      itysup = 1
      nbval  = 1
      irtyp  = 1
      rubriq = 'num_fac_par_ce_phase'//cphase
      call ecrsui(impavx,rubriq,len(rubriq),itysup,nbval,  irtyp,ifapat)
!     Pour la distance reelle, on a besoin d'un tableau provisoire
!     on ne prend que la phase 1
      allocate(w1(ncelet))
      do iel = 1, ncel
        ifac = ifapat(iel)
        w1(iel) = sqrt( (cdgfbo(1,ifac)-xyzcen(1,iel))**2 &
                      + (cdgfbo(2,ifac)-xyzcen(2,iel))**2 &
                      + (cdgfbo(3,ifac)-xyzcen(3,iel))**2 )
      enddo
      iecr   = 1
      itysup = 1
      nbval  = 1
      irtyp  = 2
      rubriq = 'dist_fac_par_ce_phase'//cphase
      call ecrsui(impavx,rubriq,len(rubriq),itysup,nbval,irtyp,w1)
      ! Free memory
      deallocate(w1)

!     Nouveau mode de calcul
    elseif (abs(icdpar).eq.1) then
      iecr   = 1
      itysup = 1
      nbval  = 1
      irtyp  = 2
      rubriq = 'dist_fac_par_ce_phase'//cphase
      call ecrsui(impavx,rubriq,len(rubriq),itysup,nbval,irtyp,dispar)
    endif
  endif

  if (iecr.ne.0) then
#if defined(_CS_LANG_FR)
    car54=' Fin de l''ecriture de la distance a la paroi         '
#else
    car54=' End writing the wall distance                        '
#endif
    write(nfecra,1110)car54
  endif

! ---> Force exterieure

  if (iphydr.eq.1) then

    itysup = 1
    nbval  = 3
    irtyp  = 2

    rubriq = 'force_ext_ce_phase'//cphase
    call ecrsui(impavx,rubriq,len(rubriq),itysup,nbval,irtyp,frcxt)

#if defined(_CS_LANG_FR)
    car54=' Fin de l''ecriture des forces exterieures            '
#else
    car54=' End writing the external forces                      '
#endif
    write(nfecra,1110)car54

  endif

! ---> Pression hydrostatique predite

  if (iphydr.eq.2) then

    itysup = 1
    nbval  = 1
    irtyp  = 2

    rubriq = 'Prhyd_pre_phase'//cphase
    call ecrsui(impavx,rubriq,len(rubriq),itysup,nbval,irtyp,prhyd(1))

#if defined(_CS_LANG_FR)
    car54=' Fin d''ecriture de la pression hydrostatique predite '
#else
    car54=' End writing the predicted hydrostatic pressure       '
#endif
    write(nfecra,1110)car54

  endif

! ---> Methode ALE

  if (iale.eq.1) then

    itysup = 4
    nbval  = 3
    irtyp  = 2

    rubriq = 'vertex_displacement'
    call ecrsui(impavx,rubriq,len(rubriq),itysup,nbval,irtyp,depale)

!     Viscosite de maillage (elle est souvent definie geometriquement sur le
!       maillage initial ... il est donc plus facile de la relire ensuite)

    rubriq = 'type_visc_mail'
    itysup = 0
    nbval  = 1
    irtyp  = 1
    call ecrsui(impavx,rubriq,len(rubriq),itysup,nbval,irtyp,iortvm)

    rubriq = 'visc_maillage_x'
    itysup = 1
    nbval  = 1
    irtyp  = 2
    call ecrsui(impavx,rubriq,len(rubriq),itysup,nbval,irtyp,     &
                propce(1,ipproc(ivisma(1))))

    if (iortvm.eq.1) then
      rubriq = 'visc_maillage_y'
      call ecrsui(impavx,rubriq,len(rubriq),itysup,nbval,irtyp,   &
                  propce(1,ipproc(ivisma(2))))
      rubriq = 'visc_maillage_z'
      call ecrsui(impavx,rubriq,len(rubriq),itysup,nbval,irtyp,   &
                  propce(1,ipproc(ivisma(3))))
    endif

#if defined(_CS_LANG_FR)
    car54=' Fin de l''ecriture des donnees ALE    '
#else
    car54=' End writing the ALE data              '
#endif
    write(nfecra,1110)car54

    ngbstr(1) = nbstru
    ngbstr(2) = nbaste

    rubriq = 'nombre_structures'
    itysup = 0
    nbval  = 2
    irtyp  = 1
    call ecrsui(impavx,rubriq,len(rubriq),itysup,nbval,irtyp,ngbstr)

    if (nbstru.gt.0) then

      nfmtst = 99
      cindst = 'XX'
!     Codage en chaine de caracteres du numero de la structure
      do istr = 1, min(nbstru ,nfmtst)
        write(cstruc(istr),'(I2.2)') istr
      enddo
      do istr = min(nbstru ,nfmtst)+1,nbstru
        cstruc(istr) = cindst
      enddo

      do istr = 1, nbstru
        rubriq = 'donnees_structure_'//cstruc(istr)
        itysup = 0
        nbval  = 27
        irtyp  = 2

        do ii = 1, 3
          tmpstr(   ii) = xstr  (ii,istr)
          tmpstr(3 +ii) = xpstr (ii,istr)
          tmpstr(6 +ii) = xppstr(ii,istr)
          tmpstr(9 +ii) = xsta  (ii,istr)
          tmpstr(12+ii) = xpsta (ii,istr)
          tmpstr(15+ii) = xppsta(ii,istr)
          tmpstr(18+ii) = xstp  (ii,istr)
          tmpstr(21+ii) = forstr(ii,istr)
          tmpstr(24+ii) = forsta(ii,istr)
        enddo

        call ecrsui(impavx,rubriq,len(rubriq),itysup,nbval,irtyp,tmpstr)
      enddo

#if defined(_CS_LANG_FR)
      car54=' Fin de l''ecriture des donnees des structures (ALE)'
#else
      car54=' End writing the structures data (ALE)              '
#endif
      write(nfecra,1110)car54

    endif
  endif


! ---> Grandeurs complementaires pour la combustion gaz

!     Modele COD3P :
!     ============

  if ( ippmod(icod3p).ge.0 ) then

    rubriq = 'hinfue_cod3p'
    itysup = 0
    nbval  = 1
    irtyp  = 2
    call ecrsui(impavx,rubriq,len(rubriq),itysup,nbval,irtyp,hinfue)

    rubriq = 'hinoxy_cod3p'
    itysup = 0
    nbval  = 1
    irtyp  = 2
    call ecrsui(impavx,rubriq,len(rubriq),itysup,nbval,irtyp,hinoxy)

    rubriq = 'tinfue_cod3p'
    itysup = 0
    nbval  = 1
    irtyp  = 2
    call ecrsui(impavx,rubriq,len(rubriq),itysup,nbval,irtyp,tinfue)

    rubriq = 'tinoxy_cod3p'
    itysup = 0
    nbval  = 1
    irtyp  = 2
    call ecrsui(impavx,rubriq,len(rubriq),itysup,nbval,irtyp,tinoxy)

!       Numero des zones
    itysup = 3
    nbval  = 1
    irtyp  = 1
    rubriq = 'num_zone_fb_cod3p'
    call ecrsui(impavx,rubriq,len(rubriq),itysup,nbval,irtyp,izfppp)

!       Entree Fuel (si ce n'est pas NOZPPM, erreur)
    itysup = 0
    nbval  = nozppm
    irtyp  = 1
    rubriq = 'ientfu_zone_bord_cod3p'
    call ecrsui(impavx,rubriq,len(rubriq),itysup,nbval,irtyp,ientfu)

!       Entree oxydant (si ce n'est pas NOZPPM, erreur)
    itysup = 0
    nbval  = nozppm
    irtyp  = 1
    rubriq = 'ientox_zone_bord_cod3p'
    call ecrsui(impavx,rubriq,len(rubriq),itysup,nbval,irtyp,ientox)

#if defined(_CS_LANG_FR)
    car54=' Fin de l''ecriture des informations combustion COD3P'
#else
    car54=' End writing combustion information (COD3P)         '
#endif
    write(nfecra,1110)car54

  endif

!      Modele EBU :
!      ==========

  if ( ippmod(icoebu).ge.0 ) then

    rubriq = 'temperature_gaz_frais_ebu'
    itysup = 0
    nbval  = 1
    irtyp  = 2
    call ecrsui(impavx,rubriq,len(rubriq),itysup,nbval,irtyp,tgf)

    rubriq = 'frmel_ebu'
    itysup = 0
    nbval  = 1
    irtyp  = 2
    call ecrsui(impavx,rubriq,len(rubriq),itysup,nbval,irtyp,frmel)

    ! Numero des zones
    itysup = 3
    nbval  = 1
    irtyp  = 1
    rubriq = 'num_zone_fb_ebu'
    call ecrsui(impavx,rubriq,len(rubriq),itysup,nbval,irtyp,izfppp)

!       Entree Gaz brule(si ce n'est pas NOZPPM, erreur)
    itysup = 0
    nbval  = nozppm
    irtyp  = 1
    rubriq = 'ientgb_zone_bord_ebu'
    call ecrsui(impavx,rubriq,len(rubriq),itysup,nbval,irtyp,ientgb)

!       Entree gaz frais (si ce n'est pas NOZPPM, erreur)
    itysup = 0
    nbval  = nozppm
    irtyp  = 1
    rubriq = 'ientgf_zone_bord_ebu'
    call ecrsui(impavx,rubriq,len(rubriq),itysup,nbval,irtyp,ientgf)

!       FMENT (si ce n'est pas NOZPPM, erreur)
    itysup = 0
    nbval  = nozppm
    irtyp  = 2
    rubriq = 'fment_zone_bord_ebu'
    call ecrsui(impavx,rubriq,len(rubriq),itysup,nbval,irtyp,fment)

!       TKENT (si ce n'est pas NOZPPM, erreur)
    itysup = 0
    nbval  = nozppm
    irtyp  = 2
    rubriq = 'tkent_zone_bord_ebu'
    call ecrsui(impavx,rubriq,len(rubriq),itysup,nbval,irtyp,tkent)

#if defined(_CS_LANG_FR)
    car54=' Fin de l''ecriture des informations combustion EBU '
#else
    car54=' End writing the combustion information (EBU)      '
#endif
    write(nfecra,1110)car54

  endif

!      Modele LWC :
!      ==========

  if ( ippmod(icolwc).ge.0 ) then

    rubriq = 'fmin_lwc'
    itysup = 0
    nbval  = 1
    irtyp  = 2
    call ecrsui(impavx,rubriq,len(rubriq),itysup,nbval,irtyp,fmin)

    rubriq = 'fmax_lwc'
    itysup = 0
    nbval  = 1
    irtyp  = 2
    call ecrsui(impavx,rubriq,len(rubriq),itysup,nbval,irtyp,fmax)

    rubriq = 'hmin_lwc'
    itysup = 0
    nbval  = 1
    irtyp  = 2
    call ecrsui(impavx,rubriq,len(rubriq),itysup,nbval,irtyp,hmin)

    rubriq = 'hmax_lwc'
    itysup = 0
    nbval  = 1
    irtyp  = 2
    call ecrsui(impavx,rubriq,len(rubriq),itysup,nbval,irtyp,hmax)

!       Numero des zones
    itysup = 3
    nbval  = 1
    irtyp  = 1
    rubriq = 'num_zone_fb_lwc'
    call ecrsui(impavx,rubriq,len(rubriq),itysup,nbval,irtyp,izfppp)

!       Entree Gaz brule(si ce n'est pas NOZPPM, erreur)
    itysup = 0
    nbval  = nozppm
    irtyp  = 1
    rubriq = 'ientgb_zone_bord_lwc'
    call ecrsui(impavx,rubriq,len(rubriq),itysup,nbval,irtyp,ientgb)

!       Entree gaz frais (si ce n'est pas NOZPPM, erreur)
    itysup = 0
    nbval  = nozppm
    irtyp  = 1
    rubriq = 'ientgf_zone_bord_lwc'
    call ecrsui(impavx,rubriq,len(rubriq),itysup,nbval,irtyp,ientgf)

!       FMENT (si ce n'est pas NOZPPM, erreur)
    itysup = 0
    nbval  = nozppm
    irtyp  = 2
    rubriq = 'fment_zone_bord_lwc'
    call ecrsui(impavx,rubriq,len(rubriq),itysup,nbval,irtyp,fment)

!       TKENT (si ce n'est pas NOZPPM, erreur)
    itysup = 0
    nbval  = nozppm
    irtyp  = 2
    rubriq = 'tkent_zone_bord_lwc'
    call ecrsui(impavx,rubriq,len(rubriq),itysup,nbval,irtyp,tkent)

#if defined(_CS_LANG_FR)
    car54=' Fin de l''ecriture des informations combustion LWC '
#else
    car54=' End writing combustion information (LWC)          '
#endif
    write(nfecra,1110)car54

  endif

! ---> Grandeurs complementaires pour la combustion CP

  if (ippmod(icpl3c).ge.0 .or. ippmod(iccoal).ge.0) then

!     Charbon PuLVerise : masse vol des charbons

    itysup = 0
    nbval  = 1
    irtyp  = 2
    do icha = 1, ncharb
      if (icha.le.nfmtch) then
        write(car2,'(I2.2)') icha
      else
        car2 = cindfc
      endif
      rubriq = 'masse_volumique_charbon'//car2
      call ecrsui(impavx,rubriq,len(rubriq),itysup,nbval,irtyp,   &
                  rhock(icha))
    enddo


!     Charbon PuLVerise : type de zones de bord, ientat, inmoxy, ientcp, timpat
!       x20, pour le calcul de rho au bord en entree

!       Numero des zones
    itysup = 3
    nbval  = 1
    irtyp  = 1
    rubriq = 'num_zone_fb_charbon_pulverise'
    call ecrsui(impavx,rubriq,len(rubriq),itysup,nbval,irtyp,izfppp)

!       Type entree air ou cp (si ce n'est pas NOZPPM, erreur)
    itysup = 0
    nbval  = nozppm
    irtyp  = 1
    rubriq = 'ientat_zone_bord_charbon_pulverise'
    call ecrsui(impavx,rubriq,len(rubriq),itysup,nbval,irtyp,ientat)

!       ientat, inmoxy et x20 ne servent pas pour le CP couple Lagrangien (cplphy)
    if (ippmod(iccoal).ge.0) then

      itysup = 0
      nbval  = nozppm
      irtyp  = 1
      rubriq = 'ientcp_zone_bord_charbon_pulverise'
      call ecrsui(impavx,rubriq,len(rubriq),itysup,nbval,irtyp,ientcp)

      itysup = 0
      nbval  = nozppm
      irtyp  = 1
      rubriq = 'inmoxy_zone_bord_charbon_pulverise'
      call ecrsui(impavx,rubriq,len(rubriq),itysup,nbval,irtyp,inmoxy)

      itysup = 0
      nbval  = nozppm
      irtyp  = 2

      idecal = 0
      do icha = 1, ncharb
        do iclapc = 1, nclpch(icha)
          icla = iclapc + idecal
          if (icha.le.nfmtch.and.iclapc.le.nfmtcl) then
            write(car2,'(I2.2)') icha
            write(car4,'(I4.4)') iclapc
          else
            car2 = cindfc
            car4 = cindfl
          endif
          rubriq = 'x20_zone_bord_charbon'//car2//'_classe'//car4
          call ecrsui(impavx,rubriq,len(rubriq),itysup,nbval,     &
                      irtyp,x20(1,icla))
        enddo
      enddo

    endif

!       Temperature
    itysup = 0
    nbval  = nozppm
    irtyp  = 2
    rubriq = 'timpat_zone_bord_charbon_pulverise'
    call ecrsui(impavx,rubriq,len(rubriq),itysup,nbval,irtyp,timpat)

#if defined(_CS_LANG_FR)
    car54=' Fin de l''ecriture des informations combustion CP    '
#else
    car54=' End writing combustion information (CP)            '
#endif
    write(nfecra,1110)car54

  endif

! ---> Grandeurs complementaires pour la FUEL

  if ( ippmod(icfuel).ge.0 ) then

!     Fioul : type de zones de bord, ientat, ientfl, timpat
!       qimpat et qimpfl  pour le calcul de rho au bord en entree

!       Numero des zones
    itysup = 3
    nbval  = 1
    irtyp  = 1
    rubriq = 'num_zone_fb_fuel'
    call ecrsui(impavx,rubriq,len(rubriq),itysup,nbval,irtyp,izfppp)

!       Type entree air ou fuel (si ce n'est pas NOZPPM, erreur)
    itysup = 0
    nbval  = nozppm
    irtyp  = 1
    rubriq = 'ientat_zone_bord_fuel'
    call ecrsui(impavx,rubriq,len(rubriq),itysup,nbval,irtyp,ientat)

    itysup = 0
    nbval  = nozppm
    irtyp  = 1
    rubriq = 'ientfl_zone_bord_fuel'
    call ecrsui(impavx,rubriq,len(rubriq),itysup,nbval,irtyp,ientfl)

    itysup = 0
    nbval  = nozppm
    irtyp  = 1
    RUBRIQ = 'inmoxy_zone_bord_fuel'
    call ecrsui(impavx,rubriq,len(rubriq),itysup,nbval,irtyp,   &
                inmoxy)

!       Timpat
    itysup = 0
    nbval  = nozppm
    irtyp  = 2
    rubriq = 'timpat_zone_bord_fuel'
    call ecrsui(impavx,rubriq,len(rubriq),itysup,nbval,irtyp,timpat)

!       Qimpat
    itysup = 0
    nbval  = nozppm
    irtyp  = 2
    rubriq = 'qimpat_zone_bord_fuel'
    call ecrsui(impavx,rubriq,len(rubriq),itysup,nbval,irtyp,qimpat)

!       Qimpfl
    itysup = 0
    nbval  = nozppm
    irtyp  = 2
    rubriq = 'qimpfl_zone_bord_fuel'
    call ecrsui(impavx,rubriq,len(rubriq),itysup,nbval,irtyp,qimpfl)

#if defined(_CS_LANG_FR)
    car54=' Fin de l''ecriture des informations combustion FUEL  '
#else
    car54=' End writing combustion information (FUEL)           '
#endif
    write(nfecra,1110)car54

  endif

! ---> Grandeurs complementaires pour les versions electriques

  iecr = 0

!     Recalage des CL pot des versions electriques

  if ( ippmod(ieljou).ge.1       ) then
    if (ielcor.eq.1) then

      iecr   = 1
      rubriq = 'coeff_recalage_joule'
      itysup = 0
      nbval  = 1
      irtyp  = 2

      call ecrsui(impavx,rubriq,len(rubriq),itysup,nbval,irtyp,coejou)

    endif
  endif

  if ( ippmod(ielarc).ge.1 .or. ippmod(ieljou).ge.1 ) then
    if (ielcor.eq.1) then

      iecr   = 1
      rubriq = 'ddpot_recalage_arc_elec'
      itysup = 0
      nbval  = 1
      irtyp  = 2

      call ecrsui(impavx,rubriq,len(rubriq),itysup,nbval,irtyp,dpot)

    endif
  endif

!     Termes sources des versions electriques

  if ( ippmod(ieljou).ge.1 .or.                                   &
       ippmod(ielarc).ge.1 .or.                                   &
       ippmod(ielion).ge.1       ) then

    iecr   = 1
    ipcefj = ipproc(iefjou)
    itysup = 1
    nbval  = 1
    irtyp  = 2

    rubriq = 'tsource_sc_ce_joule'
    call ecrsui(impavx,rubriq,len(rubriq),itysup,nbval,irtyp,     &
                propce(1,ipcefj))

  endif

  if (ippmod(ielarc).ge.1) then

    iecr   = 1
    ipcla1 = ipproc(ilapla(1))
    ipcla2 = ipproc(ilapla(2))
    ipcla3 = ipproc(ilapla(3))
    itysup = 1
    nbval  = 1
    irtyp  = 2

    rubriq = 'tsource_ns_ce_x_laplace'
    call ecrsui(impavx,rubriq,len(rubriq),itysup,nbval,irtyp,     &
                propce(1,ipcla1))

    rubriq = 'tsource_ns_ce_y_laplace'
    call ecrsui(impavx,rubriq,len(rubriq),itysup,nbval,irtyp,     &
                propce(1,ipcla2))

    rubriq = 'tsource_ns_ce_z_laplace'
    call ecrsui(impavx,rubriq,len(rubriq),itysup,nbval,irtyp,     &
                propce(1,ipcla3))

  endif

  if (iecr.ne.0) then
#if defined(_CS_LANG_FR)
    car54=' Fin de l''ecriture des informations electriques      '
#else
    car54=' End writing the electric information                '
#endif
    write(nfecra,1110)car54
  endif


!       Fermeture du fichiers suite auxiliaire
  call clssui(impavx,ierror)

  if (ierror.ne.0) then
    write(nfecra,8011) ficsui
  endif

  write(nfecra,1200)

endif
!     Fin de l'ecriture du fichier suite auxiliaire


return

!===============================================================================
! 6. FORMATS
!===============================================================================

#if defined(_CS_LANG_FR)

 1000 format(3x,'** Ecriture du fichier suite principal',/,       &
             3x,'   ----------------------------------- ',/)
 1100 format(' Debut de l''ecriture')
 1110 format('  ',A54)
 1200 format(' Fin de l''ecriture')
 2000 format(/,3x,'** Ecriture du fichier suite auxiliaire',/,    &
               3x,'   ------------------------------------ ',/)

 7001 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION :       A L''ECRITURE DU FICHIER SUITE        ',/,&
'@    =========                                               ',/,&
'@                                                            ',/,&
'@      Le nombre de scalaires maximal NSCAMX supporte par le ',/,&
'@        format d''ecriture du fichier suite est             ',/,&
'@        NFMTSC = ',i10                                       ,/,&
'@      On a ici un nombre de scalaires maximal superieur     ',/,&
'@        NSCAMX = ',i10                                       ,/,&
'@      On ne pourra pas relire les scalaires dont le numero  ',/,&
'@        est superieur                                       ',/,&
'@                                                            ',/,&
'@    Le calcul sera execute.                                 ',/,&
'@                                                            ',/,&
'@    Voir le sous-programme ecrava.                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 7002 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION :       A L''ECRITURE DU FICHIER SUITE        ',/,&
'@    =========                                               ',/,&
'@                                                            ',/,&
'@      Le nombre de flux de masse max NVARMX supporte par le ',/,&
'@        format d''ecriture du fichier suite est             ',/,&
'@        NFMTFL = ',i10                                       ,/,&
'@      On a ici un nombre de flux      maximal superieur     ',/,&
'@        NVARMX = ',i10                                       ,/,&
'@      On ne pourra pas relire les flux      dont le numero  ',/,&
'@        est superieur                                       ',/,&
'@                                                            ',/,&
'@    Le calcul sera execute.                                 ',/,&
'@                                                            ',/,&
'@    Voir le sous-programme ecrava.                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 7003 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION :       A L''ECRITURE DU FICHIER SUITE        ',/,&
'@    =========                                               ',/,&
'@                                                            ',/,&
'@      Le nombre de moments       max NBMOMX supporte par le ',/,&
'@        format d''ecriture du fichier suite est             ',/,&
'@        NFMTMO = ',i10                                       ,/,&
'@      On a ici un nombre de moments   maximal superieur     ',/,&
'@        NBMOMX = ',i10                                       ,/,&
'@      On ne pourra pas relire les moments   dont le numero  ',/,&
'@        est superieur                                       ',/,&
'@                                                            ',/,&
'@    Le calcul sera execute.                                 ',/,&
'@                                                            ',/,&
'@    Voir le sous-programme ecrava.                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 7004 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION :       A L''ECRITURE DU FICHIER SUITE        ',/,&
'@    =========                                               ',/,&
'@                                                            ',/,&
'@      Le nombre de charbons      max NCHARM supporte par le ',/,&
'@        format d''ecriture du fichier suite est             ',/,&
'@        NFMTCH = ',i10                                       ,/,&
'@      On a ici un nombre de charbons  maximal superieur     ',/,&
'@        NCHARM = ',i10                                       ,/,&
'@      On ne pourra pas relire certaines informations        ',/,&
'@        relatives aux charbons dont le numero               ',/,&
'@        est superieur                                       ',/,&
'@                                                            ',/,&
'@    Le calcul sera execute.                                 ',/,&
'@                                                            ',/,&
'@    Voir le sous-programme ecrava.                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 7005 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION :       A L''ECRITURE DU FICHIER SUITE        ',/,&
'@    =========                                               ',/,&
'@                                                            ',/,&
'@      Le nombre de classes par charbon max NCPCMX supporte  ',/,&
'@        par le format d''ecriture du fichier suite est      ',/,&
'@        NFMTCL = ',i10                                       ,/,&
'@      On a ici un nombre de classes par charbon superieur   ',/,&
'@        NCPCMX = ',i10                                       ,/,&
'@      On ne pourra pas relire certaines informations        ',/,&
'@        relatives aux classes  dont le numero               ',/,&
'@        est superieur                                       ',/,&
'@                                                            ',/,&
'@    Le calcul sera execute.                                 ',/,&
'@                                                            ',/,&
'@    Voir le sous-programme ecrava.                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 8000 format(                                                          &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ERREUR A L''OUVERTURE DU FICHIER SUITE      ',/,&
'@    =========                                 AVAL PRINCIPAL',/,&
'@                                                            ',/,&
'@    Verifier que le fichier ',a13,'peut etre                ',/,&
'@            cree dans le repertoire de travail.             ',/,&
'@                                                            ',/,&
'@    Le calcul se poursuit...                                ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8001 format(                                                          &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ERREUR A L''OUVERTURE DU FICHIER SUITE      ',/,&
'@    =========                                AVAL AUXILIAIRE',/,&
'@                                                            ',/,&
'@    Verifier que le fichier ',a13,'peut etre                ',/,&
'@            cree dans le repertoire de travail.             ',/,&
'@                                                            ',/,&
'@    Le calcul se poursuit...                                ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 8010 format(                                                          &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ERREUR A LA FERMETURE DU FICHIER SUITE      ',/,&
'@    =========                                 AVAL PRINCIPAL',/,&
'@                                                            ',/,&
'@    Probleme sur le fichier de nom (',a13,')                ',/,&
'@                                                            ',/,&
'@    Le calcul se poursuit...                                ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8011 format(                                                          &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ERREUR A LA FERMETURE DU FICHIER SUITE      ',/,&
'@    =========                                AVAL AUXILIAIRE',/,&
'@                                                            ',/,&
'@    Probleme sur le fichier de nom (',a13,')                ',/,&
'@                                                            ',/,&
'@    Le calcul se poursuit...                                ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

#else

 1000 format(3x,'** Writing the main restart file',/,             &
             3x,'   -----------------------------',/)
 1100 format(' Start writing'                                      )
 1110 format('  ',a54                                              )
 1200 format(' End writing'                                        )
 2000 format(/,3x,'** Writing the auxiliary restart file',/,      &
               3x,'   ----------------------------------',/)

 7001 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: WHILE WRITING THE RESTART FILE                 ',/,&
'@    ========                                                ',/,&
'@                                                            ',/,&
'@      The maximum number of scalars NSCAMX handled by the   ',/,&
'@        restart file writing format is                      ',/,&
'@        NFMTSC = ',i10                                       ,/,&
'@      The current maximum number of scalars is greater.     ',/,&
'@        NSCAMX = ',i10                                       ,/,&
'@      The scalars with a larger number will not be read.    ',/,&
'@                                                            ',/,&
'@    The calculation will be run.                            ',/,&
'@                                                            ',/,&
'@    Refer to the subroutine ecrava.                         ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 7002 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: WHILE WRITING THE RESTART FILE                 ',/,&
'@    ========                                                ',/,&
'@                                                            ',/,&
'@      The maximum number of mass flux NVARMX handled by the ',/,&
'@        restart file writing format is                      ',/,&
'@        NFMTFL = ',i10                                       ,/,&
'@      The current maximum number of mass fluxes is greater. ',/,&
'@        NVARMX = ',i10                                       ,/,&
'@      The fluxes with a larger number will not be read.     ',/,&
'@                                                            ',/,&
'@    The calculation will be run.                            ',/,&
'@                                                            ',/,&
'@    Refer to the subroutine ecrava.                         ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 7003 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: WHILE WRITING THE RESTART FILE                 ',/,&
'@    ========                                                ',/,&
'@                                                            ',/,&
'@      The maximum number of moments NBMOMX handled by the   ',/,&
'@        restart file writing format is                      ',/,&
'@        NFMTMO = ',i10                                       ,/,&
'@      The current maximum number of moments is greater.     ',/,&
'@        NBMOMX = ',i10                                       ,/,&
'@      The moments with a larger number will not be read.    ',/,&
'@                                                            ',/,&
'@    The calculation will be run.                            ',/,&
'@                                                            ',/,&
'@    Refer to the subroutine ecrava.                         ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 7004 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: WHILE WRITING THE RESTART FILE                 ',/,&
'@    ========                                                ',/,&
'@                                                            ',/,&
'@      The maximum number of coals NCHARM handled by the     ',/,&
'@        restart file writing format is                      ',/,&
'@        NFMTCH = ',i10                                       ,/,&
'@      The current maximum number of coals is greater.       ',/,&
'@        NCHARM = ',i10                                       ,/,&
'@      Some information relative to coals with a greater     ',/,&
'@        number will not be read.                            ',/,&
'@                                                            ',/,&
'@    The calculation will be run.                            ',/,&
'@                                                            ',/,&
'@    Refer to the subroutine ecrava.                         ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 7005 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: WHILE WRITING THE RESTART FILE                 ',/,&
'@    ========                                                ',/,&
'@                                                            ',/,&
'@      The number of coal classes NCPCMX handled by the      ',/,&
'@        restart file writing format is                      ',/,&
'@        NFMTCL = ',i10                                       ,/,&
'@      The current number of coal classes is greater.        ',/,&
'@        NCPCMX = ',i10                                       ,/,&
'@      Some information relative to classes with a greater   ',/,&
'@        number will not be read.                            ',/,&
'@                                                            ',/,&
'@    The calculation will be run.                            ',/,&
'@                                                            ',/,&
'@    Refer to the subroutine ecrava.                         ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 8000 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: ERROR WHILE OPENING THE MAIN RESTART FILE      ',/,&
'@    ========                                                ',/,&
'@                                                            ',/,&
'@    Verify that the file ',a13,'can be created              ',/,&
'@            in the working directory.                       ',/,&
'@                                                            ',/,&
'@    The calculation will be run.                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8001 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: ERROR WHILE OPENING THE AUXILIARY RESTART FILE ',/,&
'@    ========                                                ',/,&
'@                                                            ',/,&
'@    Verify that the file ',a13,'can be created              ',/,&
'@            in the working directory.                       ',/,&
'@                                                            ',/,&
'@    The calculation will be run.                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 8010 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: ERROR WHILE CLOSING THE MAIN RESTART FILE      ',/,&
'@    ========                                                ',/,&
'@                                                            ',/,&
'@    Problem with the file of name (',a13,')                 ',/,&
'@                                                            ',/,&
'@    The calculation will be run.                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8011 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: ERROR WHILE CLOSING THE AUXILIARY RESTART FILE ',/,&
'@    ========                                                ',/,&
'@                                                            ',/,&
'@    Problem with the file of name (',a13,')                 ',/,&
'@                                                            ',/,&
'@    The calculation will be run.                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

#endif


end subroutine
