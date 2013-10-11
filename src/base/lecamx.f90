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

subroutine lecamx &
!================

 ( ncelet , ncel   , nfabor , nvar   , nscal  ,                   &
   dt     , propce ,                                              &
   coefa  , coefb  , frcxt  , prhyd  )

!===============================================================================

! FONCTION :
! ----------
! LECTURE DU FICHIER SUITE AUXILIAIRE

! ON PASSE ICI SI ILEAUX = 1
! ON S'ARRETE SI ON NE PEUT PAS OUVRIR LE FICHIER
!          OU SI CE N'EST PAS UN FICHIER AUXILIAIRE
!          OU SI NCEL N'EST PAS CORRECT
!          OU SI ON NE PEUT PAS LIRE JPHAS, JTURB, JDTVAR
!          OU SI ON NE PEUT PAS LIRE UNE MOYENNE QU'ON VEUT POURSUIVRE

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
! dt(ncelet)       ! tr ! --> ! pas de temps                                   !
! rtp              ! tr ! --> ! variables de calcul au centre des              !
! (ncelet,*)       !    !     !    cellules (instant courant        )          !
! propce           ! tr ! --> ! proprietes physiques au centre des             !
! (ncelet,*)       !    !     !    cellules                                    !
! coefa, coefb     ! tr ! --> ! conditions aux limites aux                     !
!  (nfabor,*)      !    !     !    faces de bord                               !
! frcxt(3,ncelet)  ! tr ! --> ! force exterieure generant la pression          !
!                  !    !     !  hydrostatique                                 !
! prhyd(ncelet)    ! ra ! --> ! hydrostatic pressure predicted                 !
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens, only: ndimfb
use cstphy
use cstnum
use entsor
use optcal
use pointe
use numvar
use albase
use alstru
use alaste
use parall
use ppppar
use ppthch
use ppincl
use coincl
use cpincl
use cs_fuel_incl
use elincl
use ppcpfu
use mesh, only: isympa
use cs_c_bindings
use field
!===============================================================================

implicit none

! Arguments

integer          ncelet , ncel   , nfabor
integer          nvar   , nscal

double precision dt(ncelet)
double precision propce(ncelet,*)
double precision coefa(ndimfb,*), coefb(ndimfb,*)
double precision frcxt(3,ncelet), prhyd(ncelet)

! Local variables


character        rubriq*64,car4*4,car2*2
character        cnum4*4, car54*54
character        cindfp*2,cindfs*4,cindff*4,cindfm*4
character        cindfc*2,cindfl*4
character        cphase*2
character        nomflu(nvarmx)*18,nomcli(nvarmx)*18
character        cstruc(nstrmx)*2, cindst*2
character        ficsui*32
logical          lprev
integer          iel   , ifac, ii, istr
integer          ivar  , iscal , jphas , isco
integer          idecal, iclapc, icha  , icla
integer          imom  , imold
integer          jdtvar
integer          jortvm, ipcvmx, ipcvmy, ipcvmz
integer          iclvar, iclvaf
integer          idtcm
integer          iptsna, iptsta, iptsca
integer          numero, ipcefj, ipcla1, ipcla2, ipcla3
integer          iok   , inifok, iokdt
integer          ncelok, nfaiok, nfabok, nsomok
integer          ierror, irtyp,  itysup, nbval
integer          nberro, inierr, ivers
integer          ilu   , ilecec, ideblu, iannul, ierrch
integer          impamx
integer          nfmtsc, nfmtfl, nfmtmo, nfmtch, nfmtcl
integer          nfmtst
integer          jturb , jtytur, jale
integer          f_id, nfld, iflmas, iflmab
integer          ngbstr(2)
double precision d2s3  , tsrii , cdtcm
double precision tmpstr(27)

integer, allocatable, dimension(:) :: mflnum
double precision, dimension(:), pointer :: sval

!===============================================================================

!===============================================================================
! 0. INITIALISATION
!===============================================================================

!  ---> Banniere
write(nfecra,1000)

!  --->  On code en chaine le numero des phases et scalaires

!     Nombre max pour les formats choisis
nfmtsc = 9999
nfmtfl = 9999
nfmtmo = 9999
nfmtch = 99
nfmtcl = 9999

!     Indefini a 2 et 4 caracteres
CINDFP='YY'
CINDFS='YYYY'
CINDFF='YYYY'
CINDFM='YYYY'
CINDFC='YY'
CINDFL='YYYY'

!  Codage en chaine de caracteres du numero de la phase
cphase='01'

!     Avertissement
if(nscamx.gt.nfmtsc) then
  write(nfecra,8001)nfmtsc,nscamx
endif
if(nvarmx.gt.nfmtfl) then
  write(nfecra,8002)nfmtfl,nvarmx
endif
if(nbmomx.gt.nfmtmo) then
  write(nfecra,8003)nfmtmo,nbmomx
endif

!===============================================================================
! 1. OUVERTURE DU FICHIER SUITE AUXILIAIRE
!===============================================================================

!  ---> Ouverture du fichier (ILECEC = 1 : lecture)
ilecec = 1
ficsui = 'auxiliary'
call opnsui(ficsui,len(ficsui),ilecec,impamx,ierror)
!==========
if (ierror.ne.0) then
  write(nfecra,9000) ficsui
  call csexit (1)
endif

! ---> Debut de la lecture
write(nfecra,1100)

!===============================================================================
! 2. ENTETES DU FICHIER SUITE OU STOP
!===============================================================================

!  Rubrique "fichier suite aux"
!        Pourrait porter le numero de version si besoin.
!        On ne se sert pas de IVERS pour le moment.

itysup = 0
nbval  = 1
irtyp  = 1
RUBRIQ = 'version_fichier_suite_auxiliaire'
call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp,ivers,   &
            ierror)

if (ierror.ne.0) then
  write(nfecra,9100)ficsui
  call csexit (1)
endif

!     Supports

call tstsui(impamx,ncelok,nfaiok,nfabok,nsomok)
!==========

if (ncelok.eq.0) then
  write(nfecra,9101)
  call csexit (1)
endif

IF (NFAIOK.EQ.0) WRITE(NFECRA,8200)'internes','internes'

IF (NFABOK.EQ.0) WRITE(NFECRA,8200)'de bord ','de bord '


!     On n'a besoin que
!       du nombre de phases et du modele de turbulence

nberro = 0

itysup = 0
nbval  = 1
irtyp  = 1

RUBRIQ = 'nombre_phases'
call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp,jphas,   &
            ierror)
nberro=nberro+ierror

!  ---> On s'arrete si erreur
!     Si on ne peut pas relire un entier, c'est que le fichier n'est pas bon
if (nberro.ne.0) then
  CAR54 ='ERREUR A lA LECTURE DES DIMENSIONS (NPHAS)            '
  write(nfecra,9200)car54
  call csexit (1)
endif


!  ---> On ne sait relire que des calculs monophasiques
if (jphas.ne.1) then
  write(nfecra,8205) jphas
  call csexit(1)
endif

!     Modele de turbulence

nberro = 0

RUBRIQ = 'modele_turbulence_phase'//CPHASE
itysup = 0
nbval  = 1
irtyp  = 1
call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp,       &
            jturb,ierror)
nberro=nberro+ierror
jtytur=jturb/10

! --->  Stop si erreur
!     Si on ne peut pas relire un entier, c'est que le fichier n'est pas bon
if (nberro.ne.0) then
  CAR54 ='ERREUR A lA LECTURE DES MODELES DE TURBULENCE         '
  write(nfecra,9200)car54
  call csexit (1)
endif

!     Methode ALE

nberro = 0

RUBRIQ = 'methode_ALE'
itysup = 0
nbval  = 1
irtyp  = 1
call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp,jale,    &
            ierror)
nberro=nberro+ierror

! --->  Message si erreur (pas de stop pour compatibilite avec les fichiers anterieurs)
!       -> on n'affiche le message que si IALE=1 (sinon RAS)
if (nberro.ne.0) then
  if (iale.eq.1) write(nfecra,9210)
  jale = 0
endif

! ---> Pas d'iteration d'initialisation si suite de calcul ALE
if (italin.eq.-999) then
  if (iale.eq.1 .and. jale.eq.1) then
    italin = 0
  else if (iale.eq.1) then
    italin = 1
  else
    italin = 0
  endif
endif

! --->  Donnees modifiees
if (iturb .ne. jturb)                             &
     write(nfecra,8220) iturb, jturb

CAR54 =' Fin de la lecture des options                        '
write(nfecra,1110)car54

!===============================================================================
! 3. PROPRIETES PHYSIQUES
!===============================================================================

nberro = 0

!     On lit les infos des phases communes

! ---> Point de reference de pression
!     On lit les coordonnees si l'utilisateur n'a rien specifie, i.e.
!       si IXYZP0=-1, et on met IXYZP0 a 1 pour eviter de le changer ensuite.
if (ixyzp0.eq.-1) then
  RUBRIQ = 'ref_presstot'//CPHASE
  itysup = 0
  nbval  = 3
  irtyp  = 2
  call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp,     &
              xyzp0(1),ierror)
  nberro = nberro+ierror
  if (ierror.eq.0) then
    write(nfecra,7000) (xyzp0(ii),ii=1,3)
    ixyzp0 = 1
  endif
endif

! Here the physical variables below are required for the low-Mach algorithm
if (idilat.eq.3) then

  !the reference density updated with the low-Mach algorithm
  rubriq = 'ro0'//cphase
  itysup = 0
  nbval  = 1
  irtyp  = 2
  call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp,ro0,ierror)
  nberro=nberro+ierror

  ! the thermodynamic pressure for the previous time step
  rubriq = 'pther'//cphase
  itysup = 0
  nbval  = 1
  irtyp  = 2
  call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp,   &
       pther,ierror)
  nberro=nberro+ierror
endif


! ---> Masse volumique
!     On la lit, qu'elle soit extrapolee ou pas,
!       pour permettre les sous-relaxations
!     Pour les suites a rho constant, cependant, on ne la lit pas,
!       afin de ne pas ecraser la valeur RO0 fixee par l'utilisateur
!       et eventuellement modifiee.

inierr = 0

if (irovar.eq.1) then

  ! Masse volumique - cellules
  rubriq = 'rho_ce_phase01'
  itysup = 1
  nbval  = 1
  irtyp  = 2
  call field_get_val_s(icrom, sval)
  call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp,sval,ierror)
  nberro = nberro+ierror
  inierr = inierr+ierror

  ! Masse volumique - faces de bord
  if (nfabok.eq.1) then
    RUBRIQ = 'rho_fb_phase'//CPHASE
    itysup = 3
    nbval  = 1
    irtyp  = 2
    call field_get_val_s(ibrom, sval)
    call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp,sval,ierror)
    nberro = nberro+ierror
    inierr = inierr+ierror
  endif

  !     Si on a reussi a initialiser la masse volumique aux cellules ET
  !       aux faces de bord, on l'indique (pour schtmp)
  if (nfabok.eq.1.and.inierr.eq.0) then
    initro = 1
  endif

else
  !     Si la masse volumique est constante, elle est initialisee
  !       correctement
  initro = 1
endif

! ---> Viscosite moleculaire et "turbulente" ou de "sous-maille"
!     Si elle est extrapolee en temps, on la lit
!     Si on reussit, on l'indique

!     On cherche a lire uniquement si on doit extrapoler en temps
if(iviext.gt.0) then

  inierr = 0

  !         Viscosite moleculaire - cellules
  !         Uniquement si elle est variable
  if(ivivar.eq.1) then
    RUBRIQ = 'viscl_ce_phase'//CPHASE
    itysup = 1
    nbval  = 1
    irtyp  = 2
    call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp,   &
         propce(1,ipproc(iviscl)),ierror)
    nberro = nberro+ierror
    inierr = inierr+ierror
  endif

  !         Viscosite turbulente ou de sous-maille - cellules
  RUBRIQ = 'visct_ce_phase'//CPHASE
  itysup = 1
  nbval  = 1
  irtyp  = 2
  call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp,     &
       propce(1,ipproc(ivisct)),ierror)
  nberro = nberro+ierror
  inierr = inierr+ierror

  !     Si on a initialise les viscosites, on l'indique (pour schtmp)
  if (inierr.eq.0) then
    initvi = 1
  endif

endif


! ---> Chaleur massique
!     On cherche a la lire si elle est variable
!       et qu'on l'extrapole ou qu'on est en Joule
!     Si on reussit, on l'indique
!       (Ca sert quand elle est extrapolee en temps
!        et quand l'utilisateur peut s'en servir pour passer
!        de H a T, comme en effet Joule par exemple).

if((icpext.gt.0.and.icp.gt.0).or.                &
     (ippmod(ieljou).ge.1.and.icp.gt.0)) then

  inierr = 0

  !         Chaleur massique - cellules
  RUBRIQ = 'cp_ce_phase'//CPHASE
  itysup = 1
  nbval  = 1
  irtyp  = 2
  call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp,     &
       propce(1,ipproc(icp)),ierror)
  nberro = nberro+ierror
  inierr = inierr+ierror

  !     Si on a initialise Cp, on l'indique (pour schtmp)
  if (inierr.eq.0) then
    initcp = 1
  endif

endif

!     Si on a des scalaires, on lit a diffusivite
!       si le scalaire a  un correspondant et
!       si on doit extrapoler la diffusivite
!       (et qu'elle est variable, et que le scalaire n'est pas une variance)

if(nscal.gt.0) then

  ! --->  Donnees modifiees
  do iscal = 1, nscal
    isco = iscold(iscal)
    if ( isco         .gt.0.and.                                  &
         ivsext(iscal).gt.0.and.                                  &
         ivisls(iscal).gt.0.and.                                  &
        (iscavr(iscal).le.0.or.iscavr(iscal).gt.nscal) ) then

      inierr = 0

      ! Cell diffusivity
      if (isco.le.nfmtsc) then
        write(car4,'(i4.4)')isco
        rubriq = 'visls_ce_scalaire'//car4
        itysup = 1
        nbval  = 1
        irtyp  = 2
        call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp, &
                    propce(1,ipproc(ivisls(iscal))),ierror)
        nberro = nberro+ierror
        inierr = inierr+ierror
      else
        nberro = nberro-1
        inierr = inierr-1
      endif

!     Si on a initialise visls, on l'indique (pour schtmp)
      if (inierr.eq.0) then
        initvs(iscal) = 1
      endif
    endif

  enddo

!     Pour les variances, il suffit de dire qu'on a initialise ou non
!       selon ce qu'on a fait au scalaire correspondant
  do iscal = 1, nscal
    if(iscavr(iscal).gt.0.and.iscavr(iscal).le.nscal) then
      initvs(iscal) = initvs(iscavr(iscal))
    endif
  enddo

endif

!     Si erreur, on previent mais pas stop :
!       auparavant on n'avait pas stocke les prop phy
!         si on n'etait pas a l'ordre 2
!       c'est discutable pour rho

if (nberro.ne.0) then
  CAR54 = 'LECTURE DES PROPRIETES PHYSIQUES                    '
  write(nfecra,8300)car54
endif

CAR54 = ' Fin de la lecture des proprietes physiques           '
write(nfecra,1110)car54


!===============================================================================
! 4. PAS DE TEMPS
!===============================================================================


!  ---> Indicateur de pas de temps variable
RUBRIQ = 'indic_dt_variable'
itysup = 0
nbval  = 1
irtyp  = 1
call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp,jdtvar,  &
            ierror)

!  ---> On s'arrete si erreur
!     Si on ne peut pas relire un entier, c'est que le fichier n'est pas bon
if (ierror.ne.0) then
  CAR54 ='ERREUR A LA LECTURE DU MODE DE MARCHE EN TEMPS        '
  write(nfecra,9200)car54
  call csexit(1)
endif

!  ---> Pas de temps
!     Rq : si marche en temps differente, on conserve la valeur par defaut
!          DTREF imposee dans INIVA0

nberro = 0

if (idtvar.ne.jdtvar) then
  write(nfecra,8400)idtvar,jdtvar,dtref

elseif (idtvar.eq.1) then
  RUBRIQ = 'dt_variable_temps'
  itysup = 0
  nbval  = 1
  irtyp  = 2
  call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp,dt(1), &
              ierror)
  nberro=nberro+ierror
  do iel = 1, ncel
    dt(iel) = dt(1)
  enddo

elseif (idtvar.eq.2) then
  RUBRIQ = 'dt_variable_espace_ce'
  itysup = 1
  nbval  = 1
  irtyp  = 2
  call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp,dt,    &
              ierror)
  nberro=nberro+ierror
endif

!     Si erreur, on previent mais pas stop :
!       le pas de temps n'est pas une donnee majeure
!       c'est discutable

if (nberro.ne.0) then
  CAR54 = 'LECTURE DU PAS DE TEMPS                               '
  write(nfecra,8300)car54
endif

CAR54 = ' Fin de la lecture du pas de temps                    '
write(nfecra,1110)car54

!===============================================================================
! 5. FLUX DE MASSE
!===============================================================================
!     Pour retrouver la correspondance entre les variables
!     et les flux de masse, on utilise le nom de chaque variable
!     (nomflu(I)= nom de la ieme variable)
!     Ensuite, pour chaque variable, si on a deja le flux, on ne fait
!       rien, sinon on lit quel est le numero local du
!       flux qui lui est associe (en pratique 1 ou 2) et le flux lui meme.

!     Les flux de masse ne sont a lire que si les supports des faces
!     de bord ou des faces internes coincident

!     On lit d'abord le flux de masse (a l'instant n) et
!       ensuite a l'instant precedent si on est en schema en temps
!       particulier (ISTMPF NE 1)

if (nfaiok.eq.1 .or. nfabok.eq.1) then

  nberro=0

  ! Initialize work arrays

  call field_get_n_fields(nfld)

  allocate(mflnum(nfld))

  do f_id = 1, nfld
    mflnum(f_id) = 0
  enddo

  ! Name of flux associated with variable in previous calculation
  nomflu(ipr)='fm_p_phase'//cphase
  nomflu(iu)='fm_u_phase'//cphase
  nomflu(iv)='fm_v_phase'//cphase
  nomflu(iw)='fm_w_phase'//cphase
  if (itytur.eq.2.and. (jtytur.eq.2.or.jtytur.eq.5)) then
    nomflu(ik)='fm_k_phase'//cphase
    nomflu(iep)='fm_eps_phase'//cphase
  elseif (itytur.eq.2.and.jtytur.eq.3) then
    nomflu(ik)='fm_R11_phase'//cphase
    nomflu(iep)='fm_eps_phase'//cphase
  elseif (itytur.eq.2.and.jturb.eq.60) then
    nomflu(ik)='fm_k_phase'//cphase
    nomflu(iep)='fm_omega_phase'//cphase
  elseif (itytur.eq.3 .and. (jtytur.eq.2.or.jtytur.eq.50)) then
    nomflu(ir11)='fm_k_phase'//cphase
    nomflu(ir22)='fm_k_phase'//cphase
    nomflu(ir33)='fm_k_phase'//cphase
    nomflu(ir12)='fm_k_phase'//cphase
    nomflu(ir13)='fm_k_phase'//cphase
    nomflu(ir23)='fm_k_phase'//cphase
    nomflu(iep)='fm_eps_phase'//cphase
  elseif (itytur.eq.3.and.jtytur.eq.3) then
    nomflu(ir11)='fm_R11_phase'//cphase
    nomflu(ir22)='fm_R22_phase'//cphase
    nomflu(ir33)='fm_R33_phase'//cphase
    nomflu(ir12)='fm_R12_phase'//cphase
    nomflu(ir13)='fm_R13_phase'//cphase
    nomflu(ir23)='fm_R23_phase'//cphase
    nomflu(iep)='fm_eps_phase'//cphase
    if (iturb.eq.32.and.jturb.eq.32) then
      nomflu(ial)='fm_alp_phase'//cphase
    endif
    if (iturb.eq.32.and.jturb.ne.32) then
      nomflu(ial)='fm_eps_phase'//cphase
    endif
  elseif (itytur.eq.3.and.jturb.eq.60) then
    nomflu(ir11)='fm_k_phase'//cphase
    nomflu(ir22)='fm_k_phase'//cphase
    nomflu(ir33)='fm_k_phase'//cphase
    nomflu(ir12)='fm_k_phase'//cphase
    nomflu(ir13)='fm_k_phase'//cphase
    nomflu(ir23)='fm_k_phase'//cphase
    nomflu(iep)='fm_omega_phase'//cphase
  elseif (itytur.eq.5.and.jtytur.eq.2) then
    nomflu(ik)='fm_k_phase'//cphase
    nomflu(iep)='fm_eps_phase'//cphase
    nomflu(iphi)='fm_k_phase'//cphase
    if (iturb.eq.50) then
      nomflu(ifb)='fm_k_phase'//cphase
    elseif (iturb.eq.51) then
      nomflu(ial)='fm_k_phase'//cphase
    endif
  elseif (itytur.eq.5.and.jtytur.eq.3) then
    nomflu(ik)='fm_R11_phase'//cphase
    nomflu(iep)='fm_eps_phase'//cphase
    nomflu(iphi)='fm_R11_phase'//cphase
    nomflu(ifb)='fm_R11_phase'//cphase
    if (iturb.eq.50) then
      nomflu(ifb)='fm_R11_phase'//cphase
    elseif (iturb.eq.51) then
      nomflu(ial)='fm_R11_phase'//cphase
    endif
  elseif (iturb.eq.50.and.jturb.eq.50) then
    nomflu(ik)='fm_k_phase'//cphase
    nomflu(iep)='fm_eps_phase'//cphase
    nomflu(iphi)='fm_phi_phase'//cphase
    nomflu(ifb)='fm_fb_phase'//cphase
  elseif (iturb.eq.51.and.jturb.eq.51) then
    nomflu(ik)='fm_k_phase'//cphase
    nomflu(iep)='fm_eps_phase'//cphase
    nomflu(iphi)='fm_phi_phase'//cphase
    nomflu(ial)='fm_al_phase'//cphase
  elseif (iturb.eq.50.and.jturb.eq.51) then
    nomflu(ik)='fm_k_phase'//cphase
    nomflu(iep)='fm_eps_phase'//cphase
    nomflu(iphi)='fm_phi_phase'//cphase
    nomflu(ifb)='fm_al_phase'//cphase
  elseif (iturb.eq.51.and.jturb.eq.50) then
    nomflu(ik)='fm_k_phase'//cphase
    nomflu(iep)='fm_eps_phase'//cphase
    nomflu(iphi)='fm_phi_phase'//cphase
    nomflu(ial)='fm_fb_phase'//cphase
  elseif (iturb.eq.50.and.jturb.eq.60) then
    nomflu(ik)='fm_k_phase'//cphase
    nomflu(iep)='fm_omega_phase'//cphase
    nomflu(iphi)='fm_k_phase'//cphase
    nomflu(ifb)='fm_k_phase'//cphase
  elseif (iturb.eq.60.and.(jtytur.eq.2.or.jturb.eq.50)) then
    nomflu(ik  )='fm_k_phase'//cphase
    nomflu(iomg)='fm_eps_phase'//cphase
  elseif (iturb.eq.60.and.jtytur.eq.3) then
    nomflu(ik  )='fm_R11_phase'//cphase
    nomflu(iomg)='fm_eps_phase'//cphase
  elseif (iturb.eq.60.and.jturb.eq.60) then
    nomflu(ik  )='fm_k_phase'//cphase
    nomflu(iomg)='fm_omega_phase'//cphase
  elseif (iturb.eq.70.and.jturb.eq.70) then
    nomflu(inusa)='fm_nusa_phase'//cphase
  endif
  if (nscal.gt.0) then
    do iscal = 1, nscal
      if (iscold(iscal).gt.0) then
        if (iscold(iscal).le.nfmtsc) then
          write(car4,'(i4.4)')iscold(iscal)
        else
          car4 = cindfs
        endif
        nomflu(isca(iscal))='fm_scalaire'//car4
      endif
    enddo
  endif
  if (iale.eq.1) then
    nomflu(iuma)='fm_vit_maill_u'
    nomflu(ivma)='fm_vit_maill_v'
    nomflu(iwma)='fm_vit_maill_w'
  endif

  ! For variables
  do ivar = 1, nvar

    f_id = ivarfl(ivar)

    ! Is there an associated flux ?

    call field_get_key_int(f_id, kimasf, iflmas)

    if (iflmas.ge.0) then

      ! If flux has not been read yet, do it
      if (mflnum(iflmas).eq.0) then

        rubriq = nomflu(ivar)
        itysup = 0
        nbval  = 1
        irtyp  = 1
        call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp, &
                    numero,ierror)
        nberro = nberro+ierror
        ! If it exists, read it
        if (ierror.eq.0) then
          inifok = 0
          if (numero.gt.0) then
            if (numero.le.nfmtfl) then
              write(car4,'(i4.4)') numero
            else
              car4 = cindff
            endif
            if (nfaiok.eq.1) then
              call field_get_val_s(iflmas, sval)
              rubriq = 'flux_masse_fi_'//car4
              itysup = 2
              nbval  = 1
              irtyp  = 2
              call lecsui(impamx, rubriq, len(rubriq), itysup, nbval,  &
                          irtyp, sval, ierror)
              nberro = nberro+ierror
              if (ierror.eq.0) inifok = inifok+1
            endif
            if (nfabok.eq.1) then
              call field_get_key_int(f_id, kbmasf, iflmab)
              call field_get_val_s(iflmab, sval)
              rubriq = 'flux_masse_fb_'//car4
              itysup = 3
              nbval  = 1
              irtyp  = 2
              call lecsui(impamx, rubriq, len(rubriq), itysup, nbval,  &
                          irtyp, sval, ierror)
              nberro = nberro+ierror
              if (ierror.eq.0) inifok = inifok+1
            endif
          endif
          ! If everything is OK, mark this flux as read
          if (inifok.eq.2) then
            mflnum(iflmas) = 1
          endif
        endif
      endif
    endif
  enddo

  ! An now, the same for the preceding time step

  ! Initialize work arrays

  do f_id = 1, nfld
    mflnum(f_id) = 0
  enddo

  ! Name of flux associated with variable for previous calculation
  nomflu(ipr)='fm_a_p_phase'//cphase
  nomflu(iu)='fm_a_u_phase'//cphase
  nomflu(iv)='fm_a_v_phase'//cphase
  nomflu(iw)='fm_a_w_phase'//cphase
  if (itytur.eq.2 .and. (jtytur.eq.2.or.jtytur.eq.5)) then
    nomflu(ik)='fm_a_k_a_phase'//cphase
    nomflu(iep)='fm_a_eps_phase'//cphase
  elseif (itytur.eq.2.and.jtytur.eq.3) then
    nomflu(ik)='fm_a_R11_phase'//cphase
    nomflu(iep)='fm_a_eps_phase'//cphase
  elseif (itytur.eq.2.and.jturb.eq.60) then
    nomflu(ik)='fm_a_k_phase'//cphase
    nomflu(iep)='fm_a_omega_phase'//cphase
  elseif (itytur.eq.3 .and. (jtytur.eq.2.or.jtytur.eq.5)) then
    nomflu(ir11)='fm_a_k_phase'//cphase
    nomflu(ir22)='fm_a_k_phase'//cphase
    nomflu(ir33)='fm_a_k_phase'//cphase
    nomflu(ir12)='fm_a_k_phase'//cphase
    nomflu(ir13)='fm_a_k_phase'//cphase
    nomflu(ir23)='fm_a_k_phase'//cphase
    nomflu(iep)='fm_a_eps_phase'//cphase
  elseif (itytur.eq.3.and.jtytur.eq.3) then
    nomflu(ir11)='fm_a_R11_phase'//cphase
    nomflu(ir22)='fm_a_R22_phase'//cphase
    nomflu(ir33)='fm_a_R33_phase'//cphase
    nomflu(ir12)='fm_a_R12_phase'//cphase
    nomflu(ir13)='fm_a_R13_phase'//cphase
    nomflu(ir23)='fm_a_R23_phase'//cphase
    nomflu(iep)='fm_a_eps_phase'//cphase
    if (iturb.eq.32.and.jturb.eq.32) then
      nomflu(ial)='fm_a_alp_phase'//cphase
    endif
    if (iturb.eq.32.and.jturb.ne.32) then
      nomflu(ial)='fm_a_eps_phase'//cphase
    endif
  elseif (itytur.eq.3.and.jturb.eq.60) then
    nomflu(ir11)='fm_a_k_phase'//cphase
    nomflu(ir22)='fm_a_k_phase'//cphase
    nomflu(ir33)='fm_a_k_phase'//cphase
    nomflu(ir12)='fm_a_k_phase'//cphase
    nomflu(ir13)='fm_a_k_phase'//cphase
    nomflu(ir23)='fm_a_k_phase'//cphase
    nomflu(iep)='fm_a_omega_phase'//cphase
    if (iturb.eq.32) then
      nomflu(ial)='fm_a_omega_phase'//cphase
    endif
  elseif (itytur.eq.5.and.jtytur.eq.2) then
    nomflu(ik)='fm_a_k_phase'//cphase
    nomflu(iep)='fm_a_eps_phase'//cphase
    nomflu(iphi)='fm_a_k_phase'//cphase
    if (iturb.eq.50) then
      nomflu(ifb)='fm_a_k_phase'//cphase
    elseif (iturb.eq.51) then
      nomflu(ial)='fm_a_k_phase'//cphase
    endif
  elseif (itytur.eq.5.and.jtytur.eq.3) then
    nomflu(ik)='fm_a_R11_phase'//cphase
    nomflu(iep)='fm_a_eps_phase'//cphase
    nomflu(iphi)='fm_a_R11_phase'//cphase
    if (iturb.eq.50) then
      nomflu(ifb)='fm_a_R11_phase'//cphase
    elseif (iturb.eq.51) then
      nomflu(ial)='fm_a_R11_phase'//cphase
    endif
  elseif (iturb.eq.50.and.jturb.eq.50) then
    nomflu(ik)='fm_a_k_phase'//cphase
    nomflu(iep)='fm_a_eps_phase'//cphase
    nomflu(iphi)='fm_a_phi_phase'//cphase
    nomflu(ifb)='fm_a_fb_phase'//cphase
  elseif (iturb.eq.51.and.jturb.eq.51) then
    nomflu(ik)='fm_a_k_phase'//cphase
    nomflu(iep)='fm_a_eps_phase'//cphase
    nomflu(iphi)='fm_a_phi_phase'//cphase
    nomflu(ial)='fm_a_al_phase'//cphase
  elseif (iturb.eq.50.and.jturb.eq.51) then
    nomflu(ik)='fm_a_k_phase'//cphase
    nomflu(iep)='fm_a_eps_phase'//cphase
    nomflu(iphi)='fm_a_phi_phase'//cphase
    nomflu(ifb)='fm_a_al_phase'//cphase
  elseif (iturb.eq.51.and.jturb.eq.50) then
    nomflu(ik)='fm_a_k_phase'//cphase
    nomflu(iep)='fm_a_eps_phase'//cphase
    nomflu(iphi)='fm_a_phi_phase'//cphase
    nomflu(ial)='fm_a_fb_phase'//cphase
  elseif (itytur.eq.5.and.jturb.eq.60) then
    nomflu(ik)='fm_a_k_phase'//cphase
    nomflu(iep)='fm_a_omega_phase'//cphase
    nomflu(iphi)='fm_a_k_phase'//cphase
    if (iturb.eq.50) then
      nomflu(ifb)='fm_a_k_phase'//cphase
    elseif (iturb.eq.51) then
      nomflu(ial)='fm_a_k_phase'//cphase
    endif
  elseif (iturb.eq.60 .and. (jtytur.eq.2.or.jtytur.eq.5)) then
    nomflu(ik  )='fm_a_k_a_phase'//cphase
    nomflu(iomg)='fm_a_eps_phase'//cphase
  elseif (iturb.eq.60.and.jtytur.eq.3) then
    nomflu(ik  )='fm_a_R11_phase'//cphase
    nomflu(iomg)='fm_a_eps_phase'//cphase
  elseif (iturb.eq.60.and.jturb.eq.60) then
    nomflu(ik  )='fm_a_k_phase'//cphase
    nomflu(iomg)='fm_a_omega_phase'//cphase
  elseif (iturb.eq.70.and.jturb.eq.70) then
    nomflu(inusa)='fm_a_nusa_phase'//cphase
  endif
  if (nscal.gt.0) then
    do iscal = 1, nscal
      if (iscold(iscal).gt.0) then
        if (iscold(iscal).le.nfmtsc) then
          write(car4,'(i4.4)')iscold(iscal)
        else
          car4 = cindfs
        endif
        nomflu(isca(iscal))='fm_a_scalaire'//car4
      endif
    enddo
  endif
  if (iale.eq.1) then
    nomflu(iuma)='fm_a_vit_maill_u'
    nomflu(ivma)='fm_a_vit_maill_v'
    nomflu(iwma)='fm_a_vit_maill_w'
  endif

  ! For variables

  do ivar = 1, nvar

    f_id = ivarfl(ivar)

    ! If the variable is not associated with a mass flux, do nothing

    call field_get_key_int(f_id, kimasf, iflmas) ! interior mass flux

    if (iflmas.ge.0) then

      ! If flux has not been read yet, do it
      if (mflnum(iflmas).eq.0) then

        call field_have_previous(iflmas, lprev)

        if (.not. lprev) cycle ! skip to next loop variable

        ! Read local number of matching flux
        rubriq = nomflu(ivar)
        itysup = 0
        nbval  = 1
        irtyp  = 1
        call lecsui(impamx, rubriq, len(rubriq), itysup, nbval, irtyp,  &
                    numero, ierror)
        nberro = nberro+ierror
        ! If it exists, read it
        if (ierror.eq.0) then
          inifok = 0
          if (numero.gt.0) then
            if (numero.le.nfmtfl) then
              write(car4,'(i4.4)')numero
            else
              car4 = cindff
            endif
            if (nfaiok.eq.1) then
              call field_get_val_prev_s(iflmas, sval)
              rubriq = 'flux_masse_a_fi_'//car4
              itysup = 2
              nbval  = 1
              irtyp  = 2
              call lecsui(impamx, rubriq, len(rubriq), itysup, nbval,  &
                          irtyp, sval, ierror)
              nberro = nberro+ierror
              if (ierror.eq.0) inifok = inifok+1
            endif
            if (nfabok.eq.1) then
              call field_get_key_int(f_id, kbmasf, iflmab)
              call field_get_val_prev_s(iflmab, sval)
              rubriq = 'flux_masse_a_fb_'//car4
              itysup = 3
              nbval  = 1
              irtyp  = 2
              call lecsui(impamx, rubriq, len(rubriq), itysup, nbval,  &
                          irtyp, sval, ierror)
              nberro = nberro+ierror
              if (ierror.eq.0) inifok = inifok+1
            endif
          endif
          ! If everything is OK, mark this flux as read
          if (inifok.eq.2) then
            mflnum(iflmas) = 1
          endif
        endif
      endif
    endif
  enddo

  deallocate(mflnum)

  ! In case of error, warn but do not stop
  if (nberro.ne.0) then
    car54 = 'Lecture des flux de masse                             '
    write(nfecra,8300) car54
  endif

  car54 = ' Fin de la lecture des flux de masse                  '
  write(nfecra,1110) car54

endif
!     fin de "s'il faut lire les flux de masse (ie. supports coincidents)"

!===============================================================================
! 6. CONDITIONS AUX LIMITES
!===============================================================================

!     A ne relire que si les supports sont identiques (faces de bord)

ilu = 0

if (nfabok.eq.1) then

  ilu = 1

  nberro=0

!     Nom de variable associe a la variable dans le calcul precedent
  NOMCLI(IPR)='_p_phase'//CPHASE
  NOMCLI(IU)='_u_phase'//CPHASE
  NOMCLI(IV)='_v_phase'//CPHASE
  NOMCLI(IW)='_w_phase'//CPHASE
  !     Pour un calcul k-eps, on peut recuperer les CL d'un calcul k-eps
  !     ou d'un calcul v2f
  if (itytur.eq.2.and.                                   &
       (jtytur.eq.2.or.jtytur.eq.5)) then
    NOMCLI(IK)='_k_phase'//CPHASE
    NOMCLI(IEP)='_eps_phase'//CPHASE
  elseif (itytur.eq.3.and.jtytur.eq.3) then
    NOMCLI(IR11)='_R11_phase'//CPHASE
    NOMCLI(IR22)='_R22_phase'//CPHASE
    NOMCLI(IR33)='_R33_phase'//CPHASE
    NOMCLI(IR12)='_R12_phase'//CPHASE
    NOMCLI(IR13)='_R13_phase'//CPHASE
    NOMCLI(IR23)='_R23_phase'//CPHASE
    NOMCLI(IEP)='_eps_phase'//CPHASE
    if (iturb.eq.32.and.jturb.eq.32) then
      nomcli(ial)='_alp_phase'//CPHASE
    endif
    if (iturb.eq.32.and.jturb.ne.32) then
      nomcli(ial)='_eps_phase'//CPHASE
    endif
  elseif (itytur.eq.5.and.jtytur.eq.5) then
    NOMCLI(IK)='_k_phase'//CPHASE
    NOMCLI(IEP)='_eps_phase'//CPHASE
    NOMCLI(IPHI)='_phi_phase'//CPHASE
    if (iturb.eq.50 .and. jturb.eq.50) then
      NOMCLI(IFB)='_fb_phase'//CPHASE
    elseif (iturb.eq.51 .and. jturb.eq.51) then
      NOMCLI(IAL)='_al_phase'//CPHASE
    endif
  elseif (iturb.eq.60.and.jturb.eq.60) then
    NOMCLI(IK)='_k_phase'//CPHASE
    NOMCLI(IOMG)='_omega_phase'//CPHASE
  elseif (iturb.eq.70.and.jturb.eq.70) then
    nomcli(inusa)='_nusa_phase'//CPHASE
    !     On peut aussi recuperer les CL de k et eps pour un calcul v2f suite
    !     d'un calcul k-eps
  elseif (itytur.eq.5.and.jtytur.eq.2) then
    NOMCLI(IK)='_k_phase'//CPHASE
    NOMCLI(IEP)='_eps_phase'//CPHASE
    !     On peut aussi recuperer les CL de k pour un calcul k-eps, v2f ou k-omega suite
    !     d'un calcul k-eps, v2f ou k-omega et qui n'est pas deja un des cas ci-dessus.
  elseif ( (itytur.eq.2 .or. itytur.eq.50          &
       .or. iturb.eq.60)  .and. (jtytur.eq.2    &
       .or. jtytur.eq.5 .or. jturb.eq.60) ) then
    NOMCLI(IK)='_k_phase'//CPHASE
  endif
  if (nscal.gt.0) then
    do iscal = 1, nscal
      if (iscold(iscal).gt.0) then
        if(iscold(iscal).le.nfmtsc) then
          write(car4,'(i4.4)')iscold(iscal)
        else
          car4 = cindfs
        endif
        nomcli(isca(iscal))='_scalaire'//car4
      endif
    enddo
  endif
  if (iale.eq.1) then
    nomcli(iuma)='_vit_maillage_u'
    nomcli(ivma)='_vit_maillage_v'
    nomcli(iwma)='_vit_maillage_w'
  endif

!     --Pour les variables
  do ivar = 1, nvar

    itysup = 3
    nbval  = 1
    irtyp  = 2

!              Coefficients numeros 1
    iclvar = iclrtp(ivar,icoef)
    RUBRIQ = 'cla1'//NOMCLI(IVAR)
    call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp,     &
                coefa(1,iclvar),ierror)
    nberro=nberro+ierror

    RUBRIQ = 'clb1'//NOMCLI(IVAR)
    call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp,     &
                coefb(1,iclvar),ierror)
    nberro=nberro+ierror

!              Coefficients numeros 2
    iclvaf = iclrtp(ivar,icoeff)
    if (iclvar.ne.iclvaf) then

      RUBRIQ = 'cla2'//NOMCLI(IVAR)
      call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,         &
                  irtyp,coefa(1,iclvaf),ierror)
      nberro=nberro+ierror

      RUBRIQ = 'clb2'//NOMCLI(IVAR)
      call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,         &
                  irtyp,coefb(1,iclvaf),ierror)
      nberro=nberro+ierror
    endif

  enddo

!     Type symetrie (utilise pour les gradients par moindres carres
!       sur support etendu, avec extrapolation du gradient au bord).

  RUBRIQ = 'isympa_fb_phase'//CPHASE
  itysup = 3
  nbval  = 1
  irtyp  = 1
  call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp,isympa,ierror)
  nberro = nberro+ierror

endif
!     fin du test "si les supports des faces de bord sont identiques"

if (ilu.eq.1) then

!     Si erreur, on previent mais pas stop :
!       (on n'a pas forcement les coefs 2, on n'a pas forcement isympa
!        si on prend les fichiers d'une version anterieure)
  if (nberro.ne.0) then
    car54 =                                                       &
         'LECTURE DES CONDITIONS AUX LIMITES                    '
    write(nfecra,8300)car54
  endif

  CAR54 = ' Fin de la lecture des conditions aux limites         '
  write(nfecra,1110)car54

endif


!===============================================================================
! 7. TERMES SOURCES EXTRAPOLES
!===============================================================================


nberro=0
ilu = 0

! ---> Termes sources Navier-Stokes

!     Si le terme est a l'ordre 2
if(isno2t.gt.0) then

  iptsna = ipproc(itsnsa)
  itysup = 1
  nbval  = 1
  irtyp  = 2

  RUBRIQ = 'tsource_ns_ce_x_phase'//CPHASE
  call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp,     &
       propce(1,iptsna  ),ierror)
  nberro=nberro+ierror

  RUBRIQ = 'tsource_ns_ce_y_phase'//CPHASE
  call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp,     &
       propce(1,iptsna+1),ierror)
  nberro=nberro+ierror

  RUBRIQ = 'tsource_ns_ce_z_phase'//CPHASE
  call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp,     &
       propce(1,iptsna+2),ierror)
  nberro=nberro+ierror

  ilu = ilu + 1

endif

! ---> Termes sources turbulence

!     Si le terme est a l'ordre 2
if(isto2t.gt.0) then

  itysup = 1
  nbval  = 1
  irtyp  = 2

  !     Keps suite de keps ou v2f
  if(itytur.eq.2.and.                                    &
       (jtytur.eq.2.or.jtytur.eq.5)) then
    iptsta = ipproc(itstua)

    RUBRIQ = 'tsource_tu_ce_k_phase'//CPHASE
    call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp,   &
         propce(1,iptsta  ),ierror)
    nberro=nberro+ierror

    RUBRIQ = 'tsource_tu_ce_eps_phase'//CPHASE
    call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp,   &
         propce(1,iptsta+1),ierror)
    nberro=nberro+ierror

    ilu = ilu + 1

    !     Keps suite de Rij, on on utilise 1/2 de tsr11+tsr22+ts33
  elseif(itytur.eq.2.and.jtytur.eq.3) then

    !       Ici on veut lire tout ou rien, comme on fait des operations
    !          pour reconstruire le terme source sur k (discutable...)
    !          IDEBLU = 1 indique qu'on a commence a lire.
    !          Si erreur quand on a commence a lire, on remet tout a zero
    ideblu = 0
    iannul = 0

    iptsta = ipproc(itstua)

    RUBRIQ = 'tsource_tu_ce_R11_phase'//CPHASE
    call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp,   &
         propce(1,iptsta  ),ierror)
    nberro=nberro+ierror
    if(ierror.eq.0) ideblu = 1

    if(ideblu.eq.1) then
      RUBRIQ = 'tsource_tu_ce_R22_phase'//CPHASE
      call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp, &
           propce(1,iptsta+1),ierror)
      nberro=nberro+ierror
      if(ierror.ne.0) iannul = 1
    endif

    if(ideblu.eq.1.and.iannul.eq.0) then
      do iel = 1, ncel
        propce(iel,iptsta  ) =                                  &
             propce(iel,iptsta  )+propce(iel,iptsta+1)
      enddo
      RUBRIQ = 'tsource_tu_ce_R33_phase'//CPHASE
      call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp, &
           propce(1,iptsta+1),ierror)
      nberro=nberro+ierror
      if(ierror.ne.0) iannul = 1
    endif

    if(ideblu.eq.1.and.iannul.eq.0) then
      do iel = 1, ncel
        propce(iel,iptsta  ) = 0.5d0*(                          &
             propce(iel,iptsta  )+propce(iel,iptsta+1))
      enddo
      RUBRIQ = 'tsource_tu_ce_eps_phase'//CPHASE
      call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp, &
           propce(1,iptsta+1),ierror)
      nberro=nberro+ierror
      if(ierror.ne.0) iannul = 1
    endif

    if(iannul.eq.1) then
      do iel = 1, ncel
        propce(iel,iptsta  ) = 0.d0
        propce(iel,iptsta+1) = 0.d0
      enddo
    endif

    ilu = ilu + 1

    !     Rij suite de keps ou de v2f, on utilise 2/3 de tsk
  elseif(itytur.eq.3.and.                                &
       (jtytur.eq.2.or.jtytur.eq.50)) then
    iptsta = ipproc(itstua)

    RUBRIQ = 'tsource_tu_ce_k_phase'//CPHASE
    call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp,   &
         propce(1,iptsta  ),ierror)
    nberro=nberro+ierror

    d2s3 = 2.d0/3.d0
    do iel = 1, ncel
      tsrii = d2s3*propce(iel,iptsta  )
      propce(iel,iptsta  ) = tsrii
      propce(iel,iptsta+1) = tsrii
      propce(iel,iptsta+2) = tsrii
      propce(iel,iptsta+3) = 0.d0
      propce(iel,iptsta+4) = 0.d0
      propce(iel,iptsta+5) = 0.d0
    enddo

    RUBRIQ = 'tsource_tu_ce_eps_phase'//CPHASE
    call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp,   &
         propce(1,iptsta+6),ierror)
    nberro=nberro+ierror
    if (iturb.eq.32) then
      do iel = 1, ncel
        propce(iel,iptsta+7) = 0.d0
      enddo
    endif

    ilu = ilu + 1

    !     Rij suite de Rij
  elseif(itytur.eq.3.and.jtytur.eq.3) then
    iptsta = ipproc(itstua)

    RUBRIQ = 'tsource_tu_ce_R11_phase'//CPHASE
    call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp,   &
         propce(1,iptsta  ),ierror)
    nberro=nberro+ierror
    RUBRIQ = 'tsource_tu_ce_R22_phase'//CPHASE
    call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp,   &
         propce(1,iptsta+1),ierror)
    nberro=nberro+ierror
    RUBRIQ = 'tsource_tu_ce_R33_phase'//CPHASE
    call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp,   &
         propce(1,iptsta+2),ierror)
    nberro=nberro+ierror
    RUBRIQ = 'tsource_tu_ce_R12_phase'//CPHASE
    call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp,   &
         propce(1,iptsta+3),ierror)
    nberro=nberro+ierror
    RUBRIQ = 'tsource_tu_ce_R13_phase'//CPHASE
    call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp,   &
         propce(1,iptsta+4),ierror)
    nberro=nberro+ierror
    RUBRIQ = 'tsource_tu_ce_R23_phase'//CPHASE
    call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp,   &
         propce(1,iptsta+5),ierror)
    nberro=nberro+ierror

    RUBRIQ = 'tsource_tu_ce_eps_phase'//CPHASE
    call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp,   &
         propce(1,iptsta+6),ierror)
    nberro=nberro+ierror

    if (iturb.eq.32.and.jturb.eq.32) then
      RUBRIQ = 'tsource_tu_ce_alp_phase'//CPHASE
      call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp, &
                propce(1,iptsta+7),ierror)
      nberro=nberro+ierror
    endif

    ilu = ilu + 1

    !     v2f suite de keps : On ne relit que les TS de k et eps
    !     on laisse 0 pour phi et f_barre (valeur mise dans iniva0)
  elseif(itytur.eq.5.and.jtytur.eq.2) then
    iptsta = ipproc(itstua)

    RUBRIQ = 'tsource_tu_ce_k_phase'//CPHASE
    call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp,   &
         propce(1,iptsta  ),ierror)
    nberro=nberro+ierror

    RUBRIQ = 'tsource_tu_ce_eps_phase'//CPHASE
    call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp,   &
         propce(1,iptsta+1),ierror)
    nberro=nberro+ierror

    ilu = ilu + 1

    !     v2f suite de Rij : On ne relit que les TS de k et eps
    !     on laisse 0 pour phi et f_barre (valeur mise dans iniva0)
    !     on utilise 1/2 de tsr11+tsr22+ts33 pour le ts de k
  elseif(iturb.eq.50.and.jtytur.eq.3) then

    !       Ici on veut lire tout ou rien, comme on fait des operations
    !          pour reconstruire le terme source sur k (discutable...)
    !          IDEBLU = 1 indique qu'on a commence a lire.
    !          Si erreur quand on a commence a lire, on remet tout a zero
    ideblu = 0
    iannul = 0

    iptsta = ipproc(itstua)

    RUBRIQ = 'tsource_tu_ce_R11_phase'//CPHASE
    call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp,   &
         propce(1,iptsta  ),ierror)
    nberro=nberro+ierror
    if(ierror.eq.0) ideblu = 1

    if(ideblu.eq.1) then
      RUBRIQ = 'tsource_tu_ce_R22_phase'//CPHASE
      call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp, &
           propce(1,iptsta+1),ierror)
      nberro=nberro+ierror
      if(ierror.ne.0) iannul = 1
    endif

    if(ideblu.eq.1.and.iannul.eq.0) then
      do iel = 1, ncel
        propce(iel,iptsta  ) =                                  &
             propce(iel,iptsta  )+propce(iel,iptsta+1)
      enddo
      RUBRIQ = 'tsource_tu_ce_R33_phase'//CPHASE
      call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp, &
           propce(1,iptsta+1),ierror)
      nberro=nberro+ierror
      if(ierror.ne.0) iannul = 1
    endif

    if(ideblu.eq.1.and.iannul.eq.0) then
      do iel = 1, ncel
        propce(iel,iptsta  ) = 0.5d0*(                          &
             propce(iel,iptsta  )+propce(iel,iptsta+1))
      enddo
      RUBRIQ = 'tsource_tu_ce_eps_phase'//CPHASE
      call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp, &
           propce(1,iptsta+1),ierror)
      nberro=nberro+ierror
      if(ierror.ne.0) iannul = 1
    endif

    if(iannul.eq.1) then
      do iel = 1, ncel
        propce(iel,iptsta  ) = 0.d0
        propce(iel,iptsta+1) = 0.d0
      enddo
    endif

    ilu = ilu + 1

  elseif(itytur.eq.5.and.jtytur.eq.5) then
    iptsta = ipproc(itstua)

    RUBRIQ = 'tsource_tu_ce_k_phase'//CPHASE
    call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp,   &
         propce(1,iptsta  ),ierror)
    nberro=nberro+ierror

    RUBRIQ = 'tsource_tu_ce_eps_phase'//CPHASE
    call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp,   &
         propce(1,iptsta+1),ierror)
    nberro=nberro+ierror

    RUBRIQ = 'tsource_tu_ce_phi_phase'//CPHASE
    call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp,   &
         propce(1,iptsta+2),ierror)
    nberro=nberro+ierror

    RUBRIQ = 'tsource_tu_ce_fb_phase'//CPHASE
    call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp,   &
         propce(1,iptsta+3),ierror)
    nberro=nberro+ierror

    ilu = ilu + 1

    !     k-omega suite de k-omega
    !     Pour le k-omega, on ne relit les termes sources que
    !     si on est deja en k-omega.
    !     Il n'est en effet pas forcement tres judicieux de relire
    !     le TS de k et de mettre 0 pour le TS de eps.
    !     Quant a essayer de transformer le TS de omega en TS de eps... no comment           !
  elseif(iturb.eq.60.and.jturb.eq.60) then
    iptsta = ipproc(itstua)

    RUBRIQ = 'tsource_tu_ce_k_phase'//CPHASE
    call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp,   &
         propce(1,iptsta  ),ierror)
    nberro=nberro+ierror

    RUBRIQ = 'tsource_tu_ce_omega_phase'//CPHASE
    call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp,   &
         propce(1,iptsta+1),ierror)
    nberro=nberro+ierror

    ilu = ilu + 1
  elseif (iturb.eq.70.and.jturb.eq.70) then
    iptsta = ipproc(itstua)

    RUBRIQ = 'tsource_tu_ce_nusa_phase'//CPHASE
    call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp,   &
         propce(1,iptsta  ),ierror)
    nberro=nberro+ierror

    ilu = ilu + 1

  endif

endif

! ---> Termes sources scalaires

!     Boucle sur les scalaires
do iscal = 1, nscal
  isco = iscold(iscal)
!     Si il y a un correspondant
  if(isco.gt.0) then
!     Si ce correspondant est ok pour le format
    if(isco.le.nfmtsc) then
!     Si le terme est a l'ordre 2
      if(isso2t(iscal).gt.0) then
        iptsca = ipproc(itssca(iscal))
        itysup = 1
        nbval  = 1
        irtyp  = 2
        WRITE(CAR4,'(I4.4)')ISCO
        RUBRIQ = 'tsource_sc_ce_scalaire'//CAR4
        call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp, &
                    propce(1,iptsca  ),ierror)
        nberro=nberro+ierror
        ilu = ilu + 1
      endif
    else
      nberro=nberro-1
      ilu = ilu + 1
    endif
  endif
enddo

!     Si erreur, on previent mais pas stop :
!       (on n'a pas forcement les ts a l'ordre 2 avant)
if (nberro.ne.0) then
    car54 =                                                       &
         'LECTURE DES TERMES SOURCES                            '
    write(nfecra,8300)car54
endif

if(ilu.ne.0) then
  CAR54 =' Fin de la lecture des termes sources                 '
  write(nfecra,1110)car54
endif

!===============================================================================
! 8. MOYENNES
!===============================================================================

!     Indicateur ok (=0) ou non (>0)
iok = 0

ilu = 0

do imom = 1, nbmomt

  imold = imoold(imom)
!     Si on doit lire la moyenne
  if(imold.gt.0) then
    ilu = 1
!     Si ce correspondant est ok pour le format
    if(imold.le.nfmtmo) then
      WRITE(CAR4,'(I4.4)')IMOLD
    else
      car4 = cindfm
    endif

!       On la lit
    itysup = 1
    nbval  = 1
    irtyp  = 2
    RUBRIQ = 'cumul_ce_moment'//CAR4
    call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp,     &
                propce(1,ipproc(icmome(imom))),ierror)
    nberro=nberro+ierror

!       Si ca ne marche pas, on s'arrete
!         (on pourrait tenter de prendre des mesures correctives,
!          mais le cumul de la duree associee sert peut etre a d'autres
!            moyennes
!          en outre, si l'utilisateur a indique qu'il voulait relire
!            une moyenne, c'est qu'il veut faire un calcul propre
    if(ierror.ne.0) then
      write(nfecra,9300)imold,imom,imom
      iok = iok + 1
    endif



!       Si on a reussi a lire la moyenne, il faut obligatoirement disposer
!         du cumul de duree associe, sinon, on devra s'arreter.
!       IOKDT different de 0 indiquera qu'on n'a pas pu l'obtenir.
    iokdt = 0

!       On cherche le numero (local au fichier suite)
!         du cumul de temps de la moyenne du calcul precedent
!         correspondant a la moyenne IMOM du calcul courant (ouf          !)
    itysup = 0
    nbval  = 1
    irtyp  = 1
    RUBRIQ = 'numero_cumul_temps_moment'//CAR4
    numero = 0
    call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp,     &
                numero,ierror)
    nberro=nberro+ierror

!       Si on n'a pas trouve, on a echoue
    if(numero.eq.0.or.ierror.ne.0) then
      iokdt = 1

!       Sinon, on a trouve un cumul en temps correspondant
    else

!         Si NUMERO > 0, il s'agissait d'un cumul variable en espace
!             et le cumul courant est forcement variable en espace
      if(numero.gt.0) then

        if(numero.le.nfmtmo) then
          WRITE(CNUM4,'(I4.4)')NUMERO
        else
          cnum4 = cindfm
        endif

        itysup = 1
        nbval  = 1
        irtyp  = 2
        RUBRIQ = 'cumul_temps_ce_'//CNUM4
        idtcm  = ipproc(icdtmo(idtmom(imom)))
        call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp, &
                    propce(1,idtcm),ierror)
        nberro=nberro+ierror

!           Si on n'a pas lu, on a echoue
        if(ierror.ne.0) then
          iokdt = 1
        endif

!         Si NUMERO < 0, il s'agissait d'un cumul uniforme
!             et le cumul courant est variable en espace ou non
      elseif(numero.lt.0) then

        numero = -numero
        if(numero.le.nfmtmo) then
          WRITE(CNUM4,'(I4.4)')NUMERO
        else
          cnum4 = cindfm
        endif

        itysup = 0
        nbval  = 1
        irtyp  = 2
        RUBRIQ = 'cumul_temps_'//CNUM4
        cdtcm  = 0.d0
        call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp, &
                    cdtcm,ierror)
        nberro=nberro+ierror

!           Si on n'a pas lu, on a echoue
        if(ierror.ne.0) then
          iokdt = 1

!           Sinon, selon la nature du cumul actuel, on l'affecte
        else
          if(idtmom(imom).gt.0) then
            idtcm  = ipproc(icdtmo(idtmom(imom)))
            do iel = 1, ncel
              propce(iel,idtcm) = cdtcm
            enddo
          elseif(idtmom(imom).lt.0) then
            idtcm  = -idtmom(imom)
            dtcmom(idtcm) = cdtcm
          endif
        endif

      endif

    endif

    if(iokdt.ne.0) then
      write(nfecra,9310)imold,imom,imom
      iok = iok + 1
    endif

  endif
enddo

!     Si pb on s'arrete car si on cherche a faire des moyennes
!       c'est qu'on veut vraiement les relire
if(iok.ne.0) then
  call csexit(1)
endif

if(ilu.gt.0) then
  CAR54 = ' Fin de la lecture des moyennes temporelles           '
  write(nfecra,1110)car54
endif

!===============================================================================
! 9. DISTANCE A LA PAROI
!===============================================================================

ilu = 0
nberro = 0

!     MODE DE CALCUL DIRECT (NON COMPATIBLE PARALLELISME ET PERIODICITE)

if(abs(icdpar).eq.2) then

!     On la lit si on en a besoin uniquement.

!     Si l'utilisateur a force le recalcul, on ne la lit pas
!       il faudra la mettre a jour (sauf si zero pas de temps).

!     Sinon, on cherche a la lire.

!     On ne relit les numeros des faces de bord que si on a toujours
!       le meme nombre de faces de bord (sinon, la numerotation a change)
  if(ineedy.eq.1) then
    if(icdpar.eq.2.or.inpdt0.eq.1) then
      if(nfabok.eq.1) then

        itysup = 1
        nbval  = 1
        irtyp  = 1
        RUBRIQ = 'num_fac_par_ce_phase'//CPHASE
        call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,   &
             irtyp,ifapat,ierror)
        nberro=nberro+ierror
        ilu   = ilu + 1

      endif
    endif
  endif


!     MODE DE CALCUL PAR EQUATION DE DIFFUSION

elseif(abs(icdpar).eq.1) then

!     On la lit si on en a besoin uniquement.

!     Si l'utilisateur a force le recalcul, on ne la lit pas
!       il faudra la mettre a jour (sauf si zero pas de temps).

!     Sinon, on cherche a la lire.
!       On pourrait la relire aussi quand le nombre de faces a
!       change, mais il vaut mieux la recalculer au cas ou des faces de
!       paroi auraient disparu
!       Si on arrive a la lire, on note qu'elle est a jour (sauf si ALE).

  if(ineedy.eq.1) then
    if(icdpar.eq.1.or.inpdt0.eq.1) then
      if(nfabok.eq.1) then
        itysup = 1
        nbval  = 1
        irtyp  = 2
        RUBRIQ = 'dist_fac_par_ce_phase'//CPHASE
        call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp, &
                    dispar,ierror)
        nberro=nberro+ierror
        if(ierror.eq.0 .and. iale.eq.0 ) then
          imajdy = 1
        endif
        ilu   = ilu + 1
      endif
    endif
  endif

endif

if (nberro.ne.0) then
  car54 =                                                         &
         'LECTURE DE LA DISTANCE A LA PAROI                     '
  write(nfecra,8300)car54
endif

if(ilu.ne.0) then
  CAR54=' Fin de la lecture de la distance a la paroi          '
  write(nfecra,1110)car54
endif


!===============================================================================
! 10.  FORCE EXTERIEURE
!===============================================================================

if(iphydr.eq.1) then
  nberro=0

  itysup = 1
  nbval  = 3
  irtyp  = 2

  ! TODO read the old format
  rubriq = 'force_ext_ce_phase'//cphase
  call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp,     &
       frcxt,ierror)
  nberro=nberro+ierror

 if (nberro.ne.0) then
    car54 =                                                       &
         'LECTURE DES FORCES EXTERIEURES                        '
    write(nfecra,8300)car54
  endif

  CAR54 =' Fin de la lecture des forces exterieures             '
  write(nfecra,1110)car54

endif

!===============================================================================
! 11. PRESSION HYDROSTATIQUE PREDITE
!===============================================================================

if(iphydr.eq.2) then
  nberro=0

  itysup = 1
  nbval  = 1
  irtyp  = 2

  RUBRIQ = 'prhyd_pre_phase'//CPHASE
  call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp,     &
       prhyd(1),ierror)
  nberro=nberro+ierror

 if (nberro.ne.0) then
    car54 =                                                       &
         'LECTURE DE LA PRESSION HYDROSTATIQUE PREDITE          '
    write(nfecra,8300)car54
  endif

  CAR54 =' Fin de la lecture de la pression hydro. predite      '
  write(nfecra,1110)car54

endif

!===============================================================================
! 12.  DEPLACEMENT AUX NOEUDS EN ALE
!===============================================================================

if (iale.eq.1 .and. jale.eq.1) then
  nberro = 0

  itysup = 4

  call restart_read_real_3_t_compat                       &
         (impamx, 'vertex_displacement',                  &
         'deplact_x_no', 'deplact_y_no', 'deplact_z_no',  &
         itysup, depale, ierror)

  nberro=nberro+ierror

! Si JALE=1, on doit avoir le deplacement dans le fichier suite, sinon
!   les resultats relus n'ont pas de sens -> on s'arrete si pb
  if (nberro.ne.0) then
    write(nfecra,9320)
    call csexit(1)
  endif

  nberro = 0
  RUBRIQ = 'type_visc_mail'
  itysup = 0
  nbval  = 1
  irtyp  = 1
  call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp,       &
       jortvm,ierror)

  RUBRIQ = 'visc_maillage_x'
  itysup = 1
  nbval  = 1
  irtyp  = 2
  call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp,       &
       propce(1,ipproc(ivisma(1))),ierror)

  if (iortvm.eq.1) then
    if (jortvm.eq.1) then
      RUBRIQ = 'visc_maillage_y'
      call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp,   &
           propce(1,ipproc(ivisma(2))),ierror)
      RUBRIQ = 'visc_maillage_z'
      call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp,   &
           propce(1,ipproc(ivisma(3))),ierror)
    else
      do iel = 1, ncel
        ipcvmx = ipproc(ivisma(1))
        ipcvmy = ipproc(ivisma(2))
        ipcvmz = ipproc(ivisma(3))
        propce(iel,ipcvmy) = propce(iel,ipcvmx)
        propce(iel,ipcvmz) = propce(iel,ipcvmx)
      enddo
    endif
  endif

  CAR54 =' Fin de la lecture des donnees ALE                    '
  write(nfecra,1110)car54

  nberro=0
  RUBRIQ = 'nombre_structures'
  itysup = 0
  nbval  = 2
  irtyp  = 1
  call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp,       &
       ngbstr,ierror)
  nberro=nberro+ierror

  nbstru = ngbstr(1)
  nbaste = ngbstr(2)

  if (nbstru.gt.0) then

    nfmtst = 99
    CINDST='YY'
    do istr = 1, min(nbstru,nstrmx)
      WRITE(CSTRUC(ISTR),'(I2.2)') ISTR
    enddo
    do istr = min(nbstru,nfmtst)+1,nbstru
      cstruc(istr) = cindst
    enddo
    if(nstrmx.gt.nfmtst) then
      write(nfecra,8004)nfmtst,nstrmx
    endif

    do istr = 1, nbstru

      RUBRIQ = 'donnees_structure_'//CSTRUC(ISTR)
      itysup = 0
      nbval  = 27
      irtyp  = 2

      call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp,   &
           tmpstr,ierror)
      nberro=nberro+ierror

      do ii = 1, 3
        xstr  (ii,istr) = tmpstr(   ii)
        xpstr (ii,istr) = tmpstr(3 +ii)
        xppstr(ii,istr) = tmpstr(6 +ii)
        xsta  (ii,istr) = tmpstr(9 +ii)
        xpsta (ii,istr) = tmpstr(12+ii)
        xppsta(ii,istr) = tmpstr(15+ii)
        xstp  (ii,istr) = tmpstr(18+ii)
        forstr(ii,istr) = tmpstr(21+ii)
        forsta(ii,istr) = tmpstr(24+ii)
      enddo

    enddo

    CAR54 =' Fin de la lecture des donnees des structures ALE   '
    write(nfecra,1110)car54

  endif

  if (nberro.ne.0) then
    write(nfecra,9321)
    call csexit(1)
  endif

endif

!===============================================================================
! 13. LECTURE DES INFORMATIONS COMPLEMENTAIRES COMBUSTION GAZ, CP ET
!                                                                  FUEL
!===============================================================================

nberro = 0
ilu = 0

!     Modele COD3P :
!     ============

if ( ippmod(icod3p).ge.0 ) then

  RUBRIQ = 'hinfue_cod3p'
  itysup = 0
  nbval  = 1
  irtyp  = 2
  call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp,hinfue,&
              ierror)
  nberro=nberro+ierror
  ilu = ilu + 1
  if (ierror.ne.0) then
    write(nfecra,9400)
  endif

  RUBRIQ = 'hinoxy_cod3p'
  itysup = 0
  nbval  = 1
  irtyp  = 2
  ilu    = ilu + 1
  call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp,hinoxy,&
              ierror)
  nberro=nberro+ierror
  ilu = ilu + 1
  if (ierror.ne.0) then
    write(nfecra,9400)
  endif

  RUBRIQ = 'tinfue_cod3p'
  itysup = 0
  nbval  = 1
  irtyp  = 2
  ilu    = ilu + 1
  call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp,tinfue,&
              ierror)
  nberro=nberro+ierror
  ilu = ilu + 1
  if (ierror.ne.0) then
    write(nfecra,9400)
  endif

  RUBRIQ = 'tinoxy_cod3p'
  itysup = 0
  nbval  = 1
  irtyp  = 2
  ilu    = ilu + 1
  call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp,tinoxy,&
              ierror)
  nberro=nberro+ierror
  ilu = ilu + 1
  if (ierror.ne.0) then
    write(nfecra,9400)
  endif

!       Il faut le meme nbr de faces de bord, sinon on ne lit pas
  if(nfabok.eq.1) then

    ilu = ilu + 1

    ierrch = 0

!       Numero des zones
    itysup = 3
    nbval  = 1
    irtyp  = 1
    RUBRIQ = 'num_zone_fb_cod3p'
    call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp,     &
                izfppp, ierror)
    nberro=nberro+ierror

!       Type entree Fuel
    itysup = 0
    nbval  = nozppm
    irtyp  = 1
    RUBRIQ = 'ientfu_zone_bord_cod3p'
    call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp,     &
                ientfu, ierror)
    ierrch=ierrch+ierror
    nberro=nberro+ierror

!       Type entree Oxydant
    itysup = 0
    nbval  = nozppm
    irtyp  = 1
    RUBRIQ = 'ientox_zone_bord_cod3p'
    call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp,     &
                ientox, ierror)
    ierrch=ierrch+ierror
    nberro=nberro+ierror

!     Par securite, si on ne parvient pas a lire la
!       IENTCPFU ou IENTOX, on remet a zero le numero des zones IZFPPP
!       car il a peut etre ete lu.
!       Ceci permettra d'eviter de se servir des valeurs par defaut

    if(ierrch.ne.0) then
      do ifac = 1, nfabor
        izfppp(ifac) = 0
      enddo
    endif

  endif

endif

!     Modele EBU :
!     ==========

if ( ippmod(icoebu).ge.0 ) then

  RUBRIQ = 'temperature_gaz_frais_ebu'
  itysup = 0
  nbval  = 1
  irtyp  = 2
  ilu    = 1
  call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp,tgf,   &
              ierror)
  nberro=nberro+ierror
  if (ierror.ne.0) then
    write(nfecra,9500)
  endif

  RUBRIQ = 'frmel_ebu'
  itysup = 0
  nbval  = 1
  irtyp  = 2
  ilu    = 1
  call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp,frmel, &
              ierror)
  nberro=nberro+ierror
  if (ierror.ne.0) then
    write(nfecra,9500)
  endif

!       Il faut le meme nbr de faces de bord, sinon on ne lit pas
  if(nfabok.eq.1) then

    ilu = ilu + 1

    ierrch = 0

!       Numero des zones
    itysup = 3
    nbval  = 1
    irtyp  = 1
    RUBRIQ = 'num_zone_fb_ebu'
    call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp,     &
                izfppp, ierror)
    nberro=nberro+ierror

!       Type entree Gaz brulee
    itysup = 0
    nbval  = nozppm
    irtyp  = 1
    RUBRIQ = 'ientgb_zone_bord_ebu'
    call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp,     &
                ientgb, ierror)
    ierrch=ierrch+ierror
    nberro=nberro+ierror

!       Type entree gaz frais
    itysup = 0
    nbval  = nozppm
    irtyp  = 1
    RUBRIQ = 'ientgf_zone_bord_ebu'
    call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp,     &
                ientgf, ierror)
    ierrch=ierrch+ierror
    nberro=nberro+ierror

!       FMENT
    itysup = 0
    nbval  = nozppm
    irtyp  = 2
    RUBRIQ = 'fment_zone_bord_ebu'
    call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp,     &
                fment , ierror)
    ierrch=ierrch+ierror
    nberro=nberro+ierror

!       TKENT
    itysup = 0
    nbval  = nozppm
    irtyp  = 2
    RUBRIQ = 'tkent_zone_bord_ebu'
    call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp,     &
                tkent, ierror)
    ierrch=ierrch+ierror
    nberro=nberro+ierror

!     Par securite, si on ne parvient pas a lire la
!       IENTCPFU ou IENTOX, on remet a zero le numero des zones IZFPPP
!       car il a peut etre ete lu.
!       Ceci permettra d'eviter de se servir des valeurs par defaut

    if(ierrch.ne.0) then
      do ifac = 1, nfabor
        izfppp(ifac) = 0
      enddo
    endif

  endif

endif

!     Modele LWC :
!     ==========

if ( ippmod(icolwc).ge.0 ) then

  RUBRIQ = 'fmin_lwc'
  itysup = 0
  nbval  = 1
  irtyp  = 2
  ilu    = 1
  call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp,fmin,  &
              ierror)
  nberro=nberro+ierror
  if (ierror.ne.0) then
    write(nfecra,9600)
  endif

  RUBRIQ = 'fmax_lwc'
  itysup = 0
  nbval  = 1
  irtyp  = 2
  ilu    = 1
  call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp,fmax,  &
              ierror)
  nberro=nberro+ierror
  if (ierror.ne.0) then
    write(nfecra,9600)
  endif

  RUBRIQ = 'hmin_lwc'
  itysup = 0
  nbval  = 1
  irtyp  = 2
  ilu    = 1
  call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp,hmin,  &
              ierror)
  nberro=nberro+ierror
  if (ierror.ne.0) then
    write(nfecra,9600)
  endif

  RUBRIQ = 'hmax_lwc'
  itysup = 0
  nbval  = 1
  irtyp  = 2
  ilu    = 1
  call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp,hmax,  &
              ierror)
  nberro=nberro+ierror
  if (ierror.ne.0) then
    write(nfecra,9600)
  endif

!       Il faut le meme nbr de faces de bord, sinon on ne lit pas
  if(nfabok.eq.1) then

    ilu = ilu + 1

    ierrch = 0

!       Numero des zones
    itysup = 3
    nbval  = 1
    irtyp  = 1
    RUBRIQ = 'num_zone_fb_lwc'
    call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp,     &
                izfppp, ierror)
    nberro=nberro+ierror

!       Type entree Gaz brulee
    itysup = 0
    nbval  = nozppm
    irtyp  = 1
    RUBRIQ = 'ientgb_zone_bord_lwc'
    call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp,     &
                ientgb, ierror)
    ierrch=ierrch+ierror
    nberro=nberro+ierror

!       Type entree gaz frais
    itysup = 0
    nbval  = nozppm
    irtyp  = 1
    RUBRIQ = 'ientgf_zone_bord_lwc'
    call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp,     &
                ientgf, ierror)
    ierrch=ierrch+ierror
    nberro=nberro+ierror

!       FMENT
    itysup = 0
    nbval  = nozppm
    irtyp  = 2
    RUBRIQ = 'fment_zone_bord_lwc'
    call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp,     &
                fment , ierror)
    ierrch=ierrch+ierror
    nberro=nberro+ierror

!       TKENT
    itysup = 0
    nbval  = nozppm
    irtyp  = 2
    RUBRIQ = 'tkent_zone_bord_lwc'
    call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp,     &
                tkent, ierror)
    ierrch=ierrch+ierror
    nberro=nberro+ierror

!     Par securite, si on ne parvient pas a lire la
!       IENTCPFU ou IENTOX, on remet a zero le numero des zones IZFPPP
!       car il a peut etre ete lu.
!       Ceci permettra d'eviter de se servir des valeurs par defaut

    if(ierrch.ne.0) then
      do ifac = 1, nfabor
        izfppp(ifac) = 0
      enddo
    endif

  endif

endif

!     Charbon PuLVerise : masse vol des charbons
if (ippmod(icpl3c).ge.0 .or.                                      &
    ippmod(iccoal).ge.0) then
  itysup = 0
  nbval  = 1
  irtyp  = 2
  ierrch = 0
  do icha = 1, ncharb
    if(icha.le.nfmtch) then
      WRITE(CAR2,'(I2.2)')ICHA
    else
      car2 = cindfc
    endif
    RUBRIQ = 'masse_volumique_charbon'//CAR2
    call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp,     &
                rhock(icha), ierror)
    ierrch = ierrch + ierror
    nberro = nberro + ierror
    ilu = ilu + 1
  enddo
  if (ierrch.ne.0) then
    write(nfecra,8611)
    do icha = 1, ncharb
      write(nfecra,8612)icha,rhock(icha)
    enddo
    write(nfecra,8613)
  endif


!     Charbon PuLVerise : type de zones de bord, ientat, ientcp, timpat
!       et x20 pour le calcul de rho au bord en entree
!       Il faut le meme nbr de faces de bord, sinon on ne lit pas
  if(nfabok.eq.1) then

    ilu = ilu + 1

    ierrch = 0

!       Numero des zones
    itysup = 3
    nbval  = 1
    irtyp  = 1
    RUBRIQ = 'num_zone_fb_charbon_pulverise'
    call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp,     &
                izfppp, ierror)
    nberro = nberro + ierror

!       Type entree air ou cp (si ce n'est pas NOZPPM, erreur)
    itysup = 0
    nbval  = nozppm
    irtyp  = 1
    RUBRIQ = 'ientat_zone_bord_charbon_pulverise'
    call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp,     &
                ientat, ierror)
    ierrch = ierrch + ierror
    nberro = nberro + ierror

!         ientcp et x20 ne servent pas pour le CP couple Lagrangien (cplphy)
    if (ippmod(iccoal).ge.0) then

      itysup = 0
      nbval  = nozppm
      irtyp  = 1
      RUBRIQ = 'ientcp_zone_bord_charbon_pulverise'
      call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp,   &
                  ientcp, ierror)
      ierrch = ierrch + ierror
      nberro = nberro + ierror

      itysup = 0
      nbval  = nozppm
      irtyp  = 1
      RUBRIQ = 'inmoxy_zone_bord_charbon_pulverise'
      call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp,   &
                  inmoxy, ierror)
      ierrch = ierrch + ierror
      nberro = nberro + ierror

      itysup = 0
      nbval  = nozppm
      irtyp  = 2

      idecal = 0
      do icha = 1, ncharb
        do iclapc = 1, nclpch(icha)
          icla = iclapc + idecal
          if(icha.le.nfmtch.and.iclapc.le.nfmtcl) then
            WRITE(CAR2,'(I2.2)')ICHA
            WRITE(CAR4,'(I4.4)')ICLAPC
          else
            car2 = cindfc
            car4 = cindfl
          endif
          RUBRIQ = 'x20_zone_bord_charbon'//CAR2//'_classe'//CAR4
          call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,     &
               irtyp,x20(1,icla), ierror)
          ierrch = ierrch + ierror
          nberro = nberro + ierror

        enddo
      enddo

    endif

!       Temperature
    itysup = 0
    nbval  = nozppm
    irtyp  = 2
    RUBRIQ = 'timpat_zone_bord_charbon_pulverise'
    call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp,     &
                timpat, ierror)
    ierrch = ierrch + ierror
    nberro = nberro + ierror

!     Par securite, si on ne parvient pas a lire la temperature TIMPAT,
!       IENTCP ou IENTAT, on remet a zero le numero des zones IZFPPP
!       car il a peut etre ete lu.
!       Ceci permettra d'eviter de se servir des valeurs par defaut (=0)
!       de TIMPAT dans cpphyv et cplphy.
    if(ierrch.ne.0) then
      do ifac = 1, nfabor
        izfppp(ifac) = 0
      enddo
    endif

  endif

endif


!     FUEL : type de zones de bord, ientat, ientfl, timpat
!       qimpat et qimpfl pour le calcul de rho au bord en entree
if ( ippmod(icfuel).ge.0 ) then

!       Il faut le meme nbr de faces de bord, sinon on ne lit pas
  if(nfabok.eq.1) then

    ilu = ilu + 1

    ierrch = 0

!       Numero des zones
    itysup = 3
    nbval  = 1
    irtyp  = 1
    RUBRIQ = 'num_zone_fb_fuel'
    call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp,     &
                izfppp, ierror)
    nberro=nberro+ierror

!       Type entree air ou fuel (si ce n'est pas NOZPPM, erreur)
    itysup = 0
    nbval  = nozppm
    irtyp  = 1
    RUBRIQ = 'ientat_zone_bord_fuel'
    call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp,     &
                ientat, ierror)
    ierrch=ierrch+ierror
    nberro=nberro+ierror

    itysup = 0
    nbval  = nozppm
    irtyp  = 1
    RUBRIQ = 'ientfl_zone_bord_fuel'
    call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp,     &
                ientfl, ierror)
    ierrch=ierrch+ierror
    nberro=nberro+ierror

!
    itysup = 0
    nbval  = nozppm
    irtyp  = 1
    RUBRIQ = 'inmoxy_zone_bord_fuel'
    call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp,     &
                inmoxy, ierror)
    ierrch=ierrch+ierror
    nberro=nberro+ierror


!       TIMPAT
    itysup = 0
    nbval  = nozppm
    irtyp  = 2
    RUBRIQ = 'timpat_zone_bord_fuel'
    call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp,     &
                timpat, ierror)
    ierrch=ierrch+ierror
    nberro=nberro+ierror

!       QIMPAT
    itysup = 0
    nbval  = nozppm
    irtyp  = 2
    RUBRIQ = 'qimpat_zone_bord_fuel'
    call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp,     &
                qimpat, ierror)
    ierrch=ierrch+ierror
    nberro=nberro+ierror

!       QIMPFL
    itysup = 0
    nbval  = nozppm
    irtyp  = 2
    RUBRIQ = 'qimpfl_zone_bord_fuel'
    call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp,     &
                qimpfl, ierror)
    ierrch=ierrch+ierror
    nberro=nberro+ierror

!     Par securite, si on ne parvient pas a lire la temperature TIMPAT,
!       IENTCP ou IENTAT, on remet a zero le numero des zones IZFPPP
!       car il a peut etre ete lu.
!       Ceci permettra d'eviter de se servir des valeurs par defaut (=0)
!       de TIMPAT dans cpphyv et cplphy.
    if(ierrch.ne.0) then
      do ifac = 1, nfabor
        izfppp(ifac) = 0
      enddo
    endif

  endif

endif

if (nberro.ne.0) then
  car54 =                                                         &
       'LECTURE DES INFORMATIONS COMBUSTION                   '
  write(nfecra,8300)car54
endif

if(ilu.ne.0) then
  CAR54=' Fin de la lecture des informations combustion        '
  write(nfecra,1110)car54
endif

!===============================================================================
! 14. LECTURE DES INFORMATIONS COMPLEMENTAIRES ELECTRIQUES
!===============================================================================

nberro=0
ilu  = 0

!     Recalage des CL pot des versions electriques

if ( ippmod(ieljou).ge.1       ) then
  if(ielcor.eq.1) then
    ilu = ilu + 1
    RUBRIQ = 'coeff_recalage_joule'
    itysup = 0
    nbval  = 1
    irtyp  = 2
    call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp,     &
                coejou,ierror)
    nberro=nberro+ierror
  endif
endif
if ( ippmod(ielarc).ge.1  .or. ippmod(ieljou).ge.1 ) then
  if(ielcor.eq.1) then
    ilu = 1
    RUBRIQ = 'ddpot_recalage_arc_elec'
    itysup = 0
    nbval  = 1
    irtyp  = 2
    call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp,     &
                dpot  ,ierror)
    nberro=nberro+ierror
  endif
endif

! ---> Termes sources des versions electriques

if ( ippmod(ieljou).ge.1 .or.                                     &
     ippmod(ielarc).ge.1 .or.                                     &
     ippmod(ielion).ge.1       ) then

  ipcefj = ipproc(iefjou)
  itysup = 1
  nbval  = 1
  irtyp  = 2

  RUBRIQ = 'tsource_sc_ce_joule'
  call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp,       &
              propce(1,ipcefj  ),ierror)
  nberro=nberro+ierror
  ilu = ilu + 1

endif

if( ippmod(ielarc).ge.1 ) then

  ipcla1 = ipproc(ilapla(1))
  ipcla2 = ipproc(ilapla(2))
  ipcla3 = ipproc(ilapla(3))
  itysup = 1
  nbval  = 1
  irtyp  = 2

  RUBRIQ = 'tsource_ns_ce_x_laplace'
  call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp,       &
              propce(1,ipcla1  ),ierror)
  nberro=nberro+ierror
  ilu = ilu + 1

  RUBRIQ = 'tsource_ns_ce_y_laplace'
  call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp,       &
              propce(1,ipcla2  ),ierror)
  nberro=nberro+ierror
  ilu = ilu + 1

  RUBRIQ = 'tsource_ns_ce_z_laplace'
  call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp,       &
              propce(1,ipcla3  ),ierror)
  nberro=nberro+ierror
  ilu = ilu + 1

endif

if (nberro.ne.0) then
  car54 =                                                         &
       'LECTURE DES INFORMATIONS ELECTRIQUES                  '
  write(nfecra,8300)car54
endif

if(ilu.ne.0) then
  CAR54=' Fin de la lecture des informations electriques       '
  write(nfecra,1110)car54
endif

!===============================================================================
! 15.  FERMETURE DU FICHIER SUITE AUXILAIRE
!===============================================================================


!     Fermeture du fichier suite auxilaire
call clssui(impamx,ierror)

if (ierror.ne.0) then
   write(nfecra,8900) ficsui
endif

write(nfecra,1200)

!===============================================================================
! 16. SORTIE
!===============================================================================

return

!===============================================================================
! 17. FORMATS
!===============================================================================

! --- ETAPES

#if defined(_CS_LANG_FR)

 1000 format(/,                                                   &
     3X,'   LECTURE DU FICHIER SUITE AUXILIAIRE               ',/)
 1100 format(' Debut de la lecture                                    ')
 1110 format('  ',A54                                                  )
 1200 format(' Fin de la lecture                                      ')

#else

 1000 format(/,                                                   &
     3X,'      READING THE AUXILIARY RESTART FILE             ',/)
 1100 format(' Start reading                                          ')
 1110 format('  ',A54                                                  )
 1200 format(' End reading                                            ')

#endif

! --- INFORMATIONS

#if defined(_CS_LANG_FR)

 7000 format(/,                                                   &
'   Mise a jour du point de reference pour la pression totale ',/,&
'     par relecture du fichier suite                          ',/,&
'    XYZP0 = ',       E14.5,        E14.5,        E14.5        ,/)

#else

 7000 format(/,                                                   &
'   Apdatation of the reference point for the total pressure  ',/,&
'       by reading the restart file                           ',/,&
'    XYZP0 = ',       E14.5,        E14.5,        E14.5        ,/)

#endif

! --- MISES EN GARDE

#if defined(_CS_LANG_FR)

 8001 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION :       A LA LECTURE DU FICHIER SUITE         ',/,&
'@    =========                                     AUXILIAIRE',/,&
'@                                                            ',/,&
'@      Le nombre de scalaires maximal NSCAMX supporte par le ',/,&
'@        format d''ecriture du fichier suite est             ',/,&
'@        NFMTSC = ',I10                                       ,/,&
'@      On a ici un nombre de scalaires maximal superieur     ',/,&
'@        NSCAMX = ',I10                                       ,/,&
'@      On ne pourra pas relire les scalaires dont le numero  ',/,&
'@        est superieur                                       ',/,&
'@                                                            ',/,&
'@    Le calcul sera execute.                                 ',/,&
'@                                                            ',/,&
'@    Voir le sous-programme lecamx.                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8002 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION :       A LA LECTURE DU FICHIER SUITE         ',/,&
'@    =========                                     AUXILIAIRE',/,&
'@                                                            ',/,&
'@      Le nombre de flux de masse max NVARMX supporte par le ',/,&
'@        format d''ecriture du fichier suite est             ',/,&
'@        NFMTFL = ',I10                                       ,/,&
'@      On a ici un nombre de flux      maximal superieur     ',/,&
'@        NVARMX = ',I10                                       ,/,&
'@      On ne pourra pas relire les flux      dont le numero  ',/,&
'@        est superieur                                       ',/,&
'@                                                            ',/,&
'@    Le calcul sera execute.                                 ',/,&
'@                                                            ',/,&
'@    Voir le sous-programme lecamx.                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8003 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION :       A LA LECTURE DU FICHIER SUITE         ',/,&
'@    =========                                     AUXILIAIRE',/,&
'@                                                            ',/,&
'@      Le nombre de moments mx        NBMOMX supporte par le ',/,&
'@        format d''ecriture du fichier suite est             ',/,&
'@        NFMTMO = ',I10                                       ,/,&
'@      On a ici un nombre de moments   maximal superieur     ',/,&
'@        NBMOMX = ',I10                                       ,/,&
'@      On ne pourra pas relire les moments   dont le numero  ',/,&
'@        est superieur                                       ',/,&
'@                                                            ',/,&
'@    Le calcul sera execute.                                 ',/,&
'@                                                            ',/,&
'@    Voir le sous-programme lecamx.                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8004 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION :       A LA LECTURE DU FICHIER SUITE         ',/,&
'@    =========                                     AUXILIAIRE',/,&
'@                                                            ',/,&
'@      Le nombre de structures maximal NSTRMX supporte par le',/,&
'@        format d''ecriture du fichier suite est             ',/,&
'@        NFMTST = ',I10                                       ,/,&
'@      On a ici un nombre de structures maximal superieur    ',/,&
'@        NSTRMX = ',I10                                       ,/,&
'@      Si le nombre de structures effectif est superieur,    ',/,&
'@        elles ne seront pas relues.                         ',/,&
'@                                                            ',/,&
'@    Le calcul sera execute.                                 ',/,&
'@                                                            ',/,&
'@    Voir le sous-programme lecamx.                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8200 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : LECTURE DU FICHIER SUITE AUXILIAIRE         ',/,&
'@    =========                                               ',/,&
'@      DONNEES AMONT ET ACTUELLES DIFFERENTES                ',/,&
'@                                                            ',/,&
'@    Le nombre de faces ',A8  ,' a ete modifie.              ',/,&
'@                                                            ',/,&
'@    Le calcul peut etre execute mais les donnees            ',/,&
'@      sur les faces ',A8  ,' ne seront pas relues           ',/,&
'@      dans le fichier suite.                                ',/,&
'@    Elles seront initialisees par des valeurs par defaut.   ',/,&
'@                                                            ',/,&
'@    Cette situation peut se produire lorsque le fichier     ',/,&
'@      suite est issu d''un calcul realise avec des options  ',/,&
'@      de recollement differentes ou lorsque l''on modifie   ',/,&
'@      la prise en compte de periodicite.                    ',/,&
'@    Cette situation peut egalement se produire lorsque l''on',/,&
'@      realise une suite sur une machine de calcul differente',/,&
'@      et que le jeu de la precision machine modifie le      ',/,&
'@      nombre de faces issues des recollements.              ',/,&
'@                                                            ',/,&
'@    Cette situation peut enfin se produire lorsque le       ',/,&
'@      fichier suite auxiliaire ne correspond pas au cas     ',/,&
'@      traite.                                               ',/,&
'@                                                            ',/,&
'@    Verifier que le fichier suite auxiliaire utilise        ',/,&
'@      correspond bien au cas traite                         ',/,&
'@                                                            ',/,&
'@    Le calcul se poursuit...                                ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8205 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : LECTURE DU FICHIER SUITE AUXILIAIRE         ',/,&
'@    =========                                               ',/,&
'@      DONNEES AMONT MULTIPHASIQUES                          ',/,&
'@                                                            ',/,&
'@  Nombre de phases (amont) : ',I10                           ,/,&
'@                                                            ',/,&
'@    Le calcul ne peut etre execute.                         ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8220 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : LECTURE DU FICHIER SUITE AUXILIAIRE         ',/,&
'@    =========                                               ',/,&
'@      REPRISE  DE CALCUL           AVEC ITURB = ',I4         ,/,&
'@      A PARTIR D''UN CALCUL REALISE AVEC ITURB = ',I4        ,/,&
'@                                                            ',/,&
'@    Le modele de turbulence a ete modifie.                  ',/,&
'@                                                            ',/,&
'@    Il est conseille cependant de                           ',/,&
'@      verifier la valeur de ITURB                           ',/,&
'@                                                            ',/,&
'@    Verifier que le fichier suite auxiliaire utilise        ',/,&
'@      correspond bien au cas traite                         ',/,&
'@                                                            ',/,&
'@    Le calcul se poursuit...                                ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8300 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : LECTURE DU FICHIER SUITE AUXILIAIRE         ',/,&
'@    =========                                               ',/,&
'@      ', A54                                                 ,/,&
'@                                                            ',/,&
'@    Certaines grandeurs n''ont pas pu etre lues dans le     ',/,&
'@      fichier suite auxiliaire.                             ',/,&
'@    Elles seront initialisees par des valeurs par defaut.   ',/,&
'@                                                            ',/,&
'@    Cette situation peut se produire lorsque le fichier     ',/,&
'@      suite est issu d''un calcul realise avec des options  ',/,&
'@      differentes ou lorsqu''il a ete endommage.            ',/,&
'@                                                            ',/,&
'@    Le calcul se poursuit...                                ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8400 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : LECTURE DU FICHIER SUITE AUXILIAIRE         ',/,&
'@    =========                                               ',/,&
'@      REPRISE  DE CALCUL           AVEC IDTVAR = ',I10       ,/,&
'@      A PARTIR D''UN CALCUL REALISE AVEC IDTVAR = ',I10      ,/,&
'@                                                            ',/,&
'@    Le mode de marche en temps a ete modifie.               ',/,&
'@    La valeur (uniforme) du pas de temps est                ',/,&
'@      DTREF = ',E12.4   ,' fournie.                         ',/,&
'@                                                            ',/,&
'@    Il est conseille cependant de                           ',/,&
'@      verifier la valeur de IDTVAR.                         ',/,&
'@                                                            ',/,&
'@    Verifier que le fichier suite auxiliaire utilise        ',/,&
'@      correspond bien au cas traite                         ',/,&
'@                                                            ',/,&
'@    Le calcul se poursuit...                                ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8611 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : LECTURE DU FICHIER SUITE AUXILIAIRE         ',/,&
'@    =========                                               ',/,&
'@                                                            ',/,&
'@      Modele de combustion charbon pulverise.               ',/,&
'@      On ne trouve pas la masse volumique des charbons dans ',/,&
'@        le fichier suite. C''est naturel si le calcul       ',/,&
'@        precedent n''etait pas un calcul charbon pulverise. ',/,&
'@        La valeur par defaut est utilisee comme valeur      ',/,&
'@        initiale :                                          ',/,&
'@         Charbon        rho                                 '  )
 8612 format(                                                     &
'@        ',I10   ,'  ',E14.5                                    )
 8613 format(                                                     &
'@                                                            ',/,&
'@    Le calcul peut etre execute.                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8900 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ERREUR A LA FERMETURE DU FICHIER SUITE      ',/,&
'@    =========                                     AUXILIAIRE',/,&
'@                                                            ',/,&
'@    Probleme sur le fichier de nom (',A13,')                ',/,&
'@                                                            ',/,&
'@    Le calcul se poursuit...                                ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

#else

 8001 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING:       WHEN READING THE AUXILIARY RESTART FILE  ',/,&
'@    =======                                                 ',/,&
'@                                                            ',/,&
'@      The maximum number of scalars NSCAMX supported by     ',/,&
'@        the writing format of the restart file is           ',/,&
'@        NFMTSC = ',I10                                       ,/,&
'@      There is here a greater number of scalars            ',/, &
'@        NSCAMX = ',I10                                       ,/,&
'@       It is possible not to read the scalars which have    ',/,&
'@        a greater number.                                   ',/,&
'@                                                            ',/,&
'@    The run will continue.                                  ',/,&
'@                                                            ',/,&
'@    Check the subroutine lecamx.                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8002 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING:       WHEN READING THE AUXILIARY RESTART FILE  ',/,&
'@    =======                                                 ',/,&
'@                                                            ',/,&
'@      The number of the max mass flux NVARMX supported by   ',/,&
'@        the writing format of the suite file is             ',/,&
'@        NFMTFL = ',I10                                       ,/,&
'@      There is here a greater number of flux max            ',/,&
'@        NVARMX = ',I10                                       ,/,&
'@       It is not possible to read the fluxes which have     ',/,&
'@        a greater number.                                   ',/,&
'@                                                            ',/,&
'@    The run will continue.                                  ',/,&
'@                                                            ',/,&
'@    Check the subroutine lecamx.                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8003 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING:       WHEN READING THE AUXILIARY RESTART FILE  ',/,&
'@    =======                                                 ',/,&
'@                                                            ',/,&
'@      The max number of moments NBMOMX supported by         ',/,&
'@        the writing format of the suite file is             ',/,&
'@        NFMTMO = ',I10                                       ,/,&
'@      There is here a greater number of moments             ',/,&
'@        NVARMX = ',I10                                       ,/,&
'@       It is not possible to read the moments which have    ',/,&
'@        a greater number.                                   ',/,&
'@                                                            ',/,&
'@    The run will continue.                                  ',/,&
'@                                                            ',/,&
'@    Check the subroutine lecamx.                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8004 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING:       WHEN READING THE AUXILIARY RESTART FILE  ',/,&
'@    =======                                                 ',/,&
'@                                                            ',/,&
'@      The max number of structures NBMOMX supported by      ',/,&
'@        the writing format of the suite file is             ',/,&
'@        NFMTST = ',I10                                       ,/,&
'@      There is here a greater number of structures          ',/,&
'@        NSTRMX = ',I10                                       ,/,&
'@       If the effective number of structures is greater,    ',/,&
'@        these will not be reread.                           ',/,&
'@                                                            ',/,&
'@    The run will continue.                                  ',/,&
'@                                                            ',/,&
'@    Check the subroutine lecamx.                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8200 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: WHEN READING THE AUXILIARY RESTART FILE        ',/,&
'@    =======                                                 ',/,&
'@      PREVIOUS and PRESENT INPUT DATA ARE DIFFERENT         ',/,&
'@                                                            ',/,&
'@    The number of the faces ',A8  ,' has been modified      ',/,&
'@                                                            ',/,&
'@    The run can continue but the data on the                ',/,&
'@      faces ',A8  ,' will not be reread                     ',/,&
'@      in the suite file.                                    ',/,&
'@    They will be initialised by the default values.         ',/,&
'@                                                            ',/,&
'@     This situation can occur when the restart file         ',/,&
'@      originates from a run using different options         ',/,&
'@      to join the grids or when the periodicity boundary    ',/,&
'@      conditions have been modified.                        ',/,&
'@     This situation can also be generated when the          ',/,&
'@      run is conducted on a different machine               ',/,&
'@      in which case the precision of the machine modifies   ',/,&
'@      the number of faces generated when joinning the grids.',/,&
'@                                                            ',/,&
'@     Finally, this situation can be due to the fact that    ',/,&
'@      the auxiliary restart file does not correspond to     ',/,&
'@      the present case.                                     ',/,&
'@                                                            ',/,&
'@    Verify that the auxiliary restart file being used       ',/,&
'@      corresponds to the present case.                      ',/,&
'@                                                            ',/,&
'@     The run will continue...                               ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8205 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: WHEN READING THE AUXILIARY RESTART FILE        ',/,&
'@    =========                                               ',/,&
'@      CHECKPOINT DATA ARE MULTIPHASE                        ',/,&
'@                                                            ',/,&
'@  Number of phases (checkpoint) : ',I10                      ,/,&
'@                                                            ',/,&
'@    The computation cannot be executed.                     ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8220 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: WHEN READING THE AUXILIARY RESTART FILE        ',/,&
'@    =======                                                 ',/,&
'@      THE RUN RESTARTED            WITH ITURB = ',I4         ,/,&
'@      FROM RUN CONDUCTED WITH           ITURB = ',I4         ,/,&
'@                                                            ',/,&
'@    The Turbulence model has been modified.                 ',/,&
'@                                                            ',/,&
'@    It is advised however in this case to                   ',/,&
'@      verify the value of ITURB                             ',/,&
'@                                                            ',/,&
'@    Verify that the auxiliary restart file being used       ',/,&
'@      corresponds  to the present case.                     ',/,&
'@                                                            ',/,&
'@    The run will continue...                                ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8300 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: WHEN READING THE AUXILIARY RESTART FILE        ',/,&
'@    =======                                                 ',/,&
'@      ', A54                                                 ,/,&
'@                                                            ',/,&
'@    It was not possible to read some values from the        ',/,&
'@      auxiliary restart file.                               ',/,&
'@    They will be initialised by the default values.         ',/,&
'@                                                            ',/,&
'@     This situation can occur when the restart file         ',/,&
'@      originates from a run realised with different         ',/,&
'@      options or when the file is damaged.                  ',/,&
'@                                                            ',/,&
'@    The run will continue...                                ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8400 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: WHEN READING THE AUXILIARY RESTART FILE        ',/,&
'@    =======                                                 ',/,&
'@      THE RUN RESTARTED            WITH IDTVAR = ',I10       ,/,&
'@      FROM RUN CONDUCTED WITH            IDTVAR = ',I10      ,/,&
'@                                                            ',/,&
'@    The variable time step method has been modified.        ',/,&
'@    The (uniform) value of the time step is                 ',/,&
'@      DTREF = ',E12.4                                        ,/,&
'@                                                            ',/,&
'@    It is advised however in this case to                   ',/,&
'@      verify the value of IDTVAR.                           ',/,&
'@                                                            ',/,&
'@    Verify that the auxiliary restart file being used       ',/,&
'@      corresponds  to the present case.                     ',/,&
'@                                                            ',/,&
'@    The run will continue...                                ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8611 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: WHEN READING THE AUXILIARY RESTART FILE        ',/,&
'@    =======                                                 ',/,&
'@                                                            ',/,&
'@      Combustion model for pulverised coal                  ',/,&
'@      The densities of the coals can not be found           ',/,&
'@        in the restart file. This is normal if the          ',/,&
'@        previous run did not include pulverised coal.       ',/,&
'@        The default value is used as an initial             ',/,&
'@        value :                                             ',/,&
'@         Coal           rho                                 '  )
 8612 format(                                                     &
'@        ',I10   ,'  ',E14.5                                    )
 8613 format(                                                     &
'@                                                            ',/,&
'@    The run will continue...                                ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8900 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: ERROR CLOSING THE AUXILIARY RESTART FILE       ',/,&
'@    =======                                                 ',/,&
'@                                                            ',/,&
'@    Problem in the file named (',A13,')                     ',/,&
'@                                                            ',/,&
'@    The run will continue...                                ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

#endif


! --- ERREURS

#if defined(_CS_LANG_FR)

 9000 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,&
'@    =========                                     AUXILIAIRE',/,&
'@      ERREUR A L''OUVERTURE DU FICHIER SUITE                ',/,&
'@                                                            ',/,&
'@    Le calcul ne peut pas etre execute.                     ',/,&
'@                                                            ',/,&
'@    Verifier l''existence et le nom (',A13,') du            ',/,&
'@        fichier suite dans le repertoire de travail.        ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9100 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,&
'@    =========                                     AUXILIAIRE',/,&
'@      TYPE DE FICHIER INCORRECT                             ',/,&
'@                                                            ',/,&
'@    Le fichier ',A13      ,' ne semble pas etre un fichier  ',/,&
'@      suite auxiliaire.                                     ',/,&
'@                                                            ',/,&
'@    Le calcul ne peut etre execute.                         ',/,&
'@                                                            ',/,&
'@    Verifier que le fichier suite utilise correspond bien   ',/,&
'@        a un fichier suite auxiliaire.                      ',/,&
'@    Si necessaire, il est possible de desactiver la lecture ',/,&
'@      du fichier suite auxiliaire par ILEAUX = 0.           ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9101 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,&
'@    =========                                     AUXILIAIRE',/,&
'@      DONNEES AMONT ET ACTUELLES INCOHERENTES               ',/,&
'@                                                            ',/,&
'@    Le nombre de cellules a ete modifie                     ',/,&
'@                                                            ',/,&
'@    Le calcul ne peut etre execute.                         ',/,&
'@                                                            ',/,&
'@    Verifier que le fichier suite utilise correspond bien   ',/,&
'@        au cas traite.                                      ',/,&
'@    Si necessaire, il est possible de desactiver la lecture ',/,&
'@      du fichier suite auxiliaire par ILEAUX = 0.           ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9200 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION: ARRET A LA LECTURE DU FICHIER SUITE          ',/,&
'@    =========                                     AUXILIAIRE',/,&
'@                                                            ',/,&
'@      ', A54                                                 ,/,&
'@                                                            ',/,&
'@    Le calcul ne peut pas etre execute.                     ',/,&
'@                                                            ',/,&
'@    Verifier que le fichier suite utilise n''a pas ete      ',/,&
'@      endommage.                                            ',/,&
'@    Si necessaire, il est possible de desactiver la lecture ',/,&
'@      du fichier suite auxiliaire par ILEAUX = 0.           ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9210 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ERREUR A LA LECTURE DU FICHIER SUITE        ',/,&
'@    =========                                     AUXILIAIRE',/,&
'@                                                            ',/,&
'@      ERREUR A LA LECTURE DE L''INDICATEUR DE METHODE ALE   ',/,&
'@                                                            ',/,&
'@    Il se peut que le fichier suite relu corresponde a une  ',/,&
'@      version anterieure de Code_Saturne, sans methode ALE. ',/,&
'@    Le calcul sera execute en reinitialisant toutes les     ',/,&
'@      donnees ALE.                                          ',/,&
'@    Verifier neanmoins que le fichier suite utilise n''a    ',/,&
'@        pas ete endommage.                                  ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9300 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,&
'@    =========                                     AUXILIAIRE',/,&
'@                                                            ',/,&
'@    La relecture de l''ancienne moyenne ',I10  ,' qui doit  ',/,&
'@      permettre d''initialiser la nouvelle moyenne ',I10     ,/,&
'@      a echoue                                              ',/ &
'@                                                            ',/,&
'@    Le calcul ne sera pas  execute.                         ',/,&
'@                                                            ',/,&
'@    Verifier le fichier suite auxiliaire ou                 ',/,&
'@      specifier dans usipsu que la moyenne doit etre        ',/,&
'@      reinitialisee, en indiquant : IMOOLD(',I10   ,') = -1 ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9310 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,&
'@    =========                                     AUXILIAIRE',/,&
'@                                                            ',/,&
'@    La relecture du cumul de la duree associee a            ',/,&
'@                    l''ancienne moyenne ',I10  ,' qui doit  ',/,&
'@      permettre d''initialiser la nouvelle moyenne ',I10     ,/,&
'@      a echoue                                              ',/ &
'@                                                            ',/,&
'@    Le calcul ne sera pas  execute.                         ',/,&
'@                                                            ',/,&
'@    Verifier le fichier suite auxiliaire ou                 ',/,&
'@      specifier dans usipsu que la moyenne doit etre        ',/,&
'@      reinitialisee, en indiquant : IMOOLD(',I10   ,') = -1 ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9320 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,&
'@    =========                                     AUXILIAIRE',/,&
'@                                                            ',/,&
'@      ERREUR LORS DE LA LECTURE DU DEPLACEMENT AUX NOEUDS   ',/,&
'@        DU MAILLAGE (METHODE ALE)                           ',/,&
'@                                                            ',/,&
'@    Le calcul ne peut pas etre execute.                     ',/,&
'@                                                            ',/,&
'@    Verifier que le fichier suite utilise n''a pas ete      ',/,&
'@        endommage.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9321 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,&
'@    =========                                     AUXILIAIRE',/,&
'@                                                            ',/,&
'@      ERREUR LORS DE LA LECTURE DES DONNEES DES STRUCTURES  ',/,&
'@        MOBILES (METHODE ALE)                               ',/,&
'@                                                            ',/,&
'@    Le calcul ne peut pas etre execute.                     ',/,&
'@                                                            ',/,&
'@    Verifier que le fichier suite utilise n''a pas ete      ',/,&
'@        endommage.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9400 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,&
'@    =========                                     AUXILIAIRE',/,&
'@                                                            ',/,&
'@    La relecture des variable CO3DP                         ',/,&
'@                                                            ',/,&
'@    Le calcul ne sera pas  execute.                         ',/,&
'@                                                            ',/,&
'@    Verifier le fichier suite auxiliaire                    ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9500 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,&
'@    =========                                     AUXILIAIRE',/,&
'@                                                            ',/,&
'@    La relecture des variable EBU                           ',/,&
'@                                                            ',/,&
'@    Le calcul ne sera pas  execute.                         ',/,&
'@                                                            ',/,&
'@    Verifier le fichier suite auxiliaire                    ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9600 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,&
'@    =========                                     AUXILIAIRE',/,&
'@                                                            ',/,&
'@    La relecture des variable LWC                           ',/,&
'@                                                            ',/,&
'@    Le calcul ne sera pas  execute.                         ',/,&
'@                                                            ',/,&
'@    Verifier le fichier suite auxiliaire                    ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

#else

 9000 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: STOP WHILE READING THE AUXILIARY RESTART FILE  ',/,&
'@    =======                                                 ',/,&
'@      ERROR WHEN OPENING THE RESTART FILE                   ',/,&
'@                                                            ',/,&
'@    The run can not be executed.                            ',/,&
'@                                                            ',/,&
'@    Verify the existence and the name (',A13,') of the      ',/,&
'@        restart file in the work directory.                 ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9100 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: STOP WHILE READING THE AUXILIARY RESTART FILE  ',/,&
'@    =======                                                 ',/,&
'@      FILE TYPE IS INCORRECT                                ',/,&
'@                                                            ',/,&
'@    The file ',A13      ,' does not appear to be an         ',/,&
'@      auxiliary file.                                       ',/,&
'@                                                            ',/,&
'@    The run can not be executed.                            ',/,&
'@                                                            ',/,&
'@    Verify that the restart file used corresponds to        ',/,&
'@        an auxiliary restart file.                          ',/,&
'@    If necessary, it is possible to deactivate the reading  ',/,&
'@      of the auxiliary restart file by ILEAUX = 0.          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9101 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: STOP WHILE READING THE AUXILIARY RESTART FILE  ',/,&
'@    =======                                                 ',/,&
'@      INCOHERENT PREVIOUS NAD ACTUAL DATA                   ',/,&
'@                                                            ',/,&
'@    The number of cells was modified                        ',/,&
'@                                                            ',/,&
'@    The run can not be executed.                            ',/,&
'@                                                            ',/,&
'@    Verify that the restart file used corresponds to        ',/,&
'@        the present case.                                   ',/,&
'@    If necessary, it is possible to deactivate the reading  ',/,&
'@      of the auxiliary restart file by ILEAUX = 0.          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9200 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: STOP WHILE READING THE AUXILIARY RESTART FILE  ',/,&
'@    =======                                                 ',/,&
'@                                                            ',/,&
'@      ', A54                                                 ,/,&
'@                                                            ',/,&
'@    The run can not be executed.                            ',/,&
'@                                                            ',/,&
'@    Verify that the restart file used has not been damaged. ',/,&
'@                                                            ',/,&
'@    If necessary, it is possible to deactivate the reading  ',/,&
'@      of the auxiliary restart file by ILEAUX = 0.          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9210 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: ERROR WHILE READING THE AUXILIARY RESTART FILE ',/,&
'@    =======                                                 ',/,&
'@                                                            ',/,&
'@      ERROR WHEN READING THE INDICATOR OF THE ALE METHOD    ',/,&
'@                                                            ',/,&
'@    It is possible that the file read corresponds to an old ',/,&
'@      version of Code_Saturne, without the ALE method.      ',/,&
'@    The run will be executed with reinitialising all        ',/,&
'@      ALE data.                                             ',/,&
'@                                                            ',/,&
'@    Verify that the restart file used has not been damaged. ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9300 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: STOP WHILE READING THE AUXILIARY RESTART FILE  ',/,&
'@    =======                                                 ',/,&
'@                                                            ',/,&
'@    The reading of the previous average ',I10  ,' that will ',/,&
'@      allow the initialisation of the new average ',I10      ,/,&
'@      has failed                                            ',/ &
'@                                                            ',/,&
'@    The run will not be executed.                           ',/,&
'@                                                            ',/,&
'@    Verify the auxiliary restart file or specify            ',/,&
'@      in usipsu that the average has to be reinitialised,   ',/,&
'@      by indicating: IMOOLD(',I10   ,') = -1                ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9310 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: STOP WHILE READING THE AUXILIARY RESTART FILE  ',/,&
'@    =======                                                 ',/,&
'@                                                            ',/,&
'@    The reading of the accumulated period associated with   ',/,&
'@                   the previous average ',I10  ,' that will ',/,&
'@      allow to initialise the new average ',I10              ,/,&
'@      has failed                                            ',/ &
'@                                                            ',/,&
'@    The run will not be executed.                           ',/,&
'@                                                            ',/,&
'@    Verify the auxiliary restart file or specify            ',/,&
'@      in usipsu that the average has to be reinitialised,   ',/,&
'@      by indicating: IMOOLD(',I10   ,') = -1                ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9320 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: STOP WHILE READING THE AUXILIARY RESTART FILE  ',/,&
'@    =======                                                 ',/,&
'@                                                            ',/,&
'@      ERROR WHILE READING MESH VERTICES MOVEMENT DATA       ',/,&
'@        (ALE METHOD)                                        ',/,&
'@                                                            ',/,&
'@    The run can not be executed.                            ',/,&
'@                                                            ',/,&
'@    Verify that the restart file used has not been damaged  ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9321 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: STOP WHILE READING THE AUXILIARY RESTART FILE  ',/,&
'@    =======                                                 ',/,&
'@                                                            ',/,&
'@      ERROR WHILE READING MOVING STRUCTURES DATA            ',/,&
'@        (ALE METHOD)                                        ',/,&
'@                                                            ',/,&
'@    The run can not be executed.                            ',/,&
'@                                                            ',/,&
'@    Verify that the restart file used has not been damaged  ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9400 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: STOP WHILE READING THE AUXILIARY RESTART FILE  ',/,&
'@    =======                                                 ',/,&
'@                                                            ',/,&
'@    Error reading the CO3DP variables                       ',/,&
'@                                                            ',/,&
'@    The run can not be executed.                            ',/,&
'@                                                            ',/,&
'@    Verify the auxiliary restart file                       ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9500 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: STOP WHILE READING THE AUXILIARY RESTART FILE  ',/,&
'@    =======                                                 ',/,&
'@                                                            ',/,&
'@    Error reading the EBU variables                         ',/,&
'@                                                            ',/,&
'@    The run will not be executed.                           ',/,&
'@                                                            ',/,&
'@    Verify the auxiliary restart file                       ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9600 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: STOP WHILE READING THE AUXILIARY RESTART FILE  ',/,&
'@    =======                                                 ',/,&
'@                                                            ',/,&
'@    Error reading the LWC variables                         ',/,&
'@                                                            ',/,&
'@    The run will not be executed.                           ',/,&
'@                                                            ',/,&
'@    Verify the auxiliary restart file                       ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

#endif

end subroutine
