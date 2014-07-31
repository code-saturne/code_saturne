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

subroutine lecamp &
!================

 ( ncelet , ncel   ,                                              &
   nvar   , nscal  ,                                              &
   rtp    )

!===============================================================================

! FONCTION :
! ----------
! LECTURE DU FICHIER SUITE PRINCIPAL

! ON S'ARRETE SI ERREUR DE LECTURE OU
!             OU NCEL DIFFERENT
!             OU NTPABS.GE.NTMABS


!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! ncelet           ! i  ! <-- ! number of extended (real + ghost) cells        !
! ncel             ! i  ! <-- ! number of cells                                !
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! jphas            ! e  ! <-- ! nombre de phases du calcul precedent           !
! ljtu             ! e  ! <-- ! longueur de jturb                              !
! jturb            ! te ! <-- ! modeles de turb calcul precedent               !
! jturbt           ! te ! <-- ! modeles de flux turb calcul precedent          !
! rtp              ! tr ! --> ! variables de calcul au centre des              !
! (ncelet,*)       !    !     !    cellules (instant courant        )          !
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
use cstphy
use cstnum
use entsor
use optcal
use pointe
use numvar
use albase
use parall
use cplsat
use field
use atincl, only: init_at_chem
use atchem, only: ichemistry
use siream, only: iaerosol
use turbomachinery

!===============================================================================

implicit none

! Arguments

integer          ncelet , ncel
integer          nvar   , nscal


double precision rtp(ncelet,*)

! Local variables

character        rubriq*64,rubrik*64,car4*4
character        cindfp*2
character        cphase*2
character        ficsui*32
character*80     fname

integer          iel
integer          ivar  , iscal , ii    ,  ivers
integer          jphas , jvar  , jscal , jscaus, jscapp
integer          iok   , ncelok, nfaiok, nfabok, nsomok
integer          ierror, irtyp,  itysup, nbval
integer          nberro, ilecec
integer          iturph, jturph, itytph, jtytph
integer          nfmtsc, nfmtru
integer          jturb, jtytur, jale, jturbt
integer          impamo
integer          ival
integer          f_id
double precision d2s3, d2s3xk
double precision rval

double precision, dimension(:,:), pointer :: xut

!===============================================================================

!===============================================================================
! 0. INITIALISATIONS
!===============================================================================

! Initialize variables to avoid compiler warnings

ivar = 0

! Memoire


!  ---> Banniere
write(nfecra,1000)

!     Longueur pour format de print
nfmtru = 36

!  --->  On code en chaine le numero des phases et scalaires

!     Nombre de scalaires max pour les formats choisis
nfmtsc = 9999

!     Indefini a 2 caracteres
cindfp='YY'

!     Codage en chaine de caracteres du numero de la phase
write(cphase,'(i2.2)') 1

!     Avertissement
if(nscamx.gt.nfmtsc) then
  write(nfecra,8001)nfmtsc,nscamx
endif

!===============================================================================
! 1. OUVERTURE DU FICHIER OU STOP
!===============================================================================

!  ---> Ouverture du fichier (ILECEC = 1 : lecture)
ilecec = 1
ficsui = 'main'
call opnsui(ficsui,len(ficsui),ilecec,impamo,ierror)
!==========
if (ierror.ne.0) then
  write(nfecra,9100) ficsui
  call csexit (1)
endif

! ---> Debut de la lecture
write(nfecra,1100)

!===============================================================================
! 2. ENTETES DU FICHIER SUITE OU STOP
!===============================================================================

!  --->  Rubrique "fichier suite ppal"
!        Pourrait porter le numero de version si besoin.
!        On ne se sert pas de IVERS pour le moment.

itysup = 0
nbval  = 1
irtyp  = 1
rubriq = 'version_fichier_suite_principal'
call lecsui(impamo,rubriq,len(rubriq),itysup,nbval,irtyp,ivers,   &
            ierror)

if (ierror.ne.0) then
  write(nfecra,9200)ficsui
  call csexit (1)
endif

!  --->  Tests sur les supports

call tstsui(impamo,ncelok,nfaiok,nfabok,nsomok)
!==========

if (ncelok.eq.0) then
  write(nfecra,9201)
  call csexit (1)
endif

!     Inutile de tester les supports "faces internes" et "faces de bord"
!     ils ne sont pas utilises ici


!  --->  Tests sur les nombres de variables

nberro = 0

itysup = 0
nbval  = 1
irtyp  = 1

rubriq = 'nombre_variables'
call lecsui(impamo,rubriq,len(rubriq),itysup,nbval,irtyp,jvar,    &
            ierror)
nberro=nberro+ierror

rubriq = 'nombre_scalaires'
call lecsui(impamo,rubriq,len(rubriq),itysup,nbval,irtyp,jscal,   &
            ierror)
nberro=nberro+ierror

rubriq = 'nombre_scalaires_us'
call lecsui(impamo,rubriq,len(rubriq),itysup,nbval,irtyp,jscaus,  &
            ierror)
nberro=nberro+ierror

rubriq = 'nombre_scalaires_pp'
call lecsui(impamo,rubriq,len(rubriq),itysup,nbval,irtyp,jscapp,  &
            ierror)
nberro=nberro+ierror

rubriq = 'nombre_phases'
call lecsui(impamo,rubriq,len(rubriq),itysup,nbval,irtyp,jphas,   &
            ierror)
nberro=nberro+ierror


!  ---> On s'arrete si erreur

if (nberro.ne.0) then
  write(nfecra,9210)
  call csexit (1)
endif

!  ---> On ne sait relire que des calculs monophasiques

if (jphas.ne.1) then
  write(nfecra,8205) jphas
  call csexit(1)
endif

!  ---> On previent si des parametres sont differents

if ( jvar  .ne.nvar   .or. jscal .ne.nscal  .or.                  &
     jscaus.ne.nscaus .or. jscapp.ne.nscapp ) then
  write(nfecra,8210)                                              &
         jvar, jscal, jscaus, jscapp,                             &
         nvar, nscal, nscaus, nscapp
endif

! ---> Fin de la lecture des dimensions
write(nfecra,1299)

!===============================================================================
! 3. CORRESPONDANCE DES SCALAIRES OU STOP
!===============================================================================

!  ---> Correspondance avec les scalaires du calcul precedent.
!       Le tableau ISCOLD manipule des numero de scalaire "globaux"
!         i.e. de 1 a NSCAL (et non pas de 1 a NSCAUS ou 1 a NSCAPP)

!     On complete si besoin le tableau ISCOLD. Si ses valeurs sont
!       egales a -999, c'est que l'utilisateur ne les a pas modifiees,
!       et qu'il ne souhaite donc rien de particulier quant a la
!       correspondance des nouveaux et des anciens scalaires.
!       On etablit alors la correspondance
!       avec l'ancien scalaire (qui est peut etre une variance) de meme
!       numero s'il existe ou a rien du tout sinon.

!     Pour bien faire la difference entre scalaires utilisateurs et
!       physiques particulieres, on part du principe suivant :
!         par defaut
!         pour le nouveau scalaire utilisateur j
!           l'ancien scalaire correspondant est
!             l'ancien scalaire utilisateur j s'il existe
!                                             et rien sinon
!         pour le nouveau scalaire physique particuliere j
!           l'ancien scalaire correspondant est pour le moment
!             l'ancien scalaire physique particuliere j s'il existe
!                                                       et rien sinon
!             a terme, on pourra preciser les correspondances plus
!             finement, en particulier si on souhaite faire des suites
!             entre deux modeles differents.

!     Attention, on suppose ici que les scalaires sont numerotes dans
!       l'ordre scalaires utilisateurs (NSCAUS) puis scalaires
!       physiques particulieres (NSCAPP)

if(nscaus.gt.0) then
  do ii = 1, nscaus
    iscal = ii
    if(iscold(iscal).eq.-999) then
      if(ii.le.jscaus) then
        iscold(iscal) = ii
      else
        iscold(iscal) = 0
      endif
    endif
  enddo
endif

if(nscapp.gt.0) then
  do ii = 1, nscapp
    iscal = iscapp(ii)
    if(iscold(iscal).eq.-999) then
      if(ii.le.jscapp) then
        iscold(iscal) = ii+jscaus
      else
        iscold(iscal) = 0
      endif
    endif
  enddo
endif

!     On verifie les valeurs du tableau ISCOLD. Elles doivent toutes
!       etre inferieures a JSCAL et positives ou nulles.
if(nscal.gt.0) then
  iok = 0
  do iscal = 1, nscal
    if(iscold(iscal).gt.jscal) then
      write(nfecra,9320)iscal, iscold(iscal), jscal
      iok = iok+1
    endif
    if(iscold(iscal).lt.0) then
      write(nfecra,9321)iscal, iscold(iscal), jscal
      iok = iok+1
    endif
  enddo
  if(iok.ne.0) then
    write(nfecra,9330) iok
    call csexit (1)
  endif
endif

!===============================================================================
! 4. OPTIONS OU STOP
!===============================================================================

!  --->  Lecture des options

nberro = 0

!     Nombre de pas de temps, instant precedent

rubriq = 'nbre_pas_de_temps'
itysup = 0
nbval  = 1
irtyp  = 1
call lecsui(impamo,rubriq,len(rubriq),itysup,nbval,irtyp,ival,ierror)
ntpabs = ival ! no direct read to avoid pointer issue
nberro=nberro+ierror

rubriq = 'instant_precedent'
itysup = 0
nbval  = 1
irtyp  = 2
call lecsui(impamo,rubriq,len(rubriq),itysup,nbval,irtyp,rval,ierror)
ttpabs = rval ! no direct read to avoid pointer issue
nberro=nberro+ierror


!     Modeles de turbulence

rubriq = 'modele_turbulence_phase'//cphase
itysup = 0
nbval  = 1
irtyp  = 1
call lecsui(impamo,rubriq,len(rubriq),itysup,nbval,irtyp,       &
            jturb,ierror)
nberro=nberro+ierror
jtytur=jturb/10

! --->  Stop si erreur
if (nberro.ne.0) then
  write(nfecra,9400)
  call csexit (1)
endif

!     Methode ALE

nberro = 0

rubriq = 'methode_ALE'
itysup = 0
nbval  = 1
irtyp  = 1
call lecsui(impamo,rubriq,len(rubriq),itysup,nbval,irtyp,jale,    &
            ierror)
nberro=nberro+ierror

! --->  Message si erreur (pas de stop pour compatibilite avec les fichiers anterieurs)
!       -> on n'affiche le message que si IALE=1 (sinon RAS)
if (nberro.ne.0) then
  if (iale.eq.1) write(nfecra,9401)
  jale = 0
endif

! --->  Stop si pas de temps incoherent
if(ntpabs.ge.ntmabs.and.inpdt0.eq.0) then
  write(nfecra,9410)ntpabs,ntmabs
  call csexit (1)
endif

! --->  Informations
write(nfecra,2410) ntpabs
write(nfecra,2411) ttpabs

! --->  Donnees modifiees
if (iturb .ne. jturb)                            &
    write(nfecra,8410) iturb, jturb

! --->  Si le calcul precedent etait en ALE, on DOIT relire les
!         coordonnees des noeuds dans le fichier auxiliaire
if (iale.eq.1 .and. jale.eq.1) then
  if (ileaux.ne.1) then
    write(nfecra,9402)jale,iale,ileaux
    call csexit(1)
  endif
endif

!     Instant de maillage mobile precedent (rotor/stator)

nberro = 0

rubriq = 'instant_mobile_precedent'
itysup = 0
nbval  = 1
irtyp  = 2
call lecsui(impamo,rubriq,len(rubriq),itysup,nbval,irtyp,rval,ierror)
ttpmob = rval ! no direct read to avoid pointer issue
nberro=nberro+ierror

! --->  Message si erreur (pas de stop pour compatibilite avec les fichiers anterieurs)
!       -> on n'affiche le message que si imobil=1 ou iturbo=2 (sinon RAS)
if (nberro.ne.0) then
  if (imobil.eq.1 .or. iturbo.eq.2) write(nfecra,9403) ttpabs
  ttpmob = ttpabs
endif

! --->  Information (uniquement si imobil=1 ou iturbo=2
!         et pas d affichage precedent)
if (imobil.eq.1 .or. iturbo.eq.2) then
  if (nberro.eq.0)  write(nfecra,2412) ttpmob
endif

! --->  Fin de la lecture des options
write(nfecra,1499)

!================================================================
! 5. LECTURE DES VARIABLES
!================================================================

nberro = 0

! --->  Pression
!     (a priori une seule (non fonction du nombre de phases))

rubriq = 'pression_ce_phase'//cphase
itysup = 1
nbval  = 1
irtyp  = 2
call lecsui(impamo,rubriq,len(rubriq),itysup,nbval,irtyp,         &
            rtp(1,ipr),ierror)
nberro=nberro+ierror


! --->  Vitesse et turbulence
!     On boucle sur les phases communes afin de recuperer
!       les infos existantes dans le fichier suite.
!       Pour les nouvelles phases : deja fait par defaut dans INIVA0

iturph = iturb
jturph = jturb
itytph = itytur
jtytph = jtytur

!     Vitesse

itysup = 1
nbval  = 1
irtyp  = 2

rubriq = 'vitesse_u_ce_phase'//cphase
call lecsui(impamo,rubriq,len(rubriq),itysup,nbval,irtyp,      &
     rtp(1,iu),ierror)
nberro=nberro+ierror

rubriq = 'vitesse_v_ce_phase'//cphase
call lecsui(impamo,rubriq,len(rubriq),itysup,nbval,irtyp,      &
     rtp(1,iv),ierror)
nberro=nberro+ierror

rubriq = 'vitesse_w_ce_phase'//cphase
call lecsui(impamo,rubriq,len(rubriq),itysup,nbval,irtyp,      &
     rtp(1,iw),ierror)
nberro=nberro+ierror



!     Turbulence (ke -> ke, ke -> rij, rij -> ke, rij -> rij ...)
!       Quand on ne sait pas deduire les grandeurs turbulence des infos
!         disponibles dans le fichier suite, on ne fait rien :
!         les variables gardent alors les initialisations par defaut
!         faites dans INIVA0


!   -- Le nouveau calcul est en k-epsilon


if (itytph.eq.2) then

  !     * k-e ou v2f (phi-fbar ou BL-v2/k) -> k-e

  if(jtytph.eq.2 .or. jtytph.eq.5) then

    itysup = 1
    nbval  = 1
    irtyp  = 2

    rubriq = 'k_ce_phase'//cphase
    call lecsui(impamo,rubriq,len(rubriq),itysup,nbval,irtyp,   &
         rtp(1,ik),ierror)
    nberro=nberro+ierror

    rubriq = 'eps_ce_phase'//cphase
    call lecsui(impamo,rubriq,len(rubriq),itysup,nbval,irtyp,   &
         rtp(1,iep),ierror)
    nberro=nberro+ierror

    !     * rij -> k-e

  elseif(jtytph.eq.3) then

    itysup = 1
    nbval  = 1
    irtyp  = 2

    rubriq = 'R11_ce_phase'//cphase
    call lecsui(impamo,rubriq,len(rubriq),itysup,nbval,irtyp,  &
         rtp(1,ik),ierror)
    nberro=nberro+ierror

    !            La variable epsilon sert de tableau de travail
    rubriq = 'R22_ce_phase'//cphase
    call lecsui(impamo,rubriq,len(rubriq),itysup,nbval,irtyp,  &
         rtp(1,iep),ierror)
    nberro=nberro+ierror

    do iel = 1, ncel
      rtp(iel,ik) = rtp(iel,ik) + rtp(iel,iep)
    enddo

    rubriq = 'R33_ce_phase'//cphase
    call lecsui(impamo,rubriq,len(rubriq),itysup,nbval,irtyp,  &
         rtp(1,iep),ierror)
    nberro=nberro+ierror

    do iel = 1, ncel
      rtp(iel,ik) = 0.5d0*(rtp(iel,ik)+rtp(iel,iep))
    enddo

    rubriq = 'eps_ce_phase'//cphase
    call lecsui(impamo,rubriq,len(rubriq),itysup,nbval,irtyp,  &
         rtp(1,iep),ierror)
    nberro=nberro+ierror

    !     * k-omega -> k-e

  else if(jturph.eq.60) then

    itysup = 1
    nbval  = 1
    irtyp  = 2

    rubriq = 'k_ce_phase'//cphase
    call lecsui(impamo,rubriq,len(rubriq),itysup,nbval,irtyp,   &
         rtp(1,ik),ierror)
    nberro=nberro+ierror

    rubriq = 'omega_ce_phase'//cphase
    call lecsui(impamo,rubriq,len(rubriq),itysup,nbval,irtyp,   &
         rtp(1,iep),ierror)
    nberro=nberro+ierror
    !           On transforme ensuite omega en epsilon
    do iel = 1, ncel
      rtp(iel,iep) = cmu*rtp(iel,iep)*rtp(iel,ik)
    enddo

  endif
  !         Rq : laminaire -> k-e  (On ne fait rien, deja fait dans iniva0)


  !   -- Le nouveau calcul est en Rij-epsilon

elseif(itytph.eq.3) then

  !     * k-e ou v2f (phi-fbar ou BL-v2/k) -> rij

  if (jtytph.eq.2 .or. jturph.eq.50) then

    itysup = 1
    nbval  = 1
    irtyp  = 2

    rubriq = 'k_ce_phase'//cphase
    call lecsui(impamo,rubriq,len(rubriq),itysup,nbval,irtyp,  &
         rtp(1,ir11),ierror)
    nberro=nberro+ierror

    d2s3 = 2.d0/3.d0
    do iel = 1, ncel
      d2s3xk = rtp(iel,ir11)*d2s3
      rtp(iel,ir11) = d2s3xk
      rtp(iel,ir22) = d2s3xk
      rtp(iel,ir33) = d2s3xk
      rtp(iel,ir12) = 0.d0
      rtp(iel,ir13) = 0.d0
      rtp(iel,ir23) = 0.d0
    enddo

    rubriq = 'eps_ce_phase'//cphase
    call lecsui(impamo,rubriq,len(rubriq),itysup,nbval,irtyp,  &
         rtp(1,iep),ierror)
    nberro=nberro+ierror

    !     * rij -> rij

  elseif(jtytph.eq.3) then

    do ii   = 1, 7
      if (ii  .eq.1) then
        rubriq = 'R11_ce_phase'//cphase
        ivar = ir11
      elseif (ii  .eq.2) then
        rubriq = 'R22_ce_phase'//cphase
        ivar = ir22
      elseif (ii  .eq.3) then
        rubriq = 'R33_ce_phase'//cphase
        ivar = ir33
      elseif (ii  .eq.4) then
        rubriq = 'R12_ce_phase'//cphase
        ivar = ir12
      elseif (ii  .eq.5) then
        rubriq = 'R13_ce_phase'//cphase
        ivar = ir13
      elseif (ii  .eq.6) then
        rubriq = 'R23_ce_phase'//cphase
        ivar = ir23
      elseif (ii  .eq.7) then
        rubriq = 'eps_ce_phase'//cphase
        ivar = iep
      endif
      itysup = 1
      nbval  = 1
      irtyp  = 2
      call lecsui(impamo,rubriq,len(rubriq),itysup,nbval,irtyp,   &
           rtp(1,ivar),ierror)
      nberro=nberro+ierror
    enddo
    if (jturph.eq.32) then

      rubriq = 'alp_ce_phase'//cphase
      call lecsui(impamo,rubriq,len(rubriq),itysup,nbval,irtyp,   &
                rtp(1,ial),ierror)
      nberro=nberro+ierror
    endif

    !     * k-omega -> rij

  else if (jturph.eq.60) then

    itysup = 1
    nbval  = 1
    irtyp  = 2

    rubriq = 'k_ce_phase'//cphase
    call lecsui(impamo,rubriq,len(rubriq),itysup,nbval,irtyp,   &
         rtp(1,ir11),ierror)
    nberro=nberro+ierror

    rubriq = 'omega_ce_phase'//cphase
    call lecsui(impamo,rubriq,len(rubriq),itysup,nbval,irtyp,   &
         rtp(1,iep),ierror)
    nberro=nberro+ierror
    !     On transforme ensuite omega en epsilon
    do iel = 1, ncel
      rtp(iel,iep) = cmu*rtp(iel,iep)*rtp(iel,ir11)
    enddo

    d2s3 = 2.d0/3.d0
    do iel = 1, ncel
      d2s3xk = rtp(iel,ir11)*d2s3
      rtp(iel,ir11) = d2s3xk
      rtp(iel,ir22) = d2s3xk
      rtp(iel,ir33) = d2s3xk
      rtp(iel,ir12) = 0.d0
      rtp(iel,ir13) = 0.d0
      rtp(iel,ir23) = 0.d0
    enddo

  endif
  !       Rq : laminaire -> rij   (On ne fait rien, deja fait dans iniva0)

  !   -- Le nouveau calcul est en v2f

  !   -- Le nouveau calcul est en v2f (phi-fbar ou BL-v2/k)

elseif(itytph.eq.5) then

  !     * k-e -> v2f (phi-fbar ou BL-v2/k)
  if(jtytph.eq.2) then

    itysup = 1
    nbval  = 1
    irtyp  = 2

    rubriq = 'k_ce_phase'//cphase
    call lecsui(impamo,rubriq,len(rubriq),itysup,nbval,irtyp,   &
         rtp(1,ik),ierror)
    nberro=nberro+ierror

    rubriq = 'eps_ce_phase'//cphase
    call lecsui(impamo,rubriq,len(rubriq),itysup,nbval,irtyp,   &
         rtp(1,iep),ierror)
    nberro=nberro+ierror
    !     On laisse pour phi et fb les initialisations de iniva0

    !     * rij -> v2f (phi-fbar ou BL-v2/k)

  elseif(jtytph.eq.3) then

    itysup = 1
    nbval  = 1
    irtyp  = 2

    rubriq = 'R11_ce_phase'//cphase
    call lecsui(impamo,rubriq,len(rubriq),itysup,nbval,irtyp,  &
         rtp(1,ik),ierror)
    nberro=nberro+ierror

    !            La variable epsilon sert de tableau de travail
    rubriq = 'R22_ce_phase'//cphase
    call lecsui(impamo,rubriq,len(rubriq),itysup,nbval,irtyp,  &
         rtp(1,iep),ierror)
    nberro=nberro+ierror

    do iel = 1, ncel
      rtp(iel,ik) = rtp(iel,ik) + rtp(iel,iep)
    enddo

    rubriq = 'R33_ce_phase'//cphase
    call lecsui(impamo,rubriq,len(rubriq),itysup,nbval,irtyp,  &
         rtp(1,iep),ierror)
    nberro=nberro+ierror

    do iel = 1, ncel
      rtp(iel,ik) = 0.5d0*(rtp(iel,ik)+rtp(iel,iep))
    enddo

    rubriq = 'eps_ce_phase'//cphase
    call lecsui(impamo,rubriq,len(rubriq),itysup,nbval,irtyp,  &
         rtp(1,iep),ierror)
    nberro=nberro+ierror
    !     On laisse pour phi et fb l'initialisations de iniva0
    !     (le v2 qui intervient dans phi n'est pas vraiment une composante de Rij
    !      et qui plus est, on ne saurait pas quelle composante prendre ...)

    !     * v2f (phi-fbar ou BL-v2/k) -> v2f (phi-fbar ou BL-v2/k)

  elseif(jtytph.eq.5) then

    itysup = 1
    nbval  = 1
    irtyp  = 2

    rubriq = 'k_ce_phase'//cphase
    call lecsui(impamo,rubriq,len(rubriq),itysup,nbval,irtyp,   &
         rtp(1,ik),ierror)
    nberro=nberro+ierror

    rubriq = 'eps_ce_phase'//cphase
    call lecsui(impamo,rubriq,len(rubriq),itysup,nbval,irtyp,   &
         rtp(1,iep),ierror)
    nberro=nberro+ierror

    rubriq = 'phi_ce_phase'//cphase
    call lecsui(impamo,rubriq,len(rubriq),itysup,nbval,irtyp,   &
         rtp(1,iphi),ierror)
    nberro=nberro+ierror

    if(iturph.eq.50.and.jturph.eq.50) then
      rubriq = 'fb_ce_phase'//cphase
      call lecsui(impamo,rubriq,len(rubriq),itysup,nbval,irtyp,   &
           rtp(1,ifb),ierror)
      nberro=nberro+ierror
    elseif(iturph.eq.51.and.jturph.eq.51) then
      rubriq = 'al_ce_phase'//cphase
      call lecsui(impamo,rubriq,len(rubriq),itysup,nbval,irtyp,   &
           rtp(1,ial),ierror)
      nberro=nberro+ierror
    endif
    !     Si (phi-fbar -> BL-v2/k) ou (BL-v2/k -> phi-fbar)
    !      on laisse pour al ou fb l'initialisations de iniva0

    !     * k-omega -> v2f

  else if(jturph.eq.60) then

    itysup = 1
    nbval  = 1
    irtyp  = 2

    rubriq = 'k_ce_phase'//cphase
    call lecsui(impamo,rubriq,len(rubriq),itysup,nbval,irtyp,   &
         rtp(1,ik),ierror)
    nberro=nberro+ierror

    rubriq = 'omega_ce_phase'//cphase
    call lecsui(impamo,rubriq,len(rubriq),itysup,nbval,irtyp,   &
         rtp(1,iep),ierror)
    nberro=nberro+ierror
    !           On transforme ensuite omega en epsilon
    do iel = 1, ncel
      rtp(iel,iep) = cmu*rtp(iel,iep)*rtp(iel,ik)
    enddo
    !     On laisse pour phi et fb l'initialisations de iniva0


  endif
  !         Rq : laminaire -> v2f  (On ne fait rien, deja fait dans iniva0)


  !   -- Le nouveau calcul est en k-omega

else if (iturph.eq.60) then

  !     * k-e ou v2f -> k-omega

  if(jtytph.eq.2 .or. jturph.eq.50) then

    itysup = 1
    nbval  = 1
    irtyp  = 2

    rubriq = 'k_ce_phase'//cphase
    call lecsui(impamo,rubriq,len(rubriq),itysup,nbval,irtyp,   &
         rtp(1,ik),ierror)
    nberro=nberro+ierror

    rubriq = 'eps_ce_phase'//cphase
    call lecsui(impamo,rubriq,len(rubriq),itysup,nbval,irtyp,   &
         rtp(1,iomg),ierror)
    nberro=nberro+ierror
    !           On transforme ensuite epsilon en omega
    do iel = 1, ncel
      rtp(iel,iomg) = rtp(iel,iomg)/cmu/rtp(iel,ik)
    enddo

    !     * rij -> k-omega

  elseif(jtytph.eq.3) then

    itysup = 1
    nbval  = 1
    irtyp  = 2

    rubriq = 'R11_ce_phase'//cphase
    call lecsui(impamo,rubriq,len(rubriq),itysup,nbval,irtyp,  &
         rtp(1,ik),ierror)
    nberro=nberro+ierror

    !            La variable omega sert de tableau de travail
    rubriq = 'R22_ce_phase'//cphase
    call lecsui(impamo,rubriq,len(rubriq),itysup,nbval,irtyp,  &
         rtp(1,iomg),ierror)
    nberro=nberro+ierror

    do iel = 1, ncel
      rtp(iel,ik) = rtp(iel,ik) + rtp(iel,iomg)
    enddo

    rubriq = 'R33_ce_phase'//cphase
    call lecsui(impamo,rubriq,len(rubriq),itysup,nbval,irtyp,  &
         rtp(1,iomg),ierror)
    nberro=nberro+ierror

    do iel = 1, ncel
      rtp(iel,ik) = 0.5d0*(rtp(iel,ik)+rtp(iel,iomg))
    enddo

    rubriq = 'eps_ce_phase'//cphase
    call lecsui(impamo,rubriq,len(rubriq),itysup,nbval,irtyp,  &
         rtp(1,iomg),ierror)
    nberro=nberro+ierror
    !           On transforme ensuite epsilon en omega
    do iel = 1, ncel
      rtp(iel,iomg) = rtp(iel,iomg)/cmu/rtp(iel,ik)
    enddo

    !     * k-omega -> k-omega

  else if(jturph.eq.60) then

    itysup = 1
    nbval  = 1
    irtyp  = 2

    rubriq = 'k_ce_phase'//cphase
    call lecsui(impamo,rubriq,len(rubriq),itysup,nbval,irtyp,   &
         rtp(1,ik),ierror)
    nberro=nberro+ierror

    rubriq = 'omega_ce_phase'//cphase
    call lecsui(impamo,rubriq,len(rubriq),itysup,nbval,irtyp,   &
         rtp(1,iomg),ierror)
    nberro=nberro+ierror

  endif
  !         Rq : laminaire -> k-omega  (On ne fait rien, deja fait dans iniva0)

  !   -- The new computation is with the Spalart Allmaras (SA) model

else if (iturph.eq.70) then

  !     * SA -> SA

  if(jturph.eq.70) then

    itysup = 1
    nbval  = 1
    irtyp  = 2

    rubriq = 'nusa_ce_phase'//cphase
    call lecsui(impamo,rubriq,len(rubriq),itysup,nbval,irtyp,   &
         rtp(1,inusa),ierror)
    nberro=nberro+ierror

  endif

  !TODO perform the conversion from other models to SA.

  !         Rq : laminar -> SA  (We do nothing, it has already been done in iniva0)

  !   -- Le nouveau calcul est en laminaire, longueur de melange ou LES
  !           --> rien a lire

endif

if (nberro.ne.0) then
   write(nfecra,9510)
   call csexit (1)
endif

! --->  Fin de la lecture des variables pression, vitesse, turbulence
write(nfecra,1598)


!  ---> Scalaires

nberro = 0

if (nscal.gt.0) then
  do iscal = 1, nscal
    ivar = isca(iscal)
    ! Si le scalaire existait precedemment on le lit
    ! sinon on ne fait rien (init par defaut dans INIVA0)
    if (iscold(iscal).gt.0) then
      if(iscold(iscal).le.nfmtsc) then
        write(car4,'(i4.4)')iscold(iscal)
        rubriq = 'scalaire_ce_'//car4
        itysup = 1
        nbval  = 1
        irtyp  = 2
        call lecsui(impamo,rubriq,len(rubriq),itysup,nbval,irtyp, &
                    rtp(1,ivar),ierror)
        nberro=nberro+ierror
        rubriq = 'turbulent_flux_model'//car4
        itysup = 0
        nbval  = 1
        irtyp  = 1
        call lecsui(impamo,rubriq,len(rubriq),itysup,nbval,irtyp,jturbt,    &
                    ierror)
        ! If the old calculation has no turbulent flux model, set it to 0
        if (ierror.ne.0) jturbt = 0

        ! ---> Modified data
        if (iturt(iscal) .ne. jturbt) write(nfecra,8411) iturt(iscal), jturbt

      else
        ierror= -1
        nberro=nberro+ierror
      endif
      if(ierror.ne.0) then
        rubrik = rubriq(1:1)
        do ii = 2, min(len(rubriq),nfmtru)
          rubrik = rubrik(1:ii-1)//rubriq(ii:ii)
        enddo
        do ii = min(len(rubriq),nfmtru)+1,nfmtru
          RUBRIK = RUBRIK(1:II-1)//' '
        enddo
        write(nfecra,9511)rubrik
      endif
    endif
    if (ityturt(iscal).eq.2 .or. ityturt(iscal).eq.3) then

      ! Assume name of the scalar in the previous computation has not changed,
      ! so do do apply iscold() here (the restart file provides no information
      ! on the scalar variable field names/ids relation, so we cannot do better
      ! for existing restart files).
      call field_get_name(ivarfl(ivar), fname)
      rubriq = trim(fname)//'_turbulent_flux_ce'

      ! Index of the corresponding turbulent flux
      call field_get_id(trim(fname)//'_turbulent_flux', f_id)
      call field_get_val_v(f_id, xut)

      itysup = 1
      nbval  = 3
      irtyp  = 2
      call lecsui(impamo,rubriq,len(rubriq),itysup,nbval,irtyp,xut,ierror)
      nberro=nberro+ierror

    endif
  enddo
endif

if (nberro.ne.0) then
  write(nfecra,9512)
  call csexit (1)
endif

write(nfecra,1599)


!  ---> Vitesse de maillage en ALE

if (iale.eq.1 .and. jale.eq.1) then

  nberro = 0

  itysup = 1
  nbval  = 1
  irtyp  = 2

  rubriq = 'vit_maillage_u_ce'
  call lecsui(impamo,rubriq,len(rubriq),itysup,nbval,irtyp,       &
              rtp(1,iuma),ierror)
  nberro=nberro+ierror

  rubriq = 'vit_maillage_v_ce'
  call lecsui(impamo,rubriq,len(rubriq),itysup,nbval,irtyp,       &
              rtp(1,ivma),ierror)
  nberro=nberro+ierror

  rubriq = 'vit_maillage_w_ce'
  call lecsui(impamo,rubriq,len(rubriq),itysup,nbval,irtyp,       &
              rtp(1,iwma),ierror)
  nberro=nberro+ierror

  if (nberro.ne.0) then
    write(nfecra,9513)
    do iel = 1, ncel
      rtp(iel,iuma) = 0.d0
      rtp(iel,ivma) = 0.d0
      rtp(iel,iwma) = 0.d0
    enddo
  endif

  write(nfecra,1600)

endif

!================================================================
! 6. LECTURE D'INFORMATIONS COMPLEMENTAIRES LEGERES
!================================================================

if (ichemistry.gt.0.or.iaerosol.gt.0) then
  rubriq = 'atmospheric_chem'
  itysup = 0
  nbval  = 1
  irtyp  = 1
  call lecsui(impamo,rubriq,len(rubriq),itysup,nbval,irtyp,       &
              init_at_chem,ierror)

  if (ierror.eq.0.and.init_at_chem.gt.0) then
    init_at_chem = 0
  endif
endif

!================================================================
! 7. FERMETURE DU FICHIER SUITE PRINCIPAL
!================================================================

call clssui(impamo,ierror)
!==========

if (ierror.ne.0) then
  write(nfecra,8711) ficsui
endif

write(nfecra,1799)

!===============================================================================
! 8. SORTIE
!===============================================================================

return

!===============================================================================
! 9. FORMATS
!===============================================================================

! --- ETAPES

#if defined(_CS_LANG_FR)

 1000 format(/,                                                   &
     3X,'   LECTURE DU FICHIER SUITE PRINCIPAL                ',/)
 1100 format(' Debut de la lecture                                  '  )
 1299 format('   Fin de la lecture des dimensions                   '  )
 1499 format('   Fin de la lecture des options                      '  )
 1598 format('   Fin de la lecture des variables                    ',/,&
       '        pression, vitesse, turbulence                 '  )
 1599 format('   Fin de la lecture des scalaires                    '  )
 1600 format('   Fin de la lecture de la vitesse de maillage        ',/,&
       '                                  (methode ALE)       '  )
 1799 format(' Fin de la lecture                                    '  )

#else

 1000 format(/,                                                   &
     3X,'   READING THE MAIN RESTART FILE                     ',/)
 1100 format(' Start reading                                        '  )
 1299 format('   Reading dimensions complete                        '  )
 1499 format('   Reading options complete                           '  )
 1598 format('   Reading the pressure, velocity, turbulent          ',/,&
       '        variables complete                            '  )
 1599 format('   Reading scalars complete                           '  )
 1600 format('   Reading the mesh velocity (ALE method) complete    '  )
 1799 format(' Reading complete                                     '  )

#endif

! --- INFORMATIONS

#if defined(_CS_LANG_FR)

 2410 format                                                            &
 ('   Lecture du pas de temps precedent (suite) ',                &
                                                  'NTPABS = ',I10)
 2411 format                                                            &
 ('   Lecture du pas de temps precedent (suite) ',                &
                                                'TTPABS = ',E12.4)
 2412 format                                                            &
 ('   Lecture du temps de maillage mobile precedent (suite) ',    &
                                                'TTPMOB = ',E12.4)

#else

 2410 format                                                            &
 ('   Reading the previous time step number ',                    &
                      '(restarting computation)  NTPABS =   ',I10)
 2411 format                                                            &
 ('   Reading the previous time step number ',                    &
                      '(restarting computation)  TTPABS = ',E12.4)
 2412 format                                                            &
 ('   Reading the previous moving mesh moment ',                  &
                      '(restarting computation)  TTPMOB = ',E12.4)

#endif


! --- MISES EN GARDE

#if defined(_CS_LANG_FR)

 8001 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION :       A LA LECTURE DU FICHIER SUITE         ',/,&
'@    =========                                      PRINCIPAL',/,&
'@                                                            ',/,&
'@      Le nombre de scalaires maximal NSCAMX supporte pas le ',/,&
'@        format d''ecriture du fichier suite est             ',/,&
'@        NFMTSC = ',I10                                       ,/,&
'@      On a ici un nombre de scalaires maximal superieur     ',/,&
'@        NSCAMX = ',I10                                       ,/,&
'@      On ne pourra pas relire les scalaires dont le numero  ',/,&
'@        est superieur                                       ',/,&
'@                                                            ',/,&
'@    Le calcul sera execute.                                 ',/,&
'@                                                            ',/,&
'@    Voir le sous-programme lecamp.                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8205 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : LECTURE DU FICHIER SUITE PRINCIPAL          ',/,&
'@    =========                                               ',/,&
'@      DONNEES AMONT MULTIPHASIQUES                          ',/,&
'@                                                            ',/,&
'@  Nombre de phases (amont) : ',I10                           ,/,&
'@                                                            ',/,&
'@    Le calcul ne peut etre execute.                         ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8210 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : LECTURE DU FICHIER SUITE PRINCIPAL          ',/,&
'@    =========                                               ',/,&
'@      DONNEES AMONT ET ACTUELLES DIFFERENTES                ',/,&
'@                                                            ',/,&
'@    Le nombre de variables ou de scalaires a ete modifie.   ',/,&
'@    Le calcul peut etre execute.                            ',/,&
'@                                                            ',/,&
'@    Il est cependant conseille de verifier                  ',/,&
'@      les dimensions suivantes :                            ',/,&
'@                                                            ',/,&
'@                NVAR     NSCAL    NSCAUS    NSCAPP          ',/,&
'@  AMONT : ',4I10                                             ,/,&
'@  ACTUEL: ',4I10                                             ,/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8410 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : LECTURE DU FICHIER SUITE PRINCIPAL          ',/,&
'@    =========                                               ',/,&
'@                                                            ',/,&
'@      REPRISE  DE CALCUL           AVEC ITURB = ',I4         ,/,&
'@      A PARTIR D''UN CALCUL REALISE AVEC ITURB = ',I4        ,/,&
'@                                                            ',/,&
'@    Le modele de turbulence en temps a ete modifie.         ',/,&
'@    Le calcul peut etre execute.                            ',/,&
'@                                                            ',/,&
'@    Il est conseille cependant de                           ',/,&
'@      verifier la valeur de ITURB(',I2,')                   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8411 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : LECTURE DU FICHIER SUITE PRINCIPAL          ',/,&
'@    =========                                               ',/,&
'@                                                            ',/,&
'@      REPRISE  DE CALCUL           AVEC ITURB = ',I4         ,/,&
'@      A PARTIR D''UN CALCUL REALISE AVEC ITURB = ',I4        ,/,&
'@                                                            ',/,&
'@    Le modele de flux turbulent a ete modifie.              ',/,&
'@    Le calcul peut etre execute.                            ',/,&
'@                                                            ',/,&
'@    Il est conseille cependant de                           ',/,&
'@      verifier la valeur de ITURB(',I2,')                   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8711 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ERREUR A LA FERMETURE DU FICHIER SUITE      ',/,&
'@    =========                                      PRINCIPAL',/,&
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
'@ @@ WARNING :  WHEN READING THE MAIN RESTART FILE        ',/,   &
'@    =========                                               ',/,&
'@                                                            ',/,&
'@      The maximum number of scalars NSCAMX supported        ',/,&
'@        by the restart file writing format is               ',/,&
'@        NFMTSC = ',I10                                       ,/,&
'@      Here is the maximum number of scalars larger          ',/,&
'@        NSCAMX = ',I10                                       ,/,&
'@      Scalars whose number is larger wont be read.          ',/,&
'@                                                            ',/,&
'@    The calculation will be run.                            ',/,&
'@                                                            ',/,&
'@    See the sub-routine lecamp                          ',/,    &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8205 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING : WHEN READING THE MAIN RESTARTING FILE       ',/,  &
'@    =========                                               ',/,&
'@      CHECKPOINT DATA ARE MULTIPHASE                        ',/,&
'@                                                            ',/,&
'@  Number of phases (checkpoint) : ',I10                      ,/,&
'@                                                            ',/,&
'@    The computation cannot be executed.                     ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8210 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING : WHEN READING THE MAIN RESTARTING FILE       ',/,  &
'@    =========                                               ',/,&
'@      THE RESTART AND CHECKPOINT DATA ARE DIFFERENT         ',/,&
'@                                                            ',/,&
'@    The number of variables or scalars has changed.         ',/,&
'@    The calculation will be run.                            ',/,&
'@                                                            ',/,&
'@    However, it is strongly advised to check                ',/,&
'@      the following dimensions:                             ',/,&
'@                                                            ',/,&
'@                NVAR     NSCAL    NSCAUS    NSCAPP          ',/,&
'@ PREVIOUS:',4I10                                             ,/,&
'@ CURRENT :',4I10                                             ,/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8410 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING : WHEN READING THE MAIN RESTART FILE            ',/,&
'@    =========                                               ',/,&
'@                                                            ',/,&
'@    THE CURRENT CALCULATION USES ITURB = ',I4                ,/,&
'@    BUT RESTARTS FROM ANOTHER ONE USING ITURB = ',I4         ,/,&
'@                                                            ',/,&
'@    The turbulence model has changed.                       ',/,&
'@    The computation can be executed.                        ',/,&
'@                                                            ',/,&
'@    However, it is strongly advised to check                ',/,&
'@      the value of the variable ITURB(',I2,')               ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8411 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING : WHEN READING THE MAIN RESTART FILE            ',/,&
'@    =========                                               ',/,&
'@                                                            ',/,&
'@    THE CURRENT CALCULATION USES ITURB = ',I4                ,/,&
'@    BUT RESTARTS FROM ANOTHER ONE USING ITURB = ',I4         ,/,&
'@                                                            ',/,&
'@    The turbulent flux model has changed.                   ',/,&
'@    The computation can be executed.                        ',/,&
'@                                                            ',/,&
'@    However, it is strongly advised to check                ',/,&
'@      the value of the variable ITURB(',I2,')               ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8711 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING : ERROR AT CLOSING THE MAIN RESTART FILE        ',/,&
'@    =========                                               ',/,&
'@                                                            ',/,&
'@    Problem on file called (',A13,')                        ',/,&
'@                                                            ',/,&
'@    The calculation will be run.                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

#endif

! --- ERREURS

#if defined(_CS_LANG_FR)

 9100 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,&
'@    =========                                               ',/,&
'@      ERREUR A L''OUVERTURE DU FICHIER SUITE PRINCIPAL      ',/,&
'@                                                            ',/,&
'@    Le calcul ne peut pas etre execute.                     ',/,&
'@                                                            ',/,&
'@    Verifier l''existence et le nom (',A13,') du            ',/,&
'@        fichier suite dans le repertoire de travail.        ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9200 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,&
'@    =========                                      PRINCIPAL',/,&
'@      TYPE DE FICHIER INCORRECT                             ',/,&
'@                                                            ',/,&
'@    Le fichier ',A13      ,' ne semble pas etre un fichier  ',/,&
'@      suite principal.                                      ',/,&
'@                                                            ',/,&
'@    Le calcul ne peut etre execute.                         ',/,&
'@                                                            ',/,&
'@    Verifier que le fichier suite utilise correspond bien   ',/,&
'@        a un fichier suite principal.                       ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9201 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,&
'@    =========                                      PRINCIPAL',/,&
'@      DONNEES AMONT ET ACTUELLES INCOHERENTES               ',/,&
'@                                                            ',/,&
'@    Le nombre de cellules a ete modifie                     ',/,&
'@                                                            ',/,&
'@    Le calcul ne peut etre execute.                         ',/,&
'@                                                            ',/,&
'@    Verifier que le fichier suite utilise correspond bien   ',/,&
'@        au cas traite.                                      ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9210 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,&
'@    =========                                      PRINCIPAL',/,&
'@      ERREUR A LA LECTURE DES DIMENSIONS                    ',/,&
'@                                                            ',/,&
'@    Le calcul ne peut etre execute.                         ',/,&
'@                                                            ',/,&
'@    Verifier que le fichier suite utilise n''a pas ete      ',/,&
'@        endommage.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9320 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,&
'@    =========                                      PRINCIPAL',/,&
'@                                                            ',/,&
'@    On souhaite faire correspondre le scalaire  ',I10        ,/,&
'@            du present calcul avec le scalaire  ',I10        ,/,&
'@            du calcul precedent, or, le nombre maximal de   ',/,&
'@            scalaires dans le fichier suite est ',I10        ,/,&
'@                                                            ',/,&
'@    Le calcul ne peut etre execute.                         ',/,&
'@                                                            ',/,&
'@    Verifier les valeurs de ISCOLD.                         ',/,&
'@    Verifier que le fichier suite utilise correspond bien   ',/,&
'@        au cas traite.                                      ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9321 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,&
'@    =========                                      PRINCIPAL',/,&
'@                                                            ',/,&
'@    On souhaite faire correspondre le scalaire  ',I10        ,/,&
'@            du present calcul avec le scalaire  ',I10        ,/,&
'@            du calcul precedent, or, le numero des anciens  ',/,&
'@            scalaires doit etre strictement positif et      ',/,&
'@            inferieur ou egal a                 ',I10        ,/,&
'@                                                            ',/,&
'@    Le calcul ne peut etre execute.                         ',/,&
'@                                                            ',/,&
'@    Verifier les valeurs de ISCOLD.                         ',/,&
'@    Verifier que le fichier suite utilise correspond bien   ',/,&
'@        au cas traite.                                      ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9330 format (                                                          &
'@                                                            ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,&
'@    =========                                      PRINCIPAL',/,&
'@                                                            ',/,&
'@      CORRESPONDANCE AVEC LES ANCIENS SCALAIRES IMPOSSIBLE  ',/,&
'@      ',I10,   ' ERREURS REPORTEES CI-DESSUS                ',/,&
'@                                                            ',/,&
'@    Le calcul ne peut etre execute.                         ',/,&
'@                                                            ',/,&
'@      Les lignes precedentes identifient les scalaires      ',/,&
'@        definis dans le present calcul pour lesquels la     ',/,&
'@        la correspondance fournie par ISCOLD                ',/,&
'@        est incorrecte.                                     ',/,&
'@                                                            ',/,&
'@      Le calcul ne peut etre execute.                       ',/,&
'@                                                            ',/,&
'@      Verifier que le fichier suite UTILISE correspond bien ',/,&
'@        au cas traite.                                      ',/,&
'@      Verifier ISCOLD.                                      ',/,&
'@        ISCOLD(ISCAL) = 0 indique que le scalaire ISCAL     ',/,&
'@          du present calcul ne correspond a aucun scalaire  ',/,&
'@          du calcul precedent                               ',/,&
'@        ISCOLD(ISCAL) > 0 indique le numero du scalaire     ',/,&
'@          du calcul precedent auquel correspond ISCAL       ',/,&
'@      Si les correspondances avec les anciens scalaires     ',/,&
'@        ne sont pas necessaires, ne pas definir ISCOLD.     ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9400 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,&
'@    =========                                      PRINCIPAL',/,&
'@                                                            ',/,&
'@      ERREUR A LA LECTURE DES INFORMATIONS TEMPORELLES OU   ',/,&
'@        DES MODELES DE TURBULENCE                           ',/,&
'@                                                            ',/,&
'@    Le calcul ne peut pas etre execute.                     ',/,&
'@                                                            ',/,&
'@    Verifier que le fichier suite utilise n''a pas ete      ',/,&
'@        endommage.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9401 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ERREUR A LA LECTURE DU FICHIER SUITE        ',/,&
'@    =========                                      PRINCIPAL',/,&
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
 9402 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,&
'@    =========                                      PRINCIPAL',/,&
'@      INDICATEUR IALE DU CALCUL PRECEDENT = ',I10            ,/,&
'@      INDICATEUR IALE DU CALCUL ACTUEL    = ',I10            ,/,&
'@                                                            ',/,&
'@    Les coordonnees des noeuds du maillage doivent etre     ',/,&
'@      relues. Elles sont stockees dans le fichier suite     ',/,&
'@      auxiliaire.                                           ',/,&
'@    L''indicateur ILEAUX doit donc etre positionne a 1.     ',/,&
'@    Il vaut ici ILEAUX = ',I10                               ,/,&
'@                                                            ',/,&
'@    Le calcul ne peut pas etre execute.                     ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@    Verifier ILEAUX.                                        ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9403 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ERREUR A LA LECTURE DU FICHIER SUITE        ',/,&
'@    =========                                      PRINCIPAL',/,&
'@                                                            ',/,&
'@      ERREUR A LA LECTURE DE L INSTANT DE MAILLAGE MOBILE  ',/,&
'@                                                   PRECEDENT',/,&
'@    Il se peut que le fichier suite relu corresponde a une  ',/,&
'@      version anterieure de Code_Saturne, sans couplage     ',/,&
'@      rotor/stator instationnaire.                          ',/,&
'@    Le calcul sera execute en initialisant l instant de     ',/,&
'@      maillage mobile precedent a TTCMOB = ',E12.4           ,/,&
'@    Verifier neanmoins que le fichier suite utilise n''a    ',/,&
'@        pas ete endommage.                                  ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9410 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,&
'@    =========                                      PRINCIPAL',/,&
'@      NUMERO DU PAS DE TEMPS PRECEDENT NTPABS = ',I10        ,/,&
'@      NUMERO DU PAS DE TEMPS VISE      NTMABS = ',I10        ,/,&
'@                                                            ',/,&
'@    Le nombre de pas de temps (absolu) vise, NTMABS,        ',/,&
'@      doit etre superieur au nombre de pas de temps         ',/,&
'@      (absolu) deja effectues, NTPABS.                      ',/,&
'@                                                            ',/,&
'@    Le calcul ne peut etre execute.                         ',/,&
'@                                                            ',/,&
'@    Verifier (augmenter) NTMABS.                            ',/,&
'@    Verifier que le fichier suite utilise correspond bien   ',/,&
'@        au cas traite.                                      ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9510 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,&
'@    =========                                      PRINCIPAL',/,&
'@                                                            ',/,&
'@      ERREUR LORS DE LA LECTURE DES VARIABLES PRINCIPALES   ',/,&
'@        PRESSION, VITESSE, TURBULENCE                       ',/,&
'@                                                            ',/,&
'@    Le calcul ne peut etre execute.                         ',/,&
'@                                                            ',/,&
'@    Verifier que le fichier suite utilise n''a pas ete      ',/,&
'@        endommage.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
! NFMTRU = 36 pour A36
 9511 format(                                                     &
'@ Erreur a la lecture de ',A36                                  )
 9512 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,&
'@    =========                                      PRINCIPAL',/,&
'@                                                            ',/,&
'@      ERREUR LORS DE LA LECTURE DES VARIABLES PRINCIPALES   ',/,&
'@        SCALAIRES                                           ',/,&
'@                                                            ',/,&
'@    Le calcul ne peut etre execute.                         ',/,&
'@                                                            ',/,&
'@    Verifier que le fichier suite utilise n''a pas ete      ',/,&
'@        endommage.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9513 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,&
'@    =========                                      PRINCIPAL',/,&
'@                                                            ',/,&
'@      ERREUR LORS DE LA LECTURE DE LA VITESSE DE MAILLAGE   ',/,&
'@        (METHODE ALE)                                       ',/,&
'@                                                            ',/,&
'@    Verifier que le fichier suite utilise n''a pas ete      ',/,&
'@        endommage.                                          ',/,&
'@                                                            ',/,&
'@    Le calcul peut neanmoins etre execute (vitesse de       ',/,&
'@                                     maillage reinitialisee)',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

#else

 9100 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING : STOP AT THE RESTART FILE READING              ',/,&
'@    =========                                               ',/,&
'@      ERROR AT OPENING THE MAIN RESTART FILE                ',/,&
'@                                                            ',/,&
'@    The computation cannot be executed.                     ',/,&
'@                                                            ',/,&
'@    Please check the existence and the name (',A13,')       ',/,&
'@        of the restart file in the working directory.       ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9200 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING : STOP AT THE MAIN RESTART FILE READING         ',/,&
'@    =========                                               ',/,&
'@      WRONG FILE TYPE                                       ',/,&
'@                                                            ',/,&
'@    The file ',A13      ,' does not look like a proper      ',/,&
'@      main restart file.                                    ',/,&
'@                                                            ',/,&
'@    The calculation cannot be executed.                     ',/,&
'@                                                            ',/,&
'@    Please make sure the file used as a restart file        ',/,&
'@        actually is a correct main restart file.            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9201 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING : STOP AT THE MAIN RESTART FILE READING         ',/,&
'@    =========                                               ',/,&
'@      INCONSISTANT RESTART AND CHECKPOINT DATA              ',/,&
'@                                                            ',/,&
'@    The number of cells has changed                         ',/,&
'@                                                            ',/,&
'@    The calculation cannot be executed.                     ',/,&
'@                                                            ',/,&
'@    Please make sure the file used as restart file does     ',/,&
'@        correspond to your case                             ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9210 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING : STOP AT THE RESTART FILE READING              ',/,&
'@    =========                                               ',/,&
'@      ERROR AT READING DIMENSIONS                           ',/,&
'@                                                            ',/,&
'@    The calculation cannot be executed.                     ',/,&
'@                                                            ',/,&
'@    Please check the integrity of the file used as          ',/,&
'@        restart file                                        ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9320 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING : STOP AT THE MAIN RESTART FILE READING         ',/,&
'@    =========                                               ',/,&
'@                                                            ',/,&
'@    The scalar ',I10   ,' of the current calculation        ',/,&
'@        is attempted to match the scalar      ',I10          ,/,&
'@        of the previous one. But the maximum number         ',/,&
'@        of scalar read in the restart file is ',I10          ,/,&
'@                                                            ',/,&
'@    The calculation cannot be excuted.                      ',/,&
'@                                                            ',/,&
'@    Please check the value of ISCOLD.                       ',/,&
'@    Please make sure the file used as restart file does     ',/,&
'@        correspond to your case                             ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9321 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING : STOP AT THE MAIN RESTART FILE READING         ',/,&
'@    =========                                               ',/,&
'@                                                            ',/,&
'@    The scalar ',I10   ,' of the current calculation        ',/,&
'@        is attempted to match the scalar ',I10               ,/,&
'@        of the previous one. But the number                 ',/,&
'@        of the previous scalar must be stricly larger       ',/,&
'@        than zero and smaller than or equal to     ',I10     ,/,&
'@                                                            ',/,&
'@    The calculation cannot be executed.                     ',/,&
'@                                                            ',/,&
'@    Please check the value of ISCOLD.                       ',/,&
'@    Please make sure the file used as restart file does     ',/,&
'@        correspond to your case                             ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9330 format (                                                          &
'@                                                            ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING : STOP AT THE MAIN RESTART FILE READING         ',/,&
'@    =========                                               ',/,&
'@                                                            ',/,&
'@      IMPOSSIBLE MATCHING WITH PREVIOUS SCALAR              ',/,&
'@      ',I10,   ' ERRORS REPORTED IN THE ABOVE MESSAGES      ',/,&
'@                                                            ',/,&
'@    The calculation cannot be executed.                     ',/,&
'@                                                            ',/,&
'@      The previous lines identify the scalars defined      ',/, &
'@        in the current computation for which the matching   ',/,&
'@        defined by ISCOLD is not correct                    ',/,&
'@                                                            ',/,&
'@      The computation cannot be executed.                   ',/,&
'@                                                            ',/,&
'@      Please make sure the file used as restart file does   ',/,&
'@          correspond to your case                           ',/,&
'@      Please check the value of ISCOLD.                     ',/,&
'@        ISCOLD(ISCAL) = 0 means that the scalar ISCAL       ',/,&
'@          of the current calculation does not correspond to ',/,&
'@          any of the scalars of the previous calculation    ',/,&
'@        ISCOLD(ISCAL) > 0 defines the scalar number in the  ',/,&
'@          previous calculation ISCAL corresponds to         ',/,&
'@      If matching with the previous scalar is unnecessary   ',/,&
'@        please do not define ISCOLD.                        ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9400 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING : STOP AT THE MAIN RESTART FILE READING         ',/,&
'@    =========                                               ',/,&
'@                                                            ',/,&
'@      ERROR AT READING THE TEMPORAL INFORMATION OR          ',/,&
'@        TURBULENCE MODELS                                   ',/,&
'@                                                            ',/,&
'@    The computation cannot be executed.                     ',/,&
'@                                                            ',/,&
'@    Please check the integrity of the file used as          ',/,&
'@        restart file                                        ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9401 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING : ERROR AT THE MAIN RESTART FILE READING        ',/,&
'@    =========                                               ',/,&
'@                                                            ',/,&
'@      ERROR AT READING THE INDICATOR OF ALE METHOD          ',/,&
'@                                                            ',/,&
'@    The read restart file might come from a previous        ',/,&
'@      version of Code Saturne, without ALE.                 ',/,&
'@    The calculation will be executed but                    ',/,&
'@      ALE data will be reset.                               ',/,&
'@    Please check the integrity of the file used as          ',/,&
'@        restart file, however.                              ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9402 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING : STOP AT THE MAIN RESTART FILE READING         ',/,&
'@    =========                                               ',/,&
'@      IALE INDICATOR OF THE PREVIOUS CALCULATION = ',I10     ,/,&
'@      IALE INDICATOR OF THE CURRECT CALCULATION  = ',I10     ,/,&
'@                                                            ',/,&
'@    The coordinates of the mesh nodes need to be read again.',/,&
'@      They are stored in the auxiliary restart file.        ',/,&
'@    Therefore the ILEAUX indicator needs to be equal to 1.  ',/,&
'@    Its current value is ILEAUX = ',I10                     ,/, &
'@                                                            ',/,&
'@    The calculation cannot be executed.                     ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@    Please check the value of ILEAUX.                       ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9403 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING : ERROR AT THE MAIN RESTART FILE READING        ',/,&
'@    =========                                               ',/,&
'@                                                            ',/,&
'@      ERROR AT READING THE PREVIOUS MOVING MESH MOMENT      ',/,&
'@                                                            ',/,&
'@    The read restart file might come from a previous        ',/,&
'@      version of Code Saturne, without unsteady             ',/,&
'@      rotor/stator coupling method.                         ',/,&
'@    The calculation will be executed with the previous      ',/,&
'@      moving mesh moment initialized to TTCMOB = ',E12.4     ,/,&
'@    Please check the integrity of the file used as          ',/,&
'@        restart file, however.                              ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9410 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING : STOP AT THE MAIN RESTART FILE READING         ',/,&
'@    =========                                               ',/,&
'@      NUMBER OF THE PREVIOUS TIME STEP  NTPABS = ',I10       ,/,&
'@      NUMBER OF TIME STEPS WANTED       NTMABS = ',I10       ,/,&
'@                                                            ',/,&
'@    The number of time steps (absolute) wanted, NTMABS,     ',/,&
'@      has to be larger than the number of time steps        ',/,&
'@      (absolute) already done, NTPABS.                      ',/,&
'@                                                            ',/,&
'@    The calculation cannot be executed.                     ',/,&
'@                                                            ',/,&
'@    Please check (increase) NTMABS.                         ',/,&
'@    Please make sure the file used as restart file does     ',/,&
'@          correspond to your case                           ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9510 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING : STOP AT THE MAIN RESTART FILE READING         ',/,&
'@    =========                                               ',/,&
'@                                                            ',/,&
'@      ERROR AT THE READING OF MAIN VARIABLES                ',/,&
'@        PRESSURE, VELOCITY, TURBULENCE                      ',/,&
'@                                                            ',/,&
'@    The calculation cannot be executed.                     ',/,&
'@                                                            ',/,&
'@    Please check the integrity of the file used as          ',/,&
'@        restart file                                        ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9511 format(                                                     &
'@ Error at the reading of ',A36                                 )
 9512 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING : STOP AT THE MAIN RESTART FILE READING         ',/,&
'@    =========                                               ',/,&
'@                                                            ',/,&
'@      ERROR AT THE READING OF MAIN SCALAR VARIABLES         ',/,&
'@                                                            ',/,&
'@    The calculation cannot be executed.                     ',/,&
'@                                                            ',/,&
'@    Please check the integrity of the file used as          ',/,&
'@        restart file                                        ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9513 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING : STOP AT THE MAIN RESTART FILE READING         ',/,&
'@    =========                                               ',/,&
'@                                                            ',/,&
'@      ERROR AT THE READING OF MESH VELOCITY                 ',/,&
'@        (ALE METHOD)                                       ',/, &
'@                                                            ',/,&
'@    Please check the integrity of the file used as          ',/,&
'@        restart file                                        ',/,&
'@                                                            ',/,&
'@    The calculation can be executed however (mesh           ',/,&
'@                                     velocity reset)        ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

#endif

end subroutine
