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

subroutine mttsns &
!================

 ( idbia0 , idbra0 ,                                              &
   nvar   , nscal  , nphas  , ncepdp , ncesmp ,                   &
   ivar   , iphas  ,                                              &
   icepdc , icetsm , itypsm ,                                     &
   ia     ,                                                       &
   dt     , rtpa   , propce , propfa , propfb ,                   &
   coefa  , coefb  , ckupdc , smacel ,                            &
   crvexp , crvimp ,                                              &
   dam    , xam    ,                                              &
   w1     , w2     , w3     , w4     , w5     , w6     ,          &
   ra     )

!===============================================================================
!C FONCTION :
! ----------

! CALCUL DES TERMES SOURCES POUR LA COMPOSANTE DE VITESSE IVAR

!    SOUS-PROGRAMME SPECIFIQUE A MATISSE (COPIE DE USTSNS)


! LES TERMES SOURCES SONT LES PERTES QUI N'ONT PAS ETE TRAITEES
!    DANS MTKPDC (I.E. TOUT SAUF LES REGISTRES)

! IL SERAIT POSSIBLE DE REGROUPER MTTSNS ET MTKPDC EN UN SEUL
!    SOUS-PROGRAMME.


! ON RESOUT RHO*VOLUME*D(VAR)/DT = CRVIMP*VAR + CRVEXP

! ON FOURNIT ICI CRVIMP ET CRVEXP (ILS CONTIENNENT RHO*VOLUME)
!    CRVEXP en kg variable/s :
!     ex : pour la vitesse            kg m/s2
!    CRVIMP en kg /s :

! VEILLER A UTILISER UN CRVIMP NEGATIF
! (ON IMPLICITERA CRVIMP
!  IE SUR LA DIAGONALE DE LA MATRICE, LE CODE AJOUTERA :
!   MAX(-CRVIMP,0) EN SCHEMA STANDARD EN TEMPS
!       -CRVIMP    SI LES TERMES SOURCES SONT A L'ORDRE 2

! CES TABLEAUX SONT INITIALISES A ZERO AVANT APPEL A CE SOUS
!   PROGRAMME ET AJOUTES ENSUITE AUX TABLEAUX PRIS EN COMPTE
!   POUR LA RESOLUTION

! EN CAS D'ORDRE 2 DEMANDE SUR LES TERMES SOURCES, ON DOIT
!   FOURNIR CRVEXP A L'INSTANT N     (IL SERA EXTRAPOLE) ET
!           CRVIMP A L'INSTANT N+1/2 (IL EST  DANS LA MATRICE,
!                                     ON LE SUPPOSE NEGATIF)

!-------------------------------------------------------------------------------
!ARGU                             ARGUMENTS
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! nphas            ! i  ! <-- ! number of phases                               !
! ncepdp           ! i  ! <-- ! number of cells with head loss                 !
! ncesmp           ! i  ! <-- ! number of cells with mass source term          !
! ivar             ! i  ! <-- ! variable number                                !
! iphas            ! i  ! <-- ! phase number                                   !
! icepdc(ncelet    ! te ! <-- ! numero des ncepdp cellules avec pdc            !
! icetsm(ncesmp    ! te ! <-- ! numero des cellules a source de masse          !
! itypsm           ! te ! <-- ! type de source de masse pour les               !
! (ncesmp,nvar)    !    !     !  variables (cf. ustsma)                        !
! ia(*)            ! ia ! --- ! main integer work array                        !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtpa             ! tr ! <-- ! variables de calcul au centre des              !
! (ncelet,*)       !    !     !    cellules (instant            prec)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
! ckupdc           ! tr ! <-- ! tableau de travail pour pdc                    !
!  (ncepdp,6)      !    !     !                                                !
! smacel           ! tr ! <-- ! valeur des variables associee a la             !
! (ncesmp,*   )    !    !     !  source de masse                               !
!                  !    !     !  pour ivar=ipr, smacel=flux de masse           !
! crvexp(ncelet    ! tr ! --> ! tableau de travail pour part explicit          !
! crvimp(ncelet    ! tr ! --> ! tableau de travail pour part implicit          !
! dam(ncelet       ! tr ! --- ! tableau de travail pour matrice                !
! xam(nfac,*)      ! tr ! --- ! tableau de travail pour matrice                !
! w1...6(ncelet    ! tr ! --- ! tableau de travail                             !
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
use pointe
use numvar
use entsor
use optcal
use cstphy
use cstnum
use parall
use period
use matiss
use mesh

!===============================================================================

implicit none

! Arguments

integer          idbia0 , idbra0
integer          nvar   , nscal  , nphas
integer          ncepdp , ncesmp
integer          ivar   , iphas

integer          icepdc(ncepdp)
integer          icetsm(ncesmp), itypsm(ncesmp,nvar)
integer          ia(*)

double precision dt(ncelet), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision ckupdc(ncepdp,6), smacel(ncesmp,nvar)
double precision crvexp(ncelet), crvimp(ncelet)
double precision dam(ncelet ),xam(nfac ,2)
double precision w1(ncelet),w2(ncelet),w3(ncelet)
double precision w4(ncelet),w5(ncelet),w6(ncelet)
double precision ra(*)

! Local variables

character*80     chaine
integer          idebia, idebra
integer          ipcrom, ipp
integer          iel   , ifml  , icoul
integer          indixy, izone , nzones
integer          modntl
double precision xlimin, xlimax, yrgmin, yrgmax, zalmin, zalmax
double precision xzomin, xzomax, yzomin, yzomax
double precision xtzc  , xlzc  , xhzc  , xlce  , xlcs
double precision salv  , xlalv , dhplen
double precision cxsd  , cysd  , ctsd  , rs0   , gapscx
double precision rapsrf, smesr , smssr
double precision fpfl  , qlqp
double precision roe   , we0   , ve0
double precision ros   , ws0   , vs0   , ts0
double precision rom   , w0    , walv  , un0   , unh   , ureel
double precision reye0 , reys0 , reyalv, reyres, reyhor
double precision richar, betmat, gmat
double precision teta0 , sde   , sds   , xxaa  , xxem
double precision acoeff, aprim
double precision cpe   , cps   , cpc   , pdce  , pdcs
double precision vamcaa, vamhaa, pdcfaa, pdcbaa, pdcraa, cmptaa

!===============================================================================
!===============================================================================
! 1. INITIALISATIONS
!===============================================================================

! --- Gestion memoire

idebia = idbia0
idebra = idbra0

! --- Reperage des variables

ipcrom = ipproc(irom  (iphas))

! --- Reperage du pas de temps pour les impressions listing

if(ntlist.gt.0) then
  modntl = mod(ntcabs,ntlist)
elseif(ntlist.eq.-1.and.ntcabs.eq.ntmabs) then
  modntl = 0
else
  modntl = 1
endif

! --- Impressions

ipp    = ipprtp(ivar)
chaine = nomvar(ipp)
if(iwarni(ivar).ge.1) then
  write(nfecra,1000) chaine(1:8)
endif
if ((irangp.le.0).and.(modntl.eq.0)) then
  write(nfecra,2001)chaine(1:8)
endif

! --- Compteurs pour moyennes et affichages

!     Vitesse alveole moyenne corrigee (reelle)
vamcaa = 0.d0
!     Perte de charge moyenne alveole par frottement (Pa)
pdcfaa = 0.d0
!     Perte de charge moyenne alveole par bifurcation laterale (Pa)
pdcbaa = 0.d0
!     Perte de charge moyenne alveole par alimentation (reunification) (Pa)
pdcraa = 0.d0
!     Vitesse alveole moyenne homogeneisee
vamhaa = 0.d0
!     Compteur pour moyenne
cmptaa = 0.d0


!===============================================================================
! 1. DIMENSIONS GEOMETRIQUES
!===============================================================================


! --- Dimensions de la zone de stockage
!     (Transverse (x), Longitudinal (y), Hauteur (z))

xtzc = ptrres*nptran
xlzc = plgres*nplgrs
xhzc = epchel*nchest


! --- Section debitante des cheminees (SDE : Entree Amont et SDS : Sortie)
!       FRDTRA, le facteur de reduction transverse, est le facteur
!         par lequel il faut diviser les sections d'entree ou les
!         debits si l'on represente un entrepot (reel, 3D) par un
!         calcul 2D.
!       SDCHEA et SDCHES sont les sections passantes reelles en entree
!         et sortie, au dessus des convergents eventuels.
!       SDE    et SDS    sont les sections passantes correspondantes
!         a prendre en compte dans le calcul (sections reelles divisees
!         par FRDTRA, rendant compte d'une reduction eventuelle du domaine
!         dans la direction transverse par rapport a la realite).
!         Neanmoins, SDE et SDS ne sont pas necessairement les sections
!         passantes geometriquement representees dans le maillage (voir
!         les rapports SMESR et SMSSR).

sde = sdchea/frdtra
sds = sdches/frdtra


! --- Longueur cheminee
!     (XLCE : cheminee d'Entree sans le convergent,
!      XLCS : cheminee de Sortie sans le convergent ou sans le toit)

xlce = hchali - hconve
if (itypen.eq.1) then
  xlcs = hcheva - hfttoi
else
  xlcs = hcheva - hconve
endif


! --- Surface debitante d'une alveole (entre chemise et colis)

salv = pi * abs(dhalve**2 - dmcont**2)*0.25d0


! --- Rapport des sections modele/voulu : SMESR et SMSSR

!       On multiplie les vitesses issues directement du calcul
!         dans les cheminees (au dessus des convergents eventuels)
!         par SMESR ou SMSRS pour obtenir des vitesses reelles.
!       SMESR est associe a la cheminee d'entree et SMSSR a la
!         cheminee de sortie.
!       SMESR et SMSSR different de FRDTRA (voir FRDTRA plus haut).

!       SMESR (resp. SMSSR) est le rapport entre la section passante
!         de la cheminee d'entree (resp. de sortie) telle qu'elle est
!         representee geometriquement dans le modele et la section
!         passante que l'on souhaite a modeliser (les sections sont
!         prises au dessus du convergent eventuel).
!       La section representee geometriquement peut etre differente de
!         la section que l'on souhaite modeliser par exemple si les
!         cheminees reelles sont de section circulaire alors qu'on les
!         modelise de section carree.
!       La section que l'on souhaite modeliser peut etre differente de
!         la section reelle si l'on realise un calcul 2D pour
!         representer une configuration reelle 3D (rapport FRDTRA).

!      Pour un cas modelise en 3D :
!         FRDTRA = 1, i.e. SDE=SDCHEA : la section que l'on souhaite
!           modeliser est egale a la section reelle,
!         on a, pour la cheminee d'entree (et de meme pour la cheminee
!           de sortie) :
!             SMESR = EPCHEM*XTZC/(RCONVE*SDCHEA)
!         ceci est bien le rapport des sections passantes
!         "(representee geometriquement)/(que l'on souhaite modeliser)"
!         (les sections sont prises au dessus du convergent eventuel).

!       Pour un cas modelise en 2D :
!         FRDTRA non unite, i.e. SDE=SDCHEA/FRDTRA : la section que
!           l'on souhaite modeliser est differente de la section reelle
!           (en general la section que l'on souhaite modeliser est une
!            fraction de la section reelle, dans la mesure ou la
!            puissance des colis modelises en 2D est une fraction de
!            la puissance totale reelle)
!         RCONV = 1 puisqu'il n'est pas possible de representer
!           geometriquement un convergent sur un maillage 2D,
!         on a, pour la cheminee d'entree (et de meme pour la cheminee
!           de sortie) :
!             SMESR = EPCHEM*XTZC/SDE
!         ceci est bien le rapport des sections passantes
!         "(representee geometriquement)/(que l'on souhaite modeliser)"
!         (les sections etant prises au dessus du convergent eventuel).

smesr = epchem*xtzc/(rconve*sde)
smssr = epchem*xtzc/(rconve*sds)


! --- Rapport pas du reseau / diametre d'un conteneur

cxsd = ptrres/dmcont
cysd = plgres/dmcont
ctsd = sqrt((cxsd*0.5d0)**2 + cysd**2)


! --- Porosite de gap entre tubes

gapscx = (cxsd - 1.d0) / cxsd

!===============================================================================
! 2. GRANDEURS PHYSIQUES
!===============================================================================

! --- Temperature et vitesse verticale en entree et sortie
!       (la vitesse servira pour le calcul des pertes de charges)

!     W moyenne dans la cheminee d'entree et volume correspondant
we0 = 0.d0
ve0 = 0.d0
!     W moyenne dans la cheminee de sortie et volume correspondant
ws0 = 0.d0
vs0 = 0.d0
!     T moyenne dans la cheminee de sortie
ts0 = 0.d0

do iel = 1, ncel
  w0 = rtpa(iel,iw(iphas))
  ifml  = ifmcel(iel   )
  icoul = iprfml(ifml,1)
!     Cheminee d'entree
  if(icoul.eq.icmtci) then
    ve0 = ve0 + volume(iel)
    we0 = we0 + volume(iel)*w0
!     Cheminee de sortie
  elseif(icoul.eq.icmtco) then
    vs0 = vs0 + volume(iel)
    ws0 = ws0 + volume(iel)*w0
    ts0 = ts0 + volume(iel)*rtpa(iel,isca(itaamt))
  endif
enddo

we0 = abs(we0)/max(ve0,epzero)
ws0 = abs(ws0)/max(vs0,epzero)
ts0 =     ts0 /max(vs0,epzero)


! --- Masse volumique entree a TINIT et sortie a TS0

roe = (trfmat+tkelvi)*rrfmat/(tinit+tkelvi)
ros = (trfmat+tkelvi)*rrfmat/(ts0  +tkelvi)


! --- Calcul des coefficients de perte de charge des cheminees

!   - Entree

!     Nombre de Reynolds (on impose 1 au minimum)
reye0 = roe*we0*smesr*dhchea/xmumat + 1.d0

!     Coefficient de perte de charge : frottement + diffuseur + filtre
cpe = 0.3164d0 * reye0**(-0.25d0) * (xlce/dhchea)                 &
     + pdccha + pdcfch


!   - Sortie

!     Nombre de Reynolds (on impose 1 au minimum)
reys0 = ros*ws0*smssr*dhches/xmumat + 1.d0

!     Coefficient de perte de charge : frottement + diffuseur + clapet
cps = 0.3164d0 * reys0**(-0.25d0) * (xlcs/dhches)                 &
     + pdcche + pdccch

!     Pertes de charge (en Pa) et affichage
if (ivar.eq.iw(iphas)) then

  pdce = cpe * 0.5d0*roe*(we0*smesr)**2
  pdcs = cps * 0.5d0*ros*(ws0*smssr)**2

  if ((irangp.le.0).and.(ntcabs.eq.ntmabs)) then
    write(impmat,5008) pdce
    write(impmat,5009) pdcs
  endif

  if ((irangp.le.0).and.(modntl.eq.0)) then
    write(nfecra,2002)                                            &
         reye0    , reys0    ,                                    &
         we0*smesr, ws0*smssr,                                    &
         we0      , ws0      ,                                    &
         smesr    , smssr    ,                                    &
         pdce     , pdcs     ,                                    &
         tinit    , ts0      ,                                    &
         roe      , ros      ,                                    &
         cpe      , cps
  endif

endif


!  --- Perte de charge reseau Diagramme 8.11 et 8.12 IdelCik
!     Definition de XXEM et XXAA

!    - Reseau de conteneurs en ligne (reseau a pas carre)

if(iconlg.eq.1) then

  xxem = -0.2d0
  rs0 = (cxsd - 1.d0)/(cysd - 1.d0)
  xxaa = nplgrs * (cxsd - 1.d0)**(-0.5d0)
  if(ptrres.le.plgres) then
    xxaa = xxaa * 1.52d0 * rs0**(-0.2d0)
  else
    xxem = xxem * rs0**(-2.d0)
    xxaa = xxaa * 0.32d0 * (rs0-0.9d0)**(-0.68d0)
  endif

!    - Reseau de conteneurs en quinconce (reseau a pas triangulaire)

else

  xxem =-0.27d0
  rs0 = (cxsd - 1.d0)/(ctsd - 1.d0)
  xxaa = (4.6d0 - 2.7d0*rs0)*(2.d0 - cxsd) + 3.2d0
  xxaa = min(xxaa,3.2d0)
!       Pour N RANGEE de tubes
  xxaa = xxaa * (nplgrs+1)

endif

!     Le Reynolds est calcule pour un ecoulement de circulation naturelle
!       pour un diametre de DH m
!       Il sert dans les pertes de charges plus bas
reyres = (roe+ros)*0.5d0 * vitref * dmcont / (gapscx * xmumat)

!     La vitesse sert a l'affichage
un0  = vitref/(gapscx**2)

!     Affichage
if( ivar.eq.iw(iphas).and.(modntl.eq.0)                           &
                     .and.(irangp.le.0) ) then
  write(nfecra,2003)                                              &
       0.5d0*rrfmat*un0*vitref*amppdc*xxaa*(reyres**xxem),        &
       xxaa*(reyres**xxem)
endif


!===============================================================================
! 3. DEFINITION DES PERTES DE CHARGE
!===============================================================================

! === Boucle sur les cellules et selection selon leur couleur
! ===========================================================

do iel = 1, ncel

  rom = propce(iel,ipcrom)

  ifml  = ifmcel(iel   )
  icoul = iprfml(ifml,1)

! === Cheminee d'alimentation (au dessus du convergent eventuel)
! ===========================================================

  if(icoul.eq.icmtci) then
    crvimp(iel) =                                                 &
         -0.5d0*volume(iel)*rom*we0*smesr*smesr/xlce*cpe

! === Cheminee d'evacuation (au dessus du convergent eventuel)
! ===========================================================

  elseif(icoul.eq.icmtco) then
    crvimp(iel) =                                                 &
         -0.5d0*volume(iel)*rom*ws0*smssr*smssr/xlcs*cps

! === Sortie de la zone de stockage
! ===========================================================

  elseif((icoul.eq.icmtro).or.                                    &
         (icoul.eq.icmtjo.and.ialveo.eq.1)) then

!     Pour representer les murs en sortie, on impose des pertes de
!       charge "infinies". Les portes de sortie doivent rester libres.

!     La position des murs dans la direction X est definie par des
!       zones dont les bornes XLIMIN et XLIMAX sont stockees dans
!       VIZCAR (avec l'indicateur ICPDCS faisant reference aux pertes de
!       charges a la sortie de la zone de stockage)
    nzones = nzocar(iligne,icpdcs)
    do izone = 1, nzones
      xlimin = ptrres*vizcar(1,izone,iligne,icpdcs)
      xlimax = ptrres*vizcar(2,izone,iligne,icpdcs)
      if ( (xyzcen(1,iel).ge.xlimin).and.                         &
           (xyzcen(1,iel).lt.xlimax) ) then
        crvimp(iel) = crvimp(iel) - 1.d5*volume(iel)
      endif
    enddo

!     La position des murs dans la direction Z est definie par des
!       zones dont les bornes ZALMIN et ZALMAX sont stockees dans
!       VIZCAR (avec l'indicateur ICPDCS faisant reference aux pertes de
!       charges a la sortie de la zone de stockage)
    nzones = nzocar(ialtit,icpdcs)
    do izone = 1, nzones
      zalmin = epchel*vizcar(1,izone,ialtit,icpdcs)
      zalmax = epchel*vizcar(2,izone,ialtit,icpdcs)
      if ( (xyzcen(3,iel).ge.zalmin).and.                         &
           (xyzcen(3,iel).lt.zalmax) ) then
        crvimp(iel) = crvimp(iel) - 1.d5*volume(iel)
      endif
    enddo


! === Reseau dans la zone de stockage
! ===========================================================

!     On traite successivement :
!       . pertes de charge de base en reseau
!       . pertes de charge verticales dans les alveoles
!       . pertes de charge horizontales pour le stockage en alveoles
!       . pertes de charge en zone libre de colis, stockage sans alveole

!     En resume :

!     - par defaut,   partout,              pour x, y, z : pdc=    reseau

!     - en alveole,   zone de colis, pour       z : pdc=    frot., bif., reu.
!     - en alveole,   plenums,              pour x, y    : pdc=pdc+frot., bif., reu.
!     - en alveole,   zone de colis, pour x, y    : pdc=    plenum inf*1.D4

!     - sans alveole, plenum sup,    pour x, y, z : pdc=    zero
!     - sans alveole, zone de colis, pour x, y, z : pdc=pdc*cartes (*0 ou *1)


  elseif(icoul.eq.icmtst) then

! --- Vitesse de reference (corrigee : vitesse en reseau)

    un0 = sqrt(rtpa(iel,iu(iphas))**2                             &
         + rtpa(iel,iv(iphas))**2                                 &
         + rtpa(iel,iw(iphas))**2)
    un0 = max(vitref,un0)
    un0 = un0/(gapscx**2)


! --- Pertes de charge de base en reseau (cas standard)

    if(ivar.eq.iw(iphas)) then
      teta0 = 0.53d0 * xxaa * reyres**xxem
    else
      teta0 = 1.00d0 * xxaa * reyres**xxem
    endif

!       Perte de charge en reseau de base
!       Attention, cette valeur est ecrasee ensuite
!         - pour le traitement des alveoles eventuelles
!         - dans les zones libres de colis

    crvimp(iel) =                                                 &
         -0.5d0*volume(iel)*rom*teta0*un0*amppdc/xlzc


! --- Pertes de charges verticales dans les alveoles

!     On traite successivement :
!       . la perte de charge de frottement
!       . la perte de charge par bifurcation branche laterale
!       . la perte de charge d'alimentation apres alveole (reunification)

    if(ialveo.eq.1) then

      if ( xyzcen(3,iel).lt.hreso.and.                            &
           xyzcen(3,iel).gt.hplen ) then

        if(ivar.eq.iw(iphas)) then


!   - Perte de charge verticale par frottement

!     Vitesse verticale homogeneisee dans l'alveole
          walv    = abs (rtpa(iel,iw(iphas)))
!     Rapport surface d'un pas / surface passante reelle d'une alveole
          rapsrf = ptrres * plgres / salv
!     Reynolds en alveole base sur une vitesse corrigee (minimum a 1)
          reyalv  = 1.d0 +                                        &
               rom*walv*rapsrf*abs(dhalve-dmcont)/xmumat
!     Hauteur d'alveole
          xlalv   = hreso - hplen
!     Coefficient de perte de charge reguliere par frottement par metre
          cpc     =                                               &
               0.3164d0*reyalv**(-0.25d0)/abs(dhalve-dmcont)

!     Perte de charge verticale par frottement
!       On ecrase la perte de charge standard en reseau
!       calculee precedemment.
          crvimp(iel) =                                           &
               - 0.5d0*volume(iel)*rom*walv*(rapsrf**2)*cpc


!     Compteurs pour impressions

          cmptaa = cmptaa + 1.d0
          pdcfaa = pdcfaa                                         &
               + 0.5d0*rom*(walv*rapsrf)**2 * cpc*xlalv
          vamcaa = vamcaa + walv * rapsrf
          vamhaa = vamhaa + walv


!   - Perte de charge verticale par bifurcation branche laterale
!     (cf Idel'cik, page 263, diagramme 7.21)

!     Rapport de surface (FPFL) et de debit (QLQP)
!       entre la branche laterale et la branche rectiligne de
!       la bifurcation (cf Idel'cik)
          fpfl = hplen * ptrres / salv
          qlqp = 1.d0 / nplgrs

          if(qlqp*fpfl.gt.0.8d0) then
            aprim = 0.9d0
          else
            aprim = 1.0d0
          endif

!     Coefficient de perte de charge singuliere par bifurcation
          cpc = aprim*(1.d0+(qlqp*fpfl)**2)/((qlqp*fpfl)**2)

!     Perte de charge verticale par bifurcation laterale
!       On l'ajoute a la perte de charge par frottement
!       calculee precedemment.
!       Comme la perte de charge est repartie sur la hauteur de l'alveole
!       il faut diviser par XLALV
          crvimp(iel) = crvimp(iel)                               &
               - 0.5d0*volume(iel)*rom*walv*(rapsrf**2)           &
               * cpc/xlalv


!     Compteurs pour impressions
          pdcbaa = pdcbaa + 0.5d0*rom*(walv*rapsrf)**2 * cpc


!   - Perte de charge verticale par alimentation apres alveole (reunification)
!     (cf Idel'cik, page 249, diagramme 7.7)

!     Rapport de surface (FPFL) et de debit (QLQP)
!       entre la branche laterale et la branche rectiligne de
!       l'alimention (reunification) (cf Idel'cik)
          fpfl = (xhzc - hreso) * ptrres / salv
          qlqp = 1.d0 / nplgrs

          acoeff = 1.d0

!     Coefficient de perte de charge singuliere par alimentation
          cpc = acoeff                                            &
               * ( 1.d0 + (qlqp*fpfl)**2 - 2.d0*(1.d0-qlqp) )     &
               / ( (qlqp*fpfl)**2 )


!     Perte de charge verticale par alimentation (reunification)
!       On l'ajoute a la perte de charge par frottement et a la
!       perte de charge singuliere calculee precedemment.
!       Comme la perte de charge est repartie sur la hauteur de l'alveole
!       il faut diviser par XLALV
          crvimp(iel) = crvimp(iel)                               &
               - 0.5d0*volume(iel)*rom*walv*(rapsrf**2)           &
               * cpc/xlalv

!     Compteurs pour impressions
          pdcraa = pdcraa + 0.5d0*rom*(walv*rapsrf)**2 * cpc


!     Fin du traitement des pertes de charges verticales en alveole
        endif
      endif
    endif


! --- Pertes de charges horizontales pour le stockage en alveoles

!     On traite successivement :
!       . la perte de charge de reunion (au dessus des colis) et
!                            de bifurcation (au dessous des colis)
!       . la perte de charge de frottement
!                             au dessus et au dessous des colis
!       . la perte de charge dans la zone de colis

    if(ialveo.eq.1) then

      if(ivar.eq.iu(iphas).or.ivar.eq.iv(iphas)) then


!   - Perte de charge horizontale de reunion et de bifurcation

!     La perte de charge calculee dans la zone des colis sera modifiee
!       plus bas (pour simplifier les tests, les formules valables
!       au dessus et au dessous des colis sont egalement appliquees
!       dans la zone des colis).

!     Rapport de debit (QLQP) entre la branche laterale et
!       la branche rectiligne (cf Idel'cik)
        qlqp = 1.d0 / nplgrs

        if(xyzcen(3,iel).gt.hreso) then

!     Coefficient de perte de charge horizontale de reunion
!       pour le stockage en alveoles (apres alveole, branche rectiligne)
!       (cf Idel'cik, page 249, diagramme 7.7)

          cpc = (1.55d0*qlqp-qlqp**2) / (1.d0-qlqp)**2

        else

!     Coefficient de perte de charge horizontale de bifurcation
!       pour le stockage en alveoles (avant alveole, branche rectiligne)
!       (cf Idel'cik, page 265, diagramme 7.23)

!     Ici, on traite la zone XYZCEN(3,IEL).LE.HRESO et donc non
!       seulement le dessous des colis mais aussi la zone des colis.
!       Neanmoins, la perte de charge dans la zone des colis sera
!       modifiee plus bas.

          cpc = 0.4d0 * qlqp**2 / (1.d0-qlqp)**2

        endif


!     Perte de charge horizontale de reunion ou de bifurcation
!       On l'ajoute a CRVIMP
!       Comme la perte de charge est repartie sur la hauteur de l'alveole
!         il faut diviser par XLZC
!       On utilise la vitesse de gap UN0 (les puits sont supposes traverser
!         les plenums superieur et inferieur)
        crvimp(iel) = crvimp(iel)                                 &
             - 0.5d0*volume(iel)*rom*un0 * cpc/xlzc


!   - Perte de charge horizontale de frottement
!       au dessus et au dessous des colis

!     La perte de charge calculee dans la zone des colis sera modifiee
!       plus bas (pour simplifier les tests, les formules valables
!       au dessus et au dessous des colis sont egalement appliques
!       dans la zone des colis).

!     Vitesse homogeneisee locale
        unh = sqrt(rtpa(iel,iu(iphas))**2                         &
             + rtpa(iel,iv(iphas))**2                             &
             + rtpa(iel,iw(iphas))**2)

!     Reynolds au dessus et au dessous des colis (minimum 1)
!       seul le diametre hydraulique considere differe
        if(xyzcen(3,iel).gt.hreso) then
          dhplen = 2.d0*(xhzc - hreso)
        else
          dhplen = 2.d0*hplen
        endif
        reyhor = rom*unh*dhplen/xmumat + 1.d0

!     Coefficient de perte de charge reguliere par frottement par metre
        cpc =                                                     &
             0.3164d0 * reyhor**(-0.25d0) / dhplen

!     Perte de charge horizontale par frottement
!       On l'ajoute a la perte de charge par reunion et bifurcation
!       On ne doit pas diviser par XLZC car CPC est deja calcule par
!         metre : on n'a pas multiplie par XLZC dans CPC
        crvimp(iel) = crvimp(iel)                                 &
             - 0.5d0*volume(iel)*rom*un0 * cpc


!   - Perte de charge horizontale dans la zone des colis
!     On augmente la valeur de perte de charge calculee
!       precedemment, pour guider l'ecoulement sur la verticale

        if ( xyzcen(3,iel).lt.hreso.and.                          &
             xyzcen(3,iel).gt.hplen) then
          crvimp(iel) = crvimp(iel) * 1.d4
        endif

!     Fin du traitement des pertes de charges horizontales en alveole
      endif
    endif


! --- Perte de charge en zone libre de colis en stockage sans alveole
!     Dans les zones libres, on annule la perte de charge

    if(ialveo.ne.1) then

!     Au dessus des colis
      if (xyzcen(3,iel).ge.hreso) then
        crvimp(iel) = 0.d0
      endif

!     Les entrepots ne sont pas uniformement remplis.
!     Pour représenter la carte d'encombrement des rangees de colis, on
!       repere les rangees occupees (coordonnee Y, indicateur IRANGE)
!       par leurs bornes YZOMIN et YZOMAX qui sont stockees dans VIZCAR
!       (avec l'indicateur ICPDCR faisant reference aux pertes de
!       charge en reseau)

      yrgmin = 0.d0
      yrgmax = nplgrs*plgres

      if ( (xyzcen(2,iel).ge.yrgmin).and.                         &
           (xyzcen(2,iel).lt.yrgmax) ) then

        indixy = 0
        nzones = nzocar(irange,icpdcr)
        do izone = 1, nzones
          yzomin = plgres*vizcar(1,izone,irange,icpdcr)
          yzomax = plgres*vizcar(2,izone,irange,icpdcr)
          if ( (xyzcen(2,iel).ge.yzomin).and.                     &
               (xyzcen(2,iel).lt.yzomax) ) then
            indixy = 1
          endif
        enddo

        if (indixy.eq.0) then
          crvimp(iel) = 0.d0
        endif

      endif

!     Les entrepots ne sont pas uniformement remplis.
!     Pour représenter la carte d'encombrement des lignes de colis, on
!       repere les lignes occupees (coordonnee X, indicateur ILIGNE)
!       par leurs bornes XZOMIN et XZOMAX qui sont stockees dans VIZCAR
!       (avec l'indicateur ICPDCR faisant reference aux pertes de
!       charge en reseau)

      xlimin = 0.d0
      xlimax = nptran*ptrres

      if ( (xyzcen(1,iel).ge.xlimin).and.                         &
           (xyzcen(1,iel).lt.xlimax) ) then
        indixy = 0
        nzones = nzocar(iligne,icpdcr)
        do izone = 1, nzones
          xzomin = ptrres*vizcar(1,izone,iligne,icpdcr)
          xzomax = ptrres*vizcar(2,izone,iligne,icpdcr)
          if ( (xyzcen(1,iel).ge.xzomin).and.                     &
               (xyzcen(1,iel).lt.xzomax) ) then
            indixy = 1
          endif
        enddo

        if (indixy.eq.0) then
          crvimp(iel) = 0.d0
        endif

      endif

!     Fin du traitement des pertes de charges en zone libre sans alveole
    endif


! === Entree de la zone de stockage
! ===========================================================

  elseif((icoul.eq.icmtri).or.                                    &
         (icoul.eq.icmtji.and.ialveo.eq.1)) then

!     Pour représenter les murs en entree, on impose des pertes de
!       charge "infinies". Les portes d'entree doivent rester libres.

!     La position des murs dans la direction X est definie par des
!       zones dont les bornes XLIMIN et XLIMAX sont stockees dans
!       VIZCAR (avec l'indicateur ICPDCE faisant reference a la sortie)
    nzones = nzocar(iligne,icpdce)
    do izone = 1, nzones
      xlimin = ptrres*vizcar(1,izone,iligne,icpdce)
      xlimax = ptrres*vizcar(2,izone,iligne,icpdce)
      if ( (xyzcen(1,iel).ge.xlimin).and.                         &
           (xyzcen(1,iel).lt.xlimax) ) then
        crvimp(iel) = crvimp(iel) - 1.d5*volume(iel)
      endif
    enddo

!     La position des murs dans la direction Z est definie par des
!       zones dont les bornes ZALMIN et ZALMAX sont stockees dans
!       VIZCAR (avec l'indicateur ICPDCE faisant reference a la sortie)
    nzones = nzocar(ialtit,icpdce)
    do izone = 1, nzones
      zalmin = epchel*vizcar(1,izone,ialtit,icpdce)
      zalmax = epchel*vizcar(2,izone,ialtit,icpdce)
      if ( (xyzcen(3,iel).ge.zalmin).and.                         &
           (xyzcen(3,iel).lt.zalmax) ) then
        crvimp(iel) = crvimp(iel) - 1.d5*volume(iel)
      endif
    enddo

!     Fin du test sur les couleurs et de la boucle sur les elements
  endif
enddo

!===============================================================================
! 4. CALCULS ET AFFICHAGES FINALS
!===============================================================================

! --- Moyennes relatives au cas alveoles

if(cmptaa.gt.0.d0) then

  vamcaa = vamcaa/cmptaa
  pdcfaa = pdcfaa/cmptaa
  pdcbaa = pdcbaa/cmptaa
  pdcraa = pdcraa/cmptaa
  vamhaa = vamhaa/cmptaa

  if ((modntl.eq.0).and.(irangp.le.0)) then
    write(nfecra,2004) vamcaa, vamhaa,                            &
         pdcfaa, pdcbaa, pdcraa
  endif
endif


! --- En convection naturelle,
!       . actualisation et impression de la vitesse de reference
!         (utilisee dans mtphyv en particulier)
!       . calcul et impression du Richardson

if(icofor.eq.0.and.ivar.eq.iw(iphas))then

!   - Vitesse de reference

!     Calcul de la moyenne des vitesses homogenes moyennes
!       dans les cheminees d'entree et de sortie :
!       0.5D0*(WE0+WS0)
!     Cette vitesse est relative a une section de cheminee
!       sur le maillage de EPCHEM*XTZC/RCONVE (section dans
!       la partie haute, au dessus du convergent eventuel) :
!       le debit associe, represente sur le maillage, est donc
!       Q=0.5D0*(WE0+WS0)*EPCHEM*XTZC/RCONVE
!     Ce debit traverse ensuite la zone de stockage de section
!       passante XTZC*XHZC. La vitesse qui s'en deduit dans
!       cette region est donc Q/(XTZC*XHZC), soit
!       0.5D0*(WE0+WS0)*EPCHEM*XTZC/(RCONVE*XTZC*XHZC)
!     On la relaxe a partir de l'ancienne valeur de VITREF

  vitref = ( 0.5d0*(we0+ws0)*epchem/(rconve*xhzc)                 &
       + vitref )*0.5d0

!     Impressions
  if ((modntl.eq.0).and.(irangp.le.0)) then
    write(nfecra,2005) vitref, viscl0(iphas)
  endif


!   - Richardson

!     Calcul selon la formule Ri = g beta DeltaT H_ref / U**2 avec :
!       . beta = 1/((TS0+TINIT)*0.5D0 + TKELVI)
!       . DeltaT = TS0-TINIT (temperature sortie-entree)
!       . H_ref = XHZC (hauteur zone de stockage)
!       . U = VITREF/GAPSCX (U gap)

  gmat   = sqrt(gx**2+gy**2+gz**2)
  betmat = 1.d0/((ts0+tinit)*0.5d0 + tkelvi)
  ureel  = vitref/gapscx

  richar = gmat*betmat*(ts0-tinit)*xhzc/(ureel**2)

!     Impressions
  if ((modntl.eq.0).and.(irangp.le.0)) then
    write(nfecra,2006) ureel, richar
  endif
  if ((irangp.le.0).and.(ntcabs.eq.ntmabs)) then
    write(impmat,5011) richar
  endif

endif


!--------
! FORMATS
!--------

 1000 format(' TERMES SOURCES MATISSE POUR LA VARIABLE ',A8,/)
 2001 FORMAT (/,3X,'** INFORMATIONS SUR MATISSE, VARIABLE ',A8,/, &
          3X,'   -------------------------------------------')
 2002 format (                                                          &
'mati --------------------------------------------------------',/,&
'mati Cheminees (au dessus des convergents eventuels)         ',/,&
'mati ---------------------------   entree      sortie   unite',/,&
'mati --------------------------------------------------------',/,&
'mati Nombre de Reynolds        ', E12.5   , E12.5    ,'     -',/,&
'mati Vitesse reelle            ', E12.5   , E12.5    ,'   m/s',/,&
'mati Vitesse homogeneisee      ', E12.5   , E12.5    ,'   m/s',/,&
'mati Sct° maillage/a representer', E11.5  , E12.5    ,'     -',/,&
'mati Perte de charge           ', E12.5   , E12.5    ,'    Pa',/,&
'mati T° (TINIT, moyenne sortie)', E12.5   , E12.5    ,'    °C',/,&
'mati Masse volumique associee  ', E12.5   , E12.5    ,' kg/m3',/,&
'mati Coeff. de perte de charge ', E12.5   , E12.5    ,'     -',/,&
'mati --------------------------------------------------------',/)
 2003 format (                                                          &
'mati --------------------------------------------------------',/,&
'mati Reseau                                                  ',/,&
'mati ------                                                  ',/,&
'mati Perte de charge reseau    (/0.53)     ', E12.5  ,'    Pa',/,&
'mati Coeff. de perte de charge (/0.53)     ', E12.5  ,'     -',/,&
'mati --------------------------------------------------------',/)
 2004 format (                                                          &
'mati --------------------------------------------------------',/,&
'mati Alveoles                                                ',/,&
'mati --------                                                ',/,&
'mati Vitesse moyenne reelle                ', E12.5  ,'   m/s',/,&
'mati Vitesse moyenne homogeneisee          ', E12.5  ,'   m/s',/,&
'mati Pdc reguliere  moyenne de frottement  ', E12.5  ,'    Pa',/,&
'mati Pdc singuliere moyenne de bifurcation ', E12.5  ,'    Pa',/,&
'mati Pdc singuliere moyenne de reunification', E11.5 ,'    Pa',/,&
'mati --------------------------------------------------------',/)
 2005 format (                                                          &
'mati --------------------------------------------------------',/,&
'mati Vitesse de reference (zone de stockage)                 ',/,&
'mati --------------------                                    ',/,&
'mati Vitesse debitante homogeneisee       ', E12.5  ,'    m/s',/,&
'mati Viscosite dynamique totale           ', E12.5  ,' kg/m/s',/,&
'mati --------------------------------------------------------',/)
 2006 format (                                                          &
'mati --------------------------------------------------------',/,&
'mati Richardson (zone de stockage)                           ',/,&
'mati ----------                                              ',/,&
'mati Vitesse reelle de reference           ', E12.5  ,'   m/s',/,&
'mati Richardson base sur la vitesse reelle ', E12.5  ,'     -',/,&
'mati --------------------------------------------------------',/)

 5008 format(' Perte de charge de cheminee d''entree                  ',&
 ' :', E12.5, ' Pa')
 5009 format(' Perte de charge de cheminee de sortie                  ',&
 ':', E12.5, ' Pa')
 5011 format(' Nombre de Richardson                                   ',&
 ':', E12.5)


!----
! FIN
!----

return

end subroutine
