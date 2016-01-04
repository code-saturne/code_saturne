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

subroutine lagopt
!================

!===============================================================================
!  FONCTION  :
!  ---------

!   SOUS-PROGRAMME DU MODULE LAGRANGIEN :
!   -------------------------------------

!  1) Initialisation par defaut du parametrage du module
!     lagrangien

!  2) Lecture du parametrage utilisateur

!  3) Verifications du parametrage utilisateur et
!     controles de coherence

!  4) Initialisation des variables en COMMON et des pointeurs
!     sur les tableaux lies aux particules, aux statistiques,
!     aux conditions aux limites, aux variables parietales,
!     aux donnees pour le couplage retour.

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
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
use parall
use dimens
use numvar
use entsor
use optcal
use cstphy
use lagpar
use lagdim
use lagran
use ppppar
use ppthch
use ppincl
use cpincl
use radiat
use ihmpre

!===============================================================================

implicit none

integer  ii , ip , irf , icha , i1 , i2 , i3, iok

!===============================================================================
! 0. INITIALISATION
!===============================================================================


!===============================================================================
! 1. INITIALISATIONS PAR DEFAUT DU MODULE LAGRANGIEN
!                        ^^^^^^
!    TOUTES CES INITIALISATIONS PRENNENT DES VALEURS NON ADMISSIBLES
!    A LA VERIFICATION POUR ETRE SUR QUE L'UTILISATEUR LES A MODIFIEES
!===============================================================================

!     IILAGR = 0 : PAS DE CALCUL LAGRANGIEN
!            = 1 : DIPHASIQUE LAGRANGIEN SANS COUPLAGE RETOUR
!            = 2 : DIPHASIQUE LAGRANGIEN AVEC COUPLAGE RETOUR
!            = 3 : DIPHASIQUE LAGRANGIEN SUR CHAMP FIGE

iilagr = 0

!     ISUILA = 0 : PAS DE SUITE LAGRANGIENNE
!            = 1 : SUITE LAGRANGIENNE

isuila = 0

!     Suite statistiques et TS couplage retour si =1

isuist = 0

!     IPHYLA = 0 : PUREMENT DYNAMIQUE
!            = 1 : EQUATION SUR TP, DP et MP
!            = 2 : CHARBON

iphyla = 0

!     EQUATION SUR LE DIAMETRE (O : NON , 1 : OUI)

idpvar = 0

!     EQUATION SUR LA TEMPERATURE (O : NON , 1 : OUI)

itpvar = 0

!     EQUATION SUR LA MASSE (O : NON , 1 : OUI)

impvar = 0

!     TEMPERATURE D'INITIALISATION DES PARTICULES DEJA PRESENTES

tpart = -999.d0

!     CHALEUR MASSIQUE D'INITIALISATION DES PARTICULES DEJA PRESENTES

cppart = -999.d0

!     ENCRASSEMENT
!     IENCRA = 0 pas d'encrassement
!            = 1 encrassement
!     TPRENC : TEMPERATURE MINIMUM POUR l'ENCRASSEMENT
!     VISREF : VISCOSITE CRITIQUE

iencra = 0

do icha = 1 , ncharm2
  tprenc(icha) = -999.d0
  visref(icha) = -999.d0
enddo

!     NOMBRE DE PARTICULES MAXIMAL AUTORISE DANS LE DOMAINE

nbpmax = 1000

!     NOMBRE DE VARIABLES SUPPLEMENTAIRE SUR LES PARTICULES

nvls = 0

!     CARACTERE STATIONNAIRE DE L'ECOULEMENT DE LA PHASE CONTINUE

isttio = 0

!     Nombre de passages absolus (i.e. suite comprise) avant de faire
!     une moyenne en temps (CALCUL STATIONNAIRE) des termes sources
!     de couplage retour (seulement si ISTTIO = 1)

nstits = 1

!     COUPLAGE RETOUR SUR LA DYNAMIQUE (Vitesse + Turbulence)

ltsdyn = 0

!     COUPLAGE UTILE SUR LA MASSE

ltsmas = 0

!     COUPLAGE UTILE SUR LA THERMIQUE OU LES VARIABLES CHARBON

ltsthe = 0

!     STATISTIQUES

!     Calcul statistiques si  = 1 ; si 0 pas de stat ; si -1 stop

istala = 0

!     Nombre de variables statistiques supplementaires

nvlsts = 0

!     Numero du pas de temps pour debut statistiques

idstnt = 1

!     Debut calcul stationnaire

nstist = 1

!     Seuil en POIDS STAT de particules pour les stats

seuil = 0.d0

!     Nom des variables

do ii = 1,nvplmx

  write(nomlag(ii),'(A6,I4.4)') 'MoyLag',ii

  write(nomlav(ii),'(A6,I4.4)') 'VarLag',ii

  write(nombrd(ii),'(A6,I4.4)') 'BrdLag',ii

enddo

!     Historique pour la moyenne et la variance

do ii = 1,nvplmx
  ihslag(ii)  = 0
enddo

!     INJECTION CONTINUE

injcon = 0

!     ROULETTE RUSSE

iroule = 0

!     ORDRE D'INTEGRATION

nordre = 1

!     DISPERSION TURBULENTE

idistu = 1

!     DIFFUSION TURBULENTE

idiffl = 0

!     MODEL COMPLET

modcpl = 0

!     DIRECTION ASSOCIEE AU MODEL COMPLET

idirla = 0

!     CORRECTION DE PRESSION (EXPERIMENTAL)

ilapoi = 0

!     POSTPROCESSING MODE : TRAJECTOIRES

iensi1 = 0

!     POSTPROCESSING MODE : DEPLACEMENTS

iensi2 = 0

!     NOMBRE DE PARTICULES A VISUALISER MAXIMUM=NLISTE

nbvis  = 20

!     FREQUENCE D'AQUISITION DES DONNEES A VISUALISER

nvisla = 1

!     INITIALISATION  : PAR DEFAUT ON NE VISUALISE AUCUNE PARTICULE

do ii = 1, nliste
  liste(ii) = -1
  list0(ii) = -1
  nplist(ii) = 0
enddo

!     POSTPROCESSING VARIABLE : VITESSE DU FLUIDE VU

ivisv1  = 0

!     POSTPROCESSING VARIABLE : VITESSE DE LA PARTICULE

ivisv2  = 0

!     POSTPROCESSING VARIABLE : TEMPS DE SEJOUR

ivistp  = 0

!     POSTPROCESSING VARIABLE : DIAMETRE

ivisdm  = 0

!     POSTPROCESSING VARIABLE : MASSE

ivismp  = 0

!     POSTPROCESSING VARIABLE : TEMPERATURE

iviste  = 0

!     POSTPROCESSING VARIABLE : TEMPERATURE CHARBON

ivishp  = 0

!     POSTPROCESSING VARIABLE : DIAMETRE DU COEUR RETRECISSANT

ivisdk  = 0

!     POSTPROCESSING VARIABLE : MASSE CHARBON REACTIF

ivisch  = 0

!     POSTPROCESSING VARIABLE : MASSE COKE

ivisck  = 0

!     POSTPROCESSING MODE : INTERACTIONS PARTICULES/FRONTIERES

iensi3 = 0

!     POSTPROCESSING : STAT PARIETALES STATIONNAIRES

nstbor = 1

!     Seuil en POIDS STAT de particules pour les stats

seuilf = 0.d0

!     INFORMATIONS A ENREGISTRER

inbrbd = 0
iflmbd = 0
iangbd = 0
ivitbd = 0
iencbd = 0

nusbor = 0

!     Type de moyenne applicable pour affichage et post-processing

do ii = 1,nusbrd+10
  imoybr(ii) = 0
enddo

!===============================================================================
! 2.1 INITIALISATIONS UTILISATEUR DU MODULE LAGRANGIEN
!                     ^^^^^^^^^^^
!===============================================================================

if (iihmpr.eq.1) then

  do ii = 1, nvplmx
    call fclag1(nomlag(ii), len(nomlag(ii)), ii)
    call fclag2(nomlav(ii), len(nomlav(ii)), ii)
    call fclag3(nombrd(ii), len(nombrd(ii)), ii)
  enddo

  call uilag1                                                      &
  !==========
 ( iilagr, isuila, isuist, nbpmax, isttio, injcon,                 &
   iphyla, idpvar, itpvar, impvar,                                 &
   iencra, tprenc, visref, enc1, enc2,                             &
   nstits, ltsdyn, ltsmas, ltsthe,                                 &
   nordre, idistu, idiffl, modcpl, idirla,                         &
   iensi1, iensi2, ntlal,  nbvis, nvisla,                          &
   ivisv1, ivisv2, ivistp, ivisdm, iviste,                         &
   ivismp, ivishp, ivisdk, ivisch, ivisck,                         &
   istala, nbclst, seuil, idstnt,  nstist,                         &
   ihslag, iensi3, seuilf, nstbor,                                 &
   inbrbd, iflmbd, iangbd, ivitbd, iencbd, imoybr,                 &
   iactfv, iactvx, iactvy, iactvz, iactts)

  do ii = 1, nvplmx
    call cfname(1, nomlag(ii), len(nomlag(ii)), ii)
    call cfname(2, nomlav(ii), len(nomlav(ii)), ii)
    call cfname(3, nombrd(ii), len(nombrd(ii)), ii)
  enddo

  do ii = 1, nbvis
    liste(ii) = ii
  enddo

endif

call uslag1
!==========

if (iilagr.eq.0) return

!===============================================================================
! 2.2 VERIFICATION DES INITIALISATIONS UTILISATEUR DU MODULE LAGRANGIEN
!===============================================================================

! on doit verifier toutes les options entrees par l'utilisateur
!  qui est inventif
!  et si une valeur ne correspond pas, on ne corrige pas, on s'arrete.

iok = 0

!     IILAGR

if (iilagr.lt.0 .or. iilagr.gt.3) then
  write(nfecra,1010) iilagr
  iok = iok + 1
endif

!    CALCUL SUR CHAMP FIGE : SUITE OBLIGATOIRE
!    ATTENTION : LE CHAMP FIGE LAGRANGIEN N"EST PAS LE MEME QUE CELUI
!    DE L'EULERIEN. POUR LE LAGRANGIEN VITESSE PRESSION ET SCALAIRES
!    SONT CONSTANTS.

if (iilagr.eq.3 .and. isuite.ne.1) then
  write(nfecra,1012) iilagr, isuite
  iok = iok + 1
endif

if (iilagr.eq.3) iccvfg = 1

if (iilagr.ne.2 .and. ippmod(icpl3c).ge.1) then
  write(nfecra,1013) iilagr, ippmod(icpl3c)
  iok = iok + 1
endif

if (iilagr.gt.0 .and. idtvar.eq.2) then
  write(nfecra,1014) iilagr, idtvar
  iok = iok + 1
endif

!     ISUILA ISUIST

if (isuila.lt.0 .or. isuila.gt.1) then
  write(nfecra,1020) isuila
  iok = iok + 1
endif

if (isuila.eq.1 .and. isuite.eq.0) then
  write(nfecra,1021)
  iok = iok + 1
endif

if (isuila.eq.1) then
  if (isuist.lt.0 .or. isuist.gt.1) then
    write(nfecra,1022) isuist
    iok = iok + 1
  endif
else
  isuist = 0
endif

!     IPHYLA

if (iphyla.lt.0 .or. iphyla.gt.2) then
  write(nfecra,1030) iphyla
  iok = iok + 1
endif

if (iok.ne.0) call csexit (1)
              !==========

!     IDPVAR ITPVAR IMPVAR

!     Couplage-retour uniquement vers la phase continue

if (iphyla.eq.1) then
  if (idpvar.lt.0 .or. idpvar.gt.1) then
    write(nfecra,1031) idpvar
    iok = iok + 1
  endif
  if (itpvar.lt.0 .or. itpvar.gt.1) then
    write(nfecra,1032) itpvar
    iok = iok + 1
  endif
  if (impvar.lt.0 .or. impvar.gt.1) then
    write(nfecra,1033) impvar
    iok = iok + 1
  endif
  if (itpvar.eq.1 .and. iscalt.eq.-1) then
    write(nfecra,1034) itpvar, iscalt
    iok = iok + 1
  endif
else
  itpvar = 0
  impvar = 0
  idpvar = 0
endif

if (isuila.eq.1 .and. iphyla.eq.1 .and. itpvar.eq.1) then
  if (cppart.lt.0) then
    write(nfecra,1036) cppart
    iok = iok + 1
  endif
  if (tpart.lt.tkelvn) then
    write(nfecra,1037) tkelvn, tpart
    iok = iok + 1
  endif
endif

if (iok.ne.0) call csexit (1)
              !==========

!     IENCRA TPRENC VISREF

if (iphyla.eq.2) then
  if (iencra.lt.0 .or. iencra.gt.1) then
    write(nfecra,1040) iencra
    iok = iok + 1
  endif

  do icha = 1 , ncharb
    if (iencra.eq.1 .and. visref(icha).lt.0) then
      write(nfecra,1041) iencra, visref(icha), icha
      iok = iok + 1
    endif
    if (iencra.eq.1 .and. tprenc(icha).lt.tkelvn) then
      write(nfecra,1042) iencra, tkelvn, tprenc(icha), icha
      iok = iok + 1
    endif
  enddo

else
  iencra = 0
endif

if (iphyla.ne.2 .and. ippmod(icpl3c).ge.0) then
  write(nfecra,1043) iphyla, ippmod(icpl3c)
  iok = iok + 1
endif

if ( iphyla.eq.2 .and. (ippmod(icpl3c).lt.0 .and.                 &
                        ippmod(icp3pl).lt.0) ) then
  write(nfecra,1044) iphyla, ippmod(icpl3c), ippmod(icp3pl)
  iok = iok + 1
endif

if (iok.ne.0) call csexit (1)
              !==========

!     NBPMAX

if (nbpmax.lt.0) then
  write(nfecra,1050) nbpmax
  iok = iok + 1
endif

!     NVLS

if (nvls.lt.0 .or. nvls.gt.nusvar) then
  write(nfecra,1060) nusvar, nvls
  iok = iok + 1
endif

if (iok.ne.0) call csexit (1)
              !==========

!     ISTTIO NSTITS LTSDYN LTSMAS LTSTHE

!     Si champs figes alors forcement en stationnaire
if (iilagr.eq.3) isttio = 1

if (isttio.lt.0 .or. isttio.gt.1) then
  write(nfecra,1061) isttio
  iok = iok + 1
endif

if (iilagr.eq.2) then
  if (isttio.eq.1 .and. nstits.lt.1) then
    write(nfecra,1062) nstits
    iok = iok + 1
  endif
  if (ltsdyn.lt.0 .or. ltsdyn.gt.1) then
    write(nfecra,1063) ltsdyn
    iok = iok + 1
  endif
  if (iphyla.eq.1 .and. (impvar.eq.1 .or. idpvar.eq.1)) then
    if (ltsmas.lt.0 .or. ltsmas.gt.1) then
      write(nfecra,1064) ltsmas
      iok = iok + 1
    endif
  else
    ltsmas = 0
  endif
  if ((iphyla.eq.1 .and. itpvar.eq.1) .or. iphyla.eq.2) then
    if (ltsthe.lt.0 .or. ltsthe.gt.1) then
      write(nfecra,1065) ltsthe
      iok = iok + 1
    endif
  else
    ltsthe = 0
  endif
  if (ltsdyn.eq.1 .and. iccvfg.eq.1) then
    write(nfecra,1066) ltsdyn, iccvfg
    iok = iok + 1
  endif
  if (ltsdyn.ne.1 .and. ltsthe.ne.1 .and. ltsmas.ne.1) then
    write(nfecra,1067) iilagr, ltsdyn, ltsthe, ltsmas
    iok = iok + 1
  endif
else
  ltsdyn = 0
  ltsmas = 0
  ltsthe = 0
endif

if (iok.ne.0) call csexit (1)
              !==========

!     ISTALA SEUIL IDSTNT NSTIST NVLSTS

if (istala.lt.0 .or. istala.gt.1) then
  write(nfecra,1070) istala
  iok = iok + 1
endif

if (istala.eq.1) then
  if (seuil.lt.0.d0) then
    write(nfecra,1071) seuil
    iok = iok + 1
  endif
  if (idstnt.lt.1) then
    write(nfecra,1072) idstnt
    iok = iok + 1
  endif
  if (isttio.eq.1) then
    if (nstist.lt.idstnt) then
      write(nfecra,1073) idstnt, nstist
      iok = iok + 1
    endif
  endif
  if (nvlsts.lt.0 .or. nvlsts.gt.nussta) then
    write(nfecra,1074) nussta, nvlsts
    iok = iok + 1
  endif
  if ( nbclst .lt. 0 .or. nbclst.gt. nclstm ) then
    write(nfecra,1075) nclstm,nbclst
    iok = iok + 1
  endif

endif

if (iok.ne.0) call csexit (1)
              !==========

if (istala.eq.0) then
  isuist = 0
  seuil  = 0.d0
  idstnt = 0
  nstist = 0
endif

!     INJCON

if (injcon.lt.0 .or. injcon.gt.1) then
  write(nfecra,1080) injcon
  iok = iok + 1
endif

!     IROULE

if (iroule.lt.0 .or. iroule.gt.2) then
  write(nfecra,1090) iroule
  iok = iok + 1
endif

!     NORDRE

if (nordre.ne.1) then
  write(nfecra,2000) nordre
  iok = iok + 1
endif

!     IDISTU

if (idistu.lt.0 .or. idistu.gt.1) then
  write(nfecra,2010) idistu
  iok = iok + 1
endif

if (idistu.eq.1 .and. itytur.ne.2 .and. itytur.ne.3       &
     .and.  iturb.ne.50 .and. iturb.ne.60 ) then
  write(nfecra,2011) iilagr, idistu, iturb
  iok = iok + 1
else if (idistu.eq.0 .and. iturb.ne.0 .and.                   &
         itytur.ne.2 .and. itytur.ne.3                    &
     .and.  iturb.ne.50 .and. iturb.ne.60) then
  write(nfecra,2012) iilagr, idistu, iturb
  iok = iok + 1
endif

!     IDISTU

if (idiffl.lt.0 .or. idiffl.gt.1) then
  write(nfecra,2013) idiffl
  iok = iok + 1
endif

!     MODCPL IDIRLA

if (modcpl.lt.0) then
  write(nfecra,2014) modcpl
  iok = iok + 1
endif
if (modcpl.gt.0) then
  if (modcpl.lt.idstnt) then
    write(nfecra,2015) modcpl, idstnt
    iok = iok + 1
  endif
  if (istala.eq.0) then
    write(nfecra,2018) modcpl, istala
    iok = iok + 1
  endif
  if (idirla.ne.1 .and. idirla.ne.2 .and. idirla.ne.3) then
    write(nfecra,2016) idirla
    iok = iok + 1
  endif
endif

!     ILAPOI

if (ilapoi.lt.0 .or. ilapoi.gt.1) then
  write(nfecra,2017) ilapoi
  iok = iok + 1
endif

if (iok.ne.0) call csexit (1)
              !==========

!     IENSI1 IENSI2

if (iensi1.lt.0 .or. iensi1.gt.1) then
  write(nfecra,2030) iensi1
  iok = iok + 1
endif
if (iensi2.lt.0 .or. iensi2.gt.1) then
  write(nfecra,2031) iensi2
  iok = iok + 1
endif

if (irangp.ge.0) then
   if (iensi1.gt.0 .or. iensi2.gt.0) then
      write(nfecra,3015)
      iok = iok + 1
   endif
endif

!     NBVIS NVISLA

if (iensi1.eq.1 .or. iensi2.eq.1) then
  if (nbvis.gt.nbpmax .or. nbvis.gt.nliste .or.  nbvis.lt.0) then
    write(nfecra,2032) nbpmax, nliste, nbvis
    iok = iok + 1
  endif
  if (nvisla.le.0 ) then
    write(nfecra,2033) nvisla
    iok = iok + 1
  endif

!     IVISV1 IVISV2 IVISTP IVISDM IVISTE

  if (ivisv1.lt.0 .or. ivisv1.gt.1) then
    write(nfecra,2040) ivisv1
    iok = iok + 1
  endif
  if (ivisv2.lt.0 .or. ivisv2.gt.1) then
    write(nfecra,2041) ivisv2
    iok = iok + 1
  endif
  if (ivistp.lt.0 .or. ivistp.gt.1) then
    write(nfecra,2042) ivistp
    iok = iok + 1
  endif
  if (ivisdm.lt.0 .or. ivisdm.gt.1) then
    write(nfecra,2043) ivisdm
    iok = iok + 1
  endif
  if (iphyla.eq.1 .and. itpvar.eq.1) then
    if (iviste.lt.0 .or. iviste.gt.1) then
      write(nfecra,2044) iviste
      iok = iok + 1
     endif
  else
    iviste = 0
  endif

!       IVISHP IVISDK IVISCH IVISCK

  if (iphyla.eq.2) then
    if (ivishp.lt.0 .or. ivishp.gt.1) then
      write(nfecra,2045) ivishp
      iok = iok + 1
    endif
    if (ivisdk.lt.0 .or. ivisdk.gt.1) then
      write(nfecra,2046) ivisdk
      iok = iok + 1
    endif
    if (ivisch.lt.0 .or. ivisch.gt.1) then
      write(nfecra,2047) ivisch
      iok = iok + 1
    endif
    if (ivisck.lt.0 .or. ivisck.gt.1) then
      write(nfecra,2048) ivisck
      iok = iok + 1
    endif
  else
    ivishp = 0
    ivisdk = 0
    ivisch = 0
    ivisck = 0
  endif

endif

!     IENSI3 NSTBOR

if (iensi3.lt.0 .or. iensi3.gt.1) then
  write(nfecra,2050) iensi3
  iok = iok + 1
endif
if (iensi3.eq.1 .and. isttio.eq.1) then
  if (nstbor.lt.1) then
    write(nfecra,2057) nstbor
    iok = iok + 1
  endif
else
  nstbor = 1
endif
if (iensi3.eq.1) then
  if (seuilf.lt.0.d0) then
    write(nfecra,2058) seuilf
    iok = iok + 1
  endif
else
  seuilf = 0.d0
endif

!     INBRBD IFLMBD IANGBD IVITBD IENCBD NUSBOR

if (iensi3.eq.1) then

  if (iphyla.eq.2 .and. iencra.eq.1) then
    if (iencbd.lt.0 .or. iencbd.gt.1) then
      write(nfecra,2051)
      iok = iok + 1
    endif
  else
    iencbd = 0
  endif

  if (inbrbd.lt.0 .or. inbrbd.gt.1) then
    write(nfecra,2052) inbrbd
    iok = iok + 1
  endif
  if (iflmbd.lt.0 .or. iflmbd.gt.1) then
    write(nfecra,2053) iflmbd
    iok = iok + 1
  endif
  if (iangbd.lt.0 .or. iangbd.gt.1) then
    write(nfecra,2054) iangbd
    iok = iok + 1
  endif
  if (ivitbd.lt.0 .or. ivitbd.gt.1) then
    write(nfecra,2055) ivitbd
    iok = iok + 1
  endif
  if (nusbor.lt.0 .or. nusbor.gt.nusbrd) then
    write(nfecra,2056) nusbrd, nusbor
    iok = iok + 1
  endif

if (iok.ne.0) call csexit (1)
              !==========

  irf = 0

  if (inbrbd.eq.1) then
    irf = irf + 1
    if (imoybr(irf).eq.2) then
      write(nfecra,2060) inbrbd, imoybr(irf)
    endif
  endif
  if (iflmbd.eq.1) then
    irf = irf + 1
    if (imoybr(irf).eq.2 .and. inbrbd.eq.0) then
      iok = iok + 1
      WRITE(NFECRA,2061) NOMBRD(IRF), 'IFLMBD',  IFLMBD,          &
                         imoybr(irf), inbrbd
    endif
  endif
  if (iangbd.eq.1) then
    irf = irf + 1
    if (imoybr(irf).eq.2 .and. inbrbd.eq.0) then
      iok = iok + 1
      WRITE(NFECRA,2061) NOMBRD(IRF), 'IANGBD',  IANGBD,          &
                         imoybr(irf), inbrbd
    endif
  endif
  if (ivitbd.eq.1) then
    irf = irf + 1
    if (imoybr(irf).eq.2 .and. inbrbd.eq.0) then
      iok = iok + 1
      WRITE(NFECRA,2061) NOMBRD(IRF), 'IVITBD',  IVITBD,          &
                         imoybr(irf), inbrbd
    endif
  endif
  if (iphyla.eq.2 .and. iencra.eq.1 .and. iencbd.eq.1) then
    irf = irf + 1
    if (imoybr(irf).eq.2 .and. inbrbd.eq.0) then
      iok = iok + 1
      WRITE(NFECRA,2061) NOMBRD(IRF), 'IENCBD',  IENCBD,          &
                         imoybr(irf), inbrbd
    endif
  endif
  if (nusbor.gt.0) then
    do ii = 1,nusbor
      irf = irf + 1
      if (imoybr(irf).eq.2 .and. inbrbd.eq.0) then
        iok = iok + 1
        write(nfecra,2062) nombrd(irf), ii, imoybr(irf), inbrbd
      endif
    enddo
  endif

  do ii = 1, irf
    if (imoybr(ii).ne.0 .and.                                     &
        imoybr(ii).ne.1 .and. imoybr(ii).ne.2) then
    iok = iok + 1
      write(nfecra,2063) imoybr(irf), nombrd(irf)
  endif
  enddo

endif

if (iok.ne.0) call csexit (1)
              !==========

!===============================================================================
! 3. INITIALISATIONS DES VARIABLES EN COMMON

!                            ATTENTION :
!                            ^^^^^^^^^^^

!    CES INITIALISATIONS NE DOIVENT ETRE MODIFIEES PAR L'UTILISATEUR

!===============================================================================



! 3.1 GENERALITES (D'AUTRES INITIALISATIONS SONT FAITES DANS LAGLEC)


!     NOMBRE DE PASSAGES ABSOLUS DANS LE MODULE LAGRANGIEN

iplas = 0

!     NOMBRE DE PASSAGES RELATIFS DANS LE MODULE LAGRANGIEN

iplar = 0

!     PAS DE TEMPS LAGRANGIEN (LAGUNE) : Par defaut le pas de temps
!                                 de reference de la phase continue
dtp = dtref

!     TEMPS COURANT PHYSIQUE LAGRANGIEN

ttclag = 0.d0

!     INDICATEUR D'ERREUR (LAGCEL)

ierr = 0

!    NBPART/DNBPAR : NOMBRE DE PARTICULES PRESENTES DANS LE DOMAINE
!                      DE CALCUL A CHAQUE ITERATION

nbpart = 0
dnbpar = 0.d0

!     NBPERR/DNBPER : NOMBRE DE PARTICULES ELIMINES EN ERREUR

nbperr = 0
dnbper = 0.d0

!     NBPDEP/DNBDEP : NOMBRE DE PARTICULES DEPOSEES

nbpdep = 0
dnbdep = 0.d0


!     NBPERT : NOMBRE DE PARTICULES ELIMINEES EN ERREUR DANS
!              LE CALCUL DEPUIS LE DEBUT SUITE COMPRISE

nbpert = 0

!     NBPTOT : NOMBRE DE PARTICULES DU CALCUL (SUITES COMPRISES)

nbptot = 0

!     NBPOUT/DNBPOU : Contient les particules sorties de facon normale,
!                       plus les particules sorties en erreur de reperage.

nbpout = 0
dnbpou = 0.d0

!     NDEPOT : Nombre de particules deposees definitivement
!               dont on garde une trace en memoire pour le
!               post-processing en mode deplacement.

ndepot = 0

!     NPCLON/DNPCLO : NOMBRE DE NOUVELLES PARTICULES PAR CLONNAGE
!     NPKILL/DNPCSU : NOMBRE DE PARTICULES VICTIMES DE LA ROULETTE RUSSE
!     NPCSUP/DNPKIL : NOMBRE DE PARTICULES QUI ON SUBI LE CLONNAGE

npclon = 0
dnpclo = 0.d0

npcsup = 0
dnpcsu = 0.d0

npkill = 0
dnpkil = 0.d0

!     NPENCR/DNPENC : nombre de grains de charbon "encrasses"

npencr = 0
dnpenc = 0.d0


!     CONDITIONS AUX LIMITES

do ii = 1,nflagm
  ilflag(ii) = 0
  iusncl(ii) = 0
  iusclb(ii) = 0
  iusmoy(ii) = 0
  deblag(ii) = 0
enddo

!     STATISTIQUES VOLUMIQUES

!     Nombre de pas de temps DEPUIS LE DEBUT DU CALCUL STATIONNAIRES
!     des stats

npst =  0

!     Nombre de pas de temps total des stats depuis le debut
!     du calcul, partie instationnaire comprise

npstt = 0

!     Temps physique des stats

tstat = 0.d0

!     STATISTIQUES AUX FRONTIERES

!     Nombre de pas de temps DEPUIS LE DEBUT DU CALCUL STATIONNAIRES
!     des stats aux frontieres

npstf =  0

!     Nombre de pas de temps total des stats aux frontieres
!     depuis le debut du calcul, partie instationnaire comprise

npstft = 0

!     Temps physique des stats aux frontieres

tstatp = 0.d0

!     COUPLAGE RETOUR

!     Nombre de pas de temps DEPUIS LE DEBUT DU CALCUL STATIONNAIRES
!     des termes sources pour le couplage retour

npts =  0

!     Initialisation du sous-pas

nor = 0

!     NOMBRE D'ENREGISTREMENT POUR LE POST DEPLACEMENT (ENSWAF)

itlag = 0

!     TEMPS PHYSIQUE LAGRANGIEN POUR LE POST (ENSWAF)

do ii = 1,9999
  timlag(ii) = 0.d0
enddo

!     FINALISATION DE LA LISTE DE PARTICULES A VISUALISER

do ii = nbvis+1, nliste
  liste (ii) = -1
enddo

!    les trous, les repetitions dans le tableau LISTE seront
!    suprimes, et les numeros seront ranges par ordre croissant.

call lagtri
!==========

! ------------------------------------------------
! 3.2 DIMENSIONS DES TABLEAUX LIEES AUX PARTICULES
! ------------------------------------------------

!-->  NOMBRE MINIMAL DE VARIABLES LIEES AUX PARTICULES

!     NVP   : Variables sur les particules avec equation (ETTP et ETTPA)
!     NVEP  : Variables d'etat (reels)   sur les particules (TEPA)
!     NIVEP : Variables d'etat (entiers) sur les particules (ITEPA)

nvp   = 12
nvep  = 2
nivep = 2

if (nbclst.gt.0) then

!       --> Statistique par classe :
!          1 VARIABLE D'ETAT ENTIERE SUPPLEMENTAIRE : Numero de la
!          classe statique a laquelle appartient la particule

  nivep = nivep +1

endif


if (iphyla.eq.1) then

!       --> EQUATION SUR LA TEMPERATURE :
!           3 VARIABLES SUPPLEMENTAIRES Tp, Tf, Cp dans ETTP et ETTPA

  if (itpvar.eq.1) nvp = nvp + 3

!       --> EQUATION SUR LA TEMPERATURE ET RAYONNEMENT : EMISSIVITE

  if (itpvar.eq.1 .and. iirayo.gt.0) nvep = nvep + 1

else if (iphyla.eq.2) then

!    --> CHARBON :

!     ETTP et ETTPA :
!     -------------
!       5 VARIABLES SUPPLEMENTAIRES Tp, Tf, Mch, Mck, Cp
  nvp = nvp + 5

!     TEPA :
!     ----
!       3 VARIABLES D'ETATS REELLES SUPPLEMENTAIRES : Dck, D0P, R0P

  nvep = nvep + 3

!     ITEPA :
!     -----
!       1 VARIABLE D'ETAT ENTIERE SUPPLEMENTAIRE : Numero du charbon
  nivep = nivep + 1

endif

! Modele de deposition : 2 tableaux supp dans TEPA  : YPLUS , JRINPF
!                        6 tableaux supp dans ITEPA : MARKO , DIEL,DFAC, DIFEL, TRAJ, JPTDET

if ( idepst .eq. 1 ) then
  nvep  = nvep  + 2
  nivep = nivep + 7
endif

!-->  VARIABLES UTILISATEURS SUPPLEMENTAIRES : NVLS

if (nvls.gt.0) nvp = nvp + nvls

!-->  NVP1 represente le nombre de variables sur les particules en
!     enlevant position, vitesse particule et vitesse fluides (TSVAR)

nvp1 = nvp - 9


! 3.3 DEFINITION DES POINTEURS SUR LES VARIABLES LIEES AUX PARTICULES


!   3.3.1 TABLEAU ETTP
!   ~~~~~~~~~~~~~~~~~~

!   Attention il faut que
!   JMP, JDP, JXP, JYP, JZP, JUP, JVP, JWP, JUF, JVF, JWF
!   soient les derniers pointeurs pour les tableaux ETTP et ETTPA
!   a cause du remplissage et de la dimension du tableau
!   TSVAR(NBPMAX,NVP1) et de son mode de lecture/ecriture
!   (cf. LAGITG).

!    JXP,JYP,JZP  : COORDONNES DE LA POSITION DE LA PARTICULE NPT
!    JUP,JVP,JWP  : COMPOSANTES DE LA VITESSE ABSOLUE
!    JUF,JVF,JWF  : COMPOSANTES DE LA VITESSE DU FLUIDE VU

!    JMP,JDP      : MASSE, DIAMETRE
!    JTP,JTF,JCP  : TEMPERATURE PARTICULE ET FLUIDE ET CHALEUR SPECIFIQUE
!    JTAUX        : AUXILIAIRE DE CALCUL UTILE EN PARALLELE
!    JVLS(NUSVAR) : VARIABLE SUPPLEMENTAIRES

!   Charbon
!   -------
!    JHP  : Temperature en degres Celsius grain de charbon
!    JMCH : MASSE DE CHARBON REACTIF
!    JMCK : MASSE DE COKE

jtp  = 0
jtf  = 0
jcp  = 0
jhp  = 0
jmch = 0
jmck = 0
jtaux = 0
do ii = 1,nusvar
  jvls(ii) = 0
enddo

irf  = 0

if (iphyla.eq.1) then

  if (itpvar.eq.1) then
    jtp = irf + 1
    jtf = jtp + 1
    jcp = jtf + 1
    irf = jcp
  endif

else if (iphyla.eq.2) then

  jhp  = irf  + 1
  jtf  = jhp  + 1
  jmch = jtf  + 1
  jmck = jmch + 1
  jcp  = jmck + 1
  irf  = jcp

endif

if (nvls.gt.0) then
  do ii = 1,nvls
    irf = irf + 1
    jvls(ii) = irf
  enddo
  irf = jvls(nvls)
endif

jmp = irf + 1
jdp = jmp + 1
jxp = jdp + 1
jyp = jxp + 1
jzp = jyp + 1
jup = jzp + 1
jvp = jup + 1
jwp = jvp + 1
juf = jwp + 1
jvf = juf + 1
jwf = jvf + 1
jtaux = jwf + 1
irf = jtaux

if (irf.gt.nvp) then
  write(nfecra,3004) irf, nvp
  call csexit(1)
endif

!   3.3.2 TABLEAU TEPA
!   ~~~~~~~~~~~~~~~~~~

!     JRTSP       : TEMPS DE SEJOUR DES PARTICULES
!     JRPOI       : POIDS DES PARTICULES
!     JREPS       : EMISSIVITE DES PARTICULES

!   Charbon
!   -------
!     JRDCK       : DIAMETRE DU COEUR RETRECISSANT
!     JRD0P       : DIAMETRE INITIAL DES PARTICULES
!     JRR0P       : MASSE VOLUMIQUE INITIALE DES PARTICULES

jreps = 0
jrdck = 0
jrr0p = 0
jrr0p = 0

jrtsp = 1
jrpoi = 2
irf   = jrpoi

if (iphyla.eq.1 .and. itpvar.eq.1 .and. iirayo.gt.0) then
  jreps = irf   + 1
  irf   = jreps
endif

if (iphyla.eq.2) then
  jrdck = irf   + 1
  jrd0p = jrdck + 1
  jrr0p = jrd0p + 1
  irf   = jrr0p
endif

! Modele de deposition : 2 tableaux supp dans TEPA  : YPLUS et DX

if ( idepst .eq. 1 ) then
  jryplu = irf    +1
  jrinpf = jryplu +1
  irf    = jrinpf
endif


if (irf.ne.nvep) then
  write(nfecra,3005) irf, nvep
  call csexit(1)
endif



!   3.3.3 TABLEAU ITEPA
!   ~~~~~~~~~~~~~~~~~~~

!     JISOR       : MAILLE D'ARRIVEE

!   Statistique par classe
!   ----------------------

!     JCLST       : classe (statique) a laquelle la particule appartient

!   Charbon
!   -------
!     JINCH       : NUMERO DU CHARBON DE LA PARTICULE

jinch = 0

jisor = 1

jgnum = 2
irf = jgnum

if (nbclst .gt. 0) then
  jclst = irf + 1
  irf = jclst
endif

if (iphyla.eq.2) then
  jinch = irf + 1
  irf   = jinch
endif

! Modele de Deposition :  tableaux supp dans ITEPA : MARKO
!                                                    JDIEL
!                                                    JDFAC


if ( idepst .eq. 1 ) then
  jimark = irf    + 1
  jdiel  = jimark + 1
  jdfac  = jdiel  + 1
  jdifel = jdfac  + 1
  jtraj  = jdifel + 1
  jptdet = jtraj  + 1
  jinjst = jptdet + 1
  irf    = jinjst
endif



if (irf.ne.nivep) then
  write(nfecra,3006) irf, nivep
  call csexit(1)
endif


! 3.4 DEFINITION DES POINTEURS LIES A L'INJECTION DES PARTICULES


!   3.4.1 TABLEAU RUSLAG (DONNEES D'ENTREE)
!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!     IUNO  : Norme de la vitesse
!     IUPT  : U par classe et zones
!     IVPT  : V par classe et zones
!     IWPT  : W par classe et zones
!     IDEBT : Debit
!     IPOIT : Poids de la particule
!     IDPT  : Diametre
!     IVDPT : Ecart-type du diametre
!     ITPT  : Temperature
!     ICPT  : Cp
!     IEPSI : Emissivite
!     IROPT : Masse volumique
!     IHPT  : Temperature
!     IMCHT : Masse de charbon reactif
!     IMCKT : Masse de coke
!     IDCKT : Diametre du coeur retrecissant

iuno   = 1
iupt   = iuno  + 1
ivpt   = iupt  + 1
iwpt   = ivpt  + 1
itpt   = iwpt  + 1
idpt   = itpt  + 1
ivdpt  = idpt  + 1
iropt  = ivdpt + 1
icpt   = iropt + 1
iepsi  = icpt  + 1
ipoit  = iepsi + 1
idebt  = ipoit + 1
irf    = idebt

! Specifique Charbon

ihpt   = irf   + 1
imcht  = ihpt  + 1
imckt  = imcht + 1
irf    = imckt

if (irf.gt.ndlagm) then
  write(nfecra,3001) irf, ndlagm
  call csexit(1)
endif

!   3.4.2 TABLEAU IUSLAG (DONNEES D'ENTREE)
!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!     IJNBP  : nbr de part par classe et zones frontieres
!     IJFRE  : frequence d'injection
!     ICLST : numero de groupe auquel appartient la particule
!             (uniquement si on souhaite des statistiques par groupe)
!     IJUVW  : type de condition vitesse
!     IJPRTP : profil de temperature pour les particules
!     IJPRPD : profil de poids statistiques
!     IJPRDP : profil de diametre pour les particules
!     INUCHL : numero du charbon de la particule (si IPHYLA=2)

ijnbp  = 1
ijfre  = ijnbp  + 1
iclst  = ijfre  + 1
ijuvw  = iclst  + 1
ijprtp = ijuvw  + 1
ijprdp = ijprtp + 1
ijprpd = ijprdp + 1
irf    = ijprpd

! Specifique Charbon

inuchl = irf + 1
irf    = inuchl

if (irf.gt.ndlaim) then
  write(nfecra,3002) irf, ndlaim
  call csexit(1)
endif

! -----------------------------------------------
! 3.5 DIMENSIONS DU TABLEAU LIEES AUX STATISTIQUES
! -----------------------------------------------

!     NVLSTA   : NOMBRE DE VARIABLES LIEES AUX STATISTIQUES
!                DIMENSION DU TABLEAU STATIS
!     NVLSTS   : VARIABLES LIEES AUX STATISTIQUES SUPPLEMENTAIRES






! --------------------------------------------------
! 3.6 DEFINITION DES POINTEURS LIES AUX STATISTIQUES
! --------------------------------------------------

!     ILVX,ILVY,ILVZ    : Vitesse
!     ILFV              : Concentrations volumiques
!     ILPD              : Somme des poids statistiques

!     ILTP              : Temperature
!     ILDP              : Diametre
!     ILMP              : Masse

!     ILHP              : Temperature
!     ILMCH             : Masse de charbon reactif
!     ILMCK             : Masse de coke
!     ILMDK             : Diametre du coeur retrecissant

!     ILVU(NUSSTA)      : Statistiques utilisateur

ilfv = 0
ilvx = 0
ilvy = 0
ilvz = 0
ilts = 0

if (istala.eq.1) then

  irf = 0

  if (iactfv.eq.1) then
     irf = irf + 1
     ilfv  = irf
     nomlag(ilfv) = 'Part_vol_frac'
     nomlav(ilfv) = 'var_Part_vol_frac'
  endif

  if (iactvx.eq.1) then
     irf = irf + 1
     ilvx = irf
     nomlag(ilvx) = 'Part_velocity_X'
     nomlav(ilvx) = 'var_Part_velocity_X'
  endif

  if (iactvy.eq.1) then
     irf = irf + 1
     ilvy = irf
     nomlag(ilvy) = 'Part_velocity_Y'
     nomlav(ilvy) = 'var_Part_velocity_Y'
  endif

  if (iactvz.eq.1) then
     irf = irf + 1
     ilvz = irf
     nomlag(ilvz) = 'Part_velocity_Z'
     nomlav(ilvz) = 'var_Part_velocity_Z'
  endif

  if (iactts.eq.1) then
     irf = irf + 1
     ilts = irf
     nomlag(ilts) = 'Part_resid_time'
     nomlav(ilts) = 'var_Part_resid_time'
  endif

  iltp  = 0
  ildp  = 0
  ilmp  = 0
  ilhp  = 0
  ilmch = 0
  ilmck = 0
  ildck = 0

  do ii = 1,nussta
     ilvu(ii) = 0
  enddo

  if (iphyla.eq.1) then

     if (itpvar.eq.1) then
        iltp = irf + 1
        irf  = iltp
     endif

     if (idpvar.eq.1) then
        ildp = irf + 1
        irf  = ildp
     endif

     if (impvar.eq.1) then
        ilmp = irf + 1
        irf  = ilmp
     endif

  else if (iphyla.eq.2) then

     ilhp  = irf   + 1
     ilmch = ilhp  + 1
     ilmck = ilmch + 1
     ildck = ilmck + 1
     irf   = ildck

  endif

  if (nvlsts.gt.0) then
     do ii = 1,nvlsts
        ilvu(ii) = irf + ii
     enddo
     irf = irf + nvlsts
  endif

  ilpd  = irf  + 1
  nomlag(ilpd) = 'Part_statis_weight'

  nvlsta = ilpd

endif

if(nvlsta.gt.nvplmx) then
  write(nfecra,3003) nvlsta, nvplmx
  call csexit(1)
endif


! 3.7 DEFINITION DES POINTEURS LIES AUX STATISTIQUES AUX FRONTIERES


!     INBRBD : NOMBRE D'INTERACTIONS PARTICULES/FRONTIERES
!     IFLMBD : FLUX DE MASSE PARTICULAIRE
!     IANGBD : ANGLE VITESSE
!     IVITBD : VITESSE DE LA PARTICULE
!     IENCBD : MASSE DE GRAINS DE CHARBON ENCRASSES
!     NUSBOR : INFORMATIONS UTILISATEUR SUPPLEMENTAIRES
!     NVISBR : NOMBRE TOTAL D'INTERACTIONS A ENREGISTRER

if (iensi3.eq.1) then

  irf = 0

  if (inbrbd.eq.1) then
    irf = irf + 1
    inbr = irf
    nombrd(inbr) = 'Part_impact_number'
    imoybr(inbr) = 0
  endif

  if (iflmbd.eq.1) then
    irf = irf + 1
    iflm = irf
    nombrd(iflm) = 'Part_bndy_mass_flux'
    imoybr(iflm) = 1
  endif

  if (iangbd.eq.1) then
    irf = irf + 1
    iang = irf
    nombrd(iang) = 'Part_impact_angle'
    imoybr(iang) = 2
  endif

  if (ivitbd.eq.1) then
    irf = irf + 1
    ivit = irf
    nombrd(ivit) = 'Part_impact_velocity'
    imoybr(ivit) = 2
  endif

  if (iphyla.eq.2 .and. iencra.eq.1 .and. iencbd.eq.1) then
    irf = irf + 1
    ienc = irf
    nombrd(ienc) = 'foulingMass'
    imoybr(ienc) = 0
  endif

  if (nusbor.gt.0) then
    do ii = 1,nusbor
      irf = irf + 1
      iusb(ii) = irf
      write(nombrd(iusb(ii)),'(a8,i4.4)') 'addRec',II
      imoybr(iusb(ii)) = 0
    enddo
  endif

  nvisbr = irf

else

  nvisbr = 0

endif


! 3.8 DEFINITION DES POINTEURS LIES AUX TERMES SOURCES LAGRANGIEN
!     POUR COUPLAGE RETOUR


!     Nombre de termes sources de couplage-retour

ntersl = 0

irf = 0
itsvx  = 0
itsvy  = 0
itsvz  = 0
itsli  = 0
itske  = 0
itsr11 = 0
itsr12 = 0
itsr13 = 0
itsr22 = 0
itsr23 = 0
itsr33 = 0
itsmas = 0
itste  = 0
itsti  = 0
do icha = 1,ncharm2
  itsmv1(icha) = 0
  itsmv2(icha) = 0
enddo
itsco  = 0
itsfp4 = 0

! Dynamique : Vitesse + Turbulence

if (ltsdyn.eq.1) then

  ntersl = ntersl + 4

  itsvx  = irf   + 1
  itsvy  = itsvx + 1
  itsvz  = itsvy + 1
  itsli  = itsvz + 1
  irf    = itsli

  if (itytur.eq.2 .or. iturb.eq.50                  &
       .or. iturb.eq.60) then
! K-eps, v2f et k-omega
    ntersl = ntersl + 1

    itske  = irf    + 1
    irf    = itske

  else if (itytur.eq.3) then
! RIJ
    ntersl = ntersl + 6

    itsr11 = irf    + 1
    itsr12 = itsr11 + 1
    itsr13 = itsr12 + 1
    itsr22 = itsr13 + 1
    itsr23 = itsr22 + 1
    itsr33 = itsr23 + 1
    irf    = itsr33
  else
    write(nfecra,3010) iilagr, ltsdyn, iturb
    call csexit (1)
    !==========
  endif

endif

! Modele de depot

if (idepst.eq.1 .and. nordre.eq.2) then
    write(nfecra,3014)
    call csexit (1)
    !==========
endif

! Masse

if (ltsmas.eq.1) then

  ntersl = ntersl + 1

  itsmas = irf    + 1
  irf    = itsmas

endif

! Thermique

if (ltsthe.eq.1) then

  if (iphyla.eq.1) then

!     Temperature

    if (itpvar.eq.1) then

      ntersl = ntersl + 2

      itste  = irf    + 1
      itsti  = itste  + 1
      irf    = itsti

    endif

!     Charbon

  else if (iphyla.eq.2) then

    ntersl = ntersl + 4 + 2*ncharb

    itste  = irf    + 1
    itsti  = itste  + 1

    do icha = 1,ncharb
      itsmv1(icha) = itsti + icha
    enddo

    do icha = 1,ncharb
      itsmv2(icha) = itsmv1(ncharb) + icha
    enddo

    itsco  = itsmv2(ncharb) + 1
    itsfp4 = itsco          + 1
    irf    = itsfp4

  endif

endif

!===============================================================================

!--------
! FORMATS
!--------

 1010 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========                                               ',/,&
'@    L''INDICATEUR SUR LE MODULE LAGRANGIEN A UNE VALEUR     ',/,&
'@       NON PERMISE (LAGOPT).                                ',/,&
'@                                                            ',/,&
'@    IILAGR DEVRAIT ETRE UN ENTIER EGAL A 0, 1, 2 OU 3.      ',/,&
'@       IL VAUT ICI IILAGR = ', I10                           ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier la valeur de IILAGR dans la subroutine USLAG1.   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 1012 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========                                               ',/,&
'@    L''INDICATEUR SUR LE MODULE LAGRANGIEN IILAGR ET        ',/,&
'@      L''INDICATEUR DE SUITE ONT DES VALEURS INCOMPATIBLES  ',/,&
'@      (LAGOPT).                                             ',/,&
'@                                                            ',/,&
'@       IILAGR = ', I10                                       ,/,&
'@       ISUITE = ', I10                                       ,/,&
'@                                                            ',/,&
'@  Le module lagrangien est active en mode champs figes,     ',/,&
'@   alors que le calcul de la phase continue n''est pas      ',/,&
'@   une suite.                                               ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier la valeur de IILAGR dans la subroutine USLAG1.   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 1013 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========                                               ',/,&
'@    LA PHYSIQUE PARTICULIERE COMBUTION CHARBON PULVERISE    ',/,&
'@      COUPLE AU TRANSPORT LAGRANGIEN DES PARTICULES         ',/,&
'@      DE CHARBON EST ACTIVEE (USPPMO), ALORS QUE LE COUPLAGE',/,&
'@      RETOUR DE LA PHASE DISPERSEE SUR LE PHASE CONTINUE    ',/,&
'@      N''EST PAS ENCLENCHE (LAGOPT).                        ',/,&
'@                                                            ',/,&
'@       IILAGR = ', I10                                       ,/,&
'@       IPPMOD(ICPL3C) = ', I10                               ,/,&
'@                                                            ',/,&
'@  Le module lagrangien doit etre active en mode couplage    ',/,&
'@   retour pour etre couple avec la combustion d''une        ',/,&
'@   flamme de charbon pulverise en phase continue.           ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier la valeur de IILAGR dans la subroutine USLAG1 et ',/,&
'@  verifier la valeur de IPPMOD dans la subroutine USPPMO.   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 1014 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========                                               ',/,&
'@    L''INDICATEUR SUR LE MODULE LAGRANGIEN IILAGR ET        ',/,&
'@      LE CHOIX DU TYPE DE PAS DE TEMPS IDTVAR               ',/,&
'@      ONT DES VALEURS INCOMPATIBLES (LAGOPT).               ',/,&
'@                                                            ',/,&
'@       IILAGR = ', I10                                       ,/,&
'@       IDTVAR = ', I10                                       ,/,&
'@                                                            ',/,&
'@  Le module lagrangien ne peut pas etre active avec un pas  ',/,&
'@   de temps variable en temps et en espace. Seuls les pas   ',/,&
'@   de temps uniforme et constant, et variable en temps et   ',/,&
'@   uniforme en espace sont possibles.                       ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier la valeur de IILAGR dans la subroutine USLAG1 et ',/,&
'@  verifier la valeur de IDTVAR dans la subroutine USINI1.   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 1020 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========                                               ',/,&
'@    L''INDICATEUR DE SUITE DU MODULE LAGRANGIEN A UNE       ',/,&
'@       VALEUR NON PERMISE (LAGOPT).                         ',/,&
'@                                                            ',/,&
'@    ISUILA DEVRAIT ETRE UN ENTIER EGAL A 0 OU 1             ',/,&
'@       IL VAUT ICI ISUILA = ', I10                           ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier la valeur de IILAGR dans la subroutine USLAG1.   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 1021 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ALERTE A L''EXECUTION DU MODULE LAGRANGIEN  ',/,&
'@    =========                                               ',/,&
'@                                                  (LAGOPT). ',/,&
'@                                                            ',/,&
'@  Le module lagrangien est active en suite de calcul,       ',/,&
'@   alors que le calcul de la phase continue n''est pas      ',/,&
'@   une suite.                                               ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier la valeur de ISUILA dans la subroutine USLAG1.   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 1022 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========                                               ',/,&
'@    L''INDICATEUR DE SUITE DE CALCUL SUR LES STATISTIQUES   ',/,&
'@       VOLUMIQUE ET AUX FRONTIERES, AINSI QUE SUR LES       ',/,&
'@       TERMES SOURCES DE COUPLAGES RETOUR                   ',/,&
'@       A UNE VALEUR NON PERMISE (LAGOPT).                   ',/,&
'@                                                            ',/,&
'@    ISUIST DEVRAIT ETRE UN ENTIER EGAL A 0 OU 1             ',/,&
'@       IL VAUT ICI ISUIST = ', I10                           ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier la valeur de ISUIST dans la subroutine USLAG1.   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 1030 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========                                               ',/,&
'@    L''INDICATEUR DES MODELES PHYSIQUES LIES AUX PARTICULES ',/,&
'@       A UNE VALEUR NON PERMISE (LAGOPT).                   ',/,&
'@                                                            ',/,&
'@    IPHYLA DEVRAIT ETRE UN ENTIER EGAL A 0 1 OU 2           ',/,&
'@       IL VAUT ICI IPHYLA = ', I10                           ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier la valeur de IPHYLA dans la subroutine USLAG1.   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 1031 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========                                               ',/,&
'@    L''INDICATEUR SUR L''EQUATION DU DIAMETRE DES           ',/,&
'@       PARTICULES A UNE VALEUR NON PERMISE (LAGOPT).        ',/,&
'@                                                            ',/,&
'@     IDPVAR DEVRAIT ETRE UN ENTIER EGAL A 0 OU 1            ',/,&
'@       IL VAUT ICI IDPVAR = ', I10                           ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier la valeur de IDPVAR dans la subroutine USLAG1.   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 1032 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========                                               ',/,&
'@    L''INDICATEUR SUR L''EQUATION DE LA TEMPERATURE DES     ',/,&
'@       PARTICULES A UNE VALEUR NON PERMISE (LAGOPT).        ',/,&
'@                                                            ',/,&
'@     ITPVAR DEVRAIT ETRE UN ENTIER EGAL A 0 OU 1            ',/,&
'@       IL VAUT ICI ITPVAR = ', I10                           ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier la valeur de ITPVAR dans la subroutine USLAG1.   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 1033 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========                                               ',/,&
'@    L''INDICATEUR SUR L''EQUATION DE LA MASSE DES           ',/,&
'@       PARTICULES A UNE VALEUR NON PERMISE (LAGOPT).        ',/,&
'@                                                            ',/,&
'@     IMPVAR DEVRAIT ETRE UN ENTIER EGAL A 0 OU 1            ',/,&
'@       IL VAUT ICI IMPVAR = ', I10                           ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier la valeur de IMPVAR dans la subroutine USLAG1.   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 1034 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========                                               ',/,&
'@    L''INDICATEUR SUR L''EQUATION DE LA TEMPERATURE DES     ',/,&
'@       PARTICULES EST ACTIVE (ITPVAR = ',I10,'),'            ,/,&
'@       ALORS QU''AUCUN SCALAIRE THERMIQUE N''EST DISPONIBLE ',/,&
'@                                                            ',/,&
'@     ISCALT DEVRAIT ETRE UN ENTIER SUPERIEUR OU EGAL 1      ',/,&
'@       IL VAUT ICI ISCALT = ',I10                            ,/,&
'@                                                            ',/,&
'@  La valeur de ISCALT est renseignee automatiquement        ',/,&
'@    si une physique particuliere est activee dans USPPMO.   ',/,&
'@                                                            ',/,&
'@  L''utilisateur doit renseigner ISCALT dans USINI1 si      ',/,&
'@    aucune physique particuliere n''est activee dans USPPMO.',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier la valeur de ITPVAR dans la subroutine USLAG1,   ',/,&
'@  verifier la valeur de ISCALT dans la subroutine USINI1 et ',/,&
'@  verifier la valeur de IPPMOD dans la subroutine USPPMO.   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 1036 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========                                               ',/,&
'@    LA CHALEUR MASSIQUE D''INITIALISATION DES PARTICULES    ',/,&
'@       DEJA PRESENTE DANS LE DOMAINE DE CALCUL              ',/,&
'@       A UNE VALEUR NON PERMISE (LAGOPT).                   ',/,&
'@                                                            ',/,&
'@     CPPART DEVRAIT ETRE UN REEL STRICTEMENT POSITIF        ',/,&
'@       IL VAUT ICI CPPART = ', E14.5                         ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier la valeur de CPPART dans la subroutine USLAG1.   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 1037 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========                                               ',/,&
'@    LA TEMPERATURE D''INITIALISATION DES PARTICULES         ',/,&
'@       DEJA PRESENTE DANS LE DOMAINE DE CALCUL              ',/,&
'@       A UNE VALEUR NON PERMISE (LAGOPT).                   ',/,&
'@                                                            ',/,&
'@     TPART DEVRAIT ETRE UN REEL SUPERIEUR A ',E14.5          ,/,&
'@       (EN DEGRES CELSIUS)                                  ',/,&
'@       IL VAUT ICI TPART = ', E14.5                          ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier la valeur de TPART dans la subroutine USLAG1.    ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 1040 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========                                               ',/,&
'@    L''INDICATEUR SUR L''ENCRASSEMENT DES PARTICULES        ',/,&
'@       DE CHARBON A UNE VALEUR NON PERMISE (LAGOPT).        ',/,&
'@                                                            ',/,&
'@     IENCRA DEVRAIT ETRE UN ENTIER EGAL A 0 OU 1            ',/,&
'@       IL VAUT ICI IENCRA = ', I10                           ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier la valeur de IENCRA dans la subroutine USLAG1.   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 1041 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========                                               ',/,&
'@    L''INDICATEUR SUR L''ENCRASSEMENT DES PARTICULES        ',/,&
'@       DE CHARBON EST ACTIVE (IENCRA = ',I10,')'             ,/,&
'@       AVEC UNE VALEUR DE VISCOSITE CRITIQUE                ',/,&
'@       NON PERMISE (LAGOPT).                                ',/,&
'@                                                            ',/,&
'@     VISREF DEVRAIT ETRE UN REEL STRICTEMENT POSITIF (Pa.s) ',/,&
'@       IL VAUT ICI VISREF = ', E14.5                         ,/,&
'@       POUR LE CHARBON :' , I10                              ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier la valeur de VISREF dans la subroutine USLAG1.   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 1042 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========                                               ',/,&
'@    L''INDICATEUR SUR L''ENCRASSEMENT DES PARTICULES        ',/,&
'@       DE CHARBON EST ACTIVE (IENCRA = ',I10,')'             ,/,&
'@       AVEC UNE VALEUR DE TEMPERATURE SEUIL                 ',/,&
'@       NON PERMISE (LAGOPT).                                ',/,&
'@                                                            ',/,&
'@     TPRENC DEVRAIT ETRE UN REEL SUPERIEUR A ',E14.5         ,/,&
'@       (EN DEGRES CELSIUS)                                  ',/,&
'@       IL VAUT ICI TPRENC = ', E14.5                         ,/,&
'@       POUR LE CHARBON :' , I10                              ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier la valeur de TPRENC dans la subroutine USLAG1.   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 1043 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========                                               ',/,&
'@    LA PHYSIQUE PARTICULIERE COMBUTION CHARBON PULVERISE    ',/,&
'@      COUPLE AU TRANSPORT LAGRANGIEN DES PARTICULES         ',/,&
'@      DE CHARBON EST ACTIVEE (USPPMO), ALORS QUE L''OPTION  ',/,&
'@      TRANSPORT DE PARTICULE DE CHARBON                     ',/,&
'@      N''EST PAS ENCLENCHEE (LAGOPT).                       ',/,&
'@                                                            ',/,&
'@       IPHYLA = ', I10                                       ,/,&
'@       IPPMOD(ICPL3C) = ', I10                               ,/,&
'@                                                            ',/,&
'@  Le module lagrangien doit etre active en mode transport   ',/,&
'@   de particules de charbon pour etre couple avec la        ',/,&
'@   combustion d''une flamme de charbon pulverise en phase   ',/,&
'@   continue.                                                ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier la valeur de IPHYLA dans la subroutine USLAG1 et ',/,&
'@  verifier la valeur de IPPMOD dans la subroutine USPPMO.   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 1044 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========                                               ',/,&
'@    LE TRANSPORT LAGRANGIEN DE PARTICULES DE CHARBON        ',/,&
'@      EST ACTIVE (LAGOPT), ALORS QU''AUCUNE PHYSIQUE        ',/,&
'@      PARTICULIERE SUR LA COMBUSTION DU CHABON PULVERISE    ',/,&
'@      N''EST PAS ENCLENCHE (USPPMO).                        ',/,&
'@                                                            ',/,&
'@       IPHYLA = ', I10                                       ,/,&
'@       IPPMOD(ICPL3C) = ', I10                               ,/,&
'@       IPPMOD(ICP3PL) = ', I10                               ,/,&
'@                                                            ',/,&
'@  Le transport lagrangien de particule de charbon doit      ',/,&
'@   etre couple avec la combustion d''une flamme de charbon  ',/,&
'@   pulverise en phase continue.                             ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier la valeur de IPHYLA dans la subroutine USLAG1 et ',/,&
'@  verifier la valeur de IPPMOD dans la subroutine USPPMO.   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 1050 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========                                               ',/,&
'@    L''INDICATEUR SUR LE NOMBRE DE PARTIUCULES A TRAITER    ',/,&
'@      PAR LE MODULE LAGRANGIEN A UNE VALEUR                 ',/,&
'@       NON PERMISE (LAGOPT).                                ',/,&
'@                                                            ',/,&
'@    NBPMAX DEVRAIT ETRE UN ENTIER STRICTEMENT POSITIF       ',/,&
'@       IL VAUT ICI NBPMAX = ', I10                           ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier la valeur de NBPMAX dans la subroutine USLAG1.   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 1060 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========                                               ',/,&
'@    L''INDICATEUR SUR LE NOMBRE DE VARIABLES                ',/,&
'@       SUPPLEMENTAIRES LIEES AUX PARTICULES                 ',/,&
'@       A UNE VALEUR NON PERMISE (LAGOPT).                   ',/,&
'@                                                            ',/,&
'@     NVLS DEVRAIT ETRE UN ENTIER ENTRE 0 ET ',I10            ,/,&
'@       IL VAUT ICI NVLS = ', I10                             ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier la valeur de NVLS dans la subroutine USLAG1.     ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 1061 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========                                               ',/,&
'@    L''INDICATEUR SUR LE CARACTERE STATIONNAIRE DE          ',/,&
'@       L''ECOULEMENT DE LA PHASE CONTINUE                   ',/,&
'@       A UNE VALEUR NON PERMISE (LAGOPT).                   ',/,&
'@                                                            ',/,&
'@    ISTTIO DEVRAIT ETRE UN ENTIER EGAL A 0 OU 1             ',/,&
'@       IL VAUT ICI ISTTIO = ', I10                           ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier la valeur de ISTTIO dans la subroutine USLAG1.   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 1062 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========                                               ',/,&
'@    L''INDICATEUR SUR LE DECLENCHEMENT DU CALCUL            ',/,&
'@       STATIONNAIRE DES STATISTIQUES POUR UN COUPLAGE RETOUR',/,&
'@       A UNE VALEUR NON PERMISE (LAGOPT).                   ',/,&
'@                                                            ',/,&
'@    NSTITS DEVRAIT ETRE UN ENTIER SUPERIEUR OU EGAL A 1     ',/,&
'@       IL VAUT ICI NSTITS = ', I10                           ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier la valeur de NSTITS dans la subroutine USLAG1.   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 1063 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========                                               ',/,&
'@    L''INDICATEUR SUR LE COUPLAGE RETOUR SUR LA DYNAMIQUE   ',/,&
'@       A UNE VALEUR NON PERMISE (LAGOPT).                   ',/,&
'@                                                            ',/,&
'@    LTSDYN DEVRAIT ETRE UN ENTIER EGAL A 0 OU 1             ',/,&
'@       IL VAUT ICI LTSDYN = ', I10                           ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier la valeur de LTSDYN dans la subroutine USLAG1.   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 1064 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========                                               ',/,&
'@    L''INDICATEUR SUR LE COUPLAGE RETOUR SUR LA MASSE       ',/,&
'@       A UNE VALEUR NON PERMISE (LAGOPT).                   ',/,&
'@                                                            ',/,&
'@    LTSMAS DEVRAIT ETRE UN ENTIER EGAL A 0 OU 1             ',/,&
'@       IL VAUT ICI LTSMAS = ', I10                           ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier la valeur de LTSMAS dans la subroutine USLAG1.   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 1065 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========                                               ',/,&
'@    L''INDICATEUR SUR LE COUPLAGE RETOUR SUR LA THERMIQUE   ',/,&
'@       A UNE VALEUR NON PERMISE (LAGOPT).                   ',/,&
'@                                                            ',/,&
'@    LTSTHE DEVRAIT ETRE UN ENTIER EGAL A 0 OU 1             ',/,&
'@       IL VAUT ICI LTSTHE = ', I10                           ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier la valeur de LTSTHE dans la subroutine USLAG1.   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 1066 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========                                               ',/,&
'@    L''INDICATEUR SUR LE COUPLAGE RETOUR SUR LA DYNAMIQUE   ',/,&
'@       EST ACTIVE (LTSDYN = ',I10,') (LAGOPT),              ',/,&
'@       ALORS QUE LA PHASE PORTEUSE EST CALCULEE AVEC        ',/,&
'@       L''OPTION CHAMP FIGE  (ICCVFG = ',I10,') (USINI1).   ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier la valeur de LTSDYN dans la subroutine USLAG1 et ',/,&
'@  verifier la valeur de ICCVFG dans la subroutine USINI1.   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 1067 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========   (LAGOPT).                                   ',/,&
'@                                                            ',/,&
'@    L''INDICATEUR SUR LE COUPLAGE RETOUR EST ACTIVE         ',/,&
'@        IILAGR = ',I10                                       ,/,&
'@      ALORS QU''AUCUN COUPLAGE RETOUR N''EST ENCLENCHE      ',/,&
'@        DYNAMIQUE : LTSDYN = ',I10                           ,/,&
'@        THERMIQUE : LTSTHE = ',I10                           ,/,&
'@        MASSIQUE  : LTSMAS = ',I10                           ,/,&
'@                                                            ',/,&
'@    LES COUPLAGES RETOUR SUR LA THERMIQUE ET SUR LA MASSE   ',/,&
'@      NECESSITENT L''ACTIVATION D''UNE PHYSIQUE ADEQUATE    ',/,&
'@      ASSOCIEE AUX PARTICULES.                              ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier la valeur de IILAGR dans la subroutine USLAG1.   ',/,&
'@  Verifier la valeur de IPHYLA dans la subroutine USLAG1.   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 1070 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========                                               ',/,&
'@   L''INDICATEUR SUR LE LANCEMENT DU CALCUL DES STATISTIQUES',/,&
'@       A UNE VALEUR NON PERMISE (LAGOPT).                   ',/,&
'@                                                            ',/,&
'@    ISTALA DEVRAIT ETRE UN ENTIER EGAL A 0 OU 1             ',/,&
'@       IL VAUT ICI ISTALA = ', I10                           ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier la valeur de ISTALA dans la subroutine USLAG1.   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 1071 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========                                               ',/,&
'@    L''INDICATEUR DE SEUIL POUR LE CALCUL DES STATISTIQUES  ',/,&
'@       VOLUMIQUES A UNE VALEUR NON PERMISE (LAGOPT).        ',/,&
'@                                                            ',/,&
'@    SEUIL DEVRAIT ETRE UN REEL SUPERIEUR OU EGAL A 0.D0     ',/,&
'@       IL VAUT ICI SEUIL = ', E14.5                          ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier la valeur de SEUIL dans la subroutine USLAG1.    ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 1072 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========                                               ',/,&
'@    L''INDICATEUR DE SEUIL POUR LE CALCUL DES STATISTIQUES  ',/,&
'@       A UNE VALEUR NON PERMISE (LAGOPT).                   ',/,&
'@                                                            ',/,&
'@    IDSTNT DEVRAIT ETRE UN ENTIER SUPERIEUR OU EGAL A 1     ',/,&
'@       IL VAUT ICI IDSTNT = ', I10                           ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier la valeur de IDSTNT dans la subroutine USLAG1.   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 1073 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========                                               ',/,&
'@    L''INDICATEUR DE CALCUL STATIONNAIRE DES STATISTIQUES   ',/,&
'@       A UNE VALEUR NON PERMISE (LAGOPT).                   ',/,&
'@                                                            ',/,&
'@    NSTIST DEVRAIT ETRE UN ENTIER SUPERIEUR OU EGAL         ',/,&
'@       A IDSTNT  = ', I10                                    ,/,&
'@       IL VAUT ICI NSTIST = ', I10                           ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier la valeur de NSTIST dans la subroutine USLAG1.   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 1074 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========                                               ',/,&
'@    L''INDICATEUR SUR LE NOMBRE DE VARIABLES                ',/,&
'@       SUPPLEMENTAIRES LIEES AUX STATISTIQUES               ',/,&
'@       A UNE VALEUR NON PERMISE (LAGOPT).                   ',/,&
'@                                                            ',/,&
'@     NVLSTS DEVRAIT ETRE UN ENTIER ENTRE 0 ET ',I10          ,/,&
'@       IL VAUT ICI NVLSTS = ', I10                           ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier la valeur de NVLSTS dans la subroutine USLAG1.   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

1075  format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========                                               ',/,&
'@    L''INDICATEUR SUR LE NOMBRE DE CLASSE DE STATISTIQUE    ',/,&
'@       A UNE VALEUR NON PERMISE (LAGOPT).                   ',/,&
'@                                                            ',/,&
'@     NBCLST DEVRAIT ETRE UN ENTIER ENTRE 1 ET ',I10          ,/,&
'@       IL VAUT ICI NBCLST = ', I10                           ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier la valeur de NBCLST dans la subroutine USLAG1.   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 1080 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========                                               ',/,&
'@    L''INDICATEUR SUR L''UTILISATION DE L''INJECTION        ',/,&
'@       CONTINUE A UNE VALEUR NON PERMISE (LAGOPT).          ',/,&
'@                                                            ',/,&
'@    INJCON DEVRAIT ETRE UN ENTIER EGAL A 0 OU 1             ',/,&
'@       IL VAUT ICI INJCON = ', I10                           ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier la valeur de INJCON dans la subroutine USLAG1.   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 1090 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========                                               ',/,&
'@    L''INDICATEUR SUR L''UTILISATION DE LA METHODE          ',/,&
'@       DE CLONAGE/FUSION DES PARTICULES                     ',/,&
'@       A UNE VALEUR NON PERMISE (LAGOPT).                   ',/,&
'@                                                            ',/,&
'@    IROULE DEVRAIT ETRE UN ENTIER COMPRIS ENTRE 0 ET 2      ',/,&
'@       IL VAUT ICI IROULE = ', I10                           ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier la valeur de IROULE dans la subroutine USLAG1.   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 2000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========                                               ',/,&
'@    L''INDICATEUR SUR L''ORDRE D''INTEGRATION               ',/,&
'@       DES EQUATIONS DIFFERENTIELLES STOCHASTIQUES          ',/,&
'@       A UNE VALEUR NON PERMISE (LAGOPT).                   ',/,&
'@                                                            ',/,&
'@    NORDRE DEVRAIT ETRE UN ENTIER EGAL A 1                  ',/,&
'@    (L''ORDRE 2 ETANT INOPERANT DANS CETTE VERSION)         ',/,&
'@       IL VAUT ICI NORDRE = ', I10                           ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier la valeur de NORDRE dans la subroutine USLAG1.   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 2010 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========                                               ',/,&
'@    L''INDICATEUR SUR LA PRISE EN COMPTE DE LA DISPERSION   ',/,&
'@       TURBULENTE A UNE VALEUR NON PERMISE (LAGOPT).        ',/,&
'@                                                            ',/,&
'@    IDISTU DEVRAIT ETRE UN ENTIER EGAL A 0 OU 1             ',/,&
'@       IL VAUT ICI IDISTU = ', I10                           ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier la valeur de IDISTU dans la subroutine USLAG1.   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 2011 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========                                               ',/,&
'@    LE MODULE LAGRANGIEN EST INCOMPATIBLE AVEC LE MODELE    ',/,&
'@    DE TURBULENCE SELECTIONNE (LAGOPT).                     ',/,&
'@                                                            ',/,&
'@   Le module lagrangien a ete active avec IILAGR = ',I10     ,/,&
'@     et la dispersion turbulente est prise en compte        ',/,&
'@                                     avec IDISTU = ',I10     ,/,&
'@   Le modele de turbulence                                  ',/,&
'@     correspond a ITURB  = ',I10                             ,/,&
'@   Or, les seuls traitements de la turbulence compatibles   ',/,&
'@     avec le module Lagrangien et la dispersion turbulente  ',/,&
'@     sont k-epsilon et Rij-epsilon, v2f et k-omega          ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier la valeur de IILAGR et IDISTU dans la subroutine ',/,&
'@  USLAG1 et verifier la valeur de ITURB  dans la subroutine ',/,&
'@  USINI1.                                                   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 2012 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========                                               ',/,&
'@    LE MODULE LAGRANGIEN EST INCOMPATIBLE AVEC LE MODELE    ',/,&
'@    DE TURBULENCE SELECTIONNE (LAGOPT).                     ',/,&
'@                                                            ',/,&
'@   Le module lagrangien a ete active avec IILAGR = ',I10     ,/,&
'@     et la dispersion turbulente n''est pas prise en compte ',/,&
'@                                     avec IDISTU = ',I10     ,/,&
'@   Le modele de turbulence                                  ',/,&
'@     correspond a ITURB  = ',I10                             ,/,&
'@   Or, les seuls traitements de la turbulence compatibles   ',/,&
'@     avec le module lagrangien sont : calcul laminaire,     ',/,&
'@     k-epsilon, Rij-epsilon, v2f et k-omega.                ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier la valeur de IILAGR et IDISTU dans la subroutine ',/,&
'@  USLAG1 et verifier la valeur de ITURB  dans la subroutine ',/,&
'@  USINI1.                                                   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 2013 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========                                               ',/,&
'@    L''INDICATEUR SUR LA PRISE EN COMPTE DE LA DIFFUSION    ',/,&
'@       TURBULENTE A UNE VALEUR NON PERMISE (LAGOPT).        ',/,&
'@                                                            ',/,&
'@    IDIFFL DEVRAIT ETRE UN ENTIER EGAL A 0 OU 1             ',/,&
'@       IL VAUT ICI IDIFFL = ', I10                           ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier la valeur de IDIFFL dans la subroutine USLAG1.   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 2014 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========                                               ',/,&
'@    L''INDICATEUR SUR LE CHOIX DU MODELE DE DISPERSION      ',/,&
'@       TURBULENTE A UNE VALEUR NON PERMISE (LAGOPT).        ',/,&
'@                                                            ',/,&
'@    MODCPL DEVRAIT ETRE UN ENTIER SUPERIEUR OU EGAL A 0     ',/,&
'@       IL VAUT ICI MODCPL = ', I10                           ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier la valeur de MODCPL dans la subroutine USLAG1.   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 2015 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========                                               ',/,&
'@    L''INDICATEUR SUR LE CHOIX DU MODELE DE DISPERSION      ',/,&
'@       TURBULENTE EST INCOMPATIBLE AVEC CELUI DU CALCUL     ',/,&
'@       DES STATISTIQUES (LAGOPT).                           ',/,&
'@                                                            ',/,&
'@    LE MODELE COMPLET DE DISPERSION TURBULENTE EST ACTIVE   ',/,&
'@      (MODCPL = ',I10,')'                                    ,/,&
'@      AVANT LE DEBUT DU CALCUL DES STATISTIQUES             ',/,&
'@      (IDSTNT = ',I10,')'                                    ,/,&
'@                                                            ',/,&
'@  Il est necessaire d''avoir calcule des statistiques       ',/,&
'@  pour declencher le modele de dispersion turbulent complet.',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier la valeur de MODCPL dans la subroutine USLAG1.   ',/,&
'@  Verifier la valeur de IDSTNT dans la subroutine USLAG1.   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 2016 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========                                               ',/,&
'@    LE CHOIX DE LA DIRECTION DU MODELE COMPLET              ',/,&
'@       A UNE VALEUR NON PERMISE (LAGOPT).                   ',/,&
'@                                                            ',/,&
'@    IDIRLA DEVRAIT ETRE UN ENTIER EGAL A 1, 2 OU 3          ',/,&
'@       (LA VALEUR 1 POUR UN ECOULEMENT SELON L''AXE X,      ',/,&
'@        LA VALEUR 2 POUR UN ECOULEMENT SELON L''AXE Y,      ',/,&
'@        LA VALEUR 3 POUR UN ECOULEMENT SELON L''AXE Z)      ',/,&
'@       IL VAUT ICI IDIRLA = ', I10                           ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier la valeur de IDIRLA dans la subroutine USLAG1.   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 2017 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========                                               ',/,&
'@    L''INDICATEUR DE RESOLUTION DE L''EQUATION DE POISSON   ',/,&
'@      POUR LES VITESSE MOYENNES ET DE CORRECTION DES        ',/,&
'@      VITESSES INSTANTANNEES DES PARTICULES                 ',/,&
'@      A UNE VALEUR NON PERMISE (LAGOPT).                    ',/,&
'@                                                            ',/,&
'@    ILAPOI DEVRAIT ETRE UN ENTIER EGAL A 0 OU 1             ',/,&
'@       IL VAUT ICI ILAPOI = ', I10                           ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier la valeur de ILAPOI dans la subroutine USLAG1.   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 2018 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========                                               ',/,&
'@    L''INDICATEUR SUR LE CHOIX DU MODELE DE DISPERSION      ',/,&
'@      TURBULENTE A UNE VALEUR NON COMPATIBLE AVEC           ',/,&
'@       CELUI DE L''ENCLENCHEMENT DES STATISTIQUES (LAGOPT). ',/,&
'@                                                            ',/,&
'@    LE MODELE COMPLET DE DISPERSION TURBULENTE EST ACTIVE   ',/,&
'@      (MODCPL = ',I10,')'                                    ,/,&
'@      AVANT L''ENCLENCHEMENT DU CALCUL DES STATISTIQUES     ',/,&
'@      (ISTALA = ',I10,')'                                    ,/,&
'@                                                            ',/,&
'@  Il est necessaire d''activer le calcul des statistiques   ',/,&
'@  pour declencher le modele de dispersion turbulent complet.',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier la valeur de MODCPL dans la subroutine USLAG1.   ',/,&
'@  Verifier la valeur de ISTALA dans la subroutine USLAG1.   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 2030 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========                                               ',/,&
'@    L''INDICATEUR DE POST-PROCESSING EN MODE TRAJECTOIRES   ',/,&
'@      A UNE VALEUR NON PERMISE (LAGOPT).                    ',/,&
'@                                                            ',/,&
'@    IENSI1 DEVRAIT ETRE UN ENTIER EGAL A 0 OU 1             ',/,&
'@       IL VAUT ICI IENSI1 = ', I10                           ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier la valeur de IENSI1 dans la subroutine USLAG1.   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 2031 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========                                               ',/,&
'@    L''INDICATEUR DE POST-PROCESSING EN MODE DEPLACEMENTS   ',/,&
'@      PARTICULAIRES A UNE VALEUR NON PERMISE (LAGOPT).      ',/,&
'@                                                            ',/,&
'@    IENSI2 DEVRAIT ETRE UN ENTIER EGAL A 0 OU 1             ',/,&
'@       IL VAUT ICI IENSI2 = ', I10                           ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier la valeur de IENSI2 dans la subroutine USLAG1.   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 2032 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========                                               ',/,&
'@    LE NOMBRE DE PARTICULE A VISUALISER EN POST-PROCESSING  ',/,&
'@      A UNE VALEUR NON PERMISE (LAGOPT).                    ',/,&
'@                                                            ',/,&
'@    NBVIS DEVRAIT ETRE UN ENTIER SUPERIEUR OU EGAL A 0      ',/,&
'@       INFERIEUR AU NOMBRE MAX DE PARTICULES NBPMAX = ',I10  ,/,&
'@       INFERIEUR AU PARAMETRE NLISTE = ',I10                 ,/,&
'@                                                            ',/,&
'@       IL VAUT ICI NBVIS = ', I10                            ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier la valeur de NBVIS dans la subroutine USLAG1.    ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 2033 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========                                               ',/,&
'@    L''INDICATEUR DE LA FREQUENCE D''ACQUISITION DES DONNES ',/,&
'@       POUR LES SORTIES DE POST-PROCESSING EN               ',/,&
'@       MODE TRAJECTOIRES OU MODE DEPLACEMENTS               ',/,&
'@       A UNE VALEUR NON PERMISE (LAGOPT).                   ',/,&
'@                                                            ',/,&
'@    NVISLA DEVRAIT ETRE UN ENTIER SUPERIEUR OU EGAL A 1     ',/,&
'@       IL VAUT ICI NVISLA = ', I10                           ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier la valeur de NVISLA dans la subroutine USLAG1.   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 2040 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========                                               ',/,&
'@    L''INDICATEUR DE POST-PROCESSING SUR LA VARIABLE        ',/,&
'@       "VITESSE DU FLUIDE VU"                               ',/,&
'@       A UNE VALEUR NON PERMISE (LAGOPT).                   ',/,&
'@                                                            ',/,&
'@    IVISV1 DEVRAIT ETRE UN ENTIER EGAL A 0 OU 1             ',/,&
'@       IL VAUT ICI IVISV1 = ', I10                           ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier la valeur de IVISV1 dans la subroutine USLAG1.   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 2041 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========                                               ',/,&
'@    L''INDICATEUR DE POST-PROCESSING SUR LA VARIABLE        ',/,&
'@       "VITESSE DES PARTICULES"                             ',/,&
'@       A UNE VALEUR NON PERMISE (LAGOPT).                   ',/,&
'@                                                            ',/,&
'@    IVISV2 DEVRAIT ETRE UN ENTIER EGAL A 0 OU 1             ',/,&
'@       IL VAUT ICI IVISV2 = ', I10                           ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier la valeur de IVISV2 dans la subroutine USLAG1.   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 2042 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========                                               ',/,&
'@    L''INDICATEUR DE POST-PROCESSING SUR LA VARIABLE        ',/,&
'@       "TEMPS DE SEJOUR DES PARTICULES"                     ',/,&
'@       A UNE VALEUR NON PERMISE (LAGOPT).                   ',/,&
'@                                                            ',/,&
'@    IVISTP DEVRAIT ETRE UN ENTIER EGAL A 0 OU 1             ',/,&
'@       IL VAUT ICI IVISTP = ', I10                           ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier la valeur de IVISTP dans la subroutine USLAG1.   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 2043 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========                                               ',/,&
'@    L''INDICATEUR DE POST-PROCESSING SUR LA VARIABLE        ',/,&
'@       "DIAMETRE DES PARTICULES"                            ',/,&
'@       A UNE VALEUR NON PERMISE (LAGOPT).                   ',/,&
'@                                                            ',/,&
'@    IVISDM DEVRAIT ETRE UN ENTIER EGAL A 0 OU 1             ',/,&
'@       IL VAUT ICI IVISDM = ', I10                           ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier la valeur de IVISDM dans la subroutine USLAG1.   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 2044 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========                                               ',/,&
'@    L''INDICATEUR DE POST-PROCESSING SUR LA VARIABLE        ',/,&
'@       "TEMPERATURE DES PARTICULES"                         ',/,&
'@       A UNE VALEUR NON PERMISE (LAGOPT).                   ',/,&
'@                                                            ',/,&
'@    IVISTE DEVRAIT ETRE UN ENTIER EGAL A 0 OU 1             ',/,&
'@       IL VAUT ICI IVISTE = ', I10                           ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier la valeur de IVISTE dans la subroutine USLAG1.   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 2045 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========                                               ',/,&
'@    L''INDICATEUR DE POST-PROCESSING SUR LA VARIABLE        ',/,&
'@       "TEMPERATURE OU ENTHALPIE DES PARTICULES DE CHARBON" ',/,&
'@       A UNE VALEUR NON PERMISE (LAGOPT).                   ',/,&
'@                                                            ',/,&
'@    IVISHP DEVRAIT ETRE UN ENTIER EGAL A 0 OU 1             ',/,&
'@       IL VAUT ICI IVISHP = ', I10                           ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier la valeur de IVISHP dans la subroutine USLAG1.   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 2046 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========                                               ',/,&
'@    L''INDICATEUR DE POST-PROCESSING SUR LA VARIABLE        ',/,&
'@       "DIAMETRE DU COEUR RETRECISSANT DES PARTICULES       ',/,&
'@        DE CHARBON" A UNE VALEUR NON PERMISE (LAGOPT).      ',/,&
'@                                                            ',/,&
'@    IVISDK DEVRAIT ETRE UN ENTIER EGAL A 0 OU 1             ',/,&
'@       IL VAUT ICI IVISDK = ', I10                           ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier la valeur de IVISDK dans la subroutine USLAG1.   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 2047 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========                                               ',/,&
'@    L''INDICATEUR DE POST-PROCESSING SUR LA VARIABLE        ',/,&
'@       "MASSE CHARBON REACTIF DES PARTICULES DE CHARBON"    ',/,&
'@       A UNE VALEUR NON PERMISE (LAGOPT).                   ',/,&
'@                                                            ',/,&
'@    IVISCH DEVRAIT ETRE UN ENTIER EGAL A 0 OU 1             ',/,&
'@       IL VAUT ICI IVISCH = ', I10                           ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier la valeur de IVISCH dans la subroutine USLAG1.   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 2048 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========                                               ',/,&
'@    L''INDICATEUR DE POST-PROCESSING SUR LA VARIABLE        ',/,&
'@       "MASSE DE COKE DES PARTICULES DE CHARBON"            ',/,&
'@       A UNE VALEUR NON PERMISE (LAGOPT).                   ',/,&
'@                                                            ',/,&
'@    IVISCK DEVRAIT ETRE UN ENTIER EGAL A 0 OU 1             ',/,&
'@       IL VAUT ICI IVISCK = ', I10                           ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier la valeur de IVISCK dans la subroutine USLAG1.   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 2050 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========                                               ',/,&
'@    L''INDICATEUR SUR L''ACTIVATION DES STATISTIQUES        ',/,&
'@       SUR LES INTERACTIONS PARTICULES/FRONTIERES           ',/,&
'@       A UNE VALEUR NON PERMISE (LAGOPT).                   ',/,&
'@                                                            ',/,&
'@    IENSI3 DEVRAIT ETRE UN ENTIER EGAL A 0 OU 1             ',/,&
'@       IL VAUT ICI IENSI3 = ', I10                           ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier la valeur de IENSI3 dans la subroutine USLAG1.   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 2051 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========                                               ',/,&
'@    L''INDICATEUR DE CALCUL DE LA STATISTIQUE AUX FRONTIERES',/,&
'@       "MASSE DE GRAINS DE CHARBON ENCRASSES"               ',/,&
'@       A UNE VALEUR NON PERMISE (LAGOPT).                   ',/,&
'@                                                            ',/,&
'@    IENCBD DEVRAIT ETRE UN ENTIER EGAL A 0 OU 1             ',/,&
'@       IL VAUT ICI IENCBD = ', I10                           ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier la valeur de IENCBD dans la subroutine USLAG1.   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 2052 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========                                               ',/,&
'@    L''INDICATEUR DE CALCUL DE LA STATISTIQUE AUX FRONTIERES',/,&
'@       "NOMBRE D''INTERACTIONS PARTICULES/FRONTIERES"       ',/,&
'@       A UNE VALEUR NON PERMISE (LAGOPT).                   ',/,&
'@                                                            ',/,&
'@    INBRBD DEVRAIT ETRE UN ENTIER EGAL A 0 OU 1             ',/,&
'@       IL VAUT ICI INBRBD = ', I10                           ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier la valeur de INBRBD dans la subroutine USLAG1.   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 2053 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========                                               ',/,&
'@    L''INDICATEUR DE CALCUL DE LA STATISTIQUE AUX FRONTIERES',/,&
'@       "FLUX DE MASSE PARTICULAIRE LIE AUX INTERACTIONS"    ',/,&
'@       A UNE VALEUR NON PERMISE (LAGOPT).                   ',/,&
'@                                                            ',/,&
'@    IFLMBD DEVRAIT ETRE UN ENTIER EGAL A 0 OU 1             ',/,&
'@       IL VAUT ICI IFLMBD = ', I10                           ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier la valeur de IFLMBD dans la subroutine USLAG1.   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 2054 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========                                               ',/,&
'@    L''INDICATEUR DE CALCUL DE LA STATISTIQUE AUX FRONTIERES',/,&
'@       "ANGLE ENTRE LA VITESSE DE LA PARTICULE ET LE PLAN   ',/,&
'@        DE LA FACE FRONTIERE"                               ',/,&
'@       A UNE VALEUR NON PERMISE (LAGOPT).                   ',/,&
'@                                                            ',/,&
'@    IANGBD DEVRAIT ETRE UN ENTIER EGAL A 0 OU 1             ',/,&
'@       IL VAUT ICI IANGBD = ', I10                           ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier la valeur de IANGBD dans la subroutine USLAG1.   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 2055 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========                                               ',/,&
'@    L''INDICATEUR DE CALCUL DE LA STATISTIQUE AUX FRONTIERES',/,&
'@       "VITESSE DE LA PARTICULE AU MOMENT DE L''INTERACTION"',/,&
'@       A UNE VALEUR NON PERMISE (LAGOPT).                   ',/,&
'@                                                            ',/,&
'@    IVITBD DEVRAIT ETRE UN ENTIER EGAL A 0 OU 1             ',/,&
'@       IL VAUT ICI IVITBD = ', I10                           ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier la valeur de IVITBD dans la subroutine USLAG1.   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 2056 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========                                               ',/,&
'@    LE NOMBRE D''INFORMATIONS SUPPLEMENTAIRES POUR LE       ',/,&
'@       CALCUL DES STATISTIQUES AUX FRONTIERES               ',/,&
'@       A UNE VALEUR NON PERMISE (LAGOPT).                   ',/,&
'@                                                            ',/,&
'@    NUSBOR DEVRAIT ETRE UN ENTIER SUPERIEUR OU EGAL A 0     ',/,&
'@       ET INFERIEUR OU EGAL A NUSBRD = ',I10                 ,/,&
'@       IL VAUT ICI NUSBOR = ', I10                           ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier la valeur de NUSBOR dans la subroutine USLAG1.   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 2057 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========                                               ',/,&
'@    L''INDICATEUR DE CALCUL STATIONNAIRE DES STATISTIQUES   ',/,&
'@       AUX FRONTIERES A UNE VALEUR NON PERMISE (LAGOPT).    ',/,&
'@                                                            ',/,&
'@    NSTBOR DEVRAIT ETRE UN ENTIER SUPERIEUR OU EGAL A 1     ',/,&
'@       IL VAUT ICI NSTBOR = ', I10                           ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier la valeur de NSTBOR dans la subroutine USLAG1.   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 2058 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========                                               ',/,&
'@    L''INDICATEUR DE SEUIL POUR LE CALCUL DES STATISTIQUES  ',/,&
'@       AUX FRONTIERES A UNE VALEUR NON PERMISE (LAGOPT).    ',/,&
'@                                                            ',/,&
'@    SEUILF DEVRAIT ETRE UN REEL SUPERIEUR OU EGAL A 0.D0    ',/,&
'@       IL VAUT ICI SEUILF = ', E14.5                         ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier la valeur de SEUILF dans la subroutine USLAG1.   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 2060 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : A L''EXECUTION DU MODULE LAGRANGIEN         ',/,&
'@    =========                                               ',/,&
'@    LA STATISTIQUE AUX FRONTIERES                           ',/,&
'@     "NOMBRE D''INTERACTIONS PARTICULES/FRONTIERES"         ',/,&
'@     EST ACTIVE AVEC APPLICATION DE LA MOYENNE              ',/,&
'@     PARTICULAIRE (LAGOPT).                                 ',/,&
'@                                                            ',/,&
'@    LES INDICATEURS DE CALCUL DES STATISTIQUES VALENT :     ',/,&
'@       INBRBD       = ',I10                                  ,/,&
'@       IMOYBR(INBR) = ',I10                                  ,/,&
'@                                                            ',/,&
'@  Le calcul continue mais risque de donner des affichages   ',/,&
'@    et des sorties graphiques incoherentes. Une suite       ',/,&
'@    de calcul est possible sans pertes de donnees.          ',/,&
'@                                                            ',/,&
'@  Verifier la valeur de IMOYBR(INBR) dans la  USLAG1.       ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 2061 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : A L''EXECUTION DU MODULE LAGRANGIEN         ',/,&
'@    =========                                               ',/,&
'@    LA STATISTIQUE AUX FRONTIERES :                         ',/,&
'@ ',A60                                                       ,/,&
'@     EST ACTIVE AVEC APPLICATION DE LA MOYENNE              ',/,&
'@     PARTICULAIRE, ALORS QUE LE COMPTE DU NOMBRE            ',/,&
'@     D''INTERACTIONS PARTICULES/FRONTIERES N''EST PAS       ',/,&
'@     ENCLENCHE (LAGOPT).                                    ',/,&
'@                                                            ',/,&
'@    LES INDICATEURS DE CALCUL DES STATISTIQUES VALENT :     ',/,&
'@       ',A12,     ' = ',I10                                  ,/,&
'@             IMOYBR = ',I10                                  ,/,&
'@                                                            ',/,&
'@    INBRBD DEVRAIT ETRE A 1, IL VAUT ICI INBRBD = ',I10      ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier la valeur de INBRBD dans la subroutine USLAG1.   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 2062 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : A L''EXECUTION DU MODULE LAGRANGIEN         ',/,&
'@    =========                                               ',/,&
'@    LA STATISTIQUE AUX FRONTIERES UTILISATEUR :             ',/,&
'@ ',A60                                                       ,/,&
'@     EST ACTIVE AVEC APPLICATION DE LA MOYENNE              ',/,&
'@     PARTICULAIRE, ALORS QUE LE COMPTE DU NOMBRE            ',/,&
'@     D''INTERACTIONS PARTICULES/FRONTIERES N''EST PAS       ',/,&
'@     ENCLENCHE (LAGOPT).                                    ',/,&
'@                                                            ',/,&
'@    LES INDICATEURS DE CALCUL DES STATISTIQUES VALENT :     ',/,&
'@       STAT UTILISATEUR NUMERO ',I10                         ,/,&
'@       IMOYBR                = ',I10                         ,/,&
'@                                                            ',/,&
'@    INBRBD DEVRAIT ETRE A 1, IL VAUT ICI INBRBD = ',I10      ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier la valeur de INBRBD dans la subroutine USLAG1.   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 2063 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========                                               ',/,&
'@    L''INDICATEUR DE MOYENNE POUR LE CALCUL DES STATISTIQUES',/,&
'@       AUX FRONTIERES A UNE VALEUR NON PERMISE (LAGOPT).    ',/,&
'@                                                            ',/,&
'@    IMOYBR DEVRAIT ETRE UN ENTIER EGAL A 0, 1 OU 2          ',/,&
'@       IL VAUT ICI IMOYBR = ', I10                           ,/,&
'@       POUR LA STATISTIQUE AUX FRONTERES :                  ',/,&
'@ ',A60                                                       ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier la valeur de IMOYBR dans la subroutine USLAG1.   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 3001 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========                                               ',/,&
'@    LE NOMBRE DE DONNEES REELLES SUR L''INJECTION           ',/,&
'@     DES PARTICULES DANS LE DOMAINE DE CALCUL EST           ',/,&
'@     SUPERIEUR AU NOMBRE MAXIMAL PREVU PAR DEFAUT.          ',/,&
'@    (LAGOPT).                                               ',/,&
'@                                                            ',/,&
'@    LE NOMBRE DE DONNEES REELLES DEMANDEES VAUT : ',I10      ,/,&
'@                                                            ',/,&
'@    LE NOMBRE DE DONNEES MAXIMAL PREVU PAR DEFAUT VAUT      ',/,&
'@      NDLAGM = ',I10                                         ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 3002 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========                                               ',/,&
'@    LE NOMBRE DE DONNEES ENTIERES SUR L''INJECTION          ',/,&
'@     DES PARTICULES DANS LE DOMAINE DE CALCUL EST           ',/,&
'@     SUPERIEUR AU NOMBRE MAXIMAL PREVU PAR DEFAUT.          ',/,&
'@    (LAGOPT).                                               ',/,&
'@                                                            ',/,&
'@    LE NOMBRE DE DONNEES ENTIERES DEMANDEES VAUT : ',I10     ,/,&
'@                                                            ',/,&
'@    LE NOMBRE DE DONNEES MAXIMAL PREVU PAR DEFAUT VAUT      ',/,&
'@      NDLAIM = ',I10                                         ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 3003 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========                                       (LAGOPT)',/,&
'@    LE NOMBRE DE STATISTIQUES VOLUMIQUES LAGRANGIENNES      ',/,&
'@      NVLSTA = ', I10                                        ,/,&
'@    EST SUPERIEUR AU NOMBRE MAXIMAL PREVU                   ',/,&
'@      NVPLMX = ', I10                                        ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 3004 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========                                       (LAGOPT)',/,&
'@    LE NOMBRE DE POINTEURS DANS LES TABLEAUX DE             ',/,&
'@      VARIABLES PARTICULAIRES ETTP ET ETTPA                 ',/,&
'@      IRF = ', I10                                           ,/,&
'@    EST SUPERIEUR AU NOMBRE MAXIMAL CALCULE                 ',/,&
'@      NVP = ', I10                                           ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 3005 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========                                       (LAGOPT)',/,&
'@    LE NOMBRE DE POINTEURS DANS LE TABLEAU DE               ',/,&
'@      PARAMETRES PARTICULAIRES TEPA                         ',/,&
'@      IRF = ', I10                                           ,/,&
'@    EST DIFFERENT DU NOMBRE MAXIMAL CALCULE                 ',/,&
'@      NVEP = ', I10                                          ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 3006 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========                                       (LAGOPT)',/,&
'@    LE NOMBRE DE POINTEURS DANS LE TABLEAUX DE              ',/,&
'@      PARAMETRES PARTICULAIRES ITEPA                        ',/,&
'@      IRF = ', I10                                           ,/,&
'@    EST SUPERIEUR AU NOMBRE MAXIMAL CALCULE                 ',/,&
'@      NIVEP = ', I10                                         ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 3010 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========                                               ',/,&
'@    LE MODULE LAGRANGIEN EST INCOMPATIBLE AVEC LE MODELE    ',/,&
'@    DE TURBULENCE SELECTIONNE (LAGOPT).                     ',/,&
'@                                                            ',/,&
'@   Le module lagrangien a ete active avec IILAGR = ',I10     ,/,&
'@     et le couplage inverse sur la dynamique est pris en    ',/,&
'@                              compte avec LTSDYN = ',I10     ,/,&
'@   Le modele de turbulence                                  ',/,&
'@     correspond a ITURB  = ',I10                             ,/,&
'@   Or, les seuls traitements de la turbulence compatibles   ',/,&
'@     avec le module Lagrangien et le couplage inverse sur   ',/,&
'@     la dynamique sont k-epsilon, Rij-epsilon, v2f          ',/,&
'@     et k-omega                                             ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier la valeur de IILAGR et IDISTU dans la subroutine ',/,&
'@  USLAG1 et verifier la valeur de ITURB  dans la subroutine ',/,&
'@  USINI1.                                                   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)


 3014 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========                                               ',/,&
'@    LE MODELE SPECIFIQUE DE DEPOT (Guingo & Minier, 2008)   ',/,&
'@     EST ACTIVABLE UNIQUEMENT AVEC  UN SCHEMA D''ORDRE  1   ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier les valeurs de idepst et nordre dans la          ',/,&
'@  subroutine USLAG1.                                        ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)


 3015 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========                                               ',/,&
'@    LE POST TRAITEMENT DES PARTICULES EN MODE TRAJECTOIRE   ',/,&
'@    (IENSI1 = 1) OU DEPLACEMENT (IENSI2 = 1) N''EST PAS     ',/,&
'@    COMPATIBLE AVEC UN CALCUL PARALLELE DANS CETTE VERSION  ',/,&
'@    DE CODE_SATURNE                                         ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)



return

end subroutine
