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

subroutine coini1
!================


!===============================================================================
!  FONCTION  :
!  ---------

!         INIT DES OPTIONS DES VARIABLES POUR
!              POUR LA COMBUSTION
!        FLAMME DE DIFFUSION ET DE PREMELANGE
!   EN COMPLEMENT DE CE QUI A DEJA ETE FAIT DANS USINI1

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
use dimens
use numvar
use optcal
use cstphy
use entsor
use cstnum
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use radiat

!===============================================================================

implicit none

! Local variables

integer          ipp , ii , jj, iok, idirac
integer          isc
double precision wmolme

!===============================================================================
!===============================================================================
! 0. VERIFICATION ISCALT, ISCSTH
!===============================================================================
!     L'utilisateur ne doit pas y avoir touche.

if(iscalt.ne.-1) then
  write(nfecra,1000)iscalt
  call csexit (1)
  !==========
endif
do ii = 1, nscapp
  if(iscsth(iscapp(ii)).ne.-10) then
    write(nfecra,1001)ii,iscapp(ii),iscapp(ii),iscsth(iscapp(ii))
    call csexit (1)
    !==========
  endif
enddo
!===============================================================================
! 1. VARIABLES TRANSPORTEES
!===============================================================================

! 1.1 Definition des scamin et des scamax des variables transportees
! ==================================================================

! --> Flamme de diffusion : chimie 3 points

if ( ippmod(icod3p).ge.0 ) then

! ---- Taux de melange
  scamin(ifm) = 0.d0
  scamax(ifm) = 1.d0

! ---- Variance du taux de melange
!        Type de clipping superieur pour la variance
!        0 pas de clipping, 1 clipping var max de fm, 2 clipping a SCAMAX
 iclvfl(ifp2m) = 1
!        SCAMIN(IFP2M) = 0.D0
!        SCAMAX(IFP2M) = 0.25D0

! ---- Enthalpie
  if ( ippmod(icod3p).eq.1 ) then
    scamin(ihm) = -grand
    scamax(ihm) = +grand
  endif

endif


! --> Flamme de premelange : modele EBU

if ( ippmod(icoebu).ge.0 ) then

! ---- Fraction massique des gaz frais
  scamin(iygfm) = 0.d0
  scamax(iygfm) = 1.d0

! ---- Taux de melange
  if ( ippmod(icoebu).eq.2 .or.                                   &
       ippmod(icoebu).eq.3       ) then
    scamin(ifm) = 0.d0
    scamax(ifm) = 1.d0
  endif

! ---- Enthalpie
  if ( ippmod(icoebu).eq.1 .or.                                   &
       ippmod(icoebu).eq.3      ) then
    scamin(ihm) = -grand
    scamax(ihm) = +grand
  endif

endif


! --> Flamme de premelange : modele LWC

if ( ippmod(icolwc).ge.0 ) then
  scamin(ifm) = 0.d0
  scamax(ifm) = 1.d0

  iclvfl(ifp2m) = 0

  scamin(iyfm) = 0.d0
  scamax(iyfm) = 1.d0

  iclvfl(iyfp2m) = 0

  if ( ippmod(icolwc).ge.2 ) then
    scamin(icoyfp) =-0.25d0
    scamax(icoyfp) = 0.25d0
  endif

endif

! 1.2 Nature des scalaires transportes
! ====================================

do isc = 1, nscapp


! ---- Type de scalaire (0 passif, 1 temperature en K
!                                 -1 temperature en C
!                                  2 enthalpie)
!      La distinction -1/1 sert pour le rayonnement
  iscsth(iscapp(isc)) = 0

enddo


! 1.3 L'utilisation de la variable enthalpie necessite un traitement
!     particulier
! ==================================================================


! ---- On resout en enthalpie avec un CP constant (Cf. cpvarp)
if ( ippmod(icod3p).eq.1 .or.                                     &
     ippmod(icoebu).eq.1 .or.                                     &
     ippmod(icoebu).eq.3 .or.                                     &
     ippmod(icolwc).eq.1 .or.                                     &
     ippmod(icolwc).eq.3 .or.                                     &
     ippmod(icolwc).eq.5     ) then

  iscalt = ihm
  iscsth(ihm) = 2

endif


! 1.4 Donnees physiques ou numeriques propres aux scalaires COMBUSTION
! ====================================================================

do isc = 1, nscapp

  jj = iscapp(isc)

  if ( iscavr(jj).le.0 ) then

! ---- En combustion on considere que la viscosite turbulente domine
!      ON S'INTERDIT DONC LE CALCUL DES FLAMMES LAMINAIRES AVEC Le =/= 1

    visls0(jj) = viscl0

  endif

! ---- Schmidt ou Prandtl turbulent

  sigmas(jj) = 0.7d0

! ---- Coeff dissipation des fluctuations

  rvarfl(jj) = 0.8d0

  ii = isca(iscapp(isc))

! ------ Niveau de detail des impressions pour les variables et
!          donc les scalaires (valeurs 0 ou 1)
!          Si = -10000 non modifie par l'utilisateur -> niveau 1
  if(iwarni(ii).eq.-10000) then
    iwarni(ii) = 1
  endif

! ---- Informations relatives a la resolution des scalaires

!       - Facteur multiplicatif du pas de temps

  cdtvar(ii) = 1.d0

!         - Schema convectif % schema 2ieme ordre
!           = 0 : upwind
!           = 1 : second ordre
  blencv(ii) = 1.d0

!         - Type de schema convetif second ordre (utile si BLENCV > 0)
!           = 0 : Second Order Linear Upwind
!           = 1 : Centre
  ischcv(ii) = 1

!         - Test de pente pour basculer d'un schema centre vers l'upwind
!           = 0 : utilisation automatique du test de pente
!           = 1 : calcul sans test de pente
  isstpc(ii) = 0

!         - Reconstruction des flux de convetion et de diffusion aux faces
!           = 0 : pas de reconstruction
  ircflu(ii) = 1

enddo


! 1.5 Variable courante : nom, sortie chrono, suivi listing, sortie histo

!     Comme pour les autres variables,
!       si l'on n'affecte pas les tableaux suivants,
!       les valeurs par defaut seront utilisees

!     NOMVAR( ) = nom de la variable
!     ICHRVR( ) = sortie chono (oui 1/non 0)
!     ILISVR( ) = suivi listing (oui 1/non 0)
!     IHISVR( ) = sortie historique (nombre de sondes et numeros)
!     si IHISVR(.,1)  = -1 sortie sur toutes les sondes

!     NB : Les 8 premiers caracteres du noms seront repris dans le
!          listing 'developpeur'

! =======================================================================

! --> Flamme de diffusion : chimie 3 points - chmie equilibre

if ( ippmod(icod3p).ge.0 ) then

! ---- Taux de melange
  ipp = ipprtp(isca(ifm))
  NOMVAR(IPP)  = 'Fra_MEL'
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1

! ---- Variance du taux melange
  ipp = ipprtp(isca(ifp2m))
  NOMVAR(IPP)  = 'Var_FrMe'
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1

endif


! --> Flamme de premelange : Modele EBU

if ( ippmod(icoebu).ge.0 ) then

! ---- Fraction massique des gaz frais
  ipp = ipprtp(isca(iygfm))
  NOMVAR(IPP)  = 'Fra_GF'
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1
endif

! ---- Taux de melange
if ( ippmod(icoebu).ge.2 ) then
  ipp = ipprtp(isca(ifm))
  NOMVAR(IPP)  = 'Fra_MEL'
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1
endif


! --> Flamme de premelange : Modeles BML et LWC

if ( ippmod(icolwc).ne.-1 ) then
! --- Taux de melange
  ipp = ipprtp(isca(ifm))
  NOMVAR(IPP) = 'Fra_Mel'
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1
! --- Variance du taux de melange
  ipp = ipprtp(isca(ifp2m))
  NOMVAR(IPP) = 'Var_FMe'
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1
! --- Fraction massique
  ipp = ipprtp(isca(iyfm))
  NOMVAR(IPP) = 'Fra_Mas'
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1
! --- Variance de la fraction massique
  ipp = ipprtp(isca(iyfp2m))
  NOMVAR(IPP) = 'Var_FMa'
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1
! --- Covariance
  if (ippmod(icolwc).ge.2) then
    ipp = ipprtp(isca(icoyfp))
    NOMVAR(IPP) = 'COYF_PP4'
    ichrvr(ipp)  = 1
    ilisvr(ipp)  = 1
    ihisvr(ipp,1)= -1
  endif

endif
! --> Si Transport de l'enthalpie

 if ( ippmod(icod3p).eq.1 .or.                                    &
      ippmod(icoebu).eq.1 .or.                                    &
      ippmod(icoebu).eq.3 .or.                                    &
      ippmod(icolwc).eq.1 .or.                                    &
      ippmod(icolwc).eq.3 .or.                                    &
      ippmod(icolwc).eq.5   ) then
   ipp = ipprtp(isca(ihm))
   NOMVAR(IPP)  = 'Enthalpy'
   ichrvr(ipp)  = 1
   ilisvr(ipp)  = 1
   ihisvr(ipp,1)= -1
  endif


!===============================================================================
! 2. VARIABLES ALGEBRIQUES OU D'ETAT
!===============================================================================

! --> Flamme de diffusion :

if ( ippmod(icod3p).ge.0 ) then
  ipp = ipppro(ipproc(itemp))
  NOMVAR(IPP)   = 'Temperature'
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1
  ipp = ipppro(ipproc(iym(1)))
  NOMVAR(IPP)   = 'YM_Fuel'
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1
  ipp = ipppro(ipproc(iym(2)))
  NOMVAR(IPP)   = 'YM_Oxyd'
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1
  ipp = ipppro(ipproc(iym(3)))
  NOMVAR(IPP)   = 'YM_Prod'
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1
endif

! ---> Physique particuliere : Flamme de premelange (modele EBU)

if ( ippmod(icoebu).ge.0 ) then
  ipp = ipppro(ipproc(itemp))
  NOMVAR(IPP)   = 'Temperature'
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1
  ipp = ipppro(ipproc(iym(1)))
  NOMVAR(IPP)   = 'YM_Fuel'
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1
  ipp = ipppro(ipproc(iym(2)))
  NOMVAR(IPP)   = 'YM_Oxyd'
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1
  ipp = ipppro(ipproc(iym(3)))
  NOMVAR(IPP)   = 'YM_Prod'
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1
endif

! ---> Physique particuliere : Flamme de premelange (modele LWC)

if ( ippmod(icolwc).ge. 0 ) then
  ipp = ipppro(ipproc(itsc))
  NOMVAR(IPP)   = 'Source Term'
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1
  ipp = ipppro(ipproc(itemp))
  NOMVAR(IPP)   = 'Temperature'
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1
  ipp = ipppro(ipproc(iym(1)))
  NOMVAR(IPP)   = 'YM_Fuel'
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1
  ipp = ipppro(ipproc(iym(2)))
  NOMVAR(IPP)   = 'YM_Oxyd'
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1
  ipp = ipppro(ipproc(iym(3)))
  NOMVAR(IPP)   = 'YM_Prod'
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

  do idirac = 1, ndirac
    ipp = ipppro(ipproc(irhol(idirac)))
    WRITE(NOMVAR(IPP),'(A4,I1)') 'RHOL',IDIRAC
    ichrvr(ipp)   = 1
    ilisvr(ipp)   = 1
    ihisvr(ipp,1) = -1

    ipp = ipppro(ipproc(iteml(idirac)))
    WRITE(NOMVAR(IPP),'(A4,I1)') 'TEML',IDIRAC
    ichrvr(ipp)   = 1
    ilisvr(ipp)   = 1
    ihisvr(ipp,1) = -1

    ipp = ipppro(ipproc(ifmel(idirac)))
    WRITE(NOMVAR(IPP),'(A4,I1)') 'FMEL',IDIRAC
    ichrvr(ipp)   = 1
    ilisvr(ipp)   = 1
    ihisvr(ipp,1) = -1

    ipp = ipppro(ipproc(ifmal(idirac)))
    WRITE(NOMVAR(IPP),'(A4,I1)') 'FMAL',IDIRAC
    ichrvr(ipp)   = 1
    ilisvr(ipp)   = 1
    ihisvr(ipp,1) = -1

    ipp = ipppro(ipproc(iampl(idirac)))
    WRITE(NOMVAR(IPP),'(A4,I1)') 'AMPL',IDIRAC
    ichrvr(ipp)   = 1
    ilisvr(ipp)   = 1
    ihisvr(ipp,1) = -1

    ipp = ipppro(ipproc(itscl(idirac)))
    WRITE(NOMVAR(IPP),'(A4,I1)') 'TSCL',IDIRAC
    ichrvr(ipp)   = 1
    ilisvr(ipp)   = 1
    ihisvr(ipp,1) = -1

    ipp = ipppro(ipproc(imaml(idirac)))
    WRITE(NOMVAR(IPP),'(A4,I1)') 'MAML',IDIRAC
    ichrvr(ipp)   = 1
    ilisvr(ipp)   = 1
    ihisvr(ipp,1) = -1
  enddo

endif

! ---> Physique particuliere : Flamme de diffusion (chimie 3 points)
!                              Flamme de premelange (modele EBU 1 et 3)
!                              Flamme de premelange (modele LWC 1, 3 et 5)
!                              AVEC RAYONNEMENT

if ( ( ippmod(icod3p).eq.1 .or.                                   &
       ippmod(icoebu).eq.1 .or. ippmod(icoebu).eq.3 .or.          &
       ippmod(icolwc).eq.1 .or. ippmod(icolwc).eq.3 .or.          &
       ippmod(icolwc).eq.5 )                                      &
      .and. (iirayo.ge.1) ) then
  ipp = ipppro(ipproc(ickabs))
  NOMVAR(IPP)   = 'KABS'
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1
  ipp = ipppro(ipproc(it4m))
  NOMVAR(IPP)   = 'TEMP4'
  ichrvr(ipp)   = 0
  ilisvr(ipp)   = 0
  ihisvr(ipp,1) = -1
  ipp = ipppro(ipproc(it3m))
  NOMVAR(IPP)   = 'TEMP3'
  ichrvr(ipp)   = 0
  ilisvr(ipp)   = 0
  ihisvr(ipp,1) = -1
endif


!===============================================================================
! 3. INFORMATIONS COMPLEMENTAIRES
!===============================================================================

! --> Calcul de RO0 a partir de T0 et P0

 if ( ippmod(icod3p).ne.-1 .or.                                   &
      ippmod(icoebu).ne.-1 .or.                                   &
      ippmod(icolwc).ne.-1     ) then
   wmolme = wmolg(2)
   ro0 = p0*wmolme / (rr*t0)
 endif

! On met les constantes a -GRAND pour obliger l'utilisateur a les definir
!  (ou les laisser) dans usebu1, usd3p1 ou uslwc1.
! --> Constante modele EBU par defaut
cebu   =-grand

! --> Constantes modele LWC par defaut
vref  =-grand
lref  =-grand
ta    =-grand
tstar =-grand

! --> Coefficient de relaxation de la masse volumique
!      RHO(n+1) = SRROM * RHO(n) + (1-SRROM) * RHO(n+1)
srrom =-grand

! --> Viscosite laminaire associee au scalaire enthalpie
!       DIFTL0 (diffusivite dynamique en kg/(m s))
diftl0      =-grand

! --> Diffusion 3 points, tableaux HH HF THF
!           (generes dans d3pphy)

nmaxf = 9
nmaxh = 9

! ---> Masse volumique variable et viscosite constante (pour les suites)
irovar = 1
ivivar = 0

!===============================================================================
! 4. ON REDONNE LA MAIN A L'UTLISATEUR
!===============================================================================

if ( ippmod(icoebu).ge.0 ) then
  call usebu1
  !==========
else if( ippmod(icod3p).ge.0 ) then
  call usd3p1
  !==========
else if( ippmod(icolwc).ge.0 ) then
  call uslwc1
!       ==========
endif

!===============================================================================
! 5. VERIFICATION DES DONNERS FOURNIES PAR L'UTLISATEUR
!===============================================================================

iok = 0
if ( ippmod(icoebu).ge.0 ) then

  call ebuver (iok)
  !==========
  if(iok.gt.0) then
    write(nfecra,9999)iok
    call csexit (1)
    !==========
  else
    write(nfecra,9998)
  endif

else if( ippmod(icod3p).ge.0 ) then
  call d3pver (iok)
  !==========
  if(iok.gt.0) then
    write(nfecra,9991)iok
    call csexit (1)
  else
    write(nfecra,9990)
  endif

else if( ippmod(icolwc).ge.0 ) then
  call lwcver (iok)
  !==========
  if(iok.gt.0) then
    write(nfecra,9993)iok
    call csexit (1)
  else
    write(nfecra,9992)
  endif

endif

 1000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    PHYSIQUE PARTICULIERE (COMBUSTION) DEMANDEE             ',/,&
'@                                                            ',/,&
'@  La valeur de ISCALT est renseignee automatiquement.       ',/,&
'@                                                            ',/,&
'@  L''utilisateur ne doit pas la renseigner dans usini1, or  ',/,&
'@    elle a ete affectee comme suit :                        ',/,&
'@    ISCALT = ',I10                                           ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier usini1.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 1001 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    PHYSIQUE PARTICULIERE (COMBUSTION) DEMANDEE             ',/,&
'@                                                            ',/,&
'@  Les valeurs de ISCSTH sont renseignees automatiquement.   ',/,&
'@                                                            ',/,&
'@  L''utilisateur ne doit pas les renseigner dans usini1, or ',/,&
'@    pour le scalaire ',I10   ,' correspondant au scalaire   ',/,&
'@    physique particuliere ',I10   ,' on a                   ',/,&
'@    ISCSTH(',I10   ,') = ',I10                               ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier usini1.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9998 format(                                                           &
'                                                             ',/,&
' Pas d erreur detectee lors de la verification des donnees   ',/,&
'                                                    (usebu1).',/)
 9999 format(                                                           &
'@                                                            ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    LES PARAMETRES DE CALCUL SONT INCOHERENTS OU INCOMPLETS ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute (',I10,' erreurs).          ',/,&
'@                                                            ',/,&
'@  Se reporter aux impressions precedentes pour plus de      ',/,&
'@    renseignements.                                         ',/,&
'@  Verifier usebu1.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9990 format(                                                           &
'                                                             ',/,&
' Pas d erreur detectee lors de la verification des donnees   ',/,&
'                                                    (usd3p1).',/)
 9991 format(                                                           &
'@                                                            ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    LES PARAMETRES DE CALCUL SONT INCOHERENTS OU INCOMPLETS ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute (',I10,' erreurs).          ',/,&
'@                                                            ',/,&
'@  Se reporter aux impressions precedentes pour plus de      ',/,&
'@    renseignements.                                         ',/,&
'@  Verifier usd3p1.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9992 format(                                                           &
'                                                             ',/,&
' Pas d erreur detectee lors de la verification des donnees   ',/,&
'                                                    (uslwc1).',/)
 9993 format(                                                           &
'@                                                            ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    LES PARAMETRES DE CALCUL SONT INCOHERENTS OU INCOMPLETS ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute (',I10,' erreurs).          ',/,&
'@                                                            ',/,&
'@  Se reporter aux impressions precedentes pour plus de      ',/,&
'@    renseignements.                                         ',/,&
'@  Verifier uslwc1.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)


return
end subroutine
