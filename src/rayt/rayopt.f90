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

subroutine rayopt
!================

!===============================================================================
!  FONCTION  :
!  ---------

!   SOUS-PROGRAMME DU MODULE RAYONNEMENT :
!   --------------------------------------

!  1) Initialisation par defaut du parametrage du module de
!     transferts thermiques radiatifs
!  2) Lecture du parametrage utilisateur
!  3) Controle de coherence avec les physiques particulieres
!  4) Verifications du parametrage utilisateur

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
use entsor
use optcal
use cstphy
use ihmpre
use ppppar
use ppthch
use cpincl
use cs_fuel_incl
use ppincl
use radiat

!===============================================================================

implicit none

! Local variables

integer          ii, iok , iiscal, iscaok, ipp, nmodpp
integer          irphas
character        car4*4
character*3      num

!===============================================================================
! 0. REDEFINITION DU NOMBRE DE PHASES POUR LE CHARBON PULVERISE
!===============================================================================

!--> nrphas: for pulverized coal and fuel combustion:
!            nrphas = 1 (gaz) + number of classes (particles or droplets)

!--> For pulverized coal and fuel combustion:
if ( ippmod(iccoal) .ge. 0 ) then
  nrphas = 1 + nclacp
else if ( ippmod(icfuel) .ge. 0 ) then
  nrphas = 1 + nclafu
else
  nrphas = 1
endif

nmodpp = 0
do ipp = 2, nmodmx
  if (ippmod(ipp).ne.-1) then
    nmodpp = nmodpp+1
  endif
enddo

!===============================================================================
! 1. INITIALISATIONS PAR DEFAUT DU MODULE DE TRANSFERTS RADIATIFS
!                        ^^^^^^
!===============================================================================

!-->  IIRAYO = 0 : PAS DE TRANSFERTS RADIATIFS
!            = 1 : TRANSFERTS RADIATIFS, METHODE DES ORDONNEES DISCRETES
!            = 2 : TRANSFERTS RADIATIFS, APPROXIMATION P-1
!     On initialise a -1 pour montrer que ce n'est pas initialise ...
!        (on fera un test apres usray1)
iirayo = -1

!-->  CALCUL DU COEFFICIENT D'ABSORPTION
!      IMODAK = 0 : sans utiliser modak
!               1 : a l'aide modak

imodak = 0

!-->  INDICATEUR SUITE DE CALCUL (LECTURE DU FICHIER SUITE)

isuird = -1

!-->  FREQUENCE DE PASSAGE DANS LE MODULE RAYONNEMENT

nfreqr = -1

!-->  NUMERO DE LA QUADRATURE ET PARAMETRE DE LA Tn

i_quadrature = 1
ndirec = -1

!-->  POURCENTAGE DE CELLULES OU L'ON ADMET QUE LA LONGUEUR OPTIQUE DEPASSE
!       L'UNITE POUR LE MODELE P-1

xnp1mx = 10.d0

!-->  INITIALISATION DU MODE DE CALCUL DU TERME SOURCE RADIATIF EXPLICITE
!     IDIVER = 0 => CALCUL SEMI-ANALYTIQUE (OBLIGATOIRE SI TRANSPARENT)
!     IDIVER = 1 => CALCUL CONSERVATIF
!     IDIVER = 2 => CALCUL SEMI-ANALYTIQUE CORRIGE POUR ETRE CONSERVATIF
!     REMARQUE : SI TRANSPARENT IDIVER = -1 AUTOMATIQUEMENT DANS RAYDOM

idiver = -1

!--> NIVEAU D'AFFICHAGE (0,1,2) DES RENSEIGNEMENTS TEMPERATURE DE PAROI

iimpar = -1

!--> NIVEAU D'AFFICHAGE (0,1,2) DES RENSEIGNEMENTS RESOLUTION LUMINANCE

iimlum = -1

!   - Interface Code_Saturne
!     ======================

if (iihmpr.eq.1) then

  call uiray1(iirayo, isuird, i_quadrature, ndirec, nfreqr, idiver, iimpar, iimlum)
!  ==========

endif

call usray1
!==========

!===============================================================================
! 2. VERIFICATION LA COHERENCE D'UTILISATION DU MODULE DE RAYONNEMENT
!    AVEC LA THERMIQUE OU LES PHYSIQUES PARTICULIERES (COMBUSTION)
!===============================================================================

iok = 0

!--> IIRAYO = 0 (pas de rayonnement).

if(iirayo.eq.-1) then
  iirayo = 0
endif

if (iirayo.ne.0 .and. iirayo.ne.1 .and. iirayo.ne.2) then
  write(nfecra,1010) iirayo
  iok = iok + 1
endif

if (imodak.ne.0 .and. imodak.ne.1) then
  write(nfecra,1020) imodak
  iok = iok + 1
endif

!--> ISCSTH

!     Si physique particuliere avec fichier parametrique,
!       ISCSTH a ete renseigne
!       dans ppini1 ou coini1 (a verifier en elec).
!     Si physique classique et ISCSTH pas modifie dans USRAY1
!                           ou pas de variable thermique
!       STOP.

if (ippmod(iphpar).ge.2) then

  ! Il y a une seule phase ; si on rayonne
  if (iirayo.eq.1 .or. iirayo.eq.2) then
    if (iscalt.le.0) then
      write(nfecra,3001)
      iok = iok + 1
    else if (itherm.ne.2) then
      write(nfecra,3000) itherm
      iok = iok + 1
    endif
  endif

else

  ! Pour la phase qui rayonne
  if (iirayo.eq.1 .or. iirayo.eq.2) then

      ! On cherche s'il y a un scalaire thermique
      iscaok = 0
      do iiscal = 1, nscal
        if (iiscal.eq.iscalt) then
          iscaok = 1

          !           Et on regarde si on a dit temp C, K ou enthalpie
          if (itherm.ne.1 .and. itherm.ne.2) then
            write(nfecra,3010) iiscal,iiscal
            iok = iok + 1
          endif

        endif
      enddo
      if(iscaok.eq.0)then
        write(nfecra,3011)
        iok = iok + 1
      endif

    endif

endif

!--> Stop si erreur.

if(iok.ne.0) then
  call csexit (1)
  !==========
endif

 1010 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    IIRAYO NE PEUT PRENDRE POUR VALEURS QUE 0 1 OU 2        ',/,&
'@    IIRAYO vaut ',I10                                        ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Arret dans rayopt.                                        ',/,&
'@  Verifier usray1 ou l''interface graphique.                ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 1020 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    IMODAK NE PEUT PRENDRE POUR VALEURS QUE 0 OU 1          ',/,&
'@    IMODAK vaut',I10                                         ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Arret dans rayopt.                                        ',/,&
'@  Verifier usray1 ou l''interface graphique.                ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 3000 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : RAYONNEMENT : ARRET A L''ENTREE DES DONNEES ',/,&
'@    =========                                               ',/,&
'@    PHYSIQUE PARTICULIERE ACTIVEE : ENTHALPIE NECESSAIRE    ',/,&
'@                                                            ',/,&
'@  Avec rayonnement, il faut                                 ',/,&
'@    preciser le modele thermique en renseignant ITHERM      ',/,&
'@    soit :                                                  ',/,&
'@                1 temperature                               ',/,&
'@                2 enthalpie                                 ',/,&
'@                                                            ',/,&
'@  Ici, ITHERM = ',I10,'                                     ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Avec physique particuliere, cette initialisation aurait   ',/,&
'@    du etre automatique.                           ~~~~~~   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 3001 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : RAYONNEMENT : ARRET A L''ENTREE DES DONNEES ',/,&
'@    =========                                               ',/,&
'@    PHYSIQUE PARTICULIERE ACTIVEE : ENTHALPIE NECESSAIRE    ',/,&
'@                                                            ',/,&
'@  Lorsque le rayonnement est utilise, il                    ',/,&
'@    faut indiquer qu''un scalaire represente la variable    ',/,&
'@    energetique (enthalpie) en renseignant                  ',/,&
'@    ISCALT(',I10   ,') (numero du scalaire).                ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Cette initialisation aurait du etre automatique.          ',/,&
'@                       ~~~~~~                               ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 3010 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : RAYONNEMENT : ARRET A L''ENTREE DES DONNEES ',/,&
'@    =========                                               ',/,&
'@    ISCSTH DOIT ETRE RENSEIGNE OBLIGATOIREMENT              ',/,&
'@                                                            ',/,&
'@  Avec rayonnement, il faut                                 ',/,&
'@    preciser la variable energetique representee par le     ',/,&
'@    scalaire ',I10   ,' en renseignant ISCSTH(',I10   ,')   ',/,&
'@    soit :                                                  ',/,&
'@               -1 temperature en C                          ',/,&
'@                1 temperature en K                          ',/,&
'@                2 enthalpie                                 ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier les parametres.                                  ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 3011 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : RAYONNEMENT : ARRET A L''ENTREE DES DONNEES ',/,&
'@    =========                                               ',/,&
'@    IL FAUT UTILISER UNE VARIABLE ENERGETIQUE.              ',/,&
'@                                                            ',/,&
'@  Lorsque le rayonnement, il                                ',/,&
'@    faut indiquer qu''un scalaire represente la variable    ',/,&
'@    energetique (temperature ou enthalpie) en renseignant   ',/,&
'@    ISCALT(',I10   ,') (numero du scalaire).               ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier les parametres.                                  ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

!===============================================================================
! 3. VERIFICATIONS (uniquement s'il y a du rayonnement)
!===============================================================================

iok = 0

if (iirayo.gt.0) then

! Positionnement des pointeurs

  call rayprp
  !==========

! --> ISUIRD

  if (isuird.ne.0 .and. isuird.ne.1 ) then
    write(nfecra,4000) isuird
    iok = iok + 1
  endif

! --> NFREQR

  if (nfreqr.le.0) then
    write(nfecra,4010) nfreqr
    iok = iok + 1
  endif

! --> i_quadrature
!     Selection de la quadrature
  if (iirayo.eq.1) then
    if (i_quadrature.lt.1 .or. i_quadrature.gt.8) then
      write(nfecra, 4015) i_quadrature
      iok = iok + 1
    endif
  endif

! --> NDIREC
!     Parametre Quadrature Tn
  if (iirayo.eq.1 .and. i_quadrature.eq.6) then
    if (ndirec.lt.3) then
      write(nfecra, 4020) ndirec
      iok = iok + 1
    endif
  endif

! --> IDIVER
!     Choix entre 0  1 et 2
  if (idiver.ne.0 .and. idiver.ne.1 .and. idiver.ne.2) then
    write(nfecra,4030) idiver
    iok = iok + 1
  endif

! --> IIMPAR
!     Choix entre 0  1 et 2
  if (iimpar.ne.0 .and. iimpar.ne.1 .and. iimpar.ne.2) then
    write(nfecra,4040) iimpar
    iok = iok + 1
  endif

! --> IIMLUM
!     Choix entre 0  1 et 2
  if (iimlum.ne.0 .and. iimlum.ne.1 .and. iimlum.ne.2) then
    write(nfecra,4050) iimlum
    iok = iok + 1
  endif

else
  return
endif

!--> Stop si erreur.

if(iok.ne.0) then
  call csexit (1)
  !==========
endif

!===============================================================================
! 4. Quadrature initialization
!===============================================================================

call raydir
!==========

!===============================================================================
!  5. INITIALISATIONS UTILISATEURS
!                    ^^^^^^^^^^^^
!===============================================================================

!   - Code_Saturne GUI
!     ================

if (iihmpr.eq.1) then

  call uiray4(iirayo)
  !==========

  ! properties on cells
  do ii = 1, nproce
    ipp = ipppro(ii)
    call fcnmva (nomprp(ii), len(nomprp(ii)), ipp)
    !==========
  enddo

  call csenso                                                     &
  !==========
     ( nvppmx, ncapt,  nthist, frhist, ntlist, iecaux,            &
       ipstdv, ichrvr, ilisvr, ihisvr, tplfmt, isca, iscapp,      &
       ipprtp, xyzcap )

  do ii = 1, nproce
    ipp = ipppro(ii)
    call cfnmva(nomprp(ii), len(nomprp(ii)), ipp)
    !==========
  enddo

  call nvamem
  !==========

  ! take into acount user modifications
  call usipes(nmodpp)
  !==========

endif

call usray1
!==========

!--> Stop si erreur.

if(iok.ne.0) then
  call csexit (1)
  !==========
endif

 4000 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : RAYONNEMENT : ARRET A L''ENTREE DES DONNEES ',/,&
'@    =========                                               ',/,&
'@    INDICATEUR DE SUITE DE CALCUL NON ADMISSIBLE            ',/,&
'@                                                            ',/,&
'@  L''indicateur de suite de calcul doit etre 0 ou 1 (ISUIRD)',/,&
'@    Il vaut ici ISUIRD = ',I10                               ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier usray1.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 4010 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : RAYONNEMENT : ARRET A L''ENTREE DES DONNEES ',/,&
'@    =========                                               ',/,&
'@     FREQUENCE DE PASSAGE DANS LE MODULE DE RAYONNEMENT     ',/,&
'@     NON ADMISSIBLE                                         ',/,&
'@                                                            ',/,&
'@  La frequence de passage doit etre superieure ou egale a 1 ',/,&
'@    Elle vaut ici NFREQR = ',I10                             ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier usray1.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 4015 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : RAYONNEMENT : ARRET A L''ENTREE DES DONNEES ',/,&
'@    =========                                               ',/,&
'@    LE NUMERO DE LA QUADRATURE DOIT ETRE ENTRE 1 ET 6 :     ',/,&
'@        S4 = 1                                              ',/,&
'@        S6 = 2                                              ',/,&
'@        S8 = 3                                              ',/,&
'@        T2 = 4                                              ',/,&
'@        T4 = 5                                              ',/,&
'@        Tn = 6                                              ',/,&
'@                                                            ',/,&
'@    Il vaut ici i_quadrature = ',I10                         ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier usray1.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 4020 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : RAYONNEMENT : ARRET A L''ENTREE DES DONNEES ',/,&
'@    =========                                               ',/,&
'@    LE PARAMETRE n DE LA QUADRATURE Tn DOIT ETRE >= 3       ',/,&
'@                                                            ',/,&
'@    Il vaut ici ndirec = ',I10                               ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier usray1.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 4030 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : RAYONNEMENT : ARRET A L''ENTREE DES DONNEES ',/,&
'@    =========                                               ',/,&
'@    INDICATEUR DU MODE DE CALCUL DU TERME SOURCE RADIATIF   ',/,&
'@    EXPLICITE NON ADMISSIBLE                                ',/,&
'@                                                            ',/,&
'@  L''indicateur du mode de calcul doit etre 0, 1 ou 2       ',/,&
'@    Il vaut ici IDIVER = ',I10                               ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier usray1.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 4040 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : RAYONNEMENT : ERREUR A L''ENTREE DES DONNEES',/,&
'@    =========                                               ',/,&
'@    NIVEAU D''AFFICHAGE DES RENSIGNEMENTS DES               ',/,&
'@    TEMPERATURE DE PAROI NON ADMISSIBLE                     ',/,&
'@                                                            ',/,&
'@  Le niveau d''affichage doit etre 0, 1 ou 2  (IIMPAR)      ',/,&
'@    Il vaut ici IIMPAR = ',I10                               ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier usray1.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 4050 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : RAYONNEMENT : ERREUR A L''ENTREE DES DONNEES',/,&
'@    =========                                               ',/,&
'@    NIVEAU D''AFFICHAGE DES RENSIGNEMENTS SUR LA            ',/,&
'@    RESOLUTION DE LA LUMINANCE NON ADMISSIBLE               ',/,&
'@                                                            ',/,&
'@  Le niveau d''affichage doit etre 0, 1 ou 2  (IIMLUM)      ',/,&
'@    Il vaut ici IIMLUM = ',I10                               ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier usray1.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)


return

end subroutine
