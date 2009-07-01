!-------------------------------------------------------------------------------

!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2008 EDF S.A., France

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
!    nom           !type!mode !                   role                         !
!__________________!____!_____!________________________________________________!
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

implicit none

!===============================================================================
!     DONNEES EN COMMON
!===============================================================================

include "paramx.h"
include "dimens.h"
include "numvar.h"
include "entsor.h"
include "optcal.h"
include "cstphy.h"
include "radiat.h"
include "ihmpre.h"
include "ppppar.h"
include "ppthch.h"
include "cpincl.h"
include "fuincl.h"
include "ppincl.h"

!===============================================================================

! VARIABLES LOCALES

integer          ii, jj, iphas, iphaok, iok , iiscal, iscaok
character        car4*4

!===============================================================================
!===============================================================================
! 1. INITIALISATIONS PAR DEFAUT DU MODULE DE TRANSFERTS RADIATIFS
!                        ^^^^^^
!===============================================================================


!-->  IIRAYO = 0 : PAS DE TRANSFERTS RADIATIFS
!            = 1 : TRANSFERTS RADIATIFS

iirayo = 0

!-->  POUR CHAQUE PHASE:
!     IRAYON = 0 : PAS DE TRANSFERTS RADIATIFS
!            = 1 : TRANSFERTS RADIATIFS, METHODE DES ORDONNEES DISCRETES
!            = 2 : TRANSFERTS RADIATIFS, APPROXIMATION P-1
!     On initialise a -1 pour montrer que ce n'est pas initialise ...
!        (on fera un test apres usray1)
do iphas = 1,nphsmx
  irayon(iphas) = -1
enddo


!-->  CALCUL DU COEFFICIENT D'ABSORPTION
!     (REMPLI AUTOMATIQUEMENT ENSUITE POUR LA PHYSIQUE PARTICULIERE)
!      IMODAK = 0 : sans utiliser modak
!               1 : a l'aide modak

imodak = 0

!-->  INDICATEUR SUITE DE CALCUL (LECTURE DU FICHIER SUITE)

isuird = -1

!-->  FREQUENCE DE PASSAGE DANS LE MODULE RAYONNEMENT

nfreqr = -1

!-->  NOMBRE DE DIRECTIONS : 32 OU 128

ndirec = -1

!-->  NOMBRE DE BANDES SPECTRALES (PAS UTILISE)

nbande = 1

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

!--> INITIALISATION DES DONNEES POST-PROCESSING
!    NBRAYP : NOMBRE MAX DES SORTIES DE VARIABLES CELLULES
!    NBRAYF : NOMBRE MAX DES SORTIES DE VARIABLES FACETTES DE BORD

do iphas = 1,nphsmx
  do ii = 1,nbrayp
    WRITE(CAR4,'(I4.4)') II + NBRAYP * (IPHAS -1)
    NBRVAP(II,IPHAS) = 'RAYTCE'//CAR4
    irayvp(ii,iphas) = -1
  enddo
  do ii = 1,nbrayf
    WRITE(CAR4,'(I4.4)') II+ NBRAYF * (IPHAS -1)
    NBRVAF(II,IPHAS) = 'RAYTFB'//CAR4
    irayvf(ii,iphas) = -1
  enddo
enddo

!===============================================================================
! 2. INITIALISATIONS UTILISATEURS
!                    ^^^^^^^^^^^^
!===============================================================================

!   - Interface Code_Saturne
!     ======================

if (iihmpr.eq.1) then

  do iphas = 1, nphas
    do jj = 1,nbrayp
      ii = jj + nbrayp * (iphas -1)
      call fcnmra(nbrvap(jj,iphas), len(nbrvap(jj,iphas)), ii)
    enddo
  enddo
  do iphas = 1, nphas
    do jj = 1, nbrayf
      ii = jj + nphas * nbrayp + nbrayf* (iphas -1)
      call fcnmra(nbrvaf(jj,iphas), len(nbrvaf(jj,iphas)), ii)
    enddo
  enddo

  call uiray1                                                     &
  !==========
     ( nbrayp, nbrayf, nphas,                                     &
       irayon, isuird, ndirec, nfreqr, idiver, iimpar, iimlum,    &
       irayvp, irayvf )

  do iphas = 1, nphas
    do jj = 1,nbrayp
      ii = jj + nbrayp * (iphas -1)
      call cfnmra(nbrvap(jj,iphas), len(nbrvap(jj,iphas)), ii)
    enddo
  enddo
  do iphas = 1, nphas
    do jj = 1, nbrayf
      ii = jj + nphas * nbrayp + nbrayf* (iphas -1)
      call cfnmra(nbrvaf(jj,iphas), len(nbrvaf(jj,iphas)), ii)
    enddo
  enddo

endif

call usray1
!==========

!===============================================================================
! 3. VERIFICATION LA COHERENCE D'UTILISATION DU MODULE DE RAYONNEMENT
!    AVEC LA THERMIQUE OU LES PHYSIQUES PARTICULIERES (COMBUSTION)
!===============================================================================

iok = 0

!--> IRAYON

!     Si physique classique ou particuliere sans fichier paramétrique
!       et IRAYON pas modifie dans USRAY1
!       IRAYON(IPHAS) = 0 (pas de rayonnement).

if( ippmod(iphpar).le.1 ) then

  do iphas = 1, nphas
    if(irayon(iphas).eq.-1) then
      irayon(iphas) = 0
    endif
  enddo


!     Sinon, IRAYON doit etre issu du fichier
!       parametrique ; IMODAK aussi.

else

  iphas = 1

  if (irayon(iphas).eq.-1) then
    if (iraypp.eq.0) then
      irayon(iphas) = 0
      imodak        = 0
    elseif (iraypp.eq.1) then
      irayon(iphas) = 1
      imodak        = 0
    elseif (iraypp.eq.2) then
      irayon(iphas) = 1
      imodak        = 1
    elseif (iraypp.eq.3) then
      irayon(iphas) = 2
      imodak        = 0
    elseif (iraypp.eq.4) then
      irayon(iphas) = 2
      imodak        = 1
    else
      write(nfecra,1010)
      iok = iok + 1
    endif

  else
    write(nfecra,2000)
    iok = iok + 1
  endif

endif

!--> ISCSTH

!     Si physique particuliere avec fichier parametrique,
!       ISCSTH a ete renseigne
!       dans ppini1 ou coini1 (a verifier en elec).
!     Si physique classique et ISCSTH pas modifie dans USRAY1
!                           ou pas de variable thermique
!       STOP.

if(ippmod(iphpar).ge.2) then

!     Il y a une seule phase ; si on rayonne
  do iphas = 1, nphas
    if (irayon(iphas).eq.1 .or. irayon(iphas).eq.2) then
!        On cherche s'il y a un scalaire thermique
      iscaok = 0
      do iiscal = 1, nscal
        if (iiscal.eq.iscalt(iphas)) then
          iscaok = 1
!           Et on regarde si on a dit enthalpie
          if (iscsth(iiscal).ne.2) then
            write(nfecra,3000) iphas,iiscal,iiscal
            iok = iok + 1
          endif
        endif
      enddo
      if (iscaok.eq.0) then
        write(nfecra,3001) iphas,iphas
        iok = iok + 1
      endif
    endif
  enddo

else

!     Pour toutes les phases qui rayonnent
  do iphas = 1, nphas
    if (irayon(iphas).eq.1 .or. irayon(iphas).eq.2) then

!        On cherche s'il y a un scalaire thermique
      iscaok = 0
      do iiscal = 1, nscal
        if (iiscal.eq.iscalt(iphas)) then
          iscaok = 1

!           Et on regarde si on a dit temp C, K ou enthalpie
          if (abs(iscsth(iiscal)).ne.1.and.                       &
                  iscsth(iiscal) .ne.2      ) then
            write(nfecra,3010) iphas,iiscal,iiscal
            iok = iok + 1
          endif

        endif
      enddo
      if(iscaok.eq.0)then
        write(nfecra,3011)iphas,iphas
        iok = iok + 1
      endif

    endif
  enddo

endif

!--> Stop si erreur.

if(iok.ne.0) then
  call csexit (1)
  !==========
endif

 1010 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    IRAYPP NE PEUT PRENDRE POUR VALEURS QUE 0 1 2 3 OU 4    ',/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Arret dans rayopt.                                        ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    EN PHYSIQUE PARTICULIERE IRAYON EST DONNE PAR FICHIER   ',/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier usray1 (IRAYON ne doit pas etre modifie).        ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 3000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : RAYONNEMENT : ARRET A L''ENTREE DES DONNEES ',/,&
'@    =========                                               ',/,&
'@    PHYSIQUE PARTICULIERE ACTIVEE : ENTHALPIE NECESSAIRE    ',/,&
'@                                                            ',/,&
'@  Avec rayonnement, pour la phase ',I10   ,', il faut       ',/,&
'@    preciser la variable energetique representee par le     ',/,&
'@    scalaire ',I10   ,' en renseignant ISCSTH(',I10   ,')   ',/,&
'@    dans usini1 : soit                                      ',/,&
'@               -1 temperature en C                          ',/,&
'@                1 temperature en K                          ',/,&
'@                2 enthalpie                                 ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Avec physique particuliere, cette initialisation aurait   ',/,&
'@    du etre automatique.                           ~~~~~~   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 3001 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : RAYONNEMENT : ARRET A L''ENTREE DES DONNEES ',/,&
'@    =========                                               ',/,&
'@    PHYSIQUE PARTICULIERE ACTIVEE : ENTHALPIE NECESSAIRE    ',/,&
'@                                                            ',/,&
'@  Lorsque le rayonnement est utilise (phase ',I10   ,', il  ',/,&
'@    faut indiquer qu''un scalaire represente la variable    ',/,&
'@    energetique (enthalpie) en renseignant                  ',/,&
'@    ISCALT(',I10   ,') dans usini1 (numero du scalaire).    ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Avec physique particuliere, cette initialisation aurait   ',/,&
'@    du etre automatique.                           ~~~~~~   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 3010 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : RAYONNEMENT : ARRET A L''ENTREE DES DONNEES ',/,&
'@    =========                                               ',/,&
'@    ISCSTH DOIT ETRE RENSEIGNE OBLIGATOIREMENT DANS USINI1  ',/,&
'@                                                            ',/,&
'@  Avec rayonnement, pour la phase ',I10   ,', il faut       ',/,&
'@    preciser la variable energetique representee par le     ',/,&
'@    scalaire ',I10   ,' en renseignant ISCSTH(',I10   ,')   ',/,&
'@    dans usini1 : soit                                      ',/,&
'@               -1 temperature en C                          ',/,&
'@                1 temperature en K                          ',/,&
'@                2 enthalpie                                 ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier usini1.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 3011 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : RAYONNEMENT : ARRET A L''ENTREE DES DONNEES ',/,&
'@    =========                                               ',/,&
'@    IL FAUT UTILISER UNE VARIABLE ENERGETIQUE.              ',/,&
'@                                                            ',/,&
'@  Lorsque le rayonnement est utilise (phase ',I10   ,', il  ',/,&
'@    faut indiquer qu''un scalaire represente la variable    ',/,&
'@    energetique (temperature ou enthalpie) en renseignant   ',/,&
'@    ISCALT(',I10   ,') dans usini1 (numero du scalaire).    ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier usini1.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

!===============================================================================
! 4. FAIT-ON DU RAYONNEMENT SEMI-TRANSPARENT ?
!===============================================================================

iirayo = 0
do iphas = 1,nphas
  if (irayon(iphas).eq.1 .or. irayon(iphas).eq.2) iirayo = 1
enddo

!===============================================================================
! 5. VERIFICATIONS (uniquement s'il y a du rayonnement)
!===============================================================================

iok = 0

if (iirayo.eq.1) then

! --> ISUIRD

  if (isuird.ne.0 .and. isuird.ne.1 ) then
    write(nfecra,4000)isuird
    iok = iok + 1
  endif

! --> NFREQR

  if (nfreqr.le.0) then
    write(nfecra,4010)nfreqr
    iok = iok + 1
  endif

! --> NDIREC
!     Choix entre 32 et 128 directions (cf raysol)
  if (ndirec.ne.32 .and. ndirec.ne.128 ) then
    write(nfecra,4020) ndirec
    iok = iok + 1
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

! --> IRAYVP & IRAYVF
!     Choix entre -1 et 1
  do iphas = 1,nphas
    if (irayon(iphas).eq.1 .or. irayon(iphas).eq.2) then
      do ii = 1,nbrayp
        if (irayvp(ii,iphas).ne.1 .and.                           &
            irayvp(ii,iphas).ne.-1) then
          write(nfecra,4060) nbrvap(ii,iphas), irayvp(ii,iphas)
          iok = iok + 1
        endif
      enddo
      do ii = 1,nbrayf
        if (irayvf(ii,iphas).ne.1 .and.                           &
            irayvf(ii,iphas).ne.-1) then
          write(nfecra,4070) nbrvaf(ii,iphas), irayvf(ii,iphas)
          iok = iok + 1
        endif
      enddo
    endif
  enddo

! --> CONNECTIVITE
  iphaok = 0
  do iphas = 1,nphas
    do ii = 1,nbrayf
      if (irayvf(ii,iphas).eq.1) iphaok = iphaok + 1
    enddo
  enddo
  if (iphaok.ne.0) then
    if (nnod.eq.0 .or. (lndfac.eq.0 .and. lndfbr.eq.0)) then
      write(nfecra,4080) nnod, lndfac, lndfbr
      iok = iok + 1
    endif
  endif

endif

!--> Stop si erreur.

if(iok.ne.0) then
  call csexit (1)
  !==========
endif

 4000 format(                                                           &
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
 4010 format(                                                           &
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
 4020 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : RAYONNEMENT : ARRET A L''ENTREE DES DONNEES ',/,&
'@    =========                                               ',/,&
'@    NOMBRE DE DIRECTIONS NON ADMISSIBLE                     ',/,&
'@                                                            ',/,&
'@  Le nombre de directions doit etre 32 ou 128 (NDIREC)      ',/,&
'@    Il vaut ici NDIREC = ',I10                               ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier usray1.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 4030 format(                                                           &
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
 4040 format(                                                           &
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
 4050 format(                                                           &
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
 4060 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : RAYONNEMENT : ERREUR A L''ENTREE DES DONNEES',/,&
'@    =========                                               ',/,&
'@    INDICATEUR DE SORTIE EN POSTPROCESSING                  ',/,&
'@    POUR ',A40                                               ,/,&
'@    NON ADMISSIBLE                                          ',/,&
'@                                                            ',/,&
'@  L''indicateur de postprocessing doit etre -1 ou 1         ',/,&
'@    Il vaut ici IRAYVP = ',I10                               ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier usray1.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 4070 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : RAYONNEMENT : ERREUR A L''ENTREE DES DONNEES',/,&
'@    =========                                               ',/,&
'@    INDICATEUR DE SORTIE EN POSTPROCESSING                  ',/,&
'@    POUR ',A40                                               ,/,&
'@    NON ADMISSIBLE                                          ',/,&
'@                                                            ',/,&
'@  L''indicateur de postprocessing doit etre -1 ou 1         ',/,&
'@    Il vaut ici IRAYVF = ',I10                               ,/,&
'@                                                            ',/,&
'@  Le calcul sera engage.                                    ',/,&
'@                                                            ',/,&
'@  On conseille donc de verifier usray1.                     ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 4080 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : RAYONNEMENT : ERREUR A L''ENTREE DES DONNEES',/,&
'@    =========                                               ',/,&
'@    LES CONNECTIVITES FACES -> SOMMETS OU LES COORDONNEES   ',/,&
'@    DES NOEUDS SOMMETS NE SONT PAS FOURNIS                  ',/,&
'@                                                            ',/,&
'@  Le module rayonnement est active mais les noeuds sommets  ',/,&
'@   ou les connectivites faces -> sommets n''ont pas ete     ',/,&
'@   fournis au noyau.                                        ',/,&
'@                                                            ',/,&
'@    NNDOD  = ', I10                                          ,/,&
'@    LNDFAC = ', I10                                          ,/,&
'@    LNDFBR = ', I10                                          ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier si le fichier geomet contient bien les           ',/,&
'@    connectivites faces -> sommets si IFOENV = 0.           ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)


!===============================================================================
! 6. NOMBRE DE PHASES POUR LESQUELLES ON FAIT DU RAYONNEMENT
!===============================================================================

!--> NPHAST : NOMBRE DE PHASES POUR LESQUELLES ON FAIT DU RAYONNEMENT
!--> IRAPHA : NUMEROS DES PHASES POUR LESQUELLES ON FAIT DU RAYONNEMENT

!--> NPHASC : EN CP, NOMBRE DE PHASE (= 1 FORCEMENT) + NOMBRE DE CLASSES
!             SERT A STOCKER LE COEFF D'ABSORPTION PAR EXEMPLE POUR LE
!             MELANGE ET LES CLASSES DE PARTICULES (VOIR MEMRA1)

nphast = 0
do iphas = 1,nphas
  if (irayon(iphas).ge.1) then
    nphast = nphast +1
    irapha(nphast) = iphas
  endif
enddo

nphasc = nphast

!--> Pour le Charbon Pulverise :
if ( ippmod(icp3pl) .ge. 0 ) then
  nphasc = nphasc + nclacp
endif

!--> Pour le Fuel :
if ( ippmod(icfuel) .ge. 0 ) then
  nphasc = nphasc + nclafu
endif

return

end
