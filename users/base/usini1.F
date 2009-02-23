c@a
c@versb
C-----------------------------------------------------------------------
C
CVERS
C
C
C     This file is part of the Code_Saturne Kernel, element of the
C     Code_Saturne CFD tool.
C
C     Copyright (C) 1998-2008 EDF S.A., France
C
C     contact: saturne-support@edf.fr
C
C     The Code_Saturne Kernel is free software; you can redistribute it
C     and/or modify it under the terms of the GNU General Public License
C     as published by the Free Software Foundation; either version 2 of
C     the License, or (at your option) any later version.
C
C     The Code_Saturne Kernel is distributed in the hope that it will be
C     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
C     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C     GNU General Public License for more details.
C
C     You should have received a copy of the GNU General Public License
C     along with the Code_Saturne Kernel; if not, write to the
C     Free Software Foundation, Inc.,
C     51 Franklin St, Fifth Floor,
C     Boston, MA  02110-1301  USA
C
C-----------------------------------------------------------------------
c@verse
C  FONCTION  :
C  ---------
c@foncb
CFONC
CFONC  ROUTINES UTILISATEUR POUR ENTREE DES PARAMETRES DE CALCUL (COMMONS)
CFONC   ON PASSE ICI POUR TOUS LES CALCULS
CFONC
CFONC  SI L'ON UTILISE L'IHM DE CODE_SATURNE, LE FICHIER PRESENT N'EST
CFONC    PAS INDISPENSABLE (IL COMPLETE OU PREND LE PAS SUR LES
CFONC    PARAMETRES ENTRES PAR L'INTERFACE)
CFONC
CFONC  ON TROUVERA PLUSIEURS ROUTINES DANS LE FICHIER PRESENT, CHACUNE
CFONC    DESTINEE A RENSEIGNER CERTAINS PARAMETRES PARTICULIERS
CFONC
CFONC  POUR MODIFIER LA VALEUR PAR DEFAUT DE PARAMETRES N'APPARAISSANT
CFONC    PAS DANS LES EXEMPLES FOURNIS ICI, INTERVENIR CI-DESSOUS DANS
CFONC    - USIPSU   POUR LES OPTIONS NUMERIQUES ET PHYSIQUES
CFONC    - USIPES   POUR LES OPTIONS RELATIVES AUX ENTREES-SORTIES
CFONC
CFONC  PAR CONVENTION "PHYSIQUE PARTICULIERE" RENVOIE AUX MODULES
CFONC    SUIVANTS UNIQUEMENT :
CFONC            CHARBON PULVERISE, COMBUSTION GAZ, ELECTRIQUE
CFONC
CFONC
c@fonce
C-----------------------------------------------------------------------
C
C
C***********************************************************************
C
C
c@a
                        SUBROUTINE USIPPH
C                       *****************
C     --------------------------------------------------------------
     & ( NPHMAX, NPHAS , IIHMPU, NFECRA , ITURB , ICP , IVERIF )
C     --------------------------------------------------------------
C
C***********************************************************************
C  FONCTION  :
C  ---------
c@foncb
CFONC
CFONC  ROUTINE UTILISATEUR POUR ENTREE
CFONC    DES PARAMETRES DEPENDANT DU NOMBRE DE PHASES
CFONC
c@fonce
C-----------------------------------------------------------------------
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! NPHMAX       ! E  !  -> ! NOMBRE DE PHASES MAXIMUM             !
CARGU ! NPHAS        ! E  !  -> ! NOMBRE DE PHASES ACTIVES             !
CARGU ! IIHMPU       ! E  !  -> ! INDIQUE SI UN FICHIER DE DONNEES DE  !
CARGU !              !    !     ! L'IHM EST UTILISE (1: OUI, 0: NON)   !
CARGU ! NFECRA       ! E  !  -> ! NUMERO D'UNITE STD POUR IMPRESSION   !
CARGU ! ITURB(NPHMAX)! TE ! <-> ! MODELE DE TURBULENCE                 !
CARGU ! ICP  (NPHMAX)! TE ! <-> ! INDICATEUR DE CP UNIFORME OU NON     !
CARGU ! IVERIF       ! E  !  -> ! INDICATEUR DES TESTS ELEMENTAIRES    !
CARGU !______________!____!_____!______________________________________!
c@argue
C
c@commb
CCOMM                             COMMONS
CCOMM .______________.____._____.______________________________________.
CCOMM !    NOM       !TYPE!MODE !                   ROLE               !
CCOMM !______________!____!_____!______________________________________!
CCOMM !______________!____!_____!______________________________________!
c@comme
C
C     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
C            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
C     MODE : -> DONNEE, <- RESULTAT, <-> DONNEE MODIFIEE,
C            - TABLEAU DE TRAVAIL
C***********************************************************************
C
      IMPLICIT NONE
C
C***********************************************************************
C     DONNEES EN COMMON
C***********************************************************************
C
C
C     AUCUN COMMON NE DOIT APPARAITRE ICI
C
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER NPHMAX, NPHAS, IIHMPU, NFECRA
      INTEGER ITURB(NPHMAX), ICP(NPHMAX)
      INTEGER IVERIF
C
C VARIABLES LOCALES
C
      INTEGER IPHAS
C
C***********************************************************************
C
C TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_DEBUT
C=======================================================================
C 0.  CE TEST PERMET A L'UTILISATEUR D'ETRE CERTAIN QUE C'EST
C       SA VERSION DU SOUS PROGRAMME QUI EST UTILISEE
C       ET NON CELLE DE LA BIBLIOTHEQUE
C     SI UN FICHIER ISSU DE L'IHM EST UTILISE, LE SOUS-PROGRAMME N'EST
C       PAS NECESSAIREMENT INDISPENSABLE (RETURN DANS LA VERSION DE
C       LA BIBLIOTHEQUE)
C=======================================================================
C
      IF (IVERIF.EQ.0) THEN
        IF(IIHMPU.EQ.1) THEN
          RETURN
        ELSE
          WRITE(NFECRA,9000)
          CALL CSEXIT (1)
        ENDIF
      ENDIF
C
 9000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@     LE SOUS-PROGRAMME UTILISATEUR usipph DOIT ETRE COMPLETE',/,
     &'@       DANS LE FICHIER usini1.F                             ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
C TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_FIN
C
C=======================================================================
C
C
C     ON NE RENSEIGNERA DANS CE SOUS-PROGRAMME QUE LES PARAMETRES
C
C       qui y apparaissent deja, A L'EXCLUSION DE tout autre.
C                                ================
C
C
C     SI L'ON NE DISPOSE PAS DE L'INTERFACE DE CODE_SATURNE :
C
C       on renseignera TOUS les parametres qui apparaissent dans
C                      ====                          ce sous-programme.
C
C
C CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_DEBUT
C
C     SI L'ON    DISPOSE     DE L'INTERFACE DE CODE_SATURNE :
C
C       on trouvera dans les sous-programmes utilisateur
C       des exemples commentes sur le modele de la section presente.
C
C       l'utilisateur pourra, si necessaire, les decommenter et les
C       adapter a ses besoins.
C
C CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_FIN
C
C=======================================================================
C
C --- Turbulence (pour chaque phase)
C       0...Laminaire
C      10...Longueur de melange
C      20...k-epsilon
C      21...k-epsilon a production lineaire
C      30...Rij-epsilon standard (LRR)
C      31...Rij-epsilon SSG
C      40...LES (Smagorinsky)
C      41...LES (Dynamique)
C      50...v2f (phi-model)
C      60...k-omega SST
C  Pour 10, contacter l'equipe de developpement avant utilisation
C
C CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_DEBUT
C
      IPHAS = 1
      ITURB(IPHAS) = 20
C
C CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_FIN
C
C
C --- Chaleur massique variable (ICP=1) ou non (ICP=0)
C       pour chaque phase IPHAS
C
C     A renseigner uniquement si les physiques particulieres (charbon,
C       combustion, electrique) NE SONT PAS activees.
C
C     Pour ces physiques particulieres, il NE FAUT PAS modifier ICP ici
C       et les choix suivants sont imposes
C          charbon et combustion : CP constant ;
C          electrique            : CP variable.
C
C     Attention : completer usphyv avec la loi donnant Cp
C     =========   si et seulement si on a choisi CP variable ici
C                                                   (avec ICP(IPHAS)=1)
C
C CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_DEBUT
C
      IPHAS = 1
      ICP(IPHAS) = 0
C
C CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_FIN
C
C----
C FORMATS
C----
C
C
C
      RETURN
      END
c@z
C
C
C***********************************************************************
C
C
c@a
                        SUBROUTINE USINSC
C                       *****************
C     --------------------------------------------------------------
     & ( IIHMPU, NFECRA , NSCAUS , IVERIF )
C     --------------------------------------------------------------
C
C***********************************************************************
C  FONCTION  :
C  ---------
c@foncb
CFONC
CFONC  ROUTINE UTILISATEUR POUR ENTREE
CFONC    DU NOMBRE DE SCALAIRES UTILISATEUR
CFONC
c@fonce
C-----------------------------------------------------------------------
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! IIHMPU       ! E  !  -> ! INDIQUE SI UN FICHIER DE DONNEES DE  !
CARGU !              !    !     ! L'IHM EST UTILISE (1: OUI, 0: NON)   !
CARGU ! NFECRA       ! E  !  -> ! NUMERO D'UNITE STD POUR IMPRESSION   !
CARGU ! NSCAUS       ! E  ! <-> ! NOMBRE DE SCALAIRES UTILISATEUR      !
CARGU ! IVERIF       ! E  !  -> ! INDICATEUR DES TESTS ELEMENTAIRES    !
CARGU !______________!____!_____!______________________________________!
c@argue
C
c@commb
CCOMM                             COMMONS
CCOMM .______________.____._____.______________________________________.
CCOMM !    NOM       !TYPE!MODE !                   ROLE               !
CCOMM !______________!____!_____!______________________________________!
CCOMM !______________!____!_____!______________________________________!
c@comme
C
C     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
C            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
C     MODE : -> DONNEE, <- RESULTAT, <-> DONNEE MODIFIEE,
C            - TABLEAU DE TRAVAIL
C***********************************************************************
C
      IMPLICIT NONE
C
C***********************************************************************
C     DONNEES EN COMMON
C***********************************************************************
C
C
C     AUCUN COMMON NE DOIT APPARAITRE ICI
C
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER IIHMPU, NFECRA
      INTEGER NSCAUS
      INTEGER IVERIF
C
C VARIABLES LOCALES
C
C
C***********************************************************************
C
C TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_DEBUT
C=======================================================================
C 0.  CE TEST PERMET A L'UTILISATEUR D'ETRE CERTAIN QUE C'EST
C       SA VERSION DU SOUS PROGRAMME QUI EST UTILISEE
C       ET NON CELLE DE LA BIBLIOTHEQUE
C     SI UN FICHIER ISSU DE L'IHM EST UTILISE, LE SOUS-PROGRAMME N'EST
C       PAS NECESSAIREMENT INDISPENSABLE (RETURN DANS LA VERSION DE
C       LA BIBLIOTHEQUE)
C=======================================================================
C
      IF (IVERIF.EQ.0) THEN
        IF(IIHMPU.EQ.1) THEN
          RETURN
        ELSE
          WRITE(NFECRA,9000)
          CALL CSEXIT (1)
        ENDIF
      ENDIF
C
 9000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@     LE SOUS-PROGRAMME UTILISATEUR usinsc DOIT ETRE COMPLETE',/,
     &'@       DANS LE FICHIER usini1.F                             ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
C TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_FIN
C
C=======================================================================
C
C
C     ON NE RENSEIGNERA DANS CE SOUS-PROGRAMME QUE LES PARAMETRES
C
C       qui y apparaissent deja, A L'EXCLUSION DE tout autre.
C                                ================
C
C
C     SI L'ON NE DISPOSE PAS DE L'INTERFACE DE CODE_SATURNE :
C
C       on renseignera TOUS les parametres qui apparaissent dans
C                      ====                          ce sous-programme.
C
C
C CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_DEBUT
C
C     SI L'ON    DISPOSE     DE L'INTERFACE DE CODE_SATURNE :
C
C       on trouvera dans les sous-programmes utilisateur
C       des exemples commentes sur le modele de la section presente.
C
C       l'utilisateur pourra, si necessaire, les decommenter et les
C       adapter a ses besoins.
C
C CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_FIN
C
C=======================================================================
C
C --- Nombre de scalaires UTILISATEUR (thermiques ou non et quelle que
C       soit leur phase porteuse). Il s'agit des scalaires qui viennent
C       en supplement des scalaires "de base" suivants (naturellement
C       inclus dans la modelisation) :
C        - pression
C        - grandeurs turbulentes
C        - NSCAPP scalaires introduits par l'eventuel modele de
C          combustion, charbon ou electrique mis en oeuvre
C
C     Ainsi, pour un calcul sans modele de combustion, de charbon
C       pulverise ni de modele electrique, les scalaires utilisateurs
C       pourront etre par exemple :
C        - la temperature ou l'enthalpie,
C        - des fractions massiques de scalaires transportes
C        - la moyenne du carre des fluctuations d'un autre scalaire
C          utilisateur
C
C     Le nombre de scalaires maximal est donne par NSCAMX dans paramx.h :
C       il s'agit de la valeur maximale admissible pour NSCAUS + NSCAPP.
C
C
C     Imposer NSCAUS = 0 s'il n'y a pas de scalaire utilisateur.
C
C CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_DEBUT
C
      NSCAUS = 0
C
C CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_FIN
C
C----
C FORMATS
C----
C
C
C
      RETURN
      END
c@z
C
C
C***********************************************************************
C
C
c@a
                        SUBROUTINE USIPSC
C                       *****************
C     --------------------------------------------------------------
     & ( NSCMAX, NSCAUS, IIHMPU, NFECRA, ISCAVR, IVISLS , IVERIF )
C     --------------------------------------------------------------
C
C***********************************************************************
C  FONCTION  :
C  ---------
c@foncb
CFONC
CFONC  ROUTINE UTILISATEUR POUR ENTREE
CFONC    DE PARAMETRES DEPENDANT DU NOMBRE DE SCALAIRES UTILISATEUR
CFONC
c@fonce
C-----------------------------------------------------------------------
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! NSCMAX       ! E  !  -> ! NOMBRE DE SCALAIRES MAXIMAL          !
CARGU ! NSCAUS       ! E  !  -> ! NOMBRE DE SCALAIRES UTILISATEUR      !
CARGU ! IIHMPU       ! E  !  -> ! INDIQUE SI UN FICHIER DE DONNEES DE  !
CARGU !              !    !     ! L'IHM EST UTILISE (1: OUI, 0: NON)   !
CARGU ! NFECRA       ! E  !  -> ! NUMERO D'UNITE STD POUR IMPRESSION   !
CARGU ! ISCAVR(NSCMAX! TE ! <-> ! NUMERO DU SCALAIRE CORRESPONDANT     !
CARGU !              !    !     ! POUR LES SCALAIRES VARIANCE          !
CARGU ! IVISLS(NSCMAX! TE ! <-> ! DIFFUSIVITE DES SCALAIRES UNIFORME   !
CARGU !              !    !     ! OU NON                               !
CARGU ! IVERIF       ! E  !  -> ! INDICATEUR DES TESTS ELEMENTAIRES    !
CARGU !______________!____!_____!______________________________________!
c@argue
C
c@commb
CCOMM                             COMMONS
CCOMM .______________.____._____.______________________________________.
CCOMM !    NOM       !TYPE!MODE !                   ROLE               !
CCOMM !______________!____!_____!______________________________________!
CCOMM !______________!____!_____!______________________________________!
c@comme
C
C     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
C            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
C     MODE : -> DONNEE, <- RESULTAT, <-> DONNEE MODIFIEE,
C            - TABLEAU DE TRAVAIL
C***********************************************************************
C
      IMPLICIT NONE
C
C***********************************************************************
C     DONNEES EN COMMON
C***********************************************************************
C
C
C     AUCUN COMMON NE DOIT APPARAITRE ICI
C
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER NSCMAX, NSCAUS, IIHMPU, NFECRA
      INTEGER ISCAVR(NSCMAX), IVISLS(NSCMAX)
      INTEGER IVERIF
C
C VARIABLES LOCALES
C
      INTEGER IUTILE, ISCAL
C
C***********************************************************************
C
C TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_DEBUT
C=======================================================================
C 0.  CE TEST PERMET A L'UTILISATEUR D'ETRE CERTAIN QUE C'EST
C       SA VERSION DU SOUS PROGRAMME QUI EST UTILISEE
C       ET NON CELLE DE LA BIBLIOTHEQUE
C     SI UN FICHIER ISSU DE L'IHM EST UTILISE, LE SOUS-PROGRAMME N'EST
C       PAS NECESSAIREMENT INDISPENSABLE (RETURN DANS LA VERSION DE
C       LA BIBLIOTHEQUE)
C=======================================================================
C
      IF (IVERIF.EQ.0) THEN
        IF(IIHMPU.EQ.1) THEN
          RETURN
        ELSE
          WRITE(NFECRA,9000)
          CALL CSEXIT (1)
        ENDIF
      ENDIF
C
 9000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@     LE SOUS-PROGRAMME UTILISATEUR usipsc DOIT ETRE COMPLETE',/,
     &'@       DANS LE FICHIER usini1.F                             ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
C TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_FIN
C
C=======================================================================
C
C
C     ON NE RENSEIGNERA DANS CE SOUS-PROGRAMME QUE LES PARAMETRES
C
C       qui y apparaissent deja, A L'EXCLUSION DE tout autre.
C                                ================
C
C
C     SI L'ON NE DISPOSE PAS DE L'INTERFACE DE CODE_SATURNE :
C
C       on renseignera TOUS les parametres qui apparaissent dans
C                      ====                          ce sous-programme.
C
C
C CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_DEBUT
C
C     SI L'ON    DISPOSE     DE L'INTERFACE DE CODE_SATURNE :
C
C       on trouvera dans les sous-programmes utilisateur
C       des exemples commentes sur le modele de la section presente.
C
C       l'utilisateur pourra, si necessaire, les decommenter et les
C       adapter a ses besoins.
C
C CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_FIN
C
C=======================================================================
C
C --- Moyenne du carre des fluctuations d'un scalaire UTILISATEUR :
C     Si on souhaite qu'un scalaire utilisateur j represente la moyenne
C      du carre des fluctuations du scalaire utilisateur k, on indique
C      ISCAVR(j) = k.
C     Les valeurs prises par ISCAVR sont donc naturellement superieures
C      ou egales a 1 et inferieures ou egales au nombre total de scalaires.
C      Ainsi, si on impose ISCAVR(j) = k, il faut
C            0<j<NSCAUS+1, 0<k<NSCAUS+1 et j different de k.
C
C     Par exemple pour que le scalaire utilisateur 3 soit la moyenne
C      du carre des fluctuations du scalaire utilisateur 2, on impose :
C                         ISCAVR(3) = 2
C      avec NSCAUS au moins egal a 3.
C
C     Ne pas intervenir si l'on ne souhaite pas inclure explicitement a
C       la simulation la moyenne du carre des fluctuations d'un scalaire
C       utilisateur.
C
C     Pour les scalaires non utilisateur relatifs a des physiques
C       particulieres, (charbon, combustion, electrique : voir usppmo)
C       implicitement definis selon le modele, les informations sont
C       donnees automatiquement par ailleurs : on ne modifie pas ISCAVR.
C
C
C CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_DEBUT
C
C     Le test sur IUTILE permet de desactiver les instructions (qui
C       ne sont fournies qu'a titre d'exemple a adapter)
C
      IUTILE = 0
      IF(IUTILE.EQ.1) THEN
        ISCAVR(3) = 2
      ENDIF
C
C CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_FIN
C
C
C
C
C --- Diffusivite variable (IVISLS=1) ou constante (IVISLS=0) pour
C       chaque scalaire UTILISATEUR     hormis ceux qui
C       representent la moyenne du carre des fluctuations d'un autre.
C
C     Pour les scalaires utilisateur ISCAL qui representent la moyenne du
C       carre des fluctuations d'un autre scalaire utilisateur,
C                                on ne renseigne pas IVISLS(ISCAL) ici.
C       C'est l'objet du test sur ISCAVR(ISCAL) dans l'exemple ci-dessous.
C       En effet, la diffusivite de la moyenne du carre des fluctuations
C       d'un scalaire est supposee avoir le meme comportement que la
C       diffusivite de ce scalaire.
C
C     Pour les scalaires non utilisateur relatifs a des physiques
C       particulieres, (charbon, combustion, electrique : voir usppmo)
C       implicitement definis selon le modele,
C       les informations sont donnees automatiquement par ailleurs :
C                                          on ne modifie pas IVISLS ici.
C
C     Attention : completer usphyv avec la loi donnant la diffusivite
C     =========   si et seulement si on a choisi IVISLS=1 ici
C
C
C
C
C CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_DEBUT
C
      DO ISCAL = 1, NSCAUS
C
C     Pour les scalaires utilisateur qui ne representent pas
C       la moyenne du carre des fluctuations d'un autre scalaire
        IF(ISCAVR(ISCAL).LE.0) THEN
C
          IVISLS(ISCAL) = 0
C
        ENDIF
C
      ENDDO
C
C CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_FIN
C
C----
C FORMATS
C----
C
C
C
      RETURN
      END
c@z
C
C
C***********************************************************************
C
C
c@a
                        SUBROUTINE USIPGL
C                       *****************
C     --------------------------------------------------------------
     & ( NPHMAX, NESMAX,
     &   IESPRE, IESDER, IESCOR, IESTOT,
     &   NPHAS , IIHMPU, NFECRA,
     &   IDTVAR, IPUCOU, IPHYDR, IALGCE , IESCAL , IVERIF )
C     --------------------------------------------------------------
C
C***********************************************************************
C  FONCTION  :
C  ---------
c@foncb
CFONC
CFONC  ROUTINE UTILISATEUR POUR ENTREE
CFONC    DE PARAMETRES GLOBAUX
CFONC
c@fonce
C-----------------------------------------------------------------------
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! NPHMAX       ! E  !  -> ! NOMBRE DE PHASES MAXIMAL             !
CARGU ! NESMAX       ! E  !  -> ! NOMBRE D'ESTIMATEURS D'ERREUR        !
CARGU !              !    !     !   MAXIMAL POUR UNE PHASE             !
CARGU ! IESPRE       ! E  !  -> ! NUMERO DE L'ESTIMATEUR D'ERREUR      !
CARGU !              !    !     !   PREDICTION                         !
CARGU ! IESDER       ! E  !  -> ! NUMERO DE L'ESTIMATEUR D'ERREUR      !
CARGU !              !    !     !   DERIVE                             !
CARGU ! IESCOR       ! E  !  -> ! NUMERO DE L'ESTIMATEUR D'ERREUR      !
CARGU !              !    !     !   CORRECTION                         !
CARGU ! IESTOT       ! E  !  -> ! NUMERO DE L'ESTIMATEUR D'ERREUR      !
CARGU !              !    !     !   TOTAL                              !
CARGU ! NPHAS        ! E  !  -> ! NOMBRE DE PHASES ACTIVES             !
CARGU ! IIHMPU       ! E  !  -> ! INDIQUE SI UN FICHIER DE DONNEES DE  !
CARGU !              !    !     ! L'IHM EST UTILISE (1: OUI, 0: NON)   !
CARGU ! NFECRA       ! E  !  -> ! NUMERO D'UNITE STD POUR IMPRESSION   !
CARGU ! IDTVAR       ! E  ! <-  ! INDICATEUR PAS DE TEMPS VARIABLE     !
CARGU ! IPUCOU       ! E  ! <-  ! INDICATEUR COUPLAGE U-P RENFORCE     !
CARGU ! IPHYDR       ! E  ! <-  ! INDICATEUR DE PRISE EN COMPTE DE     !
CARGU !              !    !     ! L'EQUILIBRE ENTRE LE GRADIENT DE     !
CARGU !              !    !     ! PRESSION ET LES TERMES DE GRAVITE ET !
CARGU !              !    !     ! DE PERTE DE CHARGE                   !
CARGU ! IALGCE       ! E  !  -> ! INDICATEUR DE METHODE DE CALCUL DES  !
CARGU !              !    !     !  CENTRES DE GRAVITE DES CELLULES     !
CARGU ! IESCAL       ! TE ! <-> ! INDICATEUR D'ACTIVATION DES          !
CARGU !(NESMAX,NPHMAX!    !     ! ESTIMATEURS D'ERREUR POUR NAVIER     !
CARGU !              !    !     ! STOKES                               !
CARGU ! IVERIF       ! E  !  -> ! INDICATEUR DES TESTS ELEMENTAIRES    !
CARGU !______________!____!_____!______________________________________!
c@argue
C
c@commb
CCOMM                             COMMONS
CCOMM .______________.____._____.______________________________________.
CCOMM !    NOM       !TYPE!MODE !                   ROLE               !
CCOMM !______________!____!_____!______________________________________!
CCOMM !______________!____!_____!______________________________________!
c@comme
C
C     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
C            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
C     MODE : -> DONNEE, <- RESULTAT, <-> DONNEE MODIFIEE,
C            - TABLEAU DE TRAVAIL
C***********************************************************************
C
      IMPLICIT NONE
C
C***********************************************************************
C     DONNEES EN COMMON
C***********************************************************************
C
C
C     AUCUN COMMON NE DOIT APPARAITRE ICI
C
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER NPHMAX, NESMAX
      INTEGER IESPRE, IESDER, IESCOR, IESTOT
      INTEGER NPHAS , IIHMPU, NFECRA
      INTEGER IDTVAR, IPUCOU, IPHYDR
      INTEGER IESCAL(NESMAX,NPHMAX)
      INTEGER IVERIF
C
C VARIABLES LOCALES
C
      INTEGER IPHAS, IALGCE
C
C***********************************************************************
C
C TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_DEBUT
C=======================================================================
C 0.  CE TEST PERMET A L'UTILISATEUR D'ETRE CERTAIN QUE C'EST
C       SA VERSION DU SOUS PROGRAMME QUI EST UTILISEE
C       ET NON CELLE DE LA BIBLIOTHEQUE
C     SI UN FICHIER ISSU DE L'IHM EST UTILISE, LE SOUS-PROGRAMME N'EST
C       PAS NECESSAIREMENT INDISPENSABLE (RETURN DANS LA VERSION DE
C       LA BIBLIOTHEQUE)
C=======================================================================
C
      IF (IVERIF.EQ.0) THEN
        IF(IIHMPU.EQ.1) THEN
          RETURN
        ELSE
          WRITE(NFECRA,9000)
          CALL CSEXIT (1)
        ENDIF
      ENDIF
C
 9000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@     LE SOUS-PROGRAMME UTILISATEUR usipgl DOIT ETRE COMPLETE',/,
     &'@       DANS LE FICHIER usini1.F                             ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
C TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_FIN
C
C=======================================================================
C
C
C     ON NE RENSEIGNERA DANS CE SOUS-PROGRAMME QUE LES PARAMETRES
C
C       qui y apparaissent deja, A L'EXCLUSION DE tout autre.
C                                ================
C
C
C     SI L'ON NE DISPOSE PAS DE L'INTERFACE DE CODE_SATURNE :
C
C       on renseignera TOUS les parametres qui apparaissent dans
C                      ====                          ce sous-programme.
C
C
C CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_DEBUT
C
C     SI L'ON    DISPOSE     DE L'INTERFACE DE CODE_SATURNE :
C
C       on trouvera dans les sous-programmes utilisateur
C       des exemples commentes sur le modele de la section presente.
C
C       l'utilisateur pourra, si necessaire, les decommenter et les
C       adapter a ses besoins.
C
C CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_FIN
C
C=======================================================================
C
C --- Pas de temps (0 : uniforme et constant
C                   1 : variable en temps et uniforme en espace
C                   2 : variable en espace et en temps
C                  -1 : algorithme stationnaire)
C
C CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_DEBUT
C
      IDTVAR = 0
C
C CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_FIN
C
C
C --- Couplage vitesse/pression (0 : algorithme classique,
C                                1 : couplage instationnaire)
C     Uniquement en monophasique
C
C CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_DEBUT
C
      IPUCOU = 0
C
C CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_FIN
C
C
C --- Prise en compte de la pression hydrostatique
C                               (0 : algorithme usuel
C                                1 : prise en compte explicite)
C     Uniquement en monophasique
C
C CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_DEBUT
C
      IPHYDR = 0
C
C CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_FIN
C
C
C --- Estimateurs pour Navier-Stokes (champ de vitesse non fige)
C     On conseille de realiser une suite de calcul sur quelques pas de
C       temps en activant les plus parlants ci dessous
C        (=2 pour activer, =0 pour desactiver).
C
C CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_DEBUT
C
      IPHAS = 1
C       div(rho u) -Gamma
      IESCAL(IESCOR,IPHAS) = 0
C       precision de la resolution de la quantite de mouvement
      IESCAL(IESTOT,IPHAS) = 0
C
C CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_FIN
C
C
C----
C FORMATS
C----
C
C
C
      RETURN
      END
c@z
C
C
C***********************************************************************
C
C
c@a
                        SUBROUTINE USIPSU
C                       *****************
C     --------------------------------------------------------------
     & ( NMODPP , IVERIF )
C     --------------------------------------------------------------
C
C***********************************************************************
C  FONCTION  :
C  ---------
c@foncb
CFONC
CFONC  ROUTINE UTILISATEUR POUR ENTREE
CFONC    DE PARAMETRES UTILISATEUR SUPPLEMENTAIRES
CFONC
c@fonce
C-----------------------------------------------------------------------
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! NMODPP       ! E  !  -> ! NOMBRE DE MODELES PHYS.PART. ACTIVES !
CARGU ! IVERIF       ! E  !  -> ! INDICATEUR DES TESTS ELEMENTAIRES    !
CARGU !______________!____!_____!______________________________________!
c@argue
C
c@commb
CCOMM                             COMMONS
CCOMM .______________.____._____.______________________________________.
CCOMM !    NOM       !TYPE!MODE !                   ROLE               !
CCOMM !______________!____!_____!______________________________________!
CCOMM !______________!____!_____!______________________________________!
c@comme
C
C     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
C            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
C     MODE : -> DONNEE, <- RESULTAT, <-> DONNEE MODIFIEE,
C            - TABLEAU DE TRAVAIL
C***********************************************************************
C
      IMPLICIT NONE
C
C***********************************************************************
C     DONNEES EN COMMON
C***********************************************************************
C
      INCLUDE "paramx.h"
      INCLUDE "cstnum.h"
      INCLUDE "dimens.h"
      INCLUDE "numvar.h"
      INCLUDE "optcal.h"
      INCLUDE "cstphy.h"
      INCLUDE "entsor.h"
      INCLUDE "vector.h"
      INCLUDE "parall.h"
      INCLUDE "period.h"
      INCLUDE "ihmpre.h"
      INCLUDE "ppppar.h"
      INCLUDE "ppthch.h"
      INCLUDE "ppincl.h"
      INCLUDE "coincl.h"
      INCLUDE "cpincl.h"
      INCLUDE "elincl.h"
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER NMODPP
      INTEGER IVERIF
C
C VARIABLES LOCALES
C
      INTEGER IPHAS, IUTILE, II, JJ, IMOM
C
C***********************************************************************
C
C TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_DEBUT
C=======================================================================
C 0.  CE TEST PERMET A L'UTILISATEUR D'ETRE CERTAIN QUE C'EST
C       SA VERSION DU SOUS PROGRAMME QUI EST UTILISEE
C       ET NON CELLE DE LA BIBLIOTHEQUE
C     SI UN FICHIER ISSU DE L'IHM EST UTILISE, LE SOUS-PROGRAMME N'EST
C       PAS NECESSAIREMENT INDISPENSABLE (RETURN DANS LA VERSION DE
C       LA BIBLIOTHEQUE)
C=======================================================================
C
      IF (IVERIF.EQ.0) THEN
        IF(IIHMPR.EQ.1) THEN
          RETURN
        ELSE
          WRITE(NFECRA,9000)
          CALL CSEXIT (1)
        ENDIF
      ENDIF
C
 9000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@     LE SOUS-PROGRAMME UTILISATEUR usipsu DOIT ETRE COMPLETE',/,
     &'@       DANS LE FICHIER usini1.F                             ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
C TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_FIN
C
C=======================================================================
C
C
C     CE SOUS-PROGRAMME PERMET DE RENSEIGNER LES PARAMETRES
C
C       qui n'apparaissent pas deja dans les autres sous-programmes
C       du present fichier.
C
C
C     IL EST POSSIBLE D'AJOUTER OU DE RETRANCHER DES PARAMETRES
C
C
C     LE NUMERO DES PROPRIETES PHYSIQUES ET DES  VARIABLES EST CONNU
C
C
C CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_DEBUT
C
C     SI L'ON    DISPOSE     DE L'INTERFACE DE CODE_SATURNE :
C
C       on trouvera dans les sous-programmes utilisateur
C       des exemples commentes sur le modele de la section presente.
C
C       l'utilisateur pourra, si necessaire, les decommenter et les
C       adapter a ses besoins.
C
C CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_FIN
C
C=======================================================================
C
C
C Options du calcul (optcal.h)
C ============================
C
C --- Suite de calcul : ISUITE ( = 1) ou non (0)
C     Avec relecture du fichier suite auxiliaire ILEAUX ( = 1) ou non (0)
C
C CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_DEBUT
C
      ISUITE = 0
      ILEAUX = 1
C
C CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_FIN
C
C --- Duree
C       NTMABS = numero absolu du dernier pas de temps desire
C         si on a deja fait 10 pas de temps
C         et qu'on veut en faire 10 autres,
C           il faut imposer NTMABS a 10 + 10 = 20
C
C CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_DEBUT
C
      NTMABS = 10
C
C CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_FIN
C
C --- Pas de temps de reference
C     L'exemple donne ici est probablement inadapte a votre cas
C
C CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_DEBUT
C
      DTREF  = 0.01D0
C
C CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_FIN
C
C --- Pas de temps maximal DTMAX
C     Imposer une valeur a partir de grandeurs caracteristiques du cas
C      sinon, le code prendra par defaut un multiple de dtref
C     Ex. avec
C        Ld : longueur "dynamique" ie par exemple longueur du domaine
C        Ud : Vitesse debitante caracteristique
C        Lt : longueur thermique ie par exemple la hauteur du domaine
C                                                    (selon la gravite)
C        Delta_rho/rho : ecart relatif de masse volumique
C        g : acceleration de la pesanteur
C
C
C     DTMAX = Min(Ld/Ud,Sqrt(Lt/(g Delta_rho/rho)))
C
C
C
C
C
C
C
C
C --- Temperature ou enthalpie
C
C
C
C   Lorsque des physiques particulieres sont activees
C                                       (charbon, combustion, electrique)
C     on NE renseigne PAS cette section : on NE modifie ni ISCALT ni ISCSTH
C                                   (le test IF(NMODPP.EQ.0) sert a cela).
C
C
C   Par contre, si ces physiques particulieres ne sont PAS activees :
C
C     Si un scalaire UTILISATEUR represente la temperature ou l'enthalpie
C       (de la phase IPHAS) :
C          on donne le numero de ce scalaire dans ISCALT(IPHAS) et
C          on indique  ISCSTH(ISCALT(IPHAS)) = 1 si c'est la temperature
C                  ou  ISCSTH(ISCALT(IPHAS)) = 2 si c'est l'enthalpie
C
C     Si aucun scalaire ne represente la temperature ou l'enthalpie (de
C       la phase IPHAS)
C          on indique ISCALT(IPHAS) = -1
C          et on ne renseigne pas ISCSTH(ISCALT(IPHAS))
C
C
C     Pour le module rayonnement mis en oeuvre en dehors des physiques
C      charbon, combustion, electrique :
C      si l'on a choisi de resoudre en temperature (c'est-a-dire si
C      ISCSTH(ISCALT(IPHAS)) = 1), la temperature du fluide est alors
C      supposee en KELVIN (attention aux conditions aux limites et a
C      l'expression des proprietes physiques dependantes de la
C      temperature).
C      Neanmoins, bien que ce soit deconseille, si l'on souhaite
C      que le solveur fluide travaille avec une temperature en degres
C      Celsius, on doit imposer ISCSTH(ISCALT(IPHAS)) = -1.
C      Ce choix est source d'erreurs pour l'utilisateur. En effet, les
C      conditions aux limites pour la temperature du fluide seront alors
C      en degres Celsius, alors que les conditions aux limites imposees
C      pour le rayonnement dans usray2 devront, pour leur part, etre
C      obligatoirement en Kelvin.
C
C
C
C    Si les physiques particulieres ne sont pas activees
C       (charbon, combustion, electrique : voir usppmo)
C
C CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_DEBUT
C
      IF(NMODPP.EQ.0) THEN
C
        IPHAS = 1
C         Numero du scalaire representant la temperature ou l'enthalpie
C                                               ou -1 s'il n'y en a pas
C         Lorsque le choix est fait par l'interface de Code_Saturne,
C           le scalaire representant la temperature ou l'enthalpie
C           est obligatoirement le premier
        ISCALT(IPHAS) = -1
C
C         S'il y a une variable temperature ou enthalpie
        IF(ISCALT(IPHAS).GT.0) THEN
C           on indique si c'est la temperature (=1) ou l'enthalpie (=2)
          ISCSTH(ISCALT(IPHAS)) = 1
        ENDIF
C
      ENDIF
C
C CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_FIN
C
C --- Calcul (suite) a champ de vitesse fige (1 oui, 0 non)
C
C CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_DEBUT
C
      ICCVFG = 0
C
C CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_FIN
C
C --- Methode des vortex pour les conditions d'entree en L.E.S.
C              (0 : non activee,  1 : activee)
C     La methode des vortex ne concerne que les modeles L.E.S.
C       et n'est valable qu'avec une seule phase
C     Pour utiliser la methode des vortex, veuillez renseigner
C       le sous-programme usvort.F
C
C CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_DEBUT
C
      IPHAS = 1
      IF (ITYTUR(IPHAS).EQ.4) THEN
        IVRTEX = 0
      ENDIF
C
C CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_FIN
C
C --- Schema convectif
C
C       BLENCV = 0 pour upwind (ordre 1 en espace, "stable mais diffusif")
C              = 1 pour centre/second order (ordre 2 en espace)
C           on peut utiliser des valeurs reelles intermediaires
C       Ici on choisit
C         pour la vitesse de la phase 1 et les scalaires utilisateurs :
C            un schema centre upwind a 100% de centre (BLENCV=1)
C         pour les autres variables
C            la valeur par defaut du code (upwind en standard, centre en LES)
C
C     En particulier, pour les scalaires utilisateur
C       si l'on suspecte un niveau de diffusion numerique trop grand
C         sur une variable IVAR representant un scalaire utilisateur
C         ISCAL (avec IVAR=ISCA(ISCAL)), il peut etre utile d'imposer
C         BLENCV(IVAR) = 1.0D0 pour utiliser un schema d'ordre 2 en
C         espace pour la convection. Pour la temperature ou l'enthalpie,
C         en particulier, on pourra donc choisir dans ce cas :
C          BLENCV(ISCA(ISCALT(IPHAS))) = 1.0D0
C
C       Pour les scalaires non utilisateur relatifs a des physiques
C         particulieres, (charbon, combustion, electrique : voir usppmo)
C         implicitement definis selon le modele,
C         les informations sont donnees automatiquement par ailleurs :
C                                         on ne modifie pas BLENCV ici.
C
C
C CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_DEBUT
C
      IPHAS = 1
C
      BLENCV(IU(IPHAS)) = 1.0D0
      BLENCV(IV(IPHAS)) = 1.0D0
      BLENCV(IW(IPHAS)) = 1.0D0
      IF(NSCAUS.GE.1) THEN
        DO II = 1, NSCAUS
          BLENCV(ISCA(II)) = 1.0D0
        ENDDO
      ENDIF
C
C CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_FIN
C
C
C --- Parametres  multigrille algebrique
C
C     IMGR = 0 : PAS DE MULTIGRILLE
C     IMGR = 1 : MULTIGRILLE ALGEBRIQUE
C
C     Uniquement disponible pour la pression
C
C CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_DEBUT
C
      IPHAS = 1
      IMGR(IPR(IPHAS)) = 1
C
C CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_FIN
C
C
C=================================================================
C
C --- Stabilisation en regime turbulent
C
C     Pour les cas difficiles, une stabilisation peut etre obtenue en
C     ne reconstruisant pas les flux de convection et de diffusion
C     pour les variables du modele de turbulence, soit
C       en k-epsilon : IF (ITYTUR(IPHAS).EQ.2) THEN
C          IRCFLU(IK(IPHAS))   = 0 et IRCFLU(IEP(IPHAS))  = 0
C       en Rij-epsilon : IF (ITYTUR(IPHAS).EQ.3) THEN
C          IRCFLU(IR11(IPHAS)) = 0,   IRCFLU(IR22(IPHAS)) = 0,
C          IRCFLU(IR33(IPHAS)) = 0,
C          IRCFLU(IR12(IPHAS)) = 0,   IRCFLU(IR23(IPHAS)) = 0,
C          IRCFLU(IR23(IPHAS)) = 0,
C                                  et IRCFLU(IEP(IPHAS))  = 0
C     (noter que la variable ITYTUR vaut ITURB/10)
C
C CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_DEBUT
C
C     Le test sur IUTILE permet de desactiver les instructions (qui
C       ne sont fournies qu'a titre d'exemple a adapter)
C
      IUTILE = 0
      IF(IUTILE.EQ.1) THEN
C
        IPHAS = 1
        IF (ITURB(IPHAS).EQ.20) THEN
          IRCFLU(IK(IPHAS))   = 0
          IRCFLU(IEP(IPHAS))  = 0
        ENDIF
C
      ENDIF
C
C CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_FIN
C
C
C Constantes physiques (cstphy.h)
C ===============================
C
C --- gravite (g en m/s2, avec le signe dans le repere de calcul)
C
C CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_DEBUT
C
      GX = 0.D0
      GY = 0.D0
      GZ = 0.D0
C
C CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_FIN
C
C --- Proprietes  de reference du fluide (pour chaque phase)
C
C       RO0        : masse volumique en kg/m3
C       VISCL0     : viscosite dynamique en kg/(m s)
C       CP0        : chaleur specifique en J/(degres kg)
C       T0         : temperature de reference en Kelvin
C       P0         : pression totale de reference en Pascal
C                    le calcul est par contre effectue sur
C                    une pression reduite P*=Ptot-ro0*g.(x-xref)
C                    (sauf en compressible)
C       XYZP0(3,.) : coordonnees du point de reference pour
C                    la pression totale (la ou elle vaut P0)
C
C     En general, il n'est pas necessaire de fournir un point
C       de reference XYZ0. S'il y a des sorties, le code
C       prendra le centre de la face de sortie de reference.
C       Par contre, si on prevoit de fixer explicitement des
C       conditions de Dirichlet pour la pression, mieux vaut
C       indiquer a quelle reference les valeurs se rapportent
C       (pour une meilleure resolution de la pression reduite)
C
C
C     Les autres proprietes seront fournies par defaut dans tous
C       les cas.
C
C     Neanmoins, on peut noter que :
C
C       Dans les cas standard (ni combustion gaz, ni charbon, ni electrique,
C                              ni compressible) :
C       ---------------------
C         RO0, VISCL0 et CP0
C             sont utiles et representent soit les proprietes du fluide
C             si elles sont constantes, soit de simples valeurs moyennes
C             pour l'initialisation si les proprietes sont variables et
C             donnees dans usphyv.
C         T0  est inutile
C         P0  est utile mais n'est pas utilise dans une loi d'etat. P0
C             est une valeur de reference pour le solveur incompressible
C             qui servira a caler la pression a la sortie (eventuelle)
C             du domaine. On peut la prendre nulle ou egale a une valeur
C             physique en Pascal.
C
C       En modules electriques :
C       ----------------------
C         RO0, VISCL0 et CP0
C             sont utiles mais representent de simples valeurs moyennes
C             initiales ; la masse volumique, la viscosite dynamique
C             moleculaire et la chaleur massique sont donnees
C             obligatoirement dans PROPCE (qu'elles soient physiquement
C             variables ou non) : voir uselph pour le module effet Joule
C             et le fichier de donnees dp_ELE en arc electrique.
C         T0  est utile et doit etre en Kelvin (> 0) mais represente
C             une simple valeur d'initialisation.
C         P0  est utile mais n'est pas utilise dans une loi d'etat. P0
C             est une valeur de reference pour le solveur incompressible
C             qui servira a caler la pression a la sortie (eventuelle)
C             du domaine. On peut la prendre nulle ou egale a une valeur
C             physique en Pascal.
C
C       En combustion gaz :
C       -----------------
C         RO0 est inutile (il est recalcule automatiquement par la
C             loi des gaz parfaits a partir de T0 et P0).
C         VISCL0 est indispensable : c'est la viscosite dynamique
C             moleculaire, supposee constante pour le fluide
C         CP0 est indispensable : c'est la chaleur massique,
C             supposee constante (modelisation des termes sources faisant
C             intervenir un Nusselt local dans
C             le module lagrangien, valeur de reference permettant
C             de calculer un couple (temperature, coefficient d'echange)
C             en rayonnement)
C         T0  est indispensable et doit etre en Kelvin (> 0)
C         P0  est indispensable et doit etre en Pascal (> 0)
C
C       En charbon pulverise :
C       --------------------
C         RO0 est inutile (il est recalcule automatiquement par la
C             loi des gaz parfaits a partir de T0 et P0).
C         VISCL0 est indispensable : c'est la viscosite dynamique
C             moleculaire, supposee constante pour le fluide (son
C             effet est a priori faible devant les effets turbulents)
C         CP0 est indispensable : c'est la chaleur massique,
C             supposee constante (modelisation des termes sources faisant
C             intervenir un Nusselt local dans le module charbon ou
C             le module lagrangien, valeur de reference permettant
C             de calculer un couple (temperature, coefficient d'echange)
C             en rayonnement)
C         T0  est indispensable et doit etre en Kelvin (> 0)
C         P0  est indispensable et doit etre en Pascal (> 0)
C
C       En compressible :
C       --------------------
C         RO0 est inutile, stricto sensu ; neanmoins, comme l'experience
C             montre que les utilisateurs se servent souvent de cette
C             variable, on demande de lui affecter ici une valeur
C             strictement positive (par exemple, une valeur initiale)
C         VISCL0 est utile et represente la viscosite dynamique
C             moleculaire, lorsqu'elle est constante ou une valeur qui
C             servira lors des initialisations (ou dans les conditions
C             d'entree de la turbulence, selon le choix de l'utilisateur)
C         CP0 est indispensable : c'est la chaleur massique,
C             supposee constante dans la thermodynamique disponible par
C             defaut
C         T0  est indispensable et doit etre en Kelvin (> 0)
C         P0  est indispensable et doit etre en Pascal (> 0)
C             Avec la loi thermodynamique disponible par defaut,
C             T0 et P0 servent a l'initialisation de la masse volumique.
C         XYZP0 est inutile car la variable pression represente directement
C             la presion totale
C
C
C CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_DEBUT
C
      IPHAS = 1
C
      RO0   (IPHAS) = 0.235D0
      VISCL0(IPHAS) = 0.84D-6
      CP0   (IPHAS) = 1219.D0
C
C CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_FIN
C
      IPHAS = 1
C
      T0    (IPHAS) = 1000.D0 + 273.15D0
      P0    (IPHAS) = 1.013D5
C     On ne specifie XYZ0 que si on compte explicitement fixer
C       des conditions de Dirichlet sur la pression
C      XYZP0(1,IPHAS) = 0.D0
C      XYZP0(2,IPHAS) = 0.D0
C      XYZP0(3,IPHAS) = 0.D0
C
C
C --- IROVAR IVIVAR : masse volumique et viscosite constantes ou non
C
C     Lorsqu'un modele physique particuliere a ete active
C                  (charbon, combustion, modules electriques, compressible)
C       on ne renseigne PAS les variables IROVAR et IVIVAR ici
C                                  elles sont definies automatiquement.
C       Pour le compressible, cependant, IVIVAR peut être modifie dans
C                                  le sous-programme utilisateur uscfx1.
C
C     Lorsqu'aucun modele physique particuliere n'a ete active,
C         (i.e. charbon, combustion, electrique, compressible : voir usppmo)
C       il faut indiquer si la masse volumique et la viscosite moleculaire
C         sont constantes (IROVAR=0, IVIVAR=0)
C           ou variables  (IROVAR=1, IVIVAR=1)
C
C       si elles sont variables , il faut alors definir la loi dans usphyv.
C       si elles sont constantes, elles prennent les valeurs RO0 et VISCL0
C
C       a titre d'exemple, on suppose ci-dessous qu'elles sont constantes
C
C CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_DEBUT
C
      IF(NMODPP.EQ.0) THEN
        IPHAS = 1
        IROVAR(IPHAS) = 0
        IVIVAR(IPHAS) = 0
      ENDIF
C
C CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_FIN
C
C
C --- Valeurs min SCAMIN et max SCAMAX admissibles pour
C        chaque scalaire UTILISATEUR :
C
C      Les resultats sont clipppes a la fin de chaque pas de temps
C
C      Si SCAMIN > SCAMAX, on ne clippe pas.
C
C      Pour un scalaire JJ representant la moyenne du carre des
C         fluctuations d'un autre, on pourra s'abstenir de renseigner
C         ces valeurs (un clipping par defaut est mis en place).
C         C'est l'objet du test sur ISCAVR(JJ) dans l'exemple ci-dessous.
C
C
C      Pour les scalaires non utilisateur relatifs a des physiques
C         particulieres, (charbon, combustion, electrique : voir usppmo)
C         implicitement definis selon le modele, les informations sont
C         donnees automatiquement par ailleurs : on ne renseigne pas
C         SCAMIN SCAMAX.
C
C CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_DEBUT
C
C     S'il y a des scalaires utilisateur
      IF(NSCAUS.GT.0) THEN
C
C       On boucle sur les scalaires utilisateurs :
        DO JJ = 1, NSCAUS
C         Pour les scalaires qui ne sont pas des variances
          IF(ISCAVR(JJ).LE.0) THEN
C           On definit la borne min et la borne max
            SCAMIN(JJ) =-GRAND
            SCAMAX(JJ) =+GRAND
          ENDIF
        ENDDO
C
      ENDIF
C
C CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_FIN
C
C --- Diffusivite de reference VISLS0 en kg/(m s) pour
C        chaque scalaire UTILISATEUR   hormis  pour ceux qui
C        representent la moyenne du carre des fluctuations d'un autre.
C
C     Pour les scalaires non utilisateur relatifs a des physiques
C       particulieres, (charbon, combustion, electrique : voir usppmo)
C       implicitement definis selon le modele,
C       les informations sont donnees automatiquement par ailleurs :
C                                        on ne modifie pas VISLS0 ici.
C
C     Pour les scalaires utilisateur JJ qui representent la moyenne du
C       carre des fluctuations d'un autre scalaire utilisateur,
C                                   on ne renseigne pas VISLS0(JJ) ici.
C       C'est l'objet du test sur ISCAVR(JJ) dans l'exemple ci-dessous.
C       En effet la diffusivite de la moyenne du carre des fluctuations
C       d'un scalaire est supposee identique à la diffusivite de ce
C       scalaire.
C
C     Lorsqu'on n'a pas active de physique particuliere
C       (charbon, combustion, electrique) et si un scalaire utilisateur
C        represente la temperature ou l'enthalpie,
C        on definit ici pour la temperature ou l'enthalpie :
C                                     VISLS0(ISCALT(IPHAS)) = Lambda/Cp
C
C     Ici, a titre d'exemple, on affecte a la VISCL0 la viscosite de
C       la phase porteuse, ce qui convient a des traceurs passifs qui
C       suivent le fluide.
C
C
C CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_DEBUT
C
C     S'il y a des scalaires utilisateur
      IF(NSCAUS.GT.0) THEN
C
C       On boucle sur les scalaires utilisateurs :
        DO JJ = 1, NSCAUS
C         Pour les scalaires qui ne sont pas des variances
          IF(ISCAVR(JJ).LE.0) THEN
C           On definit la diffusivite
            VISLS0(JJ) = VISCL0(IPHSCA(JJ))
          ENDIF
        ENDDO
C
      ENDIF
C
C CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_FIN
C
C
C --- Vitesse de reference pour l'initialisation de la turbulence (m2/s)
C      (utile seulement en turbulence)
C
C CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_DEBUT
C
      IPHAS = 1
      UREF(IPHAS)    = 1.D0
C
C CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_FIN
C
C --- Echelle de longueur de reference en metres pour
C      initialisation de epsilon (et clipping particulier
C      de la turbulence, mais ce n'est pas l'option par defaut)
C      Donner une valeur de l'ordre de la plus grande
C      dimension du domaine physique dans lequel l'ecoulement
C      peut se developper
C      (utile seulement en turbulence)
C
C CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_DEBUT
C
      IPHAS = 1
      ALMAX(IPHAS) = -GRAND
C
C CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_FIN
C
C --- Definition des moments
C     (au maximum NBMOMX moments, correlations d'ordre maximum NDGMOX)
C
C     On calcule des moyennes temporelles du type <f1*f2*f3*...*fn>
C     Les fi sont des variables definies aux cellules (tableaux RTP et PROPCE)
C
C        IDFMOM(i,IMOM) repere la variable fi du moment IMOM
C          si IDFMOM > 0 c'est une variable resolue (RTP)
C          si IDFMOM < 0 c'est une variable auxiliaire (PROPCE)
C        IMOOLD(IMOM) donne en cas de suite le numero, dans l'ancien calcul
C                du moment a utiliser pour initialiser le moment IMOM du
C                nouveau calcul (par defaut IMOOLD(IMOM)=IMOM).
C                La valeur -1 indique qu'on doit reinitialiser le moment IMOM
C        NTDMOM(IMOM) donne le pas de temps de debut du calcul du moment
C
C     On donne ci dessous l'exemple du calcul des moments <u> et <rho u v>
C       le moment <u> est relu dans le fichier suite si on est en suite,
C         le moment <rho u v> est reinitialise a zero
C       le moment <u> est calcule a partir du pas de temps 1000
C         le moment <rho u v> est calcule a partir du pas de temps 10000
C
C
C     Le test sur IUTILE permet de desactiver les instructions (qui
C       ne sont fournies qu'a titre d'exemple a adapter)
C
      IUTILE = 0
      IF(IUTILE.EQ.1) THEN
C
C     Premier moment : <u>
        IMOM  = 1
        IPHAS = 1
        IDFMOM(1,IMOM) =  IU(IPHAS)
        NTDMOM(IMOM)   =  1000
C     Second moment : <rho u v>
        IMOM  = 2
        IPHAS = 1
        IDFMOM(1,IMOM) = -IROM(IPHAS)
        IDFMOM(2,IMOM) =  IU(IPHAS)
        IDFMOM(3,IMOM) =  IV(IPHAS)
        IMOOLD(IMOM)   = -1
        NTDMOM(IMOM)   =  10000
C
      ENDIF
C
C----
C FORMATS
C----
C
C
C
      RETURN
      END
c@z
C
C
C***********************************************************************
C
C
c@a
                        SUBROUTINE USIPES
C                       *****************
C     --------------------------------------------------------------
     & ( NMODPP , IVERIF )
C     --------------------------------------------------------------
C
C***********************************************************************
C  FONCTION  :
C  ---------
c@foncb
CFONC
CFONC  ROUTINE UTILISATEUR POUR ENTREE
CFONC    DE PARAMETRES UTILISATEUR SUPPLEMENTAIRES
CFONC
c@fonce
C-----------------------------------------------------------------------
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! NMODPP       ! E  !  -> ! NOMBRE DE MODELES PHYS.PART. ACTIVES !
CARGU ! IVERIF       ! E  !  -> ! INDICATEUR DES TESTS ELEMENTAIRES    !
CARGU !______________!____!_____!______________________________________!
c@argue
C
c@commb
CCOMM                             COMMONS
CCOMM .______________.____._____.______________________________________.
CCOMM !    NOM       !TYPE!MODE !                   ROLE               !
CCOMM !______________!____!_____!______________________________________!
CCOMM !______________!____!_____!______________________________________!
c@comme
C
C     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
C            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
C     MODE : -> DONNEE, <- RESULTAT, <-> DONNEE MODIFIEE,
C            - TABLEAU DE TRAVAIL
C***********************************************************************
C
      IMPLICIT NONE
C
C***********************************************************************
C     DONNEES EN COMMON
C***********************************************************************
C
      INCLUDE "paramx.h"
      INCLUDE "cstnum.h"
      INCLUDE "dimens.h"
      INCLUDE "numvar.h"
      INCLUDE "optcal.h"
      INCLUDE "cstphy.h"
      INCLUDE "entsor.h"
      INCLUDE "vector.h"
      INCLUDE "parall.h"
      INCLUDE "period.h"
      INCLUDE "ihmpre.h"
      INCLUDE "ppppar.h"
      INCLUDE "ppthch.h"
      INCLUDE "ppincl.h"
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER NMODPP
      INTEGER IVERIF
C
C VARIABLES LOCALES
C
      INTEGER IPHAS, IPP, IMOM
C
C***********************************************************************
C
C TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_DEBUT
C=======================================================================
C 0.  CE TEST PERMET A L'UTILISATEUR D'ETRE CERTAIN QUE C'EST
C       SA VERSION DU SOUS PROGRAMME QUI EST UTILISEE
C       ET NON CELLE DE LA BIBLIOTHEQUE
C     SI UN FICHIER ISSU DE L'IHM EST UTILISE, LE SOUS-PROGRAMME N'EST
C       PAS NECESSAIREMENT INDISPENSABLE (RETURN DANS LA VERSION DE
C       LA BIBLIOTHEQUE)
C=======================================================================
C
      IF (IVERIF.EQ.0) THEN
        IF(IIHMPR.EQ.1) THEN
          RETURN
        ELSE
          WRITE(NFECRA,9000)
          CALL CSEXIT (1)
        ENDIF
      ENDIF
C
 9000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@     LE SOUS-PROGRAMME UTILISATEUR usipes DOIT ETRE COMPLETE',/,
     &'@       DANS LE FICHIER usini1.F                             ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
C TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_FIN
C
C=======================================================================
C
C
C     CE SOUS-PROGRAMME PERMET DE RENSEIGNER LES PARAMETRES
C
C       qui n'apparaissent pas deja dans les autres sous-programmes
C       du present fichier.
C
C
C     IL EST POSSIBLE D'AJOUTER OU DE RETRANCHER DES PARAMETRES
C
C
C     LE NUMERO DES PROPRIETES PHYSIQUES ET DES  VARIABLES EST CONNU
C
C
C CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_DEBUT
C
C     SI L'ON    DISPOSE     DE L'INTERFACE DE CODE_SATURNE :
C
C       on trouvera dans les sous-programmes utilisateur
C       des exemples commentes sur le modele de la section presente.
C
C       l'utilisateur pourra, si necessaire, les decommenter et les
C       adapter a ses besoins.
C
C CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_FIN
C
C=======================================================================
C
C=======================================================================
C 1. ENTREE-SORTIE (optcal.h)
C=======================================================================
C
C --- ecriture du fichier suite auxiliaire IECAUX = 1 oui, 0 non
C
C CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_DEBUT
C
      IECAUX = 1
C
C CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_FIN
C
C
C --- pas des sorties listing
C
C CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_DEBUT
C
      NTLIST = 1
C
C CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_FIN
C
C --- post traitement
C
C     ICHRVL : post traitement du domaine fluide (oui 1/non 0)
C     ICHRBO : post traitement du bord du domaine (oui 1/non 0)
C     ICHRSY : Post traitement des zones couplees avec Syrthes
C              (1 oui, 0 non)
C     ICHRMD : indique si les maillages ecrits seront :
C               0 : fixes,
C               1 : deformables a topologie constante,
C               2 : modifiables (pourront etre completement redefinis en
C                   cours de calcul via le sous-programme USMPST).
C              10 : comme INDMOD = 0, avec champ de déplacement
C              11 : comme INDMOD = 1, avec champ de déplacement
C              12 : comme INDMOD = 2, avec champ de déplacement
C
C     FMTCHR : format de sortie, parmi
C              'EnSight Gold', 'MED_fichier', ou 'CGNS'
C     OPTCHR : options associees au format de sortie, separees
C              par des virgules, parmi
C              'text'              (format texte, pour EnSight)
C              'binary'            (format binaire, choix par defaut)
C              'big_endian'        (force les sorties EnSight binaires
C                                   en mode 'big-endian')
C              'discard_polygons'  (ignore faces de type polygone)
C              'discard_polyhedra' (ignore cellules de type polyedre)
C              'divide_polygons'   (decoupe faces de type polygone)
C              'divide_polyhedra'  (decoupe cellules de type polyedre)
C              'split_tensors'     (ecrit les tenseurs en tant que
C                                   scalaires separes)
C
C CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_DEBUT
C
      ICHRVL = 1
      ICHRBO = 0
      ICHRSY = 0
C
      ICHRMD = 0
C
      FMTCHR = 'EnSight Gold'
      OPTCHR = 'binary'
C
C CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_FIN
C
C
C --- pas des chrono (-1 : une seule sortie en fin de calcul)
C                    (valeur strictement positive : periode)
C
C CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_DEBUT
C
      NTCHR = -1
C
C CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_FIN
C
C
C --- pas des sorties historiques
C
C CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_DEBUT
C
      NTHIST = 1
C
C CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_FIN
C
C
C --- Nombre de sondes et positions (limite a NCAPTM=100)
C
C CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_DEBUT
C
      NCAPT  = 4
      XYZCAP(1,1) = 0.30D0
      XYZCAP(2,1) = 0.15D0
      XYZCAP(3,1) = 0.01D0
C
      XYZCAP(1,2) = 0.30D0
      XYZCAP(2,2) = 0.00D0
      XYZCAP(3,2) = 0.01D0
C
      XYZCAP(1,3) = 0.30D0
      XYZCAP(2,3) =-0.08D0
      XYZCAP(3,3) = 0.01D0
C
      XYZCAP(1,4) = 0.60D0
      XYZCAP(2,4) =-0.05D0
      XYZCAP(3,4) = 0.01D0
C
C CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_FIN
C
C
C --- extension historique
      EXTHIS = 'hst'
C
C
C --- variable courante
C
C     Comme pour les autres variables,
C       si l'on n'affecte pas les tableaux suivants,
C       les valeurs par defaut seront utilisees
C
C     NOMVAR( ) = nom de la variable
C     ICHRVR( ) = sortie chono (oui 1/non 0)
C     ILISVR( ) = suivi listing (oui 1/non 0)
C     IHISVR( ) = sortie historique (nombre de sondes et numeros)
C     si IHISVR(.,1)  = -1 sortie sur toutes les sondes
C
C     NB : Seuls les 8 premiers caracteres du nom seront repris dans le
C          listing le plus detaille
C
C
C
C CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_DEBUT
C
C     Variables dynamiques courantes
C
C     Exemples pour la phase 1
      IPHAS = 1
C
C     variable pression
      IPP = IPPRTP(IPR   (IPHAS))
      NOMVAR(IPP)   = 'Pression'
      ICHRVR(IPP)   = 1
      ILISVR(IPP)   = 1
      IHISVR(IPP,1) = -1
C
C     variable v1x
      IPP = IPPRTP(IU    (IPHAS))
      NOMVAR(IPP)   = 'VitesseX'
      ICHRVR(IPP)   = 1
      ILISVR(IPP)   = 1
      IHISVR(IPP,1) = -1
C
C     variable v1y
      IPP = IPPRTP(IV    (IPHAS))
      NOMVAR(IPP)   = 'VitesseY'
      ICHRVR(IPP)   = 1
      ILISVR(IPP)   = 1
      IHISVR(IPP,1) = -1
C
C     variable v1z
      IPP = IPPRTP(IW    (IPHAS))
      NOMVAR(IPP)   = 'VitesseZ'
      ICHRVR(IPP)   = 1
      ILISVR(IPP)   = 1
      IHISVR(IPP,1) = -1
C
      IF(ITYTUR(IPHAS).EQ.2) THEN
C
C     energie turbulente
        IPP = IPPRTP(IK    (IPHAS))
        NOMVAR(IPP)   = 'EnerTurb'
        ICHRVR(IPP)   = 1
        ILISVR(IPP)   = 1
        IHISVR(IPP,1) = -1
C
C     dissipation turbulente
        IPP = IPPRTP(IEP   (IPHAS))
        NOMVAR(IPP)   = 'Dissip'
        ICHRVR(IPP)   = 1
        ILISVR(IPP)   = 1
        IHISVR(IPP,1) = -1
C
      ELSEIF(ITYTUR(IPHAS).EQ.3) THEN
C
C     tensions de Reynolds
        IPP = IPPRTP(IR11  (IPHAS))
        NOMVAR(IPP)   = 'Tens.R11'
        ICHRVR(IPP)   = 1
        ILISVR(IPP)   = 1
        IHISVR(IPP,1) = -1
C
C     tensions de Reynolds
        IPP = IPPRTP(IR22  (IPHAS))
        NOMVAR(IPP)   = 'Tens.R22'
        ICHRVR(IPP)   = 1
        ILISVR(IPP)   = 1
        IHISVR(IPP,1) = -1
C
C     tensions de Reynolds
        IPP = IPPRTP(IR33  (IPHAS))
        NOMVAR(IPP)   = 'Tens.R33'
        ICHRVR(IPP)   = 1
        ILISVR(IPP)   = 1
        IHISVR(IPP,1) = -1
C
C     tensions de Reynolds
        IPP = IPPRTP(IR12  (IPHAS))
        NOMVAR(IPP)   = 'Tens.R12'
        ICHRVR(IPP)   = 1
        ILISVR(IPP)   = 1
        IHISVR(IPP,1) = -1
C
C     tensions de Reynolds
        IPP = IPPRTP(IR13  (IPHAS))
        NOMVAR(IPP)   = 'Tens.R13'
        ICHRVR(IPP)   = 1
        ILISVR(IPP)   = 1
        IHISVR(IPP,1) = -1
C
C     tensions de Reynolds
        IPP = IPPRTP(IR23  (IPHAS))
        NOMVAR(IPP)   = 'Tens.R23'
        ICHRVR(IPP)   = 1
        ILISVR(IPP)   = 1
        IHISVR(IPP,1) = -1
C
C     dissipation turbulente
        IPP = IPPRTP(IEP   (IPHAS))
        NOMVAR(IPP)   = 'Dissip'
        ICHRVR(IPP)   = 1
        ILISVR(IPP)   = 1
        IHISVR(IPP,1) = -1
C
      ELSEIF(ITURB(IPHAS).EQ.50) THEN
C
C     energie turbulente
        IPP = IPPRTP(IK    (IPHAS))
        NOMVAR(IPP)   = 'EnerTurb'
        ICHRVR(IPP)   = 1
        ILISVR(IPP)   = 1
        IHISVR(IPP,1) = -1
C
C     dissipation turbulente
        IPP = IPPRTP(IEP   (IPHAS))
        NOMVAR(IPP)   = 'Dissip'
        ICHRVR(IPP)   = 1
        ILISVR(IPP)   = 1
        IHISVR(IPP,1) = -1
C
C     phi
        IPP = IPPRTP(IPHI  (IPHAS))
        NOMVAR(IPP)   = 'Phi'
        ICHRVR(IPP)   = 1
        ILISVR(IPP)   = 1
        IHISVR(IPP,1) = -1
C
C     f_barre
        IPP = IPPRTP(IFB   (IPHAS))
        NOMVAR(IPP)   = 'f_barre'
        ICHRVR(IPP)   = 1
        ILISVR(IPP)   = 1
        IHISVR(IPP,1) = -1
C
      ELSEIF(ITURB(IPHAS).EQ.60) THEN
C
C-->  energie turbulente
        IPP = IPPRTP(IK    (IPHAS))
        NOMVAR(IPP)   = 'EnerTurb'
        ICHRVR(IPP)   = 1
        ILISVR(IPP)   = 1
        IHISVR(IPP,1) = -1
C
C-->  omega
        IPP = IPPRTP(IOMG  (IPHAS))
        NOMVAR(IPP)   = 'Omega'
        ICHRVR(IPP)   = 1
        ILISVR(IPP)   = 1
        IHISVR(IPP,1) = -1
C
      ENDIF
C
C CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_FIN
C
C     Variables scalaires utilisateur
C
C     On peut modifier ici les tableaux relatifs aux scalaires utilisateurs
C       mais les scalaires reserves pour les physiques particulieres sont
C       geres automatiquement. D'ou les tests sur NSCAUS qui assurent que
C       les scalaires vises sont effectivement des scalaires utilisateurs.
C       Par physique particuliere, on entend
C        uniquement celles definies par des modules specifiques
C        du code tels que charbon, combustion, electrique : voir usppmo)
C       Pour les scalaires non utilisateur relatifs a des physiques
C         particulieres, (charbon, combustion, electrique : voir usppmo)
C         implicitement definis selon le modele, les informations sont
C         donnees automatiquement par ailleurs : on ne modifie pas BLENCV.
C
C CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_DEBUT
C
      IF(ISCA(1).GT.0.AND.NSCAUS.GE.1) THEN
       IPP = IPPRTP(ISCA  (1))
        NOMVAR(IPP)  = 'scal 1'
        ICHRVR(IPP)  = 1
        ILISVR(IPP)  = 1
        IHISVR(IPP,1)= -1
      ENDIF
C
      IF(ISCA(2).GT.0.AND.NSCAUS.GE.2) THEN
        IPP = IPPRTP(ISCA  (2))
        NOMVAR(IPP)  = 'scal 2'
        ICHRVR(IPP)  = 1
        ILISVR(IPP)  = 1
        IHISVR(IPP,1)= -1
      ENDIF
C
C CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_FIN
C
C
C     Autres variables
C
      IPHAS = 1
C
C CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_DEBUT
C
C     variable masse volumique (sortie en post uniquement si variable
C                               ou si physique particuliere)
      IPP = IPPPRO(IPPROC(IROM  (IPHAS)))
      NOMVAR(IPP)   = 'masse vol'
      ICHRVR(IPP)   = MAX(IROVAR(IPHAS),NMODPP)
      ILISVR(IPP)   = 1
      IHISVR(IPP,1) = -1
C
C     chaleur specifique
      IF(ICP   (IPHAS).GT.0) THEN
        IPP = IPPPRO(IPPROC(ICP   (IPHAS)))
        NOMVAR(IPP)   = 'chal. spec.'
        ICHRVR(IPP)   = 0
        ILISVR(IPP)   = 0
        IHISVR(IPP,1) = 0
      ENDIF
C
C     viscosite laminaire
      IPP = IPPPRO(IPPROC(IVISCL(IPHAS)))
      NOMVAR(IPP)   = 'visc. laminaire'
      ICHRVR(IPP)   = 0
      ILISVR(IPP)   = 0
      IHISVR(IPP,1) = 0
C
C     viscosite turbulente
      IPP = IPPPRO(IPPROC(IVISCT(IPHAS)))
      NOMVAR(IPP)   = 'visc. turb1'
      ICHRVR(IPP)   = 1
      ILISVR(IPP)   = 1
      IHISVR(IPP,1) = -1
C
C     nombre de Courant
      IPP = IPPPRO(IPPROC(ICOUR(IPHAS)))
      NOMVAR(IPP)   = 'Nb Courant'
      ICHRVR(IPP)   = 1
      ILISVR(IPP)   = 0
      IHISVR(IPP,1) = -1
C
C     nombre de Fourier
      IPP = IPPPRO(IPPROC(IFOUR(IPHAS)))
      NOMVAR(IPP)   = 'Nb Fourier'
      ICHRVR(IPP)   = 1
      ILISVR(IPP)   = 0
      IHISVR(IPP,1) = -1
C
C     variable CSMAGO dans le cas des modeles L.E.S. dynamiques
C     (carre de la "constante" de Smagorinsky)
      IF(ISMAGO(IPHAS).GT.0) THEN
        IPP = IPPPRO(IPPROC(ISMAGO(IPHAS)))
        NOMVAR(IPP)   = 'Csdyn2'
        ICHRVR(IPP)   = 1
        ILISVR(IPP)   = 1
        IHISVR(IPP,1) = -1
      ENDIF
C
C     moyennes temporelles (exemple pour le moment 1)
      IF(NBMOMT.GT.0) THEN
        IMOM = 1
        IPP = IPPPRO(IPPROC(ICMOME(IMOM)))
        NOMVAR(IPP) = 'MoyTps01'
        ICHRVR(IPP) = 1
        ILISVR(IPP) = 1
        IHISVR(IPP,1) = -1
      ENDIF
C
C     pression totale (non definie en compressible)
      IF (IPPMOD(ICOMPF).LT.0) THEN
        IPP = IPPPRO(IPPROC(IPRTOT(IPHAS)))
        NOMVAR(IPP)   = 'Pression totale'
        ICHRVR(IPP)   = 1
        ILISVR(IPP)   = 1
        IHISVR(IPP,1) = -1
      ENDIF
C
C     pas de temps local
      IPP = IPPDT
      NOMVAR(IPP)   = 'pdt local'
      ICHRVR(IPP)   = 1
      ILISVR(IPP)   = 1
      IHISVR(IPP,1) = -1
C
C     temps caracteristique du couplage
C        instationnaire vitesse/pression
      IPP = IPPTX
      NOMVAR(IPP)   = 'Tx'
      ICHRVR(IPP)   = 1
      ILISVR(IPP)   = 1
      IHISVR(IPP,1) = -1
C
      IPP = IPPTY
      NOMVAR(IPP)   = 'Ty'
      ICHRVR(IPP)   = 1
      ILISVR(IPP)   = 1
      IHISVR(IPP,1) = -1
C
      IPP = IPPTZ
      NOMVAR(IPP)   = 'Tz'
      ICHRVR(IPP)   = 1
      ILISVR(IPP)   = 1
      IHISVR(IPP,1) = -1
C
C CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_FIN
C
C
C----
C FORMATS
C----
C
C
C
      RETURN
      END
c@z
C
C
C***********************************************************************
C
C
c@a
                        SUBROUTINE USTBTR
C                       *****************
C     -------------------------------------------------------------
     & ( NCEL   , NCELET , NFAC   , NFABOR , NNOD   ,
     &   LONGIA , LONGRA ,
     &   NIDEVE , NITUSE , NRDEVE , NRTUSE )
C     -------------------------------------------------------------
C***********************************************************************
C  FONCTION  :
C  ---------
c@foncb
CFONC
CFONC  ROUTINE UTILISATEUR POUR DEFINIR LES DIMENSIONS
CFONC    DES MACROS TABLEAUX IA ET RA
CFONC    DES TABLEAUX UTILISATEUR ITUSER ET RTUSER
CFONC    DES TABLEAUX DEVELOPPEUR IDEVEL ET RDEVEL
CFONC
c@fonce
C-----------------------------------------------------------------------
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! NCEL         ! E  !  -> ! NOMBRE DE CELLULES HORS HALO         !
CARGU ! NCELET       ! E  !  -> ! NOMBRE DE CELLULES HALO INCLUS       !
CARGU ! NFAC         ! E  !  -> ! NOMBRE DE FACES INTERNES             !
CARGU ! NFABOR       ! E  !  -> ! NOMBRE DE FACES DE BORD              !
CARGU ! NNOD         ! E  !  -> ! NOMBRE DE NOEUDS                     !
CARGU ! LONGIA       ! E  ! <-  ! DIMENSION DU TABLEAU IA              !
CARGU ! LONGRA       ! E  ! <-  ! DIMENSION DU TABLEAU RA              !
CARGU ! NIDEVE       ! E  ! <-  ! DIMENSION DU TABLEAU IDEVEL          !
CARGU ! NITUSE       ! E  ! <-  ! DIMENSION DU TABLEAU ITUSER          !
CARGU ! NRDEVE       ! E  ! <-  ! DIMENSION DU TABLEAU RDEVEL          !
CARGU ! NRTUSE       ! E  ! <-  ! DIMENSION DU TABLEAU RTUSER          !
CARGU !______________!____!_____!______________________________________!
c@argue
C
c@commb
CCOMM                             COMMONS
CCOMM .______________.____._____.______________________________________.
CCOMM !    NOM       !TYPE!MODE !                   ROLE               !
CCOMM !______________!____!_____!______________________________________!
CCOMM !______________!____!_____!______________________________________!
c@comme
C
C     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
C            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
C     MODE : -> DONNEE, <- RESULTAT, <-> DONNEE MODIFIEE,
C            - TABLEAU DE TRAVAIL
C***********************************************************************
C
      IMPLICIT NONE
C
C***********************************************************************
C     DONNEES EN COMMON
C***********************************************************************
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          NCEL  , NCELET, NFAC  , NFABOR, NNOD
      INTEGER          LONGIA, LONGRA
      INTEGER          NIDEVE, NITUSE, NRDEVE, NRTUSE
C
C***********************************************************************
C
C TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_DEBUT
C=======================================================================
C
      IF(1.EQ.1) RETURN
C
C=======================================================================
C TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_FIN
C
C=======================================================================
C 1. DIMENSION DES MACROS TABLEAUX IA ET RA :
C
C  L'utilisateur peut avoir a modifier la taille des tableaux
C    d'entiers et de reels ici : LONGIA et LONGRA respectivement
C
C  Le nombre d'entiers LONGIA et le nombre de reels LONGRA
C    dependent des options de calcul, du type d'elements,
C    d es caracteristiques du maillage (2d, 3d, hybride,
C    non conforme ...) et du nombre de variables, par exemple.
C    En k-epsilon, si on note NCEL le nombre de cellules du
C    maillage, on peut la plupart du temps utiliser la majoration
C    grossiere suivante : LONGIA =  45*NCEL et LONGRA = 220*NCEL
C    En Rij-epsilon, une majoration de 20% pourra etre appliquee.
c    Ces valeurs sont relativement elevees pour prendre en compte
C    les maillages 2D qui comprennent de nombreuses faces de bord :
C    une formule plus precise mais plus complexe serait necessaire.
C    Pour les gros cas 3D, une estimation plus precise est donnee
C    par : LONGIA =  25*NCEL et LONGRA = 120*NCEL.
C
C  Si LONGIA et LONGRA sont laisses a 0, alors ces valeurs sont
C    remplies de maniere automatique
C
C=======================================================================
C
C CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_DEBUT
C
C     DIMENSION DU MACRO TABLEAU D'ENTIERS IA
C
      LONGIA = 0
C
C     DIMENSION DU MACRO TABLEAU DE REELS RA
C
      LONGRA = 0
C
C CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_FIN
C
C=======================================================================
C 2. DIMENSIONS DES TABLEAUX UTILISATEUR ITUSER ET RTUSER
C=======================================================================
C
C CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_DEBUT
C
C     DIMENSION DU TABLEAU D'ENTIERS ITUSER
C
      NITUSE = 0
C
C     DIMENSION DU TABLEAU DE REELS RTUSER
C
      NRTUSE = 0
C
C CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_FIN
C
C
C----
C FORMATS
C----
C
C
C
      RETURN
      END
c@z
