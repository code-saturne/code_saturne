c@a
c@versb
C-----------------------------------------------------------------------
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
                        SUBROUTINE MTINI1
C                       *****************
C     -------------------------------------------------------------
C     -------------------------------------------------------------
C***********************************************************************
C  FONCTION  :
C  ---------
c@foncb
CFONC
CFONC          INITIALISATION DES COMMON MATISSE
c@fonce
C-----------------------------------------------------------------------
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
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
      INCLUDE "numvar.h"
      INCLUDE "optcal.h"
      INCLUDE "ihmpre.h"
      INCLUDE "matiss.h"
      INCLUDE "cstphy.h"
      INCLUDE "cstnum.h"
      INCLUDE "entsor.h"
      INCLUDE "parall.h"
C
C***********************************************************************
C
C ARGUMENTS
C
C VARIABLES LOCALES
C
      INTEGER          II    , JJ    , KK
      INTEGER          IPHAS
      INTEGER          IARRET
      DOUBLE PRECISION TSOR0
      DOUBLE PRECISION XMAPMN, XMAPMX
      DOUBLE PRECISION XMAPVL
      DOUBLE PRECISION CXSD  , CYSD  , POR0
      DOUBLE PRECISION USTMAT, DRFMAT
      DOUBLE PRECISION PRFMAT, ERFMAT, SRFMAT
C
C***********************************************************************
C
C=======================================================================
C 1. PAR DEFAUT MATISSE N'EST PAS ACTIVE
C=======================================================================
C
C     Ceci doit etre present en dehors du if defined xml :
C       on ne peut activer Matisse que avec IHM et XML
C
C     Ainsi, dans les autres sous-programmes Matisse, il est inutile de
C       tester si l'interface ou Xml sont disponibles (en effet, si
C       l'interface ou Xml ne sont pas disponibles, mtini1 est le seul
C       sous-programme Matisse dans lequel on entre).
C
      IMATIS = 0
C
C
C=======================================================================
C 2. LECTURE DE LA BALISE MATISSE EN XML (si xml present)
C=======================================================================
C
C --- Lecture de IMATIS
C     (si IHM presente ; sinon, Matisse ne sera pas actif)
C
      IF (IIHMPR.EQ.1) THEN
        CALL CSMTPR(IMATIS)
      ENDIF
C
C --- Si Matisse n'est pas actif, on sort a la fin du if
C     Sinon, on lit les donnees, on les traite, ...
C
      IF (IMATIS.EQ.1) THEN
C
C
C=======================================================================
C 3. REMPLISSAGE DES COMMONS A PARTIR DU XML
C=======================================================================
C
C On remplit les COMMON  de données géométriques et physiques
C  a partir du contenu du fichier xml
C
C --- /IMTGEO/
C
        CALL CSGEIN(NPTRAN, NPLGRS, NELGRS, NCHEST, NETRAN, ITYPEN)
C
C --- /RMTGEO/
C
        CALL CSGEDB(EPREGI, EPCHEM, HCONVE, RCONVE, HCHALI, HCHEVA,
     &              HFTTOI, PTRRES, FRDTRA, PLGRES, EPCHEL, DMCONT)
C
C --- /IMTPHY/
C
        CALL CSPHAT(IMDCNT, ICOFOR, ICONLG, IALVEO)
C
C --- /RMTPHY/
C
        CALL CSPHDB(DTDTMX, PUICON, TINIT,  TCRIT,  EMICON, EMIMUR,
     &              HEPCNT, DHPCNT, DEBMAS, PDCCHA, PDCFCH, DHCHEA,
     &              SDCHEA, PDCCHE, PDCCCH, DHCHES, SDCHES, PDCALG,
     &              PDCATV, ARGAMT, PDCSLG, PDCSTV, ARGAVL, AMPPDC,
     &              DHALVE, HRESO,  HPLEN,  DPVENT)
C
C --- Cartes 2D et 3D
C
        DO II = 1, NCARTE
          DO JJ = 1, NMTDIR
            CALL CSNBMP(JJ,II, NZOCAR(JJ,II))
            IF (NZOCAR(JJ,II) .GT. NZONMX) THEN
              WRITE(NFECRA,9011) NZOCAR(JJ,II), NZONMX
              CALL CSEXIT (1)
            ENDIF
          ENDDO
        ENDDO
C
        DO II = 1, NCARTE
          DO JJ = 1, NMTDIR
            DO KK = 1 , NZOCAR(JJ,II)
              CALL CSDFMP (KK ,JJ, II, XMAPMN, XMAPMX, XMAPVL)
              VIZCAR(1, KK, JJ, II) = XMAPMN
              VIZCAR(2, KK, JJ, II) = XMAPMX
              IF (II.EQ.ICPUIS) THEN
                VCARTH(KK, JJ) = XMAPVL
              ENDIF
            ENDDO
          ENDDO
        ENDDO
C
C
C=======================================================================
C 4. CALCULS COMPLEMENTAIRES AVANT VERIFICATIONS
C=======================================================================
C
C     Variables en common
C
C --- Hauteur d'erosion reduite a un nombre entier de mailles
        HERCNT = EPCHEL*INT(HEPCNT/EPCHEL)
C
C
C=======================================================================
C 5. VERIFICATIONS
C=======================================================================
C
C --- Par defaut, tout est suppose correct (on ne s'arretera donc pas)
        IARRET = 0
C
C
C --- Matisse n'est pas parallele
C
        IF (IRANGP.GE.0) THEN
          IARRET = 1
          WRITE(NFECRA,9001)
        ENDIF
C
C
C --- NSCAUS = 3
C     Matisse a besoin de 3 scalaires et 3 seulement
C     Des developpements supplementaires sont necessaires pour
C       que d'autres scalaires utilisateurs puissent etre pris en compte
C       (exemple : acces a ustssc)
C
        IF(NSCAUS.NE.3) THEN
          IARRET = 1
          WRITE(NFECRA,9021) NSCAUS
        ENDIF
C
C
C --- TCRIT > TINIT (cf calcul de VITREF)
C     TCRIT est la "Temperature d'air de sortie critique en degres C"
C       c'est une estimation de la temperature de l'air en sortie de
C       l'entrepot et donc naturellement strictement superieure a TINIT
C
        IF(TCRIT.LE.TINIT) THEN
          IARRET = 1
          WRITE(NFECRA,9031) TCRIT, TINIT
        ENDIF
C
C
C --- RCONVE = 1 si 2D
C     RCONVE est le rapport du convergent represente sur le maillage.
C       Attention, "respresente sur le maillage" est important :
C         il ne s'agit pas du convergent reel.
C       En 2D, il ne peut pas y avoir de convergent (le convergent
C         reduit la dimension X, transverse, des cheminees) et RCONVE
C         est donc necessairement unite.
C       Pour savoir si l'on fait une simulation 2D, on regarde si
C         FRDTRA est non unite (necessaire mais non suffisant : FRDTRA
C         non unite indique que l'on represente une configuration
C         3D par un calcul 2D, i.e. a une seule maille dans la
C         direction X).
C
        IF ( (ABS(FRDTRA-1.D0).GE.EPZERO).AND.
     &       (ABS(RCONVE-1.D0).GE.EPZERO)     ) THEN
          IARRET = 1
          WRITE(NFECRA,9041) FRDTRA, RCONVE
        ENDIF
C
C
C --- Si on modelise les panaches, la hauteur doit etre correcte
C     HERCNT est la hauteur reduite a un numbre entier de mailles
C       dans la direction verticale : c'est la hauteur que l'on
C       va reellement prendre en compte dans le calcul. Elle doit
C       etre positive et ne pas depasser la hauteur de l'espace
C       libre au dessus des colis :
C
C       HERCNT > 0 et HERCNT <= NCHEST*EPCHEL-HRESO
C
        IF(IMDCNT.EQ.1) THEN
          IF ( (HERCNT.LE.0.D0).OR.
     &         (HERCNT.GT.NCHEST*EPCHEL-HRESO) ) THEN
            IARRET = 1
            WRITE(NFECRA,9051) IMDCNT,
     &           HEPCNT, HERCNT,
     &           NCHEST*EPCHEL, NCHEST, EPCHEL, HRESO,
     &           HERCNT, NCHEST*EPCHEL-HRESO
          ENDIF
        ENDIF
C
C
C --- ARRET EVENTUEL
C
        IF(IARRET.NE.0) THEN
          CALL CSEXIT (1)
C         ===========
        ENDIF
C
C
C=======================================================================
C 6. VALEURS PAR DEFAUT
C=======================================================================
C
C === On ne traite qu'une seule et unique phase
C ===========================================================
C
C     On ne le dit qu'une fois, au debut, pour eviter d'en oublier
C
        IPHAS = 1
C
C
C === Options numeriques
C ===========================================================
C
C
C --- IPHYDR = 0 à cause des pertes de charges dans la cheminee
C     Plus exactement, a cause des pertes de charges en
C       sortie, avec presence de gravite.
C     A noter, cependant, qu'il semble que meme avec ICALHY = 0
C       l'option IPHYDR cause des problemes (et que c'était déjà le cas
C       dans les premieres versions de Matisse, puisque l'option etait
C       certes proposee par défaut =1, mais mise systematiquement =0
C       dans tous les cas de calcul).
C     Bien qu'il s'agisse de la valeur par defaut de Code_Saturne 1.2,
C       on laisse l'option ici dans la mesure ou des tests
C       supplementaires seraient necessaires pour comprendre pourquoi
C       on rencontre des difficultes avec les configurations etudiees
C       dans Matisse.
C
        IPHYDR = 0
C
C
C --- Par contre, l'extrapolation de la pression fournit des resultats
C       satisfaisants. Elle est utile et ne fait pas apparaitre de
C       probleme sur les cas testes.
C
        EXTRAG(IPR(IPHAS)) = 1.D0
C
C
C --- IPUCOU prend sa valeur par defaut de Code_Saturne 1.2
C     L'option est cependant conservee ici, dans la mesure ou elle
C       etait prise egale a 1 auparavant et ou cela fait une difference
C       (toute petite, certes, mais on garde ici la ligne pour memoire).
C
        IPUCOU = 0
C
C
C --- Pas de convection pour les temperatures solides
C     Pas de phenomene diffusif pour les parois (le rayonnement est
C       traite par terme source)
C
        ICONV(ISCA(ITPCMT)) = 0
        ICONV(ISCA(ITPPMT)) = 0
        IDIFF(ISCA(ITPPMT)) = 0
C
C
C --- Schema convectif : robuste
C
        BLENCV(IU(IPHAS)) = 0.0D0
        BLENCV(IV(IPHAS)) = 0.0D0
        BLENCV(IW(IPHAS)) = 0.0D0
        IF(NSCAUS.GE.1) THEN
          DO II = 1, NSCAUS
            BLENCV(ISCA(II)) = 0.0D0
          ENDDO
        ENDIF
C
C
C === Proprietes physiques et pas de temps
C ===========================================================
C
C
C --- Connectivite du rayonnement a calculer
C     (au premier passage dans mttssc)
C
        ICNROK = 0
C
C
C --- Masse volumique, chaleur massique
C
C     On prend comme reference :
C     - la masse volumique RRFMAT a la temperature TRFMAT en degres C
C     - on utilise la loi des gaz parfaits
C     - on fixe CP0
C
        RO0(IPHAS) = (TRFMAT + TKELVI)*RRFMAT /(TINIT + TKELVI)
        CP0(IPHAS) = CRFMAT
C
C     On fixe P0 a 1 atm et PRED0 a 0
C     La reference de pression est prise au niveau de l'alimentation
        P0(IPHAS) = 1.013D5
        PRED0(IPHAS) = 0.D0
        XYZP0(1,IPHAS) = 0.D0
        XYZP0(2,IPHAS) = 0.D0
        XYZP0(3,IPHAS) = HCHALI
        IXYZP0(IPHAS) = 1
C
C
C --- Vitesse de reference
C     pour le calcul de la viscosite turbulente et
C       pour les pertes de charge, le Richardson, ...
C     la vitesse de reference n'est pas utilisee pour le calcul du pas
C       de temps : le pas de temps est pris ici uniforme et constant ;
C       s'il etait envisage de le prendre variable en temps, pour
C       l'adapter automatiquement a la configuration finale du
C       calcul une fois l'etat stationnaire etabli, il conviendrait
C       alors d'utiliser la vitesse de reference pour calculer une
C       limite basee par exemple sur un nombre de Courant.
C
C
C   - En convection naturelle
C       l'echelle de vitesse est calculee par analyse dimensionnelle
C         comme VITREF = P / (E S), a partir :
C       . de l'energie (E, en Joule/m3) que l'air initialement a TINIT
C          peut extraire s'il atteint la temperature critique de
C          sortie TCRIT : E = RO0(IPHAS)*CP0(IPHAS)*(TSOR0-TINIT)
C       . de la puissance (P, en Joule/s) d'une ligne (x cst) de
C           NPLGRS conteneurs, P = PUICON * NPLGRS
C       . de la surface verticale (S en m2) de la coupe transverse
C           d'une ligne de conteneurs, S = PTRRES*EPCHEL*NCHEST
C
        IF(ICOFOR.EQ.0)THEN
C
          TSOR0  = TCRIT
          PRFMAT = PUICON * NPLGRS
          ERFMAT = RO0(IPHAS)*CP0(IPHAS)*(TSOR0-TINIT)
          SRFMAT = PTRRES*EPCHEL*NCHEST
C
          VITREF = PRFMAT/(ERFMAT*SRFMAT)
C
C
C   - En convection forcee
C       La vitesse de reference est calculee comme le rapport du debit
C         reel total (DEBMAS/RO0(IPHAS)) a la surface de la zone de
C         stockage (NPTRAN*PTRRES*NCHEST*EPCHEL).
C       Elle est corrigee par le facteur de reduction transverse
C         FRDTRA du maillage par rapport a la realite afin d'obtenir
C         une vitesse representative de la vitesse reelle.
C
        ELSE
C
          VITREF =
     &         DEBMAS/(RO0(IPHAS)*FRDTRA*NPTRAN*PTRRES*NCHEST*EPCHEL)
C
C TSOR0 inutile ici ; on conserve la formule pour memoire
C     TSOR0= TINIT + NPTRAN*NPLGRS*PUICON/(DEBMAS/FRDTRA)/CP0(IPHAS)
C
        ENDIF
C
C
C --- Viscosite dynamique totale (modelisee, constante et uniforme)
C
C     Le tableau VISCL0 contient la viscosite dynamique totale
C       (moleculaire + turbulente). On evalue d'abord la viscosite
C       turbulente puis on ajoute la viscosite moleculaire XMUMAT.
C
C     La viscosite dynamique turbulente est modelisee sous la forme
C       mu_t = rho * u*_ref * D_ref, formule dans laquelle :
C     . u*_ref est la vitesse de frottement deduite de la vitesse de
C              reference VITREF calculee ci-dessus et d'une intensite
C              turbulente imposee RTURB0
C     . D_ref  est une distance de reference basee sur l'encombrement
C              du milieu
C
C     La valeur de u*_ref est calculee a partir de l'energie cinetique
C       turbulente k par u*_ref = Cmu**(1/4) k**(1/2). La valeur de k
C       se deduit de l'intensite turbulente I, supposee connue, par
C       k = 3/2 (I V_ref)**2 avec V_ref la vitesse de reference
C       VITREF calculee precedemment.
C
C     La valeur de D_ref est calculee par D_ref = 0.2(Pt - d) ou Pt est
C       le pas transverse du reseau et d le diametre des conteneurs.
C
C     Ainsi, mu_t = rho * u*_ref * L_ref
C            . u*_ref = Cmu**(1/4) * (3/2)**(1/2) * I * V_ref
C            . L_ref  = 0.2 * (Pt - d)
C
C     On ajoute la viscosite moleculaire XMUMAT a mu_t pour obtenir
C       la viscosite dynamique totale VISCL0
C
C     VISCL0 est recalcule dans mtphyv (VITREF est mis a jour au cours
C       du calcul)
C
C
C   - Calcul de la viscosite turbulente
C
        USTMAT = CMU**0.25D0 * SQRT(1.5D0) * (RTURB0/100.D0) * VITREF
        DRFMAT = 0.2D0 * (PTRRES-DMCONT)
        VISCL0(IPHAS) = RO0(IPHAS) * USTMAT * DRFMAT
C
C   - Ajout de la viscosite moleculaire
C
        VISCL0(IPHAS) = VISCL0(IPHAS) + XMUMAT
C
C
C --- Pas de temps
C
C   - Pas de temps uniforme en espace et constant en temps
C       Initialement, on souhaitait adopter dans Matisse un pas de
C         temps variable en temps, mais l'option n'a jamais ete
C         effectivement utilisee.
C       Quels que soient les choix par defaut de Code_Saturne, on
C         conserve IDTVAR ici pour memoire, au cas ou l'on reprendrait
C         les tests plus tard.
C
        IDTVAR = 0
C
C   - Calcul de la porosite verticale liee a la presence d'un
C       conteneur :
C       POR0 = section horizontale passante /section horizontale totale
C       (on suppose le conteneur cylindrique à base circulaire,
C        d'axe vertical)
C       Utilise pour le calcul du pas de temps ci-dessous
C
        CXSD = PTRRES/DMCONT
        CYSD = PLGRES/DMCONT
        POR0 = 1.D0 - PI/(CXSD*CYSD*4.D0)
C
C   - On calcule le pas de temps par analyse dimensionnelle pour que
C       le volume d'air avoisinant un conteneur ne voie pas sa
C       temperature croitre de plus de DTDTMX degres en un pas de
C       temps.
C     On a donc DTREF = E * V / P avec :
C       . E = RO0(IPHAS)*CP0(IPHAS)*DTDTMX, en Joule/m3, l'energie
C           maximale que l'air peut recevoir par pas de temps,
C       . V = (PLGRES*PTRRES*EPCHEL*NCHEST)*POR0, en m3, le volume
C           d'air entourant un conteneur.
C       . P = PUICON, en Joule/s, la puissance d'un conteneur
C
        DTREF = RO0(IPHAS)*CP0(IPHAS)*DTDTMX
     &       * (PLGRES*PTRRES*EPCHEL*NCHEST)*POR0 / PUICON
C
C
C --- Diffusivite des scalaires
C     Une diffusivite de reference est fournie ici pour tous les
C       scalaires, mais la diffusivite est variable (voir mtphyv).
C     Il s'agit plutot d'une securite.
C     Pour la temperature de l'air ambiant (scalaire ITAAMT),
C       l'initialisation est ici de type "Prandtl=1" (on considere
C       un Prandtl "efficace" incluant la turbulence)
C     Pour la temperature de peau des parois (scalaires ITPPMT),
C       il n'y a pas de phenomene de diffusion pris en compte :
C       peu importe donc la valeur adoptee.
C
        IF(NSCAUS.GT.0) THEN
C
C       On boucle sur les scalaires utilisateurs :
          DO II = 1, NSCAUS
C         Pour les scalaires qui ne sont pas des variances
            IF(ISCAVR(II).LE.0) THEN
C           On definit la diffusivite
              VISLS0(II) = VISCL0(IPHSCA(II))
            ENDIF
          ENDDO
C
        ENDIF
C
C     Fin du test IF (IMATIS.EQ.1)
      ENDIF
C
C
C--------
C FORMATS
C--------
C
 9001 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : MATISSE                                     ',/,
     &'@    *********                                               ',/,
     &'@      ARRET SUR LE PARALLELISME DANS MTINI1                 ',/,
     &'@                                                            ',/,
     &'@    Le parallelisme n''est pas prevu avec Matisse           ',/,
     &'@      or le calcul semble avoir ete lance en parallele.     ',/,
     &'@                                                            ',/,
     &'@    Le calcul ne sera pas execute.                          ',/,
     &'@                                                            ',/,
     &'@    Contacter l''equipe de developpement.                   ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 9011 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET MATISSE (MTINI1)                      ',/,
     &'@    *********                                               ',/,
     &'@      ARRET SUR LE NOMBRE DE ZONES DES CARTES 2D ET 3D.     ',/,
     &'@                                                            ',/,
     &'@    Le nombre de zones demande est                 ',I10   ,  /,
     &'@    Le nombre maximal de zones autorise est NZONMX ',I10   ,  /,
     &'@                                                            ',/,
     &'@    Le calcul ne sera pas execute.                          ',/,
     &'@                                                            ',/,
     &'@    Contacter l''equipe de developpement.                   ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 9021 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET MATISSE (MTINI1)                      ',/,
     &'@    *********                                               ',/,
     &'@      NOMBRE DE SCALAIRES INCORRECT.                        ',/,
     &'@                                                            ',/,
     &'@    NSCAUS represente                                       ',/,
     &'@      le nombre de scalaires requis pour Matisse.           ',/,
     &'@                                                            ',/,
     &'@    NSCAUS doit etre exactement egal a 3.                   ',/,
     &'@                                                            ',/,
     &'@    Il vaut ici                                             ',/,
     &'@      NSCAUS = ',I10   ,'                                   ',/,
     &'@                                                            ',/,
     &'@    Le calcul ne sera pas execute.                          ',/,
     &'@                                                            ',/,
     &'@    Contacter l''equipe de developpement.                   ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 9031 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET MATISSE (MTINI1)                      ',/,
     &'@    *********                                               ',/,
     &'@      ECART TCRIT-TINIT NEGATIF OU NUL.                     ',/,
     &'@                                                            ',/,
     &'@    TCRIT represente                                        ',/,
     &'@      la temperature d''air de sortie critique en degres C. ',/,
     &'@    TINIT represente                                        ',/,
     &'@      la temperature d''air en entree          en degres C. ',/,
     &'@                                                            ',/,
     &'@    TCRIT est une estimation de la temperature de l''air    ',/,
     &'@      en sortie de l''entrepot.                             ',/,
     &'@    TCRIT doit donc naturellement etre strictement          ',/,
     &'@      superieure a TINIT.                                   ',/,
     &'@                                                            ',/,
     &'@    Or, les donnees saisies sont telles que                 ',/,
     &'@      TCRIT = ',E12.5   ,'                                  ',/,
     &'@      TINIT = ',E12.5   ,'                                  ',/,
     &'@                                                            ',/,
     &'@    Le calcul ne sera pas execute.                          ',/,
     &'@                                                            ',/,
     &'@    Modifier les donnes saisies pour que TCRIT > TINIT.     ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 9041 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET MATISSE (MTINI1)                      ',/,
     &'@    *********                                               ',/,
     &'@      RAPPORT DE CONVERGENT NON UNITE EN 2D.                ',/,
     &'@                                                            ',/,
     &'@    RCONVE represente                                       ',/,
     &'@      le rapport du convergent represente sur le maillage.  ',/,
     &'@                                                            ',/,
     &'@    RCONVE doit etre unite pour un calcul bidimensionnel.   ',/,
     &'@                                                            ',/,
     &'@    Ici, le calcul est apparemment bidimensionnel puisque   ',/,
     &'@      le facteur de reduction transverse saisi n''est pas   ',/,
     &'@      unite : FRDTRA = ',E12.5   ,'                         ',/,
     &'@    Or, la valeur saisie pour le rapport du convergent      ',/,
     &'@      n''est pas unite : RCONVE = ',E12.5   ,'              ',/,
     &'@                                                            ',/,
     &'@    Le calcul ne sera pas execute.                          ',/,
     &'@                                                            ',/,
     &'@    Modifier la valeur du facteur de reduction transverse   ',/,
     &'@      ou du rapport de convergent.                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 9051 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET MATISSE (MTINI1)                      ',/,
     &'@    *********                                               ',/,
     &'@      HAUTEUR D''EROSION DES PANACHES INADAPTEE.            ',/,
     &'@                                                            ',/,
     &'@    La modelisation des panaches est activee                ',/,
     &'@      avec IMDCNT = ',I10   ,'                              ',/,
     &'@    La hauteur d''erosion des panaches reduite a un nombre  ',/,
     &'@      entier de mailles doit etre strictement positive et   ',/,
     &'@      inferieure ou egale a la distance separant le plafond ',/,
     &'@      du sommet du reseau de colis.                         ',/,
     &'@                                                            ',/,
     &'@    La hauteur d''erosion saisie est   HEPCNT = ',E12.5      ,/,
     &'@                                                            ',/,
     &'@    La hauteur d''erosion reduite a un nombre entier de     ',/,
     &'@      mailles est                      HERCNT = ',E12.5      ,/,
     &'@    Le plafond est a la hauteur NCHEST*EPCHEL = ',E12.5      ,/,
     &'@      avec NCHEST = ',I10   ,' et EPCHEL = ',E12.5   ,'     ',/,
     &'@    La hauteur du reseau de colis est  HRESO  = ',E12.5      ,/,
     &'@                                                            ',/,
     &'@    L''inegalite suivante n''est pas verifiee :             ',/,
     &'@             0 <    HERCNT    <= (NCHEST*EPCHEL - HRESO)    ',/,
     &'@                 ',E12.5   ,'           ',E12.5   ,'        ',/,
     &'@                                                            ',/,
     &'@    Le calcul ne sera pas execute.                          ',/,
     &'@                                                            ',/,
     &'@    Descativer la modelisation des panaches ou modifier la  ',/,
     &'@      hauteur d''erosion prescrite.                         ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
C
      RETURN
      END
c@z
