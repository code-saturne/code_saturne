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
                        SUBROUTINE MTIMPI
C                       *****************
C     -------------------------------------------------------------
C     -------------------------------------------------------------
C***********************************************************************
C  FONCTION  :
C  ---------
c@foncb
CFONC
CFONC          AFFICHAGE DES DONNEES DE CALCUL MATISSE
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
      INCLUDE "ihmpre.h"
      INCLUDE "matiss.h"
      INCLUDE "cstphy.h"
      INCLUDE "entsor.h"
      INCLUDE "optcal.h"
C
C***********************************************************************
C
C ARGUMENTS
C
C VARIABLES LOCALES
C
      CHARACTER*15     NAME
      INTEGER          KECHRG, KERGRS, KECLGR, KERGCH, KECIEL
      INTEGER          II    , IPHAS
      DOUBLE PRECISION DJECHR, DJERCL, DJECLR, DJERCH
      DOUBLE PRECISION HBDTOI
      DOUBLE PRECISION TSOR0 , DBM   , CXSD
      DOUBLE PRECISION UREEL , BETMAT, HRFMAT, RICHAR
C
C***********************************************************************
C=======================================================================
C
C
C=======================================================================
C RECAPITULATIF DANS LE FICHIER resuMatisse.txt
C=======================================================================
C
C
C --- Ouverture du fichier
C       (il sera ferme au dernier pas de temps par mtproj)
C
      NAME='resuMatisse'
      OPEN(UNIT=IMPMAT,FILE=NAME,
     &         FORM='FORMATTED', STATUS='UNKNOWN', ERR=900)
      GOTO 950
C
C   - En cas d'erreur : message et stop
C
 900  WRITE (0, 9998) NAME
      WRITE (NFECRA, 9999) NAME
      CALL CSEXIT (1)
C     ===========
 950  CONTINUE
C
C
C --- Ecriture des donnees geometriques
C
C     On recupere les donnees non stockees dans les COMMON
C       et utiles seulement pour l'affichage
C
      CALL CSMHDB(DJECHR, DJERCL, DJECLR, DJERCH,
     &            KECHRG, KERGRS, KECLGR, KERGCH,
     &            HBDTOI, KECIEL)
C
C     On imprime
C
      WRITE(IMPMAT,1001)
      WRITE(IMPMAT,1002)
      WRITE(IMPMAT,1001)
C
      IF (ITYPEN.EQ.1) THEN
         WRITE(IMPMAT,2096) ' Emm  '
      ELSE
         WRITE(IMPMAT,2096) ' Vault'
      ENDIF
C
      IF (IALVEO.EQ.0) THEN
         WRITE(IMPMAT,2097) HRESO
      ELSE
         WRITE(IMPMAT,2098) HRESO
         WRITE(IMPMAT,2099) HPLEN
      ENDIF
C
      WRITE(IMPMAT,2100) EPREGI
      WRITE(IMPMAT,2101) EPCHEM
      WRITE(IMPMAT,2102) DJECHR
      WRITE(IMPMAT,2103) KECHRG
      WRITE(IMPMAT,2104) DJERCL
      WRITE(IMPMAT,2105) KERGRS
      IF (ITYPEN.EQ.0) THEN
        WRITE(IMPMAT,2106) DJECLR
        WRITE(IMPMAT,2107) KECLGR
        WRITE(IMPMAT,2108) DJERCH
        WRITE(IMPMAT,2109) KERGCH
      ENDIF
      IF (ITYPEN.EQ.1) THEN
        WRITE(IMPMAT,2124) HBDTOI
        WRITE(IMPMAT,2125) HFTTOI
        WRITE(IMPMAT,2126) KECIEL
      ENDIF
C
      WRITE(IMPMAT,2110) HCONVE
      WRITE(IMPMAT,2111) RCONVE
      WRITE(IMPMAT,2112) HCHALI
      WRITE(IMPMAT,2113) HCHEVA
      WRITE(IMPMAT,2114) PTRRES
      WRITE(IMPMAT,2115) NPTRAN
      WRITE(IMPMAT,2116) NETRAN
      WRITE(IMPMAT,2117) FRDTRA
      WRITE(IMPMAT,2118) PLGRES
      WRITE(IMPMAT,2119) NPLGRS
      WRITE(IMPMAT,2120) NELGRS
      WRITE(IMPMAT,2121) EPCHEL
      WRITE(IMPMAT,2122) NCHEST
      WRITE(IMPMAT,2123) DMCONT
C
C
C --- Ecriture des parametres de calcul
C
      WRITE(IMPMAT,*)
      WRITE(IMPMAT,1001)
      WRITE(IMPMAT,1003)
      WRITE(IMPMAT,1001)
C
      WRITE(IMPMAT,3099) NTMABS
      WRITE(IMPMAT,3100) IPHYDR
      WRITE(IMPMAT,1151)
      WRITE(IMPMAT,3101) DTDTMX
C
C
C --- Chargement thermique
C
      WRITE(IMPMAT,*)
      WRITE(IMPMAT,1001)
      WRITE(IMPMAT,1004)
      WRITE(IMPMAT,1001)
C
      WRITE(IMPMAT,3102) IMDCNT
      WRITE(IMPMAT,1151)
      WRITE(IMPMAT,3103) PUICON
      WRITE(IMPMAT,3104) TINIT
      WRITE(IMPMAT,3105) TCRIT
      WRITE(IMPMAT,3106) EMICON
      WRITE(IMPMAT,3107) EMIMUR
C
C
C --- Chargement hydraulique
C
      WRITE(IMPMAT,*)
      WRITE(IMPMAT,1001)
      WRITE(IMPMAT,1005)
      WRITE(IMPMAT,1001)
C
      WRITE(IMPMAT,4107) ICOFOR
      IF (ICONLG.EQ.1) THEN
        WRITE(IMPMAT,4108) ' En ligne        '
      ELSE
        WRITE(IMPMAT,4108) ' Pas triangulaire'
      ENDIF
      WRITE(IMPMAT,4109) IALVEO
      WRITE(IMPMAT,4110) DEBMAS
      WRITE(IMPMAT,4111) PDCCHA
      WRITE(IMPMAT,4112) PDCFCH
      WRITE(IMPMAT,4113) DHCHEA
      WRITE(IMPMAT,4114) SDCHEA
      WRITE(IMPMAT,4115) PDCCHE
      WRITE(IMPMAT,4116) PDCCCH
      WRITE(IMPMAT,4117) DHCHES
      WRITE(IMPMAT,4118) SDCHES
      WRITE(IMPMAT,4119) PDCALG
      WRITE(IMPMAT,4120) PDCATV
      WRITE(IMPMAT,4121) ARGAMT
      WRITE(IMPMAT,4122) PDCSLG
      WRITE(IMPMAT,4123) PDCSTV
      WRITE(IMPMAT,4124) ARGAVL
      WRITE(IMPMAT,4125) AMPPDC
      WRITE(IMPMAT,4126) DPVENT
      IF (IALVEO.EQ.1) THEN
        WRITE(IMPMAT,4127) DHALVE
      ENDIF
C
C
C --- Pertes de charge : cartes
C
      WRITE(IMPMAT,*)
      WRITE(IMPMAT,1001)
      WRITE(IMPMAT,1006)
      WRITE(IMPMAT,1001)
C
      WRITE(IMPMAT,5100)
      WRITE(IMPMAT,5000)
      DO II = 1, NZOCAR(IRANGE,ICPDCE)
        WRITE(IMPMAT,5003) VIZCAR(1, II, IRANGE, ICPDCE),
     &                     VIZCAR(2, II, IRANGE, ICPDCE)
      ENDDO
      WRITE(IMPMAT,5001)
      DO II = 1, NZOCAR(ILIGNE,ICPDCE)
        WRITE(IMPMAT,5003) VIZCAR(1, II, ILIGNE, ICPDCE),
     &                     VIZCAR(2, II, ILIGNE, ICPDCE)
      ENDDO
      WRITE(IMPMAT,5002)
      DO II = 1, NZOCAR(IALTIT,ICPDCE)
        WRITE(IMPMAT,5003) VIZCAR(1, II, IALTIT, ICPDCE),
     &                     VIZCAR(2, II, IALTIT, ICPDCE)
      ENDDO
C
      WRITE(IMPMAT,5101)
      WRITE(IMPMAT,5000)
      DO II = 1, NZOCAR(IRANGE,ICPDCS)
        WRITE(IMPMAT,5003) VIZCAR(1, II, IRANGE, ICPDCS),
     &                     VIZCAR(2, II, IRANGE, ICPDCS)
      ENDDO
      WRITE(IMPMAT,5001)
      DO II = 1, NZOCAR(ILIGNE,ICPDCS)
        WRITE(IMPMAT,5003) VIZCAR(1, II, ILIGNE, ICPDCS),
     &                     VIZCAR(2, II, ILIGNE, ICPDCS)
      ENDDO
      WRITE(IMPMAT,5002)
      DO II = 1, NZOCAR(IALTIT,ICPDCS)
        WRITE(IMPMAT,5003) VIZCAR(1, II, IALTIT, ICPDCS),
     &                     VIZCAR(2, II, IALTIT, ICPDCS)
      ENDDO
C
      WRITE(IMPMAT,5102)
      WRITE(IMPMAT,5000)
      DO II = 1, NZOCAR(IRANGE,ICPDCR)
        WRITE(IMPMAT,5003) VIZCAR(1, II, IRANGE, ICPDCR),
     &                     VIZCAR(2, II, IRANGE, ICPDCR)
      ENDDO
      WRITE(IMPMAT,5001)
      DO II = 1, NZOCAR(ILIGNE,ICPDCR)
        WRITE(IMPMAT,5003) VIZCAR(1, II, ILIGNE, ICPDCR),
     &                     VIZCAR(2, II, ILIGNE, ICPDCR)
      ENDDO
C
C
C --- Sources de chaleur homogeneisees : cartes
C
      WRITE(IMPMAT,*)
      WRITE(IMPMAT,1001)
      WRITE(IMPMAT,1007)
      WRITE(IMPMAT,1001)
C
      WRITE(IMPMAT,5000)
      DO II = 1, NZOCAR(IRANGE,ICPUIS)
        WRITE(IMPMAT,5004) VIZCAR(1, II, IRANGE, ICPUIS),
     &                     VIZCAR(2, II, IRANGE, ICPUIS),
     &                     VCARTH(II, IRANGE)
C
      ENDDO
      WRITE(IMPMAT,5001)
      DO II = 1, NZOCAR(ILIGNE,ICPUIS)
        WRITE(IMPMAT,5004) VIZCAR(1, II, ILIGNE, ICPUIS),
     &                     VIZCAR(2, II, ILIGNE, ICPUIS),
     &                     VCARTH(II, ILIGNE)
C
      ENDDO
      WRITE(IMPMAT,5002)
      DO II = 1, NZOCAR(IALTIT,ICPUIS)
        WRITE(IMPMAT,5004) VIZCAR(1, II, IALTIT, ICPUIS),
     &                     VIZCAR(2, II, IALTIT, ICPUIS),
     &                     VCARTH(II, IALTIT)
C
      ENDDO
C
C
C --- Affichage banniere de resultats
C
      WRITE(IMPMAT,*)
      WRITE(IMPMAT,1001)
      WRITE(IMPMAT,1008)
      WRITE(IMPMAT,1001)
C
C
C --- Affichage Richardson en convection forcee
C       selon la formule Ri = g beta DeltaT H_ref / U**2 avec :
C       . g = -GZ gravite verticale, valeur positive
C       . beta = 1/T en gaz parfait, ici 1/((TINIT+TSOR0)*0.5 + TKELVI)
C       . DeltaT = TSOR0-TINIT
C       . H_ref = EPCHEL*NCHEST hauteur de la zone de stockage
C       . U = UREEL vitesse reelle horizontale dans la zone de stockage
C
C       en outre :
C       . TSOR0 = temperature de sortie evaluee a partir de la
C           puissance des conteneurs, ie telle que on ait equilibre entre
C           debit * Cp (TSOR0-TINIT) en Joule/s et
C           NPTRAN*NPLGRS*PUICON en Joule/s
C       . le debit est calcule a partir du debit reel en corrigeant
C           par le facteur d'echelle transverse du modele par rapport au
C           cas reel FRDTRA.
C
C     En convection naturelle, le Richardson est affiche par mttsns
C
      IF (ICOFOR.EQ.1) THEN
C
        IPHAS  = 1
        DBM    = DEBMAS/FRDTRA
        TSOR0  = TINIT + NPTRAN*NPLGRS*PUICON/(DBM*CP0(IPHAS))
C
        CXSD   = PTRRES/DMCONT
        UREEL  = VITREF/(CXSD-1.D0)*CXSD
C
        BETMAT = 1.D0/((TINIT+TSOR0)*0.5D0 + TKELVI)
C
        HRFMAT = EPCHEL*NCHEST
C
        RICHAR = -GZ*BETMAT*(TSOR0-TINIT)*HRFMAT/(UREEL**2)
        WRITE(IMPMAT,1011) RICHAR
C
      ENDIF
C
C
C --- Des resultats seront affiches par mttsns et mtproj
C       le fichier n'est donc pas ferme
C
C=======================================================================
C FORMATS
C=======================================================================
C
C --- Texte
C
 1001 FORMAT(74('-'))
 1002 FORMAT(7(' '),'Geometrie :')
 1003 FORMAT(7(' '),'Parametres de calcul :')
 1004 FORMAT(7(' '),'Chargement thermique :')
 1005 FORMAT(7(' '),'Chargement hydraulique :')
 1006 FORMAT(7(' '),'Zones de pertes de charges :')
 1007 FORMAT(7(' '),'Repartition des sources de chaleurs',
     &  ' homogeneisees:')
 1008 FORMAT(7(' '),'Resultats :')
 1011 FORMAT(' Nombre de Richardson                                   ',
     & ':', E12.5)
C
 1151 FORMAT(' (1 oui, 0 non)')
C
C --- Geometrie
C
 2096 FORMAT(' Concept d''entrepot                                   ',
     & '  :',A6)
 2097 FORMAT(' Hauteur du réseau de conteneurs                       ',
     & ' :',E12.5,' m')
 2098 FORMAT(' Hauteur max des alvéoles                              ',
     & ' :',E12.5,' m')
 2099 FORMAT(' Hauteur min des alvéoles                              ',
     & ' :',E12.5,' m')
 2100 FORMAT(' Epaisseur des registres/cloisons amont et aval        ',
     & ' :',E12.5,' m')
 2101 FORMAT(' Epaisseur des cheminees                               ',
     & ' :',E12.5,' m')
 2102 FORMAT(' Jeu amont entre cheminee/registre                     ',
     & ' :',E12.5,' m')
 2103 FORMAT(' Nb d''elements entre cheminee/registre amont          ',
     & '  :',I12)
 2104 FORMAT(' Jeu entre registre amont/colis                        ',
     & ' :',E12.5,' m')
 2105 FORMAT(' Nb d''elements entre registre amont/reseau de colis   ',
     & '  :',I12)
 2106 FORMAT(' Jeu entre colis/registre aval                         ',
     & ' :',E12.5,' m')
 2107 FORMAT(' Nb d''elements entre colis/registre aval              ',
     & '  :',I12)
 2108 FORMAT(' Jeu aval entre registre/cheminee                      ',
     & ' :',E12.5,' m')
 2109 FORMAT(' Nb d''elements entre registre/cheminee aval           ',
     & '  :',I12)
 2110 FORMAT(' Hauteur du convergent                                 ',
     & ' :',E12.5,' m')
 2111 FORMAT(' Rapport du convergent                                 ',
     & ' :',E12.5,' m')
 2112 FORMAT(' Hauteur de la cheminee d''alimentation                ',
     & '  :',E12.5,' m')
 2113 FORMAT(' Hauteur de la cheminee d''evacuation                  ',
     & '  :',E12.5,' m')
 2114 FORMAT(' Pas transversal du reseau de conteneur                ',
     & ' :',E12.5,' m')
 2115 FORMAT(' Nombre de pas d''espace transversal                   ',
     & '  :',I12)
 2116 FORMAT(' Nombre d''elements par pas transversal                ',
     & '  :',I12)
 2117 FORMAT(' Facteur de reduction transversal du modele/reel       ',
     & ' :',E12.5)
 2118 FORMAT(' Pas longitudinal du reseau de conteneur               ',
     & ' :',E12.5,' m')
 2119 FORMAT(' Nombre de pas d''espace longitudinal                  ',
     & '  :',I12)
 2120 FORMAT(' Nombre d''elements par pas longitudinal               ',
     & '  :',I12)
 2121 FORMAT(' Epaisseur d''une couche d''element (zone stockage)    ',
     & '   :',E12.5,' m')
 2122 FORMAT(' Nombre de couche d''element dans la zone stockage     ',
     & '  :',I12)
 2123 FORMAT(' Diametre des conteneurs                               ',
     & ' :',E12.5,' m')
 2124 FORMAT(' Hauteur du bord de toit                               ',
     & ' :',E12.5,' m')
 2125 FORMAT(' Hauteur du faite de toit                              ',
     & ' :',E12.5,' m')
 2126 FORMAT(' Nombre de couches d''elements du ciel d''entrepot     ',
     & '   :',I12)
C
C --- Parametres de calcul
C
 3099 FORMAT(' Nombre de pas de temps                                ',
     & ' :',I12)
 3100 FORMAT(' Prise en compte de la pression hydrostatique          ',
     & ' :',I12)
 3101 FORMAT(' delta temperature max/pas de temps                    ',
     & ' :',E12.5)
 3102 FORMAT(' Modelisation des panaches de convection naturelle     ',
     & ' :',I12)
 3103 FORMAT(' Puissance d''un conteneur                             ',
     & '  :',E12.5,' W')
 3104 FORMAT(' Temperature d''air en entree                          ',
     & '  :',E12.5,' °C')
 3105 FORMAT(' Temperature d''air de sortie critique                 ',
     & '  :',E12.5,' °C')
 3106 FORMAT(' Emissivite des conteneurs                             ',
     & ' :',E12.5)
 3107 FORMAT(' Emissivite des murs                                   ',
     & ' :',E12.5)
C
C --- Chargement hydraulique
C
 4107 FORMAT(' Regime hydraulique de circulation forcee (1 oui, 0 non',
     & '):',I12)
 4108 FORMAT(' Reseau de conteneur                                   ',
     & ' :',A17)
 4109 FORMAT(' Entreposage en alveole (1 oui, 0 non)                 ',
     & ' :',I12)
 4110 FORMAT(' Debit de circulation forcee                           ',
     & ' :',E12.5, ' kg/s')
 4111 FORMAT(' Perte de charge du diffuseur de cheminee d''alimentati',
     & 'on:',E12.5)
 4112 FORMAT(' Perte de charge du filtre de cheminee d''alimentation ',
     & '  :',E12.5)
 4113 FORMAT(' Diametre hydraulique de la cheminee d''alimentation   ',
     & '  :',E12.5, ' m')
 4114 FORMAT(' Surface debitante de la cheminee d''alimentation      ',
     & '  :',E12.5, ' m²')
 4115 FORMAT(' Perte de charge du diffuseur de cheminee d''evacuation',
     & '  :',E12.5)
 4116 FORMAT(' Perte de charge du clapet de cheminee d''evacuation   ',
     & '  :',E12.5)
 4117 FORMAT(' Diametre hydraulique de la cheminee d''evacuation     ',
     & '  :',E12.5, ' m')
 4118 FORMAT(' Surface debitante de la cheminee d''evacuation        ',
     & '  :',E12.5, ' m²')
 4119 FORMAT(' Perte de charge porte d''entree AMONT longitudinale   ',
     & '  :',E12.5)
 4120 FORMAT(' Perte de charge porte d''entree AMONT transversale    ',
     & '  :',E12.5)
 4121 FORMAT(' Angle d''inclinaison du registre AMONT (degre)        ',
     & '  :',E12.5, ' °')
 4122 FORMAT(' Perte de charge porte de sortie AVAL longitudinale    ',
     & ' :',E12.5)
 4123 FORMAT(' Perte de charge porte de sortie AVAL transversale     ',
     & ' :',E12.5)
 4124 FORMAT(' Angle d''inclinaison du registre AMONT (degre)        ',
     & '  :',E12.5, ' °')
 4125 FORMAT(' Amplificateur des pertes de charge de reseau          ',
     & ' :',E12.5)
 4126 FORMAT(' Differentiel de pression entree/sortie                ',
     & ' :',E12.5, ' Pa')
 4127 FORMAT(' Diametre hydraulique de l''alveole                    ',
     & '  :',E12.5, ' m')
C
C
C --- Cartes
C
 5000 FORMAT('   - Direction longitudinale',
     &       ' (distance en nombre de rangees)')
 5001 FORMAT('   - Direction transversale',
     &       ' (distance en nombre de lignes)')
 5002 FORMAT('   - Direction verticale',
     &       ' (hauteur en nombre de mailles)')
C
 5003 FORMAT('       * min :', E12.5,' max :', E12.5)
 5004 FORMAT('       * min :', E12.5,' max :', E12.5,' val :', E12.5)
C
 5100 FORMAT(' Entree')
 5101 FORMAT(' Sortie')
 5102 FORMAT(' Réseau')
C
C --- Erreurs
C
 9998 FORMAT(/,
     &'Code_Saturne : Erreur d''initialisation :',/,
     &'Impossible d''ouvrir le fichier : ',A,/)
 9999 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET MATISSE (MTIMPI)                      ',/,
     &'@    *********                                               ',/,
     &'@      ARRET A L''OUVERTURE DU FICHIER RESUME MATISSE        ',/,
     &'@                                                            ',/,
     &'@  Le fichier ',A15,' ne peut etre ouvert.                   ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
C
      RETURN
      END
c@z
