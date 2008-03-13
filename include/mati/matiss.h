c@a
c@versb
C-----------------------------------------------------------------------
C
CVERS                  Code_Saturne version 1.3
C                      ------------------------
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
C                             matiss.h
C
C***********************************************************************
C
C            INCLUDE POUR MATISSE
C
C-----------------------------------------------------------------------
C
C
C --> PARAMETRES PHYSIQUES EN DUR
C
C     TRFMAT : Temperature de reference pour l'air en degre
C     RRFMAT : Masse volumique de reference de l'air a TRFMAT degres
C     CRFMAT : CP de reference pour l'air
C     XMUMAT : Viscosite moléculaire dynamique de l'air consideree
C     RTURB0 : Intensite turbulente (selon k = (3/2)*(V_ref*RTURB0/100)**2)
C
      DOUBLE PRECISION TRFMAT, RRFMAT, CRFMAT, XMUMAT, RTURB0
      PARAMETER(TRFMAT=20.D0)
      PARAMETER(RRFMAT=1.177D0)
      PARAMETER(CRFMAT=1004.D0)
      PARAMETER(XMUMAT=1.85D-5)
      PARAMETER(RTURB0=30.D0)
C
C
C --> INDICATEUR MATISSE
C
C     IMATIS : Indicateur permettant de savoir si on utilise Matisse
      INTEGER          IMATIS
      COMMON /IMTMAT/  IMATIS
C
C
C --> VARIABLES NUMERIQUES
C
C     IMPMAT : Unite logique du fichier de resultats
      INTEGER          IMPMAT
      COMMON /IMTFIC/  IMPMAT
C
C     IICONR : "pointeur" sur ICONRA (connectivite rayonnement et panaches)
      INTEGER          IICONR
      COMMON /IMTMEM/  IICONR
C
C     ICNROK : indique si la connectivite pour le rayonnement ICONRA
C              a ete calculee ( = 1) ou non ( = 0)
      INTEGER          ICNROK
      COMMON /IMTRAY/  ICNROK
C
C
C --> DONNEES GEOMETRIQUES
C
C   - Entiers
C
C     NPTRAN : Nombre de pas d espace transversal
C     NPLGRS : Nombre de pas d espace longitudinal
C     NELGRS : Nombre d'elements par pas longitudinal
C     NCHEST : Nombre de couche d element dans la zone stockage
C     ITYPEN : Type d'entreposage
C                * 1 : Emm
C                * 0 : Vault
C
C   - Reels
C
C     EPREGI : Epaisseur des registres/cloisons amont et aval (en y)
C     EPCHEM : Epaisseur des cheminees (en y)
C     HCONVE : Hauteur du convergent eventuel
C     RCONVE : Rapport du convergent eventuel sur le maillage (>=1)
C     HCHALI : Hauteur de la cheminee d alimentation
C     HCHEVA : Hauteur de la cheminee d evacuation
C     HFTTOI : Hauteur du faite du toit
C     PTRRES : Pas transversal du reseau de conteneur
C     FRDTRA : Facteur de reduction transversal du modele/reel
C     PLGRES : Pas longitudinal du reseau de conteneur
C     EPCHEL : Epaisseur d une couche d element (zone stockage)
C     DMCONT : Diametre des conteneurs
C     HRESO  : Hauteur du reseau  de colis
C     HPLEN  : Hauteur du plenum inferieur (cas alveole uniquement)
C
      INTEGER          NPTRAN, NPLGRS, NELGRS, NCHEST, NETRAN
      INTEGER          ITYPEN
      COMMON /IMTGEO/  NPTRAN, NPLGRS, NELGRS, NCHEST, NETRAN,
     &                 ITYPEN
C
      DOUBLE PRECISION EPREGI, EPCHEM, HCONVE, RCONVE, HCHALI, HCHEVA
      DOUBLE PRECISION HFTTOI, PTRRES, FRDTRA, PLGRES, EPCHEL, DMCONT
      DOUBLE PRECISION HRESO , HPLEN
      COMMON /RMTGEO/  EPREGI, EPCHEM, HCONVE, RCONVE, HCHALI, HCHEVA,
     &                 HFTTOI, PTRRES, FRDTRA, PLGRES, EPCHEL, DMCONT,
     &                 HRESO , HPLEN
C
C
C --> DONNEES PHYSIQUES
C
C   - Entiers
C
C     IMDCNT (0 ou 1) : Modelisation des panaches de convection naturelle
C     ICOFOR (0 ou 1) : Regime hydraulique de circulation forcee
C     ICONLG (0 ou 1) : Reseau de conteneur en ligne (pas triangulaire sinon)
C     IALVEO (0 ou 1) : Entreposage en alveole
C
C   - Reels
C
C     DTDTMX : Delta temperature max / pas de temps
C     PUICON : Puissance d'un conteneur
C     TINIT  : Temperature d'air en entree en degres C
C     TCRIT  : Temperature d'air de sortie critique en degres C
C     EMICON : Emissivite des conteneurs
C     EMIMUR : Emissivite des murs
C     HEPCNT : Hauteur d erosion des panaches de convection naturelle
C     DHPCNT : Debit enthalpique des panaches de convection naturelle
C     DEBMAS : Débit de circulation forcee
C     PDCCHA : Perte de charge du diffuseur de cheminee d'ALIMENTATION
C     PDCFCH : Perte de charge du filtre de cheminee d'ALIMENTATION
C     DHCHEA : Diametre hydraulique de cheminee d'ALIMENTATION
C     SDCHEA : Surface debitante de cheminee d'ALIMENTATION
C     PDCCHE : Perte de charge du diffuseur de cheminee d'EVACUATION
C     PDCCCH : Perte de charge du clapet de cheminee d'EVACUATION
C     DHCHES : Diametre hydraulique de cheminee d'EVACUATION
C     SDCHES : Surface debitante de cheminee d'EVACUATION
C     PDCALG : Perte de charge porte d'entree AMONT longitudinal
C     PDCATV : Perte de charge porte d'entree AMONT transversale (sur z)
C     ARGAMT : Angle d inclinaison du registre AMONT (degre)
C     PDCSLG : Perte de charge porte de sortie AVAL longitudinale
C     PDCSTV : Perte de charge porte de sortie AVAL transversale (sur z)
C     ARGAVL : Angle d inclinaison du registre AVAL (degre)
C     AMPPDC : Amplification des pertes de charge de reseau
C     DHALVE : Diametre hydraulique de l'alveole
C     VITREF : Vitesse de reference pour calculer les pertes de charge
C     PUITOT : Puissance totale de l'installation
C     DPVENT : Differentiel de pression athmospherique entree/sortie
C
      INTEGER          IMDCNT, ICOFOR, ICONLG, IALVEO
      COMMON /IMTPHY/  IMDCNT, ICOFOR, ICONLG, IALVEO
C
      DOUBLE PRECISION DTDTMX, PUICON, TINIT , TCRIT , EMICON, EMIMUR
      DOUBLE PRECISION HEPCNT, DHPCNT, DEBMAS, PDCCHA, PDCFCH, DHCHEA
      DOUBLE PRECISION SDCHEA, PDCCHE, PDCCCH, DHCHES, SDCHES, PDCALG
      DOUBLE PRECISION PDCATV, ARGAMT, PDCSLG, PDCSTV, ARGAVL, AMPPDC
      DOUBLE PRECISION DHALVE, VITREF, PUITOT, DEBCON, CFECCA, CFECMA
      DOUBLE PRECISION DPVENT
      COMMON /RMTPHY/  DTDTMX, PUICON, TINIT , TCRIT , EMICON, EMIMUR,
     &                 HEPCNT, DHPCNT, DEBMAS, PDCCHA, PDCFCH, DHCHEA,
     &                 SDCHEA, PDCCHE, PDCCCH, DHCHES, SDCHES, PDCALG,
     &                 PDCATV, ARGAMT, PDCSLG, PDCSTV, ARGAVL, AMPPDC,
     &                 DHALVE, VITREF, PUITOT, DEBCON, CFECCA, CFECMA,
     &                 DPVENT
C
C
C --> CARTES 2D ET 3D
C
C   - Dimensions
C
C     NZONMX : nombre de zones maximum pour la definition des
C               - Cartes 2D des pertes de charges de porte d'entree
C               - Cartes 2D des pertes de charges de porte de sortie
C               - Cartes 3D des pertes de charges reseau
C               - Cartes 3D des puissances
C     NCARTE : nombre de cartes 2D et 3D a definir
C     NMTDIR : nombre de directions d'espace
C
C   - Indicateurs de numero de carte (pour NZOCAR et VIZCAR)
C       de 1 a NCARTE
C
C     ICPDCE : indicateur carte 2D pdc de porte d'entree
C     ICPDCS : indicateur carte 2D pdc de porte de sortie
C     ICPDCR : indicateur carte 2D pdc de reseau
C     ICPUIS : indicateur carte 3D puissance
C
C   - Indicateurs de direction pour les cartes (NZOCAR, VIZCAR, VZOCAR)
C       de 1 a NMTDIR
C
C     ILIGNE : ligne    ( x variable, une ligne    = (y;z) constant )
C     IRANGE : rangee   ( y variable, une rangee   = (x;z) constant )
C     IALTIT : altitude ( z variable, une altitude = (x;y) constant )
C
C   - Tableaux de donnees
C
C     NZOCAR(NMTDIR, NCARTE)
C            : nombre de zones pour chaque direction de chaque carte
C     VIZCAR(2, NZONMX, NMTDIR, NCARTE)
C            : definition du debut et de la fin des zones pour
C                chaque zone dans chaque direction de chaque carte
C              c'est un reel puisque l'on peut definir une demi ligne
C                par exemple
C     VZOCAR(NZONMX, NMTDIR)
C            : valeurs associees a la carte 3D de puissance pour
C                chaque zone dans chaque direction
C              en x et y, valeurs comprises entre 0 et 1, indiquant
C                la fraction de colis representee sur une maille
C              en z, les valeurs sont renormalisees a l'unite
C
C              il n'est plus necessaire de stocker les valeurs pour les
C                cartes de pertes de charges car il n'y a que 2 valeurs
C                possibles : 1. perte de charge        (colis present)
C                            0. pas de perte de charge
C              par defaut, on considere qu'il n'y a pas de perte de
C                charge ; on repere avec VIZCAR le debut et la fin
C                des zones dans lesquelles les pertes de charges doivent
C                etre activees
C
C   - Parametres des zones
C
      INTEGER   NZONMX
      PARAMETER(NZONMX=100)
C
C   - Parametres des types de carte
C
      INTEGER   NCARTE
      PARAMETER(NCARTE=4)
C
      INTEGER   ICPDCE, ICPDCS, ICPDCR, ICPUIS
      PARAMETER(ICPDCE=1)
      PARAMETER(ICPDCS=2)
      PARAMETER(ICPDCR=3)
      PARAMETER(ICPUIS=4)
C
C   - Parametres des directions
C
      INTEGER   NMTDIR
      PARAMETER(NMTDIR=3)
C
      INTEGER   ILIGNE, IRANGE, IALTIT
      PARAMETER(ILIGNE=1)
      PARAMETER(IRANGE=2)
      PARAMETER(IALTIT=3)
C
C   - Numero des cartes
C
      INTEGER          NZOCAR(NMTDIR, NCARTE)
      COMMON /IMTCAR/  NZOCAR
C
C   - Valeur des cartes
C
      DOUBLE PRECISION VIZCAR(2, NZONMX, NMTDIR, NCARTE)
      DOUBLE PRECISION VCARTH(NZONMX, NMTDIR)
      COMMON /RMTCAR/  VIZCAR, VCARTH
C
C
C --> NUMERO DES SCALAIRES
C
C     ITAAMT : temperature air ambiant
C     ITPCMT : temperature de peau des colis
C     ITPPMT : temperature de peau des parois
C                             (murs et alveoles eventuelles)
C     A utiliser comme
C       ISCA(ITAAMT), ISCA(ITPCMT), ISCA(ITPPMT)
C
      INTEGER   ITAAMT, ITPCMT, ITPPMT
      PARAMETER(ITAAMT=1)
      PARAMETER(ITPCMT=2)
      PARAMETER(ITPPMT=3)
C
C
C --> NUMERO DES COULEURS
C
C   - Couleurs d'elements
C
C     ICMTDF : couleur par defaut      ("df" pour "defaut")
C     ICMTST : zone de stockage        ("st" pour "stockage")
C     ICMTCI : cheminee d'alimentation ("ci" pour "cheminee inlet")
C     ICMTCO : cheminee d'evacuation   ("co" pour "cheminee outlet")
C     ICMTRI : registre amont          ("ri" pour "registre inlet")
C     ICMTRO : registre aval           ("ro" pour "registre outlet")
C     ICMTJI : jeu entre colis et registre amont
C                                      ("ji" pour "jeu inlet")
C     ICMTJO : jeu entre colis et registre aval
C                                      ("jo" pour "jeu outlet")
C
C     ICMTCI et ICMTCO designent la partie des cheminees situees
C       au dessus des convergents enventuels
C     ICMTJI et ICMTJO ne sont utilises que pour les configurations
C       avec alveoles (et il n'y a pas de jeu aval en cathedrale)
C
      INTEGER   ICMTDF, ICMTST
      INTEGER   ICMTCI, ICMTCO, ICMTRI, ICMTRO, ICMTJI, ICMTJO
      PARAMETER(ICMTDF =  0)
      PARAMETER(ICMTST =  8)
      PARAMETER(ICMTCI =  3)
      PARAMETER(ICMTCO =  6)
      PARAMETER(ICMTRI =  2)
      PARAMETER(ICMTRO =  4)
      PARAMETER(ICMTJI =  7)
      PARAMETER(ICMTJO =  9)
C
C   - Couleurs de faces (conditions aux limites)
C
C     ICMTFI : entree   ("i" pour "in")
C     ICMTFO : sortie   ("o" pour "out")
C     ICMTFG : sol      ("g" pour "ground")
C     ICMTFC : plafond  ("c" pour "ceiling")
C     ICMTFS : symetrie ("s" pour "symmetry")
C     ICMTFW : paroi    ("w" pour "wall")
C
      INTEGER   ICMTFI, ICMTFO, ICMTFG, ICMTFC, ICMTFS, ICMTFW
      PARAMETER(ICMTFI =  5)
      PARAMETER(ICMTFO =  1)
      PARAMETER(ICMTFG = 10)
      PARAMETER(ICMTFC = 11)
      PARAMETER(ICMTFS = 12)
      PARAMETER(ICMTFW = 13)
C
C
C --> VARIABLES DEDUITES
C
C     HERCNT : Hauteur d'erosion HEPCNT reduite a un nombre entier de
C                mailles en altitude (ici pour eviter de faire le calcul
C                plusieurs fois, et donc de risquer des erreurs)
C
      DOUBLE PRECISION HERCNT
      COMMON /RMTBID/  HERCNT
C
C FIN
c@z
