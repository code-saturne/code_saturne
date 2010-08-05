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

! Module for Matisse

module matiss

  !===========================================================================

  ! --> PARAMETRES PHYSIQUES EN DUR

  !     TRFMAT : Temperature de reference pour l'air en degre
  !     RRFMAT : Masse volumique de reference de l'air a TRFMAT degres
  !     CRFMAT : CP de reference pour l'air
  !     XMUMAT : Viscosite moleculaire dynamique de l'air consideree
  !     RTURB0 : Intensite turbulente (selon k = (3/2)*(V_ref*RTURB0/100)**2)

  double precision trfmat, rrfmat, crfmat, xmumat, rturb0
  parameter(trfmat=20.d0)
  parameter(rrfmat=1.177d0)
  parameter(crfmat=1004.d0)
  parameter(xmumat=1.85d-5)
  parameter(rturb0=30.d0)

  ! --> INDICATEUR MATISSE

  !     IMATIS : Indicateur permettant de savoir si on utilise Matisse
  integer, save ::          imatis


  ! --> VARIABLES NUMERIQUES

  !     IMPMAT : Unite logique du fichier de resultats
  integer, save ::          impmat

  !     IICONR : "pointeur" sur ICONRA (connectivite rayonnement et panaches)
  integer, save ::          iiconr

  !     ICNROK : indique si la connectivite pour le rayonnement ICONRA
  !              a ete calculee ( = 1) ou non ( = 0)
  integer, save ::          icnrok

  ! --> DONNEES GEOMETRIQUES

  !   - Entiers

  !     NPTRAN : Nombre de pas d espace transversal
  !     NPLGRS : Nombre de pas d espace longitudinal
  !     NELGRS : Nombre d'elements par pas longitudinal
  !     NCHEST : Nombre de couche d element dans la zone stockage
  !     ITYPEN : Type d'entreposage
  !                * 1 : Emm
  !                * 0 : Vault

  !   - Reels

  !     EPREGI : Epaisseur des registres/cloisons amont et aval (en y)
  !     EPCHEM : Epaisseur des cheminees (en y)
  !     HCONVE : Hauteur du convergent eventuel
  !     RCONVE : Rapport du convergent eventuel sur le maillage (>=1)
  !     HCHALI : Hauteur de la cheminee d alimentation
  !     HCHEVA : Hauteur de la cheminee d evacuation
  !     HFTTOI : Hauteur du faite du toit
  !     PTRRES : Pas transversal du reseau de conteneur
  !     FRDTRA : Facteur de reduction transversal du modele/reel
  !     PLGRES : Pas longitudinal du reseau de conteneur
  !     EPCHEL : Epaisseur d une couche d element (zone stockage)
  !     DMCONT : Diametre des conteneurs
  !     HRESO  : Hauteur du reseau  de colis
  !     HPLEN  : Hauteur du plenum inferieur (cas alveole uniquement)

  integer, save ::          nptran, nplgrs, nelgrs, nchest, netran, itypen

  double precision, save :: epregi, epchem, hconve, rconve, hchali, hcheva
  double precision, save :: hfttoi, ptrres, frdtra, plgres, epchel, dmcont
  double precision, save :: hreso , hplen

  ! --> DONNEES PHYSIQUES

  !   - Entiers

  !     IMDCNT (0 ou 1) : Modelisation des panaches de convection naturelle
  !     ICOFOR (0 ou 1) : Regime hydraulique de circulation forcee
  !     ICONLG (0 ou 1) : Reseau de conteneur en ligne (pas triangulaire sinon)
  !     IALVEO (0 ou 1) : Entreposage en alveole

  !   - Reels

  !     DTDTMX : Delta temperature max / pas de temps
  !     PUICON : Puissance d'un conteneur
  !     TINIT  : Temperature d'air en entree en degres C
  !     TCRIT  : Temperature d'air de sortie critique en degres C
  !     EMICON : Emissivite des conteneurs
  !     EMIMUR : Emissivite des murs
  !     HEPCNT : Hauteur d erosion des panaches de convection naturelle
  !     DHPCNT : Debit enthalpique des panaches de convection naturelle
  !     DEBMAS : Debit de circulation forcee
  !     PDCCHA : Perte de charge du diffuseur de cheminee d'ALIMENTATION
  !     PDCFCH : Perte de charge du filtre de cheminee d'ALIMENTATION
  !     DHCHEA : Diametre hydraulique de cheminee d'ALIMENTATION
  !     SDCHEA : Surface debitante de cheminee d'ALIMENTATION
  !     PDCCHE : Perte de charge du diffuseur de cheminee d'EVACUATION
  !     PDCCCH : Perte de charge du clapet de cheminee d'EVACUATION
  !     DHCHES : Diametre hydraulique de cheminee d'EVACUATION
  !     SDCHES : Surface debitante de cheminee d'EVACUATION
  !     PDCALG : Perte de charge porte d'entree AMONT longitudinal
  !     PDCATV : Perte de charge porte d'entree AMONT transversale (sur z)
  !     ARGAMT : Angle d inclinaison du registre AMONT (degre)
  !     PDCSLG : Perte de charge porte de sortie AVAL longitudinale
  !     PDCSTV : Perte de charge porte de sortie AVAL transversale (sur z)
  !     ARGAVL : Angle d inclinaison du registre AVAL (degre)
  !     AMPPDC : Amplification des pertes de charge de reseau
  !     DHALVE : Diametre hydraulique de l'alveole
  !     VITREF : Vitesse de reference pour calculer les pertes de charge
  !     PUITOT : Puissance totale de l'installation
  !     DPVENT : Differentiel de pression athmospherique entree/sortie

  integer, save ::          imdcnt, icofor, iconlg, ialveo

  double precision, save :: dtdtmx, puicon, tinit , tcrit , emicon, emimur
  double precision, save :: hepcnt, dhpcnt, debmas, pdccha, pdcfch, dhchea
  double precision, save :: sdchea, pdcche, pdccch, dhches, sdches, pdcalg
  double precision, save :: pdcatv, argamt, pdcslg, pdcstv, argavl, amppdc
  double precision, save :: dhalve, vitref, puitot, debcon, cfecca, cfecma
  double precision, save :: dpvent

  ! --> CARTES 2D ET 3D

  !   - Dimensions

  !     NZONMX : nombre de zones maximum pour la definition des
  !               - Cartes 2D des pertes de charges de porte d'entree
  !               - Cartes 2D des pertes de charges de porte de sortie
  !               - Cartes 3D des pertes de charges reseau
  !               - Cartes 3D des puissances
  !     NCARTE : nombre de cartes 2D et 3D a definir
  !     NMTDIR : nombre de directions d'espace

  !   - Indicateurs de numero de carte (pour NZOCAR et VIZCAR)
  !       de 1 a NCARTE

  !     ICPDCE : indicateur carte 2D pdc de porte d'entree
  !     ICPDCS : indicateur carte 2D pdc de porte de sortie
  !     ICPDCR : indicateur carte 2D pdc de reseau
  !     ICPUIS : indicateur carte 3D puissance

  !   - Indicateurs de direction pour les cartes (NZOCAR, VIZCAR, VZOCAR)
  !       de 1 a NMTDIR

  !     ILIGNE : ligne    ( x variable, une ligne    = (y;z) constant )
  !     IRANGE : rangee   ( y variable, une rangee   = (x;z) constant )
  !     IALTIT : altitude ( z variable, une altitude = (x;y) constant )

  !   - Tableaux de donnees

  !     NZOCAR(NMTDIR, NCARTE)
  !            : nombre de zones pour chaque direction de chaque carte
  !     VIZCAR(2, NZONMX, NMTDIR, NCARTE)
  !            : definition du debut et de la fin des zones pour
  !                chaque zone dans chaque direction de chaque carte
  !              c'est un reel puisque l'on peut definir une demi ligne
  !                par exemple
  !     VZOCAR(NZONMX, NMTDIR)
  !            : valeurs associees a la carte 3D de puissance pour
  !                chaque zone dans chaque direction
  !              en x et y, valeurs comprises entre 0 et 1, indiquant
  !                la fraction de colis representee sur une maille
  !              en z, les valeurs sont renormalisees a l'unite

  !              il n'est plus necessaire de stocker les valeurs pour les
  !                cartes de pertes de charges car il n'y a que 2 valeurs
  !                possibles : 1. perte de charge        (colis present)
  !                            0. pas de perte de charge
  !              par defaut, on considere qu'il n'y a pas de perte de
  !                charge ; on repere avec VIZCAR le debut et la fin
  !                des zones dans lesquelles les pertes de charges doivent
  !                etre activees

  !   - Parametres des zones

  integer   nzonmx
  parameter(nzonmx=100)

  !   - Parametres des types de carte

  integer   ncarte
  parameter(ncarte=4)

  integer   icpdce, icpdcs, icpdcr, icpuis
  parameter(icpdce=1)
  parameter(icpdcs=2)
  parameter(icpdcr=3)
  parameter(icpuis=4)

  !   - Parametres des directions

  integer   nmtdir
  parameter(nmtdir=3)

  integer   iligne, irange, ialtit
  parameter(iligne=1)
  parameter(irange=2)
  parameter(ialtit=3)

  !   - Numero des cartes

  integer, save ::          nzocar(nmtdir, ncarte)

  !   - Valeur des cartes

  double precision, save :: vizcar(2, nzonmx, nmtdir, ncarte)
  double precision, save :: vcarth(nzonmx, nmtdir)

  ! --> NUMERO DES SCALAIRES

  !     ITAAMT : temperature air ambiant
  !     ITPCMT : temperature de peau des colis
  !     ITPPMT : temperature de peau des parois
  !                             (murs et alveoles eventuelles)
  !     A utiliser comme
  !       ISCA(ITAAMT), ISCA(ITPCMT), ISCA(ITPPMT)

  integer   itaamt, itpcmt, itppmt
  parameter(itaamt=1)
  parameter(itpcmt=2)
  parameter(itppmt=3)

  ! --> NUMERO DES COULEURS

  !   - Couleurs d'elements

  !     ICMTDF : couleur par defaut      ("df" pour "defaut")
  !     ICMTST : zone de stockage        ("st" pour "stockage")
  !     ICMTCI : cheminee d'alimentation ("ci" pour "cheminee inlet")
  !     ICMTCO : cheminee d'evacuation   ("co" pour "cheminee outlet")
  !     ICMTRI : registre amont          ("ri" pour "registre inlet")
  !     ICMTRO : registre aval           ("ro" pour "registre outlet")
  !     ICMTJI : jeu entre colis et registre amont
  !                                      ("ji" pour "jeu inlet")
  !     ICMTJO : jeu entre colis et registre aval
  !                                      ("jo" pour "jeu outlet")

  !     ICMTCI et ICMTCO designent la partie des cheminees situees
  !       au dessus des convergents enventuels
  !     ICMTJI et ICMTJO ne sont utilises que pour les configurations
  !       avec alveoles (et il n'y a pas de jeu aval en cathedrale)

  integer   icmtdf, icmtst
  integer   icmtci, icmtco, icmtri, icmtro, icmtji, icmtjo
  parameter(icmtdf =  0)
  parameter(icmtst =  8)
  parameter(icmtci =  3)
  parameter(icmtco =  6)
  parameter(icmtri =  2)
  parameter(icmtro =  4)
  parameter(icmtji =  7)
  parameter(icmtjo =  9)

  !   - Couleurs de faces (conditions aux limites)

  !     ICMTFI : entree   ("i" pour "in")
  !     ICMTFO : sortie   ("o" pour "out")
  !     ICMTFG : sol      ("g" pour "ground")
  !     ICMTFC : plafond  ("c" pour "ceiling")
  !     ICMTFS : symetrie ("s" pour "symmetry")
  !     ICMTFW : paroi    ("w" pour "wall")

  integer   icmtfi, icmtfo, icmtfg, icmtfc, icmtfs, icmtfw
  parameter(icmtfi =  5)
  parameter(icmtfo =  1)
  parameter(icmtfg = 10)
  parameter(icmtfc = 11)
  parameter(icmtfs = 12)
  parameter(icmtfw = 13)

  ! --> VARIABLES DEDUITES

  !     HERCNT : Hauteur d'erosion HEPCNT reduite a un nombre entier de
  !                mailles en altitude (ici pour eviter de faire le calcul
  !                plusieurs fois, et donc de risquer des erreurs)

  double precision, save :: hercnt

  !=============================================================================

end module matiss
