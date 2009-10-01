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

!                              ppincl.h

!===============================================================================

!            INCLUDE GENERAL PROPRE A LA PHYSIQUE PARTICULIERE

! necessite ppthch et ppppar
!-------------------------------------------------------------------------------


!--> TABLEAU INDICATEURS DU CHOIX DE LA PHYSIQUE PARTICULIERE CHOISIE

integer   nmodmx
parameter(nmodmx = 50)
integer           ippmod(nmodmx)
common / iippmd / ippmod

! ---- Indicateur global de physique particuliere
!        IPPMOD(IPHPAR) = 0 : pas de physique particuliere
!                         1 : physique particuliere enclenchee
!                         2 : physique particuliere avec pilotage du
!                             rayonnement par fichier parametrique
integer iphpar

! ---- Modeles propres a la combustion gaz ICO...
integer icod3p, icodeq, icoebu, icobml, icolwc

! ---- Modeles propres a la combustion charbon pulverise ICP...
integer icp3pl

! ---- Modeles propres aux versions effet Joule et conduction ionique
integer  ieljou, ielarc, ielion

! ---- Modeles propres a la combustion charbon pulverise couplee Lagrangien
integer icpl3c
! ---- Modeles propres a la combustion fuel
integer icfuel

! ---- Modele compressible
integer  icompf

! ---- Modele atmospherique
integer  iatmos

! ---- Modele aerorefrigerants
integer  iaeros

parameter       (iphpar = 1 , icod3p = 2 , icodeq = 3 ,           &
                 icoebu = 4 , icobml = 5 , icolwc = 6 ,           &
                 icp3pl = 7 , icpl3c = 8 , icfuel = 9 ,           &
                 ieljou = 10, ielarc = 11, ielion = 12,           &
                 icompf = 13, iatmos = 14, iaeros = 15)


!--> NOMBRE DE VARIABLES ALGEBRIQUES OU D'ETAT
!    pour la physique particuliere NSALPP
!    total NSALTO
integer nsalpp, nsalto


!--> POINTEURS VARIABLES COMBUSTION GAZ

! ---- Variables transportees
integer ifm, ifp2m, iygfm, icm, icp2m, ifpcpm
integer iyfm, iyfp2m, icoyfp

! ---- Variables d'etat
integer iym(ngazgm), itemp, ifmin, ifmax
integer ickabs, it4m, it3m

! --- Pointeurs proprietes (PROPCE)
integer itsc

!--> POINTEURS VARIABLES COMBUSTION CHARBON PULVERISE

! ---- Variables transportees
!        Phase continue (melange gazeux)
integer if1m(ncharm), if2m(ncharm)
integer if3m, if4m, if4p2m, if5m, if6m, if7m
integer if3mc2
!        Phase dispersee (classe de particules)
integer ixck(nclcpm), ixch(nclcpm), inp(nclcpm)
integer ih2(nclcpm) , ixwt(nclcpm)

! ---- Variables d'etat
!        Phase continue (melange gazeux)
integer iym1(ngazem), itemp1 , irom1 , immel
!        Phase dispersee (classes de particules)
integer itemp2(nclcpm), irom2(nclcpm), idiam2(nclcpm)
integer ix2(nclcpm)
integer igmdch(nclcpm), igmhet(nclcpm) , ighco2(nclcpm)
integer igmdv1(nclcpm), igmdv2(nclcpm)
integer igmsec(nclcpm)

!--> POINTEURS VARIABLES COMBUSTION FUEL

! ---- Variables transportees
!        Phase continue
 integer ifvap, ifhtf
!        Phase dispersee
 integer ihlf(nclcpm)
 integer ixkf(nclcpm), ixfol(nclcpm), ing(nclcpm)
 integer itemp3(nclcpm), irom3(nclcpm), idiam3(nclcpm)
 integer ix3(nclcpm)

! ---- Variables d'etat
!        Phase continue
 integer iyfol(nclcpm)
!        Phase dispersee
 integer ih1hlf(nclcpm), igmhtf(nclcpm), igmeva(nclcpm)

!--> POINTEURS VARIABLES VERSION ELECTRIQUES

!        Dimension des 'vecteurs electriques'      = NDIMVE
integer    ndimve
parameter (ndimve = 3)

! ---- Variables transportees
!        Potentiel reel       = IPOTR
!        Potentiel imaginaire = IPOTI
!        Composantes du potentiel vecteur magnetique = IPOTVA()
!        Fraction massique des constituants = IYCOEL()

integer ipotr, ipoti, ipotva(ndimve), iycoel(ngazgm)

! ---- Variables d'etat dans PROPCE
!        Puissance volumique dissipee par effet Joule W/m3 = IEFJOU
!        Forces electromagnetiques de Laplace en N/m3      = ILAPLA()
!        Charge electrique volumique C/m3                  = IQELEC
!        Densite de courant electrique reelle A/m2         = IDJR()
!        Puissance volumique rayonnee W/m3
!        Densite de courant electrique imaginaire en A/m2  = IDJI()
!         ou coeff d'absorption en m-1                     = IDRAD

!      Les autres variables deduites seront locales pour
!        economiser de la place memoire.
!        Cela n'empeche pas de les sortir en post-traitement.
!        Gradient du potentiel reel en V/m                 = IGPOTR()
!        Gradient du potentiel imaginaire en V/m           = IGPOTI()
!        Composantes du champ magnetique en Tesla          = IBMG()


integer iefjou , iqelec , ilapla(ndimve)
integer idjr(ndimve)    , idji(ndimve)   , idrad

!--> POINTEURS VARIABLES CONDUCTION IONIQUE

integer   nesiom
parameter (nesiom = 10)
integer   nespio

! ---- Variables transportees
!        par espece
integer iymion(nesiom)

! ---- Variables d'etat


!--> POINTEURS COMPRESSIBLE

! ---- Variables transportees par phase
integer irho(nphsmx), ienerg(nphsmx), itempk(nphsmx)
! ---- Proprietes supplementaires par phase
integer icv(nphsmx), iviscv(nphsmx), ieos(nphsmx)

!     COMMON complete plus bas

! ---- Aliases pour les conditions aux limites
integer irun(nphsmx), irunh(nphsmx)

common / iaclcf / irun , irunh

! ---- Proprietes supplementaires par phase
double precision cv0(nphsmx), viscv0(nphsmx)

common / rpocfp / cv0 , viscv0

! ---- Prediction de pression par une equation d'evolution
integer ippred(nphsmx)
! ---- Flux de masse specifique pour la vitesse
integer iflmau(nphsmx)
! ---- Utilisation de la pression predite pour resoudre Navier-Stokes
integer igrdpp(nphsmx)
! --- Conditions aux limites prenant en compte l'equilibre hydrostatique
integer icfgrp(nphsmx)

common / ipocfo / ippred , iflmau , igrdpp , icfgrp

! ---- Flux de bord convectifs QDM et energie (numero de PROPFB)
integer           ifbrhu(nphsmx) , ifbrhv(nphsmx) ,               &
                  ifbrhw(nphsmx) , ifbene(nphsmx)

common / ipobcf/  ifbrhu         , ifbrhv         ,               &
                  ifbrhw         , ifbene


!--> POINTEURS AEROREFRIGERANTS

! ---- Variables transportees
integer itemp4, ihumid


!--> POINTEUR RELATIF A LA VARIABLE ENTHALPIE

integer ihm


!--> REMPLISSAGE COMMON POINTEURS VARIABLES TRANSPORTEES
!                                 VARIABLES D'ETAT

common / ipovst /                                                 &

! ---- Combustion gaz
                  ifm, ifp2m, iygfm, icm, icp2m, ifpcpm,          &
                  iyfm, iyfp2m, icoyfp,                           &
! ---- Combustion charbon pulverise
                  if1m, if2m, if3m, if4m, if4p2m, if5m,           &
                  if6m, if7m, if3mc2 ,                            &
                  ixck, ixch, inp , ih2 , ixwt  ,                 &
! ---- Combustion fuel
                  ihlf, ifvap, ifhtf,                             &
                  ixkf, ixfol, ing, ix3 ,                         &

! ---- Versions electriques
                  ipoti, ipotr, ipotva, iycoel,                   &
! ---- Conduction ionique
                  nespio, iymion,                                 &

! ---- Compressible
                  irho   , ienerg , itempk ,                      &
                  icv    , iviscv , ieos   ,                      &
! ---- Enthalpie
                  ihm

common / ipovsa /                                                 &

! ---- Nb de variables d'etat ou algebriques
                  nsalpp, nsalto,                                 &
! ---- Combustion gaz
                  iym, itemp, ifmin, ifmax, ickabs, it4m,         &
                                                    it3m,         &
! ---- Combustion charbon pulverise
                  iym1, itemp1, irom1 ,immel,                     &
                  itemp2, irom2, idiam2, ix2,                     &
                  igmdch, igmdv1, igmdv2, igmhet, ighco2 ,        &
                  igmsec,                                         &
! ---- Combustion fuel
                  iyfol, itemp3, irom3 , idiam3,                  &
                  ih1hlf, igmeva, igmhtf,                         &
! ---- Versions electriques
                  iefjou, ilapla , iqelec ,                       &
                  idjr  , idji   , idrad  ,                       &
! ---- Version aerorefrigerant
                  itemp4, ihumid

!--> Modele de flamme de premelange LWC

common / ilwcpp / itsc

!--> OPTIONS NUMERIQUES

! ---- Coefficient de relaxation de la masse volumique
!      RHO(n+1) = SRROM * RHO(n) + (1-SRROM) * RHO(n+1)
double precision srrom

common / roptcp / srrom


!--> GRANDEURS FOURNIES PAR L'UTILISATEUR EN CONDITIONS AUX LIMITES
!      PERMETTANT DE CALCULER AUTOMATIQUEMENT LA VITESSE, LA TURBULENCE,
!      L'ENTHALPIE D'ENTREE.
!    LES GRANDEURS CI-DESSOUS SONT COMMUNES A LA COMBUSTION GAZ ET AU
!      CHARBON.

!    POUR LES ENTREES UNIQUEMENT , IENT ETANT LE NUMERO DE ZONE FRONTIERE

!       IQIMP (IENT) --> Indicateur zone a debit impose
!       ICALKE(IENT) --> Indicateur type de condition sur la turbulence
!         0 : Utilisateur donne les valeurs
!         1 : Automatique a partir de DH
!                                         et de la vitesse d'entree
!         2 : Automatique a partir de l'intensite turbulente
!                                         et de la vitesse d'entree
!       XINTUR(IENT) --> Intensite turbulente (k=1.5(UREF*XINTUR)**2)
!       DH    (IENT) --> Diametre hydraulique
!       QCALC (IENT) --> Debit calcule  : raf la ; direc ds le sspgm

integer          iqimp(nozppm)  , icalke(nozppm)
double precision xintur(nozppm) , dh(nozppm)
!, QCALC(NOZPPM)

common / ippcli / iqimp         , icalke
common / rppcli / xintur        , dh
!, QCALC


! Pointeur dans IA sur IZFPPP pour reperage des zones frontieres associees
!   aux faces de bord
! Peut etre serait il plus approprie de le verser dans pointe


integer           iizfpp
common / ifropp / iizfpp

! NZFPPP Nombre de zones de bord (sur le proc courant)
! ILZPPP Liste des numeros de zone de bord (du proc courant)
! NOZAPM Numero de zone de bord atteint max
!   exemple zones 1 4 2 : NZFPPP=3,NOZAPM=4

integer           nozapm, nzfppp, ilzppp(nbzppm)
common / izonpp / nozapm, nzfppp, ilzppp

! FIN
