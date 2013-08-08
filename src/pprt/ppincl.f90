!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2013 EDF S.A.
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

!> \file ppincl.f90
!> General module for specific physics

module ppincl

  !===========================================================================

  use ppppar
  use ppthch

  !=============================================================================

  implicit none

  !===========================================================================

  !--> TABLEAU INDICATEURS DU CHOIX DE LA PHYSIQUE PARTICULIERE CHOISIE

  integer   nmodmx
  parameter(nmodmx = 50)

  integer, save ::           ippmod(nmodmx)

  ! ---- Indicateur global de physique particuliere
  !        IPPMOD(IPHPAR) = 0 : pas de physique particuliere
  !                         1 : physique particuliere enclenchee
  !                         2 : physique particuliere avec pilotage du
  !                             rayonnement par fichier parametrique
  integer :: iphpar

  ! ---- Modeles propres a la combustion gaz ICO...
  integer ::  icod3p, icodeq, icoebu, icobml, icolwc, isoot

  ! ---- Modeles propres aux versions effet Joule et conduction ionique
  integer ::  ieljou, ielarc, ielion

  ! ---- Modeles propres a la combustion charbon pulverise couplee Lagrangien
  integer ::  icpl3c
  integer ::  iccoal
  ! ---- Coal with drift (0: without drift (default), 1: with)
  integer ::  i_coal_drift

  ! ---- Modeles propres a la combustion fuel
  integer ::  icfuel

  ! ---- Modele compressible
  integer ::  icompf

  ! ---- Modele atmospherique
  integer ::  iatmos

  ! ---- Modele aerorefrigerants
  integer ::  iaeros

  parameter       (iphpar = 1 , icod3p = 2 , icodeq = 3 ,           &
                   icoebu = 4 , icobml = 5 , icolwc = 6 ,           &
                   icpl3c = 7 , icfuel = 8 , ieljou = 9 ,           &
                   ielarc = 10, ielion = 11, icompf = 12,           &
                   iatmos = 13, iaeros = 14, iccoal = 15)


  !--> NOMBRE DE VARIABLES ALGEBRIQUES OU D'ETAT
  !    pour la physique particuliere NSALPP
  !    total NSALTO
  integer, save :: nsalpp, nsalto

  !--> POINTEURS VARIABLES COMBUSTION GAZ

  ! ---- Variables transportees
  integer, save :: ifm, ifp2m, iygfm, icm, icp2m, ifpcpm
  integer, save :: iyfm, iyfp2m, icoyfp

  ! ---- Variables d'etat
  integer, save :: iym(ngazgm), itemp, ifmin, ifmax
  integer, save :: ickabs, it4m, it3m

  ! --- Pointeurs proprietes (PROPCE)
  integer, save :: itsc

  ! --- Pointers for soot model
  !     INPM  : pointer for soot precursor number in isca (isoot = 1)
  !     IFSM  : pointer for soot mass fraction in isca (isoot = 1)
  integer, save :: inpm, ifsm

  !--> POINTEURS VARIABLES COMBUSTION CHARBON PULVERISE

  ! ---- Variables transportees
  !        Phase continue (melange gazeux)
  integer, save :: if1m(ncharm), if2m(ncharm)
  integer, save :: if3m, if4m, if5m, if6m, if7m , if8m, if9m
  integer, save :: if4p2m , ifvp2m , if3mc2
  !        Phase dispersee (classe de particules)
  integer, save :: ixck(nclcpm), ixch(nclcpm), inp(nclcpm)
  integer, save :: ih2(nclcpm) , ixwt(nclcpm)
  ! Pointers to x2*age(particles) and
  ! x1*age(gas phase)
  integer, save :: iagecp_temp(nclcpm), iaggas_temp

  ! ---- Variables d'etat
  !        Phase continue (melange gazeux)
  integer, save :: iym1(ngazem), itemp1 , irom1 , immel
  !        Phase dispersee (classes de particules)
  integer, save :: itemp2(nclcpm), irom2(nclcpm), idiam2(nclcpm)
  integer, save :: ix2(nclcpm)
  integer, save :: igmdch(nclcpm), igmhet(nclcpm) , ighco2(nclcpm)
  integer, save :: igmdv1(nclcpm), igmdv2(nclcpm)
  integer, save :: igmsec(nclcpm)
  !        Bilan
  integer, save :: ibcarbone,iboxygen,ibhydrogen

  !--> POINTEURS VARIABLES COMBUSTION FUEL

  ! ---- Variables transportees
  !        Phase continue
  integer, save :: ifvap
  !        Phase dispersee
  integer, save :: ihlf(nclcpm)
  integer, save :: ixkf(nclcpm), ixfol(nclcpm), ing(nclcpm)

  ! ---- Variables d'etat
  !        Phase continue
  integer, save :: iyfol(nclcpm)
  !        Phase dispersee
  integer, save :: ih1hlf(nclcpm), igmhtf(nclcpm), igmeva(nclcpm)

  !--> POINTEURS VARIABLES VERSION ELECTRIQUES

  !        Dimension des 'vecteurs electriques'      = NDIMVE
  integer    ndimve
  parameter (ndimve = 3)

  ! ---- Variables transportees
  !        Potentiel reel       = IPOTR
  !        Potentiel imaginaire = IPOTI
  !        Composantes du potentiel vecteur magnetique = IPOTVA()
  !        Fraction massique des constituants = IYCOEL()

  integer, save :: ipotr, ipoti, ipotva(ndimve), iycoel(ngazgm)

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

  integer, save :: iefjou , iqelec , ilapla(ndimve)
  integer, save :: idjr(ndimve)    , idji(ndimve)   , idrad

  !--> POINTEURS VARIABLES CONDUCTION IONIQUE

  integer   nesiom
  parameter (nesiom = 10)
  integer   nespio

  ! ---- Variables transportees
  !        par espece
  integer, save :: iymion(nesiom)

  ! ---- Variables d'etat

  !--> POINTEURS COMPRESSIBLE

  ! ---- Variables transportees par phase
  integer, save :: irho, ienerg, itempk
  ! ---- Proprietes supplementaires par phase
  integer, save :: icv, iviscv, ieos

  !     COMMON complete plus bas

  ! ---- Aliases pour les conditions aux limites
  integer, save :: irun, irunh

  ! ---- Proprietes supplementaires par phase
  double precision, save :: cv0, viscv0

  ! ---- Prediction de pression par une equation d'evolution
  integer, save :: ippred
  ! ---- Flux de masse specifique pour la vitesse
  integer, save :: iflmau
  ! ---- Utilisation de la pression predite pour resoudre Navier-Stokes
  integer, save :: igrdpp
  ! --- Conditions aux limites prenant en compte l'equilibre hydrostatique
  integer, save :: icfgrp

  ! ---- Flux de bord convectifs QDM et energie (numero de PROPFB)
  integer, save ::           ifbrhu , ifbrhv ,               &
                             ifbrhw , ifbene

  !--> POINTEURS AEROREFRIGERANTS

  ! ---- Variables transportees
  integer, save :: itemp4, ihumid

  !--> POINTEUR RELATIF A LA VARIABLE ENTHALPIE

  integer, save :: ihm

  !--> REMPLISSAGE COMMON POINTEURS VARIABLES TRANSPORTEES
  !                                 VARIABLES D'ETAT

  !--> OPTIONS NUMERIQUES

  ! ---- Coefficient de relaxation de la masse volumique
  !      RHO(n+1) = SRROM * RHO(n) + (1-SRROM) * RHO(n+1)
  double precision, save :: srrom

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

  integer, save ::          iqimp(nozppm)  , icalke(nozppm)
  double precision, save :: xintur(nozppm) , dh(nozppm)

  ! NZFPPP Nombre de zones de bord (sur le proc courant)
  ! ILZPPP Liste des numeros de zone de bord (du proc courant)
  ! NOZAPM Numero de zone de bord atteint max
  !   exemple zones 1 4 2 : NZFPPP=3,NOZAPM=4

  integer, save ::           nozapm, nzfppp, ilzppp(nbzppm)

  !===========================================================================

end module ppincl
