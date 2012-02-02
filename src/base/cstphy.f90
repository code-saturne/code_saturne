!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2012 EDF S.A.
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

! Module for physical constants

module cstphy

  !=============================================================================

  use paramx

  !=============================================================================

  ! Kelvin

  ! tkelvi       --> =  273,15
  ! tkelvn       --> = -273,15

  double precision :: tkelvi, tkelvn
  parameter(tkelvi = 273.15d0, tkelvn = -273.15d0)

  ! Calories (1 cal = xcal2j J)

  double precision :: xcal2j
  parameter(xcal2j = 4.1855d0)

  ! Stephan Boltzmann

  double precision :: stephn
  parameter(stephn = 5.6703d-8)

  ! Perfect gas constant for air (mixture)

  double precision :: rair
  parameter(rair = 287.d0)

  ! Gravity

  double precision, save :: gx, gy, gz

  ! Rotation vector

  integer, save :: icorio

  double precision, save :: omegax, omegay, omegaz
  double precision, save :: irot(3,3), prot(3,3), qrot(3,3), rrot(3,3)

  ! Constantes physiques du fluide
  !   ixyzp0 : indicateur de remplissage de xyzp0
  !   ro0    : masse volumique    de reference
  !   viscl0 : viscosite          de reference
  !   p0     : pression totale    de reference
  !   pred0  : pression reduite   de reference
  !   xyzp0  : position pression  de reference
  !   t0     : temperature        de reference
  !   cp0    : chaleur specifique de reference

  integer, save ::          ixyzp0
  double precision, save :: ro0    , viscl0,     &
                            p0     , pred0 ,     &
                            xyzp0(3), t0    ,     &
                            cp0

  ! Turbulence
  !   ivisls = 0 : viscosite laminaire constante = visls0
  !   xkappa : cst de Karman (~0.42)
  !   cstlog : cst de la loi log: 1/xkappa*log(yplus) + cstlog (~5.2)
  !   ypluli : yplus limite 1./xkappa ou 10.88 si ideuch=2
  !   *pow   : coeff Werner and Wengle
  !   cmu025 = cmu**0.25
  !   ce1, ce2, sigmak, sigmae :
  !            constantes du k-epsilon
  !   c*rij* : constantes du Rij-epsilon standard (LRR)
  !   cssg*  : constantes specifiques du Rij-epsilon SSG
  !   cebm*  : constants of the Rij-epsilon EBRSM
  !   csebm  : constant of the Rij-epsilon EBRSM
  !   cebme2 : constant of the Rij-epsilon EBRSM
  !   cebmmu : constant of the Rij-epsilon EBRSM
  !   xcl    : constant of the Rij-epsilon EBRSM
  !   sigebm : constants sigmae for the Rij-epsilon EBRSM
  !   xa1    : constant in the expression of Ce1' for the Rij-epsilon EBRSM
  !   xct    : constant of the Rij-epsilon EBRSM
  !   xceta  : constant of the Rij-epsilon EBRSM
  !   sigebm : constant sigmae for the Rij-epsilon EBRSM
  !   cv2f*  : constantes specifiques du v2f "phi-model" (f-barre)
  !   cpal*  : constantes specifiques du v2f "BL-v2k" (ou phi-alpha)
  !   ckw*   : constantes specifiques du k-omega SST
  !            (sk=sigma_k, sw=sigma_w, bt=beta, gm=gamma)
  !   csa*   : constantes specifiques de Spalart-Allmaras
  !   almax  : echelle de longueur turbulente
  !   uref   : vitesse de reference
  !   xlomlg : longueur pour longueur de melange
  !   xlesfl, ales, bles
  !       delta = xlesfl * (ales*volume)^bles (largeur du filtre utilise
  !       en fonction du volume de la cellule)
  !   csmago
  !       la constante de Smagorinsky theorique vaut 0.18
  !       pour un canal plan, on prendra cependant plutot 0.065
  !   xlesfd
  !       Dans le cas d un modele dynamique, xlesfd est le rapport entre la
  !       largeur du filtre explicite et celle du filtre implicite
  !   smagmx
  !       Constante de Smagorinsky maximale souhaitee (on peut prendre 10*csmago)
  !   idries
  !       Amortissement Van Driest active (=1) ou non (=0)
  !   CDRIES
  !       Constante de Van Driest dans (1-exp(-y+/cdries))
  !   ce4    : Coefficient du terme interfacial dans k-eps
  !            (Ce coefficient sert en Lagrangien)
  !   volmin : volume de controle minimal
  !   volmax : volume de controle maximal
  !   voltot : volume total du domaine

  double precision, save :: xkappa , cstlog , ypluli  ,                     &
                            apow   , bpow   , cpow   , dpow   ,             &
                            cmu    , cmu025 , ce1    , ce2    , ce4    ,    &
                            sigmak , sigmae ,                               &
                            crij1  , crij2  , crij3  , crijep , csrij  ,    &
                            crijp1 , crijp2 ,                               &
                            cssge2 , cssgs1 , cssgs2 ,                      &
                            cssgr1 , cssgr2 , cssgr3 , cssgr4 , cssgr5 ,    &
                            cebms1 , cebms2 , cebmr1 ,                      &
                            cebmr2 , cebmr3 , cebmr4 , cebmr5 , cebmr6 ,    &
                            csebm  , cebme2 , cebmmu , xcl    , sigebm ,    &
                            xa1    , xct    , xceta  ,                      &
                            cv2fa1 , cv2fe2 , cv2fmu , cv2fc1 , cv2fc2 ,    &
                            cv2fct , cv2fcl , cv2fet ,                      &
                            cpale1 , cpale2 , cpale3 , cpale4 , cpalse ,    &
                            cpalmu , cpalc1 , cpalc2 ,                      &
                            cpalct , cpalcl , cpalet ,                      &
                            ckwsk1 , ckwsk2 , ckwsw1 , ckwsw2 , ckwbt1 ,    &
                            ckwbt2 , ckwgm1 , ckwgm2 , ckwa1  , ckwc1  ,    &
                            csab1  , csab2  , csasig , csav1  , csaw1  ,    &
                            csaw2  , csaw3  ,                               &
                            volmin , volmax , voltot ,                      &
                            almax   , uref  ,                               &
                            xlomlg  ,                                       &
                            xlesfl  , ales  , bles,                         &
                            csmago  , cdries,                               &
                            xlesfd  , smagmx,                               &
                            cwale


  ! Constantes pour les scalaires

  ! iscsth :
  !   -1 : de type temperature en C (      Cp pour la loi de paroi)
  !    0 : scalaire passif      (ie pas de Cp pour la loi de paroi)
  !    1 : de type temperature en K (      Cp pour la loi de paroi)
  !    2 : enthalpie            (ie pas de Cp pour la loi de paroi)
  !      la distinction C/K sert en rayonnement
  ! ivisls : si positif strictement, indique que la viscosite associee
  !            au scalaire est variable, et la valeur est le numero
  !            d'ordre de la viscosite dans le tableau des viscosites
  !            variables
  ! ivissa : comme ivisls sauf que sert au stockage de la viscosite au
  !          pas de temps precedent
  ! iclvfl : 0 : clipping des variances a zero
  !          1 : clipping des variances a zero et a f(1-f)
  !          2 : clipping des variances a max(zero,scamin) et scamax
  ! iscavr : numero du scalaire associe a la variance ou zero
  !          si le scalaire n'est pas une variance
  ! scamin, scamax : min et max pour clipping des scalaires
  !                  on ne clippe que si scamin < scamax
  ! visls0 : viscosite des scalaires si constante
  ! sigmas : prandtl des scalaires
  ! rvarfl : coeff de dissipation des variances

  integer, save ::          iscsth(nscamx), ivisls(nscamx), ivissa(nscamx),  &
                            iclvfl(nscamx), iscavr(nscamx)
  double precision, save :: scamin(nscamx), scamax(nscamx),                  &
                            visls0(nscamx),sigmas(nscamx),                   &
                            rvarfl(nscamx)

  !=============================================================================

end module cstphy
