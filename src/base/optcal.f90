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

! Module for calculation options

module optcal

  !=============================================================================

  use, intrinsic :: iso_c_binding

  use paramx

  !=============================================================================

  ! Definition des equations
  !   istat
  !     = 1 prise en compte du terme instationnaire
  !     = 0 prise en compte du terme instationnaire
  !   iconv
  !     = 1 prise en compte de la convection
  !     = 0 non prise en compte de la convection
  !   idiff
  !     = 1 prise en compte de la diffusion (moleculaire et turbulente)
  !     = 0 non prise en compte de la diffusion (moleculaire et turbulente)
  !   idifft : si idiff = 1
  !     = 1 prise en compte de la diffusion turbulente
  !     = 0 non prise en compte de la diffusion turbulente
  !   idften
  !     = 1 scalar diffusivity
  !     = 3 orthotropic diffusivity
  !     = 6 symmetric tensor diffusivity
  !   iswdyn
  !     = 0 no dynamic relaxation
  !     = 1 dynamic relaxation depending on
  !       \$f \delta \varia^k \f$
  !     = 2 dynamic relaxation depending on
  !       \$f \delta \varia^k \f$  and
  !       \$f \delta \varia^{k-1} \f$

  integer, save :: istat(nvarmx), iconv(nvarmx), idiff(nvarmx), idifft(nvarmx),&
                   idften(nvarmx), iswdyn(nvarmx)

  ! Proprietes physiques rho et viscl constantes ou variables
  !    =1 variable, =0 constant
  !     sert lors des lectures de fichier suite pour eviter d'ecraser
  !     la valeur fournie par la valeur de l'ancien calcul.
  integer, save :: irovar, ivivar

  ! Algorithm to take into account the density variation in time

  !     idilat = 0 : boussinesq algorithm with constant density
  !              1 : dilatable steady algorithm (default)
  !              2 : dilatable unsteady algorithm
  !              3 : low-Mach algorithm
  !
  !           epsdp: parameter of diagonal pressure strengthening

  integer, save :: idilat

  double precision, save :: epsdp

  ! Schema en temps

  !  ischtp : indicateur de schema en temps
  !     = 2 : ordre 2
  !     = 1 : standard
  !  istmpf : indicateur de schema flux de masse
  !     = 2 theta schema avec theta > 0 (= 0.5 : ordre 2)
  !     = 0 theta schema avec theta = 0 (explicite)
  !     = 1 schema standard v1.0
  !  nterup : nombre d'iteration sur navier-stokes pour couplage vitesse/
  !           pression
  !  isno2t : indicateur d'extrapolation de termes sources Navier Stokes
  !           pour le schema en temps
  !  isto2t : indicateur d'extrapolation de termes sources des grandeurs
  !           turbulentes pour le schema en temps
  !  isso2t : indicateur d'extrapolation de termes sources des scalaires
  !           pour le theta schema en temps
  !  iroext : indicateur d'extrapolation de la masse volumique
  !           pour le schema en temps
  !  iviext : indicateur d'extrapolation de la viscosite totale
  !           pour le schema en temps
  !  ivsext : indicateur d'extrapolation de la diffusivite scalaire

  !  initvi : =1 si viscosite totale relue dans un suite

  !  initro : =1 si masse volumique relue dans un suite

  !  icpext : indicateur d'extrapolation de la masse volumique
  !           pour le schema en temps

  !  initcp : =1 si  chaleur specifique relue dans un suite
  !  initvs : =1 si  diffusivite scalaire relue dans un suite

  !  thetav : ponderation entre les pas de temps n et n+1 pour les
  !           variable principales
  !     = 1 : schema Euler implicite
  !     =1/2: schema centre en temps

  !  thetsn : schema en temps pour les termes sources de Navier Stokes
  !     = 0 : viscosite secondaire explicite
  !     =1/2: viscosite secondaire extrapolee en n+1/2
  !     = 1 : viscosite secondaire extrapolee en n+1
  !  thetst : schema en temps pour les termes sources des grandeurs turbulentes
  !     = 0 : viscosite secondaire explicite
  !     =1/2: viscosite secondaire extrapolee en n+1/2
  !     = 1 : viscosite secondaire extrapolee en n+1
  !  thetss : schema en temps pour les termes sources des scalaires
  !     = 0 : viscosite secondaire explicite
  !     =1/2: viscosite secondaire extrapolee en n+1/2
  !     = 1 : viscosite secondaire extrapolee en n+1
  !  thetfl : schema en temps pour le flux de masse
  !     = 0 : flux de masse explicite
  !     =1/2: flux de masse extrapole en n+1/2
  !     = 1 : flux de masse extrapole en n+1
  !  thetvi : schema en temps pour la viscosite totale
  !     = 0 : viscosite totale explicite
  !     =1/2: viscosite totale extrapolee en n+1/2
  !     = 1 : viscosite totale extrapolee en n+1
  !  thetro : schema en temps pour la masse volumique
  !     = 0 : masse volumique totale explicite
  !     =1/2: masse volumique totale extrapolee en n+1/2
  !     = 1 : masse volumique extrapolee en n+1
  !  thetcp : schema en temps pour la masse volumique
  !     = 0 : chaleur specifique totale explicite
  !     =1/2: chaleur specifique totale extrapolee en n+1/2
  !     = 1 : chaleur specifique extrapolee en n+1
  !  epsup  : tests de convergence du systeme vitesse/pression quand ce
  !           dernier est resolu par sous-iterations (point fixe)
  !  xnrmu  : norme de u(k+1) - u(k)
  !  xnrmu0 : norme de u(0)

  integer, save ::          nterup,                         &
                            ischtp, istmpf,                 &
                            isno2t, isto2t, isso2t(nscamx), &
                            iroext,                         &
                            iviext, icpext, ivsext(nscamx), &
                            initro, initvi,                 &
                            initcp, initvs(nscamx)
  double precision, save :: thetav(nvarmx), thetsn, thetst, &
                            thetss(nscamx),                 &
                            thetfl, thetro, thetvi,         &
                            thetcp, thetvs(nscamx), epsup , &
                            xnrmu0, xnrmu

  ! Schema convectif

  !  blencv : 100*(1-blencv) est le pourcentage d'upwind
  !     = 1 : pas d'upwind en dehors du test de pente
  !     = 0 : upwind
  !  ischcv : schema convectif centre ou second order
  !     = 1 : centre
  !     = 0 : second order
  !  isstpc : indicateur sans ou avec test de pente
  !     = 1 : sans test de pente
  !     = 0 : avec test de pente

  integer, save ::          ischcv(nvarmx), isstpc(nvarmx)
  double precision, save :: blencv(nvarmx)

  ! Reconstruction des gradients et des seconds membres
  !   imrgra : methode de recontruction des gradients
  !     = 0  : recontruction 97
  !     = 1  : moindres carres 99
  !     = 2  : moindres carres support etendu complet
  !     = 3  : moindres carres avec selection du support etendu
  !     = 4  : reconstruction 97 avec initialisation moindres carres
  !   anomax : angle de non orthogonalite des faces en radian au dela duquel
  !            on retient dans le support etendu des cellules voisines
  !            de la face les cellules dont un noeud est sur la face
  !   nswrgr : nombre de sweeps de reconstruction des gradients 97
  !   nswrsm : nombre de sweeps de reconstruction des seconds membres
  !   epsrgr : precision pour la   reconstruction des gradients 97
  !   epsrsm : precision pour la   reconstruction des seconds membres
  !   imligr : limitation des gradients
  !     < 0  : pas de limitation des gradients
  !     = 0  : premier ordre
  !     = 1  : second ordre
  !   climgr : facteur de limitation (>=1, =1 : forte limitation)
  !   ircflu : reconstruction des flux aux faces
  !     = 0  : non
  !     = 1  : oui
  !   extrag : extrapolation des gradients au bord (0 <= extrag <= 1)
  !     = 0  : non
  !     = 1  : oui

  integer, save ::          imrgra, nswrgr(nvarmx), nswrsm(nvarmx),   &
                            imligr(nvarmx)        , ircflu(nvarmx)

  double precision, save :: anomax ,                                  &
                            epsrgr(nvarmx), epsrsm(nvarmx),           &
                            climgr(nvarmx), extrag(nvarmx)

  ! Solveurs iteratifs
  !   nitmax : nombre d'iterations max
  !   epsilo : precision relative cherchee
  !   iresol
  !     =-1 : calcule automatiquement (0 si iconv=0, 1 sinon)
  !     = 0 : gradient conjugue
  !     = 1 : Jacobi
  !     = 2 : bi-CGSTAB
  !    et on ajoute ipol*1000 ou ipol est le degre du polynome de
  !       preconditionnement de Neumann
  !     en pratique, il semble que ce preconditonnement ne soit pas efficace
  !        on gagne 10% cpu sur un cas, on perd 3% sur un autre avec ipol=1
  !        on perd avec ipol=2
  !        ces valeurs ont ete obtenues sur de petits cas.
  !   idircl : decalage de la diagonale de la matrice s'il n'y a pas de Dirichlet
  !     = 0 : non
  !     = 1 : oui
  !     le code calcule automatiquement pour chaque variable ndircl, nombre de
  !        CL de Dirichlet, et en deduit s'il doit decaler ou pas la diagonale

  integer, save ::          nitmax(nvarmx),iresol(nvarmx),idircl(nvarmx),   &
                            ndircl(nvarmx)

  double precision, save :: epsilo(nvarmx)

  ! Multigrille
  !   imgr
  !     = 0 pas de multigrille
  !     = 1        multigrille algebrique
  !   ncymax : nombre max de cycles
  !   nitmgf : nombre d'iter sur maillage fin
  !   rlxp1  :

  integer, save ::          imgr(nvarmx), ncymax(nvarmx), nitmgf(nvarmx)
  double precision, save :: rlxp1

  ! Gestion du calcul
  !   isuite : suite de calcul
  !     = 0 pour sfs
  !     = 1 pour suite de calcul
  !   iscold : correspondance nouveaux-anciens scalaires
  !   iecaux : ecriture du suite auxiliaire
  !   ileaux : lecture  du suite auxiliaire
  !   isuit1 : suite du module thermique 1D en paroi
  !   isuict : suite du module aerorefrigerant
  !   isuivo : suite de la methode des vortex
  !   isuisy : suite des methodes d entree LES

  integer, save :: isuite , ileaux, iecaux, iscold(nscamx),        &
                   isuit1 , isuict, isuivo, isuisy

  ! Time step

  !> Absolute time step number for previous calculation.
  integer(c_int), pointer, save :: ntpabs

  !> Current absolute time step number.
  !> In case of restart, this is equal to ntpabs + number of new iterations.
  integer(c_int), pointer, save :: ntcabs

  !> Maximum absolute time step number.
  integer(c_int), pointer, save :: ntmabs

  !> Absolute time value for previous calculation.
  real(c_double), pointer, save :: ttpabs

  !> Current absolute time.
  !> In case of restart, this is equal to ttpabs + additional computed time.
  real(c_double), pointer, save :: ttcabs

  !> Maximum absolute time.
  real(c_double), pointer, save :: ttmabs

  ! Option pas de temps
  !   inpdt0 : indicateur "zero pas de temps"
  !     = 0 : calcul standard
  !     = 1 : pour ne faire aucun pas de temps (0 sinon)
  !             pour les calculs non suite :
  !               on saute uniquement les resolutions (Navier-Stokes,
  !               turbulence, scalaires...)
  !             pour les calculs suite :
  !               on saute les resolutions (navier-stokes,
  !               turbulence, scalaires...) et le calcul des proprietes
  !               physiques, les conditions aux limites (les grandeurs
  !               sont lues dans le fichier suite)
  !   idtvar : pas de temps variable
  !     = -1 : algorithme stationnaire
  !     =  0 : pas de temps constant
  !     =  1 : pas de temps uniforme en espace et variable en temps
  !     =  2 : pas de temps variable en espace et variable en temps
  !   iptlro : limitation du pas de temps liee aux effets de densite
  !     = 0 : non
  !     = 1 : oui
  !   coumax : nombre de Courant         maximum        (idtvar non nul)
  !   foumax : nombre de         Fourier maximum        (idtvar non nul)
  !   varrdt : variation relative permise de dt         (idtvar non nul)
  !   dtmin, dtmax : valeur limite min et max de dt     (idtvar non nul)
  !       prendre pour dtmax = max (ld/ud, sqrt(lt/(gdelta rho/rho)), ...)
  !   cdtvar : coef multiplicatif pour le pas de temps de chaque variable
  !         pour u,v,w,p il est inutilise
  !         pour k,e    on prend la meme valeur : celle de k
  !         pour Rij, e on prend la meme valeur : celle de r11
  !   relaxv : relaxation des variables (1 pas de relax)
  !   relxst : coefficient de relaxation de base stationnaire

  integer, save ::          inpdt0, idtvar,iptlro
  double precision, save :: dtref,coumax,foumax,                  &
                            dtmin,dtmax ,varrdt,cdtvar(nvarmx),   &
                            relaxv(nvarmx), relxst

  ! turbulence
  !  iturb
  !    = 0  pas de turbulence
  !    = 10 longueur de melange
  !    = 20, 21 k-epsilon
  !         * 20 modele standard
  !         * 21 modele a production lineaire
  !    = 30, 31 Rij-epsilon
  !         * 30 modele standard (LRR)
  !         * 31 modele ssg
  !    = 40, 41, 42 les
  !         * 40 modele de Smagorinsky constant
  !         * 41 modele de Smagorinsky dynamique "classique"
  !         * 42 modele de Smagorinsky dynamique de "Piomelli et Liu"
  !    = 50 v2f phi-model
  !    = 60 k-omega sst
  !    = 70 Spalart-Allmaras
  !  itytur
  !    = int(iturb/10) pour distinguer rapidement les classes de modeles
  !  irccor
  !    = 0 no rotation/curvature correction for an eddy viscosity turbulence
  !        model
  !    = 1 activate rotation/curvature correction for an eddy viscosity
  !        turbulence model
  !  itycor
  !    = 1 Cazalbou correction: default when irccor = 1 and itytur = 2 or 5
  !    = 2 Spalart-Shur correction: default when irccor = 1 and
  !        iturb = 60 or 70
  !  ideuch
  !    = 0 une echelle       (deux echelles = faux)
  !    = 1 deux echelles     (deux echelles = vrai)
  !    = 2 deux echelles limitation de yplus a ypluli (scalable wall function)
  !  ilogpo
  !    = 0 une echelle  avec loi en puissance
  !    = 1 une echelles avec loi log
  !  iclkep
  !    = 0 clipping en valeur absolue de k et epsilon
  !    = 1 clipping couple k-epsilon base sur des relations physiques
  !  igrhok
  !    = 1     prise en compte de 2/3 rho grad k dans navier stokes
  !    = 0 non prise en compte de 2/3 rho grad k dans navier stokes
  !  igrake
  !    = 1 gravite dans k-epsilon
  !    = 0 sinon
  !  igrari
  !    = 1 gravite dans Rij-epsilon
  !    = 0 sinon
  !  iscalt numero du scalaire qui tient lieu de temperature
  !    donc variable isca(iscalt)
  !  ikecou
  !    = 1 k-epsilon couple en increments
  !    = 0 sinon
  !  irijnu
  !         = 1 viscosite dans la matrice en increments de vitesse (Rij)
  !         = 0 sinon
  !  irijrb
  !         = 1 traitement precis de Rij au bord, voir condli      (Rij)
  !         = 0 sinon
  !  idifre
  !         = 1 traitement complet de la diagonale du tenseur de
  !             diffusion de Rij et epsilon (Rij)
  !         = 0 traitement simplifie
  !  iclsyr
  !         = 1 implicitation partielle de Rij dans les cl de symetrie
  !         = 0 pas d'implicitation
  !  iclptr
  !         = 1 implicitation partielle de Rij et epsilon dans les cl
  !             de paroi turbulente
  !         = 0 pas d'implicitation
  !  idries : amortissement de type Van Driest a la paroi
  !         = 0 sans amortissement
  !         = 1 avec amortissement
  !  ivrtex : utilisation de la methode des vortex
  !         = 0 sans methode des vortex
  !         = 1 avec methode des vortex

  integer, save :: iturb , itytur,                 &
                   irccor, itycor,                 &
                   ideuch, ilogpo, iclkep,         &
                   igrhok, igrake,                 &
                   iscalt, ikecou,                 &
                   irijnu, irijrb, irijec,         &
                   igrari, idifre, iclsyr,         &
                   iclptr, idries, ivrtex

  ! ivisse prise en compte de -2/3 grad(mu div(u)) + div(mu (grad_t(u)))

  integer, save :: ivisse

  ! Stokes
  !   irevmc
  !     = 2 pour reconstruction des vitesses de type rt0
  !     = 1 pour reconstruction des vitesses avec gradient de l'increment
  !           de pression par moindres carres
  !     = 0 sinon
  !   iprco
  !     = 0 pour calcul sans pression continuite
  !     = 1 pour calcul avec pression continuite
  !   arak proportion d'Arakawa (1 pour Arakawa complet)
  !   relaxv relaxation des variables (1 pas de relax)
  !   rnormp normalisation pour la convergence de resolp

  integer, save ::          irevmc, iprco , irnpnw
  double precision, save :: rnormp, arak

  !   ivelco
  !     = 1 resolution couplee des composantes de vitesse
  !     = 0 resolution decouplee des composantes de vitesse (Standard)

  integer, save :: ivelco

  !   iporos
  !     = 1 Taking porosity into account
  !     = 0 Standard algorithm (Without porosity)

  integer, save :: iporos

  ! ipucou algorithme couplage instationnaire vitesse/pression

  integer, save :: ipucou

  ! iccvfg calcul a champ de vitesse fige

  integer, save :: iccvfg

  ! Calcul de la viscosite

  integer, save :: imvisf

  ! modele de flux turbulent u'T' pour un scalaire T
  !  iturbt
  !    = 0  SGDH
  !    = 10 GGDH
  !    = 20 AFM
  !    = 30 DFM (Transport equation modelized)
  !    = 11 EB-GGDH
  !    = 21 EB-AFM
  !    = 31 EB-DFM
  !  ityturt
  !    = int(iturbt/10) pour distinguer rapidement les classes de modeles

  integer, save :: iturbt , ityturt

  ! Type des conditions limites et index min et max
  !                 des sous listes defaces de bord

  integer, save :: idebty(ntypmx), ifinty(ntypmx)

  !  itrbrb = 1 traitement precis de la temperature au bord, voir condli
  !             (utilise pour couplage syrthes)
  !         = 0 sinon
  !  icpsyr = 1 si scalaire couple a syrthes
  !    donc pour le moment vaut 1 pour iscalt uniquement

  integer, save :: itbrrb, icpsyr(nscamx)

  !   Prise en compte de l'equilibre entre le gradient de pression
  !        et les termes sources de gravite et de perte de charge

  !     iphydr = 0 algorithme sans prise en compte de l'equilibre
  !            = 1 algorithme avec prise en compte de l'equilibre
  !     icalhy = 0 pas de calcul de la pression hydrostatique pour les
  !                dirichlets de pression en sortie
  !            = 1        calcul de la pression hydrostatique pour les
  !                Dirichlets de pression en sortie

  integer, save :: iphydr, icalhy

  !   Calcul des estimateurs

  integer, save :: iescal(nestmx)

  !   Calcul des moyennes temporelles (calcul des moments)

  !  nbmomt : nombre de moyennes demandees
  !  nbdtcm : nombre de tableaux ncel pour le temps cumule
  !  ntdmom : numero du pas de temps initial pour le calcul du moment
  !  ttdmom : temps initial pour le calcul du moment
  !  imoold : numero de l'ancien moment correspondant en cas de suite
  !  icmome : pointeur pour les moments (donne un numero de propriete)
  !           s'utilise ainsi propce(iel,ipproc(icmome(imom)))
  !  idtmom : numero du temps cumule associe aux moments
  !           ce numero va de 1 a n pour les temps cumules non uniformes
  !                     et de -1 a -p pour les temps cumules uniformes
  !           s'utilise ainsi :
  !              si idtmom(imom) > 0 propce(iel,ipropc(icdtmo(idtmom(imom))))
  !              si idtmom(imom) < 0 dtcmom(-idtmom(imom))
  !  idfmom : numero des variables composant le moment idfmom(jj,imom)
  !  idgmom : degre du moment
  !  icdtmo : numero de propriete du temps cumule (voir idtmom)
  !  ippmom : repere pour le post si on doit diviser la variable
  !           par un temps cumule (voir memtri et usvpst)
  !  dtcmom : valeur du pas de temps cumule quand il est uniforme (voir idtmom).

  integer, save ::          nbmomt, nbdtcm,                                 &
                            ntdmom(nbmomx), imoold(nbmomx),                 &
                            icmome(nbmomx), idtmom(nbmomx),                 &
                            idfmom(ndgmox,nbmomx),          idgmom(nbmomx), &
                            icdtmo(nbmomx), ippmom(nvppmx)
  double precision, save :: dtcmom(nbmomx), ttdmom(nbmomx)

  ! Indicateur pertes de charge global (ie somme sur les processeurs
  !   de ncepdc)

  integer, save :: ncpdct

  ! Indicateur module thermique 1d global (ie somme sur les processeurs
  !   de nfpt1d)

  integer, save :: nfpt1t

  ! Indicateur termes sources de masse global (ie somme sur les processeurs
  !   de ncetsm)

  integer, save :: nctsmt

  ! Indicateur de passage dans l'initialisation des
  !                       variables par l'utilisateur
  !          iusini = 1 passage dans usiniv ou ppiniv
  !                   0 pas de passage (ni iusini ni ppiniv)
  !          iuscfp = 1 passage dans uscfpv
  !                   0 pas de passage

  integer, save :: iusini, iuscfp

  ! Parametres numeriques pour le calcul de la distance a la paroi

  ! ineedy : = 1 distance a la paroi est necessaire pour le calcul
  !          = 0 distance a la paroi n'est pas necessaire
  ! imajdy : = 1 distance a la paroi a ete mise a jour
  !          = 0 distance a la paroi n'a pas ete mise a jour
  ! icdpar : = 1 calcul standard (et relecture en suite de calcul)
  !          = 2 calcul ancien   (et relecture en suite de calcul)
  !          =-1 forcer le recalcul en suite (par calcul standard)
  !          =-2 forcer le recalcul en suite (par calcul ancien)
  ! nitmay : nombre max d'iterations pour les resolutions iteratives
  ! nswrsy : nombre de sweep pour reconstruction des s.m.
  ! nswrgy : nombre de sweep pour reconstruction des gradients
  ! imligy : methode de limitation du gradient
  ! ircfly : indicateur pour reconstruction des flux
  ! ischcy : indicateur du schema en espace
  ! isstpy : indicateur pour test de pente
  ! imgrpy : multigrille
  ! iwarny : niveau d'impression
  ! ntcmxy : nombre max d'iteration pour la convection de y

  integer, save :: ineedy, imajdy, icdpar,    &
                   nitmay, nswrsy, nswrgy,    &
                   imligy, ircfly, ischcy,    &
                   isstpy, imgrpy, iwarny,    &
                   ntcmxy

  ! blency : 1 - proportion d'upwind
  ! epsily : precision pour resolution iterative
  ! epsrsy : precision pour la reconstruction du second membre
  ! epsrgy : precision pour la reconstruction des gradients
  ! climgy : coef gradient*distance/ecart
  ! extray : coef d'extrapolation des gradients
  ! coumxy : valeur max   du courant pour equation convection
  ! epscvy : precision pour convergence equation convection stationnaire
  ! yplmxy : valeur max   de yplus au dessus de laquelle l'amortissement de
  !          Van Driest est sans effet et donc pour laquelle un calcul de
  !          yplus moins precis est suffisant

  double precision, save :: blency, epsily, epsrsy,    &
                            epsrgy, climgy, extray,    &
                            coumxy, epscvy, yplmxy

  ! Parametres numeriques pour le calcul des efforts aux bords

  ! ineedf : = 1 on calcule les efforts aux parois
  !          = 0 on ne calcule pas les efforts aux parois

  integer, save :: ineedf

  !=============================================================================

  interface

    !---------------------------------------------------------------------------

    !> \cond DOXYGEN_SHOULD_SKIP_THIS

    !---------------------------------------------------------------------------

    ! Interface to C function retrieving pointers to members of the
    ! global time step structure

    subroutine cs_f_time_step_get_pointers(nt_prev, nt_cur, nt_max,  &
                                           t_prev, t_cur, t_max)     &
      bind(C, name='cs_f_time_step_get_pointers')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: nt_prev, nt_cur, nt_max
      type(c_ptr), intent(out) :: t_prev, t_cur, t_max
    end subroutine cs_f_time_step_get_pointers

    !---------------------------------------------------------------------------

    !> \endcond DOXYGEN_SHOULD_SKIP_THIS

    !---------------------------------------------------------------------------

  end interface

  !=============================================================================

contains

  !=============================================================================

  !> \brief Initialize Fortran time step API.
  !> This maps Fortran pointers to global C structure members.

  subroutine time_step_init

    use, intrinsic :: iso_c_binding
    implicit none

    ! Local variables

    type(c_ptr) :: c_ntpabs, c_ntcabs, c_ntmabs
    type(c_ptr) :: c_ttpabs, c_ttcabs, c_ttmabs

    call cs_f_time_step_get_pointers(c_ntpabs, c_ntcabs, c_ntmabs, &
                                     c_ttpabs, c_ttcabs, c_ttmabs)

    call c_f_pointer(c_ntpabs, ntpabs)
    call c_f_pointer(c_ntcabs, ntcabs)
    call c_f_pointer(c_ntmabs, ntmabs)

    call c_f_pointer(c_ttpabs, ttpabs)
    call c_f_pointer(c_ttcabs, ttcabs)
    call c_f_pointer(c_ttmabs, ttmabs)

  end subroutine time_step_init

  !=============================================================================

end module optcal
