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

subroutine clptur &
!================

 ( nvar   , nscal  ,                                              &
   isvhb  ,                                                       &
   icodcl ,                                                       &
   dt     , rtp    , rtpa   , propce , propfa , propfb , rcodcl , &
   coefu  , rijipb , coefa  , coefb  , visvdr ,                   &
   hbord  , thbord )

!===============================================================================
! FONCTION :
! --------

! CONDITIONS LIMITES EN PAROI TURBULENTE POUR TOUTES LES VARIABLES

! ON SUPPOSE QUE ICODCL(IU) = 5 =>
!                     PAROI POUR TOUTES LES VARIABLES TURBULENTES
!  (A PRIORI PEU RESTRICTIF EN MONOPHASIQUE)

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! isvhb            ! e  ! <-- ! indicateur de sauvegarde des                   !
!                  !    !     !  coefficients d'echange aux bords              !
! icodcl           ! te ! --> ! code de condition limites aux faces            !
!  (nfabor,nvar    !    !     !  de bord                                       !
!                  !    !     ! = 1   -> dirichlet                             !
!                  !    !     ! = 3   -> densite de flux                       !
!                  !    !     ! = 4   -> glissemt et u.n=0 (vitesse)           !
!                  !    !     ! = 5   -> frottemt et u.n=0 (vitesse)           !
!                  !    !     ! = 6   -> rugosite et u.n=0 (vitesse)           !
!                  !    !     ! = 9   -> entree/sortie libre (vitesse          !
!                  !    !     !  entrante eventuelle     bloquee               !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! rcodcl           ! tr ! --> ! valeur des conditions aux limites              !
!  (nfabor,nvar    !    !     !  aux faces de bord                             !
!                  !    !     ! rcodcl(1) = valeur du dirichlet                !
!                  !    !     ! rcodcl(2) = valeur du coef. d'echange          !
!                  !    !     !  ext. (infinie si pas d'echange)               !
!                  !    !     ! rcodcl(3) = valeur de la densite de            !
!                  !    !     !  flux (negatif si gain) w/m2 ou                !
!                  !    !     !  hauteur de rugosite (m) si icodcl=6           !
!                  !    !     ! pour les vitesses (vistl+visct)*gradu          !
!                  !    !     ! pour la pression             dt*gradp          !
!                  !    !     ! pour les scalaires                             !
!                  !    !     !        cp*(viscls+visct/sigmas)*gradt          !
! coefu            ! tr ! <-- ! tab de trav pour valeurs en iprime             !
! (nfabor,3   )    !    !     !  des comp de la vitesse au bord                !
! rijipb           ! tr ! <-- ! tab de trav pour valeurs en iprime             !
! (nfabor,6   )    !    !     !  des rij au bord                               !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
! visvdr(ncelet)   ! tr ! <-- ! viscosite dynamique ds les cellules            !
!                  !    !     !  de bord apres amortisst de v driest           !
! hbord            ! tr ! --> ! coefficients d'echange aux bords               !
! (nfabor)         !    !     !                                                !
! thbord           ! tr ! <-- ! temperature aux bords en i'                    !
! (nfabor)         !    !     !    (plus exactmt : var. energetique)           !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use optcal
use cstphy
use cstnum
use pointe
use entsor
use albase
use parall
use ppppar
use ppthch
use ppincl
use radiat
use cplsat
use mesh
use lagran

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal
integer          isvhb

integer          icodcl(nfabor,nvar)

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision rcodcl(nfabor,nvar,3)
double precision coefu(nfabor,ndim), rijipb(nfabor,6)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision visvdr(ncelet)
double precision hbord(nfabor),thbord(nfabor)

! Local variables

integer          ifac, iel, ivar, isou, ii, jj, kk, ll, isvhbl
integer          ihcp, iscal
integer          imprim, modntl
integer          inturb, inlami, iuiptn
integer          iclu  , iclv  , iclw  , iclk  , iclep
integer          iclnu
integer          icl11 , icl22 , icl33 , icl12 , icl13 , icl23
integer          icluf , iclvf , iclwf , iclphi, iclfb , iclal , iclomg
integer          ipcrom, ipcvis, ipcvst, ipccp , ipccv
integer          iclvar, ipcvsl, iclvaf
integer          iclalp
double precision rnx, rny, rnz, rxnn
double precision tx, ty, tz, txn, txn0, t2x, t2y, t2z
double precision utau, upx, upy, upz, usn
double precision uiptn, uiptnf, uiptmn, uiptmx
double precision uetmax, uetmin, ukmax, ukmin, yplumx, yplumn
double precision uk, uet, nusury, yplus, unturb, dplus
double precision sqrcmu, clsyme, ek
double precision xnuii, xmutlm
double precision rcprod, rcflux
double precision cpp, rkl,  prdtl
double precision hflui, hredui, hint, hext, pimp
double precision und0, deuxd0
double precision eloglo(3,3), alpha(6,6)
double precision rcodcx, rcodcy, rcodcz, rcodsn
double precision visclc, visctc, romc  , distbf, srfbnf, cpscv
double precision cofimp, cofimpf, ypup
double precision ekip


integer          ntlast , iaff
data             ntlast , iaff /-1 , 0/
save             ntlast , iaff

!===============================================================================

!===============================================================================
! 1.  INITIALISATIONS
!===============================================================================

! Initialize variables to avoid compiler warnings

iclep = 0
iclk = 0
iclomg = 0
iclfb = 0
iclphi = 0
iclnu = 0
icl11 = 0
icl22 = 0
icl33 = 0
icl12 = 0
icl13 = 0
icl23 = 0
iclvar = 0
ipccv = 0

ek = 0.d0

! --- Memoire

! --- Constantes
uet = 1.d0
utau = 1.d0
sqrcmu = sqrt(cmu)

und0   = 1.d0
deuxd0 = 2.d0

! --- Conditions aux limites
iclu   = iclrtp(iu ,icoef)
iclv   = iclrtp(iv ,icoef)
iclw   = iclrtp(iw ,icoef)
if(itytur.eq.2) then
  iclk   = iclrtp(ik ,icoef)
  iclep  = iclrtp(iep,icoef)
elseif(itytur.eq.3) then
  icl11  = iclrtp(ir11,icoef)
  icl22  = iclrtp(ir22,icoef)
  icl33  = iclrtp(ir33,icoef)
  icl12  = iclrtp(ir12,icoef)
  icl13  = iclrtp(ir13,icoef)
  icl23  = iclrtp(ir23,icoef)
  iclep  = iclrtp(iep,icoef)
  if(iturb.eq.32) iclalp = iclrtp(ial, icoef)
elseif(itytur.eq.5) then
  iclk   = iclrtp(ik ,icoef)
  iclep  = iclrtp(iep,icoef)
  iclphi = iclrtp(iphi,icoef)
  if(iturb.eq.50) then
    iclfb  = iclrtp(ifb,icoef)
  elseif(iturb.eq.51) then
    iclal  = iclrtp(ial,icoef)
  endif
elseif(iturb.eq.60) then
  iclk   = iclrtp(ik ,icoef)
  iclomg = iclrtp(iomg,icoef)
elseif(iturb.eq.70) then
  iclnu  = iclrtp(inusa,icoef)
endif

icluf  = iclrtp(iu ,icoeff)
iclvf  = iclrtp(iv ,icoeff)
iclwf  = iclrtp(iw ,icoeff)

! --- Grandeurs physiques
ipcrom = ipproc(irom  )
ipcvis = ipproc(iviscl)
ipcvst = ipproc(ivisct)
if(icp.gt.0) then
  ipccp  = ipproc(icp   )
else
  ipccp = 0
endif

! --- Compressible

if ( ippmod(icompf) .ge. 0 ) then
  if(icv.gt.0) then
    ipccv  = ipproc(icv   )
  else
    ipccv = 0
  endif
endif

! MIN ET MAX DE LA VITESSE TANGENTIELLE EN PAROI
uiptmx = -grand
uiptmn =  grand

! MIN ET MAX DE LA VITESSE DE FROTTEMENT EN PAROI
uetmax = -grand
uetmin =  grand
ukmax  = -grand
ukmin  =  grand

! MIN ET MAX DE YPLUS
yplumx = -grand
yplumn =  grand

! COMPTEURS (TURBULENT, LAMINAIRE, RETOURNEMENT, CORRECTION
!            D'ECHELLE )
inturb = 0
inlami = 0
iuiptn = 0


!     En modele type v2f (phi-fbar et BL-v2/k) on met directement u=0 donc
!     UIPTMX et UIPTMN vaudront forcement 0
if (itytur.eq.5) then
  uiptmx = 0.d0
  uiptmn = 0.d0
endif

! --- Boucle sur les faces : debut
do ifac = 1, nfabor

! --- Test sur la presence d'une condition de paroi vitesse : debut
  if( icodcl(ifac,iu).eq.5 ) then

    iel = ifabor(ifac)

! --- Proprietes physiques
    visclc = propce(iel,ipcvis)
    visctc = propce(iel,ipcvst)
    romc   = propce(iel,ipcrom)

! --- Grandeurs geometriques
    distbf = distb(ifac)
    srfbnf = surfbn(ifac)

!===============================================================================
! 1. REPERE LOCAL
!===============================================================================


! ---> NORMALE UNITAIRE

    rnx = surfbo(1,ifac)/srfbnf
    rny = surfbo(2,ifac)/srfbnf
    rnz = surfbo(3,ifac)/srfbnf

! ---> PRISE EN COMPTE DE LA VITESSE DE DEFILEMENT

    rcodcx = rcodcl(ifac,iu,1)
    rcodcy = rcodcl(ifac,iv,1)
    rcodcz = rcodcl(ifac,iw,1)

!     Si on n'est pas en ALE, on force la vitesse de deplacement
!       de la face a etre tangentielle (et on met a jour rcodcl
!       pour une utilisation eventuelle)
    if (iale.eq.0.and.imobil.eq.0) then
      rcodsn = rcodcx*rnx+rcodcy*rny+rcodcz*rnz
      rcodcx = rcodcx -rcodsn*rnx
      rcodcy = rcodcy -rcodsn*rny
      rcodcz = rcodcz -rcodsn*rnz
      rcodcl(ifac,iu,1) = rcodcx
      rcodcl(ifac,iv,1) = rcodcy
      rcodcl(ifac,iw,1) = rcodcz
    endif


! ---> VITESSE TANGENTIELLE RELATIVE

    upx = coefu(ifac,1) - rcodcx
    upy = coefu(ifac,2) - rcodcy
    upz = coefu(ifac,3) - rcodcz

    usn = upx*rnx+upy*rny+upz*rnz
    tx  = upx -usn*rnx
    ty  = upy -usn*rny
    tz  = upz -usn*rnz
    txn = sqrt( tx**2 +ty**2 +tz**2 )
    utau= txn

! ---> TANGENTE UNITAIRE

    if( txn.ge.epzero) then

      txn0 = 1.d0

      tx  = tx/txn
      ty  = ty/txn
      tz  = tz/txn

    elseif(itytur.eq.3) then

!      SI LA VITESSE EST NULLE, LE VECTEUR T EST NORMAL ET QCQUE
!        ON EN A BESOIN POUR LE CHGT DE REPERE DE RIJ
!        ET ON ANNULERA LA VITESSE

      txn0 = 0.d0

      if(abs(rny).ge.epzero.or.abs(rnz).ge.epzero)then
        rxnn = sqrt(rny**2+rnz**2)
        tx  =  0.d0
        ty  =  rnz/rxnn
        tz  = -rny/rxnn
      elseif(abs(rnx).ge.epzero.or.abs(rnz).ge.epzero)then
        rxnn = sqrt(rnx**2+rnz**2)
        tx  =  rnz/rxnn
        ty  =  0.d0
        tz  = -rnx/rxnn
      else
        write(nfecra,1000)ifac,rnx,rny,rnz
        call csexit (1)
      endif

    else

!       SI LA VITESSE EST NULLE ET QU'ON N'EST PAS EN RIJ
!         TX, TY, TZ NE SERT PAS (ON ANNULE LA VITESSE)
!         ET ON LUI DONNE UNE VALEUR BIDON (NULLE PAR EXEMPLE)

      txn0 = 0.d0

      tx  = 0.d0
      ty  = 0.d0
      tz  = 0.d0

    endif

! ---> ON COMPLETE EVENTUELLEMENT POUR LE RIJ-EPSILON

    if (itytur.eq.3) then

!     --> T2 = RN X T (OU X EST LE PRODUIT VECTORIEL)


      t2x = rny*tz - rnz*ty
      t2y = rnz*tx - rnx*tz
      t2z = rnx*ty - rny*tx

!     --> MATRICE ORTHOGONALE DE CHANGEMENT DE BASE ELOGLOij
!         (DE LA BASE LOCALE VERS LA BASE GLOBALE)

!                            |TX  -RNX  T2X|
!                   ELOGLO = |TY  -RNY  T2Y|
!                            |TZ  -RNZ  T2Z|

!         SA TRANSPOSEE ELOGLOt EST SON INVERSE


      eloglo(1,1) =  tx
      eloglo(1,2) = -rnx
      eloglo(1,3) =  t2x
      eloglo(2,1) =  ty
      eloglo(2,2) = -rny
      eloglo(2,3) =  t2y
      eloglo(3,1) =  tz
      eloglo(3,2) = -rnz
      eloglo(3,3) =  t2z

!     --> ON CALCULE ALPHA(6,6)

!       SOIT f LE CENTRE DE LA FACE DE BORD ET
!            I LE CENTRE DE LA CELLULE CORRESPONDANTE

!       EN NOTE RG (RESP RL) INDICE PAR f OU PAR I
!          LE TENSEUR DE REYNOLDS DANS LA BASE GLOBALE (RESP LOCALE)

!       LA MATRICE ALPHA APPLIQUEE AU VECTEUR GLOBAL EN I'
!         (RG11,I'|RG22,I'|RG33,I'|RG12,I'|RG13,I'|RG23,I')t
!         DOIT DONNER LES VALEURS A IMPOSER A LA FACE
!         (RG11,f |RG22,f |RG33,f |RG12,f |RG13,f |RG23,f )t
!         AUX CONDITIONS LIMITES DE DIRICHLET PRES (AJOUTEES ENSUITE)

!       ON LA DEFINIT EN CALCULANT RG,f EN FONCTION DE RG,I' COMME SUIT

!         RG,f = ELOGLO.RL,f.ELOGLOt (PRODUITS MATRICIELS)

!                          | RL,I'(1,1)     B*U*.Uk     C*RL,I'(1,3) |
!           AVEC    RL,f = | B*U*.Uk       RL,I'(2,2)       0        |
!                          | C*RL,I'(1,3)     0         RL,I'(3,3)   |

!                  AVEC    RL,I = ELOGLOt.RG,I'.ELOGLO
!                          B = 0
!                    ET    C = 0 EN PAROI (1 EN SYMETRIE)



!          ON CALCULE EN FAIT   ELOGLO.PROJECTEUR.ELOGLOt


      clsyme=0.d0
      call clca66 ( clsyme , eloglo , alpha )
      !==========

    endif

!===============================================================================
! 2. VITESSES DE FROTTEMENT
!===============================================================================

! ---> ON CALCULE UET SUIVANT SI ON EST DANS LA ZONE LOG OU NON
!       EN UNE OU DEUX ECHELLES DE VITESSE
!       ET UK A PARTIR DE EK

    nusury = visclc/(distbf*romc)
! PSEUDO DECALAGE DE LA PAROI QUAND IDEUCH = 2
    dplus = 0.d0

    if (ideuch.eq.0) then

      if(ilogpo.eq.0) then
! AVEC LOI EN PUISSANCE (WERNER & WENGLE)
        uet = (utau/(apow*(1.0d0/nusury)**bpow))**dpow
      else
! AVEC LOI LOG
        imprim = max(iwarni(iu),2)
        xnuii = visclc/romc
        call causta                                               &
        !==========
      ( ifac  , imprim , xkappa , cstlog , ypluli ,        &
        apow  , bpow   , dpow   ,                                 &
        utau  , distbf , xnuii  , uet    )
      endif

! On reprend les deux lignes suivantes apres
!       l'appel eventuel a la fonction utilisateur
      uk = uet
      yplus = uk/nusury

    else
! Si IDEUCH=1 ou 2 on calcule uk et uet

      if (itytur.eq.2 .or. itytur.eq.5 .or. iturb.eq.60) then
        ek = rtp(iel,ik)
      else if (itytur.eq.3) then
        ek = 0.5d0*(rtp(iel,ir11)+rtp(iel,ir22)+rtp(iel,ir33))
      endif

      uk = cmu025*sqrt(ek)
      yplus = uk/nusury
      uet = utau/(log(yplus)/xkappa+cstlog)

    endif

    if(ideuch.eq.0) then
      uk = uet
      yplus = uk/nusury
    endif

! ---> ON TRAITE LES CAS OU YPLUS TEND VERS ZERO
! En une echelle, CAUSTA calcule d'abord u* en supposant une loi lineaire,
! puis si necessaire teste une loi log. Mais comme YPLULI est fixe a 1/kappa et pas
! 10,88 (valeur de continuite), la loi log peut donner une valeur de u* qui redonne
! un y+<YPLULI. Dans ce cas, on recalcule u* et y+ a partir de la loi lineaire : on
! obtient un y+ superieur a YPLULI, mais le frottement est sans doute correct.
! -> travail en cours sur les lois de paroi

    if (yplus.gt.ypluli) then

!       On est hors ss couche visqueuse : uet, uk et yplus sont bons
      unturb = 1.d0
      inturb = inturb + 1
    else

!       On est en sous couche visqueuse :

!       Si on utilise les "scalable wall functions", on decale la valeur de YPLUS,
!       on recalcule uet et on se considere hors de la sous-couche visqueuse
       if (ideuch.eq.2) then
          dplus = ypluli - yplus
          yplus = ypluli
          uet = utau/(log(yplus)/xkappa+cstlog)
          unturb = 1.d0
          inturb = inturb + 1
       else
!         Sinon on est reellement en sous-couche visqueuse
          unturb = 0.d0
          inlami = inlami + 1
!         On annule uk pour annuler les CL
          uk = 0.d0

!         On recalcule les valeurs fausses
!          en une  echelle  : uet et yplus sont faux
!          en deux echelles : uet est faux

!         En deux echelles :
!           On recalcule uet mais il ne sert plus a rien
!           (il intervient dans des termes multiplies par UNTURB=0)
          if(ideuch.eq.1) then
             if(yplus.gt.epzero)  then
                uet = abs(utau/yplus)
             else
                uet = 0.d0
             endif
!         En une echelle :
!           On recalcule uet : il sert pour la LES
!           On recalcule yplus (qui etait deduit de uet) : il sert pour hturbp
          else
             uet = sqrt(utau*nusury)
             yplus = uet/nusury
          endif
!       On est en ss couche visqueuse

       endif

    endif

    uetmax = max(uet,uetmax)
    uetmin = min(uet,uetmin)
    ukmax  = max(uk,ukmax)
    ukmin  = min(uk,ukmin)
    yplumx = max(yplus,yplumx)
    yplumn = min(yplus,yplumn)

! Sauvegarde de la vitesse de frottement et de la viscosite turbulente
! apres amortissement de van Driest pour la LES
! On n'amortit pas mu_t une seconde fois si on l'a deja fait
! (car une cellule peut avoir plusieurs faces de paroi)
! ou
! Sauvegarde de la vitesse de frottement et distance a la paroi yplus
! si le modele de depot de particules est active.

    if(itytur.eq.4.and.idries.eq.1) then
      uetbor(ifac) = uet
      if (visvdr(iel).lt.-900.d0) then
        propce(iel,ipcvst) = propce(iel,ipcvst)                   &
             *(1.d0-exp(-yplus/cdries))**2
        visvdr(iel) = propce(iel,ipcvst)
        visctc = propce(iel,ipcvst)
      endif
    else if (iilagr.gt.0.and.idepst.gt.0) then
      uetbor(ifac) = uet
    endif


! Sauvegarde de yplus si post traite

    if(mod(ipstdv,ipstyp).eq.0) then
      yplbr(ifac) = yplus
    endif


!===============================================================================
! 3. CONDITIONS AUX LIMITES SUR LA VITESSE
!===============================================================================

!              UIPTN  respecte la production de k
!              de facon conditionnelle    --> Coef RCPROD
!              UIPTNF respecte le flux
!               de facon conditionnelle   --> Coef RCFLUX

    if (itytur.eq.2 .or. iturb.eq.60) then

      xmutlm = xkappa*visclc*yplus

!     Si YPLUS=0, on met UIPTN et UIPTNF a 0 directement pour eviter les divisions
!     par 0. De toute facon, dans ce cas UNTURB=0
      if (yplus.gt.epzero) then
         rcprod =                                                 &
             min(xkappa , max(und0,sqrt(xmutlm/visctc))/yplus )
         rcflux = max(xmutlm,visctc)/(visclc+visctc)

         uiptn  = utau + distbf*uet*uk*romc/xkappa/visclc*(       &
              und0/(deuxd0*yplus-dplus) - deuxd0*rcprod )
         uiptnf = utau - distbf*uet*uk*romc/xmutlm*rcflux
      else
         uiptn = 0.d0
         uiptnf = 0.d0
      endif

      ! Coupled solving of the velocity components
      if (ivelco.eq.1) then
        if (yplus.ge.ypluli) then
          ! On implicite le terme de bord pour le gradient de vitesse
          ypup =  yplus/(log(yplus)/xkappa +cstlog)
          cofimp  = 1.d0 - ypup/xkappa*                        &
                           (und0/(deuxd0*yplus-dplus) - deuxd0*rcprod )
          ! On implicite le terme (rho*uet*uk)
          !NB mu/(mu+muT) est la car la viscosite turbulente au bord est non nulle
          cofimpf = visclc /(visclc+visctc)*ypup
        else
          ! Dans la sous couche visceuse : U_F=0
          cofimp  = 0.d0
          cofimpf = visclc/(visclc+visctc)
        endif
      endif

      uiptmx = max(uiptn*unturb,uiptmx)
      uiptmn = min(uiptn*unturb,uiptmn)
      if(uiptn*unturb.lt.-epzero) iuiptn = iuiptn + 1

      coefa(ifac,iclu)   = uiptn *tx*unturb *txn0
      coefa(ifac,iclv)   = uiptn *ty*unturb *txn0
      coefa(ifac,iclw)   = uiptn *tz*unturb *txn0
      coefa(ifac,icluf)  = uiptnf*tx*unturb *txn0
      coefa(ifac,iclvf)  = uiptnf*ty*unturb *txn0
      coefa(ifac,iclwf)  = uiptnf*tz*unturb *txn0

      coefa(ifac,iclu)   = coefa(ifac,iclu)  + rcodcx
      coefa(ifac,iclv)   = coefa(ifac,iclv)  + rcodcy
      coefa(ifac,iclw)   = coefa(ifac,iclw)  + rcodcz
      coefa(ifac,icluf)  = coefa(ifac,icluf) + rcodcx
      coefa(ifac,iclvf)  = coefa(ifac,iclvf) + rcodcy
      coefa(ifac,iclwf)  = coefa(ifac,iclwf) + rcodcz

      coefb(ifac,iclu)   = 0.d0
      coefb(ifac,iclv)   = 0.d0
      coefb(ifac,iclw)   = 0.d0
      coefb(ifac,icluf)  = 0.d0
      coefb(ifac,iclvf)  = 0.d0
      coefb(ifac,iclwf)  = 0.d0

      ! Coupled solving of the velocity components
      if (ivelco.eq.1) then

        coefau(1,ifac)     = rcodcx
        coefau(2,ifac)     = rcodcy
        coefau(3,ifac)     = rcodcz
        cofafu(1,ifac)     = rcodcx
        cofafu(2,ifac)     = rcodcy
        cofafu(3,ifac)     = rcodcz

        ! Projection in order to have the velocity parallel to the wall
        ! B = cofimp * ( IDENTITY - n x n )

        coefbu(1,1,ifac)   = cofimp*(1.d0-rnx**2)
        coefbu(2,2,ifac)   = cofimp*(1.d0-rny**2)
        coefbu(3,3,ifac)   = cofimp*(1.d0-rnz**2)
        coefbu(1,2,ifac)   = -cofimp*rnx*rny
        coefbu(1,3,ifac)   = -cofimp*rnx*rnz
        coefbu(2,1,ifac)   = -cofimp*rny*rnx
        coefbu(2,3,ifac)   = -cofimp*rny*rnz
        coefbu(3,1,ifac)   = -cofimp*rnz*rnx
        coefbu(3,2,ifac)   = -cofimp*rnz*rny

        ! Projection in order to have the shear stress parallel to the wall
        !  B = IDENTITY - cofimpf*( IDENTITY - n x n )

        cofbfu(1,1,ifac)   = 1.d0 - cofimpf*(1.d0-rnx**2)
        cofbfu(2,2,ifac)   = 1.d0 - cofimpf*(1.d0-rny**2)
        cofbfu(3,3,ifac)   = 1.d0 - cofimpf*(1.d0-rnz**2)

        cofbfu(1,2,ifac)   = cofimpf*rnx*rny
        cofbfu(1,3,ifac)   = cofimpf*rnx*rnz
        cofbfu(2,1,ifac)   = cofimpf*rny*rnx
        cofbfu(2,3,ifac)   = cofimpf*rny*rnz
        cofbfu(3,1,ifac)   = cofimpf*rnz*rnx
        cofbfu(3,2,ifac)   = cofimpf*rnz*rny
      endif

    elseif(iturb.eq.0 .or.iturb.eq.10.or.           &
           itytur.eq.3) then

      ! Dans le cadre de la ponderation elliptique, on ne doit pas
      ! tenir compte des lois de paroi. On fait donc un test sur le modele
      ! de turbulence :
      ! si on est en LRR ou SSG on laisse les lois de paroi, si on est en
      ! EBRSM, on impose l adherence.
      if (iturb.eq.32) then
        coefa(ifac,iclu)   = 0.d0
        coefa(ifac,iclv)   = 0.d0
        coefa(ifac,iclw)   = 0.d0

        uiptn = utau

        ! Coupled solving of the velocity components
        if (ivelco.eq.1) then
          cofimp  = 0.d0
          cofimpf = visclc/(visclc+visctc)
        endif

      else

        ! If ILOGPO=0, then IDEUCH=0
        if(ilogpo.eq.0) then
          uiptn  = utau                                             &
               + uet*apow*bpow*yplus**bpow*(2.d0**(bpow-1.d0)-2.d0)

          ! Coupled solving of the velocity components
          if (ivelco.eq.1) then
            if(yplus.ge.ypluli) then
              ! On implicite le terme de bord pour le gradient de vitesse
              ! Faute de calcul mis a part, a*yplus^b/U+ = uet^(b+1-1/d)
              ypup = utau**(2.d0*dpow-1.d0)/apow**(2.d0*dpow)
              cofimp  = 1.d0+bpow*uet**(bpow+1.d0-1.d0/dpow)*( 2.d0**(bpow-1.d0)-2.d0 )
              ! On implicite le terme (rho*uet*uk)
              !NB mu/(mu+muT) est la car la viscosite turbulente au bord est non nulle
              cofimpf = visclc /(visclc+visctc)*ypup
            else
              !Dans la sous couche visceuse : U_F=0
              cofimp  = 0.d0
              cofimpf = visclc/(visclc+visctc)
            endif
          endif

        else
          ! If YPLUS=0, UIPTN is set to 0 to avoid division by 0.
          ! By the way, in this case: UNTURB=0
          if (yplus.gt.epzero) then
            uiptn = utau - distbf*romc*uet*uk/xkappa/visclc*(     &
                  deuxd0/yplus - und0/(deuxd0*yplus-dplus) )
          else
            uiptn = 0.d0
          endif

          ! Coupled solving of the velocity components
          if (yplus.ge.ypluli) then
            ! On implicite le terme de bord pour le gradient de vitesse
            ypup   =  yplus/(log(yplus)/xkappa +cstlog)
            cofimp  = 1.d0-ypup/xkappa*(deuxd0/yplus - und0/(deuxd0*yplus-dplus) )
            ! On implicite le terme (rho*uet*uk)
            !NB mu/(mu+muT) est la car la viscosite turbulente au bord est non nulle
            cofimpf = visclc /(visclc+visctc)*ypup
          else
            !Dans la sous couche visceuse : U_F=0
            cofimp  = 0.d0
            cofimpf = visclc/(visclc+visctc)
          endif

        endif

        coefa(ifac,iclu)   = uiptn *tx*unturb *txn0
        coefa(ifac,iclv)   = uiptn *ty*unturb *txn0
        coefa(ifac,iclw)   = uiptn *tz*unturb *txn0
      endif

      uiptmx = max(uiptn*unturb,uiptmx)
      uiptmn = min(uiptn*unturb,uiptmn)
      if(uiptn*unturb.lt.-epzero) iuiptn = iuiptn + 1

      coefa(ifac,iclu)   = coefa(ifac,iclu)  + rcodcx
      coefa(ifac,iclv)   = coefa(ifac,iclv)  + rcodcy
      coefa(ifac,iclw)   = coefa(ifac,iclw)  + rcodcz

      coefb(ifac,iclu)   = 0.d0
      coefb(ifac,iclv)   = 0.d0
      coefb(ifac,iclw)   = 0.d0

      ! Coupled solving of the velocity components
      if (ivelco.eq.1) then

        coefau(1,ifac)     = rcodcx
        coefau(2,ifac)     = rcodcy
        coefau(3,ifac)     = rcodcz
        cofafu(1,ifac)     = rcodcx
        cofafu(2,ifac)     = rcodcy
        cofafu(3,ifac)     = rcodcz

        ! Projection in order to have the velocity parallel to the wall
        ! B = cofimp * ( IDENTITY - n x n )

        coefbu(1,1,ifac)   = cofimp*(1.d0-rnx**2)
        coefbu(2,2,ifac)   = cofimp*(1.d0-rny**2)
        coefbu(3,3,ifac)   = cofimp*(1.d0-rnz**2)
        coefbu(1,2,ifac)   = -cofimp*rnx*rny
        coefbu(1,3,ifac)   = -cofimp*rnx*rnz
        coefbu(2,1,ifac)   = -cofimp*rny*rnx
        coefbu(2,3,ifac)   = -cofimp*rny*rnz
        coefbu(3,1,ifac)   = -cofimp*rnz*rnx
        coefbu(3,2,ifac)   = -cofimp*rnz*rny

        ! Projection in order to have the shear stress parallel to the wall
        !  B = IDENTITY - cofimpf*( IDENTITY - n x n )

        cofbfu(1,1,ifac)   = 1.d0 - cofimpf*(1.d0-rnx**2)
        cofbfu(2,2,ifac)   = 1.d0 - cofimpf*(1.d0-rny**2)
        cofbfu(3,3,ifac)   = 1.d0 - cofimpf*(1.d0-rnz**2)

        cofbfu(1,2,ifac)   = cofimpf*rnx*rny
        cofbfu(1,3,ifac)   = cofimpf*rnx*rnz
        cofbfu(2,1,ifac)   = cofimpf*rny*rnx
        cofbfu(2,3,ifac)   = cofimpf*rny*rnz
        cofbfu(3,1,ifac)   = cofimpf*rnz*rnx
        cofbfu(3,2,ifac)   = cofimpf*rnz*rny

      endif

    ! En LES on est forcement en IDEUCH=0, pas la peine d'exprimer les flux en
    ! version "scalable wall function". Idem pour le modele de Spalart Allmaras.
    elseif(itytur.eq.4.or.iturb.eq.70) then
      if(ilogpo.eq.0) then
        uiptn  = utau                                             &
             + uet*apow*bpow*yplus**bpow*(2.d0**(bpow-1.d0)-2.d0)

        ! Coupled solving of the velocity components
        if (ivelco.eq.1) then
          ! On implicite le terme de bord pour le gradient de vitesse
          ! Faute de calcul mis a part, a*yplus^b/U+ = uet^(b+1-1/d)
          ypup = utau**(2.d0*dpow-1.d0)/apow**(2.d0*dpow)
          cofimp  = 1.d0+bpow*uet**(bpow+1.d0-1.d0/dpow)*( 2.d0**(bpow-1.d0)-2.d0 )
          ! On implicite le terme (rho*uet*uk)
          !NB mu/(mu+muT) est la car la viscosite turbulente au bord est non nulle
          cofimpf = visclc /(visclc+visctc)*ypup
        endif

      else
        uiptn  = utau - uet/xkappa*1.5d0

        ! Coupled solving of the velocity components
        if (ivelco.eq.1) then
          if (uet.gt.epzero) then
            cofimp = 1.d0 - visclc/romc/uet/(xkappa*distbf)*1.5d0
          else
            cofimp = 1.d0
          endif

        endif
      endif

      ! If (mu+mut) becomes zero (dynamic models), an arbitrary value is set
      ! (nul flux) but without any problems because teh flux is really zero at this face.
      if(visctc+visclc.le.0) then
        uiptnf = utau
        if (ivelco.eq.1) cofimpf= 1.d0
      else

        ! Coupled solving of the velocity components
        if (ivelco.eq.1) then
          if(yplus.ge.ypluli) then
            ! On implicite le terme de bord pour le gradient de vitesse
            ypup   =  yplus/(log(yplus)/xkappa +cstlog)
            ! On implicite le terme (rho*uet*uk)
            !NB mu/(mu+muT) est la car la viscosite turbulente au bord est non nulle
            cofimpf = visclc /(visclc+visctc)*ypup
          else
            !Dans la sous couche visceuse : U_F=0
            cofimpf = visclc/(visclc+visctc)
          endif
        endif

        uiptnf = utau -romc*distbf*(uet**2)/(visctc+visclc)
      endif

      uiptmx = max(uiptn*unturb,uiptmx)
      uiptmn = min(uiptn*unturb,uiptmn)
      if(uiptn*unturb.lt.-epzero) iuiptn = iuiptn + 1

      coefa(ifac,iclu)   = uiptn *tx*unturb *txn0
      coefa(ifac,iclv)   = uiptn *ty*unturb *txn0
      coefa(ifac,iclw)   = uiptn *tz*unturb *txn0
      coefa(ifac,icluf)  = uiptnf*tx*unturb *txn0
      coefa(ifac,iclvf)  = uiptnf*ty*unturb *txn0
      coefa(ifac,iclwf)  = uiptnf*tz*unturb *txn0

      coefa(ifac,iclu)   = coefa(ifac,iclu)  + rcodcx
      coefa(ifac,iclv)   = coefa(ifac,iclv)  + rcodcy
      coefa(ifac,iclw)   = coefa(ifac,iclw)  + rcodcz
      coefa(ifac,icluf)  = coefa(ifac,icluf) + rcodcx
      coefa(ifac,iclvf)  = coefa(ifac,iclvf) + rcodcy
      coefa(ifac,iclwf)  = coefa(ifac,iclwf) + rcodcz

      coefb(ifac,iclu)   = 0.d0
      coefb(ifac,iclv)   = 0.d0
      coefb(ifac,iclw)   = 0.d0
      coefb(ifac,icluf)  = 0.d0
      coefb(ifac,iclvf)  = 0.d0
      coefb(ifac,iclwf)  = 0.d0

      ! Coupled solving of the velocity components
      if (ivelco.eq.1) then

        coefau(1,ifac)     = rcodcx
        coefau(2,ifac)     = rcodcy
        coefau(3,ifac)     = rcodcz
        cofafu(1,ifac)     = rcodcx
        cofafu(2,ifac)     = rcodcy
        cofafu(3,ifac)     = rcodcz

        ! Projection in order to have the velocity parallel to the wall
        ! B = cofimp * ( IDENTITY - n x n )

        coefbu(1,1,ifac)   = cofimp*(1.d0-rnx**2)
        coefbu(2,2,ifac)   = cofimp*(1.d0-rny**2)
        coefbu(3,3,ifac)   = cofimp*(1.d0-rnz**2)
        coefbu(1,2,ifac)   = -cofimp*rnx*rny
        coefbu(1,3,ifac)   = -cofimp*rnx*rnz
        coefbu(2,1,ifac)   = -cofimp*rny*rnx
        coefbu(2,3,ifac)   = -cofimp*rny*rnz
        coefbu(3,1,ifac)   = -cofimp*rnz*rnx
        coefbu(3,2,ifac)   = -cofimp*rnz*rny

        ! Projection in order to have the shear stress parallel to the wall
        !  B = IDENTITY - cofimpf*( IDENTITY - n x n )

        cofbfu(1,1,ifac)   = 1.d0 - cofimpf*(1.d0-rnx**2)
        cofbfu(2,2,ifac)   = 1.d0 - cofimpf*(1.d0-rny**2)
        cofbfu(3,3,ifac)   = 1.d0 - cofimpf*(1.d0-rnz**2)

        cofbfu(1,2,ifac)   = cofimpf*rnx*rny
        cofbfu(1,3,ifac)   = cofimpf*rnx*rnz
        cofbfu(2,1,ifac)   = cofimpf*rny*rnx
        cofbfu(2,3,ifac)   = cofimpf*rny*rnz
        cofbfu(3,1,ifac)   = cofimpf*rnz*rnx
        cofbfu(3,2,ifac)   = cofimpf*rnz*rny

      endif

    elseif (itytur.eq.5) then

!     Avec ces conditions, pas besoin de calculer UIPTMX, UIPTMN
!     et IUIPTN qui sont nuls (valeur d'initialisation)

      coefa(ifac,iclu)   = rcodcx
      coefa(ifac,iclv)   = rcodcy
      coefa(ifac,iclw)   = rcodcz

      coefb(ifac,iclu)   = 0.d0
      coefb(ifac,iclv)   = 0.d0
      coefb(ifac,iclw)   = 0.d0

      ! Coupled solving of the velocity components
      if(ivelco.eq.1) then

        coefau(1,ifac)     = coefa(ifac,iclu)
        coefau(2,ifac)     = coefa(ifac,iclv)
        coefau(3,ifac)     = coefa(ifac,iclw)
        cofafu(1,ifac)     = coefa(ifac,iclu)
        cofafu(2,ifac)     = coefa(ifac,iclv)
        cofafu(3,ifac)     = coefa(ifac,iclw)

        coefbu(1,1,ifac)   = 0.d0
        coefbu(2,2,ifac)   = 0.d0
        coefbu(3,3,ifac)   = 0.d0

        coefbu(1,2,ifac)   = 0.d0
        coefbu(1,3,ifac)   = 0.d0
        coefbu(2,1,ifac)   = 0.d0
        coefbu(2,3,ifac)   = 0.d0
        coefbu(3,1,ifac)   = 0.d0
        coefbu(3,2,ifac)   = 0.d0

        ! Projection in order to have the shear stress parallel to the wall
        !  B = IDENTITY - cofimpf*( IDENTITY - n x n )

        cofbfu(1,1,ifac)   = rnx**2
        cofbfu(2,2,ifac)   = rny**2
        cofbfu(3,3,ifac)   = rnz**2

        cofbfu(1,2,ifac)   = rnx*rny
        cofbfu(1,3,ifac)   = rnx*rnz
        cofbfu(2,1,ifac)   = rny*rnx
        cofbfu(2,3,ifac)   = rny*rnz
        cofbfu(3,1,ifac)   = rnz*rnx
        cofbfu(3,2,ifac)   = rnz*rny

      endif

    endif

!===============================================================================
! 4. CONDITIONS AUX LIMITES SUR K ET EPSILON
!===============================================================================

    if (itytur.eq.2) then

      coefa(ifac,iclk)   = uk**2/sqrcmu
      coefb(ifac,iclk)   = 0.d0

!     Si YPLUS=0, on met COEFA a 0 directement pour eviter une division
!     par 0.
      if (yplus.gt.epzero) then
         coefa(ifac,iclep) = distbf*4.d0*uk**5*romc**2/           &
              (xkappa*visclc**2*(yplus+dplus)**2)
      else
         coefa(ifac,iclep) = 0.d0
      endif
      coefb(ifac,iclep)  = 1.d0

!===============================================================================
! 5. CONDITIONS AUX LIMITES SUR RIJ ET EPSILON
!===============================================================================

    elseif (itytur.eq.3) then

! ---> TENSEUR RIJ (PARTIELLEMENT IMPLICITE)

      do isou = 1, 6

        if(isou.eq.1) iclvar = icl11
        if(isou.eq.2) iclvar = icl22
        if(isou.eq.3) iclvar = icl33
        if(isou.eq.4) iclvar = icl12
        if(isou.eq.5) iclvar = icl13
        if(isou.eq.6) iclvar = icl23

        coefa(ifac,iclvar) = 0.0d0
        coefb(ifac,iclvar) = 0.0d0

      enddo

      ! for the LRR and the Standard SGG.
      if ((iturb.eq.30).or.(iturb.eq.31)) then

        do isou = 1,6

          if(isou.eq.1) then
            iclvar = icl11
            jj = 1
            kk = 1
          else if(isou.eq.2) then
            iclvar = icl22
            jj = 2
            kk = 2
          else if(isou.eq.3) then
            iclvar = icl33
            jj = 3
            kk = 3
          else if(isou.eq.4) then
            iclvar = icl12
            jj = 1
            kk = 2
          else if(isou.eq.5) then
            iclvar = icl13
            jj = 1
            kk = 3
          else if(isou.eq.6) then
            iclvar = icl23
            jj = 2
            kk = 3
          endif

          if (iclptr.eq.1) then
            do ii = 1, 6
              if(ii.ne.isou) then
                coefa(ifac,iclvar) = coefa(ifac,iclvar) +           &
                   alpha(isou,ii) * rijipb(ifac,ii)
              endif
            enddo
            coefb(ifac,iclvar) = alpha(isou,isou)
          else
            do ii = 1, 6
              coefa(ifac,iclvar) = coefa(ifac,iclvar) +             &
                   alpha(isou,ii) * rijipb(ifac,ii)
            enddo
            coefb(ifac,iclvar) = 0.d0
          endif

          coefa(ifac,iclvar) = coefa(ifac,iclvar)  -                &
                             (eloglo(jj,1)*eloglo(kk,2)+            &
                              eloglo(jj,2)*eloglo(kk,1))*uet*uk

          ! If laminar: zero Reynolds' stresses

          if (unturb.le.epzero) then
            coefa(ifac,iclvar) = 0.d0
            coefb(ifac,iclvar) = 0.d0
          endif

        enddo

      endif

! ---> SCALAIRE EPSILON
!      ICI AUSSI, POSSIBILITE DE FORME FLUX OU DIRICHLET ...
!      ON NE RECONSTRUIT PAS EPSILON
!      POSSIBILITE D'IMPLICITATION PARTIELLE

!     Si YPLUS=0, on met COEFA a 0 directement pour eviter une division
!     par 0.
      if ((iturb.eq.30).or.(iturb.eq.31)) then
        if (yplus.gt.epzero) then
           coefa(ifac,iclep) = distbf*4.d0*uk**5*romc**2/           &
                (xkappa*visclc**2*(yplus+dplus)**2)
        else
           coefa(ifac,iclep) = 0.d0
        endif
        if (iclptr.eq.1) then
          coefb(ifac,iclep) = 1.d0
        else
          coefa(ifac,iclep) = rtp(iel,iclep) + coefa(ifac,iclep)
          coefb(ifac,iclep) = 0.d0
        endif
      elseif (iturb.eq.32) then
        ! Use k at I'
        ekip = 0.5d0*(rijipb(ifac,1)+rijipb(ifac,2)+rijipb(ifac,3))

        coefa(ifac,iclep) = 2.d0*propce(iel,ipcvis)*ekip /       &
                            (distbf**2*propce(iel,ipcrom))
        coefb(ifac,iclep) = 0.d0

        ! Alpha
        coefa(ifac,iclalp) = 0.0d0
        coefb(ifac,iclalp) = 0.0d0

      endif

!===============================================================================
! 6a.CONDITIONS AUX LIMITES SUR K, EPSILON, F_BARRE ET PHI
!    DANS LE MODELE PHI_FBAR
!===============================================================================

    elseif (iturb.eq.50) then

      coefa(ifac,iclk) = 0.d0
      coefb(ifac,iclk) = 0.d0
      coefa(ifac,iclep) =                                         &
           2.0d0*propce(iel,ipcvis)/propce(iel,ipcrom)            &
           *rtp(iel,ik)/distbf**2
      coefb(ifac,iclep) = 0.d0
      coefa(ifac,iclphi) = 0.0d0
      coefb(ifac,iclphi) = 0.0d0
      coefa(ifac,iclfb) = 0.0d0
      coefb(ifac,iclfb) = 0.0d0

!===============================================================================
! 6b.CONDITIONS AUX LIMITES SUR K, EPSILON, PHI ET ALPHA
!    DANS LE MODELE BL-V2/K
!===============================================================================

    elseif (iturb.eq.51) then

      coefa(ifac,iclk) = 0.d0
      coefb(ifac,iclk) = 0.d0
      coefa(ifac,iclep) =                                         &
                 propce(iel,ipcvis)/propce(iel,ipcrom)            &
           *rtp(iel,ik)/distbf**2
      coefb(ifac,iclep) = 0.d0
      coefa(ifac,iclphi) = 0.0d0
      coefb(ifac,iclphi) = 0.0d0
      coefa(ifac,iclal) = 0.0d0
      coefb(ifac,iclal) = 0.0d0

!===============================================================================
! 7. CONDITIONS AUX LIMITES SUR K ET OMEGA
!===============================================================================

    elseif (iturb.eq.60) then

!     Si on est hors de la sous-couche visqueuse (reellement ou via les
!     scalable wall functions)
      if (unturb.eq.1) then
        coefa(ifac,iclk)   = uk**2/sqrcmu
        coefb(ifac,iclk)   = 0.d0
!     Comme UNTURB=1 YPLUS est forcement >0
        coefa(ifac,iclomg) = distbf*4.d0*uk**3*romc**2/           &
              (sqrcmu*xkappa*visclc**2*(yplus+dplus)**2)
        coefb(ifac,iclomg) = 1.d0

      else
!     Si on est en sous-couche visqueuse
        coefa(ifac,iclk)   = 0.d0
        coefb(ifac,iclk)   = 0.d0
        coefa(ifac,iclomg) = distbf*120.d0*8.d0*visclc/romc       &
             /(ckwbt1*distbf**3)
        coefb(ifac,iclomg)  = 1.d0
      endif

!===============================================================================
! 7.1 CONDITIONS AUX LIMITES SUR LE MODELE DE SPALART ALLMARAS
!===============================================================================

    elseif (iturb.eq.70) then

      coefa(ifac,iclnu)   = 0.d0
      coefb(ifac,iclnu)   = 0.d0

    endif

!===============================================================================
! 8. CONDITIONS AUX LIMITES SUR LES SCALAIRES
!              (AUTRES QUE PRESSION, K, EPSILON, RIJ, VARIANCES)
!    Pour les variances, pas de traitement specifique en paroi : voir
!      condli.
!===============================================================================

    if(nscal.ge.1) then

      do ll = 1, nscal

        if(iscavr(ll).le.0) then

          ivar = isca(ll)
          iclvar = iclrtp(ivar,icoef)
          iclvaf = iclrtp(ivar,icoeff)

          isvhbl = 0
          if(ll.eq.isvhb) then
            isvhbl = isvhb
          endif

          ihcp = 0
          iscal = ll
          if(iscsth(iscal).eq.0.or.iscsth(iscal).eq.2             &
                               .or.iscsth(iscal).eq.3) then
            ihcp = 0
          elseif(abs(iscsth(iscal)).eq.1) then
            if(ipccp.gt.0) then
              ihcp = 2
            else
              ihcp = 1
            endif
          endif

          cpp = 1.d0
          if(ihcp.eq.0) then
            cpp = 1.d0
          elseif(ihcp.eq.2) then
            cpp = propce(iel,ipccp )
          elseif(ihcp.eq.1) then
            cpp = cp0
          endif
          hint = cpp

          if(ivisls(ll).gt.0) then
            ipcvsl = ipproc(ivisls(ll))
          else
            ipcvsl = 0
          endif
          if (ipcvsl.le.0) then
            rkl = visls0(ll)
            prdtl = visclc/rkl
          else
            rkl = propce(iel,ipcvsl)
            prdtl = visclc/rkl
          endif

!  Compressible : On suppose que le nombre de Pr doit etre
!               defini de la meme facon que l'on resolve
!               en enthalpie ou en energie, soit Mu*Cp/Lambda.
!               Si l'on resout en energie, on a calcule ci-dessus
!               Mu*Cv/Lambda.

          if ( ippmod(icompf).ge.0 ) then
            if(iscsth(iscal).eq.3) then
              if(ipccp.gt.0) then
                prdtl = prdtl*propce(iel,ipccp )
              else
                prdtl = prdtl*cp0
              endif
              if(ipccv.gt.0) then
                prdtl = prdtl/propce(iel,ipccv )
              else
                prdtl = prdtl/cv0
              endif
            endif
          endif

!          CAS TURBULENT
          if (iturb.ne.0) then
            if ( ippmod(icompf) .ge. 0 ) then
!                 En compressible, pour l'energie LAMBDA/CV+CP/CV*(MUT/SIGMAS)
              if(ipccp.gt.0) then
                cpscv = propce(iel,ipproc(icp))
              else
                cpscv = cp0
              endif
              if(ipccv.gt.0) then
                cpscv = cpscv/propce(iel,ipproc(icv))
              else
                cpscv = cpscv/cv0
              endif
              hint = hint*(rkl+cpscv*visctc/sigmas(ll))/distbf
            else
              hint = hint*(rkl+visctc/sigmas(ll))/distbf
            endif
!          CAS LAMINAIRE
          else
            hint  = hint*rkl/distbf
          endif

          if(iturb.ne.0.and.icodcl(ifac,ivar).eq.5)then
            call hturbp (prdtl,sigmas(ll),xkappa,yplus,hflui)
            !==========
            if (ideuch.eq.2) then
              hflui = cpp*uk*romc/(yplus*prdtl) *hflui
            else
              hflui = cpp*rkl/distbf *hflui
            endif
          else
            hflui = hint
          endif

          if (isvhbl .gt. 0) hbord(ifac) = hflui



! --->  C.L DE TYPE DIRICHLET AVEC OU SANS COEFFICIENT D'ECHANGE

!     Si on a deux types de conditions aux limites (ICLVAR, ICLVAF)
!       il faut que le flux soit traduit par ICLVAF.
!     Si on n'a qu'un type de condition, peu importe (ICLVAF=ICLVAR)
!     Pour le moment, dans cette version compressible, on impose un
!       flux nul pour ICLVAR, lorsqu'il est different de ICLVAF (cette
!       condition ne sert qu'a la reconstruction des gradients et
!       s'applique a l'energie totale qui inclut l'energie cinetique :


          if( icodcl(ifac,ivar).eq.5 ) then
            hext = rcodcl(ifac,ivar,2)
            pimp = rcodcl(ifac,ivar,1)
            hredui = hint/hflui
            coefa(ifac,iclvaf) = hext*pimp/(hint+hext*hredui)
            coefb(ifac,iclvaf) = (hint-(1.d0-hredui)*hext)/       &
                                 (hint+hext*hredui)
            if(iclvar.ne.iclvaf) then
              coefa(ifac,iclvar) = 0.d0
              coefb(ifac,iclvar) = 1.d0
            endif

!--> Rayonnement :

!      On stocke le coefficient d'echange lambda/distance
!      (ou son equivalent en turbulent) quelle que soit la
!      variable thermique transportee (temperature ou enthalpie)
!      car on l'utilise pour realiser des bilans aux parois qui
!      sont faits en temperature (on cherche la temperature de
!      paroi quelle que soit la variable thermique transportee pour
!      ecrire des eps sigma T4.

!     donc :

!       lorsque la variable transportee est la temperature
!         ABS(ISCSTH(II)).EQ.1 : RA(IHCONV-1+IFAC) = HINT
!         puisque HINT = VISLS * CP / DISTBR
!                      = lambda/distance en W/(m2 K)

!       lorsque la variable transportee est l'enthalpie
!         ISCSTH(II).EQ.2 : RA(IHCONV-1+IFAC) = HINT*CPR
!         avec
!            IF(IPCCP.GT.0) THEN
!              CPR = PROPCE(IEL,IPCCP )
!            ELSE
!              CPR = CP0
!            ENDIF
!         puisque HINT = VISLS / DISTBR
!                      = lambda/(CP * distance)

!       lorsque la variable transportee est l'energie (compressible)
!         ISCSTH(II).EQ.3 :
!         on procede comme pour l'enthalpie avec CV au lieu de CP
!         (rq : il n'y a pas d'hypothese, sf en non orthogonal :
!               le flux est le bon et le coef d'echange aussi)

!      De meme dans condli.



!               Si on rayonne et que
!                  le scalaire est la variable energetique

            if (iirayo.ge.1 .and. ll.eq.iscalt) then

!                On calcule le coefficient d'echange en W/(m2 K)

!                Si on resout en enthalpie
              if(iscsth(ll).eq.2) then
!                  Si Cp variable
                if(ipccp.gt.0) then
                  propfb(ifac,ipprob(ihconv)) = hflui*propce(iel,ipccp )
                else
                  propfb(ifac,ipprob(ihconv)) = hflui*cp0
                endif

!                  Si on resout en energie (compressible)
              elseif(iscsth(ll).eq.3) then
!                    Si Cv variable
                if(ipccv.gt.0) then
                  propfb(ifac,ipprob(ihconv)) = hflui*propce(iel,ipccv )
                else
                  propfb(ifac,ipprob(ihconv)) = hflui*cv0
                endif

!                Si on resout en temperature
              elseif(abs(iscsth(ll)).eq.1) then
                propfb(ifac,ipprob(ihconv)) = hflui
              endif

!                On recupere le flux h(Ti'-Tp) (sortant ou
!                             negatif si gain pour le fluide) en W/m2

              propfb(ifac,ipprob(ifconv)) =                       &
                   hint*( (1.d0-coefb(ifac,iclvaf))*thbord(ifac)  &
                         - coefa(ifac,iclvaf))
            endif

          endif

! --->  C.L DE TYPE FLUX : VOIR CONDLI

        endif

      enddo

    endif


  endif
! --- Test sur la presence d'une condition de paroi vitesse : fin



enddo
! --- Boucle sur les faces : fin

if (irangp.ge.0) then
  call parmin (uiptmn)
  !==========
  call parmax (uiptmx)
  !==========
  call parmin (uetmin)
  !==========
  call parmax (uetmax)
  !==========
  call parmin (ukmin)
  !==========
  call parmax (ukmax)
  !==========
  call parmin (yplumn)
  !==========
  call parmax (yplumx)
  !==========
  call parcpt (inturb)
  !==========
  call parcpt (inlami)
  !==========
  call parcpt (iuiptn)
  !==========
endif

!===============================================================================
! 9.  IMPRESSIONS
!===============================================================================

!     Remarque : afin de ne pas surcharger les listings dans le cas ou
!       quelques yplus ne sont pas corrects, on ne produit le message
!       qu'aux deux premiers pas de temps ou le message apparait et
!       aux deux derniers pas de temps du calcul, ou si IWARNI est >= 2
!       On indique aussi le numero du dernier pas de temps auquel on
!       a rencontre des yplus hors bornes admissibles

if(iwarni(iu).ge.0) then
  if(ntlist.gt.0) then
    modntl = mod(ntcabs,ntlist)
  elseif(ntlist.eq.-1.and.ntcabs.eq.ntmabs) then
    modntl = 0
  else
    modntl = 1
  endif

  if ( (iturb.eq.0.and.inturb.ne.0)        .or.            &
       (itytur.eq.5.and.inturb.ne.0)       .or.            &
       ((itytur.eq.2.or.itytur.eq.3)                &
        .and.inlami.gt.0)                             )           &
       ntlast = ntcabs

  if ( (ntlast.eq.ntcabs.and.iaff.lt.2         ).or.              &
       (ntlast.ge.0     .and.ntcabs.ge.ntmabs-1).or.              &
       (ntlast.eq.ntcabs.and.iwarni(iu).ge.2)    ) then
    iaff = iaff + 1
    write(nfecra,2010) &
         uiptmn,uiptmx,uetmin,uetmax,ukmin,ukmax,yplumn,yplumx,   &
         iuiptn,inlami,inlami+inturb
    if (iturb.eq. 0)                                       &
         write(nfecra,2020)  ntlast,ypluli
    if (itytur.eq.5)                                       &
         write(nfecra,2030)  ntlast,ypluli
    ! No warnings in EBRSM
    if (itytur.eq.2.or.iturb.eq.30.or.iturb.eq.31)         &
         write(nfecra,2040)  ntlast,ypluli
    if (iwarni(iu).lt.2.and.iturb.ne.32) then
      write(nfecra,2050)
    elseif (iturb.ne.32) then
      write(nfecra,2060)
    endif

  else if (modntl.eq.0 .or. iwarni(iu).ge.2) then
    write(nfecra,2010) &
         uiptmn,uiptmx,uetmin,uetmax,ukmin,ukmax,yplumn,yplumx,   &
         iuiptn,inlami,inlami+inturb
  endif

endif

!===============================================================================
! 10.  FORMATS
!===============================================================================

#if defined(_CS_LANG_FR)

 1000 format(/,' LA NORMALE A LA FACE DE BORD DE PAROI ',I10,/,   &
         ' EST NULLE ; COORDONNEES : ',3E12.5)

 2010 format(/,                                                   &
 3X,'** CONDITIONS AUX LIMITES EN PAROI LISSE',/,           &
 '   ----------------------------------------',/,           &
 '------------------------------------------------------------',/,&
 '                                         Minimum     Maximum',/,&
 '------------------------------------------------------------',/,&
 '   Vitesse rel. en paroi    uiptn : ',2E12.5                 ,/,&
 '   Vitesse de frottement    uet   : ',2E12.5                 ,/,&
 '   Vitesse de frottement    uk    : ',2E12.5                 ,/,&
 '   Distance adimensionnelle yplus : ',2E12.5                 ,/,&
 '   ------------------------------------------------------   ',/,&
 '   Nbre de retournements de la vitesse en paroi : ',I10      ,/,&
 '   Nbre de faces en sous couche visqueuse       : ',I10      ,/,&
 '   Nbre de faces de paroi total                 : ',I10      ,/,&
 '------------------------------------------------------------',  &
 /,/)

 2020 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : RAFFINEMENT INSUFFISANT DU MAILLAGE EN PAROI',/,&
'@    =========                                               ',/,&
'@    Le maillage semble insuffisamment raffine en paroi      ',/,&
'@      pour pouvoir realiser un calcul laminaire.            ',/,&
'@                                                            ',/,&
'@    Le dernier pas de temps auquel ont ete observees de trop',/,&
'@      grandes valeurs de la distance adimensionnelle a la   ',/,&
'@      paroi (yplus) est le pas de temps ',I10                ,/,&
'@                                                            ',/,&
'@    La valeur minimale de yplus doit etre inferieure a la   ',/,&
'@      valeur limite YPLULI = ',E14.5                         ,/,&
'@                                                            ',/,&
'@    Observer la repartition de yplus en paroi (sous Ensight ',/,&
'@      par exemple) pour determiner dans quelle mesure la    ',/,&
'@      qualite des resultats est susceptible d etre affectee.')

 2030 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : RAFFINEMENT INSUFFISANT DU MAILLAGE EN PAROI',/,&
'@    =========                                               ',/,&
'@    Le maillage semble insuffisamment raffine en paroi      ',/,&
'@      pour pouvoir realiser un calcul type v2f              ',/,&
'@            (phi-fbar ou BL-v2/k)                           ',/,&
'@                                                            ',/,&
'@    Le dernier pas de temps auquel ont ete observees de trop',/,&
'@      grandes valeurs de la distance adimensionnelle a la   ',/,&
'@      paroi (yplus) est le pas de temps ',I10                ,/,&
'@                                                            ',/,&
'@    La valeur minimale de yplus doit etre inferieure a la   ',/,&
'@      valeur limite YPLULI = ',E14.5                         ,/,&
'@                                                            ',/,&
'@    Observer la repartition de yplus en paroi (sous Ensight ',/,&
'@      par exemple) pour determiner dans quelle mesure la    ',/,&
'@      qualite des resultats est susceptible d etre affectee.')

 2040 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : MAILLAGE TROP FIN EN PAROI                  ',/,&
'@    =========                                               ',/,&
'@    Le maillage semble trop raffine en paroi pour utiliser  ',/,&
'@      un modele de turbulence haut Reynolds.                ',/,&
'@                                                            ',/,&
'@    Le dernier pas de temps auquel ont ete observees des    ',/,&
'@      valeurs trop faibles de la distance adimensionnelle a ',/,&
'@      la paroi (yplus) est le pas de temps ',I10             ,/,&
'@                                                            ',/,&
'@    La valeur minimale de yplus doit etre superieure a la   ',/,&
'@      valeur limite YPLULI = ',E14.5                         ,/,&
'@                                                            ',/,&
'@    Observer la repartition de yplus en paroi (sous Ensight ',/,&
'@      par exemple) pour determiner dans quelle mesure la    ',/,&
'@      qualite des resultats est susceptible d etre affectee.')
 2050 format(                                                           &
'@                                                            ',/,&
'@    Ce message ne s''affiche qu''aux deux premieres         ',/,&
'@      occurences du probleme et aux deux derniers pas de    ',/,&
'@      temps du calcul. La disparition du message ne signifie',/,&
'@      pas forcement la disparition du probleme.             ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2060 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

#else

 1000 format(/,' THE NORMAL TO THE WALL BOUNDARY FACE ',I10,/,    &
         ' IS NULL; COORDINATES: ',3E12.5)

 2010 format(/,                                                   &
 3X,'** BOUNDARY CONDITIONS FOR SMOOTH WALLS',/,            &
 '   ---------------------------------------',/,            &
 '------------------------------------------------------------',/,&
 '                                         Minimum     Maximum',/,&
 '------------------------------------------------------------',/,&
 '   Rel velocity at the wall uiptn : ',2E12.5                 ,/,&
 '   Friction velocity        uet   : ',2E12.5                 ,/,&
 '   Friction velocity        uk    : ',2E12.5                 ,/,&
 '   Dimensionless distance   yplus : ',2E12.5                 ,/,&
 '   ------------------------------------------------------   ',/,&
 '   Nb of reversal of the velocity at the wall   : ',I10      ,/,&
 '   Nb of faces within the viscous sub-layer     : ',I10      ,/,&
 '   Total number of wall faces                   : ',I10      ,/,&
 '------------------------------------------------------------',  &
 /,/)

 2020 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: MESH NOT ENOUGH REFINED AT THE WALL            ',/,&
'@    ========                                                ',/,&
'@    The mesh does not seem to be enough refined at the wall ',/,&
'@      to be able to run a laminar simulation.               ',/,&
'@                                                            ',/,&
'@    The last time step at which too large values for the    ',/,&
'@      dimensionless distance to the wall (yplus) have been  ',/,&
'@      observed is the time step ',I10                        ,/,&
'@                                                            ',/,&
'@    The minimum value for yplus must be lower than the      ',/,&
'@      limit value YPLULI = ',E14.5                           ,/,&
'@                                                            ',/,&
'@    Have a look at the distribution of yplus at the wall    ',/,&
'@      (with EnSight for example) to conclude on the way     ',/,&
'@      the results quality might be affected.                ')

 2030 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: MESH NOT ENOUGH REFINED AT THE WALL            ',/,&
'@    ========                                                ',/,&
'@    The mesh does not seem to be enough refined at the wall ',/,&
'@      to be able to run a v2f simulation                    ',/,&
'@      (phi-fbar or BL-v2/k)                                 ',/,&
'@                                                            ',/,&
'@    The last time step at which too large values for the    ',/,&
'@      dimensionless distance to the wall (yplus) have been  ',/,&
'@      observed is the time step ',I10                        ,/,&
'@                                                            ',/,&
'@    The minimum value for yplus must be lower than the      ',/,&
'@      limit value YPLULI = ',E14.5                           ,/,&
'@                                                            ',/,&
'@    Have a look at the distribution of yplus at the wall    ',/,&
'@      (with EnSight for example) to conclude on the way     ',/,&
'@      the results quality might be affected.                ')

 2040 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: MESH TOO REFINED AT THE WALL                   ',/,&
'@    ========                                                ',/,&
'@    The mesh seems to be too refined at the wall to use     ',/,&
'@      a high-Reynolds turbulence model.                     ',/,&
'@                                                            ',/,&
'@    The last time step at which too small values for the    ',/,&
'@      dimensionless distance to the wall (yplus) have been  ',/,&
'@      observed is the time step ',I10                        ,/,&
'@                                                            ',/,&
'@    The minimum value for yplus must be greater than the    ',/,&
'@      limit value YPLULI = ',E14.5                           ,/,&
'@                                                            ',/,&
'@    Have a look at the distribution of yplus at the wall    ',/,&
'@      (with EnSight for example) to conclude on the way     ',/,&
'@      the results quality might be affected.                ')
 2050 format(                                                           &
'@                                                            ',/,&
'@    This warning is only printed at the first two           ',/,&
'@      occurences of the problem and at the last time step   ',/,&
'@      of the calculation. The vanishing of the message does ',/,&
'@      not necessarily mean the vanishing of the problem.    ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2060 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

#endif

!----
! FIN
!----

return
end subroutine
