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

subroutine clptrg &
!================

 ( idbia0 , idbra0 ,                                              &
   nvar   , nscal  , nphas  ,                                     &
   iphas  , isvhb  ,                                              &
   icodcl ,                                                       &
   ia     ,                                                       &
   dt     , rtp    , rtpa   , propce , propfa , propfb , rcodcl , &
   coefu  , rijipb , coefa  , coefb  , uetbor , visvdr ,          &
   hbord  , thbord ,                                              &
   w1     , w2     , w3     , w4     , w5     , w6     ,          &
   ra     )

!===============================================================================
! FONCTION :
! --------

! CONDITIONS LIMITES EN PAROI RUGUEUSE POUR TOUTES LES VARIABLES
!  DE LA PHASE IPHAS

! ON SUPPOSE QUE ICODCL(IU) = 6 =>
!             PAROI RUGUEUSE POUR TOUTES LES VARIABLES TURBULENTES
!  (A PRIORI PEU RESTRICTIF EN MONOPHASIQUE)

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! nphas            ! i  ! <-- ! number of phases                               !
! iphas            ! i  ! <-- ! phase number                                   !
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
! ia(*)            ! ia ! --- ! main integer work array                        !
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
! uetbor           ! tr ! --> ! vitesse de frottement au bord                  !
! (nfabor,nphas    !    !     !  pour van driest en les                        !
! visvdr(nphas)    ! tr ! <-- ! viscosite dynamique ds les cellules            !
! (ncelet,nphas    !    !     !  de bord apres amortisst de v driest           !
! hbord            ! tr ! --> ! coefficients d'echange aux bords               !
! (nfabor)         !    !     !                                                !
! thbord           ! tr ! <-- ! temperature aux bords en i'                    !
! (nfabor)         !    !     !    (plus exactmt : var. energetique)           !
! w1,2,3,4,5,6     ! ra ! --- ! work arrays                                    !
!  (ncelet)        !    !     !  (computation of pressure gradient)            !
! ra(*)            ! ra ! --- ! main real work array                           !
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

!===============================================================================

implicit none

! Arguments

integer          idbia0 , idbra0
integer          nvar   , nscal  , nphas
integer          iphas  , isvhb

integer          icodcl(nfabor,nvar)
integer          ia(*)

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision rcodcl(nfabor,nvar,3)
double precision coefu(nfabor,ndim), rijipb(nfabor,6)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision uetbor(nfabor,nphas), visvdr(ncelet,nphas)
double precision hbord(nfabor),thbord(nfabor)
double precision w1(ncelet),w2(ncelet),w3(ncelet)
double precision w4(ncelet),w5(ncelet),w6(ncelet)
double precision ra(*)

! Local variables

integer          idebia, idebra
integer          ifac, iel, ivar, isou, ii, jj, kk, ll, isvhbl
integer          ihcp, iscal
integer          modntl
integer          iuiptn
integer          iuiph , iviph , iwiph
integer          ikiph , iepiph, iphiph, ifbiph, iomgip
integer          ir11ip, ir22ip, ir33ip, ir12ip, ir13ip, ir23ip
integer          iclu  , iclv  , iclw  , iclk  , iclep
integer          icl11 , icl22 , icl33 , icl12 , icl13 , icl23
integer          icluf , iclvf , iclwf , iclphi, iclfb , iclomg
integer          ipcrom, ipcvis, ipcvst, ipccp , ipccv
integer          iclvar, ipcvsl, iclvaf
integer          iyplbp
double precision rnx, rny, rnz, rxnn
double precision tx, ty, tz, txn, txn0, t2x, t2y, t2z
double precision utau, upx, upy, upz, usn
double precision uiptn, uiptnf, uiptmn, uiptmx
double precision uetmax, uetmin, ukmax, ukmin, yplumx, yplumn
double precision uk, uet, tet, yplus
double precision gredu, rib, lmo, q0, e0
double precision cfnnu, cfnns, cfnnk, cfnne
double precision sqrcmu, clsyme, ek
double precision xmutlm
double precision rcprod, rcflux
double precision cpp, rkl,  prdtl
double precision hflui, hredui, hint, hext, pimp
double precision und0, deuxd0
double precision eloglo(3,3), alpha(6,6)
double precision rcodcx, rcodcy, rcodcz, rcodsn
double precision visclc, visctc, romc  , distbf, srfbnf, cpscv
double precision distbf0,rugd,rugt,ydep,act

!===============================================================================

!===============================================================================
! 1.  INITIALISATIONS
!===============================================================================

! Initialize variables to avoid compiler warnings

iepiph = 0
ifbiph = 0
ikiph = 0
iomgip = 0
iphiph = 0
ir11ip = 0
ir22ip = 0
ir33ip = 0
ir12ip = 0
ir13ip = 0
ir23ip = 0
iclep = 0
iclk = 0
iclomg = 0
iclfb = 0
iclphi = 0
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
idebia = idbia0
idebra = idbra0

! --- Constantes
uet = 1.d0
utau = 1.d0
sqrcmu = sqrt(cmu)

! --- Constantes de relaxation en atmosphere non neutre
cfnnu=1.d0
cfnns=1.d0
cfnnk=1.d0
cfnne=1.d0

und0   = 1.d0
deuxd0 = 2.d0

! --- Variables
iuiph  = iu  (iphas)
iviph  = iv  (iphas)
iwiph  = iw  (iphas)
if(itytur(iphas).eq.2) then
  ikiph  = ik  (iphas)
  iepiph = iep (iphas)
elseif(itytur(iphas).eq.3) then
  ir11ip = ir11(iphas)
  ir22ip = ir22(iphas)
  ir33ip = ir33(iphas)
  ir12ip = ir12(iphas)
  ir13ip = ir13(iphas)
  ir23ip = ir23(iphas)
  iepiph = iep (iphas)
elseif(iturb(iphas).eq.50) then
  ikiph  = ik  (iphas)
  iepiph = iep (iphas)
  iphiph = iphi(iphas)
  ifbiph = ifb (iphas)
elseif(iturb(iphas).eq.60) then
  ikiph  = ik  (iphas)
  iomgip = iomg(iphas)
endif

! --- Conditions aux limites
iclu   = iclrtp(iuiph ,icoef)
iclv   = iclrtp(iviph ,icoef)
iclw   = iclrtp(iwiph ,icoef)
if(itytur(iphas).eq.2) then
  iclk   = iclrtp(ikiph ,icoef)
  iclep  = iclrtp(iepiph,icoef)
elseif(itytur(iphas).eq.3) then
  icl11  = iclrtp(ir11ip,icoef)
  icl22  = iclrtp(ir22ip,icoef)
  icl33  = iclrtp(ir33ip,icoef)
  icl12  = iclrtp(ir12ip,icoef)
  icl13  = iclrtp(ir13ip,icoef)
  icl23  = iclrtp(ir23ip,icoef)
  iclep  = iclrtp(iepiph,icoef)
elseif(iturb(iphas).eq.50) then
  iclk   = iclrtp(ikiph ,icoef)
  iclep  = iclrtp(iepiph,icoef)
  iclphi = iclrtp(iphiph,icoef)
  iclfb  = iclrtp(ifbiph,icoef)
elseif(iturb(iphas).eq.60) then
  iclk   = iclrtp(ikiph ,icoef)
  iclomg = iclrtp(iomgip,icoef)
endif

icluf  = iclrtp(iuiph ,icoeff)
iclvf  = iclrtp(iviph ,icoeff)
iclwf  = iclrtp(iwiph ,icoeff)

! --- Grandeurs physiques
ipcrom = ipproc(irom  (iphas))
ipcvis = ipproc(iviscl(iphas))
ipcvst = ipproc(ivisct(iphas))
if(icp(iphas).gt.0) then
  ipccp  = ipproc(icp   (iphas))
else
  ipccp = 0
endif

! --- Compressible

if ( ippmod(icompf) .ge. 0 ) then
  if(icv(iphas).gt.0) then
    ipccv  = ipproc(icv   (iphas))
  else
    ipccv = 0
  endif
endif

! --- Post traitement de Yplus
iyplbp = iyplbr+(iphas-1)*nfabor


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

! COMPTEUR RETOURNEMENT

iuiptn = 0


!     En v2f on met directement u=0 donc UIPTMX et UIPTMN vaudront
!     forcement 0
if (iturb(iphas).eq.50) then
  uiptmx = 0.d0
  uiptmn = 0.d0
endif

! --- Boucle sur les faces : debut
do ifac = 1, nfabor

! --- Test sur la presence d'une condition de paroi rugueuse
  if( icodcl(ifac,iuiph).eq.6 ) then

    iel = ifabor(ifac)

! --- Proprietes physiques
    visclc = propce(iel,ipcvis)
    visctc = propce(iel,ipcvst)
    romc   = propce(iel,ipcrom)

! --- Grandeurs geometriques
    distbf = ra(idistb-1+ifac)
    srfbnf = surfbn(ifac)

!===============================================================================
! 1. REPERE LOCAL
!===============================================================================


! ---> NORMALE UNITAIRE

    rnx = surfbo(1,ifac)/srfbnf
    rny = surfbo(2,ifac)/srfbnf
    rnz = surfbo(3,ifac)/srfbnf

! ---> PRISE EN COMPTE DE LA VITESSE DE DEFILEMENT

    rcodcx = rcodcl(ifac,iuiph,1)
    rcodcy = rcodcl(ifac,iviph,1)
    rcodcz = rcodcl(ifac,iwiph,1)

!     Si on n'est pas en ALE, on force la vitesse de deplacement
!       de la face a etre tangentielle (et on met a jour rcodcl
!       pour une utilisation eventuelle)
    if (iale.eq.0.and.imobil.eq.0) then
      rcodsn = rcodcx*rnx+rcodcy*rny+rcodcz*rnz
      rcodcx = rcodcx -rcodsn*rnx
      rcodcy = rcodcy -rcodsn*rny
      rcodcz = rcodcz -rcodsn*rnz
      rcodcl(ifac,iuiph,1) = rcodcx
      rcodcl(ifac,iviph,1) = rcodcy
      rcodcl(ifac,iwiph,1) = rcodcz
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

    elseif(itytur(iphas).eq.3) then

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

    if (itytur(iphas).eq.3) then

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
!==========================================================================

   if (abs(utau).le.epzero) utau = epzero

! RUGD : rugosite de paroi pour les variables dynamiques
!        seule la valeur stockee pour IU est utilisee
   rugd=rcodcl(ifac,iuiph,3)
! NB: en rugueux la signification de yplus change           !
   yplus=distbf/rugd

! PSEUDO DECALAGE DE LA PAROI de RUGD ((DISTBF+RUGD)/RUGD):
   uet = utau/log(yplus+1.d0)*xkappa

    if (ideuch(iphas).eq.0) then
      uk = uet
    else
! Si IDEUCH=1 ou 2 on calcule uk et uet

      if(itytur(iphas).eq.2 .or. iturb(iphas).eq.50               &
           .or. iturb(iphas).eq.60) then
        ek = rtp(iel,ikiph)
      elseif(itytur(iphas).eq.3) then
        ek = 0.5d0*                                               &
               (rtp(iel,ir11ip)+rtp(iel,ir22ip)+rtp(iel,ir33ip))
      endif

      uk = cmu025*sqrt(ek)

    endif

  if (ippmod(iatmos).ge.1) then

    ! Compute reduced gravity for non horizontal walls :
    gredu = gx*rnx + gy*rny + gz*rnz

    call atmcls                                                   &
    !==========
 ( idebia , idebra ,                                              &
   nvar   , nscal  , nphas  ,                                     &
   iphas  , ifac   , iel    ,                                     &
   uk     , utau   , yplus  ,                                     &
   uet    ,                                                       &
   gredu  , q0     , e0     , rib    ,lmo     ,                   &
   cfnnu  , cfnns  , cfnnk  , cfnne  ,                            &
   icodcl ,                                                       &
   ia     ,                                                       &
   dt     , rtp    ,          propce , propfa , propfb , rcodcl , &
   w1     , w2     , w3     , w4     , w5     , w6     ,          &
   ra     )

  endif

    if(ideuch(iphas).eq.0) then
      uk = uet
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

    if(itytur(iphas).eq.4.and.idries(iphas).eq.1) then
      uetbor(ifac,iphas) = uet
      if (visvdr(iel,iphas).lt.-900.d0) then
        propce(iel,ipcvst) = propce(iel,ipcvst)
! NB amortissement de van Driest a revoir en rugueux :
!    &             *(1.D0-EXP(-YPLUS/CDRIES(IPHAS)))**2
        visvdr(iel,iphas) = propce(iel,ipcvst)
        visctc = propce(iel,ipcvst)
      endif
    endif


! Sauvegarde de yplus si post traite

    if(mod(ipstdv,ipstyp).eq.0) then
      ra(iyplbp+ifac-1) = yplus
    endif


!===============================================================================
! 3. CONDITIONS AUX LIMITES SUR LA VITESSE
!===============================================================================

!              UIPTN  respecte la production de k
!              de facon conditionnelle    --> Coef RCPROD
!              UIPTNF respecte le flux
!               de facon conditionnelle   --> Coef RCFLUX

! On traite tous les modeles de la meme maniere (a revisiter)
    if (itytur(iphas).eq.2 .or. iturb(iphas).eq.60                &
        .or.iturb(iphas).eq.0 .or.iturb(iphas).eq.10.or.          &
           itytur(iphas).eq.3.or.itytur(iphas).eq.4) then

! Pseudo decalage de la paroi de la distance RUGD :
! modified for non neutral boundary layer (cfnnu)
      distbf0=distbf+rugd
      xmutlm = xkappa*uk*distbf0*romc

      rcprod = max(distbf/distbf0,                                &
              deuxd0*distbf*sqrt(xmutlm/visctc/distbf0**2)        &
              -1./(2.+rugd/distbf0))

      rcflux = max(xmutlm,visctc)/(visclc+visctc)*distbf/distbf0

      uiptn  = min(utau,max(utau - uk/xkappa*rcprod*cfnnu,0.d0))
      uiptnf = utau - uet/xkappa*rcflux*cfnnu

!     Clipping :
!     On borne U_f,grad entre 0 et Utau (il y a surement mieux...)
!     - 0    : on interdit le retournement en face de bord, qui est en
!              contradiction avec l'hypothèse de loi log.
!     - Utau : la production turbulente ne peut etre nulle
!     On empeche U_f,flux d'etre negatif

      if (uiptnf.lt.epzero) uiptnf = 0.d0

      uiptmx = max(uiptn,uiptmx)
      uiptmn = min(uiptn,uiptmn)

      if(uiptn.lt.-epzero) iuiptn = iuiptn + 1

      uetmax = max(uet,uetmax)
      uetmin = min(uet,uetmin)
      ukmax = max(uk,ukmax)
      ukmin = min(uk,ukmin)

      coefa(ifac,iclu)   = uiptn *tx *txn0
      coefa(ifac,iclv)   = uiptn *ty *txn0
      coefa(ifac,iclw)   = uiptn *tz *txn0
      coefa(ifac,icluf)  = uiptnf*tx *txn0
      coefa(ifac,iclvf)  = uiptnf*ty *txn0
      coefa(ifac,iclwf)  = uiptnf*tz *txn0

      coefa(ifac,iclu)   = coefa(ifac,iclu)  +rcodcx
      coefa(ifac,iclv)   = coefa(ifac,iclv)  +rcodcy
      coefa(ifac,iclw)   = coefa(ifac,iclw)  +rcodcz
      coefa(ifac,icluf)  = coefa(ifac,icluf) +rcodcx
      coefa(ifac,iclvf)  = coefa(ifac,iclvf) +rcodcy
      coefa(ifac,iclwf)  = coefa(ifac,iclwf) +rcodcz

      coefb(ifac,iclu)   = 0.d0
      coefb(ifac,iclv)   = 0.d0
      coefb(ifac,iclw)   = 0.d0
      coefb(ifac,icluf)  = 0.d0
      coefb(ifac,iclvf)  = 0.d0
      coefb(ifac,iclwf)  = 0.d0


    elseif(iturb(iphas).eq.50) then

!     Avec ces conditions, pas besoin de calculer UIPTMX, UIPTMN
!     et IUIPTN qui sont nuls (valeur d'initialisation)

      coefa(ifac,iclu)   = rcodcx
      coefa(ifac,iclv)   = rcodcy
      coefa(ifac,iclw)   = rcodcz

      coefb(ifac,iclu)   = 0.d0
      coefb(ifac,iclv)   = 0.d0
      coefb(ifac,iclw)   = 0.d0

    endif

!===============================================================================
! 4. CONDITIONS AUX LIMITES SUR K ET EPSILON
!===============================================================================

    ydep = distbf*0.5d0+rugd
    if (itytur(iphas).eq.2) then

      coefa(ifac,iclk)   = uk**2/sqrcmu*cfnnk
      coefb(ifac,iclk)   = 0.d0

      coefa(ifac,iclep)  = uk**3/(xkappa*ydep**2)*distbf*cfnne
      coefb(ifac,iclep)  = 1.d0

!===============================================================================
! 5. CONDITIONS AUX LIMITES SUR RIJ ET EPSILON
!===============================================================================

    elseif (itytur(iphas).eq.3) then

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

        if (iclptr(iphas).eq.1) then
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
                            eloglo(jj,2)*eloglo(kk,1))*uet*uk*cfnnk

      enddo

! ---> SCALAIRE EPSILON
!      on traite comme en k-eps


      coefa(ifac,iclep)  = uk**3/(xkappa*ydep**2)*distbf*cfnne
      coefb(ifac,iclep) = 1.d0

!===============================================================================
! 6. CONDITIONS AUX LIMITES SUR K, EPSILON, F_BARRE ET PHI
!===============================================================================

    elseif (iturb(iphas).eq.50) then

      coefa(ifac,iclk) = 0.d0
      coefb(ifac,iclk) = 0.d0
      coefa(ifac,iclep) =                                         &
           2.0d0*propce(iel,ipcvis)/propce(iel,ipcrom)            &
           *rtp(iel,ikiph)/distbf**2
      coefb(ifac,iclep) = 0.d0
      coefa(ifac,iclphi) = 0.0d0
      coefb(ifac,iclphi) = 0.0d0
      coefa(ifac,iclfb) = 0.0d0
      coefb(ifac,iclfb) = 0.0d0
!===============================================================================
! 7. CONDITIONS AUX LIMITES SUR K ET OMEGA
!===============================================================================

    elseif (iturb(iphas).eq.60) then

!     on est toujours hors de la sous couche visqeuse
        coefa(ifac,iclk)   = uk**2/sqrcmu
        coefb(ifac,iclk)   = 0.d0
        coefa(ifac,iclomg) = distbf*4.d0*uk**3*romc**2/           &
              (sqrcmu*xkappa*visclc**2*yplus**2)
        coefb(ifac,iclomg) = 1.d0

    endif

!===============================================================================
! 8. CONDITIONS AUX LIMITES SUR LES SCALAIRES
!              (AUTRES QUE PRESSION, K, EPSILON, RIJ, VARIANCES)
!    Pour les variances, pas de traitement specifique en paroi : voir
!      condli.
!===============================================================================

    if(nscal.ge.1) then

      do ll = 1, nscal

        if(iphsca(ll).eq.iphas.and.iscavr(ll).le.0) then

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
            cpp = cp0(iphas)
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
!               defini de la meme façon que l'on resolve
!               en enthalpie ou en energie, soit Mu*Cp/Lambda.
!               Si l'on resout en energie, on a calcule ci-dessus
!               Mu*Cv/Lambda.

          if ( ippmod(icompf).ge.0 ) then
            if(iscsth(iscal).eq.3) then
              if(ipccp.gt.0) then
                prdtl = prdtl*propce(iel,ipccp )
              else
                prdtl = prdtl*cp0(iphas)
              endif
              if(ipccv.gt.0) then
                prdtl = prdtl/propce(iel,ipccv )
              else
                prdtl = prdtl/cv0(iphas)
              endif
            endif
          endif

!          CAS TURBULENT
          if (iturb(iphas).ne.0) then
            if ( ippmod(icompf) .ge. 0 ) then
!                 En compressible, pour l'energie LAMBDA/CV+CP/CV*(MUT/SIGMAS)
              if(ipccp.gt.0) then
                cpscv = propce(iel,ipproc(icp(iphas)))
              else
                cpscv = cp0(iphas)
              endif
              if(ipccv.gt.0) then
                cpscv = cpscv/propce(iel,ipproc(icv(iphas)))
              else
                cpscv = cpscv/cv0(iphas)
              endif
              hint = hint*(rkl+cpscv*visctc/sigmas(ll))/distbf
            else
              hint = hint*(rkl+visctc/sigmas(ll))/distbf
            endif
!          CAS LAMINAIRE
          else
            hint  = hint*rkl/distbf
          endif

          if(iturb(iphas).ne.0.and.icodcl(ifac,ivar).eq.6)then

! Loi rugueuse, on recalcule le coefficient d'echange fluide - paroi
            rugt=rcodcl(ifac,iviph,3)
            act = xkappa/log((distbf+rugt)/rugt)
            hflui = romc*cpp*uet*act*cfnns
          else
            hflui = hint
          endif

          if (isvhbl .gt. 0) hbord(ifac) = hflui



! --->  C.L DE TYPE DIRICHLET AVEC OU SANS COEFFICIENT D'ECHANGE

!     Si on a deux types de conditions aux limites (ICLVAR, ICLVAF)
!       il faut que le flux soit traduit par ICLVAF.
!     Si on n'a qu'un type de condition, peu importe (ICLVAF=ICLVAR)
!     Pour le moment, dans cette version compressible, on impose un
!       flux nul pour ICLVAR, lorsqu'il est différent de ICLVAF (cette
!       condition ne sert qu'à la reconstruction des gradients et
!       s'applique à l'energie totale qui inclut l'energie cinétique :


          if( icodcl(ifac,ivar).eq.6 ) then
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
!         ABS(ISCSTH(II)).EQ.1 : RA(IHCONV-1+IFAC+NFABOR*(IPH-1)) = HINT
!         puisque HINT = VISLS * CP / DISTBR
!                      = lambda/distance en W/(m2 K)

!       lorsque la variable transportee est l'enthalpie
!         ISCSTH(II).EQ.2 : RA(IHCONV-1+IFAC+NFABOR*(IPH-1)) = HINT*CPR
!         avec
!            IF(IPCCP.GT.0) THEN
!              CPR = PROPCE(IEL,IPCCP )
!            ELSE
!              CPR = CP0(IPHAS)
!            ENDIF
!         puisque HINT = VISLS / DISTBR
!                      = lambda/(CP * distance)

!       lorsque la variable transportee est l'energie (compressible)
!         ISCSTH(II).EQ.3 :
!         on procede comme pour l'enthalpie avec CV au lieu de CP
!         (rq : il n'y a pas d'hypothèse, sf en non orthogonal :
!               le flux est le bon et le coef d'echange aussi)

!      De meme dans condli.



!               Si on rayonne sur la phase et que
!                  le scalaire est la variable energetique

            if (iirayo.ge.1         .and.                         &
                ll.eq.iscalt(iphas) .and. iphas.eq.irapha   ) then

!                On calcule le coefficient d'echange en W/(m2 K)

!                Si on resout en enthalpie
              if(iscsth(ll).eq.2) then
!                  Si Cp variable
                if(ipccp.gt.0) then
                  propfb(ifac,ipprob(ihconv)) = hflui*propce(iel,ipccp )
                else
                  propfb(ifac,ipprob(ihconv)) = hflui*cp0(iphas)
                endif

!                  Si on resout en energie (compressible)
              elseif(iscsth(ll).eq.3) then
!                    Si Cv variable
                if(ipccv.gt.0) then
                  propfb(ifac,ipprob(ihconv)) = hflui*propce(iel,ipccv )
                else
                  propfb(ifac,ipprob(ihconv)) = hflui*cv0(iphas)
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

if(iwarni(iuiph).ge.0) then
  if(ntlist.gt.0) then
    modntl = mod(ntcabs,ntlist)
  elseif(ntlist.eq.-1.and.ntcabs.eq.ntmabs) then
    modntl = 0
  else
    modntl = 1
  endif

  if (modntl.eq.0 .or. iwarni(iuiph).ge.2) then
    write(nfecra,2010) iphas,                                     &
         uiptmn,uiptmx,uetmin,uetmax,ukmin,ukmax,yplumn,yplumx,   &
         iuiptn
  endif

endif

!===============================================================================
! 10.  FORMATS
!===============================================================================

#if defined(_CS_LANG_FR)

 1000 format(/,' LA NORMALE A LA FACE DE BORD DE PAROI ',I10,/,   &
         ' EST NULLE ; COORDONNEES : ',3E12.5)

 2010 format(/,                                                   &
 3X,'** CONDITIONS AUX LIMITES EN PAROI RUGUEUSE',/,        &
 '   -------------------------------------------',/,        &
 '------------------------------------------------------------',/,&
 ' Phase ',I6,'                            Minimum     Maximum',/,&
 '------------------------------------------------------------',/,&
 '   Vitesse rel. en paroi    uiptn : ',2E12.5                 ,/,&
 '   Vitesse de frottement    uet   : ',2E12.5                 ,/,&
 '   Vitesse de frottement    uk    : ',2E12.5                 ,/,&
 '   Distance adim. rugueuse  yplus : ',2E12.5                 ,/,&
 '   ------------------------------------------------------   ',/,&
 '   Nbre de retournements de la vitesse en paroi : ',I10      ,/,&
 '------------------------------------------------------------',  &
 /,/)


#else

 1000 format(/,' THE NORMAL TO THE WALL BOUNDARY FACE ',I10,/,    &
         ' IS NULL; COORDINATES: ',3E12.5)

 2010 format(/,                                                   &
 3X,'** BOUNDARY CONDITIONS FOR ROUGH WALLS',/,             &
 '   --------------------------------------',/,             &
 '------------------------------------------------------------',/,&
 ' Phase ',I6,'                            Minimum     Maximum',/,&
 '------------------------------------------------------------',/,&
 '   Rel velocity at the wall uiptn : ',2E12.5                 ,/,&
 '   Friction velocity        uet   : ',2E12.5                 ,/,&
 '   Friction velocity        uk    : ',2E12.5                 ,/,&
 '   Rough dimensionless dist yplus : ',2E12.5                 ,/,&
 '   ------------------------------------------------------   ',/,&
 '   Nb of reversal of the velocity at the wall   : ',I10      ,/,&
 '------------------------------------------------------------',  &
 /,/)


#endif

!----
! FIN
!----

return
end subroutine
