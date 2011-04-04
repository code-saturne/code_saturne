!-------------------------------------------------------------------------------

!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2010 EDF S.A., France

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

subroutine resrij &
!================

 ( idbia0 , idbra0 ,                                              &
   nvar   , nscal  , nphas  , ncepdp , ncesmp ,                   &
   iphas  , ivar   , isou   , ipp    ,                            &
   icepdc , icetsm , itpsmp ,                                     &
   ia     ,                                                       &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  , produc , grarox , graroy , graroz ,          &
   ckupdc , smcelp , gamma  ,                                     &
   viscf  , viscb  , coefax ,                                     &
   tslage , tslagi ,                                              &
   dam    , xam    , drtp   , smbr   , rovsdt ,                   &
   w1     , w2     , w3     , w4     ,                            &
   w5     , w6     , w7     , w8     , w9     ,                   &
   ra     )

!===============================================================================
! FONCTION :
! ----------

! RESOLUTION DES EQUATIONS CONVECTION DIFFUSION TERME SOURCE
!   POUR Rij (modele standard LRR)
! VAR  = R11 R22 R33 R12 R13 R23
! ISOU =  1   2   3   4   5   6

!-------------------------------------------------------------------------------
!ARGU                             ARGUMENTS
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! nphas            ! i  ! <-- ! number of phases                               !
! ncepdp           ! i  ! <-- ! number of cells with head loss                 !
! ncesmp           ! i  ! <-- ! number of cells with mass source term          !
! iphas            ! i  ! <-- ! phase number                                   !
! ivar             ! i  ! <-- ! variable number                                !
! isou             ! e  ! <-- ! numero de passage                              !
! ipp              ! e  ! <-- ! numero de variable pour sorties post           !
! icepdc(ncelet    ! te ! <-- ! numero des ncepdp cellules avec pdc            !
! icetsm(ncesmp    ! te ! <-- ! numero des cellules a source de masse          !
! itpsmp           ! te ! <-- ! type de source de masse pour la                !
! (ncesmp)         !    !     !  variables (cf. ustsma)                        !
! ia(*)            ! ia ! --- ! main integer work array                        !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
! produc           ! tr ! <-- ! tableau de travail pour production             !
!  (6,ncelet)      !    !     ! (sans rho volume)                              !
! grarox,y,z       ! tr ! <-- ! tableau de travail pour grad rom               !
!  (ncelet)        !    !     !                                                !
! ckupdc           ! tr ! <-- ! tableau de travail pour pdc                    !
!  (ncepdp,6)      !    !     !                                                !
! smcelp(ncesmp    ! tr ! <-- ! valeur de la variable associee a la            !
!                  !    !     !  source de masse                               !
! gamma(ncesmp)    ! tr ! <-- ! valeur du flux de masse                        !
! viscf(nfac)      ! tr ! --- ! visc*surface/dist aux faces internes           !
! viscb(nfabor     ! tr ! --- ! visc*surface/dist aux faces de bord            !
! coefax(nfabor    ! tr ! --- ! tab de trav pour cond.lim. paroi               !
!                  ! tr ! --- !   attention : uniquement avec echo             !
!                  ! tr ! --- !   de paroi et abs(icdpar) = 1                  !
! tslage(ncelet    ! tr ! <-- ! ts explicite couplage retour lagr.             !
! tslagi(ncelet    ! tr ! <-- ! ts implicite couplage retour lagr.             !
! dam(ncelet       ! tr ! --- ! tableau de travail pour matrice                !
! xam(nfac,*)      ! tr ! --- ! tableau de travail pour matrice                !
! drtp(ncelet      ! tr ! --- ! tableau de travail pour increment              !
! smbr(ncelet      ! tr ! --- ! tableau de travail pour sec mem                !
! rovsdt(ncelet    ! tr ! --- ! tableau de travail pour terme instat           !
! w1...9(ncelet    ! tr ! --- ! tableau de travail                             !
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
use dimens, only: ndimfb
use numvar
use entsor
use optcal
use cstphy
use cstnum
use pointe
use parall
use period
use lagpar
use lagran
use mesh

!===============================================================================

implicit none

! Arguments

integer          idbia0 , idbra0
integer          nvar   , nscal  , nphas
integer          ncepdp , ncesmp
integer          iphas  , ivar   , isou   , ipp

integer          icepdc(ncepdp)
integer          icetsm(ncesmp), itpsmp(ncesmp)
integer          ia(*)

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(ndimfb,*)
double precision coefa(ndimfb,*), coefb(ndimfb,*)
double precision produc(6,ncelet)
double precision grarox(ncelet), graroy(ncelet), graroz(ncelet)
double precision ckupdc(ncepdp,6)
double precision smcelp(ncesmp), gamma(ncesmp)
double precision viscf(nfac), viscb(nfabor), coefax(nfabor)
double precision tslage(ncelet),tslagi(ncelet)
double precision dam(ncelet), xam(nfac,2)
double precision drtp(ncelet), smbr(ncelet), rovsdt(ncelet)
double precision w1(ncelet), w2(ncelet), w3(ncelet)
double precision w4(ncelet), w5(ncelet), w6(ncelet)
double precision w7(ncelet), w8(ncelet), w9(ncelet)
double precision ra(*)

! Local variables

integer          idebia, idebra, ifinia
integer          init  , ifac  , iel   , inc   , iccocg
integer          ii    , jj    , iiun
integer          ir11ip, ir22ip, ir33ip, ir12ip, ir13ip, ir23ip
integer          ieiph , iuiph
integer          ipcrom, ipcvis, iflmas, iflmab, ipcroo
integer          iclvar, iclvaf
integer          nswrgp, imligp, iwarnp, iphydp
integer          iconvp, idiffp, ndircp, ireslp
integer          nitmap, nswrsp, ircflp, ischcp, isstpp, iescap
integer          imgrp , ncymxp, nitmfp
integer          iptsta
integer          isoluc
integer          maxelt, ils
double precision blencp, epsilp, epsrgp, climgp, extrap, relaxp
double precision epsrsp
double precision trprod, trrij , cstrij, rctse , deltij
double precision grdpx , grdpy , grdpz , grdsn
double precision surfn2
double precision tuexpr, thets , thetv , thetp1
double precision d1s3  , d2s3


!===============================================================================

!===============================================================================
! 1. INITIALISATION
!===============================================================================

idebia = idbia0
idebra = idbra0

if(iwarni(ivar).ge.1) then
  write(nfecra,1000) nomvar(ipp)
endif

iuiph  = iu  (iphas)
ir11ip = ir11(iphas)
ir22ip = ir22(iphas)
ir33ip = ir33(iphas)
ir12ip = ir12(iphas)
ir13ip = ir13(iphas)
ir23ip = ir23(iphas)
ieiph  = iep (iphas)

ipcrom = ipproc(irom  (iphas))
ipcvis = ipproc(iviscl(iphas))
iflmas = ipprof(ifluma(iuiph))
iflmab = ipprob(ifluma(iuiph))

iclvar = iclrtp(ivar,icoef)
iclvaf = iclrtp(ivar,icoeff)

deltij = 1.0d0
if(isou.gt.3) then
  deltij = 0.0d0
endif
d1s3 = 1.d0/3.d0
d2s3 = 2.d0/3.d0

!     S pour Source, V pour Variable
thets  = thetst(iphas)
thetv  = thetav(ivar )

ipcroo = ipcrom
if(isto2t(iphas).gt.0.and.iroext(iphas).gt.0) then
  ipcroo = ipproc(iroma(iphas))
endif
iptsta = 0
if(isto2t(iphas).gt.0) then
  iptsta = ipproc(itstua(iphas))
endif

do iel = 1, ncel
  smbr(iel) = 0.d0
enddo
do iel = 1, ncel
  rovsdt(iel) = 0.d0
enddo

!===============================================================================
! 2. TERMES SOURCES  UTILISATEURS
!===============================================================================
!(le premier argument PRODUC est lu en GRDVIT dans ustsri, mais ce
! tableau n'est dimensionne et utilise qu'en modele Rij SSG)

maxelt = max(ncelet, nfac, nfabor)
ils    = idebia
ifinia = ils + maxelt
call iasize('resrij',ifinia)

call ustsri                                                       &
!==========
 ( ifinia , idebra ,                                              &
   nvar   , nscal  , nphas  , ncepdp , ncesmp ,                   &
   iphas  , ivar   ,                                              &
   maxelt , ia(ils),                                              &
   icepdc , icetsm , itpsmp ,                                     &
   ia     ,                                                       &
   dt     , rtpa   , propce , propfa , propfb ,                   &
   coefa  , coefb  , ckupdc , smcelp , gamma  , produc , produc , &
   smbr   , rovsdt ,                                              &
!        ------   ------
   viscf  , viscb  , xam    ,                                     &
   w1     , w2     , w3     , w4     , w5     , w6     ,          &
   w7     , w8     , w9     , dam    , drtp   ,                   &
   ra     )

!     Si on extrapole les T.S.
if(isto2t(iphas).gt.0) then
  do iel = 1, ncel
!       Sauvegarde pour echange
    tuexpr = propce(iel,iptsta+isou-1)
!       Pour la suite et le pas de temps suivant
    propce(iel,iptsta+isou-1) = smbr(iel)
!       Second membre du pas de temps precedent
!       On suppose -ROVSDT > 0 : on implicite
!          le terme source utilisateur (le reste)
    smbr(iel) = rovsdt(iel)*rtpa(iel,ivar)  - thets*tuexpr
!       Diagonale
    rovsdt(iel) = - thetv*rovsdt(iel)
  enddo
else
  do iel = 1, ncel
    smbr(iel)   = rovsdt(iel)*rtpa(iel,ivar) + smbr(iel)
    rovsdt(iel) = max(-rovsdt(iel),zero)
  enddo
endif

!===============================================================================
! 2. TERMES SOURCES  LAGRANGIEN : COUPLAGE RETOUR
!===============================================================================

!     Ordre 2 non pris en compte
 if (iilagr.eq.2 .and. ltsdyn.eq.1 .and. iphas.eq.ilphas) then
   do iel = 1,ncel
     smbr(iel)   = smbr(iel)   + tslage(iel)
     rovsdt(iel) = rovsdt(iel) + max(-tslagi(iel),zero)
   enddo
 endif

!===============================================================================
! 3. TERME SOURCE DE MASSE
!===============================================================================


if (ncesmp.gt.0) then

!       Entier egal a 1 (pour navsto : nb de sur-iter)
  iiun = 1

!       On incremente SMBR par -Gamma RTPA et ROVSDT par Gamma (*theta)
  call catsma                                                     &
  !==========
 ( ncelet , ncel   , ncesmp , iiun   , isto2t(iphas) , thetv  ,   &
   icetsm , itpsmp ,                                              &
   volume , rtpa(1,ivar) , smcelp , gamma  ,                      &
   smbr   ,  rovsdt , w1 )

!       Si on extrapole les TS on met Gamma Pinj dans PROPCE
  if(isto2t(iphas).gt.0) then
    do iel = 1, ncel
      propce(iel,iptsta+isou-1) =                                 &
      propce(iel,iptsta+isou-1) + w1(iel)
    enddo
!       Sinon on le met directement dans SMBR
  else
    do iel = 1, ncel
      smbr(iel) = smbr(iel) + w1(iel)
    enddo
  endif

endif

!===============================================================================
! 4. TERME D'ACCUMULATION DE MASSE -(dRO/dt)*VOLUME
!    ET TERME INSTATIONNAIRE
!===============================================================================

! ---> Calcul de mij

init = 1
call divmas(ncelet,ncel,nfac,nfabor,init,nfecra,                  &
               ifacel,ifabor,propfa(1,iflmas),propfb(1,iflmab),w1)

! ---> Ajout au second membre

do iel = 1, ncel
  smbr(iel) = smbr(iel)                                           &
              + iconv(ivar)*w1(iel)*rtpa(iel,ivar)
enddo

! ---> Ajout dans la diagonale de la matrice
!     Extrapolation ou non, meme forme par coherence avec bilsc2

do iel=1,ncel
  rovsdt(iel) = rovsdt(iel)                                       &
            + istat(ivar)*(propce(iel,ipcrom)/dt(iel))*volume(iel)&
            - iconv(ivar)*w1(iel)*thetv
enddo


!===============================================================================
! 5. PRODUCTION, PHI1, PHI2, ET DISSIPATION
!===============================================================================


! ---> Calcul de k pour la suite du sous-programme
!       on utilise un tableau de travail puisqu'il y en a...
do iel = 1, ncel
  w8(iel) = 0.5d0 * (rtpa(iel,ir11ip) + rtpa(iel,ir22ip) +        &
                     rtpa(iel,ir33ip))
enddo

! ---> Terme source

!      (1-CRIJ2) Pij (pour toutes les composantes de Rij)

!      DELTAIJ*(2/3.CRIJ2.P+2/3.CRIJ1.EPSILON)
!                    (termes diagonaux pour R11, R22 et R33)

!      -DELTAIJ*2/3*EPSILON

!     Si on extrapole les TS
!       On modifie la partie implicite :
!         Dans PHI1, on ne prendra que RHO CRIJ1 EPSILON/K et non pas
!                                  RHO CRIJ1 EPSILON/K (1-2/3 DELTAIJ)
!         Cela permet de conserver k^n au lieu de (R11^(n+1)+R22^n+R33^n)
!         Ce choix est discutable. C'est la solution ISOLUC = 1
!       Si on veut tout prendre en implicite (comme c'est fait
!         en ordre 1 std), c'est la solution ISOLUC = 2
!       -> a tester plus precisement si necessaire


!     Si on extrapole les TS
if(isto2t(iphas).gt.0) then

  isoluc = 1

  do iel = 1, ncel

!     Demi-traces de Prod et R
    trprod = 0.5d0*(produc(1,iel)+produc(2,iel)+produc(3,iel))
    trrij  = w8(iel)

!     Calcul de Prod+Phi1+Phi2-Eps
!       = rhoPij-C1rho eps/k(Rij-2/3k dij)-C2rho(Pij-1/3Pkk dij)-2/3rho eps dij
!       Dans PROPCE :
!       = rhoPij-C1rho eps/k(   -2/3k dij)-C2rho(Pij-1/3Pkk dij)-2/3rho eps dij
!       = rho{2/3dij[C2 Pkk/2+(C1-1)eps)]+(1-C2)Pij           }
    propce(iel,iptsta+isou-1) = propce(iel,iptsta+isou-1)         &
                          + propce(iel,ipcroo) * volume(iel)      &
      *(   deltij*d2s3*                                           &
           (  crij2*trprod                                        &
            +(crij1-1.d0)* rtpa(iel,ieiph)  )                     &
         +(1.0d0-crij2)*produc(isou,iel)               )
!       Dans SMBR
!       =       -C1rho eps/k(Rij         )
!       = rho{                                     -C1eps/kRij}
    smbr(iel) = smbr(iel) + propce(iel,ipcrom) * volume(iel)      &
      *( -crij1*rtpa(iel,ieiph)/trrij * rtpa(iel,ivar)  )

!     Calcul de la partie implicite issue de Phi1
!       = C1rho eps/k(1        )
    rovsdt(iel) = rovsdt(iel) + propce(iel,ipcrom) * volume(iel)  &
                            *crij1*rtpa(iel,ieiph)/trrij*thetv

  enddo

!     Si on veut impliciter un bout de -C1rho eps/k(   -2/3k dij)
  if(isoluc.eq.2) then

    do iel = 1, ncel

      trrij  = w8(iel)

!    On enleve a PROPCE (avec IPCROO)
!       =       -C1rho eps/k(   -1/3Rij dij)
      propce(iel,iptsta+isou-1) = propce(iel,iptsta+isou-1)       &
                          - propce(iel,ipcroo) * volume(iel)      &
      *(deltij*d1s3*crij1*rtpa(iel,ieiph)/trrij * rtpa(iel,ivar))
!    On ajoute a SMBR (avec IPCROM)
!       =       -C1rho eps/k(   -1/3Rij dij)
      smbr(iel)                 = smbr(iel)                       &
                          + propce(iel,ipcrom) * volume(iel)      &
      *(deltij*d1s3*crij1*rtpa(iel,ieiph)/trrij * rtpa(iel,ivar))
!    On ajoute a ROVSDT (avec IPCROM)
!       =        C1rho eps/k(   -1/3    dij)
      rovsdt(iel) = rovsdt(iel) + propce(iel,ipcrom) * volume(iel)&
      *(deltij*d1s3*crij1*rtpa(iel,ieiph)/trrij                 )
    enddo

  endif

!     Si on n'extrapole pas les termes sources
else

  do iel = 1, ncel

!     Demi-traces de Prod et R
    trprod = 0.5d0*(produc(1,iel)+produc(2,iel)+produc(3,iel))
    trrij  = w8(iel)

!     Calcul de Prod+Phi1+Phi2-Eps
!       = rhoPij-C1rho eps/k(Rij-2/3k dij)-C2rho(Pij-1/3Pkk dij)-2/3rho eps dij
!       = rho{2/3dij[C2 Pkk/2+(C1-1)eps)]+(1-C2)Pij-C1eps/kRij}
    smbr(iel) = smbr(iel) + propce(iel,ipcrom) * volume(iel)      &
      *(   deltij*d2s3*                                           &
           (  crij2*trprod                                        &
            +(crij1-1.d0)* rtpa(iel,ieiph)  )                     &
         +(1.0d0-crij2)*produc(isou,iel)                          &
         -crij1*rtpa(iel,ieiph)/trrij * rtpa(iel,ivar)  )

!     Calcul de la partie implicite issue de Phi1
!       = C1rho eps/k(1-1/3 dij)
    rovsdt(iel) = rovsdt(iel) + propce(iel,ipcrom) * volume(iel)  &
         *(1.d0-d1s3*deltij)*crij1*rtpa(iel,ieiph)/trrij
  enddo

endif

!===============================================================================
! 6. TERMES D'ECHO DE PAROI
!===============================================================================

if(irijec(iphas).eq.1) then

  do iel = 1, ncel
    w7(iel) = 0.d0
  enddo

  call rijech                                                     &
  !==========
 ( idebia , idebra ,                                              &
   nvar   , nscal  , nphas  ,                                     &
   iphas  , ivar   , isou   , ipp    ,                            &
   ia     ,                                                       &
   rtp    , rtpa   , propce , propfa , propfb ,                   &
   coefa  , coefb  , produc , w7   ,                              &
   coefax , viscb  ,                                              &
   w1     , w2     , w3     , w4     , w5     , w6     ,          &
   ra     )

!     Si on extrapole les T.S. : PROPCE
if(isto2t(iphas).gt.0) then
  do iel = 1, ncel
     propce(iel,iptsta+isou-1) =                                  &
     propce(iel,iptsta+isou-1) + w7(iel)
   enddo
!     Sinon SMBR
 else
   do iel = 1, ncel
     smbr(iel) = smbr(iel) + w7(iel)
   enddo
 endif

endif


!===============================================================================
! 7. TERMES DE GRAVITE
!===============================================================================

if(igrari(iphas).eq.1) then

  do iel = 1, ncel
    w7(iel) = 0.d0
  enddo

  call rijthe                                                     &
  !==========
 ( idebia , idebra ,                                              &
   nvar   , nscal  , nphas  ,                                     &
   iphas  , ivar   , isou   , ipp    ,                            &
   ia     ,                                                       &
   rtp    , rtpa   , propce , propfa , propfb ,                   &
   coefa  , coefb  , grarox , graroy , graroz , w7     ,          &
   w1     , w2     , w3     , w4     , w5     , w6     ,          &
   ra     )

!     Si on extrapole les T.S. : PROPCE
if(isto2t(iphas).gt.0) then
  do iel = 1, ncel
     propce(iel,iptsta+isou-1) =                                  &
     propce(iel,iptsta+isou-1) + w7(iel)
   enddo
!     Sinon SMBR
 else
   do iel = 1, ncel
     smbr(iel) = smbr(iel) + w7(iel)
   enddo
 endif

endif


!===============================================================================
! 8. TERMES DE DIFFUSION  A.grad(Rij) : PARTIE EXTRADIAGONALE EXPLICITE
!===============================================================================

! ---> Calcul du grad(Rij)


iccocg = 1
inc = 1

nswrgp = nswrgr(ivar )
imligp = imligr(ivar )
iwarnp = iwarni(ivar )
epsrgp = epsrgr(ivar )
climgp = climgr(ivar )
extrap = extrag(ivar )
iphydp = 0

call grdcel                                                       &
!==========
 ( idebia , idebra ,                                              &
   nphas  ,                                                       &
   ivar   , imrgra , inc    , iccocg , nswrgp , imligp , iphydp , &
   iwarnp , nfecra , epsrgp , climgp , extrap ,                   &
   ia     ,                                                       &
   w1     , w1     , w1     ,                                     &
   rtpa(1,ivar )   , coefa(1,iclvar) , coefb(1,iclvar) ,          &
   w1     , w2     , w3     ,                                     &
!        ------   ------   ------
   w4     , w5     , w6     ,                                     &
   ra     )

! ---> Calcul des termes extradiagonaux de A.grad(Rij)

 do iel = 1, ncel
  trrij = w8(iel)
  cstrij = propce(iel,ipcroo) * csrij *trrij / rtpa(iel,ieiph)
  w4(iel) = cstrij * ( rtpa(iel,ir12ip) * w2(iel)                 &
                      +rtpa(iel,ir13ip) * w3(iel) )
  w5(iel) = cstrij * ( rtpa(iel,ir12ip) * w1(iel)                 &
                      +rtpa(iel,ir23ip) * w3(iel) )
  w6(iel) = cstrij * ( rtpa(iel,ir13ip) * w1(iel)                 &
                      +rtpa(iel,ir23ip) * w2(iel) )
 enddo


! ---> Assemblage de { A.grad(Rij) } .S aux faces

 call vectds                                                      &
!!==========
 (ia     ,                                                        &
  ra(ipond) ,                                                     &
  w4     , w5     , w6     ,                                      &
  viscf  , viscb  , ra     )

init = 1
call divmas(ncelet,ncel,nfac,nfabor,init,nfecra,                  &
                                   ifacel,ifabor,viscf,viscb,w4)

!     Si on extrapole les termes sources
if(isto2t(iphas).gt.0) then
  do iel = 1, ncel
    propce(iel,iptsta+isou-1) =                                   &
    propce(iel,iptsta+isou-1) + w4(iel)
  enddo
!     Sinon
else
  do iel = 1, ncel
    smbr(iel) = smbr(iel) + w4(iel)
  enddo
endif


!===============================================================================
! 9. TERMES DE DIFFUSION  A.grad(Rij) : PARTIE DIAGONALE
!===============================================================================
!     Implicitation de (grad(Rij).n)n en gradient facette
!     Si IDIFRE=1, terme correctif explicite
!        grad(Rij)-(grad(Rij).n)n calcule en gradient cellule
!     Les termes de bord sont uniquement pris en compte dans la partie
!        en (grad(Rij).n)n
!     (W1,W2,W3) contient toujours le gradient de la variable traitee

!     Attention en periodicite on traite le gradient comme si c'etait
!       un vecteur (alors que dans grdcel on l'a fait comme si c'etait
!       un tenseur ...).
!     A modifier eventuellement.

if (idifre(iphas).eq.1) then

  do iel = 1, ncel
    trrij = w8(iel)
    cstrij = propce(iel,ipcroo) * csrij *trrij / rtpa(iel,ieiph)
    w4(iel) = cstrij*rtpa(iel,ir11ip)
    w5(iel) = cstrij*rtpa(iel,ir22ip)
    w6(iel) = cstrij*rtpa(iel,ir33ip)
  enddo

! --->  TRAITEMENT DU PARALLELISME ET DE LA PERIODICITE
!        (il reste des doutes sur la periodicite)

  if (irangp.ge.0.or.iperio.eq.1) then
    call synvec(w1, w2, w3)
    !==========
    call syndia(w4, w5, w6)
    !==========
  endif


  do ifac = 1, nfac

    ii = ifacel(1,ifac)
    jj = ifacel(2,ifac)

    surfn2 = surfan(ifac)**2

    grdpx = 0.5d0*(w1(ii)+w1(jj))
    grdpy = 0.5d0*(w2(ii)+w2(jj))
    grdpz = 0.5d0*(w3(ii)+w3(jj))
    grdsn = grdpx*surfac(1,ifac)+grdpy*surfac(2,ifac)             &
           +grdpz*surfac(3,ifac)
    grdpx = grdpx-grdsn*surfac(1,ifac)/surfn2
    grdpy = grdpy-grdsn*surfac(2,ifac)/surfn2
    grdpz = grdpz-grdsn*surfac(3,ifac)/surfn2

    viscf(ifac)= 0.5d0*(                                          &
          (w4(ii)+w4(jj))*grdpx*surfac(1,ifac)                    &
         +(w5(ii)+w5(jj))*grdpy*surfac(2,ifac)                    &
         +(w6(ii)+w6(jj))*grdpz*surfac(3,ifac))

  enddo

  do ifac = 1, nfabor
    viscb(ifac) = 0.d0
  enddo

  init = 1
  call divmas(ncelet,ncel,nfac,nfabor,init,nfecra,                &
       ifacel,ifabor,viscf,viscb,w1)

!     Si on extrapole les termes sources
  if(isto2t(iphas).gt.0) then
    do iel = 1, ncel
      propce(iel,iptsta+isou-1) =                                 &
      propce(iel,iptsta+isou-1) + w1(iel)
    enddo
!     Sinon
  else
    do iel = 1, ncel
      smbr(iel) = smbr(iel) + w1(iel)
    enddo
  endif

endif


! ---> Viscosite orthotrope pour partie implicite

if( idiff(ivar).ge. 1 ) then
  do iel = 1, ncel
    trrij = w8(iel)
    rctse = propce(iel,ipcrom) * csrij * trrij / rtpa(iel,ieiph)
    w1(iel) = propce(iel,ipcvis)                                  &
                            + idifft(ivar)*rctse*rtpa(iel,ir11ip)
    w2(iel) = propce(iel,ipcvis)                                  &
                            + idifft(ivar)*rctse*rtpa(iel,ir22ip)
    w3(iel) = propce(iel,ipcvis)                                  &
                            + idifft(ivar)*rctse*rtpa(iel,ir33ip)
  enddo

  call visort                                                     &
  !==========
 ( idebia , idebra ,                                              &
   imvisf ,                                                       &
   ia     ,                                                       &
   w1     , w2     , w3     ,                                     &
   viscf  , viscb  ,                                              &
   ra     )

else

  do ifac = 1, nfac
    viscf(ifac) = 0.d0
  enddo
  do ifac = 1, nfabor
    viscb(ifac) = 0.d0
  enddo

endif


!===============================================================================
! 10. RESOLUTION
!===============================================================================

if(isto2t(iphas).gt.0) then
  thetp1 = 1.d0 + thets
  do iel = 1, ncel
    smbr(iel) = smbr(iel) + thetp1*propce(iel,iptsta+isou-1)
  enddo
endif


iconvp = iconv (ivar)
idiffp = idiff (ivar)
ndircp = ndircl(ivar)
ireslp = iresol(ivar)
nitmap = nitmax(ivar)
nswrsp = nswrsm(ivar)
nswrgp = nswrgr(ivar)
imligp = imligr(ivar)
ircflp = ircflu(ivar)
ischcp = ischcv(ivar)
isstpp = isstpc(ivar)
iescap = 0
imgrp  = imgr  (ivar)
ncymxp = ncymax(ivar)
nitmfp = nitmgf(ivar)
!MO      IPP    =
iwarnp = iwarni(ivar)
blencp = blencv(ivar)
epsilp = epsilo(ivar)
epsrsp = epsrsm(ivar)
epsrgp = epsrgr(ivar)
climgp = climgr(ivar)
extrap = extrag(ivar)
relaxp = relaxv(ivar)

call codits                                                       &
!==========
 ( idebia , idebra ,                                              &
   nvar   , nscal  , nphas  ,                                     &
   idtvar , ivar   , iconvp , idiffp , ireslp , ndircp , nitmap , &
   imrgra , nswrsp , nswrgp , imligp , ircflp ,                   &
   ischcp , isstpp , iescap ,                                     &
   imgrp  , ncymxp , nitmfp , ipp    , iwarnp ,                   &
   blencp , epsilp , epsrsp , epsrgp , climgp , extrap ,          &
   relaxp , thetv  ,                                              &
   ia     ,                                                       &
   rtpa(1,ivar)    , rtpa(1,ivar)    ,                            &
                     coefa(1,iclvar) , coefb(1,iclvar) ,          &
                     coefa(1,iclvaf) , coefb(1,iclvaf) ,          &
                     propfa(1,iflmas), propfb(1,iflmab),          &
   viscf  , viscb  , viscf  , viscb  ,                            &
   rovsdt , smbr   , rtp(1,ivar)     ,                            &
   dam    , xam    , drtp   ,                                     &
   w1     , w2     , w3     , w4     , w5     ,                   &
   w6     , w7     , w8     , w9     ,                            &
   ra     )


!===============================================================================
! 11. IMPRESSIONS
!===============================================================================


!--------
! FORMATS
!--------

#if defined(_CS_LANG_FR)

 1000 format(/,'           RESOLUTION POUR LA VARIABLE ',A8,/)

#else

 1000 format(/,'           SOLVING VARIABLE ',A8           ,/)

#endif

!12345678 : MAX: 12345678901234 MIN: 12345678901234 NORM: 12345678901234
!----
! FIN
!----

return

end subroutine
