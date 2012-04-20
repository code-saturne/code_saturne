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

subroutine turbkw &
!================

 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   icepdc , icetsm , itypsm ,                                     &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   tslagr , coefa  , coefb  , ckupdc , smacel )

!===============================================================================
! FONCTION :
! ----------

! RESOLUTION DES EQUATIONS K-OMEGA SST 1 PHASE INCOMPRESSIBLE OU
! RHO VARIABLE SUR UN PAS DE TEMPS

!-------------------------------------------------------------------------------
!ARGU                             ARGUMENTS
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! ncepdp           ! i  ! <-- ! number of cells with head loss                 !
! ncesmp           ! i  ! <-- ! number of cells with mass source term          !
! icepdc(ncelet    ! te ! <-- ! numero des ncepdp cellules avec pdc            !
! icetsm(ncesmp    ! te ! <-- ! numero des cellules a source de masse          !
! itypsm           ! te ! <-- ! type de source de masse pour les               !
! (ncesmp,nvar)    !    !     !  variables (cf. ustsma)                        !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! tslagr           ! tr ! <-- ! terme de couplage retour du                    !
!(ncelet,*)        !    !     !     lagrangien                                 !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
! ckupdc           ! tr ! <-- ! tableau de travail pour pdc                    !
!  (ncepdp,6)      !    !     !                                                !
! smacel           ! tr ! <-- ! valeur des variables associee a la             !
! (ncesmp,*   )    !    !     !  source de masse                               !
!                  !    !     !  pour ivar=ipr, smacel=flux de masse           !
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
use cstnum
use cstphy
use optcal
use lagran
use pointe, only: s2kw, divukw, ifapat, dispar
use parall
use mesh

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal
integer          ncepdp , ncesmp

integer          icepdc(ncepdp)
integer          icetsm(ncesmp), itypsm(ncesmp,nvar)

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(ndimfb,*)
double precision tslagr(ncelet,*)
double precision coefa(ndimfb,*), coefb(ndimfb,*)
double precision ckupdc(ncepdp,6), smacel(ncesmp,nvar)

! Local variables

character*80     chaine
integer          iel   , ifac  , init  , inc   , iccocg, ivar
integer          ii, iivar , iiun  , ifacpt
integer          iclipk, iclipw, isqrt
integer          nswrgp, imligp
integer          iclikp, iclomg
integer          iclvar, iclvaf
integer          iconvp, idiffp, ndircp, ireslp
integer          nitmap, nswrsp, ircflp, ischcp, isstpp, iescap
integer          imgrp , ncymxp, nitmfp
integer          ipcrom, ipbrom, ipcvst, ipcvis, iflmas, iflmab
integer          iwarnp, ipp
integer          iptsta
integer          ipcroo, ipbroo, ipcvto, ipcvlo

double precision rnorm , d2s3, divp23, epz2
double precision deltk , deltw, a11, a12, a22, a21
double precision unsdet, romvsd
double precision prdtur, xk, xw, xeps, xnu
double precision visct , rom
double precision blencp, epsilp, epsrgp, climgp, extrap, relaxp
double precision thetp1, thetak, thetaw, thets, thetap, epsrsp
double precision tuexpk, tuexpw
double precision cdkw, xarg1, xxf1, xgamma, xbeta, sigma, produc
double precision var, vrmin, vrmax

double precision rvoid(1)

double precision, allocatable, dimension(:) :: viscf, viscb
double precision, allocatable, dimension(:) :: dam
double precision, allocatable, dimension(:) :: smbrk, smbrw, rovsdt
double precision, allocatable, dimension(:) :: tinstk, tinstw, xf1
double precision, allocatable, dimension(:,:) :: gradk, grado, grad
double precision, allocatable, dimension(:) :: w1, w2, w3
double precision, allocatable, dimension(:) :: w5, w6
double precision, allocatable, dimension(:) :: w7, w8

!===============================================================================

!===============================================================================
! 1. INITIALISATION
!===============================================================================

! Allocate temporary arrays for the turbulence resolution
allocate(viscf(nfac), viscb(nfabor))
allocate(dam(ncelet))
allocate(smbrk(ncelet), smbrw(ncelet), rovsdt(ncelet))
allocate(tinstk(ncelet), tinstw(ncelet), xf1(ncelet))

! Allocate work arrays
allocate(w1(ncelet), w2(ncelet), w3(ncelet))
allocate(w5(ncelet), w6(ncelet))
allocate(w7(ncelet), w8(ncelet))


epz2 = epzero**2

iclikp = iclrtp(ik ,icoef)
iclomg = iclrtp(iomg,icoef)

ipcrom = ipproc(irom  )
ipcvst = ipproc(ivisct)
ipcvis = ipproc(iviscl)
iflmas = ipprof(ifluma(ik))
iflmab = ipprob(ifluma(ik))
ipbrom = ipprob(irom  )

thets  = thetst

ipcroo = ipcrom
ipbroo = ipbrom
ipcvto = ipcvst
ipcvlo = ipcvis
if(isto2t.gt.0) then
  if (iroext.gt.0) then
    ipcroo = ipproc(iroma)
    ipbroo = ipprob(iroma)
  endif
  if(iviext.gt.0) then
    ipcvto = ipproc(ivista)
    ipcvlo = ipproc(ivisla)
  endif
endif

if(isto2t.gt.0) then
  iptsta = ipproc(itstua)
else
  iptsta = 0
endif

if(iwarni(ik).ge.1) then
  write(nfecra,1000)
endif


!===============================================================================
! 2. CALCUL DE dk/dxj.dw/dxj
!      Le terme est stocke dans         W1
!      En sortie de l'etape on conserve W1
!===============================================================================

! Allocate temporary arrays for gradients calculation
allocate(gradk(ncelet,3), grado(ncelet,3))

iccocg = 1
inc = 1

nswrgp = nswrgr(ik)
imligp = imligr(ik)
iwarnp = iwarni(ik)
epsrgp = epsrgr(ik)
climgp = climgr(ik)
extrap = extrag(ik)

call grdcel &
!==========
 ( ik  , imrgra , inc    , iccocg , nswrgp , imligp ,             &
   iwarnp , nfecra , epsrgp , climgp , extrap ,                   &
   rtpa(1,ik)   , coefa(1,iclikp) , coefb(1,iclikp) ,             &
   gradk  )

nswrgp = nswrgr(iomg)
imligp = imligr(iomg)
iwarnp = iwarni(iomg)
epsrgp = epsrgr(iomg)
climgp = climgr(iomg)
extrap = extrag(iomg)

call grdcel &
!==========
 ( iomg , imrgra , inc    , iccocg , nswrgp , imligp ,            &
   iwarnp , nfecra , epsrgp , climgp , extrap ,                   &
   rtpa(1,iomg)  , coefa(1,iclomg) , coefb(1,iclomg) ,            &
   grado  )

do iel = 1, ncel
  w1(iel) = gradk(iel,1)*grado(iel,1) &
          + gradk(iel,2)*grado(iel,2) &
          + gradk(iel,3)*grado(iel,3)
enddo

! Free memory
deallocate(gradk, grado)

!====================================================
! 3. CALCUL DU COEFFICIENT DE PONDERATION F1
!      Le terme est stocke dans         XF1
!      En sortie de l'etape on conserve W1,XF1
!====================================================


if(abs(icdpar).eq.2) then
  do iel = 1, ncel
    ifacpt = ifapat(iel)
    w2(iel) = (cdgfbo(1,ifacpt)-xyzcen(1,iel))**2 &
            + (cdgfbo(2,ifacpt)-xyzcen(2,iel))**2 &
            + (cdgfbo(3,ifacpt)-xyzcen(3,iel))**2
    w2(iel) = sqrt(w2(iel))
  enddo
else
  do iel = 1, ncel
    w2(iel) =  max(dispar(iel),epzero)
  enddo
endif

!     En cas d'ordre 2 on utilise les valeurs en n car le terme en (1-F1)*W1
!     sera dans PROPCE. Du coup, on aura quand meme certaines "constantes"
!     intervenant dans des termes en n+1/2 (ex sigma_k pour la diffusion) calcules
!     a partir de F1 en n -> mais l'effet sur les "constantes" est faible
!     -> a garder en tete si on fait vraiment de l'ordre 2 en temps en k-omega
do iel = 1, ncel
  rom = propce(iel,ipcroo)
  xnu = propce(iel,ipcvlo)/rom
  xk = rtpa(iel,ik)
  xw  = rtpa(iel,iomg)
  cdkw = 2*rom/ckwsw2/xw*w1(iel)
  cdkw = max(cdkw,1.d-20)
  xarg1 = max(sqrt(xk)/cmu/xw/w2(iel), 500.d0*xnu/xw/w2(iel)**2)
  xarg1 = min(xarg1, 4.d0*rom*xk/ckwsw2/cdkw/w2(iel)**2)
  xf1(iel) = tanh(xarg1**4)
enddo

!===============================================================================
! 4. CALCUL DU TERME DE PRODUCTION
!      Les termes sont stockes dans     TINSTK,TINSTW
!      En sortie de l'etape on conserve W1,XF1,TINSTK,TINSTW
!===============================================================================

d2s3 = 2.d0/3.d0
do iel = 1, ncel
  xk   = rtpa(iel,ik)
  xeps = cmu*rtpa(iel,iomg)*xk
  tinstk(iel) = s2kw(iel) - d2s3*divukw(iel)*divukw(iel)
  tinstw(iel) = propce(iel,ipcvto)*tinstk(iel)                    &
       -d2s3*propce(iel,ipcroo)*xk*divukw(iel)
  tinstk(iel) = min(tinstw(iel),ckwc1*propce(iel,ipcroo)*xeps)
enddo

!===============================================================================
! 4. CALCUL DU TERME DE GRAVITE
!      Les termes sont stockes dans     TINSTK,TINSTW,W2
!      En sortie de l'etape on conserve W1,W2,XF1,TINSTK,TINSTW
!===============================================================================

if(igrake.eq.1) then

  ! Allocate a temporary array for the gradient calculation
  allocate(grad(ncelet,3))

! --- Terme de gravite G = BETA*G*GRAD(SCA)/PRDTUR/RHO
!     Ici on calcule   G =-G*GRAD(RHO)/PRDTUR/RHO

  iccocg = 1
  inc = 1

!     Le choix ci dessous a l'avantage d'etre simple

  nswrgp = nswrgr(ik)
  epsrgp = epsrgr(ik)
  imligp = imligr(ik)
  iwarnp = iwarni(ik)
  climgp = climgr(ik)
  extrap = extrag(ik)

!     Conditions aux limites sur ROM : Dirichlet ROMB
!       On utilise VISCB pour stocker le COEFB relatif a ROM
!       On impose en Dirichlet (COEFA) la valeur ROMB

  do ifac = 1, nfabor
    viscb(ifac) = 0.d0
  enddo

  iivar = 0

  call grdcel                                                     &
  !==========
 ( iivar  , imrgra , inc    , iccocg , nswrgp , imligp ,          &
   iwarnp , nfecra , epsrgp , climgp , extrap ,                   &
   propce(1,ipcroo), propfb(1,ipbroo), viscb  ,                   &
   grad   )


!      Production et terme de gravite
!        TINSTK=MIN(P,C1*EPS)+G et TINSTW=P+(1-CE3)*G
!        On conserve G dans W2 pour la phase de couplage des termes sources

  if(iscalt.gt.0.and.nscal.ge.iscalt) then
    prdtur = sigmas(iscalt)
  else
    prdtur = 1.d0
  endif

  do iel = 1, ncel
    w2(iel) = -(grad(iel,1)*gx + grad(iel,2)*gy + grad(iel,3)*gz) / &
               (propce(iel,ipcroo)*prdtur)
    tinstw(iel)=tinstw(iel)+propce(iel,ipcvto)*max(w2(iel),zero)
    tinstk(iel)=tinstk(iel)+propce(iel,ipcvto)*w2(iel)
  enddo

  ! Free memory
  deallocate(grad)

endif


!===============================================================================
! 5. PRISE EN COMPTE DES TERMES SOURCES UTILISATEURS
!      Les termes sont stockes dans     SMBRK,SMBRW,DAM,W3
!      En sortie de l'etape on conserve W1-3,XF1,TINSTK,TINSTW,SMBRK,SMBRW,DAM
!===============================================================================

do iel = 1, ncel
  smbrk(iel) = 0.d0
  smbrw(iel) = 0.d0
  dam  (iel) = 0.d0
  w3   (iel) = 0.d0
enddo

call ustskw                                                       &
!==========
 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   icepdc , icetsm , itypsm ,                                     &
   dt     , rtpa   , propce , propfa , propfb ,                   &
   coefa  , coefb  , ckupdc , smacel , s2kw   , divukw ,          &
   w1     , w2     , xf1    ,                                     &
   smbrk  , smbrw  , dam    , w3     )
!        ------   ------   ------   ------

!===============================================================================
! 7. ON FINALISE LE CALCUL DES TERMES SOURCES

!      Les termes sont stockes dans     SMBRK, SMBRW
!      En sortie de l'etape on conserve SMBRK,SMBRW,DAM,W1-4
!===============================================================================

do iel = 1, ncel

  visct  = propce(iel,ipcvto)
  rom    = propce(iel,ipcroo)
  xk     = rtpa(iel,ik)
  xw     = rtpa(iel,iomg)
  xxf1   = xf1(iel)
  xgamma = xxf1*ckwgm1 + (1.d0-xxf1)*ckwgm2
  xbeta  = xxf1*ckwbt1 + (1.d0-xxf1)*ckwbt2

  smbrk(iel) = smbrk(iel) + volume(iel)*(                         &
       tinstk(iel)                                                &
       -cmu*rom*xw*xk )

  smbrw(iel) = smbrw(iel) + volume(iel)*(                         &
       rom*xgamma/visct*tinstw(iel)                               &
       -xbeta*rom*xw**2                                           &
       +2.d0*rom/xw*(1.d0-xxf1)/ckwsw2*w1(iel) )

enddo

!===============================================================================
! 8. PRISE EN COMPTE DES TERMES D'ACCUMULATION DE MASSE ET
!         DE LA DEUXIEME PARTIE DES TS UTILISATEURS (PARTIE EXPLICITE)
!         STOCKAGE POUR EXTRAPOLATION EN TEMPS
!      On utilise                       SMBRK,SMBRW
!      En sortie de l'etape on conserve SMBRK,SMBRW,DAM,W2-4

!    Remarque : l'extrapolation telle qu'elle est ecrite n'a pas grand
!               sens si IKECOU=1
!===============================================================================

!     Si on extrapole les T.S.
if(isto2t.gt.0) then

  do iel = 1, ncel

!       Sauvegarde pour echange
    tuexpk = propce(iel,iptsta)
!       Pour la suite et le pas de temps suivant
    propce(iel,iptsta) = smbrk(iel)
!       Termes dependant de la variable resolue et theta PROPCE
    smbrk(iel) = - thets*tuexpk
!       On suppose -DAM > 0 : on implicite
!         le terme utilisateur dependant de la variable resolue
    smbrk(iel) = dam(iel)*rtpa(iel,ik) + smbrk(iel)

!       Sauvegarde pour echange
    tuexpw = propce(iel,iptsta+1)
!       Pour la suite et le pas de temps suivant
    propce(iel,iptsta+1) = smbrw(iel)
!       Termes dependant de la variable resolue et theta PROPCE
    smbrw(iel) = - thets*tuexpw
!       On suppose -W3 > 0 : on implicite
!         le terme utilisateur dependant de la variable resolue
    smbrw(iel) =  w3(iel)*rtpa(iel,iomg) + smbrw(iel)

  enddo

!     Si on n'extrapole pas les T.S.
else
  do iel = 1, ncel
    smbrk(iel) = smbrk(iel) + dam(iel)*rtpa(iel,ik)
    smbrw(iel) = smbrw(iel) + w3 (iel)*rtpa(iel,iomg)
  enddo
endif

!===============================================================================
! 8.1 PRISE EN COMPTE DES TERMES SOURCES LAGRANGIEN : PARTIE EXPLICITE
!     COUPLAGE RETOUR
!===============================================================================

!     Ordre 2 non pris en compte
if (iilagr.eq.2 .and. ltsdyn.eq.1) then

  do iel = 1,ncel

! Termes sources explicte et implicte sur k

    smbrk(iel)  = smbrk(iel) + tslagr(iel,itske)

! Termes sources explicte sur omega : on reprend la constante CE4 directement
!    du k-eps sans justification ... a creuser si necessaire           !

    smbrw(iel)  = smbrw(iel)                                      &
                + ce4 *tslagr(iel,itske) * propce(iel,ipcroo)     &
                /propce(iel,ipcvto)
  enddo

endif


!===============================================================================
! 9. PRISE EN COMPTE DES TERMES DE CONV/DIFF DANS LE SECOND MEMBRE

!      Tableaux de travail              W7, W8, W1, TINSTK, TINSTW
!      Les termes sont stockes dans     W5 et W6, puis ajoutes a SMBRK, SMBRW
!      En sortie de l'etape on conserve W2-6,SMBRK,SMBRW,DAM
!===============================================================================

!     Ceci ne sert a rien si IKECOU n'est pas egal a 1

if (ikecou.eq.1) then

  do iel = 1, ncel
    w5 (iel) = 0.d0
    w6 (iel) = 0.d0
  enddo

! ---> Traitement de k

  ivar   = ik

  ipp    = ipprtp(ivar)

  iclvar = iclrtp(ivar,icoef )
  iclvaf = iclrtp(ivar,icoeff)
  chaine = nomvar(ipp)

  if( idiff(ivar).ge. 1 ) then

    do iel = 1, ncel
      xxf1 = xf1(iel)
      sigma = xxf1*ckwsk1 + (1.d0-xxf1)*ckwsk2
      w7(iel) = propce(iel,ipcvis)                                &
           + idifft(ivar)*propce(iel,ipcvst)/sigma
    enddo
    call viscfa                                                   &
    !==========
 ( imvisf ,                                                       &
   w7     ,                                                       &
   viscf  , viscb  )

  else

    do ifac = 1, nfac
      viscf(ifac) = 0.d0
    enddo
    do ifac = 1, nfabor
      viscb(ifac) = 0.d0
    enddo

  endif

  iccocg = 1
  inc    = 1
  iconvp = iconv (ivar)
  idiffp = idiff (ivar)
  nswrgp = nswrgr(ivar)
  imligp = imligr(ivar)
  ircflp = ircflu(ivar)
  ischcp = ischcv(ivar)
  isstpp = isstpc(ivar)
  iwarnp = iwarni(ivar)
  blencp = blencv(ivar)
  epsrgp = epsrgr(ivar)
  climgp = climgr(ivar)
  extrap = extrag(ivar)
  relaxp = relaxv(ivar)
  thetap = thetav(ivar)

  call bilsc2                                                     &
  !==========
 ( nvar   , nscal  ,                                              &
   idtvar , ivar   , iconvp , idiffp , nswrgp , imligp , ircflp , &
   ischcp , isstpp , inc    , imrgra , iccocg ,                   &
   ipp    , iwarnp ,                                              &
   blencp , epsrgp , climgp , extrap , relaxp , thetap ,          &
   rtpa(1,ivar)    , rtpa(1,ivar)    ,                            &
   coefa(1,iclvar) , coefb(1,iclvar) ,                            &
   coefa(1,iclvaf) , coefb(1,iclvaf) ,                            &
   propfa(1,iflmas), propfb(1,iflmab), viscf  , viscb  ,          &
   w5     )


  if (iwarni(ivar).ge.2) then
    isqrt = 1
    call prodsc(ncel,isqrt,smbrk,smbrk,rnorm)
    write(nfecra,1100) chaine(1:8) ,rnorm
  endif


! ---> Traitement de omega

  ivar   = iomg

  ipp    = ipprtp(ivar)

  iclvar = iclrtp(ivar,icoef )
  iclvaf = iclrtp(ivar,icoeff)
  chaine = nomvar(ipp)

  if( idiff(ivar).ge. 1 ) then
    do iel = 1, ncel
      xxf1 = xf1(iel)
      sigma = xxf1*ckwsw1 + (1.d0-xxf1)*ckwsw2
      w7(iel) = propce(iel,ipcvis)                                &
           + idifft(ivar)*propce(iel,ipcvst)/sigma
    enddo
    call viscfa                                                   &
    !==========
 ( imvisf ,                                                       &
   w7     ,                                                       &
   viscf  , viscb  )

  else

    do ifac = 1, nfac
      viscf(ifac) = 0.d0
    enddo
    do ifac = 1, nfabor
      viscb(ifac) = 0.d0
    enddo

  endif

  iccocg = 1
  inc    = 1
  iconvp = iconv (ivar)
  idiffp = idiff (ivar)
  nswrgp = nswrgr(ivar)
  imligp = imligr(ivar)
  ircflp = ircflu(ivar)
  ischcp = ischcv(ivar)
  isstpp = isstpc(ivar)
  iwarnp = iwarni(ivar)
  blencp = blencv(ivar)
  epsrgp = epsrgr(ivar)
  climgp = climgr(ivar)
  extrap = extrag(ivar)
  relaxp = relaxv(ivar)
  thetap = thetav(ivar)

  call bilsc2                                                     &
  !==========
 ( nvar   , nscal  ,                                              &
   idtvar , ivar   , iconvp , idiffp , nswrgp , imligp , ircflp , &
   ischcp , isstpp , inc    , imrgra , iccocg ,                   &
   ipp    , iwarnp ,                                              &
   blencp , epsrgp , climgp , extrap , relaxp , thetap ,          &
   rtpa(1,ivar)    , rtpa(1,ivar)    ,                            &
   coefa(1,iclvar) , coefb(1,iclvar) ,                            &
   coefa(1,iclvaf) , coefb(1,iclvaf) ,                            &
   propfa(1,iflmas), propfb(1,iflmab), viscf  , viscb  ,          &
   w6  )
!        --

  if (iwarni(ivar).ge.2) then
    isqrt = 1
    call prodsc(ncel,isqrt,smbrw,smbrw,rnorm)
    write(nfecra,1100) chaine(1:8) ,rnorm
  endif

  do iel = 1,ncel
    smbrk(iel) = smbrk(iel) + w5(iel)
    smbrw(iel) = smbrw(iel) + w6(iel)
  enddo

endif

!===============================================================================
! 10. AJOUT DES TERMES SOURCES DE MASSE EXPLICITES

!       Les parties implicites eventuelles sont conservees dans W7 et W8
!         et utilisees dans la phase d'implicitation cv/diff

!       Les termes sont stockes dans     SMBRK, SMBRW, W7, W8
!       En sortie de l'etape on conserve W2-8,SMBRK,SMBRW,DAM
!===============================================================================

if (ncesmp.gt.0) then

  do iel = 1, ncel
    w7(iel) = 0.d0
    w8(iel) = 0.d0
  enddo

!       Entier egal a 1 (pour navsto : nb de sur-iter)
  iiun = 1

!       On incremente SMBRS par -Gamma RTPA et ROVSDT par Gamma (*theta)
  ivar = ik
  call catsma &
  !==========
 ( ncelet , ncel   , ncesmp , iiun   ,                            &
                                 isto2t , thetav(ivar) ,   &
   icetsm , itypsm(1,ivar) ,                                      &
   volume , rtpa(1,ivar) , smacel(1,ivar) , smacel(1,ipr) ,       &
   smbrk  , w7     , tinstk )
  ivar = iomg
  call catsma &
  !==========
 ( ncelet , ncel   , ncesmp , iiun   ,                            &
                                 isto2t , thetav(ivar) ,   &
   icetsm , itypsm(1,ivar) ,                                      &
   volume , rtpa(1,ivar) , smacel(1,ivar) , smacel(1,ipr) ,       &
   smbrw  , w8     , tinstw )

!       Si on extrapole les TS on met Gamma Pinj dans PROPCE
  if(isto2t.gt.0) then
    do iel = 1, ncel
      propce(iel,iptsta  ) = propce(iel,iptsta  ) + tinstk(iel)
      propce(iel,iptsta+1) = propce(iel,iptsta+1) + tinstw(iel)
    enddo
!       Sinon on le met directement dans SMBR
  else
    do iel = 1, ncel
      smbrk(iel) = smbrk(iel) + tinstk(iel)
      smbrw(iel) = smbrw(iel) + tinstw(iel)
    enddo
  endif

endif

!     Finalisation des termes sources
if(isto2t.gt.0) then
  thetp1 = 1.d0 + thets
  do iel = 1, ncel
    smbrk(iel) = smbrk(iel) + thetp1 * propce(iel,iptsta)
    smbrw(iel) = smbrw(iel) + thetp1 * propce(iel,iptsta+1)
  enddo
endif

!===============================================================================
! 11. INCREMENTS DES TERMES SOURCES DANS LE SECOND MEMBRE

!       Les termes sont stockes dans     SMBRK, SMBRW
!       En sortie de l'etape on conserve W3-8,SMBRK,SMBRW
!===============================================================================

!     Ordre 2 non pris en compte
if(ikecou.eq.1) then

  do iel = 1, ncel

    rom = propce(iel,ipcrom)

!   RESOLUTION COUPLEE

    romvsd     = 1.d0/(rom*volume(iel))
    smbrk(iel) = smbrk(iel)*romvsd
    smbrw(iel) = smbrw(iel)*romvsd
    divp23     = d2s3*max(divukw(iel),zero)
    produc     = s2kw(iel)-d2s3*divukw(iel)**2+w2(iel)
    xk         = rtpa(iel,ik)
    xw         = rtpa(iel,iomg)
    xxf1       = xf1(iel)
    xgamma     = xxf1*ckwgm1 + (1.d0-xxf1)*ckwgm2
    xbeta      = xxf1*ckwbt1 + (1.d0-xxf1)*ckwbt2

    a11 = 1.d0/dt(iel)                                            &
         - 1.d0/xw*min(produc,zero)+divp23+cmu*xw
    a12 = cmu*xk
    a21 = 0.d0
    a22 = 1.d0/dt(iel)+xgamma*divp23+2.d0*xbeta*xw

    unsdet = 1.d0/(a11*a22 -a12*a21)

    deltk = ( a22*smbrk(iel) -a12*smbrw(iel) )*unsdet
    deltw = (-a21*smbrk(iel) +a11*smbrw(iel) )*unsdet

!     NOUVEAU TERME SOURCE POUR CODITS

    romvsd = rom*volume(iel)/dt(iel)

    smbrk(iel) = romvsd*deltk
    smbrw(iel) = romvsd*deltw

  enddo

endif

!===============================================================================
! 12. TERMES INSTATIONNAIRES

!     Les termes sont stockes dans     TINSTK, TINSTW
!     En sortie de l'etape on conserve SMBRK, SMBRW,  TINSTK, TINSTW
!===============================================================================

! --- PARTIE EXPLICITE

!     on enleve la convection/diffusion au temps n a SMBRK et SMBRW
!     s'ils ont ete calcules
if (ikecou.eq.1) then
  do iel = 1, ncel
    smbrk(iel) = smbrk(iel) - w5(iel)
    smbrw(iel) = smbrw(iel) - w6(iel)
  enddo
endif

! --- RHO/DT et DIV

do iel = 1, ncel
  rom = propce(iel,ipcrom)
  romvsd = rom*volume(iel)/dt(iel)
  tinstk(iel) = istat(ik)*romvsd
  tinstw(iel) = istat(iomg)*romvsd
enddo

! --- Source de masse (le theta est deja inclus par catsma)
if (ncesmp.gt.0) then
  do iel = 1, ncel
    tinstk(iel) = tinstk(iel) + w7(iel)
    tinstw(iel) = tinstw(iel) + w8(iel)
  enddo
endif

! --- Termes sources utilisateurs
if(isto2t.gt.0) then
  thetak = thetav(ik)
  thetaw = thetav(iomg)
  do iel = 1, ncel
    tinstk(iel) = tinstk(iel) -dam(iel)*thetak
    tinstw(iel) = tinstw(iel) -w3 (iel)*thetaw
  enddo
else
  do iel = 1, ncel
    tinstk(iel) = tinstk(iel) + max(-dam(iel),zero)
    tinstw(iel) = tinstw(iel) + max(-w3 (iel),zero)
  enddo
endif

! --- PRISE EN COMPTE DES TERMES LAGRANGIEN : COUPLAGE RETOUR

!     Ordre 2 non pris en compte
if (iilagr.eq.2 .and. ltsdyn.eq.1) then

  do iel = 1,ncel

! Termes sources implicite sur k

    tinstk(iel) = tinstk(iel) + max(-tslagr(iel,itsli),zero)

! Termes sources implicte sur omega

    tinstw(iel) = tinstw(iel)                                     &
          + max( (-ce4*tslagr(iel,itske)/rtpa(iel,ik)) , zero)

  enddo

endif

! Si IKECOU=0, on implicite plus fortement k et omega

if(ikecou.eq.0)then
  do iel=1,ncel
    xw    = rtpa(iel,iomg)
    xxf1  = xf1(iel)
    xbeta = xxf1*ckwbt1 + (1.d0-xxf1)*ckwbt2
    rom = propce(iel,ipcrom)
      tinstk(iel) = tinstk(iel) +                                 &
         volume(iel)*cmu*rom*xw
      tinstw(iel) = tinstw(iel) +                                 &
         volume(iel)*xbeta*rom*xw
  enddo
endif


!===============================================================================
! 13. RESOLUTION

!       On utilise                      SMBRK, SMBRW,  TINSTK, TINSTW
!===============================================================================

! ---> Traitement de k

ivar = ik
iclvar = iclrtp(ivar,icoef )
iclvaf = iclrtp(ivar,icoeff)

ipp    = ipprtp(ivar)

!     "VITESSE" DE DIFFUSION FACETTE

if( idiff(ivar).ge. 1 ) then

  do iel = 1, ncel
    xxf1 = xf1(iel)
    sigma = xxf1*ckwsk1 + (1.d0-xxf1)*ckwsk2
    w1(iel) = propce(iel,ipcvis)                                  &
                        + idifft(ivar)*propce(iel,ipcvst)/sigma
  enddo
  call viscfa                                                     &
  !==========
 ( imvisf ,                                                       &
   w1     ,                                                       &
   viscf  , viscb  )

else

  do ifac = 1, nfac
    viscf(ifac) = 0.d0
  enddo
  do ifac = 1, nfabor
    viscb(ifac) = 0.d0
  enddo

endif

!     RESOLUTION POUR K

iconvp = iconv (ivar)
idiffp = idiff (ivar)
ireslp = iresol(ivar)
ndircp = ndircl(ivar)
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
thetap = thetav(ivar)

call codits                                                       &
!==========
 ( nvar   , nscal  ,                                              &
   idtvar , ivar   , iconvp , idiffp , ireslp , ndircp , nitmap , &
   imrgra , nswrsp , nswrgp , imligp , ircflp ,                   &
   ischcp , isstpp , iescap ,                                     &
   imgrp  , ncymxp , nitmfp , ipp    , iwarnp ,                   &
   blencp , epsilp , epsrsp , epsrgp , climgp , extrap ,          &
   relaxp , thetap ,                                              &
   rtpa(1,ivar)    , rtpa(1,ivar)    ,                            &
                     coefa(1,iclvar) , coefb(1,iclvar) ,          &
                     coefa(1,iclvaf) , coefb(1,iclvaf) ,          &
                     propfa(1,iflmas), propfb(1,iflmab),          &
   viscf  , viscb  , viscf  , viscb  ,                            &
   tinstk , smbrk  , rtp(1,ivar)     ,                            &
   rvoid  )


! ---> Traitement de omega

ivar = iomg
iclvar = iclrtp(ivar,icoef )
iclvaf = iclrtp(ivar,icoeff)

ipp    = ipprtp(ivar)


!     "VITESSE" DE DIFFUSION FACETTE

if( idiff(ivar).ge. 1 ) then
  do iel = 1, ncel
    xxf1 = xf1(iel)
    sigma = xxf1*ckwsw1 + (1.d0-xxf1)*ckwsw2
    w1(iel) = propce(iel,ipcvis)                                  &
                        + idifft(ivar)*propce(iel,ipcvst)/sigma
  enddo
  call viscfa                                                     &
  !==========
 ( imvisf ,                                                       &
   w1     ,                                                       &
   viscf  , viscb  )

else

  do ifac = 1, nfac
    viscf(ifac) = 0.d0
  enddo
  do ifac = 1, nfabor
    viscb(ifac) = 0.d0
  enddo

endif

! RESOLUTION POUR OMEGA

iconvp = iconv (ivar)
idiffp = idiff (ivar)
ireslp = iresol(ivar)
ndircp = ndircl(ivar)
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
thetap = thetav(ivar)

call codits                                                       &
!==========
 ( nvar   , nscal  ,                                              &
   idtvar , ivar   , iconvp , idiffp , ireslp , ndircp , nitmap , &
   imrgra , nswrsp , nswrgp , imligp , ircflp ,                   &
   ischcp , isstpp , iescap ,                                     &
   imgrp  , ncymxp , nitmfp , ipp    , iwarnp ,                   &
   blencp , epsilp , epsrsp , epsrgp , climgp , extrap ,          &
   relaxp , thetap ,                                              &
   rtpa(1,ivar)    , rtpa(1,ivar)    ,                            &
                     coefa(1,iclvar) , coefb(1,iclvar) ,          &
                     coefa(1,iclvaf) , coefb(1,iclvaf) ,          &
                     propfa(1,iflmas), propfb(1,iflmab),          &
   viscf  , viscb  , viscf  , viscb  ,                            &
   tinstw , smbrw  , rtp(1,ivar)     ,                            &
   rvoid  )


!===============================================================================
! 14. CLIPPING
!===============================================================================

!     Calcul des Min/Max avant clipping, pour affichage
do ii = 1, 2
  if(ii.eq.1) then
    ivar = ik
  elseif(ii.eq.2) then
    ivar = iomg
  endif
  ipp  = ipprtp(ivar)

  vrmin =  grand
  vrmax = -grand
  do iel = 1, ncel
    var = rtp(iel,ivar)
    vrmin = min(vrmin,var)
    vrmax = max(vrmax,var)
  enddo
  if (irangp.ge.0) then
    call parmax (vrmax)
    !==========
    call parmin (vrmin)
    !==========
  endif
  varmna(ipp) = vrmin
  varmxa(ipp) = vrmax

enddo

!     On clippe simplement k et omega par valeur absolue
iclipk = 0
iclipw = 0
do iel = 1, ncel
  xk = rtp(iel,ik )
  xw = rtp(iel,iomg)
  if (abs(xk).le.epz2) then
    iclipk = iclipk + 1
    rtp(iel,ik) = max(rtp(iel,ik),epz2)
  elseif(xk.le.0.d0) then
    iclipk = iclipk + 1
    rtp(iel,ik) = -xk
  endif
  if (abs(xw).le.epz2) then
    iclipw = iclipw + 1
    rtp(iel,iomg) = max(rtp(iel,iomg),epz2)
  elseif(xw.le.0.d0) then
    iclipw = iclipw + 1
    rtp(iel,iomg) = -xw
  endif
enddo

if (irangp.ge.0) then
  call parcpt (iclipk)
  !==========
  call parcpt (iclipw)
  !==========
endif

! ---  Stockage nb de clippings pour listing

iclpmn(ipprtp(ik )) = iclipk
iclpmn(ipprtp(iomg)) = iclipw


! Free memory
deallocate(viscf, viscb)
deallocate(dam)
deallocate(smbrk, smbrw, rovsdt)
deallocate(tinstk, tinstw, xf1)
deallocate(w1, w2, w3)
deallocate(w5, w6)
deallocate(w7, w8)

!--------
! FORMATS
!--------

#if defined(_CS_LANG_FR)

 1000 format(/,                                                   &
'   ** RESOLUTION DU K-OMEGA                     ',/,&
'      ---------------------                     ',/)
 1100 format(1X,A8,' : BILAN EXPLICITE = ',E14.5)

#else

 1000 format(/,                                                   &
'   ** SOLVING K-OMEGA'                           ,/,&
'      ---------------'                           ,/)
 1100 format(1X,A8,' : EXPLICIT BALANCE = ',E14.5)

#endif

!----
! FIN
!----

return

end subroutine
