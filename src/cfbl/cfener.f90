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

subroutine cfener &
!================

 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   iscal  ,                                                       &
   icepdc , icetsm , itypsm ,                                     &
   dt     , rtp    , rtpa   , propce , propfb ,                   &
   coefa  , coefb  , ckupdc , smacel ,                            &
   viscf  , viscb  ,                                              &
   smbrs  , rovsdt )

!===============================================================================
! FONCTION :
! ----------

! RESOLUTION DES EQUATIONS CONVECTION DIFFUSION TERME SOURCE
!   POUR L'ENERGIE TOTALE SUR UN PAS DE TEMPS

!-------------------------------------------------------------------------------
!ARGU                             ARGUMENTS
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! ncepdp           ! i  ! <-- ! number of cells with head loss                 !
! ncesmp           ! i  ! <-- ! number of cells with mass source term          !
! iscal            ! i  ! <-- ! scalar number                                  !
! icepdc(ncelet    ! te ! <-- ! numero des ncepdp cellules avec pdc            !
! icetsm(ncesmp    ! te ! <-- ! numero des cellules a source de masse          !
! itypsm           ! te ! <-- ! type de source de masse pour les               !
! (ncesmp,nvar)    !    !     !  variables (cf. ustsma)                        !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
! ckupdc           ! tr ! <-- ! tableau de travail pour pdc                    !
!  (ncepdp,6)      !    !     !                                                !
! smacel           ! tr ! <-- ! valeur des variables associee a la             !
! (ncesmp,*   )    !    !     !  source de masse                               !
!                  !    !     !  pour ivar=ipr, smacel=flux de masse           !
! viscf(nfac)      ! tr ! --- ! visc*surface/dist aux faces internes           !
! viscb(nfabor     ! tr ! --- ! visc*surface/dist aux faces de bord            !
! smbrs(ncelet     ! tr ! --- ! tableau de travail pour sec mem                !
! rovsdt(ncelet    ! tr ! --- ! tableau de travail pour terme instat           !
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
use entsor
use optcal
use cstphy
use cstnum
use parall
use period
use ppppar
use ppthch
use ppincl
use cfpoin
use mesh
use field

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal
integer          ncepdp , ncesmp
integer          iscal

integer          icepdc(ncepdp)
integer          icetsm(ncesmp), itypsm(ncesmp,nvar)

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision ckupdc(ncepdp,6), smacel(ncesmp,nvar)
double precision viscf(nfac), viscb(nfabor)
double precision smbrs(ncelet)
double precision rovsdt(ncelet)

! Local variables

character*80     chaine
integer          ivar
integer          ifac  , iel
integer          init  , isqrt , iii
integer          iclvar, iclvaf
integer          ipcrom, ipcvst, ipcvsl, iflmas, iflmab
integer          ippvar, ipp
integer          nswrgp, imligp, iwarnp
integer          iconvp, idiffp, ndircp, ireslp, nitmap
integer          nswrsp, ircflp, ischcp, isstpp, iescap
integer          imgrp , ncymxp, nitmfp
double precision epsrgp, climgp, extrap, blencp, epsilp
double precision sclnor, thetap, epsrsp

integer          inc    , iccocg
integer          ivar0  , iij , ii , jj
integer          iccfth , imodif
integer          iel1  , iel2
integer          iterns

double precision flux
double precision dijpfx, dijpfy, dijpfz, pnd  , pip   , pjp
double precision diipfx, diipfy, diipfz, djjpfx, djjpfy, djjpfz

double precision rvoid(1)

double precision, allocatable, dimension(:) :: wb
double precision, allocatable, dimension(:) :: coefap, coefbp
double precision, allocatable, dimension(:,:) :: grad
double precision, allocatable, dimension(:) :: w1
double precision, allocatable, dimension(:) :: w4, w5, w6
double precision, allocatable, dimension(:) :: w7, w8, w9
double precision, dimension(:), pointer :: imasfl, bmasfl

!===============================================================================
!===============================================================================
! 1. INITIALISATION
!===============================================================================

! Allocate a temporary array
allocate(wb(nfabor))

! Allocate work arrays
allocate(grad(ncelet,3))
allocate(w1(ncelet))
allocate(w4(ncelet), w5(ncelet), w6(ncelet))
allocate(w7(ncelet), w8(ncelet), w9(ncelet))


! --- Numero de variable de calcul et de post associe au scalaire traite
ivar   = isca(iscal)
ippvar = ipprtp(ivar)

! --- Numero des conditions aux limites
iclvar = iclrtp(ivar,icoef)
iclvaf = iclrtp(ivar,icoeff)

! --- Numero des grandeurs physiques
ipcrom = ipproc(irom  )
ipcvst = ipproc(ivisct)

call field_get_key_int(ivarfl(ivar), kimasf, iflmas)
call field_get_key_int(ivarfl(ivar), kbmasf, iflmab)
call field_get_val_s(iflmas, imasfl)
call field_get_val_s(iflmab, bmasfl)

if(ivisls(iscal).gt.0) then
  ipcvsl = ipproc(ivisls(iscal))
else
  ipcvsl = 0
endif

! --- Impressions
chaine = nomvar(ippvar)

if(iwarni(ivar).ge.1) then
  write(nfecra,1000) chaine(1:8)
endif

!===============================================================================
! 2. TERMES SOURCES
!===============================================================================

! --> Theta-schema de resolution

! Pour l'instant on prend THETA=1 et on ne code pas le theta-schema

! --> Initialisation

do iel = 1, ncel
  smbrs(iel) = 0.d0
enddo
do iel = 1, ncel
  rovsdt(iel) = 0.d0
enddo


!     TERME SOURCE VOLUMIQUE DE CHALEUR : RHO*PHI *VOLUME
!     =================================          v

call ustssc                                                       &
!==========
 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   iscal  ,                                                       &
   icepdc , icetsm , itypsm ,                                     &
   dt     , rtpa   , rtp    , propce , propfb ,                   &
   ckupdc , smacel , smbrs  , rovsdt )

do iel = 1, ncel
  smbrs(iel) = smbrs(iel) + rovsdt(iel)*rtp(iel,ivar)
  rovsdt(iel) = max(-rovsdt(iel),zero)
enddo


!     TERMES DE SOURCE DE MASSE
!     =========================

!     GAMMA(IEL) = SMACEL(IEL,IPR)

!     Terme implicite : GAMMA*VOLUME
!                                                        n
!     Terme explicite : GAMMA*VOLUME*e   - GAMMA*VOLUME*e
!                                     inj
if (ncesmp.gt.0) then
  iterns = 1
  call catsma ( ncelet , ncel , ncesmp , iterns ,                 &
                isno2t, thetav(ivar),                             &
                icetsm , itypsm(1,ivar) ,                         &
                volume , rtpa(1,ivar) , smacel(1,ivar) ,          &
                smacel(1,ipr) , smbrs , rovsdt , w1    )
endif

!                                      RHO*VOLUME
!     TERME INSTATIONNAIRE IMPLICITE : ----------
!     ==============================       DT

do iel = 1, ncel
  rovsdt(iel) = rovsdt(iel)                                       &
           + istat(ivar)*(propce(iel,ipcrom)/dt(iel))*volume(iel)
enddo

!                                       __        v
!     TERME DE DISSIPATION VISQUEUSE  : >  ((SIGMA *U).n)  *S
!     ==============================    --               ij  ij

if( idiff(iu).ge. 1 ) then
!                             ^^^
  call cfdivs                                                     &
  !==========
 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   icepdc , icetsm , itypsm ,                                     &
   rtpa   , propce , propfb ,                                     &
   coefa  , coefb  , ckupdc , smacel ,                            &
   smbrs  , rtp(1,iu), rtp(1,iv), rtp(1,iw) )
!        ------

endif


!                                         __   P        n+1
!     TERME DE TRANSPORT DE PRESSION  : - >  (---)  *(Q    .n)  *S
!     ==============================      --  RHO ij   pr     ij  ij


if(igrdpp.gt.0) then
  do iel = 1, ncel
    w9(iel) = rtp(iel,isca(irho))
  enddo
else
  do iel = 1, ncel
    w9(iel) = rtpa(iel,isca(irho))
  enddo
endif

!     Avec Reconstruction : ca pose probleme pour l'instant


!   Calcul du gradient de P/RHO

!      do iel = 1, ncel
!        w7(iel) = rtp(iel,ipr)/w9(iel)
!      enddo

! Rq : A defaut de connaitre les parametres pour P/RHO on prend ceux de P

!      iii = ipr
!      inc = 1
!      iccocg = 1
!      nswrgp = nswrgr(iii)
!      imligp = imligr(iii)
!      iwarnp = iwarni(iii)
!      epsrgp = epsrgr(iii)
!      climgp = climgr(iii)
!      extrap = extrag(iii)

!       On alloue localement 2 tableaux de NFABOR pour le calcul
!       de COEFA et COEFB de P/RHO

!      allocate(coefap(nfabor))
!      allocate(coefbp(nfabor))

!      do ifac = 1, nfabor
!        coefap(ifac) = zero
!        coefbp(ifac) = 1.d0
!      enddo

! En periodique et parallele, echange avant calcul du gradient
!      if (irangp.ge.0.or.iperio.eq.1) then
!        call synsca(w7)
!        !==========
!      endif

!  IVAR0 = 0 (indique pour la periodicite de rotation que la variable
!     n'est pas la vitesse ni Rij)
!      ivar0 = 0
!      call grdcel
!      !==========
!     & ( ivar0  , imrgra , inc    , iccocg , nswrgp , imligp ,
!     &   iwarnp , nfecra , epsrgp , climgp , extrap ,
!     &   w7     , coefap , coefbp ,
!     &   grad   )

!     Faces internes
!      do ifac = 1, nfac

!        ii = ifacel(1,ifac)
!        jj = ifacel(2,ifac)

!        dijpfx = dijpf(1,ifac)
!        dijpfy = dijpf(2,ifac)
!        dijpfz = dijpf(3,ifac)

!        pnd   = pond(ifac)

!        Calcul II' et JJ'

!        diipfx = cdgfac(1,ifac) - (xyzcen(1,ii)+ (1.d0-pnd) * dijpfx)
!        diipfy = cdgfac(2,ifac) - (xyzcen(2,ii)+ (1.d0-pnd) * dijpfy)
!        diipfz = cdgfac(3,ifac) - (xyzcen(3,ii)+ (1.d0-pnd) * dijpfz)
!        djjpfx = cdgfac(1,ifac) -  xyzcen(1,jj)+       pnd  * dijpfx
!        djjpfy = cdgfac(2,ifac) -  xyzcen(2,jj)+       pnd  * dijpfy
!        djjpfz = cdgfac(3,ifac) -  xyzcen(3,jj)+       pnd  * dijpfz

!        pip = w7(ii) +grad(ii,1)*diipfx+grad(ii,2)*diipfy+grad(ii,3)*diipfz
!        pjp = w7(jj) +grad(jj,1)*djjpfx+grad(jj,2)*djjpfy+grad(jj,3)*djjpfz

!        flui = (imasfl(ifac)+abs(imasfl(ifac)))
!        fluj = (imasfl(ifac)-abs(imasfl(ifac)))

!        viscf(ifac) = -(pnd*pip*flui+pnd*pjp*fluj)

!      enddo

!     Sans Reconstruction

!     En periodique et parallele, echange avant utilisation
!       des valeurs aux faces
if (irangp.ge.0.or.iperio.eq.1) then
  call synsca(w9)
  !==========
endif

!     Faces internes
do ifac = 1, nfac
  iel1 = ifacel(1,ifac)
  iel2 = ifacel(2,ifac)
  viscf(ifac) =                                                   &
     - rtp(iel1,ipr)/w9(iel1)                              &
     *0.5d0*( imasfl(ifac) +abs(imasfl(ifac)) )     &
     - rtp(iel2,ipr)/w9(iel2)                              &
     *0.5d0*( imasfl(ifac) -abs(imasfl(ifac)) )
enddo

!     Faces de bord : pour les faces ou on a calcule un flux de Rusanov,
!       on remplace la contribution standard par le flux de Rusanov qui
!       contient tous les flux convectifs (et il faudra donc eliminer le
!       flux convectif dans cfbsc2)

do ifac = 1, nfabor
  if (ifbrus(ifac).eq.0) then

    iel = ifabor(ifac)
    viscb(ifac) = - bmasfl(ifac)                         &
  * ( coefa(ifac,iclrtp(ipr,icoef))                      &
    + coefb(ifac,iclrtp(ipr,icoef))*rtp(iel,ipr) )&
  / ( coefa(ifac,iclrtp(isca(irho),icoef))               &
    + coefb(ifac,iclrtp(isca(irho),icoef))*w9(iel) )

  else
    viscb(ifac) = - propfb(ifac,ipprob(ifbene))
  endif
enddo

!     Divergence
init = 0
call divmas(ncelet,ncel,nfac,nfabor,init,nfecra,                  &
              ifacel,ifabor,viscf,viscb,smbrs)


!     TERME DE FORCES DE PESANTEUR : RHO*g.U *VOLUME
!     ============================

do iel = 1, ncel
  smbrs(iel) = smbrs(iel) + w9(iel)*volume(iel)                   &
                           *( gx*rtp(iel,iu)               &
                            + gy*rtp(iel,iv)               &
                            + gz*rtp(iel,iw) )
enddo

!                                       Kij*Sij           LAMBDA   Cp   MUT
!     "VITESSE" DE DIFFUSION FACETTE : --------- avec K = ------ + -- .------
!     ==============================    IJ.nij              Cv     Cv  SIGMAS

if( idiff(ivar).ge. 1 ) then

!     MUT/SIGMAS
  do iel = 1, ncel
    w1(iel) = propce(iel,ipcvst)/sigmas(iscal)
  enddo
!     CP*MUT/SIGMAS
  if(icp.gt.0) then
    do iel = 1, ncel
      w1(iel) = w1(iel)*propce(iel,ipproc(icp))
    enddo
  else
    do iel = 1, ncel
      w1(iel) = w1(iel)*cp0
    enddo
  endif
!     (CP/CV)*MUT/SIGMAS
  if(icv.gt.0) then
    do iel = 1, ncel
      w1(iel) = w1(iel)/propce(iel,ipproc(icv))
    enddo
  else
    do iel = 1, ncel
      w1(iel) = w1(iel)/cv0
    enddo
  endif
!     (CP/CV)*MUT/SIGMAS+LAMBDA/CV
  if(ipcvsl.eq.0)then
    do iel = 1, ncel
      w1(iel) = w1(iel) + visls0(iscal)
    enddo
  else
    do iel = 1, ncel
      w1(iel) = w1(iel) + propce(iel,ipcvsl)
    enddo
  endif

  call viscfa                                                     &
  !==========
 ( imvisf ,                                                       &
   w1     ,                                                       &
   viscf  , viscb  )


!     TERME DIFFUSIF COMPLEMENTAIRE : - div( K grad ( epsilon - Cv.T ) )
!     =============================                   1  2
!                                     - div( K grad ( -.u  ) )
!                                                     2

!     Terme complementaire au centre des cellules
  iccfth = 7
  imodif = 0
  call cfther                                                     &
  !==========
 ( nvar   ,                                                       &
   iccfth , imodif ,                                              &
   dt     , rtp    , rtpa   , propce ,                            &
   w9     , wb     , w8     , w4     )

!     Calcul de la divergence avec reconstruction


!   Calcul du gradient de (0.5*u*u+EPSILONsup)


do iel = 1, ncel
  w7(iel) =0.5d0*( rtp(iel,iu)**2                          &
                  +rtp(iel,iv)**2                          &
                  +rtp(iel,iw)**2 ) + w9(iel)
enddo

! Rq : A defaut de connaitre les parametres, on prend ceux de la Vitesse

iii = iu
inc = 1
iccocg = 1
nswrgp = nswrgr(iii)
imligp = imligr(iii)
iwarnp = iwarni(iii)
epsrgp = epsrgr(iii)
climgp = climgr(iii)
extrap = extrag(iii)

! Allocate temporary arrays
allocate(coefap(nfabor))
allocate(coefbp(nfabor))

do ifac = 1, nfabor
  coefap(ifac) = zero
  coefbp(ifac) = 1.d0
enddo

! En periodique et parallele, echange avant calcul du gradient
if (irangp.ge.0.or.iperio.eq.1) then
  call synsca(w7)
  !==========
endif

!  IVAR0 = 0 (indique pour la periodicite de rotation que la variable
!     n'est pas la vitesse ni Rij)
ivar0 = 0
call grdcel                                                       &
!==========
 ( ivar0  , imrgra , inc    , iccocg , nswrgp , imligp ,          &
   iwarnp , nfecra , epsrgp , climgp , extrap ,                   &
   w7     , coefap , coefbp ,                                     &
   grad   )

! Free memory
deallocate(coefap, coefbp)

!     Faces internes

do ifac = 1, nfac

  ii = ifacel(1,ifac)
  jj = ifacel(2,ifac)

  dijpfx = dijpf(1,ifac)
  dijpfy = dijpf(2,ifac)
  dijpfz = dijpf(3,ifac)

  pnd   = pond(ifac)

!        Calcul II' et JJ'

  diipfx = cdgfac(1,ifac) - (xyzcen(1,ii) + (1.d0-pnd) * dijpfx)
  diipfy = cdgfac(2,ifac) - (xyzcen(2,ii) + (1.d0-pnd) * dijpfy)
  diipfz = cdgfac(3,ifac) - (xyzcen(3,ii) + (1.d0-pnd) * dijpfz)
  djjpfx = cdgfac(1,ifac) -  xyzcen(1,jj) +  pnd  * dijpfx
  djjpfy = cdgfac(2,ifac) -  xyzcen(2,jj) +  pnd  * dijpfy
  djjpfz = cdgfac(3,ifac) -  xyzcen(3,jj) +  pnd  * dijpfz

  pip = w7(ii) + grad(ii,1)*diipfx+grad(ii,2)*diipfy+grad(ii,3)*diipfz
  pjp = w7(jj) + grad(jj,1)*djjpfx+grad(jj,2)*djjpfy+grad(jj,3)*djjpfz

  flux = viscf(ifac)*(pip-pjp)

  smbrs(ii) = smbrs(ii) - flux
  smbrs(jj) = smbrs(jj) + flux

enddo


!       Assemblage a partir des facettes de bord
!     Pour les faces à flux imposé ou temperature imposée, tout est
!       pris par le terme de diffusion de l'energie. On ne doit donc
!       pas prendre en compte la contribution des termes en u2 et e-CvT
!       quand ifbet(ifac).ne.0

  do ifac = 1, nfabor

    if (ifbet(ifac).eq.0) then

      iel = ifabor(ifac)

      flux = viscb(ifac)*(w1(iel)/distb(ifac))*            &
            ( w9(iel) - wb(ifac)                           &
           + 0.5d0*( rtp(iel,iu)**2 -                      &
   ( coefa(ifac,iclrtp(iu,icoef))                          &
   + coefb(ifac,iclrtp(iu,icoef))*rtp(iel,iu) )**2         &
                   + rtp(iel,iv)**2 -                      &
   ( coefa(ifac,iclrtp(iv,icoef))                          &
   + coefb(ifac,iclrtp(iv,icoef))*rtp(iel,iv) )**2         &
                   + rtp(iel,iw)**2 -                      &
   ( coefa(ifac,iclrtp(iw,icoef))                          &
   + coefb(ifac,iclrtp(iw,icoef))*rtp(iel,iw) )**2))

      smbrs(iel) = smbrs(iel) - flux

    endif

  enddo

else

  do ifac = 1, nfac
    viscf(ifac) = 0.d0
  enddo
  do ifac = 1, nfabor
    viscb(ifac) = 0.d0
  enddo

endif

!===============================================================================
! 4. RESOLUTION
!===============================================================================

iconvp = iconv (ivar)
idiffp = idiff (ivar)
ireslp = iresol(ivar)
nitmap = nitmax(ivar)
ndircp = ndircl(ivar)
nswrsp = nswrsm(ivar)
nswrgp = nswrgr(ivar)
imligp = imligr(ivar)
ircflp = ircflu(ivar)
ischcp = ischcv(ivar)
isstpp = isstpc(ivar)
imgrp  = imgr  (ivar)
ncymxp = ncymax(ivar)
nitmfp = nitmgf(ivar)
ipp    = ippvar
iwarnp = iwarni(ivar)
blencp = blencv(ivar)
epsilp = epsilo(ivar)
epsrsp = epsrsm(ivar)
epsrgp = epsrgr(ivar)
climgp = climgr(ivar)
extrap = extrag(ivar)
thetap = thetav(ivar)
iescap = 0

call cfcdts                                                       &
!==========
 ( nvar   , nscal  ,                                              &
   ivar  , iconvp , idiffp , ireslp , ndircp , nitmap ,           &
   imrgra , nswrsp , nswrgp , imligp , ircflp ,                   &
   ischcp , isstpp , iescap ,                                     &
   imgrp  , ncymxp , nitmfp , ipp    , iwarnp ,                   &
   blencp , epsilp , epsrsp , epsrgp , climgp , extrap , thetap , &
   rtpa(1,ivar)    , coefa(1,iclvar) , coefb(1,iclvar) ,          &
                     coefa(1,iclvaf) , coefb(1,iclvaf) ,          &
                     imasfl , bmasfl ,                            &
   viscf  , viscb  , viscf  , viscb  ,                            &
   rovsdt , smbrs  , rtp(1,ivar)     ,                            &
   rvoid  )

!===============================================================================
! 5. IMPRESSIONS ET CLIPPINGS
!===============================================================================

! Valeur bidon
  iii = 1

call clpsca(ncelet, ncel, iscal, rtp(1,iii), rtp)
!==========

! --- Traitement utilisateur pour gestion plus fine des bornes
!       et actions correctives éventuelles.
  iccfth = -4
  imodif = 0
  call cfther                                                     &
  !==========
 ( nvar   ,                                                       &
   iccfth , imodif ,                                              &
   dt     , rtp    , rtpa   , propce ,                            &
   w6     , w7     , w8     , w9     )


! --- Bilan explicite (voir codits : on enleve l'increment)

if (iwarni(ivar).ge.2) then
  do iel = 1, ncel
    smbrs(iel) = smbrs(iel)                                       &
            - istat(ivar)*(propce(iel,ipcrom)/dt(iel))*volume(iel)&
                *(rtp(iel,ivar)-rtpa(iel,ivar))                   &
                * max(0,min(nswrsm(ivar)-2,1))
  enddo
  isqrt = 1
  call prodsc(ncel,isqrt,smbrs,smbrs,sclnor)
  write(nfecra,1200)chaine(1:8) ,sclnor
endif

!===============================================================================
! 6. ACTUALISATION FINALE DE LA PRESSION (et calcul de la température)
!===============================================================================
!                               n+1      n+1  n+1
! On utilise l'equation d'etat P   =P(RHO   ,H   )

! --- Calcul de P et T au centre des cellules
  iccfth = 24
  imodif = 0
  call cfther                                                     &
  !==========
 ( nvar   ,                                                       &
   iccfth , imodif ,                                              &
   dt     , rtp    , rtpa   , propce ,                            &
   rtp(1,ipr) , rtp(1,isca(itempk)) , w8     , w9 )

!===============================================================================
! 7. COMMUNICATION DE LA PRESSION, DE L'ENERGIE ET DE LA TEMPERATURE
!===============================================================================

if (irangp.ge.0.or.iperio.eq.1) then
  call synsca(rtp(1,ipr))
  !==========
  call synsca(rtp(1,ivar))
  !==========
  call synsca(rtp(1,isca(itempk)))
  !==========
endif

! Free memory
deallocate(wb)
deallocate(grad)
deallocate(w1)
deallocate(w4, w5, w6)
deallocate(w7, w8, w9)

!--------
! FORMATS
!--------

 1000 format(/,                                                   &
'   ** RESOLUTION POUR LA VARIABLE ',A8                        ,/,&
'      ---------------------------                            ',/)
 1200 format(1X,A8,' : BILAN EXPLICITE = ',E14.5)

!----
! FIN
!----

return

end subroutine
