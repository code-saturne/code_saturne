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

subroutine navsto &
!================

 ( nvar   , nscal  , iterns , icvrge ,                            &
   isostd ,                                                       &
   dt     , tpucou , rtp    , rtpa   , propce , propfb ,          &
   tslagr , coefa  , coefb  , frcxt  , prhyd  ,                   &
   trava  , ximpa  , uvwk   )

!===============================================================================
! FONCTION :
! ----------

! RESOLUTION DES EQUATIONS N-S 1 PHASE INCOMPRESSIBLE OU RO VARIABLE
! SUR UN PAS DE TEMPS (CONVECTION/DIFFUSION - PRESSION /CONTINUITE)

!-------------------------------------------------------------------------------
!ARGU                             ARGUMENTS
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! iterns           ! e  ! <-- ! numero d'iteration sur navsto                  !
! icvrge           ! e  ! <-- ! indicateur de convergence du pnt fix           !
! isostd           ! te ! <-- ! indicateur de sortie standard                  !
!    (nfabor+1)    !    !     !  +numero de la face de reference               !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! tpucou(ncelet,3) ! ra ! <-- ! velocity-pressure coupling                     !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
! frcxt(3,ncelet)  ! tr ! <-- ! force exterieure generant la pression          !
!                  !    !     !  hydrostatique                                 !
! prhyd(ncelet)    ! tr ! <-- ! hydrostatic pressure predicted at cell centers !
! tslagr           ! tr ! <-- ! terme de couplage retour du                    !
!(ncelet,*)        !    !     !     lagrangien                                 !
! trava,ximpa      ! tr ! <-- ! tableau de travail pour couplage               !
!(ncelet,3)        !    !     ! vitesse pression par point fixe                !
! uvwk             ! tr ! <-- ! tableau de travail pour couplage u/p           !
!(ncelet,3)        !    !     ! sert a stocker la vitesse de                   !
!                  !    !     ! l'iteration precedente                         !
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
use cstphy
use cstnum
use optcal
use pointe
use albase
use parall
use period
use ppppar
use ppthch
use ppincl
use cplsat
use mesh
use field

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal  , iterns , icvrge

integer          isostd(nfabor+1)

double precision dt(ncelet), tpucou(ncelet,3), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*), propfb(ndimfb,*)
double precision tslagr(ncelet,*)
double precision coefa(ndimfb,*), coefb(ndimfb,*)
double precision frcxt(3,ncelet)
double precision prhyd(ncelet)
double precision trava(ncelet,ndim),ximpa(ncelet,ndim)
double precision uvwk(ncelet,ndim)

! Local variables

integer          iccocg, inc, iel, iel1, iel2, ifac, imax, imaxt
integer          ii    , inod
integer          isou, ivar, iitsm
integer          iclipr, iclipf
integer          icliup, iclivp, icliwp, init
integer          icluma, iclvma, iclwma
integer          iflmas, iflmab, ipcrom, ipbrom
integer          iflmb0
integer          nswrgp, imligp, iwarnp, imaspe , itypfl
integer          nbrval, iappel, iescop, idtsca
integer          ndircp, icpt  , iecrw
integer          numcpl
double precision rnorm , rnormt, rnorma, rnormi, vitnor
double precision dtsrom, unsrom, surf  , rhom
double precision epsrgp, climgp, extrap, xyzmax(3)
double precision thetap, xdu, xdv, xdw
double precision xxp0 , xyp0 , xzp0
double precision rhofac, dtfac, ddepx , ddepy, ddepz
double precision xnrdis
double precision vitbox, vitboy, vitboz

double precision rvoid(1)

double precision, allocatable, dimension(:), target :: viscf, viscb
double precision, allocatable, dimension(:), target :: wvisfi, wvisbi
double precision, allocatable, dimension(:) :: drtp, smbr, rovsdt
double precision, allocatable, dimension(:,:) :: trav, grad
double precision, allocatable, dimension(:) :: w1, w2, w3
double precision, allocatable, dimension(:) :: w4, w5, w6
double precision, allocatable, dimension(:) :: w7, w8, w9
double precision, allocatable, dimension(:) :: w10
double precision, allocatable, dimension(:,:) :: dfrcxt
double precision, allocatable, dimension(:,:) :: frchy, dfrchy
double precision, allocatable, dimension(:,:) :: grdphd
double precision, allocatable, dimension(:) :: esflum, esflub
double precision, allocatable, dimension(:) :: flint, flbrd
double precision, allocatable, dimension(:) :: coefap

double precision, pointer, dimension(:) :: viscfi => null(), viscbi => null()
double precision, dimension(:), pointer :: imasfl, bmasfl

!===============================================================================

!===============================================================================
! 0. Initialization
!===============================================================================

! Allocate temporary arrays for the velocity-pressure resolution
allocate(viscf(nfac), viscb(nfabor))
allocate(drtp(ncelet), smbr(ncelet), rovsdt(ncelet))
allocate(trav(ncelet,3))

! Array for delta p gradient boundary conditions
allocate(coefap(nfabor))

! Allocate other arrays, depending on user options
allocate(dfrcxt(3,ncelet))
if (icalhy.eq.1) allocate(frchy(ncelet,ndim), dfrchy(ncelet,ndim))
if (iphydr.eq.2) allocate(grdphd(ncelet,ndim))
if (iescal(iestot).gt.0) allocate(esflum(nfac), esflub(nfabor))
if (itytur.eq.3.and.irijnu.eq.1) then
  allocate(wvisfi(nfac), wvisbi(nfabor))
  viscfi => wvisfi(1:nfac)
  viscbi => wvisbi(1:nfabor)
else
  viscfi => viscf(1:nfac)
  viscbi => viscb(1:nfabor)
endif

! Allocate work arrays
allocate(w1(ncelet), w2(ncelet), w3(ncelet))
allocate(w4(ncelet), w5(ncelet), w6(ncelet))
allocate(w7(ncelet), w8(ncelet), w9(ncelet))
if (irnpnw.eq.1) allocate(w10(ncelet))

if (iwarni(iu).ge.1) then
  write(nfecra,1000)
endif

! Initialize variables to avoid compiler warnings

ivar = 0
iflmas = 0
ipcrom = 0
imax = 0

! Memoire


if (nterup.gt.1) then

  do isou = 1, 3
    if (isou.eq.1) ivar = iu
    if (isou.eq.2) ivar = iv
    if (isou.eq.3) ivar = iw
    !     La boucle sur NCELET est une securite au cas
    !       ou on utiliserait UVWK par erreur a ITERNS = 1
    !$omp parallel do
    do iel = 1,ncelet
      uvwk(iel,isou) = rtp(iel,ivar)
    enddo
  enddo

  ! Calcul de la norme L2 de la vitesse
  if (iterns.eq.1) then
    xnrmu0 = 0.d0
    !$omp parallel do reduction(+:xnrmu0)
    do iel = 1, ncel
      xnrmu0 =   xnrmu0 +(rtpa(iel,iu)**2 + rtpa(iel,iv)**2 + rtpa(iel,iw)**2) &
               * volume(iel)
    enddo
    if (irangp.ge.0) then
      call parsom (xnrmu0)
      !==========
    endif
    ! En cas de couplage entre deux instances de Code_Saturne, on calcule
    ! la norme totale de la vitesse
    ! Necessaire pour que l'une des instances ne stoppe pas plus tot que les autres
    ! (il faudrait quand meme verifier les options numeriques, ...)
    do numcpl = 1, nbrcpl
      call tbrcpl ( numcpl, 1, 1, xnrmu0, xnrdis )
      !==========
      xnrmu0 = xnrmu0 + xnrdis
    enddo
    xnrmu0 = sqrt(xnrmu0)
  endif

  ! On assure la periodicite ou le parallelisme de UVWK et la pression
  ! (cette derniere vaut la pression a l'iteration precedente)
  if (iterns.gt.1) then
    if (irangp.ge.0.or.iperio.eq.1) then
      call synvec(uvwk(1,1), uvwk(1,2), uvwk(1,3))
      !==========
      call synsca(rtpa(1,ipr))
      !==========
    endif
  endif

endif

!===============================================================================
! 1. Prediction of the mass flux in case of Low Mach compressible algorithm
!===============================================================================

if ((idilat.eq.2.or.idilat.eq.3).and. &
    (ntcabs.gt.1.or.isuite.gt.0)) then

  call predfl(nvar, ncetsm, icetsm, dt, propce, smacel)
  !==========

endif

!===============================================================================
! 2. Hydrostatic pressure prediction in case of Low Mach compressible algorithm
!===============================================================================

if (iphydr.eq.2) then

  call prehyd(propce, prhyd, grdphd)
  !==========

endif

!===============================================================================
! 3.  ETAPE DE PREDICTION DES VITESSES
!===============================================================================

iappel = 1

! Id of the mass flux
call field_get_key_int(ivarfl(iu), kimasf, iflmas)
call field_get_key_int(ivarfl(iu), kbmasf, iflmab)

! Pointers to the mass fluxes
call field_get_val_s(iflmas, imasfl)
call field_get_val_s(iflmab, bmasfl)

call preduv &
!==========
( iappel ,                                                       &
  nvar   , nscal  , iterns ,                                     &
  ncepdc , ncetsm ,                                              &
  icepdc , icetsm , itypsm ,                                     &
  dt     , rtp    , rtpa   , propce , propfb ,                   &
  imasfl , bmasfl ,                                              &
  tslagr , coefa  , coefb  ,                                     &
  ckupdc , smacel , frcxt  , grdphd ,                            &
  trava  , ximpa  , uvwk   , dfrcxt , tpucou ,  trav  ,          &
  viscf  , viscb  , viscfi , viscbi ,                            &
  drtp   , smbr   , rovsdt ,                                     &
  w1     , w7     , w8     , w9     , w10    )

! --- Sortie si pas de pression continuite
!       on met a jour les flux de masse, et on sort

if ( iprco.le.0 ) then

  icliup = iclrtp(iu ,icoef)
  iclivp = iclrtp(iv ,icoef)
  icliwp = iclrtp(iw ,icoef)

  ipcrom = ipproc(irom  )
  ipbrom = ipprob(irom  )

  init   = 1
  inc    = 1
  iccocg = 1
  iflmb0 = 1
  if (iale.eq.1) iflmb0 = 0
  nswrgp = nswrgr(iu)
  imligp = imligr(iu)
  iwarnp = iwarni(iu)
  epsrgp = epsrgr(iu)
  climgp = climgr(iu)
  extrap = extrag(iu)

  imaspe = 1
  itypfl = 1

  call inimas &
  !==========
( iu  , iv  , iw  , imaspe , itypfl ,                            &
  iflmb0 , init   , inc    , imrgra , iccocg , nswrgp , imligp , &
  iwarnp , nfecra ,                                              &
  epsrgp , climgp , extrap ,                                     &
  propce(1,ipcrom), propfb(1,ipbrom),                            &
  rtp(1,iu) , rtp(1,iv) , rtp(1,iw) ,                            &
  coefa(1,icliup), coefa(1,iclivp), coefa(1,icliwp),             &
  coefb(1,icliup), coefb(1,iclivp), coefb(1,icliwp),             &
  imasfl , bmasfl )


!     En ALE on doit rajouter la composante en vitesse de maillage
  if (iale.eq.1) then

    icluma = iclrtp(iuma ,icoef)
    iclvma = iclrtp(ivma ,icoef)
    iclwma = iclrtp(iwma ,icoef)

!     On change de signe car on veut l'oppose de la vitesse de maillage
!       aux faces
    !$omp parallel do
    do iel = 1, ncelet
      rtp(iel,iuma) = -rtp(iel,iuma)
      rtp(iel,ivma) = -rtp(iel,ivma)
      rtp(iel,iwma) = -rtp(iel,iwma)
    enddo
    !$omp parallel do if(nfabor > thr_n_min)
    do ifac = 1, nfabor
      coefa(ifac,icluma) = -coefa(ifac,icluma)
      coefa(ifac,iclvma) = -coefa(ifac,iclvma)
      coefa(ifac,iclwma) = -coefa(ifac,iclwma)
    enddo

!     One temporary array needed for internal faces, in case some internal vertices
!       are moved directly by the user
    allocate(flint(nfac))

    ipcrom = ipproc(irom  )
    ipbrom = ipprob(irom  )

    init   = 0
    inc    = 1
    iccocg = 1
    iflmb0 = 1
    nswrgp = nswrgr(iuma )
    imligp = imligr(iuma )
    iwarnp = iwarni(iuma )
    epsrgp = epsrgr(iuma )
    climgp = climgr(iuma )
    extrap = extrag(iuma )

    imaspe = 1
    itypfl = 1

    !$omp parallel do
    do ifac = 1, nfac
      flint(ifac) = 0.d0
    enddo

    call inimas &
    !==========
 ( iuma   , ivma   , iwma   , imaspe , itypfl ,                   &
   iflmb0 , init   , inc    , imrgra , iccocg , nswrgp , imligp , &
   iwarnp , nfecra ,                                              &
   epsrgp , climgp , extrap ,                                     &
   propce(1,ipcrom), propfb(1,ipbrom),                            &
   rtp(1,iuma )    , rtp(1,ivma )    , rtp(1,iwma)     ,          &
   coefa(1,icluma), coefa(1,iclvma), coefa(1,iclwma),             &
   coefb(1,icluma), coefb(1,iclvma), coefb(1,iclwma),             &
   flint  , bmasfl )

    !$omp parallel do
    do iel = 1, ncelet
      rtp(iel,iuma) = -rtp(iel,iuma)
      rtp(iel,ivma) = -rtp(iel,ivma)
      rtp(iel,iwma) = -rtp(iel,iwma)
    enddo
    !$omp parallel do if(nfabor > thr_n_min)
    do ifac = 1, nfabor
      coefa(ifac,icluma) = -coefa(ifac,icluma)
      coefa(ifac,iclvma) = -coefa(ifac,iclvma)
      coefa(ifac,iclwma) = -coefa(ifac,iclwma)
    enddo

    !$omp parallel do private(ii, inod, iel1, iel2, dtfac, rhofac, iecrw,  &
    !$omp                     icpt, ddepx, ddepy, ddepz)
    do ifac = 1, nfac
      iecrw = 0
      ddepx = 0.d0
      ddepy = 0.d0
      ddepz = 0.d0
      icpt  = 0
      do ii = ipnfac(ifac),ipnfac(ifac+1)-1
        inod = nodfac(ii)
        if (impale(inod).eq.0) iecrw = iecrw + 1
        icpt = icpt + 1
        ddepx = ddepx + depale(1,inod) + xyzno0(1,inod)-xyznod(1,inod)
        ddepy = ddepy + depale(2,inod) + xyzno0(2,inod)-xyznod(2,inod)
        ddepz = ddepz + depale(3,inod) + xyzno0(3,inod)-xyznod(3,inod)
      enddo
      ! If all the face vertices have prescribed displacement, w is evaluated
      ! from this displacement
      if (iecrw.eq.0) then
        iel1 = ifacel(1,ifac)
        iel2 = ifacel(2,ifac)
        dtfac = 0.5d0*(dt(iel1) + dt(iel2))
        rhofac = 0.5d0*(propce(iel1,ipcrom) + propce(iel2,ipcrom))
        imasfl(ifac) = imasfl(ifac) - rhofac*(                   &
             ddepx*surfac(1,ifac) +ddepy*surfac(2,ifac) +ddepz*surfac(3,ifac)) &
             /dtfac/icpt
        !     Else w is calculated from the cell-centre mesh velocity
      else
        imasfl(ifac) = imasfl(ifac) + flint(ifac)
      endif
    enddo

    ! Free memory
    deallocate(flint)

  endif

  ! Ajout de la vitesse du solide dans le flux convectif,
  ! si le maillage est mobile (solide rigide)
  ! En turbomachine, on connaît exactement la vitesse de maillage à ajouter
  if (imobil.eq.1) then

    ipcrom = ipproc(irom  )
    ipbrom = ipprob(irom  )

    !$omp parallel do private(iel1, iel2, dtfac, rhofac, &
    !$omp                     vitbox, vitboy, vitboz)
    do ifac = 1, nfac
      iel1 = ifacel(1,ifac)
      iel2 = ifacel(2,ifac)
      dtfac  = 0.5d0*(dt(iel1) + dt(iel2))
      rhofac = 0.5d0*(propce(iel1,ipcrom) + propce(iel2,ipcrom))
      vitbox = omegay*cdgfac(3,ifac) - omegaz*cdgfac(2,ifac)
      vitboy = omegaz*cdgfac(1,ifac) - omegax*cdgfac(3,ifac)
      vitboz = omegax*cdgfac(2,ifac) - omegay*cdgfac(1,ifac)
      imasfl(ifac) = imasfl(ifac) - rhofac*(                                       &
           vitbox*surfac(1,ifac) + vitboy*surfac(2,ifac) + vitboz*surfac(3,ifac) )
    enddo
    !$omp parallel do private(iel, dtfac, rhofac, &
    !$omp                     vitbox, vitboy, vitboz) if(nfabor > thr_n_min)
    do ifac = 1, nfabor
      iel = ifabor(ifac)
      dtfac  = dt(iel)
      rhofac = propfb(ifac,ipbrom)
      vitbox = omegay*cdgfbo(3,ifac) - omegaz*cdgfbo(2,ifac)
      vitboy = omegaz*cdgfbo(1,ifac) - omegax*cdgfbo(3,ifac)
      vitboz = omegax*cdgfbo(2,ifac) - omegay*cdgfbo(1,ifac)
      bmasfl(ifac) = bmasfl(ifac) - rhofac*(                                       &
           vitbox*surfbo(1,ifac) + vitboy*surfbo(2,ifac) + vitboz*surfbo(3,ifac) )
    enddo

  endif

  return

endif

!===============================================================================
! 4.  ETAPE DE PRESSION/CONTINUITE ( VITESSE/PRESSION )
!===============================================================================

if (iwarni(iu).ge.1) then
  write(nfecra,1200)
endif

! --- Pas de temps scalaire ou pas
idtsca = 0
if ((ipucou.eq.1).or.(ncpdct.gt.0)) idtsca = 1

call resolp &
!==========
 ( nvar   , ncetsm ,                                              &
   icepdc , isostd , idtsca ,                                     &
   dt     , rtp    , rtpa   , propce , propfb ,                   &
   coefa  , coefb  , coefap ,                                     &
   smacel ,                                                       &
   frcxt  , dfrcxt , tpucou , trav   ,                            &
   viscf  , viscb  ,                                              &
   drtp   , smbr   , rovsdt , tslagr ,                            &
   trava  )


!===============================================================================
! 5.  REACTUALISATION DU CHAMP DE VITESSE
!===============================================================================

iclipr = iclrtp(ipr,icoef)
iclipf = iclrtp(ipr,icoeff)
icliup = iclrtp(iu ,icoef)
iclivp = iclrtp(iv ,icoef)
icliwp = iclrtp(iw ,icoef)

ipcrom = ipproc(irom  )
ipbrom = ipprob(irom  )


!       IREVMC = 0 : Methode standard (pas par moindres carres) : on
!                      ajoute un gradient d'increment de pression standard
!                      a la vitesse predite opur obtenir la vitesse corrigee

!       IREVMC = 1 : On applique la methode par moindres carres a
!                      l'ecart entre le flux de masse predit et le flux
!                      de masse actualise,
!                      c'est-a-dire au gradient d'increment de pression
!                    On ajoute la grandeur obtenue aux cellules a la vitesse
!                      predite pour obtenir la vitesse actualisee
!                    Cette methode correspond a IREVMC = 0 avec
!                      gradient par moindres carres IMRGRA=1 dans la
!                      reactualisation des vitesses.

!       IREVMC = 2 : On applique la methode par moindres carres au
!                      flux de masse actualise
!                      pour obtenir la vitesse actualisee
!                    Cette methode correspond a la methode RT0.

!       La methode IREVMC = 2 semble plus "diffusive", mais semble aussi la
!         seule issue pour certains ecoulements atmospheriques de mercure.
!       La methode IREVMC = 1 semble ne pas trop "diffuser", avec un
!         gain du a l'utilisation du gradient moindres carres. Elle
!         se rapproche beaucoup de IREVMC=0.


if (irevmc.eq.1) then

  ! Allocate temporary arrays
  allocate(flint(nfac), flbrd(nfabor))

  !     on ote la partie en u-predit dans le flux de masse final,
  !     on projete des faces vers le centre, puis on rajoute u-predit.
  init   = 1
  inc    = 1
  iccocg = 1
  iflmb0 = 1
  if (iale.eq.1) iflmb0 = 0
  nswrgp = nswrgr(iu)
  imligp = imligr(iu)
  iwarnp = iwarni(iu)
  epsrgp = epsrgr(iu)
  climgp = climgr(iu)
  extrap = extrag(iu)

  imaspe = 1
  itypfl = 1

  call inimas &
  !==========
 ( iu  , iv  , iw  , imaspe , itypfl ,                            &
   iflmb0 , init   , inc    , imrgra , iccocg , nswrgp , imligp , &
   iwarnp , nfecra ,                                              &
   epsrgp , climgp , extrap ,                                     &
   propce(1,ipcrom), propfb(1,ipbrom),                            &
   rtp(1,iu)    , rtp(1,iv)    , rtp(1,iw)    ,                   &
   coefa(1,icliup), coefa(1,iclivp), coefa(1,icliwp),             &
   coefb(1,icliup), coefb(1,iclivp), coefb(1,icliwp),             &
   flint  , flbrd  )

  !$omp parallel do
  do ifac = 1, nfac
    flint(ifac) = imasfl(ifac) - flint(ifac)
  enddo
  !$omp parallel do if(nfabor > thr_n_min)
  do ifac = 1, nfabor
    flbrd(ifac) = bmasfl(ifac) - flbrd(ifac)
  enddo

  call recvmc                                                     &
  !==========
 ( propce(1,ipcrom), flint  , flbrd  ,                            &
   w1     , w2     , w3     ,                                     &
   w4     , w5     , w6     )

  !$omp parallel do
  do iel = 1, ncel
    rtp(iel,iu) = rtp(iel,iu) + w1(iel)
    rtp(iel,iv) = rtp(iel,iv) + w2(iel)
    rtp(iel,iw) = rtp(iel,iw) + w3(iel)
  enddo

  ! Free memory
  deallocate(flint, flbrd)

else if (irevmc.eq.2) then

  call recvmc &
  !==========
( propce(1,ipcrom), imasfl, bmasfl,                              &
  rtp(1,iu), rtp(1,iv), rtp(1,iw),                               &
  w4     , w5     , w6     )


else

  ! Allocate a work array for the gradient calculation
  allocate(grad(ncelet,3))

  !     On corrige la vitesse predite par le gradient cellule de
  !       l'increment de pression

  !     GRADIENT DE L'INCREMENT TOTAL DE PRESSION

  if (idtvar.lt.0) then
    !$omp parallel do
    do iel = 1, ncel
      drtp(iel) = (rtp(iel,ipr) -rtpa(iel,ipr)) / relaxv(ipr)
    enddo
  else
    !$omp parallel do
    do iel = 1, ncel
      drtp(iel) = rtp(iel,ipr) -rtpa(iel,ipr)
    enddo
  endif

  ! --->    TRAITEMENT DU PARALLELISME ET DE LA PERIODICITE

  if (irangp.ge.0.or.iperio.eq.1) then
    call synsca(drtp)
    !==========
  endif


  iccocg = 1
  inc = 0
  if (iphydr.eq.1) inc = 1
  nswrgp = nswrgr(ipr)
  imligp = imligr(ipr)
  iwarnp = iwarni(ipr)
  epsrgp = epsrgr(ipr)
  climgp = climgr(ipr)
  extrap = extrag(ipr)

  call grdpot &
  !==========
 ( ipr    , imrgra , inc    , iccocg , nswrgp , imligp , iphydr , &
   iwarnp , nfecra , epsrgp , climgp , extrap ,                   &
   rvoid  ,                                                       &
   dfrcxt ,                                                       &
   drtp   , coefap        , coefb(1,iclipr)  ,                    &
   grad   )

  !     REACTUALISATION DU CHAMP DE VITESSES

  thetap = thetav(ipr)
  if (iphydr.eq.0.or.iphydr.eq.2) then
    if (idtsca.eq.0) then
      !$omp parallel do private(dtsrom)
      do iel = 1, ncel
        dtsrom = -thetap*dt(iel)/propce(iel,ipcrom)
        rtp(iel,iu) = rtp(iel,iu)+dtsrom*grad(iel,1)
        rtp(iel,iv) = rtp(iel,iv)+dtsrom*grad(iel,2)
        rtp(iel,iw) = rtp(iel,iw)+dtsrom*grad(iel,3)
      enddo
    else
      !$omp parallel do private(unsrom)
      do iel = 1, ncel
        unsrom = -thetap/propce(iel,ipcrom)
        rtp(iel,iu) = rtp(iel,iu) + unsrom*tpucou(iel,1)*grad(iel,1)
        rtp(iel,iv) = rtp(iel,iv) + unsrom*tpucou(iel,2)*grad(iel,2)
        rtp(iel,iw) = rtp(iel,iw) + unsrom*tpucou(iel,3)*grad(iel,3)
      enddo
    endif
  else
    if (idtsca.eq.0) then
      !$omp parallel do private(dtsrom)
      do iel = 1, ncel
        dtsrom = thetap*dt(iel)/propce(iel,ipcrom)
        rtp(iel,iu) = rtp(iel,iu) + dtsrom*(dfrcxt(1 ,iel)-grad(iel,1) )
        rtp(iel,iv) = rtp(iel,iv) + dtsrom*(dfrcxt(2 ,iel)-grad(iel,2) )
        rtp(iel,iw) = rtp(iel,iw) + dtsrom*(dfrcxt(3 ,iel)-grad(iel,3) )
      enddo
    else
      !$omp parallel do private(unsrom)
      do iel = 1, ncel
        unsrom = thetap/propce(iel,ipcrom)
        rtp(iel,iu) = rtp(iel,iu) &
             +unsrom*tpucou(iel,1)*(dfrcxt(1 ,iel)-grad(iel,1))
        rtp(iel,iv) = rtp(iel,iv) &
             +unsrom*tpucou(iel,2)*(dfrcxt(2 ,iel)-grad(iel,2))
        rtp(iel,iw) = rtp(iel,iw) &
             +unsrom*tpucou(iel,3)*(dfrcxt(3 ,iel)-grad(iel,3))
      enddo
    endif
    !     mise a jour des forces exterieures pour le calcul des gradients
    !$omp parallel do
    do iel=1,ncel
      frcxt(1 ,iel) = frcxt(1 ,iel) + dfrcxt(1 ,iel)
      frcxt(2 ,iel) = frcxt(2 ,iel) + dfrcxt(2 ,iel)
      frcxt(3 ,iel) = frcxt(3 ,iel) + dfrcxt(3 ,iel)
    enddo
    if (irangp.ge.0.or.iperio.eq.1) then
      call synvin(frcxt)
      !==========
    endif
    !     mise a jour des Dirichlets de pression en sortie dans COEFA
    iclipr = iclrtp(ipr,icoef)
    iclipf = iclrtp(ipr,icoeff)
    do ifac = 1,nfabor
      if (isostd(ifac).eq.1) then
        coefa(ifac,iclipr) = coefa(ifac,iclipr) + coefap(ifac)
      endif
    enddo
  endif

  ! Free memory
  deallocate(grad)

endif

! In the ALE framework, we add the mesh velocity
if (iale.eq.1) then

  icluma = iclrtp(iuma ,icoef)
  iclvma = iclrtp(ivma ,icoef)
  iclwma = iclrtp(iwma ,icoef)

!     On change de signe car on veut l'oppose de la vitesse de maillage
!       aux faces
  !$omp parallel do
  do iel = 1, ncelet
    rtp(iel,iuma) = -rtp(iel,iuma)
    rtp(iel,ivma) = -rtp(iel,ivma)
    rtp(iel,iwma) = -rtp(iel,iwma)
  enddo
  !$omp parallel do if(nfabor > thr_n_min)
  do ifac = 1, nfabor
    coefa(ifac,icluma) = -coefa(ifac,icluma)
    coefa(ifac,iclvma) = -coefa(ifac,iclvma)
    coefa(ifac,iclwma) = -coefa(ifac,iclwma)
  enddo

  ! One temporary array needed for interior faces, in case some
  ! interior vertices are moved directly by the user
  allocate(flint(nfac))

  ipcrom = ipproc(irom  )
  ipbrom = ipprob(irom  )

  init   = 0
  inc    = 1
  iccocg = 1
  iflmb0 = 1
  nswrgp = nswrgr(iuma )
  imligp = imligr(iuma )
  iwarnp = iwarni(iuma )
  epsrgp = epsrgr(iuma )
  climgp = climgr(iuma )
  extrap = extrag(iuma )

  imaspe = 1
  itypfl = 1

  !$omp parallel do
  do ifac = 1, nfac
    flint(ifac) = 0.d0
  enddo

  call inimas &
  !==========
 ( iu  , iv  , iw  , imaspe , itypfl ,                            &
   iflmb0 , init   , inc    , imrgra , iccocg , nswrgp , imligp , &
   iwarnp , nfecra ,                                              &
   epsrgp , climgp , extrap ,                                     &
   propce(1,ipcrom), propfb(1,ipbrom),                            &
   rtp(1,iuma )    , rtp(1,ivma )    , rtp(1,iwma)     ,          &
   coefa(1,icluma), coefa(1,iclvma), coefa(1,iclwma),             &
   coefb(1,icluma), coefb(1,iclvma), coefb(1,iclwma),             &
   flint  , propfb(1,iflmab) )

  !$omp parallel do
  do iel = 1, ncelet
    rtp(iel,iuma) = -rtp(iel,iuma)
    rtp(iel,ivma) = -rtp(iel,ivma)
    rtp(iel,iwma) = -rtp(iel,iwma)
  enddo
  !$omp parallel do if(nfabor > thr_n_min)
  do ifac = 1, nfabor
    coefa(ifac,icluma) = -coefa(ifac,icluma)
    coefa(ifac,iclvma) = -coefa(ifac,iclvma)
    coefa(ifac,iclwma) = -coefa(ifac,iclwma)
  enddo

  !$omp parallel do private(iecrw, ddepx, ddepy, ddepz, icpt, ii, inod, &
  !$omp                     iel1, iel2, dtfac, rhofac)
  do ifac = 1, nfac
    iecrw = 0
    ddepx = 0.d0
    ddepy = 0.d0
    ddepz = 0.d0
    icpt  = 0
    do ii = ipnfac(ifac),ipnfac(ifac+1)-1
      inod = nodfac(ii)
      if (impale(inod).eq.0) iecrw = iecrw + 1
      icpt = icpt + 1
      ddepx = ddepx + depale(1,inod) + xyzno0(1,inod)-xyznod(1,inod)
      ddepy = ddepy + depale(2,inod) + xyzno0(2,inod)-xyznod(2,inod)
      ddepz = ddepz + depale(3,inod) + xyzno0(3,inod)-xyznod(3,inod)
    enddo
    !     If all the face vertices have imposed displacement, w is evaluated from
    !       this displacement
    if (iecrw.eq.0) then
      iel1 = ifacel(1,ifac)
      iel2 = ifacel(2,ifac)
      dtfac = 0.5d0*(dt(iel1) + dt(iel2))
      rhofac = 0.5d0*(propce(iel1,ipcrom) + propce(iel2,ipcrom))
      imasfl(ifac) = imasfl(ifac) - rhofac*(                      &
           ddepx*surfac(1,ifac)                                   &
           +ddepy*surfac(2,ifac)                                  &
           +ddepz*surfac(3,ifac) )/dtfac/icpt
      !     Else w is calculated from the cell-centre mesh velocity
    else
      imasfl(ifac) = imasfl(ifac) + flint(ifac)
    endif
  enddo

  ! Free memory
  deallocate(flint)

endif

! Ajout de la vitesse du solide dans le flux convectif,
! si le maillage est mobile (solide rigide)
! En turbomachine, on connaît exactement la vitesse de maillage à ajouter
if (imobil.eq.1) then

  ipcrom = ipproc(irom  )
  ipbrom = ipprob(irom  )

  !$omp parallel do private(iel1, iel2, rhofac, vitbox, vitboy, vitboz)
  do ifac = 1, nfac
    iel1 = ifacel(1,ifac)
    iel2 = ifacel(2,ifac)
    rhofac = 0.5d0*(propce(iel1,ipcrom) + propce(iel2,ipcrom))
    vitbox = omegay*cdgfac(3,ifac) - omegaz*cdgfac(2,ifac)
    vitboy = omegaz*cdgfac(1,ifac) - omegax*cdgfac(3,ifac)
    vitboz = omegax*cdgfac(2,ifac) - omegay*cdgfac(1,ifac)
    imasfl(ifac) = imasfl(ifac) - rhofac*(                                       &
         vitbox*surfac(1,ifac) + vitboy*surfac(2,ifac) + vitboz*surfac(3,ifac) )
  enddo
  !$omp parallel do private(iel, rhofac, vitbox, vitboy, vitboz) &
  !$omp          if(nfabor > thr_n_min)
  do ifac = 1, nfabor
    iel = ifabor(ifac)
    rhofac = propfb(ifac,ipbrom)
    vitbox = omegay*cdgfbo(3,ifac) - omegaz*cdgfbo(2,ifac)
    vitboy = omegaz*cdgfbo(1,ifac) - omegax*cdgfbo(3,ifac)
    vitboz = omegax*cdgfbo(2,ifac) - omegay*cdgfbo(1,ifac)
    bmasfl(ifac) = bmasfl(ifac) - rhofac*(                                       &
         vitbox*surfbo(1,ifac) + vitboy*surfbo(2,ifac) + vitboz*surfbo(3,ifac) )
  enddo

endif


!===============================================================================
! 6.  CALCUL D'UN ESTIMATEUR D'ERREUR DE L'ETAPE DE CORRECTION ET TOTAL
!===============================================================================


if (iescal(iescor).gt.0.or.iescal(iestot).gt.0) then

  ! ---> REPERAGE DES VARIABLES

  icliup = iclrtp(iu ,icoef)
  iclivp = iclrtp(iv ,icoef)
  icliwp = iclrtp(iw ,icoef)

  ipcrom = ipproc(irom  )
  ipbrom = ipprob(irom  )


  ! ---> ECHANGE DES VITESSES ET PRESSION EN PERIODICITE ET PARALLELISME

  !    Pour les estimateurs IESCOR et IESTOT, la vitesse doit etre echangee.

  !    Pour l'estimateur IESTOT, la pression doit etre echangee aussi.

  !    Cela ne remplace pas l'echange du debut de pas de temps
  !     a cause de cs_user_extra_operations qui vient plus tard et des calculs suite)


  ! --- Vitesse

  if (irangp.ge.0.or.iperio.eq.1) then
    call synvec(rtp(1,iu), rtp(1,iv), rtp(1,iw))
    !==========
  endif


  !  -- Pression

  if (iescal(iestot).gt.0) then

    if (irangp.ge.0.or.iperio.eq.1) then
      call synsca(rtp(1,ipr))
      !==========
    endif

  endif


  ! ---> CALCUL DU FLUX DE MASSE DEDUIT DE LA VITESSE REACTUALISEE

  init   = 1
  inc    = 1
  iccocg = 1
  iflmb0 = 1
  if (iale.eq.1) iflmb0 = 0
  nswrgp = nswrgr(iu )
  imligp = imligr(iu )
  iwarnp = iwarni(iu )
  epsrgp = epsrgr(iu )
  climgp = climgr(iu )
  extrap = extrag(iu )

  imaspe = 1
  itypfl = 1

  call inimas &
  !==========
 ( iu  , iv  , iw  , imaspe , itypfl ,                            &
   iflmb0 , init   , inc    , imrgra , iccocg , nswrgp , imligp , &
   iwarnp , nfecra ,                                              &
   epsrgp , climgp , extrap ,                                     &
   propce(1,ipcrom), propfb(1,ipbrom),                            &
   rtp(1,iu)    , rtp(1,iv)    , rtp(1,iw)    ,                   &
   coefa(1,icliup), coefa(1,iclivp), coefa(1,icliwp),             &
   coefb(1,icliup), coefb(1,iclivp), coefb(1,icliwp),             &
   esflum , esflub )


  ! ---> CALCUL DE L'ESTIMATEUR CORRECTION : DIVERGENCE DE ROM * U (N + 1)
  !                                          - GAMMA

  if (iescal(iescor).gt.0) then
    init = 1
    call divmas(ncelet,ncel,nfac,nfabor,init,nfecra,            &
                ifacel,ifabor,esflum,esflub,w1)

    if (ncetsm.gt.0) then
      !$omp parallel do private(iel) if(ncetsm > thr_n_min)
      do iitsm = 1, ncetsm
        iel = icetsm(iitsm)
        w1(iel) = w1(iel)-volume(iel)*smacel(iitsm,ipr)
      enddo
    endif

    if (iescal(iescor).eq.2) then
      iescop = ipproc(iestim(iescor))
      !$omp parallel do
      do iel = 1, ncel
        propce(iel,iescop) = abs(w1(iel))
      enddo
    elseif (iescal(iescor).eq.1) then
      iescop = ipproc(iestim(iescor))
      !$omp parallel do
      do iel = 1, ncel
        propce(iel,iescop) = abs(w1(iel)) / volume(iel)
      enddo
    endif
  endif


  ! ---> CALCUL DE L'ESTIMATEUR TOTAL

  if (iescal(iestot).gt.0) then

    !   INITIALISATION DE TRAV AVEC LE TERME INSTATIONNAIRE

    !$omp parallel do
    do iel = 1, ncel
      trav(iel,1) = propce(iel,ipcrom) * volume(iel) *          &
           (rtpa(iel,iu)- rtp(iel,iu))/dt(iel)
      trav(iel,2) = propce(iel,ipcrom) * volume(iel) *          &
           (rtpa(iel,iv)- rtp(iel,iv))/dt(iel)
      trav(iel,3) = propce(iel,ipcrom) * volume(iel) *          &
           (rtpa(iel,iw)- rtp(iel,iw))/dt(iel)
    enddo

    !   APPEL A PREDUV AVEC RTP ET RTP AU LIEU DE RTP ET RTPA
    !                  AVEC LE FLUX DE MASSE RECALCULE
    iappel = 2
    call preduv &
    !==========
 ( iappel ,                                                       &
   nvar   , nscal  , iterns ,                                     &
   ncepdc , ncetsm ,                                              &
   icepdc , icetsm , itypsm ,                                     &
   dt     , rtp    , rtp    , propce , propfb ,                   &
   esflum , esflub ,                                              &
   tslagr , coefa  , coefb  ,                                     &
   ckupdc , smacel , frcxt  , grdphd ,                            &
   trava  , ximpa  , uvwk   , dfrcxt , tpucou , trav   ,          &
   viscf  , viscb  , viscfi , viscbi ,                            &
   drtp   , smbr   , rovsdt ,                                     &
   w1     , w7     , w8     , w9     , w10    )

  endif

endif

!===============================================================================
! 7.  TRAITEMENT DU POINT FIXE SUR LE SYSTEME VITESSE/PRESSION
!===============================================================================

if (nterup.gt.1) then
! TEST DE CONVERGENCE DE L'ALGORITHME ITERATIF
! On initialise ICVRGE a 1 et on le met a 0 si une des phases n'est
! pas convergee

  icvrge = 1

  xnrmu = 0.d0
  !$omp parallel do private(xdu, xdv, xdw) reduction(+:xnrmu)
  do iel = 1, ncel
    xdu = rtp(iel,iu) - uvwk(iel,1)
    xdv = rtp(iel,iv) - uvwk(iel,2)
    xdw = rtp(iel,iw) - uvwk(iel,3)
    xnrmu = xnrmu +(xdu**2 + xdv**2 + xdw**2) * volume(iel)
  enddo
  ! --->    TRAITEMENT DU PARALLELISME

  if (irangp.ge.0) call parsom (xnrmu)
  !==========
  ! -- >    TRAITEMENT DU COUPLAGE ENTRE DEUX INSTANCES DE CODE_SATURNE
  do numcpl = 1, nbrcpl
    call tbrcpl ( numcpl, 1, 1, xnrmu, xnrdis )
    !==========
    xnrmu = xnrmu + xnrdis
  enddo
  xnrmu = sqrt(xnrmu)

  ! Indicateur de convergence du point fixe
  if (xnrmu.ge.epsup*xnrmu0) icvrge = 0

endif

! ---> RECALAGE DE LA PRESSION SUR UNE PRESSION A MOYENNE NULLE
!  On recale si on n'a pas de Dirichlet. Or le nombre de Dirichlets
!  calcule dans typecl.F est NDIRCL si IDIRCL=1 et NDIRCL-1 si IDIRCL=0
!  (ISTAT vaut toujours 0 pour la pression)

if (idircl(ipr).eq.1) then
  ndircp = ndircl(ipr)
else
  ndircp = ndircl(ipr)-1
endif
if (ndircp.le.0) then
  call prmoy0 &
  !==========
( ncelet , ncel   , volume , rtp(1,ipr) )
endif

! Calcul de la pression totale IPRTOT : (definie comme propriete )
! En compressible, la pression resolue est deja la pression totale

if (ippmod(icompf).lt.0) then
  xxp0   = xyzp0(1)
  xyp0   = xyzp0(2)
  xzp0   = xyzp0(3)
  !$omp parallel do
  do iel = 1, ncel
    propce(iel,ipproc(iprtot)) = rtp(iel,ipr)           &
         + ro0*( gx*(xyzcen(1,iel)-xxp0)               &
         + gy*(xyzcen(2,iel)-xyp0)                     &
         + gz*(xyzcen(3,iel)-xzp0) )                   &
         + p0 - pred0
  enddo
endif

!===============================================================================
! 8.  IMPRESSIONS
!===============================================================================

ipcrom = ipproc(irom  )
ipbrom = ipprob(irom  )

if (iwarni(iu).ge.1) then

  write(nfecra,2000)

  rnorm = -1.d0
  !$omp parallel do reduction(max:rnorm)
  do iel = 1, ncel
    rnorm  = max(rnorm,abs(rtp(iel,ipr)))
  enddo
  if (irangp.ge.0) call parmax (rnorm)
  !==========
  write(nfecra,2100)rnorm

  rnorm = -1.d0
  imax = 1
  !$omp parallel private(imaxt, rnormt, vitnor)
  rnormt = -1.d0
  imaxt = 1
  !$omp do
  do iel = 1, ncel
    vitnor = sqrt(rtp(iel,iu)**2+rtp(iel,iv)**2+rtp(iel,iw)**2)
    if (vitnor.ge.rnormt) then
      rnormt = vitnor
      imaxt  = iel
    endif
  enddo
  !$omp critical
  if (rnormt.ge.rnorm) then
    rnorm = rnormt
    imax  = imaxt
  endif
  !$omp end critical
  !$omp end parallel

  xyzmax(1) = xyzcen(1,imax)
  xyzmax(2) = xyzcen(2,imax)
  xyzmax(3) = xyzcen(3,imax)

  if (irangp.ge.0) then
    nbrval = 3
    call parmxl (nbrval, rnorm, xyzmax)
    !==========
  endif

  write(nfecra,2200) rnorm,xyzmax(1),xyzmax(2),xyzmax(3)


  ! Pour la periodicite et le parallelisme, rom est echange dans phyvar

  rnorma = -grand
  rnormi =  grand
  !$omp parallel do reduction(max:rnorma) reduction(min:rnormi) &
  !$omp          private(iel1, iel2, surf, rhom, rnorm)
  do ifac = 1, nfac
    iel1 = ifacel(1,ifac)
    iel2 = ifacel(2,ifac)
    surf = surfan(ifac)
    rhom = (propce(iel1,ipcrom)+propce(iel2,ipcrom))*0.5d0
    rnorm = imasfl(ifac)/(surf*rhom)
    rnorma = max(rnorma,rnorm)
    rnormi = min(rnormi,rnorm)
  enddo
  if (irangp.ge.0) then
    call parmax (rnorma)
    !==========
    call parmin (rnormi)
    !==========
  endif
  write(nfecra,2300)rnorma, rnormi

  rnorma = -grand
  rnormi =  grand
  !$omp parallel do reduction(max:rnorma) reduction(min:rnormi) &
  !$omp          private(rnorm) if(nfabor > thr_n_min)
  do ifac = 1, nfabor
    rnorm = bmasfl(ifac)/(surfbn(ifac)*propfb(ifac,ipbrom))
    rnorma = max(rnorma,rnorm)
    rnormi = min(rnormi,rnorm)
  enddo
  if (irangp.ge.0) then
    call parmax (rnorma)
    !==========
    call parmin (rnormi)
    !==========
  endif
  write(nfecra,2400)rnorma, rnormi

  rnorm = 0.d0
  !$omp parallel do reduction(+:rnorm) if(nfabor > thr_n_min)
  do ifac = 1, nfabor
    rnorm = rnorm + bmasfl(ifac)
  enddo

  if (irangp.ge.0) call parsom (rnorm)
  !==========

  write(nfecra,2500)rnorm

  write(nfecra,2001)

  if (nterup.gt.1) then
    if (icvrge.eq.0) then
      write(nfecra,2600) iterns
      write(nfecra,2601) xnrmu, xnrmu0, epsup
      write(nfecra,2001)
      if (iterns.eq.nterup) then
        write(nfecra,2603)
        write(nfecra,2001)
      endif
    else
      write(nfecra,2602) iterns
      write(nfecra,2601) xnrmu, xnrmu0, epsup
      write(nfecra,2001)
    endif
  endif

endif

! Free memory
deallocate(viscf, viscb)
deallocate(drtp, smbr, rovsdt)
deallocate(trav)
deallocate(dfrcxt)
if (allocated(frchy))  deallocate(frchy, dfrchy)
if (allocated(grdphd)) deallocate(grdphd)
if (allocated(esflum)) deallocate(esflum, esflub)
if (allocated(wvisfi)) deallocate(wvisfi, wvisbi)
deallocate(w1, w2, w3)
deallocate(w4, w5, w6)
deallocate(w7, w8, w9)
if (allocated(w10)) deallocate(w10)
deallocate(coefap)

!--------
! FORMATS
!--------
#if defined(_CS_LANG_FR)

 1000 format(/,                                                   &
'   ** RESOLUTION POUR LA VITESSE                             ',/,&
'      --------------------------                             ',/)
 1200 format(/,                                                   &
'   ** RESOLUTION POUR LA PRESSION CONTINUITE                 ',/,&
'      --------------------------------------                 ',/)
 2000 format(/,' APRES PRESSION CONTINUITE',/,                    &
'-------------------------------------------------------------'  )
 2100 format(                                                           &
' Pression max.',E12.4   ,' (max. de la valeur absolue)       ',/)
 2200 format(                                                           &
' Vitesse  max.',E12.4   ,' en',3E11.3                         ,/)
 2300 format(                                                           &
' Vitesse  en face interne max.',E12.4   ,' ; min.',E12.4        )
 2400 format(                                                           &
' Vitesse  en face de bord max.',E12.4   ,' ; min.',E12.4        )
 2500 format(                                                           &
' Bilan de masse   au bord   ',E14.6                             )
 2600 format(                                                           &
' Informations Point fixe a l''iteration :',I10                ,/)
 2601 format('norme = ',E12.4,' norme 0 = ',E12.4,' toler  = ',E12.4 ,/)
 2602 format(                                                           &
' Convergence du point fixe a l''iteration ',I10               ,/)
 2603 format(                                                           &
' Non convergence du couplage vitesse pression par point fixe  ' )
 2001 format(                                                           &
'-------------------------------------------------------------',/)

#else

 1000 format(/,                                                   &
'   ** SOLVING VELOCITY'                                       ,/,&
'      ----------------'                                       ,/)
 1200 format(/,                                                   &
'   ** SOLVING CONTINUITY PRESSURE'                            ,/,&
'      ---------------------------'                            ,/)
 2000 format(/,' AFTER CONTINUITY PRESSURE',/,                    &
'-------------------------------------------------------------'  )
 2100 format(                                                           &
' Max. pressure',E12.4   ,' (max. absolute value)'             ,/)
 2200 format(                                                           &
' Max. velocity',E12.4   ,' en',3E11.3                         ,/)
 2300 format(                                                           &
' Max. velocity at interior face',E12.4   ,' ; min.',E12.4       )
 2400 format(                                                           &
' Max. velocity at boundary face',E12.4   ,' ; min.',E12.4       )
 2500 format(                                                           &
' Mass balance  at boundary  ',E14.6                             )
 2600 format(                                                           &
' Fixed point informations at iteration:',I10                  ,/)
 2601 format('norm = ',E12.4,' norm 0 = ',E12.4,' toler  = ',E12.4   ,/)
 2602 format(                                                           &
' Fixed point convergence at iteration ',I10                   ,/)
 2603 format(                                                           &
' Non convergence of fixed point for velocity pressure coupling' )
 2001 format(                                                           &
'-------------------------------------------------------------',/)

#endif

!----
! End
!----

return

end subroutine
