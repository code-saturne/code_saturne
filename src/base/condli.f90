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

subroutine condli &
!================

 ( nvar   , nscal  ,                                              &
   isvhb  , isvtb  ,                                              &
   icodcl , isostd ,                                              &
   dt     , rtp    , rtpa   , propce , propfa , propfb , rcodcl , &
   coefa  , coefb  , visvdr , hbord  , thbord , frcxt  )

!===============================================================================
! FONCTION :
! --------

! TRADUCTION DES CONDITIONS AUX LIMITES FOURNIES PAR cs_user_boundary
! SOUS UNE FORME "SIMPLEMENT" ADMISSIBLE PAR LE SOLVEUR

! CETTE TRADUCTION SE PRESENTE SOUS LA FORME D'UNE VALEUR PFAC DE
! LA VARIABLE P CONSIDEREE A LA FACETTE :
!        PFAC = COEFA +COEFB.P(I)
! P(I) : VALEUR DE LA VARIABLE DANS LA CELLULE FLUIDE ADJACENTE

! ATTENTION : SI ON CONSIDERE L'INCREMENT DE LA VARIABLE, LA C.L SE
! REDUIT A : d(PFAC) = COEFB.d(P(I))

! CAS PARTICULIER DES VITESSES :
! --> C.L PEUVENT COUPLER LES 3 COMPOSANTES DE VITESSES
!         (POUR L'INSTANT CE N'EST PAS LE CAS)

!  UXFAC = COEFAX +COEFBX  *UX(I) +COEFU(1)*UY(I) +COEFU(2)*UZ(I)
!  UYFAC = COEFAY +COEFU(1)*UX(I) +COEFBY  *UY(I) +COEFU(3)*UZ(I)
!  UZFAC = COEFAZ +COEFU(2)*UX(I) +COEFU(3)*UY(I) +COEFBZ  *UZ(I)

! On dispose du tableau de tri des faces de bord du
!   pas de temps precedent (sauf au premier pas de temps, ou
!   ITRIFB n'a pas ete renseigne)

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! isvhb            ! e  ! <-- ! indicateur de sauvegarde des                   !
!                  !    !     !  coefficients d'echange aux bords              !
! isvtb            ! e  ! <-- ! indicateur de sauvegarde des                   !
!                  !    !     !  temperatures aux bords                        !
! icodcl           ! te ! --> ! code de condition limites aux faces            !
!  (nfabor,nvar    !    !     !  de bord                                       !
!                  !    !     ! = 1   -> dirichlet                             !
!                  !    !     ! = 3   -> densite de flux                       !
!                  !    !     ! = 4   -> glissemt et u.n=0 (vitesse)           !
!                  !    !     ! = 5   -> frottemt et u.n=0 (vitesse)           !
!                  !    !     ! = 6   -> rugosite et u.n=0 (vitesse)           !
!                  !    !     ! = 9   -> entree/sortie libre (vitesse          !
!                  !    !     !  entrante eventuelle     bloquee               !
! isostd           ! te ! --> ! indicateur de sortie standard                  !
!    (nfabor+1)    !    !     !  +numero de la face de reference               !
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
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
! visvdr(ncelet)   ! tr ! --> ! viscosite dynamique ds les cellules            !
!                  !    !     !  de bord apres amortisst de v driest           !
! hbord            ! tr ! --> ! coefficients d'echange aux bords               !
! (nfabor)         !    !     !                                                !
! thbord           ! tr ! --> ! temperature aux bords en i'                    !
! (nfabor)         !    !     !    (plus exactmt : var. energetique)           !
! frcxt(ncelet,3)  ! tr ! <-- ! force exterieure generant la pression          !
!                  !    !     !  hydrostatique                                 !
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

integer          nvar   , nscal
integer          isvhb  , isvtb

integer          icodcl(nfabor,nvar)
integer          isostd(nfabor+1)

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(ndimfb,*)
double precision rcodcl(nfabor,nvar,3)
double precision frcxt(ncelet,3)
double precision coefa(ndimfb,*), coefb(ndimfb,*)
double precision visvdr(ncelet)
double precision hbord(nfabor),thbord(nfabor)

! Local variables

integer          ifac  , iel   , ivar
integer          isou  , jsou  , ii    , iii
integer          ihcp  , iscal , iscat
integer          inc   , iccocg
integer          iok   , iok1
integer          icodcu
integer          isoent, isorti, ncpt,   isocpt(2)
integer          iclsym, ipatur, ipatrg, isvhbl
integer          ipcvis, ipcvst, ipccp , ipcvsl, ipccv
integer          iclpr , iclu  , iclv  , iclw  , iclk  , iclep
integer          iclnu
integer          icl11 , icl22 , icl33 , icl12 , icl13 , icl23
integer          icluf , iclvf , iclwf , iclphi, iclfb , iclal , iclomg
integer          iclvar, iclvaf, icluma, iclvma, iclwma
integer          iclalp, ipcrom
integer          nswrgp, imligp, iwarnp, icliva

double precision sigma , cpp   , rkl
double precision hint  , hext  , pimp  , xdis
double precision flumbf, visclc, visctc, distbf, srfbn2
double precision epsrgp, climgp, extrap
double precision xxp0, xyp0, xzp0
double precision srfbnf, rnx   , rny   , rnz
double precision upx   , upy   , upz   , vistot
double precision xk, xe, xnu
double precision xllke, xllkmg, xlldrb

logical          ilved

double precision, allocatable, dimension(:) :: w1
double precision, allocatable, dimension(:,:) :: coefu, rijipb
double precision, allocatable, dimension(:,:) :: grad
double precision, dimension(:,:,:), allocatable :: gradv

!===============================================================================

!===============================================================================
! 1.  INITIALISATIONS
!===============================================================================

! Allocate temporary arrays
allocate(coefu(nfabor,3))

!  On a besoin des COEFA et COEFB pour le calcul des gradients
!     pour les cond lim turb en paroi
!   Leur valeur en entree n'est donc pas ecrasee (au premier pas de
!     temps ils sont initialises dans INIVAR a flux nul)

!  COEFU sert a stocker la vitesse en I'
!    On l'utilise aussi pour stocker la pression en I' (dans TYPECL), etc

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
iclvaf = 0
icluf = 0
iclvf = 0
iclwf = 0
ipccv = 0

! Memoire


!  Initialisation du tableau pour stockage de yplus
!     On le remplit dans clptur

if(mod(ipstdv,ipstyp).eq.0) then
  do ifac = 1, nfabor
    yplbr(ifac) = 0.d0
  enddo
endif


!===============================================================================
! 2.  TRAITEMENT DES CONDITIONS DONNES PAR ITYPFB
!===============================================================================


if(ippmod(iphpar).ge.1) then
  call pptycl                                                     &
  !==========
 ( nvar   , nscal  ,                                              &
   icodcl , itrifb , itypfb , izfppp ,                            &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  , rcodcl )
endif

if (iale.eq.1) then
  call altycl                                                     &
  !==========
 ( nvar   , nscal  ,                                              &
   itypfb , ialtyb , icodcl , impale ,                            &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   rcodcl , xyzno0 , depale )
endif

if (imobil.eq.1) then

  call mmtycl &
  !==========
 ( nvar   , nscal  ,                                              &
   itypfb , icodcl ,                                              &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   rcodcl )

endif

call typecl                                                       &
!==========
 ( nvar   , nscal  ,                                              &
   itypfb , itrifb , icodcl , isostd ,                            &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  , rcodcl , frcxt  )

!===============================================================================
! 2.  VERIFICATION DE LA CONSISTANCE DES CL
!===============================================================================

call vericl                                                       &
!==========
 ( nvar   , nscal  ,                                              &
   icodcl ,                                                       &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  , rcodcl )

!===============================================================================
! 4. DISTANCE A LA PAROI ANCIEN MODELE
!===============================================================================
! attention, si on s'en sert pour autre chose, disons le module X
!   bien faire attention dans verini avec quelles options le module
!   X doit alors etre incompatible (perio, parall).

iok1 = 0

if(ineedy.eq.1.and.abs(icdpar).eq.2) then

  ! Allocate a temporary array
  allocate(w1(ncelet))

  ! ON FERA ATTENTION EN PARALLELISME OU PERIODICITE
  !    (UNE PAROI PEUT ETRE PLUS PROCHE EN TRAVERSANT UN BORD ...)

  do iel = 1, ncel
    w1(iel) = grand
  enddo

  do ifac = 1, nfabor
    icodcu = icodcl(ifac,iu)
    if( icodcu.eq.5 .or. icodcu.eq.6 ) then
      do iel = 1, ncel
        xdis =                                                &
             (cdgfbo(1,ifac)-xyzcen(1,iel))**2                   &
             +(cdgfbo(2,ifac)-xyzcen(2,iel))**2                   &
             +(cdgfbo(3,ifac)-xyzcen(3,iel))**2
        if(w1(iel).gt.xdis) then
          w1(iel) = xdis
          ifapat(iel) = ifac
        endif
      enddo
    endif
  enddo

  ! Free memory
  deallocate(w1)

  iok = 0
  do iel = 1, ncel
    if(ifapat(iel).le.0)then
      iok = iok + 1
    endif
  enddo
  if(iok.gt.0) then
    write(nfecra,1000) irijec, idries
    iok1 = 1
  endif

endif

!     Normalement, on ne passe pas en parallele ici,
!       mais au cas ou ...
if(irangp.ge.0) then
  call parcpt(iok1)
endif

if(iok1.ne.0) then
  call csexit (1)
  !==========
endif


!===============================================================================
! 6.  REPERAGE DES VARIABLES
!===============================================================================

! --- Variables
xxp0   = xyzp0(1)
xyp0   = xyzp0(2)
xzp0   = xyzp0(3)

! --- Conditions aux limites
iclpr  = iclrtp(ipr,icoef)
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
  if (iturb.eq.32) iclalp = iclrtp(ial,icoef)
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
ipcvis = ipproc(iviscl)
ipcvst = ipproc(ivisct)
if(icp.gt.0) then
  ipccp  = ipproc(icp   )
else
  ipccp = 0
endif
! --- Compressible
if ( ippmod(icompf).ge.0 ) then
  if(icv.gt.0) then
    ipccv  = ipproc(icv   )
  else
    ipccv = 0
  endif
endif



!===============================================================================
! 6.  CONSTRUCTION DE LA TEMPERATURE OU ENTHALPIE
!        AU CENTRE DES FACES DE BORD (OBTENUS PAR Fi + II'.GRAD(Fi))

!          POUR LE COUPLAGE SYRTHES
!            THBORD EST UTILISE PAR COUPBO EN SORTIE DE CONDLI
!          POUR LE COUPLAGE AVEC LE MODULE THERMIQUE 1D DE PAROI
!            THBORD EST UTILISE PAR COU1DO EN SORTIE DE CONDLI
!          POUR LE RAYONNEMENT
!            THBORD EST DANS LA BOUCLE POUR CONSTUIRE LE FLUX QUI SERT
!            DANS RAYPAR.


!        CECI POURRAIT EN PRATIQUE ETRE HORS DE LA BOUCLE.

!===============================================================================

! Allocate a temporary array for the gradient reconstruction
allocate(grad(ncelet,3))

!  Pour le couplage SYRTHES ou module thermique 1D
!  -----------------------------------------------
!  Ici, on fait une boucle "inutile"  (on ne fait quelque chose
!    que pour ICPSYR(ISCAL) = 1). C'est pour preparer le traitement
!    eventuel de plusieurs temperatures (ie plusieurs couplages
!    SYRTHES a la fois ; noter cependant que meme dans ce cas,
!    une seule temperature sera recue de chaque couplage. En polyph,
!    il faudrait ensuite reconstruire les enthalpies ...
!    plus tard si necessaire).
!  Ici, il ne peut y avoir qu'un seul scalaire avec ICPSYR = 1 et
!    ce uniquement s'il y a effectivement couplage avec SYRTHES
!    (sinon, on s'est arrete dans verini)
! Dans le cas du couplage avec le module 1D, on utilise le scalaire
!    couple avec Syrthes s'il y a couplage, sinon ISCALT(1).
!  La valeur de ISVTB a ete initialisee dans tridim
!    au numero du scalaire couple.


!  Pour le rayonnement
!  -------------------
!  On calcule la valeur en I' s'il y a une variable
!    thermique


!  On recherche l'unique scalaire qui convient
!     (ce peut etre T, H, ou E (en compressible))

iscat = 0

!     Si un scalaire est couple a SYRTHES ou au module 1D
if(isvtb.ne.0) then
  !         si ce n'est pas la variable thermique, ca ne va pas.
  if(isvtb.ne.iscalt) then
    write(nfecra,8000)isvtb,iscalt
    call csexit (1)
    !==========
    !         sinon, on calcule le gradient.
  else
    iscat = isvtb
  endif
endif


!     S'il y a du rayonnement
!       (il y a forcement une variable energetique)
!       on en calcule le gradient
if(iirayo.ge.1) then
  iscat = iscalt
endif

!     S'il y a un scalaire dont il faut calculer le gradient
!       ... on le calcule.
if (iscat .gt. 0) then

  ivar   = isca(iscat)

  if (ntcabs.gt.1 .and. itbrrb.eq.1) then

    inc = 1
    iccocg = 1
    nswrgp = nswrgr(ivar)
    imligp = imligr(ivar)
    iwarnp = iwarni(ivar)
    epsrgp = epsrgr(ivar)
    climgp = climgr(ivar)
    extrap = extrag(ivar)
    icliva = iclrtp(ivar,icoef)

    call grdcel                                                 &
    !==========
 ( ivar   , imrgra , inc    , iccocg , nswrgp , imligp ,          &
   iwarnp , nfecra ,                                              &
   epsrgp , climgp , extrap ,                                     &
   rtpa(1,ivar)    , coefa(1,icliva) , coefb(1,icliva) ,          &
   grad   )

    do ifac = 1 , nfabor
      iel = ifabor(ifac)
      thbord(ifac) = rtpa(iel,ivar) &
                     + grad(iel,1)*diipb(1,ifac) &
                     + grad(iel,2)*diipb(2,ifac) &
                     + grad(iel,3)*diipb(3,ifac)
    enddo

  else

    do ifac = 1 , nfabor
      iel = ifabor(ifac)
      thbord(ifac) = rtpa(iel,ivar)
    enddo

  endif

endif

!===============================================================================
! 6.  CONSTRUCTION DE LA VITESSE ET DU TENSEUR DE REYNOLDS
!        AU CENTRE DES FACES DE BORD (OBTENUS PAR Fi + II'.GRAD(Fi))
!        S'IL Y A DES SYMETRIES OU DES PAROIS TURBULENTES
!===============================================================================

! ---> INDICATEUR SYMETRIES OU PAROIS TURBULENTES

iclsym = 0
ipatur = 0
ipatrg = 0
do ifac = 1, nfabor
  if ( icodcl(ifac,iu).eq.4 ) then
    iclsym = 1
  elseif ( icodcl(ifac,iu).eq.5 ) then
    ipatur = 1
  elseif ( icodcl(ifac,iu).eq.6 ) then
    ipatrg = 1
  endif
  if (iclsym.ne.0.and.ipatur.ne.0.and.ipatrg.ne.0 ) goto 100
enddo
100 continue

if (irangp.ge.0) then
  call parcmx(iclsym)
  call parcmx(ipatur)
  call parcmx(ipatrg)
endif


! ---> CONSTRUCTION DE LA VITESSE AU CENTRE DES FACES DE BORD

if (iclsym.ne.0.or.ipatur.ne.0.or.ipatrg.ne.0) then


  if(ntcabs.gt.1) then

    ! Allocate a temporary array
    allocate(gradv(ncelet,3,3))

    iccocg = 1
    inc    = 1
    nswrgp = nswrgr(iu)
    imligp = imligr(iu)
    iwarnp = iwarni(iu)
    epsrgp = epsrgr(iu)
    climgp = climgr(iu)
    extrap = extrag(iu)
    icliva = iclrtp(iu,icoef)

    if (ivelco.eq.1) then

      ilved = .false.

      call grdvec &
      !==========
    ( iu     , imrgra , inc    , nswrgp , imligp ,                   &
      iwarnp , nfecra ,                                              &
      epsrgp , climgp , extrap ,                                     &
      ilved ,                                                        &
      rtpa(1,iu) ,  coefau , coefbu,                                 &
      gradv  )

    else

      call grdvni &
      !==========
    ( iu     , imrgra , inc    , iccocg , nswrgp , imligp ,          &
      iwarnp , nfecra ,                                              &
      epsrgp , climgp , extrap ,                                     &
      rtpa(1,iu)      , coefa(1,icliva) , coefb(1,icliva) ,          &
      gradv  )

    endif

    do isou = 1, 3
      if(isou.eq.1) ivar = iu
      if(isou.eq.2) ivar = iv
      if(isou.eq.3) ivar = iw

      do ifac = 1, nfabor
        iel = ifabor(ifac)
        coefu(ifac,isou) = gradv(iel,1,isou)*diipb(1,ifac)    &
                         + gradv(iel,2,isou)*diipb(2,ifac)    &
                         + gradv(iel,3,isou)*diipb(3,ifac)    &
                         + rtpa(iel,ivar)
      enddo
    enddo

    deallocate(gradv)

  else

    do isou = 1, 3
      if(isou.eq.1) ivar = iu
      if(isou.eq.2) ivar = iv
      if(isou.eq.3) ivar = iw

      do ifac = 1, nfabor
        iel = ifabor(ifac)
        coefu(ifac,isou) = rtpa(iel,ivar)
      enddo

    enddo

  endif

endif


! ---> CONSTRUCTION DU TENSEUR DE REYNOLDS AU CENTRE DES FACES DE BORD

if ((iclsym.ne.0.or.ipatur.ne.0.or.ipatrg.ne.0)                 &
     .and.itytur.eq.3) then

  ! Allocate a work array to store Rij values at boundary faces
  allocate(rijipb(nfabor,6))

  do isou = 1 , 6

    if(isou.eq.1) ivar = ir11
    if(isou.eq.2) ivar = ir22
    if(isou.eq.3) ivar = ir33
    if(isou.eq.4) ivar = ir12
    if(isou.eq.5) ivar = ir13
    if(isou.eq.6) ivar = ir23


    if(ntcabs.gt.1.and.irijrb.eq.1) then

      ! CALCUL DU GRADIENT CELLULE DE Rij EN I

      inc = 1
      iccocg = 1
      nswrgp = nswrgr(ivar)
      imligp = imligr(ivar)
      iwarnp = iwarni(ivar)
      epsrgp = epsrgr(ivar)
      climgp = climgr(ivar)
      extrap = extrag(ivar)
      icliva = iclrtp(ivar,icoef)

      call grdcel                                               &
      !==========
 ( ivar   , imrgra , inc    , iccocg , nswrgp , imligp ,          &
   iwarnp , nfecra ,                                              &
   epsrgp , climgp , extrap ,                                     &
   rtpa(1,ivar)    , coefa(1,icliva) , coefb(1,icliva) ,          &
   grad   )


      ! CALCUL DE LA VALEUR EN I' DE Rij

      do ifac = 1 , nfabor
        iel = ifabor(ifac)
        rijipb(ifac,isou) = rtpa(iel,ivar) &
                          + grad(iel,1)*diipb(1,ifac) &
                          + grad(iel,2)*diipb(2,ifac) &
                          + grad(iel,3)*diipb(3,ifac)
      enddo


      !   AU PREMIER PAS DE TEMPS, ON NE CONNAIT PAS COEFA ET COEFB
      !   (ILS SONT ANNULES DANS CONDLI), LE CALCUL DE RI' EST SIMPLIFIE

    else

      do ifac = 1 , nfabor
        iel = ifabor(ifac)
        rijipb(ifac,isou) = rtpa(iel,ivar)
      enddo

    endif

  enddo

endif

! Free memory
deallocate(grad)

!===============================================================================
! 7.  TURBULENCE EN PAROI : TOUTES LES VARIABLES CONCERNEES
!       (U,V,W,K,EPSILON,RIJ,TEMPERATURE)
!===============================================================================
! --- On a besoin de COEFU et de RIJIPB (et THBORD pour le rayonnement)

!     On initialise VISVDR a -999.D0.
!     Dans clptur, on amortit la viscosite turbulente sur les cellules
!     de paroi si on a active van Driest. La valeur finale est aussi
!     stockee dans VISVDR.
!     Plus loin, dans vandri, la viscosite sur les cellules
!     de paroi sera amortie une seconde fois. On se sert alors de
!     VISVDR pour lui redonner une valeur correcte.
if(itytur.eq.4.and.idries.eq.1) then
  do iel=1,ncel
    visvdr(iel) = -999.d0
  enddo
endif

if (ipatur.ne.0) then

  ! Smooth wall laws
  call clptur                                                   &
  !==========
 ( nvar   , nscal  ,                                              &
   isvhb  ,                                                       &
   icodcl ,                                                       &
   dt     , rtp    , rtpa   , propce , propfa , propfb , rcodcl , &
   coefu  , rijipb , coefa  , coefb  , visvdr ,                   &
   hbord  , thbord )

endif

if (ipatrg.ne.0) then

  ! Rough wall laws
  call clptrg                                                   &
  !==========
 ( nvar   , nscal  ,                                              &
   isvhb  ,                                                       &
   icodcl ,                                                       &
   dt     , rtp    , rtpa   , propce , propfa , propfb , rcodcl , &
   coefu  , rijipb , coefa  , coefb  , visvdr ,                   &
   hbord  , thbord )

endif

!===============================================================================
! 7.  SYMETRIES POUR LES VECTEURS ET TENSEURS
!       (U,V,W,RIJ)
!===============================================================================
!   On a besoin de COEFU et de RIJIPB

do ifac = 1, nfabor
  isympa(ifac) = 1
enddo

if (iclsym.ne.0) then

  call clsyvt                                                   &
  !==========
 ( nvar   , nscal  ,                                              &
   icodcl ,                                                       &
   dt     , rtp    , rtpa   , propce , propfa , propfb , rcodcl , &
   coefu  , rijipb , coefa  , coefb  )

endif

!===============================================================================
! 8.  VITESSE : SORTIE, DIRICHLET, NEUMANN
!===============================================================================

! ---> SORTIE : SI FLUX ENTRANT, ON "BLOQUE" A L'INFINI AVAL

isoent = 0
isorti = 0
do ifac = 1, nfabor

  flumbf = propfb(ifac,ipprob(ifluma(iu)))

  if( icodcl(ifac,iu).eq.9 ) then

    isorti = isorti + 1
    if( flumbf.lt.-epzero) then
      coefa(ifac,iclu) = 0.d0
      coefb(ifac,iclu) = 0.d0
      coefa(ifac,iclv) = 0.d0
      coefb(ifac,iclv) = 0.d0
      coefa(ifac,iclw) = 0.d0
      coefb(ifac,iclw) = 0.d0
      isoent = isoent + 1
    else
      coefa(ifac,iclu) = 0.d0
      coefb(ifac,iclu) = 1.d0
      coefa(ifac,iclv) = 0.d0
      coefb(ifac,iclv) = 1.d0
      coefa(ifac,iclw) = 0.d0
      coefb(ifac,iclw) = 1.d0
    endif

    ! Coupled solving of the velocity components
    if (ivelco.eq.1) then
      coefau(1,ifac) = coefa(ifac,iclu)
      coefau(2,ifac) = coefa(ifac,iclv)
      coefau(3,ifac) = coefa(ifac,iclw)

      coefbu(1,1,ifac) = coefb(ifac,iclu)
      coefbu(2,2,ifac) = coefb(ifac,iclu)
      coefbu(3,3,ifac) = coefb(ifac,iclu)

      coefbu(1,2,ifac) = 0.d0
      coefbu(1,3,ifac) = 0.d0
      coefbu(2,1,ifac) = 0.d0
      coefbu(2,3,ifac) = 0.d0
      coefbu(3,1,ifac) = 0.d0
      coefbu(3,2,ifac) = 0.d0
    endif


  endif

enddo

if (mod(ntcabs,ntlist).eq.0 .or. iwarni(iu).ge. 0) then
  isocpt(1) = isoent
  isocpt(2) = isorti
  if (irangp.ge.0) then
    ncpt = 2
    call parism(ncpt, isocpt)
  endif
  if (isocpt(2).gt.0 .and. (iwarni(iu).ge.2.or.isocpt(1).gt.0)) then
    write(nfecra,3010) isocpt(1), isocpt(2)
  endif
endif

! ---> DIRICHLET ET FLUX

do isou = 1, 3

  if(isou.eq.1) then
    ivar   = iu
    iclvar = iclu
  elseif(isou.eq.2) then
    ivar   = iv
    iclvar = iclv
  elseif(isou.eq.3) then
    ivar   = iw
    iclvar = iclw
  endif

  do ifac = 1, nfabor

    iel = ifabor(ifac)

    ! --- Proprietes physiques
    visclc = propce(iel,ipcvis)
    visctc = propce(iel,ipcvst)

    ! --- Grandeurs geometriques
    distbf = distb(ifac)

    if (itytur.eq.3) then
      hint =   visclc         /distbf
    else
      hint = ( visclc+visctc )/distbf
    endif

    !      C.L DE TYPE DIRICHLET
    if( icodcl(ifac,ivar).eq.1 ) then
      hext = rcodcl(ifac,ivar,2)
      if(abs(hext).gt.rinfin*0.5d0) then
        pimp = rcodcl(ifac,ivar,1)
        coefa(ifac,iclvar) = pimp
        coefb(ifac,iclvar) = 0.d0
        ! Coupled solving of the velocity components
        if (ivelco.eq.1) then
          coefau(isou,ifac) = pimp
          coefbu(isou,1,ifac) = 0.d0
          coefbu(isou,2,ifac) = 0.d0
          coefbu(isou,3,ifac) = 0.d0
        endif
      else
        pimp = rcodcl(ifac,ivar,1)
        coefa(ifac,iclvar) = hext*pimp/(hint +hext)
        coefb(ifac,iclvar) = hint     /(hint +hext)
        ! Coupled solving of the velocity components
        if (ivelco.eq.1) then
          coefau(isou,ifac) = hext*pimp/(hint +hext)
          do jsou = 1, 3
            if (jsou.eq.isou) then
              coefbu(isou,jsou,ifac) = hint/(hint +hext)
            else
              coefbu(isou,jsou,ifac) = 0.d0
            endif
          enddo
        endif
      endif

    !      C.L DE TYPE FLUX
    elseif( icodcl(ifac,ivar).eq.3 ) then
      coefa(ifac,iclvar) = -rcodcl(ifac,ivar,3)/hint
      coefb(ifac,iclvar) = 1.d0
      ! Coupled solving of the velocity components
      if (ivelco.eq.1) then
        coefau(isou,ifac) = -rcodcl(ifac,ivar,3)/hint
        do jsou = 1, 3
          if(jsou.eq.isou) then
            coefbu(isou,jsou,ifac) = 1.d0
          else
            coefbu(isou,jsou,ifac) = 0.d0
          endif
        enddo
      endif

    endif

  enddo

enddo

! ---> COEFAF ET COEFBF
!       POUR TOUS LES CODES SAUF 4, 5 ET 6 TRAITES SEPAREMENT

do isou = 1, 3

  if(isou.eq.1) then
    ivar   = iu
    iclvar = iclu
    iclvaf = icluf
  elseif(isou.eq.2) then
    ivar   = iv
    iclvar = iclv
    iclvaf = iclvf
  elseif(isou.eq.3) then
    ivar   = iw
    iclvar = iclw
    iclvaf = iclwf
  endif

  if (iclvaf.ne.iclvar) then
    do ifac = 1, nfabor
      if( icodcl(ifac,ivar).eq.1.or.icodcl(ifac,ivar).eq.3.or.  &
          icodcl(ifac,ivar).eq.9                          ) then
        coefa(ifac,iclvaf) = coefa(ifac,iclvar)
        coefb(ifac,iclvaf) = coefb(ifac,iclvar)
      endif
    enddo
  endif

  ! Coupled solving of the velocity components
  if (ivelco.eq.1) then
    do ifac = 1, nfabor
      if( icodcl(ifac,ivar).eq.1.or.icodcl(ifac,ivar).eq.3.or.  &
          icodcl(ifac,ivar).eq.9                          ) then
        cofafu(isou,ifac) = coefau(isou,ifac)
        do jsou = 1, 3
          cofbfu(isou,jsou,ifac) = coefbu(isou,jsou,ifac)
        enddo
      endif
    enddo
  endif

enddo


!===============================================================================
! 9.  PRESSION : DIRICHLET, NEUMANN
!===============================================================================

do ifac = 1, nfabor

  iel = ifabor(ifac)

  ! --- Grandeurs geometriques
  distbf = distb(ifac)

  ! ON MET UN FLUX EN DT.GRAD P (W/m2) DANS cs_user_boundary
  hint = dt(iel)/distbf

  ! On doit remodifier la valeur du  Dirichlet de pression de manière
  !  à retrouver P*. Car dans typecl.f90 on a travaillé avec la pression
  ! totale fournie par l'utilisateur :  Ptotale= P*+ rho.g.r
  ! En compressible, on laisse RCODCL tel quel

  !      C.L DE TYPE DIRICHLET AVEC OU SANS COEFFICIENT D'ECHANGE

  if( icodcl(ifac,ipr).eq.1 ) then
    hext = rcodcl(ifac,ipr,2)
    if ( ippmod(icompf).ge.0 ) then
      pimp = rcodcl(ifac,ipr,1)
    else
      pimp = rcodcl(ifac,ipr,1)                              &
           - ro0*( gx*(cdgfbo(1,ifac)-xxp0)                  &
           + gy*(cdgfbo(2,ifac)-xyp0)                  &
           + gz*(cdgfbo(3,ifac)-xzp0) )                &
           + pred0 - p0
    endif
    if( abs(hext).gt.rinfin*0.5d0 ) then
      coefa(ifac,iclpr) = pimp
      coefb(ifac,iclpr) = 0.d0
    else
      coefa(ifac,iclpr) = hext*pimp/(hint +hext)
      coefb(ifac,iclpr) = hint     /(hint +hext)
    endif
  endif

  !      C.L DE TYPE FLUX
  if( icodcl(ifac,ipr).eq.3 ) then
    coefa(ifac,iclpr) = -rcodcl(ifac,ipr,3)/hint
    coefb(ifac,iclpr) = 1.d0
  endif

enddo


!===============================================================================
! 10.  K, EPSILON, RIJ, V2F, OMEGA : DIRICHLET, NEUMANN
!===============================================================================

! ---> K-EPSILON ET K-OMEGA

if(itytur.eq.2 .or. iturb.eq.60) then

  do ii = 1, 2

    !     Pour le k-omega, on met les valeurs sigma_k2 et sigma_w2 car ce terme
    !     ne concerne en pratique que les entrees (pas de pb en paroi ou en flux
    !     nul)
    if(ii.eq.1 .and. itytur.eq.2) then
      ivar   = ik
      iclvar = iclk
      sigma  = sigmak
    elseif(ii.eq.1 .and. iturb.eq.60) then
      ivar   = ik
      iclvar = iclk
      sigma  = ckwsk2
    elseif (itytur.eq.2) then
      ivar   = iep
      iclvar = iclep
      sigma  = sigmae
    else
      ivar   = iomg
      iclvar = iclomg
      sigma  = ckwsw2
    endif

    do ifac = 1, nfabor

      iel = ifabor(ifac)

      ! --- Proprietes physiques
      visclc = propce(iel,ipcvis)
      visctc = propce(iel,ipcvst)
      flumbf = propfb(ifac,ipprob(ifluma(ik)))

      ! --- Grandeurs geometriques
      distbf = distb(ifac)

      hint = (visclc+visctc/sigma)/distbf

      !      C.L DE TYPE DIRICHLET AVEC OU SANS COEFFICIENT D'ECHANGE
      if(icodcl(ifac,ivar).eq.1) then
        hext = rcodcl(ifac,ivar,2)
        pimp = rcodcl(ifac,ivar,1)
        coefa(ifac,iclvar) = hext*pimp/(hint +hext)
        coefb(ifac,iclvar) = hint     /(hint +hext)
        !      C.L DE TYPE FLUX
      elseif(icodcl(ifac,ivar).eq.3)then
        coefa(ifac,iclvar) = -rcodcl(ifac,ivar,3)/hint
        coefb(ifac,iclvar) = 1.d0
      endif
    enddo

  enddo

  ! ---> RIJ-EPSILON
  !         (ATTENTION, PAS DE VISCT)

elseif(itytur.eq.3) then

  !   --> RIJ

  do isou = 1, 6

    if(isou.eq.1) then
      ivar   = ir11
      iclvar = icl11
    elseif(isou.eq.2) then
      ivar   = ir22
      iclvar = icl22
    elseif(isou.eq.3) then
      ivar   = ir33
      iclvar = icl33
    elseif(isou.eq.4) then
      ivar   = ir12
      iclvar = icl12
    elseif(isou.eq.5) then
      ivar   = ir13
      iclvar = icl13
    elseif(isou.eq.6) then
      ivar   = ir23
      iclvar = icl23
    endif

    do ifac = 1, nfabor

      iel = ifabor(ifac)

      ! --- Proprietes physiques
      visclc = propce(iel,ipcvis)
      flumbf = propfb(ifac,ipprob(ifluma(ir11)))

      ! --- Grandeurs geometriques
      distbf = distb(ifac)

      if(icodcl(ifac,ivar).eq.1) then
        hint = visclc/distbf

        hext = rcodcl(ifac,ivar,2)
        pimp = rcodcl(ifac,ivar,1)
        coefa(ifac,iclvar) = hext*pimp/(hint +hext)
        coefb(ifac,iclvar) = hint     /(hint +hext)

      elseif(icodcl(ifac,ivar).eq.3)then

        hint = visclc/distbf

        coefa(ifac,iclvar) = -rcodcl(ifac,ivar,3)/hint
        coefb(ifac,iclvar) = 1.d0

      endif

    enddo

  enddo


  !   --> EPSILON

  ivar   = iep
  iclvar = iclep

  do ifac = 1, nfabor

    iel = ifabor(ifac)

    ! --- Proprietes physiques
    visclc = propce(iel,ipcvis)
    flumbf = propfb(ifac,ipprob(ifluma(iep)))

    ! --- Grandeurs geometriques
    distbf = distb(ifac)

    hint = visclc/distbf

    !      C.L DE TYPE DIRICHLET AVEC OU SANS COEFFICIENT D'ECHANGE
    if( icodcl(ifac,ivar).eq.1) then
      hext = rcodcl(ifac,ivar,2)
      pimp = rcodcl(ifac,ivar,1)
      coefa(ifac,iclvar) = hext*pimp/(hint +hext)
      coefb(ifac,iclvar) = hint     /(hint +hext)
      !      C.L DE TYPE FLUX
    elseif(                                                     &
         icodcl(ifac,ivar).eq.3)then
      coefa(ifac,iclvar) = -rcodcl(ifac,ivar,3)/hint
      coefb(ifac,iclvar) = 1.d0
    endif

  enddo

  !   --> ALPHA for the EBRSM

  if(iturb.eq.32)then
    ivar   = ial
    iclvar = iclalp

    ipcrom = ipproc(irom)

    do ifac = 1, nfabor

      iel = ifabor(ifac)

      visclc = propce(iel,ipcvis)
      flumbf = propfb(ifac,ipprob(ifluma(ivar)))

      distbf = distb(ifac)

      xk = 0.5d0*(rtpa(iel,ir11)+rtpa(iel,ir22)+rtpa(iel,ir33))

      xnu  = visclc/propce(iel,ipcrom)
      ! Echelle de longueur integrale
      xllke = xk**(3.d0/2.d0)/rtpa(iel,iep)
      ! Echelle de longueur de Kolmogorov
      xllkmg = xceta*(xnu**3/rtpa(iel,iep))**(0.25d0)
      ! Echelle de longueur de Durbin
      xlldrb = xcl*max(xllke,xllkmg)

      hint = (xlldrb**2)/distbf

      ! Dirichlet boundary condition
      if (icodcl(ifac,ivar).eq.1) then
        hext = rcodcl(ifac,ivar,2)
        pimp = rcodcl(ifac,ivar,1)
        coefa(ifac,iclvar) = hext*pimp/(hint +hext)
        coefb(ifac,iclvar) = hint     /(hint +hext)
      ! Neumann boundary condition
      elseif (icodcl(ifac,ivar).eq.3) then
        coefa(ifac,iclvar) = -rcodcl(ifac,ivar,3)/hint
        coefb(ifac,iclvar) = 1.d0
      endif

    enddo
  endif


! ---> MODELES TYPE V2F (PHI_BAR et BL-V2/K)

elseif(itytur.eq.5) then

  !   --> K, EPSILON ET PHI
  do ii = 1, 3

    if(ii.eq.1) then
      ivar   = ik
      iclvar = iclk
      sigma  = sigmak
    elseif(ii.eq.2) then
      ivar   = iep
      iclvar = iclep
      sigma  = sigmae
    else
      ivar   = iphi
      iclvar = iclphi
      sigma  = sigmak
    endif

    do ifac = 1, nfabor

      iel = ifabor(ifac)

      ! --- Proprietes physiques
      visclc = propce(iel,ipcvis)
      visctc = propce(iel,ipcvst)
      flumbf = propfb(ifac,ipprob(ifluma(ik)))

      ! --- Grandeurs geometriques
      distbf = distb(ifac)

      hint = (visclc+visctc/sigma)/distbf

      !      C.L DE TYPE DIRICHLET AVEC OU SANS COEFFICIENT D'ECHANGE
      if(icodcl(ifac,ivar).eq.1) then
        hext = rcodcl(ifac,ivar,2)
        pimp = rcodcl(ifac,ivar,1)
        coefa(ifac,iclvar) = hext*pimp/(hint +hext)
        coefb(ifac,iclvar) = hint     /(hint +hext)
        !      C.L DE TYPE FLUX
      elseif(icodcl(ifac,ivar).eq.3)then
        coefa(ifac,iclvar) = -rcodcl(ifac,ivar,3)/hint
        coefb(ifac,iclvar) = 1.d0
      endif

    enddo

  enddo

  if(iturb.eq.50) then

    !   --> FB

    ivar   = ifb
    iclvar = iclfb

    do ifac = 1, nfabor

      ! --- Proprietes physiques
      visclc = 1.d0
      flumbf = propfb(ifac,ipprob(ifluma(ifb)))

      ! --- Grandeurs geometriques
      distbf = distb(ifac)

      hint = visclc/distbf

      !      C.L DE TYPE DIRICHLET AVEC OU SANS COEFFICIENT D'ECHANGE
      if( icodcl(ifac,ivar).eq.1) then
        hext = rcodcl(ifac,ivar,2)
        pimp = rcodcl(ifac,ivar,1)
        coefa(ifac,iclvar) = hext*pimp/(hint +hext)
        coefb(ifac,iclvar) = hint     /(hint +hext)
        !      C.L DE TYPE FLUX
      elseif(icodcl(ifac,ivar).eq.3)then
        coefa(ifac,iclvar) = -rcodcl(ifac,ivar,3)/hint
        coefb(ifac,iclvar) = 1.d0
      endif

    enddo

  elseif(iturb.eq.51) then

    !   --> ALPHA

    ivar   = ial
    iclvar = iclal

    do ifac = 1, nfabor

! --- Proprietes physiques
      visclc = 1.d0
      flumbf = propfb(ifac,ipprob(ifluma(ial)))

! --- Grandeurs geometriques
      distbf = distb(ifac)

      hint = visclc/distbf

!      C.L DE TYPE DIRICHLET AVEC OU SANS COEFFICIENT D'ECHANGE
      if( icodcl(ifac,ivar).eq.1) then
        hext = rcodcl(ifac,ivar,2)
        pimp = rcodcl(ifac,ivar,1)
        coefa(ifac,iclvar) = hext*pimp/(hint +hext)
        coefb(ifac,iclvar) = hint     /(hint +hext)
!      C.L DE TYPE FLUX
      elseif(icodcl(ifac,ivar).eq.3)then
        coefa(ifac,iclvar) = -rcodcl(ifac,ivar,3)/hint
        coefb(ifac,iclvar) = 1.d0
      endif

    enddo

  endif

  ! ---> SPALART ALLMARAS

elseif(iturb.eq.70) then

  ivar   = inusa
  iclvar = iclnu

  do ifac = 1, nfabor

    iel = ifabor(ifac)

    ! --- Proprietes physiques
    visclc = propce(iel,ipcvis)
    flumbf = propfb(ifac,ipprob(ifluma(inusa)))

    ! --- Grandeurs geometriques
    distbf = distb(ifac)
    hint = visclc/distbf

    !      C.L DE TYPE DIRICHLET AVEC OU SANS COEFFICIENT D'ECHANGE
    if( icodcl(ifac,ivar).eq.1) then
      hext = rcodcl(ifac,ivar,2)
      pimp = rcodcl(ifac,ivar,1)
      coefa(ifac,iclvar) = hext*pimp/(hint +hext)
      coefb(ifac,iclvar) = hint     /(hint +hext)
      !      C.L DE TYPE FLUX
    elseif(                                                     &
         icodcl(ifac,ivar).eq.3)then
      coefa(ifac,iclvar) = -rcodcl(ifac,ivar,3)/hint
      coefb(ifac,iclvar) = 1.d0
    endif

  enddo

endif

!===============================================================================
! 11. SCALAIRES (AUTRES QUE PRESSION, K, EPSILON, RIJ, OMEGA, VARIANCES)
!                     : DIRICHLET, NEUMANN
!===============================================================================

if(nscal.ge.1) then

  do ii = 1, nscal

    ivar   = isca(ii)
    iclvar = iclrtp(ivar,icoef)
    iclvaf = iclrtp(ivar,icoeff)

    isvhbl = 0
    if(ii.eq.isvhb) then
      isvhbl = isvhb
    endif

    if(ivisls(ii).gt.0) then
      ipcvsl = ipproc(ivisls(ii))
    else
      ipcvsl = 0
    endif

    ! --- Indicateur de prise en compte de Cp ou non
    !       (selon si le scalaire (scalaire associe pour une fluctuation)
    !        doit etre ou non traite comme une temperature)
    !      Si le scalaire est une variance et que le
    !        scalaire associe n'est pas resolu, on suppose alors qu'il
    !        doit etre traite comme un scalaire passif (defaut IHCP = 0)
    ihcp = 0
    if(iscavr(ii).le.nscal) then
      if(iscavr(ii).gt.0) then
        iscal = iscavr(ii)
      else
        iscal = ii
      endif
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
    endif

    ! --- Boucle sur les faces
    do ifac = 1, nfabor

      iel = ifabor(ifac)

      ! --- Proprietes physiques
      visctc = propce(iel,ipcvst)
      flumbf = propfb(ifac,ipprob(ifluma(ivar)))

      ! --- Grandeurs geometriques
      distbf = distb(ifac)

      ! --- Prise en compte de Cp ou CV
      !      (dans le Cas compressible IHCP=0)

      cpp = 1.d0
      if(ihcp.eq.0) then
        cpp = 1.d0
      elseif(ihcp.eq.2) then
        cpp = propce(iel,ipccp )
      elseif(ihcp.eq.1) then
        cpp = cp0
      endif
      hint = cpp

      ! --- Viscosite variable ou non
      if (ipcvsl.le.0) then
        rkl = visls0(ii)
      else
        rkl = propce(iel,ipcvsl)
      endif

      ! --- Cas turbulent
      if (iturb.ne.0) then
        hint = hint*(rkl+visctc/sigmas(ii))/distbf
        !     Cas laminaire
      else
        hint  = hint*rkl/distbf
      endif

      ! --->  C.L DE TYPE DIRICHLET AVEC OU SANS COEFFICIENT D'ECHANGE

      if( icodcl(ifac,ivar).eq.1) then
        hext = rcodcl(ifac,ivar,2)
        if(abs(hext).ge.rinfin*0.5d0) then
          pimp = rcodcl(ifac,ivar,1)
          coefa(ifac,iclvar) = pimp
          coefb(ifac,iclvar) = 0.d0
        else
          pimp = rcodcl(ifac,ivar,1)
          coefa(ifac,iclvar) = hext*pimp/(hint+hext)
          coefb(ifac,iclvar) = hint     /(hint+hext)
        endif

        !     On utilise le Dirichlet pour les calculs de gradients
        !       et pour les flux de bord.

        if(iclvaf.ne.iclvar) then
          coefa(ifac,iclvaf) = coefa(ifac,iclvar)
          coefb(ifac,iclvaf) = coefb(ifac,iclvar)
        endif

        ! ---> COUPLAGE : on stocke le hint (lambda/d      en temperature,
        !                                    lambda/(cp d) en enthalpie,
        !                                    lambda/(cv d) en energie)

        if (isvhbl .gt. 0) then
          hbord(ifac) = hint
        endif


        !--> Rayonnement :

        !      On stocke le coefficient d'echange lambda/distance
        !      (ou son equivalent en turbulent) quelle que soit la
        !      variable thermique transportee (temperature ou enthalpie)
        !      car on l'utilise pour realiser des bilans aux parois qui
        !      sont faits en temperature (on cherche la temperature de
        !      paroi quelle que soit la variable thermique transportee pour
        !      ecrire des eps sigma T4).

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
        !              CPR = CP0
        !            ENDIF
        !         puisque HINT = VISLS / DISTBR
        !                      = lambda/(CP * distance)

        !       lorsque la variable transportee est l'energie
        !         ISCSTH(II).EQ.3 :
        !         on procede comme pour l'enthalpie avec CV au lieu de CP
        !         (rq : il n'y a pas d'hypothèse, sf en non orthogonal :
        !               le flux est le bon et le coef d'echange aussi)

        !      De meme plus bas et de meme dans clptur.

        !               Si on rayonne et que
        !                  le scalaire est la variable energetique

        if (iirayo.ge.1 .and. ii.eq.iscalt) then

          !                On calcule le coefficient d'echange en W/(m2 K)

          !                  Si on resout en enthalpie
          if(iscsth(ii).eq.2) then
            !                    Si Cp variable
            if(ipccp.gt.0) then
              propfb(ifac,ipprob(ihconv)) = hint*propce(iel,ipccp )
            else
              propfb(ifac,ipprob(ihconv)) = hint*cp0
            endif
            !                  Si on resout en energie (compressible)
          elseif(iscsth(ii).eq.3) then
            !                    Si Cv variable
            if(ipccv.gt.0) then
              propfb(ifac,ipprob(ihconv)) = hint*propce(iel,ipccv )
            else
              propfb(ifac,ipprob(ihconv)) = hint*cv0
            endif
            !                  Si on resout en temperature
          elseif(abs(iscsth(ii)).eq.1) then
            propfb(ifac,ipprob(ihconv)) = hint
          endif

          !                On recupere le flux h(Ti'-Tp) (sortant ou
          !                             negatif si gain pour le fluide) en W/m2

          propfb(ifac,ipprob(ifconv)) =                       &
               hint*( (1.d0-coefb(ifac,iclvaf))*thbord(ifac)  &
               - coefa(ifac,iclvaf))

        endif

        ! --->  C.L DE TYPE FLUX

      elseif(icodcl(ifac,ivar).eq.3)then
        coefa(ifac,iclvaf)  = -rcodcl(ifac,ivar,3)/hint
        coefb(ifac,iclvaf)  = 1.d0
        if(iclvar.ne.iclvaf) then
          coefa(ifac,iclvar)  = 0.d0
          coefb(ifac,iclvar)  = 1.d0
        endif
        if (isvhbl .gt. 0) hbord(ifac) = hint

        !--> Rayonnement :

        if (iirayo.ge.1 .and. ii.eq.iscalt) then

          !                On calcule le coefficient d'echange en W/(m2 K)

          !                Si on resout en enthalpie
          if(iscsth(ii).eq.2) then
            !                  Si Cp variable
            if(ipccp.gt.0) then
              propfb(ifac,ipprob(ihconv)) = hint*propce(iel,ipccp )
            else
              propfb(ifac,ipprob(ihconv)) = hint*cp0
            endif
          elseif(iscsth(ii).eq.3) then
            !                    Si Cv variable
            if(ipccv.gt.0) then
              propfb(ifac,ipprob(ihconv)) = hint*propce(iel,ipccv )
            else
              propfb(ifac,ipprob(ihconv)) = hint*cv0
            endif
            !                Si on resout en temperature
          elseif(abs(iscsth(ii)).eq.1) then
            propfb(ifac,ipprob(ihconv)) = hint
          endif

          !              On recupere le flux h(Ti'-Tp) (sortant ou
          !                             negatif si gain pour le fluide)

          propfb(ifac,ipprob(ifconv)) = rcodcl(ifac,ivar,3)
        endif

      endif

    enddo

  enddo

endif

!===============================================================================
! 13.  VITESSE DE MAILLAGE EN ALE : DIRICHLET, NEUMANN
!===============================================================================
! (les conditions de glissement on ete traitees dans ALTYCL

if (iale.eq.1) then

  icluma = iclrtp(iuma ,icoef)
  iclvma = iclrtp(ivma ,icoef)
  iclwma = iclrtp(iwma ,icoef)

  do ifac = 1, nfabor

    iel = ifabor(ifac)
    distbf = distb(ifac)
    srfbn2 = surfbn(ifac)**2
    if (iortvm.eq.0) then
      hint = propce(iel,ipproc(ivisma(1)))/distbf
    else
      hint = ( propce(iel,ipproc(ivisma(1)))*surfbo(1,ifac)**2    &
             + propce(iel,ipproc(ivisma(2)))*surfbo(2,ifac)**2    &
             + propce(iel,ipproc(ivisma(3)))*surfbo(3,ifac)**2 )  &
           /distbf/srfbn2
    endif

    do isou = 1, 3
      if (isou.eq.1) ivar = iuma
      if (isou.eq.2) ivar = ivma
      if (isou.eq.3) ivar = iwma
      iclvar = iclrtp(ivar,icoef)

!      C.L DE TYPE DIRICHLET
      if( icodcl(ifac,ivar).eq.1 ) then
        hext = rcodcl(ifac,ivar,2)
        if(abs(hext).gt.rinfin*0.5d0) then
          pimp = rcodcl(ifac,ivar,1)
          coefa(ifac,iclvar) = pimp
          coefb(ifac,iclvar) = 0.d0
        else
          pimp = rcodcl(ifac,ivar,1)
          coefa(ifac,iclvar) = hext*pimp/(hint +hext)
          coefb(ifac,iclvar) = hint     /(hint +hext)
        endif

!      C.L DE TYPE FLUX
      elseif( icodcl(ifac,ivar).eq.3 ) then
        coefa(ifac,iclvar) = -rcodcl(ifac,ivar,3)/hint
        coefb(ifac,iclvar) = 1.d0
      endif

      ! Coupled solving of the velocity components
      if (ivelco.eq.1) then
        cfaale(isou,ifac) = coefa(ifac,iclvar)
        do jsou = 1, 3
          cfbale(isou,jsou,ifac) = 0.d0
          if(isou.eq.jsou) cfbale(isou,jsou,ifac) = coefb(ifac,iclvar)
        enddo
      endif

    enddo

    if (icodcl(ifac,iuma).eq.4) then
!     Face de glissement (si IUMA est de type 4, les autres aussi)
!     On force la vitesse de maillage normale a etre nulle, CL de Neumann
!       sur les autres composantes (calque sur CLSYVT). On prend directement
!       la valeur au centre de la cellule et pas reconstruite a la face (on
!       ne cherche pas une precision particuliere sur w)
      srfbnf = surfbn(ifac)
      rnx = surfbo(1,ifac)/srfbnf
      rny = surfbo(2,ifac)/srfbnf
      rnz = surfbo(3,ifac)/srfbnf
      upx = rtpa(iel,iuma)
      upy = rtpa(iel,ivma)
      upz = rtpa(iel,iwma)
      coefa(ifac,iuma) = - rnx*(rny*upy+rnz*upz)
      coefb(ifac,iuma) = 1.d0-rnx**2
      coefa(ifac,ivma) = - rny*(rnz*upz+rnx*upx)
      coefb(ifac,ivma) = 1.d0-rny**2
      coefa(ifac,iwma) = - rnz*(rnx*upx+rny*upy)
      coefb(ifac,iwma) = 1.d0-rnz**2
      ! Coupled solving of the velocity components
      if (ivelco.eq.1) then
        cfaale(1,ifac) = 0.d0
        cfaale(2,ifac) = 0.d0
        cfaale(3,ifac) = 0.d0

        cfbale(1,1,ifac) = 1.d0-rnx**2
        cfbale(2,2,ifac) = 1.d0-rny**2
        cfbale(3,3,ifac) = 1.d0-rnz**2

        cfbale(1,2,ifac) = -rnx*rny
        cfbale(2,1,ifac) = -rny*rnx
        cfbale(1,3,ifac) = -rnx*rnz
        cfbale(3,1,ifac) = -rnz*rnx
        cfbale(2,3,ifac) = -rny*rnz
        cfbale(3,2,ifac) = -rnz*rny
      endif
    endif

  enddo

endif

!===============================================================================
! 14.  CALCUL DES EFFORTS AUX BORDS (partie 1/5)
!===============================================================================

if (ineedf.eq.1) then

  if (ivelco.eq.0) then
    do ifac = 1, nfabor
      iel = ifabor(ifac)
      visclc = propce(iel,ipproc(iviscl))
      visctc = propce(iel,ipproc(ivisct))
      if (itytur.eq.3) then
        vistot = visclc
      else
        vistot = visclc + visctc
      endif
      distbf = distb(ifac)
      srfbnf = surfbn(ifac)
      forbr(1,ifac) = -vistot * ( coefa(ifac,icluf)  &
           + (coefb(ifac,icluf)-1.d0)*coefu(ifac,1) )/distbf*srfbnf
      forbr(2,ifac) = -vistot * ( coefa(ifac,iclvf)  &
           + (coefb(ifac,iclvf)-1.d0)*coefu(ifac,2) )/distbf*srfbnf
      forbr(3,ifac) = -vistot * ( coefa(ifac,iclwf)  &
           + (coefb(ifac,iclwf)-1.d0)*coefu(ifac,3) )/distbf*srfbnf
    enddo

  ! Coupled solving of the velocity components
  else
    do ifac = 1, nfabor
      iel = ifabor(ifac)
      visclc = propce(iel,ipproc(iviscl))
      visctc = propce(iel,ipproc(ivisct))
      if (itytur.eq.3) then
        vistot = visclc
      else
        vistot = visclc + visctc
      endif
      distbf = distb(ifac)
      srfbnf = surfbn(ifac)

      ! The implicit term is added after having updated the velocity
      forbr(1,ifac) = -vistot * ( cofafu(1,ifac)   &
           + (cofbfu(1,1,ifac)-1.d0)*coefu(ifac,1) &
           +  cofbfu(1,2,ifac)      *coefu(ifac,2) &
           +  cofbfu(1,3,ifac)      *coefu(ifac,3) )/distbf*srfbnf
      forbr(2,ifac) = -vistot * ( cofafu(2,ifac)   &
           +  cofbfu(2,1,ifac)      *coefu(ifac,1) &
           + (cofbfu(2,2,ifac)-1.d0)*coefu(ifac,2) &
           +  cofbfu(2,3,ifac)      *coefu(ifac,2) )/distbf*srfbnf
      forbr(3,ifac) = -vistot * ( coefau(3,ifac)   &
           +  cofbfu(3,1,ifac)      *coefu(ifac,1) &
           +  cofbfu(3,2,ifac)      *coefu(ifac,2) &
           + (cofbfu(3,3,ifac)-1.d0)*coefu(ifac,3) )/distbf*srfbnf
    enddo
  endif
endif

! Free memory
deallocate(coefu)
if (allocated(rijipb)) deallocate(rijipb)

!===============================================================================
! 15.  FORMATS
!===============================================================================

#if defined(_CS_LANG_FR)

 1000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET LORS DE L''ENTREE DES COND. LIM.      ',/,&
'@    =========                                               ',/,&
'@      INCOHERENCE ENTRE OPTIONS DE CALCUL ET COND. LIM.     ',/,&
'@                                                            ',/,&
'@      la prise en compte des termes d''echo de paroi        ',/,&
'@      du modele de turbulence Rij-epsilon est activee       ',/,&
'@      IRIJEC = ',I10,'                                      ',/,&
'@      Ou bien l amortissement de la viscosite turbulente    ',/,&
'@      est active IDRIES = ',I10,'en LES                     ',/,&
'@    mais aucune face de bord de type paroi n''est detectee. ',/,&
'@    L''incoherence indiquee ci-dessus n''est pas bloquante  ',/,&
'@      mais peut resulter d''une erreur lors de la           ',/,&
'@      specification des conditions aux limites.             ',/,&
'@                                                            ',/,&
'@    Par securite, le calcul ne sera pas execute.            ',/,&
'@                                                            ',/,&
'@    Verifier les conditions aux limites dans                ',/,&
'@      cs_user_boundary si le domaine comporte des parois.   ',/,&
'@    Eliminer l''option IRIJEC de usini1 si le domaine ne    ',/,&
'@      comporte pas de paroi (ou conditions ICODCL = 5 en    ',/,&
'@      vitesse).                                             ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 3010 format(                                                           &
 'Debit entrant retenu en ',I10   ,                  &
                                      ' faces de sortie sur ',I10)
 8000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : COND. LIM.                                  ',/,&
'@    =========                                               ',/,&
'@     Le scalaire ',I10   ,' est couple a SYRTHES            ',/,&
'@      mais n''est pas la variable energetique               ',/,&
'@         ISCALT = ',I10                               ,/,&
'@                                                            ',/,&
'@     Le calcul ne sera pas execute.                         ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

#else

 1000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: ABORT IN THE BOUNDARY CONDITIONS SPECIFICATION ',/,&
'@    =========                                               ',/,&
'@      INCOHERENCY BETWEEN CALCULATION OPTIONS AND BOUND COND',/,&
'@                                                            ',/,&
'@      The wall-echo terms of the Rij-epsilon turbulence     ',/,&
'@      model are taken into account                          ',/,&
'@      IRIJEC = ',I10,'                                      ',/,&
'@      Or the Van Driest damping of the turbulent viscosity  ',/,&
'@      is active IDRIES = ',I10,'in LES                      ',/,&
'@    but no wall boundary face is detected.                  ',/,&
'@    This incoherency is not blocking but may result from    ',/,&
'@      an error during the boundary conditions               ',/,&
'@      specification.                                        ',/,&
'@                                                            ',/,&
'@    By safety, the calculation will not be run.             ',/,&
'@                                                            ',/,&
'@    Verify the boundary conditions in cs_user_boundary in   ',/,&
'@      if the domain has any walls.                          ',/,&
'@    Remove the option IRIJEC from usini1 if the domain does ',/,&
'@      not have any wall (or conditions ICODCL = 5 for the   ',/,&
'@      velocity).                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 3010 format(                                                           &
 'Incoming flow detained for ', I10   ,              &
                                          ' outlet faces on ',I10)
 8000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: BOUNDARY CONDITIONS                            ',/,&
'@    ========                                                ',/,&
'@     The scalar ',I10   ,' is coupled with SYRTHES          ',/,&
'@      but is not the energy variable                        ',/,&
'@         ISCALT = ',I10                               ,/,&
'@                                                            ',/,&
'@     The calculation will not be run.                       ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

#endif

!----
! FIN
!----

return
end subroutine
