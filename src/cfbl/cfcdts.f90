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

subroutine cfcdts &
!================

 ( idbia0 , idbra0 ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ivar   , iconvp , idiffp , ireslp , ndircp , nitmap ,          &
   imrgra , nswrsp , nswrgp , imligp , ircflp ,                   &
   ischcp , isstpp , iescap , iifbru ,                            &
   imgrp  , ncymxp , nitmfp , ipp    , iwarnp ,                   &
   blencp , epsilp , epsrsp , epsrgp , climgp , extrap , thetap , &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml , ifrusb ,          &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   pvara  , coefap , coefbp , cofafp , cofbfp , flumas , flumab , &
   viscfm , viscbm , viscfs , viscbs ,                            &
   rovsdt , smbrp  , pvar   ,                                     &
   dam    , xam    , dpvar  ,                                     &
   w1     , w2     , w3     , w4     , w5     ,                   &
   w6     , w7     , w8     , smbini ,                            &
   rdevel , rtuser , ra     )

!===============================================================================
! FONCTION :
! ----------

! RESOLUTION SUR UN PAS DE TEMPS D'UNE EQUATION DE CONVECTION
! /DIFFUSION/TERME SOURCE POUR LA VARIABLE PVAR EN COMPRESSIBLE

! ROVSDT.( PVAR -PVARA )
!        (    ->           --->        )
!   + DIV( RO.U  PVAR -VISC GRAD( PVAR ) ).VOLUME  = SMBR
!        (                             )

! ON RESOUT EN FAIT :

! ROVSDT.DPVAR
!        (    ->            --->         )
!   + DIV( RO.U  DPVAR -VISC GRAD( DPVAR ) ).VOLUME  = SMBR
!        (                               )
! AVEC
!  SMBR = SMBR
!        (    ->n           --->      n  )
!   - DIV( RO.U  DPVAR -VISC GRAD( DPVAR ) ).VOLUME
!        (                               )
! ET DPVAR = INCREMENT VARIABLE PVAR

! Attention, on suppose qu'on arrive ici avec PVAR initialise
!  Y COMPRIS DANS LE HALO EN PERIODICITE (appel a promav pour
!  le calcul de la norme)






!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
! ndim             ! i  ! <-- ! spatial dimension                              !
! ncelet           ! i  ! <-- ! number of extended (real + ghost) cells        !
! ncel             ! i  ! <-- ! number of cells                                !
! nfac             ! i  ! <-- ! number of interior faces                       !
! nfabor           ! i  ! <-- ! number of boundary faces                       !
! nfml             ! i  ! <-- ! number of families (group classes)             !
! nprfml           ! i  ! <-- ! number of properties per family (group class)  !
! nnod             ! i  ! <-- ! number of vertices                             !
! lndfac           ! i  ! <-- ! size of nodfac indexed array                   !
! lndfbr           ! i  ! <-- ! size of nodfbr indexed array                   !
! ncelbr           ! i  ! <-- ! number of cells with faces on boundary         !
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! nphas            ! i  ! <-- ! number of phases                               !
! nideve, nrdeve   ! i  ! <-- ! sizes of idevel and rdevel arrays              !
! nituse, nrtuse   ! i  ! <-- ! sizes of ituser and rtuser arrays              !
! iconvp           ! e  ! <-- ! indicateur = 1 convection, 0 sinon             !
! idiffp           ! e  ! <-- ! indicateur = 1 diffusion , 0 sinon             !
! ndircp           ! e  ! <-- ! indicateur = 0 si decalage diagonale           !
! ireslp           ! e  ! <-- ! indicateur = 0 gradco                          !
!                  !    !     !            = 1 jacobi                          !
!                  !    !     !            = 2 bi-cgstab                       !
! imrgra           ! e  ! <-- ! indicateur = 0 gradrc 97                       !
!                  ! e  ! <-- !            = 1 gradmc 99                       !
! nswrsp           ! e  ! <-- ! nombre de sweep pour reconstruction            !
!                  !    !     !             du second membre                   !
! nswrgp           ! e  ! <-- ! nombre de sweep pour reconstruction            !
!                  !    !     !             des gradients                      !
! imligp           ! e  ! <-- ! methode de limitation du gradient              !
!                  !    !     !  < 0 pas de limitation                         !
!                  !    !     !  = 0 a partir des gradients voisins            !
!                  !    !     !  = 1 a partir du gradient moyen                !
! ircflp           ! e  ! <-- ! indicateur = 1 rec flux ; 0 sinon              !
! ischcp           ! e  ! <-- ! indicateur = 1 centre , 0 2nd order            !
! isstpp           ! e  ! <-- ! indicateur = 1 sans test de pente              !
!                  !    !     !            = 0 avec test de pente              !
! iescap           ! e  ! <-- ! =1 calcul de l'indicateur prediction           !
! iifbru           ! e  ! <-- ! pointeur flux de bord rusanov                  !
! imgrp            ! e  ! <-- ! indicateur = 0 pas de mgm                      !
!                  !    !     !            = 1 sinon                           !
! nitmap           ! e  ! <-- ! nombre max d'iter pour resol iterativ          !
! ipp              ! e  ! <-- ! numero de variable pour post                   !
! iwarnp           ! i  ! <-- ! verbosity                                      !
! nfecrl           ! e  ! <-- ! unite du fichier sortie std                    !
! blencp           ! r  ! <-- ! 1 - proportion d'upwind                        !
! epsilp           ! r  ! <-- ! precision pour resol iter                      !
! epsrgp           ! r  ! <-- ! precision relative pour la                     !
!                  !    !     !  reconstruction des gradients 97               !
! climgp           ! r  ! <-- ! coef gradient*distance/ecart                   !
! extrap           ! r  ! <-- ! coef extrap gradient                           !
! ifacel(2, nfac)  ! ia ! <-- ! interior faces -> cells connectivity           !
! ifabor(nfabor)   ! ia ! <-- ! boundary faces -> cells connectivity           !
! ifmfbr(nfabor)   ! ia ! <-- ! boundary face family numbers                   !
! ifmcel(ncelet)   ! ia ! <-- ! cell family numbers                            !
! iprfml           ! ia ! <-- ! property numbers per family                    !
!  (nfml, nprfml)  !    !     !                                                !
! ifrusb(nfabor    ! te ! <-- ! indicateur flux de rusanov                     !
! ipnfac(nfac+1)   ! ia ! <-- ! interior faces -> vertices index (optional)    !
! nodfac(lndfac)   ! ia ! <-- ! interior faces -> vertices list (optional)     !
! ipnfbr(nfabor+1) ! ia ! <-- ! boundary faces -> vertices index (optional)    !
! nodfbr(lndfbr)   ! ia ! <-- ! boundary faces -> vertices list (optional)     !
! idevel(nideve)   ! ia ! <-> ! integer work array for temporary development   !
! ituser(nituse)   ! ia ! <-> ! user-reserved integer work array               !
! ia(*)            ! ia ! --- ! main integer work array                        !
! xyzcen           ! ra ! <-- ! cell centers                                   !
!  (ndim, ncelet)  !    !     !                                                !
! surfac           ! ra ! <-- ! interior faces surface vectors                 !
!  (ndim, nfac)    !    !     !                                                !
! surfbo           ! ra ! <-- ! boundary faces surface vectors                 !
!  (ndim, nfabor)  !    !     !                                                !
! cdgfac           ! ra ! <-- ! interior faces centers of gravity              !
!  (ndim, nfac)    !    !     !                                                !
! cdgfbo           ! ra ! <-- ! boundary faces centers of gravity              !
!  (ndim, nfabor)  !    !     !                                                !
! xyznod           ! ra ! <-- ! vertex coordinates (optional)                  !
!  (ndim, nnod)    !    !     !                                                !
! volume(ncelet)   ! ra ! <-- ! cell volumes                                   !
! pvara(ncelet     ! tr ! <-- ! variable resolue (instant precedent)           !
! coefap, b        ! tr ! <-- ! tableaux des cond lim pour p                   !
!   (nfabor)       !    !     !  sur la normale a la face de bord              !
! cofafp, b        ! tr ! <-- ! tableaux des cond lim pour le flux de          !
!   (nfabor)       !    !     !  diffusion de p                                !
! flumas(nfac)     ! tr ! <-- ! flux de masse aux faces internes               !
! flumab(nfabor    ! tr ! <-- ! flux de masse aux faces de bord                !
! viscfm(nfac)     ! tr ! <-- ! visc*surface/dist aux faces internes           !
!                  !    !     !  pour la matrice                               !
! viscbm(nfabor    ! tr ! <-- ! visc*surface/dist aux faces de bord            !
!                  !    !     !  pour la matrice                               !
! viscfs(nfac)     ! tr ! <-- ! idem viscfm pour second membre                 !
! viscbs(nfabor    ! tr ! <-- ! idem viscbm pour second membre                 !
! rovsdt(ncelet    ! tr ! <-- ! rho*volume/dt                                  !
! smbrp(ncelet     ! tr ! <-- ! bilan au second membre                         !
! pvar (ncelet     ! tr ! <-- ! variable resolue                               !
! dam(ncelet       ! tr ! --> ! tableau de travail pour matrice                !
!                  !    !     !  et resultat estimateur                        !
! xam(nfac,*)      ! tr ! --- ! tableau de travail pour matrice                !
! w1...8(ncelet    ! tr ! --- ! tableau de travail                             !
! smbini(ncelet    ! tr ! --- ! tableau de travail                             !
! rdevel(nrdeve)   ! ra ! <-> ! real work array for temporary development      !
! rtuser(nrtuse)   ! ra ! <-> ! user-reserved real work array                  !
! ra(*)            ! ra ! --- ! main real work array                           !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

implicit none

!===============================================================================
! Common blocks
!===============================================================================

include "paramx.f90"
include "numvar.f90"
include "cstnum.f90"
include "entsor.f90"
include "period.f90"
include "parall.f90"

!===============================================================================

! Arguments

integer          idbia0 , idbra0
integer          ndim   , ncelet , ncel   , nfac   , nfabor
integer          nfml   , nprfml
integer          nnod   , lndfac , lndfbr , ncelbr
integer          nvar   , nscal  , nphas
integer          nideve , nrdeve , nituse , nrtuse
integer          ivar   , iconvp , idiffp , ndircp
integer          nitmap
integer          imrgra , nswrsp , nswrgp , imligp , ircflp
integer          ischcp , isstpp , iescap , iifbru , imgrp
integer          ncymxp , nitmfp
integer          ipp    , iwarnp
double precision blencp , epsilp , epsrgp , climgp , extrap
double precision thetap , epsrsp

integer          ifacel(2,nfac) , ifabor(nfabor)
integer          ifmfbr(nfabor) , ifmcel(ncelet)
integer          iprfml(nfml,nprfml)
integer          ifrusb(nfabor)
integer          ipnfac(nfac+1), nodfac(lndfac)
integer          ipnfbr(nfabor+1), nodfbr(lndfbr)
integer          idevel(nideve), ituser(nituse)
integer          ia(*)

double precision xyzcen(ndim,ncelet)
double precision surfac(ndim,nfac), surfbo(ndim,nfabor)
double precision cdgfac(ndim,nfac), cdgfbo(ndim,nfabor)
double precision xyznod(ndim,nnod), volume(ncelet)
double precision pvara(ncelet), coefap(nfabor), coefbp(nfabor)
double precision                cofafp(nfabor), cofbfp(nfabor)
double precision flumas(nfac), flumab(nfabor)
double precision viscfm(nfac), viscbm(nfabor)
double precision viscfs(nfac), viscbs(nfabor)
double precision rovsdt(ncelet), smbrp(ncelet)
double precision pvar(ncelet)
double precision dam(ncelet), xam(nfac ,2)
double precision dpvar(ncelet)
double precision w1(ncelet), w2(ncelet), w3(ncelet), w4(ncelet)
double precision w5(ncelet), w6(ncelet), w7(ncelet), w8(ncelet)
double precision smbini(ncelet)
double precision rdevel(nrdeve), rtuser(nrtuse), ra(*)

! Local variables

character*80     chaine
character*8      cnom
integer          idebia, idebra
integer          isym,ireslp,ireslq,ipol,isqrt
integer          inc,isweep,niterf,iccocg,iel,icycle,nswmod
integer          iphas,idimte,itenso,iinvpe, iinvpp
double precision residu,rnorm

!===============================================================================

!===============================================================================
! 1.  INITIALISATIONS
!===============================================================================

idebia = idbia0
idebra = idbra0

! NOMS
chaine = nomvar(ipp)
cnom   = chaine(1:8)

! MATRICE A PRIORI SYMETRIQUE ( = 1)
isym  = 1
if( iconvp.gt.0 ) isym  = 2


! METHODE DE RESOLUTION ET DEGRE DU PRECOND DE NEUMANN
!     0 SI CHOIX AUTOMATIQUE GRADCO OU BICGSTAB
!     0 SI CHOIX AUTOMATIQUE JACOBI
!     DONNE PAR IRESLP/1000 SI NON AUTOMATIQUE
if (ireslp.eq.-1) then
  ireslq = 0
  ipol   = 0
  if( iconvp.gt.0 ) then
    ireslq = 1
    ipol   = 0
  endif
else
  ireslq = mod(ireslp,1000)
  ipol   = (ireslp-ireslq)/1000
endif


! PRISE DE SQRT DANS PS
isqrt = 1

! PRISE EN COMPTE DE LA PERIODICITE

!    Initialisation pour test avant promav
iinvpe = 0

if(iperio.eq.1) then


!    Par defaut, toutes les periodicites seront traitees dans percom,
!      les variables etant assimilees a des scalaires (meme si ce sont
!      des composantes de vecteurs ou de tenseur)
  idimte = 0
  itenso = 0

  iinvpe = 1

  do iphas = 1, nphas
    if(ivar.eq.iu(iphas).or.ivar.eq.iv(iphas).or.                 &
                            ivar.eq.iw(iphas).or.                 &
       ivar.eq.ir11(iphas).or.ivar.eq.ir12(iphas).or.             &
       ivar.eq.ir13(iphas).or.ivar.eq.ir22(iphas).or.             &
       ivar.eq.ir23(iphas).or.ivar.eq.ir33(iphas)) then

!    Pour la vitesse et les tensions de Reynolds
!      seules seront echangees les informations sur les faces periodiques
!      de translation dans percom ; on ne touche pas aux informations
!      relatives aux faces de periodicite de rotation.
      idimte = 0
      itenso = 1

!      Lors de la resolution par increments, on echangera egalement les
!      informations relatives aux faces de periodicite de translation.
!      Pour les faces de periodicite de rotation, l'increment sera
!      annule dans percom (iinvpe=2).
      iinvpe = 2

    endif
  enddo

endif

!===============================================================================
! 1.  CONSTRUCTION MATRICE "SIMPLIFIEE" DE RESOLUTION
!===============================================================================

call matrix                                                       &
!==========
 ( ncelet , ncel   , nfac   , nfabor ,                            &
   iconvp , idiffp , ndircp , isym   , nfecra ,                   &
   thetap ,                                                       &
   ifacel , ifabor ,                                              &
   coefbp , rovsdt , flumas , flumab , viscfm , viscbm ,          &
   dam    , xam    )


!===============================================================================
! 2.  BOUCLES SUR LES NON ORTHOGONALITES
!       (A PARTIR DE LA SECONDE ITERATION)
!===============================================================================


!     AVANT DE BOUCLER SUR LES SWEEP, ON STOCKE LE SECOND MEMBRE SANS
!     RECONSTRUCTION DANS LE TABLEAU AUXILIAIRE SMBINI

do iel = 1, ncel
   smbini(iel) = smbrp(iel)
enddo

!     On passe toujours dans cfbsc2 avec INC=1
inc = 1
!     Sauf pour les matrices poids (NSWRSP=-1 ... a priori bloque dans
!                                   verini en compressible)
if (nswrsp.eq.-1) then
  nswrsp = 1
  inc = 0
endif


!  Ca serait bien de le simplifier plus tard aussi ce nombre de sweeps
!  Attention, pour les matrices poids il faut pouvoir ne faire
!     qu'un seul sweep
nswmod = max( nswrsp, 1 )
do 100 isweep = 1, nswmod

! ---> INCREMENTATION ET RECONSTRUCTION DU SECOND MEMBRE
!      ON NE RECALCULE COCG QU'AU PREMIER PASSAGE (PRESQUE)

  if( isweep.eq.1) then
    iccocg = 1

  else
    iccocg = 0
    do iel = 1, ncel
!     SMBINI CONTIENT LES TERMES INSTAT, EN DIV(RHO U) ET SOURCE DE MASSE
!     DU SECOND MEMBRE  MIS A JOUR A CHAQUE SWEEP
      smbini(iel) = smbini(iel) - rovsdt(iel)*dpvar(iel)
      smbrp(iel)  = smbini(iel)
    enddo
  endif

  call cfbsc2                                                     &
  !==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ivar   , iconvp , idiffp , nswrgp , imligp , ircflp ,          &
   ischcp , isstpp , inc    , imrgra , iccocg , iifbru ,          &
   ipp    , iwarnp ,                                              &
   blencp , epsrgp , climgp , extrap ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml , ifrusb ,          &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   pvar   , coefap , coefbp , cofafp , cofbfp ,                   &
   flumas , flumab , viscfs , viscbs ,                            &
   smbrp  ,                                                       &
!        ------
   w1     , w2     , w3     , w4     , w5     , w6     ,          &
   rdevel , rtuser , ra     )


  call prodsc(ncelet,ncel,isqrt,smbrp,smbrp,residu)

! ---> RESIDU DE NORMALISATION CALCULE AU PREMIER SWEEP
!    (NORME C.L +TERMES SOURCES+ TERMES DE NON ORTHOGONALITE)

!       Attention, lors de l'appel a promav, ici pour une variable qui
!         n'est pas en increments et qui est supposee initialisee
!         y compris dans le halo.
!         Pour les variables vitesse et les tensions de Reynolds
!         (IINVPE=2), il ne faudra donc pas annuler le halo
!         des periodicites de rotation, mais au contraire le laisser
!         inchange.
!         Pour les autres variables (scalaires) IINVPE=1 permettra de
!         tout echanger, meme si c'est superflu.
  if( isweep.eq.1 ) then
    if(iinvpe.eq.2) then
      iinvpp = 3
    else
      iinvpp = iinvpe
    endif
    call promav(ncelet,ncel,nfac,isym,iinvpp,ifacel,              &
                dam,xam,pvar,w1)
    do iel = 1, ncel
      w1(iel) = w1(iel) + smbrp(iel)
    enddo
    call prodsc(ncelet,ncel,isqrt,w1,w1,rnorm)
    rnsmbr(ipp) = rnorm
  endif

! ---> RESOLUTION IMPLICITE SUR L'INCREMENT DPVAR

  do iel = 1, ncel
    dpvar(iel) = 0.d0
  enddo

  call invers                                                     &
  !==========
 ( cnom   , idebia , idebra ,                                     &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nideve , nrdeve , nituse , nrtuse ,                            &
   isym   , ipol   , ireslq , nitmap , imgrp  ,                   &
   ncymxp , nitmfp ,                                              &
   iwarnp , nfecra , niterf , icycle , iinvpe ,                   &
   epsilp , rnorm  , residu ,                                     &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dam    , xam    , smbrp  , dpvar  ,                            &
   w3     , w4     , w5     , w6     , w7     , w8     ,          &
   rdevel , rtuser , ra     )


  nbivar(ipp) = niterf
  if(abs(rnorm).gt.epzero) then
    resvar(ipp) = residu/rnorm
  else
    resvar(ipp) = 0.d0
  endif

! ---> INCREMENTATION SOLUTION

  do iel = 1, ncel
    pvar(iel) = pvar(iel)+dpvar(iel)
  enddo

! ---> TRAITEMENT DU PARALLELISME

if(irangp.ge.0) call parcom (pvar)
                !==========

! ---> TRAITEMENT DE LA PERIODICITE : SEULE LA PERIODICITE IMPLICITE
!      EST ASSUREE (SCALAIRE ET TRANSLATION DE VECTEUR ET DE TENSEUR)

if(iperio.eq.1) then
  call percom                                                     &
  !==========
  ( idimte , itenso ,                                             &
    pvar   , pvar   , pvar  ,                                     &
    pvar   , pvar   , pvar  ,                                     &
    pvar   , pvar   , pvar  )
endif

! ---> TEST DE CONVERGENCE

call prodsc(ncelet,ncel,isqrt,smbrp,smbrp,residu)

if( residu.le.epsrsp*rnorm ) then
   if(iwarnp.ge.1) then
      write( nfecra,1000) cnom,isweep,residu,rnorm
   endif
   goto 200
endif

if(iwarnp.ge.3) then
   write( nfecra,1000) cnom,isweep,residu,rnorm
endif

 100  continue

if(iwarnp.ge.2) then
   write( nfecra,1100)cnom, nswmod
endif

!===============================================================================
! 3.  SORTIE OU CALCUL D'ESTIMATEURS POUR LES VITESSES
!       A L'ETAPE DE PREDICTION
!===============================================================================

 200  continue

! ---> TEST DE PASSAGE DANS LE CALCUL

if (iescap.gt.0) then

! ---> CALCUL DE LA CONTRIBUTION COMPOSANTE PAR COMPOSANTE. DE L ESTIMATEUR


!     SMBINI CONTIENT LES TERMES INSTAT ET EN DIV(U) DU SECOND MEMBRE
!     MIS A JOUR A CHAQUE SWEEP,DONC AU DERNIER, POUR KMAX +1, ON A:

  do iel = 1,ncel
    smbrp(iel) = smbini(iel) - rovsdt(iel)*dpvar(iel)
  enddo

  inc    = 1
  iccocg = 1

  call cfbsc2                                                     &
  !==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ivar   , iconvp , idiffp , nswrgp , imligp , ircflp ,          &
   ischcp , isstpp , inc    , imrgra , iccocg , iifbru ,          &
   ipp    , iwarnp ,                                              &
   blencp , epsrgp , climgp , extrap ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml , ifrusb ,          &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   pvar   , coefap , coefbp , cofafp , cofbfp ,                   &
   flumas , flumab , viscfs , viscbs ,                            &
   smbrp  ,                                                       &
!        ------
   w1     , w2     , w3     , w4     , w5     , w6     ,          &
   rdevel , rtuser , ra     )

!     CONTRIBUTION DES NORMES L2 DES DIFFERENTES COMPOSANTES
!       DANS LE TABLEAU DAM QUI EST ICI DISPONIBLE.

  do iel = 1,ncel
    dam(iel) = (smbrp(iel)/ volume(iel))**2
  enddo

endif

!--------
! FORMATS
!--------

 1000 format (                                                          &
 1X,A8,' : CV-DIF-TS',I5,' IT - RES= ',E12.5,' NORME= ', E12.5)
 1100 format (                                                          &
'@                                                            ',/,&
'@ @@ ATTENTION : ',A8 ,' CONVECTION-DIFFUSION-TERMES SOURCES ',/,&
'@    =========                                               ',/,&
'@  Nombre d''iterations maximal ',I10   ,' atteint           ',/,&
'@                                                            '  )

!12345678 : CV-DIF-TS 2000 IT - RES= 1234567890234 NORME= 12345678901234
!ATTENTION 12345678 : NON CONVERGENCE DU SYSTEME CONV-DIFF-TS
!----
! FIN
!----

end subroutine
