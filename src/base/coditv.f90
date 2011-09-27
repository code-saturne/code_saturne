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

subroutine coditv &
!================

 ( nvar   , nscal  ,                                              &
   idtvar , ivar   , iconvp , idiffp , ireslp , ndircp , nitmap , &
   imrgra , nswrsp , nswrgp , imligp , ircflp , ivisep ,          &
   ischcp , isstpp , iescap ,                                     &
   imgrp  , ncymxp , nitmfp , ippu   , ippv   , ippw   , iwarnp , &
   blencp , epsilp , epsrsp , epsrgp , climgp , extrap ,          &
   relaxp , thetap ,                                              &
   pvara  , pvark  ,                                              &
   coefav , coefbv , cofafv , cofbfv ,                            &
   flumas , flumab ,                                              &
   viscfm , viscbm , viscfs , viscbs , secvif , secvib ,          &
   fimp   ,                                                       &
   smbr   ,                                                       &
   pvar   ,                                                       &
   eswork )

!===============================================================================
! FONCTION :
! ----------

! RESOLUTION SUR UN PAS DE TEMPS DE L'EQUATION DE CONVECTION
! /DIFFUSION/TERME SOURCE POUR LA VARIABLE (UX, UY, UZ)

! ------>
! ------>  ->   -->
! FIMPuv.( U  - UA )
!                          --->
!        ( ->     ->       --->  ->  )           ->
!   + DIV( U x RO.U  -VISC GRAD( U ) ).VOLUME  = S
!        (                           )

! ON RESOUT EN FAIT :

! ------>
! ------>  -->
! FIMPuv.( DU )
!                           --->
!        ( -->     ->       --->  -->  )           ----->
!   + DIV( DU x RO.U  -VISC GRAD( DU ) ).VOLUME  = SMBR1
!        (                             )

! AVEC
!  ----->  ->
!  SMBR1 = S
!                           --->
!        ( -->     ->n      --->  -->  )
!   - DIV( UA x RO.U  -VISC GRAD( UA ) ).VOLUME
!        (                             )
!    -->                     ->  ->  -->
! ET DU = INCREMENT VARIABLE U = U - UA

! ET PLUS EXACTEMENT

! ------>
! ------>  -->
! FIMPuv.( DU )
!                           --->
!        ( -->     ->       --->  -->  )           ----->
!   + DIV( DU x RO.U  -VISC GRAD( DU ) ).VOLUME  = SMBR2
!        (                             )


! AVEC
!              ------>
!  ----->  ->  ------>  -->  -->
!  SMBR2 = S - FIMPuv.( Ui - UA )
!                           --->
!        ( -->     ->n      --->  -->  )
!   - DIV( Ui x RO.U  -VISC GRAD( Ui ) ).VOLUME
!        (                             )
!    -->                     ->  ->  ->
! ET DU = INCREMENT VARIABLE U = U - Ui


! ATTENTION : IL EST INTERDIT DE MODIFIER FIMPuv ICI
!                    ========             ======

!             il sert 32 fois dans le rayonnement (raysol).
!             et dans le calcul de y+ (distyp)

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! idtvar           ! e  ! <-- ! indicateur du schema temporel                  !
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
! ivisep           ! e  ! <-- ! indicateur = 1 pour la prise en compte         !
!                  !    !     !                div(T Grad(vel))                !
!                  !    !     !                -2/3 Grad(div(vel))             !
!                  !    !     !              0 sinon                           !
! ischcp           ! e  ! <-- ! indicateur = 1 centre , 0 2nd order            !
! isstpp           ! e  ! <-- ! indicateur = 1 sans test de pente              !
!                  !    !     !            = 0 avec test de pente              !
! iescap           ! e  ! <-- ! =1 calcul de l'indicateur prediction           !
! imgrp            ! e  ! <-- ! indicateur = 0 pas de mgm                      !
!                  !    !     !            = 1 sinon                           !
! nitmap           ! e  ! <-- ! nombre max d'iter pour resol iterativ          !
! ipp*             ! e  ! <-- ! numero de variable pour post                   !
! iwarnp           ! i  ! <-- ! verbosity                                      !
! nfecrl           ! e  ! <-- ! unite du fichier sortie std                    !
! blencp           ! r  ! <-- ! 1 - proportion d'upwind                        !
! epsilp           ! r  ! <-- ! precision pour resol iter                      !
! epsrgp           ! r  ! <-- ! precision relative pour la                     !
!                  !    !     !  reconstruction des gradients 97               !
! climgp           ! r  ! <-- ! coef gradient*distance/ecart                   !
! extrap           ! r  ! <-- ! coef extrap gradient                           !
! relaxp           ! r  ! <-- ! coefficient de relaxation                      !
! thetap           ! r  ! <-- ! coefficient de ponderation pour le             !
!                  !    !     ! theta-schema (on ne l'utilise pour le          !
!                  !    !     ! moment que pour u,v,w et les scalaire          !
!                  !    !     ! - thetap = 0.5 correspond a un schema          !
!                  !    !     !   totalement centre en temps (mixage           !
!                  !    !     !   entre crank-nicolson et adams-               !
!                  !    !     !   bashforth)                                   !
! pvara(3,ncelet)  ! tr ! <-- ! variable resolue (instant precedent)           !
! pvark(3,ncelet)  ! tr ! <-- ! variable de la sous-iteration                  !
!                  !    !     !  precedente. pour un point fixe sur            !
!                  !    !     !  navsto elle permet d'initialiser par          !
!                  !    !     !  autre chose que pvara (elle vaut              !
!                  !    !     !  pvar=pvara pour les scalaires)                !
! coefav           ! tr ! <-- ! tableaux des cond lim pour u, v, w             !
!   (3,nfabor)     !    !     !  sur la normale a la face de bord              !
! coefbv           ! tr ! <-- ! tableaux des cond lim pour u, v, w             !
!   (3,3,nfabor)   !    !     !  sur la normale a la face de bord              !
! cofafv           ! tr ! <-- ! tableaux des cond lim pour le flux de          !
!   (3,nfabor)     !    !     !  diffusion de u, v, w                          !
! cofbfv           ! tr ! <-- ! tableaux des cond lim pour le flux de          !
!   (3,3,nfabor)   !    !     !  diffusion de u, v, w                          !
! flumas(nfac)     ! tr ! <-- ! flux de masse aux faces internes               !
! flumab(nfabor    ! tr ! <-- ! flux de masse aux faces de bord                !
! viscfm(nfac)     ! tr ! <-- ! visc*surface/dist aux faces internes           !
!                  !    !     !  pour la matrice                               !
! viscbm(nfabor    ! tr ! <-- ! visc*surface/dist aux faces de bord            !
!                  !    !     !  pour la matrice                               !
! viscfs(nfac)     ! tr ! <-- ! idem viscfm pour second membre                 !
! viscbs(nfabor    ! tr ! <-- ! idem viscbm pour second membre                 !
! secvif(nfac)     ! tr ! --- ! secondary viscosity at interior faces          !
! secvib(nfabor)   ! tr ! --- ! secondary viscosity at boundary faces          !
! fimp(3,3,ncelet) ! tr ! <-- ! rho*volume/dt                                  !
! smbr(3,ncelet)   ! tr ! <-- ! bilan au second membre                         !
! pvar(3,ncelet)   ! tr ! <-- ! variable resolue                               !
! eswork(3,ncelet) ! ra ! <-- ! prediction-stage error estimator (iescap > 0)  !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use dimens, only: ndimfb
use paramx
use numvar
use cstnum
use entsor
use parall
use period
use mltgrd
use optcal, only: rlxp1
use mesh

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal
integer          idtvar , ivar   , iconvp , idiffp , ndircp
integer          nitmap
integer          imrgra , nswrsp , nswrgp , imligp , ircflp
integer          ischcp , isstpp , iescap , imgrp
integer          ncymxp , nitmfp
integer          iwarnp
integer          ippu   , ippv   , ippw   , ivisep
double precision blencp , epsilp , epsrgp , climgp , extrap
double precision relaxp , thetap , epsrsp

double precision pvara(3,ncelet)
double precision pvark(3,ncelet)
double precision pvar(3,ncelet)
double precision coefav(3,ndimfb)
double precision cofafv(3,ndimfb)
double precision coefbv(3,3,ndimfb)
double precision cofbfv(3,3,ndimfb)
double precision flumas(nfac), flumab(nfabor)
double precision viscfm(nfac), viscbm(nfabor)
double precision viscfs(nfac), viscbs(nfabor)
double precision secvif(nfac), secvib(nfabor)
double precision fimp(3,3,ncelet)
double precision smbr(3,ncelet)
double precision eswork(3,ncelet)

! Local variables


character*80     chaine
character*16     cnom(3)
integer          lchain
integer          isym,ireslp,ireslq,ipol,isqrt
integer          inc,isweep,niterf,iccocg,iel,icycle,nswmod
integer          iinvpe,iinvpp
integer          idtva0
integer          iagmax, nagmax, npstmg
double precision thetex

integer          isou , jsou, ifac
integer          ibsize
double precision residu, rnorm

double precision tps1, tps2, tempsjf

double precision, allocatable, dimension(:,:,:) :: dam
double precision, allocatable, dimension(:,:) :: xam
double precision, allocatable, dimension(:,:) :: dpvar, smbini, w1

save tempsjf
data tempsjf /0/

!===============================================================================

!===============================================================================
! 1.  INITIALISATIONS
!===============================================================================

! Allocate temporary arrays
! be carefull here, xam is interleaved
allocate(dam(3,3,ncelet), xam(2,nfac))
allocate(dpvar(3,ncelet), smbini(3,ncelet))

! Names
chaine = nomvar(ippu)
cnom(1)= chaine(1:16)
chaine = nomvar(ippv)
cnom(2)= chaine(1:16)
chaine = nomvar(ippw)
cnom(3)= chaine(1:16)

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

! LA PRISE EN COMPTE DE LA PERIODICITE EST AUTOMATIQUE

!    Initialisation pour test avant promav
iinvpe = 0

!===============================================================================
! 1.  CONSTRUCTION MATRICE "SIMPLIFIEE" DE RESOLUTION
!===============================================================================

! xam est le meme pour les 3 composantes

call matrxv                                                       &
!==========
 ( ncelet , ncel   , nfac   , nfabor ,                            &
   iconvp , idiffp , ndircp , isym   , nfecra ,                   &
   thetap ,                                                       &
   ifacel , ifabor ,                                              &
   coefbv , cofbfv , fimp   ,                                     &
   flumas , flumab , viscfm , viscbm ,                            &
   dam    , xam    )

!     En stationnaire, on relaxe la diagonale
if (idtvar.lt.0) then
  do iel = 1, ncel
    do isou=1,3
      do jsou=1,3
        dam(isou,jsou,iel) = dam(isou,jsou,iel)/relaxp
      enddo
    enddo
  enddo
endif

! PAS DE MULTIGRILLE POUR LA VITESSE

!===============================================================================
! 2.  BOUCLES SUR LES NON ORTHOGONALITES
!       (A PARTIR DE LA SECONDE ITERATION)
!===============================================================================

! Application du theta schema

! On calcule le bilan explicite total
thetex = 1.d0 - thetap


! Si THETEX=0, ce n'est pas la peine d'en rajouter
if(abs(thetex).gt.epzero) then
  inc    = 1
! ON POURRAIT METTRE ICCOCG A 0 DANS LES APPELS SUIVANT
  iccocg = 1

  call bilsc4                                                     &
  !==========
 ( nvar   , nscal  ,                                              &
   idtvar , ivar   , iconvp , idiffp , nswrgp , imligp , ircflp , &
   ischcp , isstpp , inc    , imrgra , iccocg , ivisep ,          &
   ippu   , ippv   , ippw   , iwarnp ,                            &
   blencp , epsrgp , climgp , extrap , relaxp , thetex ,          &
   pvar   , pvara  ,                                              &
   coefav , coefbv , cofafv , cofbfv ,                            &
   flumas , flumab , viscfs , viscbs , secvif , secvib ,          &
   smbr   )
endif

!     AVANT DE BOUCLER SUR LES SWEEP, ON STOCKE LE SECOND MEMBRE SANS
!     RECONSTRUCTION DANS LE TABLEAU AUXILIAIRE SMBINI

do iel = 1, ncel
  do isou=1,3
    smbini(isou,iel) = smbr(isou,iel)
  enddo
enddo

!     On initialise sur NCELET pour eviter une communication
!     u*k contient rtpa(*) initialement
do iel = 1, ncelet
  do isou=1,3
    pvar(isou,iel) = pvark(isou,iel)
  enddo
enddo

!     On passe toujours dans bilsc4 avec INC=1
inc = 1

!     Sauf pour les matrices poids (NSWRSP=-1)
if (nswrsp.eq.-1) then
  nswrsp = 1
  inc = 0
endif


!  Attention, pour les matrices poids il faut pouvoir ne faire
!     qu'un seul sweep
nswmod = max( nswrsp, 1 )
do 100 isweep = 1, nswmod

! ---> INCREMENTATION ET RECONSTRUCTION DU SECOND MEMBRE
!      ON NE RECALCULE COCG QU'AU PREMIER PASSAGE (PRESQUE)

  if( isweep.eq.1) then
    iccocg = 1

!  On est entre avec un smb explicite base sur PVARA.
!     si on initialise avec PVAR avec autre chose que PVARA
!     on doit donc corriger SMBR (c'est le cas lorsqu'on itere sur navsto)
    do iel = 1, ncel
      do isou=1,3
        smbini(isou,iel) = smbini(isou,iel)                      &
                   -fimp(isou,1,iel)*(pvar(1,iel) - pvara(1,iel))  &
                   -fimp(isou,2,iel)*(pvar(2,iel) - pvara(2,iel))  &
                   -fimp(isou,3,iel)*(pvar(3,iel) - pvara(3,iel))
        smbr(isou,iel) = smbini(isou,iel)
      enddo
    enddo

  else
    iccocg = 0
    do iel = 1, ncel
!     SMBINI CONTIENT LES TERMES INSTAT, EN DIV(RHO U) ET SOURCE DE MASSE
!     DU SECOND MEMBRE  MIS A JOUR A CHAQUE SWEEP
      do isou=1,3
        smbini(isou,iel) = smbini(isou,iel)                &
                  - fimp(isou,1,iel)*dpvar(1,iel)           &
                  - fimp(isou,2,iel)*dpvar(2,iel)           &
                  - fimp(isou,3,iel)*dpvar(3,iel)
        smbr(isou,iel) = smbini(isou,iel)
      enddo
    enddo
  endif

  call bilsc4                                                     &
  !==========
 ( nvar   , nscal  ,                                              &
   idtvar , ivar   , iconvp , idiffp , nswrgp , imligp , ircflp , &
   ischcp , isstpp , inc    , imrgra , iccocg , ivisep ,          &
   ippu   , ippv   , ippw   , iwarnp ,                            &
   blencp , epsrgp , climgp , extrap , relaxp , thetap ,          &
   pvar   , pvara  ,                                              &
   coefav , coefbv , cofafv , cofbfv ,                            &
   flumas , flumab , viscfs , viscbs , secvif , secvib ,          &
   smbr   )

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

  if (isweep.eq.1) then

    allocate(w1(3, ncelet))  ! Allocate a temporary array

    if(iinvpe.eq.2) then
      iinvpp = 3
    else
      iinvpp = iinvpe
    endif
    call promav(ncelet,ncel,nfac,isym,3,iinvpp,ifacel,dam,xam,pvar,w1)
    !==========
    do iel = 1, ncel
       do isou = 1, 3
          w1(isou,iel) = w1(isou,iel) + smbr(isou,iel)
       enddo
    enddo
    call prodsc(3*ncelet,3*ncel,isqrt,w1(1,1),w1(1,1),rnorm)
    !==========
    rnsmbr(ippu) = rnorm
    rnsmbr(ippv) = rnorm
    rnsmbr(ippw) = rnorm

    deallocate(w1)  ! Free memory

  endif


  ! ---> RESOLUTION IMPLICITE SUR L'INCREMENT DPVAR

  do iel = 1, ncelet
    do isou =1,3
      dpvar(isou,iel) = 0.d0
    enddo
  enddo

  ! iinvpe is useless in the vectorial framework
  iinvpe = 0

  ! Matrix block size
  ibsize = 3

  call invers                                                     &
  !==========
 ( cnom(1), isym   , ibsize , ipol   , ireslq , nitmap , imgrp  , &
   ncymxp , nitmfp ,                                              &
   iwarnp , nfecra , niterf , icycle , iinvpe ,                   &
   epsilp , rnorm  , residu          ,                            &
   dam    , xam    , smbr   , dpvar  )

  nbivar(ippu) = niterf
  nbivar(ippv) = niterf
  nbivar(ippw) = niterf

  if(abs(rnorm)/sqrt(3.d0).gt.epzero) then
    resvar(ippu) = residu/rnorm
    resvar(ippv) = residu/rnorm
    resvar(ippw) = residu/rnorm
  else
    resvar(ippu) = 0.d0
    resvar(ippv) = 0.d0
    resvar(ippw) = 0.d0
  endif


! ---> INCREMENTATION SOLUTION

  do iel = 1, ncelet
    do isou = 1, 3
       pvar (isou,iel) = pvar (isou,iel) + dpvar(isou,iel)
    enddo
  enddo


! ---> TRAITEMENT DU PARALLELISME ET DE LA PERIODICITE

if (irangp.ge.0.or.iperio.eq.1) then
  call synvin(pvar)
  !==========
endif

! ---> TEST DE CONVERGENCE

call prodsc(3*ncelet,3*ncel,isqrt,smbr(1,1),smbr(1,1),residu)

if( residu.le.epsrsp*rnorm         ) then
   if(iwarnp.ge.1) then
      write( nfecra,1000) cnom(1),isweep,residu,rnorm
      write( nfecra,1000) cnom(2),isweep,residu,rnorm
      write( nfecra,1000) cnom(3),isweep,residu,rnorm
   endif
   goto 200
endif

if(iwarnp.ge.3) then
   write(nfecra,1000) cnom(1),isweep,residu,rnorm
   write(nfecra,1000) cnom(2),isweep,residu,rnorm
   write(nfecra,1000) cnom(3),isweep,residu,rnorm
endif

 100  continue

if(iwarnp.ge.2) then
   write(nfecra,1100) cnom(1), nswmod
   write(nfecra,1100) cnom(2), nswmod
   write(nfecra,1100) cnom(3), nswmod
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
    do isou = 1, 3
      smbr(isou,iel) = smbini(isou,iel) - fimp(isou,1,iel)*dpvar(1,iel) &
                                        - fimp(isou,2,iel)*dpvar(2,iel) &
                                        - fimp(isou,3,iel)*dpvar(3,iel)
    enddo
  enddo

  inc    = 1
  iccocg = 1
!     On calcule sans relaxation meme en stationnaire
  idtva0 = 0

  call bilsc4                                                     &
  !==========
 ( nvar   , nscal  ,                                              &
   idtvar , iu     , iconvp , idiffp , nswrgp , imligp , ircflp , &
   ischcp , isstpp , inc    , imrgra , iccocg , ivisep ,          &
   ippu   , ippv   , ippw   , iwarnp ,                            &
   blencp , epsrgp , climgp , extrap , relaxp , thetap ,          &
   pvar   , pvara  ,                                              &
   coefav , coefbv , cofafv , cofbfv ,                            &
   flumas , flumab , viscfs , viscbs , secvif , secvib ,          &
   smbr   )

!     CONTRIBUTION DES NORMES L2 DES DIFFERENTES COMPOSANTES
!       DANS LE TABLEAU ESWORK

  do iel = 1,ncel
    do isou=1,3
      eswork(isou,iel) = (smbr(isou,iel)/ volume(iel))**2
    enddo
  enddo

endif

! SUPPRESSION DE LA HIERARCHIE DE MAILLAGES

if (imgrp.gt.0) then
  chaine = nomvar(ippu)
  lchain = 16
  call dsmlga(chaine(1:16), lchain)
  !==========
endif

!--------
! FORMATS
!--------

#if defined(_CS_LANG_FR)

 1000 format (                                                          &
 1X,A16,' : CV-DIF-TS',I5,' IT - RES= ',E12.5,' NORME= ', E12.5)
 1100 format (                                                          &
'@                                                            ',/,&
'@ @@ ATTENTION : ',A8 ,' CONVECTION-DIFFUSION-TERMES SOURCES ',/,&
'@    =========                                               ',/,&
'@  Nombre d''iterations maximal ',I10   ,' atteint           ',/,&
'@                                                            '  )

#else

 1000 format (                                                          &
 1X,A16,' : CV-DIF-TS',I5,' IT - RES= ',E12.5,' NORM= ', E12.5)
 1100 format (                                                          &
'@                                                            ',/,&
'@ @@ WARNING: ',A8 ,' CONVECTION-DIFFUSION-SOURCE TERMS      ',/,&
'@    ========                                                ',/,&
'@  Maximum number of iterations ',I10   ,' reached           ',/,&
'@                                                            '  )

#endif

!12345678 : CV-DIF-TS 2000 IT - RES= 1234567890234 NORME= 12345678901234
!ATTENTION 12345678 : NON CONVERGENCE DU SYSTEME CONV-DIFF-TS
!----
! FIN
!----

end subroutine
