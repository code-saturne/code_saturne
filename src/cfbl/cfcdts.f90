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

subroutine cfcdts &
!================

 ( nvar   , nscal  ,                                              &
   ivar   , iconvp , idiffp , ireslp , ndircp , nitmap ,          &
   imrgra , nswrsp , nswrgp , imligp , ircflp ,                   &
   ischcp , isstpp , iescap ,                                     &
   imgrp  , ncymxp , nitmfp , ipp    , iwarnp ,                   &
   blencp , epsilp , epsrsp , epsrgp , climgp , extrap , thetap , &
   pvara  , coefap , coefbp , cofafp , cofbfp , flumas , flumab , &
   viscfm , viscbm , viscfs , viscbs ,                            &
   rovsdt , smbrp  , pvar   ,                                     &
   eswork )

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
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
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
! eswork(ncelet)   ! ra ! <-- ! prediction-stage error estimator (iescap > 0)  !
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
use cstnum
use entsor
use parall
use period
use mesh

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal
integer          ivar   , iconvp , idiffp , ndircp
integer          nitmap
integer          imrgra , nswrsp , nswrgp , imligp , ircflp
integer          ischcp , isstpp , iescap , imgrp
integer          ncymxp , nitmfp
integer          ipp    , iwarnp
double precision blencp , epsilp , epsrgp , climgp , extrap
double precision thetap , epsrsp

double precision pvara(ncelet), coefap(nfabor), coefbp(nfabor)
double precision                cofafp(nfabor), cofbfp(nfabor)
double precision flumas(nfac), flumab(nfabor)
double precision viscfm(nfac), viscbm(nfabor)
double precision viscfs(nfac), viscbs(nfabor)
double precision rovsdt(ncelet), smbrp(ncelet)
double precision pvar(ncelet)
double precision eswork(ncelet)
integer          imucpp

! Local variables

character*80     chaine
character*16     cnom
integer          isym,ireslp,ireslq,ipol,isqrt
integer          inc,isweep,niterf,iccocg,iel,icycle,nswmod
integer          itenso,iinvpe, iinvpp, ibsize, iesize

double precision residu,rnorm

double precision rvoid(1)

double precision, allocatable, dimension(:) :: dam
double precision, allocatable, dimension(:,:) :: xam
double precision, allocatable, dimension(:) :: dpvar, smbini, w1

!===============================================================================

!===============================================================================
! 1.  INITIALISATIONS
!===============================================================================

! Allocate temporary arrays
allocate(dam(ncelet))
allocate(xam(nfac,2))
allocate(dpvar(ncelet), smbini(ncelet))

! NOMS
chaine = nomvar(ipp)
cnom   = chaine(1:16)

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
  ireslq = mod(ireslp+10000,1000)
  ipol   = (ireslp-ireslq)/1000
endif


! PRISE DE SQRT DANS PS
isqrt = 1

! PRISE EN COMPTE DE LA PERIODICITE

!    Initialisation pour test avant promav
iinvpe = 0

if (irangp.ge.0 .or. iperio.eq.1) then

!    Par defaut, toutes les periodicites seront traitees,
!      les variables etant assimilees a des scalaires (meme si ce sont
!      des composantes de tenseur)
  itenso = 0

  iinvpe = 1

  if (ivar.eq.ir11.or.ivar.eq.ir12.or.             &
      ivar.eq.ir13.or.ivar.eq.ir22.or.             &
      ivar.eq.ir23.or.ivar.eq.ir33) then

    !    Pour les tensions de Reynolds
    !      seules seront echangees les informations sur les faces periodiques
    !      de translation ; on ne touche pas aux informations
    !      relatives aux faces de periodicite de rotation.
    itenso = 1

    !      Lors de la resolution par increments, on echangera egalement les
    !      informations relatives aux faces de periodicite de translation.
    !      Pour les faces de periodicite de rotation, l'increment sera
    !      annule en appelant syncmp au lieu de synsca (iinvpe=2).
    iinvpe = 2

  endif

endif

!===============================================================================
! 1.  CONSTRUCTION MATRICE "SIMPLIFIEE" DE RESOLUTION
!===============================================================================
imucpp = 0

call matrix &
!==========
 ( ncelet , ncel   , nfac   , nfabor ,                            &
   iconvp , idiffp , ndircp , isym   , nfecra ,                   &
   thetap , imucpp ,                                              &
   ifacel , ifabor ,                                              &
   coefbp , cofbfp , rovsdt , flumas , flumab , viscfm , viscbm , &
   rvoid  , dam    , xam    )

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
 ( nvar   , nscal  ,                                              &
   ivar   , iconvp , idiffp , nswrgp , imligp , ircflp ,          &
   ischcp , isstpp , inc    , imrgra , iccocg ,                   &
   ipp    , iwarnp ,                                              &
   blencp , epsrgp , climgp , extrap ,                            &
   pvar   , coefap , coefbp , cofafp , cofbfp ,                   &
   flumas , flumab , viscfs , viscbs ,                            &
   smbrp  )

  call prodsc(ncel,isqrt,smbrp,smbrp,residu)

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
    ! Allocate a temporary array
    allocate(w1(ncelet))
    if(iinvpe.eq.2) then
      iinvpp = 3
    else
      iinvpp = iinvpe
    endif
    ibsize = 1
    iesize = 1
    call promav(isym,ibsize, iesize,iinvpp,dam,xam,pvar,w1)
    do iel = 1, ncel
      w1(iel) = w1(iel) + smbrp(iel)
    enddo
    call prodsc(ncel,isqrt,w1,w1,rnorm)
    rnsmbr(ipp) = rnorm
    ! Free memory
    deallocate(w1)
  endif

! ---> RESOLUTION IMPLICITE SUR L'INCREMENT DPVAR

  do iel = 1, ncel
    dpvar(iel) = 0.d0
  enddo

  ibsize = 1
  iesize = 1

  call invers &
  !==========
 ( cnom   , isym   , ibsize , iesize ,                            &
   ipol   , ireslq , nitmap , imgrp  ,                            &
   ncymxp , nitmfp ,                                              &
   iwarnp , nfecra , niterf , icycle , iinvpe ,                   &
   epsilp , rnorm  , residu ,                                     &
   dam    , xam    , smbrp  , dpvar  )


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

! ---> Handle parallelism and periodicity
!      (periodicity of rotation is not ensured here)

  if (irangp.ge.0 .or. iperio.eq.1) then
    if (itenso.eq.0) then
      call synsca (pvar)
      !==========
    else if (itenso.eq.1) then
      call syncmp (pvar)
      !==========
    endif
  endif

! ---> TEST DE CONVERGENCE

call prodsc(ncel,isqrt,smbrp,smbrp,residu)

if( residu.le.epsrsp*rnorm ) then
   if(iwarnp.ge.1) then
      write(nfecra,1000) cnom,isweep,residu,rnorm
   endif
   goto 200
endif

if(iwarnp.ge.3) then
   write(nfecra,1000) cnom,isweep,residu,rnorm
endif

 100  continue

if(iwarnp.ge.2) then
   write(nfecra,1100) cnom, nswmod
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
 ( nvar   , nscal  ,                                              &
   ivar   , iconvp , idiffp , nswrgp , imligp , ircflp ,          &
   ischcp , isstpp , inc    , imrgra , iccocg ,                   &
   ipp    , iwarnp ,                                              &
   blencp , epsrgp , climgp , extrap ,                            &
   pvar   , coefap , coefbp , cofafp , cofbfp ,                   &
   flumas , flumab , viscfs , viscbs ,                            &
   smbrp  )

!     CONTRIBUTION DES NORMES L2 DES DIFFERENTES COMPOSANTES
!       DANS LE TABLEAU ESWORK

  do iel = 1,ncel
    eswork(iel) = (smbrp(iel)/ volume(iel))**2
  enddo

endif

! Free memory
deallocate(dam, xam)
deallocate(dpvar, smbini)

!--------
! FORMATS
!--------

 1000 format (                                                          &
 1X,A16,' : CV-DIF-TS',I5,' IT - RES= ',E12.5,' NORME= ', E12.5)
 1100 format (                                                          &
'@                                                            ',/,&
'@ @@ ATTENTION : ',A16 ,' CONVECTION-DIFFUSION-TERMES SOURCES ',/,&
'@    =========                                               ',/,&
'@  Nombre d''iterations maximal ',I10   ,' atteint           ',/,&
'@                                                            '  )

!12345678 : CV-DIF-TS 2000 IT - RES= 1234567890234 NORME= 12345678901234
!ATTENTION 12345678 : NON CONVERGENCE DU SYSTEME CONV-DIFF-TS
!----
! FIN
!----

end subroutine
