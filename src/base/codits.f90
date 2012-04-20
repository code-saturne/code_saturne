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

!===============================================================================
! Function:
! ---------

!> \file codits.f90
!>
!> \brief This function solves an advection diffusion equation with source terms
!> for one time step for the variable \f$ a \f$.
!>
!> The equation reads:
!>
!> \f[
!> f_s^{imp}(a^{n+1}-a^n)
!> + \divs \left( a^{n+1} \rho \vect {u} - \mu \grad a^{n+1} \right)
!> = Rhs
!> \f]
!>
!> This equation is rewritten as:
!>
!> \f[
!> f_s^{imp} \delta a
!> + \divs \left( \delta a \rho \vect{u} - \mu \grad \delta a \right)
!> = Rhs^1
!> \f]
!>
!> where \f$ \delta a = a^{n+1} - a^n\f$ and
!> \f$ Rhs^1 = Rhs - \divs( a^n \rho \vect{u} - \mu \grad a^n)\f$
!>
!>
!> It is in fact solved with the following iterative process:
!>
!> \f[
!> f_s^{imp} \delta a^k
!> + \divs \left(\delta a^k \rho \vect{u}-\mu\grad\delta a^k \right)
!> = Rhs^k
!> \f]
!>
!> where \f$Rhs^k=Rhs-f_s^{imp}(a^k-a^n)
!> - \divs \left( a^k\rho\vect{u}-\mu\grad a^k \right)\f$
!>
!> Be careful, it is forbidden to modify \f$ f_s^{imp} \f$ here!
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nvar          total number of variables
!> \param[in]     nscal         total number of scalars
!> \param[in]     idtvar        indicateur du schema temporel
!> \param[in]     ivar          index of the current variable
!> \param[in]     iconvp        indicator
!>                               - 1 convection,
!>                               - 0 sinon
!> \param[in]     idiffp        indicator
!>                               - 1 diffusion,
!>                               - 0 sinon
!> \param[in]     ireslp        indicator
!>                               - 0 conjugate gradient
!>                               - 1 jacobi
!>                               - 2 bi-cgstab
!> \param[in]     ndircp        indicator (0 if the diagonal is stepped aside)
!> \param[in]     nitmap        maximum number of iteration to solve
!>                               the iterative process
!> \param[in]     imrgra        indicator
!>                               - 0 iterative gradient
!>                               - 1 least square gradient
!> \param[in]     nswrsp        number of reconstruction sweeps for the
!>                               Right Hand Side
!> \param[in]     nswrgp        number of reconstruction sweeps for the
!>                               gradients
!> \param[in]     imligp        clipping gradient method
!>                               - < 0 no clipping
!>                               - = 0 thank to neighbooring gradients
!>                               - = 1 thank to the mean gradient
!> \param[in]     ircflp        indicator
!>                               - 1 flux reconstruction,
!>                               - 0 otherwise
!> \param[in]     ischcp        indicator
!>                               - 1 centred
!>                               - 0 2nd order
!> \param[in]     isstpp        indicator
!>                               - 1 without slope test
!>                               - 0 with slope test
!> \param[in]     iescap        compute the predictor indicator if 1
!> \param[in]     imgrp         indicator
!>                               - 0 no multi-grid
!>                               - 1 otherwise
!> \param[in]     ipp           index of the variable for post-processing
!> \param[in]     iwarnp        verbosity
!> \param[in]     blencp        fraction of upwinding
!> \param[in]     epsilp        precision pour resol iter
!> \param[in]     epsrgp        relative precision for the gradient
!>                               reconstruction
!> \param[in]     climgp        clipping coeffecient for the computation of
!>                               the gradient
!> \param[in]     extrap        coefficient for extrapolation of the gradient
!> \param[in]     relaxp        coefficient of relaxation
!> \param[in]     thetap        weightening coefficient for the theta-schema,
!>                               - thetap = 0: explicit scheme
!>                               - thetap = 0.5: time-centred
!>                               scheme (mix between Crank-Nicolson and
!>                               Adams-Bashforth)
!>                               - thetap = 1: implicit scheme
!> \param[in]     pvara         variable at the previous time step
!>                               \f$ a^n \f$
!> \param[in]     pvark         variable at the previous sub-iteration
!>                               \f$ a^k \f$.
!>                               If you sub-iter on Navier-Stokes, then
!>                               it allows to initialize by something else than
!>                               pvara (usually pvar=pvara)
!> \param[in]     coefap        boundary condition array for the variable
!>                               (Explicit part)
!> \param[in]     coefbp        boundary condition array for the variable
!>                               (Impplicit part)
!> \param[in]     cofafp        boundary condition array for the diffusion
!>                               of the variable (Explicit part)
!> \param[in]     cofbfp        boundary condition array for the diffusion
!>                               of the variable (Implicit part)
!> \param[in]     flumas        mass flux at interior faces
!> \param[in]     flumab        mass flux at boundary faces
!> \param[in]     viscfm        \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
!>                               at interior faces for the matrix
!> \param[in]     viscbm        \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
!>                               at border faces for the matrix
!> \param[in]     viscfs        \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
!>                               at interior faces for the r.h.s.
!> \param[in]     viscbs        \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
!>                               at border faces for the r.h.s.
!> \param[in]     rovsdt        \f$ f_s^{imp} \f$
!> \param[in]     smbrp         Right hand side \f$ Rhs^k \f$
!> \param[in,out] pvar          current variable
!> \param[out]    eswork        prediction-stage error estimator
!>                              (if iescap > 0)
!_______________________________________________________________________________

subroutine codits &
!================

 ( nvar   , nscal  ,                                              &
   idtvar , ivar   , iconvp , idiffp , ireslp , ndircp , nitmap , &
   imrgra , nswrsp , nswrgp , imligp , ircflp ,                   &
   ischcp , isstpp , iescap ,                                     &
   imgrp  , ncymxp , nitmfp , ipp    , iwarnp ,                   &
   blencp , epsilp , epsrsp , epsrgp , climgp , extrap ,          &
   relaxp , thetap ,                                              &
   pvara  , pvark  ,                                              &
   coefap , coefbp , cofafp , cofbfp , flumas , flumab ,          &
   viscfm , viscbm , viscfs , viscbs ,                            &
   rovsdt , smbrp  , pvar   ,                                     &
   eswork )

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
integer          ipp    , iwarnp
double precision blencp , epsilp , epsrgp , climgp , extrap
double precision relaxp , thetap , epsrsp

double precision pvara(ncelet), pvark(ncelet)
double precision coefap(nfabor), coefbp(nfabor)
double precision cofafp(nfabor), cofbfp(nfabor)
double precision flumas(nfac), flumab(nfabor)
double precision viscfm(nfac), viscbm(nfabor)
double precision viscfs(nfac), viscbs(nfabor)
double precision rovsdt(ncelet), smbrp(ncelet)
double precision pvar(ncelet)
double precision eswork(ncelet)

! Local variables

character*80     chaine
character*16     cnom
integer          lchain
integer          isym,ireslp,ireslq,ipol,isqrt
integer          inc,isweep,niterf,iccocg,iel,icycle,nswmod
integer          idimte,itenso,iinvpe, iinvpp
integer          idtva0
integer          nagmax, npstmg
integer          ibsize

double precision residu,rnorm
double precision thetex

double precision, allocatable, dimension(:) :: dam
double precision, allocatable, dimension(:,:) :: xam
double precision, allocatable, dimension(:) :: dpvar, smbini, w1

!===============================================================================

!===============================================================================
! 0.  Initialization
!===============================================================================

! Allocate temporary arrays
allocate(dam(ncelet), xam(nfac,2))
allocate(dpvar(ncelet), smbini(ncelet))

! Names
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
  ireslq = mod(ireslp,1000)
  ipol   = (ireslp-ireslq)/1000
endif


! PRISE DE SQRT DANS PS
isqrt = 1

! PRISE EN COMPTE DE LA PERIODICITE

!    Initialisation pour test avant promav
itenso = 0
iinvpe = 0

if (iperio.eq.1) then

!    Par defaut, toutes les periodicites seront traitees,
!      les variables etant assimilees a des scalaires (meme si ce sont
!      des composantes de vecteurs ou de tenseur)

  iinvpe = 1

  if (ivar.eq.iu.or.ivar.eq.iv.or.ivar.eq.iw.or.   &
      ivar.eq.ir11.or.ivar.eq.ir12.or.             &
      ivar.eq.ir13.or.ivar.eq.ir22.or.             &
      ivar.eq.ir23.or.ivar.eq.ir33) then

    !    Pour la vitesse et les tensions de Reynolds, et les tpucou
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
! 1.  Building of the "simplified" matrix
!===============================================================================

call matrix                                                       &
!==========
 ( ncelet , ncel   , nfac   , nfabor ,                            &
   iconvp , idiffp , ndircp , isym   , nfecra ,                   &
   thetap ,                                                       &
   ifacel , ifabor ,                                              &
   coefbp , rovsdt , flumas , flumab , viscfm , viscbm ,          &
   dam    , xam    )

!     En stationnaire, on relaxe la diagonale
if (idtvar.lt.0) then
  !$omp parallel do
  do iel = 1, ncel
    dam(iel) = dam(iel)/relaxp
  enddo
endif
!      CREATION DE LA HIERARCHIE DE MAILLAGE SI MULTIGRILLE

if (imgrp.gt.0) then

! --- Creation de la hierarchie de maillages

  chaine = nomvar(ipp)
  iwarnp = iwarni(ivar)
  nagmax = nagmx0(ivar)
  npstmg = ncpmgr(ivar)
  lchain = 16

  call clmlga                                                     &
  !==========
 ( chaine(1:16) ,    lchain ,                                     &
   ncelet , ncel   , nfac   ,                                     &
   isym   , nagmax , npstmg , iwarnp ,                            &
   ngrmax , ncegrm ,                                              &
   rlxp1  ,                                                       &
   dam    , xam    )

endif


!===============================================================================
! 2. Iterative process to handle non orthogonlaities (starting rom the second
! iteration).
!===============================================================================

! Application du theta schema

! On calcule le bilan explicite total
thetex = 1.d0 - thetap


! Si THETEX=0, ce n'est pas la peine d'en rajouter
if (abs(thetex).gt.epzero) then
  inc    = 1
! ON POURRAIT METTRE ICCOCG A 0 DANS LES APPELS SUIVANT
  iccocg = 1
  call bilsc2                                                     &
  !==========
 ( nvar   , nscal  ,                                              &
   idtvar , ivar   , iconvp , idiffp , nswrgp , imligp , ircflp , &
   ischcp , isstpp , inc    , imrgra , iccocg ,                   &
   ipp    , iwarnp ,                                              &
   blencp , epsrgp , climgp , extrap , relaxp , thetex ,          &
   pvara  , pvara  ,  coefap , coefbp , cofafp , cofbfp ,         &
   flumas , flumab , viscfs , viscbs ,                            &
   smbrp  )
endif

!     AVANT DE BOUCLER SUR LES SWEEP, ON STOCKE LE SECOND MEMBRE SANS
!     RECONSTRUCTION DANS LE TABLEAU AUXILIAIRE SMBINI

!$omp parallel do
do iel = 1, ncel
  smbini(iel) = smbrp(iel)
enddo

!     On initialise sur NCELET pour eviter une communication
!$omp parallel do
do iel = 1, ncelet
  pvar(iel)   = pvark(iel)
enddo

!     On passe toujours dans bilsc2 avec INC=1
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
    !$omp parallel do
    do iel = 1, ncel
      smbini(iel) = smbini(iel) -                                 &
                    rovsdt(iel)*(pvar(iel) - pvara(iel))
      smbrp(iel)  = smbini(iel)
    enddo

  else
    iccocg = 0
    !$omp parallel do
    do iel = 1, ncel
!     SMBINI CONTIENT LES TERMES INSTAT, EN DIV(RHO U) ET SOURCE DE MASSE
!     DU SECOND MEMBRE  MIS A JOUR A CHAQUE SWEEP
      smbini(iel) = smbini(iel) - rovsdt(iel)*dpvar(iel)
      smbrp(iel)  = smbini(iel)
    enddo
  endif

  call bilsc2                                                     &
  !==========
 ( nvar   , nscal  ,                                              &
   idtvar , ivar   , iconvp , idiffp , nswrgp , imligp , ircflp , &
   ischcp , isstpp , inc    , imrgra , iccocg ,                   &
   ipp    , iwarnp ,                                              &
   blencp , epsrgp , climgp , extrap , relaxp , thetap ,          &
   pvar   , pvara  , coefap , coefbp , cofafp , cofbfp ,          &
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
    call promav(isym,1,iinvpp,dam,xam,pvar,w1)
    !$omp parallel do
    do iel = 1, ncel
      w1(iel) = w1(iel) + smbrp(iel)
    enddo
    call prodsc(ncel,isqrt,w1,w1,rnorm)
    rnsmbr(ipp) = rnorm
    ! Free memory
    deallocate(w1)
  endif

! ---> RESOLUTION IMPLICITE SUR L'INCREMENT DPVAR

  !$omp parallel do
  do iel = 1, ncel
    dpvar(iel) = 0.d0
  enddo
  ibsize = 1

  call invers                                                     &
  !==========
 ( cnom   , isym   , ibsize ,                                     &
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

  !$omp parallel do
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

if (residu.le.epsrsp*rnorm) then
   if(iwarnp.ge.1) then
      write(nfecra,1000) cnom,isweep,residu,rnorm
   endif
   goto 200
endif

if (iwarnp.ge.3) then
   write(nfecra,1000) cnom,isweep,residu,rnorm
endif

 100  continue

if(iwarnp.ge.2) then
   write(nfecra,1100) cnom, nswmod
endif

!===============================================================================
! 3. After having computed the new value, an estimator is computed for the
! prediction step of the velocity.
!===============================================================================

 200  continue

! ---> TEST DE PASSAGE DANS LE CALCUL

if (iescap.gt.0) then

! ---> CALCUL DE LA CONTRIBUTION COMPOSANTE PAR COMPOSANTE. DE L ESTIMATEUR


!     SMBINI CONTIENT LES TERMES INSTAT ET EN DIV(U) DU SECOND MEMBRE
!     MIS A JOUR A CHAQUE SWEEP,DONC AU DERNIER, POUR KMAX +1, ON A:

  !$omp parallel do
  do iel = 1,ncel
    smbrp(iel) = smbini(iel) - rovsdt(iel)*dpvar(iel)
  enddo

  inc    = 1
  iccocg = 1
!     On calcule sans relaxation meme en stationnaire
  idtva0 = 0

  call bilsc2                                                     &
  !==========
 ( nvar   , nscal  ,                                              &
   idtva0 , ivar   , iconvp , idiffp , nswrgp , imligp , ircflp , &
   ischcp , isstpp , inc    , imrgra , iccocg ,                   &
   ipp    , iwarnp ,                                              &
   blencp , epsrgp , climgp , extrap , relaxp , thetap ,          &
   pvar   , pvara  , coefap , coefbp , cofafp , cofbfp ,          &
   flumas , flumab , viscfs , viscbs ,                            &
   smbrp  )

!     CONTRIBUTION DES NORMES L2 DES DIFFERENTES COMPOSANTES
!       DANS LE TABLEAU ESWORK

  !$omp parallel do
  do iel = 1,ncel
    eswork(iel) = (smbrp(iel) / volume(iel))**2
  enddo

endif

! SUPPRESSION DE LA HIERARCHIE DE MAILLAGES

if (imgrp.gt.0) then
  chaine = nomvar(ipp)
  lchain = 16
  call dsmlga(chaine(1:16), lchain)
  !==========
endif

! Free memory
deallocate(dam, xam)
deallocate(dpvar, smbini)

!--------
! Formats
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

!----
! End
!----

end subroutine
