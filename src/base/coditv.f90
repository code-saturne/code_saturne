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

!> \file coditv.f90
!>
!> \brief This function solves an advection diffusion equation with source terms
!> for one time step for the vector variable \f$ \vect{a} \f$.
!>
!> The equation reads:
!>
!> \f[
!> \tens{f_s}^{imp}(\vect{a}^{n+1}-\vect{a}^n)
!> + \divt \left( \vect{a}^{n+1} \otimes \rho \vect {u}
!>              - \mu \gradv \vect{a}^{n+1}\right)
!> = \vect{Rhs}
!> \f]
!>
!> This equation is rewritten as:
!>
!> \f[
!> \tens{f_s}^{imp} \delta \vect{a}
!> + \divt \left( \delta \vect{a} \otimes \rho \vect{u}
!>              - \mu \gradv \delta \vect{a} \right)
!> = \vect{Rhs}_1
!> \f]
!>
!> where \f$ \delta \vect{a} = \vect{a}^{n+1} - \vect{a}^n\f$ and
!> \f$ \vect{Rhs}_1 = \vect{Rhs}
!> - \divt \left( \vect{a}^n \otimes \rho \vect{u}
!>              - \mu \gradt \vect{a}^n \right)\f$
!>
!>
!> It is in fact solved with the following iterative process:
!>
!> \f[
!> \tens{f_s}^{imp} \delta \vect{a}^k
!> + \divt \left( \delta \vect{a}^k \otimes \rho \vect{u}
!>              - \mu \gradv \delta a^k \right)
!> = \vect{Rhs}^k
!> \f]
!>
!> where \f$ \vect{Rhs}^k = \vect{Rhs}
!> - \tens{f_s}^{imp} \left(\vect{a}^k-\vect{a}^n \right)
!> - \divv \left( \vect{a}^k \otimes \rho \vect{u}
!>              - \mu \gradt \vect{a}^k \right)\f$
!>
!> Be careful, it is forbidden to modify \f$ \tens{f_s}^{imp} \f$ here!
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
!> \param[in]     imrgra        indicateur
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
!> \param[in]     ivisep        indicator to take \f$ \divv
!>                               \left(\mu \gradt \transpose{\vect{a}} \right)
!>                               -2/3 \grad\left( \mu \dive \vect{a} \right)\f$
!>                               - 1 take into account,
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
!> \param[in]     ipp*          index of the variable for post-processing
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
!>                               \f$ \vect{a}^n \f$
!> \param[in]     pvark         variable at the previous sub-iteration
!>                               \f$ \vect{a}^k \f$.
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
!> \param[in]     fimp          \f$ \tens{f_s}^{imp} \f$
!> \param[in]     smbrp         Right hand side \f$ \vect{Rhs}^k \f$
!> \param[in,out] pvar          current variable
!> \param[out]    eswork        prediction-stage error estimator
!>                              (if iescap > 0)
!_______________________________________________________________________________

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
integer          inc,isweep,niterf,iel,icycle,nswmod
integer          iinvpe,iinvpp
integer          idtva0
integer          nagmax, npstmg
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

  call bilsc4                                                     &
  !==========
 ( nvar   , nscal  ,                                              &
   idtvar , ivar   , iconvp , idiffp , nswrgp , imligp , ircflp , &
   ischcp , isstpp , inc    , imrgra , ivisep ,                   &
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
   ischcp , isstpp , inc    , imrgra , ivisep ,                   &
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
    call promav(isym,3,iinvpp,dam,xam,pvar,w1)
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
!     On calcule sans relaxation meme en stationnaire
  idtva0 = 0

  call bilsc4                                                     &
  !==========
 ( nvar   , nscal  ,                                              &
   idtvar , iu     , iconvp , idiffp , nswrgp , imligp , ircflp , &
   ischcp , isstpp , inc    , imrgra , ivisep ,                   &
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
