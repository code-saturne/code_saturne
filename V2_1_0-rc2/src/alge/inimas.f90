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

subroutine inimas &
!================

 ( nvar   , nscal  ,                                              &
   ivar1  , ivar2  , ivar3  , imaspe ,                            &
   iflmb0 , init   , inc    , imrgra , iccocg , nswrgu , imligu , &
   iwarnu , nfecra ,                                              &
   epsrgu , climgu , extrau ,                                     &
   rom    , romb   ,                                              &
   ux     , uy     , uz     ,                                     &
   coefax , coefay , coefaz , coefbx , coefby , coefbz ,          &
   flumas , flumab )

!===============================================================================
! FONCTION :
! ----------
!                                                              -->
! INCREMENTATION DU FLUX DE MASSE A PARTIR DU CHAMP VECTORIEL ROM.U
!  .     .          -->   -->
!  m   = m   +(rom* U ) . n
!   ij    ij         ij    ij


! Pour la reconstruction, grad(rho u) est calcule avec des
! conditions aux limites approchees :
!    COEFA(rho u) = ROMB * COEFA(u)
!    COEFB(rho u) = COEFB (u)

! et pour le flux de masse au bord on ecrit
!  FLUMAB = [ROMB*COEFA(u) + ROMB*COEFB(u)*Ui
!                + COEFB(u)*II'.grad(rho u) ].Sij
! ce qui utilise de petites approximations sur les
! non-orthogonalites (cf. notice)
!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! ivar1            ! e  ! <-- ! variable de la direction 1                     !
! ivar2            ! e  ! <-- ! variable de la direction 2                     !
! ivar3            ! e  ! <-- ! variable de la direction 3                     !
! imaspe           ! e  ! <-- ! suivant l'appel de inimas                      !
!                  !    !     ! = 1 si appel de navsto resolp                  !
!                  !    !     ! = 2 si appel de divrij                         !
! iflmb0           ! e  ! <-- ! =1 : flux de masse annule sym-paroi            !
! init             ! e  ! <-- ! > 0 : initialisation du flux de masse          !
! inc              ! e  ! <-- ! indicateur = 0 resol sur increment             !
!                  !    !     !              1 sinon                           !
! imrgra           ! e  ! <-- ! indicateur = 0 gradrc 97                       !
!                  ! e  ! <-- !            = 1 gradmc 99                       !
! iccocg           ! e  ! <-- ! indicateur = 1 pour recalcul de cocg           !
!                  !    !     !              0 sinon                           !
! nswrgu           ! e  ! <-- ! nombre de sweep pour reconstruction            !
!                  !    !     !             des gradients                      !
! imligu           ! e  ! <-- ! methode de limitation du gradient              !
!                  !    !     !  < 0 pas de limitation                         !
!                  !    !     !  = 0 a partir des gradients voisins            !
!                  !    !     !  = 1 a partir du gradient moyen                !
! iwarnu           ! e  ! <-- ! niveau d'impression                            !
! nfecra           ! e  ! <-- ! unite du fichier sortie std                    !
! epsrgu           ! r  ! <-- ! precision relative pour la                     !
!                  !    !     !  reconstruction des gradients 97               !
! climgu           ! r  ! <-- ! coef gradient*distance/ecart                   !
! extrau           ! r  ! <-- ! coef extrap gradient                           !
! rom(ncelet       ! tr ! <-- ! masse volumique aux cellules                   !
! romb(nfabor)     ! tr ! <-- ! masse volumique aux bords                      !
! ux,y,z(ncelet    ! tr ! <-- ! vitesse                                        !
! coefax, b        ! tr ! <-- ! tableaux des cond lim pour ux, uy, uz          !
!   (nfabor)       !    !     !  sur la normale a la face de bord              !
! flumas(nfac)     ! tr ! <-- ! flux de masse aux faces internes               !
! flumab(nfabor    ! tr ! <-- ! flux de masse aux faces de bord                !
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
use pointe
use parall
use period
use mesh

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal
integer          ivar1  , ivar2  , ivar3  , imaspe
integer          iflmb0 , init   , inc    , imrgra , iccocg
integer          nswrgu , imligu
integer          iwarnu , nfecra
double precision epsrgu , climgu , extrau


double precision rom(ncelet), romb(nfabor)
double precision ux(ncelet), uy(ncelet), uz(ncelet)
double precision coefax(nfabor), coefay(nfabor), coefaz(nfabor)
double precision coefbx(nfabor), coefby(nfabor), coefbz(nfabor)
double precision flumas(nfac), flumab(nfabor)

! Local variables

integer          ifac, ii, jj, iel, iii
integer          iappel

double precision pfac,pip,uxfac,uyfac,uzfac
double precision dofx,dofy,dofz,pnd
double precision diipbx, diipby, diipbz

double precision, allocatable, dimension(:,:) :: grad
double precision, allocatable, dimension(:) :: qdmx, qdmy, qdmz
double precision, allocatable, dimension(:,:) :: coefqa

!===============================================================================

!===============================================================================
! 1.  INITIALISATION
!===============================================================================

! Allocate temporary arrays
allocate(qdmx(ncelet), qdmy(ncelet), qdmz(ncelet))
allocate(coefqa(ndimfb,3))


! ---> CALCUL DE LA QTE DE MOUVEMENT


if( init.eq.1 ) then
  do ifac = 1, nfac
    flumas(ifac) = 0.d0
  enddo
  do ifac = 1, nfabor
    flumab(ifac) = 0.d0
  enddo

elseif(init.ne.0) then
  write(nfecra,1000) init
  call csexit (1)
endif

do iel = 1, ncel
  qdmx(iel) = rom(iel)*ux(iel)
  qdmy(iel) = rom(iel)*uy(iel)
  qdmz(iel) = rom(iel)*uz(iel)
enddo

! ---> TRAITEMENT DU PARALLELISME ET DE LA PERIODICITE

if (irangp.ge.0.or.iperio.eq.1) then
  call synvec(qdmx, qdmy, qdmz)
  !==========
endif

do ifac =1, nfabor
  coefqa(ifac,1) = romb(ifac)*coefax(ifac)
  coefqa(ifac,2) = romb(ifac)*coefay(ifac)
  coefqa(ifac,3) = romb(ifac)*coefaz(ifac)
enddo

!===============================================================================
! 2.  CALCUL DU FLUX DE MASSE SANS TECHNIQUE DE RECONSTRUCTION
!===============================================================================

if( nswrgu.le.1 ) then

!     FLUX DE MASSE SUR LES FACETTES FLUIDES

  do ifac = 1, nfac

    ii = ifacel(1,ifac)
    jj = ifacel(2,ifac)

    pnd = pond(ifac)

    flumas(ifac) =  flumas(ifac)                                  &
     +(pnd*qdmx(ii)+(1.d0-pnd)*qdmx(jj) )*surfac(1,ifac)        &
     +(pnd*qdmy(ii)+(1.d0-pnd)*qdmy(jj) )*surfac(2,ifac)        &
     +(pnd*qdmz(ii)+(1.d0-pnd)*qdmz(jj) )*surfac(3,ifac)

  enddo


!     FLUX DE MASSE SUR LES FACETTES DE BORD

  do ifac = 1, nfabor

    ii = ifabor(ifac)
    uxfac = inc*coefqa(ifac,1) +coefbx(ifac)*romb(ifac)*ux(ii)
    uyfac = inc*coefqa(ifac,2) +coefby(ifac)*romb(ifac)*uy(ii)
    uzfac = inc*coefqa(ifac,3) +coefbz(ifac)*romb(ifac)*uz(ii)

    flumab(ifac) = flumab(ifac)                                   &
     +( uxfac*surfbo(1,ifac)                                      &
                    +uyfac*surfbo(2,ifac) +uzfac*surfbo(3,ifac) )

  enddo

endif


!===============================================================================
! 4.  CALCUL DU FLUX DE MASSE AVEC TECHNIQUE DE RECONSTRUCTION
!        SI LE MAILLAGE EST NON ORTHOGONAL
!===============================================================================

if( nswrgu.gt.1 ) then


  ! Allocate a temporary array for the gradient calculation
  allocate(grad(ncelet,3))


!     TRAITEMENT DE LA PERIODICITE SPEFICIQUE A INIMAS AU DEBUT
!     =========================================================

  if(iperot.gt.0) then
    iappel = 1

    call permas &
    !==========
 ( imaspe , iappel ,                 &
   rom    ,                          &
   dudxy  , drdxy  , wdudxy , wdrdxy )

  endif

!     FLUX DE MASSE SUIVANT X
!     =======================

! ---> CALCUL DU GRADIENT

  call grdcel                                                     &
  !==========
 ( ivar1  , imrgra , inc    , iccocg , nswrgu , imligu ,          &
   iwarnu , nfecra , epsrgu , climgu , extrau ,                   &
   qdmx   , coefqa(1,1) , coefbx ,                                &
   grad   )


! ---> FLUX DE MASSE SUR LES FACETTES FLUIDES

  do ifac = 1, nfac

    ii = ifacel(1,ifac)
    jj = ifacel(2,ifac)

    pnd = pond(ifac)

    dofx = dofij(1,ifac)
    dofy = dofij(2,ifac)
    dofz = dofij(3,ifac)

    flumas(ifac) = flumas(ifac)                                   &
         +( pnd*qdmx(ii) +(1.d0-pnd)*qdmx(jj)                     &
           +0.5d0*( grad(ii,1) +grad(jj,1) )*dofx                     &
           +0.5d0*( grad(ii,2) +grad(jj,2) )*dofy                     &
           +0.5d0*( grad(ii,3) +grad(jj,3) )*dofz    )*surfac(1,ifac)

  enddo

! ---> FLUX DE MASSE SUR LES FACETTES DE BORD

  do ifac = 1, nfabor

    ii = ifabor(ifac)

    diipbx = diipb(1,ifac)
    diipby = diipb(2,ifac)
    diipbz = diipb(3,ifac)

    pip = romb(ifac) * ux(ii)                                     &
          +grad(ii,1)*diipbx                                      &
          +grad(ii,2)*diipby +grad(ii,3)*diipbz
    pfac = inc*coefqa(ifac,1) +coefbx(ifac)*pip

    flumab(ifac) = flumab(ifac) +pfac*surfbo(1,ifac)

  enddo


!     FLUX DE MASSE SUIVANT Y
!     =======================

! ---> CALCUL DU GRADIENT

  call grdcel                                                     &
  !==========
 ( ivar2  , imrgra , inc    , iccocg , nswrgu , imligu ,          &
   iwarnu , nfecra , epsrgu , climgu , extrau ,                   &
   qdmy   , coefqa(1,2) , coefby ,                                &
   grad   )


! ---> FLUX DE MASSE SUR LES FACETTES FLUIDES

  do ifac = 1, nfac

    ii = ifacel(1,ifac)
    jj = ifacel(2,ifac)

    pnd = pond(ifac)

    dofx = dofij(1,ifac)
    dofy = dofij(2,ifac)
    dofz = dofij(3,ifac)

    flumas(ifac) = flumas(ifac)                                   &
         +( pnd*qdmy(ii) +(1.d0-pnd)*qdmy(jj)                     &
           +0.5d0*( grad(ii,1) +grad(jj,1) )*dofx                     &
           +0.5d0*( grad(ii,2) +grad(jj,2) )*dofy                     &
           +0.5d0*( grad(ii,3) +grad(jj,3) )*dofz    )*surfac(2,ifac)

  enddo

! ---> FLUX DE MASSE SUR LES FACETTES DE BORD

  do ifac = 1, nfabor

    ii = ifabor(ifac)

    diipbx = diipb(1,ifac)
    diipby = diipb(2,ifac)
    diipbz = diipb(3,ifac)

    pip = romb(ifac) * uy(ii)                                     &
        +grad(ii,1)*diipbx                                        &
        +grad(ii,2)*diipby +grad(ii,3)*diipbz
    pfac = inc*coefqa(ifac,2) +coefby(ifac)*pip

    flumab(ifac) = flumab(ifac) +pfac*surfbo(2,ifac)

  enddo

!     FLUX DE MASSE SUIVANT Z
!     =======================

! ---> CALCUL DU GRADIENT

  call grdcel                                                     &
  !==========
 ( ivar3  , imrgra , inc    , iccocg , nswrgu , imligu ,          &
   iwarnu , nfecra , epsrgu , climgu , extrau ,                   &
   qdmz   , coefqa(1,3) , coefbz ,                                &
   grad   )

!     FLUX DE MASSE SUR LES FACETTES FLUIDES

  do ifac = 1, nfac

    ii = ifacel(1,ifac)
    jj = ifacel(2,ifac)

    pnd = pond(ifac)

    dofx = dofij(1,ifac)
    dofy = dofij(2,ifac)
    dofz = dofij(3,ifac)

    flumas(ifac) = flumas(ifac)                                   &
         +( pnd*qdmz(ii) +(1.d0-pnd)*qdmz(jj)                     &
           +0.5d0*( grad(ii,1) +grad(jj,1) )*dofx                     &
           +0.5d0*( grad(ii,2) +grad(jj,2) )*dofy                     &
           +0.5d0*( grad(ii,3) +grad(jj,3) )*dofz    )*surfac(3,ifac)

  enddo

! ---> FLUX DE MASSE SUR LES FACETTES DE BORD

  do ifac = 1, nfabor

    ii = ifabor(ifac)

    diipbx = diipb(1,ifac)
    diipby = diipb(2,ifac)
    diipbz = diipb(3,ifac)

    pip = romb(ifac) * uz(ii)                                     &
        +grad(ii,1)*diipbx                                        &
        +grad(ii,2)*diipby +grad(ii,3)*diipbz
    pfac = inc*coefqa(ifac,3) +coefbz(ifac)*pip

    flumab(ifac) = flumab(ifac) +pfac*surfbo(3,ifac)

  enddo

  ! Free memory
  deallocate(grad)

!     TRAITEMENT DE LA PERIODICITE SPEFICIQUE A INIMAS A LA FIN
!     =========================================================

  if(iperot.gt.0) then
    iappel = 2

    call permas &
    !==========
 ( imaspe , iappel ,                 &
   rom    ,                          &
   dudxy  , drdxy  , wdudxy , wdrdxy )

  endif




endif


!===============================================================================
! 6.  POUR S'ASSURER DE LA NULLITE DU FLUX DE MASSE AUX LIMITES
!       SYMETRIES PAROIS COUPLEES
!===============================================================================

if(iflmb0.eq.1) then
! FORCAGE DE FLUMAB a 0 pour la vitesse'
  do ifac = 1, nfabor
    if(isympa(ifac).eq.0) then
      flumab(ifac) = 0.d0
    endif
  enddo
endif

! Free memory
deallocate(qdmx, qdmy, qdmz)
deallocate(coefqa)

!--------
! FORMATS
!--------

#if defined(_CS_LANG_FR)

 1000 format('INIMAS APPELE AVEC INIT =',I10)

#else

 1000 format('INIMAS CALLED WITH INIT =',I10)

#endif

!----
! FIN
!----

return

end subroutine
