!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2011 EDF S.A.
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

subroutine inimav &
!================

 ( nvar   , nscal  ,                                              &
   ivar   ,                                                       &
   iflmb0 , init   , inc    , imrgra , iccocg , nswrgu , imligu , &
   iwarnu , nfecra ,                                              &
   epsrgu , climgu , extrau ,                                     &
   rom    , romb   ,                                              &
   vel    ,                                                       &
   coefav , coefbv ,                                              &
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
! ivar             ! e  ! <-- ! variable                                       !
! iflmb0           ! e  ! <-- ! =1 : flux de masse annule sym-paroi            !
! init             ! e  ! <-- ! > 0 : initialisation du flux de masse          !
! inc              ! e  ! <-- ! indicateur = 0 resol sur increment             !
!                  !    !     !              1 sinon                           !
! imrgra           ! e  ! <-- ! indicateur = 0 gradrc 97                       !
!                  ! e  ! <-- !            = 1 gradmc 99                       !
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
! isympa           ! te ! <-- ! zero pour annuler le flux de masse             !
! (nfabor     )    !    !     !(symetries et parois avec cl couplees)          !
!                  !    !     ! un sinon                                       !
! ia(*)            ! ia ! --- ! main integer work array                        !
! rom(ncelet       ! tr ! <-- ! masse volumique aux cellules                   !
! romb(nfabor)     ! tr ! <-- ! masse volumique aux bords                      !
! vel(3,ncelet)    ! tr ! <-- ! vitesse                                        !
! coefav, b        ! tr ! <-- ! tableaux des cond lim pour ux, uy, uz          !
!   (3,nfabor)     !    !     !  sur la normale a la face de bord              !
! coefbv, q        ! tr ! <-- ! tableaux des cond lim pour ux, uy, uz          !
!   (3,3,nfabor)   !    !     !  sur la normale a la face de bord              !
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
integer          ivar
integer          iflmb0 , init   , inc    , imrgra , iccocg
integer          nswrgu , imligu
integer          iwarnu , nfecra
double precision epsrgu , climgu , extrau


double precision rom(ncelet), romb(nfabor)
double precision vel(3,ncelet)
double precision coefav(3,ndimfb)
double precision coefbv(3,3,nfabor)
double precision flumas(nfac), flumab(nfabor)

! Local variables

integer          ifac, ii, jj, iel
integer          iappel, isou, jsou
double precision pfac, pip
double precision dofx,dofy,dofz,pnd
double precision diipbx, diipby, diipbz
logical          ilved

double precision, dimension(:,:), allocatable :: qdm, coefaq
double precision, dimension(:,:,:), allocatable :: grdqdm

allocate(qdm(3,ncelet))
allocate(coefaq(3,ndimfb))

!===============================================================================

!===============================================================================
! 1.  INITIALISATION
!===============================================================================

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
  do isou = 1, 3
    qdm(isou,iel) = rom(iel)*vel(isou,iel)
  enddo
enddo

! ---> TRAITEMENT DU PARALLELISME ET DE LA PERIODICITE

if (irangp.ge.0.or.iperio.eq.1) then
  call synvin(qdm)
  !==========
endif

do ifac =1, nfabor
  do isou = 1, 3
    coefaq(isou,ifac) = romb(ifac)*coefav(isou,ifac)
  enddo
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
    ! Components U, V, W
    do isou = 1, 3
      flumas(ifac) = flumas(ifac) +                                        &
         (pnd*qdm(isou,ii)+(1.d0-pnd)*qdm(isou,jj)) *surfac(isou,ifac)
    enddo
  enddo


!     FLUX DE MASSE SUR LES FACETTES DE BORD

  do ifac = 1, nfabor
    ii = ifabor(ifac)
    ! Components U, V, W
    do isou = 1, 3
      pfac = inc*coefaq(isou,ifac)

      ! coefbv is a matrix
      do jsou = 1, 3
        pfac = pfac + romb(ifac)*coefbv(isou,jsou,ifac)*vel(jsou,ii)
      enddo

      flumab(ifac) = flumab(ifac) + pfac*surfbo(isou,ifac)

    enddo
  enddo
endif


!===============================================================================
! 4.  CALCUL DU FLUX DE MASSE AVEC TECHNIQUE DE RECONSTRUCTION
!        SI LE MAILLAGE EST NON ORTHOGONAL
!===============================================================================

if( nswrgu.gt.1 ) then

  allocate(grdqdm(3,3,ncelet))


!     CALCUL DU GRADIENT SUIVANT de QDM
!     =================================
  ! gradient vectoriel la periodicite est deja traitee
  ilved = .true.

  call grdvec &
  !==========
( ivar   , imrgra , inc    , iccocg , nswrgu , imligu ,          &
  iwarnu , nfecra , epsrgu , climgu , extrau ,                   &
  ilved  ,                                                       &
  qdm    , coefaq , coefbv ,                                     &
  grdqdm )


! ---> FLUX DE MASSE SUR LES FACETTES FLUIDES

  do ifac = 1, nfac

    ii = ifacel(1,ifac)
    jj = ifacel(2,ifac)

    pnd = pond(ifac)

    dofx = dofij(1,ifac)
    dofy = dofij(2,ifac)
    dofz = dofij(3,ifac)

! Termes suivant U, V, W
    do isou = 1, 3
      flumas(ifac) = flumas(ifac)                                   &
! Terme non reconstruit
         +( pnd*qdm(isou,ii) +(1.d0-pnd)*qdm(isou,jj)             &
!  --->
!  --->     ->    -->      ->
! (Grad(rho U ) . OFij ) . Sij
           +0.5d0*( grdqdm(isou,1,ii) +grdqdm(isou,1,jj) )*dofx   &
           +0.5d0*( grdqdm(isou,2,ii) +grdqdm(isou,2,jj) )*dofy   &
           +0.5d0*( grdqdm(isou,3,ii) +grdqdm(isou,3,jj) )*dofz   &
           )*surfac(isou,ifac)
    enddo

  enddo

! ---> FLUX DE MASSE SUR LES FACETTES DE BORD
  do ifac = 1, nfabor

    ii = ifabor(ifac)
    diipbx = diipb(1,ifac)
    diipby = diipb(2,ifac)
    diipbz = diipb(3,ifac)

! SUIVANT U, V, W
    do isou = 1, 3

      pfac = inc*coefaq(isou,ifac)

      ! coefu is a matrix
      do jsou = 1, 3

        pip =  romb(ifac)*vel(jsou,ii)                &
              +grdqdm(jsou,1,ii)*diipbx              &
              +grdqdm(jsou,2,ii)*diipby              &
              +grdqdm(jsou,3,ii)*diipbz

        pfac = pfac +coefbv(isou,jsou,ifac)*pip
      enddo

     flumab(ifac) = flumab(ifac) +pfac*surfbo(isou,ifac)
    enddo

  enddo

! DESALOCATION
deallocate(grdqdm)
deallocate(qdm, coefaq)

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

!--------
! FORMATS
!--------

#if defined(_CS_LANG_FR)

 1000 format('INIMAV APPELE AVEC INIT =',I10)

#else

 1000 format('INIMAV CALLED WITH INIT =',I10)

#endif

!----
! FIN
!----

return

end subroutine
