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

subroutine itrmas &
!================

 ( init   , inc    , imrgra , iccocg , nswrgp , imligp ,          &
   iphydp , iwarnp , nfecra ,                                     &
   epsrgp , climgp , extrap ,                                     &
   frcxt  ,                                                       &
   pvar   , coefap , coefbp , cofafp , cofbfp , viscf  , viscb  , &
   viselx , visely , viselz ,                                     &
   flumas , flumab )

!===============================================================================
! FONCTION :
! ----------

! INCREMENTATION DU FLUX DE MASSE A PARTIR DE GRAD(P)
! P = PRESSION, INCREMENT DE PRESSION, DOUBLE INCREMENT DE PRESSION
! grad(P) = GRADIENT FACETTE POUR L'INSTANT

!  .      .                --->     -->
!  m   =  m     -  Visc   grad(P) . n
!   ij     ij          ij       ij   ij


!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! init             ! e  ! <-- ! > 0 : initialisation du flux de masse          !
! inc              ! e  ! <-- ! indicateur = 0 resol sur increment             !
!                  !    !     !              1 sinon                           !
! imrgra           ! e  ! <-- ! indicateur = 0 gradrc 97                       !
!                  ! e  ! <-- !            = 1 gradmc 99                       !
! iccocg           ! e  ! <-- ! indicateur = 1 pour recalcul de cocg           !
!                  !    !     !              0 sinon                           !
! nswrgp           ! e  ! <-- ! nombre de sweep pour reconstruction            !
!                  !    !     !             des gradients                      !
! imligp           ! e  ! <-- ! methode de limitation du gradient              !
!                  !    !     !  < 0 pas de limitation                         !
!                  !    !     !  = 0 a partir des gradients voisins            !
!                  !    !     !  = 1 a partir du gradient moyen                !
! iwarnp           ! i  ! <-- ! verbosity                                      !
! iphydp           ! e  ! <-- ! indicateur de prise en compte de la            !
!                  !    !     ! pression hydrostatique                         !
! nfecra           ! e  ! <-- ! unite du fichier sortie std                    !
! epsrgp           ! r  ! <-- ! precision relative pour la                     !
!                  !    !     !  reconstruction des gradients 97               !
! climgp           ! r  ! <-- ! coef gradient*distance/ecart                   !
! extrap           ! r  ! <-- ! coef extrap gradient                           !
! pvar  (ncelet    ! tr ! <-- ! variable (pression)                            !
! coefap, b        ! tr ! <-- ! tableaux des cond lim pour pvar                !
!   (nfabor)       !    !     !  sur la normale a la face de bord              !
! cof*fp           ! tr ! <-- ! tableaux des cond lim pour pvar                !
!   (nfabor)       !    !     !  sur la normale a la face de bord              !
! viscf (nfac)     ! tr ! <-- ! "viscosite" face interne(dt*surf/dist          !
! viscb (nfabor    ! tr ! <-- ! "viscosite" face de bord(dt*surf/dist          !
! viselx(ncelet    ! tr ! <-- ! "viscosite" par cellule  dir x                 !
! visely(ncelet    ! tr ! <-- ! "viscosite" par cellule  dir y                 !
! viselz(ncelet    ! tr ! <-- ! "viscosite" par cellule  dir z                 !
! flumas(nfac)     ! tr ! <-- ! flux de masse aux faces internes               !
! flumab(nfabor    ! tr ! <-- ! flux de masse aux faces de bord                !
! frcxt            ! tr ! <-- ! force exterieure generant la pression          !
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
use pointe
use numvar
use parall
use period
use mesh

!===============================================================================

implicit none

! Arguments

integer          init   , inc    , imrgra , iccocg
integer          nswrgp , imligp
integer          iwarnp , iphydp , nfecra
double precision epsrgp , climgp , extrap


double precision pvar(ncelet), coefap(nfabor), coefbp(nfabor)
double precision cofafp(nfabor), cofbfp(nfabor)
double precision viscf(nfac), viscb(nfabor)
double precision viselx(ncelet), visely(ncelet), viselz(ncelet)
double precision flumas(nfac), flumab(nfabor)
double precision frcxt(3,ncelet)

! Local variables

integer          ifac, ii, jj, iij, iii, ig, it
double precision pfac,pip
double precision dpxf  , dpyf  , dpzf
double precision dijpfx, dijpfy, dijpfz
double precision diipbx, diipby, diipbz
double precision dijx  , dijy  , dijz

double precision rvoid(1)

double precision, allocatable, dimension(:,:) :: grad

!===============================================================================


!     C'est penible de passer viscf et viscel qui portent une info
!                similaire



!===============================================================================
! 1. INITIALISATION
!===============================================================================


if (init.ge.1) then
  !$omp parallel do
  do ifac = 1, nfac
    flumas(ifac) = 0.d0
  enddo
  !$omp parallel do if(nfabor > thr_n_min)
  do ifac = 1, nfabor
    flumab(ifac) = 0.d0
  enddo
elseif(init.ne.0) then
  write(nfecra,1000) init
  call csexit (1)
endif

! ---> TRAITEMENT DU PARALLELISME ET DE LA PERIODICITE

if (irangp.ge.0.or.iperio.eq.1) then
  call synsca(pvar)
  !==========
endif


!===============================================================================
! 2.  INCREMENT DU FLUX DE MASSE SS TECHNIQUE DE RECONSTRUCTION
!===============================================================================

if (nswrgp.le.1) then

  ! Mass flow through interior faces

  do ig = 1, ngrpi
    !$omp parallel do private(ifac, ii, jj)
    do it = 1, nthrdi
      do ifac = iompli(1,ig,it), iompli(2,ig,it)

        ii = ifacel(1,ifac)
        jj = ifacel(2,ifac)

        flumas(ifac) = flumas(ifac) + viscf(ifac)*(pvar(ii) -pvar(jj))

      enddo
    enddo
  enddo

  ! Mass flow though boundary faces

  do ig = 1, ngrpb
    !$omp parallel do private(ifac, ii, pfac) if(nfabor > thr_n_min)
    do it = 1, nthrdb
      do ifac = iomplb(1,ig,it), iomplb(2,ig,it)

        ii = ifabor(ifac)
        pfac = inc*cofafp(ifac) + cofbfp(ifac)*pvar(ii)

        flumab(ifac) = flumab(ifac) + viscb(ifac)*pfac

      enddo
    enddo
  enddo

endif

!===============================================================================
! 3.  INCREMENTATION DU FLUX DE MASSE AVEC TECHNIQUE DE
!         RECONSTRUCTION SI LE MAILLAGE EST NON ORTHOGONAL
!===============================================================================

if (nswrgp.gt.1) then

  ! Allocate a work array for the gradient calculation
  allocate(grad(ncelet,3))

  ! Compute gradient

  call grdpot                                                     &
  !==========
 ( ipr    , imrgra , inc    , iccocg , nswrgp , imligp , iphydp , &
   iwarnp , nfecra , epsrgp , climgp , extrap ,                   &
   rvoid  ,                                                       &
   frcxt  ,                                                       &
   pvar   , coefap , coefbp ,                                     &
   grad   )

  ! Handle parallelism and periodicity

  if (irangp.ge.0.or.iperio.eq.1) then
    call synvec(viselx, visely, viselz)
    !==========
  endif

  ! Mass flow through interior faces

  do ig = 1, ngrpi
    !$omp parallel do private(ifac, ii, jj, dpxf, dpyf, dpzf, &
    !$omp          dijpfx, dijpfy, dijpfz, dijx, dijy, dijz)
    do it = 1, nthrdi
      do ifac = iompli(1,ig,it), iompli(2,ig,it)

        ii = ifacel(1,ifac)
        jj = ifacel(2,ifac)

        dpxf = 0.5d0*(viselx(ii)*grad(ii,1) + viselx(jj)*grad(jj,1))
        dpyf = 0.5d0*(visely(ii)*grad(ii,2) + visely(jj)*grad(jj,2))
        dpzf = 0.5d0*(viselz(ii)*grad(ii,3) + viselz(jj)*grad(jj,3))

        dijpfx = dijpf(1,ifac)
        dijpfy = dijpf(2,ifac)
        dijpfz = dijpf(3,ifac)

        !---> Dij = IJ - (IJ.N) N
        dijx = (xyzcen(1,jj)-xyzcen(1,ii))-dijpfx
        dijy = (xyzcen(2,jj)-xyzcen(2,ii))-dijpfy
        dijz = (xyzcen(3,jj)-xyzcen(3,ii))-dijpfz

        flumas(ifac) =    flumas(ifac) + viscf(ifac)*(pvar(ii) - pvar(jj))    &
                       +    (dpxf *dijx + dpyf*dijy + dpzf*dijz)              &
                          * surfan(ifac)/dist(ifac)

      enddo
    enddo
  enddo

  ! Mass flow though boundary faces

  do ig = 1, ngrpb
    !$omp parallel do private(ifac, ii, diipbx, diipby, diipbz, pip, pfac) &
    !$omp          if(nfabor > thr_n_min)
    do it = 1, nthrdb
      do ifac = iomplb(1,ig,it), iomplb(2,ig,it)

        ii = ifabor(ifac)

        diipbx = diipb(1,ifac)
        diipby = diipb(2,ifac)
        diipbz = diipb(3,ifac)

        pip = pvar(ii) + grad(ii,1)*diipbx + grad(ii,2)*diipby + grad(ii,3)*diipbz
        pfac = inc*cofafp(ifac) + cofbfp(ifac)*pip

        flumab(ifac) = flumab(ifac) + viscb(ifac)*pfac

      enddo
    enddo
  enddo

  ! Free memory
  deallocate(grad)

endif

!--------
! Formats
!--------

#if defined(_CS_LANG_FR)

 1000 format('ITRMAS APPELE AVEC INIT = ',I10)

#else

 1000 format('ITRMAS CALLED WITH INIT = ',I10)

#endif

!----
! FIN
!----

return

end subroutine
