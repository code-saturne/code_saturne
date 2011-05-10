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

subroutine cfbsc3 &
!================

 ( idbia0 , idbra0 ,                                              &
   nvar   , nscal  ,                                              &
   ivar   , iconvp , idiffp , nswrgp , imligp , ircflp ,          &
   ischcp , isstpp , inc    , imrgra , iccocg ,                   &
   ipp    , iwarnp ,                                              &
   blencp , epsrgp , climgp , extrap ,                            &
   ia     ,                                                       &
   pvar   , coefap , coefbp , cofafp , cofbfp ,                   &
   flumas , flumab , viscf  , viscb  ,                            &
   flvarf , flvarb ,                                              &
   dpdx   , dpdy   , dpdz   , dpdxa  , dpdya  , dpdza  ,          &
   ra     )

!===============================================================================
! FONCTION :
! ---------

! CALCUL DU FLUX DE CONVECTION-DIFFUSION D'UNE VARIABLE AUX FACES

!                    .                    ----->        -->
! FLVARF (FACEij) =  m   PVAR  - Visc   ( grad PVAR )  . n
!                     ij     ij      ij             ij    ij

!                  .                 ----->       -->
! FLVARB (FABi) =  m  PVAR - Visc  ( grad PVAR ) . n
!                   i     i      i             i    i


! CALCUL EN UPWIND

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! ivar             ! e  ! <-- ! numero de la variable                          !
! iconvp           ! e  ! <-- ! indicateur = 1 convection, 0 sinon             !
! idiffp           ! e  ! <-- ! indicateur = 1 diffusion , 0 sinon             !
! nswrgp           ! e  ! <-- ! nombre de sweep pour reconstruction            !
!                  !    !     !             des gradients                      !
! imligp           ! e  ! <-- ! methode de limitation du gradient              !
!                  !    !     !  < 0 pas de limitation                         !
!                  !    !     !  = 0 a partir des gradients voisins            !
!                  !    !     !  = 1 a partir du gradient moyen                !
! ircflp           ! e  ! <-- ! indicateur = 1 rec flux, 0 sinon               !
! ischcp           ! e  ! <-- ! indicateur = 1 centre , 0 2nd order            !
! isstpp           ! e  ! <-- ! indicateur = 1 sans test de pente              !
!                  !    !     !            = 0 avec test de pente              !
! inc              ! e  ! <-- ! indicateur = 0 resol sur increment             !
!                  !    !     !              1 sinon                           !
! imrgra           ! e  ! <-- ! indicateur = 0 gradrc 97                       !
!                  ! e  ! <-- !            = 1 gradmc 99                       !
! iccocg           ! e  ! <-- ! indicateur = 1 pour recalcul de cocg           !
!                  !    !     !              0 sinon                           !
! ipp              ! e  ! <-- ! numero de variable pour post                   !
! iwarnp           ! i  ! <-- ! verbosity                                      !
! blencp           ! r  ! <-- ! 1 - proportion d'upwind                        !
! epsrgp           ! r  ! <-- ! precision relative pour la                     !
!                  !    !     !  reconstruction des gradients 97               !
! climgp           ! r  ! <-- ! coef gradient*distance/ecart                   !
! extrap           ! r  ! <-- ! coef extrap gradient                           !
! ia(*)            ! te ! --- ! macro tableau entier                           !
! pvar (ncelet     ! tr ! <-- ! variable resolue (instant precedent)           !
! coefap, b        ! tr ! <-- ! tableaux des cond lim pour p                   !
!   (nfabor)       !    !     !  sur la normale a la face de bord              !
! cofafp, b        ! tr ! <-- ! tableaux des cond lim pour le flux de          !
!   (nfabor)       !    !     !  diffusion de p                                !
! flumas(nfac)     ! tr ! <-- ! flux de masse aux faces internes               !
! flumab(nfabor    ! tr ! <-- ! flux de masse aux faces de bord                !
! viscf (nfac)     ! tr ! <-- ! visc*surface/dist aux faces internes           !
!                  !    !     !  pour second membre                            !
! viscb (nfabor    ! tr ! <-- ! visc*surface/dist aux faces de bord            !
!                  !    !     !  pour second membre                            !
! flvarf(nfac)     ! tr ! --> ! flux de convection-diffusion                   !
!                  !    !     !  aux faces internes                            !
! flvarb(nfabor    ! tr ! --> ! flux de convection-diffusion                   !
!                  !    !     !  aux faces de bord                             !
! dpdx,y,z         ! tr ! --- ! tableau de travail pour le grad de p           !
!    (ncelet)      !    !     !                                                !
! dpdxa,ya,za      ! tr ! --- ! tableau de travail pour le grad de p           !
!    (ncelet)      !    !     !  avec decentrement amont                       !
! ra(*)            ! ra ! --- ! main real work array                           !
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
use entsor
use parall
use period
use mesh

!===============================================================================

implicit none

! Arguments

integer          idbia0 , idbra0
integer          nvar   , nscal
integer          ivar   , iconvp , idiffp , nswrgp , imligp
integer          ircflp , ischcp , isstpp
integer          inc    , imrgra , iccocg
integer          iwarnp , ipp
double precision blencp , epsrgp , climgp, extrap

integer          ia(*)

double precision pvar (ncelet), coefap(nfabor), coefbp(nfabor)
double precision                cofafp(nfabor), cofbfp(nfabor)
double precision flumas(nfac), flumab(nfabor)
double precision viscf (nfac), viscb (nfabor)
double precision flvarf(nfac), flvarb(nfabor)
double precision dpdx (ncelet),dpdy (ncelet),dpdz (ncelet)
double precision dpdxa(ncelet),dpdya(ncelet),dpdza(ncelet)
double precision ra(*)

! Local variables

character*80     chaine
character*8      cnom
integer          idebia, idebra
integer          ifac,ii,jj,infac,iel, iij, iii
integer          iphydp
double precision pfac,pfacd,pip,pjp,flui,fluj,flux
double precision pif,pjf
double precision dpxf,dpyf,dpzf
double precision dijpfx, dijpfy, dijpfz
double precision diipfx, diipfy, diipfz
double precision djjpfx, djjpfy, djjpfz
double precision diipbx, diipby, diipbz
double precision pnd

!===============================================================================

!===============================================================================
! 1.  INITIALISATION
!===============================================================================

idebia = idbia0
idebra = idbra0

chaine = nomvar(ipp)
cnom   = chaine(1:8)


!===============================================================================
! 2.  CALCUL DU BILAN AVEC TECHNIQUE DE RECONSTRUCTION
!===============================================================================

! ======================================================================
! ---> CALCUL DU GRADIENT DE PVAR
! ======================================================================
!    DPDX sert pour la reconstruction des flux de diffusion
!       (convection en upwind)
!    On doit donc le calculer uniquement si on a de la diffusion
!       et qu'on reconstruit les flux

if( idiffp.ne.0 .and. ircflp.eq.1 ) then

  iphydp = 0
  call grdcel                                                     &
  !==========
 ( idebia , idebra ,                                              &
   ivar   , imrgra , inc    , iccocg , nswrgp , imligp ,  iphydp ,&
   iwarnp , nfecra , epsrgp , climgp , extrap ,                   &
   ia     ,                                                       &
   dpdxa  , dpdxa  , dpdxa  ,                                     &
   pvar   , coefap , coefbp ,                                     &
   dpdx   , dpdy   , dpdz   ,                                     &
!        ------   ------   ------
   dpdxa  , dpdya  , dpdza  ,                                     &
   ra     )

else
  do iel = 1, ncelet
    dpdx(iel) = 0.d0
    dpdy(iel) = 0.d0
    dpdz(iel) = 0.d0
  enddo
endif


! ======================================================================
! ---> ASSEMBLAGE A PARTIR DES FACETTES FLUIDES
! ======================================================================

infac = 0

do ifac = 1, nfac
  flvarf(ifac) = 0.d0
enddo

do ifac = 1, nfabor
  flvarb(ifac) = 0.d0
enddo


!  --> FLUX UPWIND PUR
!  =====================

do ifac = 1, nfac

  ii = ifacel(1,ifac)
  jj = ifacel(2,ifac)

  dijpfx = dijpf(1,ifac)
  dijpfy = dijpf(2,ifac)
  dijpfz = dijpf(3,ifac)

  pnd   = pond(ifac)

! ON RECALCULE A CE NIVEAU II' ET JJ'

  diipfx = cdgfac(1,ifac) - (xyzcen(1,ii) + (1.d0-pnd) * dijpfx)
  diipfy = cdgfac(2,ifac) - (xyzcen(2,ii) + (1.d0-pnd) * dijpfy)
  diipfz = cdgfac(3,ifac) - (xyzcen(3,ii) + (1.d0-pnd) * dijpfz)
  djjpfx = cdgfac(1,ifac) -  xyzcen(1,jj) + pnd  * dijpfx
  djjpfy = cdgfac(2,ifac) -  xyzcen(2,jj) + pnd  * dijpfy
  djjpfz = cdgfac(3,ifac) -  xyzcen(3,jj) + pnd  * dijpfz

  dpxf = 0.5d0*(dpdx(ii) + dpdx(jj))
  dpyf = 0.5d0*(dpdy(ii) + dpdy(jj))
  dpzf = 0.5d0*(dpdz(ii) + dpdz(jj))

  pip = pvar(ii) + ircflp*(dpxf*diipfx+dpyf*diipfy+dpzf*diipfz)
  pjp = pvar(jj) + ircflp*(dpxf*djjpfx+dpyf*djjpfy+dpzf*djjpfz)

  flui = 0.5d0*( flumas(ifac) +abs(flumas(ifac)) )
  fluj = 0.5d0*( flumas(ifac) -abs(flumas(ifac)) )

  pif = pvar(ii)
  pjf = pvar(jj)
  infac = infac+1

  flux = iconvp*( flui*pif +fluj*pjf ) + idiffp*viscf(ifac)*( pip -pjp )

! --- FLVARF(IFAC) : flux de convection-diffusion de la variable
!                    a la face ij

  flvarf(ifac) = flux

enddo


! ======================================================================
! ---> ASSEMBLAGE A PARTIR DES FACETTES DE BORD
! ======================================================================

do ifac = 1, nfabor

  ii = ifabor(ifac)

  diipbx = diipb(1,ifac)
  diipby = diipb(2,ifac)
  diipbz = diipb(3,ifac)

  flui = 0.5d0*( flumab(ifac) +abs(flumab(ifac)) )
  fluj = 0.5d0*( flumab(ifac) -abs(flumab(ifac)) )

  pip = pvar(ii) +ircflp*(dpdx(ii)*diipbx+dpdy(ii)*diipby+dpdz(ii)*diipbz)

  pfac  = inc*coefap(ifac) +coefbp(ifac)*pip
  pfacd = inc*cofafp(ifac) +cofbfp(ifac)*pip

  flux = iconvp*( flui*pvar(ii) +fluj*pfac ) + idiffp*viscb(ifac)*( pip -pfacd )

! --- FLVARB(IFAC) : flux de convection-diffusion de la variable
!                    a la face de bord i

  flvarb(ifac) = flux

enddo

!--------
! FORMATS
!--------


!----
! FIN
!----

return

end subroutine
