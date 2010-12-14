!-------------------------------------------------------------------------------

!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2010 EDF S.A., France

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

subroutine itrmas &
!================

 ( idbia0 , idbra0 ,                                              &
   nvar   , nscal  , nphas  ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &
   init   , inc    , imrgra , iccocg , nswrgp , imligp ,          &
   iphydp , iwarnp , nfecra ,                                     &
   epsrgp , climgp , extrap ,                                     &
   idevel , ituser , ia     ,                                     &
   fextx  , fexty  , fextz  ,                                     &
   pvar   , coefap , coefbp , viscf  , viscb  ,                   &
   viselx , visely , viselz ,                                     &
   flumas , flumab ,                                              &
   dpdx   , dpdy   , dpdz   , dpdxa  , dpdya  , dpdza  ,          &
   rdevel , rtuser , ra     )

!===============================================================================
! FONCTION :
! ----------

! INCREMENTATION DU FLUX DE MASSE A PARTIR DE GRAD(P)
! P = PRESSION, INCREMENT DE PRESSION, DOUBLE INCREMENT DE PRESSION
! grad(P) = GRADIENT FACETTE POUR L'INSTANT

!  .          .               --->     -->
!  m   =  m     -  Visc   grad(P) . n
!   ij     ij          ij       ij   ij


!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! nphas            ! i  ! <-- ! number of phases                               !
! nideve, nrdeve   ! i  ! <-- ! sizes of idevel and rdevel arrays              !
! nituse, nrtuse   ! i  ! <-- ! sizes of ituser and rtuser arrays              !
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
! viscf (nfac)     ! tr ! <-- ! "viscosite" face interne(dt*surf/dist          !
! viscb (nfabor    ! tr ! <-- ! "viscosite" face de bord(dt*surf/dist          !
! viselx(ncelet    ! tr ! <-- ! "viscosite" par cellule  dir x                 !
! visely(ncelet    ! tr ! <-- ! "viscosite" par cellule  dir y                 !
! viselz(ncelet    ! tr ! <-- ! "viscosite" par cellule  dir z                 !
! flumas(nfac)     ! tr ! <-- ! flux de masse aux faces internes               !
! flumab(nfabor    ! tr ! <-- ! flux de masse aux faces de bord                !
! fextx,y,z        ! tr ! <-- ! force exterieure generant la pression          !
!   (ncelet)       !    !     !  hydrostatique                                 !
! dpd.(ncelet      ! tr ! --- ! tableau de travail pour le grad de p           !
! dpdxa,y,z        ! tr ! --- ! tableau de travail pour le grad de p           !
!    (ncelet       !    !     !  avec decentrement amont                       !
! rdevel(nrdeve)   ! ra ! <-> ! real work array for temporary development      !
! rtuser(nrtuse)   ! ra ! <-> ! user-reserved real work array                  !
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

integer          idbia0 , idbra0
integer          nvar   , nscal  , nphas
integer          nideve , nrdeve , nituse , nrtuse
integer          init   , inc    , imrgra , iccocg
integer          nswrgp , imligp
integer          iwarnp , iphydp , nfecra
double precision epsrgp , climgp , extrap

integer          idevel(nideve), ituser(nituse)
integer          ia(*)

double precision pvar(ncelet), coefap(nfabor), coefbp(nfabor)
double precision viscf(nfac), viscb(nfabor)
double precision viselx(ncelet), visely(ncelet), viselz(ncelet)
double precision flumas(nfac), flumab(nfabor)
double precision fextx(ncelet),fexty(ncelet),fextz(ncelet)
double precision dpdx (ncelet),dpdy (ncelet),dpdz (ncelet)
double precision dpdxa(ncelet),dpdya(ncelet),dpdza(ncelet)
double precision rdevel(nrdeve), rtuser(nrtuse), ra(*)

! Local variables

integer          idebia, idebra
integer          ifac, ii, jj, iij, iii
double precision pfac,pip
double precision dpxf  , dpyf  , dpzf  ,surfan, dist
double precision dijpfx, dijpfy, dijpfz
double precision diipbx, diipby, diipbz
double precision dijx  , dijy  , dijz

!===============================================================================


!     C'est penible de passer viscf et viscel qui portent une info
!                similaire



!===============================================================================
! 1. INITIALISATION
!===============================================================================

idebia = idbia0
idebra = idbra0

if( init.ge.1 ) then
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

! ---> TRAITEMENT DU PARALLELISME ET DE LA PERIODICITE

if (irangp.ge.0.or.iperio.eq.1) then
  call synsca(pvar)
  !==========
endif


!===============================================================================
! 2.  INCREMENT DU FLUX DE MASSE SS TECHNIQUE DE RECONSTRUCTION
!===============================================================================

if( nswrgp.le.1 ) then

!     FLUX DE MASSE SUR LES FACETTES FLUIDES

  do ifac = 1, nfac

    ii = ifacel(1,ifac)
    jj = ifacel(2,ifac)

    flumas(ifac) = flumas(ifac) + viscf(ifac)*(pvar(ii) -pvar(jj))

  enddo

!     FLUX DE MASSE SUR LES FACETTES DE BORD

  do ifac = 1, nfabor

    ii = ifabor(ifac)
    pfac = inc*coefap(ifac) +coefbp(ifac)*pvar(ii)

    flumab(ifac) = flumab(ifac) +viscb(ifac)*( pvar(ii) -pfac )

  enddo

endif


!===============================================================================
! 3.  INCREMENTATION DU FLUX DE MASSE AVEC TECHNIQUE DE
!         RECONSTRUCTION SI LE MAILLAGE EST NON ORTHOGONAL
!===============================================================================

if( nswrgp.gt.1 ) then

!     CALCUL DU GRADIENT

  call grdcel                                                     &
  !==========
 ( idebia , idebra ,                                              &
   nphas  ,                                                       &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ipr(1) , imrgra , inc    , iccocg , nswrgp , imligp , iphydp , &
   iwarnp , nfecra , epsrgp , climgp , extrap ,                   &
   idevel , ituser , ia     ,                                     &
   fextx  , fexty  , fextz  ,                                     &
   pvar   , coefap , coefbp ,                                     &
   dpdx   , dpdy   , dpdz   ,                                     &
!        ------   ------   ------
   dpdxa  , dpdya  , dpdza  ,                                     &
   rdevel , rtuser , ra     )

! ---> TRAITEMENT DU PARALLELISME ET DE LA PERIODICITE

if (irangp.ge.0.or.iperio.eq.1) then
  call synvec(viselx, visely, viselz)
  !==========
endif

!     FLUX DE MASSE SUR LES FACETTES FLUIDES

  do ifac = 1, nfac

    ii = ifacel(1,ifac)
    jj = ifacel(2,ifac)

    dpxf = 0.5d0*(viselx(ii)*dpdx(ii) + viselx(jj)*dpdx(jj))
    dpyf = 0.5d0*(visely(ii)*dpdy(ii) + visely(jj)*dpdy(jj))
    dpzf = 0.5d0*(viselz(ii)*dpdz(ii) + viselz(jj)*dpdz(jj))

    iij = idijpf-1+3*(ifac-1)
    dijpfx = ra(iij+1)
    dijpfy = ra(iij+2)
    dijpfz = ra(iij+3)

!---> DIJ = IJ - (IJ.N) N
    dijx = (xyzcen(1,jj)-xyzcen(1,ii))-dijpfx
    dijy = (xyzcen(2,jj)-xyzcen(2,ii))-dijpfy
    dijz = (xyzcen(3,jj)-xyzcen(3,ii))-dijpfz

    dist   = ra(idist-1+ifac)
    surfan = ra(isrfan-1+ifac)

    flumas(ifac) = flumas(ifac)                                   &
     + viscf(ifac)*( pvar(ii) -pvar(jj) )                         &
     + ( dpxf * dijx                                              &
     +   dpyf * dijy                                              &
     +   dpzf * dijz )*surfan/dist

  enddo

!     FLUX DE MASSE SUR LES FACETTES DE BORD

  do ifac = 1, nfabor

    ii = ifabor(ifac)
    iii = idiipb-1+3*(ifac-1)
    diipbx = ra(iii+1)
    diipby = ra(iii+2)
    diipbz = ra(iii+3)

    pip = pvar(ii)                                                &
        + dpdx(ii)*diipbx                                         &
        + dpdy(ii)*diipby + dpdz(ii)*diipbz
    pfac = inc*coefap(ifac) +coefbp(ifac)*pip

    flumab(ifac) = flumab(ifac) +viscb(ifac)*( pip -pfac )

  enddo

endif

!--------
! FORMATS
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
