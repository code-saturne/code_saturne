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

subroutine projts &
!================

 ( idbia0 , idbra0 ,                                              &
   nvar   , nscal  , nphas  ,                                     &
   init   , inc    , imrgra , iccocg , nswrgu , imligu ,          &
   iwarnu , nfecra ,                                              &
   epsrgu , climgu ,                                              &
   ia     ,                                                       &
   fextx  , fexty  , fextz  ,                                     &
   coefbp ,                                                       &
   flumas , flumab , viscf  , viscb  ,                            &
   viselx , visely , viselz ,                                     &
   ra     )

!===============================================================================
! FONCTION :
! ----------

! PROJECTION SUR LES FACES DES TERMES DE FORCE EXTERIEURE
! GENERANT UNE PRESSION HYDROSTATIQUE
! EN FAIT, LE TERME CALCULE EST : DTij FEXTij.Sij
!                                      ----   -
! ET IL EST AJOUTE AU FLUX DE MASSE.
! LE CALCUL EST FAIT DE MANIERE COMPATIBLE AVEC ITRMAS (POUR LES
! FACES INTERNES) ET DE MANIERE A CORRIGER L'ERREUR SUR LA CL
! DE PRESSION EN PAROI (dP/dn=0 n'est pas adapte en fait)

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
! ia(*)            ! ia ! --- ! main integer work array                        !
! coefbp(nfabor    ! tr ! <-- ! tableaux des cond lim de pression              !
! flumas(nfac)     ! tr ! <-- ! flux de masse aux faces internes               !
! flumab(nfabor    ! tr ! <-- ! flux de masse aux faces de bord                !
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
use mesh

!===============================================================================

implicit none

! Arguments

integer          idbia0 , idbra0
integer          nvar   , nscal  , nphas
integer          init   , inc    , imrgra , iccocg
integer          nswrgu , imligu
integer          iwarnu , nfecra
double precision epsrgu , climgu

integer          ia(*)

double precision pnd
double precision fextx(ncelet),fexty(ncelet),fextz(ncelet)
double precision viscf(nfac), viscb(nfabor)
double precision viselx(ncelet), visely(ncelet), viselz(ncelet)
double precision coefbp(nfabor)
double precision flumas(nfac), flumab(nfabor)
double precision ra(*)

! Local variables

integer          idebia, idebra
integer          ifac, ii, jj, iii
double precision dijpfx,dijpfy,dijpfz
double precision diipx,diipy,diipz
double precision djjpx,djjpy,djjpz
double precision distbf,surfn

!===============================================================================

!===============================================================================
! 1.  INITIALISATION
!===============================================================================

idebia = idbia0
idebra = idbra0



if( init.eq.1 ) then
  do ifac = 1, nfac
    flumas(ifac) = 0.d0
  enddo
  do ifac = 1, nfabor
    flumab(ifac) = 0.d0
  enddo

elseif(init.ne.0) then
  write(nfecra,1000) init
  call csexit(1)
endif

!===============================================================================
! 2.  CALCUL DU FLUX DE MASSE SANS TECHNIQUE DE RECONSTRUCTION
!===============================================================================

if( nswrgu.le.1 ) then

!     FLUX DE MASSE SUR LES FACETTES FLUIDES

  do ifac = 1, nfac

    ii = ifacel(1,ifac)
    jj = ifacel(2,ifac)

    flumas(ifac) =  flumas(ifac)                                  &
         + viscf(ifac)*(                                          &
           (cdgfac(1,ifac)-xyzcen(1,ii))*fextx(ii)                &
          +(cdgfac(2,ifac)-xyzcen(2,ii))*fexty(ii)                &
          +(cdgfac(3,ifac)-xyzcen(3,ii))*fextz(ii)                &
          -(cdgfac(1,ifac)-xyzcen(1,jj))*fextx(jj)                &
          -(cdgfac(2,ifac)-xyzcen(2,jj))*fexty(jj)                &
          -(cdgfac(3,ifac)-xyzcen(3,jj))*fextz(jj) )

  enddo


!     FLUX DE MASSE SUR LES FACETTES DE BORD

  do ifac = 1, nfabor

    ii = ifabor(ifac)
    surfn = surfbn(ifac)
    distbf = distb(ifac)

    flumab(ifac) = flumab(ifac)+viscb(ifac)*distbf/surfn          &
         *(1.d0-coefbp(ifac))*(fextx(ii)*surfbo(1,ifac)           &
         +fexty(ii)*surfbo(2,ifac)+fextz(ii)*surfbo(3,ifac) )

  enddo


else


!     FLUX DE MASSE SUR LES FACETTES FLUIDES

  do ifac = 1, nfac

    ii = ifacel(1,ifac)
    jj = ifacel(2,ifac)

    pnd = pond(ifac)

    dijpfx = dijpf(1,ifac)
    dijpfy = dijpf(2,ifac)
    dijpfz = dijpf(3,ifac)

    surfn = surfan(ifac)

!     calcul de II' et JJ'
    diipx = cdgfac(1,ifac)-xyzcen(1,ii)-(1.d0-pnd)*dijpfx
    diipy = cdgfac(2,ifac)-xyzcen(2,ii)-(1.d0-pnd)*dijpfy
    diipz = cdgfac(3,ifac)-xyzcen(3,ii)-(1.d0-pnd)*dijpfz
    djjpx = cdgfac(1,ifac)-xyzcen(1,jj)+pnd*dijpfx
    djjpy = cdgfac(2,ifac)-xyzcen(2,jj)+pnd*dijpfy
    djjpz = cdgfac(3,ifac)-xyzcen(3,jj)+pnd*dijpfz

    flumas(ifac) =  flumas(ifac)                                  &
         + viscf(ifac)*(                                          &
           (cdgfac(1,ifac)-xyzcen(1,ii))*fextx(ii)                &
          +(cdgfac(2,ifac)-xyzcen(2,ii))*fexty(ii)                &
          +(cdgfac(3,ifac)-xyzcen(3,ii))*fextz(ii)                &
          -(cdgfac(1,ifac)-xyzcen(1,jj))*fextx(jj)                &
          -(cdgfac(2,ifac)-xyzcen(2,jj))*fexty(jj)                &
          -(cdgfac(3,ifac)-xyzcen(3,jj))*fextz(jj) )              &
         +surfn/dist(ifac)*0.5d0*(                                &
       (djjpx-diipx)*(viselx(ii)*fextx(ii)+viselx(jj)*fextx(jj))  &
      +(djjpy-diipy)*(visely(ii)*fexty(ii)+visely(jj)*fexty(jj))  &
      +(djjpz-diipz)*(viselz(ii)*fextz(ii)+viselz(jj)*fextz(jj)))

  enddo


!     FLUX DE MASSE SUR LES FACETTES DE BORD

  do ifac = 1, nfabor

    ii = ifabor(ifac)
    surfn = surfbn(ifac)
    distbf = distb(ifac)

    flumab(ifac) = flumab(ifac)+viscb(ifac)*distbf/surfn          &
         *(1.d0-coefbp(ifac))*(fextx(ii)*surfbo(1,ifac)           &
         +fexty(ii)*surfbo(2,ifac)+fextz(ii)*surfbo(3,ifac) )

  enddo
endif



!--------
! FORMATS
!--------

#if defined(_CS_LANG_FR)

 1000 format('PROJTS APPELE AVEC INIT =',I10)

#else

 1000 format('PROJTS CALLED WITH INIT =',I10)

#endif


!----
! FIN
!----

return

end subroutine
