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

subroutine projts &
!================

 ( init   , inc    , imrgra , iccocg , nswrgu , imligu ,          &
   iwarnu , nfecra ,                                              &
   epsrgu , climgu ,                                              &
   frcxt  ,                                                       &
   cofbfp ,                                                       &
   flumas , flumab , viscf  , viscb  ,                            &
   viselx , visely , viselz )

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
! cofbfp(nfabor    ! tr ! <-- ! tableaux des cond lim de pression              !
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
use pointe
use mesh

!===============================================================================

implicit none

! Arguments

integer          init   , inc    , imrgra , iccocg
integer          nswrgu , imligu
integer          iwarnu , nfecra
double precision epsrgu , climgu


double precision pnd
double precision frcxt(3,ncelet)
double precision viscf(nfac), viscb(nfabor)
double precision viselx(ncelet), visely(ncelet), viselz(ncelet)
double precision cofbfp(nfabor)
double precision flumas(nfac), flumab(nfabor)

! Local variables

integer          ifac, ii, jj, iii
double precision dijpfx,dijpfy,dijpfz
double precision diipx,diipy,diipz
double precision djjpx,djjpy,djjpz
double precision distbf,surfn

!===============================================================================

!===============================================================================
! 1.  INITIALISATION
!===============================================================================

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

    flumas(ifac) =  flumas(ifac)                                     &
         + viscf(ifac)*(                                             &
           (cdgfac(1,ifac)-xyzcen(1,ii))*frcxt(1, ii)                &
          +(cdgfac(2,ifac)-xyzcen(2,ii))*frcxt(2, ii)                &
          +(cdgfac(3,ifac)-xyzcen(3,ii))*frcxt(3, ii)                &
          -(cdgfac(1,ifac)-xyzcen(1,jj))*frcxt(1, jj)                &
          -(cdgfac(2,ifac)-xyzcen(2,jj))*frcxt(2, jj)                &
          -(cdgfac(3,ifac)-xyzcen(3,jj))*frcxt(3, jj) )

  enddo

!     FLUX DE MASSE SUR LES FACETTES DE BORD

  do ifac = 1, nfabor

    ii = ifabor(ifac)
    surfn = surfbn(ifac)
    distbf = distb(ifac)

    flumab(ifac) = flumab(ifac)+viscb(ifac)*distbf/surfn          &
         *cofbfp(ifac)*(frcxt(1, ii)*surfbo(1,ifac)                  &
         +frcxt(2, ii)*surfbo(2,ifac)+frcxt(3, ii)*surfbo(3,ifac) )

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

    flumas(ifac) =  flumas(ifac)                                        &
         + viscf(ifac)*(                                                &
           (cdgfac(1,ifac)-xyzcen(1,ii))*frcxt(1, ii)                   &
          +(cdgfac(2,ifac)-xyzcen(2,ii))*frcxt(2, ii)                   &
          +(cdgfac(3,ifac)-xyzcen(3,ii))*frcxt(3, ii)                   &
          -(cdgfac(1,ifac)-xyzcen(1,jj))*frcxt(1, jj)                   &
          -(cdgfac(2,ifac)-xyzcen(2,jj))*frcxt(2, jj)                   &
          -(cdgfac(3,ifac)-xyzcen(3,jj))*frcxt(3, jj) )                 &
         +surfn/dist(ifac)*0.5d0*(                                      &
       (djjpx-diipx)*(viselx(ii)*frcxt(1, ii)+viselx(jj)*frcxt(1, jj))  &
      +(djjpy-diipy)*(visely(ii)*frcxt(2, ii)+visely(jj)*frcxt(2, jj))  &
      +(djjpz-diipz)*(viselz(ii)*frcxt(3, ii)+viselz(jj)*frcxt(3, jj)))

  enddo


!     FLUX DE MASSE SUR LES FACETTES DE BORD

  do ifac = 1, nfabor

    ii = ifabor(ifac)
    surfn = surfbn(ifac)
    distbf = distb(ifac)

    flumab(ifac) = flumab(ifac)+viscb(ifac)*distbf/surfn             &
         *cofbfp(ifac)*(frcxt(1, ii)*surfbo(1,ifac)                  &
         +frcxt(2, ii)*surfbo(2,ifac)+frcxt(3, ii)*surfbo(3,ifac) )

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
