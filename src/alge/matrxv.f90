!-------------------------------------------------------------------------------

!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2011 EDF S.A., France

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

subroutine matrxv &
!================

 ( ncelet , ncel   , nfac   , nfabor ,                            &
   iconvp , idiffp , ndircp , isym   , nfecra ,                   &
   thetap ,                                                       &
   ifacel , ifabor ,                                              &
   coefbu , cofbfu , fimp   , flumas , flumab , viscf  , viscb  , &
   da     , xa     )

!===============================================================================
! FONCTION :
! ----------

! CONSTRUCTION DE LA MATRICE DE CONVECTION UPWIND/DIFFUSION/TS DE LA VITESSE

!     IL EST INTERDIT DE MODIFIER FIMP   ICI


!-------------------------------------------------------------------------------
!ARGU                             ARGUMENTS
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! ncelet           ! i  ! <-- ! number of extended (real + ghost) cells        !
! ncel             ! i  ! <-- ! number of cells                                !
! nfac             ! i  ! <-- ! number of interior faces                       !
! nfabor           ! i  ! <-- ! number of boundary faces                       !
! iconvp           ! e  ! <-- ! indicateur = 1 convection, 0 sinon             !
! idiffp           ! e  ! <-- ! indicateur = 1 diffusion , 0 sinon             !
! ndircp           ! e  ! <-- ! indicateur = 0 si decalage diagonale           !
! isym             ! e  ! <-- ! indicateur = 1 matrice symetrique              !
!                  !    !     !              2 matrice non symetrique          !
! thetap           ! r  ! <-- ! coefficient de ponderation pour le             !
!                  !    !     ! theta-schema (on ne l'utilise pour le          !
!                  !    !     ! moment que pour u,v,w et les scalaire          !
!                  !    !     ! - thetap = 0.5 correspond a un schema          !
!                  !    !     !   totalement centre en temps (mixage           !
!                  !    !     !   entre crank-nicolson et adams-               !
!                  !    !     !   bashforth)                                   !
! ifacel(2,nfac)   ! te ! <-- ! no des elts voisins d'une face intern          !
! ifabor(nfabor)   ! te ! <-- ! no de l'elt voisin d'une face de bord          !
! coefbu(3,3,nfabor! tr ! <-- ! tab b des cl pour la vitesse                   !
! cofbfu(3,3,nfabor! tr ! <-- ! tab b des cl de flux pour la vitesse           !
! flumas(nfac)     ! tr ! <-- ! flux de masse aux faces internes               !
! flumab(nfabor    ! tr ! <-- ! flux de masse aux faces de bord                !
! viscf(nfac)      ! tr ! <-- ! visc*surface/dist aux faces internes           !
! viscb(nfabor     ! tr ! <-- ! visc*surface/dist aux faces de bord            !
! da (3,3,ncelet)  ! tr ! --> ! partie diagonale de la matrice                 !
! xa (*,nfac)      ! tr ! --> ! extra interleaved diagonale de la matrice      !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use parall

!===============================================================================

implicit none

! Arguments

integer          ncelet , ncel   , nfac   , nfabor
integer          iconvp , idiffp , ndircp , isym
integer          nfecra
double precision thetap

integer          ifacel(2,nfac), ifabor(nfabor)
double precision coefbu(3,3,nfabor), fimp(3,3,ncelet), cofbfu(3,3,nfabor)
double precision flumas(nfac), flumab(nfabor)
double precision viscf(nfac), viscb(nfabor)
double precision da(3,3,ncelet),xa(isym,nfac)

! Local variables

integer          ifac,ii,jj,iel, isou, jsou
double precision flui,fluj,epsi

!===============================================================================

!===============================================================================
! 1. INITIALISATION
!===============================================================================

if(isym.ne.1.and.isym.ne.2) then
  write(nfecra,1000) isym
  call csexit (1)
endif

epsi = 1.d-7

do iel = 1, ncel
  do isou = 1, 3
    do jsou = 1, 3
      da(isou,jsou,iel) = fimp(isou,jsou,iel)
    enddo
  enddo
enddo
if(ncelet.gt.ncel) then
  do iel = ncel+1, ncelet
    do isou = 1, 3
      do jsou = 1, 3
        da(isou,jsou,iel) = 0.d0
      enddo
    enddo
  enddo
endif

if(isym.eq.2) then
  do ifac = 1, nfac
    xa(1,ifac) = 0.d0
    xa(2,ifac) = 0.d0
  enddo
else
  do ifac = 1, nfac
    xa(1,ifac) = 0.d0
  enddo
endif

!===============================================================================
! 2.    CALCUL DES TERMES EXTRADIAGONAUX
!===============================================================================

if(isym.eq.2) then

  do ifac = 1,nfac
    flui = 0.5d0*( flumas(ifac) -abs(flumas(ifac)) )
    fluj =-0.5d0*( flumas(ifac) +abs(flumas(ifac)) )
    xa(1,ifac) = thetap*(iconvp*flui -idiffp*viscf(ifac))
    xa(2,ifac) = thetap*(iconvp*fluj -idiffp*viscf(ifac))
  enddo

else

  do ifac = 1,nfac
    flui = 0.5d0*( flumas(ifac) -abs(flumas(ifac)) )
    xa(1,ifac) = thetap*(iconvp*flui -idiffp*viscf(ifac))
  enddo

endif

!===============================================================================
! 3.     CONTRIBUTION DES TERMES X-TRADIAGONAUX A LA DIAGONALE
!===============================================================================

if(isym.eq.2) then

  do ifac = 1,nfac
    ii = ifacel(1,ifac)
    jj = ifacel(2,ifac)
    do isou = 1, 3
      da(isou,isou,ii) = da(isou,isou,ii) -xa(2,ifac)
      da(isou,isou,jj) = da(isou,isou,jj) -xa(1,ifac)
    enddo
  enddo


else

  do ifac = 1,nfac
    ii = ifacel(1,ifac)
    jj = ifacel(2,ifac)
    do isou = 1, 3
      da(isou,isou,ii) = da(isou,isou,ii) -xa(1,ifac)
      da(isou,isou,jj) = da(isou,isou,jj) -xa(1,ifac)
    enddo
  enddo

endif

!===============================================================================
! 4.     CONTRIBUTION DES FACETTES DE BORDS A LA DIAGONALE
!===============================================================================

do ifac = 1,nfabor
  ii = ifabor(ifac)
  flui = 0.5d0*( flumab(ifac) -abs(flumab(ifac)) )
  fluj =-0.5d0*( flumab(ifac) +abs(flumab(ifac)) )
  do isou = 1, 3
    do jsou = 1, 3
      if(isou.eq.jsou) then
        da(isou,jsou,ii) = da(isou,jsou,ii) + thetap*(                    &
                       iconvp*(-fluj + flui*coefbu(isou,jsou,ifac) )      &
                      +idiffp*viscb(ifac)*(1.d0-cofbfu(isou,jsou,ifac))   &
                             )
      else
        da(isou,jsou,ii) = da(isou,jsou,ii) + thetap*(                    &
                       iconvp*( flui*coefbu(isou,jsou,ifac) )             &
                      +idiffp*viscb(ifac)*(-cofbfu(isou,jsou,ifac))       &
                             )
      endif
    enddo
  enddo
enddo


!===============================================================================
! 5.  NON PRESENCE DE PTS DIRICHLET --> LEGER RENFORCEMENT DE LA
!     DIAGONALE POUR DECALER LE SPECTRE DES VALEURS PROPRES
!===============================================================================
!     (si IDIRCL=0, on a force NDIRCP a valoir au moins 1 pour ne pas
!      decaler la diagonale)

if ( ndircp.le.0 ) then
  do iel=1,ncel
    do isou = 1, 3
      da(isou,isou,iel) = (1.d0+epsi)*da(isou,isou,iel)
    enddo
  enddo
endif

!--------
! FORMATS
!--------

#if defined(_CS_LANG_FR)

 1000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET DANS matrxv                           ',/,&
'@    =========                                               ',/,&
'@     APPEL DE matrxv              AVEC ISYM   = ',I10        ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut pas etre execute.                       ',/,&
'@                                                            ',/,&
'@  Contacter l''assistance.                                  ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

#else

 1000 format(                                                           &
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/,&
'@ @@ WARNING: ABORT IN matrxv'                                ,/,&
'@    ========'                                                ,/,&
'@     matrxv CALLED                WITH ISYM   = ',I10        ,/,&
'@'                                                            ,/,&
'@  The calculation will not be run.'                          ,/,&
'@'                                                            ,/,&
'@  Contact support.'                                          ,/,&
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/)

#endif

!----
! FIN
!----

return

end subroutine
