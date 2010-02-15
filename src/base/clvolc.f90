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

subroutine clvolc &
!================

 ( ncelet , ncel   ,                                              &
   volmin , volmax , voltot , volume )

!===============================================================================

!  FONCTION  :
!  --------

!     CALCUL DU VOLUME GEOMETRIQUE DES ELEMENTS
!     FORMULE DE GREEN : 3*VOLUME = SOMME DIV(R) AVEC R = t(X,Y,Z)

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! ncelet           ! i  ! <-- ! number of extended (real + ghost) cells        !
! ncel             ! i  ! <-- ! number of cells                                !
! volmin           ! r  ! --> ! volume de controle minimal                     !
! volmax           ! r  ! --> ! volume de controle maximal                     !
! voltot           ! r  ! --> ! volume total du domaine                        !
! volume           ! tr ! --> ! volume d'un des ncelet elements                !
! (ncelet)         !    !     !                                                !
!__________________.____._____.________________________________________________.

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

implicit none

!===============================================================================
! Common blocks
!===============================================================================

include "paramx.h"
include "optcal.h"
include "entsor.h"
include "period.h"
include "parall.h"

!===============================================================================

! Arguments

integer ncelet, ncel
double precision volmin, volmax, voltot
double precision volume(ncelet)

integer iel, idimte, itenso

integer          ipass
data             ipass /0/
save             ipass

!===============================================================================
! 1. INITIALISATION EFFECTUEE DANS cs_maillage_grd.c
!===============================================================================

ipass = ipass + 1

!===============================================================================
! 2. ON PREND LES MIN ET MAX
!===============================================================================

volmin =  1.d+12
volmax = -1.d+12
voltot = 0.d0

do iel  = 1, ncel

   volmin = min(volmin, volume(iel))
   volmax = max(volmax, volume(iel))
   voltot = voltot + volume(iel)

enddo

!       On communique le volume au halo pour les filtrages des modeles
!       dynamiques de L.E.S.
!       Dans le cas d'un voisinage etendu traite separement, l'appel a
!       PARCVE est fait directement dans la routine CFLITR (ce qui d'ailleurs
!       pourrait etre optimise)
if (irangp.ge.0) then
  call parcom(volume)
  !==========
  call parmin (volmin)
  !==========
  call parmax (volmax)
  !==========
  call parsom (voltot)
  !==========
endif
if(iperio.eq.1) then

  idimte = 0
  itenso = 0

  call percom                                                     &
  !==========
       ( idimte , itenso ,                                        &
         volume , volume , volume ,                               &
         volume , volume , volume ,                               &
         volume , volume , volume )

endif

!     En ALE, on passe plusieurs fois ici.
!     Au premier passage (avant calculs) on ecrit, on teste et on s'arrete
!       si pb.
!     Aux passages suivants, on n'ecrit pas, on teste et on finit le pas
!       de temps si pb.
if (ipass.eq.1) then
  write(nfecra,1000) volmin, volmax, voltot
  if (volmin.le.0.d0) then
    write(nfecra,1002)
    call csexit (1)
  endif
else
  if (volmin.le.0.d0) then
    write(nfecra,1001) volmin, volmax, voltot
    write(nfecra,1002)
    ntmabs = ntcabs
  endif
endif
!===============================================================================
! 5. FIN
!===============================================================================

#if defined(_CS_LANG_FR)

 1000 format(                                                           &
' --- Information sur les volumes                             ',/,&
'       Volume de controle minimal = ',4X,E18.9                ,/,&
'       Volume de controle maximal = ',4X,E18.9                ,/,&
'       Volume total du domaine    = ',4X,E18.9                 )
 1001 format(/,' CLVOLC : VOLUME DE CONTROLE MINIMAL     = ',E18.9,/,   &
         '          VOLUME DE CONTROLE MAXIMAL     = ',E18.9,/,   &
         '          VOLUME TOTAL DU DOMAINE        = ',E18.9,/,/)
 1002 format(/,' CLVOLC : ARRET SUITE A LA DETECTION D''UN',/,    &
         '          VOLUME NEGATIF',/)

#else

 1000 format(                                                           &
' --- Information on the volumes                              ',/,&
'       Minimum control volume      = ',4X,E18.9               ,/,&
'       Maximum control volume      = ',4X,E18.9               ,/,&
'       Total volume for the domain = ',4X,E18.9                 )
 1001 format(/,' CLVOLC : MINIMUM CONTROL VOLUME         = ',E18.9,/,   &
         '          MAXIMUM CONTROL VOLUME         = ',E18.9,/,   &
         '          TOTAL VOLUME FOR THE DOMAIN    = ',E18.9,/,/)
 1002 format(/,' CLVOLC : ABORT DUE TO THE DETECTION OF A ',/,    &
         '          NEGATIVE VOLUME',/)

#endif

return
end subroutine
