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

subroutine invers &
!================

 ( cnom   , idbia0 , idbra0 ,                                     &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nideve , nrdeve , nituse , nrtuse ,                            &
   isym   , ipol   , ireslp , nitmap , imgrp  ,                   &
   ncymxp , nitmfp ,                                              &
   iwarnp , nfecra , niterf , icycle , iinvpe ,                   &
   epsilp , rnorm  , residu ,                                     &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dam    , xam    , smbrp  , vx     ,                            &
   w1     , w2     , w3     , w4     , w5     , w6     ,          &
   rdevel , rtuser , ra     )

!===============================================================================
! FONCTION :
! ----------

! APPEL AUX ROUTINE D'INVERSION DE SYSTEMES LINEAIRE
!  MULTIGRILLE + GRADCO OU JACOBI OU BI-CGSTAB
!  GRADCO
!  JACOBI
!  BI-CGSTAB
!  ON SUPPOSE VX INITIALISE EN ENTREE

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! cnom             ! a  ! <-- ! nom de la variable                             !
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
! ndim             ! i  ! <-- ! spatial dimension                              !
! ncelet           ! i  ! <-- ! number of extended (real + ghost) cells        !
! ncel             ! i  ! <-- ! number of cells                                !
! nfac             ! i  ! <-- ! number of interior faces                       !
! nfabor           ! i  ! <-- ! number of boundary faces                       !
! nfml             ! i  ! <-- ! number of families (group classes)             !
! nprfml           ! i  ! <-- ! number of properties per family (group class)  !
! nnod             ! i  ! <-- ! number of vertices                             !
! lndfac           ! i  ! <-- ! size of nodfac indexed array                   !
! lndfbr           ! i  ! <-- ! size of nodfbr indexed array                   !
! ncelbr           ! i  ! <-- ! number of cells with faces on boundary         !
! nideve, nrdeve   ! i  ! <-- ! sizes of idevel and rdevel arrays              !
! nituse, nrtuse   ! i  ! <-- ! sizes of ituser and rtuser arrays              !
! isym             ! e  ! <-- ! indicateur = 1 matrice sym                     !
!                  !    !     !            = 2 matrice non sym                 !
! ipol             ! e  ! <-- ! degre du polynome pour precond                 !
!                  !    !     !         (0 -> diagonal)                        !
! ireslp           ! e  ! <-- ! indicateur = 0 gradco                          !
!                  !    !     !            = 1 jacobi                          !
!                  !    !     !            = 2 cgstab                          !
! nitmap           ! e  ! <-- ! nombre max d'iter pour resol iterativ          !
! imgrp            ! e  ! <-- ! indicateur = 0 pas de multigrille              !
! ncymxp           ! e  ! <-- ! nombre de cycles max pour multigrille          !
! nitmfp           ! e  ! <-- ! nombre d iter sur maillage fin                 !
! iwarnp           ! i  ! <-- ! verbosity                                      !
! nfecra           ! e  ! <-- ! unite du fichier sortie std                    !
! niterf           ! e  ! --> ! nombre d'iterations effectuees                 !
!                  !    !     !  (non multigrille)                             !
! icycle           ! e  ! --> ! nombre de cycles mgm effectues                 !
! iinvpe           ! e  ! <-- ! indicateur pour annuler les increment          !
!                  !    !     ! en periodicite de rotation (=2) ou             !
!                  !    !     ! pour les echanger normalement de               !
!                  !    !     ! maniere scalaire (=1)                          !
! epsilp           ! r  ! <-- ! precision pour resol iter                      !
! rnorm            ! r  ! <-- ! normalisation du residu                        !
! residu           ! r  ! --> ! residu final non norme                         !
! ifacel(2, nfac)  ! ia ! <-- ! interior faces -> cells connectivity           !
! ifabor(nfabor)   ! ia ! <-- ! boundary faces -> cells connectivity           !
! ifmfbr(nfabor)   ! ia ! <-- ! boundary face family numbers                   !
! ifmcel(ncelet)   ! ia ! <-- ! cell family numbers                            !
! iprfml           ! ia ! <-- ! property numbers per family                    !
!  (nfml, nprfml)  !    !     !                                                !
! ipnfac(nfac+1)   ! ia ! <-- ! interior faces -> vertices index (optional)    !
! nodfac(lndfac)   ! ia ! <-- ! interior faces -> vertices list (optional)     !
! ipnfbr(nfabor+1) ! ia ! <-- ! boundary faces -> vertices index (optional)    !
! nodfbr(lndfbr)   ! ia ! <-- ! boundary faces -> vertices list (optional)     !
! idevel(nideve)   ! ia ! <-> ! integer work array for temporary development   !
! ituser(nituse)   ! ia ! <-> ! user-reserved integer work array               !
! ia(*)            ! ia ! --- ! main integer work array                        !
! xyzcen           ! ra ! <-- ! cell centers                                   !
!  (ndim, ncelet)  !    !     !                                                !
! surfac           ! ra ! <-- ! interior faces surface vectors                 !
!  (ndim, nfac)    !    !     !                                                !
! surfbo           ! ra ! <-- ! boundary faces surface vectors                 !
!  (ndim, nfabor)  !    !     !                                                !
! cdgfac           ! ra ! <-- ! interior faces centers of gravity              !
!  (ndim, nfac)    !    !     !                                                !
! cdgfbo           ! ra ! <-- ! boundary faces centers of gravity              !
!  (ndim, nfabor)  !    !     !                                                !
! xyznod           ! ra ! <-- ! vertex coordinates (optional)                  !
!  (ndim, nnod)    !    !     !                                                !
! volume(ncelet)   ! ra ! <-- ! cell volumes                                   !
! dam(ncelet       ! tr ! <-- ! diagonale (maillage fin si mgm)                !
! xam(nfac,isym    ! tr ! <-- ! extradiagonale (maillage fin si mgm)           !
! smbrp(ncelet     ! tr ! <-- ! second membre (maillage fin si mgm)            !
! vx   (ncelet     ! tr ! <-- ! solution du systeme                            !
! w1,2,3,4,5,6     ! tr ! --- ! auxiliaires de travail                         !
!      (ncelet     !    !     !                                                !
! rdevel(nrdeve)   ! ra ! <-> ! real work array for temporary development      !
! rtuser(nrtuse)   ! ra ! <-> ! user-reserved real work array                  !
! ra(*)            ! ra ! --- ! main real work array                           !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

implicit none

!===============================================================================
! Common blocks
!===============================================================================

include "paramx.f90"
include "optcal.f90"

!===============================================================================

! Arguments

character*8      cnom
integer          idbia0 , idbra0
integer          ndim   , ncelet , ncel   , nfac   , nfabor
integer          nfml   , nprfml
integer          nnod   , lndfac , lndfbr , ncelbr
integer          nideve , nrdeve , nituse , nrtuse
integer          isym   , ipol   , ireslp , nitmap
integer          imgrp  , ncymxp , nitmfp
integer          iwarnp , nfecra
integer          niterf , icycle , iinvpe
double precision epsilp , rnorm  , residu

integer          ifacel(2,nfac) , ifabor(nfabor)
integer          ifmfbr(nfabor) , ifmcel(ncelet)
integer          iprfml(nfml,nprfml)
integer          ipnfac(nfac+1), nodfac(lndfac)
integer          ipnfbr(nfabor+1), nodfbr(lndfbr)
integer          idevel(nideve), ituser(nituse)
integer          ia(*)

double precision xyzcen(ndim,ncelet)
double precision surfac(ndim,nfac), surfbo(ndim,nfabor)
double precision cdgfac(ndim,nfac), cdgfbo(ndim,nfabor)
double precision xyznod(ndim,nnod), volume(ncelet)
double precision dam(ncelet), xam(nfac ,2)
double precision smbrp(ncelet)
double precision vx(ncelet)
double precision w1(ncelet),w2(ncelet),w3(ncelet),w4(ncelet)
double precision w5(ncelet),w6(ncelet)
double precision rdevel(nrdeve), rtuser(nrtuse), ra(*)


! Local variables

integer          idebia, idebra
integer          lnom
integer          iresds, iresas, nitmds, nitmas

!===============================================================================

! INITIALISATIONS

lnom = len(cnom)

icycle = 0
niterf = 0

idebia = idbia0
idebra = idbra0

! RESOLUTION

if( imgrp.eq.1 ) then

  iresds = ireslp
  iresas = ireslp

  nitmds = nitmfp
  nitmas = nitmfp

  call resmgr                                                     &
  !==========
 ( cnom   , lnom   , ncelet , ncel   , nfac   ,                   &
   isym   , iresds , iresas , ireslp , ipol   ,                   &
   ncymxp , nitmds , nitmas , nitmap , iinvpe ,                   &
   iwarnp , icycle , niterf , epsilp , rnorm  , residu ,          &
   ifacel , smbrp  , vx     )

elseif(imgrp.eq.0) then

  if (ireslp.ge.0 .and. ireslp.le. 2) then

    call reslin                                                   &
    !==========
 ( cnom   , lnom   , ncelet , ncel   , nfac   ,                   &
   isym   , ireslp , ipol   , nitmap , iinvpe ,                   &
   iwarnp , niterf , epsilp , rnorm  , residu ,                   &
!                 ------                     ------
   ifacel , dam    , xam    , smbrp  , vx     )
!                                   -----

  else
    write(nfecra,1000) cnom, ireslp
    call csexit (1)
  endif

endif


#if defined(_CS_LANG_FR)

 1000 format('INVERS APPELE POUR ',A8,' AVEC IRESOL = ', I10)

#else

 1000 format('INVERS CALLED FOR ',A8,' WITH IRESOL = ', I10)

#endif

!----
! FIN
!----

return

end subroutine
