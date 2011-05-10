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

subroutine alelap &
!================

 ( idbia0 , idbra0 ,                                              &
   nvar   , nscal  ,                                                                                 &
   ia     ,                                                       &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  , viscf  , viscb  ,                            &
   dam    , xam    , drtp   , smbr   , rovsdt ,                   &
   w1     , w2     , w3     , w4     ,                            &
   w5     , w6     , w7     , w8     , w9     ,                   &
   ra     )

!===============================================================================
! FONCTION :
! ----------

! RESOLUTION DE L'EQUATION DE DIFFUSION POUR LA VITESSE DE MAILLAGE
!  EN ALE

!-------------------------------------------------------------------------------
!ARGU                             ARGUMENTS
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! ia(*)            ! ia ! --- ! main integer work array                        !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
! viscf(nfac)      ! tr ! --- ! visc*surface/dist aux faces internes           !
! viscb(nfabor     ! tr ! --- ! visc*surface/dist aux faces de bord            !
! dam(ncelet       ! tr ! --- ! tableau de travail pour matrice                !
! xam(nfac,*)      ! tr ! --- ! tableau de travail pour matrice                !
! drtp(ncelet      ! tr ! --- ! tableau de travail pour increment              !
! smbr(ncelet      ! tr ! --- ! tableau de travail pour sec mem                !
! rovsdt(ncelet    ! tr ! --- ! tableau de travail pour terme instat           !
! w1...9(ncelet    ! tr ! --- ! tableau de travail                             !
! ra(*)            ! ra ! --- ! main real work array                           !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!-------------------------------------------------------------------------------
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens, only: ndimfb
use numvar
use entsor
use optcal
use cstnum
use cstphy
use pointe
use parall
use period
use mesh

!===============================================================================

implicit none

! Arguments

integer          idbia0 , idbra0
integer          nvar   , nscal

integer          ia(*)

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(ndimfb,*)
double precision coefa(ndimfb,*), coefb(ndimfb,*)
double precision viscf(nfac), viscb(nfabor)
double precision dam(ncelet), xam(nfac,2)
double precision drtp(ncelet), smbr(ncelet), rovsdt(ncelet)
double precision w1(ncelet), w2(ncelet), w3(ncelet)
double precision w4(ncelet), w5(ncelet), w6(ncelet)
double precision w7(ncelet), w8(ncelet), w9(ncelet)
double precision ra(*)

! Local variables

integer          idebia, idebra
integer          iel
integer          iclvar, iclvaf, ivar
integer          ipcvmx, ipcvmy, ipcvmz, ii, ipp
integer          iflmas, iflmab
integer          nswrgp, imligp, iwarnp
integer          iconvp, idiffp, ndircp, ireslp
integer          nitmap, nswrsp, ircflp, ischcp, isstpp, iescap
integer          imgrp , ncymxp, nitmfp
double precision blencp, epsilp, epsrgp, climgp, extrap, thetv
double precision epsrsp
double precision relaxp

!===============================================================================

!===============================================================================
! 1. INITIALISATION
!===============================================================================

idebia = idbia0
idebra = idbra0

ipcvmx = ipproc(ivisma(1))
ipcvmy = ipproc(ivisma(2))
ipcvmz = ipproc(ivisma(3))
! Flux de masse necessaire pour l'appel mais non utilise (ICONV=0)
iflmas = ipprof(ifluma(iu))
iflmab = ipprob(ifluma(iu))

if(iwarni(iuma).ge.1) then
  write(nfecra,1000)
endif


!===============================================================================
! 2. RESOLUTION EFFECTIVE DE L'EQUATION DE VITESSE DE MAILLAGE
!===============================================================================

do ii = 1, 3

  if (ii.eq.1) ivar = iuma
  if (ii.eq.2) ivar = ivma
  if (ii.eq.3) ivar = iwma

  ipp = ipprtp(ivar)
  iclvar = iclrtp(ivar,icoef )
  iclvaf = iclrtp(ivar,icoeff)

  if(iwarni(ivar).ge.1) then
    write(nfecra,1100) nomvar(ipp)
  endif

  do iel = 1, ncel
    smbr(iel)   = 0.d0
    rovsdt(iel) = 0.d0
  enddo

  if (ipcvmx.eq.ipcvmy) then
    call viscfa                                                   &
    !==========
 ( idebia , idebra ,                                              &
   imvisf ,                                                       &
   ia     ,                                                       &
   propce(1,ipcvmx),                                              &
   viscf  , viscb  ,                                              &
   ra     )
  else
    call visort                                                   &
    !==========
 ( idebia , idebra ,                                              &
   imvisf ,                                                       &
   ia     ,                                                       &
   propce(1,ipcvmx), propce(1,ipcvmy), propce(1,ipcvmz),          &
   viscf  , viscb  ,                                              &
   ra     )
  endif

  iconvp = iconv (ivar)
  idiffp = idiff (ivar)
  ireslp = iresol(ivar)
  ndircp = ndircl(ivar)
  nitmap = nitmax(ivar)
  nswrsp = nswrsm(ivar)
  nswrgp = nswrgr(ivar)
  imligp = imligr(ivar)
  ircflp = ircflu(ivar)
  ischcp = ischcv(ivar)
  isstpp = isstpc(ivar)
  iescap = 0
  imgrp  = imgr  (ivar)
  ncymxp = ncymax(ivar)
  nitmfp = nitmgf(ivar)
  iwarnp = iwarni(ivar)
  blencp = blencv(ivar)
  epsilp = epsilo(ivar)
  epsrsp = epsrsm(ivar)
  epsrgp = epsrgr(ivar)
  climgp = climgr(ivar)
  extrap = extrag(ivar)
  relaxp = 1.d0
  thetv  = 1.d0

  call codits                                                     &
  !==========
 ( idebia , idebra ,                                              &
   nvar   , nscal  ,                                                                                 &
   idtvar , ivar   , iconvp , idiffp , ireslp , ndircp , nitmap , &
   imrgra , nswrsp , nswrgp , imligp , ircflp ,                   &
   ischcp , isstpp , iescap ,                                     &
   imgrp  , ncymxp , nitmfp , ipp    , iwarnp ,                   &
   blencp , epsilp , epsrsp , epsrgp , climgp , extrap ,          &
   relaxp , thetv  ,                                              &
   ia     ,                                                       &
   rtpa(1,ivar)    , rtpa(1,ivar)    ,                            &
                     coefa(1,iclvar) , coefb(1,iclvar) ,          &
                     coefa(1,iclvaf) , coefb(1,iclvaf) ,          &
                     propfa(1,iflmas), propfb(1,iflmab),          &
   viscf  , viscb  , viscf  , viscb  ,                            &
   rovsdt , smbr   , rtp(1,ivar)     ,                            &
   dam    , xam    , drtp   ,                                     &
   w1     , w2     , w3     , w4     , w5     ,                   &
   w6     , w7     , w8     , w9     ,                            &
   ra     )

enddo
!--------
! FORMATS
!--------

#if defined(_CS_LANG_FR)

 1000    format(/,                                                &
'   ** RESOLUTION DE LA VITESSE DE MAILLAGE                   ',/,&
'      ------------------------------------                   ',/)
 1100    format(/,'           RESOLUTION POUR LA VARIABLE ',A8,/)

#else

 1000    format(/,                                                &
'   ** SOLVING MESH VELOCITY'                                  ,/,&
'      ---------------------'                                  ,/)
 1100    format(/,'           SOLVING VARIABLE ',A8           ,/)

#endif

!----
! FIN
!----

return

end subroutine
