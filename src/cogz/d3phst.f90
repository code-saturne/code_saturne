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

subroutine d3phst &
!================

 ( ncelet , ncel   , indpdf ,                                     &
   dirmin , dirmax , fdeb   , ffin   , hrec   ,                   &
   fm     , hm     ,                                              &
   hstoe  )

!===============================================================================
!  FONCTION  :
!  ---------

! ROUTINE PHYSIQUE PARTICULIERE : FLAMME DE DIFFUSION
! CALCUL DE L'ENTHALPIE STOECHIOMETRIQUE LOCALE

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! ncelet           ! i  ! <-- ! number of extended (real + ghost) cells        !
! ncel             ! i  ! <-- ! number of cells                                !
! indpdf           ! te ! <-- ! indicateur passage ou non par les pdf          !
! dirmin           ! tr ! <-- ! pdf : dirac en fmin                            !
! dirmax           ! tr ! <-- ! pdf : dirac en fmax                            !
! fdeb             ! tr ! <-- ! pdf : abscisse debut rectangle                 !
! ffin             ! tr ! <-- ! pdf : abscisse fin rectangle                   !
! hrec             ! tr ! <-- ! pdf : hauteur rectangle                        !
! fm               ! tr ! <-- ! fraction de melange moyenne                    !
! hm               ! tr ! <-- ! enthalpie massique moyenne                     !
!                  !    !     !  si ecoulement permeatique                     !
! hstoe            ! tr ! <-- ! enthalpie stoechiometrique                     !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================
!===============================================================================

implicit none

!===============================================================================
! Common blocks
!===============================================================================

include "paramx.f90"
include "numvar.f90"
include "optcal.f90"
include "cstphy.f90"
include "cstnum.f90"
include "entsor.f90"
include "parall.f90"
include "ppppar.f90"
include "ppthch.f90"
include "coincl.f90"
include "cpincl.f90"
include "ppincl.f90"

!===============================================================================
! Arguments

integer          ncelet, ncel
integer          indpdf(ncelet)
double precision dirmin(ncelet), dirmax(ncelet)
double precision fdeb(ncelet), ffin(ncelet), hrec(ncelet)
double precision fm(ncelet), hm(ncelet), hstoe(ncelet)


! Local variables

integer          icel
double precision fsir, hhh, hct, f1, f2
double precision epsi

integer          n1, n2
double precision hsmax, hsmin


!===============================================================================

!===============================================================================
! 0. INITIALISATION
!===============================================================================

epsi  = 1.d-6
fsir =  fs(1)

n1 = 0
n2 = 0
hsmin = grand
hsmax =-grand


do icel = 1, ncel

  if ( indpdf(icel) .eq. 0 ) then

!===============================================================================
! 1. DETERMINATION DE HSTOE SANS INTEGRATION
!===============================================================================

    if ( fm(icel).le.fsir .and. fm(icel).gt.epsi ) then
      hstoe(icel) = ( fsir*hm(icel)+(fm(icel)-fsir)*hinoxy )      &
                  / fm(icel)
    elseif( fm(icel).lt.(1.d0-epsi) ) then
      hstoe(icel) = ((1.d0-fsir)*hm(icel)+(fsir-fm(icel))*hinfue) &
                 / (1.d0-fm(icel))
    endif

  else

!===============================================================================
! 2. DETERMINATION DE HSTOE AVEC INTEGRATION
!===============================================================================

    hct = dirmin(icel)*hinoxy+dirmax(icel)*hinfue
    hhh = 0.d0
    if ( hrec(icel).gt.epsi ) then

      if (fdeb(icel).le.fsir) then
        f1 = fdeb(icel)
        f2 = min(fsir,ffin(icel))
        hct = hct + hrec(icel)*                                   &
              (f2-f1)*hinoxy*(2.d0*fsir-f1-f2)/(2.d0*fsir)
        hhh = hhh + hrec(icel)*(f2**2-f1**2)/(2.d0*fsir)
      endif
      if (ffin(icel).gt.fsir) then
        f1 = max( fsir,fdeb(icel))
        f2 = ffin(icel)
        hct = hct + hrec(icel) *                                  &
             (f2-f1)*hinfue*(f2+f1-2.d0*fsir)/(2.d0*(1.d0-fsir))
        hhh = hhh +                                               &
              hrec(icel)*(f2-f1)*(2.d0-f1-f2)/(2.d0*(1.d0-fsir))
      endif
      hstoe(icel) = (hm(icel)-hct)/ hhh

! Clipping a HSTOEA = HH(1)     en max
! Clipping a HSTOEA = HH(NMAXH) em min

      if ( hstoe(icel) .gt. hh(1) ) then
        n1 = n1 + 1
        hsmax = max(hstoe(icel),hsmax)
        hstoe(icel) = hh(1)
      endif

      if ( hstoe(icel) .lt. hh(nmaxh) ) then
        n2 = n2 + 1
        hsmin = min(hstoe(icel),hsmin)
        hstoe(icel) = hh(nmaxh)
      endif

    endif

  endif

enddo

if (irangp.ge.0) then
  call parcpt (n1)
  !==========
  call parcpt (n2)
  !==========
  call parmax (hsmax)
  !==========
  call parmin (hsmin)
  !==========
endif

if ( n1.gt.0 ) then
  write(nfecra,1000) n1,hsmax,hh(1)
endif
if ( n2.gt.0 ) then
  write(nfecra,1001) n2,hsmin,hh(nmaxh)
endif

!----
! FORMATS
!----

 1000   format(1X,' Clipping de HSTOE EN MAX EN ',I8,' POINTS',/, &
         1X,'     Valeur Max : ',G15.7,/,                   &
         1X,'     Valeur De Clipping : ',G15.7,/)
 1001   format(1X,' Clipping de HSTOE EN MIN EN ',I8,' POINTS',/, &
         1X,'     Valeur Max : ',G15.7,/,                   &
         1X,'     Valeur De Clipping : ',G15.7,/)

return
end subroutine

