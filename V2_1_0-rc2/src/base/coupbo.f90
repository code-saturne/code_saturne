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

subroutine coupbo &
!================

 ( nvar   , nscal  , isvtb  ,                                     &
   ncp , ncv , ientha ,                                           &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   cpcst  , cp     , cvcst  , cv     ,                            &
   hbord  , tbord  )

!===============================================================================
! Purpose:
! --------

! Send data relative to a SYRTHES coupling

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nfabor           ! i  ! <-- ! number of boundary faces                       !
! ncelet           ! i  ! <-- ! number of extended (real + ghost) cells        !
! ncp              ! i  ! <-- ! dimension de cp (ncelet ou 1)                  !
! ncv              ! i  ! <-- ! dimension de cv (ncelet ou 1)                  !
! nfabor           ! i  ! <-- ! number of boundary faces                       !
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! ientha           ! i  ! <-- ! 1 si tparoi est une enthalpie                  !
!                  ! i  ! <-- ! 2 si tparoi est une energie                    !
!                  !    !     !    (compressible)                              !
! cpcst            ! r  ! <-- ! chaleur specifique si constante                !
! cvcst            ! r  ! <-- ! chaleur specifique si constante                !
! cp(ncp)          ! ra ! <-- ! chaleur specifique si variable                 !
! cv(ncp)          ! ra ! <-- ! chaleur specifique si variable                 !
! hbord(nfabor)    ! ra ! <-- ! coefficients d'echange aux bords               !
! tbord(nfabor)    ! ra ! <-- ! temperatures aux bords                         !
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use entsor
use cstphy
use mesh

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal
integer          isvtb
integer          ncp    , ncv    , ientha

double precision cpcst  , cvcst

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*),propfa(nfac,*),propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision cp(ncp), cv(ncv)
double precision hbord(nfabor),tbord(nfabor)

! Local variables

integer          nbccou, inbcou, inbcoo, nbfcou, ifac, iloc, iel
integer          itflui, ihparo
integer          mode
integer          iccfth, imodif
integer          iepsel, iepsfa, igamag, ixmasm, ifinwa
double precision enthal, temper, energ, cvt

integer, dimension(:), allocatable :: lfcou
double precision, dimension(:), allocatable :: tfluid, hparoi, wa

!===============================================================================

! Initialize variables to avoid compiler warnings

iepsel = 0
iepsfa = 0
igamag = 0
ixmasm = 0
ifinwa = 0

!===============================================================================
! SYRTHES coupling: compute fluid temperature and exchange coefficient
!===============================================================================

! Get number of coupling cases

call nbcsyr(nbccou)
!==========

!---> Loop on couplings

do inbcou = 1, nbccou

  ! Number of boundary faces per coupling case
  inbcoo = inbcou
  call nbfsyr(inbcoo, nbfcou)
  !==========

  ! Memory management to build arrays
  allocate(lfcou(nbfcou))
  allocate(tfluid(nbfcou))
  allocate(hparoi(nbfcou))

  ! Compressible: coupling with energy

  if (ientha .eq. 2) then
    iepsel = 1
    iepsfa = iepsel + ncelet
    igamag = iepsfa + nfabor
    ixmasm = igamag + ncelet
    ifinwa = ixmasm + ncelet
    allocate(wa(ifinwa))
  endif

  ! Loop on coupled faces to compute coefficients

  inbcoo = inbcou
  call lfasyr(inbcoo, lfcou(1))
  !==========

  do iloc = 1, nbfcou

    ifac = lfcou(iloc)

    ! Saved fluid temperatures and exchange coefficients
    tfluid(iloc) = tbord(ifac)
    hparoi(iloc) = hbord(ifac)

  enddo

  ! In enthalpy formulation, transform to temperatures for SYRTHES
  !  To conserve flux Phi = (lambda/d     ) Delta T
  !                oe Phi = (lambda/(d Cp)) Delta H
  !  we multiply hbord = lambda/(d Cp) by Cp in the adjacent cell.
  !  Conservation is not guaranteed, so we add a warning.

  if (ientha.eq.1) then

    write(nfecra,1000)
    mode = 1
    do iloc = 1, nbfcou
      ifac = lfcou(iloc)
      iel  = ifabor(ifac)
      enthal = tfluid(iloc)
      call usthht (mode   , enthal , temper  )
      !==========
      tfluid(iloc) = temper
      if (ncp.eq.ncelet) then
        hparoi(iloc) = hparoi(iloc)*cp(iel)
      else
        hparoi(iloc) = hparoi(iloc)*cpcst
      endif
    enddo

  else if (ientha.eq.2) then

    ! In energy formulation, transform to temperatures for SYRTHES
    !  To conserve flux Phi = (lambda/d     ) Delta T
    !                oe Phi = (lambda/(d Cp)) Delta H
    !  we multiply hbord = lambda/(d Cp) by Cv in the adjacent cell.
    !  Note that Ei = Cv Ti + 1/2 Ui*Ui + Epsilon_sup_i
    !  and  that Ep = Cv Tp + 1/2 Ui*Ui + Epsilon_sup_i
    !    (the difference is thus Cv Delta T)

    ! Modify temperature and exchange coefficient

    ! Compute e - CvT

    iccfth = 7
    imodif = 0

    call uscfth                                                   &
    !==========
 ( nvar   , nscal  ,                                              &
   iccfth , imodif ,                                              &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   wa(iepsel) , wa(iepsfa) , wa(igamag) , wa(ixmasm) )
   !---------   ---------

    do iloc = 1, nbfcou
      ifac  = lfcou(iloc)
      iel   = ifabor(ifac)
      energ = tfluid(iloc)
      cvt   = energ                                               &
             -(0.5d0*( rtp(iel,iu)**2                      &
                      +rtp(iel,iv)**2                      &
                      +rtp(iel,iw)**2)                     &
               + wa(iepsel+iel-1)           )
       if (ncv.eq.ncelet) then
         tfluid(iloc) = cvt/cv(iel)
         hparoi(iloc) = hparoi(iloc)*cv(iel)
       else
         tfluid(iloc) = cvt/cvcst
         hparoi(iloc) = hparoi(iloc)*cvcst
       endif
     enddo

  endif

  ! Send fluid temperature and exchange coefficient

  inbcoo = inbcou
  call varsyo (inbcoo, tfluid(1), hparoi(1))
  !==========

  ! Free memory
  if (ientha .eq. 2) deallocate(wa)
  deallocate(hparoi)
  deallocate(tfluid)
  deallocate(lfcou)

enddo

!===============================================================================
! End of boundary couplings
!===============================================================================

return

! Formats

#if defined(_CS_LANG_FR)

 1000 format(                                                     &
'@                                                            ',/,&
'@ @@ ATTENTION : COUPLAGE SYRTHES AVEC CALCUL EN ENTHALPIE   ',/,&
'@    =========                                               ',/,&
'@      OPTION NON VALIDEE - CONTACTER L''EQUIPE DE DVPT      ',/,&
'@                                                            ')

#else

 1000 format(                                                     &
'@                                                            ',/,&
'@ @@ WARNING: SYRTHES COUPLING WITH ENTHALPY CALCULATION     ',/,&
'@    ========                                                ',/,&
'@      OPTION NOT VALIDATED - CONTACT THE SUPPORT            ',/,&
'@                                                            ')

#endif

end subroutine
