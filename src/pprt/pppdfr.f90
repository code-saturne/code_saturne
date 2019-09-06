!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2019 EDF S.A.
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

!===============================================================================
! Function:
! ---------

!> \file pppdfr.f90
!>
!> \brief Specific physic subroutine: Calculation of rectangle-Dirac pdf parameters
!
! from P. Plion & A. Escaich
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     ncel          number of cells
!> \param[in]     ncelet        number of extended (real + ghost) cells
!> \param[in]     indpdf        indicator for pdf integration or mean value
!> \param[out]    tpdf          indicator for pdf shape:
!                               - 0: Dirac at mean value
!                               - 1: rectangle
!                               - 2: Dirac's peak at \f$ f_{min} \f$
!                               - 3: Dirac's peak at \f$ f_{max} \f$
!                               - 4: rectangle and 2 Dirac's pics
!> \param[in]     fm            mean mixture fraction at cell centers
!> \param[in]     fp2m          mean mixture fraction variance at cell centers
!> \param[in]     fmini         mixture fraction low boundary
!> \param[in]     fmaxi         mixture fraction high boundary
!> \param[in]     dirmin        Dirac's peak value at \f$ f_{min} \f$
!> \param[in]     dirmax        Dirac's peak value at \f$ f_{max} \f$
!> \param[in]     fdeb          abscissa of rectangle low boundary
!> \param[in]     ffin          abscissa of rectangle high boundary
!> \param[in]     hrec          rectangle height
!_______________________________________________________________________________

subroutine pppdfr &
 ( ncelet , ncel   , indpdf ,                                     &
   tpdf   ,                                                       &
   fm     , fp2m   ,                                              &
   fmini  , fmaxi  ,                                              &
   dirmin , dirmax , fdeb   , ffin   , hrec )

!===============================================================================

!==============================================================================
! Module files
!==============================================================================

use paramx
use numvar
use optcal
use cstphy
use cstnum
use entsor
use pointe
use parall
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer          ncelet, ncel
integer          indpdf(ncelet)

double precision tpdf(ncelet)
double precision fm(ncelet), fp2m(ncelet)
double precision fmini(ncelet), fmaxi(ncelet)
double precision dirmin(ncelet), dirmax(ncelet)
double precision fdeb(ncelet), ffin(ncelet)
double precision hrec(ncelet)


! Local variables

integer          iel, n1, n2, n3, n4, n5 , n6 ,nfp2 , nbspdf
double precision t1, t2, t3, t1mod, t2mod , fp2max
double precision fp2mmax1,fp2mmin1,fp2mmax2,fp2mmin2

logical(kind=c_bool) :: log_active

!===============================================================================
! 0.  INITIALISATION
!===============================================================================

log_active = cs_log_default_is_active()

do iel = 1, ncel

  indpdf(iel) = 0

  tpdf  (iel) = 0.d0
  dirmin(iel) = 0.d0
  dirmax(iel) = 0.d0
  fdeb  (iel) = 0.d0
  ffin  (iel) = 0.d0
  hrec  (iel) = 0.d0

enddo

!===============================================================================
! 1.  CALCULS PRELIMINAIRES
!===============================================================================

! Parametre relatif a la variance
t1 = 1.d-08
! Parametre relatif a la moyenne
t2 = 5.d-07

do iel = 1, ncel

! Modifs des parametres T1 et T2 afin de tenir compte du fait que
!   FMINI < FM < FMAXI
  t1mod = t1*(fmaxi(iel)-fmini(iel))**2
  t2mod = t2*(fmaxi(iel)-fmini(iel))
  if ( (fp2m(iel).gt.t1mod)                                       &
       .and.(fm(iel) .ge. (fmini(iel) + t2mod))                   &
       .and.(fm(iel) .le. (fmaxi(iel) - t2mod)) ) then
    indpdf(iel) = 1
  endif
enddo

! Clipping de la variance

fp2mmin1 =  1.D+20
fp2mmax1 = -1.D+20
do iel = 1, ncel
  fp2mmin1 = min(fp2mmin1,fp2m(iel))
  fp2mmax1 = max(fp2mmax1,fp2m(iel))
enddo
if ( irangp .ge.0 ) then
  call parmin(fp2mmin1)
  call parmax(fp2mmax1)
endif

nfp2 = 0
do iel = 1, ncel
  fp2max = (fmaxi(iel)-fm(iel))*(fm(iel)-fmini(iel))
  if ( fp2m(iel) .gt. fp2max+1.d-20 ) then
    fp2m(iel) = fp2max
    nfp2 = nfp2 + 1
  endif
enddo
if (irangp .ge. 0) then
  call parcpt(nfp2)
endif
if (log_active) then

  write(nfecra,*) ' pppdfr: variance clipping points: ', nfp2

  fp2mmin2 = 1.d+20
  fp2mmax2 =-1.d+20
  do iel = 1, ncel
    fp2mmin2 = min(fp2mmin2,fp2m(iel))
    fp2mmax2 = max(fp2mmax2,fp2m(iel))
  enddo
  if (irangp .ge.0) then
    call parmin(fp2mmin2)
    call parmax(fp2mmax2)
  endif

  if (nfp2 .gt. 0) then
    write(nfecra, *) '     Variance before clipping min and max: ', &
                     fp2mmin1, fp2mmax1
    write(nfecra, *) '     Variance after  clipping min and max: ', &
                     fp2mmin2, fp2mmax2
  endif
endif

!===============================================================================
! 2.  CALCUL DES PARAMETRES DE LA FONCTION DENSITE DE PROBABILITE
!===============================================================================

do iel = 1, ncel

  if ( indpdf(iel).eq.1 ) then

    if (    (     (fm(iel) .le.(fmini(iel) + fmaxi(iel))*0.5d0)   &
            .and.(fp2m(iel).le.(fm(iel) - fmini(iel))**2/3.d0))   &
       .or. (     (fm(iel) .gt.(fmini(iel) + fmaxi(iel))*0.5d0)   &
            .and.(fp2m(iel).le.(fmaxi(iel) -fm(iel))**2/3.d0)) )  &
      then

! --> Rectangle seul

      tpdf  (iel) = 1.d0

      hrec(iel)   = sqrt(3.d0*fp2m(iel))
      dirmin(iel) = 0.d0
      dirmax(iel) = 0.d0
      fdeb(iel)   = fm(iel) - hrec(iel)
      ffin(iel)   = fm(iel) + hrec(iel)

    elseif(      (fm(iel)  .le.(fmini(iel) + fmaxi(iel))*0.5d0)   &
           .and. (fp2m(iel).le.((fm(iel) - fmini(iel))            &
              *(2.d0*fmaxi(iel) +fmini(iel)-3.d0*fm(iel))/3.d0)) )&
      then

! --> Rectangle et un Dirac en FMINI

      tpdf  (iel) = 2.d0

      fdeb(iel)   = fmini(iel)
      dirmax(iel) = 0.d0
      ffin(iel)   = fmini(iel) +1.5d0*( (fm(iel) - fmini(iel))**2 &
                                       + fp2m(iel) )              &
                                    /(fm(iel) - fmini(iel))
      dirmin(iel) = (3.d0*fp2m(iel) -(fm(iel) - fmini(iel))**2)   &
                  / (3.d0*((fm(iel) - fmini(iel))**2 +fp2m(iel)))

    elseif(      (fm(iel)  .gt.(fmini(iel) + fmaxi(iel))*0.5d0)   &
           .and. (fp2m(iel).le.((fmaxi(iel) - fm(iel))            &
               *(3.d0*fm(iel)-fmaxi(iel)-2.d0*fmini(iel))/3.d0)) )&
      then

! --> Rectangle et un Dirac en FMAXI (c'est juste ;
!                          le HI/81/02/03/A contient une erreur  p 12)

      tpdf  (iel) = 3.d0

      ffin(iel)   = fmaxi(iel)
      dirmin(iel) = 0.d0
      fdeb(iel)   = fmini(iel)                                    &
                  + ( 3.d0*( (fm(iel)-fmini(iel))**2+fp2m(iel) )  &
                     + (fmaxi(iel) - fmini(iel))**2               &
                     - 4.d0*(fm(iel) - fmini(iel))                &
                           *(fmaxi(iel) - fmini(iel)) )           &
                    / (2.d0*(fm(iel) - fmaxi(iel)))
      dirmax(iel) = ( 3.d0*fp2m(iel) -(fm(iel) - fmaxi(iel))**2 ) &
                 / ( 3.d0*((fm(iel) - fmaxi(iel))**2 +fp2m(iel)) )

    else

! --> Rectangle et deux Diracs

      tpdf  (iel) = 4.d0

      fdeb(iel)   = fmini(iel)
      ffin(iel)   = fmaxi(iel)
      dirmax(iel) = 3.d0*((fm(iel) - fmini(iel))**2 +fp2m(iel))   &
                    /(fmaxi(iel) - fmini(iel))**2                 &
                   -2.d0*(fm(iel) - fmini(iel))                   &
                    /(fmaxi(iel) - fmini(iel))
      dirmin(iel) = dirmax(iel) + 1.d0 - 2.d0*(fm(iel)-fmini(iel))&
                                          /(fmaxi(iel)-fmini(iel))

    endif

    if ( abs(ffin(iel) - fdeb(iel)).gt.epzero ) then
      hrec(iel) = ( 1.d0-dirmin(iel)-dirmax(iel) )                &
                / ( ffin(iel)-fdeb(iel) )
    else
      t3 = sqrt(3.d0*t1*(fmaxi(iel)-fmini(iel))**2)
      fdeb(iel) = min(fmaxi(iel),max(fmini(iel),fm(iel) - t3))
      ffin(iel) = min(fmaxi(iel),max(fmini(iel),fm(iel) + t3))
      if ( abs(ffin(iel) - fdeb(iel)).gt.epzero ) then
        hrec(iel) = ( 1.d0-dirmin(iel)-dirmax(iel) )              &
                   /( ffin(iel) - fdeb(iel) )
      else
        hrec(iel) = 0.d0
      endif

    endif

  else

    tpdf  (iel) = 0.d0

    dirmin(iel) = 0.d0
    dirmax(iel) = 0.d0
    fdeb(iel)   = 0.d0
    ffin(iel)   = 0.d0
    hrec(iel)   = 0.d0
  endif

enddo

! Verification : si Hrec <= 0 on passe sans les PDF

nbspdf = 0
do iel=1,ncel
  if ( hrec(iel) .le. epzero .and. indpdf(iel).eq.1 ) then

    indpdf(iel) = 0
    nbspdf      = nbspdf + 1

  endif
enddo

if (log_active) then
  if (irangp .ge. 0) then
    call parcpt(nbspdf)
  endif
  write(nfecra,*) ' pppdfr: switch off PDF ', nbspdf
endif

!===============================================================================
! 4.  IMPRESSION
!===============================================================================

if (log_active .eqv. .false.) return

n1 = 0
n2 = 0
n3 = 0
n4 = 0
n5 = 0
n6 = ncel
do iel = 1, ncel
  if (indpdf(iel).eq.1) then
    n1 = n1+1
    if (dirmin(iel).gt.epzero .and. dirmax(iel).lt.epzero) then
      n2 = n2+1
    else if (dirmin(iel).lt.epzero .and. dirmax(iel).gt.epzero) then
      n3 = n3+1
    else if (dirmin(iel).gt.epzero .and. dirmax(iel).gt.epzero) then
      n4 = n4+1
    else if (dirmin(iel).lt.epzero .and. dirmax(iel).lt.epzero) then
      n5 = n5+1
    endif
  endif
enddo

if (irangp.ge.0) then
  call parcpt (n1)
  call parcpt (n2)
  call parcpt (n3)
  call parcpt (n4)
  call parcpt (n5)
  call parcpt (n6)
endif

write(nfecra,1000) n1, n6
write(nfecra,2000) n5, n2, n3, n4

!----
! Formats
!----

 1000 format ( /,                                           &
'Rectangle PDF - Dirac peaks', /,                           &
'Mean, variance of transported tracer',/,                   &
'Number of turbulent points (using the PDFs)   = ', i6,/,   &
'Number of computation points                  = ', i6)
 2000 format(                                                        &
' Nb points with rectangle PDF without Dirac              = ', i6,/, &
' - - - - - - - - - -- - - - and Dirac in FMINI           = ', i6,/, &
' - - - - - - - - - -- - - - - - - - - -  FMAXI           = ', i6,/, &
' - - - - - - - - - - - - - - - Diracs in FMINI and FMAXI = ', i6,/)

!----
! Fin
!----

return
end subroutine
