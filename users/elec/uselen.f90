!-------------------------------------------------------------------------------

!VERS


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

subroutine uselen &
!================

 ( nummai ,                                                       &
   nvar   , nscal  ,                                              &
   ncelps , nfacps , nfbrps ,                                     &
   lstcel , lstfac , lstfbr ,                                     &
   ia     ,                                                       &
   dt     , rtpa   , rtp    , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   tracel , trafac , trafbr ,                                     &
   ra     )

!===============================================================================
! Purpose :
! --------

! For post-processing in electric module

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
!    nom           !type!mode !                   role                         !
!__________________!____!_____!________________________________________________!
! nummai           ! ec ! <-- ! numero du maillage post                        !
! nvar             ! e  ! <-- ! nombre total de variables                      !
! nscal            ! e  ! <-- ! nombre total de scalaires                      !
! ncelps           ! e  ! <-- ! nombre de cellules du maillage post            !
! nfacps           ! e  ! <-- ! nombre de faces interieur post                 !
! nfbrps           ! e  ! <-- ! nombre de faces de bord post                   !
! lstcel(ncelps    ! te ! <-- ! liste des cellules du maillage post            !
! lstfac(nfacps    ! te ! <-- ! liste des faces interieures post               !
! lstfbr(nfbrps    ! te ! <-- ! liste des faces de bord post                   !
! ia(*)            ! tr ! --- ! macro tableau entier                           !
! dt(ncelet)       ! tr ! <-- ! pas de temps                                   !
! rtp, rtpa        ! tr ! <-- ! variables de calcul au centre des              !
! (ncelet,*)       !    !     !    cellules (instant courant ou prec)          !
! propce           ! tr ! <-- ! proprietes physiques au centre des             !
! (ncelet,*)       !    !     !    cellules                                    !
! propfa           ! tr ! <-- ! proprietes physiques au centre des             !
!  (nfac,*)        !    !     !    faces internes                              !
! propfb           ! tr ! <-- ! proprietes physiques au centre des             !
!  (nfabor,*)      !    !     !    faces de bord                               !
! coefa, coefb     ! tr ! <-- ! conditions aux limites aux                     !
!  (nfabor,*)      !    !     !    faces de bord                               !
! tracel(*)        ! tr ! <-- ! tab reel valeurs cellules post                 !
! trafac(*)        ! tr ! <-- ! tab reel valeurs faces int. post               !
! trafbr(*)        ! tr ! <-- ! tab reel valeurs faces bord post               !
! ra(*)            ! tr ! --- ! macro tableau reel                             !
!__________________!____!_____!________________________________________________!

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
use dimens, only: ndimfb
use pointe
use numvar
use optcal
use cstphy
use cstnum
use entsor
use parall
use period
use ppppar
use ppthch
use ppincl
use elincl
use mesh

!===============================================================================

implicit none

! Arguments

integer          nummai
integer          nvar   , nscal
integer          ncelps , nfacps , nfbrps
integer          idimt

integer          lstcel(ncelps), lstfac(nfacps), lstfbr(nfbrps)
integer          ia(*)

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(ndimfb,*)
double precision coefa(ndimfb,*), coefb(ndimfb,*)
double precision tracel(ncelps*3)
double precision trafac(nfacps*3), trafbr(nfbrps*3)
double precision ra(*)

! Local variables

character*32     namevr
integer          iel   , iloc
integer          ivar  , ivar0 , inc   , iccocg
integer          nswrgp, imligp, iwarnp, iclimv
integer          ipcsii
integer          ientla, ivarpr
double precision epsrgp, climgp, extrap
double precision rbid(1)

double precision, allocatable, dimension(:,:) :: grad

!===============================================================================
!===============================================================================
! 0.  PAR DEFAUT, ON CONSIDERE QUE LE SOUS PROGRAMME CI-DESSOUS CONVIENT
!       A L'UTILISATEUR, C'EST-A-DIRE QUE LA MISE EN OEUVRE DU MODULE
!       ELECTRIQUE DECLENCHE LA PRODUCTION DE CHAMPS STANDARD DANS LE
!       POST-TRAITEMENT.
!     L'UTILISATEUR N'A PAS A MODIFIER LE PRESENT SOUS-PROGRAMME DANS
!       LES CONDITIONS D'UTILISATION STANDARD.
!     DANS LE CAS OU IL SOUHAITE PRODUIRE DES VARIABLES SUPPLEMENTAIRES
!       IL PEUT LES AJOUTER A LA FIN, VOIR LA DOCUMENTATION DE USEEVO
!===============================================================================

if(nummai.eq.-1) then

  ! Allocate work arrays
  allocate(grad(ncelet,3))

!===============================================================================
! 1.   Graident of the real potential
!===============================================================================

  idimt  = 3
  NAMEVR = 'Gr_PotR'

  ivar = isca(ipotr)
  iclimv = iclrtp(ivar,icoef)

  inc = 1
  iccocg = 1
  nswrgp = nswrgr(ivar)
  imligp = imligr(ivar)
  iwarnp = iwarni(ivar)
  epsrgp = epsrgr(ivar)
  climgp = climgr(ivar)
  extrap = extrag(ivar)
  ivar0 = 0
!
  call grdcel                                                     &
  !==========
 ( ivar0  , imrgra , inc    , iccocg , nswrgp , imligp ,          &
   iwarnp , nfecra , epsrgp , climgp , extrap ,                   &
   ia     ,                                                       &
   rtp(1,ivar), coefa(1,iclimv) , coefb(1,iclimv)  ,              &
!       POTR
   grad(1,1) , grad(1,2) , grad(1,3) ,                            &
!       d POTR /dx   d POTR /dy   d POTR /dz
   ra     )

!
  ientla = 0
  ivarpr = 1

  call psteva(nummai, namevr, idimt, ientla, ivarpr,              &
  !==========
              ntcabs, ttcabs, grad, rbid, rbid)

!===============================================================================
! 2.   For Joule Heating by direct conduction :
!                           gradient of the imaginary component of the potential
!===============================================================================

  if (ippmod(ieljou).eq.2 .or. ippmod(ieljou).eq.4) then

    idimt  = 3
    NAMEVR = 'Gr_PotI'

    ivar = isca(ipoti)
    iclimv = iclrtp(ivar,icoef)

    inc = 1
    iccocg = 1
    nswrgp = nswrgr(ivar)
    imligp = imligr(ivar)
    iwarnp = iwarni(ivar)
    epsrgp = epsrgr(ivar)
    climgp = climgr(ivar)
    extrap = extrag(ivar)
!
    ivar0 = 0
!
    call grdcel                                                   &
    !==========
 ( ivar0  , imrgra , inc    , iccocg , nswrgp , imligp ,          &
   iwarnp , nfecra , epsrgp , climgp , extrap ,                   &
   ia     ,                                                       &
   rtp(1,ivar), coefa(1,iclimv) , coefb(1,iclimv)  ,              &
!       POTI
   grad(1,1) , grad(1,2) , grad(1,3) ,                            &
!       d POTI /dx   d POTI /dy   d POTI /dz
   ra     )

!
    ientla = 0
    ivarpr = 1

    call psteva(nummai, namevr, idimt, ientla, ivarpr,            &
    !==========
                ntcabs, ttcabs, grad, rbid, rbid)

  endif

!===============================================================================
! 3.  For Joule heating by direct conduction :
!                                     imaginary component of the current density
!===============================================================================

  if(ippmod(ieljou).eq.2 .or. ippmod(ieljou).eq.4 ) then

    idimt  = 3
    NAMEVR = 'Cour_Im'

    ivar = isca(ipoti)
    iclimv = iclrtp(ivar,icoef)

!    As in elflux
    ipcsii = ipproc(ivisls(ipoti))

    inc = 1
    iccocg = 1
    nswrgp = nswrgr(ivar)
    imligp = imligr(ivar)
    iwarnp = iwarni(ivar)
    epsrgp = epsrgr(ivar)
    climgp = climgr(ivar)
    extrap = extrag(ivar)
!
    ivar0 = 0

    call grdcel                                                   &
    !==========
 ( ivar0  , imrgra , inc    , iccocg , nswrgp , imligp ,          &
   iwarnp , nfecra , epsrgp , climgp , extrap ,                   &
   ia     ,                                                       &
   rtp(1,ivar), coefa(1,iclimv) , coefb(1,iclimv)  ,              &
!       POTI
   grad(1,1) , grad(1,2) , grad(1,3) ,                            &
!       d POTI /dx   d POTI /dy   d POTI /dz
   ra     )

    do iloc = 1, ncelps
      iel = lstcel(iloc)
      tracel(iloc)          = -propce(iel,ipcsii)*grad(iel,1)
      tracel(iloc+ncelps)   = -propce(iel,ipcsii)*grad(iel,2)
      tracel(iloc+2*ncelps) = -propce(iel,ipcsii)*grad(iel,3)
    enddo
!
    ientla = 0
    ivarpr = 0

    call psteva(nummai, namevr, idimt, ientla, ivarpr,            &
    !==========
                ntcabs, ttcabs, tracel, rbid, rbid)

  endif

!==========================================================
! 5.   For electric arc : electromagnetic field calculation
!==========================================================

  if( ippmod(ielarc).ge.2 ) then

    idimt  = 3
    NAMEVR = 'Ch_Mag'

!   Ax Component

    ivar = isca(ipotva(1))
    iclimv = iclrtp(ivar,icoef)

    inc = 1
    iccocg = 1
    nswrgp = nswrgr(ivar)
    imligp = imligr(ivar)
    iwarnp = iwarni(ivar)
    epsrgp = epsrgr(ivar)
    climgp = climgr(ivar)
    extrap = extrag(ivar)
!
    ivar0 = 0
!
    call grdcel                                                   &
    !==========
 ( ivar0  , imrgra , inc    , iccocg , nswrgp , imligp ,          &
   iwarnp , nfecra , epsrgp , climgp , extrap ,                   &
   ia     ,                                                       &
   rtp(1,ivar), coefa(1,iclimv) , coefb(1,iclimv)  ,              &
   grad(1,1) , grad(1,2) , grad(1,3) ,                            &
!       d Ax /dx   d Ax /dy   d Ax /dz
   ra     )

!       B = rot A ( B = curl A)

    do iloc = 1, ncelps
      iel = lstcel(iloc)
      tracel(iloc)          =  zero
      tracel(iloc+ncelps)   =  grad(iel,3)
      tracel(iloc+2*ncelps) = -grad(iel,2)
    enddo

!    Ay component

    ivar = isca(ipotva(2))
    iclimv = iclrtp(ivar,icoef)

    inc = 1
    iccocg = 1
    nswrgp = nswrgr(ivar)
    imligp = imligr(ivar)
    iwarnp = iwarni(ivar)
    epsrgp = epsrgr(ivar)
    climgp = climgr(ivar)
    extrap = extrag(ivar)
!
    ivar0 = 0
!
    call grdcel                                                   &
    !==========
  ( ivar0  , imrgra , inc    , iccocg , nswrgp , imligp ,         &
    iwarnp , nfecra , epsrgp , climgp , extrap ,                  &
    ia     ,                                                      &
    rtp(1,ivar), coefa(1,iclimv) , coefb(1,iclimv) ,              &
    grad(1,1) , grad(1,2) , grad(1,3) ,                           &
!       d Ay /dx   d Ay /dy   d Ay /dz
    ra     )

!       B = rot A (B = curl A)

    do iloc = 1, ncelps
      iel = lstcel(iloc)
      tracel(iloc)          = tracel(iloc)          - grad(iel,3)
      tracel(iloc+ncelps)   = tracel(iloc + ncelps) + zero
      tracel(iloc+2*ncelps) = tracel(iloc+2*ncelps) + grad(iel,1)
    enddo

!    Az component

    ivar = isca(ipotva(3))
    iclimv = iclrtp(ivar,icoef)

    inc = 1
    iccocg = 1
    nswrgp = nswrgr(ivar)
    imligp = imligr(ivar)
    iwarnp = iwarni(ivar)
    epsrgp = epsrgr(ivar)
    climgp = climgr(ivar)
    extrap = extrag(ivar)
!
    ivar0 = 0
!
    call grdcel                                                   &
    !==========
  ( ivar0  , imrgra , inc    , iccocg , nswrgp , imligp ,         &
    iwarnp , nfecra , epsrgp , climgp , extrap ,                  &
    ia     ,                                                      &
    rtp(1,ivar), coefa(1,iclimv) , coefb(1,iclimv) ,              &
    grad(1,1) , grad(1,2) , grad(1,3) ,                           &
!       d Az /dx   d Az /dy   d Az /dz
    ra     )

!       B = rot A (B = curl A)

    do iloc = 1, ncelps
      iel = lstcel(iloc)
      tracel(iloc)          = tracel(iloc)          + grad(iel,2)
      tracel(iloc+ncelps)   = tracel(iloc+ncelps)   - grad(iel,1)
      tracel(iloc+2*ncelps) = tracel(iloc+2*ncelps) + zero
    enddo
!
    ientla = 0
    ivarpr = 0

    call psteva(nummai, namevr, idimt, ientla, ivarpr,            &
    !==========
                ntcabs, ttcabs, tracel, rbid, rbid)

  endif


!===============================================================================
! 4.   Calculation of Module and Argument of the complex potential if IELJOU = 4
!===============================================================================

  if (ippmod(ieljou).eq.4) then

    idimt  = 1
    NAMEVR = 'ModPot'

    ivar = isca(ipotr)

    do iloc = 1, ncelps
      iel = lstcel(iloc)
      tracel(iloc) =                                              &
        sqrt( rtp(iel,isca(ipotr))*rtp(iel,isca(ipotr))           &
             +rtp(iel,isca(ipoti))*rtp(iel,isca(ipoti)) )
    enddo

    ientla = 0
    ivarpr = 0

    call psteva(nummai, namevr, idimt, ientla, ivarpr,            &
    !==========
                ntcabs, ttcabs, tracel, rbid, rbid)

    idimt  = 1
    NAMEVR = 'ArgPot'

    ivar = isca(ipotr)

    do iloc = 1, ncelps
      iel = lstcel(iloc)
      if ( rtp(iel,isca(ipotr)) .ne. 0.d0 ) then
        if ( rtp(iel,isca(ipotr)) .ge. 0.d0 ) then
          tracel(iloc) =                                          &
           atan( rtp(iel,isca(ipoti))/rtp(iel,isca(ipotr)))
        else
          if ( rtp(iel,isca(ipoti)) .gt. 0.d0 ) then
            tracel(iloc) =                                        &
              4.d0*atan(1.d0)                                     &
             +atan( rtp(iel,isca(ipoti))                          &
                   /rtp(iel,isca(ipotr)))
          else
            tracel(iloc) =                                        &
             -4.d0*atan(1.d0)                                     &
             +atan( rtp(iel,isca(ipoti))                          &
                   /rtp(iel,isca(ipotr)))
          endif
        endif
      else
        tracel(iloc) = 2.d0*atan(1.d0)
      endif

      if (tracel(iloc) .lt. 0.d0) then
        tracel(iloc) = tracel(iloc) + 8.d0**atan(1.d0)
      endif

    enddo

    ientla = 0
    ivarpr = 0

    call psteva(nummai, namevr, idimt, ientla, ivarpr,            &
    !==========
                ntcabs, ttcabs, tracel, rbid, rbid)

  endif

  ! Free memory
  deallocate(grad)

endif

return

end subroutine
