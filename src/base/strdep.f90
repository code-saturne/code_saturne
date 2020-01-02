!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2020 EDF S.A.
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

subroutine strdep &
!================

 ( itrale , italim , itrfin ,                                     &
   nvar   ,                                                       &
   dt     ,                                                       &
   cofale , xprale )

!===============================================================================
! FONCTION :
! ----------

! DEPLACEMENT DES STRUCTURES MOBILES EN ALE EN COUPLAGE INTERNE

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! itrale           ! e  ! <-- ! numero d'iteration pour l'ale                  !
! italim           ! e  ! <-- ! numero d'iteration couplage implicite          !
! itrfin           ! e  ! <-- ! indicateur de derniere iteration de            !
!                  !    !     !                    couplage implicite          !
! nvar             ! i  ! <-- ! total number of variables                      !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! cofale           ! tr ! --> ! sauvegarde des cl de p et u                    !
!    (nfabor,8)    !    !     !                                                !
! xprale(ncelet    ! tr ! --> ! sauvegarde de la pression, si nterup           !
!                  !    !     !    est >1                                      !
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
use cstphy
use numvar
use optcal
use entsor
use pointe
use albase
use alstru
use alaste
use parall
use period
use mesh
use field
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer          itrale , italim , itrfin
integer          nvar

double precision dt(ncelet)
double precision xprale(ncelet)
double precision cofale(nfabor,11)

! Local variables

integer          istr, ii, iel, ifac, ntab
integer          iflmas, iflmab
integer          indast
integer          icvext, icvint, icved
integer          f_dim

double precision delta

double precision, allocatable, dimension(:,:) :: forast
double precision, dimension(:), pointer :: imasfl, bmasfl
double precision, dimension(:), pointer :: imasfl_pre, bmasfl_pre
double precision, dimension(:,:), pointer :: forbr
double precision, dimension(:,:), pointer :: coefau
double precision, dimension(:,:,:), pointer :: coefbu
double precision, dimension(:), pointer :: coefap, coefbp
double precision, dimension(:), pointer :: cvar_var, cvara_var
double precision, dimension(:,:), pointer :: cvar_varv, cvara_varv

type(var_cal_opt) :: vcopt

!===============================================================================

!===============================================================================
! Interfaces
!===============================================================================

interface

  function cs_ast_coupling_get_ext_cvg() result (icvext) &
    bind(C, name='cs_ast_coupling_get_ext_cvg')
    use, intrinsic :: iso_c_binding
    implicit none
    integer(c_int) :: icvext
  end function cs_ast_coupling_get_ext_cvg

  subroutine cs_ast_coupling_send_cvg(icv) &
    bind(C, name='cs_ast_coupling_send_cvg')
    use, intrinsic :: iso_c_binding
    implicit none
    integer(c_int), value :: icv
  end subroutine cs_ast_coupling_send_cvg

end interface

!===============================================================================
! 1. Initialization
!===============================================================================

call field_get_key_int(ivarfl(iu), kimasf, iflmas)
call field_get_key_int(ivarfl(iu), kbmasf, iflmab)
call field_get_val_s(iflmas, imasfl)
call field_get_val_prev_s(iflmas, imasfl_pre)
call field_get_val_s(iflmab, bmasfl)
call field_get_val_prev_s(iflmab, bmasfl_pre)

call field_get_val_v(iforbr, forbr)

call field_get_coefa_v(ivarfl(iu), coefau)
call field_get_coefb_v(ivarfl(iu), coefbu)

call field_get_coefa_s(ivarfl(ipr), coefap)
call field_get_coefb_s(ivarfl(ipr), coefbp)

call field_get_key_struct_var_cal_opt(ivarfl(iuma), vcopt)

!===============================================================================
! 2. Computue forces on the structures
!===============================================================================

do istr = 1, nbstru
  do ii = 1, ndim
    forsta(ii,istr) = forstr(ii,istr)
!-a tester          forsta(ii,istr) = forstP(ii,istr)
    forstr(ii,istr) = 0.d0
  enddo
enddo

! Allocate a temporary array
allocate(forast(3,nbfast))

indast = 0
do ifac = 1, nfabor
  istr = idfstr(ifac)
  if (istr.gt.0) then
    do ii = 1, 3
      forstr(ii,istr) = forstr(ii,istr) + forbr(ii,ifac)
    enddo
  else if (istr.lt.0) then
    indast = indast + 1
    do ii = 1, 3
      forast(ii,indast) = asddlf(ii,-istr)*forbr(ii,ifac)
    enddo
  endif
enddo

if (irangp.ge.0) then
  ntab = ndim*nbstru
  call parrsm(ntab,forstr)
endif

!     Calcul de l'effort envoye au structures internes
do istr = 1, nbstru
  do ii = 1, ndim
    forstp(ii,istr) = cfopre*forstr(ii,istr)+                     &
         (1.d0-cfopre)*forsta(ii,istr)
  enddo
enddo

!     Envoi de l'effort applique aux structures externes
if (nbaste.gt.0) then
  call astfor(ntcast, nbfast, forast)
endif

! Free memory
deallocate(forast)

!===============================================================================
! 3. Structure characteristic given by the user
!===============================================================================

if (nbstru.gt.0) then

    call uistr2 &
 ( xmstru, xcstru, xkstru,     &
   forstp,                     &
   dtref , ttcabs, ntcabs   )

  call usstr2                                                     &
 ( nbstru ,                                                       &
   idfstr ,                                                       &
   dt     ,                                                       &
   xmstru , xcstru , xkstru , xstreq , xstr   , xpstr  , forstp , &
   dtstr  )

endif

! If the fluid is initailizing, we do read structures
if (itrale.le.nalinf) then
  itrfin = -1
  return
endif

!===============================================================================
! 4. Internal structure displacement
!===============================================================================

do istr = 1, nbstru

  call newmrk                                                     &
 ( istr  , alpnmk  , betnmk          , gamnmk          ,          &
   xmstru(1,1,istr), xcstru(1,1,istr), xkstru(1,1,istr),          &
   xstreq(1,istr)  ,                                              &
   xstr(1,istr)    , xpstr(1,istr)   , xppstr(1,istr)  ,          &
   xsta(1,istr)    , xpsta(1,istr)   , xppsta(1,istr)  ,          &
   forstp(1,istr)  , forsta(1,istr)  , dtstr(istr)     )

enddo

!===============================================================================
! 5. Convergence test
!===============================================================================

icvext = 0
icvint = 0
icved  = 0

delta = 0.d0
do istr = 1, nbstru
  do ii = 1, 3
    delta = delta + (xstr(ii,istr)-xstp(ii,istr))**2
  enddo
enddo
if (nbstru.gt.0) then
  delta = sqrt(delta)/almax/nbstru
  if (delta.lt.epalim) icvint = 1
endif

if (nbaste.gt.0) icvext = cs_ast_coupling_get_ext_cvg()

if (nbstru.gt.0.and.nbaste.gt.0) then
  icved = icvext*icvint
elseif (nbstru.gt.0.and.nbaste.eq.0) then
  icved = icvint
elseif (nbaste.gt.0.and.nbstru.eq.0) then
  icved = icvext
endif

if (vcopt%iwarni.ge.2) write(nfecra,1000) italim, delta

! if converged
if (icved.eq.1) then
  if (itrfin.eq.1) then
    ! si itrfin=1 on sort
    if (vcopt%iwarni.ge.1) write(nfecra,1001) italim, delta
    itrfin = -1
  else
    ! otherwise one last iteration for SYRTHES/T1D/radiation
    ! and reset icved to 0 so code_aster also runs an iteration.
    itrfin = 1
    icved = 0
  endif
elseif (itrfin.eq.0 .and. italim.eq.nalimx-1) then
  ! this will be the lst iteration
  itrfin = 1
elseif (italim.eq.nalimx) then
  ! we have itrfin=1 and are finished
  if (nalimx.gt.1) write(nfecra,1100) italim, delta
  itrfin = -1
  ! Set icved to 1 so code_aster also stops
  icved = 1
endif

! Return the final convergence indicator to code_aster
call cs_ast_coupling_send_cvg(icved)

!===============================================================================
! 6. Re set previous values if required
!===============================================================================

! If nterup .gt. 1, values at previous time step have been modified after navstv
! We must then go back to a previous value
if (itrfin.ne.-1) then
  do ii = 1, nvar
    call field_get_dim (ivarfl(ii), f_dim)

    ! Fields of dimension one
    if (f_dim.eq.1) then
      call field_get_val_s(ivarfl(ii), cvar_var)
      call field_get_val_prev_s(ivarfl(ii), cvara_var)
      if (ii.eq.ipr .and. nterup.gt.1) then
        do iel = 1, ncelet
          cvara_var(iel) = xprale(iel)
        enddo
      endif
      do iel = 1, ncelet
        cvar_var(iel) = cvara_var(iel)
      enddo

    ! Vector fields
    else if (f_dim.eq.3) then
      call field_get_val_v(ivarfl(ii), cvar_varv)
      call field_get_val_prev_v(ivarfl(ii), cvara_varv)
      do iel = 1, ncelet
        cvar_varv(1, iel) = cvara_varv(1, iel)
        cvar_varv(2, iel) = cvara_varv(2, iel)
        cvar_varv(3, iel) = cvara_varv(3, iel)
      enddo
    else
      call csexit(1)
    endif
  enddo
  do ifac = 1, nfac
    imasfl(ifac) = imasfl_pre(ifac)
  enddo
  do ifac = 1, nfabor
    bmasfl(ifac) = bmasfl_pre(ifac)
    coefap(ifac) = cofale(ifac,1)
    coefau(1, ifac) = cofale(ifac,2)
    coefau(2, ifac) = cofale(ifac,3)
    coefau(3, ifac) = cofale(ifac,4)
    coefbp(ifac) = cofale(ifac,5)
    coefbu(1, 1, ifac) = cofale(ifac,6)
    coefbu(2, 2, ifac) = cofale(ifac,7)
    coefbu(3, 3, ifac) = cofale(ifac,8)
    coefbu(1, 2, ifac) = cofale(ifac,9)
    coefbu(2, 3, ifac) = cofale(ifac,10)
    coefbu(1, 3, ifac) = cofale(ifac,11)
    ! the coefficient B is supposed to be symmetric
    coefbu(2, 1, ifac) = cofale(ifac,9)
    coefbu(3, 2, ifac) = cofale(ifac,10)
    coefbu(3, 1, ifac) = cofale(ifac,11)
  enddo
endif

!--------
! Formats
!--------

#if defined(_CS_LANG_FR)

 1000 format (                                                          &
 '            ALE IMPLICITE : ITER=',I5,' DERIVE=',E12.5     )
 1001 format (                                                          &
 'CONVERGENCE ALE IMPLICITE : ITER=',I5,' DERIVE=',E12.5     )
 1100 format (                                                          &
'@                                                            ',/,&
'@ @@ ATTENTION : COUPLAGE IMPLICITE ALE                      ',/,&
'@    =========                                               ',/,&
'@  Nombre d''iterations maximal ',I10   ,' atteint           ',/,&
'@  Derive normee :',E12.5                                     ,/,&
'@                                                            '  )

#else

 1000 format (                                                          &
 '            IMPLICIT ALE: ITER=',I5,' DRIFT=',E12.5     )
 1001 format (                                                          &
 'CONVERGENCE IMPLICIT ALE: ITER=',I5,' DRIFT=',E12.5     )
 1100 format (                                                          &
'@'                                                            ,/,&
'@ @@ WARNING: IMPLICIT ALE'                                   ,/,&
'@    ========'                                                ,/,&
'@  Maximum number of iterations ',I10   ,' reached'           ,/,&
'@  Normed drift:',E12.5                                     ,/,&
'@'                                                              )

#endif

!----
! End
!----

end subroutine
