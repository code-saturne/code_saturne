!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2016 EDF S.A.
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

subroutine cs_coal_fp2st &
!=======================

 ( iscal  ,                                                        &
   propce ,                                                        &
   smbrs  , rovsdt )

!===============================================================================
! FONCTION :
! ----------

! ROUTINE PHYSIQUE PARTICULIERE : FLAMME CHARBON PULVERISE
!   TERMES SOURCES DE PRODUCTION ET DE DISSIPATION POUR
!   LA VARIANCE (BILANS EXPLICITE ET IMPLICITE)

!-------------------------------------------------------------------------------
!ARGU                             ARGUMENTS
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! smbrs(ncelet)    ! tr ! --> ! second membre explicite                        !
! rovsdt(ncelet    ! tr ! --> ! partie diagonale implicite                     !
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
use numvar
use entsor
use optcal
use cstphy
use cstnum
use parall
use period
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use ppcpfu
use cs_coal_incl
use mesh
use field

!===============================================================================

implicit none

! Arguments

integer          iscal

double precision propce(ncelet,*)
double precision smbrs(ncelet), rovsdt(ncelet)

! Local variables

integer           iel    , ifac   , ivar0
integer           icla   , icha   , numcha
integer           inc    , iccocg , nswrgp , imligp , iwarnp
integer           ipcx2c
integer           ipcgd1 , ipcgd2
integer           iold

double precision xk     , xe     , rhovst
double precision epsrgp , climgp , extrap
double precision aux
double precision gdev1 , gdev2
double precision fsd   , fdev  , diamdv , gdev

integer           iok1,iok2
double precision, dimension(:) ,allocatable :: x1,f1f2
double precision, dimension(:) ,allocatable :: coefap , coefbp
double precision, allocatable, dimension(:,:) :: grad
double precision, dimension(:), pointer ::  crom
double precision, dimension(:), pointer :: visct
double precision, dimension(:), pointer :: cvara_k, cvara_ep, cvara_omg
double precision, dimension(:), pointer :: cvara_r11, cvara_r22, cvara_r33
double precision, dimension(:), pointer :: cvar_xchcl, cvar_xckcl, cvar_xnpcl
double precision, dimension(:), pointer :: cvar_xwtcl
double precision, dimension(:), pointer :: cvar_f1m, cvar_f2m
double precision, dimension(:), pointer :: cvara_scal

!===============================================================================
! 1. Initialization
!===============================================================================

!===============================================================================
! Deallocation dynamic arrays
!----
allocate(x1(1:ncelet) , f1f2(1:ncelet),                STAT=iok1)
allocate(grad(ncelet,3), STAT=iok1)
if ( iok1 > 0 ) then
  write(nfecra,*) ' Memory allocation error inside: '
  write(nfecra,*) '     cs_coal_fp2st               '
  call csexit(1)
endif
!===============================================================================

! --- La variance n'est pas associé a un scalaire mais a f1+f2
call field_get_val_prev_s(ivarfl(isca(iscal)), cvara_scal)

! --- Numero des grandeurs physiques
call field_get_val_s(icrom, crom)
call field_get_val_s(iprpfl(ivisct), visct)

if ( itytur.eq.2 .or. iturb.eq.50 ) then
  call field_get_val_prev_s(ivarfl(ik), cvara_k)
  call field_get_val_prev_s(ivarfl(iep), cvara_ep)
elseif ( itytur.eq.3 ) then
  call field_get_val_prev_s(ivarfl(ir11), cvara_r11)
  call field_get_val_prev_s(ivarfl(ir22), cvara_r22)
  call field_get_val_prev_s(ivarfl(ir33), cvara_r33)
  call field_get_val_prev_s(ivarfl(iep), cvara_ep)
elseif ( iturb.eq.60 ) then
  call field_get_val_prev_s(ivarfl(ik), cvara_k)
  call field_get_val_prev_s(ivarfl(iomg), cvara_omg)
endif

!===============================================================================
! 2. PRISE EN COMPTE DES TERMES SOURCES DE PRODUCTION PAR LES GRADIENTS
!    ET DE DISSIPATION
!===============================================================================
if ( itytur.eq.2 .or. iturb.eq.50 .or.             &
     itytur.eq.3 .or. iturb.eq.60      ) then
  inc = 1
  iccocg = 1
! A defaut de savoir pour F1M+F2M on prend comme pour F1M(1)
  nswrgp = nswrgr(isca(if1m(1)))
  imligp = imligr(isca(if1m(1)))
  iwarnp = iwarni(isca(if1m(1)))
  epsrgp = epsrgr(isca(if1m(1)))
  climgp = climgr(isca(if1m(1)))
  extrap = extrag(isca(if1m(1)))

! --> calcul de X1

  x1( : ) = 1.d0
  do icla = 1, nclacp
    call field_get_val_s(ivarfl(isca(ixch(icla))), cvar_xchcl)
    call field_get_val_s(ivarfl(isca(ixck(icla))), cvar_xckcl)
    call field_get_val_s(ivarfl(isca(inp(icla))), cvar_xnpcl)
    if ( ippmod(iccoal) .ge. 1 ) then
      call field_get_val_s(ivarfl(isca(ixwt(icla))), cvar_xwtcl)
    endif
    do iel = 1, ncel
      x1(iel) =   x1(iel)                                        &
               -( cvar_xchcl(iel)                                &
                 +cvar_xckcl(iel)                                &
                 +cvar_xnpcl(iel)*xmash(icla) )
      if ( ippmod(iccoal) .ge. 1 ) then
        x1(iel) = x1(iel) - cvar_xwtcl(iel)
      endif
    enddo
  enddo

! --> calcul de F=F1+F2
  f1f2( : ) = zero
  do icha = 1, ncharb
    call field_get_val_s(ivarfl(isca(if1m(icha))), cvar_f1m)
    call field_get_val_s(ivarfl(isca(if2m(icha))), cvar_f2m)
    do iel = 1, ncel
      f1f2(iel) =  f1f2(iel)                                     &
                 + cvar_f1m(iel)                                 &
                 + cvar_f2m(iel)
    enddo
  enddo
  do iel = 1, ncel
    f1f2(iel) = f1f2(iel)/x1(iel)
  enddo

! --> Calcul du gradient de f1f2

  ! Allocate temporary arrays
  allocate(coefap(nfabor), coefbp(nfabor), STAT=iok1)

  if (iok1 > 0) then
    write(nfecra,*) ' Memory allocation error inside : '
    write(nfecra,*) '     cs_coal_fp2st                '
    call csexit(1)
  endif

  do ifac = 1, nfabor
    ! Homogenous Neumann on the gradient
    coefap(ifac) = zero
    coefbp(ifac) = 1.d0
  enddo

! En periodique et parallele, echange avant calcul du gradient

  if (irangp.ge.0.or.iperio.eq.1) then
    call synsca(f1f2)
    !==========
  endif

!  IVAR0 = 0 (indique pour la periodicite de rotation que la variable
!     n'est pas la vitesse ni Rij)
  ivar0  = 0
  call grdcel &
  !==========
 ( ivar0  , imrgra , inc    , iccocg , nswrgp , imligp ,          &
   iwarnp , nfecra , epsrgp , climgp , extrap ,                   &
   f1f2   , coefap , coefbp ,                                     &
   grad   )

  do iel = 1, ncel
    if ( itytur.eq.2 .or. iturb.eq.50 ) then
      xk = cvara_k(iel)
      xe = cvara_ep(iel)
    elseif ( itytur.eq.3 ) then
      xk = 0.5d0*(cvara_r11(iel)+cvara_r22(iel)+cvara_r33(iel))
      xe = cvara_ep(iel)
    elseif ( iturb.eq.60 ) then
      xk = cvara_k(iel)
      xe = cmu*xk*cvara_omg(iel)
    endif

    rhovst = propce(iel,ipproc(irom1))*xe/(xk*rvarfl(iscal))*volume(iel)
    rovsdt(iel) = rovsdt(iel) + max(zero,rhovst)
    smbrs(iel) = smbrs(iel)                                            &
                + 2.d0*visct(iel)*volume(iel)/sigmas(iscal)            &
                 *( grad(iel,1)**2.d0 + grad(iel,2)**2.d0              &
                  + grad(iel,3)**2.d0 )*x1(iel) - rhovst*cvara_scal(iel)

  enddo

endif

!===============================================================================
! 3. PRISE EN COMPTE DES TERMES SOURCES RELATIF AUX ECHANGES INTERFACIAUX
!==============================================================================
!
! 2 versions disponible
!   iold = 1 ===> c'est l'ancienne version
!   iold = 2 ===> c'est la nouvelle
!
iold = 1
!
if ( iold .eq. 1 ) then
!
  do icla=1,nclacp
    numcha = ichcor(icla)
    ipcx2c = ipproc(ix2(icla))
    call field_get_val_s(ivarfl(isca(ixch(icla))), cvar_xchcl)
    call field_get_val_s(ivarfl(isca(inp(icla))), cvar_xnpcl)
    ipcgd1 = ipproc(igmdv1(icla))
    ipcgd2 = ipproc(igmdv2(icla))
    do iel = 1, ncel
      gdev1 = -crom(iel)*propce(iel,ipcgd1)               &
                                 *cvar_xchcl(iel)
      gdev2 = -crom(iel)*propce(iel,ipcgd2)               &
                                 *cvar_xchcl(iel)
      gdev  = gdev1 + gdev2
!
      if ( cvar_xnpcl(iel) .gt. epsicp ) then
        diamdv = diam20(icla)
        fsd  =  1.d0 - (1.d0-f1f2(iel))                            &
               * exp( ( cvar_xchcl(iel)                            &
                       *(propce(iel,ipcgd1)+propce(iel,ipcgd2))  ) &
                     /( 2.d0*pi*2.77d-4*diamdv                     &
                        *cvar_xnpcl(iel)*crom(iel) ) )
        fdev = 1.d0
!
! ts explicite
!
        if ( (fsd-f1f2(iel))*(2.d0*fdev-fsd-f1f2(iel)) .gt. epsicp ) then
          smbrs(iel) = smbrs(iel)                                          &
                     + volume(iel)*(gdev1+gdev2)                           &
                      *(fsd-f1f2(iel))*(2.d0*fdev-fsd-f1f2(iel))
        endif
      endif
!
    enddo
  enddo
!
else
!
  do icla=1,nclacp
    numcha = ichcor(icla)
    ipcx2c = ipproc(ix2(icla))
    call field_get_val_s(ivarfl(isca(ixch(icla))), cvar_xchcl)
    do iel = 1, ncel
      gdev1 = -crom(iel)*propce(iel,ipproc(igmdv1(icla)))       &
                                 *cvar_xchcl(iel)
      gdev2 = -crom(iel)*propce(iel,ipproc(igmdv2(icla)))       &
                                 *cvar_xchcl(iel)
!
      aux  = (gdev1+gdev2)               *(1.d0-f1f2(iel))**2.d0
!
! ts implicite : pour l'instant on implicite de facon simple
!
      if ( abs(f1f2(iel)*(1.d0-f1f2(iel))) .GT. epsicp ) then
        rhovst = aux*cvara_scal(iel)/((f1f2(iel)*(1-f1f2(iel)))**2.d0)     &
                    *volume(iel)
      else
        rhovst = 0.d0
      endif
      rovsdt(iel) = rovsdt(iel) + max(zero,rhovst)
! ts explicite
      smbrs(iel) = smbrs(iel)+aux*volume(iel)-rhovst*cvara_scal(iel)

    enddo
  enddo
!
endif
!

!--------
! Formats
!--------

! Free memory
deallocate(x1,f1f2,grad,STAT=iok1)
deallocate(coefap, coefbp, STAT=iok2)

if ( iok1 > 0 .or. iok2 > 0) then
  write(nfecra,*) ' Memory deallocation error inside: '
  write(nfecra,*) '     cs_coal_fp2st                 '
  call csexit(1)
endif
!===============================================================================

!----
! End
!----

return
end subroutine
