!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2018 EDF S.A.
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
! --------
!> \file cs_fuel_bcond.f90
!>
!> \brief Automatic boundary conditions
!>          Fuel combustion
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role
!______________________________________________________________________________!
!> \param[in]     itypfb        boundary face types
!> \param[in]     izfppp        zone number for the boundary face for
!>                                      the specific physic module
!> \param[in,out] icodcl        face boundary condition code:
!>                               - 1 Dirichlet
!>                               - 2 Radiative outlet
!>                               - 3 Neumann
!>                               - 4 sliding and
!>                                 \f$ \vect{u} \cdot \vect{n} = 0 \f$
!>                               - 5 smooth wall and
!>                                 \f$ \vect{u} \cdot \vect{n} = 0 \f$
!>                               - 6 rough wall and
!>                                 \f$ \vect{u} \cdot \vect{n} = 0 \f$
!>                               - 9 free inlet/outlet
!>                                 (input mass flux blocked to 0)
!>                               - 13 Dirichlet for the advection operator and
!>                                    Neumann for the diffusion operator
!> \param[in,out] rcodcl        boundary conditions value on edge faces
!>                               rcodcl(1) = value of the Dirichlet
!>                               rcodcl(2) = value of the extern exchange coef.
!>                                (infinit if no exchange)
!>                               rcodcl(3) = value of the flux density
!>                                (negative if gain) \f$w \cdot m^{-2}\f$ or
!>                                the rugosity high \f$m\f$ if  \c icodcl=6
!>                               for velocity  \f$(vistl+visct)\cdot\grad{u}\f$
!>                               for pressure  \f$dt \cdot \grad{p}\f$
!>                               for scalar
!>                                      \f$C_p(viscls+visct/turb_schmidt) \grad{t}\f$
!______________________________________________________________________________!

subroutine cs_fuel_bcond &
 ( itypfb , izfppp ,                                              &
   icodcl , rcodcl )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use optcal
use cstphy
use cstnum
use dimens, only : nvar
use entsor
use parall
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use ppcpfu
use cs_fuel_incl
use mesh
use field
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer          itypfb(nfabor)
integer          izfppp(nfabor)
integer          icodcl(nfabor,nvar)

double precision rcodcl(nfabor,nvar,3)

! Local variables

integer          ii, ifac, izone, mode, iel, ige, iok
integer          icla , ioxy
integer          icke
integer          nbrval

double precision qisqc, viscla, d2s3, uref2, rhomoy, dhy, xiturb
double precision t1, t2
double precision h1(nozppm) , h2(nozppm)
double precision x20t(nozppm)
double precision xmg0(nozppm,nclcpm)
double precision x2h20t(nozppm)
double precision qimpc(nozppm) , qcalc(nozppm)
double precision coefe(ngazem)
double precision xsolid(2)
double precision hlf, totfu
double precision dmas
double precision, dimension(:), pointer :: brom, b_x1
double precision, dimension(:), pointer :: viscl

!===============================================================================

!===============================================================================
! 0. Initializations
!===============================================================================
!
call field_get_val_s(ibrom, brom)
call field_get_val_s(iviscl, viscl)
call field_get_val_s_by_name("b_x_c", b_x1)

d2s3   = 2.d0 / 3.d0

!===============================================================================
! 1.  Parallel exchanges for the user data
!===============================================================================

!  In reality we can avoid this exchange by modifying uspcl and by
!  asking the user to provide the bulks which depend of the zone
!  out of the loop on the boundary faces: the bulks
!  would be available on all processors. However, it makes the subroutine
!  more complicated and mainly if the user modified it in a wrong way
!  it will not work.
!  We asume that all provided bulks are positive, it allows to use a max that
!  all processors know. If it is not the case, it is more complicated but
!  we can still find a max anyway.

if(irangp.ge.0) then
  call parimx(nozapm,iqimp )
  call parimx(nozapm,ientat)
  call parimx(nozapm,ientfl)
  call parimx(nozapm,inmoxy)
  call parrmx(nozapm,qimpat)
  call parrmx(nozapm,timpat)
  nbrval = nozppm
  call parrmx(nbrval,qimpfl)
  nbrval = nozppm
  call parrmx(nbrval,timpfl)
  nbrval = nozppm*nclcpm
  call parrmx(nbrval,distfu)
endif

!===============================================================================
! 2.  Velocity correction (in norm) to control the imposed flows
!       Loop over all inlet faces
!                     =========================
!===============================================================================
! --- Calculated flow
do izone = 1, nozppm
  qcalc(izone) = 0.d0
  h1(izone)    = 0.d0
enddo
do ifac = 1, nfabor
  izone = izfppp(ifac)
  qcalc(izone) = qcalc(izone) - brom(ifac) *             &
                ( rcodcl(ifac,iu,1)*surfbo(1,ifac) +       &
                  rcodcl(ifac,iv,1)*surfbo(2,ifac) +       &
                  rcodcl(ifac,iw,1)*surfbo(3,ifac) )
enddo

if (irangp.ge.0) then
  call parrsm(nozapm,qcalc )
endif

do izone = 1, nozapm
  if (iqimp(izone).eq.0) then
    qimpc(izone) = qcalc(izone)
  endif
enddo

if ( ntcabs .gt. 1 ) then
  !
  ! --- Velocity correction in norm: we do it only at the
  !     second iteration because the first one the mass density is not known yet


  iok = 0
  do ii = 1, nzfppp
    izone = ilzppp(ii)
    if ( iqimp(izone).eq.1 ) then
      if(abs(qcalc(izone)).lt.epzero) then
        write(nfecra,2001)izone,iqimp(izone),qcalc(izone)
        iok = iok + 1
      endif
    endif
  enddo
  if(iok.ne.0) then
    call csexit (1)
  endif
  do ifac = 1, nfabor
    izone = izfppp(ifac)
    if ( iqimp(izone).eq.1 ) then
      qimpc(izone) = qimpat(izone) + qimpfl(izone)
      qisqc = qimpc(izone) / qcalc(izone)
      rcodcl(ifac,iu,1) = rcodcl(ifac,iu,1)*qisqc
      rcodcl(ifac,iv,1) = rcodcl(ifac,iv,1)*qisqc
      rcodcl(ifac,iw,1) = rcodcl(ifac,iw,1)*qisqc
    endif
  enddo

else

  do izone = 1, nozapm
    if (iqimp(izone).eq.1) then
      qimpc(izone) = qimpat(izone) + qimpfl(izone)
    endif
  enddo

endif


 2001 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: Specific physics module                        ',/,&
'@    =========                       fuel                    ',/,&
'@    issue in boundary conditions                            ',/,&
'@                                                            ',/,&
'@  The outflow is imposed on the zone izone =  ', I10         ,/,&
'@    because                iqimp(izone) =     ', I10         ,/,&
'@  However, on this zone, the                                ',/,&
'@      integred product rho D S is zero:                     ',/,&
'@    it worths                           = ',E14.5            ,/,&
'@    (D is the direction in which the outflow is imposed).   ',/,&
'@                                                            ',/,&
'@  The calcultaion can not run.                              ',/,&
'@                                                            ',/,&
'@  Check boundary condition definitions, in particular that  ',/,&
'@    - the vector rcodcl(ifac,IU,1),                         ',/,&
'@                 rcodcl(ifac,IV,1),                         ',/,&
'@                 rcodcl(ifac,IW,1) which determines         ',/,&
'@      the velocity direction is not zero and is not         ',/,&
'@      uniformly perpendicular to the input face             ',/,&
'@    - the input surface is not zero (or that the number     ',/,&
'@      of edge faces in the zone is not zero)                ',/,&
'@    - the mass density is not zero                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

!===============================================================================
! 3. Verifications
!        Sum of fuel distributions = 100% for the zones ientfl = 1
!===============================================================================

iok = 0
do ii = 1, nzfppp
  izone = ilzppp(ii)
  if ( ientfl(izone).eq.1 ) then
    totfu = 0.d0
    do icla = 1, nclafu
      totfu = totfu + distfu(izone,icla)
    enddo
    if(abs(totfu-100.d0).gt.epzero) then
      write(nfecra,2010)
      do icla = 1, nclafu
        write(nfecra,2011)izone,icla,                             &
               distfu(izone,icla)
      enddo
      write(nfecra,2012)izone,ientfl(izone),                      &
             totfu,totfu-100.d0
      iok = iok + 1
    endif
  endif
enddo

if(iok.ne.0) then
  call csexit (1)
endif


 2010 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: Specific physics modul                         ',/,&
'@    =========                       fuel                    ',/,&
'@    Issue in boundary conditions                            ',/,&
'@                                                            ',/,&
'@        Zone    Class          Distfu(%)                    '  )
 2011 format(                                                           &
'@  ',I10   ,' ',I10   ,'    ',E14.5                             )
 2012 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING; Specific physics modul                         ',/,&
'@    =========                        fuel                   ',/,&
'@    Issue in boundary conditions                            ',/,&
'@                                                            ',/,&
'@  We impose a fuel inflow in izone = ', I10                  ,/,&
'@    because               ientfl(izone) = ', I10            ,/, &
'@  However, on this zone, the distribution sum               ',/,&
'@    in percentage for the fuel ifol = ', I1 0                ,/,&
'@    is different from 100%: it worths totfol = ', E14.5      ,/,&
'@    with                           totfol-100 = ', E14.5     ,/,&
'@                                                            ',/,&
'@  The calculation can not run.                              ',/,&
'@                                                            ',/,&
'@  Check user_fuel_bconds.                                   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

!===============================================================================
! 4.  Filling the boundary conditions table
!       Loop on all input faces
!                     =========================
!         We determine the family and its properties
!           We impose the boundary conditions
!           for the turbulence

!===============================================================================
do ifac = 1, nfabor

  izone = izfppp(ifac)

  if ( itypfb(ifac).eq.ientre ) then

    ! ----  Automatic treatement for turbulence

    if ( icalke(izone).ne.0 ) then

      !       The turbulence is calculated by default if icalke different from 0
      !          - either from the hydrolic diameter, an reference velocity
      !            adapted to current input if icalke = 1
      !          - or from the hydrolic diameter, a reference velocity and
      !            the turbulent intensity adapted to the current input
      !            if icalke = 2

      uref2 = rcodcl(ifac,iu,1)**2                         &
            + rcodcl(ifac,iv,1)**2                         &
            + rcodcl(ifac,iw,1)**2
      uref2 = max(uref2,1.d-12)
      rhomoy = brom(ifac)
      iel    = ifabor(ifac)
      viscla = viscl(iel)
      icke   = icalke(izone)
      dhy    = dh(izone)
      xiturb = xintur(izone)

      if (icke.eq.1) then

        !   Calculation of turbulent inlet conditions using
        !     standard laws for a circular pipe
        !     (their initialization is not needed here but is good practice).
        call turbulence_bc_inlet_hyd_diam(ifac, uref2, dhy, rhomoy, viscla,  &
                                          rcodcl)
      else if (icke.eq.2) then

        ! Calculation of turbulent inlet conditions using
        !   the turbulence intensity and standard laws for a circular pipe
        !   (their initialization is not needed here but is good practice)

        call turbulence_bc_inlet_turb_intensity(ifac, uref2, xiturb, dhy,  &
                                                rcodcl)


      endif

    endif

  endif

enddo

!===============================================================================
! 5.  Filling the boundary conditions table
!     Loop on all input faces
!                     =========================
!     We determine the family and its properties
!     We impose the boundary conditions
!     for the scalars
!===============================================================================

do ii = 1, nzfppp

  izone = ilzppp(ii)

  ! An input ientre must be of type
  ! ientat = 1 or ientfl = 1
  if (ientat(izone).eq.1 .or. ientfl(izone).eq.1) then

    x20t  (izone) = zero
    x2h20t(izone) = zero

    do icla = 1, nclafu

      ! ------ Calculation of total X2 per zone
      !        Small correction in case of an closed inlet
      if(abs(qimpc(izone)).le.epzero) then
        x20(izone,icla) = 0.d0
      else
        x20(izone,icla) = qimpfl(izone)/qimpc(izone)              &
                         *distfu(izone,icla)*1.d-2
      endif
      x20t(izone)     = x20t(izone) +  x20(izone,icla)
    enddo
    ! ------ Calculation of H2, XMG0
    if ( ientfl(izone) .eq. 1 ) then
      t2        = timpfl(izone)
      xsolid(1) = 1.d0-fkc
      xsolid(2) = fkc
      mode      = -1
      call cs_fuel_htconvers2 (mode, h2(izone) , xsolid , t2)
!     =======================

      do icla = 1, nclafu
        xmg0(izone,icla) = pi/6.d0*(dinifl(icla)**3)*rho0fl
      enddo
    else
      h2(izone) = zero
      do icla = 1, nclafu
        xmg0(izone,icla) = 1.d0
      enddo
    endif
    x2h20t(izone) = x20t(izone)*h2(izone)


    ! ------ Calculation of H1(izone)
    do ige = 1, ngazem
      coefe(ige) = zero
    enddo

    ioxy = inmoxy(izone)
    dmas = wmole(io2) *oxyo2(ioxy) +wmole(in2) *oxyn2(ioxy)    &
          +wmole(ih2o)*oxyh2o(ioxy)+wmole(ico2)*oxyco2(ioxy)

    coefe(io2)  = wmole(io2 )*oxyo2(ioxy )/dmas
    coefe(ih2o) = wmole(ih2o)*oxyh2o(ioxy)/dmas
    coefe(ico2) = wmole(ico2)*oxyco2(ioxy)/dmas
    coefe(in2)  = wmole(in2 )*oxyn2(ioxy )/dmas

    hlf = zero
    t1   = timpat(izone)
    mode = -1
    call cs_fuel_htconvers1 (mode, h1(izone) , coefe , t1)
!   =======================

  endif
enddo

do ifac = 1, nfabor

  izone = izfppp(ifac)


  if ( itypfb(ifac).eq.ientre ) then

    ! ----  Automatic treatment for specific physics scalars

    do icla = 1, nclafu
      ! ------ Boundary conditions for Xfol
      rcodcl(ifac,isca(iyfol(icla)),1) = x20(izone,icla)
      ! ------ Boundary conditions for Ng
      rcodcl(ifac,isca(ing(icla)),1) = x20(izone,icla)            &
                                      /xmg0(izone,icla)
      ! ------ Boundary conditions for X2HLF
      rcodcl(ifac,isca(ih2(icla)),1) = x20(izone,icla)*h2(izone)

      if (i_comb_drift.eq.1) then
        rcodcl(ifac, isca(iv_p_x(icla)), 1) = rcodcl(ifac,iu,1)
        rcodcl(ifac, isca(iv_p_y(icla)), 1) = rcodcl(ifac,iv,1)
        rcodcl(ifac, isca(iv_p_z(icla)), 1) = rcodcl(ifac,iw,1)
      endif

    enddo
    ! ------ Boundary conditions for X1.FVAP
    rcodcl(ifac,isca(ifvap),1) = zero
    ! ------ Boundary conditions for X1.F7M
    rcodcl(ifac,isca(if7m),1) = zero
    ! ------ Boundary conditions for X1.Variance
    rcodcl(ifac,isca(ifvp2m),1)   = zero
    ! ------ Boundary conditions for HM
    rcodcl(ifac,isca(iscalt),1) = (1.d0-x20t(izone))*h1(izone)+x2h20t(izone)
    rcodcl(ifac,isca(ihgas),1) = (1.d0-x20t(izone))*h1(izone)

    ! Store the Boundary value of X1
    b_x1(ifac) = (1.d0-x20t(izone))
    ! ------ Boundary conditions for X1.F4M (Oxyd 2)
    if ( noxyd .ge. 2 ) then
      if ( inmoxy(izone) .eq. 2 ) then
        rcodcl(ifac,isca(if4m),1)   = (1.d0-x20t(izone))
      else
        rcodcl(ifac,isca(if4m),1)   = zero
      endif
    endif

    ! ------ Boundary conditions for X1.F5M (Oxyd3)

    if ( noxyd .eq. 3 ) then
      if ( inmoxy(izone) .eq. 3 ) then
        rcodcl(ifac,isca(if5m),1)   = (1.d0-x20t(izone))
      else
        rcodcl(ifac,isca(if5m),1)   = zero
      endif
    endif

    ! ------ Boundary conditions for X1.YCO2
    if ( ieqco2 .ge. 1 ) then
      rcodcl(ifac,isca(iyco2),1)   = zero
    endif

    ! ------ Boundary conditions for X1.HCN and X1.NO
    if ( ieqnox .eq. 1 ) then
      rcodcl(ifac,isca(iyhcn),1)   = zero
      rcodcl(ifac,isca(iyno ),1)   = zero
      rcodcl(ifac,isca(ihox ),1)   = (1.d0-x20t(izone))*h1(izone)
    endif

  endif

  ! Wall BCs on the particle velocity: zero Dirichlet
  if (i_comb_drift.eq.1) then
    if (itypfb(ifac).eq.iparoi.or.itypfb(ifac).eq.iparug) then

      do icla = 1, nclafu

        icodcl(ifac, isca(iv_p_x(icla))) = 1
        icodcl(ifac, isca(iv_p_y(icla))) = 1
        icodcl(ifac, isca(iv_p_z(icla))) = 1
        rcodcl(ifac, isca(iv_p_x(icla)), 1) = 0.d0
        rcodcl(ifac, isca(iv_p_y(icla)), 1) = 0.d0
        rcodcl(ifac, isca(iv_p_z(icla)), 1) = 0.d0
      enddo

    endif
  endif

enddo

!--------
! Formats
!--------

!----
! End
!----

return
end subroutine
