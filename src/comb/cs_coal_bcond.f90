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

!===============================================================================
! Function:
! --------
!> \file cs_coal_bcond.f90
!>
!> \brief Automatic boundary condition for pulverized coal combution
!
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
!> \param[in,out] rcodcl        value of the boundary conditions to edge faces
!>
!>                              boundary condition values:
!>                               - rcodcl(1) value of the dirichlet
!>                               - rcodcl(2) value of the exterior exchange
!>                               -  coefficient (infinite if no exchange)
!>                               -  rcodcl(3) value flux density
!>                               -  (negative if gain) \f$w.m^{-2} \f$ or
!>                               -  roughness in \f$m\f$ if  icodcl=6
!>                                -# for velocity:
!>                                           \f$(\mu+\mu_T)\gradv \vect{u}\f$
!>                                -# for pressure: \f$ \Delta \grad P
!>                                                 \cdot \vect{n} \f$
!>                                -# for scalar:   \f$ C_p \left ( K +
!>                                                 \dfrac{K_T}{\sigma_T} \right)
!>                                                 \grad T \cdot \vect{n} \f$
!______________________________________________________________________________!

subroutine cs_coal_bcond &
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
use cs_coal_incl
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

character(len=80) :: name

integer          ii, ifac, izone, mode, iel, ige, iok
integer          icha, iclapc, isol, icla
integer          icke, idecal
integer          nbrval, ioxy
integer          f_id, keyvar

double precision qisqc, viscla, d2s3, uref2, rhomoy, dhy, xiturb
double precision t1, t2, totcp , dmas
double precision h1    (nozppm) , h2   (nozppm,nclcpm)
double precision x2h20t(nozppm) , x20t (nozppm)
double precision qimpc (nozppm) , qcalc(nozppm)
double precision coefe (ngazem)
double precision xsolid(nsolim)
double precision f1mc  (ncharm) , f2mc (ncharm)
double precision wmh2o,wmco2,wmn2,wmo2
double precision, dimension(:), pointer :: brom, b_x1
integer, dimension (:), allocatable :: iagecp
double precision, dimension(:), pointer :: viscl
!===============================================================================
! 0. Initializations
!===============================================================================
call field_get_val_s(ibrom, brom)
call field_get_val_s(iviscl, viscl)
call field_get_val_s_by_name("b_x_c", b_x1)


d2s3 = 2.d0/3.d0

call field_get_key_id("variable_id", keyvar)

if (iage.ge.1) then

  allocate (iagecp(nclacp))

  do icla = 1, nclacp
    write(name,'(a,i2.2)')'n_p_age_', icla
    call field_get_id(name, f_id)
    call field_get_key_int(f_id, keyvar, iagecp(icla))
  enddo
endif
!===============================================================================
! 1.  Parallel exchanges for the user data
!===============================================================================
!  In fact this exchange could be avoided by changing uscpcl and by asking
!    the user to give the variables which depend of the area out of the loop
!    on the boundary faces: the variables would be available on all processors.
!  However, it makes the user subroutine a bit more complicated and especially
!    if the user modifies it through, it does not work.
!  We assume that all the provided variables are positive,
!    which allows to use a max for the proceedings know them.
!  If this is not the case, it is more complicated but we can get a max anyway.
if(irangp.ge.0) then
  call parimx(nozapm,iqimp )
  call parimx(nozapm,ientat)
  call parimx(nozapm,ientcp)
  call parimx(nozapm,inmoxy)
  call parrmx(nozapm,qimpat)
  call parrmx(nozapm,timpat)
  nbrval = nozppm*ncharm
  call parrmx(nbrval,qimpcp)
  nbrval = nozppm*ncharm
  call parrmx(nbrval,timpcp)
  nbrval = nozppm*ncharm*ncpcmx
  call parrmx(nbrval,distch)
endif

!===============================================================================
! 2.  Correction of the velocities (in norm) for controlling the imposed flow
!       Loop over all inlet faces
!                     =========================
!===============================================================================
! --- Calculated flow
do izone = 1, nozppm
  qcalc(izone) = 0.d0
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

! --- Correction of the velocities (in norm) from the second iteration,
!       otherwise we do not know the density at the edge
if ( ntcabs .gt. 1 ) then
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
      qimpc(izone) = qimpat(izone)
      do icha = 1, ncharb
        qimpc(izone) = qimpc(izone) + qimpcp(izone,icha)
      enddo
      qisqc = qimpc(izone)/qcalc(izone)
      rcodcl(ifac,iu,1) = rcodcl(ifac,iu,1)*qisqc
      rcodcl(ifac,iv,1) = rcodcl(ifac,iv,1)*qisqc
      rcodcl(ifac,iw,1) = rcodcl(ifac,iw,1)*qisqc
    endif
  enddo

else

  do izone = 1, nozapm
    qimpc(izone) = qimpat(izone)
    do icha = 1, ncharb
      qimpc(izone) = qimpc(izone) + qimpcp(izone,icha)
    enddo
  enddo

endif


 2001 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING : SPECIFIC PHYSIC MODULE                        ',/,&
'@    =========                        pulverized coal        ',/,&
'@    problem in the boundary conditions                      ',/,&
'@                                                            ',/,&
'@  The flow rate is imposed on the area izone =  ', I10       ,/,&
'@    because                iqimp(izone) =     ', I10         ,/,&
'@  However, on this area, the integrated product rho D S     ',/,&
'@    is zero                             :                   ',/,&
'@    it is                               = ',E14.5            ,/,&
'@    (D is the direction in which the flow is imposed).      ',/,&
'@                                                            ',/,&
'@  The calculation can not be executed                       ',/,&
'@                                                            ',/,&
'@  Check uscpcl, and in particular                           ',/,&
'@    - that the vector  rcodcl(ifac,iu,1)                    ',/,&
'@                       rcodcl(ifac,iv,1),                   ',/,&
'@                       rcodcl ifac,iw,1) which determines   ',/,&
'@      the direction of the velocity is not zero and is not  ',/,&
'@      uniformly perpendicular to the imput faces            ',/,&
'@    - that the surface of the imput is not zero (or the     ',/,&
'@      number of edge faces in the area is non-zero)         ',/,&
'@    - that the density is not zero                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

!===============================================================================
! 3. Verifications
!        Sum coal distribution = 100% for area ientcp = 1
!===============================================================================

iok = 0
do ii = 1, nzfppp
  izone = ilzppp(ii)
  if ( ientcp(izone).eq.1 ) then
    do icha = 1, ncharb
      totcp = 0.d0
      do iclapc = 1, nclpch(icha)
        totcp = totcp + distch(izone,icha,iclapc)
      enddo
      if(abs(totcp-100.d0).gt.epzero) then
        write(nfecra,2010)
        do iclapc = 1, nclpch(icha)
          write(nfecra,2011)izone,icha,iclapc,                    &
               distch(izone,icha,iclapc)
        enddo
        write(nfecra,2012)izone,ientcp(izone),icha,               &
             totcp,totcp-100.d0
        iok = iok + 1
      endif
    enddo
  endif
enddo

if(iok.ne.0) then
  call csexit (1)
endif


 2010 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING : SPECIFIC PHYSIC MODULE                        ',/,&
'@    =========                        pulverized coal        ',/,&
'@    probleme in the boundary conditions                     ',/,&
'@                                                            ',/,&
'@        Zone    Coal     Class         Distch(%)        '  )
 2011 format(                                                           &
'@  ',I10   ,' ',I10   ,' ',I10   ,'    ',E14.5                  )
 2012 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING : SPECIFIC PHYSIC MODULE                        ',/,&
'@    =========                        pulverized coal        ',/,&
'@    probleme in the boundary conditions                     ',/,&
'@                                                            ',/,&
'@  A coal input is imposed in izone = ', I10                  ,/,&
'@    because               ientcp(izone) = ', I10             ,/,&
'@  However, on this area, the sum of distributions by class  ',/,&
'@    in percentage for coal         icha = ', I10             ,/,&
'@    is different from 100% : it is     totcp = ', E14.5      ,/,&
'@    with                           totcp-100 = ', E14.5      ,/,&
'@                                                            ',/,&
'@  The calcul will not run                                   ',/,&
'@                                                            ',/,&
'@  Check    uscpcl.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

!===============================================================================
! 4.  Filling the table of the boundary conditions
!       Loop on all input faces
!                     =========================
!         Determining the family and its properties
!         Imposing boundary conditions for the turbulence

!===============================================================================
do ifac = 1, nfabor

  izone = izfppp(ifac)

  if ( itypfb(ifac).eq.ientre ) then

    !       The turbulence is calculated by default if icalke different from 0
    !          - or from hydraulic diameter and a reference velocity adapted
    !            for the current input if icalke = 1
    !          - either from the hydraulic diameter, a reference velocity and
    !            a turbulence intensity adapted to the current input if icalke = 2
    if ( icalke(izone).ne.0 ) then

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
  ! ientat = 1 or ientcp = 1
  if ( ientat(izone).eq.1 .or. ientcp(izone).eq.1) then

    x20t  (izone) = zero
    x2h20t(izone) = zero

    idecal = 0

    do icha = 1, ncharb

      do iclapc = 1, nclpch(icha)

        icla = iclapc + idecal
        ! ------ Calculation of total X2 per zone
        !        Small correction in case of an closed inlet
        if(abs(qimpc(izone)).lt.epzero) then
          x20(izone,icla) = 0.d0
        else
          x20(izone,icla) = qimpcp(izone,icha)/qimpc(izone)       &
                          * distch(izone,icha,iclapc)*1.d-2
        endif
        x20t(izone)     = x20t(izone) +  x20(izone,icla)
        ! ------ Calculating H2 of class icla
        do isol = 1, nsolim
          xsolid(isol) = zero
        enddo
        if ( ientcp(izone).eq.1 ) then
          t2  = timpcp(izone,icha)
          xsolid(ich(icha)) = 1.d0-xashch(icha)
          xsolid(ick(icha)) = zero
          xsolid(iash(icha)) = xashch(icha)
          !------- Taking into account humidity
          if ( ippmod(iccoal) .eq. 1 ) then
            xsolid(ich(icha)) = xsolid(ich(icha))-xwatch(icha)
            xsolid(iwat(icha)) = xwatch(icha)
          else
            xsolid(iwat(icha)) = 0.d0
          endif

        else
          t2  = timpat(izone)

          xsolid(ich(icha))  = (1.d0-xashch(icha)-xwatch(icha))
          xsolid(ick(icha))  = 0.d0
          xsolid(iash(icha)) = xashch(icha)
          xsolid(iwat(icha)) = xwatch(icha)

        endif
        mode = -1
        t1 = t2
        call cs_coal_htconvers2(mode,icla,h2(izone,icla),xsolid,t2,t1)
        x2h20t(izone) = x2h20t(izone)+x20(izone,icla)*h2(izone,icla)

      enddo

      idecal = idecal + nclpch(icha)

    enddo

    ! ------ Calculating H1(izone)
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

    do icha = 1, ncharm
      f1mc(icha) = zero
      f2mc(icha) = zero
    enddo
    t1   = timpat(izone)
    mode = -1
    call cs_coal_htconvers1(mode,h1(izone),coefe,f1mc,f2mc,t1)

  endif
enddo

do ifac = 1, nfabor

  izone = izfppp(ifac)

  if ( itypfb(ifac).eq.ientre ) then

    ! ----  Automatic processing of specific physic scalar

    idecal = 0

    do icha = 1, ncharb

      do iclapc = 1, nclpch(icha)

        icla = iclapc + idecal

        ! ------ Boundary conditions for Xch of class icla
        rcodcl(ifac,isca(ixch(icla)),1) = x20(izone,icla)         &
                                        * (1.d0-xashch(icha))
        !             Taking into account humidity
        if ( ippmod(iccoal) .eq. 1 ) then
          rcodcl(ifac,isca(ixch(icla)),1) = x20(izone,icla)       &
                                          *(1.d0-xashch(icha)     &
                                                -xwatch(icha))
        endif
        ! ------ Boundary conditions for Xck of class icla
        rcodcl(ifac,isca(ixck(icla)),1) = 0.d0

        ! ------ Boundary conditions for Np of class icla

        rcodcl(ifac,isca(inp(icla)),1) = x20(izone,icla)          &
                                        / xmp0(icla)

        ! ------ Boundary conditions for Xwater of class icla

        if ( ippmod(iccoal) .eq. 1 ) then
          rcodcl(ifac,isca(ixwt(icla)),1) = x20(izone,icla)       &
                                           *xwatch(icha)
        endif

        ! ------ Boundary conditions for H2 of class icla

        rcodcl(ifac,isca(ih2(icla)),1) = x20(izone,icla)          &
                                        *h2(izone,icla)

        ! Boundary conditions for particle age
        if (i_comb_drift.ge.1) then
          rcodcl(ifac, iagecp(icla), 1) = 0.d0
        endif
        if (i_comb_drift.eq.1) then
          rcodcl(ifac, isca(iv_p_x(icla)), 1) = rcodcl(ifac,iu,1)
          rcodcl(ifac, isca(iv_p_y(icla)), 1) = rcodcl(ifac,iv,1)
          rcodcl(ifac, isca(iv_p_z(icla)), 1) = rcodcl(ifac,iw,1)
        endif
      enddo

      idecal = idecal + nclpch(icha)

      ! ------ Boundary conditions for X1F1M and X1F2M from coal icha

      rcodcl(ifac,isca(if1m(icha)),1) = zero
      rcodcl(ifac,isca(if2m(icha)),1) = zero

    enddo
    ! Boundary condition for the age
    if (iage.ge.1) then
      rcodcl(ifac, isca(iage), 1) = zero
    endif

    ! ------ Boundary conditions for HM
    rcodcl(ifac,isca(iscalt),1) = (1.d0-x20t(izone))*h1(izone)    &
                                 +x2h20t(izone)
    ! Boundary condition for x1*h1
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

    ! ------ Boundary conditions for X1.F6M (Water)

    if ( ippmod(iccoal) .ge. 1 ) then
      rcodcl(ifac,isca(if6m),1) = zero
    endif

    ! ------ Boundary conditions for X1.F7M_O2

    rcodcl(ifac,isca(if7m),1)   = zero

    ! ------ Boundary conditions for X1.FM8_CO2

    if ( ihtco2 .eq. 1 ) then
      rcodcl(ifac,isca(if8m),1) = zero
    endif
    ! ------ Boundary conditions for X1.FM9_H2O
    if ( ihth2o .eq. 1 ) then
      rcodcl(ifac,isca(if9m),1) = zero
    endif
    ! ------ Boundary conditions for X1.Variance
    rcodcl(ifac,isca(ifvp2m),1) = zero

    ! ------ Boundary conditions for X1.YCO2
    if ( ieqco2 .eq. 1 ) then
      ioxy =  inmoxy(izone)
      wmo2   = wmole(io2)
      wmco2  = wmole(ico2)
      wmh2o  = wmole(ih2o)
      wmn2   = wmole(in2)
      dmas = ( oxyo2 (ioxy)*wmo2 +oxyn2 (ioxy)*wmn2               &
              +oxyh2o(ioxy)*wmh2o+oxyco2(ioxy)*wmco2 )
      xco2 = oxyco2(ioxy)*wmco2/dmas
      rcodcl(ifac,isca(iyco2),1)   = xco2*(1.d0-x20t(izone))
    endif

    ! ------ Boundary conditions for X1.HCN, X1.NO, T air
    if( ieqnox .eq. 1 ) then
      rcodcl(ifac,isca(iyhcn ),1)  = zero
      rcodcl(ifac,isca(iyno  ),1)  = zero
      rcodcl(ifac,isca(iynh3 ),1)  = zero
      rcodcl(ifac,isca(ihox  ),1)  = (1.d0-x20t(izone))*h1(izone)
    endif

  endif

  ! Wall BCs on the particle velocity: zero Dirichlet
  if (itypfb(ifac).eq.iparoi.or.itypfb(ifac).eq.iparug) then

    idecal = 0

    do icha = 1, ncharb

      do iclapc = 1, nclpch(icha)

        icla = iclapc + idecal

        if (i_comb_drift.eq.1) then
          icodcl(ifac, isca(iv_p_x(icla))) = 1
          icodcl(ifac, isca(iv_p_y(icla))) = 1
          icodcl(ifac, isca(iv_p_z(icla))) = 1
          rcodcl(ifac, isca(iv_p_x(icla)), 1) = 0.d0
          rcodcl(ifac, isca(iv_p_y(icla)), 1) = 0.d0
          rcodcl(ifac, isca(iv_p_z(icla)), 1) = 0.d0
        endif
      enddo

    enddo
  endif
enddo

! Free memory
if (allocated(iagecp)) deallocate(iagecp)
!--------
! Formats
!--------

!----
! End
!----

return
end subroutine
