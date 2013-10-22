!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2013 EDF S.A.
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
! Function :
! --------

!> \file clsyvt.f90
!>
!> \brief Symmetry boundary conditions for vectors and tensors.
!>
!> Correspond to the code icodcl(ivar) = 4.
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nscal         total number of scalars
!> \param[in,out] icodcl        face boundary condition code:
!>                               - 1 Dirichlet
!>                               - 3 Neumann
!>                               - 4 sliding and
!>                                 \f$ \vect{u} \cdot \vect{n} = 0 \f$
!>                               - 5 smooth wall and
!>                                 \f$ \vect{u} \cdot \vect{n} = 0 \f$
!>                               - 6 rought wall and
!>                                 \f$ \vect{u} \cdot \vect{n} = 0 \f$
!>                               - 9 free inlet/outlet
!>                                 (input mass flux blocked to 0)
!> \param[in]     propce        physical properties at cell centers
!> \param[in,out] rcodcl        boundary condition values:
!>                               - rcodcl(1) value of the dirichlet
!>                               - rcodcl(2) value of the exterior exchange
!>                                 coefficient (infinite if no exchange)
!>                               - rcodcl(3) value flux density
!>                                 (negative if gain) in w/m2 or roughtness
!>                                 in m if icodcl=6
!>                                 -# for the velocity \f$ (\mu+\mu_T)
!>                                    \gradv \vect{u} \cdot \vect{n}  \f$
!>                                 -# for the pressure \f$ \Delta t
!>                                    \grad P \cdot \vect{n}  \f$
!>                                 -# for a scalar \f$ cp \left( K +
!>                                     \dfrac{K_T}{\sigma_T} \right)
!>                                     \grad T \cdot \vect{n} \f$
!> \param[in]     velipb        value of the velocity at \f$ \centip \f$
!>                               of boundary cells
!> \param[in]     rijipb        value of \f$ R_{ij} \f$ at \f$ \centip \f$
!>                               of boundary cells
!> \param[out]    coefa         explicit boundary condition coefficient
!> \param[out]    coefb         implicit boundary condition coefficient
!_______________________________________________________________________________


subroutine clsyvt &
 ( nscal  , icodcl ,                                     &
   propce , rcodcl ,                                     &
   velipb , rijipb , coefa  , coefb  )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use optcal
use cstphy
use cstnum
use pointe
use entsor
use albase
use field
use mesh

!===============================================================================

implicit none

! Arguments

integer          nscal

integer          icodcl(nfabor,nvarcl)

double precision propce(ncelet,*)
double precision rcodcl(nfabor,nvarcl,3)
double precision velipb(nfabor,ndim), rijipb(nfabor,6)
double precision coefa(nfabor,*), coefb(nfabor,*)

! Local variables

integer          ifac, ii, isou, jsou
integer          icl11 , icl22 , icl33 , icl12 , icl13 , icl23
integer          icl11r, icl22r, icl33r, icl12r, icl13r, icl23r
integer          icl11f, icl22f, icl33f, icl12f, icl13f, icl23f
integer          iclvar, iel   , iclvrr, iclvaf
integer          iscal , ipccp , ivar
integer          f_id

double precision rnx, rny, rnz, rxnn
double precision upx, upy, upz, usn
double precision tx, ty, tz, txn, t2x, t2y, t2z
double precision clsyme
double precision eloglo(3,3), alpha(6,6)
double precision srfbnf, rcodcn, hint, visclc, visctc, distbf
double precision cpp,rkl
double precision vismsh(3), hintv(3)
double precision hintt(6)
double precision visci(3,3), fikis, viscis, distfi

character*80     fname

double precision, dimension(:,:), pointer :: coefaut, cofafut, cofarut
double precision, dimension(:,:,:), pointer :: coefbut, cofbfut, cofbrut

!===============================================================================

!===============================================================================
! 1. Initializations
!===============================================================================

! Initialize variables to avoid compiler warnings

icl11 = 0
icl22 = 0
icl33 = 0
icl12 = 0
icl13 = 0
icl23 = 0
icl11r = 0
icl12r = 0
icl13r = 0
icl22r = 0
icl23r = 0
icl33r = 0
iclvar = 0

cpp = 0.d0

! --- Gradient Boundary Conditions
if (itytur.eq.3) then
  icl11  = iclrtp(ir11,icoef)
  icl22  = iclrtp(ir22,icoef)
  icl33  = iclrtp(ir33,icoef)
  icl12  = iclrtp(ir12,icoef)
  icl13  = iclrtp(ir13,icoef)
  icl23  = iclrtp(ir23,icoef)
  icl11r = iclrtp(ir11,icoefr)
  icl22r = iclrtp(ir22,icoefr)
  icl33r = iclrtp(ir33,icoefr)
  icl12r = iclrtp(ir12,icoefr)
  icl13r = iclrtp(ir13,icoefr)
  icl23r = iclrtp(ir23,icoefr)
endif

! --- Flux Boundary Conditions
if (itytur.eq.3) then
  icl11f = iclrtp(ir11,icoeff)
  icl22f = iclrtp(ir22,icoeff)
  icl33f = iclrtp(ir33,icoeff)
  icl12f = iclrtp(ir12,icoeff)
  icl13f = iclrtp(ir13,icoeff)
  icl23f = iclrtp(ir23,icoeff)
endif

! --- Begin the loop over boundary faces
do ifac = 1, nfabor

  ! --- Test sur la presence d'une condition de symetrie vitesse : debut
  if (icodcl(ifac,iu).eq.4) then

    ! --- To cancel the mass flux
    isympa(ifac) = 0

    ! Geometric quantities
    srfbnf = surfbn(ifac)

    !===========================================================================
    ! 1. Local framework
    !===========================================================================

    ! Unit normal

    rnx = surfbo(1,ifac)/srfbnf
    rny = surfbo(2,ifac)/srfbnf
    rnz = surfbo(3,ifac)/srfbnf

    ! En ALE, on a eventuellement une vitesse de deplacement de la face
    !   donc seule la composante normale importe (on continue a determiner
    !   TX a partir de la vitesse tangentielle absolue car l'orientation
    !   de TX et T2X est sans importance pour les symetries)
    rcodcn = 0.d0
    if (iale.eq.1) then
      rcodcn = rcodcl(ifac,iu,1)*rnx                           &
             + rcodcl(ifac,iv,1)*rny                           &
             + rcodcl(ifac,iw,1)*rnz
    endif

    upx = velipb(ifac,1)
    upy = velipb(ifac,2)
    upz = velipb(ifac,3)

    if (itytur.eq.3) then

      ! Relative tangential velocity

      usn = upx*rnx+upy*rny+upz*rnz
      tx  = upx -usn*rnx
      ty  = upy -usn*rny
      tz  = upz -usn*rnz
      txn = sqrt( tx**2 +ty**2 +tz**2 )

      ! Unit tangent

      if( txn.ge.epzero) then

        tx = tx/txn
        ty = ty/txn
        tz = tz/txn

      else

        ! If the velocity is zero, vector T is normal and random;
        !   we need it for the reference change for Rij, and we cancel the velocity.

        if(abs(rny).ge.epzero.or.abs(rnz).ge.epzero)then
          rxnn = sqrt(rny**2+rnz**2)
          tx  =  0.d0
          ty  =  rnz/rxnn
          tz  = -rny/rxnn
        elseif(abs(rnx).ge.epzero.or.abs(rnz).ge.epzero)then
          rxnn = sqrt(rnx**2+rnz**2)
          tx  =  rnz/rxnn
          ty  =  0.d0
          tz  = -rnx/rxnn
        else
          write(nfecra,1000)ifac,rnx,rny,rnz
          call csexit (1)
        endif

      endif


      ! --> T2 = RN X T (where X is the cross product)

      t2x = rny*tz - rnz*ty
      t2y = rnz*tx - rnx*tz
      t2z = rnx*ty - rny*tx

      ! --> Orthogonal matrix for change of reference frame ELOGLOij
      !     (from local to global reference frame)

      !                      |TX  -RNX  T2X|
      !             ELOGLO = |TY  -RNY  T2Y|
      !                      |TZ  -RNZ  T2Z|

      !    Its transpose ELOGLOt is its inverse

      eloglo(1,1) =  tx
      eloglo(1,2) = -rnx
      eloglo(1,3) =  t2x
      eloglo(2,1) =  ty
      eloglo(2,2) = -rny
      eloglo(2,3) =  t2y
      eloglo(3,1) =  tz
      eloglo(3,2) = -rnz
      eloglo(3,3) =  t2z

      ! --> Commpute alpha(6,6)

      ! Let f be the center of the boundary faces and
      !   I the center of the matching cell

      ! We noteE Rg (resp. Rl) indexed by f or by I
      !   the Reynolds Stress tensor in the global basis (resp. local)

      ! The alpha matrix applied to the global vector in I'
      !   (Rg11,I'|Rg22,I'|Rg33,I'|Rg12,I'|Rg13,I'|Rg23,I')t
      !    must provide the values to prescribe to the face
      !   (Rg11,f |Rg22,f |Rg33,f |Rg12,f |Rg13,f |Rg23,f )t
      !    except for the Dirichlet boundary conditions (added later)

      ! We define it by computing Rg,f as a function of Rg,I' as follows

      !   RG,f = ELOGLO.RL,f.ELOGLOt (matrix products)

      !                     | RL,I'(1,1)     B*U*.Uk     C*RL,I'(1,3) |
      !      with    RL,f = | B*U*.Uk       RL,I'(2,2)       0        |
      !                     | C*RL,I'(1,3)     0         RL,I'(3,3)   |

      !             with    RL,I = ELOGLOt.RG,I'.ELOGLO
      !                     B = 0
      !              and    C = 0 at the wall (1 with symmetry)

      ! We compute in fact  ELOGLO.projector.ELOGLOt

      clsyme=1.d0
      call clca66 ( clsyme , eloglo , alpha )
      !==========

    endif

    !===========================================================================
    ! 2. Boundary conditions on the velocity
    !    (totaly implicit)
    !    The condition is a zero (except in ALE) Dirichlet on the normal component
    !                     a homogenous Neumann on the other components
    !===========================================================================

    ! Coupled solving of the velocity components

    iel = ifabor(ifac)
    ! --- Physical properties
    visclc = propce(iel,ipproc(iviscl))
    visctc = propce(iel,ipproc(ivisct))

    ! --- Geometrical quantity
    distbf = distb(ifac)

    if (itytur.eq.3) then
      hint =   visclc         /distbf
    else
      hint = ( visclc+visctc )/distbf
    endif

    ! Gradient BCs
    coefau(1,ifac) = rcodcn*rnx
    coefau(2,ifac) = rcodcn*rny
    coefau(3,ifac) = rcodcn*rnz

    coefbu(1,1,ifac) = 1.d0-rnx**2
    coefbu(2,2,ifac) = 1.d0-rny**2
    coefbu(3,3,ifac) = 1.d0-rnz**2

    coefbu(1,2,ifac) = -rnx*rny
    coefbu(1,3,ifac) = -rnx*rnz
    coefbu(2,1,ifac) = -rny*rnx
    coefbu(2,3,ifac) = -rny*rnz
    coefbu(3,1,ifac) = -rnz*rnx
    coefbu(3,2,ifac) = -rnz*rny

    ! Flux BCs
    cofafu(1,ifac) = -hint*rcodcn*rnx
    cofafu(2,ifac) = -hint*rcodcn*rny
    cofafu(3,ifac) = -hint*rcodcn*rnz

    cofbfu(1,1,ifac) = hint*rnx**2
    cofbfu(2,2,ifac) = hint*rny**2
    cofbfu(3,3,ifac) = hint*rnz**2

    cofbfu(1,2,ifac) = hint*rnx*rny
    cofbfu(1,3,ifac) = hint*rnx*rnz
    cofbfu(2,1,ifac) = hint*rny*rnx
    cofbfu(2,3,ifac) = hint*rny*rnz
    cofbfu(3,1,ifac) = hint*rnz*rnx
    cofbfu(3,2,ifac) = hint*rnz*rny

    !===========================================================================
    ! 3. Boundary conditions on Rij (partially implicited)
    !===========================================================================

    if (itytur.eq.3) then

      do isou = 1, 6

        if (isou.eq.1) then
          iclvar = icl11
          iclvrr = icl11r
        elseif (isou.eq.2) then
          iclvar = icl22
          iclvrr = icl22r
        elseif (isou.eq.3) then
          iclvar = icl33
          iclvrr = icl33r
        elseif (isou.eq.4) then
          iclvar = icl12
          iclvrr = icl12r
        elseif (isou.eq.5) then
          iclvar = icl13
          iclvrr = icl13r
        elseif (isou.eq.6) then
          iclvar = icl23
          iclvrr = icl23r
        endif

        coefa(ifac,iclvar) = 0.0d0
        coefa(ifac,iclvrr) = 0.0d0
        coefb(ifac,iclvar) = 0.0d0
        coefb(ifac,iclvrr) = 0.0d0

      enddo

      do isou = 1,6

        if (isou.eq.1) then
          ivar = ir11
          iclvar = icl11
          iclvaf = icl11f
          iclvrr = icl11r
        elseif (isou.eq.2) then
          ivar = ir22
          iclvar = icl22
          iclvaf = icl22f
          iclvrr = icl22r
        elseif (isou.eq.3) then
          ivar = ir33
          iclvar = icl33
          iclvaf = icl33f
          iclvrr = icl33r
        elseif (isou.eq.4) then
          ivar = ir12
          iclvar = icl12
          iclvaf = icl12f
          iclvrr = icl12r
        elseif (isou.eq.5) then
          ivar = ir13
          iclvar = icl13
          iclvaf = icl13f
          iclvrr = icl13r
        elseif (isou.eq.6) then
          ivar = ir23
          iclvar = icl23
          iclvaf = icl23f
          iclvrr = icl23r
        endif

        ! IMPLICITATION PARTIELLE EVENTUELLE DES CL
        if (iclsyr.eq.1) then
          do ii = 1, 6
            if (ii.ne.isou) then
              coefa(ifac,iclvar) = coefa(ifac,iclvar) +           &
                   alpha(isou,ii) * rijipb(ifac,ii)
            endif
          enddo
          coefb(ifac,iclvar) = alpha(isou,isou)
        else
          do ii = 1, 6
            coefa(ifac,iclvar) = coefa(ifac,iclvar) +             &
                 alpha(isou,ii) * rijipb(ifac,ii)
          enddo
          coefb(ifac,iclvar) = 0.d0
        endif

        ! Boundary conditions for the momentum equation
        coefa(ifac,iclvrr) = coefa(ifac,iclvar)
        coefb(ifac,iclvrr) = coefb(ifac,iclvar)

        ! Symmetric tensor diffusivity (Daly Harlow -- GGDH)
        if (idften(ivar).eq.6) then

          visci(1,1) = visclc + visten(1,iel)
          visci(2,2) = visclc + visten(2,iel)
          visci(3,3) = visclc + visten(3,iel)
          visci(1,2) =          visten(4,iel)
          visci(2,1) =          visten(4,iel)
          visci(2,3) =          visten(5,iel)
          visci(3,2) =          visten(5,iel)
          visci(1,3) =          visten(6,iel)
          visci(3,1) =          visten(6,iel)

          ! ||Ki.S||^2
          viscis = ( visci(1,1)*surfbo(1,ifac)       &
                   + visci(1,2)*surfbo(2,ifac)       &
                   + visci(1,3)*surfbo(3,ifac))**2   &
                 + ( visci(2,1)*surfbo(1,ifac)       &
                   + visci(2,2)*surfbo(2,ifac)       &
                   + visci(2,3)*surfbo(3,ifac))**2   &
                 + ( visci(3,1)*surfbo(1,ifac)       &
                   + visci(3,2)*surfbo(2,ifac)       &
                   + visci(3,3)*surfbo(3,ifac))**2

          ! IF.Ki.S
          fikis = ( (cdgfbo(1,ifac)-xyzcen(1,iel))*visci(1,1)   &
                  + (cdgfbo(2,ifac)-xyzcen(2,iel))*visci(2,1)   &
                  + (cdgfbo(3,ifac)-xyzcen(3,iel))*visci(3,1)   &
                  )*surfbo(1,ifac)                              &
                + ( (cdgfbo(1,ifac)-xyzcen(1,iel))*visci(1,2)   &
                  + (cdgfbo(2,ifac)-xyzcen(2,iel))*visci(2,2)   &
                  + (cdgfbo(3,ifac)-xyzcen(3,iel))*visci(3,2)   &
                  )*surfbo(2,ifac)                              &
                + ( (cdgfbo(1,ifac)-xyzcen(1,iel))*visci(1,3)   &
                  + (cdgfbo(2,ifac)-xyzcen(2,iel))*visci(2,3)   &
                  + (cdgfbo(3,ifac)-xyzcen(3,iel))*visci(3,3)   &
                  )*surfbo(3,ifac)

          distfi = distb(ifac)

          ! Take I" so that I"F= eps*||FI||*Ki.n when J" is in cell rji
          ! NB: eps =1.d-1 must be consistent with vitens.f90
          fikis = max(fikis, 1.d-1*sqrt(viscis)*distfi)

          hint = viscis/surfbn(ifac)/fikis

        else
          call csexit(1)
        endif

        ! Translate into Diffusive flux BCs
        coefa(ifac,iclvaf) = -hint*coefa(ifac,iclvar)
        coefb(ifac,iclvaf) = hint*(1.d0-coefb(ifac,iclvar))

      enddo

    endif

    !===========================================================================
    ! 3.bis Boundary conditions on u'T'
    !===========================================================================

    do iscal = 1, nscal

      if (ityturt(iscal).eq.3) then
        ivar = isca(iscal)
        ! Name of the scalar ivar !TODO move outside of the loop
        call field_get_name(ivarfl(ivar), fname)

        ! Index of the corresponding turbulent flux
        call field_get_id(trim(fname)//'_turbulent_flux', f_id)

        call field_get_coefa_v(f_id,coefaut)
        call field_get_coefb_v(f_id,coefbut)
        call field_get_coefaf_v(f_id,cofafut)
        call field_get_coefbf_v(f_id,cofbfut)
        call field_get_coefad_v(f_id,cofarut)
        call field_get_coefbd_v(f_id,cofbrut)

        iel = ifabor(ifac)
        ! --- Physical Propreties
        visclc = propce(iel,ipproc(iviscl))
        if (icp.gt.0) then
          ipccp  = ipproc(icp)
        else
          ipccp = 0
        endif
        cpp = 1.d0
        if (iscacp(iscal).eq.1) then
          if (ipccp.gt.0) then
            cpp = propce(iel,ipccp)
          else
            cpp = cp0
          endif
        endif

        ! --- Geometrical quantities
        distbf = distb(ifac)

        if (ivisls(iscal).le.0) then
          rkl = visls0(iscal)/cpp
        else
          rkl = propce(iel,ipproc(ivisls(iscal)))/cpp
        endif

        hintt(1) = 0.5d0*(visclc+rkl)/distbf                        &
                 + visten(1,iel)*ctheta(iscal)/distbf
        hintt(2) = 0.5d0*(visclc+rkl)/distbf                        &
                 + visten(2,iel)*ctheta(iscal)/distbf
        hintt(3) = 0.5d0*(visclc+rkl)/distbf                        &
                 + visten(3,iel)*ctheta(iscal)/distbf
        hintt(4) = visten(4,iel)*ctheta(iscal)/distbf
        hintt(5) = visten(5,iel)*ctheta(iscal)/distbf
        hintt(6) = visten(6,iel)*ctheta(iscal)/distbf

        ! Gradient BCs
        coefaut(1,ifac) = 0.d0
        coefaut(2,ifac) = 0.d0
        coefaut(3,ifac) = 0.d0

        coefbut(1,1,ifac) = 1.d0-rnx**2
        coefbut(2,2,ifac) = 1.d0-rny**2
        coefbut(3,3,ifac) = 1.d0-rnz**2

        coefbut(1,2,ifac) = -rnx*rny
        coefbut(1,3,ifac) = -rnx*rnz
        coefbut(2,1,ifac) = -rny*rnx
        coefbut(2,3,ifac) = -rny*rnz
        coefbut(3,1,ifac) = -rnz*rnx
        coefbut(3,2,ifac) = -rnz*rny

        ! Flux BCs
        cofafut(1,ifac) = 0.d0
        cofafut(2,ifac) = 0.d0
        cofafut(3,ifac) = 0.d0

        cofbfut(1,1,ifac) = hintt(1)*rnx**2  + hintt(4)*rnx*rny + hintt(6)*rnx*rnz
        cofbfut(2,2,ifac) = hintt(4)*rnx*rny + hintt(2)*rny**2  + hintt(5)*rny*rnz
        cofbfut(3,3,ifac) = hintt(6)*rnx*rnz + hintt(5)*rny*rnz + hintt(3)*rnz**2

        cofbfut(1,2,ifac) = hintt(1)*rnx*rny + hintt(4)*rny**2  + hintt(6)*rny*rnz
        cofbfut(2,1,ifac) = hintt(1)*rnx*rny + hintt(4)*rny**2  + hintt(6)*rny*rnz
        cofbfut(1,3,ifac) = hintt(1)*rnx*rnz + hintt(4)*rny*rnz + hintt(6)*rnz**2
        cofbfut(3,1,ifac) = hintt(1)*rnx*rnz + hintt(4)*rny*rnz + hintt(6)*rnz**2
        cofbfut(2,3,ifac) = hintt(4)*rnx*rnz + hintt(2)*rny*rnz + hintt(5)*rnz**2
        cofbfut(3,2,ifac) = hintt(4)*rnx*rnz + hintt(2)*rny*rnz + hintt(5)*rnz**2

        ! Boundary conditions for thermal transport equation
        do isou = 1, 3
          cofarut(isou,ifac) = coefaut(isou,ifac)
          do jsou =1, 3
            cofbrut(isou,jsou,ifac) = coefbut(isou,jsou,ifac)
          enddo
        enddo

      endif

    enddo

  endif
! --- Test sur la presence d'une condition de symetrie vitesse : fin

enddo
! ---  End of loop over boundary faces

!===============================================================================
! 4. Symmetry boundary conditions for mesh velocity (ALE module)
!===============================================================================

if (iale.eq.1) then

  do ifac = 1, nfabor
    if (icodcl(ifac,iuma).eq.4) then

      iel = ifabor(ifac)

      ! For a sliding boundary, the normal velocity is enforced to zero
      ! whereas the other components have an Homogenous Neumann
      ! NB: no recontruction in I' here

      ! --- Geometrical quantity
      distbf = distb(ifac)
      srfbnf = surfbn(ifac)

      ! --- Physical properties
      vismsh(1) = propce(iel, ipproc(ivisma(1)))
      vismsh(2) = propce(iel, ipproc(ivisma(2)))
      vismsh(3) = propce(iel, ipproc(ivisma(3)))

      hintv(1) = vismsh(1)/distbf
      hintv(2) = vismsh(2)/distbf
      hintv(3) = vismsh(3)/distbf

      ! Unit normal
      rnx = surfbo(1,ifac)/srfbnf
      rny = surfbo(2,ifac)/srfbnf
      rnz = surfbo(3,ifac)/srfbnf

      ! Coupled solving of the velocity components

      ! Gradient BCs
      claale(1,ifac) = 0.d0
      claale(2,ifac) = 0.d0
      claale(3,ifac) = 0.d0

      clbale(1,1,ifac) = 1.d0-rnx**2
      clbale(2,2,ifac) = 1.d0-rny**2
      clbale(3,3,ifac) = 1.d0-rnz**2

      clbale(1,2,ifac) = -rnx*rny
      clbale(2,1,ifac) = -rny*rnx
      clbale(1,3,ifac) = -rnx*rnz
      clbale(3,1,ifac) = -rnz*rnx
      clbale(2,3,ifac) = -rny*rnz
      clbale(3,2,ifac) = -rnz*rny

      ! Flux BCs
      cfaale(1,ifac) = 0.d0
      cfaale(2,ifac) = 0.d0
      cfaale(3,ifac) = 0.d0

      cfbale(1,1,ifac) = hintv(1)*rnx**2
      cfbale(2,2,ifac) = hintv(2)*rny**2
      cfbale(3,3,ifac) = hintv(3)*rnz**2

      cfbale(1,2,ifac) = hintv(1)*rnx*rny
      cfbale(2,1,ifac) = hintv(2)*rny*rnx
      cfbale(1,3,ifac) = hintv(1)*rnx*rnz
      cfbale(3,1,ifac) = hintv(3)*rnz*rnx
      cfbale(2,3,ifac) = hintv(2)*rny*rnz
      cfbale(3,2,ifac) = hintv(3)*rnz*rny

    endif
  enddo

endif

!===============================================================================
! 5. Formats
!===============================================================================

#if defined(_CS_LANG_FR)

 1000 format(/,' LA NORMALE A LA FACE DE BORD DE SYMETRIE ',I10,/,&
         ' EST NULLE ; COORDONNEES : ',3E12.5)

#else

 1000 format(/,' THE NORMAL TO THE SYMMETRY BOUNDARY FACE ',I10,/,&
         ' IS NULL; COORDINATES: ',3E12.5)

#endif

!----
! End
!----

return
end subroutine
