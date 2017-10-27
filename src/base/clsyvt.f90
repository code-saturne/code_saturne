!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2017 EDF S.A.
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
!>
!> Please refer to the
!> <a href="../../theory.pdf#clsyvt"><b>clsyvt</b></a> section of the
!> theory guide for more informations.
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
!>                               - 6 rough wall and
!>                                 \f$ \vect{u} \cdot \vect{n} = 0 \f$
!>                               - 9 free inlet/outlet
!>                                 (input mass flux blocked to 0)
!> \param[in,out] rcodcl        boundary condition values:
!>                               - rcodcl(1) value of the dirichlet
!>                               - rcodcl(2) value of the exterior exchange
!>                                 coefficient (infinite if no exchange)
!>                               - rcodcl(3) value flux density
!>                                 (negative if gain) in w/m2 or roughness
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
!_______________________________________________________________________________


subroutine clsyvt &
 ( nscal  , icodcl ,                                     &
   rcodcl ,                                              &
   velipb , rijipb )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use optcal
use cstphy
use cstnum
use dimens, only: nvar
use pointe
use entsor
use albase
use field
use mesh
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer          nscal

integer          icodcl(nfabor,nvar)

double precision rcodcl(nfabor,nvar,3)
double precision velipb(nfabor,ndim), rijipb(nfabor,6)

! Local variables

integer          ifac, ii, isou
integer          iel, f_dim
integer          iscal , ivar

double precision rnx, rny, rnz, rxnn
double precision upx, upy, upz, usn
double precision tx, ty, tz, txn, t2x, t2y, t2z
double precision clsyme
double precision eloglo(3,3), alpha(6,6)
double precision srfbnf, rcodcn, hint, visclc, visctc, distbf
double precision cpp
double precision vismsh(3), hintv(3)
double precision visci(3,3), fikis, viscis, distfi
double precision fcoefa(6), fcoefb(6), fcofaf(6), fcofbf(6), fcofad(6), fcofbd(6)

double precision, dimension(:,:), pointer :: coefau, cofafu, cfaale, claale
double precision, dimension(:,:), pointer :: visten
double precision, dimension(:,:,:), pointer :: coefbu, cofbfu, cfbale, clbale
double precision, dimension(:), pointer :: coefa_r11, coefaf_r11, coefad_r11
double precision, dimension(:), pointer :: coefb_r11, coefbf_r11, coefbd_r11
double precision, dimension(:), pointer :: coefa_r22, coefaf_r22, coefad_r22
double precision, dimension(:), pointer :: coefb_r22, coefbf_r22, coefbd_r22
double precision, dimension(:), pointer :: coefa_r33, coefaf_r33, coefad_r33
double precision, dimension(:), pointer :: coefb_r33, coefbf_r33, coefbd_r33
double precision, dimension(:), pointer :: coefa_r12, coefaf_r12, coefad_r12
double precision, dimension(:), pointer :: coefb_r12, coefbf_r12, coefbd_r12
double precision, dimension(:), pointer :: coefa_r13, coefaf_r13, coefad_r13
double precision, dimension(:), pointer :: coefb_r13, coefbf_r13, coefbd_r13
double precision, dimension(:), pointer :: coefa_r23, coefaf_r23, coefad_r23
double precision, dimension(:), pointer :: coefb_r23, coefbf_r23, coefbd_r23
double precision, dimension(:,:), pointer :: coefa_rij, coefaf_rij, coefad_rij
double precision, dimension(:,:,:), pointer :: coefb_rij, coefbf_rij, coefbd_rij

double precision, dimension(:), pointer :: viscl, visct
double precision, dimension(:), pointer :: cpro_visma_s
double precision, dimension(:,:), pointer :: cpro_visma_v

type(var_cal_opt) :: vcopt_rij

!===============================================================================

!===============================================================================
! 1. Initializations
!===============================================================================

cpp = 0.d0

if (itytur.eq.3 .and. idirsm.eq.1) call field_get_val_v(ivsten, visten)

call field_get_val_s(iviscl, viscl)
call field_get_val_s(ivisct, visct)

! Boundary Conditions

call field_get_coefa_v(ivarfl(iu), coefau)
call field_get_coefb_v(ivarfl(iu), coefbu)
call field_get_coefaf_v(ivarfl(iu), cofafu)
call field_get_coefbf_v(ivarfl(iu), cofbfu)

if (itytur.eq.3) then
  call field_get_key_struct_var_cal_opt(ivarfl(ir11), vcopt_rij)
  if (irijco.eq.1) then
    call field_get_coefa_v(ivarfl(irij), coefa_rij)
    call field_get_coefb_v(ivarfl(irij), coefb_rij)
    call field_get_coefaf_v(ivarfl(irij), coefaf_rij)
    call field_get_coefbf_v(ivarfl(irij), coefbf_rij)
    call field_get_coefad_v(ivarfl(irij), coefad_rij)
    call field_get_coefbd_v(ivarfl(irij), coefbd_rij)

    coefb_r11 => null()
    coefaf_r11 => null()
    coefbf_r11 => null()
    coefad_r11 => null()
    coefbd_r11 => null()

    coefa_r22 => null()
    coefb_r22 => null()
    coefaf_r22 => null()
    coefbf_r22 => null()
    coefad_r22 => null()
    coefbd_r22 => null()

    coefa_r33 => null()
    coefb_r33 => null()
    coefaf_r33 => null()
    coefbf_r33 => null()
    coefad_r33 => null()
    coefbd_r33 => null()

    coefa_r12 => null()
    coefb_r12 => null()
    coefaf_r12 => null()
    coefbf_r12 => null()
    coefad_r12 => null()
    coefbd_r12 => null()

    coefa_r13 => null()
    coefb_r13 => null()
    coefaf_r13 => null()
    coefbf_r13 => null()
    coefad_r13 => null()
    coefbd_r13 => null()

    coefa_r23 => null()
    coefb_r23 => null()
    coefaf_r23 => null()
    coefbf_r23 => null()
    coefad_r23 => null()
    coefbd_r23 => null()


  else

    call field_get_coefa_s(ivarfl(ir11), coefa_r11)
    call field_get_coefb_s(ivarfl(ir11), coefb_r11)
    call field_get_coefaf_s(ivarfl(ir11), coefaf_r11)
    call field_get_coefbf_s(ivarfl(ir11), coefbf_r11)
    call field_get_coefad_s(ivarfl(ir11), coefad_r11)
    call field_get_coefbd_s(ivarfl(ir11), coefbd_r11)

    call field_get_coefa_s(ivarfl(ir22), coefa_r22)
    call field_get_coefb_s(ivarfl(ir22), coefb_r22)
    call field_get_coefaf_s(ivarfl(ir22), coefaf_r22)
    call field_get_coefbf_s(ivarfl(ir22), coefbf_r22)
    call field_get_coefad_s(ivarfl(ir22), coefad_r22)
    call field_get_coefbd_s(ivarfl(ir22), coefbd_r22)

    call field_get_coefa_s(ivarfl(ir33), coefa_r33)
    call field_get_coefb_s(ivarfl(ir33), coefb_r33)
    call field_get_coefaf_s(ivarfl(ir33), coefaf_r33)
    call field_get_coefbf_s(ivarfl(ir33), coefbf_r33)
    call field_get_coefad_s(ivarfl(ir33), coefad_r33)
    call field_get_coefbd_s(ivarfl(ir33), coefbd_r33)

    call field_get_coefa_s(ivarfl(ir12), coefa_r12)
    call field_get_coefb_s(ivarfl(ir12), coefb_r12)
    call field_get_coefaf_s(ivarfl(ir12), coefaf_r12)
    call field_get_coefbf_s(ivarfl(ir12), coefbf_r12)
    call field_get_coefad_s(ivarfl(ir12), coefad_r12)
    call field_get_coefbd_s(ivarfl(ir12), coefbd_r12)

    call field_get_coefa_s(ivarfl(ir13), coefa_r13)
    call field_get_coefb_s(ivarfl(ir13), coefb_r13)
    call field_get_coefaf_s(ivarfl(ir13), coefaf_r13)
    call field_get_coefbf_s(ivarfl(ir13), coefbf_r13)
    call field_get_coefad_s(ivarfl(ir13), coefad_r13)
    call field_get_coefbd_s(ivarfl(ir13), coefbd_r13)

    call field_get_coefa_s(ivarfl(ir23), coefa_r23)
    call field_get_coefb_s(ivarfl(ir23), coefb_r23)
    call field_get_coefaf_s(ivarfl(ir23), coefaf_r23)
    call field_get_coefbf_s(ivarfl(ir23), coefbf_r23)
    call field_get_coefad_s(ivarfl(ir23), coefad_r23)
    call field_get_coefbd_s(ivarfl(ir23), coefbd_r23)

    coefa_rij => null()
    coefb_rij => null()
    coefaf_rij => null()
    coefbf_rij => null()
    coefad_rij => null()
    coefbd_rij => null()
  endif
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
    visclc = viscl(iel)
    visctc = visct(iel)

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

      ! Symmetric tensor diffusivity (Daly Harlow -- GGDH)
      if (iand(vcopt_rij%idften, ANISOTROPIC_RIGHT_DIFFUSION).ne.0) then

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

      ! Scalar diffusivity
      else
        hint = (visclc+visctc*csrij/cmu)/distbf
      endif

      do isou = 1, 6
        fcoefa(isou) = 0.d0
        fcofad(isou) = 0.d0
        fcoefb(isou) = 0.d0
        fcofbd(isou) = 0.d0
      enddo

      do isou = 1,6

        if (isou.eq.1) then
          ivar = ir11
        elseif (isou.eq.2) then
          ivar = ir22
        elseif (isou.eq.3) then
          ivar = ir33
        elseif (isou.eq.4) then
          ivar = ir12
        elseif (isou.eq.5) then
          ivar = ir23
        elseif (isou.eq.6) then
          ivar = ir13
        endif

        ! Partial (or total if coupled) implicitation
        if (irijco.eq.1) then
          coefa_rij(isou, ifac)  = 0.d0
          coefaf_rij(isou, ifac) = 0.d0
          coefad_rij(isou, ifac) = 0.d0
          do ii = 1, 6
            coefb_rij(isou,ii, ifac)  = alpha(ii, isou)
            if (ii.eq.isou) then
              coefbf_rij(isou,ii, ifac) = hint * (1.d0 - coefb_rij(isou,ii, ifac))
            else
              coefbf_rij(isou,ii, ifac) = - hint * coefb_rij(isou,ii, ifac)
            endif
            coefbd_rij(isou,ii, ifac) = coefb_rij(isou,ii, ifac)
          enddo

        else if (iclsyr.eq.1) then
          do ii = 1, 6
            if (ii.ne.isou) then
              fcoefa(isou) = fcoefa(isou) + alpha(isou,ii) * rijipb(ifac,ii)
            endif
          enddo
          fcoefb(isou) = alpha(isou,isou)
        else
          do ii = 1, 6
            fcoefa(isou) = fcoefa(isou) + alpha(isou,ii) * rijipb(ifac,ii)
          enddo
          fcoefb(isou) = 0.d0
        endif

        ! Boundary conditions for the momentum equation
        fcofad(isou) = fcoefa(isou)
        fcofbd(isou) = fcoefb(isou)

        ! Translate into Diffusive flux BCs
        fcofaf(isou) = -hint*fcoefa(isou)
        fcofbf(isou) = hint*(1.d0-fcoefb(isou))
      enddo
      if (irijco.ne.1) then
        do isou = 1, 6
          if (isou.eq.1) then
            coefa_r11(ifac) = fcoefa(isou)
            coefb_r11(ifac) = fcoefb(isou)
            coefaf_r11(ifac) = fcofaf(isou)
            coefbf_r11(ifac) = fcofbf(isou)
            coefad_r11(ifac) = fcofad(isou)
            coefbd_r11(ifac) = fcofbd(isou)
          else if (isou.eq.2) then
            coefa_r22(ifac) = fcoefa(isou)
            coefb_r22(ifac) = fcoefb(isou)
            coefaf_r22(ifac) = fcofaf(isou)
            coefbf_r22(ifac) = fcofbf(isou)
            coefad_r22(ifac) = fcofad(isou)
            coefbd_r22(ifac) = fcofbd(isou)
          else if (isou.eq.3) then
            coefa_r33(ifac) = fcoefa(isou)
            coefb_r33(ifac) = fcoefb(isou)
            coefaf_r33(ifac) = fcofaf(isou)
            coefbf_r33(ifac) = fcofbf(isou)
            coefad_r33(ifac) = fcofad(isou)
            coefbd_r33(ifac) = fcofbd(isou)
          else if (isou.eq.4) then
            coefa_r12(ifac) = fcoefa(isou)
            coefb_r12(ifac) = fcoefb(isou)
            coefaf_r12(ifac) = fcofaf(isou)
            coefbf_r12(ifac) = fcofbf(isou)
            coefad_r12(ifac) = fcofad(isou)
            coefbd_r12(ifac) = fcofbd(isou)
          else if (isou.eq.5) then
            coefa_r23(ifac) = fcoefa(isou)
            coefb_r23(ifac) = fcoefb(isou)
            coefaf_r23(ifac) = fcofaf(isou)
            coefbf_r23(ifac) = fcofbf(isou)
            coefad_r23(ifac) = fcofad(isou)
            coefbd_r23(ifac) = fcofbd(isou)
          else if (isou.eq.6) then
            coefa_r13(ifac) = fcoefa(isou)
            coefb_r13(ifac) = fcoefb(isou)
            coefaf_r13(ifac) = fcofaf(isou)
            coefbf_r13(ifac) = fcofbf(isou)
            coefad_r13(ifac) = fcofad(isou)
            coefbd_r13(ifac) = fcofbd(isou)
          endif
        enddo
      endif

    endif

  endif
! --- Test sur la presence d'une condition de symetrie vitesse : fin

enddo
! ---  End of loop over boundary faces

!===========================================================================
! 3.bis Boundary conditions on transported vectors
!===========================================================================

do iscal = 1, nscal
  ! u'T'
  if (ityturt(iscal).eq.3) then
    call clsyvt_scalar(iscal, icodcl)
  endif

  ! additional transported vectors
  call field_get_dim(ivarfl(isca(iscal)), f_dim)
  if (f_dim.gt.1) then
    call clsyvt_vector(iscal, icodcl)
  endif
enddo

!===============================================================================
! 4. Symmetry boundary conditions for mesh velocity (ALE module)
!===============================================================================

if (iale.eq.1) then

  call field_get_coefa_v(ivarfl(iuma), claale)
  call field_get_coefb_v(ivarfl(iuma), clbale)
  call field_get_coefaf_v(ivarfl(iuma), cfaale)
  call field_get_coefbf_v(ivarfl(iuma), cfbale)

  if (iortvm.eq.0) then
    call field_get_val_s(ivisma, cpro_visma_s)
  else
    call field_get_val_v(ivisma, cpro_visma_v)
  endif

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
      if (iortvm.eq.0) then
        vismsh(1) = cpro_visma_s(iel)
        vismsh(2) = cpro_visma_s(iel)
        vismsh(3) = cpro_visma_s(iel)
      else
        vismsh(1) = cpro_visma_v(1,iel)
        vismsh(2) = cpro_visma_v(2,iel)
        vismsh(3) = cpro_visma_v(3,iel)
      endif

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

!===============================================================================
! Local functions
!===============================================================================

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     iscal         scalar id
!> \param[in,out] icodcl        face boundary condition code:
!>                               - 1 Dirichlet
!>                               - 3 Neumann
!>                               - 4 sliding and
!>                                 \f$ \vect{u} \cdot \vect{n} = 0 \f$
!>                               - 5 smooth wall and
!>                                 \f$ \vect{u} \cdot \vect{n} = 0 \f$
!>                               - 6 rough wall and
!>                                 \f$ \vect{u} \cdot \vect{n} = 0 \f$
!>                               - 9 free inlet/outlet
!>                                 (input mass flux blocked to 0)
!_______________________________________________________________________________

subroutine clsyvt_scalar &
 ( iscal  , icodcl )

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use optcal
use cstphy
use cstnum
use dimens, only: nvar
use pointe
use entsor
use albase
use parall
use ppppar
use ppthch
use ppincl
use cplsat
use mesh
use field

!===============================================================================

implicit none

! Arguments

integer          iscal

integer          icodcl(nfabor,nvar)

! Local variables

integer          ivar, f_id
integer          ifac, iel, isou, jsou
integer          ifcvsl

double precision cpp, rkl, visclc
double precision distbf, srfbnf
double precision rnx, rny, rnz
double precision hintt(6)
double precision hint, qimp

character(len=80) :: fname

double precision, dimension(:), pointer :: val_s, crom
double precision, dimension(:), pointer :: coefap, coefbp, cofafp, cofbfp
double precision, dimension(:,:), pointer :: coefaut, cofafut, cofarut, visten
double precision, dimension(:,:,:), pointer :: coefbut, cofbfut, cofbrut
double precision, dimension(:), pointer :: a_al, b_al, af_al, bf_al

double precision, dimension(:), pointer :: viscl, viscls, cpro_cp

!===============================================================================

if (ityturt(iscal).ne.3) return

ivar = isca(iscal)
f_id = ivarfl(ivar)

call field_get_val_s(ivarfl(ivar), val_s)
call field_get_val_v(ivsten, visten)
call field_get_val_s(iviscl, viscl)

call field_get_coefa_s(f_id, coefap)
call field_get_coefb_s(f_id, coefbp)
call field_get_coefaf_s(f_id, cofafp)
call field_get_coefbf_s(f_id, cofbfp)

call field_get_val_s(icrom, crom)

if (icp.ge.0) then
  call field_get_val_s(icp, cpro_cp)
endif

call field_get_key_int (f_id, kivisl, ifcvsl)
if (ifcvsl .ge. 0) then
  call field_get_val_s(ifcvsl, viscls)
endif

! Turbulent diffusive flux of the scalar T
! (blending factor so that the component v'T' have only
!  mu_T/(mu+mu_T)* Phi_T)

! Name of the scalar ivar
call field_get_name(ivarfl(ivar), fname)

! Index of the corresponding turbulent flux
call field_get_id(trim(fname)//'_turbulent_flux', f_id)

call field_get_coefa_v(f_id,coefaut)
call field_get_coefb_v(f_id,coefbut)
call field_get_coefaf_v(f_id,cofafut)
call field_get_coefbf_v(f_id,cofbfut)
call field_get_coefad_v(f_id,cofarut)
call field_get_coefbd_v(f_id,cofbrut)

! EB-GGDH/AFM/DFM alpha boundary conditions
if (iturt(iscal).eq.11 .or. iturt(iscal).eq.21 .or. iturt(iscal).eq.31) then

  ! Name of the scalar ivar
  call field_get_name(ivarfl(ivar), fname)

  ! Index of the corresponding turbulent flux
  call field_get_id(trim(fname)//'_alpha', f_id)

  call field_get_coefa_s (f_id, a_al)
  call field_get_coefb_s (f_id, b_al)
  call field_get_coefaf_s(f_id, af_al)
  call field_get_coefbf_s(f_id, bf_al)
endif

! --- Loop on boundary faces
do ifac = 1, nfabor

  ! Test on symmetry boundary condition: start
  if (icodcl(ifac,iu).eq.4) then

    iel = ifabor(ifac)
    ! --- Physical Propreties
    visclc = viscl(iel)
    cpp = 1.d0
    if (iscacp(iscal).eq.1) then
      if (icp.ge.0) then
        cpp = cpro_cp(iel)
      else
        cpp = cp0
      endif
    endif

    ! --- Geometrical quantities
    distbf = distb(ifac)
    srfbnf = surfbn(ifac)

    rnx = surfbo(1,ifac)/srfbnf
    rny = surfbo(2,ifac)/srfbnf
    rnz = surfbo(3,ifac)/srfbnf

    if (ifcvsl .lt. 0) then
      rkl = visls0(iscal)/cpp
    else
      rkl = viscls(iel)/cpp
    endif

    !FIXME for EB DFM
    do isou = 1, 6
      if (isou.le.3) then
        hintt(isou) = (0.5d0*(visclc + rkl)    &
          + ctheta(iscal)*visten(isou,iel)/csrij) / distbf
      else
        hintt(isou) = ctheta(iscal)*visten(isou,iel) / csrij / distbf
      endif
    enddo

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
    cofbfut(2,1,ifac) = hintt(4)*rnx**2  + hintt(2)*rny*rnx + hintt(5)*rnx*rnz
    cofbfut(1,3,ifac) = hintt(1)*rnx*rnz + hintt(4)*rny*rnz + hintt(6)*rnz**2
    cofbfut(3,1,ifac) = hintt(6)*rnx**2  + hintt(5)*rny*rnx + hintt(3)*rnz*rnx
    cofbfut(2,3,ifac) = hintt(4)*rnx*rnz + hintt(2)*rny*rnz + hintt(5)*rnz**2
    cofbfut(3,2,ifac) = hintt(6)*rnx*rny + hintt(5)*rny**2  + hintt(3)*rnz*rny

    ! Boundary conditions for thermal transport equation
    do isou = 1, 3
      cofarut(isou,ifac) = coefaut(isou,ifac)
      do jsou =1, 3
        cofbrut(isou,jsou,ifac) = coefbut(isou,jsou,ifac)
      enddo
    enddo

    ! EB-GGDH/AFM/DFM alpha boundary conditions
    if (iturt(iscal).eq.11 .or. iturt(iscal).eq.21 .or. iturt(iscal).eq.31) then

      ! Dirichlet Boundary Condition
      !-----------------------------

      qimp = 0.d0

      hint = 1.d0/distbf

      call set_neumann_scalar &
        !====================
      ( a_al(ifac), af_al(ifac),             &
        b_al(ifac), bf_al(ifac),             &
        qimp      , hint       )

    endif


  endif
  ! Test on velocity symmetry condition: end

enddo

!----
! End
!----

return
end subroutine

!===============================================================================
! Local functions
!===============================================================================

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     iscal         scalar id
!> \param[in,out] icodcl        face boundary condition code:
!>                               - 1 Dirichlet
!>                               - 3 Neumann
!>                               - 4 sliding and
!>                                 \f$ \vect{u} \cdot \vect{n} = 0 \f$
!>                               - 5 smooth wall and
!>                                 \f$ \vect{u} \cdot \vect{n} = 0 \f$
!>                               - 6 rough wall and
!>                                 \f$ \vect{u} \cdot \vect{n} = 0 \f$
!>                               - 9 free inlet/outlet
!>                                 (input mass flux blocked to 0)
!_______________________________________________________________________________

subroutine clsyvt_vector &
 ( iscal  , icodcl )

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use optcal
use cstphy
use cstnum
use dimens, only: nvar
use pointe
use entsor
use albase
use parall
use ppppar
use ppthch
use ppincl
use cplsat
use mesh
use field
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer          iscal

integer          icodcl(nfabor,nvar)

! Local variables

integer          ivar, f_id
integer          ifac, iel
integer          ifcvsl

double precision rkl, visclc
double precision distbf, srfbnf
double precision rnx, rny, rnz, temp
double precision hintt(6)
double precision turb_schmidt

double precision, dimension(:), pointer :: crom
double precision, dimension(:,:), pointer :: coefav, cofafv
double precision, dimension(:,:,:), pointer :: coefbv, cofbfv
double precision, dimension(:,:), pointer :: visten

double precision, dimension(:), pointer :: viscl, visct, viscls

type(var_cal_opt) :: vcopt

!===============================================================================

ivar = isca(iscal)
f_id = ivarfl(ivar)

call field_get_key_struct_var_cal_opt(f_id, vcopt)

if (iand(vcopt%idften, ANISOTROPIC_DIFFUSION).ne.0) then
  if (iturb.ne.32) then
    call field_get_val_v(ivsten, visten)
  else ! EBRSM and (GGDH or AFM)
    call field_get_val_v(ivstes, visten)
  endif
endif

call field_get_val_s(iviscl, viscl)
call field_get_val_s(ivisct, visct)

call field_get_coefa_v(f_id, coefav)
call field_get_coefb_v(f_id, coefbv)
call field_get_coefaf_v(f_id, cofafv)
call field_get_coefbf_v(f_id, cofbfv)

call field_get_val_s(icrom, crom)

call field_get_key_int (f_id, kivisl, ifcvsl)
if (ifcvsl .ge. 0) then
  call field_get_val_s(ifcvsl, viscls)
endif

! retrieve turbulent Schmidt value for current scalar
call field_get_key_double(ivarfl(ivar), ksigmas, turb_schmidt)

! --- Loop on boundary faces
do ifac = 1, nfabor

  ! Test on symmetry boundary condition: start
  if (icodcl(ifac,ivar).eq.4) then

    iel = ifabor(ifac)
    ! --- Physical Propreties
    visclc = viscl(iel)

    ! --- Geometrical quantities
    distbf = distb(ifac)
    srfbnf = surfbn(ifac)

    rnx = surfbo(1,ifac)/srfbnf
    rny = surfbo(2,ifac)/srfbnf
    rnz = surfbo(3,ifac)/srfbnf

    if (ifcvsl .lt. 0) then
      rkl = visls0(iscal)
    else
      rkl = viscls(iel)
    endif

    ! Isotropic diffusivity
    if (iand(vcopt%idften, ISOTROPIC_DIFFUSION).ne.0) then
      hintt(1) = (vcopt%idifft*max(visct(iel),zero)/turb_schmidt + rkl)/distbf
      hintt(2) = hintt(1)
      hintt(3) = hintt(1)
      hintt(4) = hintt(1)
      hintt(5) = hintt(1)
      hintt(6) = hintt(1)
    ! Symmetric tensor diffusivity
    elseif (iand(vcopt%idften, ANISOTROPIC_DIFFUSION).ne.0) then
      temp = vcopt%idifft*ctheta(iscal)/csrij
      hintt(1) = (temp*visten(1,iel) + rkl)/distbf
      hintt(2) = (temp*visten(2,iel) + rkl)/distbf
      hintt(3) = (temp*visten(3,iel) + rkl)/distbf
      hintt(4) =  temp*visten(4,iel)       /distbf
      hintt(5) =  temp*visten(5,iel)       /distbf
      hintt(6) =  temp*visten(6,iel)       /distbf
    endif

    ! Gradient BCs
    coefav(1,ifac) = 0.d0
    coefav(2,ifac) = 0.d0
    coefav(3,ifac) = 0.d0

    coefbv(1,1,ifac) = 1.d0-rnx**2
    coefbv(2,2,ifac) = 1.d0-rny**2
    coefbv(3,3,ifac) = 1.d0-rnz**2

    coefbv(1,2,ifac) = -rnx*rny
    coefbv(1,3,ifac) = -rnx*rnz
    coefbv(2,1,ifac) = -rny*rnx
    coefbv(2,3,ifac) = -rny*rnz
    coefbv(3,1,ifac) = -rnz*rnx
    coefbv(3,2,ifac) = -rnz*rny

    ! Flux BCs
    cofafv(1,ifac) = 0.d0
    cofafv(2,ifac) = 0.d0
    cofafv(3,ifac) = 0.d0

    cofbfv(1,1,ifac) = hintt(1)*rnx**2  + hintt(4)*rnx*rny + hintt(6)*rnx*rnz
    cofbfv(2,2,ifac) = hintt(4)*rnx*rny + hintt(2)*rny**2  + hintt(5)*rny*rnz
    cofbfv(3,3,ifac) = hintt(6)*rnx*rnz + hintt(5)*rny*rnz + hintt(3)*rnz**2

    cofbfv(1,2,ifac) = hintt(1)*rnx*rny + hintt(4)*rny**2  + hintt(6)*rny*rnz
    cofbfv(2,1,ifac) = hintt(1)*rnx*rny + hintt(4)*rny**2  + hintt(6)*rny*rnz
    cofbfv(1,3,ifac) = hintt(1)*rnx*rnz + hintt(4)*rny*rnz + hintt(6)*rnz**2
    cofbfv(3,1,ifac) = hintt(1)*rnx*rnz + hintt(4)*rny*rnz + hintt(6)*rnz**2
    cofbfv(2,3,ifac) = hintt(4)*rnx*rnz + hintt(2)*rny*rnz + hintt(5)*rnz**2
    cofbfv(3,2,ifac) = hintt(4)*rnx*rnz + hintt(2)*rny*rnz + hintt(5)*rnz**2

  endif
  ! Test on velocity symmetry condition: end

enddo

!----
! End
!----

return
end subroutine
