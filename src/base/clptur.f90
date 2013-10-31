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

!> \file clptur.f90
!>
!> \brief Boundary conditions for smooth walls (icodcl = 5).
!>
!> The wall functions may change the value of the diffusive flux.
!>
!> The values at a border face \f$ \fib \f$ stored in the face center
!> \f$ \centf \f$ of the variable \f$ P \f$ and its diffusive flux \f$ Q \f$
!> are written as:
!> \f[
!> P_\centf = A_P^g + B_P^g P_\centi
!> \f]
!> and
!> \f[
!> Q_\centf = A_P^f + B_P^f P_\centi
!> \f]
!> where \f$ P_\centi \f$ is the value of the variable \f$ P \f$ at the
!> neighbooring cell.
!>
!> Warning:
!>
!> - for a vector field such as the veclocity \f$ \vect{u} \f$ the boundary
!>   conditions may read:
!>   \f[
!>   \vect{u}_\centf = \vect{A}_u^g + \tens{B}_u^g \vect{u}_\centi
!>   \f]
!>   and
!>   \f[
!>   \vect{Q}_\centf = \vect{A}_u^f + \tens{B}_u^f \vect{u}_\centi
!>   \f]
!>   where \f$ \tens{B}_u^g \f$ and \f$ \tens{B}_u^f \f$ are 3x3 tensor matrix
!>   which coupled veclocity components next to a boundary.
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nscal         total number of scalars
!> \param[in]     isvhb         indicator to save exchange coeffient
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
!> \param[in]     dt            time step (per cell)
!> \param[in]     rtp           calculated variables at cell centers
!>                               (at current time step)
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
!> \param[out]    visvdr        viscosite dynamique ds les cellules
!>                               de bord apres amortisst de v driest
!> \param[out]    hbord         coefficients d'echange aux bords
!>
!> \param[in]     theipb        boundary temperature in \f$ \centip \f$
!>                               (more exaclty the energetic variable)
!_______________________________________________________________________________

subroutine clptur &
 ( nscal  , isvhb  , icodcl ,                                     &
   rtp    , propce , rcodcl ,                                     &
   velipb , rijipb , coefa  , coefb  , visvdr ,                   &
   hbord  , theipb )

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
use parall
use ppppar
use ppthch
use ppincl
use radiat
use cplsat
use mesh
use field
use lagran
use field
use turbomachinery

!===============================================================================

implicit none

! Arguments

integer          nscal, isvhb

integer          icodcl(nfabor,nvarcl)

double precision rtp(ncelet,*)
double precision propce(ncelet,*)
double precision rcodcl(nfabor,nvarcl,3)
double precision velipb(nfabor,ndim), rijipb(nfabor,6)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision visvdr(ncelet)
double precision hbord(nfabor),theipb(nfabor)

! Local variables

integer          ifac, iel, ivar, isou, jsou, ii, jj, kk, isvhbl
integer          iscal
integer          modntl
integer          iuntur
integer          nlogla, nsubla, iuiptn
integer          iclk  , iclep
integer          iclnu
integer          icl11 , icl22 , icl33 , icl12 , icl13 , icl23
integer          icl11r, icl22r, icl33r, icl12r, icl13r, icl23r
integer          iclphi, iclfb , iclal , iclomg
integer          iclvar, iclvrr
integer          iclkf , iclepf
integer          iclnuf
integer          icl11f, icl22f, icl33f, icl12f, icl13f, icl23f
integer          iclphf, iclfbf, iclalf, iclomf
integer          iclvaf
integer          ipcvis, ipcvst, ipccp , ipccv
integer          ipcvsl, itplus, itstar
integer          f_id
integer          iwallf

double precision rnx, rny, rnz, rxnn
double precision tx, ty, tz, txn, txn0, t2x, t2y, t2z
double precision utau, upx, upy, upz, usn
double precision uiptn, uiptmn, uiptmx
double precision uetmax, uetmin, ukmax, ukmin, yplumx, yplumn
double precision tetmax, tetmin, tplumx, tplumn
double precision uk, uet, nusury, yplus, dplus, tet, phit
double precision sqrcmu, clsyme, ek
double precision xnuii, xnuit, xmutlm
double precision rcprod
double precision cpp, rkl,  prdtl
double precision hflui, hint, hext, pimp, heq, qimp
double precision und0, deuxd0
double precision eloglo(3,3), alpha(6,6)
double precision rcodcx, rcodcy, rcodcz, rcodcn
double precision visclc, visctc, romc  , distbf, srfbnf, cpscv
double precision cofimp, ypup, yptp, ypth
double precision bldr12
double precision xkip
double precision tplus
double precision rinfiv(3), pimpv(3), hintt(6), pfac
double precision visci(3,3), fikis, viscis, distfi
double precision temp

character*80     fname

double precision, dimension(:), pointer :: tplusp, tstarp
double precision, dimension(:,:), pointer :: coefaut, cofafut, cofarut
double precision, dimension(:,:,:), pointer :: coefbut, cofbfut, cofbrut
double precision, dimension(:), pointer :: crom
double precision, dimension(:), pointer :: bfconv, bhconv
integer          ntlast , iaff
data             ntlast , iaff /-1 , 0/
save             ntlast , iaff

!===============================================================================

!===============================================================================
! 1. Initializations
!===============================================================================

! Initialize variables to avoid compiler warnings
iclep = 0
iclk = 0
iclomg = 0
iclfb = 0
iclphi = 0
iclnu = 0
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
ipccv = 0
iclal = 0
iclalf = 0
iclepf = 0
iclfbf = 0
iclkf = 0
iclnuf = 0
iclomf = 0
iclphf = 0

cofimp  = 0.d0
ek = 0.d0
rcprod = 0.d0
uiptn = 0.d0
ypth = 0.d0

rinfiv(1) = rinfin
rinfiv(2) = rinfin
rinfiv(3) = rinfin

! --- Constants
uet = 1.d0
utau = 1.d0
sqrcmu = sqrt(cmu)

und0   = 1.d0
deuxd0 = 2.d0

! --- Gradient Boundary Conditions
if (itytur.eq.2) then
  iclk   = iclrtp(ik ,icoef)
  iclep  = iclrtp(iep,icoef)
elseif (itytur.eq.3) then
  icl11  = iclrtp(ir11,icoef)
  icl22  = iclrtp(ir22,icoef)
  icl33  = iclrtp(ir33,icoef)
  icl12  = iclrtp(ir12,icoef)
  icl13  = iclrtp(ir13,icoef)
  icl23  = iclrtp(ir23,icoef)
  ! Boundary conditions for the momentum equation
  icl11r = iclrtp(ir11,icoefr)
  icl22r = iclrtp(ir22,icoefr)
  icl33r = iclrtp(ir33,icoefr)
  icl12r = iclrtp(ir12,icoefr)
  icl13r = iclrtp(ir13,icoefr)
  icl23r = iclrtp(ir23,icoefr)
  iclep  = iclrtp(iep,icoef)
  if (iturb.eq.32) iclal = iclrtp(ial,icoef)
elseif (itytur.eq.5) then
  iclk   = iclrtp(ik ,icoef)
  iclep  = iclrtp(iep,icoef)
  iclphi = iclrtp(iphi,icoef)
  if (iturb.eq.50) then
    iclfb  = iclrtp(ifb,icoef)
  elseif (iturb.eq.51) then
    iclal  = iclrtp(ial,icoef)
  endif
elseif (iturb.eq.60) then
  iclk   = iclrtp(ik ,icoef)
  iclomg = iclrtp(iomg,icoef)
elseif (iturb.eq.70) then
  iclnu  = iclrtp(inusa,icoef)
endif

! --- Flux Boundary Conditions
if (itytur.eq.2) then
  iclkf  = iclrtp(ik ,icoeff)
  iclepf = iclrtp(iep,icoeff)
elseif (itytur.eq.3) then
  icl11f = iclrtp(ir11,icoeff)
  icl22f = iclrtp(ir22,icoeff)
  icl33f = iclrtp(ir33,icoeff)
  icl12f = iclrtp(ir12,icoeff)
  icl13f = iclrtp(ir13,icoeff)
  icl23f = iclrtp(ir23,icoeff)
  iclepf = iclrtp(iep,icoeff)
  if (iturb.eq.32) iclalf = iclrtp(ial,icoeff)
elseif (itytur.eq.5) then
  iclkf  = iclrtp(ik ,icoeff)
  iclepf = iclrtp(iep,icoeff)
  iclphf = iclrtp(iphi,icoeff)
  if (iturb.eq.50) then
    iclfbf = iclrtp(ifb,icoeff)
  elseif (iturb.eq.51) then
    iclalf = iclrtp(ial,icoeff)
  endif
elseif (iturb.eq.60) then
  iclkf  = iclrtp(ik ,icoeff)
  iclomf = iclrtp(iomg,icoeff)
elseif (iturb.eq.70) then
  iclnuf = iclrtp(inusa,icoeff)
endif

! --- Physical quantities
call field_get_val_s(icrom, crom)
ipcvis = ipproc(iviscl)
ipcvst = ipproc(ivisct)
if (icp.gt.0) then
  ipccp  = ipproc(icp)
else
  ipccp = 0
endif

! --- Compressible
if (ippmod(icompf) .ge. 0) then
  if (icv.gt.0) then
    ipccv  = ipproc(icv)
  else
    ipccv = 0
  endif
endif

! min. and max. of wall tangential velocity
uiptmx = -grand
uiptmn =  grand

! min. and max. of wall friction velocity
uetmax = -grand
uetmin =  grand
ukmax  = -grand
ukmin  =  grand

! min. and max. of y+
yplumx = -grand
yplumn =  grand

! min. and max. of wall friction of the thermal scalar
tetmax = -grand
tetmin =  grand

! min. and max. of T+
tplumx = -grand
tplumn =  grand

! Counters (turbulent, laminar, reversal, scale correction)
nlogla = 0
nsubla = 0
iuiptn = 0


! With v2f type model, (phi-fbar et BL-v2/k) u=0 is set directly, so
! uiptmx and uiptmn are necessarily 0
if (itytur.eq.5) then
  uiptmx = 0.d0
  uiptmn = 0.d0
endif

! pointers to T+ and T* if saved

tplusp => null()
tstarp => null()

call field_get_id_try('tplus', itplus)
if (itplus.ge.0) then
  call field_get_val_s (itplus, tplusp)
endif

call field_get_id_try('tstar', itstar)
if (itstar.ge.0) then
  call field_get_val_s (itstar, tstarp)
endif

! Pointers to specific fields
if (ifconv.ge.0) call field_get_val_s(ifconv, bfconv)
if (ihconv.ge.0) call field_get_val_s(ihconv, bhconv)

! --- Loop on boundary faces
do ifac = 1, nfabor

  ! Test on the presence of a smooth wall condition (start)
  if (icodcl(ifac,iu).eq.5) then

    iel = ifabor(ifac)

    ! Physical properties
    visclc = propce(iel,ipcvis)
    visctc = propce(iel,ipcvst)
    romc   = crom(iel)

    ! Geometric quantities
    distbf = distb(ifac)
    srfbnf = surfbn(ifac)

    !===========================================================================
    ! 1. Local framework
    !===========================================================================

    ! Unit normal

    rnx = surfbo(1,ifac)/srfbnf
    rny = surfbo(2,ifac)/srfbnf
    rnz = surfbo(3,ifac)/srfbnf

    ! Handle displacement velocity

    rcodcx = rcodcl(ifac,iu,1)
    rcodcy = rcodcl(ifac,iv,1)
    rcodcz = rcodcl(ifac,iw,1)

    ! If we are not using ALE, force the displacement velocity for the face
    !  to be tangential (and update rcodcl for possible use)
    ! In frozen rotor (iturbo = 1), the velocity is neither tangential to the
    !  wall (absolute velocity solved in a relative frame of reference)
    if (iale.eq.0.and.imobil.eq.0.and.iturbo.eq.0) then
      rcodcn = rcodcx*rnx+rcodcy*rny+rcodcz*rnz
      rcodcx = rcodcx -rcodcn*rnx
      rcodcy = rcodcy -rcodcn*rny
      rcodcz = rcodcz -rcodcn*rnz
      rcodcl(ifac,iu,1) = rcodcx
      rcodcl(ifac,iv,1) = rcodcy
      rcodcl(ifac,iw,1) = rcodcz
    endif

    ! Relative tangential velocity

    upx = velipb(ifac,1) - rcodcx
    upy = velipb(ifac,2) - rcodcy
    upz = velipb(ifac,3) - rcodcz

    usn = upx*rnx+upy*rny+upz*rnz
    tx  = upx -usn*rnx
    ty  = upy -usn*rny
    tz  = upz -usn*rnz
    txn = sqrt(tx**2 +ty**2 +tz**2)
    utau= txn

    ! Unit tangent

    if (txn.ge.epzero) then

      txn0 = 1.d0

      tx  = tx/txn
      ty  = ty/txn
      tz  = tz/txn

    elseif (itytur.eq.3) then

      ! If the velocity is zero, vector T is normal and random;
      !   we need it for the reference change for Rij, and we cancel the velocity.

      txn0 = 0.d0

      if (abs(rny).ge.epzero.or.abs(rnz).ge.epzero)then
        rxnn = sqrt(rny**2+rnz**2)
        tx  =  0.d0
        ty  =  rnz/rxnn
        tz  = -rny/rxnn
      elseif (abs(rnx).ge.epzero.or.abs(rnz).ge.epzero)then
        rxnn = sqrt(rnx**2+rnz**2)
        tx  =  rnz/rxnn
        ty  =  0.d0
        tz  = -rnx/rxnn
      else
        write(nfecra,1000)ifac,rnx,rny,rnz
        call csexit (1)
      endif

    else

      ! If the velocity is zero, and we are not using Reynolds Stresses,
      !  Tx, Ty, Tz is not used (we cancel the velocity), so we assign any
      !  value (zero for example)

      txn0 = 0.d0

      tx  = 0.d0
      ty  = 0.d0
      tz  = 0.d0

    endif

    ! Complete if necessary for Rij-Epsilon

    if (itytur.eq.3) then

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

      clsyme=0.d0
      call clca66 (clsyme , eloglo , alpha)
      !==========

    endif

    !===========================================================================
    ! 2. Friction velocities
    !===========================================================================

    ! ---> Compute Uet depending if we are in the log zone or not
    !      in 1 or 2 velocity scales
    !      and uk based on ek

    nusury = visclc/(distbf*romc)
    xnuii = visclc/romc
    xnuit = visctc/romc

    if (itytur.eq.2 .or. itytur.eq.5 .or. iturb.eq.60) then
      ek = rtp(iel,ik)
    else if (itytur.eq.3) then
      ek = 0.5d0*(rtp(iel,ir11)+rtp(iel,ir22)+rtp(iel,ir33))
    endif


    ! With power law (Werner & Wengle)
    if (ideuch.eq.0.and.ilogpo.eq.0) then
      iwallf = 3
    else
      iwallf = ideuch
    endif

    ! No wall functions for low Reynolds models
    if (iturb.eq.0.or.itytur.eq.5.or.iturb.eq.32) then
      iwallf = 4
    endif

    call wallfunctions &
  ( iwallf, ifac  , xkappa, cstlog, cmu025, ypluli,        &
    apow  , bpow  , dpow  ,                                &
    xnuii , xnuit , utau  , distbf, ek    ,                &
    iuntur, nsubla, nlogla,                                &
    uet   , uk    , yplus , ypup  , cofimp, dplus )

    uetmax = max(uet,uetmax)
    uetmin = min(uet,uetmin)
    ukmax  = max(uk,ukmax)
    ukmin  = min(uk,ukmin)
    yplumx = max(yplus-dplus,yplumx)
    yplumn = min(yplus-dplus,yplumn)

    ! Sauvegarde de la vitesse de frottement et de la viscosite turbulente
    ! apres amortissement de van Driest pour la LES
    ! On n'amortit pas mu_t une seconde fois si on l'a deja fait
    ! (car une cellule peut avoir plusieurs faces de paroi)
    ! ou
    ! Sauvegarde de la vitesse de frottement et distance a la paroi yplus
    ! si le modele de depot de particules est active.

    if (itytur.eq.4.and.idries.eq.1) then
      uetbor(ifac) = uet !TODO remove, this information is in cofaf cofbf
      if (visvdr(iel).lt.-900.d0) then
        propce(iel,ipcvst) = propce(iel,ipcvst)                   &
             *(1.d0-exp(-(yplus-dplus)/cdries))**2
        visvdr(iel) = propce(iel,ipcvst)
        visctc = propce(iel,ipcvst)
      endif
    else if (iilagr.gt.0.and.idepst.gt.0) then
      uetbor(ifac) = uet
    endif

    ! Save yplus if post-processed
    if (ipstdv(ipstyp).ne.0) then
      yplbr(ifac) = yplus-dplus
    endif

    !===========================================================================
    ! 3. Velocity boundary conditions
    !===========================================================================

    ! Deprecated power law (Werner & Wengle)
    ! If ilogpo=0, then ideuch=0
    if (ilogpo.eq.0) then
      uiptn  = utau + uet*apow*bpow*yplus**bpow*(2.d0**(bpow-1.d0)-2.d0)

    ! Dependant on the turbulence Model
    else

      ! uiptn respecte la production de k
      !  de facon conditionnelle --> Coef RCPROD

      ! --> k-epsilon and k-omega
      !--------------------------
      if (itytur.eq.2.or.iturb.eq.60) then

        xmutlm = xkappa*visclc*yplus

        ! If yplus=0, uiptn is set to 0 to avoid division by 0.
        ! By the way, in this case: iuntur=0
        if (yplus.gt.epzero) then !TODO use iuntur.eq.1
          rcprod = min(xkappa , max(und0,sqrt(xmutlm/visctc))/yplus)

          uiptn  = utau + distbf*uet*uk*romc/xkappa/visclc*(       &
               und0/(deuxd0*yplus-dplus) - deuxd0*rcprod )
        else
          uiptn = 0.d0
        endif

      ! --> No turbulence, mixing length or Rij-espilon
      !------------------------------------------------
      elseif (iturb.eq.0.or.iturb.eq.10.or.itytur.eq.3) then

        ! Dans le cadre de la ponderation elliptique, on ne doit pas
        ! tenir compte des lois de paroi. On fait donc un test sur le modele
        ! de turbulence :
        ! si on est en LRR ou SSG on laisse les lois de paroi, si on est en
        ! EBRSM, on impose l adherence.
        if (iturb.eq.32.or.iturb.eq.0) then
          uiptn = 0.d0
        else

          ! If yplus=0, uiptn is set to 0 to avoid division by 0.
          ! By the way, in this case: iuntur=0
          if (yplus.gt.epzero) then !FIXME use iuntur
            uiptn = utau - distbf*romc*uet*uk/xkappa/visclc                    &
                                 *(deuxd0/yplus - und0/(deuxd0*yplus-dplus))
          else
            uiptn = 0.d0
          endif

        endif

      ! --> LES and Spalart Allmaras
      !-----------------------------
      elseif (itytur.eq.4.or.iturb.eq.70) then

        uiptn  = utau - uet/xkappa*1.5d0

        ! If (mu+mut) becomes zero (dynamic models), an arbitrary value is set
        ! (nul flux) but without any problems because the flux
        ! is really zero at this face.
        if (visctc+visclc.le.0) then
          hflui = 0.d0!FIXME

        else

        endif

      ! --> v2f
      !--------
      elseif (itytur.eq.5) then

        ! Avec ces conditions, pas besoin de calculer uiptmx, uiptmn
        ! et iuiptn qui sont nuls (valeur d'initialisation)
        uiptn  = 0.d0

      endif
    endif

    ! Min and Max and counter of reversal layer
    uiptmn = min(uiptn*iuntur,uiptmn)
    uiptmx = max(uiptn*iuntur,uiptmx)
    if (uiptn*iuntur.lt.-epzero) iuiptn = iuiptn + 1

    ! To be coherent with a wall function, clip it to 0
    cofimp = max(cofimp, 0.d0)

    ! On implicite le terme (rho*uet*uk)
    hflui = visclc / distbf * ypup

    if (itytur.eq.3) then
      hint =  visclc          /distbf
    else
      hint = (visclc + visctc)/distbf
    endif

    ! Coupled solving of the velocity components

    ! Gradient boundary conditions
    !-----------------------------
    rcodcn = rcodcx*rnx+rcodcy*rny+rcodcz*rnz

    coefau(1,ifac) = (1.d0-cofimp)*(rcodcx - rcodcn*rnx) + rcodcn*rnx
    coefau(2,ifac) = (1.d0-cofimp)*(rcodcy - rcodcn*rny) + rcodcn*rny
    coefau(3,ifac) = (1.d0-cofimp)*(rcodcz - rcodcn*rnz) + rcodcn*rnz

    ! Projection in order to have the velocity parallel to the wall
    ! B = cofimp * ( IDENTITY - n x n )

    coefbu(1,1,ifac) = cofimp*(1.d0-rnx**2)
    coefbu(2,2,ifac) = cofimp*(1.d0-rny**2)
    coefbu(3,3,ifac) = cofimp*(1.d0-rnz**2)
    coefbu(1,2,ifac) = -cofimp*rnx*rny
    coefbu(1,3,ifac) = -cofimp*rnx*rnz
    coefbu(2,1,ifac) = -cofimp*rny*rnx
    coefbu(2,3,ifac) = -cofimp*rny*rnz
    coefbu(3,1,ifac) = -cofimp*rnz*rnx
    coefbu(3,2,ifac) = -cofimp*rnz*rny

    ! Flux boundary conditions
    !-------------------------

    cofafu(1,ifac)   = -hflui*(rcodcx - rcodcn*rnx) - hint*rcodcn*rnx
    cofafu(2,ifac)   = -hflui*(rcodcy - rcodcn*rny) - hint*rcodcn*rny
    cofafu(3,ifac)   = -hflui*(rcodcz - rcodcn*rnz) - hint*rcodcn*rnz

    ! Projection in order to have the shear stress parallel to the wall
    !  B = hflui*( IDENTITY - n x n )

    cofbfu(1,1,ifac) = hflui*(1.d0-rnx**2) + hint*rnx**2
    cofbfu(2,2,ifac) = hflui*(1.d0-rny**2) + hint*rny**2
    cofbfu(3,3,ifac) = hflui*(1.d0-rnz**2) + hint*rnz**2

    cofbfu(1,2,ifac) = (hint - hflui)*rnx*rny
    cofbfu(1,3,ifac) = (hint - hflui)*rnx*rnz
    cofbfu(2,1,ifac) = (hint - hflui)*rny*rnx
    cofbfu(2,3,ifac) = (hint - hflui)*rny*rnz
    cofbfu(3,1,ifac) = (hint - hflui)*rnz*rnx
    cofbfu(3,2,ifac) = (hint - hflui)*rnz*rny

    !===========================================================================
    ! 4. Boundary conditions on k and espilon
    !===========================================================================

    if (itytur.eq.2) then

      ! Dirichlet Boundary Condition on k
      !----------------------------------

      iclvar = iclk
      iclvaf = iclkf

      if (iuntur.eq.1) then
        pimp = uk**2/sqrcmu
      else
        pimp = 0.d0
      endif
      hint = (visclc+visctc/sigmak)/distbf

      call set_dirichlet_scalar &
           !====================
         ( coefa(ifac,iclvar), coefa(ifac,iclvaf),             &
           coefb(ifac,iclvar), coefb(ifac,iclvaf),             &
           pimp              , hint              , rinfin )


      ! Neumann Boundary Condition on epsilon
      !--------------------------------------

      iclvar = iclep
      iclvaf = iclepf

      hint = (visclc+visctc/sigmae)/distbf

      ! If yplus=0, uiptn is set to 0 to avoid division by 0.
      ! By the way, in this case: iuntur=0
      if (yplus.gt.epzero.and.iuntur.eq.1) then !FIXME use only iuntur
        pimp = distbf*4.d0*uk**5*romc**2/           &
            (xkappa*visclc**2*(yplus+dplus)**2)

        qimp = -pimp*hint !TODO transform it, it is only to be fully equivalent
      else
        qimp = 0.d0
      endif

      call set_neumann_scalar &
           !==================
         ( coefa(ifac,iclvar), coefa(ifac,iclvaf),             &
           coefb(ifac,iclvar), coefb(ifac,iclvaf),             &
           qimp              , hint )

    !===========================================================================
    ! 5. Boundary conditions on Rij-epsilon
    !===========================================================================

    elseif (itytur.eq.3) then

      ! ---> Tensor Rij (Partially implicited)

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
        coefb(ifac,iclvar) = 0.0d0
        coefa(ifac,iclvrr) = 0.0d0
        coefb(ifac,iclvrr) = 0.0d0

      enddo

      ! blending factor so that the component R(n,tau) have only
      ! -mu_T/(mu+mu_T)*uet*uk
      bldr12 = visctc/(visclc + visctc)

      do isou = 1,6

        if (isou.eq.1) then
          ivar = ir11
          iclvar = icl11
          iclvaf = icl11f
          iclvrr = icl11r
          jj = 1
          kk = 1
        else if (isou.eq.2) then
          ivar = ir22
          iclvar = icl22
          iclvaf = icl22f
          iclvrr = icl22r
          jj = 2
          kk = 2
        else if (isou.eq.3) then
          ivar = ir33
          iclvar = icl33
          iclvaf = icl33f
          iclvrr = icl33r
          jj = 3
          kk = 3
        else if (isou.eq.4) then
          ivar = ir12
          iclvar = icl12
          iclvaf = icl12f
          iclvrr = icl12r
          jj = 1
          kk = 2
        else if (isou.eq.5) then
          ivar = ir13
          iclvar = icl13
          iclvaf = icl13f
          iclvrr = icl13r
          jj = 1
          kk = 3
        else if (isou.eq.6) then
          ivar = ir23
          iclvar = icl23
          iclvaf = icl23f
          iclvrr = icl23r
          jj = 2
          kk = 3
        endif

        ! LRR and the Standard SGG.
        if ((iturb.eq.30).or.(iturb.eq.31).and.iuntur.eq.1) then

          if (iclptr.eq.1) then
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

          coefa(ifac,iclvar) = coefa(ifac,iclvar)                    &
                             - ( eloglo(jj,1)*eloglo(kk,2)           &
                               + eloglo(jj,2)*eloglo(kk,1)           &
                               )*bldr12*uet*uk

        ! In the viscous sublayer or for EBRSM: zero Reynolds' stresses
        else
          coefa(ifac,iclvar) = 0.d0
          coefa(ifac,iclvrr) = 0.d0
          coefb(ifac,iclvar) = 0.d0
          coefb(ifac,iclvrr) = 0.d0
        endif

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

      ! ---> Epsilon
      !      NB: no reconstruction, possibility of partial implicitation

      ivar = iep
      iclvar = iclep
      iclvaf = iclepf

      ! Symmetric tensor diffusivity (Daly Harlow -- GGDH)
      if (idften(ivar).eq.6) then

        visci(1,1) = visclc + visten(1,iel)/sigmae
        visci(2,2) = visclc + visten(2,iel)/sigmae
        visci(3,3) = visclc + visten(3,iel)/sigmae
        visci(1,2) =          visten(4,iel)/sigmae
        visci(2,1) =          visten(4,iel)/sigmae
        visci(2,3) =          visten(5,iel)/sigmae
        visci(3,2) =          visten(5,iel)/sigmae
        visci(1,3) =          visten(6,iel)/sigmae
        visci(3,1) =          visten(6,iel)/sigmae

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

      if ((iturb.eq.30).or.(iturb.eq.31)) then

        ! Si yplus=0, on met coefa a 0 directement pour eviter une division
        ! par 0.
        if (yplus.gt.epzero.and.iuntur.eq.1) then
          pimp = distbf*4.d0*uk**5*romc**2/           &
                (xkappa*visclc**2*(yplus+dplus)**2)
        else
          pimp = 0.d0
        endif

        ! Neumann Boundary Condition
        !---------------------------

        if (iclptr.eq.1) then !TODO not available for k-eps

          qimp = -pimp*hint !TODO transform it, it is only to be fully equivalent

           call set_neumann_scalar &
           !==================
         ( coefa(ifac,iclvar), coefa(ifac,iclvaf),             &
           coefb(ifac,iclvar), coefb(ifac,iclvaf),             &
           qimp              , hint )



        ! Dirichlet Boundary Condition
        !-----------------------------

        else

          pimp = pimp + rtp(iel,iep)

          call set_dirichlet_scalar &
               !====================
             ( coefa(ifac,iclvar), coefa(ifac,iclvaf),             &
               coefb(ifac,iclvar), coefb(ifac,iclvaf),             &
               pimp              , hint              , rinfin )



        endif
      elseif (iturb.eq.32) then
        ! Use k at I'
        xkip = 0.5d0*(rijipb(ifac,1)+rijipb(ifac,2)+rijipb(ifac,3))

        ! Dirichlet Boundary Condition
        !-----------------------------

        pimp = 2.d0*visclc*xkip/(distbf**2*romc)

        call set_dirichlet_scalar &
             !====================
           ( coefa(ifac,iclvar), coefa(ifac,iclvaf),             &
             coefb(ifac,iclvar), coefb(ifac,iclvaf),             &
             pimp              , hint              , rinfin )


        ! ---> Alpha
        iclvar = iclal
        iclvaf = iclalf

        ! Dirichlet Boundary Condition
        !-----------------------------

        pimp = 0.d0

        hint = 1.d0/distbf

        call set_dirichlet_scalar &
             !====================
           ( coefa(ifac,iclvar), coefa(ifac,iclvaf),             &
             coefb(ifac,iclvar), coefb(ifac,iclvaf),             &
             pimp              , hint              , rinfin )

      endif

    !===========================================================================
    ! 6a.Boundary conditions on k, epsilon, f_bar and phi in the phi_Fbar model
    !===========================================================================

    elseif (iturb.eq.50) then

      ! Dirichlet Boundary Condition on k
      !----------------------------------

      iclvar = iclk
      iclvaf = iclkf

      pimp = 0.d0
      hint = (visclc+visctc/sigmak)/distbf

      call set_dirichlet_scalar &
           !====================
         ( coefa(ifac,iclvar), coefa(ifac,iclvaf),             &
           coefb(ifac,iclvar), coefb(ifac,iclvaf),             &
           pimp              , hint              , rinfin )

      ! Dirichlet Boundary Condition on epsilon
      !----------------------------------------

      iclvar = iclep
      iclvaf = iclepf

      pimp = 2.0d0*visclc/romc*rtp(iel,ik)/distbf**2
      hint = (visclc+visctc/sigmae)/distbf

      call set_dirichlet_scalar &
           !====================
         ( coefa(ifac,iclvar), coefa(ifac,iclvaf),             &
           coefb(ifac,iclvar), coefb(ifac,iclvaf),             &
           pimp              , hint              , rinfin )

      ! Dirichlet Boundary Condition on Phi
      !------------------------------------

      iclvar = iclphi
      iclvaf = iclphf

      pimp = 0.d0
      hint = (visclc+visctc/sigmak)/distbf

      call set_dirichlet_scalar &
           !====================
         ( coefa(ifac,iclvar), coefa(ifac,iclvaf),             &
           coefb(ifac,iclvar), coefb(ifac,iclvaf),             &
           pimp              , hint              , rinfin )

      ! Dirichlet Boundary Condition on Fb
      !-----------------------------------

      iclvar = iclfb
      iclvaf = iclfbf

      pimp = 0.d0
      hint = 1.d0/distbf

      call set_dirichlet_scalar &
           !====================
         ( coefa(ifac,iclvar), coefa(ifac,iclvaf),             &
           coefb(ifac,iclvar), coefb(ifac,iclvaf),             &
           pimp              , hint              , rinfin )


    !===========================================================================
    ! 6b.Boundary conditions on k, epsilon, phi and alpha in the Bl-v2/k model
    !===========================================================================

    elseif (iturb.eq.51) then

      ! Dirichlet Boundary Condition on k
      !----------------------------------

      iclvar = iclk
      iclvaf = iclkf

      pimp = 0.d0
      hint = (visclc+visctc/sigmak)/distbf

      call set_dirichlet_scalar &
           !====================
         ( coefa(ifac,iclvar), coefa(ifac,iclvaf),             &
           coefb(ifac,iclvar), coefb(ifac,iclvaf),             &
           pimp              , hint              , rinfin )

      ! Dirichlet Boundary Condition on epsilon
      !----------------------------------------

      iclvar = iclep
      iclvaf = iclepf
      pimp = visclc/romc*rtp(iel,ik)/distbf**2
      hint = (visclc+visctc/sigmae)/distbf

      call set_dirichlet_scalar &
           !====================
         ( coefa(ifac,iclvar), coefa(ifac,iclvaf),             &
           coefb(ifac,iclvar), coefb(ifac,iclvaf),             &
           pimp              , hint              , rinfin )

      ! Dirichlet Boundary Condition on Phi
      !------------------------------------

      iclvar = iclphi
      iclvaf = iclphf

      pimp = 0.d0
      hint = (visclc+visctc/sigmak)/distbf

      call set_dirichlet_scalar &
           !====================
         ( coefa(ifac,iclvar), coefa(ifac,iclvaf),             &
           coefb(ifac,iclvar), coefb(ifac,iclvaf),             &
           pimp              , hint              , rinfin )

      ! Dirichlet Boundary Condition on alpha
      !--------------------------------------

      iclvar = iclal
      iclvaf = iclalf

      pimp = 0.d0
      hint = 1.d0/distbf

      call set_dirichlet_scalar &
           !====================
         ( coefa(ifac,iclvar), coefa(ifac,iclvaf),             &
           coefb(ifac,iclvar), coefb(ifac,iclvaf),             &
           pimp              , hint              , rinfin )


    !===========================================================================
    ! 7. Boundary conditions on k and omega
    !===========================================================================

    elseif (iturb.eq.60) then

      ! Dirichlet Boundary Condition on k
      !----------------------------------

      iclvar = iclk
      iclvaf = iclkf

      ! Si on est hors de la sous-couche visqueuse (reellement ou via les
      ! scalable wall functions)
      if (iuntur.eq.1) then
        pimp = uk**2/sqrcmu

      ! Si on est en sous-couche visqueuse
      else
        pimp = 0.d0
      endif

      !FIXME it is wrong because sigma is computed within the model
      ! see turbkw.f90
      hint = (visclc+visctc/ckwsk2)/distbf

      call set_dirichlet_scalar &
           !====================
         ( coefa(ifac,iclvar), coefa(ifac,iclvaf),             &
           coefb(ifac,iclvar), coefb(ifac,iclvaf),             &
           pimp              , hint              , rinfin )

      ! Neumann Boundary Condition on omega
      !------------------------------------

      iclvar = iclomg
      iclvaf = iclomf

      !FIXME it is wrong because sigma is computed within the model
      ! see turbkw.f90 (So the flux is not the one we impose!)
      hint = (visclc+visctc/ckwsw2)/distbf

      ! Si on est hors de la sous-couche visqueuse (reellement ou via les
      ! scalable wall functions)
      ! Comme iuntur=1 yplus est forcement >0
      if (iuntur.eq.1) then
        pimp = distbf*4.d0*uk**3*romc**2/           &
              (sqrcmu*xkappa*visclc**2*(yplus+dplus)**2)
        qimp = -pimp*hint !TODO transform it, it is only to be fully equivalent

      ! Si on est en sous-couche visqueuse
      else
        pimp = 120.d0*8.d0*visclc/(romc*ckwbt1*distbf**2)
        qimp = -pimp*hint !TODO transform it, it is only to be fully equivalent
      endif

      call set_neumann_scalar &
           !==================
         ( coefa(ifac,iclvar), coefa(ifac,iclvaf),             &
           coefb(ifac,iclvar), coefb(ifac,iclvaf),             &
           qimp              , hint )

    !===========================================================================
    ! 7.1 Boundary conditions on the Spalart Allmaras turbulence model
    !===========================================================================

    elseif (iturb.eq.70) then

      ! Dirichlet Boundary Condition on nusa
      !-------------------------------------

      iclvar = iclnu
      iclvaf = iclnuf

      pimp = 0.d0

      hint = visclc/distbf !FIXME nusa

      call set_dirichlet_scalar &
           !====================
         ( coefa(ifac,iclvar), coefa(ifac,iclvaf),             &
           coefb(ifac,iclvar), coefb(ifac,iclvaf),             &
           pimp              , hint              , rinfin )

    endif

    !===========================================================================
    ! 8. Boundary conditions on the other scalars
    !    (Specific treatment for the variances of the scalars next to walls:
    !     see condli)
    !===========================================================================

    if (nscal.ge.1) then

      do iscal = 1, nscal

        if (iscavr(iscal).le.0) then

          ivar = isca(iscal)
          iclvar = iclrtp(ivar,icoef)
          iclvaf = iclrtp(ivar,icoeff)

          isvhbl = 0
          if (iscal.eq.isvhb) then
            isvhbl = isvhb
          endif

          cpp = 1.d0
          if (iscacp(iscal).eq.1) then
            if (ipccp.gt.0) then
              cpp = propce(iel,ipccp)
            else
              cpp = cp0
            endif
          endif

          if (ivisls(iscal).gt.0) then
            ipcvsl = ipproc(ivisls(iscal))
          else
            ipcvsl = 0
          endif
          if (ipcvsl.le.0) then
            rkl = visls0(iscal)
            prdtl = cpp*visclc/rkl
          else
            rkl = propce(iel,ipcvsl)
            prdtl = cpp*visclc/rkl
          endif

          ! --> Compressible module:
          !  On suppose que le nombre de Pr doit etre
          !  defini de la meme facon que l'on resolve
          !  en enthalpie ou en energie, soit Mu*Cp/Lambda.
          !  Si l'on resout en energie, on a calcule ci-dessus
          !  Mu*Cv/Lambda.

          if (iscal.eq.iscalt .and. itherm.eq.3) then
            if (ipccp.gt.0) then
              prdtl = prdtl*propce(iel,ipccp )
            else
              prdtl = prdtl*cp0
            endif
            if (ipccv.gt.0) then
              prdtl = prdtl/propce(iel,ipccv )
            else
              prdtl = prdtl/cv0
            endif
          endif

          ! Scalar diffusivity
          if (idften(ivar).eq.1) then

            ! En compressible, pour l'energie LAMBDA/CV+CP/CV*(MUT/SIGMAS)
            if (ippmod(icompf) .ge. 0) then
              if (ipccp.gt.0) then
                cpscv = propce(iel,ipproc(icp))
              else
                cpscv = cp0
              endif
              if (ipccv.gt.0) then
                cpscv = cpscv/propce(iel,ipproc(icv))
              else
                cpscv = cpscv/cv0
              endif
              hint = (rkl+idifft(ivar)*cpp*cpscv*visctc/sigmas(iscal))/distbf
            else
              hint = (rkl+idifft(ivar)*cpp*visctc/sigmas(iscal))/distbf
            endif

          ! Symmetric tensor diffusivity (GGDH or AFM)
          elseif (idften(ivar).eq.6) then
            ! En compressible, pour l'energie LAMBDA/CV+CP/CV*(MUT/SIGMAS)
            if (ippmod(icompf) .ge. 0) then
              if (ipccp.gt.0) then
                cpscv = propce(iel,ipproc(icp))
              else
                cpscv = cp0
              endif
              if (ipccv.gt.0) then
                cpscv = cpscv/propce(iel,ipproc(icv))
              else
                cpscv = cpscv/cv0
              endif
              temp = idifft(ivar)*cpp*cpscv*ctheta(iscal)/csrij
            else
              temp = idifft(ivar)*cpp*ctheta(iscal)/csrij
            endif
            visci(1,1) = temp*visten(1,iel) + rkl
            visci(2,2) = temp*visten(2,iel) + rkl
            visci(3,3) = temp*visten(3,iel) + rkl
            visci(1,2) = temp*visten(4,iel)
            visci(2,1) = temp*visten(4,iel)
            visci(2,3) = temp*visten(5,iel)
            visci(3,2) = temp*visten(5,iel)
            visci(1,3) = temp*visten(6,iel)
            visci(3,1) = temp*visten(6,iel)

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

            ! Take I" so that I"F= eps*||FI||*Ki.n when I" is not in cell i
            ! NB: eps =1.d-1 must be consistent with vitens.f90
            fikis = max(fikis, 1.d-1*sqrt(viscis)*distfi)

            hint = viscis/surfbn(ifac)/fikis
          endif

          ! Dirichlet on the scalar, with wall function
          if (iturb.ne.0.and.icodcl(ifac,ivar).eq.5) then

            call hturbp (prdtl,sigmas(iscal),xkappa,yplus,dplus,hflui,ypth)

            ! Compute (y+-d+)/T+ *PrT
            yptp = hflui/prdtl
            ! Compute lambda/y * (y+-d+)/T+
            hflui = rkl/distbf *hflui

          ! Neumann on the scalar, with wall function (for post-processing)
          elseif (iturb.ne.0.and.icodcl(ifac,ivar).eq.3) then
            call hturbp (prdtl,sigmas(iscal),xkappa,yplus,dplus,hflui,ypth)
            ! y+/T+ *PrT
            yptp = hflui/prdtl
            hflui = hint

          else
            ! y+/T+ *PrT
            yptp = 1.d0/prdtl
            hflui = hint
          endif

          if (isvhbl.gt.0) hbord(ifac) = hflui

          ! ---> Dirichlet Boundary condition with a wall function correction
          !      with or without an additional exchange coefficient hext

          if (icodcl(ifac,ivar).eq.5) then
            hext = rcodcl(ifac,ivar,2)
            pimp = rcodcl(ifac,ivar,1)

            if (abs(hext).gt.rinfin*0.5d0) then
              heq = hflui
            else
              heq = hflui*hext/(hflui+hext)
            endif

            ! DFM: the gradient BCs are so that the production term
            !      of u'T' is correcty computed
            if (ityturt(iscal).ge.1) then
              ! In the log layer
              if (yplus.ge.ypth.and.iturb.ne.0) then
                xmutlm = xkappa*visclc*yplus
                rcprod = min(xkappa , max(und0,sqrt(xmutlm/visctc))/yplus)

                cofimp = 1.d0 - yptp*sigmas(iscal)/xkappa*                        &
                                (deuxd0*rcprod - und0/(deuxd0*yplus-dplus))

                ! the term (rho*tet*uk) is partially implicit
                pfac = (1.d0 -cofimp)*pimp

              ! In the viscous sub-layer
              else
                cofimp = 0.d0
                pfac = pimp
              endif
            else
              pfac = pimp * heq/hint
              cofimp = 1.d0 - heq/hint
            endif

            ! Gradient BCs
            coefa(ifac,iclvar) = pfac
            coefb(ifac,iclvar) = cofimp

            ! Flux BCs
            coefa(ifac,iclvaf) = -heq*pimp
            coefb(ifac,iclvaf) =  heq

            !--> Turbulent heat flux
            if (ityturt(iscal).eq.3) then

              ! Turbulent diffusive flux of the scalar T
              ! (blending factor so that the component v'T' have only
              !  mu_T/(mu+mu_T)* Phi_T)
              phit = (coefa(ifac,iclvaf) + coefb(ifac,iclvaf)*rtp(iel,ivar))

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

              hintt(1) = 0.5d0*(visclc+rkl)/distbf                        &
                       + visten(1,iel)*ctheta(iscal)/distbf/csrij
              hintt(2) = 0.5d0*(visclc+rkl)/distbf                        &
                       + visten(2,iel)*ctheta(iscal)/distbf/csrij
              hintt(3) = 0.5d0*(visclc+rkl)/distbf                        &
                       + visten(3,iel)*ctheta(iscal)/distbf/csrij
              hintt(4) = visten(4,iel)*ctheta(iscal)/distbf/csrij
              hintt(5) = visten(5,iel)*ctheta(iscal)/distbf/csrij
              hintt(6) = visten(6,iel)*ctheta(iscal)/distbf/csrij

              ! Dirichlet Boundary Condition
              !-----------------------------

              ! Add rho*uk*Tet to T'v' in High Reynolds
              if (yplus.ge.ypth) then
                do isou = 1, 3
                  pimpv(isou) = surfbo(isou,ifac)*phit/(surfbn(ifac)*cpp*romc)
                enddo
              else
                do isou = 1, 3
                  pimpv(isou) = 0.d0
                enddo
              endif

              call set_dirichlet_vector_ggdh &
                   !========================
                 ( coefaut(:,ifac)  , cofafut(:,ifac)  ,           &
                   coefbut(:,:,ifac), cofbfut(:,:,ifac),           &
                   pimpv            , hintt            , rinfiv )

              ! Boundary conditions used in the temperature equation
              do isou = 1, 3
                cofarut(isou,ifac) = 0.d0
                do jsou = 1, 3
                  cofbrut(isou,jsou,ifac) = 0.d0
                enddo
              enddo

            endif

            !--> Radiative module:

            ! On stocke le coefficient d'echange lambda/distance
            ! (ou son equivalent en turbulent) quelle que soit la
            ! variable thermique transportee (temperature ou enthalpie)
            ! car on l'utilise pour realiser des bilans aux parois qui
            ! sont faits en temperature (on cherche la temperature de
            ! paroi quelle que soit la variable thermique transportee pour
            ! ecrire des eps sigma T4.

!     donc :

!       lorsque la variable transportee est la temperature
!         ABS(ISCSTH(II)).EQ.1 : RA(IHCONV-1+IFAC) = HINT
!         puisque HINT = VISLS * CP / DISTBR
!                      = lambda/distance en W/(m2 K)

!       lorsque la variable transportee est l'enthalpie
!         ISCSTH(II).EQ.2 : RA(IHCONV-1+IFAC) = HINT*CPR
!         avec
!            IF (IPCCP.GT.0) THEN
!              CPR = PROPCE(IEL,IPCCP )
!            ELSE
!              CPR = CP0
!            ENDIF
!         puisque HINT = VISLS / DISTBR
!                      = lambda/(CP * distance)

!       lorsque la variable transportee est l'energie (compressible)
!         ISCSTH(II).EQ.3 :
!         on procede comme pour l'enthalpie avec CV au lieu de CP
!         (rq : il n'y a pas d'hypothese, sf en non orthogonal :
!               le flux est le bon et le coef d'echange aussi)

!      De meme dans condli.

!               Si on rayonne et que
!                  le scalaire est la variable energetique

            if (iirayo.ge.1 .and. iscal.eq.iscalt) then

              ! We compute the exchange coefficient in W/(m2 K)

              ! Enthalpy
              if (itherm.eq.2) then
                ! If Cp is variable
                if (ipccp.gt.0) then
                  bhconv(ifac) = hflui*propce(iel, ipccp)
                else
                  bhconv(ifac) = hflui*cp0
                endif

              ! Energie (compressible module)
              elseif (itherm.eq.3) then
                ! If Cv is variable
                if (ipccv.gt.0) then
                  bhconv(ifac) = hflui*propce(iel,ipccv)
                else
                  bhconv(ifac) = hflui*cv0
                endif

              ! Temperature
              elseif (itherm.eq.1) then
                bhconv(ifac) = hflui
              endif

              ! The outgoing flux is stored (Q = h(Ti'-Tp): negative if
              !  gain for the fluid) in W/m2
              bfconv(ifac) = coefa(ifac,iclvaf)              &
                                          + coefb(ifac,iclvaf)*theipb(ifac)
            endif

          endif ! End of icodcl.eq.5

          ! Save the value of T^star and T^+ for post-processing
          if (iscal.eq.iscalt) then

            ! Wall function
            if (icodcl(ifac,ivar).eq.5) then
              phit = coefa(ifac,iclvaf)+coefb(ifac,iclvaf)*theipb(ifac)

            ! Imposed flux with wall function for post-processing
            elseif (icodcl(ifac,ivar).eq.3) then
              phit = rcodcl(ifac,ivar,3)
            else
              phit = 0.d0
            endif

            tet = phit/(romc*cpp*max(uk,epzero))
            ! T+ = (T_I - T_w) / Tet
            tplus = (yplus-dplus)/yptp

            if (itplus .ge. 0) tplusp(ifac) = tplus
            if (itstar .ge. 0) tstarp(ifac) = tet

            tetmax = max(tet, tetmax)
            tetmin = min(tet, tetmin)
            tplumx = max(tplus,tplumx)
            tplumn = min(tplus,tplumn)
          endif

        endif

      enddo

    endif

  endif
  ! Test on the presence of a smooth wall (End)

enddo
! --- End of loop over faces

if (irangp.ge.0) then
  call parmin (uiptmn)
  call parmax (uiptmx)
  call parmin (uetmin)
  call parmax (uetmax)
  call parmin (ukmin)
  call parmax (ukmax)
  call parmin (yplumn)
  call parmax (yplumx)
  call parcpt (nlogla)
  call parcpt (nsubla)
  call parcpt (iuiptn)
  if (iscalt.gt.0) then
    call parmin (tetmin)
    call parmax (tetmax)
    call parmin (tplumn)
    call parmax (tplumx)
  endif
endif

!===============================================================================
! 9. Writings
!===============================================================================

!     Remarque : afin de ne pas surcharger les listings dans le cas ou
!       quelques yplus ne sont pas corrects, on ne produit le message
!       qu'aux deux premiers pas de temps ou le message apparait et
!       aux deux derniers pas de temps du calcul, ou si IWARNI est >= 2
!       On indique aussi le numero du dernier pas de temps auquel on
!       a rencontre des yplus hors bornes admissibles

if (iwarni(iu).ge.0) then
  if (ntlist.gt.0) then
    modntl = mod(ntcabs,ntlist)
  elseif (ntlist.eq.-1.and.ntcabs.eq.ntmabs) then
    modntl = 0
  else
    modntl = 1
  endif

  if ( (iturb.eq.0.and.nlogla.ne.0)                      .or.     &
       (itytur.eq.5.and.nlogla.ne.0)                     .or.     &
       ((itytur.eq.2.or.itytur.eq.3) .and. nsubla.gt.0)      )    &
       ntlast = ntcabs

  if ( (ntlast.eq.ntcabs.and.iaff.lt.2         ) .or.             &
       (ntlast.ge.0     .and.ntcabs.ge.ntmabs-1) .or.             &
       (ntlast.eq.ntcabs.and.iwarni(iu).ge.2) ) then
    iaff = iaff + 1

    if (iscalt.gt.0) then
      write(nfecra,2011) &
           uiptmn,uiptmx,uetmin,uetmax,ukmin,ukmax,yplumn,yplumx,      &
           tetmin, tetmax, tplumn, tplumx, iuiptn,nsubla,nsubla+nlogla
    else
      write(nfecra,2010) &
           uiptmn,uiptmx,uetmin,uetmax,ukmin,ukmax,yplumn,yplumx,   &
           iuiptn,nsubla,nsubla+nlogla
    endif

    if (iturb.eq. 0) write(nfecra,2020)  ntlast,ypluli
    if (itytur.eq.5) write(nfecra,2030)  ntlast,ypluli
    ! No warnings in EBRSM
    if (itytur.eq.2.or.iturb.eq.30.or.iturb.eq.31)                &
      write(nfecra,2040)  ntlast,ypluli
    if (iwarni(iu).lt.2.and.iturb.ne.32) then
      write(nfecra,2050)
    elseif (iturb.ne.32) then
      write(nfecra,2060)
    endif

  else if (modntl.eq.0 .or. iwarni(iu).ge.2) then

    if (iscalt.gt.0) then
      write(nfecra,2011) &
           uiptmn,uiptmx,uetmin,uetmax,ukmin,ukmax,yplumn,yplumx,      &
           tetmin, tetmax, tplumn, tplumx, iuiptn,nsubla,nsubla+nlogla
    else
      write(nfecra,2010) &
           uiptmn,uiptmx,uetmin,uetmax,ukmin,ukmax,yplumn,yplumx,   &
           iuiptn,nsubla,nsubla+nlogla
    endif
  endif

endif

!===============================================================================
! 10. Formats
!===============================================================================

#if defined(_CS_LANG_FR)

 1000 format(/,' LA NORMALE A LA FACE DE BORD DE PAROI ',I10,/,   &
         ' EST NULLE ; COORDONNEES : ',3E12.5)

 2010 format(/,                                                   &
 3X,'** CONDITIONS AUX LIMITES EN PAROI LISSE',/,                 &
 '   ----------------------------------------',/,                 &
 '------------------------------------------------------------',/,&
 '                                         Minimum     Maximum',/,&
 '------------------------------------------------------------',/,&
 '   Vitesse rel. en paroi    uiptn : ',2E12.5                 ,/,&
 '   Vitesse de frottement    uet   : ',2E12.5                 ,/,&
 '   Vitesse de frottement    uk    : ',2E12.5                 ,/,&
 '   Distance adimensionnelle yplus : ',2E12.5                 ,/,&
 '   ------------------------------------------------------   ',/,&
 '   Nbre de retournements de la vitesse en paroi : ',I10      ,/,&
 '   Nbre de faces en sous couche visqueuse       : ',I10      ,/,&
 '   Nbre de faces de paroi total                 : ',I10      ,/,&
 '------------------------------------------------------------',  &
 /,/)

 2011 format(/,                                                   &
 3X,'** CONDITIONS AUX LIMITES EN PAROI LISSE',/,                 &
 '   ----------------------------------------',/,                 &
 '------------------------------------------------------------',/,&
 '                                         Minimum     Maximum',/,&
 '------------------------------------------------------------',/,&
 '   Vitesse rel. en paroi    uiptn : ',2E12.5                 ,/,&
 '   Vitesse de frottement    uet   : ',2E12.5                 ,/,&
 '   Vitesse de frottement    uk    : ',2E12.5                 ,/,&
 '   Distance adimensionnelle yplus : ',2E12.5                 ,/,&
 '   Sca. thermal de frott.   tstar : ',2E12.5                 ,/,&
 '   Sca. thermal adim. rug.  tplus : ',2E12.5                 ,/,&
 '   ------------------------------------------------------'   ,/,&
 '   Nbre de retournements de la vitesse en paroi : ',I10      ,/,&
 '   Nbre de faces en sous couche visqueuse       : ',I10      ,/,&
 '   Nbre de faces de paroi total                 : ',I10      ,/,&
 '------------------------------------------------------------',  &
 /,/)

 2020 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : RAFFINEMENT INSUFFISANT DU MAILLAGE EN PAROI',/,&
'@    =========                                               ',/,&
'@    Le maillage semble insuffisamment raffine en paroi      ',/,&
'@      pour pouvoir realiser un calcul laminaire.            ',/,&
'@                                                            ',/,&
'@    Le dernier pas de temps auquel ont ete observees de trop',/,&
'@      grandes valeurs de la distance adimensionnelle a la   ',/,&
'@      paroi (yplus) est le pas de temps ',I10                ,/,&
'@                                                            ',/,&
'@    La valeur minimale de yplus doit etre inferieure a la   ',/,&
'@      valeur limite YPLULI = ',E14.5                         ,/,&
'@                                                            ',/,&
'@    Observer la repartition de yplus en paroi (sous EnSight ',/,&
'@      ou ParaView par exemple) pour determiner dans quelle  ',/,&
'@      mesure la qualite des resultats est susceptible d etre',/,&
'@      affectee.')

 2030 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : RAFFINEMENT INSUFFISANT DU MAILLAGE EN PAROI',/,&
'@    =========                                               ',/,&
'@    Le maillage semble insuffisamment raffine en paroi      ',/,&
'@      pour pouvoir realiser un calcul type v2f              ',/,&
'@            (phi-fbar ou BL-v2/k)                           ',/,&
'@                                                            ',/,&
'@    Le dernier pas de temps auquel ont ete observees de trop',/,&
'@      grandes valeurs de la distance adimensionnelle a la   ',/,&
'@      paroi (yplus) est le pas de temps ',I10                ,/,&
'@                                                            ',/,&
'@    La valeur minimale de yplus doit etre inferieure a la   ',/,&
'@      valeur limite YPLULI = ',E14.5                         ,/,&
'@                                                            ',/,&
'@    Observer la repartition de yplus en paroi (sous EnSight ',/,&
'@      ou ParaView par exemple) pour determiner dans quelle  ',/,&
'@      mesure la qualite des resultats est susceptible d etre',/,&
'@      affectee.')

 2040 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : MAILLAGE TROP FIN EN PAROI                  ',/,&
'@    =========                                               ',/,&
'@    Le maillage semble trop raffine en paroi pour utiliser  ',/,&
'@      un modele de turbulence haut Reynolds.                ',/,&
'@                                                            ',/,&
'@    Le dernier pas de temps auquel ont ete observees des    ',/,&
'@      valeurs trop faibles de la distance adimensionnelle a ',/,&
'@      la paroi (yplus) est le pas de temps ',I10             ,/,&
'@                                                            ',/,&
'@    La valeur minimale de yplus doit etre superieure a la   ',/,&
'@      valeur limite YPLULI = ',E14.5                         ,/,&
'@                                                            ',/,&
'@    Observer la repartition de yplus en paroi (sous EnSight ',/,&
'@      ou ParaView par exemple) pour determiner dans quelle  ',/,&
'@      mesure la qualite des resultats est susceptible d etre',/,&
'@      affectee.')
 2050 format(                                                     &
'@                                                            ',/,&
'@    Ce message ne s''affiche qu''aux deux premieres         ',/,&
'@      occurences du probleme et aux deux derniers pas de    ',/,&
'@      temps du calcul. La disparition du message ne signifie',/,&
'@      pas forcement la disparition du probleme.             ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2060 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

#else

 1000 format(/,' THE NORMAL TO THE WALL BOUNDARY FACE ',I10,/,    &
         ' IS NULL; COORDINATES: ',3E12.5)

 2010 format(/,                                                   &
 3X,'** BOUNDARY CONDITIONS FOR SMOOTH WALLS',/,                  &
 '   ---------------------------------------',/,                  &
 '------------------------------------------------------------',/,&
 '                                         Minimum     Maximum',/,&
 '------------------------------------------------------------',/,&
 '   Rel velocity at the wall uiptn : ',2E12.5                 ,/,&
 '   Friction velocity        uet   : ',2E12.5                 ,/,&
 '   Friction velocity        uk    : ',2E12.5                 ,/,&
 '   Dimensionless distance   yplus : ',2E12.5                 ,/,&
 '   ------------------------------------------------------'   ,/,&
 '   Nb of reversal of the velocity at the wall   : ',I10      ,/,&
 '   Nb of faces within the viscous sub-layer     : ',I10      ,/,&
 '   Total number of wall faces                   : ',I10      ,/,&
 '------------------------------------------------------------',  &
 /,/)

 2011 format(/,                                                   &
 3X,'** BOUNDARY CONDITIONS FOR SMOOTH WALLS',/,                  &
 '   ---------------------------------------',/,                  &
 '------------------------------------------------------------',/,&
 '                                         Minimum     Maximum',/,&
 '------------------------------------------------------------',/,&
 '   Rel velocity at the wall uiptn : ',2E12.5                 ,/,&
 '   Friction velocity        uet   : ',2E12.5                 ,/,&
 '   Friction velocity        uk    : ',2E12.5                 ,/,&
 '   Dimensionless distance   yplus : ',2E12.5                 ,/,&
 '   Friction thermal sca.    tstar : ',2E12.5                 ,/,&
 '   Rough dim-less th. sca.  tplus : ',2E12.5                 ,/,&
 '   ------------------------------------------------------'   ,/,&
 '   Nb of reversal of the velocity at the wall   : ',I10      ,/,&
 '   Nb of faces within the viscous sub-layer     : ',I10      ,/,&
 '   Total number of wall faces                   : ',I10      ,/,&
 '------------------------------------------------------------',  &
 /,/)


 2020 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: MESH NOT ENOUGH REFINED AT THE WALL            ',/,&
'@    ========                                                ',/,&
'@    The mesh does not seem to be enough refined at the wall ',/,&
'@      to be able to run a laminar simulation.               ',/,&
'@                                                            ',/,&
'@    The last time step at which too large values for the    ',/,&
'@      dimensionless distance to the wall (yplus) have been  ',/,&
'@      observed is the time step ',I10                        ,/,&
'@                                                            ',/,&
'@    The minimum value for yplus must be lower than the      ',/,&
'@      limit value YPLULI = ',E14.5                           ,/,&
'@                                                            ',/,&
'@    Have a look at the distribution of yplus at the wall    ',/,&
'@      (with EnSight or ParaView for example) to conclude on ',/,&
'@      the way the results quality might be affected.        ')

 2030 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: MESH NOT ENOUGH REFINED AT THE WALL            ',/,&
'@    ========                                                ',/,&
'@    The mesh does not seem to be enough refined at the wall ',/,&
'@      to be able to run a v2f simulation                    ',/,&
'@      (phi-fbar or BL-v2/k)                                 ',/,&
'@                                                            ',/,&
'@    The last time step at which too large values for the    ',/,&
'@      dimensionless distance to the wall (yplus) have been  ',/,&
'@      observed is the time step ',I10                        ,/,&
'@                                                            ',/,&
'@    The minimum value for yplus must be lower than the      ',/,&
'@      limit value YPLULI = ',E14.5                           ,/,&
'@                                                            ',/,&
'@    Have a look at the distribution of yplus at the wall    ',/,&
'@      (with EnSight or ParaView for example) to conclude on ',/,&
'@      the way the results quality might be affected.        ')

 2040 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: MESH TOO REFINED AT THE WALL                   ',/,&
'@    ========                                                ',/,&
'@    The mesh seems to be too refined at the wall to use     ',/,&
'@      a high-Reynolds turbulence model.                     ',/,&
'@                                                            ',/,&
'@    The last time step at which too small values for the    ',/,&
'@      dimensionless distance to the wall (yplus) have been  ',/,&
'@      observed is the time step ',I10                        ,/,&
'@                                                            ',/,&
'@    The minimum value for yplus must be greater than the    ',/,&
'@      limit value YPLULI = ',E14.5                           ,/,&
'@                                                            ',/,&
'@    Have a look at the distribution of yplus at the wall    ',/,&
'@      (with EnSight or ParaView for example) to conclude on ',/,&
'@      the way the results quality might be affected.        ')
 2050 format(                                                     &
'@                                                            ',/,&
'@    This warning is only printed at the first two           ',/,&
'@      occurences of the problem and at the last time step   ',/,&
'@      of the calculation. The vanishing of the message does ',/,&
'@      not necessarily mean the vanishing of the problem.    ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2060 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

#endif

!----
! End
!----

return
end subroutine
