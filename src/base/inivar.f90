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

!> \file inivar.f90
!> \brief Initialization of calculation variables, time step
!> and table that stores distance to the wall
!> by the user (after reading a restart file).
!>
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Arguments
!------------------------------------------------------------------------------
!   mode          name          role
!------------------------------------------------------------------------------
!> \param[in]     nvar          total number of variables
!> \param[in]     nscal         total number of scalars
!______________________________________________________________________________

subroutine inivar &
 ( nvar   , nscal )

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
use parall
use period
use ihmpre
use ppppar
use ppthch
use ppincl
use mesh
use field
use cfpoin, only:ithvar
use cs_c_bindings
use cs_cf_bindings
use cs_f_interfaces

use, intrinsic :: iso_c_binding

use darcy_module

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal

! Local variables

character(len=80) :: chaine
integer          ivar  , iscal
integer          iel
integer          iclip , iok   , ii
integer          kscmin, kscmax, keyvar, n_fields
integer          f_id, f_id_prv, c_id, f_dim
integer          iflid, iflidp
integer          idimf
integer          ivoid, uprtot

double precision valmax, valmin, vfmin , vfmax
double precision vdtmax, vdtmin
double precision xekmin, xepmin, xomgmn, xphmin, xphmax
double precision xnumin
double precision x11min, x22min, x33min
double precision xxp0, xyp0, xzp0
double precision xalmin, xalmax
double precision scmaxp, scminp

double precision rvoid(1)
double precision vvoid(3)

double precision, dimension(:), pointer :: dt
double precision, dimension(:), pointer :: field_s_v
double precision, dimension(:,:), pointer :: field_v_v
double precision, dimension(:), pointer :: cvar_pr
double precision, dimension(:), pointer :: cvar_k, cvar_ep, cvar_al
double precision, dimension(:), pointer :: cvar_phi, cvar_omg, cvar_nusa
double precision, dimension(:), pointer :: cvar_r11, cvar_r22, cvar_r33
double precision, dimension(:,:), pointer :: cvar_rij
double precision, dimension(:), pointer :: cvar_var
double precision, dimension(:), pointer :: cpro_prtot

type(var_cal_opt) :: vcopt

!===============================================================================

!===============================================================================
! 1. Initialization
!===============================================================================

call field_get_val_s_by_name('dt', dt)

call field_get_n_fields(n_fields)

call field_get_key_id("variable_id", keyvar)

! Key ids for clipping
call field_get_key_id("min_scalar_clipping", kscmin)
call field_get_key_id("max_scalar_clipping", kscmax)

iok = 0

!===============================================================================
! 2. ON REPASSE LA MAIN A L'UTILISATEUR POUR LA PROGRAMMATION DES
!    INITIALISATIONS QUI LUI SONT PROPRES
!===============================================================================

! Initialized thermodynamic variables indicator
ithvar = 10000

! Indicateur d'initialisation des scalaires par l'utilisateur
! (mis a 1 si passage dans USINIV ou PPINIV ou dans l'IHM ; a 0 sinon)

iusini = 1

iflidp = -1

do ivar = 1, nvar
  call field_get_key_struct_var_cal_opt(ivarfl(ivar), vcopt)
  if (vcopt%iwgrec.eq.1) then

    if (vcopt%idiff.lt.1) cycle
    iflid = ivarfl(ivar)
    if (iflid.eq.iflidp) cycle
    iflidp = iflid

    if (iand(vcopt%idften, ISOTROPIC_DIFFUSION).ne.0) then
      idimf = 1
    elseif (iand(vcopt%idften, ANISOTROPIC_DIFFUSION).ne.0) then
      idimf = 6
    endif
    call field_get_key_int(iflid, kwgrec, f_id)
    if (idimf.eq.6) then
      call field_get_val_v(f_id, field_v_v)
      do iel = 1, ncelet
        field_v_v(1,iel) = 1.d0
        field_v_v(2,iel) = 1.d0
        field_v_v(3,iel) = 1.d0
        field_v_v(4,iel) = 0.d0
        field_v_v(5,iel) = 0.d0
        field_v_v(6,iel) = 0.d0
      enddo
    else if (idimf.eq.1) then
      call field_get_val_s(f_id, field_s_v)
      do iel = 1, ncelet
        field_s_v(iel) = 1.d0
      enddo
    endif

  endif
enddo

! - Interface Code_Saturne
!   ======================

if (iihmpr.eq.1) then

  call uiiniv (isuite, ippmod(idarcy), ithvar)

endif

!   - Sous-programme utilisateur
!     ==========================

if (ippmod(iphpar).eq.0) then

  call cs_user_f_initialization(nvar, nscal, dt)

  !     Avec l'interface, il peut y avoir eu initialisation,
  !       meme si usiniv n'est pas utilise.
  if (isuite.eq.0 .and. iihmpr.eq.1) then
    iusini = 1
  endif

else

  ! ON FAIT DE LA PHYSIQUE PARTICULIERE
  !   On pourrait remonter la partie init non utilisateur de ppiniv avant lecamo
  !     dans iniva0, mais il faudrait quand meme conserver ici l'appel a
  !     ppiniv car il encapsule les appels aux ss pgm utilisateur similaires a
  !     usiniv.

  iusini = 1

  call ppiniv(nvar, nscal, dt)

  if (ippmod(icompf).ge.0.and.(    isuite.eq.0                 &
                               .or.isuite.eq.1.and.ileaux.eq.0)) then

    if (     ithvar.ne. 60000.and.ithvar.ne.100000                    &
        .and.ithvar.ne.140000.and.ithvar.ne.150000.and.ithvar.ne.210000) then
        write(nfecra,1000) ithvar
        iok = iok + 1
    endif

    ivoid = -1
    call cs_cf_thermo(ithvar, ivoid,  rvoid, rvoid, rvoid, vvoid)

  endif

endif

call user_initialization()


! Pressure / Total pressure initialisation

! Standard:
! If the user has initialized the total pressure Ptot, P* is initialized
! accordingly, only if the user has speficied the reference point.
! (all values of the total pressure have to be initialized).
! Otherwise, the total pressure is initialized using P*,
! Ptot = P* + P0 + rho.g.r

! In case of restart without auxiliary, Ptot is recomputed with P*.
! (For EVM models, the shift by 2/3*rho*k is missing)
! In case of restart with auxiliary, nothing is done.

! Compressible:
! The total pressure field does not need to be defined. The solved pressure is
! the total pressure.

! Ground water flow:
! The field of index iprtot is the pressure head (h = H - z),
! h is only used when gravity is taken into account.

if (ippmod(icompf).lt.0.and.ippmod(idarcy).lt.0) then

  call field_get_val_s(ivarfl(ipr), cvar_pr)
  call field_get_val_s(iprtot, cpro_prtot)

  uprtot = 0

  if (ixyzp0.gt.-1.and.(isuite.eq.0.or.ileaux.eq.0)) then
    uprtot = 1
    do iel = 1, ncel
      if (cpro_prtot(iel).le.-0.5d0*rinfin) then
        uprtot = 0
        exit
      endif
    enddo
  endif

  xxp0   = xyzp0(1)
  xyp0   = xyzp0(2)
  xzp0   = xyzp0(3)

  if (uprtot.gt.0) then
    do iel = 1, ncel
      cvar_pr(iel) =  cpro_prtot(iel)               &
                    - ro0*( gx*(xyzcen(1,iel)-xxp0) &
                    + gy*(xyzcen(2,iel)-xyp0)       &
                    + gz*(xyzcen(3,iel)-xzp0) )     &
                    + pred0 - p0
    enddo
  elseif (isuite.eq.0.or.ileaux.eq.0) then
    do iel = 1, ncel
      cpro_prtot(iel) =  cvar_pr(iel)                  &
                       + ro0*( gx*(xyzcen(1,iel)-xxp0) &
                       + gy*(xyzcen(2,iel)-xyp0)       &
                       + gz*(xyzcen(3,iel)-xzp0) )     &
                       + p0 - pred0
    enddo
  endif

else if ((ippmod(idarcy).ge.0).and.(darcy_gravity.ge.1)) then

  call field_get_val_s(ivarfl(ipr), cvar_pr)
  call field_get_val_s(iprtot, cpro_prtot)

  do iel = 1, ncel
    cpro_prtot(iel) = cvar_pr(iel) - xyzcen(1,iel)*darcy_gravity_x &
                                   - xyzcen(2,iel)*darcy_gravity_y &
                                   - xyzcen(3,iel)*darcy_gravity_z
  enddo

endif

!===============================================================================
! 3.  CLIPPING DES GRANDEURS TURBULENTES (UTILISATEUR OU SUITE)
!     (pour ITYTUR=2, 3, 5 ou 6)
!     Si l'utilisateur est intervenu dans USINIV, PPINIV ou via l'interface
!         et a impose des valeurs "correctes" (au sens k, eps, Rii > 0)
!         on considere qu'il s'agit d'une initialisation admissible,
!         on la clippe pour la rendre coherente avec le clipping du code
!         et on continue le calcul
!     Si l'utilisateur est intervenu dans USINIV, PPINIV ou via l'interface
!         et a impose des valeurs visiblement erronees
!         (k, eps ou Rii < 0), on s'arrete (il s'est sans doute trompe).
!     On adopte le meme traitement en suite de calcul
!       pour assurer un comportement identique en suite entre un calcul
!       ou l'utilisateur modifie une variable avec usiniv (mais pas la
!       turbulence) et un calcul ou l'utilisateur ne modifie pas usiniv.
!     S'il n'y a ni suite ni intervention dans USINIV ou PPINIV ou via l'interface,
!       les grandeurs ont deja ete clippees par defaut, sauf si UREF n'a pas
!       (ou a mal) ete initialise. Dans ce cas on avertit aussi l'utilisateur et on
!       stoppe le calcul.

!     Pour resumer :
!      -en   suite  avec des valeurs positives pour k, eps, Rii : on clippe
!      -avec usiniv ou ppiniv ou interface
!                   avec des valeurs positives pour k, eps, Rii : on clippe
!      -non suite sans usiniv ni ppiniv ni interface avec UREF positif :
!                                      grandeurs par defaut (deja clippees)
!      -non suite sans usiniv ni ppiniv ni interface avec UREF negatif : stop
!      -suite ou usiniv ou ppiniv ou interface
!                   avec une valeur negative de k, eps ou Rii : stop
!                   avec une valeur hors de [0;2] pour phi : stop
!         (on souhaite indiquer a l'utilisateur que son fichier suite est
!          bizarre ou que son initialisation est fausse et qu'il a donc
!          fait au moins une erreur qui peut en cacher d'autres)
!===============================================================================

if (iusini.eq.1.or.isuite.eq.1) then

  if(itytur.eq.2 .or. itytur.eq.5) then

    call field_get_val_s(ivarfl(ik), cvar_k)
    call field_get_val_s(ivarfl(iep), cvar_ep)

    xekmin = cvar_k(1)
    xepmin = cvar_ep(1)
    do iel = 1, ncel
      xekmin = min(xekmin,cvar_k(iel) )
      xepmin = min(xepmin,cvar_ep(iel))
    enddo
    if (irangp.ge.0) then
      call parmin (xekmin)
      call parmin (xepmin)
    endif

    if(xekmin.ge.0.d0.and.xepmin.ge.0.d0) then
      iclip = 1
      call clipke(ncelet, ncel, iclip)
    else
      write(nfecra,3020) xekmin,xepmin
      iok = iok + 1
    endif

    !     En v2-f, phi-fbar ou BL-v2/k, on verifie aussi que phi est
    !     compris entre 0 et 2
    if (itytur.eq.5) then

      call field_get_val_s(ivarfl(iphi), cvar_phi)

      xphmin = cvar_phi(1)
      xphmax = cvar_phi(1)
      do iel = 1, ncel
        xphmin = min(xphmin,cvar_phi(iel) )
        xphmax = max(xphmax,cvar_phi(iel) )
      enddo
      if (irangp.ge.0) then
        call parmin (xphmin)
        call parmax (xphmax)
      endif

      !     Par coherence avec clpv2f, on ne clippe qu'a zero et pas a 2
      !              IF(XPHMIN.LT.0.D0 .OR. XPHMAX.GT.2.D0) THEN
      if(xphmin.lt.0.d0) then
        write(nfecra,3021) xphmin,xphmax
        iok = iok + 1
      endif

      !     En v2-f, BL-v2/k, on verifie aussi que alpha est
      !     compris entre 0 et 1
      if (iturb.eq.51) then
        call field_get_val_s(ivarfl(ial), cvar_al)
        xalmin = cvar_al(1)
        xalmax = cvar_al(1)
        do iel = 1, ncel
          xalmin = min(xalmin,cvar_al(iel) )
          xalmax = max(xalmax,cvar_al(iel) )
        enddo
        if (irangp.ge.0) then
          call parmin (xalmin)
          call parmax (xalmax)
        endif

        if(xalmin.lt.0.d0 .or. xalmax.gt.1.d0) then
          write(nfecra,3022) xalmin,xalmax
          iok = iok + 1
        endif

      endif

    endif

  elseif(itytur.eq.3) then

    call field_get_val_s(ivarfl(iep), cvar_ep)
    if (irijco.eq.1) then
      call field_get_val_v(ivarfl(irij), cvar_rij)

      x11min = cvar_rij(1,1)
      x22min = cvar_rij(2,1)
      x33min = cvar_rij(3,1)
      xepmin = cvar_ep(1)
      do iel = 1, ncel
        x11min = min(x11min,cvar_rij(1,iel))
        x22min = min(x22min,cvar_rij(2,iel))
        x33min = min(x33min,cvar_rij(3,iel))
        xepmin = min(xepmin,cvar_ep(iel) )
      enddo
    else
      call field_get_val_s(ivarfl(ir11), cvar_r11)
      call field_get_val_s(ivarfl(ir22), cvar_r22)
      call field_get_val_s(ivarfl(ir33), cvar_r33)

      x11min = cvar_r11(1)
      x22min = cvar_r22(1)
      x33min = cvar_r33(1)
      xepmin = cvar_ep(1)
      do iel = 1, ncel
        x11min = min(x11min,cvar_r11(iel))
        x22min = min(x22min,cvar_r22(iel))
        x33min = min(x33min,cvar_r33(iel))
        xepmin = min(xepmin,cvar_ep(iel) )
      enddo
    endif
    if (irangp.ge.0) then
      call parmin (x11min)
      call parmin (x22min)
      call parmin (x33min)
      call parmin (xepmin)
    endif
    if (x11min.ge.0.d0.and.x22min.ge.0.d0.and.                  &
         x33min.ge.0.d0.and.xepmin.ge.0.d0 ) then
      iclip = 1
      if (irijco.eq.0) then
        call clprij(ncelet, ncel, iclip)
      endif
    else
      write(nfecra,3030) x11min,x22min,x33min,xepmin
      iok = iok + 1
    endif
    if (iturb.eq.32) then
      call field_get_val_s(ivarfl(ial), cvar_al)
      xalmin = cvar_al(1)
      xalmax = cvar_al(1)
      do iel = 1, ncel
        xalmin = min(xalmin, cvar_al(iel))
        xalmax = max(xalmax, cvar_al(iel))
      enddo
      if (irangp.ge.0) then
        call parmin (xalmin)
        call parmax (xalmax)
      endif
      if (xalmin.lt.0.or.xalmax.gt.1.d0) then
        write(nfecra,3033) xalmin, xalmax
        iok = iok + 1
      endif
    endif

  elseif(iturb.eq.60) then

    call field_get_val_s(ivarfl(ik), cvar_k)
    call field_get_val_s(ivarfl(iomg), cvar_omg)

    xekmin = cvar_k(1)
    xomgmn = cvar_omg(1)
    do iel = 1, ncel
      xekmin = min(xekmin,cvar_k(iel))
      xomgmn = min(xomgmn,cvar_omg(iel))
    enddo
    if (irangp.ge.0) then
      call parmin (xekmin)
      call parmin (xomgmn)
    endif

    !     En k-omega on clippe seulement a 0
    if(xekmin.lt.0.d0 .or. xomgmn.le.0.d0) then
      write(nfecra,3031) xekmin,xomgmn
      iok = iok + 1
    endif

  elseif(iturb.eq.70) then

    call field_get_val_s(ivarfl(inusa), cvar_nusa)

    xnumin = cvar_nusa(1)
    do iel = 1, ncel
      xnumin = min(xnumin,cvar_nusa(iel))
    enddo
    if (irangp.ge.0) then
      call parmin (xnumin)
    endif

    !     En Spalart-Allmaras on clippe seulement a 0
    if(xnumin.lt.0.d0 ) then
      write(nfecra,3032) xnumin
      iok = iok + 1
    endif

  endif

else

  if (iturb.ne.0 .and. iturb.ne.10 .and. itytur.ne.4) then
    if (uref.lt.0.d0) then
      write(nfecra,3039) uref
      iok = iok + 1
    endif
  endif

endif

!===============================================================================
! 4.  CLIPPING DES SCALAIRES (UTILISATEUR OU SUITE)
!     Si l'utilisateur est intervenu dans USINIV ou PPINIV et
!       a impose des valeurs "correctes" (au sens comprises dans des bornes
!         simplifiees a base de 0, scamin, scamax)
!         on considere qu'il s'agit d'une initialisation admissible,
!         on la clippe pour la rendre coherente avec le clipping du code
!         et on continue le calcul
!       si l'utilisateur a impose des valeurs visiblement erronees
!         (au sens comprises dans des bornes simplifiees a base de 0, scamin,
!          scamax), on s'arrete (il s'est sans doute trompe).
!     On adopte le meme traitement en suite de calcul
!       pour assurer un comportement identique en suite entre un calcul
!       ou l'utilisateur modifie une variable avec usiniv (mais pas un
!       scalaire) et un calcul ou l'utilisateur ne modifie pas usiniv.
!     Sinon, les grandeurs ont deja ete clippees apres les init par defaut

!     Pour resumer :
!      -en   suite  avec des valeurs grossierement admissibles : on clippe
!      -avec usiniv ou ppiniv
!                   avec des valeurs grossierement admissibles : on clippe
!      -non suite sans usiniv ni ppiniv :
!                                      grandeurs par defaut (deja clippees)
!      -suite ou usiniv ou ppiniv
!                   avec une valeur grossierement non admissible : stop
!         (on souhaite indiquer a l'utilisateur que son fichier suite est
!          bizarre ou que son initialisation est fausse et qu'il a donc
!          fait au moins une erreur qui peut en cacher d'autres)
!===============================================================================

! On traite tous les scalaires d'abord, car ils peuvent etre necessaires
!     pour clipper les variances

if(nscal.gt.0.and.(iusini.eq.1.or.isuite.eq.1)) then

!     Scalaires non variance

  do ii = 1, nscal
    f_id = ivarfl(isca(ii))
    call field_get_dim(f_id, f_dim)
    if((iscavr(ii).le.0.or.iscavr(ii).gt.nscal).and.f_dim.eq.1) then

      ! Get the min clipping
      call field_get_key_double(f_id, kscmin, scminp)
      call field_get_key_double(f_id, kscmax, scmaxp)

      if (scminp.le.scmaxp) then
        ivar = isca(ii)
        call field_get_val_s(ivarfl(ivar), cvar_var)
        valmax = cvar_var(1)
        valmin = cvar_var(1)
        do iel = 1, ncel
          valmax = max(valmax,cvar_var(iel))
          valmin = min(valmin,cvar_var(iel))
        enddo
        if (irangp.ge.0) then
          call parmax (valmax)
          call parmin (valmin)
        endif

!     Verification de la coherence pour les clippings
!                                           des scalaires non variance.
        if (valmin.ge.scminp.and.valmax.le.scmaxp) then
          iscal = ii
          call clpsca(iscal)
        else
          call field_get_label(ivarfl(isca(ii)), chaine)
          write(nfecra,3040) ii,chaine(1:16),                     &
                             valmin,scminp,valmax,scmaxp
          iok = iok + 1
        endif
      endif

    endif
  enddo


!     Variances

  do ii = 1, nscal
    if(iscavr(ii).gt.0.and.iscavr(ii).le.nscal) then

      ! Get the min clipping
      f_id = ivarfl(isca(ii))
      call field_get_key_double(f_id, kscmin, scminp)
      call field_get_key_double(f_id, kscmax, scmaxp)


      if (scminp.le.scmaxp) then
        ivar = isca(ii)
        call field_get_val_s(ivarfl(ivar), cvar_var)
        valmax = cvar_var(1)
        valmin = cvar_var(1)
        do iel = 1, ncel
          valmax = max(valmax,cvar_var(iel))
          valmin = min(valmin,cvar_var(iel))
        enddo
        if (irangp.ge.0) then
          call parmax (valmax)
          call parmin (valmin)
        endif

!     Verification de la coherence pour les clippings de variance.
!     Pour iclvfl = 1 on ne verifie que > 0 sinon ca va devenir difficile
!     de faire une initialisation correcte.

        if(iclvfl(ii).eq.0) then
!       On pourrait clipper dans le cas ou VALMIN.GE.0, mais ca
!       n'apporterait rien, par definition
          if(valmin.lt.0.d0) then
            call field_get_name(ivarfl(isca(ii)), chaine)
            write(nfecra,3050)ii,chaine(1:16),                     &
                              valmin,scminp,valmax,scmaxp
            iok = iok + 1
          endif
        elseif(iclvfl(ii).eq.1) then
! Ici on clippe pour etre coherent avec la valeur du scalaire
          if(valmin.ge.0.d0) then
            iscal = ii
            call clpsca(iscal)
          else
            call field_get_name(ivarfl(isca(ii)), chaine)
            write(nfecra,3050)ii,chaine(1:16),valmin,scminp,valmax,scmaxp
            iok = iok + 1
          endif
        elseif(iclvfl(ii).eq.2) then
          vfmin = 0.d0
          vfmin = max(scminp, vfmin)
          vfmax = scmaxp
! On pourrait clipper dans le cas ou VALMIN.GE.VFMIN.AND.VALMAX.LE.VFMAX
!     mais ca n'apporterait rien, par definition
          if(valmin.lt.vfmin.or.valmax.gt.vfmax) then
            call field_get_name(ivarfl(isca(ii)), chaine)
            write(nfecra,3051)ii,chaine(1:16),                     &
                              valmin,scminp,valmax,scmaxp,         &
                              ii,iclvfl(ii)
            iok = iok + 1
          endif
        endif
      endif

    endif
  enddo

endif

!===============================================================================
! 5.  IMPRESSIONS DE CONTROLE POUR LES INCONNUES, LE PAS DE TEMPS
!        LE CUMUL DES DUREE POUR LES MOYENNES
!===============================================================================

write(nfecra,2000)

!     Inconnues de calcul : on affiche les bornes
f_id = -1
c_id = 1

do ivar = 1, nvar
  f_id_prv = f_id
  f_id = ivarfl(ivar)
  if (f_id.eq.f_id_prv) then
    c_id = c_id + 1
  else
    c_id = 1
  endif

  call field_get_dim(f_id, f_dim)

  if (f_dim.gt.1) then
    call field_get_val_v(f_id, field_v_v)
  else if (f_dim.eq.1) then
    call field_get_val_s(f_id, field_s_v)
  endif

  if (f_dim.gt.1) then
    valmax = -grand
    valmin =  grand
    do iel = 1, ncel
      valmax = max(valmax, field_v_v(c_id,iel))
      valmin = min(valmin, field_v_v(c_id,iel))
    enddo
  else
    valmax = -grand
    valmin =  grand
    do iel = 1, ncel
      valmax = max(valmax, field_s_v(iel))
      valmin = min(valmin, field_s_v(iel))
    enddo
  endif

  if (irangp.ge.0) then
    call parmax (valmax)
    call parmin (valmin)
  endif
  call field_get_label(f_id, chaine)
  write(nfecra,2010)chaine(1:16),valmin,valmax
enddo
write(nfecra,2020)

if (idtvar.ge.0) then
!     Pas de temps : on affiche les bornes
!                    si < 0 on s'arrete
  vdtmax = -grand
  vdtmin =  grand
  do iel = 1, ncel
    vdtmax = max(vdtmax,dt(iel))
    vdtmin = min(vdtmin,dt(iel))
  enddo
  if (irangp.ge.0) then
    call parmax (vdtmax)
    call parmin (vdtmin)
  endif
  write(nfecra,2010) 'dt', vdtmin, vdtmax
  write(nfecra,2020)

  if (vdtmin.le.zero) then
    write(nfecra,3010) vdtmin
    iok = iok + 1
  endif

endif

!     Cumul du temps associe aux moments : on affiche les bornes
!                                          si < 0 on s'arrete

call time_moment_log_iteration

!===============================================================================
! 6.  ARRET GENERAL SI PB
!===============================================================================

if (iok.gt.0) then
  write(nfecra,3090) iok
  call csexit (1)
endif

write(nfecra,3000)

!----
! Formats
!----


#if defined(_CS_LANG_FR)

 1000 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION:   ARRET A L''INITIALISATION DES VARIABLES    ',/,&
'@    ==========   DANS LA THERMODYNAMIQUE COMPRESSIBLE.      ',/,&
'@                                                            ',/,&
'@    Le calcul ne peut etre execute.                         ',/,&
'@                                                            ',/,&
'@    Valeur de l''indicateur ithvar non prevue (',i10,').    ',/,&
'@                                                            ',/,&
'@    Deux et seulement deux variables parmi                  ',/,&
'@    P, rho, T et E (sauf T et E) doivent etre imposees      ',/,&
'@    a l''initialisation dans le GUI ou dans la routine      ',/,&
'@    d''initialisation (cs_user_initialisation.f90) ou       ',/,&
'@    encore via les deux.                                    ',/,&
'@                                                            ',/,&
'@    Verifier que ithvar n''ait pas ete partiellement        ',/,&
'@    positionnee par le GUI de maniere non consistante       ',/,&
'@    avec cs_user_initialisation.f90.                        ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 2000 format(                                                     &
                                                                /,&
' -----------------------------------------------------------', /,&
                                                                /,&
                                                                /,&
' ** INITIALISATION DES VARIABLES',                             /,&
'    ----------------------------',                             /,&
                                                                /,&
' -----------------------------------------',                   /,&
'  Variable          Valeur min  Valeur max',                   /,&
' -----------------------------------------                   '  )
 2010 format(                                                     &
 2x,     a16,      e12.4,      e12.4                             )
 2020 format(                                                     &
' ---------------------------------',                           /)

 3000 format(/,/,                                                 &
'-------------------------------------------------------------',/)
 3010 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''INITIALISATION DES VARIABLES',     /,&
'@    =========',                                               /,&
'@    PAS DE TEMPS NEGATIF OU NUL',                             /,&
'@',                                                            /,&
'@  Le calcul ne peut etre execute.',                           /,&
'@',                                                            /,&
'@  La valeur minimale du pas de temps DT est ',e14.5,          /,&
'@  Verifier l''initialisation ou le fichier de reprise.',      /,&
'@',                                                            /,&
'@  Dans le cas ou les valeurs lues dans le fichier',           /,&
'@    de reprise sont incorrectes, on peut les modifier via',   /,&
'@    cs_user_initialization.f90 ou via l''interface.',         /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 3020 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''INITIALISATION DES VARIABLES',     /,&
'@    =========',                                               /,&
'@     TURBULENCE NEGATIVE OU NULLE',                           /,&
'@',                                                            /,&
'@  Le calcul ne peut etre execute.',                           /,&
'@',                                                            /,&
'@   Valeur minimale de k       = ',e14.5,                      /,&
'@   Valeur minimale de epsilon = ',e14.5,                      /,&
'@',                                                            /,&
'@  Verifier l''initialisation, le fichier de reprise,',        /,&
'@    ou bien la valeur de UREF.',                              /,&
'@  Dans le cas ou les valeurs lues dans le fichier',           /,&
'@    de reprise sont incorrectes, on peut les modifier via',   /,&
'@    cs_user_initialization.f90 ou via l''interface.',         /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 3021 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''INITIALISATION DES VARIABLES',     /,&
'@    =========',                                               /,&
'@     VARIABLE PHI DU V2F (PHI_FBAR ou BL-V2/K)',              /,&
'@     HORS DES BORNES [0;2]',                                  /,&
'@',                                                            /,&
'@  Le calcul ne peut etre execute.',                           /,&
'@',                                                            /,&
'@   Valeur minimale de phi     = ',e14.5,                      /,&
'@   Valeur maximale de phi     = ',e14.5,                      /,&
'@',                                                            /,&
'@  Verifier l''initialisation ou le fichier de reprise.',      /,&
'@  Dans le cas ou les valeurs lues dans le fichier',           /,&
'@    de reprise sont incorrectes, on peut les modifier via',   /,&
'@    cs_user_initialization.f90 ou via l''interface.',         /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 3022 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''INITIALISATION DES VARIABLES',     /,&
'@    =========',                                               /,&
'@     VARIABLE ALPHA DU V2F (BL-V2/K)',                        /,&
'@     HORS DES BORNES [0;1]',                                  /,&
'@',                                                            /,&
'@  Le calcul ne peut etre execute.',                           /,&
'@',                                                            /,&
'@   Valeur minimale de alpha   = ',e14.5,                      /,&
'@   Valeur maximale de alpha   = ',e14.5,                      /,&
'@',                                                            /,&
'@  Verifier l''initialisation ou le fichier de reprise.',      /,&
'@  Dans le cas ou les valeurs lues dans le fichier',           /,&
'@    de reprise sont incorrectes, on peut les modifier via',   /,&
'@    cs_user_initialization.f90 ou via l''interface.',         /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 3030 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''INITIALISATION DES VARIABLES',     /,&
'@    =========',                                               /,&
'@     TURBULENCE NEGATIVE OU NULLE',                           /,&
'@',                                                            /,&
'@  Le calcul ne peut etre execute.',                           /,&
'@',                                                            /,&
'@   Valeur minimale de R11     = ',e14.5,                      /,&
'@   Valeur minimale de R22     = ',e14.5,                      /,&
'@   Valeur minimale de R33     = ',e14.5,                      /,&
'@   Valeur minimale de epsilon = ',e14.5,                      /,&
'@',                                                            /,&
'@  Verifier l''initialisation, le fichier de reprise,',        /,&
'@    ou bien la valeur de UREF.',                              /,&
'@  Dans le cas ou les valeurs lues dans le fichier',           /,&
'@    de reprise sont incorrectes, on peut les modifier via',   /,&
'@    cs_user_initialization.f90 ou via l''interface.',         /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 3031 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''INITIALISATION DES VARIABLES',     /,&
'@    =========',                                               /,&
'@    TURBULENCE NEGATIVE OU NULLE',                            /,&
'@',                                                            /,&
'@  Le calcul ne peut etre execute.',                           /,&
'@',                                                            /,&
'@   Valeur minimale de k       = ',e14.5,                      /,&
'@   Valeur minimale de omega   = ',e14.5,                      /,&
'@',                                                            /,&
'@  Verifier l''initialisation, le fichier de reprise,',        /,&
'@    ou bien la valeur de UREF.',                              /,&
'@  Dans le cas ou les valeurs lues dans le fichier',           /,&
'@    de reprise sont incorrectes, on peut les modifier via',   /,&
'@    cs_user_initialization.f90 ou via l''interface.',         /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 3032 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''INITIALISATION DES VARIABLES',     /,&
'@    =========',                                               /,&
'@    TURBULENCE NEGATIVE OU NULLE',                            /,&
'@',                                                            /,&
'@  Le calcul ne peut etre execute.',                           /,&
'@',                                                            /,&
'@   Valeur minimale de nu      = ',e14.5,                      /,&
'@',                                                            /,&
'@  Verifier l''initialisation, le fichier de reprise,',        /,&
'@    ou bien la valeur de UREF.',                              /,&
'@  Dans le cas ou les valeurs lues dans le fichier',           /,&
'@    de reprise sont incorrectes, on peut les modifier via',   /,&
'@    cs_user_initialization.f90 ou via l''interface.',         /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 3039 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''INITIALISATION DES VARIABLES',     /,&
'@    =========',                                               /,&
'@    LA VITESSE DE REFERENCE UREF N''A PAS ETE INITIALISEE',   /,&
'@    OU A ETE MAL INITIALISEE (VALEUR NEGATIVE).',             /,&
'@    ELLE VAUT ICI ',e14.5,                                    /,&
'@',                                                            /,&
'@  La turbulence n''a pas pu etre initialisee',                /,&
'@  Corriger la valeur de UREF ou bien initialiser',            /,&
'@    directement la turbulence.',                              /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 3033 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''INITIALISATION DES VARIABLES',     /,&
'@    =========',                                               /,&
'@     EBRSM ALPHA<0 OU ALPHA>1',                               /,&
'@',                                                            /,&
'@  Le calcul ne peut etre execute.',                           /,&
'@',                                                            /,&
'@   Valeur minimale de alpha   = ',e14.5,                      /,&
'@   Valeur maximale de alpha   = ',e14.5,                      /,&
'@',                                                            /,&
'@  Verifier l''initialisation ou le fichier de reprise.',      /,&
'@  Dans le cas ou les valeurs lues dans le fichier',           /,&
'@    de reprise sont incorrectes, on peut les modifier via',   /,&
'@    cs_user_initialization.f90 ou via l''interface.',         /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)

 3040 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''INITIALISATION DES VARIABLES',     /,&
'@    =========',                                               /,&
'@     GRANDEUR SCALAIRE HORS BORNES',                          /,&
'@',                                                            /,&
'@  Le calcul ne peut etre execute.',                           /,&
'@',                                                            /,&
'@  Scalaire numero ',i10,' : ',a16,                            /,&
'@  Valeur minimale             = ',e14.5,                      /,&
'@    Clipping demande a SCAMIN = ',e14.5,                      /,&
'@  Valeur maximale             = ',e14.5,                      /,&
'@    Clipping demande a SCAMAX = ',e14.5,                      /,&
'@  Les valeurs extremes ne sont pas coherentes avec les',      /,&
'@    limites SCAMIN et SCAMAX imposees.',                      /,&
'@',                                                            /,&
'@  Verifier l''initialisation ou le fichier de reprise.',      /,&
'@  Dans le cas ou les valeurs lues dans le fichier',           /,&
'@    de reprise sont incorrectes, on peut les modifier via',   /,&
'@    cs_user_initialization.f90 ou via l''interface.',         /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 3050 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''INITIALISATION DES VARIABLES',     /,&
'@    =========',                                               /,&
'@     VARIANCE NEGATIVE',                                      /,&
'@',                                                            /,&
'@  Le calcul ne peut etre execute.',                           /,&
'@',                                                            /,&
'@  Scalaire numero ',i10,' : ',a16,                            /,&
'@  Valeur minimale             = ',e14.5,                      /,&
'@  Le scalaire indique ci-dessus est une variance (ISCAVR est',/,&
'@    positif) mais l''initialisation imposee comporte',        /,&
'@    des valeurs negatives.',                                  /,&
'@',                                                            /,&
'@  Verifier l''initialisation ou le fichier de reprise.',      /,&
'@  Dans le cas ou les valeurs lues dans le fichier',           /,&
'@    de reprise sont incorrectes, on peut les modifier via',   /,&
'@    cs_user_initialization.f90 ou via l''interface.',         /,&
'@  Verifier la definition des variances.',                     /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 3051 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''INITIALISATION DES VARIABLES',     /,&
'@    =========',                                               /,&
'@     VARIANCE HORS BORNES',                                   /,&
'@',                                                            /,&
'@  Le calcul ne peut etre execute.',                           /,&
'@',                                                            /,&
'@  Scalaire numero ',i10,' : ',a16,                            /,&
'@  Valeur minimale             = ',e14.5,                      /,&
'@    Clipping demande a SCAMIN = ',e14.5,                      /,&
'@  Valeur maximale             = ',e14.5,                      /,&
'@    Clipping demande a SCAMAX = ',e14.5,                      /,&
'@  Le scalaire indique ci-dessus est une variance (ISCAVR est',/,&
'@    postif) mais l initialisation imposee comporte',          /,&
'@    des valeurs situees hors des bornes',                     /,&
'@    SCAMIN, SCAMAX ou inferieures a 0 et le mode de clipping',/,&
'@    demande est ICLVFL(',i10,') = ',i10,                      /,&
'@',                                                            /,&
'@  Verifier l''initialisation ou le fichier de reprise.',      /,&
'@  Dans le cas ou les valeurs lues dans le fichier',           /,&
'@    de reprise sont incorrectes, on peut les modifier via',   /,&
'@    cs_user_initialization.f90 ou via l''interface.',         /,&
'@  Verifier la definition des variances et le mode de',        /,&
'@    clipping demande.',                                       /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 3090 format(                                                     &
'@',                                                            /,&
'@',                                                            /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''INITIALISATION DES VARIABLES',     /,&
'@    =========',                                               /,&
'@',                                                            /,&
'@    L INITIALISATION DES VARIABLES EST INCOMPLETE OU',        /,&
'@      INCOHERENTE AVEC LES VALEURS DES PARAMETRES DE CALCUL', /,&
'@',                                                            /,&
'@  Le calcul ne sera pas execute (',i10,' erreurs).',          /,&
'@',                                                            /,&
'@  Se reporter aux impressions precedentes pour plus de',      /,&
'@    renseignements.',                                         /,&
'@  Attention a l''initialisation du pas de temps',             /,&
'@                                de la turbulence',            /,&
'@                                des scalaires et variances',  /,&
'@                                des moyennes temporelles',    /,&
'@',                                                            /,&
'@  Verifier l''initialisation ou le fichier de reprise.',      /,&
'@  Dans le cas ou les valeurs lues dans le fichier',           /,&
'@    de reprise sont incorrectes, on peut les modifier via',   /,&
'@    cs_user_initialization.f90 ou via l''interface.',         /,&
'@  Verifier les parametres de calcul.',                        /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)

#else

 1000 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING :     stop in compressible thermodynamics at    ',/,&
'@    =========     initialisation.                           ',/,&
'@                                                            ',/,&
'@    The computation will stop.                              ',/,&
'@                                                            ',/,&
'@    Unexpected value of the indicator ithvar (',i10,').     ',/,&
'@                                                            ',/,&
'@    Two and only two independant variables among            ',/,&
'@    P, rho, T and E (except T and E) should be imposed at   ',/,&
'@    the initialisation in the GUI or in the user subroutine ',/,&
'@    of initialisation (cs_user_initialisation.f90) or       ',/,&
'@    in both.                                                ',/,&
'@                                                            ',/,&
'@    Check if iccfth has not been partially set through the  ',/,&
'@    GUI in a non consistant manner with                     ',/,&
'@    cs_user_initialisation.f90.                             ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 2000 format(                                                     &
                                                                /,&
' -----------------------------------------------------------', /,&
                                                                /,&
                                                                /,&
' ** VARIABLES INITIALIZATION',                                 /,&
'    ------------------------',                                 /,&
                                                                /,&
' -----------------------------------------',                   /,&
'  Variable          Min. value  Max. value',                   /,&
' -----------------------------------------                   '  )
 2010 format(                                                     &
 2x,     a16,      e12.4,      e12.4                             )
 2020 format(                                                     &
' ---------------------------------',                           /)

 3000 format(/,/,                                                 &
'-------------------------------------------------------------',/)
 3010 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ WARNING: ABORT IN THE VARIABLES INITIALIZATION',          /,&
'@    ========',                                                /,&
'@    NEGATIVE OR NULL TIME STEP',                              /,&
'@',                                                            /,&
'@  The calculation will not be run.',                          /,&
'@',                                                            /,&
'@  The minimum value of the time-step dt is ',e14.5,           /,&
'@  Verify the initialization or the restart file.',            /,&
'@  In the case where the values read in the restart file',     /,&
'@    are incorrect, they may be modified with',                /,&
'@    cs_user_initialization.f90 or with the interface.',       /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 3020 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ WARNING: ABORT IN THE VARIABLES INITIALIZATION',          /,&
'@    ========',                                                /,&
'@     NEGATIVE OR NULL TURBULENCE',                            /,&
'@',                                                            /,&
'@  The calculation will not be run.',                          /,&
'@',                                                            /,&
'@   Minimum value of k       = ',e14.5,                        /,&
'@   Minimum value of epsilon = ',e14.5,                        /,&
'@',                                                            /,&
'@  Verify the initialization, the restart file,',              /,&
'@    and the value of UREF.',                                  /,&
'@  In the case where the values read in the restart file',     /,&
'@    are incorrect, they may be modified with',                /,&
'@    cs_user_initialization.f90 or with the interface.',       /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 3021 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ WARNING: ABORT IN THE VARIABLES INITIALIZATION',          /,&
'@    ========',                                                /,&
'@     PHI VARIABLE OF V2F (PHI_FBAR or BL-V2/K)',              /,&
'@     OUT OF BOUNDS [0;2]',                                    /,&
'@',                                                            /,&
'@  The calculation will not be run.',                          /,&
'@',                                                            /,&
'@   Minimum value of phi = ',e14.5,                            /,&
'@   Maximum value of phi = ',e14.5,                            /,&
'@',                                                            /,&
'@  Verify the initialization or the restart file.',            /,&
'@  In the case where the values read in the restart file',     /,&
'@    are incorrect, they may be modified with',                /,&
'@    cs_user_initialization.f90 or with the interface.',       /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 3022 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ WARNING: ABORT IN THE VARIABLES INITIALIZATION',          /,&
'@    ========',                                                /,&
'@     ALPHA VARIABLE OF V2F (BL-V2/K)',                        /,&
'@     OUT OF BOUNDS [0;1]',                                    /,&
'@',                                                            /,&
'@  The calculation will not be run.',                          /,&
'@',                                                            /,&
'@   Minimum value of alpha = ',e14.5,                          /,&
'@   Maximum value of alpha = ',e14.5,                          /,&
'@',                                                            /,&
'@  Verify the initialization or the restart file.',            /,&
'@  In the case where the values read in the restart file',     /,&
'@    are incorrect, they may be modified with',                /,&
'@    cs_user_initialization.f90 or with the interface.',       /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 3030 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ WARNING: ABORT IN THE VARIABLES INITIALIZATION',          /,&
'@    ========',                                                /,&
'@     NEGATIVE OR NULL TURBULENCE',                            /,&
'@',                                                            /,&
'@  The calculation will not be run.',                          /,&
'@',                                                            /,&
'@   Minimum value of R11     = ',e14.5,                        /,&
'@   Minimum value of R22     = ',e14.5,                        /,&
'@   Minimum value of R33     = ',e14.5,                        /,&
'@   Minimum value of epsilon = ',e14.5,                        /,&
'@',                                                            /,&
'@  Verify the initialization, the restart file,',              /,&
'@    and the value of UREF.',                                  /,&
'@  In the case where the values read in the restart file',     /,&
'@    are incorrect, they may be modified with',                /,&
'@    cs_user_initialization.f90 or with the interface.',       /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 3031 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ WARNING: ABORT IN THE VARIABLES INITIALIZATION',          /,&
'@    ========',                                                /,&
'@     NEGATIVE OR NULL TURBULENCE',                            /,&
'@',                                                            /,&
'@  The calculation will not be run.',                          /,&
'@',                                                            /,&
'@   Minimum value of k       = ',e14.5,                        /,&
'@   Minimum value of omega   = ',e14.5,                        /,&
'@',                                                            /,&
'@  Verify the initialization, the restart file,',              /,&
'@    and the value of UREF.',                                  /,&
'@  In the case where the values read in the restart file',     /,&
'@    are incorrect, they may be modified with',                /,&
'@    cs_user_initialization.f90 or with the interface.',       /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 3032 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ WARNING: ABORT IN THE VARIABLES INITIALIZATION',          /,&
'@    ========',                                                /,&
'@     NEGATIVE OR NULL TURBULENCE',                            /,&
'@',                                                            /,&
'@  The calculation will not be run.',                          /,&
'@',                                                            /,&
'@   Minimum value of nu      = ',e14.5,                        /,&
'@',                                                            /,&
'@  Verify the initialization, the restart file,',              /,&
'@    and the value of UREF.',                                  /,&
'@  In the case where the values read in the restart file',     /,&
'@    are incorrect, they may be modified with',                /,&
'@    cs_user_initialization.f90 or with the interface.',       /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 3039 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ WARNING: ABORT IN THE VARIABLES INITIALIZATION',          /,&
'@    ========',                                                /,&
'@    THE REFERENCE VELOCITY UREF HAS NOT BEEN INITIALIZED',    /,&
'@    OR HAS NOT BEEN CORRECTLY INITIALIZED (NEGATIVE VALUE)',  /,&
'@    ITS VALUE IS ',e14.5,                                     /,&
'@',                                                            /,&
'@  The turbulence cannot be initialized',                      /,&
'@  Correct the value of UREF or initialize the turbulence',    /,&
'@    directly with cs_user_initialization.f90',                /,&
'@    or with the interface.',                                  /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 3033 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ WARNING: ABORT IN THE VARIABLES INITIALIZATION',          /,&
'@    ========',                                                /,&
'@     EBRSM ALPHA<0 OU ALPHA>1',                               /,&
'@',                                                            /,&
'@  The calculation will not be run.',                          /,&
'@',                                                            /,&
'@   Minimum value of alpha   = ',e14.5,                        /,&
'@   Maximum value of alpha   = ',e14.5,                        /,&
'@',                                                            /,&
'@  Verify the initialization or the restart file',             /,&
'@  In the case where the values read in the restart file',     /,&
'@    are incorrect, they may be modified with',                /,&
'@    cs_user_initialization.f90 or with the interface.',       /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 3040 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ WARNING: ABORT IN THE VARIABLES INITIALIZATION',          /,&
'@    ========',                                                /,&
'@     SCALAR QUANTITIES OUT OF BOUNDS',                        /,&
'@',                                                            /,&
'@  The calculation will not be run.',                          /,&
'@',                                                            /,&
'@  Scalar number ',i10,': ',a16,                               /,&
'@  Minimum value                = ',e14.5,                     /,&
'@    Desired clipping at SCAMIN = ',e14.5,                     /,&
'@  Maximum value                = ',e14.5,                     /,&
'@    Desired clipping at SCAMAX = ',e14.5,                     /,&
'@  The bounds are not coherent with the limits SCAMIN and',    /,&
'@    SCAMAX.',                                                 /,&
'@',                                                            /,&
'@  Verify the initialization and the restart file.',           /,&
'@  In the case where the values read in the restart file',     /,&
'@    are incorrect, they may be modified with',                /,&
'@    cs_user_initialization.f90 or with the interface.',       /,&
'@  Verify the clipping values.',                               /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 3050 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ WARNING: ABORT IN THE VARIABLES INITIALIZATION',          /,&
'@    ========',                                                /,&
'@     NEGATIVE VARIANCE',                                      /,&
'@',                                                            /,&
'@  The calculation will not be run.',                          /,&
'@',                                                            /,&
'@  Scalar number ',i10,': ',a16,                               /,&
'@  Minimum value               = ',e14.5,                      /,&
'@  This scalar is a variance (ISCAVR is positive)',            /,&
'@    but the initialization has some negative values.',        /,&
'@',                                                            /,&
'@  Verify the initialization and the restart file.',           /,&
'@  In the case where the values read in the restart file',     /,&
'@    are incorrect, they may be modified with',                /,&
'@    cs_user_initialization.f90 or with the interface.',       /,&
'@  Verify the variance definition.',                           /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 3051 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ WARNING: ABORT IN THE VARIABLES INITIALIZATION',          /,&
'@    ========',                                                /,&
'@     VARIANCE OUT OF BOUNDS',                                 /,&
'@',                                                            /,&
'@  The calculation will not be run.',                          /,&
'@',                                                            /,&
'@  Scalar number ',i10,': ',a16,                               /,&
'@  Minimum value                = ',e14.5,                     /,&
'@    Desired clipping at SCAMIN = ',e14.5,                     /,&
'@  Maximum value                = ',e14.5,                     /,&
'@    Desired clipping at SCAMAX = ',e14.5,                     /,&
'@  This scalar is a variance (ISCAVR is positive)',            /,&
'@    but the initialization has some values out',              /,&
'@    of the bounds SCAMIN, SCAMAX or lower than 0 and the',    /,&
'@    desired clipping mode is ICLVFL(',i10,') = ',i10,         /,&
'@',                                                            /,&
'@  Verify the initialization and the restart file.',           /,&
'@  In the case where the values read in the restart file',     /,&
'@    are incorrect, they may be modified with',                /,&
'@    cs_user_initialization.f90 or with the interface.',       /,&
'@  Verify the variance definition and the clipping mode.',     /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 3090 format(                                                     &
'@',                                                            /,&
'@',                                                            /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ WARNING: ABORT IN THE VARIABLES INITIALIZATION',          /,&
'@    ========',                                                /,&
'@',                                                            /,&
'@    THE VARIABLES INITIALIZATION IS INCOMPLETE OR',           /,&
'@    INCOHERENT WITH THE PARAMETERS VALUE OF THE CALCULATION', /,&
'@',                                                            /,&
'@  The calculation will not be run (',i10,' errors).',         /,&
'@',                                                            /,&
'@  Refer to the previous warnings for further information.',   /,&
'@  Pay attention to the initialization of',                    /,&
'@                                the time-step',               /,&
'@                                the turbulence',              /,&
'@                                the scalars and variances',   /,&
'@                                the time averages',           /,&
'@',                                                            /,&
'@  Verify the initialization and the restart file.',           /,&
'@  In the case where the values read in the restart file',     /,&
'@    are incorrect, they may be modified with',                /,&
'@    cs_user_initialization.f90 or with the interface.',       /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)

#endif

!----
! End
!----

return
end subroutine
