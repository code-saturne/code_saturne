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

!> \file tridim.f90
!> \brief Resolution of incompressible Navier Stokes and scalar transport
!> equations for a time step.
!>
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Arguments
!------------------------------------------------------------------------------
!   mode          name          role
!------------------------------------------------------------------------------
!> \param[in]     itrale        ALE iteration number
!> \param[in]     nvar          total number of variables
!> \param[in]     nscal         total number of scalars
!> \param[in]     dt            time step (per cell)
!______________________________________________________________________________

subroutine tridim &
 ( itrale ,                                                       &
   nvar   , nscal  ,                                              &
   dt     )

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use optcal
use entsor
use cstphy
use cstnum
use pointe
use albase
use alstru
use alaste
use parall
use period
use ppppar
use ppthch
use ppincl
use cpincl
use coincl
use atincl
use ctincl
use atsoil
use lagran
use vorinc
use ihmpre
use radiat
use cplsat
use ppcpfu
use cs_fuel_incl
use mesh
use field
use rotation
use turbomachinery
use darcy_module
use cs_f_interfaces
use cs_c_bindings
use cs_tagmr, only: rob, condb, cpb, hext, text
use cs_tagms, only: t_metal, tmet0
use cs_nz_tagmr
use cs_nz_condensation
use turbomachinery

! les " use pp* " ne servent que pour recuperer le pointeur IIZFPP

!===============================================================================

implicit none

! Arguments

integer          itrale
integer          nvar   , nscal

double precision, pointer, dimension(:)   :: dt
double precision, pointer, dimension(:,:) :: frcxt => null()

! Local variables

logical          must_return

integer          iel   , ifac  , ivar  , iscal , iappel, n_fans
integer          iok   , nfld  , f_id  , f_dim  , f_type
integer          nbccou
integer          ntrela
integer          icmst
integer          st_id

integer          isvhb, iz
integer          ii
integer          iterns, inslst, icvrge
integer          italim, itrfin, itrfup, ineefl
integer          ielpdc, iflmas, iflmab
integer          kcpsyr, icpsyr

double precision cvcst
double precision xxp0, xyp0, xzp0
double precision relaxk, relaxe, relaxw, relaxn
double precision hdls(6)

double precision, save :: tpar, tmet

integer          ipass
data             ipass /0/
save             ipass

integer, pointer, dimension(:,:) :: icodcl
integer, allocatable, dimension(:) :: isostd

double precision, pointer, dimension(:) :: flmalf, flmalb, xprale
double precision, pointer, dimension(:,:) :: cofale
double precision, pointer, dimension(:,:) :: dttens
double precision, pointer, dimension(:,:,:) :: rcodcl
double precision, pointer, dimension(:) :: hbord, theipb
double precision, pointer, dimension(:) :: visvdr
double precision, allocatable, dimension(:) :: prdv2f
double precision, allocatable, dimension(:) :: mass_source
double precision, dimension(:), pointer :: brom, crom

double precision, pointer, dimension(:,:) :: trava
double precision, dimension(:,:), pointer :: vel
double precision, dimension(:,:), pointer :: cvar_vec
double precision, dimension(:), pointer :: cvar_sca
double precision, dimension(:), pointer :: cvar_pr
double precision, dimension(:), pointer :: cvar_k, cvara_k, cvar_ep, cvara_ep
double precision, dimension(:), pointer :: cvar_omg, cvara_omg
double precision, dimension(:), pointer :: cvar_nusa, cvara_nusa
double precision, dimension(:), pointer :: cvar_r11, cvar_r22, cvar_r33
double precision, dimension(:), pointer :: cvar_r12, cvar_r13, cvar_r23
double precision, dimension(:), pointer :: cpro_prtot
double precision, dimension(:), pointer :: cvar_scalt, cvar_totwt

! Darcy
integer mbrom
double precision, dimension(:), pointer :: cpro_delay, cpro_capacity, cpro_sat
double precision, dimension(:), pointer :: cproa_delay, cproa_capacity
double precision, dimension(:), pointer :: cproa_sat
double precision, dimension(:), pointer :: i_mass_flux, b_mass_flux

double precision, dimension(:), pointer :: coefap, cofafp, cofbfp

type(gwf_soilwater_partition) :: sorption_scal

type(var_cal_opt) :: vcopt, vcopt_u, vcopt_p

!===============================================================================
! Interfaces
!===============================================================================

interface

  subroutine condli &
  ( nvar   , nscal  , iterns ,                                     &
    isvhb  ,                                                       &
    itrale , italim , itrfin , ineefl , itrfup ,                   &
    flmalf , flmalb , cofale , xprale ,                            &
    icodcl , isostd ,                                              &
    dt     , rcodcl ,                                              &
    visvdr , hbord  , theipb )

    use mesh, only: nfac, nfabor

    implicit none

    integer          nvar, nscal, iterns, isvhb
    integer          itrale , italim , itrfin , ineefl , itrfup

    double precision, pointer, dimension(:) :: flmalf, flmalb, xprale
    double precision, pointer, dimension(:,:) :: cofale
    integer, pointer, dimension(:,:) :: icodcl
    integer, dimension(nfabor+1) :: isostd
    double precision, pointer, dimension(:) :: dt
    double precision, pointer, dimension(:,:,:) :: rcodcl
    double precision, pointer, dimension(:) :: visvdr, hbord, theipb

  end subroutine condli

  subroutine navstv &
  ( nvar   , nscal  , iterns , icvrge , itrale ,                   &
    isostd ,                                                       &
    dt     ,                                                       &
    frcxt  ,                                                       &
    trava  )

    use mesh, only: nfabor

    implicit none

    integer          nvar   , nscal  , iterns , icvrge , itrale

    integer          isostd(nfabor+1)

    double precision, pointer, dimension(:)   :: dt
    double precision, pointer, dimension(:,:) :: frcxt
    double precision, pointer, dimension(:,:) :: trava

  end subroutine navstv

  subroutine richards(icvrge, dt)

    implicit none

    integer  icvrge
    double precision, pointer, dimension(:)   :: dt

  end subroutine richards

  subroutine cs_lagr_head_losses(n_hl_cells, cell_ids, bc_type, cku) &
    bind(C, name='cs_lagr_head_losses')
    use, intrinsic :: iso_c_binding
    implicit none
    integer(c_int), value :: n_hl_cells
    integer(c_int), dimension(*), intent(in) :: cell_ids
    integer(c_int), dimension(*) :: bc_type
    real(kind=c_double), dimension(*) :: cku
  end subroutine cs_lagr_head_losses

end interface

!===============================================================================

! Map field arrays
call field_get_val_v(ivarfl(iu), vel)

! Number of fields
call field_get_n_fields(nfld)

call field_get_key_struct_var_cal_opt(ivarfl(iu), vcopt_u)
call field_get_key_struct_var_cal_opt(ivarfl(ipr), vcopt_p)

!===============================================================================
! 1. Initialisation
!===============================================================================

allocate(isostd(nfabor+1))

must_return = .false.

if (vcopt_u%iwarni.ge.1) then
  write(nfecra,1000)
endif

ipass = ipass + 1

! --- Indicateur de stockage d'un scalaire et de son coef
!     d'echange associe.
!     Pour le moment, on stocke uniquement dans le cas couplage SYRTHES.
!     ISVHB donne le numero du scalaire (on suppose qu'il n'y en a qu'un).
!     Dans le cas ou on a un couplage avec le module thermique 1D en paroi,
!     on utilise le meme scalaire que celui qui sert a Syrthes (s'il y a
!     couplage Syrthes), sinon on stocke le scalaire thermique.

call field_get_key_id("syrthes_coupling", kcpsyr)

call nbcsyr (nbccou)
isvhb = 0
if (nbccou .ge. 1) then
  do iscal = 1, nscal
    call field_get_key_int(ivarfl(isca(iscal)), kcpsyr, icpsyr)
    if(icpsyr.eq.1) then
      isvhb = iscal
    endif
  enddo
endif

if ((nfpt1t.gt.0).and.(nbccou.le.0)) then
  isvhb = iscalt
endif

call field_get_val_s(ivarfl(ipr), cvar_pr)

if (iphydr.eq.1) then
  call field_get_val_v_by_name('volume_forces', frcxt)
else
  frcxt => rvoid2
endif

!===============================================================================
! 2.  AU DEBUT DU CALCUL ON REINITIALISE LA PRESSION
!===============================================================================

if (ippmod(idarcy).eq.-1) then

! On le fait sur ntinit pas de temps, car souvent, le champ de flux de masse
!   initial n'est pas a divergence nulle (CL incluses) et l'obtention
!   d'un flux a divergence nulle coherent avec la contrainte stationnaire
!   peut prendre quelques pas de temps.
! Remarquer que la pression est rattrapee dans l'etape de stokes.
! On ne le fait pas dans le cas de la prise en compte de la pression
!   hydrostatique, ni dans le cas du compressible

  if( ntcabs.le.ntinit .and. isuite.eq.0 .and. (iphydr.eq.0.or.iphydr.eq.2)    &
                  .and. ippmod(icompf).lt.0                               &
                  .and. idilat .le.1                 ) then

    if(vcopt_p%iwarni.ge.2) then
      write(nfecra,2000) ntcabs
    endif
    call field_get_val_s(iprtot, cpro_prtot)
    xxp0   = xyzp0(1)
    xyp0   = xyzp0(2)
    xzp0   = xyzp0(3)
    do iel = 1, ncel
      cvar_pr(iel) = pred0
      cpro_prtot(iel) = p0                                   &
           + ro0*( gx*(xyzcen(1,iel)-xxp0)                   &
           +       gy*(xyzcen(2,iel)-xyp0)                   &
           +       gz*(xyzcen(3,iel)-xzp0) )
    enddo
  endif

endif

 2000 format(                                                           &
  ' REINITIALISATION DE LA PRESSION A L''ITERATION ',I10)

!===============================================================================
! 3.  COMMUNICATIONS
!===============================================================================

! Halo synchronization (only variables require this)

if (irangp.ge.0 .or. iperio.eq.1) then

  do f_id = 0, nfld-1

    call field_get_dim(f_id, f_dim)
    call field_get_type(f_id, f_type)

    ! Is the field of type FIELD_VARIABLE?
    if (iand(f_type, FIELD_VARIABLE).eq.FIELD_VARIABLE) then
      if (f_dim.eq.1) then

        call field_get_val_s(f_id, cvar_sca)
        call synsce (cvar_sca)

      else if (f_dim.eq.3) then

        call field_get_val_v(f_id, cvar_vec)
        call synvie(cvar_vec)

      else if (f_dim.eq.6) then

        call field_get_val_v(f_id, cvar_vec)
        call syntis(cvar_vec)
      else
        call csexit(1)
      endif

    endif
  enddo

endif

! ---> Periodicity of rotation for related fields

if (iperio.eq.1) then

  !  -- Reynolds stress tensor

  if (itytur.eq.3 .and. irijco.eq.0) then
    call field_get_val_s(ivarfl(ir11), cvar_r11)
    call field_get_val_s(ivarfl(ir22), cvar_r22)
    call field_get_val_s(ivarfl(ir33), cvar_r33)
    call field_get_val_s(ivarfl(ir12), cvar_r12)
    call field_get_val_s(ivarfl(ir13), cvar_r13)
    call field_get_val_s(ivarfl(ir23), cvar_r23)

    call perrte &
    ( cvar_r11, cvar_r12, cvar_r13,           &
      cvar_r12, cvar_r22, cvar_r23,           &
      cvar_r13, cvar_r23, cvar_r33 )
  endif

  !  -- Note for v2f:
  !     v2 (thus phi) is oriented locally, and is handled as a scalar
  !     regarding periodicity of rotation

endif

!===============================================================================
! 4.  POUR IPHYDR ON DOIT COMMUNIQUER FRCXT AU PREMIER PASSAGE
!     (FRCXT SERT DANS TYPECL)
!     If icalhy=1, rho must be synchronized before copying to previous values
!===============================================================================

if (ipass.eq.1) then

! --- Communication de FRCXT
  if (iphydr.eq.1) then

    if (irangp.ge.0 .or. iperio.eq.1) then
      call synvin (frcxt)
    endif

  endif

! --- Communication de RHO
  if (icalhy.eq.1 .or. idilat.eq.3) then

    call field_get_val_s(icrom, crom)
    if (irangp.ge.0 .or. iperio.eq.1) then
      call synsce (crom)
    endif

  endif

endif

!===============================================================================
! 5. Temporal update of previous values (mass flux, density, ...)
!===============================================================================
!  We exit before SCHTMP as otherwise for 2nd order in time the value of the
!  mass flux at the previous time is overwritten by the value at the current
!  time step. When ntmabs = 0, there is no issue since all mass fluxes are 0.

if (ntmabs.eq.ntpabs .and. isuite.eq.1) return

!  If itrale=0, we are initializing ALE, we do not touch the mas flux either.

if (itrale.gt.0) then
  iappel = 1
  call schtmp(nscal, iappel)
endif

!===============================================================================
! 6.  MISE A JOUR DE LA LOCALISATION DES INTERFACES DE COUPLAGE CS/CS
!===============================================================================

! Localisation des interfaces de couplage via la librairie FVM

! On fait cette mise a jour des interfaces de localisation juste apres
! les changements de geometries dus :
!   - soit a la methode ALE (en fin de pas de temps precedent)
!   - soit a un deplacement impose (cf ci-dessus)

if (nbrcpl.gt.0) call cscloc

!===============================================================================
! 7.  CALCUL DES PROPRIETES PHYSIQUES VARIABLES
!      SOIT VARIABLES AU COURS DU TEMPS
!      SOIT VARIABLES LORS D'UNE REPRISE DE CALCUL
!        (VISCOSITES ET MASSE VOLUMIQUE)
!===============================================================================

if (vcopt_u%iwarni.ge.1) then
  write(nfecra,1010)
endif

iterns = -1
call phyvar(nvar, nscal, iterns, dt)

if (itrale.gt.0) then
  iappel = 2
  call schtmp(nscal, iappel)
endif


! REMPLISSAGE DES COEFS DE PDC
!    ON Y PASSE MEME S'IL N'Y A PAS DE PDC SUR LE PROC COURANT AU CAS OU
!    UN UTILISATEUR DECIDERAIT D'AVOIR UN COEFF DE PDC DEPENDANT DE
!    LA VITESSE MOYENNE OU MAX.

if (ncpdct.gt.0) then

  call cs_head_losses_compute(ckupdc)

 if (iflow .eq.1) then
   call cs_lagr_head_losses(ncepdc, icepdc, itypfb, ckupdc)
 endif

endif

! REMPLISSAGE DES COEFS DE TERME SOURCE DE MASSE

!    ON Y PASSE MEME S'IL N'Y A PAS DE TSM SUR LE PROC COURANT AU CAS OU
!    UN UTILISATEUR DECIDERAIT D'AVOIR UN TSM DEPENDANT DE
!    VALEURS GLOBALES OU MAX.
if(nctsmt.gt.0) then

  !     Mise a zero du tableau de type de TS masse et source
  do ii = 1, ncetsm
    do ivar = 1, nvar
      itypsm(ii,ivar) = 0
      smacel(ii,ivar) = 0.d0
    enddo
  enddo

  iappel = 3
  call  cs_user_mass_source_terms &
( nvar   , nscal  , ncepdc ,                                     &
  ncetsm , iappel ,                                              &
  icepdc ,                                                       &
  icetsm , itypsm , izctsm ,                                     &
  dt     ,                                                       &
  ckupdc , smacel )

  if (ippmod(iaeros).gt.0) then

    allocate(mass_source(ncelet))
    ! Cooling tower model
    ! Evaporation mass exchange term
    call cs_ctwr_bulk_mass_source_term &
      (p0   , molmass_rat, mass_source)

    do ii = 1, ncetsm
      iel = icetsm(ii)
      smacel(ii, ipr) = smacel(ii, ipr) + mass_source(iel)
    enddo

    deallocate(mass_source)
  endif

endif

!------------------------------------------------------------------------
!-- Fill the condensation arrays spcond for the sink term of condensation
!-- and hpcond the thermal exchange coefficient associated to the phase
!-- change (gas phase to liquid phase)
!------------------------------------------------------------------------

if (nftcdt.gt.0) then

  iappel = 3

  ! Condensation source terms arrays initialized
  do ii = 1, nfbpcd
    do ivar = 1, nvar
      itypcd(ii,ivar) = 0
      spcond(ii,ivar) = 0.d0
      hpcond(ii)      = 0.d0
    enddo
  enddo

  call cs_user_boundary_mass_source_terms &
( nvar   , nscal  ,                                              &
  nfbpcd , iappel ,                                              &
  ifbpcd , itypcd , izftcd ,                                     &
  spcond , tpar)

  if (nzones.eq.1) then
    do ii = 1, nfbpcd
      iz = izzftcd(ii)
      zrob  (iz) = rob
      zcondb(iz) = condb
      zcpb  (iz) = cpb
      zhext (iz) = hext
      ztext (iz) = text
      ztpar (iz) = tpar
    enddo
  endif

  ! Empiric laws used by COPAIN condensation model to
  ! the computation of the condensation source term and
  ! exchange coefficient of the heat transfer imposed
  ! as boundary condition.

  call condensation_copain_model &
( nvar   , nfbpcd , ifbpcd , izzftcd ,                           &
  tpar   ,                                                       &
  spcond , hpcond )

endif

!----------------------------------------------------------
!-- Fill the condensation arrays (svcond) for the sink term
!-- of condensation and source term type (itypst) of each
!-- variable solved associated to the metal structures
!-- condensation modelling.
!----------------------------------------------------------

if (icondv.eq.0) then

  !-- Condensation source terms arrays initialized
  do iel = 1, ncelet
    ltmast(iel) = 0
    do ivar = 1, nvar
      itypst(iel, ivar) = 0
      svcond(iel, ivar) = 0.d0
    enddo
    flxmst(iel) = 0.d0
  enddo

  call cs_user_metal_structures_source_terms &
( nvar   , nscal  ,                                              &
  ncmast , ltmast,                                               &
  itypst , izmast ,                                              &
  svcond , tmet)

  ! Condensation model to compute the sink source term
  ! (svcond) and the  heat transfer flux (flxmst) imposed
  ! in the cells associated to the metal  structures
  ! volume where this phenomenon occurs.

  call metal_structures_copain_model &
( ncmast , ltmast ,                                          &
  tmet   ,                                                   &
  svcond(:, ipr)  , flxmst )

  ! array initialization if the metal structures
  ! condensation model is coupled with
  ! a 0-D thermal model
  ! FIXME add restart file later
  if (itagms.eq.1) then
    do icmst = 1, ncmast
      iel = ltmast(icmst)
      t_metal(iel,1) = tmet0
      t_metal(iel,2) = tmet0
    enddo
  endif

endif

!===============================================================================
! 7.bis Current to previous for variables and GWF module
!===============================================================================

do f_id = 0, nfld - 1
  call field_get_type(f_id, f_type)
  ! Is the field of type FIELD_VARIABLE?
  if (iand(f_type, FIELD_VARIABLE).eq.FIELD_VARIABLE) then
    call field_current_to_previous(f_id)

    ! For buoyant scalar with source termes, current to previous for them
    call field_get_key_int(f_id, kst, st_id)
    if (st_id .ge.0) then
      call field_current_to_previous(st_id)
    endif
  endif
enddo



if (ippmod(idarcy).eq.1) then

  ! Index of the corresponding field
  call field_get_val_prev_s_by_name('capacity', cproa_capacity)
  call field_get_val_prev_s_by_name('saturation', cproa_sat)
  call field_get_val_s_by_name('capacity', cpro_capacity)
  call field_get_val_s_by_name('saturation', cpro_sat)

  do iel = 1, ncel
    cproa_capacity(iel) = cpro_capacity(iel)
    cproa_sat(iel) = cpro_sat(iel)
  enddo

  do ii = 1, nscal
    ivar = ivarfl(isca(ii))
    call field_get_key_struct_gwf_soilwater_partition(ivarfl(isca(ii)), &
                                                      sorption_scal)
    call field_get_val_s(sorption_scal%idel, cpro_delay)
    call field_get_val_prev_s(sorption_scal%idel, cproa_delay)
    do iel = 1, ncel
      cproa_delay(iel) = cpro_delay(iel)
    enddo
  enddo

endif

!===============================================================================
! 8.  CALCUL DU NOMBRE DE COURANT ET DE FOURIER
!     CALCUL DU PAS DE TEMPS SI VARIABLE
!===============================================================================

if(vcopt_u%iwarni.ge.1) then
  write(nfecra,1020)
endif

call dttvar &
 ( nvar   , nscal  , ncepdc , ncetsm ,                            &
   vcopt_u%iwarni   ,                                              &
   icepdc , icetsm , itypsm ,                                     &
   dt     ,                                                       &
   ckupdc , smacel )

if (nbaste.gt.0.and.itrale.gt.nalinf) then
  ntrela = ntcabs - ntpabs
  call astpdt(dt, ncelet, ntrela)
endif

! Compute the pseudo tensorial time step if needed for the pressure solving

if (iand(vcopt_p%idften, ANISOTROPIC_DIFFUSION).ne.0               &
    .and.(ippmod(idarcy).eq.-1)) then

  call field_get_val_v(idtten, dttens)

  do iel = 1, ncel
    dttens(1, iel) = dt(iel)
    dttens(2, iel) = dt(iel)
    dttens(3, iel) = dt(iel)
    dttens(4, iel) = 0.d0
    dttens(5, iel) = 0.d0
    dttens(6, iel) = 0.d0
  enddo

  do ielpdc = 1, ncepdc
    iel = icepdc(ielpdc)

    ! dttens = (1/dt + Kpdc)^-1
    hdls(1) = ckupdc(1, ielpdc) + 1.d0/dt(iel)
    hdls(2) = ckupdc(2, ielpdc) + 1.d0/dt(iel)
    hdls(3) = ckupdc(3, ielpdc) + 1.d0/dt(iel)
    hdls(4) = ckupdc(4, ielpdc)
    hdls(5) = ckupdc(5, ielpdc)
    hdls(6) = ckupdc(6, ielpdc)

    call symmetric_matrix_inverse(hdls, dttens(:, iel))
  enddo

  if (irangp.ge.0) then
    call syntis(dttens)
  endif

endif

!===============================================================================
!     RECALAGE DE LA PRESSION Pth ET MASSE VOLUMIQUE rho
!     POUR L'ALGORITHME A MASSE VOLUMIQUE VARIABLE.
!===============================================================================

if (idilat.eq.3.or.ipthrm.eq.1) then
  call pthrbm &
 ( nvar   , ncetsm , nfbpcd , ncmast,                             &
   dt     , smacel , spcond , svcond )

endif

!===============================================================================
! 9.  CHARGEMENT ET TRADUCTION DES CONDITIONS AUX LIMITES
!===============================================================================

if(vcopt_u%iwarni.ge.1) then
  write(nfecra,1030)
endif

!  -- Methode des vortex en LES :
!     Definition ou modification eventuelle des parametres
!     Mise a jour des vortex

if (ivrtex.eq.1) then

  iappel = 2
  call usvort &
 ( nvar   , nscal  ,                                              &
   iappel ,                                                       &
   dt     )

!     Verification des donnees entrees par l'utilisateur
!       (au premier passage seulement)
  if (ipass.eq.1) then
    call vorver ( nfabor , iappel )
  endif

  if(irangp.le.0) then
    call vortex
  endif

! -- Fin de zone Methode des vortex

endif

! --- Methode ALE : debut de boucle d'implicitation du deplacement des
!       structures. ITRFIN=0 indique qu'on a besoin de refaire une iteration
!       pour Syrthes, T1D ou rayonnement.
italim = 1
itrfin = 1
ineefl = 0
if (iale.ge.1 .and. nalimx.gt.1 .and. itrale.gt.nalinf) then
!     On reserve certains tableaux pour permettre le retour a l'etat
!       initial en fin d'iteration ALE
!       - flux de masse
!       - conditions aux limites de gradient de P et U (car on a un appel
!         a GDRCEL pour les non orthogonalites pour calculer les CL reelles)
!         -> n'est peut-etre pas reellement necessaire
!       - la pression initiale (car RTPA est aussi ecrase dans le cas
!         ou NTERUP>1) -> on pourrait optimiser en ne reservant que si
!         necessaire ...
  allocate(flmalf(nfac))
  allocate(flmalb(nfabor))
  allocate(cofale(nfabor,11))
  allocate(xprale(ncelet))
  ineefl = 1

  if (nbccou.gt.0 .or. nfpt1t.gt.0 .or. iirayo.gt.0) itrfin = 0

else
  flmalf => null()
  flmalb => null()
  cofale => null()
  xprale => null()
endif

300 continue


! --- Boucle sur navstv pour couplage vitesse/pression
!     on s'arrete a NTERUP ou quand on a converge
!     ITRFUP=0 indique qu'on a besoin de refaire une iteration
!     pour Syrthes, T1D ou rayonnement.
itrfup = 1

if (nterup.gt.1) then
  allocate(trava(ndim,ncelet))
else
  trava => rvoid2
endif

if (nterup.gt.1.or.isno2t.gt.0) then
  if (nbccou.gt.0 .or. nfpt1t.gt.0 .or. iirayo.gt.0) itrfup = 0
endif

icodcl => null()
rcodcl => null()
hbord => null()
theipb => null()
visvdr => null()

! Allocate temporary arrays for boundary conditions
if (italim .eq. 1) then
  allocate(icodcl(nfabor,nvar))
  allocate(rcodcl(nfabor,nvar,3))
endif
if (isvhb.gt.0) then
  allocate(hbord(nfabor))
endif
! Boundary value of the thermal scalar in I'
if (iscalt.gt.0) then
  allocate(theipb(nfabor))
endif
if (itytur.eq.4 .and. idries.eq.1) then
  allocate(visvdr(ncelet))
endif

icvrge = 0
inslst = 0

! Darcy : in case of a steady flow, we resolve Richards only once,
! at the first time step.
if (ippmod(idarcy).eq.1) then
  if ((darcy_unsteady.eq.0).and.(ntcabs.gt.1)) goto 100
endif

iterns = 1
do while (iterns.le.nterup)

  ! Calls user BCs and computes BC coefficients
  call condli &
    (nvar   , nscal  , iterns ,                                    &
     isvhb  ,                                                      &
     itrale , italim , itrfin , ineefl , itrfup ,                  &
     flmalf , flmalb , cofale , xprale ,                           &
     icodcl , isostd ,                                             &
     dt     , rcodcl ,                                             &
     visvdr , hbord  , theipb )

  if (nftcdt.gt.0) then
    ! Coefficient exchange of the enthalpy scalar
    ivar = isca(iscalt)
    call field_get_coefa_s(ivarfl(ivar) , coefap)
    call field_get_coefaf_s(ivarfl(ivar), cofafp)
    call field_get_coefbf_s(ivarfl(ivar), cofbfp)

    ! Pass the heat transfer computed by the Empiric laws
    ! of the COPAIN condensation to impose the heat transfer
    ! at the wall due to condensation for the enthalpy scalar.
    do ii = 1, nfbpcd

      ifac= ifbpcd(ii)
      iel = ifabor(ifac)

      ! Enthalpy Boundary condition associated
      ! to the heat transfer due to condensation.
      cofafp(ifac) = -hpcond(ii)*coefap(ifac)
      cofbfp(ifac) =  hpcond(ii)

    enddo

  endif

!     ==============================================
!     Appel de l'interface sol-atmosphere
!     ==============================================

  if (ippmod(iatmos).eq.2.and.iatsoil.eq.1.and.nfmodsol.gt.0) then
    call field_get_val_s(icrom, crom)
    call field_get_val_s(ivarfl(isca(iscalt)), cvar_scalt)
    call field_get_val_s(ivarfl(isca(itotwt)), cvar_totwt)
    call solvar(cvar_scalt , cvar_totwt ,                &
                crom   , dt ,                            &
                rcodcl )
  endif

  !     UNE FOIS LES COEFFICIENTS CALCULES, ON PEUT EN DEDUIRE PLUS
  !     FACILEMENT (I.E. SANS RECALCULS INUTILES) LES TERMES A
  !     ENVOYER POUR LES COUPLAGES AUX BORDS (TYPE SYRTHES)

  ! En compressible et si on couple ave l'energie
  ! on recupere le Cv de la phase couplee

  if (itherm .eq. 3) then

    if(icv.ge.0) then
      cvcst = 0.d0
    else
      cvcst = cv0
    endif
  else
    cvcst = 0.d0
  endif

  ! On envoie le tout vers SYRTHES, en distinguant CP
  !  constant ou variable
  if (itrfin.eq.1 .and. itrfup.eq.1) then

    if (isvhb.gt.0) then
      call coupbo(itherm, cvcst, hbord, theipb)
    endif

    if (iscalt.gt.0 .and. nfpt1t.gt.0) then
      call cou1do(cvcst, hbord, theipb)

      if (iirayo.ge.1) call cou1di(nfabor, iscalt, icodcl, rcodcl)

    endif

    ! 1-D thermal model coupling with condensation
    ! on a surface region
    if (nftcdt.gt.0.and.nztag1d.eq.1) then
      call cs_tagmro &
     ( nfbpcd , ifbpcd , izzftcd ,                  &
       dt     )
    endif

     ! 0-D thermal model coupling with condensation
     ! on a volume region associated to metal structures
    if (icondv.eq.0.and.itagms.eq.1) then
      call cs_metal_structures_tag &
     ( ncmast , ltmast ,                          &
       dt     )
    endif

  endif

  !     ON N'A PLUS BESOIN DE ISVHB OU ISVHT (POUR HBORD ET TBORD)
  !     A PARTIR D'ICI

  ! Compute wall distance
  ! TODO it has to be moved before phyvar, for that bc types have to be known
  ! (itypfb)

  !       (Nouvel algorithme. L'ancien est dans condli)
  ! In ALE, this computation is done only for the first step.
  if (italim.eq.1) then

    ! Pour le moment, on suppose que l'on peut se contenter de faire
    !  cela au premier passage, sauf avec les maillages mobiles. Attention donc
    !  aux conditions aux limites variables (faces qui deviennent des parois ou
    !  parois qui deviennent autre chose)

    ! Wall distance is computed if:
    !   - it has to be updated
    !   - we need it
    ! In case there is no wall, distance is a big value.
    if (imajdy.eq.0 .and. ineedy.eq.1) then

      if (abs(icdpar).eq.1) then
        call distpr(itypfb)
      ! Deprecated algorithm
      else if (abs(icdpar).eq.2) then
        call distpr2(itypfb)
      endif
      ! Wall distance is not updated except if ALE is switched on
      if (iale.eq.0) imajdy = 1
    endif

  endif

  ! Compute y+ if needed
  ! and Van Driest "amortissement"
  if (itytur.eq.4 .and. idries.eq.1) then
    call distyp(itypfb, visvdr)
  endif

!===============================================================================
! 10. DANS LE CAS  "zero pas de temps" EN "NON SUITE" DE CALCUL
!      ON SORT ICI
!===============================================================================

  if (ntmabs.eq.ntpabs .and. isuite.eq.0) must_return = .true.

  if (iilagr.eq.3) must_return = .true.

!===============================================================================
! 11. RESOLUTION DE LA VITESSE DE MAILLAGE EN ALE
!===============================================================================

  if (iale.ge.1) then

    ! otherwise it is done in navstv.f90
    if (itrale.eq.0) then

      call cs_ale_solve_mesh_velocity(iterns, impale, ialtyb)
      must_return = .true.

    endif

  endif

!===============================================================================
! Return now if required
!===============================================================================

  if (must_return) then

    if (associated(hbord)) deallocate(hbord)
    if (associated(theipb)) deallocate(theipb)
    if (associated(visvdr)) deallocate(visvdr)

    if (nterup.gt.1) then
      deallocate(trava)
    endif

    deallocate(icodcl, rcodcl)

    deallocate(isostd)

    return

  endif

!===============================================================================

!===============================================================================
! 11. CALCUL A CHAMP DE VITESSE NON FIGE :
!      ON RESOUT VITESSE ET TURBULENCE
!    ON SUPPOSE QUE TOUTES LES PHASES SONT FIGEES OU AUCUNE
!===============================================================================

  ! En cas de champ de vitesse fige, on ne boucle pas sur U/P
  if (iccvfg.eq.0) then

!===============================================================================
! 12. Solve momentum and mass equation
!===============================================================================

    ! In case of buoyancy, scalars and momentum are coupled
    if (n_buoyant_scal.ge.1) then

      if(vcopt_u%iwarni.ge.1) then
        write(nfecra,1060)
      endif

      ! Update buoyant scalar(s)
      call scalai(nvar, nscal , iterns , dt)

      ! Diffusion terms for weakly compressible algorithm
      if (idilat.ge.4) then
        call diffst(nscal, iterns)
      endif

      ! Update the density
      !-------------------

      call phyvar(nvar, nscal, iterns, dt)

    endif

    if (vcopt_u%iwarni.ge.1) then
      write(nfecra,1040) iterns
    endif

    ! Coupled solving of the velocity components

    if (ippmod(idarcy).eq.-1) then

      call navstv &
      ( nvar   , nscal  , iterns , icvrge , itrale ,                   &
        isostd ,                                                       &
        dt     ,                                                       &
        frcxt  ,                                                       &
        trava  )

      ! Update local pointer arrays for transient turbomachinery computations
      if (iturbo.eq.2) then
        call field_get_val_s(ivarfl(ipr), cvar_pr)
      endif

    else

      call richards (icvrge, dt)

      if (iihmpr.eq.1) then
        call uidapp                                                    &
         ( darcy_anisotropic_permeability,                             &
           darcy_anisotropic_dispersion,                               &
           darcy_gravity,                                              &
           darcy_gravity_x, darcy_gravity_y, darcy_gravity_z,          &
           darcy_unsaturated)
      endif

      ! Darcy : update data specific to underground flow
      mbrom = 0
      call usphyv(nvar, nscal, mbrom, dt)

      if (darcy_unsteady.eq.0) then

        do iel = 1, ncel
          cproa_capacity(iel) = cpro_capacity(iel)
          cproa_sat(iel) = cpro_sat(iel)
        enddo

        do ii = 1, nscal
          call field_get_key_struct_gwf_soilwater_partition(ivarfl(isca(ii)), &
                                                            sorption_scal)
          call field_get_val_s(sorption_scal%idel, cpro_delay)
          call field_get_val_prev_s(sorption_scal%idel, cproa_delay)
          do iel = 1, ncel
            cproa_delay(iel) = cpro_delay(iel)
          enddo
        enddo

      endif

    endif

    !     Si c'est la derniere iteration : INSLST = 1
    if ((icvrge.eq.1).or.(iterns.eq.nterup)) then

      ! Si on a besoin de refaire une nouvelle iteration pour SYRTHES,
      ! rayonnement, paroi thermique 1D...
      ! ET que l'on est a la derniere iteration en ALE !

      ! ...alors, on remet a zero les indicateurs de convergence
      if (itrfup.eq.0.and.itrfin.eq.1) then
        itrfup = 1
        icvrge = 0
        iterns = iterns - 1

        ! ...sinon, on termine
      else
        inslst = 1
      endif

      ! For explicit mass flux
      if (istmpf.eq.0.and.inslst.eq.0) then
        iappel = 3
        call schtmp(nscal, iappel)
      endif

      if (inslst.eq.1) goto 100

    endif

  endif ! Fin si calcul sur champ de vitesse figee

  iterns = iterns + 1

enddo

100 continue

! DARCY : the hydraulic head, identified with the pressure,
! has been updated by the call to Richards.
! As diffusion of scalars depends on hydraulic head in the
! general case, in order to compute the exact
! values of the boundary faces coefficients, we have to
! call boundary conditions routine again.
! Moreover, we need an update of the boundary
! conditions in the cases where they vary in time.

if (ippmod(idarcy).eq.1) then

  ! Calls user BCs and computes BC coefficients
  call condli &
    (nvar   , nscal  , iterns ,                                    &
     isvhb  ,                                                      &
     itrale , italim , itrfin , ineefl , itrfup ,                  &
     flmalf , flmalb , cofale , xprale ,                           &
     icodcl , isostd ,                                             &
     dt     , rcodcl ,                                             &
     visvdr , hbord  , theipb )

endif

! Free memory
if (associated(hbord)) deallocate(hbord)
if (associated(theipb)) deallocate(theipb)
if (associated(visvdr)) deallocate(visvdr)

if (nterup.gt.1) then
  deallocate(trava)
endif

! Calcul sur champ de vitesse NON fige SUITE (a cause de la boucle U/P)
if (iccvfg.eq.0) then

!===============================================================================
! 13.  DEPLACEMENT DES STRUCTURES EN ALE ET TEST DE BOUCLAGE IMPLICITE
!===============================================================================

  if (nbstru.gt.0.or.nbaste.gt.0) then

    call strdep &
  ( itrale , italim , itrfin ,                                     &
    nvar   ,                                                       &
    dt     ,                                                       &
    flmalf , flmalb , cofale , xprale )

    !     On boucle eventuellement sur de deplacement des structures
    if (itrfin.ne.-1) then
      italim = italim + 1
      goto 300
    endif

    ! Free memory
    if (associated(flmalf)) then
      deallocate(flmalf, flmalb)
      deallocate(cofale)
      deallocate(xprale)
    endif

  endif

  !     On ne passe dans SCHTMP que si ISTMPF.EQ.0 (explicite)
  if (istmpf.eq.0) then
    iappel = 4
    call schtmp(nscal, iappel)
  endif

!===============================================================================
! 14. RESOLUTION TURBULENCE
!===============================================================================

  iok = 0
  if(vcopt_u%iwarni.ge.1) then
    if( itytur.eq.2 .or. itytur.eq.3              &
         .or. itytur.eq.5 .or. iturb.eq.60 ) then
      iok = 1
    endif
    if(iok.eq.1) then
      write(nfecra,1050)
    endif
  endif

  ! Si on est en v2f (phi-fbar ou BL-v2/k), on reserve un tableau
  ! de taille NCELET pour eviter de recalculer la production dans RESV2F
  if (itytur.eq.5) then
    allocate(prdv2f(ncelet))
  endif

  if ((itytur.eq.2) .or. (itytur.eq.5)) then

    call turbke &
  ( nvar   , nscal  ,                                              &
    ncepdc , ncetsm ,                                              &
    icepdc , icetsm , itypsm ,                                     &
    dt     ,                                                       &
    tslagr ,                                                       &
    ckupdc , smacel ,                                              &
    prdv2f )

    if( itytur.eq.5 )  then

      call resv2f &
    ( nvar   , nscal  ,                                              &
      ncepdc , ncetsm ,                                              &
      icepdc , icetsm , itypsm ,                                     &
      dt     ,                                                       &
      ckupdc , smacel ,                                              &
      prdv2f )

      ! Free memory
      deallocate(prdv2f)

    endif

    call field_get_val_s(ivarfl(ik), cvar_k)
    call field_get_val_prev_s(ivarfl(ik), cvara_k)
    call field_get_val_s(ivarfl(iep), cvar_ep)
    call field_get_val_prev_s(ivarfl(iep), cvara_ep)

    !  RELAXATION DE K ET EPSILON SI IKECOU=0 EN INSTATIONNAIRE
    if (ikecou.eq.0 .and. idtvar.ge.0) then
      call field_get_key_struct_var_cal_opt(ivarfl(ik), vcopt)
      relaxk = vcopt%relaxv
      call field_get_key_struct_var_cal_opt(ivarfl(iep), vcopt)
      relaxe = vcopt%relaxv
      do iel = 1,ncel
        cvar_k(iel) = relaxk*cvar_k(iel) + (1.d0-relaxk)*cvara_k(iel)
        cvar_ep(iel) = relaxe*cvar_ep(iel) + (1.d0-relaxe)*cvara_ep(iel)
      enddo
    endif

  else if(itytur.eq.3) then

    ! Compute Alpha for EBRSM
    if (iturb.eq.32) then

      call resalp(ivarfl(ial), xcl)

    endif

    call turrij &
  ( nvar   , nscal  ,                                              &
    ncepdc , ncetsm ,                                              &
    icepdc , icetsm , itypsm ,                                     &
    dt     ,                                                       &
    tslagr ,                                                       &
    ckupdc , smacel )

  else if (iturb.eq.60) then

    call turbkw &
  ( nvar   , nscal  ,                                              &
    ncepdc , ncetsm ,                                              &
    icepdc , icetsm , itypsm ,                                     &
    dt     ,                                                       &
    tslagr ,                                                       &
    ckupdc , smacel )

    call field_get_val_s(ivarfl(ik), cvar_k)
    call field_get_val_prev_s(ivarfl(ik), cvara_k)
    call field_get_val_s(ivarfl(iomg), cvar_omg)
    call field_get_val_prev_s(ivarfl(iomg), cvara_omg)

    !  RELAXATION DE K ET OMEGA SI IKECOU=0
    if (ikecou.eq.0 .and. idtvar.ge.0) then
      call field_get_key_struct_var_cal_opt(ivarfl(ik), vcopt)
      relaxk = vcopt%relaxv
      call field_get_key_struct_var_cal_opt(ivarfl(iomg), vcopt)
      relaxw = vcopt%relaxv
      do iel = 1,ncel
        cvar_k(iel)   = relaxk*cvar_k(iel)   + (1.d0-relaxk)*cvara_k(iel)
        cvar_omg(iel) = relaxw*cvar_omg(iel) + (1.d0-relaxw)*cvara_omg(iel)
      enddo
    end if

  else if( iturb.eq.70 ) then

    call turbsa &
  ( nvar   , nscal  ,                                              &
    ncepdc , ncetsm ,                                              &
    icepdc , icetsm , itypsm ,                                     &
    dt     ,                                                       &
    ckupdc , smacel ,                                              &
    itypfb )

    call field_get_val_s(ivarfl(inusa), cvar_nusa)
    call field_get_val_prev_s(ivarfl(inusa), cvara_nusa)

    !  RELAXATION DE NUSA
    if (idtvar.ge.0) then
      call field_get_key_struct_var_cal_opt(ivarfl(inusa), vcopt)
      relaxn = vcopt%relaxv
      do iel = 1,ncel
        cvar_nusa(iel) = relaxn*cvar_nusa(iel)+(1.d0-relaxn)*cvara_nusa(iel)
      enddo
    endif

  endif

endif  ! Fin si calcul sur champ de vitesse fige SUITE


!     Ici on peut liberer les eventuels tableaux SKW et DIVUKW

!===============================================================================
! 15.  RESOLUTION DES SCALAIRES
!===============================================================================

if (nscal.ge.1 .and. iirayo.gt.0) then

  if (vcopt_u%iwarni.ge.1 .and. mod(ntcabs,nfreqr).eq.0) then
    write(nfecra,1070)
  endif

  call cs_rad_transfer_solve(itypfb, nclacp, nclafu, &
                             dt, cp2fol, cp2ch, ichcor)
endif

if (nscal.ge.1) then

  if(vcopt_u%iwarni.ge.1) then
    write(nfecra,1060)
  endif

  ! Update non-buoyant scalar(s)
  iterns = -1
  call scalai(nvar, nscal, iterns, dt)

  ! Diffusion terms for weakly compressible algorithm
  if (idilat.ge.4) then
    call diffst(nscal, iterns)
  endif

endif

! Free memory
deallocate(icodcl, rcodcl)

deallocate(isostd)

!===============================================================================
! 16.  TRAITEMENT DU FLUX DE MASSE, DE LA VISCOSITE,
!      DE LA MASSE VOLUMIQUE ET DE LA CHALEUR SPECIFIQUE POUR
!      UN THETA SCHEMA
!===============================================================================


iappel = 5
call schtmp(nscal, iappel)

!===============================================================================
! Update flow through fans
!===============================================================================

n_fans = cs_fan_n_fans()
if (n_fans .gt. 0) then
  call field_get_key_int(ivarfl(iu), kimasf, iflmas)
  call field_get_key_int(ivarfl(iu), kbmasf, iflmab)
  call field_get_val_s(iflmas, i_mass_flux)
  call field_get_val_s(iflmab, b_mass_flux)
  call field_get_val_s(icrom, crom)
  call field_get_val_s(ibrom, brom)
  call debvtl(i_mass_flux, b_mass_flux, crom, brom)
endif

!===============================================================================

!--------
! Formats
!--------

#if defined(_CS_LANG_FR)

 1000 format(/,                                                   &
' ------------------------------------------------------------',/,&
                                                              /,/,&
'  INITIALISATIONS                                            ',/,&
'  ===============                                            ',/)
 1010 format(/,                                                   &
' ------------------------------------------------------------',/,&
                                                              /,/,&
'  CALCUL DES GRANDEURS PHYSIQUES                             ',/,&
'  ==============================                             ',/)
 1020 format(/,                                                   &
' ------------------------------------------------------------',/,&
                                                              /,/,&
'  CALCUL DU CFL, DU FOURIER ET DU DT VARIABLE                ',/,&
'  ===========================================                ',/)
 1030 format(/,                                                   &
' ------------------------------------------------------------',/,&
                                                              /,/,&
'  MISE EN PLACE DES CONDITIONS AUX LIMITES                   ',/,&
'  ========================================                   ',/)
 1040 format(/,                                                   &
' ------------------------------------------------------------',/,&
                                                              /,/,&
'  RESOLUTION DES EQUATIONS DE NAVIER-STOKES (sub iter: ',i3,')',/,&
'  =========================================                  ',/)
 1050 format(/,                                                   &
' ------------------------------------------------------------',/,&
                                                              /,/,&
'  RESOLUTION DES EQUATIONS DES VARIABLES TURBULENTES         ',/,&
'  ==================================================         ',/)
 1060 format(/,                                                   &
' ------------------------------------------------------------',/,&
                                                              /,/,&
'  RESOLUTION DES EQUATIONS SUR L''ENERGIE ET LES SCALAIRES   ',/,&
'  ========================================================   ',/)
 1070 format(/,                                                   &
 '------------------------------------------------------------',/,&
                                                              /,/,&
 ' RESOLUTION DES TRANSFERTS THERMIQUES RADIATIFS             ',/,&
'  ==============================================             ',/)

#else

 1000 format(/,                                                   &
' ------------------------------------------------------------',/,&
                                                              /,/,&
'  INITIALISATIONS                                            ',/,&
'  ===============                                            ',/)
 1010 format(/,                                                   &
' ------------------------------------------------------------',/,&
                                                              /,/,&
'  COMPUTATION OF PHYSICAL QUANTITIES                         ',/,&
'  ==================================                         ',/)
 1020 format(/,                                                   &
' ------------------------------------------------------------',/,&
                                                              /,/,&
'  COMPUTATION OF CFL, FOURIER AND VARIABLE DT                ',/,&
'  ===========================================                ',/)
 1030 format(/,                                                   &
' ------------------------------------------------------------',/,&
                                                              /,/,&
'  SETTING UP THE BOUNDARY CONDITIONS                         ',/,&
'  ==================================                         ',/)
 1040 format(/,                                                   &
' ------------------------------------------------------------',/,&
                                                              /,/,&
'  SOLVING NAVIER-STOKES EQUATIONS (sub iter: ',i3,')'         ,/,&
'  ===============================                             ',/)
 1050 format(/,                                                   &
' ------------------------------------------------------------',/,&
                                                              /,/,&
'  SOLVING TURBULENT VARIABLES EQUATIONS                      ',/,&
'  =====================================                      ',/)
 1060 format(/,                                                   &
' ------------------------------------------------------------',/,&
                                                              /,/,&
'  SOLVING ENERGY AND SCALARS EQUATIONS                       ',/,&
'  ====================================                       ',/)
 1070 format(/,                                                   &
 '------------------------------------------------------------',/,&
                                                              /,/,&
 ' SOLVING THERMAL RADIATIVE TRANSFER                         ',/,&
'  ==================================                         ',/)

#endif

!----
! End
!----

end subroutine
