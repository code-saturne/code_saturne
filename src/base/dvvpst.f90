!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2012 EDF S.A.
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

subroutine dvvpst &
!================

 ( nummai , numtyp ,                                              &
   nvar   , nscal  , nvlsta , nvisbr ,                            &
   ncelps , nfacps , nfbrps ,                                     &
   itypps ,                                                       &
   lstcel , lstfac , lstfbr ,                                     &
   dt     , rtpa   , rtp    , propce , propfa , propfb ,          &
   coefa  , coefb  , statce , stativ , statfb ,                   &
   tracel , trafac , trafbr )

!===============================================================================
! Purpose:
! --------

! Standard output of variables on post-processing meshes
!   (called after cs_user_extra_operations)
!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nummai           ! ec ! <-- ! numero du maillage post                        !
! numtyp           ! ec ! <-- ! numero de type de post-traitement              !
!                  !    !     ! (-1: volume, -2: bord, nummai par defaut)      !
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! nvlsta           ! e  ! <-- ! nombre de variables stat. lagrangien           !
! nvisbr           ! e  ! <-- ! nombre de statistiques aux frontieres          !
! ncelps           ! e  ! <-- ! nombre de cellules du maillage post            !
! nfacps           ! e  ! <-- ! nombre de faces interieur post                 !
! nfbrps           ! e  ! <-- ! nombre de faces de bord post                   !
! itypps(3)        ! te ! <-- ! indicateur de presence (0 ou 1) de             !
!                  !    !     ! cellules (1), faces (2), ou faces de           !
!                  !    !     ! de bord (3) dans le maillage post              !
! lstcel(ncelps    ! te ! <-- ! liste des cellules du maillage post            !
! lstfac(nfacps    ! te ! <-- ! liste des faces interieures post               !
! lstfbr(nfbrps    ! te ! <-- ! liste des faces de bord post                   !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
! statce           ! tr ! <-- ! statistiques cellules (lagrangien)             !
!  (ncelet,nvlsta) !    !     !                                                !
! stativ           ! tr ! <-- ! statistiques variance cellules (lagrangien)    !
!  (ncelet,nvlsta) !    !     !                                                !
! statfb           ! tr ! <-- ! statistiques faces bord (lagrangien)           !
!  (nfabor,nvisbr) !    !     !                                                !
! tracel(*)        ! tr ! <-- ! tab reel valeurs cellules post                 !
! trafac(*)        ! tr ! <-- ! tab reel valeurs faces int. post               !
! trafbr(*)        ! tr ! <-- ! tab reel valeurs faces bord post               !
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens, only: ndimfb
use pointe
use entsor
use cstnum
use cstphy
use optcal
use numvar
use parall
use period
use lagpar
use lagran
use ppppar
use ppthch
use ppincl
use radiat
use cplsat
use mesh

!===============================================================================

implicit none

! Arguments

integer          nummai , numtyp
integer          nvar   , nscal  , nvlsta , nvisbr
integer          ncelps , nfacps , nfbrps

integer          itypps(3)
integer          lstcel(ncelps), lstfac(nfacps), lstfbr(nfbrps)

double precision dt(ncelet), rtpa(ncelet,*), rtp(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(ndimfb,*)
double precision coefa(ndimfb,*), coefb(ndimfb,*)
double precision statce(ncelet,nvlsta), statfb(nfabor,nvisbr)
double precision stativ(ncelet,nvlsta)
double precision tracel(ncelps*3)
double precision trafac(nfacps*3), trafbr(nfbrps*3)

! Local variables

character*32     namevr, namev1, namev2
character*80     name80

integer          inc   , iccocg, nswrgp, imligp, iwarnp
integer          isorva, isaut
integer          ifac  , iloc  , ivar , iclvar, iclvaf
integer          ira   , idivdt, ineeyp
integer          ipp   , idimt , ii    , kk   , iel
integer          ivarl , iip
integer          iii, ivarl1 , ivarlm , iflu   , ilpd1  , icla
integer          iscal , ipcvsl, ipcvst, iflmab
integer          ientla, ivarpr
integer          ipccp , ipcrom

double precision xcp   , xvsl  , srfbn
double precision visct , flumab, diipbx, diipby, diipbz
double precision epsrgp, climgp, extrap
double precision pcentr

double precision rbid(1)

double precision, allocatable, dimension(:,:) :: grad
double precision, allocatable, dimension(:) :: treco

!===============================================================================

! Initialize variables to avoid compiler warnings

ipp = 0

!===============================================================================
! 1.1. Fluid domain
!===============================================================================

if (numtyp .eq. -1) then

  !  1.1.2 Variables supplementaires automatiques
  !  --------------------------------------------

  ! Distance a la paroi (si LES+VanDriest ou Rij+Echo ou K-w SST)

  if (ineedy.eq.1 .and. abs(icdpar).eq.1) then

    namevr = 'DistWall'
    idimt = 1
    ientla = 0
    ivarpr = 1

    call psteva(nummai, namevr, idimt, ientla, ivarpr,  &
    !==========
                ntcabs, ttcabs, dispar, rbid, rbid)

  endif

  ! Yplus (si LES+VanDriest)

  if (ineedy.eq.1 .and. abs(icdpar).eq.1) then

    ineeyp = 0
    if (itytur.eq.4.and.idries.eq.1) then
      ineeyp = 1
    endif

    if (ineeyp.eq.1) then

      namevr = 'Yplus'
      idimt = 1
      ientla = 0
      ivarpr = 1

      call psteva(nummai, namevr, idimt, ientla, ivarpr,          &
      !==========
                  ntcabs, ttcabs, yplpar, rbid, rbid)

    endif

  endif

  ! Vitesse et pression absolues en cas de calcul en repère relatif

  if (icorio.eq.1) then

    ipcrom = ipproc(irom)

    namevr = 'Pressure'
    idimt = 1
    ientla = 0
    ivarpr = 0

    do iloc = 1, ncelps

      iel = lstcel(iloc)

      pcentr =   0.5d0*((omegay*xyzcen(3,iel) - omegaz*xyzcen(2,iel))**2 &
                      + (omegaz*xyzcen(1,iel) - omegax*xyzcen(3,iel))**2 &
                      + (omegax*xyzcen(2,iel) - omegay*xyzcen(1,iel))**2)

      tracel(iloc) = rtp(iel,ipr) + propce(iel,ipcrom)*pcentr

    enddo

    call psteva(nummai, namevr, idimt, ientla, ivarpr,        &
    !==========
                ntcabs, ttcabs, tracel, rbid, rbid)


    namevr = 'Velocity'
    idimt = 3
    ientla = 1
    ivarpr = 0

    do iloc = 1, ncelps

      iel = lstcel(iloc)

      tracel(1 + (iloc-1)*idimt) = rtp(iel,iu) &
          + (omegay*xyzcen(3,iel) - omegaz*xyzcen(2,iel))

      tracel(2 + (iloc-1)*idimt) = rtp(iel,iv) &
          + (omegaz*xyzcen(1,iel) - omegax*xyzcen(3,iel))

      tracel(3 + (iloc-1)*idimt) = rtp(iel,iw) &
          + (omegax*xyzcen(2,iel) - omegay*xyzcen(1,iel))

    enddo

    call psteva(nummai, namevr, idimt, ientla, ivarpr,        &
    !==========
                ntcabs, ttcabs, tracel, rbid, rbid)
  endif

  ! Vitesse et pression relatives en cas de calcul en repère fixe

  if (imobil.eq.1) then

    ipcrom = ipproc(irom)

    namevr = 'Rel Pressure'
    idimt = 1
    ientla = 0
    ivarpr = 0

    do iloc = 1, ncelps

      iel = lstcel(iloc)

      pcentr =   0.5d0*((omegay*xyzcen(3,iel) - omegaz*xyzcen(2,iel))**2 &
                      + (omegaz*xyzcen(1,iel) - omegax*xyzcen(3,iel))**2 &
                      + (omegax*xyzcen(2,iel) - omegay*xyzcen(1,iel))**2)

      tracel(iloc) = rtp(iel,ipr) - propce(iel,ipcrom)*pcentr

    enddo

    call psteva(nummai, namevr, idimt, ientla, ivarpr,        &
    !==========
                ntcabs, ttcabs, tracel, rbid, rbid)


    namevr = 'Rel Velocity'
    idimt = 3
    ientla = 1
    ivarpr = 0

    do iloc = 1, ncelps

      iel = lstcel(iloc)

      tracel(1 + (iloc-1)*idimt) = rtp(iel,iu) &
          - (omegay*xyzcen(3,iel) - omegaz*xyzcen(2,iel))

      tracel(2 + (iloc-1)*idimt) = rtp(iel,iv) &
          - (omegaz*xyzcen(1,iel) - omegax*xyzcen(3,iel))

      tracel(3 + (iloc-1)*idimt) = rtp(iel,iw) &
          - (omegax*xyzcen(2,iel) - omegay*xyzcen(1,iel))

    enddo

    call psteva(nummai, namevr, idimt, ientla, ivarpr,        &
    !==========
                ntcabs, ttcabs, tracel, rbid, rbid)

  endif


!===============================================================================
! 1.2. Boundary
!===============================================================================

else if (numtyp .eq. -2) then


  ! 1.2.1 output y+ at the boundary
  ! -------------------------------

  if (mod(ipstdv,ipstyp).eq.0) then

    ! Variable name

    do ii = 1, 32
      namevr (ii:ii) = ' '
    enddo

    namevr = 'Yplus'

    idimt = 1  ! variable dimension (3: vector, 1: scalar)

    ! Compute variable on boundary faces

    do iloc = 1, nfbrps
      ifac = lstfbr(iloc)
      trafbr(1 + (iloc-1)*idimt) = yplbr(ifac)
    enddo

    ! Non interleaved values, defined in work array

    ientla = 0
    ivarpr = 0

    call psteva(nummai, namevr, idimt, ientla, ivarpr,          &
    !==========
                ntcabs, ttcabs, rbid, rbid, trafbr)

  endif ! end of test on output of y+

  !  1.2.2 Projection of variables at boundary with no reconstruction
  !  ----------------------------------------------------------------

  if (mod(ipstdv, ipstcl).eq.0) then

    ! Loop on main cell-based variables
    !----------------------------------

    isaut = 0

    do ivar = 1, nvar

      ! isaut used here to avoid multiple output for vector components 2 and 3

      ipp = ipprtp(ivar)

      isorva = 0
      if (isaut .gt. 0) then
        isaut = isaut - 1
        isorva = 0
      else if (ichrvr(ipp).eq.1) then
        isorva = 1
        name80 = nomvar(ipp)
        namevr = name80(1:32)
      endif

      if (isorva .eq. 1) then

        ! if the sign of the index in ra is negative, we have a vector

        idimt = 1
        if (ipp2ra(ipp).lt.0) then
          idimt = 3
          isaut = 2
        endif

        ! For vectors, remove X, Y, or Z at the end of the name

        if (idimt.eq.3) then
          name80 = nomvar(ipp+1)
          namev1 = name80(1:32)
          name80 = nomvar(ipp+2)
          namev2 = name80(1:32)
          call pstsnv ( namevr , namev1 , namev2 )
          !==========
        endif

        !  Compute non-reconstructed values at boundary faces

        if (     ivelco.eq.0 &
            .or. (ivar.ne.iu .and. ivar.ne.iv .and. ivar.ne.iw .and.        &
                  ivar.ne.iuma .and. ivar.ne.ivma .and. ivar.ne.iwma)) then

          do kk = 0, idimt-1

            iclvar = iclrtp(ivar+kk,icoef)
            do iloc = 1, nfbrps

              ifac = lstfbr(iloc)
              iel = ifabor(ifac)

              trafbr(kk + (iloc-1)*idimt + 1)                       &
                   =   coefa(ifac,iclvar)                           &
                     + coefb(ifac,iclvar)*rtp(iel,ivar+kk)

            enddo

          enddo

        else if (ivar.eq.iu) then

          do kk = 0, idimt-1

            do iloc = 1, nfbrps

              ifac = lstfbr(iloc)
              iel = ifabor(ifac)

              trafbr(kk + (iloc-1)*idimt + 1)                       &
                   =   coefau(kk+1,ifac)                            &
                     + coefbu(kk+1,1,ifac)*rtp(iel,ivar)            &
                     + coefbu(kk+1,2,ifac)*rtp(iel,ivar+1)          &
                     + coefbu(kk+1,3,ifac)*rtp(iel,ivar+2)

            enddo

          enddo

        else if (ivar.eq.iuma) then

          do kk = 0, idimt-1

            do iloc = 1, nfbrps

              ifac = lstfbr(iloc)
              iel = ifabor(ifac)

              trafbr(kk + (iloc-1)*idimt + 1)                       &
                   =   claale(kk+1,ifac)                            &
                     + clbale(kk+1,1,ifac)*rtp(iel,ivar)            &
                     + clbale(kk+1,2,ifac)*rtp(iel,ivar+1)          &
                     + clbale(kk+1,3,ifac)*rtp(iel,ivar+2)

            enddo

          enddo

        endif

        ! Interleaved values, defined on work array
        ientla = 1
        ivarpr = 0

        call psteva(nummai, namevr, idimt, ientla, ivarpr,        &
        !==========
                    ntcabs, ttcabs, rbid, rbid, trafbr)

      endif ! End of variable output

    enddo ! End of loop on variables

  endif ! End of test on variable boundary values output

  ! Output thermal flux at boundary
  ! -------------------------------
  !  If working with enthalpy, compute an enthalpy flux

  if (mod(ipstdv,ipstft).eq.0) then

    if (iscalt.gt.0 .and. nscal.gt.0 .and. iscalt.le.nscal) then

      !       Initialisation
      do ii = 1, 32
        namevr (ii:ii) = ' '
      enddo

      !       Nom de la variable
      namevr = 'Input thermal flux W.m-2'

      !       Dimension de la variable (3 = vecteur, 1=scalaire)
      idimt = 1

      !       Numero de la variable

      iscal  = iscalt
      ivar   = isca(iscal)
      iclvar = iclrtp(ivar,icoef)
      iclvaf = iclrtp(ivar,icoeff)

      !       Calcul des valeurs de la variable sur les faces de bord

      !          Reservation de la memoire pour reconstruction

      allocate(treco(nfabor))

      !          Calcul du gradient de la temperature / enthalpie


      !      Pour calculer le gradient de Temperature
      !        - dans les calculs paralleles, il est necessaire que
      !          les cellules situees sur un bord de sous-domaine connaissent
      !          la valeur de temperature dans les cellules situees en
      !          vis-a-vis sur le sous-domaine voisin.
      !        - dans les calculs periodiques, il est necessaire que
      !          les cellules periodiques aient acces a la valeur de la
      !          temperature des cellules periodiques correspondantes

      !      Pour cela, il est necessaire d'appeler les routines de
      !        de synchronisation des halos pour echanger les valeurs de temperature
      !        avant de calculer le gradient.
      !      En effet, on se situe ici a la fin du pas de temps n. Or,
      !        les variables RTP ne seront echangees qu'en debut du pas de
      !        temps n+1. Ici, seules les variables RTPA (obtenues a la fin
      !        du pas de temps n-1) ont deja ete echangees.

      !      Si le calcul n'est ni periodique, ni parallele, on peut conserver
      !        appels (les tests sur IPERIO et IRANGP assurent la generalite)


      !          Echange pour le parallelisme et la periodicite

      if (irangp.ge.0.or.iperio.eq.1) then
        call synsca(rtp(1,ivar))
        !==========
      endif

      ! Allocate a temporary array for the gradient calculation
      allocate(grad(ncelet,3))

      !          Calcul du gradient

      inc = 1
      iccocg = 1
      nswrgp = nswrgr(ivar)
      imligp = imligr(ivar)
      iwarnp = iwarni(ivar)
      epsrgp = epsrgr(ivar)
      climgp = climgr(ivar)
      extrap = extrag(ivar)

      call grdcel &
      !==========
 ( ivar   , imrgra , inc    , iccocg , nswrgp , imligp ,          &
   iwarnp , nfecra ,                                              &
   epsrgp , climgp , extrap ,                                     &
   rtp(1,ivar) , coefa(1,iclvar) , coefb(1,iclvar) ,              &
   grad   )


      !          Calcul de la valeur reconstruite dans les cellules de bord

      do ifac = 1, nfabor
        iel = ifabor(ifac)
        diipbx = diipb(1,ifac)
        diipby = diipb(2,ifac)
        diipbz = diipb(3,ifac)
        treco(ifac) = rtp(iel,ivar)                  &
             + diipbx*grad(iel,1)                          &
             + diipby*grad(iel,2)                          &
             + diipbz*grad(iel,3)
      enddo

      ! Free memory
      deallocate(grad)
      deallocate(treco)


      ! Calcul du flux convectif et diffusif

      if (ivisls(iscal).gt.0) then
        ipcvsl = ipproc(ivisls(iscal))
      else
        ipcvsl = 0
      endif
      ipcvst = ipproc(ivisct)
      iflmab = ipprob(ifluma(ivar))

      do iloc = 1, nfbrps
        ifac = lstfbr(iloc)
        iel = ifabor(ifac)

        if (ipcvsl.gt.0) then
          xvsl = propce(iel,ipcvsl)
        else
          xvsl = visls0(iscal)
        endif
        srfbn = max(surfbn(ifac), epzero**2)
        visct  = propce(iel,ipcvst)
        flumab = propfb(ifac,iflmab)

        trafbr(1 + (iloc-1)*idimt) =                                   &
             (coefa(ifac,iclvaf) + coefb(ifac,iclvaf)*rtp(iel,ivar))   &
             - flumab/srfbn*                                           &
             (coefa(ifac,iclvar) + coefb(ifac,iclvar)*rtp(iel,ivar))

      enddo

      !             Valeurs entrelacées, définies sur tableau de travail
      ientla = 1
      ivarpr = 0

      call psteva(nummai, namevr, idimt, ientla, ivarpr,        &
      !==========
                  ntcabs, ttcabs, rbid, rbid, trafbr)

    endif
    !         Fin du test sur variable thermique

  endif
  !      Fin du test sur sortie des flux thermiques

! --    1.2.4 TRAITEMENT DES EFFORTS AUX BORDS
!       --------------------------------------

  if (mod(ipstdv,ipstfo).eq.0) then

!       Initialisation
    do ii = 1, 32
      NAMEVR (II:II) = ' '
    enddo

!       Nom de la variable
    NAMEVR = 'Efforts'

!       Dimension de la variable (3 = vecteur, 1=scalaire)
    idimt = 3

!       Calcul des valeurs de la variable sur les faces de bord

    do iloc = 1, nfbrps
      ifac = lstfbr(iloc)
      srfbn = surfbn(ifac)
      trafbr(1 + (iloc-1)*idimt ) = forbr(1,ifac)/srfbn
      trafbr(2 + (iloc-1)*idimt ) = forbr(2,ifac)/srfbn
      trafbr(3 + (iloc-1)*idimt ) = forbr(3,ifac)/srfbn
    enddo

!           Valeurs entrelacées, définies sur tableau de travail
      ientla = 1
      ivarpr = 0

      call psteva(nummai, namevr, idimt, ientla, ivarpr,          &
      !==========
                  ntcabs, ttcabs, rbid, rbid, trafbr)

  endif
! fin du test sur sortie des efforts

endif
!     Fin du test sur le numero de maillage post.

!===============================================================================
!     2.1. VARIABLES LAGRANGIENNES
!===============================================================================

if (nummai .eq. -1) then

  if (iilagr.gt.0 .and. istala.ge.1) then

!         Toutes les statistiques standard sont de dimension 1,
!         et sont definies ou calculees sur tableau de travail
!         de maniere non entrelacee (sans importance de toutes
!         manieres pour une variable scalaire)

    idimt  = 1
    ientla = 0

    iii = nvlsta-nvlsts

    do icla  = 0, nbclst

!     -> si ICLA = 0 : statistiques globales
!        si 0 < ICLA =< NBCLST : statistiques par groupe

      do ivarl = 1, nvlsta

        ivarl1 = icla*nvlsta +ivarl
        ivarlm = ivarl1
        ilpd1  = icla*nvlsta +ilpd
        iflu   = 0

        if (ivarl.le.iii) then
          if (ivarl.eq.ivarl1) then
            name80 = nomlag(ivarl)
          else
            write(name80,'(a8,a4,i3)') nomlag(ivarl),'_grp',icla
          endif
        else if (nvlsts.gt.0) then
          if (ivarl.eq.ivarl1) then
            name80 = nomlag(ilvu(ivarl-iii))
          else
            write(name80,'(a8,a4,i3)')                            &
                 nomlag(ilvu(ivarl-iii)),'_grp',icla
          endif
        endif

        namevr = name80(1:32)

        call uslaen                                               &
        !==========
 ( nvar   , nscal  , nvlsta ,                                     &
   ivarl  , ivarl1 , ivarlm , iflu   , ilpd1  , icla   ,          &
   dt     , rtpa   , rtp    , propce , propfa , propfb ,          &
   coefa  , coefb  , statce , stativ , tracel )

!           La variable est deja definie sur le maillage volumique
!           global ; on utilise donc l'indirection  (donc IVARPR = 1)
        ivarpr = 1

        call psteva(nummai, namevr, idimt, ientla, ivarpr,        &
                    ntcabs, ttcabs, tracel, rbid, rbid)
      enddo

      do ivarl = 1, nvlsta-1

        ivarl1 = icla*(nvlsta-1)+ivarl
        ivarlm = icla*nvlsta+ivarl
        ilpd1  = icla*nvlsta +ilpd
        iflu   = 1

        if (ivarl.le.iii) then
          if (ivarl.eq.ivarl1) then
            name80 = nomlav(ivarl)
          else
            write(name80,'(a8,a4,i3)') nomlav(ivarl),'_grp',icla
          endif
        else if (nvlsts.gt.0) then
          if (ivarl.eq.ivarl1) then
            name80 = nomlav(ilvu(ivarl-iii))
          else
            write(name80,'(a8,a4,i3)')                            &
                 nomlav(ilvu(ivarl-iii)),'_grp',icla
          endif
        endif

        namevr = name80(1:32)

        call uslaen                                               &
        !==========
 ( nvar   , nscal  , nvlsta ,                                     &
   ivarl  , ivarl1 , ivarlm , iflu   , ilpd1  , icla   ,          &
   dt     , rtpa   , rtp    , propce , propfa , propfb ,          &
   coefa  , coefb  , statce , stativ , tracel )

!           La variable est deja definie sur le maillage volumique
!           global ; on utilise donc l'indirection  (donc IVARPR = 1)
        ivarpr = 1

        call psteva(nummai, namevr, idimt, ientla, ivarpr,        &
                    ntcabs, ttcabs, tracel, rbid, rbid)
      enddo

    enddo

  endif

endif

if (nummai.eq.-2) then

  if (iilagr.gt.0 .and. iensi3.eq.1) then

    iii = nvisbr-nusbor

    do ivarl = 1,nvisbr

      if (ivarl.le.iii) then
        name80 = nombrd(ivarl)
      else if (nusbor.gt.0) then
        name80 = nombrd(iusb(ivarl-iii))
      endif
      namevr = name80(1:32)

      if (imoybr(ivarl).eq.2) then

        do iloc = 1, nfbrps
          ifac = lstfbr(iloc)
          if (statfb(ifac,inbr).gt.seuilf) then
            trafbr(iloc) = statfb(ifac,ivarl)/statfb(ifac,inbr)
          else
            trafbr(iloc) = 0.d0
          endif
        enddo

      else if (imoybr(ivarl).eq.1) then

        do iloc = 1, nfbrps
          ifac = lstfbr(iloc)
          if (statfb(ifac,inbr).gt.seuilf) then
            trafbr(iloc) = statfb(ifac,ivarl) / tstatp
          else
            trafbr(iloc) = 0.d0
          endif
        enddo

      else

        do iloc = 1, nfbrps
          ifac = lstfbr(iloc)
          if (statfb(ifac,inbr).gt.seuilf) then
            trafbr(iloc) = statfb(ifac,ivarl)
          else
            trafbr(iloc) = 0.d0
          endif
        enddo

      endif

      idimt  = 1
      ientla = 0
      ivarpr = 0

      call psteva(nummai, namevr, idimt, ientla, ivarpr,          &
                  ntcabs, ttcabs, rbid, rbid, trafbr)

    enddo

    !NAME80 = 'lagrangian_boundary_zones'
    !namevr = name80(1:32)

    !do iloc = 1, nfbrps
    !  ifac = lstfbr(iloc)
    !  trafbr(iloc) = ia(iifrla+ifac-1) !! TODO: ifrlag (cf caltri)
    !enddo

    !idimt  = 1
    !ientla = 0
    !ivarpr = 0

    !call psteva(nummai, namevr, idimt, ientla, ivarpr,            &
    !     ntcabs, ttcabs, rbid, rbid, trafbr)

  endif
endif
!     Fin du test sur le numero de maillage post.

!===============================================================================
!     2.2. VARIABLES RADIATIVES AUX FRONTIERES
!===============================================================================


if (nummai.eq.-2) then

  if (iirayo.gt.0) then

    do ivarl = 1,nbrayf

      if (irayvf(ivarl).eq.1) then

        name80 = nbrvaf(ivarl)
        namevr = name80(1:32)

        if (ivarl .eq. itparp)      then
          ipp =  ipprob(itparo)
        else if (ivarl .eq. iqincp) then
          ipp = ipprob(iqinci)
        else if (ivarl .eq. ixlamp)  then
          ipp = ipprob(ixlam)
        else if (ivarl .eq. iepap)   then
          ipp = ipprob(iepa)
        else if (ivarl .eq. iepsp)   then
          ipp = ipprob(ieps)
        else if (ivarl .eq. ifnetp)  then
          ipp = ipprob(ifnet)
        else if (ivarl .eq. ifconp) then
          ipp = ipprob(ifconv)
        else if (ivarl .eq. ihconp) then
          ipp = ipprob(ihconv)
        endif

        do iloc = 1, nfbrps
          ifac = lstfbr(iloc)
          trafbr(iloc) = propfb(ifac,ipp)
        enddo

        idimt  = 1
        ientla = 0
        ivarpr = 0

        call psteva(nummai, namevr, idimt, ientla, ivarpr,        &
                    ntcabs, ttcabs, rbid, rbid, trafbr)

      endif
    enddo

    name80 = 'radiative_boundary_zones'
    namevr = name80(1:32)

    do iloc = 1, nfbrps
      ifac = lstfbr(iloc)
      trafbr(iloc) = izfrad(ifac)
    enddo

    idimt  = 1
    ientla = 0
    ivarpr = 0
!
    call psteva(nummai, namevr, idimt, ientla, ivarpr,            &
         ntcabs, ttcabs, rbid, rbid, trafbr)

  endif
endif

!===============================================================================
!     2.3. VARIABLES ELECTRIQUES
!===============================================================================

if (     ippmod(ieljou).ge.1                                      &
    .or. ippmod(ielarc).ge.1                                      &
    .or. ippmod(ielion).ge.1) then

  call uselen                                                     &
  !==========
 ( nummai ,                                                       &
   nvar   , nscal  ,                                              &
   ncelps , nfacps , nfbrps ,                                     &
   lstcel , lstfac , lstfbr ,                                     &
   dt     , rtpa   , rtp    , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   tracel , trafac , trafbr )

endif

!--------
! Formats
!--------

!----
! End
!----

return
end subroutine
