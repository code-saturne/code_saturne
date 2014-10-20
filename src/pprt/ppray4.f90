!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2014 EDF S.A.
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

subroutine ppray4 &
!================

 ( mode   ,                                                       &
   itypefb ,                                                      &
   propce ,                                                       &
   tparop , hparop , tempk  )

!===============================================================================
! FONCTION :
! ----------


!   SOUS-PROGRAMME PHYSIQUES PARTICULIERES
!   POUR LES CONVERSIONS TEMPERATURE <-> ENTHALPIE
!      POUR LE MELANGE GAZEUX



!     A) L'ENTHALPIE DOIT ETRE CONVERTIE EN TEMPERATURE EN KELVIN
!        IL FAUT ECRIRE UNE LOI DE CONVERSION :
!        ENTHALPIE   -> TEMPERATURE (MODE =  1)
!        POUR REMPLIR TEMPK


!     B) SI DE PLUS ON CALCULE LES TEMPERATURES DE PAROI
!        ALORS IL FAUT FOURNIR UNE LOI DE CONVERSION :
!        TEMPERATURE -> ENTHALPIE   (MODE = -1)
!        POUR REMPLIR HPAROP


!-------------------------------------------------------------------------------
!ARGU                             ARGUMENTS
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! mode             ! e  ! <-- ! type de conversion enthal<->tempk              !
! itypefb(nfabor   ! te ! <-- ! type des faces de bord                         !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! cofrua,cofrub    ! tr ! --> ! conditions aux limites aux                     !
!(nfabor)          !    !     !    faces de bord pour la luminances            !
! tempk(ncelet)    ! tr ! --> ! temperature en kelvin                          !
! hparop(nfabor    ! tr ! --> ! enthalpie massique de paroi en j/kg            !
!                  !    !     ! (en degres celsius ou kelvin)                  !
! tparop(nfabor    ! tr ! <-- ! temperature de paroi en kelvin                 !
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
use pointe
use ppppar
use ppthch
use coincl
use cpincl
use cs_fuel_incl
use ppincl
use radiat
use mesh
use field

!===============================================================================

implicit none

! Arguments

integer          mode

integer          itypefb(nfabor)

double precision propce(ncelet,*)

double precision tempk(ncelet)
double precision tparop(nfabor), hparop(nfabor)

! Local variables

integer          iel , ifac , icla , icha , isol , ige
integer          ipcx2c , ixchcl, ixckcl, ixnpcl, igg, iii
integer          iesp
double precision coefg(ngazgm), coefe(ngazem)
double precision f1mc(ncharm) , f2mc(ncharm)
double precision x2t , h2 , x2h2 , hf , xsolid(nsolim),t1
double precision ym    (ngazgm)
double precision diamgt,masgut,mkgout,mfgout,mkfini,rhofol
double precision, dimension(:), pointer :: bym1, bym2, bym3

type(pmapper_double_r1), dimension(:), allocatable :: cvar_xch, cvar_xck
type(pmapper_double_r1), dimension(:), allocatable :: cvar_xnp, cvar_xwt
type(pmapper_double_r1), dimension(:), allocatable :: cvar_f1m, cvar_f2m
type(pmapper_double_r1), dimension(:), allocatable :: cvara_yfol
type(pmapper_double_r1), dimension(:), allocatable :: cvar_ycoel

!===============================================================================

!===============================================================================
! 1 - INITIALISATIONS GENERALES
!===============================================================================

!===============================================================================
!  2.1 - CALCUL DE LA TEMPERATURE EN KELVIN AUX CELLULES
!===============================================================================

!---> CONVERSION ENTHALPIE -> TEMPERATURE (MODE =  1)
!     AUX CELLULES FLUIDES
!     (UTILE SI NON TRANSPARENT)
!     -----------------------------------------------


if (mode.eq.1) then

! ---- Combustion gaz : Flamme de Premelange ou Flamme de Diffusion

  if ( ippmod(icoebu).ge.0 .or.                                   &
       ippmod(icod3p).ge.0      ) then

    do iel = 1,ncel
      tempk(iel) = propce(iel, ipproc(itemp))
    enddo

! ---- Combustion charbon pulverise

  else if ( ippmod(iccoal).ge.0 ) then

    do iel = 1,ncel
      tempk(iel) = propce(iel, ipproc(itemp1))
    enddo

! ---- Combustion fuel

  else if ( ippmod(icfuel).ge.0 ) then

    do iel = 1,ncel
      tempk(iel) = propce(iel, ipproc(itemp1))
    enddo

! ---- Module Electrique

  else if ( ippmod(ieljou).ge.1 .or.                              &
            ippmod(ielarc).ge.1 .or.                              &
            ippmod(ielion).ge.1       ) then

    do iel = 1,ncel
      tempk(iel) = propce(iel, ipproc(itemp))
    enddo

  endif

endif


!===============================================================================
!  2.2 - CALCUL DE L'ENTHALPIE AU FACES DE BORD
!===============================================================================

!---> CONVERSION TEMPERATURE -> ENTHALPIE (MODE = -1)
!     AUX FACES DE BORD DE PAROI
!     (INUTILE POUR LES PAROIS ISOTH(IFAC)=3 ET EPS(IFAC)=0)
!     ------------------------------------------------------


if (mode.eq.-1) then

  call field_get_val_s(ibym(1), bym1)
  call field_get_val_s(ibym(2), bym2)
  call field_get_val_s(ibym(3), bym3)

  ! Arrays of pointers containing the fields values for each class
  ! (loop on cells outside loop on classes)
  if ( ippmod(iccoal).ge.0 ) then

    allocate(cvar_xch(nclacp))
    allocate(cvar_xck(nclacp))
    allocate(cvar_xnp(nclacp))
    allocate(cvar_xwt(nclacp))

    do icla = 1, nclacp
      call field_get_val_s(ivarfl(isca(ixch(icla))), cvar_xch(icla)%p)
      call field_get_val_s(ivarfl(isca(ixck(icla))), cvar_xck(icla)%p)
      call field_get_val_s(ivarfl(isca(inp(icla))), cvar_xnp(icla)%p)
      if ( ippmod(iccoal).eq.1 ) then
        call field_get_val_s(ivarfl(isca(ixwt(icla))), cvar_xwt(icla)%p)
      endif
    enddo

    allocate(cvar_f1m(ncharb))
    allocate(cvar_f2m(ncharb))

    do icha = 1, ncharb
      call field_get_val_s(ivarfl(isca(if1m(icha))), cvar_f1m(icha)%p)
      call field_get_val_s(ivarfl(isca(if2m(icha))), cvar_f2m(icha)%p)
    enddo

  else if ( ippmod(icfuel).ge.0 ) then

    allocate(cvara_yfol(nclafu))

    do icla = 1,nclafu
      call field_get_val_prev_s(ivarfl(isca(iyfol(icla))), cvara_yfol(icla)%p)
    enddo

  else if ( ippmod(ielarc).ge.1  ) then

    allocate(cvar_ycoel(ngazg-1))

    do iesp = 1, ngazg-1
      call field_get_val_s(ivarfl(isca(iycoel(iesp))), cvar_ycoel(iesp)%p)
    enddo

  endif

  do ifac = 1,nfabor

    if (itypefb(ifac).eq.iparoi .or. itypefb(ifac).eq.iparug) then

      iel = ifabor(ifac)

      ! ---- Combustion gaz : Flamme de Premelange ou Flamme de Diffusion

      if (ippmod(icoebu).ge.0 .or. ippmod(icod3p).ge.0) then

        do igg = 1, ngazgm
          coefg(igg) = zero
        enddo
        coefg(1) = bym1(ifac)
        coefg(2) = bym2(ifac)
        coefg(3) = bym3(ifac)
        mode     = -1
        call cothht                                               &
        !==========
        ( mode   , ngazg , ngazgm  , coefg  ,                     &
          npo    , npot   , th     , ehgazg ,                     &
          hparop(ifac) , tparop(ifac) )

      ! ---- Combustion charbon pulverise : nouveau modele

      else if ( ippmod(iccoal).ge.0 ) then

        x2t  = zero
        x2h2 = zero
        do icla = 1, nclacp
          icha   = ichcor(icla)
          ixchcl = isca(ixch(icla))
          ixckcl = isca(ixck(icla))
          ixnpcl = isca(inp(icla ))
          ipcx2c = ipproc(ix2(icla))
          x2t = x2t + propce(iel,ipcx2c)
          h2 = zero
          do isol = 1, nsolim
            xsolid(isol) = zero
          enddo
          if (propce(iel,ipcx2c).gt.epsicp) then
            xsolid(ich(icha) ) = cvar_xch(icla)%p(iel)             &
                               / propce(iel,ipcx2c)
            xsolid(ick(icha) ) = cvar_xck(icla)%p(iel)             &
                               / propce(iel,ipcx2c)
            xsolid(iash(icha)) = cvar_xnp(icla)%p(iel)*xmash(icla) &
                               / propce(iel,ipcx2c)
            if ( ippmod(iccoal).eq.1 ) then
              xsolid(iwat(icha)) = cvar_xwt(icla)%p(iel)           &
                                  /propce(iel,ipcx2c)
            endif
            iii = icla
            t1 = tparop(ifac)

            call cs_coal_htconvers2                            &
           !============================
           ( mode , iii , h2 , xsolid ,  tparop(ifac) , t1 )

          endif
          x2h2 = x2h2 + propce(iel,ipcx2c)*h2
        enddo

        do ige = 1, ngazem
          coefe(ige) = zero
        enddo
        do icha = 1, ncharb
          f1mc(icha) = cvar_f1m(icha)%p(iel)                       &
                     / (1.d0-x2t)
          f2mc(icha) =  cvar_f2m(icha)%p(iel)                      &
                     / (1.d0-x2t)
        enddo
        do icha = (ncharb+1), ncharm
          f1mc(icha) = zero
          f2mc(icha) = zero
        enddo
        do ige = 1, ngazg
          coefe(ige) = propce(iel,ipproc(iym1(ige)))
        enddo

        call cs_coal_htconvers1(mode,hf,coefe,f1mc,f2mc,tparop(ifac))
       !============================

        hparop(ifac) = (1.d0-x2t)*hf+x2h2

! ---- Combustion fuel

      else if ( ippmod(icfuel).ge.0 ) then

        x2t  = zero
        x2h2 = zero

        do icla=1,nclafu

          x2t = x2t+cvara_yfol(icla)%p(iel)

          mkfini = rho0fl*pi/6.d0*dinikf(icla)**3
          rhofol = propce(iel,ipproc(irom2(icla)))
          diamgt = propce(iel,ipproc(idiam2(icla)))
          masgut = rhofol*pi/6.d0*diamgt**3
          if (diamgt.le.dinikf(icla)) then
            mkgout = masgut
          else
            mkgout = mkfini
          endif
          mfgout = masgut - mkgout
          xsolid(1) = 1.d0-fkc
          xsolid(2) = fkc
          if(masgut.gt.epzero) then
            xsolid(1) = mfgout / masgut
            xsolid(2) = mkgout / masgut
          endif
          xsolid(1) = min(1.d0,max(0.d0,xsolid(1)))
          xsolid(2) = min(1.d0,max(0.d0,xsolid(2)))

          call cs_fuel_htconvers2                                  &
         !=======================
        ( mode   , h2     , xsolid , tparop(ifac) )

          x2h2 =  x2h2 + cvara_yfol(icla)%p(iel)*h2

        enddo

        do ige = 1,ngaze
          coefe(ige) = propce(iel,ipproc(iym1(ige)))
        enddo
        call cs_fuel_htconvers1                                   &
       !=======================
        ( mode   , hf     , coefe  , tparop(ifac) )
          hparop(ifac) = (1.d0-x2t)*hf+x2h2

! ---- Module Electrique


      else if ( ippmod(ieljou).ge.1  ) then

        call usthht (mode,hparop(ifac),tparop(ifac))

      else if ( ippmod(ielarc).ge.1  ) then

        if ( ngazg .eq. 1 ) then
          ym(1) = 1.d0
          call elthht(mode,ngazg,ym,hparop(ifac),tparop(ifac))
        else
          ym(ngazg) = 1.d0
          do iesp = 1, ngazg-1
            ym(iesp) = cvar_ycoel(iesp)%p(iel)
            ym(ngazg) = ym(ngazg) - ym(iesp)
          enddo
          call elthht(mode,ngazg,ym,hparop(ifac),tparop(ifac))
        endif

!     ELSE IF ( IPPMOD(IELION).GE.1 ) THEN
!     ... to be implemented ...

      endif

    endif

  enddo

  if ( ippmod(iccoal).ge.0 ) then
    deallocate(cvar_xch, cvar_xck, cvar_xnp)
    if ( ippmod(iccoal).eq.1 ) then
      deallocate(cvar_xwt)
    endif
    deallocate(cvar_f1m, cvar_f2m)
  else if ( ippmod(icfuel).ge.0 ) then
    deallocate(cvara_yfol)
  else if ( ippmod(ielarc).ge.1  ) then
    deallocate(cvar_ycoel)
  endif

endif

!----
! FIN
!----

return
end subroutine
