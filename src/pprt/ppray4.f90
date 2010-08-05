!-------------------------------------------------------------------------------

!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2009 EDF S.A., France

!     contact: saturne-support@edf.fr

!     The Code_Saturne Kernel is free software; you can redistribute it
!     and/or modify it under the terms of the GNU General Public License
!     as published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.

!     The Code_Saturne Kernel is distributed in the hope that it will be
!     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
!     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.

!     You should have received a copy of the GNU General Public License
!     along with the Code_Saturne Kernel; if not, write to the
!     Free Software Foundation, Inc.,
!     51 Franklin St, Fifth Floor,
!     Boston, MA  02110-1301  USA

!-------------------------------------------------------------------------------

subroutine ppray4 &
!================

 ( idbia0 , idbra0 ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , iphas  ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &

   mode   ,                                                       &

   ifacel , ifabor , ifmfbr , ifmcel , iprfml , itypfb ,          &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   w1     , w2     , w3     , w4     , w5     , w6     ,          &
   tparop , hparop , tempk  ,                                     &
   rdevel , rtuser ,                                              &
   ra     )

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
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
! ndim             ! i  ! <-- ! spatial dimension                              !
! ncelet           ! i  ! <-- ! number of extended (real + ghost) cells        !
! ncel             ! i  ! <-- ! number of cells                                !
! nfac             ! i  ! <-- ! number of interior faces                       !
! nfabor           ! i  ! <-- ! number of boundary faces                       !
! nfml             ! i  ! <-- ! number of families (group classes)             !
! nprfml           ! i  ! <-- ! number of properties per family (group class)  !
! nnod             ! i  ! <-- ! number of vertices                             !
! lndfac           ! i  ! <-- ! size of nodfac indexed array                   !
! lndfbr           ! i  ! <-- ! size of nodfbr indexed array                   !
! ncelbr           ! i  ! <-- ! number of cells with faces on boundary         !
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! iphas            ! i  ! <-- ! phase number                                   !
! nideve, nrdeve   ! i  ! <-- ! sizes of idevel and rdevel arrays              !
! nituse, nrtuse   ! i  ! <-- ! sizes of ituser and rtuser arrays              !
! mode             ! e  ! <-- ! type de conversion enthal<->tempk              !
! ifacel(2, nfac)  ! ia ! <-- ! interior faces -> cells connectivity           !
! ifabor(nfabor)   ! ia ! <-- ! boundary faces -> cells connectivity           !
! ifmfbr(nfabor)   ! ia ! <-- ! boundary face family numbers                   !
! ifmcel(ncelet)   ! ia ! <-- ! cell family numbers                            !
! iprfml           ! ia ! <-- ! property numbers per family                    !
!  (nfml, nprfml)  !    !     !                                                !
! itypfb(nfabor    ! te ! <-- ! type des faces de bord                         !
! ipnfac(nfac+1)   ! ia ! <-- ! interior faces -> vertices index (optional)    !
! nodfac(lndfac)   ! ia ! <-- ! interior faces -> vertices list (optional)     !
! ipnfbr(nfabor+1) ! ia ! <-- ! boundary faces -> vertices index (optional)    !
! nodfbr(lndfbr)   ! ia ! <-- ! boundary faces -> vertices list (optional)     !
! idevel(nideve)   ! ia ! <-> ! integer work array for temporary development   !
! ituser(nituse)   ! ia ! <-> ! user-reserved integer work array               !
! ia(*)            ! ia ! --- ! main integer work array                        !
! xyzcen           ! ra ! <-- ! cell centers                                   !
!  (ndim, ncelet)  !    !     !                                                !
! surfac           ! ra ! <-- ! interior faces surface vectors                 !
!  (ndim, nfac)    !    !     !                                                !
! surfbo           ! ra ! <-- ! boundary faces surface vectors                 !
!  (ndim, nfabor)  !    !     !                                                !
! cdgfac           ! ra ! <-- ! interior faces centers of gravity              !
!  (ndim, nfac)    !    !     !                                                !
! cdgfbo           ! ra ! <-- ! boundary faces centers of gravity              !
!  (ndim, nfabor)  !    !     !                                                !
! xyznod           ! ra ! <-- ! vertex coordinates (optional)                  !
!  (ndim, nnod)    !    !     !                                                !
! volume(ncelet)   ! ra ! <-- ! cell volumes                                   !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
! cofrua,cofrub    ! tr ! --> ! conditions aux limites aux                     !
!(nfabor)          !    !     !    faces de bord pour la luminances            !
! w1...6(ncelet    ! tr ! --- ! tableau de travail                             !
! tempk(ncelet)    ! tr ! --> ! temperature en kelvin                          !
! hparop(nfabor    ! tr ! --> ! enthalpie massique de paroi en j/kg            !
!                  !    !     ! (en degres celsius ou kelvin)                  !
! tparop(nfabor    ! tr ! <-- ! temperature de paroi en kelvin                 !
! rdevel(nrdeve)   ! ra ! <-> ! real work array for temporary development      !
! rtuser(nrtuse)   ! ra ! <-> ! user-reserved real work array                  !
! ra(*)            ! ra ! --- ! main real work array                           !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!-------------------------------------------------------------------------------
!===============================================================================

implicit none

!===============================================================================
! Common blocks
!===============================================================================

include "paramx.f90"
include "numvar.f90"
include "entsor.f90"
include "optcal.f90"
include "cstphy.f90"
include "cstnum.f90"
include "pointe.f90"
include "ppppar.f90"
include "ppthch.f90"
include "coincl.f90"
include "cpincl.f90"
include "fuincl.f90"
include "ppincl.f90"
include "radiat.f90"


!===============================================================================

! Arguments

integer          idbia0 , idbra0
integer          ndim   , ncelet , ncel   , nfac   , nfabor
integer          nfml   , nprfml
integer          nnod   , lndfac , lndfbr , ncelbr
integer          nvar   , nscal  , iphas
integer          nideve , nrdeve , nituse , nrtuse

integer          mode

integer          ifacel(2,nfac) , ifabor(nfabor)
integer          ifmfbr(nfabor) , ifmcel(ncelet)
integer          iprfml(nfml,nprfml) , itypfb(nfabor)
integer          ipnfac(nfac+1), nodfac(lndfac)
integer          ipnfbr(nfabor+1), nodfbr(lndfbr)
integer          idevel(nideve), ituser(nituse), ia(*)
double precision xyzcen(ndim,ncelet)
double precision surfac(ndim,nfac), surfbo(ndim,nfabor)
double precision cdgfac(ndim,nfac), cdgfbo(ndim,nfabor)
double precision xyznod(ndim,nnod), volume(ncelet)
double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision w1(ncelet), w2(ncelet), w3(ncelet)
double precision w4(ncelet), w5(ncelet), w6(ncelet)

double precision tempk(ncelet)
double precision tparop(nfabor), hparop(nfabor)

double precision rdevel(nrdeve), rtuser(nrtuse), ra(*)


! Local variables

integer          idebia , idebra
integer          iel , ifac , icla , icha , isol , ige
integer          ipcx2c , ixchcl, ixckcl, ixnpcl, igg, iii
integer          iesp
double precision coefg(ngazgm), coefe(ngazem)
double precision f1mc(ncharm) , f2mc(ncharm)
double precision x2t , h2 , x2h2 , hf , xsolid(nsolim),t1
double precision ym    (ngazgm)
double precision diamgt,masgut,mkgout,mfgout,mkfini,rhofol

!===============================================================================

!===============================================================================
! 0 - GESTION MEMOIRE
!===============================================================================

idebia = idbia0
idebra = idbra0

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

  else if ( ippmod(icp3pl).ge.0 ) then

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

  do ifac = 1,nfabor

    if (itypfb(ifac).eq.iparoi .or.                               &
        itypfb(ifac).eq.iparug) then

!   Numero de la cellule en regard

      iel = ifabor(ifac)

! ---- Combustion gaz : Flamme de Premelange ou Flamme de Diffusion

      if ( ippmod(icoebu).ge.0 .or.                               &
           ippmod(icod3p).ge.0      ) then

        do igg = 1, ngazgm
          coefg(igg) = zero
        enddo
        coefg(1) = propfb(ifac,ipprob(iym(1)))
        coefg(2) = propfb(ifac,ipprob(iym(2)))
        coefg(3) = propfb(ifac,ipprob(iym(3)))
        mode     = -1
        call cothht                                               &
        !==========
        ( mode   , ngazg , ngazgm  , coefg  ,                     &
          npo    , npot   , th     , ehgazg ,                     &
          hparop(ifac) , tparop(ifac) )

! ---- Combustion charbon pulverise

      else if ( ippmod(icp3pl).ge.0 ) then

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
            xsolid(ich(icha) ) = rtp(iel,ixchcl)                  &
                               / propce(iel,ipcx2c)
            xsolid(ick(icha) ) = rtp(iel,ixckcl)                  &
                               / propce(iel,ipcx2c)
            xsolid(iash(icha)) = rtp(iel,ixnpcl)*xmash(icla)      &
                               / propce(iel,ipcx2c)
            if ( ippmod(icp3pl).eq.1 ) then
              xsolid(iwat(icha)) = rtp(iel,isca(ixwt(icla)))      &
                                  /propce(iel,ipcx2c)
            endif
            iii = icla
            t1 = tparop(ifac)
            call cpthp2                                           &
            !==========
            ( mode   , iii   , h2     , xsolid ,                  &
              tparop(ifac) , t1 )
          endif
          x2h2 = x2h2 + propce(iel,ipcx2c)*h2
        enddo

        do ige = 1, ngazem
          coefe(ige) = zero
        enddo
        do icha = 1, ncharb
          f1mc(icha) = rtp(iel,isca(if1m(icha)))                  &
                     / (1.d0-x2t)
          f2mc(icha) =  rtp(iel,isca(if2m(icha)))                 &
                     / (1.d0-x2t)
        enddo
        do icha = (ncharb+1), ncharm
          f1mc(icha) = zero
          f2mc(icha) = zero
        enddo
        do ige = 1, (ngaze-2*ncharb)
          coefe(ige) = propce(iel,ipproc(iym1(ige)))
        enddo
        call cpthp1                                               &
        !==========
        ( mode   , hf     , coefe       ,                         &
          f1mc   , f2mc   ,tparop(ifac) )
          hparop(ifac) = (1.d0-x2t)*hf+x2h2

! ---- Combustion fuel

      else if ( ippmod(icfuel).ge.0 ) then

        x2t  = zero
        x2h2 = zero

        do icla=1,nclafu

          x2t = x2t+rtpa(iel,isca(iyfol(icla)))

          mkfini = rho0fl*pi/6.d0*dinikf(icla)**3
          rhofol = propce(iel,ipproc(irom3(icla)))
          diamgt = propce(iel,ipproc(idiam3(icla)))
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

          call futhp2                                             &
          !==========
        ( mode   , h2     , xsolid , tparop(ifac) )

          x2h2 =  x2h2 + rtpa(iel,isca(iyfol(icla)))*h2

        enddo

        do ige = 1,ngaze
          coefe(ige) = propce(iel,ipproc(iym1(ige)))
        enddo
        call futhp1                                               &
        !==========
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
            ym(iesp) = rtp(iel,isca(iycoel(iesp)))
            ym(ngazg) = ym(ngazg) - ym(iesp)
          enddo
          call elthht(mode,ngazg,ym,hparop(ifac),tparop(ifac))
        endif

!     ELSE IF ( IPPMOD(IELION).GE.1 ) THEN
!     ... to be implemented ...

      endif

    endif

  enddo

endif

!----
! FIN
!----

return
end subroutine
