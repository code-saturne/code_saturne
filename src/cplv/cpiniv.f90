!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2011 EDF S.A.
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

subroutine cpiniv &
!================

 ( nvar   , nscal  ,                                              &
   dt     , rtp    , propce , propfa , propfb , coefa  , coefb  )

!===============================================================================
! FONCTION :
! --------

! INITIALISATION DES VARIABLES DE CALCUL
!    POUR LA PHYSIQUE PARTICULIERE : COMBUSTION CP
!    PENDANT DE USINIV.F

! Cette routine est appelee en debut de calcul (suite ou non)
!     avant le debut de la boucle en temps

! Elle permet d'INITIALISER ou de MODIFIER (pour les calculs suite)
!     les variables de calcul,
!     les valeurs du pas de temps


! On dispose ici de ROM et VISCL initialises par RO0 et VISCL0
!     ou relues d'un fichier suite
! On ne dispose des variables VISCLS, CP (quand elles sont
!     definies) que si elles ont pu etre relues dans un fichier
!     suite de calcul

! Les proprietes physiaues sont accessibles dans le tableau
!     PROPCE (prop au centre), PROPFA (aux faces internes),
!     PROPFB (prop aux faces de bord)
!     Ainsi,
!      PROPCE(IEL,IPPROC(IROM  )) designe ROM   (IEL)
!      PROPCE(IEL,IPPROC(IVISCL)) designe VISCL (IEL)
!      PROPCE(IEL,IPPROC(ICP   )) designe CP    (IEL)
!      PROPCE(IEL,IPPROC(IVISLS(ISCAL))) designe VISLS (IEL ,ISCAL)

!      PROPFA(IFAC,IPPROF(IFLUMA(IVAR ))) designe FLUMAS(IFAC,IVAR)

!      PROPFB(IFAC,IPPROB(IROM  )) designe ROMB  (IFAC)
!      PROPFB(IFAC,IPPROB(IFLUMA(IVAR ))) designe FLUMAB(IFAC,IVAR)

! LA MODIFICATION DES PROPRIETES PHYSIQUES (ROM, VISCL, VISCLS, CP)
!     SE FERA EN STANDARD DANS LE SOUS PROGRAMME PPPHYV
!     ET PAS ICI

! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! dt(ncelet)       ! tr ! <-- ! valeur du pas de temps                         !
! rtp              ! tr ! <-- ! variables de calcul au centre des              !
! (ncelet,*)       !    !     !    cellules                                    !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! coefa coefb      ! tr ! <-- ! conditions aux limites aux                     !
!  (nfabor,*)      !    !     !    faces de bord                               !
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
use optcal
use cstphy
use cstnum
use entsor
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use ppcpfu
use mesh

!===============================================================================

implicit none

integer          nvar   , nscal


double precision dt(ncelet), rtp(ncelet,*), propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)

! Local variables

integer          iel, ige, mode, icla, icha

double precision t1init, h1init, coefe(ngazem)
double precision t2init, h2init
double precision f1mc(ncharm), f2mc(ncharm)
double precision xkent, xeent, d2s3

integer imode , isol , ioxy
double precision wmh2o,wmco2,wmn2,wmo2,dmas


! NOMBRE DE PASSAGES DANS LA ROUTINE

integer          ipass
data             ipass /0/
save             ipass

!===============================================================================
!===============================================================================
! 1.  INITIALISATION VARIABLES LOCALES
!===============================================================================

ipass = ipass + 1


d2s3 = 2.d0/3.d0

!===============================================================================
! 2. INITIALISATION DES INCONNUES :
!      UNIQUEMENT SI ON NE FAIT PAS UNE SUITE
!===============================================================================

! RQ IMPORTANTE : pour la combustion CP, 1 seul passage suffit

if ( isuite.eq.0 .and. ipass.eq.1 ) then

! --> Initialisation de k et epsilon

  xkent = 1.d-10
  xeent = 1.d-10

! ---- TURBULENCE

  if (itytur.eq.2) then

    do iel = 1, ncel
      rtp(iel,ik)  = xkent
      rtp(iel,iep) = xeent
    enddo

  elseif (itytur.eq.3) then

    do iel = 1, ncel
      rtp(iel,ir11) = d2s3*xkent
      rtp(iel,ir22) = d2s3*xkent
      rtp(iel,ir33) = d2s3*xkent
      rtp(iel,ir12) = 0.d0
      rtp(iel,ir13) = 0.d0
      rtp(iel,ir23) = 0.d0
      rtp(iel,iep)  = xeent
    enddo

  elseif (iturb.eq.50) then

    do iel = 1, ncel
      rtp(iel,ik)   = xkent
      rtp(iel,iep)  = xeent
      rtp(iel,iphi) = d2s3
      rtp(iel,ifb)  = 0.d0
    enddo

  elseif (iturb.eq.60) then

    do iel = 1, ncel
      rtp(iel,ik)   = xkent
      rtp(iel,iomg) = xeent/cmu/xkent
    enddo

  elseif (iturb.eq.70) then

    do iel = 1, ncel
      rtp(iel,inusa) = cmu*xkent**2/xeent
    enddo

  endif

! --> On initialise tout le domaine de calcul avec de l'air a TINITK
!                   ================================================

! ---- Calculs de H1INIT et H2INIT

  t1init = t0
  t2init = t0

! ------ Variables de transport relatives a la phase solide :
!        calcul de H

  do icla = 1, nclacp
    icha = ichcor(icla)

    do iel = 1, ncel

      rtp(iel,isca(ixch(icla))) = zero
      rtp(iel,isca(ixck(icla))) = zero
      rtp(iel,isca(inp(icla) )) = zero
      rtp(iel,isca(ih2(icla) )) = zero
      if ( ippmod(icp3pl) .eq. 1 ) then
        rtp(iel,isca(ixwt(icla))) = zero
      endif

    enddo
  enddo

! ------ Variables de transport relatives au melange

  do ige = 1, ngazem
    coefe(ige) = zero
  enddo

!       On considere l'oxydant 1

  coefe(io2) = wmole(io2)*oxyo2(1)                                &
              /( wmole(io2) *oxyo2(1) +wmole(in2) *oxyn2(1)       &
                +wmole(ih2o)*oxyh2o(1)+wmole(ico2)*oxyco2(1))
  coefe(ih2o) = wmole(ih2o)*oxyh2o(1)                             &
              /( wmole(io2) *oxyo2(1) +wmole(in2) *oxyn2(1)       &
                +wmole(ih2o)*oxyh2o(1)+wmole(ico2)*oxyco2(1))
  coefe(ico2) = wmole(ico2)*oxyco2(1)                             &
              /( wmole(io2) *oxyo2(1) +wmole(in2) *oxyn2(1)       &
                +wmole(ih2o)*oxyh2o(1)+wmole(ico2)*oxyco2(1))
  coefe(in2) = 1.d0-coefe(io2)-coefe(ih2o)-coefe(ico2)

  do icha = 1, ncharm
    f1mc(icha) = zero
    f2mc(icha) = zero
  enddo
  mode = -1
  call cpthp1                                                     &
  !==========
 ( mode   , h1init , coefe  , f1mc   , f2mc   ,                   &
   t1init )

  do iel = 1, ncel
    rtp(iel,isca(ihm)) = h1init
  enddo

! ------ Variables de transport relatives au melange gazeux
!        (scalaires passifs et variances associees)

  do icha = 1, ncharb
    do iel = 1, ncel
      rtp(iel,isca(if1m(icha))) = zero
      rtp(iel,isca(if2m(icha))) = zero
    enddo
  enddo
  do iel = 1, ncel
    rtp(iel,isca(if3m)) = zero

    if ( ihtco2 .eq. 1 ) then
      rtp(iel,isca(if3mc2)) = zero
    endif

    rtp(iel,isca(if4p2m)) = zero
    if ( ippmod(icp3pl) .eq. 1 ) then
      rtp(iel,isca(if5m)) = zero
    endif

    if ( noxyd .ge. 2 ) then
      rtp(iel,isca(if6m)) = zero
    endif
    if ( noxyd .eq. 3 ) then
      rtp(iel,isca(if7m)) = zero
    endif

    if ( ieqco2.ge.1 ) then

      ioxy   =  1
      wmo2   = wmole(io2)
      wmco2  = wmole(ico2)
      wmh2o  = wmole(ih2o)
      wmn2   = wmole(in2)

      dmas = ( oxyo2 (ioxy)*wmo2 +oxyn2 (ioxy)*wmn2               &
              +oxyh2o(ioxy)*wmh2o+oxyco2(ioxy)*wmco2 )
      xco2 = oxyco2(ioxy)*wmco2/dmas

      rtp(iel,isca(iyco2)) = oxyco2(ioxy)*wmco2/dmas
    endif
  enddo

endif

!===============================================================================
! 2.  ON DONNE LA MAIN A L'UTILISATEUR
!===============================================================================

if (ipass.eq.1) then

  call cs_user_initialization &
  !==========================
( nvar   , nscal  ,                                            &
  dt     , rtp    , propce , propfa , propfb , coefa  , coefb  )

endif

!----
! FORMATS
!----

!----
! FIN
!----

return
end subroutine
