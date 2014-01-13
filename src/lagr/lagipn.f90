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

subroutine lagipn &
!================

 ( nbpmax ,                                                       &
   npar1  , npar2  ,                                              &
   rtp    ,                                                       &
   vagaus ,                                                       &
   propce )

!===============================================================================
! FONCTION :
! ----------

!   SOUS-PROGRAMME DU MODULE LAGRANGIEN :
!   -------------------------------------

!     Initialisation de la vitesse fluide vu pour les nouvelles
!     particules.

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nbpmax           ! e  ! <-- ! nombre max de particulies autorise             !
! npar1 ,npar2     ! e  ! <-- ! borne min et max des particules                !
!                  !    !     !    a initialiser                               !
! rtp              ! tr ! <-- ! variables de calcul au centre des              !
! (ncelet,*)       !    !     !    cellules (instant courant ou prec)          !
! vagaus           ! tr ! --> ! variables aleatoires gaussiennes               !
!(nbpmax,nvgaus    !    !     !                                                !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
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
use cstnum
use numvar
use optcal
use entsor
use cstphy
use pointe
use parall
use period
use lagpar
use lagran
use ppincl
use mesh
use field

!===============================================================================

implicit none

! Arguments

integer          nbpmax
integer          npar1 , npar2

double precision rtp(ncelet,*)
double precision propce(ncelet,*)
double precision vagaus(nbpmax,*)

! Local variables

integer          iel , npt , nomb , nbfac , ifac , il , izone
double precision tu , d2s3

double precision, allocatable, dimension(:) :: w1

double precision  lvisq
double precision  px , py , pz , distp , d1
double precision  dismin,dismax, ustar, visccf, romf
double precision  unif1(1)

double precision, dimension(:), pointer :: cromf

!===============================================================================

!===============================================================================
! 1. INITIALISATION
!===============================================================================

d2s3 = 2.d0 / 3.d0

! Allocate a work array
allocate(w1(ncelet))

!===============================================================================
! 2. SIMULATION DES VITESSES TURBULENTES FLUIDES INSTANTANNEES VUES
!    PAR LES PARTICULES SOLIDES LE LONG DE LEUR TRAJECTOIRE.
!===============================================================================

if (idistu.eq.1) then

  if (itytur.eq.2 .or. iturb.eq.50                  &
       .or. iturb.eq.60) then
    do iel = 1,ncel
      w1(iel) = rtp(iel,ik)
    enddo
  else if (itytur.eq.3) then
    do iel = 1,ncel
      w1(iel) = 0.5d0 * ( rtp(iel,ir11)                    &
                        + rtp(iel,ir22)                    &
                        + rtp(iel,ir33) )
    enddo
  else
    write(nfecra,9000) iilagr, idistu, iturb
    call csexit (1)
    !==========
  endif
else
  do iel = 1,ncel
    w1(iel) = 0.d0
  enddo
endif


!---> CALCUL DES TIRAGES ALEATOIRES
!     CALCUL DU TEMPS CARACTERISTIQUE DES PARTICULES
!     remarque : NORMALEN est dans le fichier ZUFALL.F
!     ^^^^^^^^

if (idistu.eq.1) then
  nomb = npar2-npar1+1
  call normalen (nomb,vagaus(npar1,1))
  call normalen (nomb,vagaus(npar1,2))
  call normalen (nomb,vagaus(npar1,3))
else
  do npt = npar1,npar2
    vagaus(npt,1) = 0.d0
    vagaus(npt,2) = 0.d0
    vagaus(npt,3) = 0.d0
  enddo
endif

do npt = npar1,npar2

  iel = itepa(npt,jisor)

  tu = sqrt( d2s3*w1(iel) )

  ettp(npt,juf) = rtp(iel,iu) + vagaus(npt,1)*tu
  ettp(npt,jvf) = rtp(iel,iv) + vagaus(npt,2)*tu
  ettp(npt,jwf) = rtp(iel,iw) + vagaus(npt,3)*tu

enddo

!Calcul de la fluctuation de vitesse si le modèle de dépôt est activé

if (idepst.eq.1) then

  do npt = npar1,npar2
    iel = itepa(npt,jisor)
    !Calculation of the normalized wall-normal particle distance (y^+)

    dismin = 1.d+20
    dismax = -1.d+20

    distp = 1.0d+20
    tepa(npt,jryplu) = 1.0d3
    nbfac = 0

    do il = itycel(iel)+1,itycel(iel+1)
      ifac = icocel(il)
      if (ifac.lt.0) then
        ifac = -ifac
        izone = ifrlag(ifac)

        ! Test if the particle is located in a boundary cell

        if ( iusclb(izone) .eq. idepo1 .or.                   &
             iusclb(izone) .eq. idepo2 .or.                   &
             iusclb(izone) .eq. idepfa .or.                   &
             iusclb(izone) .eq. irebol     ) then


          ! Calculation of the wall units

          if (ippmod(iccoal).ge.0 .or. ippmod(icfuel).ge.0) then
            call field_get_val_s(iprpfl(ipproc(irom1)), cromf)
          else
            call field_get_val_s(icrom, cromf)
          endif

          romf = cromf(iel)
          visccf = propce(iel,ipproc(iviscl)) / romf

          ustar = uetbor(ifac)
          lvisq = visccf / ustar

          px = ettp(npt,jxp)
          py = ettp(npt,jyp)
          pz = ettp(npt,jzp)

          d1 = abs(px*dlgeo(ifac,1)+py*dlgeo(ifac,2)         &
                  +pz*dlgeo(ifac,3)+dlgeo(ifac,4))           &
               /sqrt( dlgeo(ifac,1)*dlgeo(ifac,1)            &
                  +dlgeo(ifac,2)*dlgeo(ifac,2)               &
                  +dlgeo(ifac,3)*dlgeo(ifac,3))

          if (d1.lt.distp) then
            distp = d1
            tepa(npt,jryplu) = distp/lvisq
            itepa(npt,jdfac) = ifac
          endif
        endif
      endif
    enddo

    if (tepa(npt,jryplu).lt.tepa(npt,jrinpf)) then
      itepa(npt,jimark) = 10
    elseif (tepa(npt,jryplu).gt.100.d0) then
      itepa(npt,jimark) = -1
    else

      call zufall(1, unif1(1))
      !==========
      if (unif1(1).lt.0.25d0) then
        itepa(npt,jimark) = 12
      elseif (unif1(1).gt.0.625d0) then
        itepa(npt,jimark) = 1
      elseif ((unif1(1).gt.0.25d0).and.(unif1(1).lt.0.625d0)) then
        itepa(npt,jimark) = 3
      endif
    endif

    if (tepa(npt,jryplu).le.tepa(npt,jrinpf)) then
      ettp(npt,juf) = rtp(iel,iu)
      ettp(npt,jvf) = rtp(iel,iv)
      ettp(npt,jwf) = rtp(iel,iw)
    endif

    ! No deposited particles at the injection
    itepa(npt,jdepo) = 0

  enddo

  ! Initialization of the additional "pointers"
  ! for the resuspension model

  if (ireent.gt.0) then

    tepa(npt,jfadh) = 0.d0
    tepa(npt,jmfadh) = 0.d0

    itepa(npt,jnbasg) = 0
    itepa(npt,jnbasp) = 0

    tepa(npt,jndisp) = 0.d0

  endif

endif

! Free memory
deallocate(w1)

!--------
! FORMATS
!--------

 9000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========   (LAGIPN)                                    ',/,&
'@                                                            ',/,&
'@    LE MODULE LAGRANGIEN EST INCOMPATIBLE AVEC LE MODELE    ',/,&
'@    DE TURBULENCE SELECTIONNE.                              ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@   Le module Lagrangien a ete active avec IILAGR = ',I10     ,/,&
'@     et la dispersion turbulente est prise en compte        ',/,&
'@                                     avec IDISTU = ',I10     ,/,&
'@   Le modele de turbulence                                  ',/,&
'@     correspond a ITURB  = ',I10                             ,/,&
'@   Or, les seuls traitements de la turbulence compatibles   ',/,&
'@     avec le module Lagrangien et la dispersion turbulente  ',/,&
'@     sont k-epsilon et Rij-epsilon, v2f et k-omega.         ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier la valeur de IILAGR et IDISTU dans la subroutine ',/,&
'@  USLAG1 et verifier la valeur de ITURB  dans la subroutine ',/,&
'@  USINI1.                                                   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

!----
! FIN
!----

end subroutine
