!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2015 EDF S.A.
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

 ( npar1  , npar2  ,                                              &
   iprev  )

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
! npar1 ,npar2     ! e  ! <-- ! borne min et max des particules                !
!                  !    !     !    a initialiser                               !
! iprev            ! e  ! <-- ! time step indicator for fields                 !
!                  !    !     !   0: use fields at current time step           !
!                  !    !     !   1: use fields at previous time step          !
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

integer          npar1 , npar2
integer          iprev

! Local variables

integer          iel , npt , nomb , nbfac , ifac , il , izone
double precision tu , d2s3

double precision, allocatable, dimension(:) :: w1

double precision  lvisq
double precision  px , py , pz , distp , d1
double precision  dismin,dismax, ustar, visccf, romf
double precision  unif1(1)

double precision, dimension(:), pointer :: cromf
double precision, dimension(:,:), pointer :: vel
double precision, dimension(:), pointer :: cvar_k
double precision, dimension(:), pointer :: cvar_r11, cvar_r22, cvar_r33
double precision, dimension(:), pointer :: viscl
double precision, dimension(:,:), allocatable :: vagaus

!===============================================================================

! Map field arrays
if (iprev.eq.0) then
  call field_get_val_v(ivarfl(iu), vel)

  if (itytur.eq.2 .or. iturb.eq.50                              &
       .or. iturb.eq.60) call field_get_val_s(ivarfl(ik), cvar_k)
  if (itytur.eq.3) then
    call field_get_val_s(ivarfl(ir11), cvar_r11)
    call field_get_val_s(ivarfl(ir22), cvar_r22)
    call field_get_val_s(ivarfl(ir33), cvar_r33)
  endif

else if (iprev.eq.1) then
  call field_get_val_prev_v(ivarfl(iu), vel)

  if (itytur.eq.2 .or. iturb.eq.50                              &
       .or. iturb.eq.60) call field_get_val_prev_s(ivarfl(ik), cvar_k)
  if (itytur.eq.3) then
    call field_get_val_prev_s(ivarfl(ir11), cvar_r11)
    call field_get_val_prev_s(ivarfl(ir22), cvar_r22)
    call field_get_val_prev_s(ivarfl(ir33), cvar_r33)
  endif
endif

call field_get_val_s(iprpfl(iviscl), viscl)

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

  if (itytur.eq.2 .or. iturb.eq.50                     &
       .or. iturb.eq.60) then
    do iel = 1,ncel
      w1(iel) = cvar_k(iel)
    enddo
  else if (itytur.eq.3) then
    do iel = 1,ncel
      w1(iel) = 0.5d0 * ( cvar_r11(iel)                    &
                        + cvar_r22(iel)                    &
                        + cvar_r33(iel) )
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

nomb = npar2-npar1+1
allocate(vagaus(nomb,3))

if (idistu.eq.1 .and. nomb.gt.0) then
  call normalen (nomb,vagaus(:,1))
  call normalen (nomb,vagaus(:,2))
  call normalen (nomb,vagaus(:,3))
else
  do npt = 1,nomb
    vagaus(npt,1) = 0.d0
    vagaus(npt,2) = 0.d0
    vagaus(npt,3) = 0.d0
  enddo
endif

do npt = npar1,npar2

  iel = ipepa(jisor,npt)

  tu = sqrt(d2s3*w1(iel))

  eptp(juf,npt) = vel(1,iel) + vagaus(npt-npar1+1,1)*tu
  eptp(jvf,npt) = vel(2,iel) + vagaus(npt-npar1+1,2)*tu
  eptp(jwf,npt) = vel(3,iel) + vagaus(npt-npar1+1,3)*tu

enddo

deallocate(vagaus)

!Calcul de la fluctuation de vitesse si le modèle de dépôt est activé

if (idepst.eq.1) then

  do npt = npar1,npar2
    iel = ipepa(jisor,npt)
    !Calculation of the normalized wall-normal particle distance (y^+)

    dismin = 1.d+20
    dismax = -1.d+20

    distp = 1.0d+20
    pepa(jryplu,npt) = 1.0d3
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
          visccf = viscl(iel) / romf

          ustar = uetbor(ifac)
          lvisq = visccf / ustar

          px = eptp(jxp,npt)
          py = eptp(jyp,npt)
          pz = eptp(jzp,npt)

          d1 = abs(px*dlgeo(ifac,1)+py*dlgeo(ifac,2)         &
                  +pz*dlgeo(ifac,3)+dlgeo(ifac,4))           &
               /sqrt( dlgeo(ifac,1)*dlgeo(ifac,1)            &
                  +dlgeo(ifac,2)*dlgeo(ifac,2)               &
                  +dlgeo(ifac,3)*dlgeo(ifac,3))

          if (d1.lt.distp) then
            distp = d1
            pepa(jryplu,npt) = distp/lvisq
            ipepa(jdfac,npt) = ifac
          endif
        endif
      endif
    enddo

    if (pepa(jryplu,npt).lt.pepa(jrinpf,npt)) then
      ipepa(jimark,npt) = 10
    elseif (pepa(jryplu,npt).gt.100.d0) then
      ipepa(jimark,npt) = -1
    else

      call zufall(1, unif1(1))
      !==========
      if (unif1(1).lt.0.25d0) then
        ipepa(jimark,npt) = 12
      elseif (unif1(1).gt.0.625d0) then
        ipepa(jimark,npt) = 1
      elseif ((unif1(1).gt.0.25d0).and.(unif1(1).lt.0.625d0)) then
        ipepa(jimark,npt) = 3
      endif
    endif

    if (pepa(jryplu,npt).le.pepa(jrinpf,npt)) then
      eptp(juf,npt) = vel(1,iel)
      eptp(jvf,npt) = vel(2,iel)
      eptp(jwf,npt) = vel(3,iel)
    endif

    ! No deposited particles at the injection
    ipepa(jdepo,npt) = 0

    ! Initialization of the additional "pointers"
    ! for the resuspension model

    if (ireent.gt.0) then

      pepa(jfadh,npt) = 0.d0
      pepa(jmfadh,npt) = 0.d0

      ipepa(jnbasg,npt) = 0
      ipepa(jnbasp,npt) = 0

      pepa(jndisp,npt) = 0.d0

   endif

 enddo

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

return
end subroutine lagipn
