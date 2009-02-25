!-------------------------------------------------------------------------------

!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2008 EDF S.A., France

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

subroutine fuflux &
!================

 ( idbia0 , idbra0 ,                                              &
   ncelet , ncel   ,                                              &
   rtpa   , propce , volume ,                                     &
   w1     , w2     , w3     ,                                     &
   ra     )

!===============================================================================
! FONCTION :
! --------

! CALCUL DES TERMES DE TRANSFERT DE MASSE ENTRE LA PHASE CONTINUE
! ET LA PHASE DISPERSEE


! Arguments
!__________________.____._____.________________________________________________.
!    nom           !type!mode !                   role                         !
!__________________!____!_____!________________________________________________!
! idbia0           ! e  ! <-- ! numero de la 1ere case libre dans ia           !
! idbra0           ! e  ! <-- ! numero de la 1ere case libre dans ra           !
! ncelet           ! e  ! <-- ! nombre d'elements halo compris                 !
! ncel             ! e  ! <-- ! nombre d'elements actifs                       !
! rtpa             ! tr ! <-- ! variables de calcul au centre des              !
! (ncelet,*)       !    !     !    cellules (instant precedent)                !
! propce           ! tr ! <-- ! proprietes physiques au centre des             !
! (ncelet,*)       !    !     !    cellules                                    !
! volume           ! tr ! <-- ! volume d'un des ncelet elements                !
! (ncelet          !    !     !                                                !
! w1, w2, w3       ! tr ! --- ! tableaux de travail                            !
! ra(*)            ! tr ! --- ! macro tableau reel                             !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

implicit none

!===============================================================================
!     DONNEES EN COMMON
!===============================================================================

include "paramx.h"
include "numvar.h"
include "optcal.h"
include "cstphy.h"
include "cstnum.h"
include "entsor.h"
include "parall.h"
include "ppppar.h"
include "ppthch.h"
include "coincl.h"
include "cpincl.h"
include "fuincl.h"
include "ppincl.h"

!===============================================================================

! Arguments

integer          idbia0 , idbra0
integer          ncelet , ncel

double precision rtpa(ncelet,*), propce(ncelet,*)
double precision w1(ncelet), w2(ncelet), w3(ncelet)
double precision volume(ncelet)
double precision ra(*)

! VARIABLES LOCALES

integer          idebia , idebra
integer          iel    , iphas  , icla
integer          ipcrom , ipcte1 , ipcte2 , ipcro2 , ipcdia
integer          ipcgev , ipcght , ipcyox
integer          ipcvst,ipcvsl,ipccp,ipchgl

double precision xng,xnuss
double precision pparo2 , xdffli , xdfext , xdftot0 , xdftot1
double precision diacka, xuash
double precision dcoke , surf

!===============================================================================
! 1. INITIALISATIONS ET CALCULS PRELIMINAIRES
!===============================================================================

! --- Initialisation memoire

idebia = idbia0
idebra = idbra0

! --- Initialisation des termes de transfert de masse

do icla = 1, nclafu
  ipcgev = ipproc(igmeva(icla))
  ipcght = ipproc(igmhtf(icla))
  do iel = 1, ncel
    propce(iel,ipcgev) = zero
    propce(iel,ipcght) = zero
  enddo
enddo

! --- Initialisation des tableaux de travail

do iel = 1, ncel
  w1(iel) = zero
  w2(iel) = zero
  w3(iel) = zero
enddo

! --- Calcul de la masse volumique du melange gazeux

iphas  = 1
ipcrom = ipproc(irom(iphas))
ipcte1 = ipproc(itemp1)
ipcyox = ipproc(iym1(io2))

! --- Numero des grandeurs physiques (voir usclim)
ipcrom = ipproc(irom(iphas))
ipcvst = ipproc(ivisct(iphas))

! --> Terme source pour l'enthalpie du liquide

do icla = 1, nclafu

  ipcro2 = ipproc(irom3 (icla))
  ipcdia = ipproc(idiam3(icla))
  ipcte2 = ipproc(itemp3(icla))
  ipcght = ipproc(igmhtf(icla))
  ipchgl = ipproc(ih1hlf(icla))


! ---- Contribution aux bilans explicite et implicite
!        des echanges par diffusion moleculaire
!        6 Lambda Nu / diam**2 / Rho2 * Rho * (T1-T2)
! ------ Calcul de lambda dans W1

  xnuss = 2.d0
  do iel = 1, ncel
    if ( ivisls(ihm).gt.0 ) then
      ipcvsl = ipproc(ivisls(ihm))
      if ( icp(iphas).gt.0 ) then
        ipccp   = ipproc(icp(iphas))
        w1(iel) = propce(iel,ipcvsl) * propce(iel,ipccp)
      else
        w1(iel) = propce(iel,ipcvsl) * cp0(iphas)
      endif
    else
      if ( icp(iphas).gt.0 ) then
        ipccp   = ipproc(icp(iphas))
        w1(iel) = visls0(ihm) * propce(iel,ipccp)
      else
        w1(iel) = visls0(ihm) * cp0(iphas)
      endif
    endif
  enddo

!----Contribution aux bilans explicite et implicite des
!    echanges par diffusion moleculaire
!        6 Lambda Nu / diam**2 / Rho2 * Rho
!    le diametre est en mm  donc on multiplie par 1.D-3
!    pour l'avoir en m

  do iel = 1, ncel


    if ( rtpa(iel,isca(ing(icla))) .gt. epsifl .and.              &
           propce(iel,ipcte1).gt. propce(iel,ipcte2) ) then


!          PROPCE(IEL,IPCHGL) = 6.D0 * W1(IEL) * XNUSS
!     &                         /((PROPCE(IEL,IPCDIA)*1.D-3)**2)
!     &                         /PROPCE(IEL,IPCRO2)

      propce(iel,ipchgl) = w1(iel)*xnuss*rtpa(iel,isca(ing(icla)))&
                          *pi*propce(iel,ipcdia)*1.d6

    else
      propce(iel,ipchgl) =0.d0
    endif

  enddo

enddo

!===============================================================================
! 2. TRANSFERTS DE MASSE PAR EVAPORATION
!===============================================================================

do icla = 1, nclafu

  ipcro2 = ipproc(irom3 (icla))
  ipcdia = ipproc(idiam3(icla))
  ipcte2 = ipproc(itemp3(icla))
  ipcgev = ipproc(igmeva(icla))
  ipchgl = ipproc(ih1hlf(icla))

! FO & PPl 16/09/05  Version avec parametres en dur, en attendant les
!                    lectures et inclusion en include
  tevap1 = 150.d0 + tkelvi
  tevap2 = 450.d0 + tkelvi


! PPl 161205 On répartit le flux interfacial d'enthalpie entre
!            l'échauffement de la goutte et l'évaporation

  do iel = 1, ncel

! --- Transfert de masse du a l'evaporation

    propce(iel,ipcgev) = zero

!         IF (PROPCE(IEL,IPCTE2)    .GE. TEVAP1 .AND.
!    &       PROPCE(IEL,IPCTE2)    .LE. TEVAP2 .AND.
!    &       RTPA(IEL,ISCA(IYFOL(ICLA))) .GE. EPSIFL      ) THEN
!          PROPCE(IEL,IPCGEV) = PROPCE(IEL,IPPROC(IH1HLF)) /
!     &    ( RTPA(IEL,ISCA(IYFOL)) * CP2FOL * (TEVAP2-TEVAP1)  /
!     &      ( ( RTPA(IEL,ISCA(IYFOL)) + RTPA(IEL,ISCA(IFVAP)) )
!     &    * (1    .d0-FKC))                                     + HRFVAP )
! PPl 161205

! PPl 090106 On teste plutôt les diamètres car désormais les
!        gouttes finissent trop fines

    if ( propce(iel,ipcte2)    .gt. tevap1 .and.                  &
         propce(iel,ipcdia)    .gt. dinikf(icla) .and.            &
         rtpa(iel,isca(iyfol(icla))) .gt. epsifl       ) then
      propce(iel,ipcgev) = propce(iel,ipchgl)                     &
                        /( hrfvap + cp2fol*(tevap2-tevap1) )
    endif

  enddo

enddo

!===============================================================================
! 3. CALCUL DE RHO_COKE MOYEN
!    On suppose pour le calcul de la masse volumique du coke que
!    l'evaporation a lieu a volume de coke constant
!    Par la suite, on suppose (pour commencer ?) que la
!    combustion hétérogène a lieu à volume constant =>
!    à masse volumique décroissante
!===============================================================================

! --- Initialisation

rhokf  = rho0fl

! --- Calcul de la masse volumique moyenne du coke

!===============================================================================
! 4. TRANSFERTS DE MASSE PAR COMBUSTION HETEROGENE
!===============================================================================

do icla = 1, nclafu

  ipcro2 = ipproc(irom3 (icla))
  ipcdia = ipproc(idiam3(icla))
  ipcte2 = ipproc(itemp3(icla))
  ipcgev = ipproc(igmeva(icla))
  ipcght = ipproc(igmhtf(icla))

  do iel = 1, ncel

    if ( propce(iel,ipcdia) .le. dinikf(icla) .and.               &
         propce(iel,ipcdia) .gt. diniin(icla) .and.               &
         rtpa(iel,isca(iyfol(icla))) .gt. epsifl ) then

      xng   = rtpa(iel,isca(ing(icla)))*1.d9

      xuash = xng*(1.d0-xinfol)                                   &
             *(rho0fl*pi*((dinifl(icla)*1.d-3)**2))/6.d0

! --- Calcul de la pression partielle en oxygene (atm)
!                                                 ---
!       PO2 = RHO1*RR*T*YO2/MO2

      pparo2 = propce(iel,ipproc(irom1))*rr*propce(iel,ipcte1)    &
              *propce(iel,ipcyox)/wmole(io2)
      pparo2 = pparo2 / prefth

! --- Calcul de Dcoke en metres

      dcoke = ( ( rtpa(iel,isca(iyfol(icla)))                     &
              /(rtpa(iel,isca(ing(icla)))*rho0fl)                 &
              -pi*(dinikf(icla)**3)*xinkf/6.d0  )                 &
               *6.d0/(pi*(1.d0-xinkf)) )**(1.d0/3.d0)             &
            *1.d-3
     if ( dcoke .lt. 0.d0 ) then
       WRITE(NFECRA,*) 'erreur Dcoke = ',Dcoke,IEL
       call csexit(1)
     endif

! --- Coefficient de cinetique chimique de formation de CO
!       en (kg.m-2.s-1.atm(-n))

      xdffli = ahetfl*exp(-ehetfl*4185.d0                         &
              /(rr*propce(iel,ipcte1)))

! --- Coefficient de diffusion en  (Kg/m2/s/atm) : XDFEXT
!     Coefficient global pour n=0.5 en (kg/m2/s) : XDFTOT0
!     Coefficient global pour n=1   en (Kg/m2/s) : XDFTOT1

      diacka = dcoke/(dinikf(icla)*1.d-3)
      if ( diacka .gt. epsifl ) then
        xdfext = 2.53d-7*((propce(iel,ipcte1))**0.75d0)           &
                / dcoke*2.d0
        xdftot1 = pparo2 / ( 1.d0/xdffli + 1.d0/xdfext )
        xdftot0 = -(xdffli**2)/(2.d0*xdfext**2)+(pparo2*xdffli**2 &
                  +(xdffli**4)/(2.d0*xdfext**2))**0.5d0
      else
        xdftot1 = xdffli*pparo2
        xdftot0 = xdffli*pparo2**0.5d0
      endif

!     Surface

      surf = pi*(dcoke**2)

! --- Calcul de PROPCE(IEL,IPCGHT) = - COXCK*XDFTOT0*PPARO2*XNP < 0
! --- ou        PROPCE(IEL,IPCGHT) = - COXCK*XDFTOT1*PPARO2*XNP < 0

      if (iofhet.eq.1) then
!             PROPCE(IEL,IPCGHT) = - XDFTOT1*COXCK*XNG
        propce(iel,ipcght) = - xdftot1*surf*xng
      else
!             PROPCE(IEL,IPCGHT) = - XDFTOT0*COXCK*XNG
        propce(iel,ipcght) = - xdftot0*surf*xng
      endif

    else
      propce(iel,ipcght) = 0.d0
    endif

  enddo

enddo

!===============================================================================
! FORMATS
!--------
!----
! FIN
!----

return
end
