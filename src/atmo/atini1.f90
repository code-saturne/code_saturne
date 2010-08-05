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

subroutine atini1
!================


!===============================================================================
!  FONCTION  :
!  ---------

!   INIT DES OPTIONS DES VARIABLES POUR LA VERSION ATMOSPHERIQUE
!      EN COMPLEMENT DE CE QUI A DEJA ETE FAIT DANS USINI1

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

implicit none

!===============================================================================
! Common blocks
!===============================================================================

include "paramx.f90"
include "dimens.f90"
include "ihmpre.f90"
include "numvar.f90"
include "optcal.f90"
include "cstphy.f90"
include "entsor.f90"
include "cstnum.f90"
include "ppppar.f90"
include "ppthch.f90"
include "ppincl.f90"
include "atincl.f90"

!===============================================================================

! Local variables

integer          iphas, ii, isc, jj, ipp

!===============================================================================


!===============================================================================
! 0. VERIFICATION ISCALT, ISCSTH pour IPPMOD(IATMOS) = 1 or 2
!===============================================================================
!     L'utilisateur ne doit pas y avoir touche.

if ( ippmod(iatmos).ge.1 ) then
  do iphas = 1, nphas
    if(iscalt(iphas).ne.-1) then
      write(nfecra,1000)iphas,iscalt(iphas)
      call csexit (1)
      !==========
    endif
  enddo

  do ii = 1, nscapp
    if(iscsth(iscapp(ii)).ne.-10) then
    write(nfecra,1001)ii,iscapp(ii),iscapp(ii),iscsth(iscapp(ii))
     call csexit (1)
     !==========
    endif
  enddo
endif


!===============================================================================
! 1. INFORMATIONS GENERALES
!===============================================================================

! ---> Initialisation
iphas = 1
rair = 287.d0

! ---> Masse volumique et viscosite
irovar(iphas) = 0
ivivar(iphas) = 0

!===============================================================================
! 2. VARIABLES TRANSPORTEES pour IPPMOD(IATMOS) = 1 or 2
!===============================================================================

! 2.1  Dry atmosphere
! =====================

if ( ippmod(iatmos).eq.1 ) then

  iscsth(itempp) = 1
  iscalt(iphas) = itempp
  scamin(itempp)   = 0.d0
  scamax(itempp)   = +grand

!  for the dry atmosphere case, non constant density
  irovar(iphas) = 1

! --> Donnees physiques ou numeriques propres aux scalaires

  do isc = 1, nscapp

    jj = iscapp(isc)

    if (iscavr(jj).le.0) then
      visls0(jj) = viscl0(iphsca(jj))
    endif

    blencv(isca(jj)) = 1.d0

  enddo

  ipp = ipprtp(isca(itempp))
  nomvar(IPP)  = 'PotTemp'
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1

endif

! 2.2  Humid atmosphere
! =====================

if ( ippmod(iatmos).eq.2 ) then

  iscsth(itempl) = 1
  iscsth(itotwt) = 0
  iscsth(intdrp) = 0

  iscalt(iphas) = itempl

!  for the humid atmosphere case, non constant density
  irovar(iphas) = 1


  ipp = ipprtp(isca(itempl))
  nomvar(IPP)  = 'LqPotTmp'
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1

  ipp = ipprtp(isca(itotwt))
  nomvar(IPP)  = 'TotWater'
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1

  ipp = ipprtp(isca(intdrp))
  nomvar(IPP)  = 'TotDrop'
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1

endif

!===============================================================================
! 3. VARIABLES D'ETAT pour IPPMOD(IATMOS) = 1 or 2
!===============================================================================

! 3.1  Dry or humid atmosphere
! =============================

if (ippmod(iatmos).eq.1 .or. ippmod(iatmos).eq.2) then

  ipp = ipppro(ipproc(itempc))
  nomvar(IPP)   = 'TempC'
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

endif


!===============================================================================
! 4. One scale turbulent model for k-eps closure for IPPMOD(IATMOS) = 1 or 2
!===============================================================================

if (ippmod(iatmos).eq.1 .or. ippmod(iatmos).eq.2) then

  do iphas = 1, nphas
    if (itytur(iphas).eq.2) then
      ideuch(iphas) = 0
    endif
  enddo

endif

!===============================================================================
! 5. Turbulent Schmidt number for atmospheric flows
!===============================================================================

if (nscal.gt.0) then
  do ii = 1, nscal
    do iphas = 1, nphas
        sigmas(ii) = 0.7d0
    enddo
  enddo
endif

!===============================================================================
! 6. ON DONNE LA MAIN A L'UTLISATEUR
!===============================================================================

!   - Interface Code_Saturne
!     ======================

if (iihmpr.eq.1) then

  call uiati1 (imeteo)
  !==========

endif

call usati1
!==========

!--------
! FORMATS
!--------

 1000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    PHYSIQUE PARTICULIERE (ATMOSPHERIQUE) DEMANDEE          ',/,&
'@                                                            ',/,&
'@  La valeur de ISCALT est renseignee automatiquement.       ',/,&
'@                                                            ',/,&
'@  L''utilisateur ne doit pas la renseigner dans usini1, or  ',/,&
'@    elle a ete affectee comme suit :                        ',/,&
'@    ISCALT(',I10   ,') = ',I10                               ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier usini1.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 1001 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    PHYSIQUE PARTICULIERE (ATMOSPHERIQUE) DEMANDEE          ',/,&
'@                                                            ',/,&
'@  Les valeurs de ISCSTH sont renseignees automatiquement.   ',/,&
'@                                                            ',/,&
'@  L''utilisateur ne doit pas les renseigner dans usini1, or ',/,&
'@    pour le scalaire ',I10   ,' correspondant au scalaire   ',/,&
'@    physique particuliere ',I10   ,' on a                   ',/,&
'@    ISCSTH(',I10   ,') = ',I10                               ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier usini1.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

return
end subroutine
