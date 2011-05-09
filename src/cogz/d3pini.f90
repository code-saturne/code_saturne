!-------------------------------------------------------------------------------

!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2010 EDF S.A., France

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

subroutine d3pini &
!================

 ( idbia0 , idbra0 ,                                              &
   nvar   , nscal  , nphas  ,                                     &
   maxelt , lstelt ,                                              &
   ia     ,                                                       &
   dt     , rtp    , propce , propfa , propfb , coefa  , coefb  , &
   ra     )

!===============================================================================
! FONCTION :
! --------

! INITIALISATION DES VARIABLES DE CALCUL
!    POUR LA PHYSIQUE PARTICULIERE :
!        COMBUSTION GAZ - FLAMME DE DIFFUSION CHIMIE 3 POINTS
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

! Les proprietes physiques sont accessibles dans le tableau
!     PROPCE (prop au centre), PROPFA (aux faces internes),
!     PROPFB (prop aux faces de bord)
!     Ainsi,
!      PROPCE(IEL,IPPROC(IROM  )) designe ROM   (IEL ,IPHAS)
!      PROPCE(IEL,IPPROC(IVISCL)) designe VISCL (IEL ,IPHAS)
!      PROPCE(IEL,IPPROC(ICP   )) designe CP    (IEL ,IPHAS)
!      PROPCE(IEL,IPPROC(IVISLS(ISCAL))) designe VISLS (IEL ,ISCAL)

!      PROPFA(IFAC,IPPROF(IFLUMA(IVAR ))) designe FLUMAS(IFAC,IVAR)

!      PROPFB(IFAC,IPPROB(IROM  )) designe ROMB  (IFAC,IPHAS)
!      PROPFB(IFAC,IPPROB(IFLUMA(IVAR ))) designe FLUMAB(IFAC,IVAR)

! LA MODIFICATION DES PROPRIETES PHYSIQUES (ROM, VISCL, VISCLS, CP)
!     SE FERA EN STANDARD DANS LE SOUS PROGRAMME PPPHYV
!     ET PAS ICI

! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! nphas            ! i  ! <-- ! number of phases                               !
! maxelt           ! i  ! <-- ! max number of cells and faces (int/boundary)   !
! lstelt(maxelt)   ! ia ! --- ! work array                                     !
! ia(*)            ! ia ! --- ! main integer work array                        !
! dt(ncelet)       ! tr ! <-- ! valeur du pas de temps                         !
! rtp              ! tr ! <-- ! variables de calcul au centre des              !
! (ncelet,*)       !    !     !    cellules                                    !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! coefa coefb      ! tr ! <-- ! conditions aux limites aux                     !
!  (nfabor,*)      !    !     !    faces de bord                               !
! ra(*)            ! ra ! --- ! main real work array                           !
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
use parall
use period
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use mesh

!===============================================================================

implicit none

integer          idbia0 , idbra0
integer          nvar   , nscal  , nphas

integer          maxelt, lstelt(maxelt)
integer          ia(*)

double precision dt(ncelet), rtp(ncelet,*), propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision ra(*)

! Local variables

character*80     chaine
integer          idebia, idebra
integer          iel, igg, iphas, mode
integer          iscal, ivar, ii
double precision coefg(ngazgm), hair, tinitk
double precision valmax, valmin

! NOMBRE DE PASSAGES DANS LA ROUTINE

integer          ipass
data             ipass /0/
save             ipass

!===============================================================================
!===============================================================================
! 1.  INITIALISATION VARIABLES LOCALES
!===============================================================================

ipass = ipass + 1

idebia = idbia0
idebra = idbra0

do igg = 1, ngazgm
  coefg(igg) = zero
enddo

iphas    = 1


!===============================================================================
! 2. INITIALISATION DES INCONNUES :
!      UNIQUEMENT SI ON NE FAIT PAS UNE SUITE
!===============================================================================

if ( isuite.eq.0 ) then

! ---> Initialisation au 1er passage avec de l'air a TINITK
!                                    ======================

  if ( ipass.eq.1 ) then

! ----- Calcul de l'enthalpie de l'air HAIR a TINITK

    tinitk   = t0
    coefg(1) = zero
    coefg(2) = 1.d0
    coefg(3) = zero
    mode     = -1
    call cothht                                                   &
    !==========
  ( mode   , ngazg , ngazgm  , coefg  ,                           &
    npo    , npot   , th     , ehgazg ,                           &
    hair   , tinitk )

! ----- On en profite pour initialiser TINFUE et TINOXY
!       et calculer HINFUE et HINOXY a partir de TINITK et HAIR
!       CAR on n'a pas encore vu usd3pc.F

    tinoxy = tinitk
    tinfue = tinitk
    hinoxy = hair
    hinfue = hair

    do iel = 1, ncel

! ----- Moyenne et variance du taux de melange

      rtp(iel,isca(ifm))   = zero
      rtp(iel,isca(ifp2m)) = zero

! ----- Enthalpie

      if ( ippmod(icod3p).eq.1 ) then
        rtp(iel,isca(ihm)) = hair
      endif

    enddo

! ---> Initialisation au 2eme passage

  else if ( ipass.eq.2 ) then

    do iel = 1, ncel

! ----- Moyenne et variance du taux de melange

      rtp(iel,isca(ifm))   = fs(1)
      rtp(iel,isca(ifp2m)) = zero

! ----- Enthalpie

      if ( ippmod(icod3p).eq.1 ) then
        rtp(iel,isca(ihm)) = hinfue*fs(1)+hinoxy*(1.d0-fs(1))
      endif

    enddo

! ----- On donne la main a l'utilisateur

    call usd3pi                                                   &
    !==========
 ( idebia , idebra ,                                              &
   nvar   , nscal  , nphas  ,                                     &
   maxelt , lstelt ,                                              &
   ia     ,                                                       &
   dt     , rtp    , propce , propfa , propfb , coefa  , coefb  , &
   ra     )

! ----- En periodique et en parallele,
!       il faut echanger ces initialisations

    if (irangp.ge.0.or.iperio.eq.1) then
      call synsca(rtp(1,isca(ifm)))
      !==========
      call synsca(rtp(1,isca(ifp2m)))
      !==========
      if ( ippmod(icod3p).eq.1 ) then
        call synsca(rtp(1,isca(ihm)))
        !==========
      endif
    endif


!      Impressions de controle

    write(nfecra,2000)

    do ii  = 1, nscapp
      iscal = iscapp(ii)
      ivar  = isca(iscal)
      valmax = -grand
      valmin =  grand
      do iel = 1, ncel
        valmax = max(valmax,rtp(iel,ivar))
        valmin = min(valmin,rtp(iel,ivar))
      enddo
      chaine = nomvar(ipprtp(ivar))
      if (irangp.ge.0) then
        call parmin(valmin)
        !==========
        call parmax(valmax)
        !==========
      endif
      write(nfecra,2010)chaine(1:8),valmin,valmax
    enddo

    write(nfecra,2020)

  endif

endif


!----
! FORMATS
!----

 2000 format(                                                           &
'                                                             ',/,&
' ----------------------------------------------------------- ',/,&
'                                                             ',/,&
'                                                             ',/,&
' ** INITIALISATION DES VARIABLES PROPRES AU GAZ (FL DIF 3PT) ',/,&
'    -------------------------------------------------------- ',/,&
'           2eme PASSAGE                                      ',/,&
' ---------------------------------                           ',/,&
'  Variable  Valeur min  Valeur max                           ',/,&
' ---------------------------------                           '  )

 2010 format(                                                           &
 2x,     a8,      e12.4,      e12.4                              )

 2020 format(                                                           &
' ---------------------------------                           ',/)


!----
! FIN
!----

return
end subroutine
