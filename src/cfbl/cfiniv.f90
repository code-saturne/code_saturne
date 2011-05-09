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

subroutine cfiniv &
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
!    POUR LA PHYSIQUE PARTICULIERE : COMPRESSIBLE SANS CHOC
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

integer          idebia, idebra
integer          ifinia, ifinra
integer          iwcel1, iwcel2, iwcel3, iwcel4
integer          iel   , iphas , iccfth, imodif
integer          iirom , iiromb, ifac


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

! --- Reservation de la memoire pour appel à uscfth

call memcfv                                                       &
!==========
 ( idebia , idebra ,                                              &
   iwcel1 , iwcel2 , iwcel3 , iwcel4 ,                            &
   ifinia , ifinra )

idebia = ifinia
idebra = ifinra

!===============================================================================
! 2. INITIALISATION DES INCONNUES :
!      UNIQUEMENT SI ON NE FAIT PAS UNE SUITE
!===============================================================================

if ( isuite.eq.0 ) then

  if ( ipass.eq.1 ) then

! ----- Initialisations par defaut

!     ON MET LA TEMPERATURE A T0
    do iel = 1, ncel
      rtp(iel,isca(itempk)) = t0
    enddo

!     On initialise Cv, rho et l'energie
    iccfth = 0
    imodif = 1

    call uscfth                                                 &
    !==========
 ( nvar   , nscal  , nphas  ,                                     &
   iccfth , imodif , iphas  ,                                     &
   dt     , rtp    , rtp    , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   ra(iwcel1), ra(iwcel2), ra(iwcel3), ra(iwcel4) )

!     On initialise la diffusivite thermique
    visls0(ienerg) = visls0(itempk)/cv0

    if(ivisls(ienerg).gt.0) then
      if(ivisls(itempk).gt.0) then
        if(icv.gt.0) then
          do iel = 1, ncel
            propce(iel,ipproc(ivisls(ienerg))) =         &
                 propce(iel,ipproc(ivisls(itempk)))        &
                 / propce(iel,ipproc(icv))
          enddo
        else
          do iel = 1, ncel
            propce(iel,ipproc(ivisls(ienerg))) =         &
                 propce(iel,ipproc(ivisls(itempk))) / cv0
          enddo
        endif
      else
        do iel = 1, ncel
          propce(iel,ipproc(ivisls(ienerg))) =         &
               visls0(itempk) / propce(iel,ipproc(icv))
        enddo
      endif
    endif

! ----- On donne la main a l'utilisateur

    call uscfxi                                                   &
    !==========
 ( idebia , idebra ,                                              &
   nvar   , nscal  , nphas  ,                                     &
   maxelt , lstelt ,                                              &
   ia     ,                                                       &
   dt     , rtp    , propce , propfa , propfb , coefa  , coefb  , &
   ra(iwcel1), ra(iwcel2), ra(iwcel3), ra(iwcel4),                &
   ra     )

! ----- Initialisation des proprietes physiques ROM et ROMB

    iirom  = ipproc(irom  )
    iiromb = ipprob(irom  )

    do iel = 1, ncel
      propce(iel,iirom)  = rtp(iel,isca(irho))
    enddo

    do ifac = 1, nfabor
      iel = ifabor(ifac)
      propfb(ifac,iiromb) =                                     &
           coefa(ifac,iclrtp(isca(irho),icoef))           &
           + coefb(ifac,iclrtp(isca(irho),icoef))           &
           * rtp(iel,isca(irho))
    enddo

! ----- Initialisation de la viscosite en volume

    if(iviscv.gt.0) then
      do iel = 1, ncel
        propce(iel,ipproc(iviscv)) = viscv0
      enddo
    endif

  endif

else

  if ( ipass.eq.1 ) then

! ----- Initialisations par defaut

    !     On initialise Cv

    iccfth = 0
    imodif = 1

    call uscfth                                                 &
    !==========
 ( nvar   , nscal  , nphas  ,                                     &
   iccfth , imodif , iphas  ,                                     &
   dt     , rtp    , rtp    , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   ra(iwcel1), ra(iwcel2), ra(iwcel3), ra(iwcel4) )

  endif

endif

!----
! FORMATS
!----


!----
! FIN
!----

return
end subroutine
