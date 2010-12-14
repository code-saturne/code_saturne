!-------------------------------------------------------------------------------

!VERS


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

subroutine usativ &
!================

 ( idbia0 , idbra0 ,                                              &
   nvar   , nscal  , nphas  ,                                     &
   nbmetd , nbmett , nbmetm ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &
   maxelt , lstelt ,                                              &
   idevel , ituser , ia     ,                                     &
   dt     , rtp    , propce , propfa , propfb , coefa  , coefb  , &
   tmprom , ztprom , zdprom , xmet   , ymet   , pmer   ,          &
   ttprom , qvprom , uprom  , vprom  , ekprom , epprom ,          &
   rprom  , tpprom , phprom ,                                     &
   rdevel , rtuser , ra     )

!===============================================================================
! FONCTION :
! --------

! ROUTINE UTILISATEUR : INITIALISATION DES VARIABLES DE CALCUL
!     POUR LA PHYSIQUE PARTICULIERE: VERSION ATMOSPHERIQUE
!     PENDANT DE USINIV

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
!      PROPCE(IEL,IPPROC(IROM  (IPHAS))) designe ROM   (IEL ,IPHAS)
!      PROPCE(IEL,IPPROC(IVISCL(IPHAS))) designe VISCL (IEL ,IPHAS)
!      PROPCE(IEL,IPPROC(ICP   (IPHAS))) designe CP    (IEL ,IPHAS)
!      PROPCE(IEL,IPPROC(IVISLS(ISCAL))) designe VISLS (IEL ,ISCAL)

!      PROPFA(IFAC,IPPROF(IFLUMA(IVAR ))) designe FLUMAS(IFAC,IVAR)

!      PROPFB(IFAC,IPPROB(IROM  (IPHAS))) designe ROMB  (IFAC,IPHAS)
!      PROPFB(IFAC,IPPROB(IFLUMA(IVAR ))) designe FLUMAB(IFAC,IVAR)





! LA MODIFICATION DES PROPRIETES PHYSIQUES (ROM, VISCL, VISCLS, CP)
!     SE FERA EN STANDARD DANS LE SOUS PROGRAMME USPHYV
!     ET PAS ICI


! Cells identification
! ====================

! Cells may be identified using the 'getcel' subroutine.
! The syntax of this subroutine is described in the 'usclim' subroutine,
! but a more thorough description can be found in the user guide.


! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! nphas            ! i  ! <-- ! number of phases                               !
! nideve, nrdeve   ! i  ! <-- ! sizes of idevel and rdevel arrays              !
! nituse, nrtuse   ! i  ! <-- ! sizes of ituser and rtuser arrays              !
! maxelt           ! i  ! <-- ! max number of cells and faces (int/boundary)   !
! lstelt(maxelt)   ! ia ! --- ! work array                                     !
! idevel(nideve)   ! ia ! <-> ! integer work array for temporary development   !
! ituser(nituse)   ! ia ! <-> ! user-reserved integer work array               !
! ia(*)            ! ia ! --- ! main integer work array                        !
! dt(ncelet)       ! tr ! <-- ! valeur du pas de temps                         !
! rtp              ! tr ! <-- ! variables de calcul au centre des              !
! (ncelet,*)       !    !     !    cellules                                    !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! coefa coefb      ! tr ! <-- ! conditions aux limites aux                     !
!  (nfabor,*)      !    !     !    faces de bord                               !
! rdevel(nrdeve)   ! ra ! <-> ! real work array for temporary development      !
! rtuser(nrtuse)   ! ra ! <-> ! user-reserved real work array                  !
! ra(*)            ! ra ! --- ! main real work array                           !
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use pointe
use numvar
use optcal
use cstphy
use entsor
use parall
use period
use mesh

!===============================================================================

implicit none

integer          idbia0 , idbra0
integer          nvar   , nscal  , nphas
integer          nbmetd , nbmett , nbmetm
integer          nideve , nrdeve , nituse , nrtuse

integer          maxelt, lstelt(maxelt)
integer          idevel(nideve), ituser(nituse), ia(*)

double precision dt(ncelet), rtp(ncelet,*), propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision tmprom(nbmetm)
double precision ztprom(nbmett) , zdprom(nbmetd)
double precision xmet(nbmetm)   , ymet(nbmetm)  , pmer(nbmetm)
double precision ttprom(nbmett,nbmetm) , qvprom(nbmett,nbmetm)
double precision uprom(nbmetd,nbmetm)  , vprom(nbmetd,nbmetm)
double precision ekprom(nbmetd,nbmetm) , epprom(nbmetd,nbmetm)
double precision rprom(nbmett,nbmetm)  , tpprom(nbmett,nbmetm)
double precision phprom(nbmett,nbmetm)
double precision rdevel(nrdeve), rtuser(nrtuse), ra(*)

! Local variables

integer          idebia, idebra
integer          iel, iutile,iphas
double precision d2s3
double precision zent,xuent,xvent,xkent,xeent,tpent

!===============================================================================

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================

if(1.eq.1) then
!       Indicateur de non passage dans le sous-programme
  iusini = 0
  return
endif

!===============================================================================
! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END


!===============================================================================
! 1.  INITIALISATION VARIABLES LOCALES
!===============================================================================

idebia = idbia0
idebra = idbra0
d2s3 = 2.d0/3.d0

!===============================================================================
! 2. INITIALISATION DES INCONNUES :
!      UNIQUEMENT SI ON NE FAIT PAS UNE SUITE
!===============================================================================

! --- INCONNUES
!     Exemple : Initialisation a partir des profils meteo


if (isuite.eq.0) then

! Initialisation with the input meteo profile

  iphas=1

  do iel = 1, ncel

    zent=xyzcen(3,iel)

    call intprf                                                   &
    !==========
   (nbmetd, nbmetm,                                               &
    zdprom, tmprom, uprom , zent  , ttcabs, xuent )

    call intprf                                                   &
    !==========
   (nbmetd, nbmetm,                                               &
    zdprom, tmprom, vprom , zent  , ttcabs, xvent )

    call intprf                                                   &
    !==========
   (nbmetd, nbmetm,                                               &
    zdprom, tmprom, ekprom, zent  , ttcabs, xkent )

    call intprf                                                   &
    !==========
   (nbmetd, nbmetm,                                               &
    zdprom, tmprom, epprom, zent  , ttcabs, xeent )

    rtp(iel,iu(iphas))=xuent
    rtp(iel,iv(iphas))=xvent
    rtp(iel,iw(iphas))=0.d0

!     ITYTUR est un indicateur qui vaut ITURB/10
    if    (itytur(iphas).eq.2) then

      rtp(iel,ik(iphas))  = xkent
      rtp(iel,iep(iphas)) = xeent

    elseif(itytur(iphas).eq.3) then

      rtp(iel,ir11(iphas)) = d2s3*xkent
      rtp(iel,ir22(iphas)) = d2s3*xkent
      rtp(iel,ir33(iphas)) = d2s3*xkent
      rtp(iel,ir12(iphas)) = 0.d0
      rtp(iel,ir13(iphas)) = 0.d0
      rtp(iel,ir23(iphas)) = 0.d0
      rtp(iel,iep(iphas))  = xeent

    elseif(iturb(iphas).eq.50) then

      rtp(iel,ik(iphas))   = xkent
      rtp(iel,iep(iphas))  = xeent
      rtp(iel,iphi(iphas)) = d2s3
      rtp(iel,ifb(iphas))  = 0.d0

    elseif(iturb(iphas).eq.60) then

      rtp(iel,ik(iphas))   = xkent
      rtp(iel,iomg(iphas)) = xeent/cmu/xkent

    endif

    if (iscalt(iphas).ge.0) then
! On suppose que le scalaire est la temperature potentielle :
      call intprf                                                 &
      !==========
   (nbmett, nbmetm,                                               &
    ztprom, tmprom, tpprom, zent  , ttcabs, tpent )

      rtp(iel,isca(iscalt(iphas))) = tpent

    endif
  enddo

endif

!----
! FORMATS
!----

!----
! FIN
!----

return
end subroutine
