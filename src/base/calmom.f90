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

subroutine calmom &
!================

 ( idbia0 , idbra0 , ncel   , ncelet ,                            &
   nideve , nrdeve , nituse , nrtuse ,                            &
   idevel , ituser , ia     ,                                     &
   rtp    , dt     , propce ,                                     &
   rdevel , rtuser , ra )

!===============================================================================
! FONCTION :
! ----------

! CALCUL DE LA MOYENNE TEMPORELLE DE CORRELATIONS

!-------------------------------------------------------------------------------
!ARGU                             ARGUMENTS
!__________________.____._____.________________________________________________.
!    nom           !type!mode !                   role                         !
!__________________!____!_____!________________________________________________!
! idbia0           ! e  ! <-- ! numero de la 1ere case libre dans ia           !
! idbra0           ! e  ! <-- ! numero de la 1ere case libre dans ra           !
! ncelet           ! e  ! <-- ! nombre d'elements halo compris                 !
! ncel             ! e  ! <-- ! nombre de cellules                             !
! nideve nrdeve    ! e  ! <-- ! longueur de idevel rdevel                      !
! nituse nrtuse    ! e  ! <-- ! longueur de ituser rtuser                      !
! idevel(nideve    ! te ! <-- ! tab entier complementaire developemt           !
! ituser(nituse    ! te ! <-- ! tab entier complementaire utilisateur          !
! ia(*)            ! tr ! --- ! macro tableau entier                           !
! rtp              ! tr ! <-- ! tableaux des variables au pdt courant          !
! (ncelet,*   )    !    !     !                                                !
! dt(ncelet)       ! tr ! <-- ! pas de temps                                   !
! propce           ! tr ! <-- ! proprietes physiques au centre des             !
! (ncelet,*)       !    !     !    cellules                                    !
! rdevel(nrdeve    ! tr ! <-- ! tab reel complementaire developemt             !
! rtuser(nrtuse    ! tr ! <-- ! tab reel complementaire utilisateur            !
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
include "optcal.h"
include "numvar.h"

!===============================================================================

! Arguments

integer          idbia0 , idbra0 , ncel   , ncelet
integer          nideve , nrdeve , nituse , nrtuse
integer          idevel(nideve), ituser(nituse)
integer          ia(*)
double precision rtp(ncelet,*) , dt(ncelet), propce(ncelet,*)
double precision rdevel(nrdeve), rtuser(nrtuse), ra(*)

! VARIABLES LOCALES

integer          ii, iel, imom, icmom, idtcm, itravm, ifinra
integer          jmom, icumut

!===============================================================================

!     PHASE DE CUMUL DES TERMES
!     =========================

! Pour un tableau, on reserve la memoire ici.

itravm = idbra0
ifinra  = itravm + ncelet
CALL RASIZE('CALMOM',IFINRA)

! Boucle et test sur les moments a calculer
do imom = 1, nbmomt
  if(ntcabs.ge.ntdmom(imom)) then

!   Position dans PROPCE du tableau de cumul des moments
    icmom = ipproc(icmome(imom))

!   Tableau de travail pour recueillir la correlation par cellule
    do iel = 1, ncel
      ra(itravm+iel-1) = 1.d0
    enddo

!   Calcul de la correlation (on suppose qu'il y a au moins une variable)
    do ii = 1, idgmom(imom)

!     Si la variable est dans RTP
      if(idfmom(ii,imom).gt.0) then
        do iel = 1, ncel
          ra(itravm+iel-1) =                                      &
          ra(itravm+iel-1) * rtp(iel,idfmom(ii,imom))
        enddo
!     Si la variable est dans PROPCE
      elseif(idfmom(ii,imom).lt.0) then
        do iel = 1, ncel
          ra(itravm+iel-1) =                                      &
          ra(itravm+iel-1) * propce(iel,ipproc(-idfmom(ii,imom)))
        enddo
      endif

    enddo

!   Incrementation du tableau de cumul des correlations
    do iel = 1, ncel
      propce(iel,icmom) =                                         &
      propce(iel,icmom) + dt(iel)*ra(itravm+iel-1)
    enddo

!   Incrementation du tableau de cumul du temps (ncel ou constante)

!     Si plusieurs moyennes partagent le meme cumul temporel,
!     il ne faut incrementer que pour une seule (la premiere)

!     On decide si on doit incrementer
    icumut = 1
    if (imom.gt.1) then
      do jmom = 1, imom-1
        if(idtmom(jmom).eq.idtmom(imom)) then
          icumut = 0
        endif
      enddo
    endif

!     On incremente
    if(icumut.eq.1) then
      if(idtmom(imom).gt.0) then
        idtcm = ipproc(icdtmo(idtmom(imom)))
        do iel = 1, ncel
          propce(iel,idtcm) = propce(iel,idtcm) + dt(iel)
        enddo
      elseif(idtmom(imom).lt.0) then
        idtcm = -idtmom(imom)
        dtcmom(idtcm) = dtcmom(idtcm) + dt(1)
      endif
    endif

  endif
enddo

return

end subroutine
