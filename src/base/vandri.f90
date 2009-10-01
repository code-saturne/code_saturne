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

subroutine vandri &
!================

 (  ndim   , ncelet , ncel   , nfac   , nfabor , nphas ,          &
    nideve , nrdeve , nituse , nrtuse , iphas  ,                  &
    itypfb , ifabor , ifapat , idevel , ituser , ia    ,          &
    xyzcen , cdgfbo , uetbor , visvdr , yplusc , propce ,         &
    rdevel , rtuser , ra    )

!===============================================================================
! FONCTION :
! ----------

! IMPOSITION D'UN AMORTISSEMENT DE TYPE VAN DRIEST POUR LA LES
! nut est amortie par (1-exp(-y+/d+))**2 ou d+ est mis par defaut a 26

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
!    nom           !type!mode !                   role                         !
!__________________!____!_____!________________________________________________!
! ndim             ! e  ! <-- ! dimension de l'espace                          !
! ncelet           ! e  ! <-- ! nombre d'elements halo compris                 !
! ncel             ! e  ! <-- ! nombre d'elements actifs                       !
! nfac             ! e  ! <-- ! nombre de faces internes                       !
! nfabor           ! e  ! <-- ! nombre de faces de bord                        !
! nphas            ! e  ! <-- ! nombre de phases                               !
! nideve nrdeve    ! e  ! <-- ! longueur de idevel rdevel                      !
! nituse nrtuse    ! e  ! <-- ! longueur de ituser rtuser                      !
! iphas            ! e  ! <-- ! numero de la phase traitee                     !
! itypfb(nfabor    ! te ! <-- ! type des faces de bord                         !
!  nphas)          !    !     !                                                !
! ifabor           ! te ! <-- ! element  voisin  d'une face de bord            !
! (nfabor)         !    !     !                                                !
! ifapat           ! te ! <-- ! no de face de brd code 5 la + proche           !
! (ncelet)         !    !     !    (rij et echo de paroi      )                !
! idevel(nideve    ! te ! <-- ! tab entier complementaire developemt           !
! ituser(nituse    ! te ! <-- ! tab entier complementaire utilisateur          !
! ia(*)            ! tr ! --- ! macro tableau entier                           !
! xyzcen           ! tr ! <-- ! point associes aux volumes de control          !
! (ndim,ncelet     !    !     !                                                !
! cdgfbo           ! tr ! <-- ! centre de gravite des faces de bord            !
! (ndim,nfabor)    !    !     !                                                !
! uetbor           ! tr ! <-- ! vitesse de frottement au bord                  !
! (nfabor,nphas    !    !     !  pour van driest en les                        !
! visvdr(nphas)    ! tr ! <-- ! viscosite dynamique ds les cellules            !
! (ncelet,nphas    !    !     !  de bord apres amortisst de v driest           !
! yplusc           ! tr ! <-- ! valeur de yplus aux cellules                   !
! (ncelet  )       !    !     !    dans le cas abs(icdpar).eq.1                !
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
include "numvar.h"
include "optcal.h"
include "entsor.h"
include "cstphy.h"
include "parall.h"

!===============================================================================

! Arguments

integer          ndim, ncelet , ncel   , nfac   , nfabor, nphas
integer          nideve , nrdeve , nituse , nrtuse
integer          iphas
integer          itypfb(nfabor,nphas),ifabor(nfabor)
integer          ifapat(ncelet)
integer          idevel(nideve), ituser(nituse)
integer          ia(*)

double precision xyzcen(ndim,ncelet),cdgfbo(ndim,nfabor)
double precision uetbor(nfabor,nphas), visvdr(ncelet,nphas)
double precision yplusc(ncelet)
double precision propce(ncelet,*)
double precision rdevel(nrdeve), rtuser(nrtuse), ra(*)

! VARIABLES LOCALES

integer          iel   , ifac  , ipcvis, ipcvst, ipcrom
double precision yplus , yminpa, viscos

!===============================================================================

ipcrom = ipproc(irom  (iphas))
ipcvis = ipproc(iviscl(iphas))
ipcvst = ipproc(ivisct(iphas))

!     Calcul direct de la distance a la paroi (non compatible parall/perio)
if(abs(icdpar).eq.2) then

!     En sequentiel, RAS
  if(irangp.lt.0) then
    do iel = 1, ncel
      ifac = ifapat(iel)
      viscos = propce(iel,ipcvis)/propce(iel,ipcrom)
      yminpa = sqrt((cdgfbo(1,ifac)-xyzcen(1,iel))**2             &
           +        (cdgfbo(2,ifac)-xyzcen(2,iel))**2             &
           +        (cdgfbo(3,ifac)-xyzcen(3,iel))**2)
      yplus = uetbor(ifac,iphas) * yminpa/ viscos
      propce(iel,ipcvst) = propce(iel,ipcvst)*                    &
           (1.0d0-exp(-yplus/cdries(iphas)))**2
    enddo
!     En parallele, on n'amortit que la premiere maille de paroi :
!     dangereux mais a priori inutile (car l'utilisation de
!     ICDPAR=+/-2 en parallele est bloque dans verini)
  else
    write(nfecra,1000)
    do ifac = 1, nfabor
      if(itypfb(ifac,iphas).eq.iparoi .or.                        &
         itypfb(ifac,iphas).eq.iparug ) then
        iel = ifabor(ifac)
        viscos = propce(iel,ipcvis)/propce(iel,ipcrom)
        yminpa = sqrt((cdgfbo(1,ifac)-xyzcen(1,iel))**2           &
             +        (cdgfbo(2,ifac)-xyzcen(2,iel))**2           &
             +        (cdgfbo(3,ifac)-xyzcen(3,iel))**2)
        yplus = uetbor(ifac,iphas) * yminpa/ viscos
        propce(iel,ipcvst) = propce(iel,ipcvst)*                  &
             (1.0d0-exp(-yplus/cdries(iphas)))**2
      endif
    enddo
endif

!     Nouveau mode de calcul : c'est plus simple
elseif(abs(icdpar).eq.1) then
  do iel = 1, ncel
    yplus = yplusc(iel)
    propce(iel,ipcvst) = propce(iel,ipcvst)*                      &
           (1.0d0-exp(-yplus/cdries(iphas)))**2
  enddo
endif

!     Pour les cellules de paroi on remet la viscosite turbulente
!     qui avait ete amortie dans clptur et qui a servi a calculer
!     les conditions aux limites
do iel = 1, ncel
  if (visvdr(iel,iphas).gt.-900.d0)                               &
       propce(iel,ipcvst) = visvdr(iel,iphas)
enddo

!--------
! FORMATS
!--------
#if defined(_CS_LANG_FR)

 1000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : DANS LE CAS DE LA LES AVEC AMORTISSEMENT    ',/,&
'@    =========                                               ',/,&
'@    L''AMORTISSEMENT DE VAN DRIEST N''EST FAIT QUE SUR LA   ',/,&
'@    PREMIERE CELLULE A LA PAROI EN CAS DE PARALLELISME      ',/,&
'@                                                            ',/,&
'@  Le calcul se poursuit.                                    ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

#else

 1000 format(                                                           &
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/,&
'@ @@ WARNING: IN CASE OF LES WITH DAMPING'                    ,/,&
'@    ========'                                                ,/,&
'@    VAN DRIEST DAMPING IS ONLY EFFECTIVE ON THE FIRST CELL'  ,/,&
'@    OFF-WALL IN CASE OF PARALLELISM'                         ,/,&
'@'                                                            ,/,&
'@  The calculation will be run.'                              ,/,&
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/)

#endif
!----
! FIN
!----


return
end subroutine
