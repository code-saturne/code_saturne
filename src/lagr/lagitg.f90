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

subroutine lagitg &
!================

 ( nbpmax , nvp    , nvp1   ,                                     &
   ivar   ,                                                       &
   isorti , ibord  ,                                              &
   ettp   , ettpa  , tcarac , pip    , tsvar  )

!===============================================================================
! FONCTION :
! ----------

!   SOUS-PROGRAMME DU MODULE LAGRANGIEN :
!   -------------------------------------

!             INTEGRATION DE L'EDS POUR LA VARIABLE IVAR


!                 d V     V - PIP
!         EDS :   --- = - -------
!                 d t      TCARAC

!               Lorsqu'il y a eu interaction avec une face de bord,
!                 l'integration est degeneree a l'ordre 1
!                 (meme si on est a l'ordre 2 par ailleurs).

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nbpmax           ! e  ! <-- ! nombre max de particulies autorise             !
! nvp              ! e  ! <-- ! nombre de variables particulaires              !
! nvp1             ! e  ! <-- ! nvp sans position, vfluide, vpart              !
! ivar             ! e  ! <-- ! numero de la variable a integrer               !
!                  !    !     ! dans le tableau ettp                           !
! isorti(nbpmax    ! te ! <-- ! pour chaque particule :                        !
!                  !    !     !    * numero de sa cellule                      !
!                  !    !     !    * 0 si sortie du domaine                    !
! ibord            ! te ! <-- ! contient le numero de la                       !
!   (nbpmax)       !    !     !   face d'interaction part/frontiere            !
! ettp             ! tr ! --> ! tableaux des variables liees                   !
!  (nbpmax,nvp)    !    !     !   aux particules etape courante                !
! ettpa            ! tr ! <-- ! tableaux des variables liees                   !
!  (nbpmax,nvp)    !    !     !   aux particules etape precedente              !
! tcarac(nbpmax    ! tr ! <-- ! temps caracteristique associe a l'eds          !
! pip(nbpmax)      ! tr ! <-- ! second membre associe a l'eds                  !
! tsvar            ! tr ! --> ! prediction 1er sous-pas pour la                !
! (nbpmax,nvp1)    !    !     !   variable ivar, utilise pour la               !
!                  !    !     !   correction au 2eme sous-pas                  !
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
use cstphy
use cstnum
use optcal
use entsor
use lagpar
use lagran

!===============================================================================

implicit none

! Arguments

integer          nbpmax , nvp , nvp1
integer          ivar
integer          isorti(nbpmax) , ibord(nbpmax)

double precision ettp(nbpmax,nvp) ,  ettpa(nbpmax,nvp)
double precision tcarac(nbpmax) , pip(nbpmax)
double precision tsvar(nbpmax,nvp1)

! Local variables

integer          npt
double precision aux1 , aux2 , ter1 , ter2 , ter3

!===============================================================================


  if (nor.eq.1) then

    do npt = 1, nbpart
      if (isorti(npt).gt.0) then

        if (tcarac(npt).le.0.d0) then
          write(nfecra,2000) ivar, tcarac(npt), npt
          call csexit (1)
        endif

        aux1 = dtp/tcarac(npt)
        aux2 = exp(-aux1)

        ter1 = ettpa(npt,ivar) * aux2
        ter2 = pip(npt) * (1.d0-aux2)

!        Pour le cas NORDRE= 1 ou s'il y a rebond,
!          le ETTP suivant est le resultat final

        ettp(npt,ivar) = ter1 + ter2

!        Pour le cas NORDRE= 2, on calcule en plus TSVAR pour NOR= 2

        ter3 = ( -aux2 + (1.d0-aux2) / aux1 ) * pip(npt)
        tsvar(npt,ivar) = 0.5d0 * ter1 + ter3
      endif
    enddo

  else if (nor.eq.2) then

    do npt = 1, nbpart
      if (isorti(npt).gt.0 .and. ibord(npt).eq.0)  then

        if (tcarac(npt).le.0.d0) then
          write(nfecra,2000) ivar, tcarac(npt), npt
          call csexit (1)
        endif

        aux1 = dtp/tcarac(npt)
        aux2 = exp(-aux1)

        ter1 = 0.5d0 * ettpa(npt,ivar) * aux2
        ter2 = pip(npt) * ( 1.d0 - (1.d0-aux2) / aux1 )

!        Pour le cas NORDRE= 2, le ETTP suivant est le resultat final

        ettp(npt,ivar) = tsvar(npt,ivar) + ter1 + ter2
      endif
    enddo

  else
    write(nfecra,1000) nor
    call csexit (1)
    !==========
  endif

!===============================================================================

!-------
! FORMAT
!-------

 1000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========                                               ',/,&
'@                                                            ',/,&
'@    L''INDICATEUR SUR L''ORDRE D''INTEGRATION               ',/,&
'@       DES EQUATIONS DIFFERENTIELLES STOCHASTIQUES          ',/,&
'@       A UNE VALEUR NON PERMISE (LAGITG).                   ',/,&
'@                                                            ',/,&
'@    NORDRE DEVRAIT ETRE UN ENTIER EGAL A 1 OU 2             ',/,&
'@       IL VAUT ICI NORDRE = ', I10                           ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier la valeur de NORDRE dans la subroutine USLAG1.   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 2000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========                                               ',/,&
'@                                                            ',/,&
'@    LE TEMPS CARACTERISTIQUE LIE A L''EQUATION              ',/,&
'@      DIFFERENTIELLE STOCHASTIQUE DE LA VARIABLE            ',/,&
'@      NUMERO ',I10 ,'UNE VALEUR NON PERMISE (LAGITG).       ',/,&
'@                                                            ',/,&
'@    TCARAC DEVRAIT ETRE UN ENTIER STRICTEMENT POSITIF       ',/,&
'@       IL VAUT ICI TCARAC = ', E10.4                         ,/,&
'@       POUR LA PARTICULE NUMERO ',I10                        ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

!----
! FIN
!----

end subroutine
