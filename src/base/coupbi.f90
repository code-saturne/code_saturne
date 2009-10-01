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

subroutine coupbi &
!================

 ( idbia0 , idbra0 ,                                              &
   nfabor ,                                                       &
   nvar   , nscal  , nphas  ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &
   icodcl ,                                                       &
   idevel , ituser , ia     ,                                     &
   rcodcl ,                                                       &
   rdevel , rtuser , ra     )

!===============================================================================
! FONCTION :
! ---------

! LECTURE DE DONNEES RELATIVES A UN COUPLAGE AVEC SYRTHES

!-------------------------------------------------------------------------------
!ARGU                             ARGUMENTS
!__________________.____._____.________________________________________________.
!    nom           !type!mode !                   role                         !
!__________________!____!_____!________________________________________________!
! idbia0           ! e  ! <-- ! numero de la 1ere case libre dans ia           !
! idbra0           ! e  ! <-- ! numero de la 1ere case libre dans ra           !
! nfabor           ! e  ! <-- ! nombre de faces de bord                        !
! nvar             ! e  ! <-- ! nombre total de variables                      !
! nscal            ! e  ! <-- ! nombre total de scalaires                      !
! nphas            ! e  ! <-- ! nombre de phases                               !
! nideve nrdeve    ! e  ! <-- ! longueur de idevel rdevel                      !
! nituse nrtuse    ! e  ! <-- ! longueur de ituser rtuser                      !
! icodcl           ! te ! --> ! code de condition limites aux faces            !
!  (nfabor,nvar    !    !     !  de bord                                       !
!                  !    !     ! = 1   -> dirichlet                             !
!                  !    !     ! = 3   -> densite de flux                       !
!                  !    !     ! = 4   -> glissemt et u.n=0 (vitesse)           !
!                  !    !     ! = 5   -> frottemt et u.n=0 (vitesse)           !
!                  !    !     ! = 6   -> rugosite et u.n=0 (vitesse)           !
!                  !    !     ! = 9   -> entree/sortie libre (vitesse          !
!                  !    !     !  entrante eventuelle     bloquee               !
! idevel(nideve    ! te ! <-- ! tab entier complementaire developemt           !
! ituser(nituse    ! te ! <-- ! tab entier complementaire utilisateur          !
! ia(*)            ! tr ! --- ! macro tableau entier                           !
! rcodcl           ! tr ! --> ! valeur des conditions aux limites              !
!  (nfabor,nvar    !    !     !  aux faces de bord                             !
!                  !    !     ! rcodcl(1) = valeur du dirichlet                !
!                  !    !     ! rcodcl(2) = valeur du coef. d'echange          !
!                  !    !     !  ext. (infinie si pas d'echange)               !
!                  !    !     ! rcodcl(3) = valeur de la densite de            !
!                  !    !     !  flux (negatif si gain) w/m2 ou                !
!                  !    !     !  hauteur de rugosite (m) si icodcl=6           !
!                  !    !     ! pour les vitesses (vistl+visct)*gradu          !
!                  !    !     ! pour la pression             dt*gradp          !
!                  !    !     ! pour les scalaires                             !
!                  !    !     !        cp*(viscls+visct/sigmas)*gradt          !
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
include "cstnum.h"
include "cstphy.h"
include "entsor.h"
include "ppppar.h"
include "ppthch.h"
include "ppincl.h"

!===============================================================================

! Arguments

integer          idbia0, idbra0
integer          nfabor
integer          nvar , nscal , nphas
integer          nideve , nrdeve , nituse , nrtuse
integer          icodcl(nfabor,nvar)
integer          idevel(nideve), ituser(nituse), ia(*)
double precision rcodcl(nfabor,nvar,3)
double precision rdevel(nrdeve), rtuser(nrtuse), ra(*)

! VARIABLES LOCALES

integer          ll, nbccou, inbcou, inbcoo, nbfcou
integer          ifac, iloc, iscal , iphas
integer          idebia, idebra, ifinia, ifinra, ipfcou, ithpar
integer          icldef
integer          mode
double precision temper, enthal


!===============================================================================

idebia = idbia0
idebra = idbra0

!===============================================================================
!     COUPLAGE SYRTHES : RECUPERATION DE LA TEMPERATURE DE PAROI
!===============================================================================

!     RECUPERATION DU NOMBRE DE CAS DE COUPLAGE

call nbcsyr (nbccou)
!==========

!---> BOUCLE SUR LES COUPLAGES : RECEPTION DES TABLEAUX TPAROI POUR
!                                CHAQUE COUPLAGE ET APPLICATION DE
!                                LA CONDITION LIMITE CORRESPONDANTE

do inbcou = 1, nbccou

!       NOMBRE DE FACES DE BORD PAR CAS DE COUPLAGE
  inbcoo = inbcou
  call nbfsyr (inbcoo, nbfcou)
  !==========

!        GESTION MEMOIRE POUR RECEVOIR LE TABLEAU

  ipfcou = idebia
  ifinia = ipfcou + nbfcou

  ithpar = idebra
  ifinra = ithpar + nbfcou

  CALL IASIZE('COUPBI',IFINIA)
  !==========
  CALL RASIZE('COUPBI',IFINRA)
  !==========

!        LECTURE DU MESSAGE (TEMPERATURE PAROI) ET
!        INTERPOLATION DES VALEURS AUX SOMMETS AUX FACES
  call varsyi (inbcou, ra(ithpar))
  !==========

!        ON IMPOSE LA TEMPERATURE A LA PAROI

  inbcoo = inbcou
  call lfasyr(inbcoo, ia(ipfcou))
  !==========

!        TYPE DE CONDITION PAR DEFAUT
  icldef = 5


  do iscal = 1, nscal

    if(icpsyr(iscal).eq.1) then

! --- Pour les scalaires qui sont couples a SYRTHES
!       on impose une condition de Dirichlet aux faces couplees
!     Pour le moment, on ne passe ici qu'une seule fois,
!       etant entendu que un seul scalaire est couple a SYRTHES
!     Pour le module compressible, on resout en energie, mais on
!       stocke la temperature a part, pour que ce soit plus clair
!       dans les conditions aux limites


      ll = isca(iscal)
      if(ippmod(icompf).ge.0) then
        iphas = iphsca(iscal)
        if(iscal.eq.ienerg(iphas)) then
          ll = isca(itempk(iphas))
        else
          write(nfecra,9000)ienerg(iphas),iscal
          call csexit (1)
        endif
      endif


      do iloc = 1, nbfcou

        ifac = ia(ipfcou+iloc-1)

        if ((icodcl(ifac,ll) .ne. 1) .and.                        &
            (icodcl(ifac,ll) .ne. 5) .and.                        &
            (icodcl(ifac,ll) .ne. 6)) icodcl(ifac,ll) = icldef

        rcodcl(ifac,ll,1) = ra(ithpar+iloc-1)
        rcodcl(ifac,ll,2) = rinfin
        rcodcl(ifac,ll,3) = 0.d0

      enddo

      ! Conversion eventuelle temperature -> enthalpie

      if(iscsth(iscal).eq.2) then

        do iloc = 1, nbfcou

          ifac = ia(ipfcou+iloc-1)

          temper = rcodcl(ifac,ll,1)
          mode   = -1
          call usthht(mode,enthal,temper)
          !==========
          rcodcl(ifac,ll,1) = enthal

        enddo

      endif

    endif

  enddo

enddo

!===============================================================================
!     FIN DES COUPLAGES DE BORD
!===============================================================================

return

! FORMATS

#if defined(_CS_LANG_FR)

 9000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET LORS DU COUPLAGE SYRTHES              ',/,&
'@    =========                                               ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Avec le module compressible, seul le scalaire ',I10        ,/,&
'@    peut etre couple a SYRTHES. Ici, on cherche a coupler   ',/,&
'@    le scalaire ',I10                                        ,/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

#else

 9000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: ABORT IN THE SYRTHES COUPLING                  ',/,&
'@    ========                                                ',/,&
'@                                                            ',/,&
'@  The calculation will not be run.                          ',/,&
'@                                                            ',/,&
'@  With the compressible module, only the scalar ',I10        ,/,&
'@    may be coupled with SYRTHES. Here, one tries to couple  ',/,&
'@    with the scalar ',I10                                    ,/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

#endif

end subroutine
