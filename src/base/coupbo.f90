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

subroutine coupbo &
!================

 ( idbia0 , idbra0 ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  , isvtb  ,                            &
   nideve , nrdeve , nituse , nrtuse , ncp , ncv , ientha ,       &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   cpcst  , cp     , cvcst  , cv     ,                            &
   hbord  , tbord  ,                                              &
   rdevel , rtuser , ra     )

!===============================================================================
! FONCTION :
! ---------

! ECRITURE DE DONNEES RELATIVES A UN COUPLAGE AVEC SYRTHES

!-------------------------------------------------------------------------------
!ARGU                             ARGUMENTS
!__________________.____._____.________________________________________________.
!    nom           !type!mode !                   role                         !
!__________________!____!_____!________________________________________________!
! idbia0           ! e  ! <-- ! numero de la 1ere case libre dans ia           !
! idbra0           ! e  ! <-- ! numero de la 1ere case libre dans ra           !
! nfabor           ! e  ! <-- ! nombre de faces de bord                        !
! ncelet           ! e  ! <-- ! nombre d'elements halo compris                 !
! ncp              ! e  ! <-- ! dimension de cp (ncelet ou 1)                  !
! ncv              ! e  ! <-- ! dimension de cv (ncelet ou 1)                  !
! nfabor           ! e  ! <-- ! nombre de faces de bord                        !
! nvar             ! e  ! <-- ! nombre total de variables                      !
! nscal            ! e  ! <-- ! nombre total de scalaires                      !
! nphas            ! e  ! <-- ! nombre de phases                               !
! nideve nrdeve    ! e  ! <-- ! longueur de idevel rdevel                      !
! nituse nrtuse    ! e  ! <-- ! longueur de ituser rtuser                      !
! ientha           ! e  ! <-- ! 1 si tparoi est une enthalpie                  !
!                  ! e  ! <-- ! 2 si tparoi est une energie                    !
!                  !    !     !    (compressible)                              !
! ifabor           ! te ! <-- ! element  voisin  d'une face de bord            !
! (nfabor)         !    !     !                                                !
! idevel(nideve    ! te ! <-- ! tab entier complementaire developemt           !
! ituser(nituse    ! te ! <-- ! tab entier complementaire utilisateur          !
! ia(*)            ! tr ! --- ! macro tableau entier                           !
! cpcst            ! r  ! <-- ! chaleur specifique si constante                !
! cvcst            ! r  ! <-- ! chaleur specifique si constante                !
! cp(ncp)          ! tr ! <-- ! chaleur specifique si variable                 !
! cv(ncp)          ! tr ! <-- ! chaleur specifique si variable                 !
! hbord            ! tr ! <-- ! coefficients d'echange aux bords               !
! (nfabor)         !    !     !                                                !
! tbord            ! tr ! <-- ! temperatures aux bords                         !
! (nfabor)         !    !     !                                                !
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
include "entsor.h"
include "cstphy.h"

!===============================================================================

! Arguments

integer          idbia0, idbra0
integer          ndim   , ncelet , ncel   , nfac   , nfabor
integer          nfml   , nprfml
integer          nnod   , lndfac , lndfbr , ncelbr
integer          nvar   , nscal  , nphas
integer          isvtb
integer          nideve , nrdeve , nituse , nrtuse
integer          ncp    , ncv    , ientha

integer          ifacel(2,nfac) , ifabor(nfabor)
integer          ifmfbr(nfabor) , ifmcel(ncelet)
integer          iprfml(nfml,nprfml)
integer          ipnfac(nfac+1), nodfac(lndfac)
integer          ipnfbr(nfabor+1), nodfbr(lndfbr)
integer          idevel(nideve), ituser(nituse), ia(*)

double precision cpcst  , cvcst

double precision xyzcen(ndim,ncelet)
double precision surfac(ndim,nfac), surfbo(ndim,nfabor)
double precision cdgfac(ndim,nfac), cdgfbo(ndim,nfabor)
double precision xyznod(ndim,nnod), volume(ncelet)
double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*),propfa(nfac,*),propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision cp(ncp), cv(ncv)
double precision hbord(nfabor),tbord(nfabor)
double precision rdevel(nrdeve), rtuser(nrtuse), ra(*)

! VARIABLES LOCALES

integer          nbccou, inbcou, inbcoo, nbfcou, ifac, iloc, iel
integer          itflui, ihparo
integer          idebia, idebra, ifinia, ifinra, ipfcou, mode
integer          iccfth, imodif, iphas
integer          iepsel, iepsfa, igamag, ixmasm
double precision enthal, temper, energ, cvt

!===============================================================================

idebia = idbia0
idebra = idbra0

!===============================================================================
!     COUPLAGE SYRTHES : CALCUL DE LA TEMPERATURE FLUIDE ET DU
!                        COEFFICIENT D'ECHANGE
!===============================================================================

!     RECUPERATION DU NOMBRE DE CAS DE COUPLAGE

call nbcsyr (nbccou)
!==========

!---> BOUCLE SUR LES CAS DE COUPLAGE DE TYPE SYRTHES

do inbcou = 1, nbccou

!        NOMBRE DE FACES DE BORD PAR CAS DE COUPLAGE
  inbcoo = inbcou
  call nbfsyr (inbcoo, nbfcou)
  !==========

!        GESTION MEMOIRE POUR CONSTRUIRE LES TABLEAUX

  ipfcou = idebia
  ifinia = ipfcou + nbfcou

  itflui = idebra
  ihparo = itflui + nbfcou
  ifinra = ihparo + nbfcou

! Compressible : couplage avec l'energie

  if(ientha .eq. 2) then
    iepsel = ifinra
    iepsfa = iepsel + ncelet
    igamag = iepsfa + nfabor
    ixmasm = igamag + ncelet
    ifinra = ixmasm + ncelet
  endif

! Fin Compressible

  CALL RASIZE('COUPBO',IFINIA)
  !==========
  CALL RASIZE('COUPBO',IFINRA)
  !==========

!       BOUCLE SUR LES FACES DE COUPLAGE ET CALCUL DES COEFFICIENTS

  inbcoo = inbcou
  call lfasyr(inbcoo, ia(ipfcou))
  !==========

  do iloc = 1, nbfcou

    ifac = ia(ipfcou+iloc-1)

!           TEMPERATURES FLUIDES SAUVEGARDEES
    ra(itflui+iloc-1) = tbord(ifac)

!           COEFFICIENTS D'ECHANGE SAUVEGARDES
    ra(ihparo+iloc-1) = hbord(ifac)

  enddo

!        SI ENTHALPIE, ON TRANSFORME EN TEMPERATURE
!          Il est necessaire de transmettre a SYRTHES des Temperatures
!          Afin de conserver le flux Phi = (lambda/d     ) Delta T
!                                 ou Phi = (lambda/(d Cp)) Delta H
!            on multiplie HBORD = lambda/(d Cp) par Cp pris dans la
!              cellule adjacente.
!          Le resultat n'est pas garanti (conservation en particulier),
!             on ajoute donc un avertissement.


  if(ientha.eq.1) then

    write(nfecra,1000)
    mode = 1
    do iloc = 1, nbfcou
      ifac = ia(ipfcou+iloc-1)
      iel  = ifabor(ifac)
      enthal = ra(itflui+iloc-1)
      call usthht (mode   , enthal , temper  )
      !==========
      ra(itflui+iloc-1) = temper
      if(ncp.eq.ncelet) then
        ra(ihparo+iloc-1) = ra(ihparo+iloc-1)*cp(iel)
      else
        ra(ihparo+iloc-1) = ra(ihparo+iloc-1)*cpcst
      endif
    enddo

  else if(ientha.eq.2) then

!        SI ENERGIE, ON TRANSFORME EN TEMPERATURE
!          Il est necessaire de transmettre a SYRTHES des Temperatures
!          Afin de conserver le flux Phi = (lambda/d     ) Delta T
!                                 ou Phi = (lambda/(d Cv)) (Ei - Ep)
!            on multiplie HBORD = lambda/(d Cv) par Cv pris dans la
!              cellule adjacente.
!            noter que Ei = Cv Ti + 1/2 Ui*Ui + Epsilon_sup_i
!               et que Ep = Cv Tp + 1/2 Ui*Ui + Epsilon_sup_i
!               (l'ecart est donc bien Cv Delta T)

!     Modif temperature et coeff d'echange

!       Calcul de e - CvT

    iccfth = 7
    imodif = 0
    iphas  = iphsca(isvtb)

    call uscfth                                                   &
    !==========
 ( ifinia , ifinra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   iccfth , imodif , iphas  ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   ra(iepsel) , ra(iepsfa) , ra(igamag) , ra(ixmasm) ,            &
!        ----------   ---------
   rdevel , rtuser , ra     )

    do iloc = 1, nbfcou
      ifac  = ia(ipfcou+iloc-1)
      iel   = ifabor(ifac)
      energ = ra(itflui+iloc-1)
      cvt   = energ                                               &
             -(0.5d0*( rtp(iel,iu(iphas))**2                      &
                      +rtp(iel,iv(iphas))**2                      &
                      +rtp(iel,iw(iphas))**2)                     &
               + ra(iepsel+iel-1)           )
       if(ncv.eq.ncelet) then
         ra(itflui+iloc-1) = cvt/cv(iel)
         ra(ihparo+iloc-1) = ra(ihparo+iloc-1)*cv(iel)
       else
         ra(itflui+iloc-1) = cvt/cvcst
         ra(ihparo+iloc-1) = ra(ihparo+iloc-1)*cvcst
       endif
     enddo

  endif

!       ENVOI DE LA TEMPERATURE FLUIDE ET DU COEFFICIENT D'ECHANGE

  inbcoo = inbcou
  call varsyo (inbcoo, ra(itflui), ra(ihparo))
  !==========

enddo

!===============================================================================
!     FIN DES COUPLAGES DE BORD
!===============================================================================

return

! FORMATS

#if defined(_CS_LANG_FR)

 1000 format(                                                           &
'@                                                            ',/,&
'@ @@ ATTENTION : COUPLAGE SYRTHES AVEC CALCUL EN ENTHALPIE   ',/,&
'@    =========                                               ',/,&
'@      OPTION NON VALIDEE - CONTACTER L''EQUIPE DE DVPT      ',/,&
'@                                                            ')

#else

 1000 format(                                                           &
'@                                                            ',/,&
'@ @@ WARNING: SYRTHES COUPLING WITH ENTHALPY CALCULATION     ',/,&
'@    ========                                                ',/,&
'@      OPTION NOT VALIDATED - CONTACT THE SUPPORT            ',/,&
'@                                                            ')

#endif

end
