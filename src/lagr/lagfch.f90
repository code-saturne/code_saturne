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

subroutine lagfch &
!================

 ( idbia0 , idbra0 ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   nbpmax , nvp    , nvp1   , nvep   , nivep  ,                   &
   ntersl , nvlsta , nvisbr ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &
   itepa  , ibord  , idevel , ituser , ia     ,                   &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtpa   , rtp    , propce , propfa , propfb ,          &
   ettp   , ettpa  , tepa   , statis , stativ ,                   &
   taup   , tlag   , piil   ,                                     &
   tsuf   , tsup   , bx     , tsfext ,                            &
   vagaus , gradpr , gradvf ,                                     &
   romp   , fextla ,                                              &
   rdevel , rtuser , ra     )

!===============================================================================
! FONCTION :
! ----------

!   SOUS-PROGRAMME DU MODULE LAGRANGIEN :
!   -------------------------------------

!    CALCUL DES FORCES DLVO

!       - FORCES DE VAN DER WAALS
!       - FORCES ELECTROSTATIQUES

!    ELLES DOIVENT ETRE CONNUES EN CHAQUE CELLULE
!      ET ETRE HOMOGENES A LA GRAVITE (M/S2)


!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
! ndim             ! i  ! <-- ! spatial dimension                              !
! ncelet           ! i  ! <-- ! number of extended (real + ghost) cells        !
! ncel             ! i  ! <-- ! number of cells                                !
! nfac             ! i  ! <-- ! number of interior faces                       !
! nfabor           ! i  ! <-- ! number of boundary faces                       !
! nfml             ! i  ! <-- ! number of families (group classes)             !
! nprfml           ! i  ! <-- ! number of properties per family (group class)  !
! nnod             ! i  ! <-- ! number of vertices                             !
! lndfac           ! i  ! <-- ! size of nodfac indexed array                   !
! lndfbr           ! i  ! <-- ! size of nodfbr indexed array                   !
! ncelbr           ! i  ! <-- ! number of cells with faces on boundary         !
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! nphas            ! i  ! <-- ! number of phases                               !
! nbpmax           ! e  ! <-- ! nombre max de particulies autorise             !
! nvp              ! e  ! <-- ! nombre de variables particulaires              !
! nvp1             ! e  ! <-- ! nvp sans position, vfluide, vpart              !
! nvep             ! e  ! <-- ! nombre info particulaires (reels)              !
! nivep            ! e  ! <-- ! nombre info particulaires (entiers)            !
! ntersl           ! e  ! <-- ! nbr termes sources de couplage retour          !
! nvlsta           ! e  ! <-- ! nombre de var statistiques lagrangien          !
! nvisbr           ! e  ! <-- ! nombre de statistiques aux frontieres          !
! nideve, nrdeve   ! i  ! <-- ! sizes of idevel and rdevel arrays              !
! nituse, nrtuse   ! i  ! <-- ! sizes of ituser and rtuser arrays              !
! itepa            ! te ! <-- ! info particulaires (entiers)                   !
! (nbpmax,nivep    !    !     !   (cellule de la particule,...)                !
! ibord            ! te ! --> ! si nordre=2, contient le numero de la          !
!   (nbpmax)       !    !     !   face d'interaction part/frontiere            !
! idevel(nideve)   ! ia ! <-> ! integer work array for temporary development   !
! ituser(nituse)   ! ia ! <-> ! user-reserved integer work array               !
! ia(*)            ! ia ! --- ! main integer work array                        !
! xyzcen           ! ra ! <-- ! cell centers                                   !
!  (ndim, ncelet)  !    !     !                                                !
! surfac           ! ra ! <-- ! interior faces surface vectors                 !
!  (ndim, nfac)    !    !     !                                                !
! surfbo           ! ra ! <-- ! boundary faces surface vectors                 !
!  (ndim, nfabor)  !    !     !                                                !
! cdgfac           ! ra ! <-- ! interior faces centers of gravity              !
!  (ndim, nfac)    !    !     !                                                !
! cdgfbo           ! ra ! <-- ! boundary faces centers of gravity              !
!  (ndim, nfabor)  !    !     !                                                !
! xyznod           ! tr ! <-- ! coordonnes des noeuds                          !
! (ndim,nnod)      !    !     !                                                !
! volume(ncelet    ! tr ! <-- ! volume d'un des ncelet elements                !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! ettp             ! tr ! <-- ! tableaux des variables liees                   !
!  (nbpmax,nvp)    !    !     !   aux particules etape courante                !
! ettpa            ! tr ! <-- ! tableaux des variables liees                   !
!  (nbpmax,nvp)    !    !     !   aux particules etape precedente              !
! tepa             ! tr ! <-- ! info particulaires (reels)                     !
! (nbpmax,nvep)    !    !     !   (poids statistiques,...)                     !
! statis           ! tr ! <-- ! cumul pour les moyennes des                    !
!(ncelet,nvlsta    !    !     !    statistiques volumiques                     !
! stativ           ! tr ! <-- ! cumul pour les variances des                   !
!(ncelet,          !    !     !    statistiques volumiques                     !
!   nvlsta-1)      !    !     !                                                !
! taup(nbpmax)     ! tr ! <-- ! temps caracteristique dynamique                !
! tlag(nbpmax)     ! tr ! <-- ! temps caracteristique fluide                   !
! piil(nbpmax,3    ! tr ! <-- ! terme dans l'integration des eds up            !
! tsup(nbpmax,3    ! tr ! <-- ! prediction 1er sous-pas pour                   !
!                  !    !     !   la vitesse des particules                    !
! tsuf(nbpmax,3    ! tr ! <-- ! prediction 1er sous-pas pour                   !
!                  !    !     !   la vitesse du fluide vu                      !
! bx(nbpmax,3,2    ! tr ! <-- ! caracteristiques de la turbulence              !
! tsfext(nbpmax    ! tr ! <-- ! infos pour le couplage retour                  !
! vagaus           ! tr ! <-- ! variables aleatoires gaussiennes               !
!(nbpmax,nvgaus    !    !     !                                                !
! gradpr(ncel,3    ! tr ! <-- ! gradient de pression                           !
! gradvf(ncel,3    ! tr ! <-- ! gradient de la vitesse du fluide               !
! romp             ! tr ! --- ! masse volumique des particules                 !
! fextla           ! tr ! --> ! champ de forces exterieur                      !
!(ncelet,3)        !    !     !    utilisateur (m/s2)                          !
! rdevel(nrdeve)   ! ra ! <-> ! real work array for temporary development      !
! rtuser(nrtuse)   ! ra ! <-> ! user-reserved real work array                  !
! ra(*)            ! ra ! --- ! main real work array                           !
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
include "numvar.f90"
include "cstnum.f90"
include "cstphy.f90"
include "optcal.f90"
include "entsor.f90"
include "lagpar.f90"
include "lagran.f90"
include "ppppar.f90"
include "ppthch.f90"
include "ppincl.f90"

!===============================================================================

! Arguments

integer          idbia0 , idbra0
integer          ndim   , ncelet , ncel   , nfac   , nfabor
integer          nfml   , nprfml
integer          nnod   , lndfac , lndfbr , ncelbr
integer          nvar   , nscal  , nphas
integer          nbpmax , nvp    , nvp1   , nvep  , nivep
integer          ntersl , nvlsta , nvisbr
integer          nideve , nrdeve , nituse , nrtuse
integer          itepa(nbpmax,nivep) , ibord(nbpmax)
integer          idevel(nideve), ituser(nituse)
integer          ia(*)

double precision xyzcen(ndim,ncelet)
double precision surfac(ndim,nfac) , surfbo(ndim,nfabor)
double precision cdgfac(ndim,nfac) , cdgfbo(ndim,nfabor)
double precision xyznod(ndim,nnod) , volume(ncelet)
double precision dt(ncelet) , rtp(ncelet,*) , rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*) , propfb(nfabor,*)
double precision ettp(nbpmax,nvp) , ettpa(nbpmax,nvp)
double precision tepa(nbpmax,nvep)
double precision statis(ncelet,*),stativ(ncelet,*)
double precision taup(nbpmax) , tlag(nbpmax,3)
double precision piil(nbpmax,3) , bx(nbpmax,3,2)
double precision tsuf(nbpmax,3) , tsup(nbpmax,3)
double precision tsfext(nbpmax)
double precision vagaus(nbpmax,*)
double precision gradpr(ncelet,3) , gradvf(ncelet,9)
double precision romp(nbpmax)
double precision fextla(nbpmax,3)
double precision rdevel(nrdeve), rtuser(nrtuse)
double precision ra(*)

! Local variables

integer          idebia, idebra,ifinia , ifinra
integer          idppar, inxpar, inypar, inzpar
integer          ip , iel , iphas , mode

double precision val , tempf , dnorm
double precision debye, aa

!===============================================================================

!===============================================================================
! 0.  GESTION MEMOIRE ET INITIALISATION
!===============================================================================

idebia = idbia0
idebra = idbra0

iphas  = 1

!===============================================================================
! 1. CALCUL DE LA DISTANCE A LA PAROI + NORMAL A LA PAROI
!===============================================================================

ifinia = idebia
idppar = idebra
inxpar = idppar + nbpart
inypar = inxpar + nbpart
inzpar = inypar + nbpart
ifinra = inzpar + nbpart
CALL RASIZE('LAGFCH',IFINRA)
!==========

do ip = 1,nbpart
  ra(idppar+ip-1) = 0.d0
  ra(inxpar+ip-1) = 0.d0
  ra(inypar+ip-1) = 0.d0
  ra(inzpar+ip-1) = 0.d0
enddo

call usladp                                                       &
!==========
 ( ifinia , ifinra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   nbpmax , nvp    , nvp1   , nvep   , nivep  ,                   &
   ntersl , nvlsta , nvisbr ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &
   itepa  ,                                                       &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtpa   , propce , propfa , propfb ,                   &
   ettp   , ettpa  , tepa   , statis ,                            &
   taup   , tlag   , piil   ,                                     &
   vagaus , gradpr , gradvf ,                                     &
   romp   ,                                                       &
   ra(idppar) ,ra(inxpar) , ra(inypar) , ra(inzpar) ,             &
   rdevel , rtuser , ra     )

!===============================================================================
! 2. FORCES DE VAN DER WAALS
!    Pour etre homogene a des m/s2 on divise par la masse
!===============================================================================

do ip = 1,nbpart

! Force = -A/6 dp/2 /D**2

  if ( ra(idppar+ip-1) .gt. dparmn ) then

    val = (cstham*ettp(ip,jdp)/2.d0)                              &
         /(6.d0*ra(idppar+ip-1)*ra(idppar+ip-1))

    dnorm = sqrt( ra(inxpar+ip-1)*ra(inxpar+ip-1)                 &
                 +ra(inypar+ip-1)*ra(inypar+ip-1)                 &
                 +ra(inzpar+ip-1)*ra(inzpar+ip-1) )

! Attention la normale est oriente du fluide vers l'exterieur

    aa = dnorm*ettp(ip,jmp)

    fextla(ip,1) = fextla(ip,1) + val*ra(inxpar+ip-1) /aa
    fextla(ip,2) = fextla(ip,2) + val*ra(inypar+ip-1) /aa
    fextla(ip,3) = fextla(ip,3) + val*ra(inzpar+ip-1) /aa

  endif

enddo

!===============================================================================
! 3. FORCES ELECTROSTATIQUES
!    Pour etre homogene a des m/s2 on divise par la masse
!===============================================================================

do ip = 1,nbpart

  iel = itepa(ip,jisor)

! Calcul de la temperature du fluide en fonction du type
! d'ecoulement

  if ( ra(idppar+ip-1) .gt. dparmn ) then

    if ( ippmod(icp3pl).ge.0 .or.                                 &
         ippmod(icpl3c).ge.0      ) then

      tempf = propce(iel,ipproc(itemp1)) - tkelvi

    else if ( ippmod(icod3p).ge.0 .or.                            &
              ippmod(icoebu).ge.0 .or.                            &
              ippmod(ielarc).ge.0 .or.                            &
              ippmod(ieljou).ge.0      ) then

      tempf = propce(iel,ipproc(itemp)) - tkelvi

    else if ( iscsth(iscalt(iphas)).eq.-1 ) then
      tempf = rtp(iel,isca(iscalt(iphas)))

    else if ( iscsth(iscalt(iphas)).eq.1 ) then
      tempf = rtp(iel,isca(iscalt(iphas))) - tkelvi

    else if ( iscsth(iscalt(iphas)).eq.2 ) then
      mode = 1
      call usthht (mode, rtp(iel,isca(iscalt(iphas))), tempf)
      !==========
    else
      tempf = t0(iphas)
    endif

! FORCE :

! Longueur de Debye

    if (fion .ne. 0 .and. cstfar .gt. 0.d0) then
      debye  = sqrt( (epseau*epsvid*rr*tempf)                     &
              /(2000.d0*cstfar*cstfar*fion) )
    else
      write(nfecra,9001) fion,cstfar
      call csexit(1)
    endif

    if ( debye .gt. 0.d0 ) then
      debye = sqrt(debye)
    else
      write(nfecra,9002) ip,debye,tempf,fion,epseau,epsvid
      call csexit(1)
    endif

    val = -4.d0*pi*epseau*epsvid*phi1*phi2*(ettp(ip,jdp)/2.d0)    &
            *exp(-ra(idppar+ip-1)/debye)                          &
            /debye

    dnorm = sqrt( ra(inxpar+ip-1)*ra(inxpar+ip-1)                 &
                 +ra(inypar+ip-1)*ra(inypar+ip-1)                 &
                 +ra(inzpar+ip-1)*ra(inzpar+ip-1) )

! Attention la normale est oriente du fluide vers l'exterieur

    fextla(ip,1)= fextla(ip,1)+val*ra(inxpar+ip-1)                &
                 /(dnorm*ettp(ip,jmp))
    fextla(ip,2)= fextla(ip,2)+val*ra(inypar+ip-1)                &
                 /(dnorm*ettp(ip,jmp))
    fextla(ip,3)= fextla(ip,3)+val*ra(inzpar+ip-1)                &
                 /(dnorm*ettp(ip,jmp))

  endif

enddo

!===============================================================================
! 4. FORCES D'ADHESION
!===============================================================================

if ( dcoup .gt. 0.d0 ) then
  gamasv = cstham/(24.d0*pi*dcoup*dcoup)
else
  write(nfecra,9010) dcoup
  call csexit(1)
endif

do ip = 1,nbpart

! Force = 3*PI*(Dp/2)*Gamma_SV + SIG2*PI*(Dp/2)/Eps0

  if ( ra(idppar+ip-1) .le. dparmn ) then

    val = 3.d0*pi*(ettp(ip,jdp)/2.d0)*gamasv                      &
        + sigch*sigch*pi*(ettp(ip,jdp)/2.d0)/epsvid

    dnorm = sqrt( ra(inxpar+ip-1)*ra(inxpar+ip-1)                 &
                 +ra(inypar+ip-1)*ra(inypar+ip-1)                 &
                 +ra(inzpar+ip-1)*ra(inzpar+ip-1) )

! Attention la normale est oriente du fluide vers l'exterieur

    aa = dnorm*ettp(ip,jmp)

    fextla(ip,1)= fextla(ip,1)+val*ra(inxpar+ip-1) /aa
    fextla(ip,2)= fextla(ip,2)+val*ra(inypar+ip-1) /aa
    fextla(ip,3)= fextla(ip,3)+val*ra(inzpar+ip-1) /aa

  endif

enddo

!==============================================================================

!--------
! FORMATS
!--------

 9001 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========                                               ',/,&
'@    LA VALEUR DE LA FORCE IONIQUE EST NULLE                 ',/,&
'@ OU LA VALEUR DE LA CONSTANTE DE FARADET EST NEGATIVE       ',/,&
'@                                          OU NULLEE         ',/,&
'@                                                            ',/,&
'@       FORCE IONIQUE                : ',G15.7                ,/,&
'@       CSTE DE FARADET              : ',G15.7                ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier les valeurs de FION et de CSTFAR                 ',/,&
'@                              dans la subroutine USLAG1.    ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 9002 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========                                               ',/,&
'@    LA VALEUR DE L''EPAISSEUR DE LA DOUBLE COUCHE           ',/,&
'@    EST NEGATIVE OU NULLE :                                 ',/,&
'@       NUMERO DE PARTICULE          : ',I10                  ,/,&
'@       EPAISSEUR                    : ',G15.7                ,/,&
'@       TEMPERATURE                  : ',G15.7                ,/,&
'@       FORCE IONIQUE                : ',G15.7                ,/,&
'@       CSTE DIELECTIQUE DU VIDE     : ',G15.7                ,/,&
'@       CSTE DIELECTIQUE DE L''EAU   : ',G15.7                ,/,&
'@       CSTE DE FARADET              : ',G15.7                ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier les valeurs des CSTES dans la subroutine USLAG1. ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 9010 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========                                               ',/,&
'@    LA VALEUR DE DISTANCE DE COUPURE EST                    ',/,&
'@    EST NEGATIVE OU NULLE               : ',G15.7            ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier les valeurs des CSTES dans la subroutine USLAG1. ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

!----
! FIN
!----

end subroutine
