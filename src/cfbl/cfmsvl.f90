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

subroutine cfmsvl &
!================

 ( idbia0 , idbra0 ,                                              &
   nvar   , nscal  , ncepdp , ncesmp ,                            &
   iscal  ,                                                       &
   icepdc , icetsm , itypsm ,                                     &
   ia     ,                                                       &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  , ckupdc , smacel ,                            &
   viscf  , viscb  ,                                              &
   dam    , xam    ,                                              &
   drtp   , smbrs  , rovsdt ,                                     &
   w1     , w2     , w3     , w4     , w5     , w6     ,          &
   w7     , w8     , w9     , w10    , w11    , w12    ,          &
   wflmas , wflmab ,                                              &
   coefu  ,                                                       &
   ra     )

!===============================================================================
! FONCTION :
! ----------

! RESOLUTION DES EQUATIONS CONVECTION DIFFUSION TERME SOURCE
!   POUR LA MASSE VOLUMIQUE SUR UN PAS DE TEMPS
!   (ALGORITHME COMPRESSIBLE EN RHO, U, E)

!-------------------------------------------------------------------------------
!ARGU                             ARGUMENTS
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! ncepdp           ! i  ! <-- ! number of cells with head loss                 !
! ncesmp           ! i  ! <-- ! number of cells with mass source term          !
! iscal            ! i  ! <-- ! scalar number                                  !
! itspdv           ! e  ! <-- ! calcul termes sources prod et dissip           !
!                  !    !     !  (0 : non , 1 : oui)                           !
! icepdc(ncelet    ! te ! <-- ! numero des ncepdp cellules avec pdc            !
! icetsm(ncesmp    ! te ! <-- ! numero des cellules a source de masse          !
! itypsm           ! te ! <-- ! type de source de masse pour les               !
! (ncesmp,nvar)    !    !     !  variables (cf. ustsma)                        !
! ia(*)            ! ia ! --- ! main integer work array                        !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! tslagr           ! tr ! <-- ! terme de couplage retour du                    !
!(ncelet,*)        !    !     !     lagrangien                                 !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
! ckupdc           ! tr ! <-- ! tableau de travail pour pdc                    !
!  (ncepdp,6)      !    !     !                                                !
! smacel           ! tr ! <-- ! valeur des variables associee a la             !
! (ncesmp,*   )    !    !     !  source de masse                               !
!                  !    !     !  pour ivar=ipr, smacel=flux de masse           !
! viscf(nfac)      ! tr ! --- ! visc*surface/dist aux faces internes           !
! viscb(nfabor     ! tr ! --- ! visc*surface/dist aux faces de bord            !
! dam(ncelet       ! tr ! --- ! tableau de travail pour matrice                !
! xam(nfac,*)      ! tr ! --- ! tableau de travail pour matrice                !
! drtp(ncelet      ! tr ! --- ! tableau de travail pour increment              !
! smbrs(ncelet     ! tr ! --- ! tableau de travail pour sec mem                !
! rovsdt(ncelet    ! tr ! --- ! tableau de travail pour terme instat           !
! w1..12(ncelet    ! tr ! --- ! tableau de travail                             !
! wflmas(nfac)     ! tr ! --- ! tableau de w flux de masse aux faces           !
! wflmab(nfabor    ! tr ! --- ! tableau de w flux de masse aux bords           !
! coefu(nfabo,3    ! tr ! --- ! tableau de travail cl de la qdm                !
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
use entsor
use optcal
use cstphy
use cstnum
use pointe
use parall
use period
use ppppar
use ppthch
use ppincl
use mesh

!===============================================================================

implicit none

! Arguments

integer          idbia0 , idbra0
integer          nvar   , nscal
integer          ncepdp , ncesmp
integer          iscal

integer          icepdc(ncepdp)
integer          icetsm(ncesmp), itypsm(ncesmp,nvar)
integer          ia(*)

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision ckupdc(ncepdp,6), smacel(ncesmp,nvar)
double precision viscf(nfac), viscb(nfabor)
double precision dam(ncelet), xam(nfac,2)
double precision drtp(ncelet), smbrs(ncelet)
double precision rovsdt(ncelet)
double precision w1(ncelet) , w2(ncelet) , w3(ncelet)
double precision w4(ncelet) , w5(ncelet) , w6(ncelet)
double precision w7(ncelet) , w8(ncelet) , w9(ncelet)
double precision w10(ncelet), w11(ncelet), w12(ncelet)
double precision wflmas(nfac), wflmab(nfabor)
double precision coefu(nfabor,3)
double precision ra(*)

! Local variables

character*80     chaine
integer          idebia, idebra
integer          ivar
integer          ifac  , iel
integer          init  , inc   , iccocg, isqrt , ii, jj, iii
integer          iclvar, iclvaf
integer          iflmas, iflmab
integer          ippvar, ipp   , iphydp
integer          nswrgp, imligp, iwarnp
integer          istatp, iconvp, idiffp, ireslp, ndircp, nitmap
integer          nswrsp, ircflp, ischcp, isstpp, iescap
integer          imgrp , ncymxp, nitmfp
double precision epsrgp, climgp, extrap, blencp, epsilp
double precision epsrsp
double precision sclnor

integer          iccfth, imodif
integer          iij
integer          iwfabg, iwfbbg
double precision dijpfx, dijpfy, dijpfz, pnd
double precision diipfx, diipfy, diipfz, djjpfx, djjpfy, djjpfz
double precision diipbx, diipby, diipbz
double precision pip   , pjp   , thetv, relaxp

!===============================================================================
!===============================================================================
! 1. INITIALISATION
!===============================================================================

idebia = idbia0
idebra = idbra0

iwfabg = idebra
iwfbbg = iwfabg+nfac
idebra = iwfbbg+nfabor
call rasize('cfmsvl',idebra)

! --- Numero de variable de calcul et de post associe au scalaire traite
ivar   = isca(iscal)
ippvar = ipprtp(ivar)

! --- Numero des conditions aux limites
iclvar = iclrtp(ivar,icoef)
iclvaf = iclrtp(ivar,icoeff)

! --- Flux de masse associe a l'energie
iflmas = ipprof(ifluma(isca(ienerg)))
iflmab = ipprob(ifluma(isca(ienerg)))

chaine = nomvar(ippvar)

if(iwarni(ivar).ge.1) then
  write(nfecra,1000) chaine(1:8)
endif

!===============================================================================
! 2. TERMES SOURCES
!===============================================================================

! --> Initialisation

do iel = 1, ncel
  smbrs(iel) = 0.d0
enddo
do iel = 1, ncel
  rovsdt(iel) = 0.d0
enddo


!     TERME SOURCE DE MASSE
!     =====================

if (ncesmp.gt.0) then
  do ii = 1, ncesmp
    iel = icetsm(ii)
    smbrs(iel) = smbrs(iel) + smacel(iel,ipr)*volume(iel)
  enddo
endif


!     TERME INSTATIONNAIRE
!     ====================

do iel = 1, ncel
  rovsdt(iel) = rovsdt(iel) + istat(ivar)*(volume(iel)/dt(iel))
enddo

!===============================================================================
! 3. CALCUL DU "FLUX DE MASSE" ET DE LA "VISCOSITE" AUX FACES
!===============================================================================

!     Ici VISCF et VISCB sont deux tableaux de travail.
!     On calcule WFLMAS et WFLMAB, WFABGS , WFBBGS

call cfmsgs                                                       &
!==========
 ( idebia , idebra ,                                              &
   nvar   , nscal  , ncepdp , ncesmp ,                            &
   iscal  ,                                                       &
   icepdc , icetsm , itypsm ,                                     &
   ia     ,                                                       &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  , ckupdc , smacel ,                            &
   wflmas , wflmab , ra(iwfabg) , ra(iwfbbg) ,                    &
!        ------   ------   ------   ------
   w1     , w2     , w3     , w4     , w5     , w6     ,          &
   w7     , w8     , w9     , w10    , w11    , w12    ,          &
   viscf  , viscb  , coefu  , xam    ,                            &
   ra     )


!     Calcul du gradient de rho pour la reconstruction de rho
!       (contribution du terme convectif)

ircflp = ircflu(ivar)

if(ircflp.gt.0) then

  inc    = 1
  iccocg = 1
  nswrgp = nswrgr(ivar)
  imligp = imligr(ivar)
  iphydp = 0
  iwarnp = iwarni(ivar)
  epsrgp = epsrgr(ivar)
  climgp = climgr(ivar)
  extrap = extrag(ivar)

  call grdcel                                                     &
  !==========
 ( idebia , idebra ,                                              &
   ivar   , imrgra , inc    , iccocg , nswrgp , imligp , iphydp , &
   iwarnp , nfecra , epsrgp , climgp , extrap ,                   &
   ia     ,                                                       &
   w1     , w1     , w1     ,                                     &
   rtpa(1,ivar)    ,                                              &
   coefa(1,iclrtp(ivar,icoef)) , coefb(1,iclrtp(ivar,icoef)) ,    &
   w1     , w2     , w3     ,                                     &
!        ------   ------   ------
   w4     , w5     , w6     ,                                     &
   ra     )

  else
    do ii = 1, ncelet
      w1(ii) = 0.d0
      w2(ii) = 0.d0
      w3(ii) = 0.d0
    enddo
  endif


!     Au bord, on écrase WFLMAB pour imposer le débit souhaité.
!     Si on explicite le terme de convection (choix retenu par défaut,
!       seul testé), pas de problème.
!     Si on implicite le terme de convection, on n'impose pas
!       nécessairement le bon débit (cela dépend du traitement de
!       la convection au bord dans bilsc2 et de la condition à la
!       limite sur rho : ainsi, si l'on se place sur une sortie pour
!       laquelle la condition n'est pas valeur bord = valeur interne,
!       la convection étant traitée en upwind, c'est la valeur interne
!       de rho qui intervient dans le flux de bord et non pas
!       la valeur de bord comme supposé ci-dessous.

if(iconv(ivar).le.0) then

  if(ircflp.gt.0) then

    do ifac = 1, nfac
      ii = ifacel(1,ifac)
      jj = ifacel(2,ifac)

      dijpfx = dijpf(1,ifac)
      dijpfy = dijpf(2,ifac)
      dijpfz = dijpf(3,ifac)

      pnd   = pond(ifac)

      diipfx = cdgfac(1,ifac) - (xyzcen(1,ii)+                    &
               (1.d0-pnd) * dijpfx)
      diipfy = cdgfac(2,ifac) - (xyzcen(2,ii)+                    &
               (1.d0-pnd) * dijpfy)
      diipfz = cdgfac(3,ifac) - (xyzcen(3,ii)+                    &
               (1.d0-pnd) * dijpfz)
      djjpfx = cdgfac(1,ifac) -  xyzcen(1,jj)+                    &
                   pnd  * dijpfx
      djjpfy = cdgfac(2,ifac) -  xyzcen(2,jj)+                    &
                   pnd  * dijpfy
      djjpfz = cdgfac(3,ifac) -  xyzcen(3,jj)+                    &
                   pnd  * dijpfz

      pip = rtpa(ii,ivar)                                         &
           + ircflp*(w1(ii)*diipfx+w2(ii)*diipfy+w3(ii)*diipfz)
      pjp = rtpa(jj,ivar)                                         &
           + ircflp*(w1(jj)*djjpfx+w2(jj)*djjpfy+w3(jj)*djjpfz)

      wflmas(ifac) = -0.5d0*                                      &
           ( pip          *(wflmas(ifac)+abs(wflmas(ifac)))       &
           + pjp          *(wflmas(ifac)-abs(wflmas(ifac))))
    enddo

  else
    do ifac = 1, nfac
      ii = ifacel(1,ifac)
      jj = ifacel(2,ifac)
      wflmas(ifac) = -0.5d0*                                      &
           ( rtpa(ii,ivar)*(wflmas(ifac)+abs(wflmas(ifac)))       &
           + rtpa(jj,ivar)*(wflmas(ifac)-abs(wflmas(ifac))))
    enddo
  endif

  do ifac = 1, nfabor
    wflmab(ifac) = -propfb(ifac,iflmab)
  enddo

  init = 0
  call divmas(ncelet,ncel,nfac,nfabor,init,nfecra,                &
              ifacel,ifabor,wflmas,wflmab,smbrs)

  do ifac = 1, nfac
    ra(iwfabg+ifac-1) = - ra(iwfabg+ifac-1)
  enddo
  do ifac = 1, nfabor
    ra(iwfbbg+ifac-1) = - ra(iwfbbg+ifac-1)
  enddo
  init = 0
  call divmas(ncelet,ncel,nfac,nfabor,init,nfecra,                &
              ifacel,ifabor,ra(iwfabg),ra(iwfbbg),smbrs)

else

  if(ircflp.gt.0) then

    do ifac = 1, nfabor
      ii = ifabor(ifac)

      diipbx = diipb(1,ifac)
      diipby = diipb(2,ifac)
      diipbz = diipb(3,ifac)

      pip = rtpa(ii,ivar)                                         &
           +ircflp*(w1(ii)*diipbx+w2(ii)*diipby+w3(ii)*diipbz)

      wflmab(ifac) = -propfb(ifac,iflmab)/                        &
           ( coefa(ifac,iclrtp(ivar,icoef))                       &
           + coefb(ifac,iclrtp(ivar,icoef))*pip           )
    enddo

  else
    do ifac = 1, nfabor
      ii = ifabor(ifac)
      wflmab(ifac) = -propfb(ifac,iflmab)/                        &
           ( coefa(ifac,iclrtp(ivar,icoef))                       &
           + coefb(ifac,iclrtp(ivar,icoef))*rtpa(ii,ivar) )
    enddo
  endif

endif


!     On calcule VISCF et VISCB

call cfmsvs                                                       &
!==========
 ( idebia , idebra ,                                              &
   nvar   , nscal  ,                                                                                 &
   iscal  ,                                                       &
   ia     ,                                                       &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   viscf  , viscb  ,                                              &
!        ------   ------
   w1     , w9     , w10    ,                                     &
   ra     )


!     On annule la viscosité au bord afin que la contribution
!       au flux de bord soit exactement WFLMAB (dans lequel on a mis
!       le flux de masse souhaité). Si on n'annule pas, on risque
!       d'obtenir des contributions non nulles de la partie diffusive,
!       sauf si on a imposé une condition de Neumann homogene sur rho.
!       Pour le moment, on prefere prendre des precautions (la
!       modification de VISCB se traduit simplement par la modification
!       de la matrice portant sur les incréments, ou encore par une
!       une condition de Neumann homogène sur les incréments).

do ifac = 1, nfabor
  viscb (ifac) = 0.d0
enddo

!===============================================================================
! 4. RESOLUTION
!===============================================================================

istatp = istat (ivar)
iconvp = iconv (ivar)
idiffp = idiff (ivar)
ireslp = iresol(ivar)
ndircp = ndircl(ivar)
nitmap = nitmax(ivar)
nswrsp = nswrsm(ivar)
nswrgp = nswrgr(ivar)
imligp = imligr(ivar)
ircflp = ircflu(ivar)
ischcp = ischcv(ivar)
isstpp = isstpc(ivar)
iescap = 0
imgrp  = imgr  (ivar)
ncymxp = ncymax(ivar)
nitmfp = nitmgf(ivar)
ipp    = ippvar
iwarnp = iwarni(ivar)
blencp = blencv(ivar)
epsilp = epsilo(ivar)
epsrsp = epsrsm(ivar)
epsrgp = epsrgr(ivar)
climgp = climgr(ivar)
extrap = extrag(ivar)
relaxp = relaxv(ivar)
thetv  = thetav(ivar)

call codits                                                       &
!==========
 ( idebia , idebra ,                                              &
   nvar   , nscal  ,                                                                                 &
   idtvar , ivar   , iconvp , idiffp , ireslp , ndircp , nitmap , &
   imrgra , nswrsp , nswrgp , imligp , ircflp ,                   &
   ischcp , isstpp , iescap ,                                     &
   imgrp  , ncymxp , nitmfp , ipp    , iwarnp ,                   &
   blencp , epsilp , epsrsp , epsrgp , climgp , extrap ,          &
   relaxp , thetv  ,                                              &
   ia     ,                                                       &
   rtpa(1,ivar)    , rtpa(1,ivar)    ,                            &
   coefa(1,iclvar) , coefb(1,iclvar) ,                            &
                     coefa(1,iclvaf) , coefb(1,iclvaf) ,          &
                     wflmas          , wflmab          ,          &
   viscf  , viscb  , viscf  , viscb  ,                            &
   rovsdt , smbrs  , rtp(1,ivar)     ,                            &
   dam    , xam    , drtp   ,                                     &
   w1     , w2     , w3     , w4     , w5     ,                   &
   w6     , w7     , w8     , w9     ,                            &
   ra     )

!===============================================================================
! 5. IMPRESSIONS ET CLIPPINGS
!===============================================================================

! --- Clipping aux bornes définies par l'utilisateur ou par defaut
!       (par défaut, pas de borne contraignate)

!     Valeur bidon
iii = 1

call clpsca                                                       &
!==========
 ( ncelet , ncel   , nvar   , nscal  , iscal  ,                   &
   propce , rtp(1,iii)      , rtp    )

! --- Traitement utilisateur pour gestion plus fine des bornes
!       et actions correctives éventuelles.
  iccfth = -2
  imodif = 0
  call uscfth                                                     &
  !==========
 ( nvar   , nscal  ,                                                                                 &
   iccfth , imodif ,                                              &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   w7     , w8     , w9     , w10    )


! --- Bilan explicite (voir codits : on enleve l'increment)

if (iwarni(ivar).ge.2) then
  do iel = 1, ncel
    smbrs(iel) = smbrs(iel)                                       &
            - istat(ivar)*(volume(iel)/dt(iel))                   &
                *(rtp(iel,ivar)-rtpa(iel,ivar))                   &
                * max(0,min(nswrsm(ivar)-2,1))
  enddo
  isqrt = 1
  call prodsc(ncelet,ncel,isqrt,smbrs,smbrs,sclnor)
  write(nfecra,1200)chaine(1:8) ,sclnor
endif

!===============================================================================
! 6. COMMUNICATION DE RHO
!===============================================================================

if (irangp.ge.0.or.iperio.eq.1) then
  call synsca(rtp(1,ivar))
  !==========
endif

!     On ne remplit pas PROPCE et PROPFB ici, car on veut disposer de
!       rho à l'instant précédent pour résoudre en qdm et energie
!     On modifiera PROPCE et PROPFB après resolution de l'energie.


!===============================================================================
! 7. CALCUL DU FLUX DE MASSE ACOUSTIQUE AUX FACES
!===============================================================================

! Ce flux est stocke en tant que flux de masse associe a l'energie
inc    = 1
iccocg = 1

call cfbsc3                                                       &
!==========
 ( idebia , idebra ,                                              &
   nvar   , nscal  ,                                                                                 &
   ivar   , iconvp , idiffp , nswrgp , imligp , ircflp ,          &
   ischcp , isstpp , inc    , imrgra , iccocg ,                   &
   ipp    , iwarnp ,                                              &
   blencp , epsrgp , climgp , extrap ,                            &
   ia     ,                                                       &
   rtp(1,ivar)     , coefa(1,iclvar) , coefb(1,iclvar) ,          &
                     coefa(1,iclvaf) , coefb(1,iclvaf) ,          &
                     wflmas          , wflmab          ,          &
   viscf  , viscb  ,                                              &
   propfa(1,iflmas), propfb(1,iflmab),                            &
!        ----------------  ----------------
   w1     , w2     , w3     , w4     , w5     , w6     ,          &
   ra     )


!     Si ICONV = 0, le terme convectif n'est pas traite par CFBSC3
!       il faut donc le rajouter a la main
!     Noter egalement que si ICONV = 0, PROPFB contient zero au sortir
!       de cfbsc3 et qu'on lui ajoute donc le flux de masse WFLMAB
!       (cohérent avec le flux imposé au bord et avec un signe négatif,
!        car WFLMAB etait, ci-dessus, utilise au second membre)
if(iconv(ivar).le.0) then
  do ifac = 1, nfac
    propfa(ifac,iflmas) = propfa(ifac,iflmas) - wflmas(ifac)      &
                                              - ra(iwfabg+ifac-1)
  enddo
  do ifac = 1, nfabor
    propfb(ifac,iflmab) = propfb(ifac,iflmab) - wflmab(ifac)      &
                                              - ra(iwfbbg+ifac-1)
  enddo
endif


!===============================================================================
! 8. ACTUALISATION DE LA PRESSION
!===============================================================================
!                               Pred      n+1  n
! On utilise l'equation d'etat P    =P(rho   ,e )

! --- Calcul de P au centre des cellules et actualisation de RTP
if(igrdpp.gt.0) then

  iccfth = 24
  imodif = 0
  call uscfth                                                     &
  !==========
 ( nvar   , nscal  ,                                                                                 &
   iccfth , imodif ,                                              &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   rtp(1,ipr)        , w8     , w9     , w10    )

!===============================================================================
! 9. COMMUNICATION DE LA PRESSION
!===============================================================================

  if (irangp.ge.0.or.iperio.eq.1) then
    call synsca(rtp(1,ipr))
    !==========
  endif

endif

!     Pas de vérification ni de clipping de la pression, car on
!       considere que la masse volumique et l'énergie ont été vérifiees,
!       sont correctes et donc que la pression l'est aussi (elle
!       vient d'être calculée avec un sous-pgm utilisateur).

!--------
! FORMATS
!--------

 1000 format(/,                                                   &
'   ** RESOLUTION POUR LA VARIABLE ',A8                        ,/,&
'      ---------------------------                            ',/)
 1200 format(1X,A8,' : BILAN EXPLICITE = ',E14.5)

!----
! FIN
!----

return
end subroutine
