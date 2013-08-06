!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2013 EDF S.A.
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

subroutine cfmsvl &
!================

 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   iscal  ,                                                       &
   icepdc , icetsm , itypsm ,                                     &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  , ckupdc , smacel )

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
use parall
use period
use ppppar
use ppthch
use ppincl
use mesh
use field

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal
integer          ncepdp , ncesmp
integer          iscal

integer          icepdc(ncepdp)
integer          icetsm(ncesmp), itypsm(ncesmp,nvar)

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision ckupdc(ncepdp,6), smacel(ncesmp,nvar)

! Local variables

character*80     chaine
integer          ivar
integer          ifac  , iel
integer          init  , inc   , iccocg, isqrt , ii, jj, iii
integer          iclvar, iclvaf
integer          iflmas, iflmab
integer          ippvar, ipp
integer          nswrgp, imligp, iwarnp
integer          istatp, iconvp, idiffp, ireslp, ndircp, nitmap
integer          nswrsp, ircflp, ischcp, isstpp, iescap
integer          imgrp , ncymxp, nitmfp
double precision epsrgp, climgp, extrap, blencp, epsilp
double precision epsrsp
double precision sclnor

integer          iccfth, imodif
integer          iij
integer          imucpp, idftnp, iswdyp
double precision dijpfx, dijpfy, dijpfz, pnd
double precision diipfx, diipfy, diipfz, djjpfx, djjpfy, djjpfz
double precision diipbx, diipby, diipbz
double precision pip   , pjp   , thetv, relaxp

double precision rvoid(1)

double precision, allocatable, dimension(:) :: viscf, viscb
double precision, allocatable, dimension(:) :: smbrs, rovsdt
double precision, allocatable, dimension(:,:) :: grad
double precision, allocatable, dimension(:) :: w1
double precision, allocatable, dimension(:) :: w7, w8, w9
double precision, allocatable, dimension(:) :: w10
double precision, allocatable, dimension(:) :: wflmas, wflmab
double precision, allocatable, dimension(:) :: wfabg, wfbbg
double precision, allocatable, dimension(:) :: dpvar
double precision, dimension(:), pointer :: imasfl, bmasfl

!===============================================================================
!===============================================================================
! 1. INITIALISATION
!===============================================================================

! Allocate temporary arrays for the mass resolution
allocate(viscf(nfac), viscb(nfabor))
allocate(smbrs(ncelet), rovsdt(ncelet))
allocate(wflmas(nfac), wflmab(nfabor))
allocate(wfabg(nfac), wfbbg(nfabor))

! Allocate work arrays
allocate(w1(ncelet))
allocate(w7(ncelet), w8(ncelet), w9(ncelet))
allocate(w10(ncelet))
allocate(dpvar(ncelet))

! --- Numero de variable de calcul et de post associe au scalaire traite
ivar   = isca(iscal)
ippvar = ipprtp(ivar)

! --- Numero des conditions aux limites
iclvar = iclrtp(ivar,icoef)
iclvaf = iclrtp(ivar,icoeff)

! --- Flux de masse associe a l'energie
call field_get_key_int(ivarfl(isca(ienerg)), kimasf, iflmas)
call field_get_key_int(ivarfl(isca(ienerg)), kbmasf, iflmab)
call field_get_val_s(iflmas, imasfl)
call field_get_val_s(iflmab, bmasfl)

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
 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   iscal  ,                                                       &
   icepdc , icetsm , itypsm ,                                     &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  , ckupdc , smacel ,                            &
   wflmas , wflmab , wfabg  , wfbbg  ,                            &
!        ------   ------   ------   ------
   viscf  , viscb  )


!     Calcul du gradient de rho pour la reconstruction de rho
!       (contribution du terme convectif)

! Allocate a work array
allocate(grad(ncelet,3))

ircflp = ircflu(ivar)

if(ircflp.gt.0) then

  inc    = 1
  iccocg = 1
  nswrgp = nswrgr(ivar)
  imligp = imligr(ivar)
  iwarnp = iwarni(ivar)
  epsrgp = epsrgr(ivar)
  climgp = climgr(ivar)
  extrap = extrag(ivar)

  call grdcel                                                     &
  !==========
 ( ivar   , imrgra , inc    , iccocg , nswrgp , imligp ,          &
   iwarnp , nfecra , epsrgp , climgp , extrap ,                   &
   rtpa(1,ivar)    ,                                              &
   coefa(1,iclrtp(ivar,icoef)) , coefb(1,iclrtp(ivar,icoef)) ,    &
   grad   )

  else
    do ii = 1, ncelet
      grad(ii,1) = 0.d0
      grad(ii,2) = 0.d0
      grad(ii,3) = 0.d0
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

      diipfx = cdgfac(1,ifac) - (xyzcen(1,ii) + (1.d0-pnd) * dijpfx)
      diipfy = cdgfac(2,ifac) - (xyzcen(2,ii) + (1.d0-pnd) * dijpfy)
      diipfz = cdgfac(3,ifac) - (xyzcen(3,ii) + (1.d0-pnd) * dijpfz)
      djjpfx = cdgfac(1,ifac) -  xyzcen(1,jj) + pnd  * dijpfx
      djjpfy = cdgfac(2,ifac) -  xyzcen(2,jj) + pnd  * dijpfy
      djjpfz = cdgfac(3,ifac) -  xyzcen(3,jj) + pnd  * dijpfz

      pip = rtpa(ii,ivar)                                         &
           + ircflp*(grad(ii,1)*diipfx+grad(ii,2)*diipfy+grad(ii,3)*diipfz)
      pjp = rtpa(jj,ivar)                                         &
           + ircflp*(grad(jj,1)*djjpfx+grad(jj,2)*djjpfy+grad(jj,3)*djjpfz)

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
    wflmab(ifac) = -bmasfl(ifac)
  enddo

  init = 0
  call divmas(ncelet,ncel,nfac,nfabor,init,nfecra,                &
              ifacel,ifabor,wflmas,wflmab,smbrs)

  do ifac = 1, nfac
    wfabg(ifac) = - wfabg(ifac)
  enddo
  do ifac = 1, nfabor
    wfbbg(ifac) = - wfbbg(ifac)
  enddo
  init = 0
  call divmas(ncelet,ncel,nfac,nfabor,init,nfecra,                &
              ifacel,ifabor,wfabg,wfbbg,smbrs)

else

  if(ircflp.gt.0) then

    do ifac = 1, nfabor
      ii = ifabor(ifac)

      diipbx = diipb(1,ifac)
      diipby = diipb(2,ifac)
      diipbz = diipb(3,ifac)

      pip = rtpa(ii,ivar)                                         &
           +ircflp*(grad(ii,1)*diipbx+grad(ii,2)*diipby+grad(ii,3)*diipbz)

      wflmab(ifac) = -bmasfl(ifac)/                               &
           ( coefa(ifac,iclrtp(ivar,icoef))                       &
           + coefb(ifac,iclrtp(ivar,icoef))*pip           )
    enddo

  else
    do ifac = 1, nfabor
      ii = ifabor(ifac)
      wflmab(ifac) = -bmasfl(ifac)/                               &
           ( coefa(ifac,iclrtp(ivar,icoef))                       &
           + coefb(ifac,iclrtp(ivar,icoef))*rtpa(ii,ivar) )
    enddo
  endif

endif

! Free memory
deallocate(grad)

!     On calcule VISCF et VISCB

call cfmsvs                                                       &
!==========
 ( nvar   , nscal  ,                                              &
   iscal  ,                                                       &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   viscf  , viscb  ,                                              &
!        ------   ------
   w1     , w9     , w10    )


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
imucpp = 0
idftnp = idften(ivar)
iswdyp = iswdyn(ivar)
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

call codits &
!==========
 ( nvar   , nscal  ,                                              &
   idtvar , ivar   , iconvp , idiffp , ireslp , ndircp , nitmap , &
   imrgra , nswrsp , nswrgp , imligp , ircflp ,                   &
   ischcp , isstpp , iescap , imucpp , idftnp , iswdyp ,          &
   imgrp  , ncymxp , nitmfp , ipp    , iwarnp ,                   &
   blencp , epsilp , epsrsp , epsrgp , climgp , extrap ,          &
   relaxp , thetv  ,                                              &
   rtpa(1,ivar)    , rtpa(1,ivar)    ,                            &
   coefa(1,iclvar) , coefb(1,iclvar) ,                            &
   coefa(1,iclvaf) , coefb(1,iclvaf) ,                            &
   wflmas          , wflmab          ,                            &
   viscf  , viscb  , rvoid  , viscf  , viscb  , rvoid  ,          &
   rvoid  , rvoid  ,                                              &
   rovsdt , smbrs  , rtp(1,ivar)     , dpvar  ,                   &
   rvoid  , rvoid  )

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
  call cfther                                                     &
  !==========
 ( nvar   , nscal  ,                                              &
   iccfth , imodif ,                                              &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
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
  call prodsc(ncel,isqrt,smbrs,smbrs,sclnor)
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
 ( nvar   , nscal  ,                                              &
   ivar   , iconvp , idiffp , nswrgp , imligp , ircflp ,          &
   ischcp , isstpp , inc    , imrgra , iccocg ,                   &
   ipp    , iwarnp ,                                              &
   blencp , epsrgp , climgp , extrap ,                            &
   rtp(1,ivar)     , coefa(1,iclvar) , coefb(1,iclvar) ,          &
                     coefa(1,iclvaf) , coefb(1,iclvaf) ,          &
                     wflmas          , wflmab          ,          &
   viscf  , viscb  ,                                              &
   imasfl , bmasfl)


!     Si ICONV = 0, le terme convectif n'est pas traite par CFBSC3
!       il faut donc le rajouter a la main
!     Noter egalement que si ICONV = 0, PROPFB contient zero au sortir
!       de cfbsc3 et qu'on lui ajoute donc le flux de masse WFLMAB
!       (cohérent avec le flux imposé au bord et avec un signe négatif,
!        car WFLMAB etait, ci-dessus, utilise au second membre)
if(iconv(ivar).le.0) then
  do ifac = 1, nfac
    imasfl(ifac) = imasfl(ifac) - wflmas(ifac) - wfabg(ifac)
  enddo
  do ifac = 1, nfabor
    bmasfl(ifac) = bmasfl(ifac) - wflmab(ifac) - wfbbg(ifac)
  enddo
endif

! Free memory
deallocate(wfabg, wfbbg)

!===============================================================================
! 8. ACTUALISATION DE LA PRESSION
!===============================================================================
!                               Pred      n+1  n
! On utilise l'equation d'etat P    =P(rho   ,e )

! --- Calcul de P au centre des cellules et actualisation de RTP
if(igrdpp.gt.0) then

  iccfth = 24
  imodif = 0
  call cfther                                                     &
  !==========
 ( nvar   , nscal  ,                                              &
   iccfth , imodif ,                                              &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
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

! Free memory
deallocate(viscf, viscb)
deallocate(smbrs, rovsdt)
deallocate(wflmas, wflmab)
deallocate(w1)
deallocate(w7, w8, w9)
deallocate(w10)
deallocate(dpvar)

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
