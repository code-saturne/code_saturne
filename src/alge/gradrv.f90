!-------------------------------------------------------------------------------

!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2011 EDF S.A., France

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

subroutine gradrv &
!================
 ( ncelet , ncel   , nfac   , nfabor ,                            &
   imrgra , inc    , nswrgp ,                                     &
   iwarnp , nfecra , epsrgp , extrap ,                            &
   ifacel , ifabor , ivar   ,                                     &
   volume , surfac , surfbo , pond   , xyzcen , cdgfac , cdgfbo , &
   dijpf  , diipb  , dofij  ,                                     &
   coefav , coefbv , pvar   ,                                     &
   cocg   ,                                                       &
   gradv  )

!===============================================================================
! FONCTION :
! ----------

! COMPUTE THE GRADIENT OF A VECTOR WITH AN ITERATIVE TECHNIC IN ORDER TO HANDLE
! WITH NONE ORTHOGANALITIES (NSWRGP>1). WE DO NOT TAKE INTO ACCOUNT ANY VOLUMIC
! FORCE HERE.

! COMPUTE COCG AT THE FISRT CALL AND IF NEEDED.

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! ncelet           ! i  ! <-- ! number of extended (real + ghost) cells        !
! ncel             ! i  ! <-- ! number of cells                                !
! nfac             ! i  ! <-- ! number of interior faces                       !
! nfabor           ! i  ! <-- ! number of boundary faces                       !
! ivar             ! e  ! <-- ! numero de la variable                          !
!                  !    !     !   destine a etre utilise pour la               !
!                  !    !     !   periodicite uniquement (pering)              !
!                  !    !     !   on pourra donner ivar=0 si la                !
!                  !    !     !   variable n'est ni une composante de          !
!                  !    !     !   la vitesse, ni une composante du             !
!                  !    !     !   tenseur des contraintes rij                  !
! imrgra           ! e  ! <-- ! methode de reconstruction du gradient          !
!                  !    !     !  0 reconstruction 97                           !
!                  !    !     !  1 moindres carres                             !
!                  !    !     !  2 moindres carres support etendu              !
!                  !    !     !    complet                                     !
!                  !    !     !  3 moindres carres avec selection du           !
!                  !    !     !    support etendu                              !
! inc              ! e  ! <-- ! indicateur = 0 resol sur increment             !
!                  !    !     !              1 sinon                           !
! nswrgp           ! e  ! <-- ! nombre de sweep pour reconstruction            !
!                  !    !     !             des gradients                      !
! imligp           ! e  ! <-- ! methode de limitation du gradient              !
!                  !    !     !  < 0 pas de limitation                         !
!                  !    !     !  = 0 a partir des gradients voisins            !
!                  !    !     !  = 1 a partir du gradient moyen                !
! iwarnp           ! i  ! <-- ! verbosity                                      !
! nfecra           ! e  ! <-- ! unite du fichier sortie std                    !
! epsrgp           ! r  ! <-- ! precision relative pour la                     !
!                  !    !     !  reconstruction des gradients 97               !
! climgp           ! r  ! <-- ! coef gradient*distance/ecart                   !
! extrap           ! r  ! <-- ! coef extrap gradient                           !
! volume(ncelet)   ! tr ! <-- ! volume des elements                            !
! surfac(3,nfac)   ! tr ! <-- ! surf vectorielle des surfaces interne          !
! surfbo           ! tr ! <-- ! surf vectorielle des surfaces de bord          !
!   (3,nfabor)     !    !     !                                                !
! pond(nfac)       ! tr ! <-- ! ponderation geometrique (entre 0 et 1          !
! xyzcen           ! tr ! <-- ! point associes aux volumes de control          !
! (3,ncelet)       !    !     !                                                !
! cdgfac           ! tr ! <-- ! centre de gravite des faces internes           !
! (3,nfac)         !    !     !                                                !
! cdgfbo           ! tr ! <-- ! centre de gravite des faces de bord            !
! (3,nfabor)       !    !     !                                                !
! dijpf(3,nfac)    ! tr ! <-- ! vect i'j', i' (resp. j') projection            !
!                  !    !     !  du centre i (resp. j) sur la normale          !
!                  !    !     !  a la face interne                             !
! diipb            ! tr ! <-- ! vect ii', ii projection du centre i            !
!   (3,nfabor)     !    !     !  sur la normale a la face de bord              !
! cocg(ncelet,3    ! tr ! <-- ! couplage des composantes du gradient           !
!     ,3)          !    !     ! lors de la reconstruction                      !
! dofij(3,nfac)    ! tr ! --> ! vecteur of pour les faces internes             !
!                  !    !     ! o : intersection de ij et la face              !
!                  !    !     ! f : centre de la face                          !
! pond(nfac)       ! tr ! <-- ! ponderation geometrique (entre 0 et 1)         !
! pvar(3,ncelet)    ! tr ! <-- ! variable (vitesse)                             !
! coefav,coefbv    ! tr ! <-- ! tableaux des cond lim pour pvar                 !
!   (3,nfabor)     !    !     !  sur la normale a la face de bord              !
! gradv            ! tr ! --> ! gradient de vitesse                            !
!   (3,3,ncelet)   !    !     !                                                !
! gradva           ! tr ! --- ! tableau de travail pour le grad de vitesse     !
!   (3,3,ncelet)   !    !     !                                                !
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
use cstphy
use cstnum
use albase
use cplsat
use parall
use period

!===============================================================================

implicit none

! Arguments

integer          ncelet , ncel   , nfac   , nfabor
integer          ivar   , imrgra , inc    , nswrgp
integer          imligp , iwarnp , nfecra
double precision epsrgp , climgp , extrap

integer          ifacel(2,nfac), ifabor(nfabor)
double precision volume(ncelet), surfac(3,nfac)
double precision surfbo(3,nfabor)
double precision pond(nfac)
double precision xyzcen(3,ncelet)
double precision cdgfac(3,nfac),cdgfbo(3,nfabor)
double precision dijpf(3,nfac), diipb(3,nfabor)
double precision cocg(3,3,ncelet)
double precision dofij(3,nfac)

double precision pvar(3,ncelet), coefav(3,nfabor), coefbv(3,3,nfabor)
double precision gradv (3,3,ncelet)

! Local variables

integer          imlini

double precision climin

integer          lbloc
parameter       (lbloc = 1024)
integer          iel, ifac, ii, jj, kk, ll, mm
integer          isou, jsou
integer          isqrt, isweep, nswmax
integer          ibloc,nbloc,irel
double precision pfac,pip,deltpx,deltpy,deltpz
double precision rnorx,rnory,rnorz,rnorm,residu
double precision pfaci
double precision dof,fmoyx,fmoyy,fmoyz
double precision aa(3,3,lbloc),unsdet,pfsx,pfsy,pfsz
double precision a11,a12,a13,a21,a22,a23,a31,a32,a33
double precision cocg11,cocg12,cocg13,cocg21,cocg22,cocg23
double precision cocg31,cocg32,cocg33
double precision unsvol, vecfac

double precision, dimension(:,:,:), allocatable :: gradva

integer          ipassv
data             ipassv /0/
save             ipassv

!===============================================================================

!===============================================================================
! 0. PREPARATION POUR PERIODICITE DE ROTATION
!===============================================================================

! Par defaut, on traitera le gradient comme un tenseur ...
!   (i.e. on suppose que c'est le gradient d'une grandeurs vectorielle)

if(iperio.eq.1) then
  call synvin(pvar)
  !==========
endif


!===============================================================================
! 1. INITIALISATION
!===============================================================================

isqrt = 1
nswmax = nswrgp
isweep = 1

! Allocate temporary array
allocate(gradva(3,3,ncelet))

!===============================================================================
! 2. CALCUL SANS RECONSTRUCTION
!===============================================================================

!     SI INITIALISATION PAR MOINDRES CARRES (IMRGRA=4), B EST DEJA REMPLI
!     SINON (IMRGRA=0) ON CALCULE UN GRADIENT SANS RECONSTRUCTION

if (imrgra.eq.0) then


  do iel = 1, ncelet
  ! BOUCLE SUR LES COMPOSANTES X, Y, Z
    do jsou = 1, 3
    ! BOUCLE SUR LES COMPOSANTES U, V, W
      do isou = 1 ,3
        gradva(isou,jsou,iel) = 0.d0
      enddo
    enddo
  enddo

  !     ASSEMBLAGE A PARTIR DES FACETTES FLUIDES

  do ifac = 1,nfac
  ! BOUCLE SUR LES COMPOSANTES U, V, W
    do isou = 1 ,3
      ii = ifacel(1,ifac)
      jj = ifacel(2,ifac)
      pfac   = pond(ifac)*pvar(isou,ii) +(1.d0-pond(ifac))*pvar(isou,jj)
      do jsou = 1, 3
        gradva(isou,jsou,ii) = gradva(isou,jsou,ii) + pfac*surfac(jsou,ifac)
        gradva(isou,jsou,ii) = gradva(isou,jsou,jj) - pfac*surfac(jsou,ifac)
      enddo
    enddo
  enddo

  !     ASSEMBLAGE A PARTIR DES FACETTES DE BORD

  do ifac = 1,nfabor
  ! BOUCLE SUR LES COMPOSANTES U, V, W
    do isou = 1 ,3
      ii = ifabor(ifac)
      pfac = inc*coefav(isou,ifac) + coefbv(isou,1,ifac)*pvar(1,ii)    &
                                   + coefbv(isou,2,ifac)*pvar(2,ii)    &
                                   + coefbv(isou,3,ifac)*pvar(3,ii)
      do jsou = 1, 3
        gradva(isou,jsou,ii) = gradva(isou,jsou,ii) +pfac*surfbo(jsou,ifac)
      enddo
    enddo
  enddo


  ! GRAVEL = GRADIENT

  do iel = 1, ncel
    do jsou = 1, 3
      do isou = 1, 3
        unsvol = 1.d0/volume(iel)
        gradv(isou,jsou,iel) = gradva(isou,jsou,iel)*unsvol
      enddo
    enddo
  enddo

!     TRAITEMENT DU PARALLELISME ET DE LA PERIODICITE

  if(irangp.ge.0.or.iperio.eq.1) then
    call syntin (gradv)
    !==========
  endif

endif

if( nswrgp.le.1 ) then
  ! Free memory
  deallocate(gradva)

  return
endif


!     On incremente IPASS quand on calcule COCG pour la premiere fois
ipassv = ipassv + 1

!===============================================================================
! 3. RECONSTRUCTION DES GRADIENTS POUR LES MAILLAGES TORDUS
!===============================================================================

!     RESOLUTION SEMI-IMPLICITE SUR TOUT LE MAILLAGE
!     gradv = GRADIENT


if(ipassv.eq.1 .or. iale.eq.1 .or. imobil.eq.1) then

! ---> CALCUL DE COCG

  do iel =1,ncelet
    do jj = 1, 3
      do ii = 1, 3
        cocg(ii,jj,iel)= 0.d0
      enddo
    enddo
  enddo
  do iel=1,ncel
    cocg(1,1,iel) = volume(iel)
    cocg(2,2,iel) = volume(iel)
    cocg(3,3,iel) = volume(iel)
  enddo

! ---> AJOUT DES CONTRIBUTIONS DES FACES INTERNES
! ATTENTION, LA MATRICE 3X3 COCG(iel) EST ICI LA TRANSPOSEE
!            DE CELLE UTILISEE DANS LE GRADIENT SCALAIRE ET NE PREND PAS EN
!            COMPTE LES TERMES DE BORD

  do ifac=1,nfac
    ii = ifacel(1,ifac)
    jj = ifacel(2,ifac)

    do isou =1,3
      !---> DOF = OF
      dof = dofij(isou,ifac)

      pfaci  = -dof*0.5d0
      do jsou =1,3
        vecfac = pfaci*surfac(jsou,ifac)
        cocg(isou,jsou,ii) = cocg(isou,jsou,ii) + vecfac
        cocg(isou,jsou,jj) = cocg(isou,jsou,jj) - vecfac
      enddo
    enddo

  enddo


! prise en compte explicite des faces de bords: pas de recalcul de COCG

! ---> ON INVERSE POUR TOUTE LES CELLULES : LE COCG POUR TOUTES LES CELLULES
!      RESTE ENSUITE LE MEME TANT QUE LE MAILLAGE NE CHANGE PAS
  nbloc = ncel/lbloc
  if (nbloc.gt.0) then
    do ibloc = 1, nbloc
      do ii =1, lbloc
        iel = (ibloc-1)*lbloc+ii

        cocg11 = cocg(1,1,iel)
        cocg12 = cocg(1,2,iel)
        cocg13 = cocg(1,3,iel)
        cocg21 = cocg(2,1,iel)
        cocg22 = cocg(2,2,iel)
        cocg23 = cocg(2,3,iel)
        cocg31 = cocg(3,1,iel)
        cocg32 = cocg(3,2,iel)
        cocg33 = cocg(3,3,iel)

        a11=cocg22*cocg33-cocg32*cocg23
        a12=cocg32*cocg13-cocg12*cocg33
        a13=cocg12*cocg23-cocg22*cocg13
        a21=cocg31*cocg23-cocg21*cocg33
        a22=cocg11*cocg33-cocg31*cocg13
        a23=cocg21*cocg13-cocg11*cocg23
        a31=cocg21*cocg32-cocg31*cocg22
        a32=cocg31*cocg12-cocg11*cocg32
        a33=cocg11*cocg22-cocg21*cocg12

        unsdet = 1.d0/(cocg11*a11 +cocg21*a12+cocg31*a13)

        aa(1,1,ii) = a11*unsdet
        aa(1,2,ii) = a12*unsdet
        aa(1,3,ii) = a13*unsdet
        aa(2,1,ii) = a21*unsdet
        aa(2,2,ii) = a22*unsdet
        aa(2,3,ii) = a23*unsdet
        aa(3,1,ii) = a31*unsdet
        aa(3,2,ii) = a32*unsdet
        aa(3,3,ii) = a33*unsdet

      enddo

      do ii = 1, lbloc
        iel = (ibloc-1)*lbloc+ii
        cocg(1,1,iel) = aa(1,1,ii)
        cocg(1,2,iel) = aa(1,2,ii)
        cocg(1,3,iel) = aa(1,3,ii)
        cocg(2,1,iel) = aa(2,1,ii)
        cocg(2,2,iel) = aa(2,2,ii)
        cocg(2,3,iel) = aa(2,3,ii)
        cocg(3,1,iel) = aa(3,1,ii)
        cocg(3,2,iel) = aa(3,2,ii)
        cocg(3,3,iel) = aa(3,3,ii)
      enddo

    enddo

  endif

  irel = mod(ncel,lbloc)
  if (irel.gt.0) then
    ibloc = nbloc + 1
    do ii = 1, irel
      iel = (ibloc-1)*lbloc+ii

      cocg11 = cocg(1,1,iel)
      cocg12 = cocg(1,2,iel)
      cocg13 = cocg(1,3,iel)
      cocg21 = cocg(2,1,iel)
      cocg22 = cocg(2,2,iel)
      cocg23 = cocg(2,3,iel)
      cocg31 = cocg(3,1,iel)
      cocg32 = cocg(3,2,iel)
      cocg33 = cocg(3,3,iel)

      a11=cocg22*cocg33-cocg32*cocg23
      a12=cocg32*cocg13-cocg12*cocg33
      a13=cocg12*cocg23-cocg22*cocg13
      a21=cocg31*cocg23-cocg21*cocg33
      a22=cocg11*cocg33-cocg31*cocg13
      a23=cocg21*cocg13-cocg11*cocg23
      a31=cocg21*cocg32-cocg31*cocg22
      a32=cocg31*cocg12-cocg11*cocg32
      a33=cocg11*cocg22-cocg21*cocg12

      unsdet = 1.d0/(cocg11*a11 +cocg21*a12+cocg31*a13)

      aa(1,1,ii) = a11*unsdet
      aa(1,2,ii) = a12*unsdet
      aa(1,3,ii) = a13*unsdet
      aa(2,1,ii) = a21*unsdet
      aa(2,2,ii) = a22*unsdet
      aa(2,3,ii) = a23*unsdet
      aa(3,1,ii) = a31*unsdet
      aa(3,2,ii) = a32*unsdet
      aa(3,3,ii) = a33*unsdet

    enddo

    do ii = 1, irel
      iel = (ibloc-1)*lbloc+ii
      cocg(1,1,iel) = aa(1,1,ii)
      cocg(1,2,iel) = aa(1,2,ii)
      cocg(1,3,iel) = aa(1,3,ii)
      cocg(2,1,iel) = aa(2,1,ii)
      cocg(2,2,iel) = aa(2,2,ii)
      cocg(2,3,iel) = aa(2,3,ii)
      cocg(3,1,iel) = aa(3,1,ii)
      cocg(3,2,iel) = aa(3,2,ii)
      cocg(3,3,iel) = aa(3,3,ii)
    enddo

  endif

endif

! ---> CALCUL DU RESIDU DE NORMALISATION

! TODO a opimized form prodcsc9
! Euclidian norm
call prodsc(9*ncelet,9*ncel,isqrt,gradva(1,1,1),gradva(1,1,1),       &
                              rnorm )


if (volmax.gt.1.d0) rnorm = rnorm / volmax
if( rnorm.le.epzero ) then
  ! Free memory
  deallocate(gradva)
  return
endif

!  LE VECTEUR OijFij EST CALCULE DANS CLDIJP

! ---> DEBUT DES ITERATIONS

 100  continue

isweep = isweep +1


!     CALCUL DU SECOND MEMBRE


do iel = 1, ncel
  do jsou =1, 3
    do isou = 1, 3
      gradva(isou,jsou,iel) = -gradv(isou,jsou,iel)*volume(iel)
    enddo
  enddo
enddo


!     ASSEMBLAGE A PARTIR DES FACETTES FLUIDES

do ifac = 1,nfac
  ii = ifacel(1,ifac)
  jj = ifacel(2,ifac)
  do isou = 1, 3
    vecfac = pond(ifac)*pvar(isou,ii) +(1.d0-pond(ifac))*pvar(isou,jj)         &
           +0.5d0*(  (gradv(isou,1,ii)+gradv(isou,1,jj))*dofij(1,ifac)   &
                    +(gradv(isou,2,ii)+gradv(isou,2,jj))*dofij(2,ifac)   &
                    +(gradv(isou,3,ii)+gradv(isou,3,jj))*dofij(3,ifac))
    do jsou = 1, 3
      gradva(isou,jsou,ii) = gradva(isou,jsou,ii)             &
                             + vecfac*surfac(jsou,ifac)
      gradva(isou,jsou,jj) = gradva(isou,jsou,jj)             &
                             - vecfac*surfac(jsou,ifac)
    enddo
  enddo
enddo


!     ASSEMBLAGE A PARTIR DES FACETTES DE BORD

do ifac = 1,nfabor
  ii = ifabor(ifac)

  do isou = 1,3
    pfac = inc*coefav(isou,ifac)

    ! coefbv is a matrix
    do jsou = 1, 3
      pip =  pvar(jsou,ii)                                   &
         + gradv(jsou,1,ii)*diipb(1,ifac)                 &
         + gradv(jsou,2,ii)*diipb(2,ifac)                 &
         + gradv(jsou,3,ii)*diipb(3,ifac)
      pfac = pfac + coefbv(isou,jsou,ifac)*pip
    enddo

    do jsou = 1, 3
      gradva(isou,jsou,ii) = gradva(isou,jsou,ii) +pfac*surfbo(jsou,ifac)
    enddo
  enddo
enddo


!     INCREMENTATION DU GRADIENT

do iel =1,ncel

  do jsou = 1, 3
    do isou = 1, 3

      gradv(isou,jsou,iel) = gradv(isou,jsou,iel)             &
                        + gradva(isou,1,iel) *cocg(1,jsou,iel)  &
                        + gradva(isou,2,iel) *cocg(2,jsou,iel)  &
                        + gradva(isou,3,iel) *cocg(3,jsou,iel)

    enddo
  enddo

enddo


!     TRAITEMENT DU PARALLELISME

if(irangp.ge.0.or.iperio.eq.1) then
  call syntin (gradv)
endif

! ---> TEST DE CONVERGENCE

! Norme Euclidienne
call prodsc(9*ncelet,9*ncel,isqrt,gradva(1,1,1),gradva(1,1,1),       &
                              residu         )

if (volmax.gt.1.d0) residu = residu / volmax

if( residu.le.epsrgp*rnorm) then
  if( iwarnp.ge.2 ) then
    write (nfecra,1000) isweep,residu/rnorm,rnorm,ivar
  endif
  goto 101
elseif( isweep.ge.nswmax ) then
  if( iwarnp.ge.0) then
     write (nfecra,1000)isweep,residu/rnorm,rnorm,ivar
     write (nfecra,1100)
  endif
  goto 101
else
  goto 100
endif

 101  continue


! Free memory
deallocate(gradva)

!--------
! FORMATS
!--------
#if defined(_CS_LANG_FR)

 1000 format(1X,'GRDVEC ISWEEP = ',I4,' RESIDU NORME: ',E11.4,          &
         ' NORME: ',E11.4,/,1X,'PARAMETRE IVAR = ',I4 )
 1100 format(                                                           &
'@                                                            ',/,&
'@ @@ ATTENTION :         NON CONVERGENCE DE GRDVEC           ',/,&
'@    =========                                               ',/,&
'@                                                            '  )

#else

 1000 format(1X,'GRDVEC ISWEEP = ',I4,' NORMED RESIDUAL: ',E11.4,       &
         ' NORM: ',E11.4,/,1X,'PARAMETER IVAR = ',I4 )
 1100 format(                                                           &
'@'                                                            ,/,&
'@ @@ WARNING:            NON CONVERGENCE OF GRDVEC'           ,/,&
'@    ========'                                                ,/,&
'@'                                                              )

#endif

!----
! FIN
!----

return

end subroutine
