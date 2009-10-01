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

subroutine vorini &
!================

 ( ncevor , nvor   , ient   ,                                     &
   ivorce , xyz    , yzcel  , xu     , xv     , xw     ,          &
   yzvor  , signv  , temps  , tpslim )

!===============================================================================
!  FONCTION  :
!  ----------

! INITIALISATION DE LA METHODE DES VORTEX POUR LES ENTREES EN L.E.S.

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
!    nom           !type!mode !                   role                         !
!__________________!____!_____!________________________________________________!
! ncevor           ! e  ! <-- ! nombre de face a l'entree ou est               !
!                  !    !     ! utilise la methode                             !
! nvor             ! e  ! <-- ! nombre de vortex a l'entree                    !
! ient             ! e  ! <-- ! numero de l'entree                             !
! ivorce           ! te ! <-- ! numero du vortex le plus proche d'une          !
!     (nvomax)     !    !     ! face donnee                                    !
! xyz(icvmax,3)    !    ! <-- ! coordonnees des faces d'entree dans            !
!                  !    !     ! le calcul                                      !
! yzcel            ! tr ! <-- ! coordonnees des faces d'entree dans            !
!   (icvmax ,2)    !    !     ! le referentiel local                           !
! xu(icvmax)       ! tr ! --- ! composante de vitesse principale               !
! xv(icvmax)       ! tr ! <-- ! composantes de vitesse transverses             !
! xw(icvmax)       ! tr ! <-- !                                                !
! yzvor            ! tr ! <-- ! coordonnees du centre des vortex               !
!   (nvomax,2)     !    !     !                                                !
! signv(nvomax)    ! tr ! <-- ! sens de rotation des vortex                    !
! temps            ! tr ! <-- ! temps ecoule depuis la creation                !
!     (nvomax)     !    !     ! du vortex                                      !
! tpslim           ! tr ! <-- ! duree de vie du vortex                         !
!     (nvomax)     !    !     !                                                !
!__________________.____._____.________________________________________________.

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
include "cstphy.h"
include "cstnum.h"
include "optcal.h"
include "entsor.h"
include "vortex.h"

!===============================================================================

! Arguments

integer          ncevor , nvor   , ient
integer          ivorce(nvomax)

double precision xyz(icvmax,3)   , yzcel(icvmax,2)
double precision xu(icvmax)      , xv(icvmax)      , xw(icvmax)
double precision yzvor(nvomax,2) , signv(nvomax)
double precision temps(nvomax)   , tpslim(nvomax)


! VARIABLES LOCALES

integer          ii, jj, kk, iii, iun, ivort, iient, iok

double precision dd
double precision drand(1), phidat, xx, yy, zz
double precision uu, vv, ww
double precision u_vor, ek_vor, ee_vor

integer          ilect
data             ilect /0/
save             ilect
!===============================================================================
! 1. CALCUL DU REPERE LOCAL ET CHANGEMENT DE REPERE
!===============================================================================

ilect = ilect + 1
if(ilect.eq.1) then
  write(nfecra,1000) nnent, isuivo
endif

if(icas(ient).eq.1.or.icas(ient).eq.2.or.icas(ient).eq.3) then

  dir3(1,ient)=dir1(2,ient)*dir2(3,ient)-dir1(3,ient)*dir2(2,ient)
  dir3(2,ient)=dir1(3,ient)*dir2(1,ient)-dir1(1,ient)*dir2(3,ient)
  dir3(3,ient)=dir1(1,ient)*dir2(2,ient)-dir1(2,ient)*dir2(1,ient)

elseif(icas(ient).eq.4) then

!     On s'aide du vecteur surface d'une face de l'entree (supposee plane)
!     pour definir le repere local

  vv = sqrt(surf(1,ient)**2 + surf(2,ient)**2 + surf(3,ient)**2)

  dir3(1,ient) = -surf(1,ient)/vv
  dir3(2,ient) = -surf(2,ient)/vv
  dir3(3,ient) = -surf(3,ient)/vv

!     On se fixe, par exemple, x1 = 0 et y1 = 1 et on norme le vecteur

  dir1(1,ient) = 0.d0
  dir1(2,ient) = 1.d0
  if (abs(dir3(3,ient)).gt.epzero) then
    dir1(3,ient) = -dir3(2,ient)/dir3(3,ient)
  else
    dir1(3,ient) = 0.d0
  endif

  vv = sqrt(dir1(1,ient)**2 + dir1(2,ient)**2 + dir1(3,ient)**2)

  dir1(1,ient) = dir1(1,ient)/vv
  dir1(2,ient) = dir1(2,ient)/vv
  dir1(3,ient) = dir1(3,ient)/vv

!     On obtient le dernier vecteur par produit vectoriel des deux autres

  dir2(1,ient) =                                                  &
    dir3(2,ient)*dir1(3,ient) - dir3(3,ient)*dir1(2,ient)
  dir2(2,ient) =                                                  &
    dir3(3,ient)*dir1(1,ient) - dir3(1,ient)*dir1(3,ient)
  dir2(3,ient) =                                                  &
    dir3(1,ient)*dir1(2,ient) - dir3(2,ient)*dir1(1,ient)

endif

! - Changement de repere (on suppose que les vecteurs sont normes)

do ii = 1, ncevor
  xx = xyz(ii,1) - cen(1,ient)
  yy = xyz(ii,2) - cen(2,ient)
  zz = xyz(ii,3) - cen(3,ient)
  yzcel(ii,1) = dir1(1,ient)*xx+dir1(2,ient)*yy+dir1(3,ient)*zz
  yzcel(ii,2) = dir2(1,ient)*xx+dir2(2,ient)*yy+dir2(3,ient)*zz
enddo

! - Dimensions min et max de l entree

ymax(ient) = -grand
ymin(ient) =  grand
zmax(ient) = -grand
zmin(ient) =  grand
do ii = 1, ncevor
  ymax(ient) = max(ymax(ient),yzcel(ii,1))
  ymin(ient) = min(ymin(ient),yzcel(ii,1))
  zmax(ient) = max(zmax(ient),yzcel(ii,2))
  zmin(ient) = min(zmin(ient),yzcel(ii,2))
enddo

! - Verification

if(icas(ient).eq.1) then
  if(lly(ient).lt.ymax(ient)-ymin(ient).or.                       &
       llz(ient).lt.zmax(ient)-zmin(ient)) then
    write(nfecra,2000) ient
    call csexit(1)
  endif
elseif(icas(ient).eq.2) then
  if(lld(ient).lt.ymax(ient)-ymin(ient).or.                       &
       lld(ient).lt.zmax(ient)-zmin(ient)) then
    write(nfecra,2000) ient
    call csexit(1)
  endif
endif
!===============================================================================
! 2. IMPRESSIONS DES PARAMETES A CHAQUE ENTREE
!===============================================================================

call vorimp(ient)
!==========
!===============================================================================
! 3. REMPLISSAGE DES TABLEAUX DE DONNEES EN ENTREE
!===============================================================================

iun = 1

if(icas(ient).eq.1.or.icas(ient).eq.2.or.icas(ient).eq.3) then
  open(file=ficvor(ient),unit=impdvo)
  rewind(impdvo)
  do ii = 1, ndat(ient)
    read(impdvo,*)                                                &
         xdat(ii,ient) , ydat(ii,ient) , zdat(ii,ient)  ,         &
         udat(ii,ient) , vdat(ii,ient) , wdat(ii,ient)  ,         &
         dudat(ii,ient), kdat(ii,ient) , epsdat(ii,ient)
  enddo
  close(impdvo)
  write(nfecra,3000)
elseif(icas(ient).eq.4) then
  xdat(1,ient)   = cen(1,ient)
  ydat(1,ient)   = cen(2,ient)
  zdat(1,ient)   = cen(3,ient)
  udat(1,ient)   = udebit(ient)
  vdat(1,ient)   = 0.d0
  wdat(1,ient)   = 0.d0
  dudat(1,ient)  = 0.d0
  kdat(1,ient)   = kdebit(ient)
  epsdat(1,ient) = edebit(ient)
endif

! On suppose que les donnees sont fournies
! dans le repere du calcul (et non dans le repere local).
! C'est plus simple pour l'utilisateur, et plus
! naturel pour faire du couplage

do ii = 1, ndat(ient)
  xx = xdat(ii,ient) - cen(1,ient)
  yy = ydat(ii,ient) - cen(2,ient)
  zz = zdat(ii,ient) - cen(3,ient)
  uu = udat(ii,ient)
  vv = vdat(ii,ient)
  ww = wdat(ii,ient)
  ydat(ii,ient) = dir1(1,ient)*xx+dir1(2,ient)*yy+dir1(3,ient)*zz
  zdat(ii,ient) = dir2(1,ient)*xx+dir2(2,ient)*yy+dir2(3,ient)*zz
  udat(ii,ient) = dir3(1,ient)*uu+dir3(2,ient)*vv+dir3(3,ient)*ww
  vdat(ii,ient) = dir1(1,ient)*uu+dir1(2,ient)*vv+dir1(3,ient)*ww
  wdat(ii,ient) = dir2(1,ient)*uu+dir2(2,ient)*vv+dir2(3,ient)*ww
enddo

! --- Verfication des donnees

iok = 0
do ii = 1, ndat(ient)
  if(udat(ii,ient).le.0.d0.or.kdat(ii,ient).le.0.d0.or.           &
       epsdat(ii,ient).le.0.d0) then
    write(nfecra,3100) ient
    call csexit (1)
  endif
  if(icas(ient).eq.1) then
    if(ydat(ii,ient).lt.-lly(ient)/2.d0.or.                       &
         ydat(ii,ient).gt.lly(ient)/2.d0.or.                      &
         zdat(ii,ient).lt.-llz(ient)/2.d0.or.                     &
         zdat(ii,ient).gt.llz(ient)/2.d0) then
      iok = iok + 1
    endif
  elseif(icas(ient).eq.2) then
    if(ydat(ii,ient).lt.-lld(ient)/2.d0.or.                       &
         ydat(ii,ient).gt.lld(ient)/2.d0.or.                      &
         zdat(ii,ient).lt.-lld(ient)/2.d0.or.                     &
         zdat(ii,ient).gt.lld(ient)/2.d0) then
      iok = iok + 1
    endif
  endif
enddo

if(iok.gt.0) then
  write(nfecra,3200) ient
endif

!===============================================================================
! 4. LECTURE DU FICHIER SUITE / INITIALISATION DU CHAMP DE VORTEX
!===============================================================================

if(isuivo.eq.1) then

  if(ient.eq.1) then
    open(unit=impmvo,file=ficmvo)
    rewind(impmvo)
  endif

  read(impmvo,100) iient
  read(impmvo,100) ivort
  if(ivort.ne.nvor.or.iient.ne.ient) then
    write(nfecra,4500) ient, ivort, nvor
    initvo(ient) = 1
  else
    do ii = 1, nvor
      read(impmvo,200) yzvor(ii,1), yzvor(ii,2),                  &
             temps(ii), tpslim(ii), signv(ii)
    enddo
    initvo(ient) = 0
    write(nfecra,4000)
  endif

  if(ient.eq.nnent) then
    close(impmvo)
  endif

endif

if(isuivo.eq.0.or.initvo(ient).eq.1) then

!-------------------------------
!  Tirage des positions
!-------------------------------
  iun = 1
  if(icas(ient).eq.1)then
    do ii = 1, nvor
      call zufall(iun,drand(1))
      yzvor(ii,1) = lly(ient) * drand(1) - lly(ient)/2.d0
      call zufall(iun,drand(1))
      yzvor(ii,2) = llz(ient) * drand(1) - llz(ient)/2.d0
    enddo
  elseif(icas(ient).eq.2) then
    do ii = 1, nvor
 15         continue
      call zufall(iun,drand(1))
      yzvor(ii,1) = lld(ient) * drand(1) - lld(ient)/2.0d0
      call zufall(iun,drand(1))
      yzvor(ii,2) = lld(ient) * drand(1) - lld(ient)/2.0d0
      if ((yzvor(ii,1)**2+yzvor(ii,2)**2).gt.                     &
           (lld(ient)/2.d0)**2) then
        goto 15
      endif
    enddo
  elseif(icas(ient).eq.3.or.icas(ient).eq.4) then
    do ii = 1, nvor
      call zufall(iun,drand(1))
      yzvor(ii,1) = ymin(ient) + lly(ient) * drand(1)
      call zufall(iun,drand(1))
      yzvor(ii,2) = zmin(ient) + llz(ient) * drand(1)
    enddo
  endif
!--------------
! Duree de vie
!--------------
  if(itlivo(ient).eq.1) then
    do ii = 1, nvor
      call zufall(iun,drand(1))
      temps(ii) = drand(1)*tlimvo(ient)

! on fait cela pour que les vortex ne disparaissent pas tous
! en meme temps    .

      tpslim(ii) = tlimvo(ient)
    enddo
  elseif(itlivo(ient).eq.2) then
    do ii = 1, nvor
      yy = yzvor(ii,1)
      zz = yzvor(ii,2)
      iii = 0
      u_vor  =  phidat(nfecra,icas(ient),ndat(ient),yy,zz,        &
                ydat(1,ient),zdat(1,ient),udat(1,ient),iii)
      ek_vor =  phidat(nfecra,icas(ient),ndat(ient),yy,zz,        &
                ydat(1,ient),zdat(1,ient),kdat(1,ient),iii)
      ee_vor =  phidat(nfecra,icas(ient),ndat(ient),yy,zz,        &
                ydat(1,ient),zdat(1,ient),epsdat(1,ient),iii)
      tpslim(ii) = 5.d0*cmu*ek_vor**(3.d0/2.d0)/ee_vor
      tpslim(ii) = tpslim(ii)/u_vor
      temps(ii) = 0.d0
    enddo
  endif
!------------------
! Sens de rotation
!------------------
  do ii = 1, nvor
    signv(ii) = 1.d0
    call zufall(iun,drand(1))
    if (drand(1).lt.0.5d0) signv(ii) = -1.d0
  enddo
endif

!===============================================================================
! 5. AJOUT DE LA VITESSE MOYENNE POUR U
!===============================================================================

do ii = 1, ncevor
  yy = yzcel(ii,1)
  zz = yzcel(ii,2)
  iii = 0
  u_vor = phidat(nfecra,icas(ient),ndat(ient),yy,zz,              &
          ydat(1,ient),zdat(1,ient),udat(1,ient),iii)
  xu(ii)= u_vor
enddo

!===============================================================================
! 6. RECHERCHE DE LA FACE LA PLUS PROCHE DE CHAQUE VORTEX
!===============================================================================
do ii = 1, nvor
  kk = 0
  dd = grand
  do jj = 1, ncevor
    if(((yzcel(jj,1)-yzvor(ii,1))**2+                             &
        (yzcel(jj,2)-yzvor(ii,2))**2).lt.dd)then
      dd = (yzcel(jj,1)-yzvor(ii,1))**2                           &
          +(yzcel(jj,2)-yzvor(ii,2))**2
      kk = jj
    endif
  enddo
  ivorce(ii) = kk
enddo

! FORMATS

#if defined(_CS_LANG_FR)

 1000 format(                                                           &
'                                                             ',/,&
' ** METHODE DES VORTEX                                       ',/,&
'    ------------------                                       ',/,&
'       NNENT  = ',4X,I10,    ' (Nombre d entrees            )',/,&
'       ISUIVO = ',4X,I10,    ' (1 : suite de calcul         )'  )

 2000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    LES DIMENSIONS MAX DE L''ENTREE SONT INCOMPATIBLES AVEC ',/,&
'@    LES DONNEES A L''ENTREE ',I10                            ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier usvort.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 3000 format(                                                           &
'                                                             ',/,&
' --  Fin de la lecture du fichier de donnees                 ',/)
 3100 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    U, K ET EPSILON SONT DES GRANDEURS QUI SONT DEFINIES    ',/,&
'@    POSITIVES DANS LE REPERE LOCAL DE L''ENTREE             ',/,&
'@                                                            ',/,&
'@    VERIFIER LE FICHIER DE DONNEE DE L''ENTREE ',I10         ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier usvort.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 3200 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION :       A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    LES DIMENSIONS MAX DE L''ENTREE SONT INCOMPATIBLES AVEC ',/,&
'@    CELLES DU FICHIER DE DONNEES                            ',/,&
'@                                                            ',/,&
'@    VERIFIER LE FICHIER DE DONNEE DE L''ENTREE ',I10         ,/,&
'@                                                            ',/,&
'@  Le calcul sera execute.                                   ',/,&
'@                                                            ',/,&
'@  Verifier usvort.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 4000 format(                                                           &
' --  Fin de la lecture du fichier suite                      ',/)

 4500   format(                                                         &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION :      A L''ENTREE DES DONNEES                ',/,&
'@    =========                                               ',/,&
'@  LE NOMBRE DE VORTEX A CHANGE A L''ENTREE',I10              ,/,&
'@    NVORT VALAIT PRECECEDEMENT ',I10                         ,/,&
'@    A L''ENTREE ',I10                                        ,/,&
'@    ET VAUT MAINTENANT  ',I10                                ,/,&
'@                                                            ',/,&
'@  LA METHODE EST REINITIALISE A CETTE ENTREE                ',/,&
'@                                                            ',/,&
'@  Le calcul sera execute                                    ',/,&
'@                                                            ',/,&
'@  Verifier usvort.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 100  format(i10)
 200  format(5e13.5)

#else

 1000 format(                                                           &
'                                                             ',/,&
' ** VORTEX METHOD                                            ',/,&
'    -------------                                            ',/,&
'       NNENT  = ',4X,I10,    ' (Number of inlets            )',/,&
'       ISUIVO = ',4X,I10,    ' (1: calculation restart      )'  )

 2000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: ABORT IN THE DATA SPECIFICATION                ',/,&
'@    ========                                                ',/,&
'@    THE MAX DIMENSIONS OF THE INLET ARE INCOMPATIBLE WITH   ',/,&
'@    THE DATA AT THE INLET ',I10                              ,/,&
'@                                                            ',/,&
'@  The calculation will not be run.                          ',/,&
'@                                                            ',/,&
'@  Verify usvort.                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 3000 format(                                                           &
'                                                             ',/,&
' --  End reading the data file                               ',/)
 3100 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: ABORT IN THE DATA SPECIFICATION                ',/,&
'@    ========                                                ',/,&
'@    U, K AND EPSILON ARE QUANTITIES WHICH MUST BE POSITIVE  ',/,&
'@    IN THE LOCAL FRAME OF THE INLET                         ',/,&
'@                                                            ',/,&
'@    VERIFY THE DATA FILE FOR THE INLET ',I10                 ,/,&
'@                                                            ',/,&
'@  The calculation will not be run.                          ',/,&
'@                                                            ',/,&
'@  Verify usvort.                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 3200 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING:       IN THE DATA SPECIFICATION                ',/,&
'@    ========                                                ',/,&
'@    THE MAX DIMENSIONS OF THE INLET ARE INCOMPATIBLE WITH   ',/,&
'@    THE ONES FROM THE DATA FILE                             ',/,&
'@                                                            ',/,&
'@    VERIFY THE DATA FILE FOR THE INLET ',I10                 ,/,&
'@                                                            ',/,&
'@  The calculation will be run.                              ',/,&
'@                                                            ',/,&
'@  Verify usvort.                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 4000 format(                                                           &
' --  End reading the restart file                            ',/)

 4500   format(                                                         &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING:       IN THE DATA SPECIFICATION                ',/,&
'@    ========                                                ',/,&
'@  THE NUMBER OF VORTICES HAS CHANGED AT THE INLET ',I10      ,/,&
'@    NVORT WAS PREVIOUSLY ',I10                               ,/,&
'@    AT THE INLET ',I10                                       ,/,&
'@    IT IS CURRENTLY ',I10                                    ,/,&
'@                                                            ',/,&
'@  THE METHOD IS RE-INITIALIZED AT THIS INLET                ',/,&
'@                                                            ',/,&
'@  The calculation will be run.                              ',/,&
'@                                                            ',/,&
'@  Verify usvort.                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 100  format(i10)
 200  format(5e13.5)

#endif

return
end subroutine
