!-------------------------------------------------------------------------------

!VERS


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

subroutine uslabo &
!================

 ( idbia0 , idbra0 ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndnod , lndfac , lndfbr , ncelbr ,                   &
   nvar   , nscal  , nphas  ,                                     &
   nbpmax , nvp    , nvp1   , nvep   , nivep  ,                   &
   ntersl , nvlsta , nvisbr ,                                     &
   kface  , nbpt   , isuivi ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   itypfb , itrifb , ifrlag , itepa  , indep  ,                   &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   surfbn , dt     , rtpa   , rtp    , propce , propfa , propfb , &
   coefa  , coefb  ,                                              &
   ettp   , ettpa  , tepa   , parbor , vitpar , vitflu , auxl   , &
   rdevel , rtuser , ra     )

!===============================================================================
! FONCTION :
! ----------

!   SOUS-PROGRAMME DU MODULE LAGRANGIEN :
!   -------------------------------------

!  SOUS-PROGRAMME UTILISATEUR (INTERVENTION NON OBLIGATOIRE)

!  SOUS-PROGRAMME UTILISATEUR DE GESTION DU COMPORTEMENT
!    DES PARTICULES A UNE INTERACTION PARTICULE/FRONTIERE
!    ET D'ENREGISTREMENT DES STATISTIQUES AUX FRONTIERES.

!  L'UTILISATEUR N'A PAS A MODIFIER LE PRESENT SOUS-PROGRAMME DANS
!    LES CONDITIONS D'UTILISATION STANDARD. DANS LE CAS OU IL
!    SOUHAITE AVOIR DES INTERACTIONS NON PREVU, IL DOIT INTERVENIR
!    DANS LES RUBRIQUES NUMERO 8 et 10.

!  Gestion de l'interaction particule/face de frontiere
!    en fonction des informations donnees par l'utilisateur
!    (valeur de IUSCLB par zone) dans la routine USLAG2.

!  En fonction du nom stocke dans IUSCLB et affecte
!    a la face de bord KFACE, on donne le type de comportement
!    des particules.

!  En standard, on a :

!    IUSCLB conditions au bord pour les particules
!    = IENTRL -> zone d'injection de particules
!    = ISORTL -> sortie du domaine
!    = IREBOL -> rebond des particules
!    = IDEPO1 -> deposition definitive
!    = IDEPO2 -> deposition definitive mais la particule reste en
!                memoire (utile si IENSI2 = 1 uniquement)
!    = IDEPO3 -> deposition et remise en suspension possible
!                suivant les condition de l'ecoulement
!    = IENCRL -> encrassement (Charbon uniquement IPHYLA = 2)


!  De plus, si on souhaite un autre type d'interaction non prevue
!    pour zone de faces de bord, dans USLAG2 on affecte dans
!    IUSCLB(KZONE) un des noms suivants :

!               JBORD1, JBORD2, JBORD3, JBORD4, JBORD5

!    et dans cette routine on code le comportement
!    des particules pour cette zone frontiere.

!  ATTENTION : en debut de routine la variable ISUIVI est
!   initialisee a une valeur aberrante et DOIT etre modifiee
!   avant la fin de la routine.
!   Les vitesses de la particule et du fluide vu doivent etre
!   modifiees en fonction de l'interaction via les tableaux
!   VITPAR et VITFLU, elles NE DOIVENT PAS etre modifiees
!   via ETTP et ETTPA dans ce sous-programme.

!  Regle de modification de ISUIVI :
!  =================================
! 1) Mettre ISUIVI = 0 si la particule ne doit pas etre
!    suivi dans le maillage apres l'interaction de sa trajectoire
!    et de la face de bord (ex : IENTRL, ISORTL, IDEPO1, IDEPO2).

! 2) mettre ISUIVI = 1 pour continuer a suivre la particule
!    dans le maillage apres son interaction (ex : IREBOL, IDEPO3).
!    On peut tres bien avoir ISUIVI = 0 ou ISUIVI = 1 pour un type
!    d'interaction en fonction de comportement de la particule
!    (ex : IENCRL).

! REMARQUE : Lorsqu'il y a interaction, les calculs de la vitesse
! ^^^^^^^^   de la particule et de la vitesse du fluide vu sont
!            forcement a l'ordre 1 (meme si on est a l'ordre 2
!            par ailleurs).

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
! lndnod           ! e  ! <-- ! dim. connectivite cellules->faces              !
! lndfac           ! e  ! <-- ! longueur du tableau nodfac                     !
! lndfbr           ! e  ! <-- ! longueur du tableau nodfbr                     !
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
! kface            ! e  ! <-- ! numero de la face d'interaction                !
! nbpt             ! e  ! <-- ! numero de la particule traitee                 !
! isuivi           ! e  ! --> ! indicateur de suivi de la particule            !
! nideve, nrdeve   ! i  ! <-- ! sizes of idevel and rdevel arrays              !
! nituse, nrtuse   ! i  ! <-- ! sizes of ituser and rtuser arrays              !
! ifacel(2, nfac)  ! ia ! <-- ! interior faces -> cells connectivity           !
! ifabor(nfabor)   ! ia ! <-- ! boundary faces -> cells connectivity           !
! ifmfbr(nfabor)   ! ia ! <-- ! boundary face family numbers                   !
! ifmcel(ncelet)   ! ia ! <-- ! cell family numbers                            !
! iprfml           ! te ! <-- ! proprietes d'une famille                       !
!  (nfml,nprfml    !    !     !                                                !
! ipnfac           ! te ! <-- ! position du premier noeud de chaque            !
!   (lndfac)       !    !     !  face interne dans nodfac                      !
! nodfac           ! te ! <-- ! connectivite faces internes/noeuds             !
!   (nfac+1)       !    !     !                                                !
! ipnfbr           ! te ! <-- ! position du premier noeud de chaque            !
!   (lndfbr)       !    !     !  face de bord dans nodfbr                      !
! nodfbr           ! te ! <-- ! connectivite faces de bord/noeuds              !
!   (nfabor+1)     !    !     !                                                !
! itypfb(nfabor    ! te ! <-- ! type des faces de bord                         !
!  nphas)          !    !     !                                                !
! itrifb(nfabor    ! te ! <-- ! tab d'indirection pour tri des faces           !
!   nphas)         !    !     !                                                !
! ifrlag           ! te ! <-- ! numero de zone de la face de bord              !
!   (nfabor)       !    !     !  pour le module lagrangien                     !
! itepa            ! te ! --> ! info particulaires (entiers)                   !
! (nbpmax,nivep    !    !     !   (cellule de la particule,...)                !
! indep            ! te ! --> ! pour chaque particule :                        !
!   (nbpmax)       !    !     !   numero de la cellule de depart               !
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
! surfbn(nfabor    ! tr ! <-- ! surface des faces de bord                      !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! tr ! <-- ! variables de calcul au centre des              !
! (ncelet,*)       !    !     !    cellules (instant courant et prec)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
! ettp             ! tr ! <-- ! tableaux des variables liees                   !
!  (nbpmax,nvp)    !    !     !   aux particules etape courante                !
! ettpa            ! tr ! <-- ! tableaux des variables liees                   !
!  (nbpmax,nvp)    !    !     !   aux particules etape precedente              !
! tepa             ! tr ! <-- ! info particulaires (reels)                     !
! (nbpmax,nvep)    !    !     !   (poids statistiques,...)                     !
! parbor(nfabor    ! tr ! <-- ! cumul des statistiques aux frontieres          !
!    nvisbr)       !    !     !                                                !
! vitpar           ! tr ! <-- ! vitesse particule pour le traitement           !
!   (nbpmax,3)     !    !     !   interactions particules/frontieres           !
! vitflu           ! tr ! <-- ! vitesse fluide vu pour le traitement           !
!   (nbpmax,3)     !    !     !   interactions particules/frontieres           !
! auxl(nbpmax,3    ! tr ! --- ! tableau de travail                             !
! rdevel(nrdeve)   ! ra ! <-> ! real work array for temporary development      !
! rtuser(nrtuse)   ! ra ! <-> ! user-reserved real work array                  !
! ra(*)            ! ra ! --- ! main real work array                           !
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================

implicit none

!===============================================================================
! Common blocks
!===============================================================================

include "paramx.h"
include "cstnum.h"
include "numvar.h"
include "optcal.h"
include "entsor.h"
include "cstphy.h"
include "pointe.h"
include "period.h"
include "parall.h"
include "lagpar.h"
include "lagran.h"
include "ppppar.h"
include "ppthch.h"
include "cpincl.h"

!===============================================================================

! Arguments

integer          idbia0 , idbra0
integer          ndim   , ncelet , ncel   , nfac   , nfabor
integer          nfml   , nprfml
integer          nnod   , lndnod , lndfac , lndfbr , ncelbr
integer          nvar   , nscal  , nphas
integer          nbpmax , nvp    , nvp1   , nvep  , nivep
integer          ntersl , nvlsta , nvisbr
integer          kface  , nbpt   , isuivi
integer          nideve , nrdeve , nituse , nrtuse
integer          ifacel(2,nfac) , ifabor(nfabor)
integer          ifmfbr(nfabor) , ifmcel(ncelet)
integer          iprfml(nfml,nprfml)
integer          ipnfac(nfac+1) , nodfac(lndfac)
integer          ipnfbr(nfabor+1) , nodfbr(lndfbr)
integer          itypfb(nfabor,nphas) , itrifb(nfabor,nphas)
integer          ifrlag(nfabor) , itepa(nbpmax,nivep)
integer          indep(nbpmax)
integer          idevel(nideve) , ituser(nituse)
integer          ia(*)

double precision xyzcen(ndim,ncelet)
double precision surfac(ndim,nfac) , surfbo(ndim,nfabor)
double precision cdgfac(ndim,nfac) , cdgfbo(ndim,nfabor)
double precision xyznod(ndim,nnod) , volume(ncelet)
double precision surfbn(nfabor)
double precision dt(ncelet) , rtp(ncelet,*) , rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*) , propfb(nfabor,*)
double precision coefa(nfabor,*) , coefb(nfabor,*)
double precision ettp(nbpmax,nvp) , ettpa(nbpmax,nvp)
double precision tepa(nbpmax,nvep)
double precision parbor(nfabor,nvisbr) , auxl(nbpmax,3)
double precision vitpar(nbpmax,3) , vitflu(nbpmax,3)
double precision rdevel(nrdeve) , rtuser(nrtuse)
double precision ra(*)

! Local variables

integer          idebia, idebra
integer          ip , nfin , kzone , n1 , icha, iok

double precision aa
double precision xp , yp , zp
double precision xq , yq , zq
double precision xk , yk , zk
double precision xpq , ypq , zpq
double precision xnn , ynn , znn
double precision vnorl(1)  , enc3 , viscp , masse
double precision dpinit , dp03 , mp0 , trap , vnorm , ang
double precision energ , energt
double precision uxn   , vyn    , wzn

!===============================================================================

!===============================================================================
! 0.  GESTION MEMOIRE
!===============================================================================

idebia = idbia0
idebra = idbra0

!===============================================================================
! 1. Traitement en fonction du type de frontiere
!===============================================================================

iok = 0

!--> numero de la particule traitee

ip = nbpt

!--> Zone de la face de bord a traitee

kzone = ifrlag(kface)

!-->normale normee sortante de la face KFACE

aa = 1.d0 / surfbn(kface)
xnn = surfbo(1,kface) * aa
ynn = surfbo(2,kface) * aa
znn = surfbo(3,kface) * aa

!===============================================================================
! 2. Recherche du point d'intersection entre la face de bord et
!    le rayon. Les coordonnees sont stockees dans XK YK ZK
!===============================================================================

! Petit rappel de geometrie 3D :
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!            1)  equation d'un plan de normal (a,b,c) :
!                a*x + b*y + c*z + d = 0
!            2)  equation d'une droite qui passe par P et Q :
!                x = XP + (XQ-XP) * AA
!                y = YP + (YQ-YP) * AA
!                z = ZP + (ZQ-ZP) * AA
!                ou AA est un parametre qui varie dans l'ensemble des reels

!-->on determine le vecteur PQ

xp = ettpa(ip,jxp)
yp = ettpa(ip,jyp)
zp = ettpa(ip,jzp)

xq = ettp(ip,jxp)
yq = ettp(ip,jyp)
zq = ettp(ip,jzp)

xpq = xq - xp
ypq = yq - yp
zpq = zq - zp

!-->si la particule n'a pas bougee ( si elle est deposee
!   sur la face de bord), elle n'est pas traitee.

if (xpq.eq.0.d0 .and. ypq.eq.0.d0 .and. zpq.eq.0.d0) return

!-->a partir de l'equation du plan de la face et de l'equation
!   parametre du rayon, on determine le parametre du point d'intersection

aa = xpq * surfbo(1,kface)                                        &
   + ypq * surfbo(2,kface)                                        &
   + zpq * surfbo(3,kface)

if ( aa.eq.0.d0 ) then
  write (nfecra,9010) ip
  nbperr = nbperr + 1
  dnbper = dnbper + tepa(ip,jrpoi)
  isuivi = 0
  itepa(ip,jisor) = 0
  return
endif

aa =                                                              &
     ( surfbo(1,kface) * cdgfbo(1,kface)                          &
     + surfbo(2,kface) * cdgfbo(2,kface)                          &
     + surfbo(3,kface) * cdgfbo(3,kface)                          &
     - surfbo(1,kface) * xp                                       &
     - surfbo(2,kface) * yp                                       &
     - surfbo(3,kface) * zp )                                     &
     / aa

!-->on reporte ce parametre dans l'equation de la droite du rayon pour
!   avoir le point d'intersection (XK YK ZK)

xk = xp + xpq * aa
yk = yp + ypq * aa
zk = zp + zpq * aa

!===============================================================================
! 3. Depart de la particule du domaine de calcul,
!    ou deposition de la particule sur la frontiere
!===============================================================================

if ( iusclb(kzone).eq.isortl .or.                                 &
     iusclb(kzone).eq.ientrl .or.                                 &
     iusclb(kzone).eq.idepo1      ) then

  isuivi = 0
  itepa(ip,jisor) = 0

!      mise a jour du debit

  deblag(kzone) = deblag(kzone)-tepa(ip,jrpoi)*ettp(ip,jmp)

!--> La particule sort mais pour la visualisation ensight, on la place
!    correctement au point d'intersection.
!    Possibilite que le numero IP soit reutilise pour une autre
!    particule a entrer.

  ettp(ip,jxp) = xk
  ettp(ip,jyp) = yk
  ettp(ip,jzp) = zk

!===============================================================================
! 4. Deposition de la particule, celle-ci reste en memoire
!===============================================================================

else if (iusclb(kzone).eq.idepo2) then

!--> La particule ne sort pas, elle n'est plus traitee,
!    mais toujours visualisee. Le numero IP n'est pas reutilisable.

    isuivi = 0
    itepa(ip,jisor) = -itepa(ip,jisor)
  ettp(ip,jxp) = xk
  ettp(ip,jyp) = yk
  ettp(ip,jzp) = zk

  do n1 = 1,3
    vitpar(ip,n1) = 0.d0
    vitflu(ip,n1) = 0.d0
  enddo

!===============================================================================
! 5. Deposition de la particule, la remise en suspension est possible
!===============================================================================

else if (iusclb(kzone).eq.idepo3) then

  isuivi = 0
  itepa(ip,jisor) = ifabor(kface)
  ettp(ip,jxp) = xk
  ettp(ip,jyp) = yk
  ettp(ip,jzp) = zk

  do n1 = 1,3
    vitpar(ip,n1) = 0.d0
    vitflu(ip,n1) = 0.d0
  enddo
  do n1 = jup,jwf
    ettpa(ip,n1) = 0.d0
  enddo

!===============================================================================
! 6. Deposition de la particule avec force d'attachement,
!    vitesse conservee et re-entrainement possible
!===============================================================================

else if (iusclb(kzone).eq.idepfa) then


! Calcul du critere

  uxn = ettp(ip,jup)*xnn
  vyn = ettp(ip,jvp)*ynn
  wzn = ettp(ip,jwp)*znn

  energ = 0.5d0*ettp(ip,jmp)*(uxn*uxn+vyn*vyn+wzn*wzn)

  energt   = 3.34d-12*ettp(ip,jdp)

  if ( energ .ge. energt )then

    isuivi = 0
    itepa(ip,jisor) = ifabor(kface)
    ettp(ip,jxp) = xk
    ettp(ip,jyp) = yk
    ettp(ip,jzp) = zk

    vitpar(ip,1) = 0.d0
    vitpar(ip,2) = 0.d0
    vitpar(ip,3) = 0.d0

    vitflu(ip,1) = 0.d0
    vitflu(ip,2) = 0.d0
    vitflu(ip,3) = 0.d0

  else

    isuivi = 0
    itepa(ip,jisor) = indep(ip)
    ettp(ip,jxp) = ettpa(ip,jxp)
    ettp(ip,jyp) = ettpa(ip,jyp)
    ettp(ip,jzp) = ettpa(ip,jzp)

!-->changement de la vitesse de la particule au point d'arrive
!   (comme pour une rebond elastique)

    aa = abs(( vitpar(ip,1)*xnn                                   &
              +vitpar(ip,2)*ynn                                   &
              +vitpar(ip,3)*znn) )*2.d0

    vitpar(ip,1) = vitpar(ip,1) - aa*xnn
    vitpar(ip,2) = vitpar(ip,2) - aa*ynn
    vitpar(ip,3) = vitpar(ip,3) - aa*znn

!-->Annule  la vitesse du fluide vu au point d'arrive
!   (comme pour une rebond elastique)

    aa = abs( (vitflu(ip,1)*xnn                                   &
             + vitflu(ip,2)*ynn                                   &
             + vitflu(ip,3)*znn) ) * 2.d0

    vitflu(ip,1) = 0.d0
    vitflu(ip,2) = 0.d0
    vitflu(ip,3) = 0.d0
  endif

!===============================================================================
! 7. Rebond elastique de la particule sur la frontiere
!===============================================================================

else if (iusclb(kzone).eq.irebol) then

  isuivi = 1
  itepa(ip,jisor) = ifabor(kface)

!-->changement du point de depart

  ettpa(ip,jxp) = xk
  ettpa(ip,jyp) = yk
  ettpa(ip,jzp) = zk

    if (iensi1.eq.1) then
      nfin = 0
      call enslag                                                 &
      !==========
       ( idebia, idebra  ,                                        &
         nbpmax , nvp    , nvp1   , nvep   , nivep  ,             &
         nfin   , ip     ,                                        &
         itepa  ,                                                 &
         ettpa  , tepa   , ra)
    endif

!-->changement du point d'arrive
!   (la valeur absolue est pour eviter les p. scalaires
!    negatifs impossibles qui interviennent avec des erreurs d'arrondi
!    machines)

    aa = 2.d0 * abs( (xq-xk)*xnn + (yq-yk)*ynn + (zq-zk)*znn )

  ettp(ip,jxp) = xq - aa*xnn
  ettp(ip,jyp) = yq - aa*ynn
  ettp(ip,jzp) = zq - aa*znn

!-->changement de la vitesse de la particule au point d'arrive


  aa = abs( (vitpar(ip,1)*xnn                                     &
          +  vitpar(ip,2)*ynn                                     &
          +  vitpar(ip,3)*znn) ) * 2.d0

  vitpar(ip,1) = vitpar(ip,1) - aa*xnn
  vitpar(ip,2) = vitpar(ip,2) - aa*ynn
  vitpar(ip,3) = vitpar(ip,3) - aa*znn

!-->changement de la vitesse du fluide vu au point d'arrive

  aa = abs( (vitflu(ip,1)*xnn                                     &
          +  vitflu(ip,2)*ynn                                     &
          +  vitflu(ip,3)*znn) ) * 2.d0

  vitflu(ip,1) = vitflu(ip,1) - aa*xnn
  vitflu(ip,2) = vitflu(ip,2) - aa*ynn
  vitflu(ip,3) = vitflu(ip,3) - aa*znn

!===============================================================================
! 8. Encrassement des grains de charbon
!===============================================================================

else if (iusclb(kzone).eq.iencrl) then

!--> Encrassement de la particules si ses proprietes le permettent
!    et en fonction d'une probabilite
!      ICI c'est si Tp     > TPENC
!                si VISCP  > VISCREF

  icha = itepa(ip,jinch)

  if ( ettp(ip,jhp).gt.tprenc(icha) ) then

    enc3 = ( (1.d+7 * enc1(icha))/((ettp(ip,jhp)-150.d0)**2) )    &
           + enc2(icha)
      viscp = exp( log(10.d0)*enc3 ) * 0.1d0

      if ( viscp.gt.visref(icha) ) then
        n1 = 1
        call zufall(n1,vnorl(1))
        trap = 1.d0-visref(icha) / viscp
      endif

!  Si VISCP <= VISREF ===> 100% de chance d'encrasser
!  Si VISCP  > VISREF ===> probabilte TRAP = 1-VISREF/VISCP
!                          d'encrasser
!                     ===> On encrasse si VNORL est compris
!                          entre TRAP et 1.

      if ( viscp.le.visref(icha) .or.                             &
         (viscp.gt.visref(icha) .and. vnorl(1).ge.trap) ) then

! La calcul de la masse de grains de charbon encrasses est faite plus bas.

        npencr = npencr + 1
        isuivi = 0
        itepa(ip,jisor)  =  0
      ettp(ip,jxp) = xk
      ettp(ip,jyp) = yk
      ettp(ip,jzp) = zk

      endif
    endif

!--> Si pas encrassement alors rebond elastique

    if ( itepa(ip,jisor).ne.0 ) then

      isuivi = 1
      itepa(ip,jisor) = ifabor(kface)

!-->changement du point de depart

    ettpa(ip,jxp) = xk
    ettpa(ip,jyp) = yk
    ettpa(ip,jzp) = zk

      if (iensi1.eq.1) then
        nfin = 0
        call enslag                                               &
        !==========
       ( idebia, idebra  ,                                        &
         nbpmax , nvp    , nvp1   , nvep   , nivep  ,             &
         nfin   , ip     ,                                        &
         itepa  ,                                                 &
         ettpa  , tepa   , ra)
      endif

!-->changement du point d'arrive

      aa = 2.d0 * abs((xq-xk)*xnn + (yq-yk)*ynn + (zq-zk)*znn)

    ettp(ip,jxp) = xq - aa*xnn
    ettp(ip,jyp) = yq - aa*ynn
    ettp(ip,jzp) = zq - aa*znn

    endif

  if (itepa(ip,jisor).gt.0) then

!-->changement de la vitesse de la particule au point d'arrive


    aa = abs( (vitpar(ip,1)*xnn                                   &
            +  vitpar(ip,2)*ynn                                   &
            +  vitpar(ip,3)*znn) ) * 2.d0

    vitpar(ip,1) = vitpar(ip,1) - aa*xnn
    vitpar(ip,2) = vitpar(ip,2) - aa*ynn
    vitpar(ip,3) = vitpar(ip,3) - aa*znn

!-->changement de la vitesse du fluide vu au point d'arrive

    aa = abs( (vitflu(ip,1)*xnn                                   &
            +  vitflu(ip,2)*ynn                                   &
            +  vitflu(ip,3)*znn) ) * 2.d0

    vitflu(ip,1) = vitflu(ip,1) - aa*xnn
    vitflu(ip,2) = vitflu(ip,2) - aa*ynn
    vitflu(ip,3) = vitflu(ip,3) - aa*znn

    endif

!===============================================================================
! 9. Interaction utilisateur numero 1 : JBORD1
!===============================================================================

!  ON PEUT FAIRE DE MEME AVEC JBORD2, JBORD3, JBORD4 et JBORD5
!    On ne donnne l'exemple que pour JBORD1

!     On regarde si on est dans la zone cherchee
!      ELSE IF (IUSCLB(KZONE).EQ.JBORD1) THEN

!       si on doit continuer a suivre la particule
!         ISUIVI = 0 OU 1

!       l'element de maillage concerne
!         ITEPA(IP,JISOR) =

!       on change la localisation de depart de la particule
!         ETTP(IP,JXP) =
!         ETTP(IP,JYP) =
!         ETTP(IP,JZP) =

!       la vitesse de la particule au point d'arrivee
!         VITPAR(IP,1) =
!         VITPAR(IP,2) =
!         VITPAR(IP,3) =

!       la vitesse du fluide vu au point d'arrivee
!         VITFLU(IP,1) =
!         VITFLU(IP,2) =
!         VITFLU(IP,3) =


else if (iusclb(kzone).eq.jbord1                                  &
! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!    en standard, sans intervention de l'utilisateur,
!      on ne souhaite pas passer ici
!    mais on desire
!      que le test IUSCLB(KZONE).EQ.JBORD1 soit dans le us* exemple
!      que le source ci-dessous soit compile pour verifier
!         qu'il n y a pas d'erreur
         .and.(0.eq.1)                                            &
! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END
                                ) then

!     ----------------------------------------------------
!     EXEMPLE 1 : LA PARTICULE A UNE CHANCE SUR DEUX (TRAP)
!                 DE SE DEPOSER DEFINITIVEMENT A LA PAROI
!                 OU DE REBONDIR VERS L'ECOULEMENT.
!     ----------------------------------------------------


    n1 = 1
    call zufall(n1,vnorl(1))
    trap = 0.5d0

    if (vnorl(1).ge.trap) then

      isuivi = 0
      itepa(ip,jisor)  =  0
      ettp(ip,jxp) = xk
      ettp(ip,jyp) = yk
      ettp(ip,jzp) = zk

    else

      isuivi = 1
      itepa(ip,jisor) = ifabor(kface)

!-->changement du point de depart

      ettpa(ip,jxp) = xk
      ettpa(ip,jyp) = yk
      ettpa(ip,jzp) = zk

!-->changement du point d'arrive

      aa = 2.d0 * abs((xq-xk)*xnn + (yq-yk)*ynn + (zq-zk)*znn)

      ettp(ip,jxp) = xq - aa*xnn
      ettp(ip,jyp) = yq - aa*ynn
      ettp(ip,jzp) = zq - aa*znn

    endif

!-->Inutile de traiter les particule ITEPA(IP,JISOR)=0 car
!   Elle vont etre eliminees de la listes des particules.

    if (itepa(ip,jisor).gt.0) then

!-->changement de la vitesse de la particule au point d'arrive


      aa = abs( (vitpar(ip,1)*xnn                                 &
              +  vitpar(ip,2)*ynn                                 &
              +  vitpar(ip,3)*znn) ) * 2.d0

      vitpar(ip,1) = vitpar(ip,1) - aa*xnn
      vitpar(ip,2) = vitpar(ip,2) - aa*ynn
      vitpar(ip,3) = vitpar(ip,3) - aa*znn

!-->changement de la vitesse du fluide vu au point d'arrive

      aa = abs( (vitflu(ip,1)*xnn                                 &
              +  vitflu(ip,2)*ynn                                 &
              +  vitflu(ip,3)*znn) ) * 2.d0

      vitflu(ip,1) = vitflu(ip,1) - aa*xnn
      vitflu(ip,2) = vitflu(ip,2) - aa*ynn
      vitflu(ip,3) = vitflu(ip,3) - aa*znn

    endif


!===============================================================================
! 10. Verification et sortie si erreur
!===============================================================================

else
  write (nfecra,9020) kzone
  iok = iok + 1
endif

if (iok.ne.0) then
  call csexit (1)
  !==========
endif

!===============================================================================
! 11. Enregistrement de l'interaction particule-frontiere s'il y a lieu
!===============================================================================

!     L'enregistrement des statistiques parietales debute des que
!     l'indicateur est a IENSI3 = 1. Cependant tant que le numero
!     de l'iteration Lagrangienne absolue est inferieur a NSTBOR,
!     ou que l'ecoulement est instationnaire (ISTTIO = 0), alors
!     le tableau PARBOR est remis a zero avant d'entrer dans ce
!     sous-programme.

!     NPSTF   : nombre d'iterations de calcul de stat aux frontieres
!               stationnaires

!     NPSTFT  : nombre d'iterations total des stats iaux frontieres
!               depuis le du calcul, partie instationnaire comprise
!               (a utiliser que pour les affichages listing)

!     TSTATP  : Temps physique d'enregistrement des statistiques
!               des interactions particules/frontiere stationnaires,
!               si instatinnaire alors contient DTP le dernier
!               pas de temps Lagrangien

!     Avant impression dans le listing, ou sortie pour le
!     post-processing, on applique aux statistiques sur les frontieres
!     une moyenne en fonction des infos fournies dans USLAG1 de la
!     maniere suivante :



!     CES LIGNES NE SONT DONNEES QU'A TITRE INDICATIF


!     DO IVAR = 1,NVISBR

!      IF (IMOYBR(IVAR).EQ.2) THEN

!         DO IFAC = 1,NFABOR
!           IF (PARBOR(IFAC,INBR).GT.SEUILF) THEN
!             PARBOR(IFAC,IVAR) = PARBOR(IFAC,IVAR) /PARBOR(IFAC,INBR)
!           ELSE
!             PARBOR(IFAC,IVAR) = 0.D0
!           ENDIF
!         ENDDO

!       ELSE IF (IMOYBR(IVAR).EQ.1) THEN

!         DO IFAC = 1,NFABOR
!           IF (PARBOR(IFAC,INBR).GT.SEUILF) THEN
!             PARBOR(IFAC,IVAR) = PARBOR(IFAC,IVAR) / TSTATP
!           ELSE
!             PARBOR(IFAC,IVAR) = 0.D0
!           ENDIF
!         ENDDO
!       ENDIF
!     ENDDO




if ( iensi3.eq.1 ) then

!--> exemple de types d'interaction sur lesquels on souhaite
!    enregistrer des informations

  if ( iusclb(kzone).eq.irebol .or.                               &
       iusclb(kzone).eq.idepo1 .or.                               &
       iusclb(kzone).eq.idepo2 .or.                               &
       iusclb(kzone).eq.idepo3      ) then

    if (inbrbd.eq.1) then
      parbor(kface,inbr) = parbor(kface,inbr) + tepa(ip,jrpoi)
    endif

    if (iflmbd.eq.1) then
        parbor(kface,iflm) = parbor(kface,iflm)                   &
     + ( tepa(ip,jrpoi) * ettp(ip,jmp) /surfbn(kface) )
    endif

    if (iangbd.eq.1) then
      vnorm = ettp(ip,jup) * ettp(ip,jup)                         &
            + ettp(ip,jvp) * ettp(ip,jvp)                         &
            + ettp(ip,jwp) * ettp(ip,jwp)
      vnorm = sqrt( vnorm )
      ang =  ettp(ip,jup) * surfbo(1,kface)                       &
           + ettp(ip,jvp) * surfbo(2,kface)                       &
           + ettp(ip,jwp) * surfbo(3,kface)                       &
           / surfbn(kface)                                        &
           / vnorm
      ang = acos(ang)
      parbor(kface,iang) = parbor(kface,iang) + ang*tepa(ip,jrpoi)
    endif

    if (ivitbd.eq.1) then
      vnorm = ettp(ip,jup) * ettp(ip,jup)                         &
            + ettp(ip,jvp) * ettp(ip,jvp)                         &
            + ettp(ip,jwp) * ettp(ip,jwp)
      vnorm = sqrt( vnorm )
      parbor(kface,ivit) =parbor(kface,ivit) +vnorm*tepa(ip,jrpoi)
    endif

    if (nusbor.gt.0) then
      do n1 = 1,nusbor
        parbor(kface,iusb(n1)) = 0.d0
      enddo
    endif

!--> Cas particulier de la masse de grains de charbon encrasses

  else if ( iusclb(kzone).eq.iencrl .and. isuivi.eq.0 ) then

    parbor(kface,inbr) = parbor(kface,inbr) + tepa(ip,jrpoi)

    if (iencbd.eq.1) then

      icha =  itepa(ip,jinch)
      dpinit = tepa(ip,jrd0p)

      dp03 = dpinit * dpinit * dpinit
      mp0  = pi * dp03 * rho0ch(icha) / 6.d0

      masse = ettp(ip,jmch) + ettp(ip,jmck)                       &
                                + xashch(icha) * mp0

      parbor(kface,ienc) = parbor(kface,ienc)                     &
                         + tepa(ip,jrpoi)*masse

    endif

  endif

endif


!===============================================================================
! Archives. Je laisse cette partie au cas ou...
!           Creation d'un repere local associe a la face frontiere
!===============================================================================

! Le repere local (T1,T2,N) est construit de facon a ce que N soit
! la normale normee de la face et que T1 et T2 appartiennent a la face

!-->1. Je connais les vecteurs N et PK, je definis T1 tel que
!      T1 = PK vectoriel N

!     XPK = XK - ETTPA(IP,JXP)
!     YPK = YK - ETTPA(IP,JYP)
!     ZPK = ZK - ETTPA(IP,JZP)

!     XT1 = YPK*ZNN - ZPK*YNN
!     YT1 = ZPK*XNN - XPK*ZNN
!     ZT1 = XPK*YNN - YPK*XNN

!     AA = SQRT(XT1*XT1 + YT1*YT1 + ZT1*ZT1)
!     XT1 = XT1 / AA
!     YT1 = YT1 / AA
!     ZT1 = ZT1 / AA

!-->2. Maintenant je peux construire T2 = - T1 vectoriel N

!     XT2 = YT1*ZNN - ZT1*YNN
!     YT2 = ZT1*XNN - XT1*ZNN
!     ZT2 = XT1*YNN - YT1*XNN

!     AA = SQRT(XT2*XT2 + YT2*YT2 + ZT2*ZT2)
!     XT2 = -XT2 / AA
!     YT2 = -YT2 / AA
!     ZT2 = -ZT2 / AA

!--------
! FORMATS
!--------

 9010 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========   (USLABO)                                    ',/,&
'@                                                            ',/,&
'@  La normale d''une face de frontiere est perpendiculaire   ',/,&
'@   a un rayon PQ : impossible                               ',/,&
'@                                                            ',/,&
'@  La particle ',I10,' est ELIMINEE                          ',/,&
'@                                                            ',/,&
'@  Contacter l''equipe de developpement.                     ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 9020 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========   (USLABO)                                    ',/,&
'@                                                            ',/,&
'@  Le type de condition aux limites IUSCLB                   ',/,&
'@    n''est pas renseigne pour la frontiere NB = ',I10        ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier USLAG2 et USLABO.                                ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

!----
! FIN
!----

return
end subroutine
