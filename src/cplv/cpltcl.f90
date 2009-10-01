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

subroutine cpltcl &
!================

 ( idbia0 , idbra0 ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   icodcl , itrifb , itypfb , izfppp ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  , rcodcl ,                                     &
   w1     , w2     , w3     , w4     , w5     , w6     , coefu  , &
   rdevel , rtuser , ra     )

!===============================================================================
! FONCTION :
! --------

!   SOUS-PROGRAMME DU MODULE LAGRANGIEN COUPLE CHARBON PULVERISE :
!   --------------------------------------------------------------

!    ROUTINE UTILISATEUR POUR PHYSIQUE PARTICULIERE

!      COMBUSTION EULERIENNE DE CHARBON PULVERISE ET
!      TRANSPORT LAGRANGIEN DES PARTICULES DE CHARBON

!      CONDITIONS AUX LIMITES AUTOMATIQUES

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
!    nom           !type!mode !                   role                         !
!__________________!____!_____!________________________________________________!
! idbia0           ! e  ! <-- ! numero de la 1ere case libre dans ia           !
! idbra0           ! e  ! <-- ! numero de la 1ere case libre dans ra           !
! ndim             ! e  ! <-- ! dimension de l'espace                          !
! ncelet           ! e  ! <-- ! nombre d'elements halo compris                 !
! ncel             ! e  ! <-- ! nombre d'elements actifs                       !
! nfac             ! e  ! <-- ! nombre de faces internes                       !
! nfabor           ! e  ! <-- ! nombre de faces de bord                        !
! nfml             ! e  ! <-- ! nombre de familles d entites                   !
! nprfml           ! e  ! <-- ! nombre de proprietese des familles             !
! nnod             ! e  ! <-- ! nombre de sommets                              !
! lndfac           ! e  ! <-- ! longueur du tableau nodfac                     !
! lndfbr           ! e  ! <-- ! longueur du tableau nodfbr                     !
! ncelbr           ! e  ! <-- ! nombre d'elements ayant au moins une           !
!                  !    !     ! face de bord                                   !
! nvar             ! e  ! <-- ! nombre total de variables                      !
! nscal            ! e  ! <-- ! nombre total de scalaires                      !
! nphas            ! e  ! <-- ! nombre de phases                               !
! nideve nrdeve    ! e  ! <-- ! longueur de idevel rdevel                      !
! nituse nrtuse    ! e  ! <-- ! longueur de ituser rtuser                      !
! ifacel           ! te ! <-- ! elements voisins d'une face interne            !
! (2, nfac)        !    !     !                                                !
! ifabor           ! te ! <-- ! element  voisin  d'une face de bord            !
! (nfabor)         !    !     !                                                !
! ifmfbr           ! te ! <-- ! numero de famille d'une face de bord           !
! (nfabor)         !    !     !                                                !
! ifmcel           ! te ! <-- ! numero de famille d'une cellule                !
! (ncelet)         !    !     !                                                !
! iprfml           ! te ! <-- ! proprietes d'une famille                       !
! nfml  ,nprfml    !    !     !                                                !
! ipnfac           ! te ! <-- ! position du premier noeud de chaque            !
!   (lndfac)       !    !     !  face interne dans nodfac                      !
! nodfac           ! te ! <-- ! connectivite faces internes/noeuds             !
!   (nfac+1)       !    !     !                                                !
! ipnfbr           ! te ! <-- ! position du premier noeud de chaque            !
!   (lndfbr)       !    !     !  face de bord dans nodfbr                      !
! nodfbr           ! te ! <-- ! connectivite faces de bord/noeuds              !
!   (nfabor+1)     !    !     !                                                !
! icodcl           ! te ! --> ! code de condition limites aux faces            !
!  (nfabor,nvar    !    !     !  de bord                                       !
!                  !    !     ! = 1   -> dirichlet                             !
!                  !    !     ! = 3   -> densite de flux                       !
!                  !    !     ! = 4   -> glissemt et u.n=0 (vitesse)           !
!                  !    !     ! = 5   -> frottemt et u.n=0 (vitesse)           !
!                  !    !     ! = 6   -> rugosite et u.n=0 (vitesse)           !
!                  !    !     ! = 9   -> entree/sortie libre (vitesse          !
!                  !    !     !  entrante eventuelle     bloquee               !
! itrifb(nfabor    ! te ! <-- ! indirection pour tri des faces de brd          !
!  nphas      )    !    !     !                                                !
! itypfb(nfabor    ! te ! <-- ! type des faces de bord                         !
!  nphas      )    !    !     !                                                !
! izfppp           ! te ! <-- ! numero de zone de la face de bord              !
! (nfabor)         !    !     !  pour le module phys. part.                    !
! idevel(nideve    ! te ! <-- ! tab entier complementaire developemt           !
! ituser(nituse    ! te ! <-- ! tab entier complementaire utilisateur          !
! ia(*)            ! tr ! --- ! macro tableau entier                           !
! xyzcen           ! tr ! <-- ! point associes aux volumes de control          !
! (ndim,ncelet     !    !     !                                                !
! surfac           ! tr ! <-- ! vecteur surface des faces internes             !
! (ndim,nfac)      !    !     !                                                !
! surfbo           ! tr ! <-- ! vecteur surface des faces de bord              !
! (ndim,nfabor)    !    !     !                                                !
! cdgfac           ! tr ! <-- ! centre de gravite des faces internes           !
! (ndim,nfac)      !    !     !                                                !
! cdgfbo           ! tr ! <-- ! centre de gravite des faces de bord            !
! (ndim,nfabor)    !    !     !                                                !
! xyznod           ! tr ! <-- ! coordonnes des noeuds                          !
! (ndim,nnod)      !    !     !                                                !
! volume           ! tr ! <-- ! volume d'un des ncelet elements                !
! (ncelet          !    !     !                                                !
! dt(ncelet)       ! tr ! <-- ! pas de temps                                   !
! rtp, rtpa        ! tr ! <-- ! variables de calcul au centre des              !
! (ncelet,*)       !    !     !    cellules (instant courant ou prec)          !
! propce           ! tr ! <-- ! proprietes physiques au centre des             !
! (ncelet,*)       !    !     !    cellules                                    !
! propfa           ! tr ! <-- ! proprietes physiques au centre des             !
!  (nfac,*)        !    !     !    faces internes                              !
! propfb           ! tr ! <-- ! proprietes physiques au centre des             !
!  (nfabor,*)      !    !     !    faces de bord                               !
! coefa, coefb     ! tr ! <-- ! conditions aux limites aux                     !
!  (nfabor,*)      !    !     !    faces de bord                               !
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
! w1,2,3,4,5,6     ! tr ! --- ! tableaux de travail                            !
!  (ncelet         !    !     !  (calcul du gradient de pression)              !
! coefu            ! tr ! --- ! tab de trav                                    !
!  (nfabor,3)      !    !     !  (calcul du gradient de pression)              !
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

! Arguments

include "paramx.h"
include "numvar.h"
include "optcal.h"
include "cstphy.h"
include "cstnum.h"
include "entsor.h"
include "parall.h"
include "ppppar.h"
include "ppthch.h"
include "coincl.h"
include "cpincl.h"
include "ppincl.h"
include "ppcpfu.h"

!===============================================================================

! Arguments

integer          idbia0 , idbra0
integer          ndim   , ncelet , ncel   , nfac   , nfabor
integer          nfml   , nprfml
integer          nnod   , lndfac , lndfbr , ncelbr
integer          nvar   , nscal  , nphas
integer          nideve , nrdeve , nituse , nrtuse

integer          ifacel(2,nfac) , ifabor(nfabor)
integer          ifmfbr(nfabor) , ifmcel(ncelet)
integer          iprfml(nfml,nprfml)
integer          ipnfac(nfac+1), nodfac(lndfac)
integer          ipnfbr(nfabor+1), nodfbr(lndfbr)
integer          icodcl(nfabor,nvar)
integer          itrifb(nfabor,nphas), itypfb(nfabor,nphas)
integer          izfppp(nfabor)
integer          idevel(nideve), ituser(nituse), ia(*)

double precision xyzcen(ndim,ncelet)
double precision surfac(ndim,nfac), surfbo(ndim,nfabor)
double precision cdgfac(ndim,nfac), cdgfbo(ndim,nfabor)
double precision xyznod(ndim,nnod), volume(ncelet)
double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision rcodcl(nfabor,nvar,3)
double precision w1(ncelet),w2(ncelet),w3(ncelet)
double precision w4(ncelet),w5(ncelet),w6(ncelet)
double precision coefu(nfabor,ndim)
double precision rdevel(nrdeve), rtuser(nrtuse), ra(*)

! VARIABLES LOCALES

integer          idebia, idebra
integer          ii, iphas, ifac, izone, mode, iel, ige, iok
integer          icha, ipbrom, icke, ipcvis
integer          nbrval
double precision qisqc, viscla, d2s3, uref2, rhomoy, dhy, xiturb
double precision ustar2, xkent, xeent, t1
double precision h1    (nozppm)
double precision qimpc (nozppm) , qcalc(nozppm)
double precision coefe (ngazem)
double precision f1mc  (ncharm) , f2mc (ncharm)

!===============================================================================
!===============================================================================
! 1.  INITIALISATIONS
!===============================================================================

idebia = idbia0
idebra = idbra0

iphas = 1
ipbrom = ipprob(irom  (iphas))
ipcvis = ipproc(iviscl(iphas))

d2s3 = 2.d0/3.d0

!===============================================================================
! 1.  ECHANGES EN PARALLELE POUR LES DONNEES UTILISATEUR
!===============================================================================

!  En realite on pourrait eviter cet echange en modifiant uscpcl et en
!    demandant a l'utilisateur de donner les grandeurs dependant de la
!    zone hors de la boucle sur les faces de bord : les grandeurs
!    seraient ainsi disponibles sur tous les processeurs. Cependant,
!    ca rend le sous programme utilisateur un peu plus complique et
!    surtout, si l'utilisateur le modifie de travers, ca ne marche pas.
!  On suppose que toutes les grandeurs fournies sont positives, ce qui
!    permet d'utiliser un max pour que tous les procs les connaissent.
!    Si ce n'est pas le cas, c'est plus complique mais on peut s'en tirer
!    avec un max quand meme.

if(irangp.ge.0) then
  call parimx(nozapm,iqimp )
  !==========
  call parimx(nozapm,ientat)
  !==========
  call parimx(nozapm,ientcp)
  !==========
  call parrmx(nozapm,qimpat)
  !==========
  call parrmx(nozapm,timpat)
  !==========
  nbrval = nozppm*ncharm
  call parrmx(nbrval,qimpcp)
  !==========
  nbrval = nozppm*ncharm
  call parrmx(nbrval,timpcp)
  !==========
  nbrval = nozppm*ncharm*nclcpm
  call parrmx(nbrval,distch)
  !==========
endif


!===============================================================================
! 2.  CORRECTION DES VITESSES (EN NORME) POUR CONTROLER LES DEBITS
!     IMPOSES
!       ON BOUCLE SUR TOUTES LES FACES D'ENTREE
!                     =========================
!===============================================================================

! --- Debit calcule

do izone = 1, nozppm
  qcalc(izone) = 0.d0
enddo
do ifac = 1, nfabor
  izone = izfppp(ifac)
  if (izone .gt. 0) then
    if ( iqimp(izone).eq.2) then
      qcalc(izone) = qcalc(izone) -                               &
         ( rcodcl(ifac,iu(iphas),1)*surfbo(1,ifac) +              &
           rcodcl(ifac,iv(iphas),1)*surfbo(2,ifac) +              &
           rcodcl(ifac,iw(iphas),1)*surfbo(3,ifac) )
    else
      qcalc(izone) = qcalc(izone) - propfb(ifac,ipbrom) *         &
         ( rcodcl(ifac,iu(iphas),1)*surfbo(1,ifac) +              &
           rcodcl(ifac,iv(iphas),1)*surfbo(2,ifac) +              &
           rcodcl(ifac,iw(iphas),1)*surfbo(3,ifac) )
    endif
  endif
enddo

if(irangp.ge.0) then
  call parrsm(nozapm,qcalc )
endif

do izone = 1, nozapm
  if ( iqimp(izone).eq.0 ) then
    qimpc(izone) = qcalc(izone)
  endif
enddo

! --- Correction des vitesses en norme

iok = 0
do ii = 1, nzfppp
  izone = ilzppp(ii)
  if ( iqimp(izone).eq.1  .or.  iqimp(izone).eq.2) then
    if(abs(qcalc(izone)).lt.epzero) then
      write(nfecra,2001)izone,iqimp(izone),qcalc(izone)
      iok = iok + 1
    endif
  endif
enddo
if(iok.ne.0) then
  call csexit (1)
  !==========
endif
do ifac = 1, nfabor
  izone = izfppp(ifac)
   if (izone .gt. 0) then
     if ( iqimp(izone).eq.1 .or.  iqimp(izone).eq.2) then
       qimpc(izone) = qimpat(izone)
       qisqc = qimpc(izone)/qcalc(izone)
       rcodcl(ifac,iu(iphas),1) = rcodcl(ifac,iu(iphas),1)*qisqc
       rcodcl(ifac,iv(iphas),1) = rcodcl(ifac,iv(iphas),1)*qisqc
       rcodcl(ifac,iw(iphas),1) = rcodcl(ifac,iw(iphas),1)*qisqc
    endif
  endif

enddo



 2001 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : MODULE PHYSIQUES PARTICULIERES              ',/,&
'@    =========                                               ',/,&
'@     COMBUSTION CHARBON PULVERISE COUPLE AU                 ',/,&
'@     TRANSPORT LAGRANGIEN DES PARTICULES DE CHARBON :       ',/,&
'@     PROBLEME DANS LES CONDITIONS AUX LIMITES                ',/&
'@                                                            ',/,&
'@  Le debit est impose sur la zone IZONE =     ', I10         ,/,&
'@    puisque                IQIMP(IZONE) =     ', I10         ,/,&
'@  Or, sur cette zone, le produit RHO D S integre est nul :  ',/,&
'@    il vaut                             = ',E14.5            ,/,&
'@    (D est la direction selon laquelle est impose le debit).',/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier uscpcl, et en particulier                        ',/,&
'@    - que le vecteur  RCODCL(IFAC,IU(IPHAS),1),             ',/,&
'@                      RCODCL(IFAC,IV(IPHAS),1),             ',/,&
'@                      RCODCL(IFAC,IW(IPHAS),1) qui determine',/,&
'@      la direction de la vitesse est non nul et n''est pas  ',/,&
'@      uniformement perpendiculaire aux face d''entree       ',/,&
'@    - que la surface de l''entree n''est pas nulle (ou que  ',/,&
'@      le nombre de faces de bord dans la zone est non nul)  ',/,&
'@    - que la masse volumique n''est pas nulle               ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

!===============================================================================
! 3.  REMPLISSAGE DU TABLEAU DES CONDITIONS LIMITES
!       ON BOUCLE SUR TOUTES LES FACES D'ENTREE
!                     =========================
!         ON DETERMINE LA FAMILLE ET SES PROPRIETES
!           ON IMPOSE LES CONDITIONS AUX LIMITES
!           POUR LA TURBULENCE

!===============================================================================

do ifac = 1, nfabor

  izone = izfppp(ifac)

!      ELEMENT ADJACENT A LA FACE DE BORD

  if ( itypfb(ifac,iphas).eq.ientre ) then

! ----  Traitement automatique de la turbulence

    if ( icalke(izone).ne.0 ) then

!       La turbulence est calculee par defaut si ICALKE different de 0
!          - soit a partir du diametre hydraulique, d'une vitesse
!            de reference adaptes a l'entree courante si ICALKE = 1
!          - soit a partir du diametre hydraulique, d'une vitesse
!            de reference et de l'intensite turvulente
!            adaptes a l'entree courante si ICALKE = 2

      uref2 = rcodcl(ifac,iu(iphas),1)**2                         &
            + rcodcl(ifac,iv(iphas),1)**2                         &
            + rcodcl(ifac,iw(iphas),1)**2
      uref2 = max(uref2,1.d-12)
      rhomoy = propfb(ifac,ipbrom)
      iel    = ifabor(ifac)
      viscla = propce(iel,ipcvis)
      icke   = icalke(izone)
      dhy    = dh(izone)
      xiturb = xintur(izone)
      ustar2 = 0.d0
      xkent = epzero
      xeent = epzero
      if (icke.eq.1) then
        call keendb                                               &
        !==========
        ( uref2, dhy, rhomoy, viscla, cmu, xkappa,                &
          ustar2, xkent, xeent )
      else if (icke.eq.2) then
        call keenin                                               &
        !==========
        ( uref2, xiturb, dhy, cmu, xkappa, xkent, xeent )
      endif

      if (itytur(iphas).eq.2) then

        rcodcl(ifac,ik(iphas),1)  = xkent
        rcodcl(ifac,iep(iphas),1) = xeent

      elseif (itytur(iphas).eq.3) then

        rcodcl(ifac,ir11(iphas),1) = d2s3*xkent
        rcodcl(ifac,ir22(iphas),1) = d2s3*xkent
        rcodcl(ifac,ir33(iphas),1) = d2s3*xkent
        rcodcl(ifac,ir12(iphas),1) = 0.d0
        rcodcl(ifac,ir13(iphas),1) = 0.d0
        rcodcl(ifac,ir23(iphas),1) = 0.d0
        rcodcl(ifac,iep(iphas),1)  = xeent

      elseif (iturb(iphas).eq.50) then

        rcodcl(ifac,ik(iphas),1)   = xkent
        rcodcl(ifac,iep(iphas),1)  = xeent
        rcodcl(ifac,iphi(iphas),1) = d2s3
        rcodcl(ifac,ifb(iphas),1)  = 0.d0

      elseif (iturb(iphas).eq.60) then

        rcodcl(ifac,ik(iphas),1)   = xkent
        rcodcl(ifac,iomg(iphas),1) = xeent/cmu/xkent

      endif

    endif

  endif

enddo

!===============================================================================
! 2.  REMPLISSAGE DU TABLEAU DES CONDITIONS LIMITES
!       ON BOUCLE SUR TOUTES LES FACES D'ENTREE
!                     =========================
!         ON DETERMINE LA FAMILLE ET SES PROPRIETES
!           ON IMPOSE LES CONDITIONS AUX LIMITES
!           POUR LES SCALAIRES
!===============================================================================

do ii = 1, nzfppp

  izone = ilzppp(ii)

! Une entree IENTRE est forcement du type
!            IENTAT = 1 ou IENTCP = 1

  if ( ientat(izone).eq.1 ) then

! ------ Calcul de H1(IZONE)

    do ige = 1, ngazem
      coefe(ige) = zero
    enddo
    coefe(io2) = wmole(io2) / (wmole(io2)+xsi*wmole(in2))
    coefe(in2) = 1.d0 - coefe(io2)
    do icha = 1, ncharm
      f1mc(icha) = zero
      f2mc(icha) = zero
    enddo
    t1   = timpat(izone)
    mode = -1
    call cpthp1                                                   &
    !==========
    ( mode  , h1(izone) , coefe  , f1mc   , f2mc   ,              &
      t1    )

    endif

enddo

do ifac = 1, nfabor

  izone = izfppp(ifac)

!      ELEMENT ADJACENT A LA FACE DE BORD

  if ( itypfb(ifac,iphas).eq.ientre ) then

! ----  Traitement automatique des scalaires physiques particulieres

    do icha = 1, ncharb

! ------ CL pour F1M et F2M du charbon ICHA
      rcodcl(ifac,isca(if1m(icha)),1) = zero
      rcodcl(ifac,isca(if2m(icha)),1) = zero

    enddo

! ------ CL pour F3M
    rcodcl(ifac,isca(if3m),1) = zero
! ------ CL pour FP4M
    rcodcl(ifac,isca(if4p2m),1)   = zero
! ------ CL pour HM
    do ige = 1, ngazem
      coefe(ige) = zero
    enddo
    coefe(io2) = wmole(io2) / (wmole(io2)+xsi*wmole(in2))
    coefe(in2) = 1.d0 - coefe(io2)
    do icha = 1, ncharm
      f1mc(icha) = zero
      f2mc(icha) = zero
    enddo
    t1   = timpat(izone)
    mode = -1
    call cpthp1                                                   &
    !==========
    ( mode  , h1(izone) , coefe  , f1mc   , f2mc   ,              &
      t1    )
    rcodcl(ifac,isca(ihm),1) = h1(izone)

  endif

enddo


!----
! FORMATS
!----


!----
! FIN
!----

return
end subroutine
