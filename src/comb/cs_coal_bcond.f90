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

subroutine cs_coal_bcond &
!=======================

 ( nvar   , nscal  ,                                              &
   icodcl , itrifb , itypfb , izfppp ,                            &
   dt     , rtp    , rtpa   , propce , propfb ,                   &
   coefa  , coefb  , rcodcl )

!===============================================================================
! FONCTION :
! --------
!    BOUNDARY CONDITION AUTOMATIC FOR PULVERISED COAL COMBUSTION
!
!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! icodcl           ! te ! --> ! code de condition limites aux faces            !
!  (nfabor,nvar    !    !     !  de bord                                       !
!                  !    !     ! = 1   -> dirichlet                             !
!                  !    !     ! = 3   -> densite de flux                       !
!                  !    !     ! = 4   -> glissemt et u.n=0 (vitesse)           !
!                  !    !     ! = 5   -> frottemt et u.n=0 (vitesse)           !
!                  !    !     ! = 6   -> rugosite et u.n=0 (vitesse)           !
!                  !    !     ! = 9   -> entree/sortie libre (vitesse          !
!                  !    !     !  entrante eventuelle     bloquee               !
! itrifb(nfabor)   ! ia ! <-- ! indirection for boundary faces ordering        !
! itypfb(nfabor)   ! ia ! <-- ! boundary face types                            !
! izfppp(nfabor)   ! te ! <-- ! numero de zone de la face de bord              !
!                  !    !     !  pour le module phys. part.                    !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
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
use optcal
use cstphy
use cstnum
use entsor
use parall
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use ppcpfu
use cs_coal_incl
use mesh
use field

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal

integer          icodcl(nfabor,nvarcl)
integer          itrifb(nfabor), itypfb(nfabor)
integer          izfppp(nfabor)

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision rcodcl(nfabor,nvarcl,3)

! Local variables

character*80     name

integer          idebia, idebra
integer          ii, ifac, izone, mode, iel, ige, iok
integer          icha, iclapc, isol, icla
integer          ipbrom, icke, idecal, ipcvis
integer          nbrval, ioxy
integer          f_id, iaggas, keyvar

double precision qisqc, viscla, d2s3, uref2, rhomoy, dhy, xiturb
double precision ustar2, xkent, xeent, t1, t2, totcp , dmas
double precision h1    (nozppm) , h2   (nozppm,nclcpm)
double precision x2h20t(nozppm) , x20t (nozppm)
double precision qimpc (nozppm) , qcalc(nozppm)
double precision coefe (ngazem)
double precision xsolid(nsolim)
double precision f1mc  (ncharm) , f2mc (ncharm)
double precision wmh2o,wmco2,wmn2,wmo2

integer, dimension (:), allocatable :: iagecp

!===============================================================================
! 1. Initialization
!===============================================================================

ipbrom = ipprob(irom  )
ipcvis = ipproc(iviscl)

d2s3 = 2.d0/3.d0

call field_get_key_id("variable_id", keyvar)

call field_get_id('X_Age_Gas', f_id)

if (f_id.ne.-1) then
  call field_get_key_int(f_id, keyvar, iaggas)

  allocate (iagecp(nclacp))

  do icla = 1, nclacp
    write(name,'(a8,i2.2)')'X_Age_CP', icla
    call field_get_id(name, f_id)
    call field_get_key_int(f_id, keyvar, iagecp(icla))
  enddo
endif

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
!
if(irangp.ge.0) then
  call parimx(nozapm,iqimp )
  !==========
  call parimx(nozapm,ientat)
  !==========
  call parimx(nozapm,ientcp)
  !==========
  call parimx(nozapm,inmoxy)
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
  nbrval = nozppm*ncharm*ncpcmx
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
  qcalc(izone) = qcalc(izone) - propfb(ifac,ipbrom) *             &
                ( rcodcl(ifac,iu,1)*surfbo(1,ifac) +       &
                  rcodcl(ifac,iv,1)*surfbo(2,ifac) +       &
                  rcodcl(ifac,iw,1)*surfbo(3,ifac) )
enddo

if(irangp.ge.0) then
  call parrsm(nozapm,qcalc )
endif

do izone = 1, nozapm
  if ( iqimp(izone).eq.0 ) then
    qimpc(izone) = qcalc(izone)
  endif
enddo

! --- Correction des vitesses en norme a partir de la 2eme iter
!     car sinon on ne connait pas la masse volumique au bord

if ( ntcabs .gt. 1 ) then
  iok = 0
  do ii = 1, nzfppp
    izone = ilzppp(ii)
    if ( iqimp(izone).eq.1 ) then
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
    if ( iqimp(izone).eq.1 ) then
      qimpc(izone) = qimpat(izone)
      do icha = 1, ncharb
        qimpc(izone) = qimpc(izone) + qimpcp(izone,icha)
      enddo
      qisqc = qimpc(izone)/qcalc(izone)
      rcodcl(ifac,iu,1) = rcodcl(ifac,iu,1)*qisqc
      rcodcl(ifac,iv,1) = rcodcl(ifac,iv,1)*qisqc
      rcodcl(ifac,iw,1) = rcodcl(ifac,iw,1)*qisqc
    endif
  enddo
!
else
!
  do izone = 1, nozapm
    qimpc(izone) = qimpat(izone)
    do icha = 1, ncharb
      qimpc(izone) = qimpc(izone) + qimpcp(izone,icha)
    enddo
  enddo
!
endif
!
 2001 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : MODULE PHYSIQUES PARTICULIERES              ',/,&
'@    =========                        CHARBON PULVERISE      ',/,&
'@    PROBLEME DANS LES CONDITIONS AUX LIMITES                ',/,&
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
'@    - que le vecteur  RCODCL(IFAC,IU,1)                     ',/,&
'@                      RCODCL(IFAC,IV,1),                    ',/,&
'@                      RCODCL(IFAC,IW,1) qui determine       ',/,&
'@      la direction de la vitesse est non nul et n''est pas  ',/,&
'@      uniformement perpendiculaire aux face d''entree       ',/,&
'@    - que la surface de l''entree n''est pas nulle (ou que  ',/,&
'@      le nombre de faces de bord dans la zone est non nul)  ',/,&
'@    - que la masse volumique n''est pas nulle               ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

!===============================================================================
! 3. VERIFICATIONS
!        Somme des DISTributions CHarbon = 100% pour les zones IENTCP =1
!===============================================================================

iok = 0
do ii = 1, nzfppp
  izone = ilzppp(ii)
  if ( ientcp(izone).eq.1 ) then
    do icha = 1, ncharb
      totcp = 0.d0
      do iclapc = 1, nclpch(icha)
        totcp = totcp + distch(izone,icha,iclapc)
      enddo
      if(abs(totcp-100.d0).gt.epzero) then
        write(nfecra,2010)
        do iclapc = 1, nclpch(icha)
          write(nfecra,2011)izone,icha,iclapc,                    &
               distch(izone,icha,iclapc)
        enddo
        write(nfecra,2012)izone,ientcp(izone),icha,               &
             totcp,totcp-100.d0
        iok = iok + 1
      endif
    enddo
  endif
enddo

if(iok.ne.0) then
  call csexit (1)
  !==========
endif


 2010 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : MODULE PHYSIQUES PARTICULIERES              ',/,&
'@    =========                        CHARBON PULVERISE      ',/,&
'@    PROBLEME DANS LES CONDITIONS AUX LIMITES                ',/,&
'@                                                            ',/,&
'@        Zone    Charbon     Classe         Distch(%)        '  )
 2011 format(                                                           &
'@  ',I10   ,' ',I10   ,' ',I10   ,'    ',E14.5                  )
 2012 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : MODULE PHYSIQUES PARTICULIERES              ',/,&
'@    =========                        CHARBON PULVERISE      ',/,&
'@    PROBLEME DANS LES CONDITIONS AUX LIMITES                ',/,&
'@                                                            ',/,&
'@  On impose une entree charbon en IZONE = ', I10             ,/,&
'@    puisque               IENTCP(IZONE) = ', I10             ,/,&
'@  Or, sur cette zone, la somme des distributions par classe ',/,&
'@    en pourcentage pour le charbon ICHA = ', I10             ,/,&
'@    est differente de 100% : elle vaut TOTCP = ', E14.5      ,/,&
'@    avec                           TOTCP-100 = ', E14.5      ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier uscpcl.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

!===============================================================================
! 4.  REMPLISSAGE DU TABLEAU DES CONDITIONS LIMITES
!       ON BOUCLE SUR TOUTES LES FACES D'ENTREE
!                     =========================
!         ON DETERMINE LA FAMILLE ET SES PROPRIETES
!           ON IMPOSE LES CONDITIONS AUX LIMITES
!           POUR LA TURBULENCE

!===============================================================================
do ifac = 1, nfabor

  izone = izfppp(ifac)

!      ELEMENT ADJACENT A LA FACE DE BORD

  if ( itypfb(ifac).eq.ientre ) then

! ----  Traitement automatique de la turbulence

    if ( icalke(izone).ne.0 ) then

!       La turbulence est calculee par defaut si ICALKE different de 0
!          - soit a partir du diametre hydraulique, d'une vitesse
!            de reference adaptes a l'entree courante si ICALKE = 1
!          - soit a partir du diametre hydraulique, d'une vitesse
!            de reference et de l'intensite turvulente
!            adaptes a l'entree courante si ICALKE = 2

      uref2 = rcodcl(ifac,iu,1)**2                         &
            + rcodcl(ifac,iv,1)**2                         &
            + rcodcl(ifac,iw,1)**2
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

      if (itytur.eq.2) then

        rcodcl(ifac,ik,1)  = xkent
        rcodcl(ifac,iep,1) = xeent

      elseif (itytur.eq.3) then

        rcodcl(ifac,ir11,1) = d2s3*xkent
        rcodcl(ifac,ir22,1) = d2s3*xkent
        rcodcl(ifac,ir33,1) = d2s3*xkent
        rcodcl(ifac,ir12,1) = 0.d0
        rcodcl(ifac,ir13,1) = 0.d0
        rcodcl(ifac,ir23,1) = 0.d0
        rcodcl(ifac,iep,1)  = xeent

      elseif (iturb.eq.50) then

        rcodcl(ifac,ik,1)   = xkent
        rcodcl(ifac,iep,1)  = xeent
        rcodcl(ifac,iphi,1) = d2s3
        rcodcl(ifac,ifb,1)  = 0.d0

      elseif (iturb.eq.60) then

        rcodcl(ifac,ik,1)   = xkent
        rcodcl(ifac,iomg,1) = xeent/cmu/xkent

      endif

    endif

  endif

enddo

!===============================================================================
! 5.  REMPLISSAGE DU TABLEAU DES CONDITIONS LIMITES
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
  if ( ientat(izone).eq.1 .or. ientcp(izone).eq.1) then

    x20t  (izone) = zero
    x2h20t(izone) = zero

    idecal = 0

    do icha = 1, ncharb

      do iclapc = 1, nclpch(icha)

        icla = iclapc + idecal
! ------ Calcul de X2 total par zone
!         Petite retouche au cas ou l'entree est fermee
        if(abs(qimpc(izone)).lt.epzero) then
          x20(izone,icla) = 0.d0
        else
          x20(izone,icla) = qimpcp(izone,icha)/qimpc(izone)       &
                          * distch(izone,icha,iclapc)*1.d-2
        endif
        x20t(izone)     = x20t(izone) +  x20(izone,icla)
! ------ Calcul de H2 de la classe ICLA
        do isol = 1, nsolim
          xsolid(isol) = zero
        enddo
        if ( ientcp(izone).eq.1 ) then
          t2  = timpcp(izone,icha)
          xsolid(ich(icha)) = 1.d0-xashch(icha)
          xsolid(ick(icha)) = zero
          xsolid(iash(icha)) = xashch(icha)

!------- Prise en compte de l'humidite
          if ( ippmod(iccoal) .eq. 1 ) then
            xsolid(ich(icha)) = xsolid(ich(icha))-xwatch(icha)
            xsolid(iwat(icha)) = xwatch(icha)
          else
            xsolid(iwat(icha)) = 0.d0
          endif

        else
          t2  = timpat(izone)

          xsolid(ich(icha))  = (1.d0-xashch(icha)-xwatch(icha))
          xsolid(ick(icha))  = 0.d0
          xsolid(iash(icha)) = xashch(icha)
          xsolid(iwat(icha)) = xwatch(icha)

        endif
        mode = -1
        t1 = t2
        call cs_coal_htconvers2(mode,icla,h2(izone,icla),xsolid,t2,t1)
        !======================

        x2h20t(izone) = x2h20t(izone)+x20(izone,icla)*h2(izone,icla)

      enddo

      idecal = idecal + nclpch(icha)

    enddo

! ------ Calcul de H1(IZONE)
    do ige = 1, ngazem
      coefe(ige) = zero
    enddo

    ioxy = inmoxy(izone)
    dmas = wmole(io2) *oxyo2(ioxy) +wmole(in2) *oxyn2(ioxy)       &
          +wmole(ih2o)*oxyh2o(ioxy)+wmole(ico2)*oxyco2(ioxy)

    coefe(io2)  = wmole(io2 )*oxyo2(ioxy )/dmas
    coefe(ih2o) = wmole(ih2o)*oxyh2o(ioxy)/dmas
    coefe(ico2) = wmole(ico2)*oxyco2(ioxy)/dmas
    coefe(in2)  = wmole(in2 )*oxyn2(ioxy )/dmas

    do icha = 1, ncharm
      f1mc(icha) = zero
      f2mc(icha) = zero
    enddo
    t1   = timpat(izone)
    mode = -1
    call cs_coal_htconvers1(mode,h1(izone),coefe,f1mc,f2mc,t1)
   !=======================

  endif

enddo

do ifac = 1, nfabor

  izone = izfppp(ifac)

!      ELEMENT ADJACENT A LA FACE DE BORD

  if ( itypfb(ifac).eq.ientre ) then

! ----  Traitement automatique des scalaires physiques particulieres

    idecal = 0

    do icha = 1, ncharb

      do iclapc = 1, nclpch(icha)

        icla = iclapc + idecal
! ------ CL pour Xch de la classe ICLA
        rcodcl(ifac,isca(ixch(icla)),1) = x20(izone,icla)         &
                                        * (1.d0-xashch(icha))
!             Prise en compte de l'humidite
        if ( ippmod(iccoal) .eq. 1 ) then
          rcodcl(ifac,isca(ixch(icla)),1) = x20(izone,icla)       &
                                          *(1.d0-xashch(icha)     &
                                                -xwatch(icha))
        endif
! ------ CL pour Xck de la classe ICLA
        rcodcl(ifac,isca(ixck(icla)),1) = 0.d0
! ------ CL pour Np de la classe ICLA
        rcodcl(ifac,isca(inp(icla)),1) = x20(izone,icla)          &
                                        / xmp0(icla)
! ------ CL pour Xwater de la classe ICLA
        if ( ippmod(iccoal) .eq. 1 ) then
          rcodcl(ifac,isca(ixwt(icla)),1) = x20(izone,icla)       &
                                           *xwatch(icha)
        endif
! ------ CL pour H2 de la classe ICLA
        rcodcl(ifac,isca(ih2(icla)),1) = x20(izone,icla)          &
                                        *h2(izone,icla)
        if (i_coal_drift.eq.1) then
          rcodcl(ifac, iagecp(icla), 1) = zero
        endif
      enddo

      idecal = idecal + nclpch(icha)

! ------ CL pour X1F1M et X1F2M du charbon ICHA
      rcodcl(ifac,isca(if1m(icha)),1) = zero
      rcodcl(ifac,isca(if2m(icha)),1) = zero

    enddo
    if (i_coal_drift.eq.1) then
      rcodcl(ifac, iaggas, 1) = zero
    endif
! ------ CL pour HM
    rcodcl(ifac,isca(ihm),1) = (1.d0-x20t(izone))*h1(izone)    &
                              +x2h20t(izone)
! ------ CL pour X1.F4M (Oxyd 2)
    if ( noxyd .ge. 2 ) then
      if ( inmoxy(izone) .eq. 2 ) then
        rcodcl(ifac,isca(if4m),1)   = (1.d0-x20t(izone))
      else
        rcodcl(ifac,isca(if4m),1)   = zero
      endif
    endif
! ------ CL pour X1.F5M (Oxyd3)
    if ( noxyd .eq. 3 ) then
      if ( inmoxy(izone) .eq. 3 ) then
        rcodcl(ifac,isca(if5m),1)   = (1.d0-x20t(izone))
      else
        rcodcl(ifac,isca(if5m),1)   = zero
      endif
    endif
! ------ CL pour X1.F6M (Water)
    if ( ippmod(iccoal) .ge. 1 ) then
      rcodcl(ifac,isca(if6m),1) = zero
    endif
! ------ CL pour X1.F7M_O2
    rcodcl(ifac,isca(if7m),1)   = zero
! ------ CL pour X1.FM8_CO2
    if ( ihtco2 .eq. 1 ) then
      rcodcl(ifac,isca(if8m),1) = zero
    endif
! ------ CL pour X1.FM9_H2O
    if ( ihth2o .eq. 1 ) then
      rcodcl(ifac,isca(if9m),1) = zero
    endif
! ------ CL pour X1.Variance
    rcodcl(ifac,isca(ifvp2m),1) = zero

! ------ CL pour X1.YCO2
    if ( ieqco2 .eq. 1 ) then
      ioxy =  inmoxy(izone)
      wmo2   = wmole(io2)
      wmco2  = wmole(ico2)
      wmh2o  = wmole(ih2o)
      wmn2   = wmole(in2)
      dmas = ( oxyo2 (ioxy)*wmo2 +oxyn2 (ioxy)*wmn2               &
              +oxyh2o(ioxy)*wmh2o+oxyco2(ioxy)*wmco2 )
      xco2 = oxyco2(ioxy)*wmco2/dmas
      rcodcl(ifac,isca(iyco2),1)   = xco2*(1.d0-x20t(izone))
    endif
! ------ CL pour X1.HCN, X1.NO, Taire
    if( ieqnox .eq. 1 ) then
      rcodcl(ifac,isca(iyhcn ),1)  = zero
      rcodcl(ifac,isca(iyno  ),1)  = zero
      rcodcl(ifac,isca(iynh3 ),1)  = zero
      rcodcl(ifac,isca(ihox  ),1)  = (1.d0-x20t(izone))*h1(izone)
    endif

  endif

enddo

! Free memory
if (allocated(iagecp)) deallocate(iagecp)

!--------
! Formats
!--------

!----
! End
!----

return
end subroutine
