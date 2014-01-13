!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2014 EDF S.A.
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

subroutine cs_fuel_bcond &
!=======================

 ( itypfb , izfppp ,                                              &
   propce ,                                                       &
   rcodcl )

!===============================================================================
! FONCTION :
! --------
!    CONDITIONS AUX LIMITES AUTOMATIQUES
!           COMBUSTION FUEL
!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! itypfb(nfabor)   ! ia ! <-- ! boundary face types                            !
! izfppp(nfabor)   ! te ! <-- ! numero de zone de la face de bord              !
!                  !    !     !  pour le module phys. part.                    !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! rcodcl           ! tr ! --> ! valeur des conditions aux limites              !
!  (nfabor,nvarcl) !    !     !  aux faces de bord                             !
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
use cs_fuel_incl
use ppincl
use ppcpfu
use mesh
use field
!===============================================================================

implicit none

! Arguments

integer          itypfb(nfabor)
integer          izfppp(nfabor)

double precision propce(ncelet,*)
double precision rcodcl(nfabor,nvarcl,3)

! Local variables

integer          ii, ifac, izone, mode, iel, ige, iok
integer          icla , ioxy
integer          icke, ipcvis
integer          nbrval
double precision qisqc, viscla, d2s3, uref2, rhomoy, dhy, xiturb
double precision xkent, xeent, t1, t2, ustar2
double precision h1(nozppm) , h2(nozppm)
double precision x20t(nozppm)
double precision xmg0(nozppm,nclcpm)
double precision x2h20t(nozppm)
double precision qimpc(nozppm) , qcalc(nozppm)
double precision coefe(ngazem)
double precision xsolid(2)
double precision hlf, totfu
double precision dmas
double precision, dimension(:), pointer ::  brom

!===============================================================================

!===============================================================================
! 1.  INITIALISATIONS
!===============================================================================
!
call field_get_val_s(ibrom, brom)
ipcvis = ipproc(iviscl)
!
d2s3   = 2.d0 / 3.d0
!
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
  call parimx(nozapm,ientfl)
  !==========
  call parimx(nozapm,inmoxy)
  !==========
  call parrmx(nozapm,qimpat)
  !==========
  call parrmx(nozapm,timpat)
  !==========
  nbrval = nozppm
  call parrmx(nbrval,qimpfl)
  !==========
  nbrval = nozppm
  call parrmx(nbrval,timpfl)
  !==========
  nbrval = nozppm*nclcpm
  call parrmx(nbrval,distfu)
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
  h1(izone)    = 0.d0
enddo
do ifac = 1, nfabor
  izone = izfppp(ifac)
  qcalc(izone) = qcalc(izone) - brom(ifac) *             &
      ( rcodcl(ifac,iu,1)*surfbo(1,ifac) +                 &
        rcodcl(ifac,iv,1)*surfbo(2,ifac) +                 &
        rcodcl(ifac,iw,1)*surfbo(3,ifac) )
enddo

if(irangp .ge. 0) then
  call parrsm(nozapm,qcalc )
endif

do izone = 1, nozapm
  if ( iqimp(izone) .eq. 0 ) then
    qimpc(izone) = qcalc(izone)
  endif
enddo

if ( ntcabs .gt. 1 ) then
!
! --- Correction des vitesses en norme :  on ne le fait qu'a la
!     2eme iteration car pour la 1ere la masse vol au bord n'est
!     pas encore connue
!
  iok = 0
  do ii = 1, nzfppp
    izone = ilzppp(ii)
    if ( iqimp(izone) .eq. 1 ) then
      if(abs(qcalc(izone)) .lt. epzero) then
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
    if ( iqimp(izone) .eq. 1 ) then
      qimpc(izone) = qimpat(izone) + qimpfl(izone)
      qisqc = qimpc(izone) / qcalc(izone)
      rcodcl(ifac,iu,1) = rcodcl(ifac,iu,1)*qisqc
      rcodcl(ifac,iv,1) = rcodcl(ifac,iv,1)*qisqc
      rcodcl(ifac,iw,1) = rcodcl(ifac,iw,1)*qisqc
    endif
  enddo
!
else
!
  do izone = 1, nozapm
    if ( iqimp(izone) .eq. 1 ) then
      qimpc(izone) = qimpat(izone) + qimpfl(izone)
    endif
  enddo
!
endif
!
 2001 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : MODULE PHYSIQUES PARTICULIERES              ',/,&
'@    =========                       FUEL                    ',/,&
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
'@  Verifier user_fuel_bconds, et en particulier  ',/,&
'@    - que le vecteur  RCODCL(IFAC,IU,1),                    ',/,&
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
!        Somme des DISTributions FUel = 100% pour les zones IENTFL =1
!===============================================================================

iok = 0
do ii = 1, nzfppp
  izone = ilzppp(ii)
  if ( ientfl(izone).eq.1 ) then
    totfu = 0.d0
    do icla = 1, nclafu
      totfu = totfu + distfu(izone,icla)
    enddo
    if(abs(totfu-100.d0).gt.epzero) then
      write(nfecra,2010)
      do icla = 1, nclafu
        write(nfecra,2011)izone,icla,                             &
               distfu(izone,icla)
      enddo
      write(nfecra,2012)izone,ientfl(izone),                      &
             totfu,totfu-100.d0
      iok = iok + 1
    endif
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
'@    =========                       FUEL                    ',/,&
'@    PROBLEME DANS LES CONDITIONS AUX LIMITES                ',/,&
'@                                                            ',/,&
'@        Zone    Classe         Distfu(%)                    '  )
 2011 format(                                                           &
'@  ',I10   ,' ',I10   ,'    ',E14.5                             )
 2012 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : MODULE PHYSIQUES PARTICULIERES              ',/,&
'@    =========                        FUEL                   ',/,&
'@    PROBLEME DANS LES CONDITIONS AUX LIMITES                ',/,&
'@                                                            ',/,&
'@  On impose une entree fuel en IZONE = ', I10                ,/,&
'@    puisque               IENTFL(IZONE) = ', I10            ,/, &
'@  Or, sur cette zone, la somme des distributions            ',/,&
'@    en pourcentage pour le fuel IFOL = ', I10                ,/,&
'@    est differente de 100% : elle vaut TOTFOL = ', E14.5     ,/,&
'@    avec                           TOTFOL-100 = ', E14.5     ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier user_fuel_bconds.                    ',/,&
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
      rhomoy = brom(ifac)
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
!         IENTAT = 1 ou IENTFL = 1
  if ( ientat(izone).eq.1 .or. ientfl(izone).eq.1) then

    x20t  (izone) = zero
    x2h20t(izone) = zero

    do icla = 1, nclafu

! ------ Calcul de X2 total par zone
!         Petite retouche au cas ou l'entree est fermee
      if(abs(qimpc(izone)).le.epzero) then
        x20(izone,icla) = 0.d0
      else
        x20(izone,icla) = qimpfl(izone)/qimpc(izone)              &
                         *distfu(izone,icla)*1.d-2
      endif
      x20t(izone)     = x20t(izone) +  x20(izone,icla)
    enddo
! ------ Calcul de H2 , XMG0
    if ( ientfl(izone) .eq. 1 ) then
      t2        = timpfl(izone)
      xsolid(1) = 1.d0-fkc
      xsolid(2) = fkc
      mode      = -1
      call cs_fuel_htconvers2 (mode, h2(izone) , xsolid , t2)
!     =======================
!
      do icla = 1, nclafu
        xmg0(izone,icla) = pi/6.d0*(dinifl(icla)**3)*rho0fl
      enddo
    else
      h2(izone) = zero
      do icla = 1, nclafu
        xmg0(izone,icla) = 1.d0
      enddo
    endif
    x2h20t(izone) = x20t(izone)*h2(izone)


! ------ Calcul de H1(IZONE)
    do ige = 1, ngazem
      coefe(ige) = zero
    enddo
!
    ioxy = inmoxy(izone)
    dmas = wmole(io2) *oxyo2(ioxy) +wmole(in2) *oxyn2(ioxy)    &
          +wmole(ih2o)*oxyh2o(ioxy)+wmole(ico2)*oxyco2(ioxy)
!
    coefe(io2)  = wmole(io2 )*oxyo2(ioxy )/dmas
    coefe(ih2o) = wmole(ih2o)*oxyh2o(ioxy)/dmas
    coefe(ico2) = wmole(ico2)*oxyco2(ioxy)/dmas
    coefe(in2)  = wmole(in2 )*oxyn2(ioxy )/dmas
!
    hlf = zero
    t1   = timpat(izone)
    mode = -1
    call cs_fuel_htconvers1 (mode, h1(izone) , coefe , t1)
!   =======================
!
  endif
enddo
!
!
do ifac = 1, nfabor

  izone = izfppp(ifac)

!      ELEMENT ADJACENT A LA FACE DE BORD

  if ( itypfb(ifac).eq.ientre ) then

! ----  Traitement automatique des scalaires physiques particulieres

    do icla = 1, nclafu
! ------ CL pour Xfol
      rcodcl(ifac,isca(iyfol(icla)),1) = x20(izone,icla)
! ------ CL pour Ng
      rcodcl(ifac,isca(ing(icla)),1) = x20(izone,icla)            &
                                      /xmg0(izone,icla)
! ------ CL pour X2HLF
      rcodcl(ifac,isca(ih2(icla)),1) = x20(izone,icla)*h2(izone)
    enddo
! ------ CL pour X1.FVAP
    rcodcl(ifac,isca(ifvap),1) = zero
! ------ CL pour X1.F7M
    rcodcl(ifac,isca(if7m),1) = zero
! ------ CL pour X1.Variance
    rcodcl(ifac,isca(ifvp2m),1)   = zero
! ------ CL pour HM
    rcodcl(ifac,isca(iscalt),1) = (1.d0-x20t(izone))*h1(izone)+x2h20t(izone)
!
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

! ------ CL pour X1.YCO2
    if ( ieqco2 .ge. 1 ) then
      rcodcl(ifac,isca(iyco2),1)   = zero
    endif

! ------ CL pour X1.HCN et X1.NO
    if ( ieqnox .eq. 1 ) then
      rcodcl(ifac,isca(iyhcn),1)   = zero
      rcodcl(ifac,isca(iyno ),1)   = zero
      rcodcl(ifac,isca(ihox ),1)   = (1.d0-x20t(izone))*h1(izone)
    endif

  endif

enddo

!----
! Formats
!----


!----
! End
!----


return
end subroutine
