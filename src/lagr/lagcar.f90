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

subroutine lagcar &
!================

 ( nvar   , nscal  ,                                              &
   nbpmax ,                                                       &
   nvlsta , iprev  ,                                              &
   dt     ,                                                       &
   taup   , tlag   ,                                              &
   piil   , bx     , tempct , statis ,                            &
   gradpr , gradvf , energi , dissip , romp   )

!===============================================================================
! FONCTION :
! ----------

!   SOUS-PROGRAMME DU MODULE LAGRANGIEN :
!   -------------------------------------

!    CALCUL DES CARACTERISTIQUES DES PARTICULES : Tp, TL et PI

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! nbpmax           ! e  ! <-- ! nombre max de particulies autorise             !
! nvlsta           ! e  ! <-- ! nombre de var statistiques lagrangien          !
! iprev            ! e  ! <-- ! time step indicator for fields                 !
!                  !    !     !   0: use fields at current time step           !
!                  !    !     !   1: use fields at previous time step          !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! taup(nbpmax)     ! tr ! --> ! temps caracteristiques dynamique               !
! tlag(nbpmax)     ! tr ! --> ! temps caracteristiques fluide                  !
! piil(nbpmax,3    ! tr ! --> ! terme dans l'integration des eds up            !
! bx(nbpmax,3,2    ! tr ! --> ! caracteristiques de la turbulence              !
! tempct           ! tr ! --> ! temps caracteristique thermique                !
!   (nbpmax,2)     !    !     !                                                !
! statis(ncelet    ! tr ! <-- ! cumul des statistiques volumiques              !
!    nvlsta)       !    !     !                                                !
! gradpr           ! tr ! <-- ! gradient de pression                           !
!   (3,ncelet)     !    !     !                                                !
! gradvf           ! tr ! <-- ! gradient de la vitesse du fluide               !
!   (3,3,ncelet)   !    !     !                                                !
! energi(ncelet    ! tr ! --- ! tableau de travail                             !
! dissip(ncelet    ! tr ! --- ! tableau de travail                             !
! romp(nbpmax)     ! tr ! --- ! tableau de travail                             !
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
use cstnum
use cstphy
use optcal
use entsor
use lagpar
use lagran
use ppppar
use ppthch
use ppincl
use mesh
use field

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal
integer          nbpmax
integer          nvlsta
integer          iprev

double precision dt(ncelet)
double precision taup(nbpmax) , tlag(nbpmax,3)
double precision piil(nbpmax,3) , bx(nbpmax,3,2)
double precision tempct(nbpmax,2)
double precision statis(ncelet,nvlsta)
double precision gradpr(3,ncelet) , gradvf(3,3,ncelet)
double precision energi(ncelet) , dissip(ncelet), romp(nbpmax)

! Local variables

integer          iel , ip , id , iscath

double precision cd1 , cd2 , rec , cl , c0 , cb , cbcb
double precision upart , vpart , wpart
double precision uflui , vflui , wflui
double precision uvwdif , tl , uvwr
double precision rep , d2 , d3 , fdr , d1s3 , d3s444 , d6spi
double precision bb1 , bb2 , bb3 , ktil , bx1 , bx2 , bx3
double precision vpmx , vpmy , vpmz
double precision r11 , r22 , r33
double precision xnul , rom , prt , fnus , xrkl , xcp

double precision, dimension(:), pointer :: cromf
double precision, dimension(:,:), pointer :: vel
double precision, dimension(:), pointer :: cvar_k, cvar_ep, cvar_omg
double precision, dimension(:), pointer :: cvar_r11, cvar_r22, cvar_r33
double precision, dimension(:), pointer :: viscl, visls, cpro_cp

!===============================================================================

! Map field arrays
if (iprev.eq.0) then
  call field_get_val_v(ivarfl(iu), vel)

  if (itytur.eq.2 .or. iturb.eq.50) then
    call field_get_val_s(ivarfl(ik), cvar_k)
    call field_get_val_s(ivarfl(iep), cvar_ep)
  else if (itytur.eq.3) then
    call field_get_val_s(ivarfl(ir11), cvar_r11)
    call field_get_val_s(ivarfl(ir22), cvar_r22)
    call field_get_val_s(ivarfl(ir33), cvar_r33)
    call field_get_val_s(ivarfl(iep), cvar_ep)
  else if (iturb.eq.60) then
    call field_get_val_s(ivarfl(ik), cvar_k)
    call field_get_val_s(ivarfl(iomg), cvar_omg)
  endif

else if (iprev.eq.1) then
  call field_get_val_prev_v(ivarfl(iu), vel)

  if (itytur.eq.2 .or. iturb.eq.50) then
    call field_get_val_prev_s(ivarfl(ik), cvar_k)
    call field_get_val_prev_s(ivarfl(iep), cvar_ep)
  else if (itytur.eq.3) then
    call field_get_val_prev_s(ivarfl(ir11), cvar_r11)
    call field_get_val_prev_s(ivarfl(ir22), cvar_r22)
    call field_get_val_prev_s(ivarfl(ir33), cvar_r33)
    call field_get_val_prev_s(ivarfl(iep), cvar_ep)
  else if (iturb.eq.60) then
    call field_get_val_prev_s(ivarfl(ik), cvar_k)
    call field_get_val_prev_s(ivarfl(iomg), cvar_omg)
  endif
endif

call field_get_val_s(iprpfl(iviscl), viscl)
if (icp.gt.0) call field_get_val_s(iprpfl(icp), cpro_cp)

!===============================================================================
! 0.  GESTION MEMOIRE
!===============================================================================

! Initialize variables to avoid compiler warnings

bb1 = 0.d0
bb2 = 0.d0
bb3 = 0.d0
ktil = 0.d0

!===============================================================================
! 1. INITIALISATIONS
!===============================================================================

cd1  = 0.15d0
cd2  = 0.687d0
rec  = 1000.d0
c0   = 2.1d0
cl   = 1.d0 / (0.5d0 + (3.d0/4.d0)*c0 )
cb   = 0.8d0
cbcb = 0.64d0

d6spi = 6.d0 / pi
d1s3 = 1.d0 / 3.d0
d3s444 = 0.44d0 * 3.d0 / 4.d0

! Pointeur sur la masse volumique en fonction de l'ecoulement

if (ippmod(iccoal).ge.0 .or. ippmod(icfuel).ge.0) then
  call field_get_val_s(iprpfl(ipproc(irom1)), cromf)
else
  call field_get_val_s(icrom, cromf)
endif

! Calcul de la masse volumique

do ip = 1,nbpart
  if ( ipepa(jisor,ip).gt.0 ) then
    d3 = eptp(jdp,ip) * eptp(jdp,ip) * eptp(jdp,ip)
    romp(ip) = eptp(jmp,ip) * d6spi / d3
  endif
enddo

!===============================================================================
! 2. CALCUL DE Tp ET DE Tc SI THERMIQUE
!===============================================================================

do ip = 1,nbpart

  if ( ipepa(jisor,ip) .gt.0 ) then

    iel = ipepa(jisor,ip)

    rom  = cromf(iel)
    xnul = viscl(iel) / rom

    uvwr = sqrt( ( eptp(juf,ip) -eptp(jup,ip) )*                  &
                 ( eptp(juf,ip) -eptp(jup,ip) )                   &
               + ( eptp(jvf,ip) -eptp(jvp,ip) )*                  &
                 ( eptp(jvf,ip) -eptp(jvp,ip) )                   &
               + ( eptp(jwf,ip) -eptp(jwp,ip) )*                  &
                 ( eptp(jwf,ip) -eptp(jwp,ip) )  )

!--->  CALCUL DU REYNOLDS LOCAL

    rep  = uvwr * eptp(jdp,ip) / xnul

!--->  CALCUL DU COEFFICIENT DE TRAINEE

    d2 = eptp(jdp,ip) * eptp(jdp,ip)

    if (rep.le.rec) then
      fdr = 18.d0 * xnul * (1.d0 + cd1 * rep**cd2) / d2
    else
      fdr = d3s444 * uvwr / eptp(jdp,ip)
    endif

!--->  CALCUL DE Tp

    taup(ip) = romp(ip) / rom / fdr

!--->  CALCUL UTILISATEUR DE Tp

    call uslatp                                                   &
    !==========
     ( nvar   , nscal  ,                                          &
       nbpmax ,                                                   &
       ip     ,                                                   &
       rep    , uvwr   , rom    , romp(ip) , xnul , taup(ip) ,    &
       dt     )

!--->  CALCUL DE Tc

    if ( (iphyla.eq.1 .and. itpvar.eq.1) .or.                     &
         (iphyla.eq.2)                        ) then

!     CP fluide

      if (icp.gt.0) then
        xcp = cpro_cp(1)
      else
        xcp = cp0
      endif

!     CALCUL DU NUSSELT LOCAL

      iscath = iscalt

! a priori en combustion gaz ou CP, la diffusvite est toujours constante

      if (ippmod(icoebu).eq.0 .or. ippmod(icoebu).eq.2) then
        xrkl = diftl0 / rom
      else if (ivisls(iscath).ge.1) then
        call field_get_val_s(iprpfl(ivisls(iscath)), visls)
        xrkl = visls(iel) / rom
      else
        xrkl = visls0(iscath) / rom
      endif

      prt  = xnul / xrkl
      fnus = 2.d0 + 0.55d0 * rep**0.5d0 * prt**(d1s3)

! Calcul du temps caracteristique thermique Tc

      tempct(ip,1) = d2 * romp(ip) * eptp(jcp,ip)                 &
                   / ( fnus * 6.d0 * rom * xcp * xrkl )

!--->  CALCUL UTILISATEUR DE Tc

    call uslatc                                                   &
    !==========
     ( nvar   , nscal  ,                                          &
       nbpmax ,                                                   &
       ip     ,                                                   &
       rep    , uvwr   , rom    , romp(ip) , xnul ,               &
       xcp    , xrkl   , tempct(ip,1) ,                           &
       dt     )

! Terme source implicite pour le couplage retour thermique

      tempct(ip,2) = fnus * pi * eptp(jdp,ip) * xrkl * rom

    endif

  endif

enddo


!===============================================================================
! 3. CALCUL DE TL
!===============================================================================

!--> Calcul de l'energie turbulente et de la dissipation
!      en fonction du modele de turbulence

if (idistu.eq.1) then

  if (itytur.eq.2 .or. iturb.eq.50) then
    do iel = 1,ncel
      energi(iel) = cvar_k(iel)
      dissip(iel) = cvar_ep(iel)
    enddo
  else if (itytur.eq.3) then
    do iel = 1,ncel
      energi(iel) = 0.5d0*( cvar_r11(iel)                  &
                           +cvar_r22(iel)                  &
                           +cvar_r33(iel) )
      dissip(iel) = cvar_ep(iel)
    enddo
  else if (iturb.eq.60) then
    do iel = 1,ncel
      energi(iel) = cvar_k(iel)
      dissip(iel) = cmu*energi(iel)*cvar_omg(iel)
    enddo
  else
    write(nfecra,2000) iilagr, idistu, iturb
    call csexit (1)
!              ======
  endif

!--> Calcul de TL et BX

  do ip = 1,nbpart

    if (ipepa(jisor,ip).gt.0) then

      iel = ipepa(jisor,ip)


      if (dissip(iel).gt.0.d0 .and.                               &
          energi(iel).gt.0.d0              ) then

      tl = cl * energi(iel) / dissip(iel)
      tl = max(tl,epzero)

      upart = eptp(jup,ip)
      vpart = eptp(jvp,ip)
      wpart = eptp(jwp,ip)
      uflui = eptp(juf,ip)
      vflui = eptp(jvf,ip)
      wflui = eptp(jwf,ip)

      if (modcpl.gt.0 .and. iplas.gt.modcpl) then
        if (statis(iel,ilpd).gt.seuil) then
          upart = statis(iel,ilvx) / statis(iel,ilpd)
          vpart = statis(iel,ilvy) / statis(iel,ilpd)
          wpart = statis(iel,ilvz) / statis(iel,ilpd)
          uflui = vel(1,iel)
          vflui = vel(2,iel)
          wflui = vel(3,iel)
        endif
      endif

      uvwdif = (uflui-upart) * (uflui-upart)                      &
             + (vflui-vpart) * (vflui-vpart)                      &
             + (wflui-wpart) * (wflui-wpart)

      uvwdif = (3.d0 * uvwdif) / (2.d0 *energi(iel))

      if (modcpl.gt.0 .and. iplas.gt.modcpl) then

        if (idirla.eq.1) then
          bb1 = sqrt( 1.d0 + cbcb *uvwdif )
          tlag(ip,1)= tl / bb1
          bb2  = sqrt( 1.d0 + 4.d0 *cbcb *uvwdif )
          tlag(ip,2)= tl / bb2
          bb3  = sqrt( 1.d0 + 4.d0 *cbcb *uvwdif )
          tlag(ip,3)= tl / bb3

        else if (idirla.eq.2) then
          bb1  = sqrt( 1.d0 + 4.d0 *cbcb *uvwdif )
          tlag(ip,1)= tl / bb1
          bb2  = sqrt( 1.d0 + cbcb *uvwdif )
          tlag(ip,2)= tl / bb2
          bb3  = sqrt( 1.d0 + 4.d0 *cbcb *uvwdif )
          tlag(ip,3)= tl / bb3

        else if (idirla.eq.3) then
          bb1  = sqrt( 1.d0 + 4.d0 *cbcb *uvwdif )
          tlag(ip,1)= tl / bb1
          bb2  = sqrt( 1.d0 + 4.d0 *cbcb *uvwdif )
          tlag(ip,2)= tl / bb2
          bb3  = sqrt( 1.d0 + cbcb *uvwdif )
          tlag(ip,3)= tl / bb3
        else
          write(nfecra,1000) idirla
          call csexit(1)
!                    ======
        endif

        if (itytur.eq.3) then
          r11 = cvar_r11(iel)
          r22 = cvar_r22(iel)
          r33 = cvar_r33(iel)
          ktil = 3.d0 * ( r11*bb1 + r22*bb2 + r33*bb3  )          &
               / (2.d0 * (bb1+bb2+bb3) )
        else if (itytur.eq.2 .or. iturb.eq.50       &
                                    .or. iturb.eq.60) then
          ktil = energi(iel)
        endif

        bx1 = dissip(iel) * ( (c0*bb1*ktil/energi(iel))           &
                 +(2.d0 *(bb1*ktil/energi(iel) -1.d0)/3.d0) )
        bx2 = dissip(iel) * ( (c0*bb2*ktil/energi(iel))           &
                 +(2.d0 *(bb2*ktil/energi(iel) -1.d0)/3.d0) )
        bx3 = dissip(iel) * ( (c0*bb3*ktil/energi(iel))           &
                 +(2.d0 *(bb3*ktil/energi(iel) -1.d0)/3.d0) )

        if (bx1.gt.0.d0) then
          bx(ip,1,nor) = sqrt(bx1)
        else
          bx(ip,1,nor) = 0.d0
        endif

        if (bx2.gt.0.d0) then
          bx(ip,2,nor) = sqrt(bx2)
        else
          bx(ip,2,nor) = 0.d0
        endif

        if (bx3.gt.0.d0) then
          bx(ip,3,nor) = sqrt(bx3)
        else
          bx(ip,3,nor) = 0.d0
        endif

      else

        tlag(ip,1) = tl
        tlag(ip,2) = tl
        tlag(ip,3) = tl

        if (idiffl.eq.0) then
          uvwdif = sqrt(uvwdif)
          tlag(ip,1) = tl/(1.d0 + cb*uvwdif)
          tlag(ip,2) = tlag(ip,1)
          tlag(ip,3) = tlag(ip,1)
        endif

        bx(ip,1,nor) = sqrt(c0*dissip(iel))
        bx(ip,2,nor) = bx(ip,1,nor)
        bx(ip,3,nor) = bx(ip,1,nor)
      endif

      else

        tlag(ip,1) = epzero
        tlag(ip,2) = epzero
        tlag(ip,3) = epzero
        bx(ip,1,nor) = zero
        bx(ip,2,nor) = zero
        bx(ip,3,nor) = zero

      endif

    endif

  enddo

else

  do ip = 1,nbpart

    if ( ipepa(jisor,ip) .gt.0 ) then
      tlag(ip,1) = epzero
      tlag(ip,2) = epzero
      tlag(ip,3) = epzero
      bx(ip,1,nor) = zero
      bx(ip,2,nor) = zero
      bx(ip,3,nor) = zero
    endif

  enddo

endif

!===============================================================================
! 4. CALCUL DE PII
!===============================================================================

do id = 1,3

  do ip = 1,nbpart

    if (ipepa(jisor,ip).gt.0) then

!--->   Calcul de II = ( -grad(P)/Rom(f)+grad(<Vf>)*(<Up>-<Uf>) )
!       ou
!       Calcul de II = ( -grad(P)/Rom(f) )

      iel = ipepa(jisor,ip)

      piil(ip,id) = gradpr(id,iel)

      if (modcpl.gt.0 .and. iplas.gt.modcpl) then
        if (statis(iel,ilpd) .gt. seuil) then
          vpmx = statis(iel,ilvx) / statis(iel,ilpd)
          vpmy = statis(iel,ilvy) / statis(iel,ilpd)
          vpmz = statis(iel,ilvz) / statis(iel,ilpd)

          uflui = vel(1,iel)
          vflui = vel(2,iel)
          wflui = vel(3,iel)

          piil(ip,id) = gradpr(id,iel)                            &
                       +gradvf(1,id,iel) * (vpmx-uflui)           &
                       +gradvf(2,id,iel) * (vpmy-vflui)           &
                       +gradvf(3,id,iel) * (vpmz-wflui)
        endif
      endif

!--->  Terme purement explicite : probleme avec petit diametre
!      ne pas effacer svp           !

!            if (iilagr.eq.2 .and. istala.ge.1) then

!              if (statis(iel,ilpd).gt.seuil) then

!                rom   = cromf(iel)

!                ff = romp(ip) / rom
!     &             *( statis(iel,ilfv) / (dble(npst)*volume(iel)) )
!     &             *( eptp(juf+(id-1),ip) - eptp(jup+(id-1),ip))
!     &                   /taup(ip)

!                piil(ip,id) = piil(ip,id) - ff

!              endif

!            endif

    endif

  enddo

enddo

!==============================================================================

!--------
! FORMATS
!--------

 1000   format(                                                         &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========                                               ',/,&
'@    LE CHOIX DE LA DIRECTION DU MODELE COMPLET              ',/,&
'@       A UNE VALEUR NON PERMISE (LAGCAR).                   ',/,&
'@                                                            ',/,&
'@    IDIRLA DEVRAIT ETRE UN ENTIER EGAL A 1 2 OU 3           ',/,&
'@       (LA VALEUR 1 POUR UN ECOULEMENT SELON L''AXE X,      ',/,&
'@        LA VALEUR 2 POUR UN ECOULEMENT SELON L''AXE Y,      ',/,&
'@        LA VALEUR 3 POUR UN ECOULEMENT SELON L''AXE Z)      ',/,&
'@       IL VAUT ICI IDIRLA = ', I10                           ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier la valeur de IDIRLA dans la subroutine USLAG1.   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 2000   format(                                                         &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========                                               ',/,&
'@    LE MODULE LAGRANGIEN EST INCOMPATIBLE AVEC LE MODELE    ',/,&
'@    DE TURBULENCE SELECTIONNE.                              ',/,&
'@                                                            ',/,&
'@   Le module Lagrangien a ete active avec IILAGR = ',I10     ,/,&
'@     et la dispersion turbulente est prise en compte        ',/,&
'@                                     avec IDISTU = ',I10     ,/,&
'@   Le modele de turbulence active                           ',/,&
'@     correspond a ITURB  = ',I10                             ,/,&
'@   Or, les seuls traitements de la turbulence compatibles   ',/,&
'@     avec le module Lagrangien et la dispersion turbulente  ',/,&
'@     sont k-epsilon et Rij-epsilon, v2f et k-omega.         ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier la valeur de IILAGR et IDISTU dans la subroutine ',/,&
'@  USLAG1 et verifier la valeur de ITURB  dans la subroutine ',/,&
'@  USINI1.                                                   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

!----
! FIN
!----

return
end subroutine lagcar
