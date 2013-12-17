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

subroutine cs_coal_scast &
!=======================

 ( iscal  ,                                                       &
   rtpa   , rtp    , propce ,                                     &
   smbrs  , rovsdt )

!===============================================================================
! FONCTION :
! ----------

! ROUTINE PHYSIQUE PARTICULIERE : FLAMME CHARBON PULVERISE
!   ON PRECISE LES TERMES SOURCES POUR UN SCALAIRE PP
!   SUR UN PAS DE TEMPS

! ATTENTION : LE TRAITEMENT DES TERMES SOURCES EST DIFFERENT
! ---------   DE CELUI DE USTSSC.F

! ON RESOUT ROVSDT*D(VAR) = SMBRS

! ROVSDT ET SMBRS CONTIENNENT DEJA D'EVENTUELS TERMES SOURCES
!  UTILISATEUR. IL FAUT DONC LES INCREMENTER ET PAS LES
!  ECRASER

! POUR DES QUESTIONS DE STABILITE, ON NE RAJOUTE DANS ROVSDT
!  QUE DES TERMES POSITIFS. IL N'Y A PAS DE CONTRAINTE POUR
!  SMBRS

! DANS LE CAS D'UN TERME SOURCE EN CEXP + CIMP*VAR ON DOIT
! ECRIRE :
!          SMBRS  = SMBRS  + CEXP + CIMP*VAR
!          ROVSDT = ROVSDT + MAX(-CIMP,ZERO)

! ON FOURNIT ICI ROVSDT ET SMBRS (ILS CONTIENNENT RHO*VOLUME)
!    SMBRS en kg variable/s :
!     ex : pour la vitesse            kg m/s2
!          pour les temperatures      kg degres/s
!          pour les enthalpies        Joules/s
!    ROVSDT en kg /s

!-------------------------------------------------------------------------------
!ARGU                             ARGUMENTS
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! iscal            ! i  ! <-- ! scalar number                                  !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! smbrs(ncelet)    ! tr ! --> ! second membre explicite                        !
! rovsdt(ncelet    ! tr ! --> ! partie diagonale implicite                     !
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

integer          iscal

double precision rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision smbrs(ncelet), rovsdt(ncelet)

! Local variables

character*80     chaine
character*80     fname
character*80     name
integer          ivar , iel
integer          numcla , numcha , icla
integer          ipcgch , ipcgd1 , ipcgd2 , ipcght , ipcsec
integer          ipghco2 , ipghh2o
integer          ixchcl , ixckcl
integer          ipcro2 , ipcte1 , ipcte2 , ipcvsl , ipccp
integer          ipcdia , ipcvst
integer          mode, ige
integer          ipcx2c , icha , ii, jj
integer          itermx,nbpauv,nbrich,nbepau,nberic,ipghc2
integer          iterch,nbpass,nbarre,nbimax
integer          iexp1 , iexp2 , iexp3
integer          iok1
integer          keyccl

double precision xnuss
double precision aux, rhovst
double precision coefe(ngazem)
double precision t1, t2, hh2ov
double precision f1mc(ncharm), f2mc(ncharm)
double precision xhdev1 , xhdev2 , xhco , xho2 , gamdv1 , gamdv2
double precision xhco2  , xhh2   , xhh2o

double precision gamhet , den
double precision xxco,xxo2,xxco2,xxh2o
double precision xkp,xkm,t0p,t0m
double precision anmr,tauchi,tautur
double precision sqh2o , x2
double precision err1mx,err2mx
double precision errch,xden
double precision fn0,fn1,fn2,anmr0,anmr1,anmr2
double precision lnk0p,l10k0e,lnk0m,t0e,xco2eq,xcoeq,xo2eq
double precision xcom,xo2m,xkcequ,xkpequ

double precision  xw1,xw2,xw3,xw4
double precision  xo2,wmel,wmhcn,wmno,wmo2,wmnh3
double precision gmdev1(ncharm),gmdev2(ncharm),gmhet(ncharm)
double precision aux1 , aux2 , aux3
double precision xch,xck,xash,xmx2
double precision tfuelmin,tfuelmax
double precision auxdev,auxht3,auxco2,auxh2o,auxwat

double precision , dimension ( : )     , allocatable :: w1,w2,w3,w4,w5
double precision , dimension ( : )     , allocatable :: tfuel
double precision, dimension(:), pointer :: gamvlei, gamvloi, xchcpi, gaheto2i, xckcpi
double precision, dimension(:), pointer :: gaseci, frmcpi, agei, gahetco2i
double precision, dimension(:), pointer :: gaheth2oi, xwtcpi, xacpip
double precision, dimension(:), pointer ::  crom

!LOCAL VARIABLES
!===============
!
! Pointers de variables d'etat
! ----------------------------
! (Temperatur d'une particule de iclas,
!  Pointeur sur CHx1,
!  Pointeur sur CHx2,
!  Masse volumique de la pahse gaz,
!  Masse volumique du melange)
integer ipctem, ipcyf1, ipcyf2, idgaz
!
! Constante cinetiques laminaires
! -------------------------------
! (NH3 + O2,
!  NH3 + NO,
!  reburning)
integer iexp4,iexp5,iexprb
! Variables auxiliaire
! -------------------------------------
double precision aux4,aux5,auxhet,auxrb1,auxrb2
double precision core1,core2,core3,para2
double precision ychx
! Combustion heterogene
! ---------------------
! (Taux de reaction de la comb. het. du char 1,
!  Taux de reaction de la comb. het. du char 2)
double precision mckcl1, mckcl2
!
!===============================================================================
!
!===============================================================================
! 1. INITIALISATION
!===============================================================================
! --- Numero du scalaire a traiter : ISCAL
! --- Numero de la variable associee au scalaire a traiter ISCAL
ivar = isca(iscal)

! --- Nom de la variable associee au scalaire a traiter ISCAL
call field_get_label(ivarfl(ivar), chaine)

! --- Numero des grandeurs physiques
call field_get_val_s(icrom, crom)
ipcvst = ipproc(ivisct)
ipcte1 = ipproc(itemp1)

! Key id of the coal scalar class
call field_get_key_id("scalar_class", keyccl)

!===============================================================================
! Deallocation dynamic arrays
!----
allocate(w1(1:ncel),w2(1:ncel),w3(1:ncel),w4(1:ncel),w5(1:ncel),stat=iok1)
if (iok1 > 0) then
  write(nfecra,*) ' Memory allocation error inside: '
  write(nfecra,*) '     cs_coal_scast               '
  call csexit(1)
endif
allocate(tfuel(1:ncel),stat=iok1)
if (iok1 > 0) then
  write(nfecra,*) ' Memory allocation error inside: '
  write(nfecra,*) '     cs_coal_scast               '
  call csexit(1)
endif
!===============================================================================

!===============================================================================
! 2. PRISE EN COMPTE DES TERMES SOURCES POUR LES VARIABLES RELATIVES
!    AUX CLASSES DE PARTICULES
!===============================================================================

! --> Terme source pour la fraction massique de charbon reactif

if ( ivar.ge.isca(ixch(1)) .and. ivar.le.isca(ixch(nclacp)) ) then

  if (iwarni(ivar).ge.1) then
    write(nfecra,1000) chaine(1:8)
  endif
  numcla = ivar-isca(ixch(1))+1
  ipcgch = ipproc(igmdch(numcla))

  do iel = 1, ncel

! ---- Calcul de  W1 = - rho.GMDCH > 0

    xw1 = - crom(iel)*propce(iel,ipcgch)*volume(iel)

! ---- Calcul des parties explicite et implicite du TS

    rovsdt(iel) = rovsdt(iel) + max(xw1,zero)
    smbrs (iel) = smbrs(iel)  - xw1*rtpa(iel,ivar)

  enddo

endif


! --> Terme source pour la fraction massique de coke
if ( ivar.ge.isca(ixck(1)) .and. ivar.le.isca(ixck(nclacp)) ) then

  if (iwarni(ivar).ge.1) then
    write(nfecra,1000) chaine(1:8)
  endif

  numcla = ivar-isca(ixck(1))+1
  ipcgch = ipproc(igmdch(numcla))
  ipcgd1 = ipproc(igmdv1(numcla))
  ipcgd2 = ipproc(igmdv2(numcla))
  ixchcl = isca(ixch(numcla))
  ixckcl = isca(ixck(numcla))
  ipcght = ipproc(igmhet(numcla))
  if ( ihtco2 .eq. 1 ) then
    ipghco2 = ipproc(ighco2(numcla))
  endif
  if ( ihth2o .eq. 1 ) then
    ipghh2o = ipproc(ighh2o(numcla))
  endif

  do iel = 1, ncel

! ---- Calcul de W1 = - rho.Xch.GMDCH.Volume > 0

    xw1=-crom(iel)*rtp(iel,ixchcl)*propce(iel,ipcgch)     &
                                           *volume(iel)

! AE : On prend RTP(IEL,IXCHCL) et pas RTPA(IEL,IXCHCL) afin
!      d'etre conservatif sur la masse

! ---- Calcul de W2 = rho.Xch.(GMDV1+GMDV2)Volume < 0

    xw2 = crom(iel)*rtp(iel,ixchcl)                  &
         *(propce(iel,ipcgd1)+propce(iel,ipcgd2))*volume(iel)

    if ( rtpa(iel,ixckcl) .gt. epsicp ) then

! Reaction C(s) + O2 ---> 0.5CO
! =============================

! ---- Calcul de la partie implicite  > 0 du TS relatif a GMHET

      xw3 = -2.d0/3.d0*crom(iel)*propce(iel,ipcght) &
           /(rtpa(iel,ixckcl))**(1.d0/3.d0)*volume(iel)

! ---- Calcul de la partie explicite < 0 du TS relatif a GMHET

      xw4 = crom(iel)*propce(iel,ipcght)             &
                * (rtpa(iel,ixckcl))**(2.d0/3.d0)*volume(iel)

    else
      xw3 = 0.d0
      xw4 = 0.d0
    endif

! ---- Calcul des parties explicite et implicite du TS

    rovsdt(iel) = rovsdt(iel) + max(xw3,zero)
    smbrs(iel)  = smbrs(iel)  + xw1 + xw2 + xw4

  enddo

  if ( ihtco2 .eq. 1 ) then

    do iel = 1, ncel

      if ( rtpa(iel,ixckcl) .gt. epsicp ) then

! Reaction C(s) + CO2 ---> 2CO
! =============================

! ---- Calcul de la partie implicite  > 0 du TS relatif a GMHET

        xw3 = -2.d0/3.d0*crom(iel)*propce(iel,ipghco2)     &
              /(rtpa(iel,ixckcl))**(1.d0/3.d0)*volume(iel)

! ---- Calcul de la partie explicite < 0 du TS relatif a GMHET

        xw4 = crom(iel)*propce(iel,ipghco2)                 &
             *(rtpa(iel,ixckcl))**(2.d0/3.d0)*volume(iel)

      else
        xw3 = 0.d0
        xw4 = 0.d0
      endif

! ---- Calcul des parties explicite et implicite du TS

      rovsdt(iel) = rovsdt(iel) + max(xw3,zero)
      smbrs(iel)  = smbrs(iel)  + xw4

    enddo

  endif
!
  if ( ihth2o .eq. 1 ) then

    do iel = 1, ncel

      if ( rtpa(iel,ixckcl) .gt. epsicp ) then

! Reaction C(s) + CO2 ---> 2CO
! =============================

! ---- Calcul de la partie implicite  > 0 du TS relatif a GMHET

        xw3 = -2.d0/3.d0*crom(iel)*propce(iel,ipghh2o)     &
              /(rtpa(iel,ixckcl))**(1.d0/3.d0)*volume(iel)

! ---- Calcul de la partie explicite < 0 du TS relatif a GMHET

        xw4 = crom(iel)*propce(iel,ipghh2o)                &
             *(rtpa(iel,ixckcl))**(2.d0/3.d0)*volume(iel)

      else
        xw3 = 0.d0
        xw4 = 0.d0
      endif

! ---- Calcul des parties explicite et implicite du TS

      rovsdt(iel) = rovsdt(iel) + max(xw3,zero)
      smbrs(iel)  = smbrs(iel)  + xw4

    enddo

  endif

endif


! --> Terme source pour la fraction massique de d'eau

if ( ippmod(iccoal) .eq. 1 ) then

  if ( ivar.ge.isca(ixwt(1)) .and.                                &
       ivar.le.isca(ixwt(nclacp)) ) then

    if (iwarni(ivar).ge.1) then
      write(nfecra,1000) chaine(1:8)
    endif
    numcla = ivar-isca(ixwt(1))+1
    numcha = ichcor(numcla)

    ipcsec = ipproc(igmsec(numcla))
    ipcx2c = ipproc(ix2(numcla))

    do iel = 1, ncel

! ---- Calcul des parties explicite et implicite du TS

     if ( rtpa(iel,ivar).gt. epsicp .and.                         &
          xwatch(numcha).gt. epsicp       ) then
       xw1 = crom(iel)*propce(iel,ipcsec)*volume(iel)    &
            *(1.d0/propce(iel,ipcx2c))*(1.d0/xwatch(numcha))

       rovsdt(iel) = rovsdt(iel) + max(xw1,zero)
       smbrs(iel)  = smbrs(iel)  - xw1*rtpa(iel,ivar)
     endif

    enddo

  endif

endif

if (i_coal_drift.eq.1) then

  call field_get_name(ivarfl(ivar), fname)

  ! Particle age source term
  if (fname(1:8).eq.'X_Age_CP') then

    ! index of the coal particle class
    call field_get_key_int(ivarfl(ivar), keyccl, icla)

    ! Array values at previous time step
    call field_get_val_prev_s_by_name(fname, xacpip)

    ! Light volatile's source term
    write(name,'(a6,i2.2)')'Ga_DV1' ,icla
    call field_get_val_s_by_name(name, gamvlei)

    ! Heavy volatile's source term
    write(name,'(a6,i2.2)')'Ga_DV2' ,icla
    call field_get_val_s_by_name(name, gamvloi)

    ! Fraction massique du charbon de icla
    write(name,'(a6,i2.2)')'Xch_CP' ,icla
    call field_get_val_s_by_name(name, xchcpi)

    ! Echelle temporelle de la combustion heter.
    write(name,'(a9,i2.2)')'Ga_HET_O2' ,icla
    call field_get_val_s_by_name(name, gaheto2i)

    ! Fraction massique du char de icla
    write(name,'(a6,i2.2)')'xck_cp' ,icla
    call field_get_val_s_by_name(name, xckcpi)

    ! Indicateur du charbon de classe icla
    numcha = ichcor(icla)

    ! Echelle temporelle de la sechage
    if (ippmod(iccoal) .eq. 1) then
      write(name,'(a6,i2.2)')'Ga_SEC' ,icla
      call field_get_val_s_by_name(name, gaseci)
    endif

    ! Fraction massique de la phase solide
    write(name,'(a6,i2.2)')'Frm_CP' ,icla
    call field_get_val_s_by_name(name, frmcpi)

    ! L'age des particules par cellule
    write(name,'(a6,i2.2)')'Age_CP' ,icla
    call field_get_val_s_by_name(name, agei)

    ! Echelle temporelle de la gazefication par CO2
    if (ihtco2 .eq. 1) then
      write(name,'(a10,i2.2)')'Ga_HET_CO2' ,icla
      call field_get_val_s_by_name(name, gahetco2i)
    endif

    ! Echelle temporelle de la gazefication par H2O
    if (ihth2o .eq. 1) then
      write(name,'(a10,i2.2)')'Ga_HET_H2O' ,icla
      call field_get_val_s_by_name(name, gaheth2oi)
    endif

    if (ippmod(iccoal).eq.1) then
      write(name,'(a6,i2.2)')'Xwt_CP' ,icla
      call field_get_val_s_by_name(name, xwtcpi)
    endif

    do iel = 1, ncel
      ! Flux de masse: Devolatilisation
      auxdev =  -(gamvlei(iel)+gamvloi(iel))*xchcpi(iel)
      ! Consommation de char par combustion heterogene
      if (xckcpi(iel) .gt. epsicp) then
        auxht3 = -gaheto2i(iel) * (xckcpi(iel))**(2.d0/3.d0)
      else
        auxht3 = 0.d0
      endif
      ! Flux de masse: Gazefication par CO2
      if (ihtco2 .eq. 1) then
        if (xckcpi(iel) .gt. epsicp ) then
          auxco2 = -gahetco2i(iel)* (xckcpi(iel))**(2.d0/3.d0)
        else
          auxco2 = 0.d0
        endif
      else
        auxco2 = 0.d0
      endif
      ! Flux de masse: Gazefication par H2O
      if (ihth2o .eq. 1) then
        if (xckcpi(iel) .gt. epsicp) then
          auxh2o = -gaheth2oi(iel)*(xckcpi(iel))**(2.d0/3.d0)
        else
          auxh2o = 0.d0
        endif
      else
        auxh2o = 0.d0
      endif
      !  Flux de masse: Sechage
      if (ippmod(iccoal) .eq. 1) then
        if (xwtcpi(iel).gt.epsicp .and.xwatch(numcha).gt.epsicp) then
          auxwat = gaseci(iel)/(frmcpi(iel)*xwatch(numcha))
        else
          auxwat = 0.d0
        endif
      else
        auxwat = 0.d0
      endif

      if (frmcpi(iel).gt.epsicp) then
         smbrs(iel) =  smbrs(iel) + (crom(iel) * volume(iel) *          &
                       (frmcpi(iel)                                    -          &
                       ((auxdev+auxht3+auxco2+auxh2o+auxwat)*xacpip(iel))) )
         rovsdt(iel) = rovsdt(iel) + max(((auxdev+auxht3+auxco2+auxh2o+auxwat)/   &
                                         frmcpi(iel)), zero)
      endif

    enddo
  endif

  ! Gas age source term
  if (fname(1:9).eq.'X_Age_Gas') then

    ! Loop over particle classes
    do icla = 1, nclacp
      ! Light volatile's source term
      write(name,'(a6,i2.2)')'Ga_DV1' ,icla
      call field_get_val_s_by_name(name, gamvlei)

      ! Heavy volatile's source term
      write(name,'(a6,i2.2)')'Ga_DV2' ,icla
      call field_get_val_s_by_name(name, gamvloi)

      ! Fraction massique du charbon de icla
      write(name,'(a6,i2.2)')'Xch_CP' ,icla
      call field_get_val_s_by_name(name, xchcpi)

      ! Echelle temporelle de la combustion heter.
      write(name,'(a9,i2.2)')'Ga_HET_O2' ,icla
      call field_get_val_s_by_name(name, gaheto2i)

      ! Fraction massique du char de icla
      write(name,'(a6,i2.2)')'xck_cp' ,icla
      call field_get_val_s_by_name(name, xckcpi)

      ! Indicateur du charbon de classe icla
      numcha = ichcor(icla)

      ! Echelle temporelle de la sechage
      if (ippmod(iccoal) .eq. 1) then
        write(name,'(a6,i2.2)')'Ga_SEC' ,icla
        call field_get_val_s_by_name(name, gaseci)
      endif

      ! Fraction massique de la phase solide
      write(name,'(a6,i2.2)')'Frm_CP' ,icla
      call field_get_val_s_by_name(name, frmcpi)

      ! L'age des particules par cellule
      write(name,'(a6,i2.2)')'Age_CP' ,icla
      call field_get_val_s_by_name(name, agei)

      ! Echelle temporelle de la gazefication par CO2
      if (ihtco2 .eq. 1) then
        write(name,'(a10,i2.2)')'Ga_HET_CO2' ,icla
        call field_get_val_s_by_name(name, gahetco2i)
      endif

      ! Echelle temporelle de la gazefication par H2O
      if (ihth2o .eq. 1) then
        write(name,'(a10,i2.2)')'Ga_HET_H2O' ,icla
        call field_get_val_s_by_name(name, gaheth2oi)
      endif

      if (ippmod(iccoal).eq.1) then
        write(name,'(a6,i2.2)')'Xwt_CP' ,icla
        call field_get_val_s_by_name(name, xwtcpi)
      endif

      do iel = 1, ncel
        ! Flux de masse: Devolatilisation
        auxdev = -(gamvlei(iel)+gamvloi(iel))*xchcpi(iel)
        ! Consommation de char par combustion heterogene
        if (xckcpi(iel) .gt. epsicp) then
          auxht3 = -gaheto2i(iel) * (xckcpi(iel))**(2.d0/3.d0)
        else
          auxht3 = 0.d0
        endif
        ! Flux de masse: Gazefication par CO2
        if (ihtco2 .eq. 1) then
          if (xckcpi(iel) .gt. epsicp ) then
            auxco2 = -gahetco2i(iel)* (xckcpi(iel))**(2.d0/3.d0)
          else
            auxco2 = 0.d0
          endif
        else
          auxco2 = 0.d0
        endif
        ! Flux de masse: Gazefication par H2O
        if (ihth2o .eq. 1) then
          if (xckcpi(iel) .gt. epsicp) then
            auxh2o = -gaheth2oi(iel)*(xckcpi(iel))**(2.d0/3.d0)
          else
            auxh2o = 0.d0
          endif
        else
          auxh2o = 0.d0
        endif
        !  Flux de masse: Sechage
        if (ippmod(iccoal) .eq. 1) then
          if (xwtcpi(iel).gt.epsicp .and.xwatch(numcha).gt.epsicp) then
            auxwat = gaseci(iel)/(frmcpi(iel)*xwatch(numcha))
          else
            auxwat = 0.d0
          endif
        else
          auxwat = 0.d0
        endif

        smbrs(iel) = smbrs(iel) + crom(iel)*volume(iel)               &
                                *( ( auxdev+auxht3+auxco2+auxh2o+auxwat)       &
                                   * agei(iel)                                 &
                                 - frmcpi(iel)                                 &
                                 )
      enddo

    enddo

    ! Finalization
    do iel = 1, ncel
      ! The formula is (1- SUM X2)
      smbrs(iel) =  smbrs(iel) + crom(iel) * volume(iel)
    enddo

  endif

endif

! --> Terme source pour l'enthalpie du solide

if ( ivar.ge.isca(ih2(1)) .and. ivar.le.isca(ih2(nclacp)) ) then

  if (iwarni(ivar).ge.1) then
    write(nfecra,1000) chaine(1:8)
  endif

  numcla = ivar-isca(ih2(1))+1
  numcha = ichcor(numcla)
  ixchcl = isca(ixch(numcla))
  ixckcl = isca(ixck(numcla))
  ipcx2c = ipproc(ix2(numcla))
  ipcro2 = ipproc(irom2(numcla ))
  ipcdia = ipproc(idiam2(numcla))
  ipcte2 = ipproc(itemp2(numcla))
  ipcght = ipproc(igmhet(numcla))
  if ( ihtco2 .eq. 1 ) then
    ipghco2 = ipproc(ighco2(numcla))
  endif
  if ( ihth2o .eq. 1 ) then
    ipghh2o = ipproc(ighh2o(numcla))
  endif
  ipcgd1 = ipproc(igmdv1(numcla))
  ipcgd2 = ipproc(igmdv2(numcla))
  ipcgch = ipproc(igmdch(numcla))
  ixchcl = isca(ixch(numcla))
  ixckcl = isca(ixck(numcla))

! ---- Contribution aux bilans explicite et implicite
!        des echanges par diffusion moleculaire
!        6 Lambda Nu / diam**2 / Rho2 * Rho * (T1-T2)

! ------ Calcul de lambda dans W1

  xnuss = 2.d0
  do iel = 1, ncel
    if ( ivisls(iscalt).gt.0 ) then
      ipcvsl = ipproc(ivisls(iscalt))
      if ( icp.gt.0 ) then
        ipccp   = ipproc(icp)
        w1(iel) = propce(iel,ipcvsl) * propce(iel,ipccp)
      else
        w1(iel) = propce(iel,ipcvsl) * cp0
      endif
    else
      if ( icp.gt.0 ) then
        ipccp   = ipproc(icp)
        w1(iel) = visls0(iscalt) * propce(iel,ipccp)
      else
        w1(iel) = visls0(iscalt) * cp0
      endif
    endif
  enddo

! ------ Calcul du diametre des particules dans W2
!        On calcule le d20 = (A0.D0**2+(1-A0)*DCK**2)**0.5

  do iel = 1, ncel
    w2(iel) = ( xashch(numcha)*diam20(numcla)**2 +                &
                (1.d0-xashch(numcha))*propce(iel,ipcdia)**2       &
              )**0.5
  enddo

! ------ Contribution aux bilans explicite et implicite de
!        des echanges par diffusion moleculaire
!      Rq : on utilise PROPCE(IEL,IPCX2C) car on veut X2 a l'iteration n

  do iel = 1, ncel
    if ( propce(iel,ipcx2c) .gt. epsicp ) then
      aux         = 6.d0 * w1(iel) * xnuss / w2(iel)**2           &
                  / propce(iel,ipcro2) * crom(iel)       &
                  * propce(iel,ipcx2c) * volume(iel)
      rhovst      = aux / cp2ch(numcha) /propce(iel,ipcx2c)

      smbrs(iel)  = smbrs(iel)-aux*(propce(iel,ipcte2)-propce(iel,ipcte1))
      rovsdt(iel) = rovsdt(iel)+ max(zero,rhovst)
    endif

  enddo


! ---- Contribution aux bilans explicite et implicite
!        du terme echange d'energie entre les phases :
!        GAMA(dev1) H(mv1,T2)+GAMA(dev2) H(mv2,T2)

  do iel = 1, ncel

!        Gama Dev1 et Gama Dev2

    gamdv1 = crom(iel)*rtp(iel,ixchcl)                   &
            *propce(iel,ipcgd1)

    gamdv2 = crom(iel)*rtp(iel,ixchcl)                   &
            *propce(iel,ipcgd2)

!        H(mv1,T2)

    do ige = 1, ngazem
      coefe(ige) = zero
    enddo
!
    den = a1(numcha)*wmole(ichx1c(numcha))+b1(numcha)*wmole(ico)  &
         +c1(numcha)*wmole(ih2o)          +d1(numcha)*wmole(ih2s) &
         +e1(numcha)*wmole(ihcn)          +f1(numcha)*wmole(inh3)
    coefe(ichx1) = a1(numcha)*wmole(ichx1c(numcha)) /den
    coefe(ico  ) = b1(numcha)*wmole(ico)            /den
    coefe(ih2o ) = c1(numcha)*wmole(ih2o)           /den
    coefe(ih2s ) = d1(numcha)*wmole(ih2s)           /den
    coefe(ihcn ) = e1(numcha)*wmole(ihcn)           /den
    coefe(inh3 ) = f1(numcha)*wmole(inh3)           /den

    t2         = propce(iel,ipcte2)
    do icha = 1, ncharm
      f1mc(icha) = zero
      f2mc(icha) = zero
    enddo
    f1mc(numcha) = 1.d0

    mode      = -1
    call cs_coal_htconvers1 &
    !======================
    ( mode , xhdev1 , coefe , f1mc , f2mc , t2 )

!        H(mv2,T2)

    do ige = 1, ngazem
      coefe(ige) = zero
    enddo
    den = a2(numcha)*wmole(ichx2c(numcha))+b2(numcha)*wmole(ico)  &
         +c2(numcha)*wmole(ih2o)          +d2(numcha)*wmole(ih2s) &
         +e2(numcha)*wmole(ihcn)          +f2(numcha)*wmole(inh3)
    coefe(ichx2) = a2(numcha)*wmole(ichx2c(numcha)) /den
    coefe(ico  ) = b2(numcha)*wmole(ico)            /den
    coefe(ih2o ) = c2(numcha)*wmole(ih2o)           /den
    coefe(ih2s ) = d2(numcha)*wmole(ih2s)           /den
    coefe(ihcn ) = e2(numcha)*wmole(ihcn)           /den
    coefe(inh3 ) = f2(numcha)*wmole(inh3)           /den

    t2         = propce(iel,ipcte2)
    do icha = 1, ncharm
      f1mc(icha) = zero
      f2mc(icha) = zero
    enddo
    f2mc(numcha) = 1.d0

    mode      = -1
    call cs_coal_htconvers1 &
    !======================
    ( mode , xhdev2 , coefe , f1mc , f2mc , t2 )

!         Contribution aux bilans explicite et implicite

    smbrs(iel) = smbrs(iel)+(gamdv1*xhdev1+gamdv2*xhdev2)*volume(iel)

  enddo

! ------ combustion heterogene : C(s) + 02 ---> 0.5 C0
!        GamHET * (28/12 H(CO,T2)-16/12 H(O2,T1) )

  do iel = 1, ncel

!        Calcul de HCO(T2)

    do ige = 1, ngazem
      coefe(ige) = zero
    enddo
    coefe(ico) = 1.d0
    do icha = 1, ncharm
      f1mc(icha) = zero
      f2mc(icha) = zero
    enddo

    t2        = propce(iel,ipcte2)
    mode      = -1
    call cs_coal_htconvers1 &
    !======================
    ( mode , xhco , coefe , f1mc , f2mc , t2 )

!        Calcul de HO2(T1)

    do ige = 1, ngazem
      coefe(ige) = zero
    enddo
    coefe(io2) = 1.d0
    do icha = 1, ncharm
      f1mc(icha) = zero
      f2mc(icha) = zero
    enddo

    t1        = propce(iel,ipcte1)
    mode      = -1
    call cs_coal_htconvers1 &
    !======================
    ( mode , xho2 , coefe , f1mc , f2mc , t1 )

!         Contribution aux bilans explicite et implicite

    if ( rtpa(iel,ixckcl) .gt. epsicp ) then

      gamhet = crom(iel)*propce(iel,ipcght)              &
               * ( (rtpa(iel,ixckcl))**(2.d0/3.d0) +              &
                2.d0/3.d0*(rtp(iel,ixckcl)-rtpa(iel,ixckcl))      &
                 /(rtpa(iel,ixckcl))**(1.d0/3.d0) )

    else
      gamhet = 0.d0
    endif

   smbrs(iel) = smbrs(iel)                                        &
                 +gamhet                                          &
                  *(28.d0/12.d0*xhco-16.d0/12.d0*xho2)            &
                  *volume(iel)

  enddo

! ------ combustion heterogene : C(s) + C02 ---> 2 C0
!        GamHET * (56/12 H(CO,T2)-44/12 H(CO2,T1) )

  if ( ihtco2 .eq. 1 ) then
    do iel = 1, ncel

!        Calcul de HCO(T2)

      do ige = 1, ngazem
        coefe(ige) = zero
      enddo
      coefe(ico) = 1.d0
      do icha = 1, ncharm
        f1mc(icha) = zero
        f2mc(icha) = zero
      enddo

      t2        = propce(iel,ipcte2)
      mode      = -1
      call cs_coal_htconvers1 &
      !======================
      ( mode , xhco , coefe , f1mc , f2mc , t2  )

!        Calcul de HCO2(T1)

      do ige = 1, ngazem
        coefe(ige) = zero
      enddo
      coefe(ico2) = 1.d0
      do icha = 1, ncharm
        f1mc(icha) = zero
        f2mc(icha) = zero
      enddo

      t1        = propce(iel,ipcte1)
      mode      = -1
      call cs_coal_htconvers1 &
      !======================
      ( mode , xhco2 , coefe , f1mc , f2mc , t1    )

!         Contribution aux bilans explicite et implicite

      if ( rtpa(iel,ixckcl) .gt. epsicp ) then

        gamhet = crom(iel)*propce(iel,ipghco2)           &
                 * ( (rtpa(iel,ixckcl))**(2.d0/3.d0) +            &
                2.d0/3.d0*(rtp(iel,ixckcl)-rtpa(iel,ixckcl))      &
                 /(rtpa(iel,ixckcl))**(1.d0/3.d0) )

      else
        gamhet = 0.d0
      endif

     smbrs(iel) = smbrs(iel)                                      &
                   +gamhet                                        &
                    *(56.d0/12.d0*xhco-44.d0/12.d0*xhco2)         &
                    *volume(iel)

    enddo

  endif
! ------ combustion heterogene : C(s) + H2O ---> CO + H2
!        GamHET * (28/12 H(CO,T2)+2/12 H(HY,T2) -18/12 H(H2O,T1) )

  if ( ihth2o .eq. 1 ) then
    do iel = 1, ncel

!        Calcul de HCO(T2)

      do ige = 1, ngazem
        coefe(ige) = zero
      enddo
      coefe(ico) = 1.d0
      do icha = 1, ncharm
        f1mc(icha) = zero
        f2mc(icha) = zero
      enddo

      t2        = propce(iel,ipcte2)
      mode      = -1
      call cs_coal_htconvers1 &
      !======================
      ( mode , xhco , coefe , f1mc , f2mc , t2 )

!        Calcul de HH2(T2)

      do ige = 1, ngazem
        coefe(ige) = zero
      enddo
      coefe(ihy) = 1.d0
      do icha = 1, ncharm
        f1mc(icha) = zero
        f2mc(icha) = zero
      enddo

      t2        = propce(iel,ipcte2)
      mode      = -1
      call cs_coal_htconvers1 &
      !======================
      ( mode , xhh2 , coefe , f1mc , f2mc , t2    )

!        Calcul de HH2O(T1)

      do ige = 1, ngazem
        coefe(ige) = zero
      enddo
      coefe(ih2o) = 1.d0
      do icha = 1, ncharm
        f1mc(icha) = zero
        f2mc(icha) = zero
      enddo

      t1        = propce(iel,ipcte1)
      mode      = -1
      call cs_coal_htconvers1 &
      !======================
      ( mode , xhh2o , coefe , f1mc , f2mc , t1    )

!         Contribution aux bilans explicite et implicite

      if ( rtpa(iel,ixckcl) .gt. epsicp ) then

        gamhet = crom(iel)*propce(iel,ipghh2o)           &
                 * ( (rtpa(iel,ixckcl))**(2.d0/3.d0) +            &
                2.d0/3.d0*(rtp(iel,ixckcl)-rtpa(iel,ixckcl))      &
                 /(rtpa(iel,ixckcl))**(1.d0/3.d0) )

      else
        gamhet = 0.d0
      endif

     smbrs(iel) = smbrs(iel)                                    &
                 +gamhet                                        &
                  *(28.d0/12.d0*xhco+ 2.d0/12.d0*xhh2           &
                                    -18.d0/12.d0*xhh2o )        &
                    *volume(iel)

    enddo

  endif

!       --> Terme source sur H2 issu du sechage)

  if ( ippmod(iccoal) .eq. 1 ) then

! ---- Contribution du TS interfacial aux bilans explicite et implicite


    ipcsec = ipproc(igmsec(numcla))
    ipcte2 = ipproc(itemp2(numcla))

    do iel = 1, ncel

!          Calcul de H(H2O) a T2

      do ige = 1, ngazem
        coefe(ige) = zero
      enddo
      coefe(ih2o) = 1.d0
      do icha = 1, ncharm
        f1mc(icha) = zero
        f2mc(icha) = zero
      enddo

      t2 = propce(iel,ipcte2)
      if ( t2 .gt. 100.d0+tkelvi ) then
        t2 = 100.d0+tkelvi
      endif
      mode      = -1
      call cpthp1 &
      !==========
      ( mode , hh2ov , coefe , f1mc , f2mc ,  t2  )

!         Contribution aux bilans explicite

      if ( rtpa(iel,isca(ixwt(numcla))).gt. epsicp .and.          &
           xwatch(numcha) .gt. epsicp       ) then

        aux = -crom(iel)*propce(iel,ipcsec)              &
       *(rtp(iel,isca(ixwt(numcla)))/propce(iel,ipcx2c))          &
       *(1.d0                    /xwatch(numcha))                 &
             *hh2ov

      else
        aux = 0.d0
      endif

      smbrs(iel) = smbrs(iel) + aux*volume(iel)

    enddo

  endif
endif

!===============================================================================
! 3. PRISE EN COMPTE DES TERMES SOURCES POUR LES VARIABLES RELATIVES
!    AU MELANGE GAZEUX
!===============================================================================

! --> Terme source pour les matieres volatiles legeres

if ( ivar.ge.isca(if1m(1)) .and. ivar.le.isca(if1m(ncharb)) ) then

  if (iwarni(ivar).ge.1) then
    write(nfecra,1000) chaine(1:8)
  endif

! ---- Calcul de GMDEV1 = - SOMME (rho.XCH.GMDV1) > 0  --> W1

  numcha = ivar-isca(if1m(1))+1
  do iel = 1, ncel
    w1(iel) = zero
  enddo
  do icla = 1, nclacp
    ipcgd1 = ipproc(igmdv1(icla))
    ixchcl = isca(ixch(icla))
    if ( ichcor(icla).eq.numcha ) then
      do iel = 1, ncel
        w1(iel) = w1(iel) - crom(iel)*rtp(iel,ixchcl)    &
                * propce(iel,ipcgd1)
      enddo
    endif
  enddo

! ---- Contribution du TS interfacial aux bilans explicite et implicite

  do iel = 1, ncel
    smbrs(iel)  = smbrs(iel)  + volume(iel) * w1(iel)
  enddo

endif


! --> Terme source pour les matieres volatiles lourdes

if ( ivar.ge.isca(if2m(1)) .and. ivar.le.isca(if2m(ncharb)) ) then

  if (iwarni(ivar).ge.1) then
    write(nfecra,1000) chaine(1:8)
  endif

! ---- Calcul de GMDEV2 = - SOMME (rho.XCH.GMDV2) >0 --> W1

  numcha = ivar-isca(if2m(1))+1
  do iel = 1, ncel
    w1(iel) = zero
  enddo
  do icla = 1, nclacp
    ipcgd2 = ipproc(igmdv2(icla))
    ixchcl = isca(ixch(icla))
    if ( ichcor(icla).eq.numcha ) then
      do iel = 1, ncel
        w1(iel) = w1(iel) - crom(iel)*rtp(iel,ixchcl)    &
                * propce(iel,ipcgd2)
      enddo
    endif
  enddo

! ---- Contribution du TS interfacial pour le bilan explicite

  do iel = 1, ncel
    smbrs(iel)  = smbrs(iel)  + volume(iel) * w1(iel)
  enddo

endif


! --> Terme source pour le traceur 7 (O2) (C de la comb. het.)

if ( ivar.eq.isca(if7m) ) then

! RQ IMPORTANTE :  On prend les meme TS que pour Xck
!                  afin d'etre conservatif

  if (iwarni(ivar).ge.1) then
    write(nfecra,1000) chaine(1:8)
  endif

  do iel = 1, ncel
    w1(iel) = zero
  enddo

  do icla = 1, nclacp
    ipcght = ipproc(igmhet(icla))
    ixckcl = isca(ixck(icla))
    do iel = 1, ncel
      if ( rtpa(iel,ixckcl) .gt. epsicp ) then
        w1(iel) = w1(iel)                                         &
                 - crom(iel)*propce(iel,ipcght)          &
                 * ( (rtpa(iel,ixckcl))**(2.d0/3.d0) +            &
                    2.d0/3.d0*(rtp(iel,ixckcl)-rtpa(iel,ixckcl))  &
                    /(rtpa(iel,ixckcl))**(1.d0/3.d0) )
      endif
    enddo
  enddo

! ---- Contribution du TS interfacial aux bilans explicite et implicite

  do iel = 1, ncel
    smbrs(iel)  = smbrs(iel)  + volume(iel) * w1(iel)
  enddo

endif



! --> Terme source pour le traceur 8 (CO2) (C de la comb. het.)

if ( ihtco2 .eq. 1 ) then
  if ( ivar.eq.isca(if8m) ) then

! RQ IMPORTANTE :  On prend les meme TS que pour Xck
!                  afin d'etre conservatif

    if (iwarni(ivar).ge.1) then
      write(nfecra,1000) chaine(1:8)
    endif

    do iel = 1, ncel
      w1(iel) = zero
    enddo

    do icla = 1, nclacp
      ipcght = ipproc(ighco2(icla))
      ixckcl = isca(ixck(icla))
      do iel = 1, ncel
        if ( rtpa(iel,ixckcl) .gt. epsicp ) then
          w1(iel) = w1(iel)                                        &
                   - crom(iel)*propce(iel,ipcght)         &
                   * ( (rtpa(iel,ixckcl))**(2.d0/3.d0) +           &
                      2.d0/3.d0*(rtp(iel,ixckcl)-rtpa(iel,ixckcl)) &
                      /(rtpa(iel,ixckcl))**(1.d0/3.d0) )
        endif
      enddo
    enddo

! ---- Contribution du TS interfacial aux bilans explicite et implicite

    do iel = 1, ncel
      smbrs(iel)  = smbrs(iel)  + volume(iel) * w1(iel)
    enddo

  endif

endif

! --> Terme source pour le traceur 9 (H2O) (comb. het. par h2O)

if ( ihth2o .eq. 1 ) then
  if ( ivar.eq.isca(if9m) ) then

! RQ IMPORTANTE :  On prend les meme TS que pour Xck
!                  afin d'etre conservatif

    if (iwarni(ivar).ge.1) then
      write(nfecra,1000) chaine(1:8)
    endif
!
    do iel = 1, ncel
      w1(iel) = zero
    enddo

    do icla = 1, nclacp
      ipcght = ipproc(ighh2o(icla))
      ixckcl = isca(ixck(icla))
      do iel = 1, ncel
        if ( rtpa(iel,ixckcl) .gt. epsicp ) then
          w1(iel) = w1(iel)                                        &
                   - crom(iel)*propce(iel,ipcght)        &
                   * ( (rtpa(iel,ixckcl))**(2.d0/3.d0) +           &
                      2.d0/3.d0*(rtp(iel,ixckcl)-rtpa(iel,ixckcl)) &
                      /(rtpa(iel,ixckcl))**(1.d0/3.d0) )
        endif
      enddo
    enddo

! ---- Contribution du TS interfacial aux bilans explicite et implicite

    do iel = 1, ncel
      smbrs(iel)  = smbrs(iel)  + volume(iel) * w1(iel)
    enddo

  endif

endif


! --> Terme source pour la variance du combustible

if ( ivar.eq.isca(ifvp2m) ) then

  if (iwarni(ivar).ge.1) then
    write(nfecra,1000) chaine(1:8)
  endif

  call cs_coal_fp2st                                               &
 !==================
 ( iscal  ,                                                        &
   rtpa   , rtp    , propce ,                                      &
   smbrs  , rovsdt )

endif

! --> Terme source pour le traceur 6 (Eau issue du séchage)

if ( ippmod(iccoal) .eq. 1 ) then

  if ( ivar.eq.isca(if6m) ) then


    if (iwarni(ivar).ge.1) then
      write(nfecra,1000) chaine(1:8)
    endif

! ---- Contribution du TS interfacial aux bilans explicite et implicite

    do iel = 1, ncel
      w1(iel) = zero
    enddo

    do icla = 1, nclacp

      ipcsec = ipproc(igmsec(icla))
      ipcx2c = ipproc(ix2(icla))
      numcha = ichcor(icla)

      do iel = 1, ncel

        if (  rtpa(iel,isca(ixwt(icla))).gt. epsicp               &
            .and.                                                 &
              xwatch(numcha) .gt. epsicp       ) then

          w1(iel) = w1(iel)                                       &
      + crom(iel)*propce(iel,ipcsec)                     &
       *(rtp(iel,isca(ixwt(icla)))/propce(iel,ipcx2c))            &
       *(1.d0                  /xwatch(numcha))

        endif

      enddo

    enddo

    do iel = 1, ncel
      smbrs(iel)  = smbrs(iel)  + volume(iel) * w1(iel)
    enddo

  endif

endif

! --> Terme source pour CO2

if ( ieqco2 .eq. 1 ) then

  if ( ivar.eq.isca(iyco2) ) then


    if (iwarni(ivar).ge.1) then
      write(nfecra,1000) chaine(1:8)
    endif

! ---- Contribution du TS interfacial aux bilans explicite et implicite

! Oxydation du CO
! ===============

!  Dryer Glassman : XK0P en (moles/m3)**(-0.75) s-1
!          XK0P = 1.26D10
!           XK0P = 1.26D7 * (1.1)**(NTCABS)
!           IF ( XK0P .GT. 1.26D10 ) XK0P=1.26D10
!           T0P  = 4807.D0
!  Howard : XK0P en (moles/m3)**(-0.75) s-1
!             XK0P = 4.11D9
!             T0P  = 15090.D0
!  Westbrook & Dryer

    lnk0p = 23.256d0
    t0p  = 20096.d0
!
!  Hawkin et Smith Purdue University Engeneering Bulletin, i
!  Research series 108 vol 33, n 3n 1949
!  Kp = 10**(4.6-14833/T)
!  Constante d'equilibre en pression partielle (atm           !)
!  XKOE est le log decimal de la constante pre-exponentielle
!  TOE  n'est PAS une temerature d'activation  ... il reste un lg(e)
!  pour repasser en Kc et utiliser des concetrations (moles/m3)
!  Kc = (1/RT)**variation nb moles * Kp
!  ici Kc = sqrt(0.082*T)*Kp

    l10k0e = 4.6d0
    t0e  = 14833.d0
! Dissociation du CO2 (Trinh Minh Chinh)
! ===================
!          XK0M = 5.D8
!          T0M  = 4807.D0
!          XK0M = 0.D0
!  Westbrook & Dryer

    lnk0m = 20.03d0
    t0m  = 20096.d0

    err1mx = 0.d0
    err2mx = 0.d0

! Nombre d'iterations
    itermx = 500
! Nombre de points converges

   nbpauv = 0
   nbepau = 0
   nbrich = 0
   nberic = 0
   nbpass = 0
   nbarre = 0
   nbimax = 0
! Precision pour la convergence
   errch = 1.d-8
!
   do iel = 1, ncel
!
     xxco  = propce(iel,ipproc(iym1(ico  )))/wmole(ico)           &
            *propce(iel,ipproc(irom1))
     xxo2  = propce(iel,ipproc(iym1(io2  )))/wmole(io2)           &
            *propce(iel,ipproc(irom1))
     xxco2 = propce(iel,ipproc(iym1(ico2 )))/wmole(ico2)          &
            *propce(iel,ipproc(irom1))
     xxh2o = propce(iel,ipproc(iym1(ih2o )))/wmole(ih2o)          &
            *propce(iel,ipproc(irom1))
!
     xxco  = max(xxco ,zero)
     xxo2  = max(xxo2 ,zero)
     xxco2 = max(xxco2,zero)
     xxh2o = max(xxh2o,zero)
     sqh2o = sqrt(xxh2o)
!
     xkp = exp(lnk0p-t0p/propce(iel,ipproc(itemp1)))
     xkm = exp(lnk0m-t0m/propce(iel,ipproc(itemp1)))
!
     xkpequ = 10.d0**(l10k0e-t0e/propce(iel,ipproc(itemp1)))
     xkcequ = xkpequ                                              &
             /sqrt(8.32d0*propce(iel,ipproc(itemp1))/1.015d5)

!        initialisation par l'état transporté

     anmr  = xxco2
     xcom  = xxco + xxco2
     xo2m  = xxo2 + 0.5d0*xxco2

     if ( propce(iel,ipproc(itemp1)) .gt. 1200.d0 ) then

!           Recherche de l'état d'équilibre
!           Recerche itérative sans controle de convergence
!            (pour conserver la parallelisation sur les mailles)
!           sur le nombre de moles de reaction séparant
!           l'etat avant réaction (tel que calculé par Cpcym)
!           de l'état d'équilibre
!          ANMR doit etre borne entre 0 et Min(XCOM,2.*XO2M)
!          on recherche la solution par dichotomie

       anmr0 = 0.d0
       anmr1 = min(xcom,2.d0*xo2m)
       iterch = 0
       fn2    = 1.d0
       fn0  = -0.5d0                           * anmr0**3         &
            + (     xcom    + xo2m - xkcequ**2) * anmr0**2        &
            - (.5d0*xcom    +2.d0*xo2m)*xcom   * anmr0            &
            +       xcom**2 * xo2m
       fn1  = -0.5d0                           * anmr1**3         &
            + (     xcom    + xo2m - xkcequ**2) * anmr1**2        &
            - (.5d0*xcom    +2.d0*xo2m)*xcom   * anmr1            &
            +       xcom**2 * xo2m

       if ( xo2m.gt.1.d-6) then
         do while ( iterch.lt.itermx .and. fn2.gt.errch )
           anmr2 = 0.5d0*(anmr0+anmr1)
           fn2  = -0.5d0                            * anmr2**3    &
                + (     xcom    + xo2m - xkcequ**2) * anmr2**2    &
                - (.5d0*xcom    +2.d0*xo2m)*xcom    * anmr2       &
                +       xcom**2 * xo2m
           if(fn0*fn2 .gt. 0.d0) then
             anmr0 = anmr2
             fn0 = fn2
           elseif(fn1*fn2 .gt. 0.d0) then
             anmr1 = anmr2
             fn1 = fn2
           elseif(fn0*fn1 .gt. 0.d0) then
             iterch = itermx
             anmr2 = min(xcom,2.d0*xo2m)
             nbarre = nbarre + 1
           endif
           iterch = iterch + 1
         enddo
!
         if ( iterch .ge. itermx) then
           nberic = nberic + 1
         else
           nbimax = max(nbimax,iterch)
         endif
         err1mx = max(err1mx,fn2)

         xco2eq = anmr2
         xcoeq  = xcom - anmr2
         xo2eq  = xo2m - 0.5d0 * anmr2
       else
         xo2eq  = 0.d0
         xcoeq  = xxco
         xco2eq = 0.d0
       endif

     else

       xco2eq = min(xcom,2.d0*xo2m)
       xo2eq  = xo2m - 0.5d0*xco2eq
       xcoeq  = xcom - xco2eq

     endif
!
     if ( xco2eq.gt.xxco2 ) then
!           oxydation
       xden = xkp*sqh2o*(xxo2)**0.25d0
     else
!           dissociation
       xden = xkm
     endif
     if ( xden .ne. 0.d0 ) then

       tauchi = 1.d0/xden
       tautur = rtpa(iel,ik)/rtpa(iel,iep)

       x2 = 0.d0
       do icla = 1, nclacp
         x2 = x2 + propce(iel,ipproc(ix2(icla)))
       enddo

!    On transporte CO2

       smbrs(iel)  = smbrs(iel)                                   &
                    +wmole(ico2)/propce(iel,ipproc(irom1))        &
         * (xco2eq-xxco2)/(tauchi+tautur)                         &
         * (1.d0-x2)                                              &
         * volume(iel) * crom(iel)
!
       w1(iel) = volume(iel)*crom(iel)/(tauchi+tautur)
       rovsdt(iel) = rovsdt(iel) +   max(w1(iel),zero)

     else
       rovsdt(iel) = rovsdt(iel) + 0.d0
       smbrs(iel)  = smbrs(iel)  + 0.d0
     endif

   enddo

   if(irangp.ge.0) then
     call parcpt(nberic)
     call parmax(err1mx)
     call parcpt(nbpass)
     call parcpt(nbarre)
     call parcpt(nbarre)
     call parcmx(nbimax)
   endif

   write(nfecra,*) ' Max Error = ',ERR1MX
   write(nfecra,*) ' no Points   ',NBERIC,NBARRE,NBPASS
   write(nfecra,*) ' Iter max number ',NBIMAX

!     Terme source : combustion heterogene par le CO2

   if ( ihtco2 .eq. 1) then

     do iel = 1, ncel

       aux = 0.d0
       do icla = 1,nclacp

         ixckcl = isca(ixck(icla))
         ipghc2 = ipproc(ighco2(icla))

         aux = aux                                                &
              + crom(iel)*propce(iel,ipghc2)             &
               *(rtpa(iel,ixckcl))**(2.d0/3.d0)*volume(iel)

       enddo

       rovsdt(iel) = rovsdt(iel) - aux*(wmole(ico2)/0.012)

     enddo

   endif

  endif

endif
!
! --> Terme source pour Enth_Ox
!                       HCN et NO : uniquement a partir de la 2eme iters
!
if ( ieqnox .eq. 1 .and. ntcabs .gt. 1) then
!
! Termes sur l'enthalpie Oxydant
!
  if ( ivar .eq. isca(ihox) ) then
!
!  Calcul de T2 moy sur la particules
!
    tfuelmin = 1.d+20
    tfuelmax =-1.d+20
    do iel=1,ncel
!
      xmx2 = 0.d0
      do icla = 1, nclacp
        xck  = rtp(iel,isca(ixck(icla)))
        xch  = rtp(iel,isca(ixch(icla)))
        xash = rtp(iel,isca(inp (icla)))*xmash(icla)
        xmx2   = xmx2 + xch + xck + xash
!
!        Prise en compte de l'humidite
!
        if ( ippmod(iccoal) .eq. 1 ) then
          xmx2 = xmx2+rtp(iel,isca(ixwt(icla)))
        endif
      enddo
!
      if ( xmx2 .gt. 0.d0 ) then
        tfuel(iel) = 0.d0
        do icla=1,nclacp
          ipcte2=ipproc(itemp2(icla))
          tfuel(iel) = tfuel(iel)                                              &
                      +( rtp(iel,isca(ixck(icla)))                             &
                        +rtp(iel,isca(ixch(icla)))                             &
                        +rtp(iel,isca(inp (icla)))*xmash(icla) )               &
                        *propce(iel,ipcte2)
!
!         Prise en compte de l'humidite
!
          if ( ippmod(iccoal) .eq. 1 ) then
            tfuel(iel) = tfuel(iel) + (rtp(iel,isca(ixwt(icla))))              &
                         *propce(iel,ipcte2)
          endif
        enddo
!
        tfuel(iel) = tfuel(iel)/xmx2
!
      else
        tfuel(iel) = propce(iel,ipcte1)
      endif
      tfuelmin = min(tfuel(iel),tfuelmin)
      tfuelmax = max(tfuel(iel),tfuelmax)
    enddo
    if (irangp .ge. 0 ) then
      call parmin(tfuelmin)
      call parmax(tfuelmax)
    endif
    write(nfecra,*) ' Min max de Tfuel pour Hoxy ',tfuelmin,tfuelmax
!
! Combustion heterogene : C + O2 ---> 0.5CO
!
    do iel=1,ncel
!
!   Calcul de HCO(T2)
!
      do ige = 1, ngazem
        coefe(ige) = zero
      enddo
      coefe(ico) = 1.d0
      do icha = 1, ncharm
        f1mc(icha) = zero
        f2mc(icha) = zero
      enddo
        t2        = tfuel(iel)
      mode      = -1
      call cs_coal_htconvers1 &
      !======================
    ( mode  , xhco    , coefe  , f1mc   , f2mc   ,  t2    )
!
!  Calcul de HO2(T1)
!
      do ige = 1, ngazem
        coefe(ige) = zero
      enddo
      coefe(io2) = 1.d0
      do icha = 1, ncharm
        f1mc(icha) = zero
        f2mc(icha) = zero
      enddo
      t1        = propce(iel,ipcte1)
      mode      = -1
      call cs_coal_htconvers1 &
      !======================
      ( mode  , xho2    , coefe  , f1mc   , f2mc   , t1    )
!
      do icla=1,nclacp
        ixckcl = isca(ixck(icla))
        ipcght = ipproc(igmhet(icla))
        if ( rtpa(iel,ixckcl) .gt. epsicp ) then
          gamhet = crom(iel)*propce(iel,ipcght)              &
                   * ( (rtpa(iel,ixckcl))**(2.d0/3.d0) +              &
                  2.d0/3.d0*(rtp(iel,ixckcl)-rtpa(iel,ixckcl))        &
                    /(rtpa(iel,ixckcl))**(1.d0/3.d0) )
        else
          gamhet = 0.d0
        endif
        smbrs(iel) = smbrs(iel)                                        &
                    -gamhet                                            &
                     *(28.d0/12.d0*xhco-16.d0/12.d0*xho2)*volume(iel)
!
      enddo
!
    enddo
!
!   Combustion heterogene : C + CO2 --->2 CO
!
!     Calcul de HO2(T1)
!
    if ( ihtco2 .eq. 1 ) then
      do iel=1,ncel
!
!      Calcul de HCO(T2)
!
        do ige = 1, ngazem
          coefe(ige) = zero
        enddo
        coefe(ico) = 1.d0
        do icha = 1, ncharm
          f1mc(icha) = zero
          f2mc(icha) = zero
        enddo

        t2        = tfuel(iel)
        mode      = -1
        call cs_coal_htconvers1 &
        !======================
        ( mode  , xhco    , coefe  , f1mc   , f2mc   , t2    )

!
!       Calcul de HCO2(T1)
!
        do ige = 1, ngazem
          coefe(ige) = zero
        enddo
        coefe(ico2) = 1.d0
        do icha = 1, ncharm
          f1mc(icha) = zero
          f2mc(icha) = zero
        enddo
        t1        = propce(iel,ipcte1)
        mode      = -1
        call cs_coal_htconvers1 &
        !======================
        ( mode  , xhco2    , coefe  , f1mc   , f2mc   , t1    )
!
        do icla=1,nclacp
          ixckcl  = isca(ixck(icla))
          ipghco2 = ipproc(ighco2(icla))
          if ( rtpa(iel,ixckcl) .gt. epsicp ) then
            gamhet = crom(iel)*propce(iel,ipghco2)              &
                     * ( (rtpa(iel,ixckcl))**(2.d0/3.d0) +               &
                    2.d0/3.d0*(rtp(iel,ixckcl)-rtpa(iel,ixckcl))         &
                      /(rtpa(iel,ixckcl))**(1.d0/3.d0) )
          else
            gamhet = 0.d0
          endif
          smbrs(iel) = smbrs(iel)                                          &
                      -gamhet                                              &
                       *(56.d0/12.d0*xhco-44.d0/12.d0*xhco2) *volume(iel)
!
        enddo
!
      enddo
!
    endif
!
!   Combustion heterogene : C + H2O ---> CO + H2
!
!     Calcul de HO2(T1)
!
    if ( ihth2o .eq. 1 ) then
!
      do iel=1,ncel
!
!      Calcul de HCO(T2)
!
        do ige = 1, ngazem
          coefe(ige) = zero
        enddo
        coefe(ico) = 1.d0
        do icha = 1, ncharm
          f1mc(icha) = zero
          f2mc(icha) = zero
        enddo

        t2        = tfuel(iel)
        mode      = -1
        call cs_coal_htconvers1 &
        !======================
        ( mode  , xhco    , coefe  , f1mc   , f2mc   ,  t2    )

!      Calcul de HH2(T2)

        do ige = 1, ngazem
          coefe(ige) = zero
        enddo
        coefe(ihy) = 1.d0
        do icha = 1, ncharm
          f1mc(icha) = zero
          f2mc(icha) = zero
        enddo

        t2        = tfuel(iel)
        mode      = -1
        call cs_coal_htconvers1 &
        !======================
        ( mode  , xhh2    , coefe  , f1mc   , f2mc   , t2    )
!
!       Calcul de HH2O(T1)
!
        do ige = 1, ngazem
          coefe(ige) = zero
        enddo
        coefe(ih2o) = 1.d0
        do icha = 1, ncharm
          f1mc(icha) = zero
          f2mc(icha) = zero
        enddo
        t1        = propce(iel,ipcte1)
        mode      = -1
        call cs_coal_htconvers1 &
        !======================
      ( mode  , xhh2o    , coefe  , f1mc   , f2mc   , t1    )
!
        do icla=1,nclacp
          ixckcl  = isca(ixck(icla))
          ipghh2o = ipproc(ighh2o(icla))
          if ( rtpa(iel,ixckcl) .gt. epsicp ) then
            gamhet = crom(iel)*propce(iel,ipghh2o)           &
                     * ( (rtpa(iel,ixckcl))**(2.d0/3.d0) +            &
                    2.d0/3.d0*(rtp(iel,ixckcl)-rtpa(iel,ixckcl))      &
                      /(rtpa(iel,ixckcl))**(1.d0/3.d0) )
          else
            gamhet = 0.d0
          endif
          smbrs(iel) = smbrs(iel)                                           &
                      -gamhet                                               &
                       *(28.d0/12.d0*xhco+ 2.d0/12.d0*xhh2                  &
                                         -18.d0/12.d0*xhh2o ) *volume(iel)
!
        enddo
!
      enddo
!
    endif
!
!   Sechage
!
    if ( ippmod(iccoal) .eq. 1 ) then
!
      do icla=1,nclacp
!
        numcha = ichcor(icla)
!
        ipcsec = ipproc(igmsec(icla))
        ipcte2 = ipproc(itemp2(icla))

        do iel = 1, ncel

!          Calcul de H(H2O) a T2

          do ige = 1, ngazem
            coefe(ige) = zero
          enddo
          coefe(ih2o) = 1.d0
          do icha = 1, ncharm
            f1mc(icha) = zero
            f2mc(icha) = zero
          enddo

          t2 = propce(iel,ipcte2)
          mode      = -1
          call cpthp1 &
          !==========
          ( mode  , hh2ov    , coefe  , f1mc   , f2mc   ,t2    )

!         Contribution aux bilans explicite

          if ( rtpa(iel,isca(ixwt(icla))).gt. epsicp .and.          &
               xwatch(numcha) .gt. epsicp       ) then

            aux = crom(iel)*propce(iel,ipcsec)             &
                 *(rtp(iel,isca(ixwt(icla)))/propce(iel,ipcx2c))    &
                 *(1.d0                    /xwatch(numcha))         &
                 *hh2ov

          else
            aux = 0.d0
          endif

          smbrs(iel) = smbrs(iel) - aux*volume(iel)
!
        enddo

      enddo

    endif
!
  endif
endif
!
if ( ieqnox .eq. 1 .and. imdnox.eq.0 .and. ntcabs .gt. 1) then
!
! Termes soures sur Y_HCN et Y_NO
!
  if ( ivar.eq.isca(iyhcn) .or. ivar.eq.isca(iyno) ) then
!
! Pointeur Termes sources
!
    iexp1  = ipproc(ighcn1)
    iexp2  = ipproc(ighcn2)
    iexp3  = ipproc(ignoth)

!
! Masse molaire
!
    wmhcn = wmole(ihcn)
    wmno  = 0.030d0
    wmo2  = wmole(io2  )

    if ( ivar.eq.isca(iyhcn) ) then
!
!        Terme source HCN
!
      if (iwarni(ivar).ge.1) then
        write(nfecra,1000) chaine(1:8)
      endif
!
      do iel=1,ncel
        wmel=propce(iel,ipproc(immel))
        xo2= propce(iel,ipproc(iym1(io2)))*wmel/wmo2
!
        aux = volume(iel)*crom(iel)                      &
             *(propce(iel,iexp2)+propce(iel,iexp1)                &
             *rtpa(iel,isca(iyno))                                &
             *propce(iel,ipproc(immel))                           &
             /wmno)
!
        smbrs(iel)  = smbrs(iel)  - aux*rtpa(iel,ivar)
        rovsdt(iel) = rovsdt(iel) + aux

        do icha=1,ncharb
          gmdev1(icha)=0.d0
          gmdev2(icha)=0.d0
          gmhet (icha)=0.d0
        enddo
!
        do icla=1,nclacp

          icha = ichcor(icla)

          gmdev1(icha) = gmdev1(icha)                              &
               +propce(iel,ipproc(igmdv1(icla)))                   &
               *crom(iel)                                 &
               *rtpa(iel,isca(ixch(icla)))
          gmdev2(icha) = gmdev2(icha)                              &
               +propce(iel,ipproc(igmdv2(icla)))                   &
               *crom(iel)                                 &
               *rtpa(iel,isca(ixch(icla)))
          gmhet(icha) = gmhet(icha)                                &
               +propce(iel,ipproc(igmhet(icla)))                   &
               *crom(iel)                                 &
               *rtpa(iel,isca(ixck(icla)))**(2.d0/3.d0)

        enddo

        do icha=1,ncharb
!        % d'azote sur pur dans le charbon
!
          aux = -volume(iel)*fn(icha)*wmhcn/(wmole(in2)/2.d0)        &
                            *(qpr(icha)*(gmdev1(icha)+gmdev2(icha)))
          if(xo2.gt.0.03d0) then
            aux=aux-volume(iel)*fn(icha)*wmhcn/(wmole(in2)/2.d0)     &
                               * (1.d0-qpr(icha)*y2ch(icha))         &
                                / (1-y2ch(icha))*gmhet(icha)         &
                                * (1.d0-xashch(icha))
          endif
          smbrs(iel)  = smbrs(iel) + aux
        enddo
!
      enddo
!
    endif

    if ( ivar.eq.isca(iyno) ) then
!
!        Terme source NO
!
      if (iwarni(ivar).ge.1) then
        write(nfecra,1000) chaine(1:8)
      endif

      do iel=1,ncel
!
        wmel=propce(iel,ipproc(immel))

        aux1 = volume(iel)*crom(iel)                     &
              *propce(iel,iexp1)*rtpa(iel,isca(iyhcn))            &
              *propce(iel,ipproc(immel))/wmhcn
        aux2 = volume(iel)*crom(iel)                     &
              *propce(iel,iexp2)*rtpa(iel,isca(iyhcn))            &
              *wmno/wmhcn
        aux3 = volume(iel)*crom(iel)**1.5d0              &
              *propce(iel,iexp3)                                  &
              *propce(iel,ipproc(iym1(in2)))

        smbrs(iel)  = smbrs(iel) - aux1*rtpa(iel,ivar)            &
                                 + aux2 + aux3
        rovsdt(iel) = rovsdt(iel) + aux1
      enddo

    endif

  endif

endif
!
if ( ieqnox .eq. 1 .and. imdnox.eq.1 .and. ntcabs .gt. 1) then
!
! Termes soures sur Y_HCN et Y_NO
!
  if ( ivar.eq.isca(iyhcn) .or. ivar.eq.isca(iyno) .or. ivar.eq.isca(iynh3)    &
     ) then
!
! Pointeur Termes sources NO phase gaz
    iexp1  = ipproc(ighcn1)
    iexp2  = ipproc(ighcn2)
    iexp3  = ipproc(ignoth)
    iexp4  = ipproc(ignh31)
    iexp5  = ipproc(ignh32)
    iexprb = ipproc(igrb)
! Pointeur sur CHx1 et CHx2
    ipcyf1 = ipproc(iym1(1))
    ipcyf2 = ipproc(iym1(2))
    idgaz  = ipproc(irom1)
!
! Masse molaire
!
    wmhcn = wmole(ihcn)
    wmno  = 0.030d0
    wmnh3 = wmole(inh3)

    if ( ivar.eq.isca(iyhcn) ) then
!
    aux = 0.d0
!
      do iel = 1,ncel

         propce(iel,ipproc(ifhcnr))  = zero
         propce(iel,ipproc(ifhcnd))  = zero
         propce(iel,ipproc(ifhcnc)) = zero

      enddo
!
!        Terme source HCN
!
      if (iwarni(ivar).ge.1) then
        write(nfecra,1000) chaine(1:8)
      endif
!
      do iel=1,ncel
!
!       Masse molaire de la melange gazeuse
        wmel=propce(iel,ipproc(immel))

!
!       Coefficient des reaction HCN + O2 et HCN + NO
        aux = volume(iel)*crom(iel)                                   &
             *(propce(iel,iexp2)+propce(iel,iexp1)                             &
             *rtpa(iel,isca(iyno))                                             &
             *wmel/wmno)
!
        smbrs(iel)  = smbrs(iel)  - aux*rtpa(iel,ivar)
        rovsdt(iel) = rovsdt(iel) + aux
!
!       Reburning ?
!       Model de Chen
        if(irb.eq.1) then
!
          do icha = 1,ncharb

             ychx = ( propce(iel,ipcyf1) * wmel/wmchx1 )                       &
                  + ( propce(iel,ipcyf2) * wmel/wmchx2 )

             aux = volume(iel)*wmhcn*propce(iel,iexprb)                        &
                 * rtpa(iel,isca(iyno))*wmel/wmno                              &
                 * ychx

             smbrs(iel)  = smbrs(iel)  + aux
!
             propce(iel,ipproc(ifhcnr)) = propce(iel,ipproc(ifhcnr)) + aux
!
          enddo
!
!       Model de Dimitiou
        elseif(irb.eq.2) then
!
          do icha = 1,ncharb

!           Reburning par CHx1
            if(propce(iel,ipcyf1).gt.0.d0) then
!
!             Nombre de point de la discretisation de la temperature
              do ii = 1,7
!
!               On cherche l'intervall teno(ii)<Tgaz<teno(ii+1)
                if(propce(iel,ipcte1).ge.teno(ii).and.propce(iel,ipcte1).lt.   &
                   teno(ii+1)) then

!                 JJ indique le rapport H/C du combustible (4=CH4;3=CH3,etc.)
                  do jj = 1,4
!
!                   On cherche l'intervall jj<chx1(icha)<jj+1
                    if(chx1(icha).ge.4.d0) then

                    core1 = ka(4,ii) + ( (ka(4,ii+1)- ka(4,ii))/(teno(ii+1)-   &
                            teno(ii)) ) * (propce(iel,ipcte1) - teno(ii))
                    core2 = kb(4,ii) + ( (kb(4,ii+1)- kb(4,ii))/(teno(ii+1)-   &
                            teno(ii)) ) * (propce(iel,ipcte1) - teno(ii))
                    para2 = chi2(ii) + ( (chi2(ii+1)-chi2(ii))/(teno(ii+1) -   &
                            teno(ii)) ) * (propce(iel,ipcte1) - teno(ii))

                    elseif(chx1(icha).le.1.d0) then

                    core1 = ka(1,ii) + ( (ka(1,ii+1)- ka(1,ii))/(teno(ii+1)-   &
                            teno(ii)) ) * (propce(iel,ipcte1) - teno(ii))
                    core2 = kb(1,ii) + ( (kb(1,ii+1)- kb(1,ii))/(teno(ii+1)-   &
                            teno(ii)) ) * (propce(iel,ipcte1) - teno(ii))
                    para2 = chi2(ii) + ( (chi2(ii+1)-chi2(ii))/(teno(ii+1) -   &
                            teno(ii)) ) * (propce(iel,ipcte1) - teno(ii))

                    elseif (chx1(icha).ge.jj.and.chx1(icha).lt.jj+1) then

                    core1 = ka(jj,ii) + ( (ka(jj+1,ii+1)- ka(jj,ii))/          &
                            (teno(ii+1)-teno(ii)) ) * (propce(iel,ipcte1)      &
                            - teno(ii))
                    core2 = kb(jj,ii) + ( (kb(jj+1,ii+1)- kb(jj,ii))/          &
                            (teno(ii+1)-teno(ii)) ) * (propce(iel,ipcte1) -    &
                            teno(ii))
                    para2 = chi2(ii) + ( (chi2(ii+1)-chi2(ii))/(teno(ii+1)-    &
                            teno(ii)) ) * (propce(iel,ipcte1) - teno(ii))

                    endif

                  enddo

                endif

              enddo

              if(chx1(icha).ge.3.d0) then

              aux = ( volume(iel)*wmhcn )                                      &
                  * ( (core1 + core2) * para2 )                                &
                  * ( rtp(iel,isca(iyno))*crom(iel)/wmno )                     &
                  * ( propce(iel,ipcyf1)*propce(iel,idgaz)/wmchx1 )

              else

              aux = ( volume(iel)*wmhcn )                                      &
                  * ( core1 + core2 )                                          &
                  * ( rtp(iel,isca(iyno))*crom(iel)/wmno )                     &
                  * ( propce(iel,ipcyf1)*propce(iel,idgaz)/wmchx1 )

              endif
!
              smbrs(iel)  = smbrs(iel)  + aux
!
              propce(iel,ipproc(ifhcnr)) = propce(iel,ipproc(ifhcnr)) + aux
!
            elseif(propce(iel,ipcyf2).gt.0.d0) then

!           Reburning par CHx2
!
!             Nombre de point de la discretisation de la temperature
              do ii = 1,7
!
!               On cherche l'intervall teno(ii)<Tgaz<teno(ii+1)
                if(propce(iel,ipcte1).ge.teno(ii).and.propce(iel,ipcte1).lt.   &
                   teno(ii+1)) then

!                 JJ indique le rapport H/C du combustible (4=CH4;3=CH3,etc.)
                  do jj = 1,4
!
!                   On cherche l'intervall jj<chx1(icha)<jj+1
                    if(chx2(icha).ge.4.d0) then

                    core1 = ka(4,ii) + ( (ka(4,ii+1)- ka(4,ii))/(teno(ii+1)-   &
                            teno(ii)) ) * (propce(iel,ipcte1) - teno(ii))
                    core2 = kb(4,ii) + ( (kb(4,ii+1)- kb(4,ii))/(teno(ii+1)-   &
                            teno(ii)) ) * (propce(iel,ipcte1) - teno(ii))
                    para2 = chi2(ii) + ( (chi2(ii+1)-chi2(ii))/(teno(ii+1) -   &
                            teno(ii)) ) * (propce(iel,ipcte1) - teno(ii))

                    elseif(chx2(icha).le.1.d0) then

                    core1 = ka(1,ii) + ( (ka(1,ii+1)- ka(1,ii))/(teno(ii+1)-   &
                            teno(ii)) ) * (propce(iel,ipcte1) - teno(ii))
                    core2 = kb(1,ii) + ( (kb(1,ii+1)- kb(1,ii))/(teno(ii+1)-   &
                            teno(ii)) ) * (propce(iel,ipcte1) - teno(ii))
                    para2 = chi2(ii) + ( (chi2(ii+1)-chi2(ii))/(teno(ii+1) -   &
                            teno(ii)) ) * (propce(iel,ipcte1) - teno(ii))

                    elseif (chx2(icha).ge.jj.and.chx2(icha).lt.jj+1) then

                    core1 = ka(jj,ii) + ( (ka(jj+1,ii+1)- ka(jj,ii))/          &
                            (teno(ii+1)-teno(ii)) ) * (propce(iel,ipcte1) -    &
                            teno(ii))
                    core2 = kb(jj,ii) + ( (kb(jj+1,ii+1)- kb(jj,ii))/          &
                            (teno(ii+1)-teno(ii)) ) * (propce(iel,ipcte1) -    &
                            teno(ii))
                    para2 = chi2(ii) + ( (chi2(ii+1)-chi2(ii))/(teno(ii+1) -   &
                            teno(ii)) ) * (propce(iel,ipcte1) - teno(ii))

                    endif

                  enddo

                endif

              enddo

              if(chx2(icha).ge.3.d0) then

              aux = ( volume(iel)*wmhcn )                                      &
                  * ( (core1 + core2) * para2 )                                &
                  * ( rtp(iel,isca(iyno))*crom(iel)/wmno )                     &
                  * ( propce(iel,ipcyf2)*propce(iel,idgaz)/wmchx2 )

              else

              aux = ( volume(iel)*wmhcn )                                      &
                  * ( core1 + core2 )                                          &
                  * ( rtp(iel,isca(iyno))*crom(iel)/wmno )                     &
                  * ( propce(iel,ipcyf2)*propce(iel,idgaz)/wmchx2 )

              endif

              smbrs(iel)  = smbrs(iel)  + aux
!
              propce(iel,ipproc(ifhcnr)) = propce(iel,ipproc(ifhcnr)) + aux
!
            endif
!
          enddo
!
        endif
!
!       Initialisation des variables
        do icha=1,ncharb
!
          gmdev1(icha)=0.d0
          gmdev2(icha)=0.d0
          gmhet (icha)=0.d0
!
        enddo
!
        do icla=1,nclacp

          icha   = ichcor(icla)
          ipctem = ipproc(itemp2(icla))
          ixckcl = isca(ixck(icla))
          ipcght = ipproc(igmhet(icla))
!
          mckcl1 = (1.d0-y1ch(icha))*a1ch(icha)                                &
                   *exp(-e1ch(icha)/(rr*propce(iel,ipctem)))

          mckcl2 = (1.d0-y2ch(icha))*a2ch(icha)                                &
                   *exp(-e2ch(icha)/(rr*propce(iel,ipctem)))
!
!         Taux de formation de la premiere reaction de pyrolyse
          gmdev1(icha) = gmdev1(icha)                                          &
               +propce(iel,ipproc(igmdv1(icla)))                               &
               *crom(iel)                                             &
               *rtpa(iel,isca(ixch(icla)))
!
!         Taux de formation de la deuxieme reaction de pyrolyse
          gmdev2(icha) = gmdev2(icha)                                          &
               +propce(iel,ipproc(igmdv2(icla)))                               &
               *crom(iel)                                             &
               *rtpa(iel,isca(ixch(icla)))
!
          if ( rtpa(iel,ixckcl) .gt. epsicp ) then
!           Taux de reaction de la combustion heterogene
            gmhet(icha) = gmhet(icha)                                          &
               +propce(iel,ipcght)                                             &
               *crom(iel)                                             &
               *( rtpa(iel,ixckcl)*((1.d0/(mckcl2/mckcl1+1.d0))*yhcnc1(icha)   &
                   +(1.d0/(mckcl1/mckcl2+1.d0))*yhcnc2(icha)) )**(2.d0/3.d0)
          endif
!
        enddo
!
!       Terme source modifie (nouveau model de NOx)
!
        do icha=1,ncharb

!          Liberation du HCN lors de la devotalisation
           aux = -volume(iel)*(gmdev1(icha)*yhcnle(icha)                       &
                 +gmdev2(icha)*yhcnlo(icha))
!
!          Liberation du HCN lors de la reaction heterogene selon la valeur
!          repnck(icha)

           aux = aux-volume(iel)*gmhet(icha)
!
           smbrs(iel)  = smbrs(iel) + aux

!          Affichage des termes sources
           propce(iel,ipproc(ifhcnd)) = propce(iel,ipproc(ifhcnd))             &
           -volume(iel)*(gmdev1(icha)*yhcnle(icha)+gmdev2(icha)*yhcnlo(icha))
!
           propce(iel,ipproc(ifhcnc))= propce(iel,ipproc(ifhcnc))              &
           -volume(iel)*gmhet(icha)

        enddo

      enddo
!
    endif
!
!
    if ( ivar.eq.isca(iynh3)) then
!
       aux = 0.d0
!
      do iel = 1, ncel

         propce(iel,ipproc(ifnh3d)) = zero

      enddo
!
!     Terme source NH3
!
      if (iwarni(ivar).ge.1) then
        write(nfecra,1000) chaine(1:8)
      endif
!
      do iel=1,ncel
!
!       Masse molaire de la melange gazeuse
        wmel = propce(iel,ipproc(immel))
!
!       Coefficient des reaction NH3 + O2 et NH3 + NO
        aux  =   volume(iel)*crom(iel)                                &
             * ( propce(iel,iexp4) + propce(iel,iexp5)                         &
             *   rtpa(iel,isca(iyno))*wmel/wmno )
!
        smbrs(iel)  = smbrs(iel)  - aux*rtpa(iel,ivar)
        rovsdt(iel) = rovsdt(iel) + aux
!
!       Initialisation des variables
        do icha=1,ncharb

          gmdev1(icha)=0.d0
          gmdev2(icha)=0.d0
          gmhet (icha)=0.d0

        enddo
!
        do icla=1,nclacp

          icha = ichcor(icla)
!
!         Taux de formation de la premiere reaction de pyrolyse
          gmdev1(icha) = gmdev1(icha)                                          &
               +propce(iel,ipproc(igmdv1(icla)))                               &
               *crom(iel)                                             &
               *rtpa(iel,isca(ixch(icla)))
!
!         Taux de formation de la deuxieme reaction de pyrolyse
          gmdev2(icha) = gmdev2(icha)                                          &
               +propce(iel,ipproc(igmdv2(icla)))                               &
               *crom(iel)                                             &
               *rtpa(iel,isca(ixch(icla)))
!
        enddo
!
        do icha=1,ncharb
!
!        Liberation du NH3 lors de la devolatisation.
           aux = -volume(iel)*(gmdev1(icha)*ynh3le(icha)                       &
                 +gmdev2(icha)*ynh3lo(icha))
!
           smbrs(iel)  = smbrs(iel) + aux

!        Affichage des termes source
           propce(iel,ipproc(ifnh3d)) = propce(iel,ipproc(ifnh3d))             &
           -volume(iel)*(gmdev1(icha)*ynh3le(icha)+gmdev2(icha)*ynh3lo(icha))
!
           propce(iel,ipproc(ifnh3c)) = zero
!
        enddo



      enddo
!
    endif
!
    if ( ivar.eq.isca(iyno) ) then
!
      do iel = 1 ,ncel

        propce(iel,ipproc(icnorb))  = zero
        propce(iel,ipproc(ifnoch)) = zero

      enddo
!
!     Terme source NO
!
      if (iwarni(ivar).ge.1) then
        write(nfecra,1000) chaine(1:8)
      endif

      do iel=1,ncel
!
!       Masse molaire de la melange gazeuse
        wmel=propce(iel,ipproc(immel))
!
!       Coefficient de la reaction HCN + NO
        aux1 = volume(iel)*crom(iel)                                  &
              *propce(iel,iexp1)*rtpa(iel,isca(iyhcn))                         &
              *wmel/wmhcn
!
        propce(iel,ipproc(icnohc)) = aux1*rtpa(iel,ivar)
!
!       Coefficient de la reaction HCN + O2
        aux2 = volume(iel)*crom(iel)                                  &
              *propce(iel,iexp2)*rtpa(iel,isca(iyhcn))                         &
              *wmno/wmhcn
!
        propce(iel,ipproc(ifnohc)) = aux2
!
!       Coefficient du NO thermique
        aux3 = volume(iel)*crom(iel)**1.5d0                           &
              *propce(iel,iexp3)                                               &
!       Pourquoi la fraction massique d'azote n'a ete pas transforme dans une
!       fraction molaire ?
              *propce(iel,ipproc(iym1(in2)))
!
        propce(iel,ipproc(ifnoth)) = aux3
!
!       Coefficient de la reaction NH3 + O2 --> NO + ...
        aux4 = volume(iel)*crom(iel)                                  &
              *propce(iel,iexp4)*rtpa(iel,isca(iynh3))                         &
              *wmno/wmnh3
!
!
        propce(iel,ipproc(ifnonh)) = aux4
!
!       Coefficient de la reaction NH3 + NO --> N2 + ...
        aux5 = volume(iel)*crom(iel)                                  &
              *propce(iel,iexp5)*rtpa(iel,isca(iynh3))                         &
              *wmel/wmnh3
!
!
        propce(iel,ipproc(icnonh)) = aux5*rtpa(iel,ivar)

!       Reburning ?
!       Model de Chen
        if(irb.eq.1) then
!
          do icha = 1,ncharb

             ychx = ( propce(iel,ipcyf1) * wmel/wmchx1 )                       &
                  + ( propce(iel,ipcyf2) * wmel/wmchx2 )

             aux = volume(iel)*wmhcn*propce(iel,iexprb)                        &
                 * rtpa(iel,isca(iyno)) * wmel/wmno  * ychx

             smbrs(iel)  = smbrs(iel)  - aux
!
             propce(iel,ipproc(icnorb)) = propce(iel,ipproc(icnorb)) + aux

          enddo

!       Model de Dimitiou
        elseif(irb.eq.2) then

          do icha = 1,ncharb
!
!           Reburning par CHx1
            if (propce(iel,ipcyf1).gt.0.d0) then

!             Nombre de point de la discretisation de la temperature
              do ii = 1,7
!
!               On cherche l'intervall teno(ii)<Tgaz<teno(ii+1)
                if(propce(iel,ipcte1).ge.teno(ii).and.propce(iel,ipcte1).lt.   &
                   teno(ii+1)) then

!                 JJ indique le rapport H/C du combustible (4=CH4;3=CH3,etc.)
                  do jj = 1,4
!
!                   On cherche l'intervall jj<chx1(icha)<jj+1
                    if(chx1(icha).ge.4.d0) then

                    core1 = ka(4,ii) + ( (ka(4,ii+1)- ka(4,ii))/(teno(ii+1)-   &
                            teno(ii)) ) * (propce(iel,ipcte1) - teno(ii))
                    core2 = kb(4,ii) + ( (kb(4,ii+1)- kb(4,ii))/(teno(ii+1)-   &
                            teno(ii)) ) * (propce(iel,ipcte1) - teno(ii))
                    core3 = kc(4,ii) + ( (kc(4,ii+1)- kc(4,ii))/(teno(ii+1)-   &
                            teno(ii)) ) * (propce(iel,ipcte1) - teno(ii))
                    para2 = chi2(ii) + ( (chi2(ii+1)-chi2(ii))/(teno(ii+1) -   &
                            teno(ii)) ) * (propce(iel,ipcte1) - teno(ii))

                    elseif(chx1(icha).le.1.d0) then

                    core1 = ka(1,ii) + ( (ka(1,ii+1)- ka(1,ii))/               &
                            (teno(ii+1)-teno(ii)) ) * (propce(iel,ipcte1)  -   &
                            teno(ii))
                    core2 = kb(1,ii) + ( (kb(1,ii+1)- kb(1,ii))/(teno(ii+1)-   &
                            teno(ii)) ) * (propce(iel,ipcte1) - teno(ii))
                    core3 = kc(1,ii) + ( (kc(1,ii+1)- kc(1,ii))/               &
                            (teno(ii+1)-teno(ii)) ) * (propce(iel,ipcte1) -    &
                            teno(ii))
                    para2 = chi2(ii) + ( (chi2(ii+1)-chi2(ii)) /               &
                            (teno(ii+1) - teno(ii)) ) * (propce(iel,ipcte1) -  &
                            teno(ii))

                    elseif (chx1(icha).ge.jj.and.chx1(icha).lt.jj+1) then

                    core1 = ka(jj,ii) + ( (ka(jj+1,ii+1)- ka(jj,ii))/          &
                            (teno(ii+1)-teno(ii)) ) * (propce(iel,ipcte1) -    &
                            teno(ii))
                    core2 = kb(jj,ii) + ( (kb(jj+1,ii+1)- kb(jj,ii))/          &
                            (teno(ii+1)-teno(ii)) ) * (propce(iel,ipcte1) -    &
                            teno(ii))
                    core3 = kc(jj,ii) + ( (kc(jj+1,ii+1)- kc(jj,ii))/          &
                            (teno(ii+1)-teno(ii)) ) * (propce(iel,ipcte1) -    &
                            teno(ii))
                    para2 = chi2(ii) + ( (chi2(ii+1)-chi2(ii))/                &
                            (teno(ii+1) - teno(ii)) ) * (propce(iel,ipcte1) -  &
                            teno(ii))

                    endif

                  enddo

                endif

              enddo

              if(chx1(icha).ge.3.d0) then

              auxrb1 = ( volume(iel)*wmno )                                    &
                     * ( (core1 + core2 + core3) * para2 )                     &
                     * ( propce(iel,ipcyf1)*propce(iel,idgaz)/wmchx1 )         &
                     * ( rtp(iel,isca(iyno)) * crom(iel)/wmno )

              else

              auxrb1 = ( volume(iel)*wmno )                                    &
                     * ( core1 + core2 + core3 )                               &
                     * ( propce(iel,ipcyf1)*propce(iel,idgaz)/wmchx1 )         &
                     * ( rtp(iel,isca(iyno)) * crom(iel)/wmno )

              endif

              smbrs(iel)  = smbrs(iel)  - auxrb1
!
              propce(iel,ipproc(icnorb)) = propce(iel,ipproc(icnorb)) + auxrb1

!           Reburning par CHx2
            elseif(propce(iel,ipcyf2).gt.0.d0) then
!
!             Nombre de point de la discretisation de la temperature
              do ii = 1,7
!
!               On cherche l'intervall teno(ii)<Tgaz<teno(ii+1)
                if(propce(iel,ipcte1).ge.teno(ii).and.propce(iel,ipcte1).lt.   &
                   teno(ii+1)) then

!                 JJ indique le rapport H/C du combustible (4=CH4;3=CH3,etc.)
                  do jj = 1,4
!
!                   On cherche l'intervall jj<chx1(icha)<jj+1
                    if(chx2(icha).ge.4.d0) then

                    core1 = ka(4,ii) + ( (ka(4,ii+1)- ka(4,ii))/               &
                            (teno(ii+1)-teno(ii)) ) * (propce(iel,ipcte1) -    &
                            teno(ii))
                    core2 = kb(4,ii) + ( (kb(4,ii+1)- kb(4,ii))/               &
                            (teno(ii+1)-teno(ii)) ) * (propce(iel,ipcte1) -    &
                            teno(ii))
                    core3 = kc(4,ii) + ( (kc(4,ii+1)- kc(4,ii))/               &
                            (teno(ii+1)-teno(ii)) ) * (propce(iel,ipcte1) -    &
                            teno(ii))
                    para2 = chi2(ii) + ( (chi2(ii+1)-chi2(ii))/                &
                            (teno(ii+1) - teno(ii)) ) * (propce(iel,ipcte1) -  &
                            teno(ii))

                    elseif(chx2(icha).le.1.d0) then

                    core1 = ka(1,ii) + ( (ka(1,ii+1)- ka(1,ii))/(teno(ii+1)-   &
                            teno(ii)) ) * (propce(iel,ipcte1) -                &
                            teno(ii))
                    core2 = kb(1,ii) + ( (kb(1,ii+1)- kb(1,ii))/               &
                            (teno(ii+1)-teno(ii)) ) * (propce(iel,ipcte1) -    &
                            teno(ii))
                    core3 = kc(1,ii) + ( (kc(1,ii+1)- kc(1,ii))/               &
                            (teno(ii+1)-teno(ii)) ) * (propce(iel,ipcte1) -    &
                            teno(ii))
                    para2 = chi2(ii) + ( (chi2(ii+1)-chi2(ii)) /               &
                            (teno(ii+1) - teno(ii)) ) * (propce(iel,ipcte1) -  &
                            teno(ii))
                    elseif (chx2(icha).ge.jj.and.chx2(icha).lt.jj+1) then

                    core1 = ka(jj,ii) + ( (ka(jj+1,ii+1)- ka(jj,ii))/          &
                            (teno(ii+1)-teno(ii)) ) * (propce(iel,ipcte1) -    &
                            teno(ii))
                    core2 = kb(jj,ii) + ( (kb(jj+1,ii+1)- kb(jj,ii))/          &
                            (teno(ii+1)-teno(ii)) ) * (propce(iel,ipcte1) -    &
                            teno(ii))
                    core3 = kc(jj,ii) + ( (kc(jj+1,ii+1)- kc(jj,ii))/          &
                            (teno(ii+1)-teno(ii)) ) * (propce(iel,ipcte1) -    &
                            teno(ii))
                    para2 = chi2(ii) + ( (chi2(ii+1)-chi2(ii))/                &
                            (teno(ii+1) - teno(ii)) ) * (propce(iel,ipcte1) -  &
                            teno(ii))

                    endif

                  enddo

                endif

              enddo

              if(chx2(icha).ge.3.d0) then

              auxrb2 = ( volume(iel)*wmno )                                    &
                     * ( (core1 + core2 + core3) * para2 )                     &
                     * ( propce(iel,ipcyf2)*propce(iel,idgaz)/wmchx2 )         &
                     * ( rtp(iel,isca(iyno)) * crom(iel)/wmno )

              else

              auxrb2 = ( volume(iel)*wmno )                                    &
                     * ( core1 + core2 + core3 )                               &
                     * ( propce(iel,ipcyf2)*propce(iel,idgaz)/wmchx2 )         &
                     * ( rtp(iel,isca(iyno)) * crom(iel)/wmno )

              endif

              smbrs(iel)  = smbrs(iel)  - auxrb2
!
              propce(iel,ipproc(icnorb)) = propce(iel,ipproc(icnorb)) + auxrb2
!
            endif
!
          enddo

        endif

!
!       Initialisation
        do icha=1,ncharb
!
          gmhet (icha)=0.d0
!
        enddo
!
        do icla=1,nclacp
!
          icha   = ichcor(icla)
          ipctem = ipproc(itemp2(icla))
          ixckcl = isca(ixck(icla))
!
          mckcl1 = (1.d0-y1ch(icha))*a1ch(icha)                                &
                   *exp(-e1ch(icha)/(rr*propce(iel,ipctem)))

          mckcl2 = (1.d0-y2ch(icha))*a2ch(icha)                                &
                   *exp(-e2ch(icha)/(rr*propce(iel,ipctem)))
!
!         Taux de reaction de la combustion heterogene
          if ( rtpa(iel,ixckcl) .gt. epsicp ) then
          gmhet(icha) = gmhet(icha)                                            &
               +propce(iel,ipproc(igmhet(icla)))                               &
               *crom(iel)                                             &
               *( rtpa(iel,ixckcl)*((1.d0/(mckcl2/mckcl1+1.d0))*ynoch1(icha)   &
                   +(1.d0/(mckcl1/mckcl2+1.d0))*ynoch2(icha)) )**(2.d0/3.d0)
          endif
!
        enddo
!
!       Coefficient de NO libere lors de la combustion heterogene
!
        do icha=1,ncharb

          auxhet = -volume(iel)*gmhet(icha)
          propce(iel,ipproc(ifnoch)) = propce(iel,ipproc(ifnoch)) + auxhet


          smbrs(iel)  = smbrs(iel) - aux1*rtpa(iel,ivar)                       &
                                   - aux5*rtpa(iel,ivar)                       &
                                   + aux2 + aux3 + aux4 + auxhet

          rovsdt(iel) = rovsdt(iel) + aux1 + aux5

        enddo

      enddo

    endif

  endif

endif

!--------
! Formats
!--------

 1000 format(' TERMES SOURCES PHYSIQUE PARTICULIERE POUR LA VARIABLE '  &
       ,a8,/)

!----
! End
!----

!===============================================================================
! Deallocation dynamic arrays
!----
deallocate(w1,w2,w3,w4,w5,STAT=iok1)
!----
if ( iok1 > 0 ) THEN
  write(nfecra,*) ' Memory deallocation error inside: '
  write(nfecra,*) '      cs_coal_scast                '
  call csexit(1)
endif
deallocate(tfuel,STAT=iok1)
if ( iok1 > 0 ) THEN
  write(nfecra,*) ' Memory deallocation error inside: '
  write(nfecra,*) '      cs_coal_scast                '
  call csexit(1)
endif
!===============================================================================

return

end subroutine
