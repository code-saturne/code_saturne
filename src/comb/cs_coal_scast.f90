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

subroutine cs_coal_scast &
!=======================

 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   iscal  ,                                                       &
   itypfb ,                                                       &
   icepdc , icetsm , itypsm ,                                     &
   izfppp ,                                                       &
   dt     , rtpa   , rtp    , propce , propfa , propfb ,          &
   coefa  , coefb  , ckupdc , smacel ,                            &
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
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! ncepdp           ! i  ! <-- ! number of cells with head loss                 !
! ncesmp           ! i  ! <-- ! number of cells with mass source term          !
! iscal            ! i  ! <-- ! scalar number                                  !
! itypfb(nfabor    ! te ! --> ! type des faces de bord                         !
! icepdc(ncelet    ! te ! <-- ! numero des ncepdp cellules avec pdc            !
! icetsm(ncesmp    ! te ! <-- ! numero des cellules a source de masse          !
! itypsm           ! te ! <-- ! type de source de masse pour les               !
! (ncesmp,nvar)    !    !     !  variables (cf. ustsma)                        !
! izfppp           ! te ! --> ! numero de zone de la face de bord              !
! (nfabor)         !    !     !  pour le module phys. part.                    !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
! ckupdc           ! tr ! <-- ! tableau de travail pour pdc                    !
!  (ncepdp,6)      !    !     !                                                !
! smacel           ! tr ! <-- ! valeur des variables associee a la             !
! (ncesmp,*   )    !    !     !  source de masse                               !
!                  !    !     !  pour ivar=ipr, smacel=flux de masse           !
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

use dimens, only: ndimfb
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

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal
integer          ncepdp , ncesmp
integer          iscal

integer          itypfb(nfabor)
integer          icepdc(ncepdp)
integer          icetsm(ncesmp), itypsm(ncesmp,nvar)
integer          izfppp(nfabor)

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(ndimfb,*)
double precision coefa(ndimfb,*), coefb(ndimfb,*)
double precision ckupdc(ncepdp,6), smacel(ncesmp,nvar)
double precision smbrs(ncelet), rovsdt(ncelet)

! Local variables

character*80     chaine
integer          ivar , ipcrom, iel
integer          numcla , numcha , icla , numtra
integer          ipcgch , ipcgd1 , ipcgd2 , ipcght , ipcsec
integer          ipghco2 , ipghh2o
integer          ixchcl , ixckcl , iscala
integer          ipcro2 , ipcte1 , ipcte2 , ipcvsl , ipccp
integer          ipcdia , ipcvst
integer          mode, ige
integer          ifac   , ifinra , icoefa , icoefb
integer          ipcx2c , icha , imode , ii
integer          itermx,nbpauv,nbrich,nbepau,nberic,ipghc2
integer          iterch,nbpass,nbarre,nbimax
integer          iexp1 , iexp2 , iexp3

double precision epsrgp , climgp , extrap , xnuss
double precision aux, rhovst
double precision coefe(ngazem)
double precision t1, t2, hh2ov
double precision f1mc(ncharm), f2mc(ncharm)
double precision xhdev1 , xhdev2 , xhco , xho2 , gamdv1 , gamdv2
double precision xhco2  , xhh2   , xhh2o

double precision gamhet , den
double precision xxco,xxo2,xxco2,xxh2o,xco2mx
double precision xkp,xk0p,xkm,xk0m,wplus,wmoins,t0p,t0m
double precision auxp,auxm
double precision anmr,xcot,xo2t,xco2e,xo2e,xcoe,tauchi,tautur
double precision sqh2o , x2
double precision err1mx,err2mx
double precision errch,ter1,ddelta,xden
double precision fn0,fn1,fn2,anmr0,anmr1,anmr2
double precision lnk0p,l10k0e,lnk0m,t0e,xco2eq,xcoeq,xo2eq
double precision xcom,xo2m,xkcequ,xkpequ

double precision  xw1,xw2,xw3,xw4
double precision  xo2,wmel,wmhcn,wmno,wmo2
double precision gmdev1(ncharm),gmdev2(ncharm),gmhet(ncharm)
double precision aux1 , aux2 , aux3
double precision xch,xck,xash,xmx2
double precision tfuelmin,tfuelmax

integer           iok1
double precision , dimension ( : )     , allocatable :: w1,w2,w3,w4,w5
double precision , dimension ( : )     , allocatable :: tfuel

!===============================================================================
!
!===============================================================================
! 1. INITIALISATION
!===============================================================================
! --- Numero du scalaire a traiter : ISCAL
! --- Numero de la variable associee au scalaire a traiter ISCAL
ivar = isca(iscal)

! --- Nom de la variable associee au scalaire a traiter ISCAL
chaine = nomvar(ipprtp(ivar))

! --- Numero des grandeurs physiques
ipcrom = ipproc(irom)
ipcvst = ipproc(ivisct)
ipcte1 = ipproc(itemp1)
!
!===============================================================================
! Deallocation dynamic arrays
!----
allocate(w1(1:ncel),w2(1:ncel),w3(1:ncel),w4(1:ncel),w5(1:ncel),STAT=iok1)
if ( iok1 > 0 ) THEN
  write(nfecra,*) ' Memory allocation error inside: '
  write(nfecra,*) '     cs_coal_scast               '
  call csexit(1)
endif
allocate(tfuel(1:ncel),STAT=iok1)
if ( iok1 > 0 ) THEN
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

    xw1 = - propce(iel,ipcrom)*propce(iel,ipcgch)*volume(iel)

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

    xw1=-propce(iel,ipcrom)*rtp(iel,ixchcl)*propce(iel,ipcgch)     &
                                           *volume(iel)

! AE : On prend RTP(IEL,IXCHCL) et pas RTPA(IEL,IXCHCL) afin
!      d'etre conservatif sur la masse

! ---- Calcul de W2 = rho.Xch.(GMDV1+GMDV2)Volume < 0

    xw2 = propce(iel,ipcrom)*rtp(iel,ixchcl)                  &
         *(propce(iel,ipcgd1)+propce(iel,ipcgd2))*volume(iel)

    if ( rtpa(iel,ixckcl) .gt. epsicp ) then

! Reaction C(s) + O2 ---> 0.5CO
! =============================

! ---- Calcul de la partie implicite  > 0 du TS relatif a GMHET

      xw3 = -2.d0/3.d0*propce(iel,ipcrom)*propce(iel,ipcght) &
           /(rtpa(iel,ixckcl))**(1.d0/3.d0)*volume(iel)

! ---- Calcul de la partie explicite < 0 du TS relatif a GMHET

      xw4 = propce(iel,ipcrom)*propce(iel,ipcght)             &
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

        xw3 = -2.d0/3.d0*propce(iel,ipcrom)*propce(iel,ipghco2)     &
              /(rtpa(iel,ixckcl))**(1.d0/3.d0)*volume(iel)

! ---- Calcul de la partie explicite < 0 du TS relatif a GMHET

        xw4 = propce(iel,ipcrom)*propce(iel,ipghco2)                 &
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

        xw3 = -2.d0/3.d0*propce(iel,ipcrom)*propce(iel,ipghh2o)     &
              /(rtpa(iel,ixckcl))**(1.d0/3.d0)*volume(iel)

! ---- Calcul de la partie explicite < 0 du TS relatif a GMHET

        xw4 = propce(iel,ipcrom)*propce(iel,ipghh2o)                &
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
       xw1 = propce(iel,ipcrom)*propce(iel,ipcsec)*volume(iel)    &
            *(1.d0/propce(iel,ipcx2c))*(1.d0/xwatch(numcha))

       rovsdt(iel) = rovsdt(iel) + max(xw1,zero)
       smbrs(iel)  = smbrs(iel)  - xw1*rtpa(iel,ivar)
     endif

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
    if ( ivisls(ihm).gt.0 ) then
      ipcvsl = ipproc(ivisls(ihm))
      if ( icp.gt.0 ) then
        ipccp   = ipproc(icp)
        w1(iel) = propce(iel,ipcvsl) * propce(iel,ipccp)
      else
        w1(iel) = propce(iel,ipcvsl) * cp0
      endif
    else
      if ( icp.gt.0 ) then
        ipccp   = ipproc(icp)
        w1(iel) = visls0(ihm) * propce(iel,ipccp)
      else
        w1(iel) = visls0(ihm) * cp0
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
                  / propce(iel,ipcro2) * propce(iel,ipcrom)       &
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

    gamdv1 = propce(iel,ipcrom)*rtp(iel,ixchcl)                   &
            *propce(iel,ipcgd1)

    gamdv2 = propce(iel,ipcrom)*rtp(iel,ixchcl)                   &
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

      gamhet = propce(iel,ipcrom)*propce(iel,ipcght)              &
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

        gamhet = propce(iel,ipcrom)*propce(iel,ipghco2)           &
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

        gamhet = propce(iel,ipcrom)*propce(iel,ipghh2o)           &
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

        aux = -propce(iel,ipcrom)*propce(iel,ipcsec)              &
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
        w1(iel) = w1(iel) - propce(iel,ipcrom)*rtp(iel,ixchcl)    &
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
        w1(iel) = w1(iel) - propce(iel,ipcrom)*rtp(iel,ixchcl)    &
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
                 - propce(iel,ipcrom)*propce(iel,ipcght)          &
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
                   - propce(iel,ipcrom)*propce(iel,ipcght)         &
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
                   - propce(iel,ipcrom)*propce(iel,ipcght)        &
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
 ( nvar   , nscal  , ncepdp , ncesmp ,                             &
   iscal  ,                                                        &
   itypfb ,                                                        &
   icepdc , icetsm , itypsm ,                                      &
   dt     , rtpa   , rtp    , propce , propfa , propfb ,           &
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
      + propce(iel,ipcrom)*propce(iel,ipcsec)                     &
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
         * volume(iel) * propce(iel,ipcrom)
!
       w1(iel) = volume(iel)*propce(iel,ipcrom)/(tauchi+tautur)
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
              + propce(iel,ipcrom)*propce(iel,ipghc2)             &
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
          tfuel(iel) = tfuel(iel)                                                  &
                      +( rtp(iel,isca(ixck(icla)))                                 &
                        +rtp(iel,isca(ixch(icla)))                                 &
                        +rtp(iel,isca(inp (icla)))*xmash(icla) )*propce(iel,ipcte2)
!
!         Prise en compte de l'humidite
!
          if ( ippmod(iccoal) .eq. 1 ) then
            tfuel(iel) = tfuel(iel) +(rtp(iel,isca(ixwt(icla))))*propce(iel,ipcte2)
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
          gamhet = propce(iel,ipcrom)*propce(iel,ipcght)              &
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
            gamhet = propce(iel,ipcrom)*propce(iel,ipghco2)              &
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

        t2        = propce(iel,ipcte2)
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
            gamhet = propce(iel,ipcrom)*propce(iel,ipghh2o)           &
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
        ipcx2c = ipproc(ix2(icla))

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

            aux = propce(iel,ipcrom)*propce(iel,ipcsec)             &
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
if ( ieqnox .eq. 1 .and. ntcabs .gt. 1) then
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
        aux = volume(iel)*propce(iel,ipcrom)                      &
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
               *propce(iel,ipcrom)                                 &
               *rtpa(iel,isca(ixch(icla)))
          gmdev2(icha) = gmdev2(icha)                              &
               +propce(iel,ipproc(igmdv2(icla)))                   &
               *propce(iel,ipcrom)                                 &
               *rtpa(iel,isca(ixch(icla)))
          gmhet(icha) = gmhet(icha)                                &
               +propce(iel,ipproc(igmhet(icla)))                   &
               *propce(iel,ipcrom)                                 &
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

        aux1 = volume(iel)*propce(iel,ipcrom)                     &
              *propce(iel,iexp1)*rtpa(iel,isca(iyhcn))            &
              *propce(iel,ipproc(immel))/wmhcn
        aux2 = volume(iel)*propce(iel,ipcrom)                     &
              *propce(iel,iexp2)*rtpa(iel,isca(iyhcn))            &
              *wmno/wmhcn
        aux3 = volume(iel)*propce(iel,ipcrom)**1.5d0              &
              *propce(iel,iexp3)                                  &
              *propce(iel,ipproc(iym1(in2)))

        smbrs(iel)  = smbrs(iel) - aux1*rtpa(iel,ivar)            &
                                 + aux2 + aux3
        rovsdt(iel) = rovsdt(iel) + aux1
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
