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

subroutine cs_fuel_scast &
!=======================

 ( iscal  ,                                                       &
   rtpa   , rtp    , propce ,                                     &
   smbrs  , rovsdt )

!===============================================================================
! FONCTION :
! ----------

! ROUTINE PHYSIQUE PARTICULIERE : FLAMME FUEL
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
use dimens, only: nvar
use cstphy
use cstnum
use parall
use period
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

integer          iscal

double precision rtp(ncelet,nflown:nvar), rtpa(ncelet,nflown:nvar)
double precision propce(ncelet,*)
double precision smbrs(ncelet), rovsdt(ncelet)

! Local variables

character*80     chaine
integer          ivar , iel, icla , numcla
integer          iexp1 , iexp2 , iexp3
integer          ipcro2 , ipcte1 , ipcte2
integer          ipcdia
integer          imode  , iesp
integer          ipcgev , ipcght , ipchgl
integer          itermx,nbpauv,nbrich,nbepau,nberic
integer          nbarre,nbimax,nbpass
integer          iterch

double precision aux, rhovst
double precision rom
double precision hfov
double precision ho2,hco,xesp(ngazem),t2mt1
double precision gmech,gmvap,gmhet
double precision xxco,xxo2,xxco2,xxh2o
double precision xkp,xkm,t0p,t0m
double precision aux1 , aux2 , aux3 , w1
double precision anmr,tauchi,tautur
double precision sqh2o , x2 , wmhcn , wmno ,wmo2
double precision err1mx,err2mx
double precision errch,fn,qpr
double precision auxmax,auxmin
double precision ymoy
double precision fn0,fn1,fn2,anmr0,anmr1,anmr2
double precision lnk0p,l10k0e,lnk0m,t0e,xco2eq,xcoeq,xo2eq
double precision xcom,xo2m,xkcequ,xkpequ,xden
double precision, dimension(:), pointer ::  crom
double precision, dimension(:), pointer :: cka, cvara_ep

!===============================================================================
! 1. INITIALISATION
!===============================================================================

! --- Numero du scalaire a traiter : ISCAL

! --- Numero de la variable associee au scalaire a traiter ISCAL
ivar = isca(iscal)

! --- Nom de la variable associee au scalaire a traiter ISCAL
call field_get_label(ivarfl(ivar), chaine)

! --- Numero des grandeurs physiques (voir cs_user_boundary_conditions)
call field_get_val_s(icrom, crom)

call field_get_val_prev_s(ivarfl(ik), cka)
call field_get_val_prev_s(ivarfl(iep), cvara_ep)

! --- Temperature phase gaz

ipcte1 = ipproc(itemp1)

!===============================================================================
! 2. PRISE EN COMPTE DES TERMES SOURCES POUR LES VARIABLES RELATIVES
!    AUX CLASSES DE PARTICULES
!===============================================================================

! --> Terme source pour l'enthalpie du liquide

if ( ivar .ge. isca(ih2(1))     .and.                            &
     ivar .le. isca(ih2(nclafu))      ) then

  if (iwarni(ivar).ge.1) then
    write(nfecra,1000) chaine(1:8)
  endif

  numcla = ivar-isca(ih2(1))+1

  ipcro2 = ipproc(irom2 (numcla))
  ipcdia = ipproc(idiam2(numcla))
  ipcte2 = ipproc(itemp2(numcla))
  ipcgev = ipproc(igmeva(numcla))
  ipcght = ipproc(igmhtf(numcla))
  ipchgl = ipproc(ih1hlf(numcla))

!       La variable est l'enthalpie du liquide ramenée à
!       la masse de mélange
!       Les flux interfaciaux contribuent à faire varier l'enthalpie
!       du liquide
!       La vapeur emporte son enthalpie
!       flux = PROPCE(IEL,IPPROC(IGMEVA))
!       enthalpie massique reconstituée à partir de EHGAZE(IFOV )
!       à la température de la goutte
!       L'oxydation héterogène comporte un flux "entrant" d'O2
!       un flux sortant de CO
!       le flux net est celui du carbone
!       fluxIN  = 16/12 * PROPCE(IEL,IPPROC(IGMHTF))
!       fluxOUT = 28/12 * PROPCE(IEL,IPPROC(IGMHTF))
!       Enthalpie entrante reconstituée à partir de EHGAZE(IO2 )
!       à la température du gaz environnant
!       Enthalpie sortante reconstituée à partir de EHGAZE(ICO )
!       à la température du grain

  imode = -1
  do iel = 1, ncel
!
    if ( rtpa(iel,isca(iyfol(numcla))) .gt. epsifl ) then
!
      rom = crom(iel)

      do iesp = 1, ngazem
        xesp(iesp) = zero
      enddo

      xesp(ifov) = 1.d0
      call cs_fuel_htconvers1(imode,hfov,xesp,propce(iel,ipcte2))

      xesp(ifov) = zero
      xesp(io2)  = 1.d0
      call cs_fuel_htconvers1(imode,ho2 ,xesp,propce(iel,ipcte1))

      xesp(io2)  = zero
      xesp(ico)  = 1.d0
      call cs_fuel_htconvers1(imode,hco,xesp,propce(iel,ipcte2))

      t2mt1 = propce(iel,ipcte2)-propce(iel,ipcte1)

      gmech = -propce(iel,ipchgl)*t2mt1
      gmvap = propce(iel,ipcgev)*hfov*t2mt1
      gmhet = 16.d0/12.d0*propce(iel,ipcght)*ho2                    &
             -28.d0/12.d0*propce(iel,ipcght)*hco

      smbrs(iel) = smbrs(iel) +                                     &
           ( gmech+gmvap+gmhet )*rom*volume(iel)
      rhovst = ( propce(iel,ipchgl)                                 &
                -propce(iel,ipcgev)*hfov )/cp2fol                   &
              *rom*volume(iel)
      rovsdt(iel) = rovsdt(iel) +  max(zero,rhovst)

    endif

  enddo

! --> T.S. pour la masse de liquide

elseif ( ivar .ge. isca(iyfol(1))     .and.                       &
         ivar .le. isca(iyfol(nclafu))        ) then

  if (iwarni(ivar).ge.1) then
    write(nfecra,1000) chaine(1:8)
  endif

  numcla = ivar-isca(iyfol(1))+1

  ipcro2 = ipproc(irom2 (numcla))
  ipcdia = ipproc(idiam2(numcla))
  ipcte2 = ipproc(itemp2(numcla))
  ipcgev = ipproc(igmeva(numcla))
  ipcght = ipproc(igmhtf(numcla))
  ipchgl = ipproc(ih1hlf(numcla))

  do iel = 1, ncel

    t2mt1 =  propce(iel,ipcte2)-propce(iel,ipcte1)
    gmvap = -propce(iel,ipcgev)*t2mt1
    gmhet = -propce(iel,ipcght)

    smbrs(iel) = smbrs(iel)                                       &
         - crom(iel)*volume(iel)*(gmvap+gmhet)
    if ( rtpa(iel,ivar).gt.epsifl ) then
      rhovst = crom(iel)*volume(iel)*(gmvap + gmhet)       &
              / rtpa(iel,ivar)
    else
      rhovst = 0.d0
    endif
    rovsdt(iel) = rovsdt(iel) + max(zero,rhovst)

  enddo

! --> T.S. pour le traceur de la vapeur

elseif ( ivar .eq. isca(ifvap) ) then

  if (iwarni(ivar).ge.1) then
    write(nfecra,1000) chaine(1:8)
  endif

  do icla = 1, nclafu

    ipcte2 = ipproc(itemp2(icla))
    ipcgev = ipproc(igmeva(icla))

    do iel = 1, ncel

      t2mt1 = propce(iel,ipcte2)-propce(iel,ipcte1)
      if ( rtpa(iel,isca(iyfol(icla))) .gt. epsifl ) then
        gmvap = -propce(iel,ipcgev)*t2mt1*rtp(iel,isca(iyfol(icla)))    &
                / rtpa(iel,isca(iyfol(icla)))
      else
        gmvap = -propce(iel,ipcgev)*t2mt1
      endif

      smbrs(iel) = smbrs(iel)                                     &
                 + gmvap*crom(iel)*volume(iel)
    enddo

  enddo

! --> T.S. pour le traceur du C ex réaction heterogene

elseif ( ivar .eq. isca(if7m) ) then

  if (iwarni(ivar).ge.1) then
    write(nfecra,1000) chaine(1:8)
  endif

  do icla = 1, nclafu

    ipcght = ipproc(igmhtf(icla))

    do iel = 1, ncel
      if (rtpa(iel,isca(iyfol(icla))) .gt. epsifl) then
        smbrs(iel) = smbrs(iel)                                        &
             -crom(iel)*propce(iel,ipcght)*volume(iel)        &
                                *rtp(iel,isca(iyfol(icla)))            &
                                /rtpa(iel,isca(iyfol(icla)))
      else
        smbrs(iel) = smbrs(iel)                                        &
                    -crom(iel)*propce(iel,ipcght)*volume(iel)
      endif

    enddo

  enddo

endif

! --> Terme source pour la variance du traceur 4 (Air)

if ( ivar.eq.isca(ifvp2m) ) then

  if (iwarni(ivar).ge.1) then
    write(nfecra,1000) chaine(1:8)
  endif

! ---- Calcul des termes sources explicite et implicite
!      relatif aux echanges interfaciaux entre phases

  call cs_fuel_fp2st &
 !==================
 ( iscal  ,                                                        &
   rtpa   , rtp    , propce ,                                      &
   smbrs  , rovsdt )

endif


! --> Terme source pour CO2

if ( ieqco2 .ge. 1 ) then

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
!
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
       tautur = cka(iel)/cvara_ep(iel)

       x2 = 0.d0
       do icla = 1, nclafu
         x2 = x2 + rtpa(iel,isca(iyfol(icla)))
       enddo

!    On transporte CO2

       smbrs(iel)  = smbrs(iel)                                   &
                    +wmole(ico2)/propce(iel,ipproc(irom1))        &
         * (xco2eq-xxco2)/(tauchi+tautur)                         &
         * (1.d0-x2)                                              &
         * volume(iel) * crom(iel)
!
       w1 = volume(iel)*crom(iel)/(tauchi+tautur)
       rovsdt(iel) = rovsdt(iel) +   max(w1,zero)

     else
       rovsdt(iel) = rovsdt(iel) + 0.d0
       smbrs(iel)  = smbrs(iel)  + 0.d0
     endif

   enddo
!
   if(irangp.ge.0) then
     call parcpt(nberic)
     call parmax(err1mx)
     call parcpt(nbpass)
     call parcpt(nbarre)
     call parcpt(nbarre)
     call parcmx(nbimax)
   endif

   write(nfecra,*) ' Max Error = ', err1mx
   write(nfecra,*) ' no Points   ', nberic, nbarre, nbpass
   write(nfecra,*) ' Iter max number ', nbimax


  endif

endif


! --> Terme source pour HCN et NO : uniquement a partir de la 2eme
!                                   iter

if ( ieqnox .eq. 1 .and. ntcabs .gt. 1) then

  if ( ivar.eq.isca(iyhcn) .or. ivar.eq.isca(iyno) ) then

    iexp1  = ipproc(ighcn1)
    iexp2  = ipproc(ighcn2)
    iexp3  = ipproc(ignoth)

! QPR= %N libéré pendant l'evaporation/taux de matieres volatiles
!          moyen

    qpr = 1.3d0

! YMOY = % vapeur en sorties

    ymoy = 0.7d0

! Azote dans le fuel

    fn = 0.015

! Masse molaire

    wmhcn = wmole(ihcn)
    wmno  = 0.030d0
    wmo2  = wmole(io2)

    if ( ivar.eq.isca(iyhcn) ) then

!        Terme source HCN

      if (iwarni(ivar).ge.1) then
        write(nfecra,1000) chaine(1:8)
      endif

      auxmin = 1.d+20
      auxmax =-1.d+20

      do iel=1,ncel
!
        xxo2 = propce(iel,ipproc(iym1(io2)))                        &
              *propce(iel,ipproc(immel))/wmo2

        aux = volume(iel)*crom(iel)                       &
             *( propce(iel,iexp2)                                  &
               +propce(iel,iexp1)*rtpa(iel,isca(iyno))             &
                                 *propce(iel,ipproc(immel))/wmno )

        smbrs(iel)  = smbrs(iel)  - aux*rtpa(iel,ivar)
        rovsdt(iel) = rovsdt(iel) + aux

        gmvap = 0.d0
        gmhet = 0.d0
        do icla=1,nclafu

          ipcgev = ipproc(igmeva(icla))
          ipcght = ipproc(igmhtf(icla))
          ipcte2 = ipproc(itemp2(icla))
          ipcte1 = ipproc(itemp1)

          gmvap = gmvap                                           &
                 + crom(iel)*propce(iel,ipcgev)          &
                  *(propce(iel,ipcte2)-propce(iel,ipcte1))

          gmhet = gmhet                                           &
                 +crom(iel)*propce(iel,ipcght)

        enddo
        if ( xxo2 .gt. 0.03d0 ) then
          aux = -volume(iel)*fn*wmhcn/(wmole(in2)/2.d0)           &
                *( qpr*gmvap+(1.d0-qpr*ymoy)/(1.d0-ymoy)*gmhet )
        else
          aux = -volume(iel)*fn*wmhcn/(wmole(in2)/2.d0)           &
                            *(qpr*gmvap)
        endif
        smbrs(iel)  = smbrs(iel) + aux

      enddo

    endif

    if ( ivar.eq.isca(iyno) ) then

!        Terme source NO

      if (iwarni(ivar).ge.1) then
        write(nfecra,1000) chaine(1:8)
      endif

      do iel=1,ncel

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

!--------
! Formats
!--------

 1000 format(' TERMES SOURCES PHYSIQUE PARTICULIERE POUR LA VARIABLE '  &
       ,a8,/)

!----
! End
!----

return

end subroutine
