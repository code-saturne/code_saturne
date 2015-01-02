!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2015 EDF S.A.
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

!                   LIBRAIRIE DE SOUS-PROGRAMMES RMODAK
!==========================================


!===============================================================================
! FONCTION :
! --------

!   SOUS-PROGRAMME DU MODULE DE RAYONNEMENT :
!   -----------------------------------------

!   SOUS-PROGRAMMES DE LA PHYSIQUE PARTICULIERE
!    RELATIFS AU CALCUL DU COEFFICIENT D'ABSORPTION AVEC MODAK

! References : MODAK A.T.,
!             "Radiation from products of combustion",
!             Fire Research, 1 pp. 339-361, 1978.

!             MECHITOUA N.
!             "Modelisation numerique du rayonnment dans les
!              milieux semi-transoparents",
!             Raport EDF, HE/44/87-15, 1987.


!==============================================================================

subroutine absorb &
!================

 ( ts    , te     , path   , sootk  ,                             &
   pco2  , ph2o   , alpha  )



!===============================================================================
! FONCTION :
! --------

!    ON CALCULE LES ABSORPTIVITES (PAR RAPPORT
!    A UNE SOURCE DE CORPS NOIR) D'UN MELANGE GAZEUX ISOTHERME,
!    HOMOGENE DE SUIE, CO2 ET H2O A LA PRESSION TOTALE D'1 ATM.

!    SI LA TEMPERATURE DU CORPS NOIR EST EGALE A LA TEMPERATURE DU
!    MELANGE, L'ABSORPTIVITE EST EGALE A L'EMISSIVITE.
!    LES EMISSIVITES AINSI CALCULEES SONT EN BON ACCORD AVEC LES
!    CALCULS SPECTRAUX ET LES MESURES EXPERIMENTALES

!    TS ET TE DOIVENT ETRE COMPRIS ENTRE 300 ET 2000 KELVIN

!    LA LONGUEUR D'ONDE 0.94 MICRON.
!    SOOTK EST LIEE A LA FRACTION VOLUMIQUE DE SUIE FV SUIVANT
!    LA FORMULE :
!                   SOOTK=7FV/0.94E-6

! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! ts               ! r  ! <-- ! temperature du corps noir (k)                  !
! te               ! r  ! <-- ! temperature du melange (k)                     !
! path             ! r  ! <-- ! penetration du rayonnement dans le             !
!                  !    !     ! melange (m)                                    !
! sootk            ! r  ! <-- ! coefficient d'absorption des suies             !
! pco2             ! r  ! <-- ! pression partielle de co2 dans un              !
!                  !    !     !  melange de presssion totale 1 atm.            !
! ph2o             ! r  ! <-- ! pression partielle de h2o dans un              !
!                  !    !     !  melange de presssion totale 1 atm.            !
! alpha            ! r  ! <-            ! absorptivite
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
use entsor

!===============================================================================

implicit none

! Arguments

double precision ts, te, path, sootk, pco2, ph2o, alpha

! Local variables

double precision tmax, tmin, ptotal, ratio, pathl, pcl, pwl
double precision as, taus, ag, power, zeta
double precision emigas

!===============================================================================
! CALCUL
!===============================================================================

tmax = 3000.d0
tmin = 298.d0
if (ts.lt.tmin  .or.  ts.gt.tmax ) goto 1
if (te.lt.tmin  .or.  te.gt.tmax ) goto 2

! --- Pression totale : PTOTAL

ptotal = pco2 + ph2o

if ( ptotal.gt.1.d0 ) goto 3

! --- Rapport temeperature melange et temperature source : RATIO

ratio = te/ts

! --- Loncueur de penetration du rayonment effectif : PATHL

pathl = path/ratio
pcl   = pco2*pathl
pwl   = ph2o*pathl
if (pcl.gt.5.98d0 .or. pwl.gt.5.98d0) goto 4

! --- Calcul de l'absortivite des suies : AS

as = 0.d0
if (sootk.le.0.d0) goto 51
call tasoot                                                       &
!==========
 ( sootk , path , ts , taus )

as = 1.d0-taus

 51   continue

! --- Calcul de l'absorptivite du gaz : AG
!                 = emissivite du gaz

ag = 0.d0
if (pco2.lt.0.0011d0  .and.   ph2o.lt.0.0011d0) goto 52
if (pcl .lt.0.0011d0  .and.   pwl .lt.0.0011d0) goto 52
ag = emigas( pathl, pco2, ph2o, ts )

! --- Calcul de la fraction de vapeur d'eau : ZETA

zeta  = ph2o/ptotal
power = 0.65d0-0.2d0*zeta
ag    = ag*(ratio**power)

 52   continue
alpha = as + ag -as*ag
!     ALPHA = ABS(ALPHA)
if (alpha.le.1.d-8) goto 8

return

 4    continue
write(nfecra,1000)
goto 8

 3    continue
write(nfecra,1001)
goto 8

 2    continue
write(nfecra,1002)
goto 8

 1    continue
write(nfecra,1003)

 8    continue
alpha= 1.d-8


!========
! FORMATS
!========

 1000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ERREUR MODAK :                                          ',/,&
'@    ============                                            ',/,&
'@ LE PRODUIT PATH*TS/T*PCO2 OU PATH*TS/T*PH2O                ',/,&
'@ DEPASSE LA VALEUR 5.98 ATM.METRE.                          ',/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 1001 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ERREUR MODAK :                                          ',/,&
'@    ============                                            ',/,&
'@ LA SOMME DES PRESSIONS PARTIELLES DES GAZ CO2 ET H2O       ',/,&
'@ DEPASSE UN ATMOSPHERE.                                     ',/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 1002 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ERREUR MODAK :                                          ',/,&
'@    ============                                            ',/,&
'@ LA TEMPERATURE DU MELANGE TE                               ',/,&
'@ SORT DES LIMITES DU DOMAINE.                               ',/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 1003 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ERREUR MODAK :                                          ',/,&
'@    ============                                            ',/,&
'@ LA TEMPERATURE DU CORPS NOIR TS                            ',/,&
'@ SORT DES LIMITES DU DOMAINE.                               ',/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

return
end subroutine

!==============================================================================

subroutine chebyc &
!================

 ( norpol , argpol , valpol )


!===============================================================================
! FONCTION :
! --------

!     CALCUL DU POLYNOME DE CHEBYCHEV D'ORDRE NORPOL

! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! norpol           ! e  ! <-- ! ordre du polynome de  chebychev                !
! argpol           ! r  ! <-- ! argument du polynome de chebychev              !
! valpol           ! r  ! --> ! valeur du polynome de chebychev                !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

!===============================================================================

implicit none

! Arguments

integer          norpol
double precision argpol , valpol

! Local variables

integer          ict
double precision f, vm2,vm1

!===============================================================================
! CALCUL
!===============================================================================

valpol = 1.d0
if (norpol.le.0) then
  goto 1
else
  goto 2
endif

 1    return

 2    valpol = argpol

if ( (norpol-1).le.0) then
  goto 1
else
  goto 3
endif
 3    f   = argpol+argpol
vm1 = argpol
vm2 = 1.d0
do ict = 2, norpol
  valpol = f*vm1-vm2
  vm2    = vm1
  vm1    = valpol
enddo

return
end subroutine

!==============================================================================

subroutine asympt &
!================

 ( zz    , zzv    )


!===============================================================================
! FONCTION :
! --------

!FONCC   CALCUL DE L'EXPANSION ASYMPTOTIQUE POUR LA FONCTION PENTAGAMMA
!FONCC
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! zz               ! r  ! <-- !                                                !
! zzv              ! r  ! --> !                                                !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

!===============================================================================

implicit none

! Arguments

double precision zz, zzv

! Local variables

double precision zi1, zi2, zi3, d1s3

!===============================================================================
! CALCUL
!===============================================================================

d1s3 = 1.d0/3.d0
zi1 = 1.d0/zz
zi2 = zi1*zi1
zi3 = zi1*zi2
zzv = zi3 * (                                                     &
              (2.d0+3.d0*zi1) +                                   &
              zi2*(2.d0 + zi2*(-1.d0 +                            &
                               zi2*( 1.d0+d1s3                    &
                                     + zi2*(-3.d0+10.d0*zi2) )    &
                              )                                   &
                  )                                               &
            )

return
end subroutine

!==============================================================================

subroutine tasoot &
!================

 ( zkled  , pathl  , tblack , taus   )


!===============================================================================
! FONCTION :
! --------

!  CALCUL DE LA TRANSMISSIVITE(TAUS) DE PATH
!  A UNE TEMPERATURE DONNEE.

! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! zkled            ! r  ! <-- !                                                !
! pathl            ! r  ! <-- ! penetration du rayonnement dans le             !
!                  !    !     ! melange                                        !
! tblack           ! r  ! <-- ! temperature source ou temperature              !
!                  !    !     !  du gaz                                        !
! taus             ! r  ! --> ! transmissivite de path                         !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

!===============================================================================

implicit none

! Arguments

double precision zkled, pathl, tblack, taus

! Local variables

double precision arg , val

!===============================================================================
! CALCUL
!===============================================================================

if ( zkled.le.0.d0 ) goto 1

arg = 1.d0 + zkled*pathl*tblack*6.5333d-5

call pentag                                                       &
!==========
 ( arg    , val )

taus = val*.1539897336d0
return

 1    taus = 1.d0

return
end subroutine

!==============================================================================

subroutine pentag &
!================

 ( argfpe , valfpe )


!===============================================================================
! FONCTION :
! --------

!   CALCUL DE LA VALEUR VAL DE LA
!   FONCTION PENTAGAMMA D ARGUMENT X.
!   ON UTILISE LES FORMULES ASYMPTOTIQUES ET DE RECURRENCE
!   D'ABRAMOWITZ ET STEGUN.

! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! argfpe           ! r  ! <-- ! argument de la fonction pentagamma             !
! valfpe           ! r  ! --> ! valeur de la fonction pentagamma               !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

!===============================================================================

implicit none

! Arguments

double precision argfpe, valfpe

! Local variables

double precision zz, zzv, zs

!===============================================================================
!   CALCUL
!===============================================================================

if (argfpe.ge.4.d0) goto 1
if (argfpe.ge.3.d0) goto 2
if (argfpe.ge.2.d0) goto 3

zs = ( 1.d0/(argfpe+2.d0)**4 + 1.d0/(argfpe+1.d0)**4              &
                             + 1.d0/argfpe**4          )*6.d0
zz = argfpe+3.d0
call asympt                                                       &
!==========
 ( zz, zzv )

goto 4

 3    continue
zs = (1.d0/(argfpe+1.d0)**4+1.d0/argfpe**4)*6.d0
zz = argfpe+2.d0
call asympt                                                       &
!==========
 ( zz, zzv )
goto 4

 2    continue
zs= 6.d0/argfpe**4
zz= argfpe+1.d0
call asympt                                                       &
!==========
 ( zz, zzv )

goto 4

 1    continue
zs = 0.d0
! La ligne suivante est ajoutee simplement pour eviter un warning
!    sous Foresys (var non initialisee) en attendant que la routine
!    soit revue de maniere globale (pour l'instant, l'utilisation
!    de rmodak est bloquee en amont).
zz= argfpe

call asympt                                                       &
!==========
 ( zz, zzv )

 4    continue
valfpe = zzv+zs

return
end subroutine

!==============================================================================

function fdleck &
!==============

 ( val    , pl     , te     )


!===============================================================================
! FONCTION :
! --------

!  CALCUL DE LA CORRECTION APPORTEE POUR LE MELANGE
!  DE CO2 ET H2O LORSQUE LES LONGUEURS D ONDES SONT AU DELA DE 2.7
!  ET 15 MICRONS

! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! val              ! r  ! ->  !                                                !
! pl               ! r  ! ->  !                                                !
! te               ! r  ! ->  !                                                !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

!===============================================================================

implicit none

! Arguments

double precision val, pl, te, fdleck

! Local variables

double precision term, term2, term3, tt, tt2, aa, bb, cc

!===============================================================================
!   CALCUL
!===============================================================================

if (pl.lt.0.1d0) goto 1

term  = val/(10.7d0+101.d0*val) - val**10.4d0/111.7d0
term2 = log10(101.325d0*pl)
term2 = term2**2.76d0
tt    = te/1000.d0
tt2   = tt*tt
aa    = -1.0204082d0
bb    =  2.2448979d0
cc    = -0.23469386d0
term3 = aa*tt2+bb*tt+cc

! --- TERM3 represente l'ajustement de temperature

fdleck = term*term2*term3

return

 1    fdleck= 0.d0

return
end function

!==============================================================================

function emigas &
!==============
!      -------------------------------------------------------------
 ( pathl  , pc     , pw     , te     )
!      -------------------------------------------------------------

!===============================================================================
! FONCTION :
! --------

!    CALCUL DE L EMISSIVITE A UN PATH DONNE
!    D'UN MELANGE DE CO2 ET H2O A LA TEMPERATURE TE

! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! pathl            ! r  ! ->  ! valeur du path                                 !
! pc               ! r  ! ->  ! pression partielle de co2                      !
! pw               ! r  ! ->  ! pression partielle de h2o                      !
! te               ! r  ! ->  ! temperature                                    !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

!===============================================================================

implicit none

! Arguments

double precision pathl, pc, pw, te, emigas

! Local variables

double precision tmin, tmax, pcl, ec
double precision pwl, pcwl, dels, ew, pcpw, xi, fdleck

!===============================================================================
!   CALCUL
!===============================================================================

tmin   = 298.d0
tmax   = 3000.d0
emigas = 0.d0

if ( te.lt.tmin .or. te.gt.tmax ) return

ec = 0.d0

if ( pc.lt.0.0011d0 .or. pc.gt.1.d0 ) goto 1

pcl = pc*pathl

if ( pcl.lt.0.0011d0 .or. pcl.gt.5.98d0 ) goto 1

call scrtch                                                       &
!==========
 ( pc , pcl , te , 1 , ec)

 1    continue

if ( pw.lt.0.0011d0 .or. pw.gt.1.d0 ) goto 2

pwl = pw*pathl

if ( pwl.lt.0.0011d0 .or. pwl.gt.5.98d0 ) goto 2

call scrtch                                                       &
!==========
 ( pw , pwl , te , 2 , ew )

emigas = ec + ew

if ( ec.le.0.0d0 ) return

pcpw = pc + pw
xi   = pw / pcpw
if ( xi.lt.0.01d0 ) return
pcwl = pcpw*pathl

if ( pcwl.lt.0.1d0 ) return

dels   = fdleck(xi,pcwl,te)
emigas = emigas - dels

return

2     continue
emigas = ec

return
end function

!==============================================================================

subroutine scrtch &
!================

 ( pp     , pl     , te     , index  , val    )


!===============================================================================
! FONCTION :
! --------



! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! pp               ! r  ! ->  !                                                !
! pl               ! r  ! ->  !                                                !
! te               ! r  ! ->  !                                                !
! index            ! e  ! ->  !                                                !
! val              ! r  ! ->  !                                                !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

!===============================================================================

implicit none

! Arguments

integer          index
double precision pp, pl, te, val

! Local variables

integer          ii, jj, kk, iii, jjj, kkk
double precision cc(3,4,4), cw(3,4,4), sc(3,4,4)
double precision xx, yy, zz, v6, v7
double precision tix, tjy, tkz

!===============================================================================
!   CALCUL
!===============================================================================

if ( index.eq.2 ) goto 2

! --- CC represente untebleau de 48 elements pour CO2

cc(1,1,1) = -.2754568d1
cc(1,1,2) = -.2997857d0
cc(1,1,3) = -.1232494d0
cc(1,1,4) =  .1279287d-1
cc(1,2,1) =  .1503051d1
cc(1,2,2) =  .3156449d0
cc(1,2,3) =  .1058126d-1
cc(1,2,4) = -.3729625d-1
cc(1,3,1) = -.247411d0
cc(1,3,2) = -.3323846d-1
cc(1,3,3) = -.1819471d-1
cc(1,3,4) =  .2289789d-1
cc(1,4,1) =  .4994029d-1
cc(1,4,2) = -.1986786d-2
cc(1,4,3) =  .3007898d-2
cc(1,4,4) = -.1175598d-2
cc(2,1,1) =  .5737722d-2
cc(2,1,2) = -.9328458d-2
cc(2,1,3) =  .2906286d-2
cc(2,1,4) =  .422752d-3
cc(2,2,1) = -.3151784d-2
cc(2,2,2) =  .5632821d-2
cc(2,2,3) = -.3260295d-2
cc(2,2,4) =  .7065884d-3
cc(2,3,1) =  .1668751d-3
cc(2,3,2) = -.7326533d-3
cc(2,3,3) =  .3639855d-3
cc(2,3,4) =  .3228318d-3
cc(2,4,1) =  .7386638d-3
cc(2,4,2) = -.7277073d-3
cc(2,4,3) =  .5925968d-3
cc(2,4,4) = -.2021413d-3
cc(3,1,1) =  .3385611d-2
cc(3,1,2) = -.5439185d-2
cc(3,1,3) =  .176456d-2
cc(3,1,4) =  .3036031d-3
cc(3,2,1) = -.18627d-2
cc(3,2,2) =  .3236275d-2
cc(3,2,3) = -.195225d-2
cc(3,2,4) =  .3474022d-3
cc(3,3,1) =  .1204807d-3
cc(3,3,2) = -.4479927d-3
cc(3,3,3) =  .2497521d-3
cc(3,3,4) =  .1812996d-3
cc(3,4,1) =  .4218169d-3
cc(3,4,2) = -.4046608d-3
cc(3,4,3) =  .3256861d-3
cc(3,4,4) = -.9514981d-4

goto 4

 2    continue

! --- CW represente untebleau de 48 elements pour H2O

cw(1,1,1) = -.2594279d1
cw(1,1,2) = -.7118472d0
cw(1,1,3) = -.9956839d-3
cw(1,1,4) =  .1226560d-1
cw(1,2,1) =  .2510331d1
cw(1,2,2) =  .6481808d0
cw(1,2,3) = -.3330587d-1
cw(1,2,4) = -.5524345d-2
cw(1,3,1) = -.4191636d0
cw(1,3,2) = -.1375180d0
cw(1,3,3) =  .3877930d-1
cw(1,3,4) =  .8862328d-3
cw(1,4,1) = -.322912d-1
cw(1,4,2) = -.1820241d-1
cw(1,4,3) = -.2223133d-1
cw(1,4,4) = -.5940781d-3
cw(2,1,1) =  .1126869d0
cw(2,1,2) = -.8133829d-1
cw(2,1,3) =  .1514940d-1
cw(2,1,4) =  .1393980d-2
cw(2,2,1) = -.9298805d-2
cw(2,2,2) =  .4550660d-1
cw(2,2,3) = -.2082008d-1
cw(2,2,4) =  .2013361d-2
cw(2,3,1) = -.4375032d-1
cw(2,3,2) =  .1924597d-1
cw(2,3,3) =  .8859877d-2
cw(2,3,4) = -.4618414d-2
cw(2,4,1) =  .7077876d-2
cw(2,4,2) = -.2096188d-1
cw(2,4,3) =  .1458262d-2
cw(2,4,4) =  .3851421d-2
cw(3,1,1) =  .5341517d-1
cw(3,1,2) = -.3407693d-1
cw(3,1,3) =  .4354611d-2
cw(3,1,4) =  .1492038d-2
cw(3,2,1) = -.4708178d-2
cw(3,2,2) =  .2086896d-1
cw(3,2,3) = -.9477533d-2
cw(3,2,4) =  .6153272d-3
cw(3,3,1) = -.2104622d-1
cw(3,3,2) =  .7515796d-2
cw(3,3,3) =  .5965509d-2
cw(3,3,4) = -.2756144d-2
cw(3,4,1) =  .4318975d-2
cw(3,4,2) = -.1005744d-1
cw(3,4,3) =  .4091084d-3
cw(3,4,4) =  .2550435d-2

 4    continue

xx = log(pp)/3.45d0 + 1.d0
yy  = (log(pl)+2.555d0) / 4.345d0
zz  = (te-1150.d0) / 850.d0
val = 0.d0

do ii =1, 3
  iii = ii-1
  call chebyc                                                     &
  !==========
  ( iii , xx , tix )

  v6 = 0.d0
  do jj = 1, 4
    jjj = jj-1
    call chebyc                                                   &
    !==========
    ( jjj, yy , tjy )

    v7 = 0.d0
    do kk = 1, 4
      kkk = kk-1
      call chebyc                                                 &
      !==========
      ( kkk , zz , tkz )

      if (index.eq.1) sc(ii,jj,kk) = cc(ii,jj,kk)
      if (index.eq.2) sc(ii,jj,kk) = cw(ii,jj,kk)
      v7 = v7 + tkz*sc(ii,jj,kk)
    enddo

      v6 = v6 + v7*tjy
  enddo

  val = val + v6*tix

enddo

val = exp(val)

return
end subroutine
