!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2020 EDF S.A.
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
!> \file nuclea.f90
!> \brief Compute aerosol cloud droplets nucleation when using the atmospheric
!> humid model using a microphysical model.
!>
!> It is taken into account as an additional step split from advection-diffusion
!> equation, hence the droplet number is first clipped if necessary.
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[out]    nc            Droplet number (scalar)in 1/cm**3
!> \param[in]     w             Vertical wind speed in m/s
!> \param[in]     rom           density of air in kg/m**3
!> \param[in]     tempc         Temperature (scalar) in degre celsius
!> \param[in]     qldia         masse fraction of liquid water in kg/kg
!> \param[in]     pphy          true pressure in pascals
!> \param[in]     refrad        radiative cooling
!-------------------------------------------------------------------------------

subroutine nuclea (nc, w, rom, tempc, qldia, pphy, refrad)

!===============================================================================
! Module files
!===============================================================================

use paramx
use entsor
use optcal
use cstphy
use cstnum, only: pi, epzero
use ppppar
use ppthch
use ppincl
use atincl
use atsoil
use spefun
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

double precision nc(ncelet),w(3,ncelet),qldia(ncelet)
double precision tempc(ncelet),rom(ncelet),pphy(ncelet)
double precision refrad(ncel)

! Local variables

integer          iel, ii
double precision nuc,constc,constk,esat,tempk
double precision fbeta,aa1,aa2,aa3,ddv,kka
double precision sursat,constmu,constbeta,yy,tmpsur
double precision aa4, cp

double precision nbion(4),frasolu(4),mmaero(4)
double precision rhoaero(4),fraero(4),coefosm(4)
double precision rayonm1,ntotal1,sigaero1,mmh2o,smax
double precision rayonm2,ntotal2,sigaero2
double precision rayonm3,ntotal3,sigaero3
double precision tauvw
double precision coefa,coefb,coefavg,coefzeta
double precision foncf1,foncg1,coefeta1
double precision foncf2,foncg2,coefeta2
double precision foncf3,foncg3,coefeta3
double precision ugauss1,sursat1,smax1,nuc1
double precision ugauss2,sursat2,smax2,nuc2
double precision ugauss3,sursat3,smax3,nuc3
double precision nuc1snc,nuc2snc,nuc3snc
double precision numcb,dencb

double precision , parameter :: rhowater=1000. ! kg/m**3
! FIXME should be set somewhere else

! ===============================================================================
! 0.  Initialisation
! ===============================================================================

cp     = cp0
constc = 0.d0
constk = 0.d0
fbeta  = 0.d0
nuc    = 0.d0

if (modnuc.eq.1) then
  !  Constants for the model of Pruppacher and Klett (1997)
  !  Case of Buffalo, New-York, Kocmond [1965]
  constc = 3500d0
  constk = 0.9d0

  !  fbeta : beta function ( constk/2 , 3/2 )
  fbeta = beta(constk/2.d0, 1.5d+00)

else if (modnuc.eq.2) then
  !  Constants for model of Cohard and Pinty (1998)
  !  (general case)
  constc    = 3270d0
  constk    = 1.56d0
  constmu   = 0.70d0
  constbeta = 136.d0

  !  Polluted Air Case. See Hudson and Li (1995)
  !  constc    = 1865. ! 8700.d0 if Cohard and Pinty used on PARISFOG
  !  constk    = 0.86
  !  constmu   = 1.50
  !  constbeta = 6.80

  fbeta = beta(constk/2.d+00, 1.5d+00)
endif

! ===============================================================================
! 1.  Computation new nc field
! ===============================================================================

do iel = 1, ncel

  if (qldia(iel).le.epzero) then

     nc(iel) = 0.d0

  else ! qldia(iel) > epzero

    ! New droplets are created by nucleation if w > 0 or refrad < 0,
    ! then if number of new droplets is superior to nc, the difference is added

    if (w(3,iel).gt.epzero.or.refrad(iel).lt.epzero) then

      tempk = tempc(iel) + tkelvi
      aa1   = 0.622d0*clatev*9.81d0/(rair*cp*tempk**2) - 9.81d0/(rair*tempk)
      esat  = cs_air_pwv_sat(tempc(iel))
      aa2   = rair*tempk/(0.622d0 *esat) + (0.622d0 *clatev**2)/(tempk*pphy(iel)*cp)
      ddv   = 0.211d0 * ( tempk/tkelvi )**(1.94d0)*(101325d0/pphy(iel))*1.d-4
      kka   = ( (5.69d0 + 0.017d0*tempc(iel))/0.239d0 ) * 1.d-3
      aa3   = 1.d0/((1000.d0*rvap*tempk)/(esat*ddv)                             &
            + clatev*1000.d0*(clatev/(tempk*rvap) - 1.d0)/(kka*tempk) )
      aa4   = -0.622d0*clatev/(rair*tempk**2)

      if ((aa1*w(3,iel)+aa4*refrad(iel)).gt.epzero) then ! alphaV/G > 0

       ! 1.1  Model of Pruppacher and Klett (1997)

       if (modnuc.eq.1) then

         nuc = (constc)**(2./(constk+2.)) * ( 0.01d0 * (aa1*w(3,iel)                 &
              + aa4*refrad(iel))**(3./2.)                                        &
              / (2.*pi*rhowater*aa2*aa3**(3.d0/2.d0)*constk*fbeta))                &
              **(constk/(constk+2.d0))

       ! 1.2  Modele of Cohard and Pinty (1998)

       else if (modnuc.eq.2) then

         ! compute max supersaturation: sursat
         sursat = 0.d0
         yy    = 1.d0
         do ii = 1, 20
           tmpsur = sursat
           yy     = hypgeo(constmu,constk/2.d0,constk/2.d0 + 3.d0/2.d0,          &
                - constbeta*sursat**2 )
           sursat = ((0.01 * (aa1*w(3,iel) + aa4*refrad(iel))                      &
                **(3.d0/2.d0) / (2.d0*constk*constc*pi*rhowater                    &
                *aa2*fbeta*aa3**(3.d0/2.d0)) )/yy )                              &
                **(1.d0/(constk + 2.d0))
         enddo
         if (abs(tmpsur - sursat).gt.1.d-2) then
           write(nfecra,1100) abs(tmpsur - sursat)
         endif

         nuc = constc * (sursat**constk) *                                       &
              hypgeo(constmu,constk/2.d0,constk/2.d0 + 1.d0,                    &
              - constbeta*sursat**2)

         if (nuc.lt.0.d0) then
           write(nfecra,1101)
           call csexit(1)
         endif

       ! 1.3  Model of Abdul-Razzak and al. (1998)
       ! This model requires a fine knowledge of aerosols chemistry.

       else if (modnuc.eq.3) then

         ! Warning: this set of constants fits for PARISFOG case.
         ! They were determined using fine data on aerosols measured
         ! during PARISFOG POI 13
         ntotal1  = 8700.d0
         rayonm1  = 0.0165d-6
         sigaero1 = 1.57d0
         ntotal2  = 8300.d0
         rayonm2  = 0.055d-6
         sigaero2 = 1.59d0
         ntotal3  = 1000.d0
         rayonm3  = 0.40d-6
         sigaero3 = 1.3d0

         ! surface tension water-steam [n/m], tempk [k]
         tempk = tempc(iel)+tkelvi
         tauvw = 0.075d0

         ! fractions in %
         ! 1: sulfate
         ! 2: nitrate
         ! 3: black carbon
         ! 4: organic carbon
         fraero(1) = 0.25d0
         fraero(2) = 0.20d0
         fraero(3) = 0.164d0
         fraero(4) = 0.386d0

         ! number of ions produced by dissociation of a salt molecule in water
         ! (if NaCl, nbion=2)
         nbion(1)= 3.d0
         nbion(2)= 2.d0
         nbion(3)= 0.d0
         nbion(4)= 1.d0

         ! osmotic coefficient
         coefosm(1) = 1.d0
         coefosm(2) = 1.d0
         coefosm(3) = 1.d0
         coefosm(4) = 1.d0

         ! massic fraction of of soluble substance in aerosol mix: msolu/mmix
         frasolu(1) = 1.d0
         frasolu(2) = 1.d0
         frasolu(3) = 0.d0
         frasolu(4) = 0.1d0

         ! molecular mass of water
         mmh2o = 18.d-3

         ! molecular mass of aerosol components
         mmaero(1) = 132.d-3
         mmaero(2) = 80.d-3
         mmaero(3) = 250.d-3
         mmaero(4) = 250.d-3

         ! density of aerosol !FIXME to check
         rhoaero(1) = 1.77d3
         rhoaero(2) = 1.77d3
         rhoaero(3) = 1.77d3
         rhoaero(4) = 1.77d3

         ! coefficient for Kelvin effect: coefa [/m]
         coefa = 2.d0*tauvw/(rhowater*rvap*tempk)

         ! FIXME
         ! cf pruppacher and klett 1997 (reprinted correction 2000) (6-28)
         ! coefa = 3.3d-7 / tempk

         ! coefficient for Raoult effect: coefb [-]
         numcb = 0.d0
         dencb = 0.d0
         do ii = 1, 4
           numcb = numcb + fraero(ii)*nbion(ii)*coefosm(ii)*frasolu(ii)   &
                   / mmaero(ii)
           dencb = dencb + fraero(ii)/rhoaero(ii)
         enddo
         coefb = mmh2o*numcb/(dencb*rhowater)

         ! supersaturation [-]
         sursat1 = (2./sqrt(coefb))*((coefa/3./rayonm1)**(1.5))
         sursat2 = (2./sqrt(coefb))*((coefa/3./rayonm2)**(1.5))
         sursat3 = (2./sqrt(coefb))*((coefa/3./rayonm3)**(1.5))

         ! standard deviation function
         foncf1 = 0.5d0 * exp(2.5d0*(log(sigaero1))**2.)
         foncf2 = 0.5d0 * exp(2.5d0*(log(sigaero2))**2.)
         foncf3 = 0.5d0 * exp(2.5d0*(log(sigaero3))**2.)

         foncg1 = 1.d0 + 0.25d0 * (log(sigaero1))
         foncg2 = 1.d0 + 0.25d0 * (log(sigaero2))
         foncg3 = 1.d0 + 0.25d0 * (log(sigaero3))

         ! coefficent corresponding to vertical velocity |aa1|
         aa1 = 0.622d0*clatev*9.81d0/(rair*cp*tempk**2.d0)-          &
               9.81d0/(rair*tempk)

         ! coefficient corresponding to liquid water condensation |aa2|
         esat = cs_air_pwv_sat(tempc(iel))
         aa2  = rair*tempk/(0.622d0*esat)+(0.622d0*clatev**2.d0) /   &
                (tempk*pphy(iel)*cp)

         ! coefficient corresponding to the droplet growth |aa3|
         ddv = 0.211d0 * (tempk/tkelvi)**(1.94d0)*                   &
               (101325.d0/pphy(iel))*1.d-4
         rvap = rvsra*rair
         kka  = ((5.69d0 + 0.017d0 * tempc(iel)) / 0.239d0) * 1.d-3
         aa3  = 1.d0 / ( (1000.d0*rvap*tempk)/(esat*ddv) + clatev *  &
                  1000.d0*(clatev/(tempk*rvap)-1.d0)/(kka*tempk) )

         ! coefficient corresponding to infrared radiation |aa4|
         aa4 = -0.622d0*clatev/(rair*tempk**2)

         ! alphaV/G
         coefavg = (aa1*w(3,iel)+aa4*refrad(iel))/aa3

         ! coefficent zeta
         coefzeta = 2.*coefa/3.*sqrt(coefavg)

         ! coefficient eta - Ntot [m^(-3)]
         coefeta1 = (coefavg**1.5d0)/(2.d0*pi*rhowater*aa2*ntotal1*1.d6)
         coefeta2 = (coefavg**1.5d0)/(2.d0*pi*rhowater*aa2*ntotal2*1.d6)
         coefeta3 = (coefavg**1.5d0)/(2.d0*pi*rhowater*aa2*ntotal3*1.d6)

         ! max supersaturation (smax)
         smax1 = (foncf1*((coefzeta/coefeta1)**1.5d0) +              &
                 foncg1*((sursat1*sursat1)/(coefeta1+3.d0*coefzeta)) &
                 **0.75d0)/sursat1/sursat1
         smax2 = (foncf2*((coefzeta/coefeta2)**1.5d0) +              &
                 foncg2*((sursat2*sursat2)/(coefeta2+3.d0*coefzeta)) &
                 **0.75d0)/sursat2/sursat2
         smax3 = (foncf3*((coefzeta/coefeta3)**1.5d0) +              &
                 foncg3*((sursat3*sursat3)/(coefeta3+3.d0*coefzeta)) &
                 **0.75d0)/sursat3/sursat3
         smax = 1.d0/sqrt(smax3+smax2+smax1)

         if (    smax.gt.1.d0.or.sursat1.gt.1.d0 &
                             .or.sursat2.gt.1.d0 &
                             .or.sursat3.gt.1.d0) then
           write(nfecra,1102)
           call csexit(1)
         endif

         ugauss1 = (2.d0*log(sursat1/smax))/(3.d0*sqrt(2.d0)*log(sigaero1))
         ugauss2 = (2.d0*log(sursat2/smax))/(3.d0*sqrt(2.d0)*log(sigaero2))
         ugauss3 = (2.d0*log(sursat3/smax))/(3.d0*sqrt(2.d0)*log(sigaero3))

         nuc1   = 0.5d0 * ntotal1 * (1.d0-ferf(ugauss1))
         nuc2   = 0.5d0 * ntotal2 * (1.d0-ferf(ugauss2))
         nuc3   = 0.5d0 * ntotal3 * (1.d0-ferf(ugauss3))
         ! nuc=nuc1+nuc2
         nuc=nuc1+nuc2+nuc3
         ! nuc=nuc3

         nuc1snc = nuc1/ntotal1
         nuc2snc = nuc2/ntotal2
         nuc3snc = nuc3/ntotal3
       endif ! modnuc

       ! 2.  Compute difference

       nuc = max(nuc - nc(iel),0.d0)

     else
       nuc = 0.d0
     endif

     nc(iel) = nc(iel) + nuc

     ! 3. if qc > 0, w <= 0 and nc = 0,
     ! we impose nc > 0 so that the mean volume radius = 10 microns

   elseif (nc(iel).lt.epzero) then

     nc(iel) = 1.d-6*(3.d0*rom(iel)*qldia(iel))/(4.d0*pi*rhowater*(10.d-6)**3)

   endif ! end w > 0

  endif ! end qldia > 0
enddo

!--------
! Formats
!--------

#if defined(_CS_LANG_FR)

 1100 format (                                                          &
'@                                                            ',/,&
'@ @@ ATTENTION : la surtaturation maximale n''a pas convergee',/,&
'@    =========                                               ',/,&
'@  Residu = ',E12.5                                           ,/,&
'@                                                            '  )

 1101 format (                                                          &
'@                                                            ',/,&
'@ @@ ERREUR : modele Cohard et Pindy (1998).                 ',/,&
'@    ======                                                  ',/,&
'@  La nucleation  est negative.                              ',/,&
'@                                                            '  )

 1102 format (                                                          &
'@                                                            ',/,&
'@ @@ ERREUR : modele Abdul-Razzak et al. (1998).             ',/,&
'@    ======                                                  ',/,&
'@  Sursaturation negative.                                   ',/,&
'@                                                            '  )

#else

 1100 format (                                                          &
'@                                                            ',/,&
'@ @@ WARNING: ',A8 ,' Maximum saturation has not converged   ',/,&
'@    ========                                                ',/,&
'@  Residue = ',E12.5                                          ,/,&
'@                                                            '  )

 1101 format (                                                          &
'@                                                            ',/,&
'@ @@ ERROR: Cohard and Pindy model (1998).                   ',/,&
'@    =====                                                   ',/,&
'@  The nucleation is negative.                               ',/,&
'@                                                            '  )

 1102 format (                                                          &
'@                                                            ',/,&
'@ @@ ERROR : Abdul-Razzak and al. model (1998).              ',/,&
'@    ======                                                  ',/,&
'@  Negative sursaturation.                                   ',/,&
'@                                                            '  )

#endif

! ----
!  End
! ----

return
end subroutine nuclea
