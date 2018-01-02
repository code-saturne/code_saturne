!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2018 EDF S.A.
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
!> \brief Atmospheric module - humid atmosphere - cloud droplet nucleation

!> \brief Compute aerosols nucleation source term
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[out]      nc       Dropplet number (scalar)in 1/cm**3
!> \param[in]       w        Vertical wind speed in m/s
!> \param[in]       rom      density of air in kg/m**3
!> \param[in]       tempc    Temperature (scalar) in degre celsius
!> \param[in]       qldia    masse fraction of liquid water in kg/kg
!> \param[in]       pphy     true pressure in pascals
!> \param[in]       refrad   radiative cooling
!-------------------------------------------------------------------------------
subroutine nuclea (nc,w,rom,tempc,qldia,pphy,refrad)

!===============================================================================
! FONCTION :
! ----------

!   Routine qui traite le terme source
!      du a la nucleation des aerosols en gouttes d'eau liquide
!      ce terme est pris en compte explicitement pour nc
!
!   C'est la premiere routine microphysique par laquelle on
!      passe apres l'advection et la diffusion, on corrige donc nc
!      si necessaire.
!
!===============================================================================
! Module files
!===============================================================================

use paramx
use entsor
use optcal
use cstphy
use cstnum, only: pi
use ppppar
use ppthch
use ppincl
use atincl
use atsoil
use spefun

!===============================================================================

implicit none

! Arguments

double precision nc(ncelet),w(3,ncelet),qldia(ncelet)
double precision tempc(ncelet),rom(ncelet),pphy(ncelet)
double precision refrad(ncelet)

! Local variables

integer          iel, ii
double precision nuc,constc,constk,esat,tempk
double precision fbeta,aa1,aa2,aa3,ddv,kka
double precision sursat,constmu,constbeta,yy,tmpsur
double precision aa4, cp

! Declaration des fonctions

external         esatliq
double precision esatliq

double precision , parameter :: rhoeau=1000. ! kg/m**3

! ===============================================================================
! 0.  Initialisation
! ===============================================================================

cp     = cp0
constc = 0.d0
constk = 0.d0
fbeta  = 0.d0
nuc    = 0.d0

if (modnuc.eq.1) then
  !  Const. pour le modele de Pruppacher et Klett (1997)
  !  Cas du Buffalo, New-York, Kocmond [1965]
  constc = 3500d0
  constk = 0.9

  !  fbeta : fonction beta ( constk/2 , 3/2 )
  fbeta = beta(constk/2.d0, 1.5d+00)

else if (modnuc.eq.2) then
  !  Const. pour le modele de Cohard et Pinty (1998)
  !  Cas general
  constc    = 3270d0
  constk    = 1.56d0
  constmu   = 0.70d0
  constbeta = 136.d0

  !  Pollued Air Case. See Hudson and Li (1995)
  !  constc    = 1865.
  !  constk    = 0.86
  !  constmu   = 1.50
  !  constbeta = 6.80

  fbeta = beta(constk/2.d+00, 1.5d+00)
endif

! ===============================================================================
! 1.  calcul du nouveau champ Nc (en cm^-3)
! ===============================================================================

do iel = 1, ncel
  if (qldia(iel).gt.0.d0) then

    ! Calcul du nombre de gouttelettes creees par nucleation
    ! si W > 0 or Refrad < 0.
    ! si ce nombre est superieur a nc, on ajoute la difference.

    if (w(3,iel).gt.0.d0) then

      tempk = tempc(iel) + tkelvi
      aa1   = 0.622*clatev*9.81/(rair*cp*tempk**2) - 9.81/(rair*tempk)
      esat  = esatliq(tempk)
      aa2   = rair*tempk/(0.622*esat) + (0.622*clatev**2)/(tempk*pphy(iel)*cp)
      ddv   = 0.211*( tempk/tkelvi )**(1.94)*(101325d0/pphy(iel))*1.d-4
      kka   = ( (5.69 + 0.017*tempc(iel))/0.239 ) * 1.d-3
      aa3   = 1.d0/((1000.d0*rvap*tempk)/(esat*ddv)                             &
            + clatev*1000.d0*(clatev/(tempk*rvap) - 1.d0)/(kka*tempk) )
      aa4   = -0.622*clatev/(rair*tempk**2)

! -------------------------------------------------------------------------------
! 1.1  Modele Pruppacher et Klett (1997)
! -------------------------------------------------------------------------------

      if (modnuc.eq.1) then

        nuc = (constc)**(2./(constk+2.)) * ( 0.01 * (aa1*w(3,iel)                 &
             + aa4*refrad(iel))**(3./2.)                                        &
             / (2.*pi*rhoeau*aa2*aa3**(3.d0/2.d0)*constk*fbeta))                &
             **(constk/(constk+2.d0))

! -------------------------------------------------------------------------------
! 1.2  Modele Cohard et Pindy (1998)
! -------------------------------------------------------------------------------

      else if (modnuc.eq.2) then

        ! calcul de la sursaturation max : sursat
        sursat = 0.d0
        yy    = 1.d0
        do ii = 1, 20
          tmpsur = sursat
          yy     = hypgeo(constmu,constk/2.d0,constk/2.d0 + 3.d0/2.d0,          &
                 - constbeta*sursat**2 )
          sursat = ((0.01 * (aa1*w(3,iel) + aa4*refrad(iel))                      &
               **(3.d0/2.d0) / (2.d0*constk*constc*pi*rhoeau                    &
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

! -------------------------------------------------------------------------------
! 1.3  Modele Abdul-Razzak et al. (1998)
! -------------------------------------------------------------------------------

      else if (modnuc.eq.3) then
        write(nfecra,1102)
        call csexit(1)
      endif ! modnuc.eq.1 , 2 , ...

! ===============================================================================
! 2.  Calcul de la difference
! ===============================================================================

      nuc     = max(nuc - nc(iel),0.d0)
      nc(iel) = nc(iel) + nuc

! ===============================================================================
! 3.  End subroutine
! ===============================================================================
      ! si qc > 0, w <= 0 & nc = 0, on cree nc > 0 tel que le rayon
      ! volumique moyen = 10 microns

    elseif (nc(iel).eq.0) then
      nc(iel) = 1.d-6*(3.*rom(iel)*qldia(iel))/(4*pi*rhoeau*(10.d-6)**3)

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
'@ @@ ERREUR : Le modele Abdul-Razzak et al. (1998) n''est    ',/,&
'@    ======                                                  ',/,&
'@             pas encore implemente.                         ',/,&
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
'@ @@ ERROR: The Abdul-Razzak et al. model (1998) is not      ',/,&
'@    =====                                                   ',/,&
'@           implemented yet.                                 ',/,&
'@                                                            '  )

#endif

! ----
!  End
! ----

return
end subroutine nuclea
