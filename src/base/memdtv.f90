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

subroutine memdtv &
!================

 ( idbia0 , idbra0 ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   iviscf , iviscb , idam   , icofbd , iw1   , iw2    , iw3     , &
   icofbr , igrarx , igrary , igrarz , iwcf  ,                    &
   iptlro , ippmcf ,                                              &
   ifinia , ifinra )

!===============================================================================
!  FONCTION
!  --------

!  GESTION MEMOIRE DTTVAR

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
! ndim             ! i  ! <-- ! spatial dimension                              !
! ncelet           ! i  ! <-- ! number of extended (real + ghost) cells        !
! ncel             ! i  ! <-- ! number of cells                                !
! nfac             ! i  ! <-- ! number of interior faces                       !
! nfabor           ! i  ! <-- ! number of boundary faces                       !
! nfml             ! i  ! <-- ! number of families (group classes)             !
! nprfml           ! i  ! <-- ! number of properties per family (group class)  !
! nnod             ! i  ! <-- ! number of vertices                             !
! lndfac           ! i  ! <-- ! size of nodfac indexed array                   !
! lndfbr           ! i  ! <-- ! size of nodfbr indexed array                   !
! ncelbr           ! i  ! <-- ! number of cells with faces on boundary         !
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! nphas            ! i  ! <-- ! number of phases                               !
! iviscf, b        ! e  ! --> ! "pointeur" sur viscf, viscb                    !
! idam             ! e  ! --> ! "pointeur" sur dam                             !
! icofbd           ! e  ! --> ! "pointeur" sur cofbdt                          !
! iw1              ! e  ! --> ! "pointeur" sur w1                              !
! iw2              ! e  ! --> ! "pointeur" sur w2                              !
! iw3              ! e  ! --> ! "pointeur" sur w3                              !
! igrarx,y,z       ! e  ! --> ! "pointeurs" sur grarox,y,z (iptlro=1)          !
! iwcf             ! e  ! --> ! "pointeur" sur wcf                             !
! ifinia           ! i  ! --> ! number of first free position in ia (at exit)  !
! ifinra           ! i  ! --> ! number of first free position in ra (at exit)  !
!__________________.____._____.________________________________________________.

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

integer          idbia0 , idbra0
integer          ndim   , ncelet , ncel  , nfac  , nfabor
integer          nfml   , nprfml
integer          nnod   , lndfac , lndfbr, ncelbr
integer          nvar   , nscal  , nphas
integer          iviscf , iviscb , idam  , icofbd, icofbr
integer          iw1    , iw2    , iw3
integer          igrarx, igrary, igrarz, iwcf
integer          ifinia , ifinra
integer          iptlro , ippmcf

integer          idebia, idebra
integer          iipmcf, iirocf

!===============================================================================

!---> INITIALISATION

idebia = idbia0
idebra = idbra0


!---> PLACE MEMOIRE RESERVEE AVEC DEFINITION DE IFINIA IFINRA

ifinia =       idebia

!  Ceci, pour le compressible, n'est pas tres satisfaisant,
!  mais on se calque sur IPTLRO (certains tableaux n'existent
!  pas toujours)

iipmcf = 0
if(ippmcf.ge.0) then
  iipmcf = 1
endif

iirocf = 0
if(iptlro.eq.1.or.ippmcf.ge.0) then
  iirocf = 1
endif

iviscf =       idebra
iviscb =       iviscf + nfac
idam   =       iviscb + nfabor
icofbd =       idam   + ncelet
iw1    =       icofbd + nfabor
iw2    =       iw1    + ncelet
iw3    =       iw2    + ncelet
icofbr =       iw3    + ncelet
igrarx =       icofbr + nfabor*iirocf
igrary =       igrarx + ncelet*iirocf
igrarz =       igrary + ncelet*iirocf
iwcf   =       igrarz + nfabor*iirocf
ifinra =       igrarz + ncelet*iipmcf

!---> VERIFICATION

call iasize('memdtv',ifinia)
!==========

call rasize('memdtv',ifinra)
!==========

return
end subroutine
