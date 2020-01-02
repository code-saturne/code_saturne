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

subroutine vorimp &
 ( ient   )

!===============================================================================
!  FONCTION  :
!  ---------

! IMPRESSION DES PARAMETRES DE LA METHODE DES VORTEX
! APRES INTERVENTION UTILISATEUR DANS usvort
!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! ient             ! e  ! <-- ! numero de l'entree                             !
! nnent            ! e  ! <-- ! nombre d'entrees utilisant des vortex          !
! nvort(nentmx)    ! te ! ->  ! nombre de vortex a chaque entree               !
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
use optcal
use vorinc

!===============================================================================

implicit none

! Arguments

integer          ient

! Local variables

integer          ipass
data             ipass /0/
save             ipass

!===============================================================================
! On ne peut pas faire ces impressions dans impini car la geometrie
!  des entrees n'est pas encore connue.

!===============================================================================
! 1. Methode vortex
!===============================================================================

ipass = ipass + 1

write(nfecra,9010) ient, nvort(ient), icvor(ient), icas(ient)
write(nfecra,9020) dir1(1,ient),dir1(2,ient),dir1(3,ient),        &
                   dir2(1,ient),dir2(2,ient),dir2(3,ient),        &
                   dir3(1,ient),dir3(2,ient),dir3(3,ient),        &
                   cen(1,ient) ,cen(2,ient) ,cen(3,ient)
if(icas(ient).eq.1) then
  write(nfecra,9030) iclvor(1,ient),iclvor(2,ient),               &
                     iclvor(3,ient),iclvor(4,ient)
endif
write(nfecra,9040) ymin(ient),ymax(ient),zmin(ient),zmax(ient)
if(icas(ient).eq.1) then
  write(nfecra,9050) lly(ient),llz(ient)
endif
if(icas(ient).eq.2) then
  write(nfecra,9060) lld(ient)
endif
write(nfecra,9070) itlivo(ient)
if(itlivo(ient).eq.1) then
  write(nfecra,9080) tlimvo(ient)
endif
write(nfecra,9090) isgmvo(ient)
if(isgmvo(ient).eq.1) then
  write(nfecra,9100) xsgmvo(ient)
endif
write(nfecra,9110) idepvo(ient)
if(idepvo(ient).eq.1) then
  write(nfecra,9120) ud(ient)
endif
if(icas(ient).eq.1.or.icas(ient).eq.2.or.icas(ient).eq.3) then
  write(nfecra,9130) ndat(ient)
elseif(icas(ient).eq.4) then
  write(nfecra,9140) udebit(ient),kdebit(ient),edebit(ient)
endif

#if defined(_CS_LANG_FR)

 9010 format(                                                           &
'-------------                                                ',/,&
'  -- Entree : ',I10                                            /,&
'-------------                                                ',/,&
'       NVORT  = ',4X,I10,    ' (Nombre de vortex            )',/,&
'       ICVOR  = ',4X,I10,   ' (Nombre de faces a l''entree  )',/,&
'       ICAS   = ',4X,I10,    ' (1 : conduite rectangulaire   ',/,&
'                ',14X,       '  2 : conduite circulaire      ',/,&
'                ',14X,       '  3 : sans CL mais avec fichier',/,&
'                ',14X,       '  4 : sans CL ni fichier      )',/)
 9020 format(                                                           &
' --- Directions principales du repere local                  ',/,&
'                  ---- X ----    ---- Y ----    ---- Z ----  ',/,&
'       DIR1  = ', E14.5,' ',E14.5,' ',E14.5                    /,&
'       DIR2  = ', E14.5,' ',E14.5,' ',E14.5                    /,&
'       DIR3  = ', E14.5,' ',E14.5,' ',E14.5                    /,&
'                                                             ',/,&
' --- Coordonnees du centre de l''entree                      ',/,&
'       CEN   = ', E14.5,' ',E14.5,' ',E14.5,                   /)
 9030 format(                                                           &
' --- Conditions aux limites dans le repere local             ',/,&
' Plan y = -LLY/2 ',4X,I10,    ' (1 : paroi                   ',/,&
' Plan z =  LLZ/2 ',4X,I10,    '  2 : symetrie                ',/,&
' Plan y =  LLY/2 ',4X,I10,    '  3 : periodicite            )',/,&
' Plan z = -LLZ/2 ',4X,I10,    '                              ',/)
 9040 format(                                                           &
' --- Dimensions de l''entree dans le repere local            ',/,&
'                  ---- min ----    ---- max ----             ',/,&
'       Y       = ',E14.5,'  ',E14.5,'                        ',/,&
'       Z       = ',E14.5,'  ',E14.5,'                        ',/)
 9050 format(                                                           &
'       LLY     = ',E14.5,    ' (longueur de la conduite dans ',/,&
'       LLZ     = ',E14.5,    '  les directions DIR1 et DIR2) ',/)
 9060 format(                                                           &
'       LLD     = ',E14.5,    ' (diametre de la conduite    ) ',/)
 9070 format(                                                           &
' --- Duree de vie des vortex                                 ',/,&
'       ITLIVO  = ',4X,I10,   ' (1 : constante                ',/,&
'                 ',14X,      '  2 : en k^(3/2).U/epsilon   ) ',/)
 9080 format(                                                           &
'       TLIMVO  = ',E14.5,    ' (1 : duree de vie imposee   ) ',/)
 9090 format(                                                           &
' --- Taille des vortex                                       ',/,&
'       ISGMVO  = ',4X,I10,   ' (1 : taille constante         ',/,&
'                 ',14X,      '  2 : en k^(3/2)/epsilon       ',/,&
'                 ',14X,      '  2 : en max[nu.k/eps,200.Lk]) ',/)
 9100 format(                                                           &
'       XSGMVO  = ',E14.5,    ' (1 : taille imposee         ) ',/)
 9110 format(                                                           &
' --- Marche en temps                                         ',/,&
'       IDEPVO  = ',4X,I10,   ' (1 : deplacement aleatoire    ',/,&
'                 ',14X,      '  2 : convection des vortex    ',/)
 9120 format(                                                           &
'       UD      = ',E14.5,    ' (1 : vit. de deplacement max) ',/)
 9130 format(                                                           &
' --- Fichier de donnees                                      ',/,&
'       NDAT    = ',4X,I10,   ' (Nombre de lignes du fichier )',/)
 9140 format(                                                           &
' --- Donnees a l''entree                                     ',/,&
'       UDEBIT  = ',E14.5,    ' (vitesse debitante imposee)   ',/,&
'       KDEBIT  = ',E14.5,    ' (energie cinetique imposee)   ',/,&
'       EDEBIT  = ',E14.5,    ' (dissipation imposee)         ',/)

#else

 9010 format(                                                           &
'-----------                                                  ',/,&
'  -- Inlet: ',I10                                              /,&
'-----------                                                  ',/,&
'       NVORT  = ',4X,I10,    ' (Number of vortices          )',/,&
'       ICVOR  = ',4X,I10,    ' (Number of faces at the inlet)',/,&
'       ICAS   = ',4X,I10,    ' (1 : rectangular duct         ',/,&
'                ',14X,       '  2 : pipe                     ',/,&
'                ',14X,       '  3 : wo BC but with file      ',/,&
'                ',14X,       '  4 : wo BC neither file      )',/)
 9020 format(                                                           &
' --- Main directions for the local frame                     ',/,&
'                  ---- X ----    ---- Y ----    ---- Z ----  ',/,&
'       DIR1  = ', E14.5,' ',E14.5,' ',E14.5                    /,&
'       DIR2  = ', E14.5,' ',E14.5,' ',E14.5                    /,&
'       DIR3  = ', E14.5,' ',E14.5,' ',E14.5                    /,&
'                                                             ',/,&
' --- Inlet center coordinates                                ',/,&
'       CEN   = ', E14.5,' ',E14.5,' ',E14.5,                   /)
 9030 format(                                                           &
' --- Boundary conditions in the local frame                  ',/,&
' Y plane = -LLY/2 ',4X,I10,    ' (1 : wall                   ',/,&
' Z plane =  LLZ/2 ',4X,I10,    '  2 : symmetry               ',/,&
' Y plane =  LLY/2 ',4X,I10,    '  3 : periodicity           )',/,&
' Z plane = -LLZ/2 ',4X,I10,    '                             ',/)
 9040 format(                                                           &
' --- Inlet dimensions in the local framae                    ',/,&
'                  ---- min ----    ---- max ----             ',/,&
'       Y       = ',E14.5,'  ',E14.5,'                        ',/,&
'       Z       = ',E14.5,'  ',E14.5,'                        ',/)
 9050 format(                                                           &
'       LLY     = ',E14.5,    ' (duct length in the directions',/,&
'       LLZ     = ',E14.5,    '  DIR1 and DIR2              ) ',/)
 9060 format(                                                           &
'       LLD     = ',E14.5,    ' (pipe diameter              ) ',/)
 9070 format(                                                           &
' --- Vortices life time                                      ',/,&
'       ITLIVO  = ',4X,I10,   ' (1 : constant                 ',/,&
'                 ',14X,      '  2 : in k^(3/2).U/epsilon   ) ',/)
 9080 format(                                                           &
'       TLIMVO  = ',E14.5,    ' (1 : given life time        ) ',/)
 9090 format(                                                           &
' --- Vortices size                                           ',/,&
'       ISGMVO  = ',4X,I10,   ' (1 : constant                 ',/,&
'                 ',14X,      '  2 : in k^(3/2)/epsilon       ',/,&
'                 ',14X,      '  2 : in max[nu.k/eps,200.Lk]) ',/)
 9100 format(                                                           &
'       XSGMVO  = ',E14.5,    ' (1 : given size             ) ',/)
 9110 format(                                                           &
' --- Time marching                                           ',/,&
'       IDEPVO  = ',4X,I10,   ' (1 : random displacement      ',/,&
'                 ',14X,      '  2 : vortices convection      ',/)
 9120 format(                                                           &
'       UD      = ',E14.5,    ' (1 : max displacement vel.  ) ',/)
 9130 format(                                                           &
' --- Data file                                               ',/,&
'       NDAT    = ',4X,I10,   ' (Number of lines in the file )',/)
 9140 format(                                                           &
' --- Inlet date                                              ',/,&
'       UDEBIT  = ',E14.5,    ' (given velocity      )        ',/,&
'       KDEBIT  = ',E14.5,    ' (given kinetic energy)        ',/,&
'       EDEBIT  = ',E14.5,    ' (given dissipation   )        ',/)

#endif

!===============================================================================
! 2. FIN
!===============================================================================

return
end subroutine
