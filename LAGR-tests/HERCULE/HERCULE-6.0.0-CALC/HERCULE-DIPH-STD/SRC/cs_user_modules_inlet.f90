!-------------------------------------------------------------------------------

!                      Code_Saturne version 6.0.0-patch
!                      --------------------------
! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2019 EDF S.A.
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

!===============================================================================
! Purpose:
! -------

!> \file cs_user_modules.f90
!>
!> \brief User-defined module: it allows to create any user array.
!>
!> See \subpage cs_user_modules for examples.
!>
!> This file is compiled before all other user Fortran files.
!> To ensure this, it must not be renamed.
!>
!> The user may define an arbitrary number of modules here, even though
!> only one is defined in the example.
!
!> \cond DOXYGEN_SHOULD_SKIP_THIS

!-------------------------------------------------------------------------------

module user_module

implicit none

contains

      SUBROUTINE INLET(U1,V1,W1,AK,EP,RAY,R11,R12,R13,R22,R23,R33)                          
!                                                                               
! *********************************************************************         
! *  INITIALISATION DES CONDITIONS D'ENTREE DE LA PHASE GAZEUSE       *         
! *    D'UNE INJECTION AXISYMETRIQUE DE PARTICULES (DP= 60 MICRONS)   *         
! *         MESURES EXPERIMENTALES HERCULE REALISEES PAR TSU          *
! *                    VERSION DU  10/12/96                           *         
! *********************************************************************         
!                                                                               
      IMPLICIT NONE                                              
!                         
      INTEGER IT,ITMX                                             
      DOUBLE PRECISION U1,V1,W1,AK,EP,RAY                       
      DOUBLE PRECISION UI(20),VI(20),WI(20),RUI(20),RVI(20),RWI(20),XI(20)  
      DOUBLE PRECISION UU,VV,WW,XLM
      DOUBLE PRECISION REYNOL,DH,RLAMBD,UETOIL   
      DOUBLE PRECISION D1,D2,D3,D4      
      DOUBLE PRECISION R11,R12,R13,R22,R23,R33
!
      DOUBLE PRECISION PY,RNU,CMU,CKARM                                   
!                                                                               
!  DISTANCE RADIALE M                                                           
!                                                                                             
!                                                                               
      DATA XI/  0.E-3,   2.E-3,   4.E-3,   6.E-3,   8.E-3,  10.E-3,   &                   
               76.E-3,  80.E-3,  84.E-3,  88.E-3,  92.E-3,  96.E-3,   &                  
              100.E-3, 104.E-3, 108.E-3, 112.E-3, 116.E-3, 120.E-3,   &                
              124.E-3, 127.E-3/                                                         
!                                                                               
!  VITESSE RADIALE EN M/S                                                       
!                                                                               
      DATA UI/ 0.043071,  0.056412,  0.071609,  0.056304,  0.067233,  &                    
               0.057133, -0.078698, -0.062107, -0.043911, -0.060412,  &                    
              -0.039152, -0.041603, -0.047923, -0.034075, -0.043817,  &                    
              -0.028827, -0.029241, -0.015086,  0.017707,  0.045791/    
!
!!    DATA UI/ 0., 0., 0., 0., 0.,                      
!!   *         0., 0., 0., 0., 0.,                      
!!   *         0., 0., 0., 0., 0.,                      
!!   *         0., 0., 0., 0., 0./                   
!                                                                               
!  VITESSE AXIALE EN M/S                                                        
!                                                                               
      DATA VI/ 3.989888, 3.955918, 3.818365, 3.615492, 3.235091,   &                   
               1.15302 , 4.609382, 5.050922, 5.349239, 5.663382,   &                   
               5.774033, 5.852311, 5.957135, 6.0781  , 6.075214,   &                   
               5.910082, 5.932579, 5.929479, 5.824071, 5.773876/                      
!                                                                               
!                                                                               
!  VITESSE CIRCONFERENTIELLE EN M/S                                             
!                                                                               
      DATA WI/ 0., 0., 0., 0., 0.,   &                   
               0., 0., 0., 0., 0.,   &                  
               0., 0., 0., 0., 0.,   &                   
               0., 0., 0., 0., 0./                     
!                                                                               
!                                                                               
!  ECART TYPE DE LA VITESSE RADIALE EN M/S                                      
!                                                                               
      DATA RUI/ 0.165755, 0.167099, 0.188935, 0.196318, 0.199175,   &                  
                0.192527, 0.297634, 0.26486 , 0.253143, 0.227446,   &                  
                0.198263, 0.180068, 0.16104 , 0.157479, 0.145587,   &                  
                0.16021 , 0.15885 , 0.186978, 0.201831, 0.214657/                    
!                                                                               
!  ECART TYPE DE LA VITESSE AXIALE EN M/S                                       
!                                                                               
      DATA RVI/ 0.2206  , 0.235237, 0.280546, 0.345076, 0.480505,   &                  
                0.462576, 0.620767, 0.412637, 0.378357, 0.33724 ,   &                  
                0.309257, 0.268617, 0.26992 , 0.248661, 0.239869,   &                  
                0.24758 , 0.238192, 0.31793 , 0.337416, 0.338541/                    
!                                                                               
!                                                                               
!  ECART TYPE DE LA VITESSE TANGENTIELLE EN M/S                                 
!                                                                               
      DATA RWI/ 0.165755, 0.167099, 0.188935, 0.196318, 0.199175,  &                   
                0.192527, 0.297634, 0.26486 , 0.253143, 0.227446,  &                   
                0.198263, 0.180068, 0.16104 , 0.157479, 0.145587,  &                   
                0.16021 , 0.15885 , 0.186978, 0.201831, 0.214657/                
!                                                                               
!  INITIALISATION                                                               
! ****************                                                              
!                                                                                       
        U1= 0.                                                               
        V1= 0.                                                               
        W1= 0. 
        UU= 0.
        VV= 0.
        WW= 0.                                                              
!                                                                               
        AK= 1.E-9   
!
        PY = 3.14156
        RNU = 1.6e-05
        CMU   = 0.09
        CKARM = 0.41                                     
!                                                                               
!  INTERPOLATION                                                                
! ***************                                                               
!                                                                               
      ITMX= 19                                                                  
!                                                   
        DO 200 IT=1,ITMX                                                        
          IF( RAY.GE.XI(IT) .AND. RAY.LE.XI(IT+1) ) THEN                    
            U1= UI(IT)+(RAY-XI(IT))*(UI(IT+1)-UI(IT))/(XI(IT+1)-XI(IT))           
            V1= VI(IT)+(RAY-XI(IT))*(VI(IT+1)-VI(IT))/(XI(IT+1)-XI(IT))           
            W1= WI(IT)+(RAY-XI(IT))*(WI(IT+1)-WI(IT))/(XI(IT+1)-XI(IT))           
!                                                                               
            UU= RUI(IT)+(RAY-XI(IT))*(RUI(IT+1)-RUI(IT))/(XI(IT+1)-XI(IT))         
            VV= RVI(IT)+(RAY-XI(IT))*(RVI(IT+1)-RVI(IT))/(XI(IT+1)-XI(IT))         
            WW= RWI(IT)+(RAY-XI(IT))*(RWI(IT+1)-RWI(IT))/(XI(IT+1)-XI(IT))         
!    
            UU= UU*UU                                                  
            VV= VV*VV                                                  
            WW= WW*WW                                                  
!                                                                               
            AK= ( UU +VV +WW )/2.                                                                                                            
          ENDIF 
 200  CONTINUE
          IF ( RAY.GT.XI(ITMX) ) THEN
            U1= UI(ITMX)+(RAY-XI(ITMX))*(UI(ITMX+1)-UI(ITMX))/(XI(ITMX+1)-XI(ITMX))           
            V1= VI(ITMX)+(RAY-XI(ITMX))*(VI(ITMX+1)-VI(ITMX))/(XI(ITMX+1)-XI(ITMX))       
            W1= WI(ITMX)+(RAY-XI(ITMX))*(WI(ITMX+1)-WI(ITMX))/(XI(ITMX+1)-XI(ITMX))           
!                                                                               
            UU= RUI(ITMX)+(RAY-XI(ITMX))*(RUI(ITMX+1)-RUI(ITMX))/(XI(ITMX+1)-XI(ITMX))         
            VV= RVI(ITMX)+(RAY-XI(ITMX))*(RVI(ITMX+1)-RVI(ITMX))/(XI(ITMX+1)-XI(ITMX))         
            WW= RWI(ITMX)+(RAY-XI(ITMX))*(RWI(ITMX+1)-RWI(ITMX))/(XI(ITMX+1)-XI(ITMX))         
!    
            UU= UU*UU                                                  
            VV= VV*VV                                                  
            WW= WW*WW                                                  
!                                                                               
            AK= ( UU +VV +WW )/2.                                              
            R11 = UU
            R12 = 0.
            R13 = 0.
            R22 = VV
            R23 = 0.
            R33 = WW
          ENDIF                           
!    
        IF ( RAY .LE. 0.01   ) THEN
          DH = 0.02
        ELSE
          DH = 0.150
        ENDIF
        VV = SQRT ( U1*U1+V1*V1+W1*W1 )
        IF ( VV .LT. 1.E-10 ) THEN
          WRITE(*,*) ' NORME VITESSE NUL  RAYON = ',RAY
          STOP
        ELSE 
          REYNOL = VV * DH / RNU
          IF (  REYNOL .LT. 3.E+4  ) THEN
            RLAMBD = 0.3164 / ( REYNOL ** 0.25 )
          ELSE
            RLAMBD = 0.184  / ( REYNOL ** 0.2 )
          ENDIF
          UETOIL = VV * SQRT ( RLAMBD / 8. )
          EP = CMU * AK * AK / (CKARM*UETOIL*DH*0.1)
        ENDIF      
!                                                           
!! on met en commentaire car on ne sait pas d'ou ca vient AD-NP(3/2/04)
!!        D1=  20.E-3 
!!        D2= 150.E-3                                                               
!!        D3= 294.E-3        
!!        IF ( RAY.LE.(D1+D2)/4. ) THEN
!!          XLM = 0.03*D1
!!        ELSE
!!          XLM = 0.03*(D3-D2)
!!        ENDIF
!!        EP=  CMU*AK**1.5/XLM   
!! on met en commentaire car on ne sait pas d'ou ca vient AD-NP(3/2/04)
!                                                           
!                                                 
      END SUBROUTINE INLET                                                                       


end module user_module

!-------------------------------------------------------------------------------

!> (DOXYGEN_SHOULD_SKIP_THIS) \endcond
