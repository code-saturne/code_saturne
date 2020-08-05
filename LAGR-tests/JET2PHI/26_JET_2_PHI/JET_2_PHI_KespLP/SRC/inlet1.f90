subroutine inlet1                               &
!================

 ( u11, v11, w11, q11, ep11, zz )

!===============================================================================
! Purpose:
! ----------

! Subroutine inlet1
! -----------------
!
!   User-defined subroutine for the validation test-case
!   "Particle-laden wall jet" (JET_2_PHI) for the validation
!   of Code_Saturne
!
!   Initialization of the inlet conditions of the fluid phase
!


!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
!    name          !type!mode !                   role                         !
!__________________!____!_____!________________________________________________!
!  u11             ! r  ! --> ! vertical component of the fluid velocity       !
!  v11             ! r  ! --> ! tangential component of the fluid velocity (=0)!
!  w11             ! r  ! --> ! transverse component of the fluid velocity     !
!  q11             ! r  ! --> ! turbulent energy                               !
!  ep11            ! r  ! --> ! turbulent dissipation                          !
!  zz              ! r  ! <-- ! transverse coordinate                          !
!__________________!____!_____!________________________________________________!
!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================


!===============================================================================
!     Modules files
!===============================================================================

use paramx
use cstnum
use cstphy
use entsor

!===============================================================================
implicit none

! Arguments

double precision u11, v11, w11, q11, ep11, zz

! Local variables

double precision di
double precision zi(1:40), ui(1:40), wi(1:40)
double precision ufi(1:40), wfi(1:40), uvi(1:40)
double precision rlambd, reynol, uetoil, vv, dh
double precision vx,vz,vpx,vpz
double precision rnu, uvx
double precision ckarm

integer i, it, itmx, itrouv

!===============================================================================

!===============================================================================
! 1.  Data initializations with experimental measurements
 !===============================================================================

!     Transverse coordinate

data zi/  0.8e-3, 1.0e-3, 1.2e-3, 1.4e-3, 1.6e-3, 1.8e-3, 2.0e-3, &
          2.2e-3, 2.4e-3, 2.6e-3, 2.8e-3, 3.0e-3, 3.2e-3, 3.4e-3, &
          3.6e-3, 3.8e-3, 4.0e-3, 4.2e-3, 4.4e-3, 4.6e-3, 4.8e-3, &
          5.0e-3, 5.2e-3, 5.4e-3, 5.6e-3, 5.8e-3, 6.0e-3, 6.2e-3, &
          6.5e-3, 7.0e-3, 7.5e-3, 8.0e-3, 8.5e-3, 9.0e-3, 9.5e-3, &
          10.e-3, 11.e-3, 16.e-3, 21.e-3, 31.e-3 /

!     Vertical mean velocity

data ui/  8.568, 9.069, 9.521, 9.828, 9.850, 9.953, 9.944, 9.990, &
          9.958, 9.906, 9.798, 9.622, 9.494, 9.119, 8.578, 7.896, &
          7.081, 5.973, 4.821, 3.639, 2.001, 1.136, 0.821, 0.723, &
          0.736, 0.860, 1.037, 1.113, 1.309, 1.572, 1.711, 1.799, &
          1.882, 1.920, 1.948, 1.954, 1.995, 2.020, 2.011, 1.995 /

!     Transverse mean velocity

data wi/  0.290, 0.290, 0.324, 0.337, 0.345, 0.353, 0.361, 0.371, &
          0.365, 0.360, 0.355, 0.349, 0.349, 0.348, 0.328, 0.310, &
          0.278, 0.237, 0.179, 0.150, 0.069,-0.034,-0.024,-0.053, &
         -0.063,-0.077,-0.086,-0.095,-0.060,-0.052,-0.054,-0.066, &
         -0.045,-0.048,-0.054,-0.056,-0.046,-0.037,-0.041,-0.031 /

!     Fluctuation of the vertical velocity

data ufi/ 0.733, 0.640, 0.455, 0.271, 0.211, 0.179, 0.170, 0.166, &
          0.189, 0.207, 0.254, 0.329, 0.402, 0.512, 0.531, 0.579, &
          0.592, 0.532, 0.513, 0.418, 0.296, 0.180, 0.140, 0.169, &
          0.162, 0.189, 0.184, 0.179, 0.183, 0.184, 0.169, 0.162, &
          0.137, 0.145, 0.131, 0.142, 0.136, 0.147, 0.141, 0.151 /

!     Fluctuation of the transverse velocity

data wfi/ 0.168, 0.171, 0.139, 0.125, 0.122, 0.121, 0.119, 0.118, &
          0.120, 0.118, 0.121, 0.115, 0.136, 0.143, 0.143, 0.142, &
          0.156, 0.149, 0.159, 0.142, 0.174, 0.146, 0.116, 0.149, &
          0.136, 0.143, 0.142, 0.132, 0.117, 0.114, 0.107, 0.114, &
          0.097, 0.112, 0.111, 0.130, 0.124, 0.140, 0.140, 0.155 /

!  Shear-stress (currently not used)
!
data uvi/ 0.0084,  0.0048,  -0.0037,  0.0017,  0.0045,            &
          0.0043,  0.0049,   0.0043,  0.0058,  0.0077,            &
          0.0096,  0.0131,   0.0151,  0.0199,  0.0240,            &
          0.0244,  0.0224,   0.0307,  0.0224,  0.0167,            &
          0.0286,  0.0191,   0.0098,  0.0164,  0.0127,            &
          0.0144,  0.0160,   0.0139,  0.0103,  0.0079,            &
          0.0078,  0.0103,   0.0074,  0.0110,  0.0103,            &
          0.0139,  0.0126,   0.0170,  0.0163,  0.0203 /

!===============================================================================
! 2.  Initializations of the constants
!===============================================================================

itmx  = 40
di    = 5.d-3
rnu   = 1.6d-05
cmu   = 0.09d0
ckarm = 0.41d0
it = -1

!===============================================================================
! 3.  Interpolation
!===============================================================================

itrouv = 0

if ( zz.lt.zi(itmx) .and. zz.gt.zi(1) )then
  do i = 1,39
    if ( zz.ge.zi(i) .and. zz.lt.zi(i+1) ) then
      if ( itrouv.eq.1 ) then
        write(nfecra,9001)
        call csexit(1)
      else
        it = i
        itrouv = 1
      endif
    endif
  enddo
else if ( zz.ge.zi(itmx) ) then
  it = itmx
else if ( zz.le.zi(1) )then
  it = 1
else
  write(nfecra,9002)
  call csexit(1)
endif

!===============================================================================
! 4.  Calculation of the velocity components and of the turbulent energy
!===============================================================================

if (  zz.ge.zi(itmx) .or. zz.le.zi(1) ) then
   vx  = ui(it)
   vz  = wi(it)
   vpx = ufi(it)
   vpz = wfi(it)
   uvx = uvi(it)
else
  vx = ui(it)                                                     &
      +(zz-zi(it))*(ui(it+1)-ui(it))/(zi(it+1)-zi(it))
  vz = wi(it)                                                     &
      +(zz-zi(it))*(wi(it+1)-wi(it))/(zi(it+1)-zi(it))
  vpx = ufi(it)                                                   &
      +(zz-zi(it))*(ufi(it+1)-ufi(it))/(zi(it+1)-zi(it))
  vpz = wfi(it)                                                   &
      +(zz-zi(it))*(wfi(it+1)-wfi(it))/(zi(it+1)-zi(it))
  uvx = uvi(it)                                                   &
      +(zz-zi(it))*(uvi(it+1)-uvi(it))/(zi(it+1)-zi(it))
endif

u11 = vx
v11 = 0.d0
w11 = vz
q11 = 0.75d0 * ( vpx**2 + vpz**2 )

!===============================================================================
! 5.  Calculation of the dissipation by means of a correlation
!===============================================================================

if ( zz.lt.di ) then

  vv = sqrt( u11*u11 + v11*v11 + w11*w11 )
  dh = 5.d-3
  if ( vv .lt. 1.d-10 ) then
    write(nfecra,9003)
    call csexit(1)
  else
    reynol = vv * dh / rnu
    if (  reynol .lt. 3.d+4  ) then
      rlambd = 0.3164d0 / ( reynol ** 0.25d0 )
    else
      rlambd = 0.184d0  / ( reynol ** 0.2d0 )
    endif
    uetoil = vv * sqrt ( rlambd / 8.d0 )
    ep11 = cmu * q11 * q11 / (ckarm*uetoil*dh*0.1)
  endif

else

  vv = sqrt ( u11 * u11 + v11 * v11 + w11 * w11 )
  dh = 9.667d-3
  if ( vv.lt.1.d-10 ) then
    write(nfecra,9001)
    call csexit(1)
  else
    reynol = vv * dh / rnu
    if (  reynol .lt. 3.d+4  ) then
      rlambd = 0.3164d0 / ( reynol ** 0.25d0 )
    else
      rlambd = 0.184d0  / ( reynol ** 0.2d0 )
    endif
    uetoil = vv * sqrt ( rlambd / 8.d0 )
    ep11 = cmu * q11 * q11 / (ckarm*uetoil*dh*0.1d0)
  endif

endif

!===============================================================================

!--------
! Formats
!--------


 9001 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ CAUTION: STOP OF THE VALIDATION TEST-CASE               ',/,&
'@    =========   (INLET1)                                    ',/,&
'@                                                            ',/,&
'@    PROBLEM IN THE BOUNDARY CONDITIONS                      ',/,&
'@                                                            ',/,&
'@  The given coordinate is too large                         ',/,&
'@                                                            ',/,&
'@  The calculation cannot be run.                            ',/,&
'@                                                            ',/,&
'@  Please check the "inlet1" user-defined subroutine         ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 9002 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ CAUTION: STOP OF THE VALIDATION TEST-CASE               ',/,&
'@    =========   (INLET1)                                    ',/,&
'@                                                            ',/,&
'@    PROBLEM IN THE BOUNDARY CONDITIONS                      ',/,&
'@                                                            ',/,&
'@  The interpolation has failed                              ',/,&
'@                                                            ',/,&
'@  The calculation cannot be run.                            ',/,&
'@                                                            ',/,&
'@  Please check the "inlet1" user-defined subroutine         ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 9003 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ CAUTION: STOP OF THE VALIDATION TEST-CASE               ',/,&
'@    =========   (INLET1)                                    ',/,&
'@                                                            ',/,&
'@    PROBLEM IN THE BOUNDARY CONDITIONS                      ',/,&
'@                                                            ',/,&
'@  The calculated velocity is equal to zero.                 ',/,&
'@                                                            ',/,&
'@  The calculation cannot be run.                            ',/,&
'@                                                            ',/,&
'@  Please check the "inlet1" user-defined subroutine         ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)



!----
! End
!----

return

end subroutine inlet1
