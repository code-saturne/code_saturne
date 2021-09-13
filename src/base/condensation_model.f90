subroutine condensation_model(nvar, nfbpcd, ifbpcd, izzftcd, &
                              spcond, hpcond) 

use ppincl, only: icondb_model, icondb_regime

!===============================================================================
implicit none

!> Arguments

integer          nvar, nfbpcd, ifbpcd(nfbpcd), izzftcd(nfbpcd)

double precision spcond(nfbpcd,nvar)
double precision hpcond(nfbpcd)

if (icondb_model == 0) then
  call condensation_copain_model &
( nvar, nfbpcd, ifbpcd, izzftcd, &
  spcond ,hpcond, icondb_regime )
else if (icondb_model == 1) then
  call condensation_copain_benteboula_dabbene_model&
( nvar, nfbpcd, ifbpcd, izzftcd, &
  spcond ,hpcond, icondb_regime )
else if (icondb_model == 2) then
  call condensation_uchida_model&
( nvar, nfbpcd, ifbpcd, izzftcd, &
  spcond ,hpcond, icondb_regime )
else if (icondb_model == 3) then
  call condensation_dehbi_model&
( nvar, nfbpcd, ifbpcd, izzftcd, &
  spcond ,hpcond, icondb_regime )
endif

end subroutine condensation_model
