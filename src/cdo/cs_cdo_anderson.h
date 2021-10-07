#ifndef __CS_CDO_ANDERSON_H__
#define __CS_CDO_ANDERSON_H__

/*============================================================================
 * Build an algebraic CDO face-based system for the Navier-Stokes equations
 * and solved it as one block (monolithic approach of the velocity-pressure
 * coupling), Anderson
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2021 EDF S.A.

  This program is free software; you can redistribute it and/or modify it under
  the terms of the GNU General Public License as published by the Free Software
  Foundation; either version 2 of the License, or (at your option) any later
  version.

  This program is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
  details.

  You should have received a copy of the GNU General Public License along with
  this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
  Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

/*----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

//   typedef struct {
// 
//   int                              mMax_aa;
//   int                              AAStart;
//   /* Set of tolerances of matrix conditioning number */
//   double                           droptol_;
//   /* Set the relaxation parameter */
//   double                           beta_;
//   /* Number of Anderson accelerations */
//   int                              n_aa_iter;
//   
//   bool                             aa_activated;
// 
//   cs_real_t                        fval_;
//   cs_real_t                        gold_;
//   cs_real_t                        fold_;
//   cs_real_t                        df_;
//
// 
//   Eigen::MatrixXd R_;
// 
// 
// } cs_cdo_anderson_t; 

static const int                        mMax_aa = 4;
static const int                        AAStart = 0;
static const int                        aa_activated = 0;
static const double                     droptol_ = 50;
static const double                     beta_ = 0.0;
int                              n_aa_iter;
int                              mAA_;
cs_real_t                        *fval_;
cs_real_t                        *gval;
cs_real_t                        *fold_;
cs_real_t                        *gold_;
cs_real_t                        *df_;
cs_real_t                        *DG_;
cs_real_t                        *Q_;
cs_sdm_t                         *R_;
cs_real_t                        *gamma_;
cs_real_t                        *tempGamma_;

/*============================================================================
 * Static inline public function prototypes
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
inline void aa_init(const int           mMax,
                    const cs_lnum_t     n_faces) {

  n_aa_iter = 0;
  mAA_ = 0;
  fval_ = NULL;
  fold_ = NULL;
  gold_ = NULL;
  df_   = NULL;
  DG_   = NULL;
  Q_    = NULL;
  R_    = NULL;
  gamma_ = NULL;
  tempGamma_ = NULL;
  gval  = NULL;

  BFT_MALLOC(fval_, n_faces, cs_real_t);
  BFT_MALLOC(fold_, n_faces, cs_real_t);
  BFT_MALLOC(gold_, n_faces, cs_real_t);
  BFT_MALLOC(df_, n_faces, cs_real_t);
  BFT_MALLOC(DG_, n_faces*mMax, cs_real_t);
  BFT_MALLOC(Q_, n_faces*mMax, cs_real_t);
  R_ = cs_sdm_square_create(mMax);
  cs_sdm_square_init(mMax, R_);
  BFT_MALLOC(gamma_, mMax, cs_real_t);
  BFT_MALLOC(tempGamma_, mMax, cs_real_t);
  cs_log_printf(CS_LOG_DEFAULT, " anderson   %d %d \n", mMax, n_faces);
  cs_log_printf(CS_LOG_DEFAULT, " anderson R_  %d %lf \n", R_->n_rows, R_->val[1]);
}

/*----------------------------------------------------------------------------*/

inline void aa_term() {

  BFT_FREE(fval_);
  BFT_FREE(fold_);
  BFT_FREE(gold_);
  BFT_FREE(df_);
  BFT_FREE(DG_);
  BFT_FREE(Q_);
  R_ = cs_sdm_free(R_);
  BFT_FREE(gamma_);
  BFT_FREE(tempGamma_);
}

/*----------------------------------------------------------------------------*/

 inline void qrdelete(cs_real_t *Q, cs_sdm_t *R, int n_faces, int mMax){
   
   int n_rows = R->n_rows;
   for (int i=0 ; i<mMax-1 ; i++)
   {  
     cs_log_printf(CS_LOG_DEFAULT, " n_rows   %d \n", n_rows);
     //int n_cols = R->n_cols;
     assert(R->n_cols == R->n_rows);
     const double temp = sqrt(R->val[i*n_rows+i+1]*R->val[i*n_rows+i+1]+R->val[(i+1)*n_rows+i+1]*R->val[(i+1)*n_rows+i+1]);
        
     const double c = R->val[i*n_rows+i+1]/temp;
     const double s = R->val[(i+1)*n_rows+i+1]/temp;
     
     R->val[i*n_rows+i+1] = temp;
     R->val[(i+1)*n_rows+i+1]=0.0;
     if (i<mMax-2)
     { //diff % matlab
       for (int j=i+2 ; j<mMax ; j++)
       {
	 const double temp0 = c*R->val[i*n_rows+j]+s*R->val[(i+1)*n_rows+j];
	 R->val[(i+1)*n_rows+j] = -s*R->val[i*n_rows+j]+c*R->val[(i+1)*n_rows+j];
	 R->val[i*n_rows+j] = temp0;
       }
     }
    cs_real_t *Q_i = Q + i*n_faces;
    cs_real_t *Q_ip1 = Q + (i+1)*n_faces;
    for (int l=0 ; l<n_faces ; l++){
      const double temp1 = c*Q_i[l] + s*Q_ip1[l];
      Q_ip1[l] = -s*Q_i[l] + c*Q_ip1[l];
      Q_i[l] = temp1;
    }
   }

   cs_sdm_t *temp2 = cs_sdm_square_create(R->n_rows) ;
   cs_sdm_square_init(R->n_rows, temp2);

   //diff % matlab the Q and R matrices are not resized
   // A block version should be used
   for (int i=0 ; i<mMax-1 ; i++){
     for (int j=i ; j<mMax-1 ; j++){
       temp2->val[i*n_rows+j] = R->val[i*n_rows+j+1];
     }
   }

   cs_sdm_copy(R, temp2);
   cs_real_t *Q_imax = Q + (mMax-1)*n_faces;
   for (int i=0 ; i<n_faces ; i++){
     Q_imax[i] = 0.0;
   }
}

/*----------------------------------------------------------------------------*/

// inline void shiftColumns( void ){
//     VECTOR * first=data_[0];
//     for (int j=0 ; j<colNumber_-1 ; j++) data_[j]=data_[j+1];
//     data_[colNumber_-1]=first;
//     colNumber_--;
//  }

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_CDO_ANDERSON_H__ */
