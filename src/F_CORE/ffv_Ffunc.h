// #################################################################
//
// FFV : Frontflow / violet Cartesian
//
// Copyright (c) 2012-2014 Institute of Industrial Science, The University of Tokyo.
// All right reserved.
//
// Copyright (c) 2012-2014 Advanced Institute for Computational Science, RIKEN.
// All right reserved.
//
// #################################################################

/** 
 * @file   ffv_Ffunc.h
 * @brief  FFV Fortran function Header
 * @author kero
 */

#include "FB_Define.h"

#ifndef _FFV_F_FUNC_H_
#define _FFV_F_FUNC_H_


extern "C" {
  
  
  //***********************************************************************************************
  // ffv_poisson.f90
  void poi_residual_ (double* res,
                      int* sz,
                      int* g,
                      REAL_TYPE* p,
                      REAL_TYPE* b,
                      int* bp,
                      double* flop);
  
  void poi_rhs_ (double* rhs,
                 REAL_TYPE* b,
                 int* sz,
                 int* g,
                 REAL_TYPE* s_0,
                 REAL_TYPE* s_1,
                 int* bp,
                 REAL_TYPE* dh,
                 REAL_TYPE* dt,
                 double* flop);
  
  void psor2sma_core_ (REAL_TYPE* p,
                       int* sz,
                       int* g,
                       int* ip,
                       int* color,
                       REAL_TYPE* omg,
                       double* res,
                       REAL_TYPE* b,
                       int* bp,
                       double* flop);
  
  void sma_comm_      (REAL_TYPE* p,
                       int* sz,
                       int* g,
                       int* col,
                       int* ip,
                       int* cf_sz,
                       REAL_TYPE* cf_x,
                       REAL_TYPE* cf_y,
                       REAL_TYPE* cf_z,
                       int* key,
                       int* nID);
  
  void sma_comm_wait_ (REAL_TYPE* p,
                       int* sz,
                       int* g,
                       int* col,
                       int* ip,
                       int* cf_sz,
                       REAL_TYPE* cf_x,
                       REAL_TYPE* cf_y,
                       REAL_TYPE* cf_z,
                       int* key);
  
  
  //***********************************************************************************************
  // ffv_vbc_outer.f90
  void pvec_vobc_wall_    (REAL_TYPE* wv, int* sz, int* g, REAL_TYPE* dh, REAL_TYPE* v00, REAL_TYPE* rei, REAL_TYPE* v,
                           int* bv, REAL_TYPE* vec, int* face, double* flop);
  
  void div_obc_drchlt_    (REAL_TYPE* div,
                           int* sz,
                           int* g,
                           int* face,
                           REAL_TYPE* v00,
                           int* bv,
                           REAL_TYPE* vec,
                           double* flop);
  
  
  //***********************************************************************************************
  // ffv_velocity_binary.f90
  void divergence_        (REAL_TYPE* dv, int* sz, int* g, REAL_TYPE* vc, int* bv, REAL_TYPE* v00, double* flop);
  void euler_explicit_    (REAL_TYPE* vc, int* sz, int* g, REAL_TYPE* dt, REAL_TYPE* v, int* bd, double* flop);
  void pvec_muscl_        (REAL_TYPE* wv, int* sz, int* g, REAL_TYPE* dh, int* c_scheme, REAL_TYPE* v00, REAL_TYPE* rei, REAL_TYPE* v, 
                           int* bv, int* bp, int* v_mode, REAL_TYPE* ut, int* wall_type, int* bd, REAL_TYPE* vcs_coef, double* flop);
  void update_vec_        (REAL_TYPE* v, REAL_TYPE* div, int* sz, int* g, REAL_TYPE* dt, REAL_TYPE* dh, REAL_TYPE* vc, REAL_TYPE* p, int* bp, int* bv, 
                           REAL_TYPE* v00, double* flop);

  //***********************************************************************************************
  // ffv_utility.f90
  void norm_v_div_max_ (double* ds,
                        int* sz,
                        int* g,
                        REAL_TYPE* div,
                        REAL_TYPE* coef,
                        int* bp,
                        double* flop);
  void find_vmax_         (REAL_TYPE* v_max, int* sz, int* g, REAL_TYPE* v00, REAL_TYPE* v, double* flop);
  
}

#endif // _FFV_F_FUNC_H_
