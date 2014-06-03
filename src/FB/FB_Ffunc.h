#ifndef _FB_F_FUNC_H_
#define _FB_F_FUNC_H_

// #################################################################
//
// CAERU Library
//
// Copyright (c) 2012-2014 Institute of Industrial Science, The University of Tokyo.
// All right reserved.
//
// Copyright (c) 2012-2014 Advanced Institute for Computational Science, RIKEN.
// All right reserved.
//
// #################################################################

/**
 * @file   FB_Ffunc.h
 * @brief  FlowBase Fortran function Header
 * @author kero
 */


extern "C" {
	// FB_util.f90
  
  void fb_delta_s_        (double* d,
                           int* sz,
                           int* g,
                           REAL_TYPE* sn,
                           REAL_TYPE* so,
                           int* bx,
                           double* flop);
  
  void fb_delta_v_        (double* d,
                           int* sz,
                           int* g,
                           REAL_TYPE* vn,
                           REAL_TYPE* vo,
                           int* bx,
                           double* flop);

  void fb_set_vector_     (REAL_TYPE* var, 
                           int* sz, 
                           int* g, 
                           REAL_TYPE* val, 
                           int* bv);

  void fb_write_sph_s_    (REAL_TYPE* s,
                           int* sz,
                           int* g,
                           char* fname,
                           int* step,
                           REAL_TYPE* time,
                           REAL_TYPE* org,
                           REAL_TYPE* pit,
                           int* d_type,
                           int* gs,
                           int* avs,
                           int* step_avr,
                           REAL_TYPE* time_avr);

  void fb_write_sph_v_    (REAL_TYPE* v,
                           int* sz,
                           int* g,
                           char* fname,
                           int* step,
                           REAL_TYPE* time,
                           REAL_TYPE* org,
                           REAL_TYPE* pit,
                           int* d_type,
                           int* gs,
                           int* avs,
                           int* step_avr,
                           REAL_TYPE* time_avr);
}

#endif // _FB_F_FUNC_H_
