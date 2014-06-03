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
 * @file   ffv_Alloc.C
 * @brief  FFV Class
 * @author kero
 */

#include "ffv.h"

#ifdef _REAL_IS_DOUBLE_
#define AllocRealS3D(guide) AllocDoubleS3D(guide)
#define AllocRealV3D(guide) AllocDoubleV3D(guide)
#else
#define AllocRealS3D(guide) AllocFloatS3D(guide)
#define AllocRealV3D(guide) AllocFloatV3D(guide)
#endif


// 主計算部分に用いる配列のアロケーション
void FFV::allocArray_Main(double &total)
{
  double mc = (double)(size[0] * size[1] * size[2]);
  size_t mt = (size[0]+2*guide) * (size[1]+2*guide) *(size[2]+2*guide);

  // d_v
  if ( !(d_v = paraMngr->AllocRealV3D(guide)) ) Exit(0);
  total+= mc * (double)sizeof(REAL_TYPE) * 3.0;
  memset(d_v, 0, sizeof(REAL_TYPE)*mt*3);
  
  // d_vc
  if ( !(d_vc = paraMngr->AllocRealV3D(guide)) ) Exit(0);
  total+= mc * (double)sizeof(REAL_TYPE) * 3.0;
  memset(d_vc, 0, sizeof(REAL_TYPE)*mt*3);
  
  // d_v0
  if ( !(d_v0 = paraMngr->AllocRealV3D(guide)) ) Exit(0);
  total+= mc * (double)sizeof(REAL_TYPE) * 3.0;
  memset(d_v0, 0, sizeof(REAL_TYPE)*mt*3);
  
  // d_wv
  if ( !(d_wv = paraMngr->AllocRealV3D(guide)) ) Exit(0);
  total+= mc * (double)sizeof(REAL_TYPE) * 3.0;
  memset(d_wv, 0, sizeof(REAL_TYPE)*mt*3);
  
  // d_wo
  if ( !(d_wo = paraMngr->AllocRealV3D(guide)) ) Exit(0);
  total+= mc * (double)sizeof(REAL_TYPE) * 3.0;
  memset(d_wo, 0, sizeof(REAL_TYPE)*mt*3);
  
  // d_p
  if ( !(d_p = paraMngr->AllocRealS3D(guide)) ) Exit(0);
  total+= mc * (double)sizeof(REAL_TYPE);
  memset(d_p, 0, sizeof(REAL_TYPE)*mt);
  
  // d_p0
  if ( !(d_p0 = paraMngr->AllocRealS3D(guide)) ) Exit(0);
  total+= mc * (double)sizeof(REAL_TYPE);
  memset(d_p0, 0, sizeof(REAL_TYPE)*mt);
  
  // d_sq
  if ( !(d_sq = paraMngr->AllocRealS3D(guide)) ) Exit(0);
  total+= mc * (double)sizeof(REAL_TYPE);
  memset(d_sq, 0, sizeof(REAL_TYPE)*mt);
  
  // d_dv
  if ( !(d_dv = paraMngr->AllocRealS3D(guide)) ) Exit(0);
  total+= mc * (double)sizeof(REAL_TYPE);
  memset(d_dv, 0, sizeof(REAL_TYPE)*mt);
  
  // d_b
  if ( !(d_b = paraMngr->AllocRealS3D(guide)) ) Exit(0);
  total+= mc * (double)sizeof(REAL_TYPE);
  memset(d_b, 0, sizeof(REAL_TYPE)*mt);
}



// 前処理に用いる配列のアロケーション
void FFV::allocArray_Prep(double &prep, double &total)
{
  double mc = (double)(size[0] * size[1] * size[2]);
  
  // d_ws
  if ( !(d_ws = paraMngr->AllocRealS3D(guide)) ) Exit(0);
  prep += mc * (double)sizeof(REAL_TYPE);
  total+= mc * (double)sizeof(REAL_TYPE);
  
  // d_mid
  if ( !(d_mid = paraMngr->AllocIntS3D(guide)) ) Exit(0);
  prep += mc * (double)sizeof(int);
  
  // d_bcd
  if ( !(d_bcd = paraMngr->AllocIntS3D(guide)) ) Exit(0);
  prep += mc * (double)sizeof(int);
  total+= mc * (double)sizeof(int);
  
  // d_bcp
  if ( !(d_bcp = paraMngr->AllocIntS3D(guide)) ) Exit(0);
  prep += mc * (double)sizeof(int);
  total+= mc * (double)sizeof(int);
  
  // d_bcv
  if ( !(d_bcv = paraMngr->AllocIntS3D(guide)) ) Exit(0);
  prep += mc * (double)sizeof(int);
  total+= mc * (double)sizeof(int);
}


// SOR2SMAのバッファ確保
void FFV::allocate_SOR2SMA_buffer(double &total)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  
  cf_sz[0] = (jx+1) * (kx+1) / 2; // バッファサイズ
  cf_sz[1] = (kx+1) * (ix+1) / 2; // +1はマージン
  cf_sz[2] = (ix+1) * (jx+1) / 2;
  
  size_t n1 = cf_sz[0]*4;
  size_t n2 = cf_sz[1]*4;
  size_t n3 = cf_sz[2]*4;
  
  if( (cf_x = new REAL_TYPE[n1]) == NULL ) Exit(0);
  if( (cf_y = new REAL_TYPE[n2]) == NULL ) Exit(0);
  if( (cf_z = new REAL_TYPE[n3]) == NULL ) Exit(0);
  
  memset(cf_x, 0, sizeof(REAL_TYPE)*n1);
  memset(cf_y, 0, sizeof(REAL_TYPE)*n2);
  memset(cf_z, 0, sizeof(REAL_TYPE)*n3);
  
  total += (double)( (n1+n2+n3)*sizeof(REAL_TYPE) );
}


// 主計算に用いる配列の確保
void FFV::allocate_Main(double &total)
{
  allocArray_Main(total);
}

