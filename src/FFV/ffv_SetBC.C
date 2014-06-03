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
 * @file   ffv_SetBC.C
 * @brief  FFV BC Class
 * @author kero
 */

#include "ffv_SetBC.h"


// #################################################################
// 外部境界条件リストから速度境界条件の成分を取り出す
REAL_TYPE SetBC3D::extractVel_OBC(const int n, REAL_TYPE* vec, const REAL_TYPE tm, const REAL_TYPE* v00, double& flop)
{
  REAL_TYPE vel = obc[n].ca[CompoList::bias] * v00[0];
  vec[0] = obc[n].nv[0] * vel;
  vec[1] = obc[n].nv[1] * vel;
  vec[2] = obc[n].nv[2] * vel;
  
  return vel;
}


// #################################################################
// 速度境界条件による速度の発散の修正ほか
// 外部境界面のdiv(u)の修正時に領域境界の流量などのモニタ値を計算し，BoundaryOuterクラスに保持 > 反復後にDomainMonitor()で集約
// avr[]のインデクスに注意 (Fortran <-> C)
void SetBC3D::mod_div(REAL_TYPE* dv, int* bv, REAL_TYPE tm, REAL_TYPE* v00, double& flop)
{
  REAL_TYPE vec[3];
  int st[3], ed[3];
  int typ=0;
  int gd = guide;
  double fcount = 0.0;
  REAL_TYPE aa[2];
  
  // 外部境界条件による修正
  for (int face=0; face<NOFACE; face++) {
    typ = obc[face].get_Class();
    
    // 計算領域の最外郭領域でないときに，境界処理をスキップ
    if( nID[face] >= 0 ) 
    {
      vec[0] =  0.0;   // sum
      vec[1] =  1.0e6; // min
      vec[2] = -1.0e6; // max
      continue;
    }
    
    switch (typ) 
    {
      case OBC_WALL:
        extractVel_OBC(face, vec, tm, v00, fcount);
        div_obc_drchlt_(dv, size, &gd, &face, v00, bv, vec, &fcount);
        break;
    }
  }
  
  flop += fcount;
}


// #################################################################
// 速度境界条件による流束の修正
void SetBC3D::mod_Pvec_Flux(REAL_TYPE* wv, REAL_TYPE* v, int* bv, REAL_TYPE tm, int v_mode, REAL_TYPE* v00, double& flop)
{
  REAL_TYPE vec[3];
  int st[3], ed[3];
  int typ;
  REAL_TYPE dh = deltaX;
  int gd = guide;
  
  // 流束形式の外部境界条件
  for (int face=0; face<NOFACE; face++) {
    typ = obc[face].get_Class();
    
    // 計算領域の最外郭領域でないときに，境界処理をスキップ，次のface面を評価
    if( nID[face] >= 0 ) continue;
    
    switch ( typ )
    {
      case OBC_WALL:
        extractVel_OBC(face, vec, tm, v00, flop);
        pvec_vobc_wall_(wv, size, &gd, &dh, v00, &rei, v, bv, vec, &face, &flop);
        break;
    }
  }
  
}



// #################################################################
// 速度境界条件によるPoisosn式のソース項の修正
void SetBC3D::mod_Psrc_VBC(REAL_TYPE* s_0, REAL_TYPE* vc, REAL_TYPE* v0, int* bv, REAL_TYPE tm, REAL_TYPE dt, REAL_TYPE* v00, double &flop)
{
  int st[3], ed[3];
  REAL_TYPE vec[3], vel;
  REAL_TYPE dh = deltaX;
  int typ;
  int gd = guide;
  double fcount = 0.0;
  
  
  // 外部境界条件による修正
  for (int face=0; face<NOFACE; face++) {
    typ = obc[face].get_Class();
    
    // 計算領域の最外郭領域でないときに，境界処理をスキップ，次のface面を評価
    if( nID[face] >= 0 ) continue;
    
    switch ( typ )
    {
      case OBC_WALL:
      {
        extractVel_OBC(face, vec, tm, v00, fcount);
        div_obc_drchlt_(s_0, size, &gd, &face, v00, bv, vec, &fcount);
        break;
      }
    }
  }
  
  flop += fcount;
}


/**
 @brief 圧力の外部境界条件
 @param d_p 圧力のデータクラス
 */
void SetBC3D::OuterPBC(REAL_TYPE* d_p)
{
  int uod, F;
  REAL_TYPE pv=0.0;
  int gd = guide;
  
  for (int face=0; face<NOFACE; face++) {
    F = obc[face].get_Class();
    
    // 周期境界条件
    if ( F == OBC_PERIODIC )
    {
    }
    else // 周期境界条件以外の処理
    {
    }
  }
}


// #################################################################
/**
 @brief 速度の外部境界条件処理(タイムステップに一度)
 @param[out] v 速度ベクトル v^{n+1}
 @param vc 速度ベクトル v^*
 @param bv BCindex V
 @param tm
 @param dt 
 @param v00
 @param flop
 */
void SetBC3D::OuterVBC(REAL_TYPE* d_v, REAL_TYPE* d_vc, int* d_bv, REAL_TYPE tm, REAL_TYPE dt, REAL_TYPE* v00, double& flop)
{
  REAL_TYPE vec[3];
  int gd = guide;
  
  for (int face=0; face<NOFACE; face++) {

    if( nID[face] >= 0 ) continue;
    
    switch ( obc[face].get_Class() )
    {
    }
    
  }  
}


// #################################################################
/**
 @brief 疑似速度の外部境界条件処理
 @param[out] vc 疑似速度ベクトル v^*
 @param v0 速度ベクトル v^n
 @param bv BCindex V
 @param tm
 @param dt 
 @param v00
 @param flop
 */
void SetBC3D::OuterVBC_Pseudo(REAL_TYPE* d_vc, REAL_TYPE* d_v0, int* d_bv, REAL_TYPE tm, REAL_TYPE dt, REAL_TYPE* v00, double& flop)
{
  REAL_TYPE v_cnv;
  REAL_TYPE dh = deltaX;
  int gd = guide;
  
  for (int face=0; face<NOFACE; face++) {
    
    if( nID[face] >= 0 ) continue; 
    
    switch ( obc[face].get_Class() )
    {
    }
    
  }  
}

