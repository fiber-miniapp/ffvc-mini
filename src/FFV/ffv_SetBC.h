#ifndef _FFV_SETBC_H_
#define _FFV_SETBC_H_

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
 * @file   ffv_SetBC.h
 * @brief  FFV BC Class Header
 * @author kero
 */

#include <math.h>

#include "DomainInfo.h"
#include "BndOuter.h"
#include "Component.h"
#include "ffv_Define.h"
#include "ffv_Ffunc.h"


class SetBC3D : public DomainInfo {
protected:

  REAL_TYPE RefV, RefL, rei;

  BoundaryOuter   obc[NOFACE];

public:
  
  /** コンストラクタ */
  SetBC3D() {
    RefV = RefL = rei = 0.0;
  }
  
  /**　デストラクタ */
  virtual ~SetBC3D() {}
  
protected:
  
  /**
   * @brief コンポーネントから速度境界条件の成分を取り出す
   * @param [in]     n    コンポーネントのインデクス
   * @param [out]    vec  ベクトル成分
   * @param [in]     tm   時刻
   * @param [in]     v00  格子速度
   * @param [in,out] flop 浮動小数点演算数
   */
  REAL_TYPE extractVel_OBC (const int n, REAL_TYPE* vec, const REAL_TYPE tm, const REAL_TYPE* v00, double& flop);
  
  
public:
  /** 
   * @brief クラスに必要な変数のコピー
   */
  void setControlVars(REAL_TYPE Reynolds, REAL_TYPE RefVelocity, REAL_TYPE RefLength)
  {
    rei       = 1.0 / Reynolds;
    RefV      = RefVelocity;
    RefL      = RefLength;
  }
  
  
  /** 外部境界リストのポインタを返す */
  BoundaryOuter* export_OBC()
  { 
    return obc;
  }


  
  /** 引数の外部境界面の外部境界リストのポインタを返す
   * @param [in] face 面番号
   */
  BoundaryOuter* export_OBC(const int face)
  { 
    return &obc[face];
  }
  
  /**
   * @brief 速度境界条件による速度の発散の修正ほか
   * @param [in,out] dv     \sum{u}
   * @param [in]     bv     BCindex V
   * @param [in]     tm     無次元時刻
   * @param [in]     v00    基準速度
   * @param [in]     flop   flop count
   * @note 外部境界面のdiv(u)の修正時に領域境界の流量などのモニタ値を計算し，BoundaryOuterクラスに保持 > 反復後にDomainMonitor()で集約
   */
  void mod_div (REAL_TYPE* dv,
                int* bv,
                REAL_TYPE tm,
                REAL_TYPE* v00,
                double& flop);
  
  
  /**
   @brief 速度境界条件によるPoisosn式のソース項の修正
   @param [out] s_0   \sum{u^*}
   @param [in]  vc    セルセンタ疑似速度
   @param [in]  v0    セルセンタ速度 u^n
   @param [in]  bv    BCindex V
   @param [in]  tm    無次元時刻
   @param [in]  dt    時間積分幅
   @param [in]  v00   基準速度
   @param [out] flop  flop count
   */
  void mod_Psrc_VBC (REAL_TYPE* s_0,
                     REAL_TYPE* vc,
                     REAL_TYPE* v0,
                     int* bv,
                     REAL_TYPE tm,
                     REAL_TYPE dt,
                     REAL_TYPE* v00,
                     double& flop);
  
  
  
  /**
   @brief 速度境界条件による流束の修正
   @param [in,out] wv     疑似速度ベクトル u^*
   @param [in]     v      速度ベクトル u^n
   @param [in]     bv     BCindex V
   @param [in]     tm     無次元時刻
   @param [in]     v_mode 粘性項のモード (0=粘性項を計算しない, 1=粘性項を計算する, 2=壁法則)
   @param [in]     v00    基準速度
   @param [out]    flop   flop count
   */
  void mod_Pvec_Flux (REAL_TYPE* wv,
                      REAL_TYPE* v,
                      int* bv,
                      REAL_TYPE tm,
                      int v_mode,
                      REAL_TYPE* v00,
                      double& flop);
  
  void OuterPBC             (REAL_TYPE* d_p);
  void OuterVBC             (REAL_TYPE* d_v, REAL_TYPE* d_vc, int* d_bv, REAL_TYPE tm, REAL_TYPE dt, REAL_TYPE* v00, double& flop);
  void OuterVBC_Pseudo      (REAL_TYPE* d_vc, REAL_TYPE* v0, int* d_bv, REAL_TYPE tm, REAL_TYPE dt, REAL_TYPE* v00, double& flop);
  void updateOuter          (REAL_TYPE* d_v, REAL_TYPE* vc);
  
  
};

#endif // _FFV_SETBC_H_
