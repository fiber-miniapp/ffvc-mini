#ifndef _FB_ITR_CTL_H_
#define _FB_ITR_CTL_H_

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
 @file   ItrCtl.h
 @brief  FlowBase ItrCtl class Header
 @author kero
 */

#include <string>
#include "FB_Define.h"

class ItrCtl {
private:
  double NormValue;      ///< ノルムの値
  double eps;            ///< 収束閾値
  REAL_TYPE omg;         ///< 加速/緩和係数
  int NormType;          ///< ノルムの種類
  int SubType;           ///< SKIP LOOP or MASK LOOP
  int ItrMax;            ///< 最大反復数
  int LinearSolver;      ///< 線形ソルバーの種類
  int SyncMode;          ///< 同期モード (comm_sync, comm_async)
  
public:
  int LoopCount;  ///< 反復回数
  
  /** 反復制御リスト */
  enum itr_cntl_key 
  {
    ic_prs_pr,
    ic_prs_cr,
    ic_vis_cn,
    ic_tdf_ei,
    ic_tdf_cn,
    ic_div,
    ic_END
  };
  
  
  /** 反復法の収束基準種別 */
  enum norm_type 
  { 
    v_div_max=1,
    v_div_dbg,
    dx_b,
    r_b,
    r_r0
  };

  
  /** コンストラクタ */
  ItrCtl() {
    NormType = 0;
    ItrMax = LoopCount = LinearSolver = SubType = 0;
    SyncMode = 0;
    eps = omg = 0.0;
    NormValue = 0.0;
  }
  
  /**　デストラクタ */
  ~ItrCtl() {}
  
  
public:
  /** @brief 線形ソルバの種類を返す */
  int get_LS() const 
  { 
    return LinearSolver; 
  }
  
  /** @brief 最大反復回数を返す */
  int get_ItrMax() const 
  { 
    return ItrMax; 
  }
  
  /** @brief ループ実装の種類を返す */
  int get_LoopType() const 
  { 
    return SubType; 
  }
  
  /** @brief 緩和/加速係数を返す */
  REAL_TYPE get_omg() const 
  { 
    return omg; 
  }
  
  /** @brief 収束閾値を返す */
  double get_eps() const
  { 
    return eps; 
  }
  
  /** @brief ノルムのタイプを返す */
  int get_normType() const 
  { 
    return NormType; 
  }
  
  /** @brief keyに対応するノルムの値を返す */
  double get_normValue() const
  { 
    return NormValue; 
  }
  
  /** @brief 同期モードを返す */
  int get_SyncMode() const 
  { 
    return SyncMode; 
  }
  
  
  /** @brief 線形ソルバの種類を設定する */
  void set_LS(const int key) 
  { 
    LinearSolver=key; 
  }
  
  /** @brief 最大反復回数を設定する */
  void set_ItrMax(const int key) 
  { 
    ItrMax=key; 
  }
  
  /** @brief ループ実装の種類を設定する */
  void set_LoopType(const int key) 
  { 
    SubType=key; 
  }

  /** @brief 緩和/加速係数を保持 */
  void set_omg(const REAL_TYPE r) 
  { 
    omg = r; 
  }
  
  /** @brief 収束閾値を保持 */
  void set_eps(const double r)
  { 
    eps = r; 
  }
  
  /** @brief ノルムのタイプを保持 */
  void set_normType(const int n) 
  { 
    NormType = n; 
  }
  
  /** @brief ノルム値を保持 */
  void set_normValue(const double r)
  { 
    NormValue = r; 
  }
  
  /** @brief 同期モードを保持 */
  void set_SyncMode(const Synch_Mode r) 
  { 
    SyncMode= r; 
  }

  /** @brief ノルムのラベルを返す */
  static std::string getNormString(const int d)
  {
    switch (d)
    {
      case v_div_dbg:
        return "Max. Norm : Divergence of velocity with Monitoring  ### Forced to be selected since Iteration Log is specified ###";
      case v_div_max:
        return "Max. Norm : Divergence of velocity";
      case dx_b:
        return "dx_b : Increment of vector x divided by RHS vector b";
      case r_b:
        return "r_b  : Residual vector divided by RHS vector b";
      case r_r0:
        return "r_r0 : Residual vector divided by initial residual vector";
      default:
        Exit(1);
    }
    /* NOTREACHED */
  }

};


#endif // _FB_ITR_CTL_H_
