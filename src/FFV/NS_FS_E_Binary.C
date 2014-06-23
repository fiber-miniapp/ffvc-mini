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
 * @file   NS_FS_E_Binary.C
 * @brief  FFV Class
 * @author kero
 */

#include "ffv.h"
#include "timing.h"


// Fractional Step法でNavier-Stokes方程式を解く．バイナリ近似．
void FFV::NS_FS_E_Binary()
{
  // local variables
  double flop;                         /// 浮動小数演算数
  double rhs_nrm = 0.0;                /// 反復解法での定数項ベクトルのL2ノルム
  double res_init = 0.0;               /// 反復解法での初期残差ベクトルのL2ノルム
  double convergence=0.0;              /// 定常収束モニター量
  
  REAL_TYPE dt = deltaT;               /// 時間積分幅
  REAL_TYPE dh = (REAL_TYPE)deltaX;    /// 空間格子幅
  REAL_TYPE coef = deltaX/dt;          /// Poissonソース項の係数
  REAL_TYPE Re = Reynolds;             /// レイノルズ数
  REAL_TYPE rei = 1.0 / Reynolds;      /// レイノルズ数の逆数
  REAL_TYPE half = 0.5;                /// 定数
  REAL_TYPE one = 1.0;                 /// 定数
  REAL_TYPE zero = 0.0;                /// 定数
  int wall_prof = 0;                   /// 壁面条件（noslip）
  int cnv_scheme = CnvScheme;          /// 対流項スキーム
  
  int v_mode=0;
  
  ItrCtl* ICp = &IC[ItrCtl::ic_prs_pr];  /// 圧力のPoisson反復
  ItrCtl* ICv = &IC[ItrCtl::ic_vis_cn];  /// 粘性項のCrank-Nicolson反復
  ItrCtl* ICd = &IC[ItrCtl::ic_div];     /// 圧力-速度反復

  // point Data
  // d_v   セルセンタ速度 v^n -> v^{n+1}
  // d_v0  セルセンタ速度 v^nの保持
  // d_vc  疑似速度ベクトル
  // d_wv  ワーク　陰解法の時の疑似速度ベクトル，射影ステップの境界条件
  // d_p   圧力 p^n -> p^{n+1}
  // d_p0  圧力 p^nの保持
  // d_ws  Poissonのソース項0　速度境界を考慮
  // d_sq  Poissonのソース項1　反復毎に変化するソース項，摩擦速度，
  // d_dv  発散値, div(u)の値を保持
  // d_b   反復の右辺ベクトル
  // d_bcd IDのビットフラグ
  // d_bcp 圧力のビットフラグ
  // d_bcv 速度のビットフラグ
  // d_wo  ワーク　壁関数利用時のWSS，ベクトル出力時のテンポラリ
  // d_t0  温度 t^n 
  // d_ab0 Adams-Bashforth用のワーク

  
  // n stepの値を保持 >> In use (d_v0, d_p0)
  flop = 0.0;
  FBUtility::xcopy(d_p0, size, guide, d_p, one, kind_scalar, flop);
  FBUtility::xcopy(d_v0, size, guide, d_v, one, kind_vector, flop);
  

  // 対流項と粘性項の評価 >> In use (d_vc, d_wv)
  TIME_START(tm_PVEC_MUSCL);
  flop = 0.0;
  v_mode = 1; // No_Slip
  pvec_muscl_(d_vc, size, &guide, &dh, &cnv_scheme, v00, &rei, d_v0, d_bcv, d_bcp, &v_mode, d_ws, &wall_prof, d_bcd, &one, &flop);
  TIME_STOP(tm_PVEC_MUSCL);
  ADD_FLOPS(tm_PVEC_MUSCL, flop);

  flop = 0.0;
  BC.mod_Pvec_Flux(d_vc, d_v0, d_bcv, CurrentTime, v_mode, v00, flop);


  // 時間積分
  flop = 0.0;
  euler_explicit_ (d_vc, size, &guide, &dt, d_v0, d_bcd, &flop);

  
  // 疑似ベクトルの境界条件
  flop = 0.0;
  BC.OuterVBC_Pseudo(d_vc, d_v0, d_bcv, CurrentTime, dt, v00, flop);

  // 疑似ベクトルの同期
  if ( numProc > 1 ) 
  {
    if ( paraMngr->BndCommV3DEx(d_vc, size[0], size[1], size[2], guide, 1) != CPM_SUCCESS ) Exit(0);
  }

  
  // Poissonのソース部分
  
  // vの初期値をvcにしておく
  flop = 0.0;
  FBUtility::xcopy(d_v, size, guide, d_vc, one, kind_vector, flop);

  
  // 非VBC面に対してのみ，セルセンターの値から \sum{u^*} を計算
  flop = 0.0;
  divergence_(d_ws, size, &guide, d_vc, d_bcv, v00, &flop);
  
  
  // Poissonソース項の速度境界条件（VBC）面による修正
  flop = 0.0;
  BC.mod_Psrc_VBC(d_ws, d_vc, d_v0, d_bcv, CurrentTime, dt, v00, flop);
  
  
  // 反復ソースの初期化
  FBUtility::xset(d_sq, size, guide, zero, kind_scalar);
  
  
  // 定数項のL2ノルム　rhs_nrm
  rhs_nrm = 0.0;
  flop = 0.0;
  poi_rhs_(&rhs_nrm, d_b, size, &guide, d_ws, d_sq, d_bcp, &dh, &dt, &flop);
  
  if ( numProc > 1 )
  {
    double m_tmp = rhs_nrm;
    if ( paraMngr->Allreduce(&m_tmp, &rhs_nrm, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
  }
  
  rhs_nrm = sqrt(rhs_nrm);
  
  
  // Initial residual
  if ( ICp->get_normType() == ItrCtl::r_r0 )
  {
    res_init = 0.0;
    flop = 0.0;
    poi_residual_(&res_init, size, &guide, d_p, d_b, d_bcp, &flop);
    
    if ( numProc > 1 )
    {
      double m_tmp = res_init;
      if ( paraMngr->Allreduce(&m_tmp, &res_init, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
    }

    res_init = sqrt(res_init);
  }
  
  
  // 反復回数の積算
  int lc = 0;
  ICp->LoopCount = 0;
  
  TIME_START(tm_VP_ITR);
  for (ICd->LoopCount=0; ICd->LoopCount< ICd->get_ItrMax(); ICd->LoopCount++)
  {
    // 線形ソルバー
    switch (ICp->get_LS())
    {
      case SOR2SMA:
        TIME_START(tm_SOR_2_SMA);
        lc += SOR_2_SMA(ICp, d_p, d_b, rhs_nrm, res_init);
        TIME_STOP(tm_SOR_2_SMA);
        break;

      default:
        printf("\tInvalid Linear Solver for Pressure\n");
        Exit(0);
        break;
    }

    
    // 速度のスカラポテンシャルによる射影と速度の発散の計算 d_dvはdiv(u)のテンポラリ保持に利用
    TIME_START(tm_UPDATE_VEC);
    flop = 0.0;
    update_vec_(d_v, d_dv, size, &guide, &dt, &dh, d_vc, d_p, d_bcp, d_bcv, v00, &flop);
    TIME_STOP(tm_UPDATE_VEC);
    ADD_FLOPS(tm_UPDATE_VEC, flop);

    
    // セルフェイス速度の境界条件による修正
    flop=0.0;
    BC.mod_div(d_dv, d_bcv, CurrentTime, v00, flop);
    
    
    // div(u^{n+1})の計算
    Norm_Div(ICd);
    
   
    // 収束判定　性能測定モードのときは収束判定を行わない
    if ( (PM_Test == OFF) && (ICd->get_normValue() < ICd->get_eps()) ) break;
  } // end of iteration
  TIME_STOP(tm_VP_ITR);

  
  ICp->LoopCount = lc;
  
  
  // 同期
  if ( numProc > 1 )
  {
    if ( paraMngr->BndCommV3DEx(d_v, size[0], size[1], size[2], guide, guide) != CPM_SUCCESS ) Exit(0);
  }
  
  
  // 流出境界のガイドセル値の更新と速度境界条件
  flop=0.0;
  BC.OuterVBC(d_v, d_vc, d_bcv, CurrentTime, dt, v00, flop);
  
  
  // ノルムの増加率が規定値をこえたら，終了
  if (convergence_prev != 0.0 ) 
  {
    convergence_rate = convergence / convergence_prev;
  }
  else 
  {
    convergence_rate = 1.0;
  }
  convergence_prev = convergence;
  
}
