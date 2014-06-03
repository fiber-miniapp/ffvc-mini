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
 * @file   ffv_Loop.C
 * @brief  FFV Class
 * @author kero
 */

#include "ffv.h"
#include "timing.h"


// タイムステップループの処理
int FFV::Loop(const unsigned step) 
{
  // 1 step elapse
  step_start = cpm_Base::GetWTime();
  double step_end;
  
  double flop_count=0.0;   /// 浮動小数演算数
  double avr_Var[3];       /// 平均値（速度、圧力、温度）
  double rms_Var[3];       /// 変動値
  REAL_TYPE vMax=0.0;      /// 最大速度成分

  TIME_START(tm_SETUP);

  
  // 時間進行
  CurrentTime += deltaT;
  CurrentStep++;
  
  
  // 参照座標速度をv00に保持する
  setV00(CurrentTime);

  
  // 速度成分の最大値
  flop_count = 0.0;
  find_vmax_(&vMax, size, &guide, v00, d_v, &flop_count);
  
  if ( numProc > 1 ) 
  {
    REAL_TYPE vMax_tmp = vMax;
    if ( paraMngr->Allreduce(&vMax_tmp, &vMax, 1, MPI_MAX) != CPM_SUCCESS ) Exit(0);
  }

  TIME_STOP(tm_SETUP);
  
  
  // Flow
  TIME_START(tm_NS_FS_E_BINARY);
  NS_FS_E_Binary();
  TIME_STOP(tm_NS_FS_E_BINARY);


  TIME_START(tm_MONITOR);
  
  // 空間平均値操作と変動量
  flop_count=0.0;
  for (int i=0; i<3; i++) 
  {
    avr_Var[i] = 0.0;
    rms_Var[i] = 0.0;
  }
  Variation_Space(avr_Var, rms_Var, flop_count);

  
  if ( numProc > 1 ) 
  {
    /// var_Velocity=0,  > FB_Define.h
    /// var_Pressure,
    /// var_Temperature,
    double src[6], dst[6]; // Vel, Prs, Tempで3*2
    
    for (int n=0; n<3; n++) {
      src[n]   = avr_Var[n];
      src[n+3] = rms_Var[n];
    }
    
    if ( paraMngr->Allreduce(src, dst, 6, MPI_SUM) != CPM_SUCCESS) Exit(0); // 変数 x (平均値+変動値)
    
    for (int n=0; n<3; n++) {
      avr_Var[n] = dst[n];
      rms_Var[n] = dst[n+3];
    }
  }

  avr_Var[var_Velocity] /= (double)G_Acell;  // 速度の空間平均
  avr_Var[var_Pressure] /= (double)G_Acell;  // 圧力の空間平均
  
  rms_Var[var_Velocity] /= (double)G_Acell;  // 速度の変動量
  rms_Var[var_Pressure] /= (double)G_Acell;  // 圧力の変動量
  rms_Var[var_Velocity] = sqrt(rms_Var[var_Velocity]);
  rms_Var[var_Pressure] = sqrt(rms_Var[var_Pressure]);
  
  
  // Historyクラスのタイムスタンプを更新
  H->updateTimeStamp(CurrentStep, (REAL_TYPE)CurrentTime, vMax);
  

  // 1 step elapse
  step_end = cpm_Base::GetWTime() - step_start;

  

  // 瞬時値のデータ出力
  if ( OutputInterval > 0 && CurrentStep % OutputInterval == 0 )
  {
    TIME_START(tm_OUTPUT);
    FileOutput();
    TIME_STOP(tm_OUTPUT);
  }
  
  
  // 基本履歴情報
  {
    // コンソール出力
    Hostonly_ H->printHistory(stdout, avr_Var, rms_Var, IC, step_end, true);

    // ファイル出力
    Hostonly_ H->printHistory(fp_b, avr_Var, rms_Var, IC, step_end, true);
  }
  
  
  // 発散時の打ち切り
  if ( CurrentStep > 1 ) 
  {
    if ( (convergence_rate > 100.0) )
    {
      Hostonly_ {
        printf      ("\tForced termination : converegence rate >> 100.0\n");
        fprintf(fp_b,"\tForced termination : converegence rate >> 100.0\n");
      }
      return -1;
    }
  }
  
  TIME_STOP(tm_MONITOR);
  
  return 1;
}
