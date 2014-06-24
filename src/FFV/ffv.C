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
 * @file   ffv.C
 * @brief  FFV Class
 * @author kero
 */

#include "ffv.h"

// コンストラクタ
FFV::FFV()
{
  ffv_procGrp = 0;
  G_Acell = 0;
  CurrentTime = 0.0;
  Session_StartTime = 0.0;
  Session_CurrentTime = 0.0;
  
  Session_LastStep = 0;
  Session_CurrentStep = 0;
  Session_StartStep = 0;
  CurrentStep = 0;

  convergence_prev = 0;
  
  for (int i=0; i<3; i++) 
  {
    G_size[i]= 0;
    G_origin[i] = 0.0;
    G_region[i] = 0.0;
  }
  
  fp_b = NULL;
  mat = NULL;
  cmp = NULL;
  paraMngr = NULL;

  // dfi管理
  for (int i=0; i<var_END; i++) dfi_mng[i]=0;
  
}


// デストラクタ
FFV::~FFV()
{
  
}


// #################################################################
// タイムステップループ
int FFV::MainLoop()
{
  int ret = 1;
  
  for (int i=1; i<=Session_LastStep; i++)
  {
    Session_CurrentStep = i;
    
    int loop_ret = Loop(i);
    
    switch (loop_ret) 
    {
      case -1: // error
        return -1;
        break;
        
      case 0: // forced terminated
        ret = 1;
        break;
        
      case 1: // normal
        ret = 1;
        break;
        
      default:
        return -1;
    }
    
    if( loop_ret == 0 ) break;
  }

  return ret;
}


// #################################################################
// V-P反復のdiv(u)ノルムを計算する
void FFV::Norm_Div(ItrCtl* IC)
{
  double nrm;
  double flop_count;
  REAL_TYPE coef = 1.0/deltaX; /// 発散値を計算するための係数

  
  switch (IC->get_normType())
  {

    case ItrCtl::v_div_max:
      flop_count=0.0;
      norm_v_div_max_(&nrm, size, &guide, d_dv, &coef, d_bcp, &flop_count);
      
      if ( numProc > 1 )
      {
        double tmp;
        tmp = nrm;
        if ( paraMngr->Allreduce(&tmp, &nrm, 1, MPI_MAX) != CPM_SUCCESS ) Exit(0); // 最大値
      }
      IC->set_normValue(nrm);
      break;
      
    default:
      stamped_printf("\tInvalid convergence type\n");
      Exit(0);
  }
}


// #################################################################
// 空間平均操作と変動量の計算を行う
// スカラ値は算術平均，ベクトル値は自乗和
void FFV::Variation_Space(double* avr, double* rms, double& flop)
{
  double m_var[2];
  
  // 速度
  fb_delta_v_(m_var, size, &guide, d_v, d_v0, d_bcd, &flop); // 速度反復でV_res_L2_を計算している場合はスキップすること
  rms[var_Velocity] = m_var[0];
  avr[var_Velocity] = m_var[1];
  
  // 圧力
  fb_delta_s_(m_var, size, &guide, d_p, d_p0, d_bcd, &flop);
  rms[var_Pressure] = m_var[0];
  avr[var_Pressure] = m_var[1];
  
}


// #################################################################
// 参照速度を設定する
void FFV::setV00(double m_time) 
{
  if ( TimeAccel == 0.0 )
  {
    v00[0] = 1.0;
  }
  else
  {
    const double c_pai = (double)(2.0*asin(1.0));
    v00[0] = 0.5*(1.0-cos(c_pai*m_time/(TimeAccel)));
    if ( m_time > (TimeAccel) ) v00[0] = 1.0;
  }

  v00[1] = 0.0;
  v00[2] = 0.0;
  v00[3] = 0.0;
}


// #################################################################
// ファイル出力
void FFV::FileOutput()
{

  // ファイルプレフィックス
  const char* f_Pressure = "prs_";
  const char* f_Velocity = "vel_";
  
  // ステップ数
  int step = (int)CurrentStep;
  
  // 時間
  REAL_TYPE time = (REAL_TYPE)CurrentTime;
  
  // ガイドセル出力
  int gc_out = 0;
  
  // 出力ファイル名
  std::string fname;
  
  // 並列出力モード
  bool pout = true;
  
  // Pressure
  fname = "./" + DFI_.Generate_FileName(f_Pressure, step, myRank, pout);
  F.writeScalar(fname, size, guide, d_p, step, time, origin, pitch, gc_out);
  Hostonly_ if ( !DFI_.Write_DFI_File(f_Pressure, step, dfi_mng[var_Pressure], pout) ) Exit(0);
  
  
  // Velocity
  fname = "./" + DFI_.Generate_FileName(f_Velocity, step, myRank, pout);
  F.writeVector(fname, size, guide, d_v, step, time, origin, pitch, gc_out);
  Hostonly_ if ( !DFI_.Write_DFI_File(f_Velocity, step, dfi_mng[var_Velocity], pout) ) Exit(0);

}
