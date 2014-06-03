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
 * @file   History.C
 * @brief  FlowBase History class
 * @author kero
 */

#include "History.h"


// #################################################################
// タイムスタンプの更新
void History::updateTimeStamp(const int m_stp, const REAL_TYPE m_tm, const REAL_TYPE vMax)
{
  step  = m_stp;
  time  = m_tm;
  v_max = vMax;
}


// #################################################################
// 標準履歴モニタのヘッダー出力
void History::printHistoryTitle(FILE* fp, const ItrCtl* IC, const bool disp)
{
  const ItrCtl* ICp1 = &IC[ItrCtl::ic_prs_pr];  /// 圧力のPoisson反復
  const ItrCtl* ICv  = &IC[ItrCtl::ic_vis_cn];  /// 粘性項のCrank-Nicolson反復
  
  fprintf(fp, "    step        time[-]    v_max[-] ItrVP v_div_max[-]");
  
  fprintf(fp, "  ItrP");
  if      (ICp1->get_normType() == ItrCtl::dx_b)       fprintf(fp, "        dx_b");
  else if (ICp1->get_normType() == ItrCtl::r_b)        fprintf(fp, "         r_b");
  else if (ICp1->get_normType() == ItrCtl::r_r0)       fprintf(fp, "        r_r0");
  
  fprintf(fp, "     deltaP       avrP     deltaV       avrV");
  
  if ( disp )
  {
    fprintf(fp, "     time[sec]");
  }
  
  fprintf(fp, "\n");
}


// #################################################################
// 標準履歴の出力
void History::printHistory(FILE* fp, const double* avr, const double* rms, const ItrCtl* IC, const double stptm, const bool disp)
{
  const ItrCtl* ICp1 = &IC[ItrCtl::ic_prs_pr];  ///< 圧力のPoisson反復
  const ItrCtl* ICd  = &IC[ItrCtl::ic_div];     ///< 圧力-速度反復
  
  fprintf(fp, "%8d %14.6e %11.4e %5d  %11.4e",
          step, time, v_max, ICd->LoopCount, ICd->get_normValue() );

  if ( (IC->get_normType() != ItrCtl::v_div_max) && (IC->get_normType() != ItrCtl::v_div_dbg) )
  {
    fprintf(fp, " %5d %11.4e", ICp1->LoopCount, ICp1->get_normValue());
  }
  else
  {
    fprintf(fp, " %5d %11.4e", ICp1->LoopCount, ICp1->get_normValue());
  }
  

  fprintf(fp, " %10.3e %10.3e %10.3e %10.3e", rms[var_Pressure], avr[var_Pressure], rms[var_Velocity], avr[var_Velocity]);
  
  if ( disp )
  {
    fprintf(fp, "%14.6e", stptm);
  }
  fprintf(fp, "\n");
  fflush(fp);
}
