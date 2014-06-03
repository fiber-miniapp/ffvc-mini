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
 * @file   ffv_LS.C
 * @brief  FFV Class
 * @author kero
 */

#include "ffv.h"
#include "timing.h"

// #################################################################
// SOR2SMA
int FFV::SOR_2_SMA(ItrCtl* IC, REAL_TYPE* x, REAL_TYPE* b, const double rhs_nrm, const double r0)
{
  int ip;                        /// ローカルノードの基点(1,1,1)のカラーを示すインデクス
                                 /// ip=0 > R, ip=1 > B
  double flop_count=0.0;         /// 浮動小数点演算数
  REAL_TYPE omg = IC->get_omg(); /// 加速係数
	double res = 0.0;              /// 残差
  int lc=0;                      /// ループカウント
  
  // x     圧力 p^{n+1}
  // b     RHS vector
	// d_bcp ビットフラグ
  
  
  for (lc=0; lc<IC->get_ItrMax(); lc++)
  {
    // 2色のマルチカラー(Red&Black)のセットアップ
    
    MPI_Request req[12]; /// 送信ID
    
    for (int i=0; i<12; i++) req[i] = MPI_REQUEST_NULL;
    
    // ip = 0 基点(1,1,1)がRからスタート
    //    = 1 基点(1,1,1)がBからスタート
    if ( numProc > 1 )
    {
      ip = (head[0]+head[1]+head[2]+1) % 2;
    }
    else
    {
      ip = 0;
    }
    
    
    // 各カラー毎の間に同期, 残差は色間で積算する
    // R - color=0 / B - color=1
    for (int color=0; color<2; color++) {
      
      TIME_START(tm_PSOR2SMA_CORE);
      flop_count = 0.0; // 色間で積算しない
      psor2sma_core_(x, size, &guide, &ip, &color, &omg, &res, b, d_bcp, &flop_count);
      TIME_STOP(tm_PSOR2SMA_CORE);
      ADD_FLOPS(tm_PSOR2SMA_CORE, flop_count);
      
      // 境界条件
      BC.OuterPBC(x);
      
      
      // 同期処理
      TIME_START(tm_POI_COMM);
      if ( numProc > 1 )
      {
        /// 通信面1面あたりの通信量
        double comm_size = count_comm_size(size, guide);
        
        if (IC->get_SyncMode() == comm_sync )
        {
          if ( paraMngr->BndCommS3D(x, size[0], size[1], size[2], guide, 1) != CPM_SUCCESS ) Exit(0); // 1 layer communication
        }
        else
        {
          int ireq[12];
          sma_comm_     (x, size, &guide, &color, &ip, cf_sz, cf_x, cf_y, cf_z, ireq, nID);
          sma_comm_wait_(x, size, &guide, &color, &ip, cf_sz, cf_x, cf_y, cf_z, ireq);
        }
      }
      TIME_STOP(tm_POI_COMM);
    }
    
    
    // Residual
    switch ( IC->get_normType() )
    {
      case ItrCtl::r_b:
      case ItrCtl::r_r0:
        
        TIME_START(tm_POI_RESIDUAL);
        res = 0.0;
        flop_count = 0.0;
        poi_residual_(&res, size, &guide, x, b, d_bcp, &flop_count);
        TIME_STOP(tm_POI_RESIDUAL);
        ADD_FLOPS(tm_POI_RESIDUAL, flop_count);
        break;
        
      case ItrCtl::dx_b:
        // nothing to do, dx is already obtained in psor_(&res,...)
        break;
    }
    
    TIME_START(tm_POI_SRC_COMM);
    if ( numProc > 1 )
    {
      double m_tmp = res;
      if ( paraMngr->Allreduce(&m_tmp, &res, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
    }
    TIME_STOP(tm_POI_SRC_COMM);
    
    res = sqrt(res);
    
    
    // 残差の保存
    switch ( IC->get_normType() )
    {
      case ItrCtl::dx_b:
        IC->set_normValue( res/rhs_nrm );
        break;
        
      case ItrCtl::r_b:
        IC->set_normValue( res/rhs_nrm );
        break;
        
      case ItrCtl::r_r0:
        IC->set_normValue( res/r0 );
        break;
        
      default:
        printf("\tInvalid Linear Solver for Pressure\n");
        Exit(0);
        break;
    }
    
    // 収束判定　性能測定モードのときは収束判定を行わない
    if ( (PM_Test == OFF) && (IC->get_normValue() < IC->get_eps()) ) break;
    
  }
  
  return lc;
}

