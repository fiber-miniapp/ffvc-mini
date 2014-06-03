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
 * @file   main.C
 * @brief  ffvのmain関数
 * @author kero
 */

#include "ffv.h"
#include "timing.h"


// return; 0 - normal
//         1 - others
int main( int argc, char **argv )
{
  // FFV classのインスタンス
  FFV ffv;

  TIME_START(tm_TOTAL);

  // ##################################################################
  // 初期化
  TIME_START(tm_INIT);
  
  // 並列管理クラスのインスタンスと初期化
  // ここでMPI_Initも行う
  if ( !ffv.importCPM(cpm_ParaManager::get_instance(argc, argv)) )
  {
    return -1;
  }
  
  int init_ret = ffv.Initialize(argc, argv);
  
  switch( init_ret )
  {
    case -1:
      if ( ffv.IsMaster() ) printf("\n\tSolver initialize error.\n\n");
      return -1;

    case 0:
      if ( ffv.IsMaster() ) printf("\n\tForced termination during initialization.\n\n");
      return -1;
      
    case 1:
      // keep going on processing
      break;
      
    default:
      if ( ffv.IsMaster() ) printf("\n\tSolver initialize error.\n\n");
      return -1;
  }
  
  TIME_STOP(tm_INIT);
  
  // ##################################################################
  // タイムステップループ

  TIME_START(tm_MAIN);
  int loop_ret = ffv.MainLoop();
  TIME_STOP(tm_MAIN);
  
  switch (loop_ret) 
  {
    case -1:
      if ( ffv.IsMaster() ) printf("\n\tSolver error.\n\n");
      break;

    case 0:
      if ( ffv.IsMaster() ) printf("\n\tSolver forced termination time-step loop.\n\n");
      break;
      
    case 1:
      if ( ffv.IsMaster() ) printf("\n\tSolver finished.\n\n");
      break;
  }
  
  TIME_STOP(tm_TOTAL);

  // ##################################################################
  // ポスト処理
  if( !ffv.Post() )
  {
    printf("solver post error.\n");
    return -1;
  }

  if ( loop_ret != 1 ) return -1;
  
  return 0;
}
