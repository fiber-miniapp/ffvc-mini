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
 * @file   ffv_Post.C
 * @brief  FFV Class
 * @author kero
 */

#include "ffv.h"
#include "timing.h"


// 終了時の処理
bool FFV::Post() 
{
#ifdef PROF_MAPROF
  if (IsMaster()) printf("\n>> Timings\n\n");
  TIME_PRINT(tm_TOTAL,          "Total                     ");
  TIME_PRINT(tm_INIT,           "  Initial                 ");
  TIME_PRINT(tm_MAIN,           "  Main Loop               ");
  TIME_PRINT(tm_SETUP,          "    Setup                 ");
  TIME_PRINT(tm_NS_FS_E_BINARY, "    NS_FS_E_Binary        ");
  TIME_PRINT(tm_PVEC_MUSCL,     "      pvec_muscl_         ");
  TIME_PRINT(tm_VP_ITR,         "      V-P Iteration       ");
  TIME_PRINT(tm_UPDATE_VEC,     "        update_vec_       ");
  TIME_PRINT(tm_SOR_2_SMA,      "        SOR_2_SMA         ");
  TIME_PRINT(tm_PSOR2SMA_CORE,  "          psor2sma_core_  ");
  TIME_PRINT(tm_POI_COMM,       "          Poisson Comm    ");
  TIME_PRINT(tm_POI_RESIDUAL,   "          poi_residual_   ");
  TIME_PRINT(tm_POI_SRC_COMM,   "          Poisson Src Comm");
  TIME_PRINT(tm_MONITOR,        "    Monitor               ");
  TIME_PRINT(tm_OUTPUT,         "    Output                ");

  if (IsMaster()) printf("\n>> FLOPS/node\n\n");
  FLOPS_PRINT(tm_PVEC_MUSCL,     "pvec_muscl_   ");
  FLOPS_PRINT(tm_UPDATE_VEC,     "update_vec_   ");
  FLOPS_PRINT(tm_PSOR2SMA_CORE,  "psor2sma_core_");
  FLOPS_PRINT(tm_POI_RESIDUAL,   "poi_residual_ ");

  maprof_setup("FFVC MINI", FFVC_MINI_VERSION);
  maprof_app_add_int("max_vp_itr", IC[ItrCtl::ic_div].get_ItrMax());
  maprof_app_add_int("max_poisson_itr", IC[ItrCtl::ic_prs_pr].get_ItrMax());
  maprof_profile_add_problem_size("size_x", size[0]);
  maprof_profile_add_problem_size("size_y", size[1]);
  maprof_profile_add_problem_size("size_z", size[2]);

  maprof_add_section("main_loop", tm_MAIN);
  maprof_add_section("pvec_muscl", tm_PVEC_MUSCL);
  maprof_add_section("update_vec", tm_UPDATE_VEC);
  maprof_add_section("psor2sma_core", tm_PSOR2SMA_CORE);
  maprof_add_section("poi_residual", tm_POI_RESIDUAL);
  maprof_add_section("poisson_comm", tm_POI_COMM);

  maprof_output();
#endif

  return true;
}

