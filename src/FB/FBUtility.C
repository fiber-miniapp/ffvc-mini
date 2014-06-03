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
 * @file   FBUtility.C
 * @brief  FlowBase FBUtility class
 * @author kero
 */


#include "FBUtility.h"

// #################################################################
// メモリ使用量を表示する
void FBUtility::MemoryRequirement(const char* mode, const double Memory, const double l_memory, FILE* fp)
{
  const double mem = Memory;
  const double lmem= l_memory;
  const double KB = 1024.0;
  const double MB = 1024.0*KB;
  const double GB = 1024.0*MB;
  const double TB = 1024.0*GB;
  const double PB = 1024.0*TB;
  const double factor = 1.05; // estimate 5% for addtional

  fprintf (fp,">> Memory required for %s : ", mode);

  // Global memory
  fprintf (fp," Global=");
  if ( mem > PB ) {
    fprintf (fp,"%6.2f (PB)", mem / PB *factor);
  }
  else if ( mem > TB ) {
    fprintf (fp,"%6.2f (TB)", mem / TB *factor);
  }
  else if ( mem > GB ) {
    fprintf (fp,"%6.2f (GB)", mem / GB *factor);
  }
  else if ( mem > MB ) {
    fprintf (fp,"%6.2f (MB)", mem / MB *factor);
  }
  else if ( mem > KB ) {
    fprintf (fp,"%6.2f (KB)", mem / KB *factor);
  }
  else if ( mem <= KB ){
    fprintf (fp,"%6.2f (B)", mem *factor);
  }
  else {
    fprintf (fp,"Caution! Memory required : %d (Byte)", (int)(mem *factor) );
  }
  
  // Local memory
  fprintf (fp," : Local=");
  if ( lmem > PB ) {
    fprintf (fp,"%6.2f (PB)\n", lmem / PB *factor);
  }
  else if ( lmem > TB ) {
    fprintf (fp,"%6.2f (TB)\n", lmem / TB *factor);
  }
  else if ( lmem > GB ) {
    fprintf (fp,"%6.2f (GB)\n", lmem / GB *factor);
  }
  else if ( lmem > MB ) {
    fprintf (fp,"%6.2f (MB)\n", lmem / MB *factor);
  }
  else if ( lmem > KB ) {
    fprintf (fp,"%6.2f (KB)\n", lmem / KB *factor);
  }
  else if ( lmem <= KB ){
    fprintf (fp,"%6.2f (B)\n", lmem *factor);
  }
  else {
    fprintf (fp,"Caution! Memory required : %d (Byte)\n", (int)(lmem *factor) );
  }

  fflush(fp);
}


// #################################################################
// スカラー倍してコピー
void FBUtility::xcopy(REAL_TYPE* dst, const int* size, const int guide, const REAL_TYPE* src, const REAL_TYPE scale, const int mode, double& flop)
{
  REAL_TYPE s = scale;
  size_t n = (size_t)(size[0]+2*guide) * (size_t)(size[1]+2*guide) * (size_t)(size[2]+2*guide);
  
  if ( mode == kind_vector ) n *= 3;
  
  flop += (double)n;
  
#pragma omp parallel for firstprivate(n, s) schedule(static)
  for (size_t i=0; i<n; i++) {
    dst[i] = s * src[i];
  }
}



// #################################################################
// 初期化
void FBUtility::xset(REAL_TYPE* dst, const int* size, const int guide, const REAL_TYPE init, const int mode)
{
  REAL_TYPE s = init;
  size_t n = (size_t)(size[0]+2*guide) * (size_t)(size[1]+2*guide) * (size_t)(size[2]+2*guide);
  
  if ( mode == kind_vector ) n *= 3;
  
#pragma omp parallel for firstprivate(n, s) schedule(static)
  for (size_t i=0; i<n; i++) {
    dst[i] = s;
  }
}
