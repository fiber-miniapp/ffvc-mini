// #################################################################
//
// FFV : Frontflow / violet
//
// Copyright (c) 2012-2014 Institute of Industrial Science, The University of Tokyo.
// All right reserved.
//
// Copyright (c) 2012-2014 Advanced Institute for Computational Science, RIKEN.
// All right reserved.
//
// #################################################################

/** 
 * @file   FileIO.C
 * @brief  FlowBase FileIO class Header
 * @author kero
 */

#include "FileIO.h"

#include <string.h>

// スカラーファイルを出力する
void FileIO::writeScalar(const std::string fname, 
                         int* sz, 
                         int gc,
                         REAL_TYPE* s, 
                         const unsigned step, 
                         const double time, 
                         const REAL_TYPE* org, 
                         const REAL_TYPE* pit, 
                         const int guide_out,
                         const bool mode,
                         const unsigned step_avr,
                         const double time_avr)
{
  if ( fname.empty() ) Exit(0);
  
  if ( fname.size() > FB_FILE_PATH_LENGTH ) {
    Hostonly_ printf ("Error : Length of file path is greater than %d\n", FB_FILE_PATH_LENGTH);
    Exit(0);
  }
  
  
  char tmp[FB_FILE_PATH_LENGTH];
  memset(tmp, 0, sizeof(char)*FB_FILE_PATH_LENGTH);
  strcpy(tmp, fname.c_str());

  int stp = (int)step;
  REAL_TYPE tm = (REAL_TYPE)time;
  int g = guide_out;
  REAL_TYPE o[3], p[3];
  o[0] = org[0];
  o[1] = org[1];
  o[2] = org[2];
  p[0] = pit[0];
  p[1] = pit[1];
  p[2] = pit[2];
  
  int avs = (mode == true) ? 1 : 0;
  int stp_a = (int)step_avr;
  REAL_TYPE tm_a = (REAL_TYPE)time_avr;
  int d_type = (sizeof(REAL_TYPE) == 4) ? 1 : 2;  // 1-float / 2-double
  
  fb_write_sph_s_ (s, sz, &gc, tmp, &stp, &tm, o, p, &d_type, &g, &avs, &stp_a, &tm_a);
  
}


// ベクトルファイルを出力する
void FileIO::writeVector(const std::string fname, 
                         int* sz, 
                         int gc, 
                         REAL_TYPE* v, 
                         const unsigned step, 
                         const double time, 
                         const REAL_TYPE* org, 
                         const REAL_TYPE* pit, 
                         const int guide_out,
                         const bool mode,
                         const unsigned step_avr,
                         const double time_avr)
{
  if ( fname.empty() ) Exit(0);
  
  if ( fname.size() > FB_FILE_PATH_LENGTH ) {
    Hostonly_ printf ("Error : Length of file path is greater than %d\n", FB_FILE_PATH_LENGTH);
    Exit(0);
  }
  
  char tmp[FB_FILE_PATH_LENGTH];
  memset(tmp, 0, sizeof(char)*FB_FILE_PATH_LENGTH);
  strcpy(tmp, fname.c_str());
  
  int stp = (int)step;
  REAL_TYPE tm = (REAL_TYPE)time;
  int g = guide_out;
  REAL_TYPE o[3], p[3];
  o[0] = org[0];
  o[1] = org[1];
  o[2] = org[2];
  p[0] = pit[0];
  p[1] = pit[1];
  p[2] = pit[2];
  
  int avs = (mode == true) ? 1 : 0;
  int stp_a = (int)step_avr;
  REAL_TYPE tm_a = (REAL_TYPE)time_avr;
  int d_type = (sizeof(REAL_TYPE) == 4) ? 1 : 2;  // 1-float / 2-double
  
  fb_write_sph_v_ (v, sz, &gc, tmp, &stp, &tm, o, p, &d_type, &g, &avs, &stp_a, &tm_a);

}

