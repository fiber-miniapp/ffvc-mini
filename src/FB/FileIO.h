#ifndef _FB_FILE_IO_H_
#define _FB_FILE_IO_H_

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
 * @file   FileIO.h
 * @brief  FlowBase FileIO class Header
 * @author kero
 */

#include <string>

#include "cpm_ParaManager.h"

#include "DomainInfo.h"
#include "FB_Define.h"
#include "FB_Ffunc.h"

class FileIO : public DomainInfo {
  
public:
  /** コンストラクタ */
  FileIO() {}
  
  /**　デストラクタ */
  ~FileIO() {}
  

  /** 
   * @brief スカラーファイルを出力する
   * @param [in] fname     ファイル名
   * @param [in] sz        分割数
   * @param [in] gc        ガイドセル数
   * @param [in] s         スカラー場
   * @param [in] step      ステップ
   * @param [in] time      時刻
   * @param [in] org       領域の基点
   * @param [in] pit       セル幅
   * @param [in] guide_out ガイドセル数
   * @param [in] mode      平均値出力指示（瞬時値のときtrue，平均値のときfalse）
   * @param [in] step_avr  平均操作したステップ数
   * @param [in] time_avr  平均操作した時間
   */
  void writeScalar(const std::string fname, 
                   int* sz, 
                   int gc,
                   REAL_TYPE* s, 
                   const unsigned step, 
                   const double time, 
                   const REAL_TYPE* org, 
                   const REAL_TYPE* pit, 
                   const int guide_out,
                   const bool mode=true,
                   const unsigned step_avr=0,
                   const double time_avr=0.0);
  
  
  /** 
   * @brief ベクトルファイルを出力する
   * @param [in] fname     ファイル名
   * @param [in] sz        分割数
   * @param [in] gc        ガイドセル数
   * @param [in] v         ベクトル場
   * @param [in] step      ステップ
   * @param [in] time      時刻
   * @param [in] org       領域の基点
   * @param [in] pit       セル幅
   * @param [in] guide_out ガイドセル数
   * @param [in] mode      平均値出力指示（瞬時値のときtrue，平均値のときfalse）
   * @param [in] step_avr  平均操作したステップ数
   * @param [in] time_avr  平均操作した時間
   */
  void writeVector(const std::string fname, 
                   int* sz, 
                   int gc, 
                   REAL_TYPE* v, 
                   const unsigned step, 
                   const double time, 
                   const REAL_TYPE* org, 
                   const REAL_TYPE* pit, 
                   const int guide_out,
                   const bool mode=true,
                   const unsigned step_avr=0,
                   const double time_avr=0.0);
  
};
#endif // _FB_FILE_IO_H_
