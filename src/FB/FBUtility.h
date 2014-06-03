#ifndef _FB_UTY_H_
#define _FB_UTY_H_

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
 * @file   FBUtility.h
 * @brief  FlowBase FBUtility class Header
 * @author kero
 */

#include <math.h>
#include <string>
#include <stdio.h>

#include "FB_Define.h"
#include "cpm_Define.h"

// #################################################################
class FBUtility {

public:
  /** コンストラクタ */
  FBUtility() {}
  
  /** デストラクタ */
  ~FBUtility() {}
 
  
public:
  
  /**
   * @brief dirの方向ラベルを返す
   * @param [in] dir 方向コード
   * @return 方向ラベル
   */
  static std::string getDirection(const int dir)
  {
    std::string face;
    if      (dir == X_MINUS) face = "X-";
    else if (dir == X_PLUS)  face = "X+";
    else if (dir == Y_MINUS) face = "Y-";
    else if (dir == Y_PLUS)  face = "Y+";
    else if (dir == Z_MINUS) face = "Z-";
    else if (dir == Z_PLUS)  face = "Z+";
    return face;
  }
  
  
  /**
   * @brief メモリ使用量を表示する
   * @param [in] mode     処理モード
   * @param [in] Memory   必要メモリ量
   * @param [in] l_memory local
   * @param [in] fp       ファイルポインタ
   */
  static void MemoryRequirement(const char* mode, const double Memory, const double l_memory, FILE* fp);
  
  
  /**
   * @brief スカラー倍コピー
   * @param [out]    dst   出力
   * @param [in]     size  配列サイズ
   * @param [in]     guide ガイドセルサイズ
   * @param [in]     src   入力
   * @param [in]     scale スカラー倍数
   * @param [in]     mode  スカラー or ベクトル
   * @param [in,out] flop  浮動小数点演算
   */
  static void xcopy (REAL_TYPE* dst,
                     const int* size,
                     const int guide,
                     const REAL_TYPE* src,
                     const REAL_TYPE scale,
                     const int mode,
                     double& flop);

  
  /**
   * @brief 初期化
   * @param [out]    dst   出力
   * @param [in]     size  配列サイズ
   * @param [in]     guide ガイドセルサイズ
   * @param [in]     init  定数
   * @param [in]     mode  スカラー or ベクトル
   */
  static void xset (REAL_TYPE* dst,
                    const int* size,
                    const int guide,
                    const REAL_TYPE init,
                    const int mode);
  
};

#endif // _FB_UTY_H_
