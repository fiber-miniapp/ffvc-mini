#ifndef _FB_HISTORY_H_
#define _FB_HISTORY_H_

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
 * @file   History.h
 * @brief  FlowBase History class Header
 * @author kero
 */

#include "ItrCtl.h"
#include "Component.h"
#include "FBUtility.h"

class History {
protected:
  REAL_TYPE time;            ///< 時刻
  REAL_TYPE v_max;           ///< 速度成分の最大値
  int step;                  ///< ステップ数
  
public:
  
  /** コンストラクタ */
  History()
  {
    time  = 0.0;
    v_max = 0.0;
    step  = 0;
  }
  
  /**　デストラクタ */
  ~History() {}
  
public:
  
  /**
   * @brief 標準履歴の出力
   * @param [in] fp    出力ファイルポインタ
   * @param [in] avr   1タイムステップの平均値　（0-pressure, 1-velocity, 2-temperature)
   * @param [in] rms   1タイムステップの変化量　（0-pressure, 1-velocity, 2-temperature)
   * @param [in] IC    ItrCtlクラスのポインタ
   * @param [in] stptm 1タイムステップの計算時間
   * @param [in] disp  計算時間表示の有無
   */
  void printHistory(FILE* fp, const double* avr, const double* rms, const ItrCtl* IC, const double stptm, const bool disp);
  
  /**
   * @brief 標準履歴の出力のヘッダー出力
   * @param [in] fp 出力ファイルポインタ
   * @param [in] IC 反復管理クラス
   * @param [in] disp  計算時間表示の有無
   */
  void printHistoryTitle(FILE* fp, const ItrCtl* IC, const bool disp);
  
  /**
   * @brief タイムスタンプの更新
   * @param [in] m_stp ステップ数
   * @param [in] m_tm  時刻
   * @param [in] vMax  速度最大値成分
   */
  void updateTimeStamp(const int m_stp, const REAL_TYPE m_tm, const REAL_TYPE vMax);
};

#endif // _FB_HISTORY_H_
