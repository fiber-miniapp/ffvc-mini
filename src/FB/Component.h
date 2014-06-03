#ifndef _FB_COMPO_H_
#define _FB_COMPO_H_

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
 * @file   Component.h
 * @brief  FlowBase CompoList class Header
 * @author kero
 */

#include <string>
#include <string.h>
#include "FB_Define.h"
#include "cpm_Define.h"


class CompoList {
public:
  
  /** 変動速度指定時の変数指定子 */
  enum coef_ID 
  {
    amplitude = 0,
    frequency,
    initphase,
    bias
  };
  
  /** 速度変化のモード */
  enum velocity_variation 
  {
    vel_constant,
    vel_harmonic,
    vel_zero
  };
  
  
private:
  int state;             /// Fluid(1) or Solid(0)
  int mat_odr;           /// 媒質リストのインデクス
  std::string name;      ///< ラベル
  
  
public:
  
  /** コンストラクタ */
  CompoList() {
    mat_odr = 0;
    state = -1;
  }
  
  /**　デストラクタ */
  ~CompoList() {}
  
public:
  
  
  /**
   * @brief BCのラベル名を返す
   */
  std::string getBCstr() const
  {
    return "";
  }
  
  
  /**
   * @brief ラベル名を返す
   */
  std::string getLabel() const
  { 
    return name; 
  }
  

  //@brief mat_odrをセットする
  //@param key MediumListのエントリ番号
  void setMatOdr(const int key)
  { 
    mat_odr = key; 
  }


  //@brief ラベル名をセットする
  void setLabel(const std::string pnt)
  { 
    name = pnt; 
  }


  //@brief stateをセットする
  //@param key セルの状態 SOLID/FLUID
  void setState(const int key)
  { 
    state = key; 
  }

  
  int getState() const
  { 
    return state; 
  }

  
  int getMatOdr() const
  { 
    return mat_odr;
  }
  
};

#endif // _FB_COMPO_H_
