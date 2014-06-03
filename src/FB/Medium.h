#ifndef _FB_MEDIUM_H_
#define _FB_MEDIUM_H_

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
 * @file   Medium.h
 * @brief  FlowBase Medium class Header
 * @author kero
 */

#include <string>
#include "FB_Define.h"
#include "cpm_Define.h"


class MediumList {
private:
  int  state;         ///< solid or fluid
  std::string name;   ///< ラベル

public:
  
  MediumList() {
    state  = -1;
  }
  ~MediumList() {}
  
public:

  /**
   * @brief 媒質の属性を取得
   * @return Fluid or Solid
   */
  int getState() const
  { 
    return state; 
  }
  
  
  /**
   * @brief ラベルを取得
   * @return 文字列
   */
  std::string getLabel() const
  { 
    return name; 
  }
  
  
  /**
   * @brief 状態をセット
   * @param[in] key fluid(1) or solid(0)
   */
  void setState(const int key) 
  { 
    state = key; 
  }
  
  
  /**
   * @brief ラベルをセット
   * @param[in] key 文字列
   */
  void setLabel(const std::string key) 
  {
    name = key;
  }

};

#endif // _FB_MEDIUM_H_
