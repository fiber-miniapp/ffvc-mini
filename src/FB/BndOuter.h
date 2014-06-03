#ifndef _FB_BND_OUTER_H_
#define _FB_BND_OUTER_H_

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
 * @file   BndOuter.h
 * @brief  FlowBase BoundaryOuter Class Header
 * @author kero
 * @note メンバ変数に追加したら，dataCopy()処理にも加えること
 */

#include <string>
#include "FB_Define.h"
#include "cpm_Define.h"

class BoundaryOuter {
private:
  int BCclass;       ///< 境界条件の種類
  int wallType;      ///< wall >> (fixed, slide)
  int gc_medium;     ///< ガイドセルの媒質インデクス
  int v_profile;     ///< 速度プロファイル（constant, harmonic, zero）
  int pType;         ///< 外部境界の圧力指定(ディリクレ，勾配ゼロ)
  std::string label; ///< ラベル
  std::string alias; ///< 別名
  
public: 
  REAL_TYPE nv[3];   ///< 法線
  REAL_TYPE ca[5];   ///< 係数
  REAL_TYPE cb[5];   ///< 係数
  REAL_TYPE p;       ///< ワーク
  
  /** 壁面の種類 */
  enum wall_kind 
  {
    fixed,
    slide
  };
  
  /** コンストラクタ */
  BoundaryOuter() 
  {
    BCclass =  wallType = 0;
    pType = v_profile = gc_medium = 0;
    p = 0.0;
    for (int i=0; i<5; i++) ca[i] = cb[i] = 0.0;
    for (int i=0; i<3; i++) nv[i] = 0.0;
  }
  
  /**　デストラクタ */
  ~BoundaryOuter() {}
  
  
public:
  int get_Class() const
  { 
    return BCclass; 
  }
  
  int get_GuideMedium() const
  { 
    return gc_medium; 
  }

  int get_pType() const 
  { 
    return pType; 
  }
  
  int get_wallType() const 
  { 
    return wallType; 
  }
  
  int get_V_Profile() const 
  { 
    return v_profile; 
  }
  
  std::string get_Label() const 
  { 
    return label; 
  }
  
  std::string get_Alias() const 
  { 
    return alias;
  }
  
  void addVec         (REAL_TYPE* vec);
  void dataCopy       (BoundaryOuter* src);
  void set_Alias      (std::string key);
  void set_Class      (const int key);
  
  void set_GuideMedium(int key);
  void set_Label      (std::string key);
  void set_pType      (int key);
  void set_wallType   (const int key);
  void set_V_Profile  (const int key);
  
};

#endif // _FB_BND_OUTER_H_
