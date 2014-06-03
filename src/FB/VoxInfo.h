#ifndef _FB_BINVOX_H_
#define _FB_BINVOX_H_

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
 * @file   VoxInfo.h
 * @brief  FlowBase VoxInfo class Header
 * @author kero
 */

#include <math.h>
#include <cstdlib>

#include "cpm_ParaManager.h"

#include "DomainInfo.h"
#include "FBUtility.h"
#include "Component.h"
#include "Medium.h"
#include "BndOuter.h"


class VoxInfo : public DomainInfo {

private:
  int NoBC;                      ///< 境界条件数
  int NoCompo;                   ///< コンポーネントの総数

public:
  /** コンストラクタ */
  VoxInfo() {
    NoBC = 0;
    NoCompo = 0;
  }
  
  /**　デストラクタ */
  ~VoxInfo() {}
  
private:

  void encPbit_N_Binary(int* bx);
  
  void encActive             (int* bx);
  void encPbit               (int* bx);
  void encPbit_OBC           (int face, int* bx, std::string key, bool dir);
  void encVbit_OBC           (int face, int* bv, std::string key, const bool enc_sw, std::string chk, int* bp, const bool enc_uwd);
  
  
  //@fn inline int offBit(int idx, const int shift)
  //@brief idxの第shiftビットをOFFにする
  inline int offBit(int idx, const int shift) {
    return ( idx & (~(0x1<<shift)) );
  }
  
  //@fn inline int onBit(int idx, const int shift)
  //@brief idxの第shiftビットをONにする
  inline int onBit(int idx, const int shift) {
    return ( idx | (0x1<<shift) );
  }
  
  
public:

  /**
   @fn void VoxInfo::copyBCIbase(int* dst, int* src)
   @brief dst[]にsrc[]のstate, activeビットの情報をコピーする
   @param dst
   @param src
   */
  void copyBCIbase           (int* dst, int* src);
  
  
  /**
   * @brief bx[]に各境界条件の共通のビット情報をエンコードする（その1）
   * @param [in,out] bx   BCindex ID
   * @param [in,out] mid  ID配列
   * @param [in]     mat  MediumList
   * @param [in,out] cmp  CompoList
   * @note 事前に，cmp[]へMediumListへのエントリ番号をエンコードしておく -> cmp[].setMatOdr()
   */
  void setBCIndex_base1(int* bd, int* mid, const MediumList* mat, CompoList* cmp);
  
  
  /**
   * @brief bx[]に各境界条件の共通のビット情報をエンコードする（その2）
   * @param [out]    bx    BCindex ID
   * @param [in,out] mid   ID配列
   * @param [in,out] cmp   CompoList
   */
  void setBCIndex_base2(int* bx, int* mid, CompoList* cmp);
  
  
  /**
   * @brief 圧力境界条件のビット情報をエンコードする
   * @param [in,out] bcd   BCindex ID
   * @param [in,out] bcp   BCindex P
   * @param [in,out] mid   ID配列
   * @param [in]     obc   外部境界リスト
   * @param [in,out] cmp   CompoList
   */
  void setBCIndexP(int* bcd, int* bcp, int* mid, BoundaryOuter* obc, CompoList* cmp);
  
  
  /**
   * @brief bv[]に境界条件のビット情報をエンコードする
   * @param [in,out] bv BCindex V
   * @param [in] mid ID配列
   * @param [in,out] bp BCindex P
   * @param [in]     obc   外部境界リスト
   * @param [in,out] cmp CompoList
   */
  void setBCIndexV(int* bv, const int* mid, int* bp, BoundaryOuter* obc, CompoList* cmp);
  

  void setNoCompo_BC         (int m_NoBC, int m_NoCompo);

  
};

#endif // _FB_BINVOX_H_
