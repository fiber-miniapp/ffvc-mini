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

//@file   BndOuter.C
//@brief  FlowBase BoundaryOuter class
//@author kero

#include "BndOuter.h"

// #################################################################
// ラベルを設定
void BoundaryOuter::set_Label(std::string key)
{
  label = key;
}


// #################################################################
// aliasラベルを設定する
void BoundaryOuter::set_Alias(std::string key)
{
  alias = key;
}


// #################################################################
// ガイドセルの媒質IDをセットする
void BoundaryOuter::set_GuideMedium(int key)     
{ 
  gc_medium = key;
}


// #################################################################
// 境界条件の種類をセットする
void BoundaryOuter::set_Class(const int key)     
{ 
  BCclass = key;
}


// #################################################################
// 壁面境界のモードをセットする
void BoundaryOuter::set_wallType(const int key)     
{ 
  wallType = key;
}



// #################################################################
// 外部境界の圧力指定
void BoundaryOuter::set_pType(int key)
{
  pType  = key;
}


// #################################################################
// 速度プロファイルの指定
void BoundaryOuter::set_V_Profile(const int key)
{
  v_profile  = key;
}


// #################################################################
// ベクトルのコピー
void BoundaryOuter::addVec(REAL_TYPE* vec) 
{
  nv[0] = vec[0];
  nv[1] = vec[1];
  nv[2] = vec[2];
}


// #################################################################
// メンバー変数のコピー
void BoundaryOuter::dataCopy(BoundaryOuter* src)
{
  BCclass   = src->BCclass;
  gc_medium = src->gc_medium;
  pType     = src->pType;
  v_profile = src->v_profile;
  p         = src->p;
  wallType  = src->wallType;
  
  label     = src->label;
  alias     = src->alias;
  
  for (int i=0; i<3; i++) {
    nv[i] = src->nv[i];
  }
  for (int i=0; i<5; i++) {
    ca[i] = src->ca[i];
    cb[i] = src->cb[i];
  }
}
