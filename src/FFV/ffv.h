#ifndef _FFV_H_
#define _FFV_H_

// #################################################################
//
// FFV : Frontflow / violet Cartesian
//
// Copyright (c) 2012-2014 Institute of Industrial Science, The University of Tokyo.
// All right reserved.
//
// Copyright (c) 2012-2014 Advanced Institute for Computational Science, RIKEN.
// All right reserved.
//
// #################################################################
//
// 以下のマクロはcpm_Define.hで定義されている
//   REAL_TYPE
//   X_MINUS, Y_MINUS, Z_MINUS, X_PLUS, Y_PLUS, Z_PLUS
//   X_DIR, Y_DIR, Z_DIR
//   PLUS2MINUS, MINUS2PLUS, BOTH

/** 
 * @file   ffv.h
 * @brief  FFV Class Header
 * @author kero
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <iostream>
#include <fstream>

#include "cpm_ParaManager.h"

#include "DomainInfo.h"

#include "FB_Define.h"
#include "ffv_Define.h"
#include "FBUtility.h"
#include "ItrCtl.h"
#include "BndOuter.h"
#include "Medium.h"
#include "VoxInfo.h"
#include "ffv_SetBC.h"
#include "History.h"
#include "FileIO.h"
#include "dfi.h"
#include "ffv_Ffunc.h"
#include "FB_Ffunc.h"

#include "CommandLine.h"

#ifdef _OPENMP
#include "omp.h"
#endif


using namespace std;


class FFV : public DomainInfo {
private:

  /** 対流項スキーム */
  enum convection_scheme 
  {
    O1_upwind=1,
    O2_central,
    O3_muscl,
    O4_central
  };

  int ffv_procGrp;         ///< プロセスグループ番号 => 0
  
  unsigned long G_Acell;   ///< グローバルなActive cell
  
  double CurrentTime;           ///< 計算開始からの積算時刻（ケース）
  double Session_StartTime;     ///< セッションの開始時間
  double Session_CurrentTime;   ///< セッションの現在時間
  
  double step_start;            ///< 1stepのelapse time(sec)
  
  unsigned CurrentStep;         ///< 計算開始からの積算ステップ（ケース）
  unsigned Session_StartStep;   ///< セッションの開始ステップ
  unsigned Session_CurrentStep; ///< セッションの現在のステップ
  unsigned Session_LastStep;    ///< セッションで計算するステップ数
  
  REAL_TYPE convergence_prev;  ///< 前回の反復の収束値
  REAL_TYPE convergence_rate;  ///< 収束値の増減比
  
  REAL_TYPE deltaT; ///< 時間積分幅（無次元）

  int cf_sz[3];     ///< SOR2SMAの反復の場合のバッファサイズ
  REAL_TYPE *cf_x;  ///< i方向のバッファ
  REAL_TYPE *cf_y;  ///< j方向のバッファ
  REAL_TYPE *cf_z;  ///< k方向のバッファ

  // dfi ファイル管理
  int dfi_mng[var_END];

  
  REAL_TYPE v00[4];      ///< 参照速度
  
  
  // データ領域ポインタ
  
  // Vector3D
  REAL_TYPE *d_v;
  REAL_TYPE *d_vc;
  REAL_TYPE *d_v0;
  REAL_TYPE *d_wv;
  REAL_TYPE *d_vf0;
  REAL_TYPE *d_wo;
  
  // Scalar3D
  int *d_mid;
  int *d_bcd;
  int *d_bcp;
  int *d_bcv;
  
  REAL_TYPE *d_p;   ///< 圧力
  REAL_TYPE *d_p0;  ///< 圧力（1ステップ前）
  REAL_TYPE *d_ws;  ///< 反復中に固定のソース
  REAL_TYPE *d_sq;  ///< 反復中に変化するソース
  REAL_TYPE *d_dv;  ///< div(u)の保存
  REAL_TYPE *d_b;   ///< Ax=bの右辺ベクトル
  
  
  FILE *fp_b;  ///< 基本情報
  
  ItrCtl IC[ItrCtl::ic_END]; ///< 反復情報管理クラス
  MediumList* mat;           ///< 媒質リスト
  CompoList* cmp;            ///< コンポーネントリスト
  VoxInfo V;                 ///< ボクセル前処理クラス
  SetBC3D BC;                ///< BCクラス
  History* H;                ///< 履歴クラス
  FileIO F;                  ///< ファイル入出力クラス
  DFI DFI_;                  ///< 分散ファイルインデクス管理クラス

  CommandLine CL;

  int NoBC;
  int NoCompo;
  int NoMedium;

  REAL_TYPE RefLength;
  REAL_TYPE RefVelocity;
  REAL_TYPE Reynolds;
  REAL_TYPE Tscale;
  REAL_TYPE CFL;
  
  int PM_Test;

  int CnvScheme;

  double TimeAccel;   ///< 加速時間（無次元）

  const static bool Log_Base = true;

  int OutputInterval;
  bool OutputGather;


public:
  /** コンストラクタ */
  FFV();
  
  /**　デストラクタ */
  ~FFV();
  
  
private:

  /**
   * @brief 並列分散時のファイル名の管理を行う
   */
  void setDFI();

  /**
   * @brief 主計算部分に用いる配列のアロケーション
   * @param [in,out] total ソルバーに使用するメモリ量
   */
  void allocArray_Main(double &total);
  
  
  /**
   * @brief 前処理に用いる配列のアロケーション
   * @param [in,out] prep  前処理に使用するメモリ量
   * @param [in,out] total ソルバーに使用するメモリ量
   */
  void allocArray_Prep(double &prep, double &total);
  
  
  /**
   * @brief 主計算部分に用いる配列のアロケーション
   * @param [in,out] total ソルバーに使用するメモリ量
   */
  void allocate_Main(double &total);
  
  
  /**
   * @brief SOR2SMAのバッファ確保
   * @param [in,out] total ソルバーに使用するメモリ量
   */
  void allocate_SOR2SMA_buffer(double &total);
  
  
  /**
   * @brief 参照速度を設定する
   * @param [in] time 設定する時刻
   */
  void setV00(double m_time);
  
   
  /**
   * @brief 制御パラメータ，物理パラメータの表示
   * @param [in]  fp   ファイルポインタ
   */
  void display_Parameters(FILE* fp);
  
  
  /**
   * @brief メモリ消費情報を表示
   * @param [in]     fp    ファイルポインタ
   * @param [in]     L_mem ローカルメモリサイズ
   * @param [in]     str   表示用文字列
   */
  void display_memory_info(FILE* fp, double L_mem, const char* str);
  
  
  /** 
   * @brief 計算領域情報を設定する
   */
  void DomainInitialize();
  
  /**
   * @brief グローバルな領域情報を表示する
   * @param [in] fp      ファイルポインタ
   */
  void printGlobalDomain(FILE* fp);
  
  
  /** 1ステップのコアの処理
   * @param [in] m_step   現在のステップ数
   */
  int Loop(const unsigned m_step);
  
  
  /**
   * @brief 線形ソルバーの選択実行
   * @param [in]  IC       ItrCtlクラス
   * @param [in]  rhs_nrm  Poisson定数項ベクトルの自乗和ノルム
   * @param [in]  res_init 初期残差ベクトル
   */
  void LS_Binary(ItrCtl* IC, const double rhs_nrm, const double rhs_init);
  
  
  /**
   * @brief VP反復の発散値を計算する
   * @param [in] IC ItrCtlクラス
   */
  void Norm_Div(ItrCtl* IC);
  
  
  /**
   * @brief Fractional Step法でNavier-Stokes方程式を解く．バイナリ近似．
   */
  void NS_FS_E_Binary();
  
  
  /**
   * @brief 履歴の出力準備
   */
  void prep_HistoryOutput();

  
  /**
   * @brief 外部境界条件の各面の情報を表示する
   * @param [in] fp    ファイルポインタ
   */
  void printFaceOBC(FILE* fp);
  
  
  /**
   * @brief 単媒質に対する熱伝導方程式を陰解法で解く
   * @param [in]  IC       IterationCtlクラス
   * @param [in]  rhs_nrm  Poisson定数項ベクトルの自乗和ノルム
   * @param [in]  r0       初期残差ベクトル
   */
  void ps_LS(ItrCtl* IC, const double rhs_nrm, const double r0);
  
  
  /**
   * @brief 初期条件の設定
   */
  void setInitialCondition();
  
  
  /** モデルをセットアップ
   */
  void setModel();


  /**
   * @brief パラメータの設定
   */
  void setParameters();

  void displayParameters();
  

  /** 2色オーダリングSORのストライドメモリアクセス版
   * @retval 反復数
   * @param [in]     IC      IterationCtlクラス
   * @param [in,out] x       解ベクトル
   * @param [in]     b  RHS  vector
   * @param [in]     rhs_nrm RHS vector
   * @param [in]     r0      初期残差ベクトル
   */
  int SOR_2_SMA(ItrCtl* IC, REAL_TYPE* x, REAL_TYPE* b, const double rhs_nrm, const double r0);
    
  
  /**
   * @brief 空間平均操作と変動量の計算を行う
   * @param [out]    avr  平均値
   * @param [out]    rms  変動値
   * @param [in,out] flop 浮動小数演算数
   */
  void Variation_Space(double* avr, double* rms, double& flop);
  
  
  /**
   * @brief BCIndexにビット情報をエンコードする
   */
  void VoxEncode();


  /**
   * @brief ファイル出力
   */
  void FileOutput();
  
  
public:
  
  /**
   * @brief CPMのポインタをコピーし、ランク情報を設定
   * @param [in] m_paraMngr  cpm_ParaManagerクラス
   * @return  エラーコード
   */
  bool importCPM(cpm_ParaManager* m_paraMngr)
  {
    if ( !m_paraMngr ) return false;
    paraMngr = m_paraMngr;
    
    setRankInfo(paraMngr, procGrp);
    
    return true;
  }
  
  
  /**  
   * @brief 初期化格子生成、ビットフラグ処理ほか
   * @param [in] argc  main関数の引数の個数
   * @param [in] argv  main関数の引数リスト
   */
  int Initialize(int argc, char **argv);
  
  
  /** 
   * @brief マスターノードのみ trueを返す
   * @return true(Rank==0) / false(Rank!=0)
   */
  bool IsMaster() const
  {
    return ( paraMngr->GetMyRankID() == 0 ) ? true : false;
  }
  
  
  /** 
   * @brief シミュレーションの1ステップの処理
   *  Loop() + stepPost()
   */
  int MainLoop();
  
  
  /** 
   * @brief シミュレーションの終了時の処理
   * プロファイルの統計処理ほか
   */
  bool Post();
  
  
};

#endif // _FFV_H_
