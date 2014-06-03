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
 * @file   VoxInfo.C
 * @brief  FlowBase VoxInfo class
 * @author kero
 */

#include "VoxInfo.h"


// #################################################################
/**
 @fn void VoxInfo::copyBCIbase(int* dst, int* src)
 @brief dst[]にsrc[]のstate, activeビットの情報をコピーする
 @param dst
 @param src
 */
void VoxInfo::copyBCIbase(int* dst, int* src)
{
  size_t nx = (size_t)(size[0]+2*guide) * (size_t)(size[1]+2*guide) * (size_t)(size[2]+2*guide);
  
  // 30, 31ビットのみコピー
  for (size_t m=0; m<nx; m++) {
    dst[m] = src[m] & 0xc0000000;
  }  
}


// #################################################################
/**
 @brief BCindexにそのセルが計算に有効(active)かどうかをエンコードする
 @param bx BCindex ID
 @note
 - IS_FLUID returns true if FLUID
 - 設定は内部領域のみ
 */
void VoxInfo::encActive(int* bx)
{
  size_t m;
  int s;
  
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        s = bx[m];
        if ( IS_FLUID( s ) )
        {
          s = onBit( s, ACTIVE_BIT );
        }
        else
        {
          s = offBit( s, ACTIVE_BIT );
        }
        bx[m] = s;
      }
    }
  }
}


// #################################################################
/**
 @brief ディリクレ条件とノイマン条件の排他性をチェックし，反復行列の非対角要素/対角要素の係数をエンコードする
 @param bx BCindex P
 @note
 - ディリクレ条件とノイマン条件の排他性のチェック
 - 非対角要素と対角要素の係数をエンコードする
 */
void VoxInfo::encPbit(int* bx)
{
  size_t m, flag;
  int coef;
  int s_e, s_w, s_n, s_s, s_t, s_b, ss;
  int d_e, d_w, d_n, d_s, d_t, d_b;
  int s;
  bool exclusive;
  
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  // ディリクレ条件とノイマン条件の排他性のチェック
  exclusive = true;
  
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        s = bx[m];
        flag = 0;
        
        // ノイマン条件：値がゼロのとき，BCがセットされている
        s_e = BIT_SHIFT(s, BC_N_E); 
        s_w = BIT_SHIFT(s, BC_N_W);
        s_n = BIT_SHIFT(s, BC_N_N);
        s_s = BIT_SHIFT(s, BC_N_S);
        s_t = BIT_SHIFT(s, BC_N_T);
        s_b = BIT_SHIFT(s, BC_N_B);
        
        // ディリクレ条件：値がゼロのとき，BCがセットされている
        d_e = BIT_SHIFT(s, BC_D_E); 
        d_w = BIT_SHIFT(s, BC_D_W);
        d_n = BIT_SHIFT(s, BC_D_N);
        d_s = BIT_SHIFT(s, BC_D_S);
        d_t = BIT_SHIFT(s, BC_D_T);
        d_b = BIT_SHIFT(s, BC_D_B);
        
        // 非対角要素の係数をエンコード
        if ( (s_e * d_e) == 0 ) s = offBit( s, BC_NDAG_E );
        if ( (s_w * d_w) == 0 ) s = offBit( s, BC_NDAG_W );
        if ( (s_n * d_n) == 0 ) s = offBit( s, BC_NDAG_N );
        if ( (s_s * d_s) == 0 ) s = offBit( s, BC_NDAG_S );
        if ( (s_t * d_t) == 0 ) s = offBit( s, BC_NDAG_T );
        if ( (s_b * d_b) == 0 ) s = offBit( s, BC_NDAG_B );
        
        bx[m] = s;
        
        if ( (s_e==0) && (d_e==0) ) flag++;
        if ( (s_w==0) && (d_w==0) ) flag++;
        if ( (s_n==0) && (d_n==0) ) flag++;
        if ( (s_s==0) && (d_s==0) ) flag++;
        if ( (s_t==0) && (d_t==0) ) flag++;
        if ( (s_b==0) && (d_b==0) ) flag++;
        
        if ( flag != 0)
        {
          Hostonly_ printf("\tDirichlet and Neumann BC are specified on the same face in cell (%d,%d,%d)\n", i,j,k);
          exclusive = false;
        }
      }
    }
  }
  if ( !exclusive ) Exit(0);
  
  // 対角要素の係数のチェックとエンコード
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        s = bx[m];
        
        s_e = BIT_SHIFT(s, BC_N_E); // 0-Neumann / 1-normal
        s_w = BIT_SHIFT(s, BC_N_W);
        s_n = BIT_SHIFT(s, BC_N_N);
        s_s = BIT_SHIFT(s, BC_N_S);
        s_t = BIT_SHIFT(s, BC_N_T);
        s_b = BIT_SHIFT(s, BC_N_B);
        
        ss = s_e + s_w + s_n + s_s + s_t + s_b;
        bx[m] = s | (ss<<BC_DIAG);
        
        if ( (ss == 0) && (BIT_IS_SHIFT(s,ACTIVE_BIT)) )
        {
          Hostonly_ printf("\tError : Coefficient of diagonal element is zero at (%d,%d,%d) : (wesnbt)[%1d %1d %1d %1d %1d %1d]\n", i,j,k,
                           IS_FLUID(bx[_F_IDX_S3D(i-1, j,   k,   ix, jx, kx, gd)]),
                           IS_FLUID(bx[_F_IDX_S3D(i+1, j,   k,   ix, jx, kx, gd)]),
                           IS_FLUID(bx[_F_IDX_S3D(i,   j-1, k,   ix, jx, kx, gd)]),
                           IS_FLUID(bx[_F_IDX_S3D(i,   j+1, k,   ix, jx, kx, gd)]),
                           IS_FLUID(bx[_F_IDX_S3D(i,   j,   k-1, ix, jx, kx, gd)]),
                           IS_FLUID(bx[_F_IDX_S3D(i,   j,   k+1, ix, jx, kx, gd)]) );
        }
        
      }
    }
  }
  
  // ゼロ割防止のためのダミー係数 >> 全領域
  size_t nx = (size[0]+2*guide) * (size[1]+2*guide) * (size[2]+2*guide);
  
  for (m=0; m<nx; m++) {
    s = bx[m];
    if ( ((s>>BC_DIAG) & 0x7) == 0 ) { // 0x7 = 3 bit
      bx[m] = s | (0x6<<BC_DIAG);
    }
  }
}



// #################################################################
/**
 @brief bcp[]に壁面境界の圧力ノイマン条件のビットフラグと固体に隣接するFセルに方向フラグ，収束判定の有効フラグをエンコードする
 @param[in,out] bx BCindex P
 @note 
 - 流体セルのうち，固体セルに隣接する面のノイマンフラグをゼロにする．ただし，内部領域のみ．
 - 固体セルに隣接する流体セルに方向フラグを付与する．全内部領域．
 */
void VoxInfo::encPbit_N_Binary(int* bx)
{
  size_t m_p, m_e, m_w, m_n, m_s, m_t, m_b;
  int s;
  
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;

  // ノイマンフラグ
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        m_p = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        s = bx[m_p];
        
        if ( IS_FLUID( s ) )
        {
          m_e = _F_IDX_S3D(i+1, j,   k,   ix, jx, kx, gd);
          m_w = _F_IDX_S3D(i-1, j,   k,   ix, jx, kx, gd);
          m_n = _F_IDX_S3D(i,   j+1, k,   ix, jx, kx, gd);
          m_s = _F_IDX_S3D(i,   j-1, k,   ix, jx, kx, gd);
          m_t = _F_IDX_S3D(i,   j,   k+1, ix, jx, kx, gd);
          m_b = _F_IDX_S3D(i,   j,   k-1, ix, jx, kx, gd);

          // X_MINUS
          if ( !IS_FLUID(bx[m_w]) )
          {
            s = offBit( s, BC_N_W );    // ノイマンフラグ
            s = onBit( s, FACING_W );   // wall locationフラグ
          }
          
          // X_PLUS
          if ( !IS_FLUID(bx[m_e]) )
          {
            s = offBit( s, BC_N_E );
            s = onBit( s, FACING_E );
          }
          
          // Y_MINUS
          if ( !IS_FLUID(bx[m_s]) )
          {
            s = offBit( s, BC_N_S );
            s = onBit( s, FACING_S );
          }
          
          // Y_PLUS
          if ( !IS_FLUID(bx[m_n]) )
          {
            s = offBit( s, BC_N_N );
            s = onBit( s, FACING_N );
          }
          
          // Z_MINUS
          if ( !IS_FLUID(bx[m_b]) )
          {
            s = offBit( s, BC_N_B );
            s = onBit( s, FACING_B );
          }
          
          // Z_PLUS
          if ( !IS_FLUID(bx[m_t]) )
          {
            s = offBit( s, BC_N_T );
            s = onBit( s, FACING_T );
          }
          
          // 収束判定の有効フラグ，計算内部領域の流体セルのみ
          s = onBit(s, VLD_CNVG);
          
          bx[m_p] = s;
        }
      }
    }
  }
}


// #################################################################
/**
 @brief 外部境界に接するセルにおいて，bx[]に圧力境界条件keyに対応するビットフラグを設定する
 @param face 外部境界面番号
 @param bx BCindex P
 @param key Dirichlet or Neumann
 @param dir 壁面の場合(true)，方向フラグをON
 @note 
 - 流体セルに対してのみ，1-Normal, 0-BC
 - 固体セルに隣接する面のノイマンフラグをゼロにし，方向フラグを立てる
 */
void VoxInfo::encPbit_OBC(const int face, int* bx, const std::string key, const bool dir)
{
  size_t m;
  int s;
  
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  switch (face)
  {
    case X_MINUS:
      if( nID[X_MINUS] < 0 )
      {
        int i = 1;
        int bit;
        if ("Neumann"==key) bit = BC_N_W;
        else                bit = BC_D_W;

        for (int k=1; k<=kx; k++) {
          for (int j=1; j<=jx; j++) {
            m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
            s = bx[m];
            if ( IS_FLUID(s) )
            {
              if (dir) s = onBit( s, FACING_W );
              bx[m] = offBit( s, bit );
            }
          }
        }
      }
      break;
      
    case X_PLUS:
      if( nID[X_PLUS] < 0 )
      {
        int i = ix;
        int bit;
        if ("Neumann"==key) bit = BC_N_E;
        else                bit = BC_D_E;

        for (int k=1; k<=kx; k++) {
          for (int j=1; j<=jx; j++) {
            m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
            s = bx[m];
            if ( IS_FLUID(s) )
            {
              if (dir) s = onBit( s, FACING_E );
              bx[m] = offBit( s, bit );
            }
          }
        }
      }
      break;
      
    case Y_MINUS:
      if( nID[Y_MINUS] < 0 )
      {
        int j = 1;
        int bit;
        if ("Neumann"==key) bit = BC_N_S;
        else                bit = BC_D_S;

        for (int k=1; k<=kx; k++) {
          for (int i=1; i<=ix; i++) {
            m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
            s = bx[m];
            if ( IS_FLUID(s) )
            {
              if (dir) s = onBit( s, FACING_S );
              bx[m] = offBit( s, bit );
            }
          }
        }
      }
      break;
      
    case Y_PLUS:
      if( nID[Y_PLUS] < 0 )
      {
        int j = jx;
        int bit;
        if ("Neumann"==key) bit = BC_N_N;
        else                bit = BC_D_N;

        for (int k=1; k<=kx; k++) {
          for (int i=1; i<=ix; i++) {
            m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
            s = bx[m];
            if ( IS_FLUID(s) )
            {
              if (dir) s = onBit( s, FACING_N );
              bx[m] = offBit( s, bit );
            }
          }
        }
      }
      break;
      
    case Z_MINUS:
      if( nID[Z_MINUS] < 0 )
      {
        int k = 1;
        int bit;
        if ("Neumann"==key) bit = BC_N_B;
        else                bit = BC_D_B;

        for (int j=1; j<=jx; j++) {
          for (int i=1; i<=ix; i++) {
            m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
            s = bx[m];
            if ( IS_FLUID(s) )
            {
              if (dir) s = onBit( s, FACING_B );
              bx[m] = offBit( s, bit );
            }
          }
        }
      }
      break;
      
    case Z_PLUS:
      if( nID[Z_PLUS] < 0 )
      {
        int k = kx;
        int bit;
        if ("Neumann"==key) bit = BC_N_T;
        else                bit = BC_D_T;

        for (int j=1; j<=jx; j++) {
          for (int i=1; i<=ix; i++) {
            m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
            s = bx[m];
            if ( IS_FLUID(s) )
            {
              if (dir) s = onBit( s, FACING_T );
              bx[m] = offBit( s, bit );
            }
          }
        }
      }
      break;

  } // end of switch
}


// #################################################################
/**
 @brief 外部境界に接するセルにおいて，各種速度境界条件に対応する媒質をチェックし，bv[]にビットフラグをセットする
 @param face 外部境界面番号
 @param bv BCindex V
 @param key fluid or solid　指定するBCが要求するガイドセルの状態 >> エラーチェックに使う
 @param enc_sw trueのとき，エンコードする．falseの場合にはガイドセルの状態チェックのみ
 @param chk ガイドセルの状態をチェックするかどうかを指定
 @param bp BCindex P
 @param enc_uwd trueのとき，1次精度のスイッチオン
 @note 
 - 外部境界条件の実装には，流束型とディリクレ型の2種類がある．
 - adjMedium_on_GC()でガイドセル上のIDを指定済み．指定BCとの適合性をチェックする
 */
void VoxInfo::encVbit_OBC(const int face, int* bv, const std::string key, const bool enc_sw, const std::string chk, int* bp, const bool enc_uwd)
{
  size_t m, mt;
  int s, q;
  
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  switch (face) {
    case X_MINUS:
      if( nID[X_MINUS] < 0 ){ // 外部境界をもつノードのみ
        for (int k=1; k<=kx; k++) {
          for (int j=1; j<=jx; j++) {
            m = _F_IDX_S3D(1, j, k, ix, jx, kx, gd);
            mt= _F_IDX_S3D(0, j, k, ix, jx, kx, gd);
            
            s = bv[m];
            q = bp[mt];
            
            if ( IS_FLUID(s) ) {
              // エンコード処理
              if ( enc_sw ) {
                bv[m]  = s | (OBC_MASK << BC_FACE_W); // OBC_MASK==31 外部境界条件のフラグ
              }
              
              // 外部境界で安定化のため，スキームを1次精度にする
              if ( enc_uwd ) bp[mt] = onBit(q, VBC_UWD);
              
            }// IS_FLUID
          }
        }        
      }
      break;
      
    case X_PLUS:
      if( nID[X_PLUS] < 0 ){
        for (int k=1; k<=kx; k++) {
          for (int j=1; j<=jx; j++) {
            m = _F_IDX_S3D(ix,   j, k, ix, jx, kx, gd);
            mt= _F_IDX_S3D(ix+1, j, k, ix, jx, kx, gd);
            
            s = bv[m];
            q = bp[mt];
            
            if ( IS_FLUID(s) ) {
              
              // エンコード処理
              if ( enc_sw ) {
                bv[m]  = s | (OBC_MASK << BC_FACE_E);
              }
              
              // 外部境界で安定化のため，スキームを1次精度にする
              if ( enc_uwd ) bp[mt] = onBit(q, VBC_UWD);
              
            }// IS_FLUID
          }
        }
      }
      break;
      
    case Y_MINUS:
      if( nID[Y_MINUS] < 0 ){
        for (int k=1; k<=kx; k++) {
          for (int i=1; i<=ix; i++) {
            m = _F_IDX_S3D(i, 1, k, ix, jx, kx, gd);
            mt= _F_IDX_S3D(i, 0, k, ix, jx, kx, gd);
            
            s = bv[m];
            q = bp[mt];
            
            if ( IS_FLUID(s) ) {
              
              // エンコード処理
              if ( enc_sw ) {
                bv[m]  = s | (OBC_MASK << BC_FACE_S);
              }
              
              // 外部境界で安定化のため，スキームを1次精度にする
              if ( enc_uwd ) bp[mt] = onBit(q, VBC_UWD);
              
            }// IS_FLUID
          }
        }
      }
      break;
      
    case Y_PLUS:
      if( nID[Y_PLUS] < 0 ){
        for (int k=1; k<=kx; k++) {
          for (int i=1; i<=ix; i++) {
            m = _F_IDX_S3D(i, jx,   k, ix, jx, kx, gd);
            mt= _F_IDX_S3D(i, jx+1, k, ix, jx, kx, gd);
            
            s = bv[m];
            q = bp[mt];
            
            if ( IS_FLUID(s) ) {
              
              // エンコード処理
              if ( enc_sw ) {
                bv[m]  = s | (OBC_MASK << BC_FACE_N);
              }
              
              // 外部境界で安定化のため，スキームを1次精度にする
              if ( enc_uwd ) bp[mt] = onBit(q, VBC_UWD);
              
            }// IS_FLUID
          }
        }
      }
      break;
      
    case Z_MINUS:
      if( nID[Z_MINUS] < 0 ){
        for (int j=1; j<=jx; j++) {
          for (int i=1; i<=ix; i++) {
            m = _F_IDX_S3D(i, j, 1, ix, jx, kx, gd);
            mt= _F_IDX_S3D(i, j, 0, ix, jx, kx, gd);
            
            s = bv[m];
            q = bp[mt];
            
            if ( IS_FLUID(s) ) {
              
              // エンコード処理
              if ( enc_sw ) {
                bv[m]  = s | (OBC_MASK << BC_FACE_B);
              }
              
              // 外部境界で安定化のため，スキームを1次精度にする
              if ( enc_uwd ) bp[mt] = onBit(q, VBC_UWD);
              
            }// IS_FLUID
          }
        }
      }
      break;
      
    case Z_PLUS:
      if( nID[Z_PLUS] < 0 ){
        for (int j=1; j<=jx; j++) {
          for (int i=1; i<=ix; i++) {
            m = _F_IDX_S3D(i, j, kx,   ix, jx, kx, gd);
            mt= _F_IDX_S3D(i, j, kx+1, ix, jx, kx, gd);
            
            s = bv[m];
            q = bp[mt];
            
            if ( IS_FLUID(s) ) {
              
              // エンコード処理
              if ( enc_sw ) {
                bv[m]  = s | (OBC_MASK << BC_FACE_T);
              }
              
              // 外部境界で安定化のため，スキームを1次精度にする
              if ( enc_uwd ) bp[mt] = onBit(q, VBC_UWD);
              
            }// IS_FLUID
          }
        }
      }
      break;
  } // end of switch
}


// #################################################################
/**
 @brief コンポーネントの操作に必要な定数の設定
 @param m_NoBC 境界条件数
 @param m_NoCompo コンポーネントの総数
 */
void VoxInfo::setNoCompo_BC(int m_NoBC, int m_NoCompo)
{
  NoBC    = m_NoBC;
  NoCompo = m_NoCompo;
}

// #################################################################
// bx[]に各境界条件の共通のビット情報をエンコードする（その1）
// 事前に，cmp[]へMediumListへのエントリ番号をエンコードしておく -> cmp[].setMatOdr()
void VoxInfo::setBCIndex_base1(int* bx, int* mid, const MediumList* mat, CompoList* cmp)
{
  int odr;
  int id;
  int s;
  int md;
  
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  size_t nx = (ix+2*gd) * (jx+2*gd) * (kx+2*gd); // ガイドセルを含む全領域を対象にする
  memset(bx, 0, sizeof(int)*nx);
  
  // セルの状態を流体で初期化
  for (size_t m=0; m<nx; m++) {
    bx[m] = onBit( bx[m], STATE_BIT );
  }
  
  // セルIDをエンコード
  for (size_t m=0; m<nx; m++) {
    md = mid[m];
    bx[m] |= (md << TOP_CELL_ID) ;
  }

  // MediumListのエントリをエンコードする
  for (int n=1; n<=NoCompo; n++) {
    id  = cmp[n].getMatOdr();

    for (size_t m=0; m<nx; m++) {
      if ( mid[m] == id ) bx[m] |= (id << TOP_MATERIAL);
      
    }
  }

  // CompoListのエントリ　セル要素のみ
  for (int n=1; n<=NoBC; n++) {
  }
  
  // 状態のエンコード
  for (size_t m=0; m<nx; m++) {
    s = bx[m];
    odr = DECODE_MAT( s );

    if ( mat[odr].getState() == FLUID )
    {
      s = onBit( s, STATE_BIT );
    }
    else  // SOLID
    {
      s = offBit( s, STATE_BIT );
    }
    bx[m] = s;
  }
  
}



// #################################################################
// bx[]に各境界条件の共通のビット情報をエンコードする（その2）
void VoxInfo::setBCIndex_base2(int* bx, int* mid, CompoList* cmp)
{
  
  // BCIndexにそのセルが計算に有効(active)かどうかをエンコードする．KindOfSolverによって異なる
  encActive(bx);
}


// #################################################################
// 圧力境界条件のビット情報をエンコードする
void VoxInfo::setBCIndexP(int* bcd, int* bcp, int* mid, BoundaryOuter* obc, CompoList* cmp)
{
  // 初期化 @note ビットを1に初期化する．初期化範囲はガイドセルを含む全領域．セルフェイスの射影処理で必要．
  size_t mx = (size[0]+2*guide) * (size[1]+2*guide) * (size[2]+2*guide);
  for (size_t m=0; m<mx; m++) {
    bcp[m] |= ( 0x3ffff << BC_NDAG_W ); // BC_NDAG_W〜BC_D_Tまで18bitまとめて1に初期化
  }

  // 計算領域内の壁面のNeumannBCのマスク処理と固体に隣接するFセルに方向フラグをエンコード
  encPbit_N_Binary(bcp);    // Binary

  // 外部境界のビットフラグをエンコード
  for (int face=0; face<NOFACE; face++) {
    int F = obc[face].get_Class();
    
    switch ( F )
    {
      case OBC_WALL:
        encPbit_OBC(face, bcp, "Neumann", true);
        break;
      default:
        Exit(0);
    }
  }

  // 内部境界のコンポーネントのエンコード
  for (int n=1; n<=NoBC; n++) {
  }

  // 全周Neumannフラグのセルと排他性をチェックし，反復行列の非対角要素/対角要素をエンコードする
  encPbit(bcp);
}



// #################################################################
// bv[]に境界条件のビット情報をエンコードする
void VoxInfo::setBCIndexV(int* bv, const int* mid, int* bp, BoundaryOuter* obc, CompoList* cmp)
{
  // ガイドセルの媒質情報をチェックし，流束形式のBCの場合にビットフラグをセット
  
  // 外部境界
  for (int face=0; face<NOFACE; face++) {
    int F = obc[face].get_Class();
    
    switch ( F )
    {
      case OBC_WALL:
        encVbit_OBC(face, bv, "solid", true, "check", bp, false); // 流束形式
        break;
      default:
        Exit(0);
    }
    
  }
  
  // 内部境界のコンポーネントのエンコード
  for (int n=1; n<=NoBC; n++) {
  }
  
}

