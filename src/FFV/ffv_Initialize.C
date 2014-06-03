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

/** 
 * @file   ffv_Initialize.C
 * @brief  FFV Class
 * @author kero
 */

#include "ffv.h"
#include <climits>   // for UINT_MAX
#include <ctime>   // for time, ctime

int FFV::Initialize(int argc, char **argv)
{
  // コマンドライン引数解析
  CL.parse(argc, argv);
  
  // メッセージ表示
  Hostonly_
  {
    printf("**************************************************************\n");
    printf("FFVC_MINI %s\n", FFVC_MINI_VERSION);
    printf("   based on FFV-C %s (FFV %d.%d.%d, FB %d.%d.%d), CPM %s\n",
           FFVC_VERSION,
           FFV_VERS/100, (FFV_VERS/10)%10, FFV_VERS%10,
           FB_VERS/100, (FB_VERS/10)%10, FB_VERS%10,
           CPM_VERSION_NO);
    printf("\n");
    printf("  Command line:\n  ");
    for(int i = 0; i < argc; i++) printf(" %s", argv[i]);
    printf("\n");
    printf("  Processes: %d\n", paraMngr->GetNumRank());
    int num_thread = 1;
#ifdef _OPENMP
    num_thread = omp_get_max_threads();
#endif
    printf("  Threads: %d\n", num_thread);
    printf("  Host: %s\n", paraMngr->GetHostName().c_str());
    time_t t = time(0);
    printf("  Date: %s", ctime(&t));
    printf("**************************************************************\n");
  }

  double TotalMemory   = 0.0;  ///< 計算に必要なメモリ量（ローカル）
  double PrepMemory    = 0.0;  ///< 初期化に必要なメモリ量（ローカル）
  
  // ランク情報をセット
  V.setRankInfo(paraMngr, procGrp);
  BC.setRankInfo(paraMngr, procGrp);
  F.setRankInfo(paraMngr, procGrp);
  
  // 領域設定
  DomainInitialize();
  
  // 各クラスで領域情報を保持
  V.setNeighborInfo(guide);
  BC.setNeighborInfo(guide);
  F.setNeighborInfo(guide);

  
  // 前処理に用いるデータクラスのアロケート -------------------------------------
  allocArray_Prep(PrepMemory, TotalMemory);

  // メモリ消費量の情報を表示
  Hostonly_
  {
    printf("\n");
  }
  display_memory_info(stdout, PrepMemory, "Preprocessor");
  
  // モデルの設定
  setModel();
  
  // BCIndexにビット情報をエンコードとコンポーネントインデクスの再構築
  VoxEncode();
  
  // パラメータの設定
  setParameters();
  
  // 必要なパラメータをSetBC3Dクラスオブジェクトにコピー
  BC.setControlVars(Reynolds, RefVelocity, RefLength);

  // mid[]を解放する  -----------------------------------------------------------
  if ( d_mid ) delete [] d_mid;
  
  // ここまでがボクセル準備の時間セクション
  
  // 計算に用いる配列のアロケート -----------------------------------------------
  allocate_Main(TotalMemory);

  // 分散時のインデクスファイル生成
  setDFI();
  
  // スタート処理
  CurrentTime = 0.0;
  CurrentStep = 0;
  setV00(CurrentTime);
  
  // 制御パラメータ，物理パラメータの表示
  Hostonly_ 
  {
    display_Parameters(stdout);
  }
  
  // 初期条件の条件設定
  setInitialCondition();

  // SOR2SMA
  allocate_SOR2SMA_buffer(TotalMemory);
  
  // メモリ使用量の表示
  display_memory_info(stdout, TotalMemory, "Solver");

  // 履歴出力準備
  Hostonly_
  {
    printf("\n");
    printf(">> Start time step loop\n\n");
  }
  prep_HistoryOutput();
  
  return 1;
}


// #################################################################
// 計算領域情報を設定する
void FFV::DomainInitialize()
{
  // ガイドセル数
  const int guide = 2;

  // 強スケーリング
  if (CL.scale == "strong") {
    // 計算領域
    G_region[0] = G_region[1] = G_region[2] = 1.0;

    // 原点座標
    G_origin[0] = G_origin[1] = G_origin[2] = - 0.5;

    // ボクセル数
    G_size[0] = G_size[1] = G_size[2] = CL.size;
  }
  // 弱スケーリング
  else if (CL.scale == "weak") {
    for (int i = 0; i < 3; i++) {
      // 計算領域
      G_region[i] = 1.0 * CL.division[i];

      // 原点座標
      G_origin[i] = - 0.5 * G_region[i];

      // ボクセル数
      G_size[i] = CL.size * CL.division[i];
    }
  }
  else {
    Exit(0);
  }
  
  pitch[0] = pitch[1] = pitch[2] = 1.0 / CL.size;
  
  G_division[0] = CL.division[0];
  G_division[1] = CL.division[1];
  G_division[2] = CL.division[2];

  // 袖通信の最大数
  size_t n_vc  = (size_t)guide;
  size_t n_cmp = 3;
  
  int sz[3]  = {G_size[0], G_size[1], G_size[2]};
  int div[3] = {G_division[0], G_division[1], G_division[2]};
  
  double org[3] = {(double)G_origin[0], (double)G_origin[1], (double)G_origin[2]};
  double reg[3] = {(double)G_region[0], (double)G_region[1], (double)G_region[2]};
  
  int ret = paraMngr->VoxelInit(div, sz, org, reg, n_vc, n_cmp);
  if ( ret != CPM_SUCCESS )
  {
    cout << "Domain decomposition error : " << ret << endl;
    Exit(0);
  }

  // 分割後のパラメータをDomainInfoクラスメンバ変数に保持
  setNeighborInfo(guide);

  // チェック
  unsigned long tz = (unsigned long)size[0] * (unsigned long)size[1] * (unsigned long)size[2];
  if ( tz >= UINT_MAX )
  {
    Hostonly_ stamped_printf("\n\tError : Product of size[] exceeds UINT_MAX\n\n");
    Exit(0);
  }

  // 領域情報の表示
  Hostonly_
  {
    printf("\n");
    printf(">> Global Domain Information\n\n");
    printGlobalDomain(stdout);
  }

}


// #################################################################
// モデルの設定
void FFV::setModel()
{
  const int id_FLUID = 1;
  const int id_SOLID = 2;
  NoMedium = 2;
 
  // 媒質リスト
  mat = new MediumList[NoMedium+1];

  mat[id_FLUID].setState(FLUID);
  mat[id_FLUID].setLabel("Air");

  mat[id_SOLID].setState(SOLID);
  mat[id_SOLID].setLabel("Fe");

  // 内部境界条件の個数
  NoBC = 0;

  // コンポーネントの数
  NoCompo = NoBC + NoMedium;

  V.setNoCompo_BC(NoBC, NoCompo);

  // CompoListクラスをインスタンス．
  // [0]はダミーとして利用しないので，配列の大きさはプラス１する
  cmp = new CompoList[NoCompo+1];

  // d_midをゼロで初期化
  size_t mt = (size[0]+2*guide) * (size[1]+2*guide) *(size[2]+2*guide) * sizeof(int);
  memset(d_mid, 0, mt);
  
  for (int k=1; k<=size[2]; k++) {
    for (int j=1; j<=size[1]; j++) {
      for (int i=1; i<=size[0]; i++) {
        size_t m = _F_IDX_S3D(i, j, k, size[2], size[1], size[0], guide);
        d_mid[m] = id_FLUID;
      }
    }
  }
  
  // midのガイドセル同期
  if ( paraMngr->BndCommS3D(d_mid, size[0], size[1], size[2], guide, 1) != CPM_SUCCESS ) Exit(0);

  G_Acell = (unsigned long)size[0] * (unsigned long)size[1] * (unsigned long)size[2];
  if ( paraMngr->IsParallel() )
  {
    unsigned long tmp = G_Acell;
    if ( paraMngr->Allreduce(&tmp, &G_Acell, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
  }
  Hostonly_
  {
    printf("\n");
    printf(">> Medium List\n\n");
    printf("\tID = %d, Medium = %s, Label = %-4s, Element = %ld\n", id_FLUID, "FLUID", mat[id_FLUID].getLabel().c_str(), G_Acell);
    printf("\tID = %d, Medium = %s, Label = %-4s, Element = %ld\n", id_SOLID, "SOLID", mat[id_SOLID].getLabel().c_str(), 0UL);
  }

  // 外部境界条件
  // Basic Outer BCリスト
  BoundaryOuter BaseBc[2];

  BaseBc[0].set_Alias("outer_wall");
  BaseBc[0].set_Class(OBC_WALL);
  BaseBc[0].set_Label("Wall");
  BaseBc[0].set_wallType(BoundaryOuter::fixed);
  BaseBc[0].set_V_Profile(CompoList::vel_zero);

  BaseBc[1].set_Alias("slide_wall");
  BaseBc[1].set_Class(OBC_WALL);
  BaseBc[1].set_Label("wall");
  BaseBc[1].set_wallType(BoundaryOuter::slide);
  BaseBc[1].set_V_Profile(CompoList::vel_constant);
  REAL_TYPE v[3] = {1.0, 0.0, 0.0};
  BaseBc[1].addVec(v);
  BaseBc[1].ca[CompoList::bias] = 1.0;
  
  // 各フェイスに境界条件を設定する
  BoundaryOuter* bc = BC.export_OBC();

  bc[X_MINUS].dataCopy(&BaseBc[0]);  // outer_wall
  bc[X_MINUS].set_GuideMedium(id_SOLID);
  if ( nID[X_MINUS] < 0 ) {
    for (int k=1; k<=size[2]; k++) {
      for (int j=1; j<=size[1]; j++) {
        size_t m = _F_IDX_S3D(0, j, k, size[0], size[1], size[2], guide);
        d_mid[m] = id_SOLID;
      }
    }
  }

  bc[X_PLUS ].dataCopy(&BaseBc[0]);  // outer_wall
  bc[X_PLUS ].set_GuideMedium(id_SOLID);
  if ( nID[X_PLUS] < 0 ) {
    for (int k=1; k<=size[2]; k++) {
      for (int j=1; j<=size[1]; j++) {
        size_t m = _F_IDX_S3D(size[0]+1, j, k, size[0], size[1], size[2], guide);
        d_mid[m] = id_SOLID;
      }
    }
  }

  bc[Y_MINUS].dataCopy(&BaseBc[0]);  // outer_wall
  bc[Y_MINUS].set_GuideMedium(id_SOLID);
  if ( nID[Y_MINUS] < 0 ) {
    for (int k=1; k<=size[2]; k++) {
      for (int i=1; i<=size[0]; i++) {
        size_t m = _F_IDX_S3D(i, 0, k, size[0], size[1], size[2], guide);
        d_mid[m] = id_SOLID;
      }
    }
  }

  bc[Y_PLUS ].dataCopy(&BaseBc[0]);  // outer_wall
  bc[Y_PLUS ].set_GuideMedium(id_SOLID);
  if ( nID[Y_PLUS] < 0 ) {
    for (int k=1; k<=size[2]; k++) {
      for (int i=1; i<=size[0]; i++) {
        size_t m = _F_IDX_S3D(i, size[1]+1, k, size[0], size[1], size[2], guide);
        d_mid[m] = id_SOLID;
      }
    }
  }

  bc[Z_MINUS].dataCopy(&BaseBc[0]);  // outer_wall
  bc[Z_MINUS].set_GuideMedium(id_SOLID);
  if ( nID[Z_MINUS] < 0 ) {
    for (int j=1; j<=size[1]; j++) {
      for (int i=1; i<=size[0]; i++) {
        size_t m = _F_IDX_S3D(i, j, 0, size[0], size[1], size[2], guide);
        d_mid[m] = id_SOLID;
      }
    }
  }

  bc[Z_PLUS ].dataCopy(&BaseBc[1]);  // slide_wall
  bc[Z_PLUS ].set_GuideMedium(id_SOLID);
  if ( nID[Z_PLUS] < 0 ) {
    for (int j=1; j<=size[1]; j++) {
      for (int i=1; i<=size[0]; i++) {
        size_t m = _F_IDX_S3D(i, j, size[2]+1, size[0], size[1], size[2], guide);
        d_mid[m] = id_SOLID;
      }
    }
  }
 
  // セルIDのノード間同期
  if ( paraMngr->BndCommS3D(d_mid, size[0], size[1], size[2], guide, 1) != CPM_SUCCESS ) Exit(0);
  
  // 媒質情報の登録
  for (int i=1; i<=NoMedium; i++) {
    cmp[NoBC+i].setState( mat[i].getState() );
    cmp[NoBC+i].setLabel( mat[i].getLabel() );
    cmp[NoBC+i].setMatOdr(i);
  }

  // 外部境界面の表示
  Hostonly_
  {
    printf("\n");
    printf(">> Outer Boundary Conditions\n\n");
    printFaceOBC(stdout);
  }

}


// #################################################################
// BCIndexにビット情報をエンコードする
void FFV::VoxEncode()
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  // 基本ビット情報（Active, State, コンポ，媒質情報）を全領域についてエンコードする
  V.setBCIndex_base1(d_bcd, d_mid, mat, cmp);

  // bcdの同期
  if ( paraMngr->BndCommS3D(d_bcd, ix, jx, kx, gd, 1) != CPM_SUCCESS ) Exit(0);
  
  V.setBCIndex_base2(d_bcd, d_mid, cmp);

  // STATEとACTIVEビットのコピー
  V.copyBCIbase(d_bcp, d_bcd);
  V.copyBCIbase(d_bcv, d_bcd);

  // BCIndexP に圧力計算のビット情報をエンコードする -----
  V.setBCIndexP(d_bcd, d_bcp, d_mid, BC.export_OBC(), cmp);

  // BCIndexV に速度計算のビット情報をエンコードする -----
  V.setBCIndexV(d_bcv, d_mid, d_bcp, BC.export_OBC(), cmp);

  // bcd/bcp/bcvの同期
  if ( paraMngr->BndCommS3D(d_bcd, ix, jx, kx, gd, 1) != CPM_SUCCESS ) Exit(0);
  if ( paraMngr->BndCommS3D(d_bcp, ix, jx, kx, gd, 1) != CPM_SUCCESS ) Exit(0);
  if ( paraMngr->BndCommS3D(d_bcv, ix, jx, kx, gd, 1) != CPM_SUCCESS ) Exit(0);

}


// #################################################################
// パラメータの設定
void FFV::setParameters()
{
  if (CL.practical) {
    PM_Test = OFF;
  }
  else {
    PM_Test = ON;
  }

  CnvScheme = O3_muscl;

  RefLength = 1.0;
  RefVelocity = 1.0;
  Reynolds = 1000.0;
  Tscale = RefLength / RefVelocity;

  Session_LastStep = CL.step;
  CFL = CL.dt;

  TimeAccel = 1.0;

  // 無次元時間積分幅
  REAL_TYPE vRef = 1.0;
  deltaT = deltaX * CFL / vRef;

  // Poisson Iteration
  IC[ItrCtl::ic_prs_pr].set_ItrMax(CL.p_itr);
  IC[ItrCtl::ic_prs_pr].set_eps(1.0e-4);
  IC[ItrCtl::ic_prs_pr].set_normType(ItrCtl::r_b);
  IC[ItrCtl::ic_prs_pr].set_omg(1.2);
  IC[ItrCtl::ic_prs_pr].set_LS(SOR2SMA);
  if (CL.comm_mode == "sync") {
    IC[ItrCtl::ic_prs_pr].set_SyncMode(comm_sync);
  }
  else if (CL.comm_mode == "async") {
    IC[ItrCtl::ic_prs_pr].set_SyncMode(comm_async);
  }
  else {
    Exit(0);
  }

  // V-P Iteration
  IC[ItrCtl::ic_div].set_ItrMax(CL.vp_itr);
  IC[ItrCtl::ic_div].set_eps(1.0e-4);
  IC[ItrCtl::ic_div].set_normType(ItrCtl::v_div_max);

  // 出力
  OutputInterval = CL.output_interval;

}


// #################################################################
// 制御パラメータ，物理パラメータの表示
void FFV::display_Parameters(FILE* fp)
{
  fprintf(fp,">> Solver Control Parameters\n\n");
  
  fprintf(fp,"\tSolver Properties\n");
  
  // Precision
  if ( sizeof(REAL_TYPE) == sizeof(float) )
  {
    fprintf(fp,"\t     Precision                :   Single Precision \n");
  }
  else
  {
    fprintf(fp,"\t     Precision                :   Double Precision \n");
  }
  
  // Flow Algorithm
  fprintf(fp,"\t     Flow Algorithm           :   Fractional Step\n");
  fprintf(fp,"\t        Time marching scheme  :   Euler Explicit O(dt1)\n");

  // Convection scheme
  switch (CnvScheme)
  {
    case O1_upwind:
      fprintf(fp,"\t     Convective flux scheme   :   Upwind O(dx1)\n");
      break;
    case O3_muscl:
      fprintf(fp,"\t     Convective flux scheme   :   MUSCL O(dx3)\n");
      break;
    default:
      stamped_printf("Error: Convection scheme section\n");
      Exit(1);
  }

  // 時間制御 ------------------
  fprintf(fp,"\n\tTime Control\n");
  
  // 加速時間
  fprintf(fp,"\t     Acceleration Time        :   %12.5e [-]\n", TimeAccel);
  
  // Time Increment
  fprintf(fp,"\t     Time Increment dt        :   %12.5e [-]: CFL [%8.5f] with Reference velocity\n",  deltaT, CFL);
  
  // Calculation time/step
  fprintf(fp,"\t     Calculation Step         :   %12d\n", Session_LastStep);


  // Criteria ------------------
  fprintf(fp,"\n\tParameter of Linear Equation\n");
  const ItrCtl* ICp1= &IC[ItrCtl::ic_prs_pr];  /// 圧力のPoisson反復
  const ItrCtl* ICd = &IC[ItrCtl::ic_div];     /// V-P反復

  if ( PM_Test == ON)
  {
    fprintf(fp,"\t ### Performance Test Mode >> The iteration number is fixed by Iteration max.\n\n");
  }
  
  // V-P iteration
  fprintf(fp,"\t     V-P Iteration \n");
  fprintf(fp,"\t       Iteration max          :   %d\n"  ,  ICd->get_ItrMax());
  fprintf(fp,"\t       Convergence eps        :   %9.3e\n", ICd->get_eps());
  fprintf(fp,"\t       Norm type              :   %s\n", ItrCtl::getNormString(ICd->get_normType()).c_str() );
  
  // 1st iteration
  fprintf(fp,"\t     1st Pressure Iteration \n");
  fprintf(fp,"\t       Iteration max          :   %d\n"  ,  ICp1->get_ItrMax());
  fprintf(fp,"\t       Convergence eps        :   %9.3e\n", ICp1->get_eps());
  fprintf(fp,"\t       Coef. of Relax./Accel. :   %9.3e\n", ICp1->get_omg());
  fprintf(fp,"\t       Norm type              :   %s\n", ItrCtl::getNormString(ICp1->get_normType()).c_str() );
  fprintf(fp,"\t       Communication Mode     :   %s\n", (ICp1->get_SyncMode()==comm_sync) ? "SYNC" : "ASYNC");
  switch (ICp1->get_LS()) 
  {
    case SOR2SMA:
      fprintf(fp,"\t       Linear Solver          :   2-colored SOR SMA (Stride Memory Access)\n");
      break;
      
    default:
      stamped_printf("Error: Linear Solver section\n");
      Exit(0);
  }

  // Output
  fprintf(fp,"\n\tOutput\n");
  if (OutputInterval > 0) {
    fprintf(fp,"\t       Interval Steps         :   %d\n", OutputInterval);
  } else {
    fprintf(fp,"\t       Interval Steps         :   %d (no output)\n", OutputInterval);
  }

  fprintf(fp,"\n");
  fprintf(fp,">> Simulation Parameters\n\n");
  
  // Reference values
  fprintf(fp,"\tRef. Length               [m]         : %12.5e\n", RefLength);
  fprintf(fp,"\tRef. Velocity             [m/s]       : %12.5e\n", RefVelocity);
  fprintf(fp,"\n");
  
  fprintf(fp,"\tSpacing                   [m] / [-]   : %12.5e / %12.5e\n", deltaX*RefLength, deltaX);
  fprintf(fp,"\tTime Scale                [sec]       : %12.5e\n", Tscale);
  fprintf(fp,"\n");

  fprintf(fp,"\tReynolds number           [-]         : %12.5e\n", Reynolds);
  fprintf(fp,"\n");

  fflush(fp);

}


// #################################################################
// メモリ使用量の表示
void FFV::display_memory_info(FILE* fp, double L_mem, const char* str)
{
  double G_mem;

  if ( numProc > 1 )
  {
    if ( paraMngr->Allreduce(&L_mem, &G_mem, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
  }
  else 
  {
    G_mem = L_mem;
  }
  
  Hostonly_
  {
    FBUtility::MemoryRequirement(str, G_mem, L_mem, fp);
  }
}


// #################################################################
// 外部境界条件の各面の情報を表示する
void FFV::printFaceOBC(FILE* fp)
{
  BoundaryOuter* bc = BC.export_OBC();

  for (int i=0; i<NOFACE; i++) {
    REAL_TYPE c;
    switch ( bc[i].get_Class() ) {
      case OBC_WALL:
        fprintf(fp,"\tSet %s up as %s : < %s >\n", 
                FBUtility::getDirection(i).c_str(), 
                "Wall", bc[i].get_Alias().c_str());
        fprintf(fp,"\t\tGuide Cell Medium = %s\n",
                mat[bc[i].get_GuideMedium()].getLabel().c_str());
        c = bc[i].ca[CompoList::bias];
        fprintf(fp,"\t\tVelocity V(%10.3e, %10.3e, %10.3e) [-]\n", 
                bc[i].nv[0]*c, bc[i].nv[1]*c, bc[i].nv[2]*c);
        break;

      default:
        printf("\n\tError : OuterBC\n");
        Exit(0);
        break;
   
    }
    fprintf(fp,"\n");
  }
  fflush(fp);
}


// #################################################################
// 初期条件の設定
void FFV::setInitialCondition()
{
  double flop_task;
  
  REAL_TYPE tm = CurrentTime * Tscale;
  
  // 速度の初期条件の設定
  REAL_TYPE U0[3] = { 0.0, 0.0, 0.0 };
  fb_set_vector_(d_v, size, &guide, U0, d_bcd);
  
  
  // 外部境界面の移流速度を計算し，外部境界条件を設定
  BC.OuterVBC(d_v, d_v, d_bcv, tm, deltaT, v00, flop_task);
  
  // 圧力
  REAL_TYPE ip = 0.0;;
  FBUtility::xset(d_p, size, guide, ip, kind_scalar);
  BC.OuterPBC(d_p);
    
  
  // 初期解の同期
  if ( numProc > 1 )
  {
    if ( paraMngr->BndCommV3DEx(d_v, size[0], size[1], size[2], guide, guide) != CPM_SUCCESS ) Exit(0);
    if ( paraMngr->BndCommS3D  (d_p, size[0], size[1], size[2], guide, 1    ) != CPM_SUCCESS ) Exit(0);
  }
  
}


// #################################################################
// グローバルな領域情報を表示する
void FFV::printGlobalDomain(FILE* fp)
{
  REAL_TYPE KB = 1000.0;
  REAL_TYPE MB = 1000.0*KB;
  REAL_TYPE GB = 1000.0*MB;
  REAL_TYPE TB = 1000.0*GB;
  REAL_TYPE PB = 1000.0*TB;

  REAL_TYPE total = (REAL_TYPE)G_size[0] * (REAL_TYPE)G_size[1] * (REAL_TYPE)G_size[2];
  
  fprintf(fp,"\timax, jmax, kmax    = %13d %13d %13d\n\t\t     >> ", 
          G_size[0], G_size[1], G_size[2]);

  
  if ( total > PB ) {
    fprintf (fp,"%6.2f (P cells)\n", total / PB);
  }
  else if ( total > TB ) {
    fprintf (fp,"%6.2f (T cells)\n", total / TB);
  }
  else if ( total > GB ) {
    fprintf (fp,"%6.2f (G cells)\n", total / GB);
  }
  else if ( total > MB ) {
    fprintf (fp,"%6.2f (M cells)\n", total / MB);
  }
  else if ( total > KB ) {
    fprintf (fp,"%6.2f (K cells)\n", total / KB);
  }
  else if ( total <= KB ){
    fprintf (fp,"%6.2f (cells)\n", total);
  }
  fprintf(fp,"\n");
  
  fprintf(fp,"\t(dx, dy, dz) [-] = (%13.6e %13.6e %13.6e)\n", pitch[0], pitch[1], pitch[2]);
  fprintf(fp,"\t(ox, oy, oz) [-] = (%13.6e %13.6e %13.6e)\n", G_origin[0], G_origin[1], G_origin[2]);
  fprintf(fp,"\t(Lx, Ly, Lz) [-] = (%13.6e %13.6e %13.6e)\n", G_region[0], G_region[1], G_region[2]);
  fprintf(fp,"\n");
  
  fflush(fp);
}


// #################################################################
// 履歴の出力準備
void FFV::prep_HistoryOutput()
{
  // マスターノードでの履歴出力準備
  H = new History();

  Hostonly_ {
    H->printHistoryTitle(stdout, IC, true);
    
    // コンポーネント情報
    {
      // 基本情報　history.log, history_compo.log, history_domfx.log
      if ( !(fp_b=fopen("history_base.txt", "w")) ) 
      {
        stamped_printf("\tSorry, can't open '%s' file. Write failed.\n", "history_base.txt");
        Exit(0);
      }
      H->printHistoryTitle(fp_b, IC, false);

    }
    
  }
}


// #################################################################
// 並列分散時のファイル名の管理を行う
void FFV::setDFI()
{
  int* g_bbox_st = new int[3*numProc];
  int* g_bbox_ed = new int[3*numProc];
  
  // host nameの取得
  string host = paraMngr->GetHostName();
  if ( host.empty() ) Exit(0);
  
  
  // 並列時のみ
  if ( numProc > 1 )
  {
    const int* m_tail = paraMngr->GetVoxelTailIndex();
    int tail[3];
    
    tail[0] = m_tail[0] + 1;
    tail[1] = m_tail[1] + 1;
    tail[2] = m_tail[2] + 1;
    
    // 集約
    if ( paraMngr->Gather(head, 3, g_bbox_st, 3, 0) != CPM_SUCCESS ) Exit(0);
    if ( paraMngr->Gather(tail, 3, g_bbox_ed, 3, 0) != CPM_SUCCESS ) Exit(0);
    
  }
  else // serial
  {
    g_bbox_st[0] = 1;
    g_bbox_st[1] = 1;
    g_bbox_st[2] = 1;
    g_bbox_ed[0] = size[0];
    g_bbox_ed[1] = size[1];
    g_bbox_ed[2] = size[2];
  }
  
  // DFIクラスの初期化
  int GuideOut = 0;  // without Guide cells
  int Start = initial_start;
  if ( !DFI_.init(G_size, paraMngr->GetDivNum(), GuideOut, Start, g_bbox_st, g_bbox_ed, host) ) Exit(0);
  
  // 後始末
  delete [] g_bbox_st;
  delete [] g_bbox_ed;
  
}
