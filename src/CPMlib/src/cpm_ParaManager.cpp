/*
 * CPM - Cartesian Partition Manager
 *
 * Copyright (C) 2012 University of Tokyo.
 *
 */

/**
 * @file   cpm_ParaManager.cpp
 * パラレルマネージャクラスのソースファイル
 * @author University of Tokyo
 * @date   2012/05/31
 */
#include "cpm_ParaManager.h"

////////////////////////////////////////////////////////////////////////////////
// 唯一のインスタンスの取得
cpm_ParaManager*
cpm_ParaManager::get_instance()
{
  // 宣言
  static cpm_ParaManager instance;

  // ポインタ
  return &instance;
}

////////////////////////////////////////////////////////////////////////////////
// 唯一のインスタンスの取得(initialize処理も実行)
cpm_ParaManager*
cpm_ParaManager::get_instance(int &argc, char**& argv)
{
  // 宣言
  static cpm_ParaManager instance;

  // MPI_Init
  if( instance.Initialize(argc,argv) != CPM_SUCCESS )
  {
    return NULL;
  }

  // ポインタ
  return &instance;
}

////////////////////////////////////////////////////////////////////////////////
// コンストラクタ
cpm_ParaManager::cpm_ParaManager()
  : cpm_Base()
{
  // 並列数、ランク番号
  m_nRank  = 1;
  m_rankNo = 0;

  // プロセスグループリストをクリア
  m_procGrpList.clear();
  m_procGrpList.push_back(MPI_COMM_WORLD);

  // 領域管理情報マップのクリア
  m_voxelInfoMap.clear();
  m_rankNoMap.clear();

  // 袖通信バッファ情報のクリア
  m_bndCommInfoMap.clear();
}

////////////////////////////////////////////////////////////////////////////////
// デストラクタ
cpm_ParaManager::~cpm_ParaManager()
{
  // プロセスグループリストのクリア
  m_procGrpList.clear();

  // 領域管理情報マップの削除、クリア
  {
    VoxelInfoMap::iterator it  = m_voxelInfoMap.begin();
    VoxelInfoMap::iterator ite = m_voxelInfoMap.end();
    for( ; it!=ite; it++ )
    {
      if( it->second ) delete it->second;
    }
    m_voxelInfoMap.clear();
  }

  // ランク番号マップの削除、クリア
  {
    RankNoMap::iterator it  = m_rankNoMap.begin();
    RankNoMap::iterator ite = m_rankNoMap.end();
    for( ; it!=ite; it++ )
    {
      if( it->second ) delete [] it->second;
    }
    m_rankNoMap.clear();
  }

  // 袖通信バッファ情報の削除、クリア
  {
    BndCommInfoMap::iterator it  = m_bndCommInfoMap.begin();
    BndCommInfoMap::iterator ite = m_bndCommInfoMap.end();
    for( ; it!=ite; it++ )
    {
      if( it->second ) delete it->second;
    }
    m_bndCommInfoMap.clear();
  }

  // MPIのFinalize
  int flag1, flag2;
  MPI_Initialized(&flag1);
  MPI_Finalized(&flag2);
  if( flag1==true && flag2==false ) MPI_Finalize();
}

////////////////////////////////////////////////////////////////////////////////
// 初期化処理(MPI_Initは実行済みである必要がある)
cpm_ErrorCode
cpm_ParaManager::Initialize()
{
  m_nRank  = 1;
  m_rankNo = 0;

  // MPI_Init実行済みかチェック
  int flag1;
  MPI_Initialized(&flag1);
  if( flag1==false )
  {
    return CPM_ERROR_NO_MPI_INIT;
  }

  // プロセス並列数の取得
  if( MPI_Comm_size(MPI_COMM_WORLD, &m_nRank) != MPI_SUCCESS )
  {
    std::cerr << "MPI_Comm_size error." << std::endl;
    return CPM_ERROR_MPI;
  }
  if( m_nRank < 1 ) m_nRank = 1;

  // 自ランク番号の取得
  if( MPI_Comm_rank(MPI_COMM_WORLD, &m_rankNo) != MPI_SUCCESS )
  {
    std::cerr << "MPI_Comm_rank error." << std::endl;
    return CPM_ERROR_MPI;
  }
  if( m_rankNo < 0 ) m_rankNo = 0;

#ifdef _DEBUG
  Barrier();
  std::cout << std::flush;
  Barrier();
  for( int i=0;i<m_nRank;i++ )
  {
    if( i==m_rankNo )
      std::cout << "[" << m_rankNo << "] #rank=" << m_nRank << " parallel=" << IsParallel() << std::endl;
    std::cout << std::flush;
    Barrier();
  }
  std::cout << std::flush;
  Barrier();
#endif

  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// 初期化処理(MPI_Initも実行)
cpm_ErrorCode
cpm_ParaManager::Initialize( int &argc, char**& argv )
{
  m_nRank  = 1;
  m_rankNo = 0;

  // MPI_Init
  int flag1;
  MPI_Initialized(&flag1);
  if( flag1==false )
  {
    if( MPI_Init(&argc,&argv) != MPI_SUCCESS )
    {
      std::cerr << "MPI_Init error." << std::endl;
      return CPM_ERROR_MPI;
    }
  }

  // initialize
  return Initialize();
}

////////////////////////////////////////////////////////////////////////////////
// 並列実行であるかチェックする
bool
cpm_ParaManager::IsParallel()
{
  if( m_nRank <= 1 )
  {
    return false;
  }
  return true;
}

////////////////////////////////////////////////////////////////////////////////
// 並列実行であるかチェックする(const)
bool
cpm_ParaManager::IsParallel() const
{
  if( m_nRank <= 1 )
  {
    return false;
  }
  return true;
}

////////////////////////////////////////////////////////////////////////////////
// 領域分割
cpm_ErrorCode
cpm_ParaManager::VoxelInit( cpm_GlobalDomainInfo* domainInfo, size_t maxVC, size_t maxN, int procGrpNo )
{
  cpm_ErrorCode ret;
  cpm_VoxelInfo *voxelInfo = NULL;

  // 入力引数のチェック
  if( !domainInfo )
  {
    return CPM_ERROR_INVALID_PTR;
  }

  // 既に領域分割済みか
  if( m_voxelInfoMap.find(procGrpNo) != m_voxelInfoMap.end() )
  {
    return CPM_ERROR_ALREADY_VOXELINIIT;
  }

  // コミュニケータを取得
  MPI_Comm comm = GetMPI_Comm( procGrpNo );
  if( IsCommNull( comm ) )
  {
    return CPM_ERROR_MPI_INVALID_COMM;
  }

  // commの並列数を取得
  int nRank = 0;
  MPI_Comm_size( comm, &nRank );

  // 領域情報のチェック
  // 活性サブドメイン配列が空のとき、全領域が活性サブドメインになる
  // その場合、CheckData関数内で活性サブドメイン情報を生成する
  if( (ret = domainInfo->CheckData( nRank )) != CPM_SUCCESS )
  {
    return ret;
  }

  // プロセス数と活性サブドメイン数のチェック
  if( nRank != domainInfo->GetSubdomainNum() )
  {
    Abort(CPM_ERROR_MISMATCH_NP_SUBDOMAIN);
    return CPM_ERROR_MISMATCH_NP_SUBDOMAIN;
  }

  // インスタンス
  voxelInfo = new cpm_VoxelInfo();
  if( !voxelInfo )
  {
    Abort(CPM_ERROR_INVALID_PTR);
    return CPM_ERROR_INVALID_PTR;
  }

  // 領域分割情報の生成
  if( (ret = voxelInfo->Init( comm, domainInfo )) != CPM_SUCCESS )
  {
    delete voxelInfo;
    Abort(ret);
    return ret;
  }

  // VOXEL空間マップに登録
  if( !m_voxelInfoMap.insert(std::make_pair(procGrpNo, voxelInfo)).second )
  {
    m_procGrpList.pop_back();
    delete voxelInfo;
    return CPM_ERROR_INSERT_VOXELMAP;
  }

  // 袖通信バッファの設定
  SetBndCommBuffer( maxVC, maxN, procGrpNo );

  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// 領域分割(領域分割数を指定)
cpm_ErrorCode
cpm_ParaManager::VoxelInit( int div[3], int vox[3], double origin[3], double region[3]
                          , size_t maxVC, size_t maxN, int procGrpNo )
{
  // 入力値のチェック
  if( vox[0] <= 0 || vox[1] <= 0 || vox[2] <= 0 )
  {
    return CPM_ERROR_INVALID_VOXELSIZE;
  }
  if( region[0] <= 0.0 || region[1] <= 0.0 || region[2] <= 0.0 )
  {
    return CPM_ERROR_INVALID_REGION;
  }

  // 領域分割数
  if( div[0]<=0 || div[1]<=0 || div[2]<=0 )
  {
    // ランク数
    int nrank = GetNumRank( procGrpNo );

    // 領域分割数の決定
    div[0] = div[1] = div[2] = 0;
    cpm_ErrorCode ret;
    if( (ret = DecideDivPattern( nrank, vox, div )) != CPM_SUCCESS )
    {
      return ret;
    }
  }

  //ピッチを計算
  double pitch[3];
  for( int i=0;i<3;i++ )
  {
    pitch[i] = region[i] / double(vox[i]);
  }

  // DomainInfoを生成
  cpm_GlobalDomainInfo dInfo;
  dInfo.SetOrigin( origin );
  dInfo.SetPitch ( pitch );
  dInfo.SetRegion( region );
  dInfo.SetVoxNum( vox );
  dInfo.SetDivNum( div );

  // 共通の処理
  return VoxelInit( &dInfo, maxVC, maxN, procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
// 領域分割(プロセスグループのランク数で自動領域分割)
cpm_ErrorCode
cpm_ParaManager::VoxelInit( int vox[3], double origin[3], double pitch[3]
                          , size_t maxVC, size_t maxN, int procGrpNo )
{
  // 領域分割
  int div[3] = {0,0,0};
  return VoxelInit( div, vox, origin, pitch, maxVC, maxN, procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
// 並列プロセス数からI,J,K方向の分割数を取得する
cpm_ErrorCode
cpm_ParaManager::DecideDivPattern( int divNum
                                 , int voxSize[3]
                                 , int divPttn[3]
                                 ) const
{
  if( !voxSize || !divPttn )
  {
    return CPM_ERROR_INVALID_PTR;
  }
  if( (voxSize[0]==0) || (voxSize[1]==0) || (voxSize[2]==0) )
  {
    return CPM_ERROR_INVALID_VOXELSIZE;
  }
  if( divNum <= 1 ){
    divPttn[0] = divPttn[1] = divPttn[2] = 1;
    return CPM_SUCCESS;
  }

  divPttn[0] = divPttn[1] = divPttn[2] = 0;

  unsigned long long minCommSize = 0;

  unsigned long long divNumll = divNum;
  unsigned long long voxSizell[3] = {0, 0, 0};
  unsigned long long divPttnll[3] = {0, 0, 0};
  voxSizell[0] = voxSize[0];
  voxSizell[1] = voxSize[1];
  voxSizell[2] = voxSize[2];

#if 0 //2012.08.03.mod.start
  bool flag = false;
  unsigned long long n;
  for(n=1; n<=divNumll; n++){
    unsigned long long i, j, k;
    for(i=1; i<=n; i++){
      if( (voxSizell[0] / i) < 1 ) break;
      unsigned long long jmax = n/i;
      for(j=1; j<=jmax; j++){
        if( (voxSizell[1] / j) < 1 ) break;
        unsigned long long kmax = jmax/j;
        for(k=1; k<=kmax; k++){
          if( (voxSizell[2] / k) < 1 ) break;
          if( i*j*k == divNumll){
            unsigned long long commSize;
            if( (commSize=CalcCommSize(i, j, k, voxSizell)) == 0 ) break;

            if( !flag ) {
              minCommSize = commSize;
              flag = true;
            }
            else if( commSize < minCommSize ){
              divPttnll[0] = i; divPttnll[1] = j; divPttnll[2] = k;
              minCommSize = commSize;
            }
          }
        }
      }
    }
  }
#else
  bool flag = false;
  unsigned long long i, j, k;
  for(i=1; i<=divNumll; i++)
  {
    if( divNumll%i != 0 ) continue;
    unsigned long long jmax = divNumll/i;
    for(j=1; j<=jmax; j++)
    {
      if( (divNumll/i)%j != 0 ) continue;
      unsigned long long kmax = divNumll/(i*j);
      for(k=1; k<=kmax; k++)
      {
        if( (divNumll/(i*j))%k != 0 ) continue;
        if( (voxSizell[2] / k) < 1 ) break;
        if( i*j*k == divNumll){
          unsigned long long commSize;
          if( (commSize=CalcCommSize(i, j, k, voxSizell)) == 0 ) break;

          if( !flag )
          {
            minCommSize = commSize;
            flag = true;
          }
          else if( commSize < minCommSize )
          {
            divPttnll[0] = i; divPttnll[1] = j; divPttnll[2] = k;
            minCommSize = commSize;
          }
        }
      }
    }
  }
#endif //2012.08.03.mod.end

  divPttn[0] = divPttnll[0];
  divPttn[1] = divPttnll[1];
  divPttn[2] = divPttnll[2];

#if 0
  if( (divPttn[0]==0) || (divPttn[1]==0) || (divPttn[2]==0) )
  {
    divPttn[0] = divPttn[1] = divPttn[2] = 1;
    int maxLenElem = 2;
    if( voxSize[maxLenElem] < voxSize[1] ) maxLenElem = 1;
    if( voxSize[maxLenElem] < voxSize[0] ) maxLenElem = 0;
    divPttn[maxLenElem] = divNum;
  }
#else
  // 決定しなかったとき、Dims_createを試す
  if( (divPttn[0]==0) || (divPttn[1]==0) || (divPttn[2]==0) )
  {
    if( MPI_Dims_create( divNum, 3, divPttn ) != MPI_SUCCESS )
    {
      return CPM_ERROR_MPI_DIMSCREATE;
    }
  }
#endif

  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// I,J,K分割を行った時の通信点数の総数を取得する
unsigned long long
cpm_ParaManager::CalcCommSize( unsigned long long iDiv
                             , unsigned long long jDiv
                             , unsigned long long kDiv
                             , unsigned long long voxSize[3]
                             ) const
{
  if( (iDiv==0) || (jDiv==0) || (kDiv==0) ) return 0;
  if( !voxSize ) return 0;

  unsigned long long Len[3];
  Len[0] = voxSize[0] / iDiv; if( Len[0] == 0 ) return 0;
  Len[1] = voxSize[1] / jDiv; if( Len[1] == 0 ) return 0;
  Len[2] = voxSize[2] / kDiv; if( Len[2] == 0 ) return 0;

  unsigned long long commFace[3];
  if( iDiv != 1) commFace[0] = Len[1]*Len[2]*iDiv;
  else commFace[0] = 0;
  if( jDiv != 1) commFace[1] = Len[2]*Len[0]*jDiv;
  else commFace[1] = 0;
  if( kDiv != 1) commFace[2] = Len[0]*Len[1]*kDiv;
  else commFace[2] = 0;

  return (commFace[0] + commFace[1] + commFace[2]);
}

////////////////////////////////////////////////////////////////////////////////
// プロセスグループの作成
int
cpm_ParaManager::CreateProcessGroup( int nproc, int *proclist, int parentProcGrpNo )
{
  // 入力チェック
  if( parentProcGrpNo < 0 || nproc < 1 || !proclist )
  {
    return -1;
  }

  // 親のMPIコミュニケータを取得
  MPI_Comm parentComm = GetMPI_Comm( parentProcGrpNo );
  if( IsCommNull( parentComm ) )
  {
    return -1;
  }

  // 親のグループを取得
  MPI_Group parentGrp;
  MPI_Comm_group(parentComm, &parentGrp);

  // 新しいグループを生成
  MPI_Group newGrp;
  MPI_Group_incl(parentGrp, nproc, proclist, &newGrp);

  // 新しいコミュニケータを生成
  MPI_Comm  newComm;
  MPI_Comm_create(parentComm, newGrp, &newComm);
  if( IsCommNull(newComm) )
  {
    // 時ランクが含まれない
    return -1;
  }

  // コミュニケータを登録
  m_procGrpList.push_back(newComm);

  // プロセスグループの番号
  return m_procGrpList.size()-1;
}

////////////////////////////////////////////////////////////////////////////////
// VOXEL空間マップを検索
const cpm_VoxelInfo*
cpm_ParaManager::FindVoxelInfo( int procGrpNo )
{
  VoxelInfoMap::iterator it = m_voxelInfoMap.find(procGrpNo);
  if( it == m_voxelInfoMap.end() ) return NULL;
  return it->second;
}

////////////////////////////////////////////////////////////////////////////////
// 領域分割数を取得
const int*
cpm_ParaManager::GetDivNum( int procGrpNo )
{
  //VOXEL空間マップを検索
  const cpm_VoxelInfo *pVoxelInfo = FindVoxelInfo( procGrpNo );
  if( !pVoxelInfo ) return NULL;

  return pVoxelInfo->GetDivNum();
}

////////////////////////////////////////////////////////////////////////////////
// ピッチを取得
const double*
cpm_ParaManager::GetPitch( int procGrpNo )
{
  //VOXEL空間マップを検索
  const cpm_VoxelInfo *pVoxelInfo = FindVoxelInfo( procGrpNo );
  if( !pVoxelInfo ) return NULL;

  return pVoxelInfo->GetPitch();
}

////////////////////////////////////////////////////////////////////////////////
// 全体ボクセル数を取得
const int*
cpm_ParaManager::GetGlobalVoxelSize( int procGrpNo )
{
  //VOXEL空間マップを検索
  const cpm_VoxelInfo *pVoxelInfo = FindVoxelInfo( procGrpNo );
  if( !pVoxelInfo ) return NULL;

  return pVoxelInfo->GetGlobalVoxelSize();
}

////////////////////////////////////////////////////////////////////////////////
// 全体空間の原点を取得
const double*
cpm_ParaManager::GetGlobalOrigin( int procGrpNo )
{
  //VOXEL空間マップを検索
  const cpm_VoxelInfo *pVoxelInfo = FindVoxelInfo( procGrpNo );
  if( !pVoxelInfo ) return NULL;

  return pVoxelInfo->GetGlobalOrigin();
}

////////////////////////////////////////////////////////////////////////////////
// 全体空間サイズを取得
const double*
cpm_ParaManager::GetGlobalRegion( int procGrpNo )
{
  //VOXEL空間マップを検索
  const cpm_VoxelInfo *pVoxelInfo = FindVoxelInfo( procGrpNo );
  if( !pVoxelInfo ) return NULL;

  return pVoxelInfo->GetGlobalRegion();
}

////////////////////////////////////////////////////////////////////////////////
// 自ランクのボクセル数を取得
const int*
cpm_ParaManager::GetLocalVoxelSize( int procGrpNo )
{
  //VOXEL空間マップを検索
  const cpm_VoxelInfo *pVoxelInfo = FindVoxelInfo( procGrpNo );
  if( !pVoxelInfo ) return NULL;

  return pVoxelInfo->GetLocalVoxelSize();
}

////////////////////////////////////////////////////////////////////////////////
// 自ランクの空間原点を取得
const double*
cpm_ParaManager::GetLocalOrigin( int procGrpNo )
{
  //VOXEL空間マップを検索
  const cpm_VoxelInfo *pVoxelInfo = FindVoxelInfo( procGrpNo );
  if( !pVoxelInfo ) return NULL;

  return pVoxelInfo->GetLocalOrigin();
}

////////////////////////////////////////////////////////////////////////////////
// 自ランクの空間サイズを取得
const double*
cpm_ParaManager::GetLocalRegion( int procGrpNo )
{
  //VOXEL空間マップを検索
  const cpm_VoxelInfo *pVoxelInfo = FindVoxelInfo( procGrpNo );
  if( !pVoxelInfo ) return NULL;

  return pVoxelInfo->GetLocalRegion();
}

////////////////////////////////////////////////////////////////////////////////
// 自ランクの領域分割位置を取得
const int*
cpm_ParaManager::GetDivPos( int procGrpNo )
{
  //VOXEL空間マップを検索
  const cpm_VoxelInfo *pVoxelInfo = FindVoxelInfo( procGrpNo );
  if( !pVoxelInfo ) return NULL;

  return pVoxelInfo->GetDivPos();
}

////////////////////////////////////////////////////////////////////////////////
// 自ランクの始点VOXELの全体空間でのインデクスを取得
const int*
cpm_ParaManager::GetVoxelHeadIndex( int procGrpNo )
{
  //VOXEL空間マップを検索
  const cpm_VoxelInfo *pVoxelInfo = FindVoxelInfo( procGrpNo );
  if( !pVoxelInfo ) return NULL;

  return pVoxelInfo->GetVoxelHeadIndex();
}

////////////////////////////////////////////////////////////////////////////////
// 自ランクの終点VOXELの全体空間でのインデクスを取得
const int*
cpm_ParaManager::GetVoxelTailIndex( int procGrpNo )
{
  //VOXEL空間マップを検索
  const cpm_VoxelInfo *pVoxelInfo = FindVoxelInfo( procGrpNo );
  if( !pVoxelInfo ) return NULL;

  return pVoxelInfo->GetVoxelTailIndex();
}

////////////////////////////////////////////////////////////////////////////////
// 自ランクの隣接ランク番号を取得
const int*
cpm_ParaManager::GetNeighborRankID( int procGrpNo )
{
  //VOXEL空間マップを検索
  const cpm_VoxelInfo *pVoxelInfo = FindVoxelInfo( procGrpNo );
  if( !pVoxelInfo ) return NULL;

  return pVoxelInfo->GetNeighborRankID();
}

////////////////////////////////////////////////////////////////////////////////
// 自ランクの周期境界の隣接ランク番号を取得
const int*
cpm_ParaManager::GetPeriodicRankID( int procGrpNo )
{
  //VOXEL空間マップを検索
  const cpm_VoxelInfo *pVoxelInfo = FindVoxelInfo( procGrpNo );
  if( !pVoxelInfo ) return NULL;

  return pVoxelInfo->GetPeriodicRankID();
}

////////////////////////////////////////////////////////////////////////////////
// 袖通信バッファのセット
cpm_ErrorCode
cpm_ParaManager::SetBndCommBuffer( size_t maxVC, size_t maxN, int procGrpNo )
{
  if( maxVC==0 || maxN==0 )
  {
    return CPM_ERROR_BNDCOMM;
  }

  // local voxel size
  const int *sz = GetLocalVoxelSize(procGrpNo);
  if( !sz ) return CPM_ERROR_BNDCOMM_VOXELSIZE;

  // buffer size
  size_t nwX = size_t(sz[1]+2*maxVC) * size_t(sz[2]+2*maxVC) * maxVC * maxN;
  size_t nwY = size_t(sz[2]+2*maxVC) * size_t(sz[0]+2*maxVC) * maxVC * maxN;
  size_t nwZ = size_t(sz[0]+2*maxVC) * size_t(sz[1]+2*maxVC) * maxVC * maxN;

  // 送受信バッファ情報のインスタンス
  S_BNDCOMM_BUFFER *bufInfo = new S_BNDCOMM_BUFFER();
  if( !bufInfo )
  {
    return CPM_ERROR_INVALID_PTR;
  }

  bufInfo->m_maxVC = maxVC;
  bufInfo->m_maxN  = maxN;
  bufInfo->m_nwX   = nwX;
  bufInfo->m_nwY   = nwY;
  bufInfo->m_nwZ   = nwZ;

  // buffer
  for( int i=0;i<4;i++ )
  {
    bufInfo->m_bufX[i] = new REAL_BUF_TYPE[nwX];
    if( !bufInfo->m_bufX[i] )
    {
      delete bufInfo;
      return CPM_ERROR_INVALID_PTR;
    }
    bufInfo->m_bufY[i] = new REAL_BUF_TYPE[nwY];
    if( !bufInfo->m_bufY[i] )
    {
      delete bufInfo;
      return CPM_ERROR_INVALID_PTR;
    }
    bufInfo->m_bufZ[i] = new REAL_BUF_TYPE[nwZ];
    if( !bufInfo->m_bufZ[i] )
    {
      delete bufInfo;
      return CPM_ERROR_INVALID_PTR;
    }
  }

  //マップにセット
  BndCommInfoMap::iterator it = m_bndCommInfoMap.find(procGrpNo);
  if( it != m_bndCommInfoMap.end() )
  {
    delete it->second;
    m_bndCommInfoMap.erase( it );
  }
  m_bndCommInfoMap.insert( std::make_pair(procGrpNo, bufInfo) );

  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// 袖通信バッファサイズの取得
size_t
cpm_ParaManager::GetBndCommBufferSize( int procGrpNo )
{
  size_t mem = 0;

  // マップを検索
  if( procGrpNo < 0 )
  {
    //全体
    BndCommInfoMap::iterator it = m_bndCommInfoMap.begin();
    for( ; it!=m_bndCommInfoMap.end(); it++ )
    {
      if( !it->second ) continue;
      size_t nwX   = it->second->m_nwX;
      size_t nwY   = it->second->m_nwY;
      size_t nwZ   = it->second->m_nwZ;
      mem += ((nwX*4 + nwY*4 + nwZ*4) * sizeof(REAL_BUF_TYPE));
    }
  }
  else
  {
    // 指定プロセスグループを検索
    S_BNDCOMM_BUFFER *bInfo = GetBndCommBuffer( procGrpNo );
    if( !bInfo ) return 0;
    size_t nwX   = bInfo->m_nwX;
    size_t nwY   = bInfo->m_nwY;
    size_t nwZ   = bInfo->m_nwZ;
    mem = (nwX*4 + nwY*4 + nwZ*4) * sizeof(REAL_BUF_TYPE);
  }

  return mem;
}

