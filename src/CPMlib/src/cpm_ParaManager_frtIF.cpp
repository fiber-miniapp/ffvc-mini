/*
 * CPM - Cartesian Partition Manager
 *
 * Copyright (C) 2012 University of Tokyo.
 *
 */

/**
 * @file   cpm_ParaManager_frtIF.cpp
 * パラレルマネージャクラスのFortranインターフェイスのソースファイル
 * @author University of Tokyo
 * @date   2012/05/31
 */
#include "cpm_ParaManager.h"

/** extern宣言 */
#define CPM_EXTERN extern "C"

#if 0
/** S3D,V3D通信でS4D版を使う */
#define _USE_S4D_
#endif

#ifndef CPM_WINDOWS
  #define cpm_Initialize_            cpm_initialize_
  #define cpm_VoxelInit_             cpm_voxelinit_
  #define cpm_VoxelInit_nodiv_       cpm_voxelinit_nodiv_
  #define cpm_IsParallel_            cpm_isparallel_
  #define cpm_GetDivNum_             cpm_getdivnum_
  #define cpm_GetPitch_              cpm_getpitch_
  #define cpm_GetGlobalVoxelSize_    cpm_getglobalvoxelsize_
  #define cpm_GetGlobalOrigin_       cpm_getglobalorigin_
  #define cpm_GetGlobalRegion_       cpm_getglobalregion_
  #define cpm_GetLocalVoxelSize_     cpm_getlocalvoxelsize_
  #define cpm_GetLocalOrigin_        cpm_getlocalorigin_
  #define cpm_GetLocalRegion_        cpm_getlocalregion_
  #define cpm_GetDivPos_             cpm_getdivpos_
  #define cpm_GetVoxelHeadIndex_     cpm_getvoxelheadindex_
  #define cpm_GetVoxelTailIndex_     cpm_getvoxeltailindex_
  #define cpm_GetNeighborRankID_     cpm_getneighborrankid_
  #define cpm_GetPeriodicRankID_     cpm_getperiodicrankid_
  #define cpm_GetMyRankID_           cpm_getmyrankid_
  #define cpm_GetNumRank_            cpm_getnumrank_
  #define cpm_Abort_                 cpm_abort_
  #define cpm_Barrier_               cpm_barrier_
  #define cpm_Wait_                  cpm_wait_
  #define cpm_Waitall_               cpm_waitall_
  #define cpm_Bcast_                 cpm_bcast_
  #define cpm_Send_                  cpm_send_
  #define cpm_Recv_                  cpm_recv_
  #define cpm_Isend_                 cpm_isend_
  #define cpm_Irecv_                 cpm_irecv_
  #define cpm_Allreduce_             cpm_allreduce_
  #define cpm_Gather_                cpm_gather_
  #define cpm_Allgather_             cpm_allgather_
  #define cpm_Gatherv_               cpm_gatherv_
  #define cpm_Allgatherv_            cpm_allgatherv_
  #define cpm_SetBndCommBuffer_      cpm_setbndcommbuffer_
  #define cpm_BndCommS3D_            cpm_bndcomms3d_
  #define cpm_BndCommV3D_            cpm_bndcommv3d_
  #define cpm_BndCommS4D_            cpm_bndcomms4d_
  #define cpm_BndCommS3D_nowait_     cpm_bndcomms3d_nowait_
  #define cpm_BndCommV3D_nowait_     cpm_bndcommv3d_nowait_
  #define cpm_BndCommS4D_nowait_     cpm_bndcomms4d_nowait_
  #define cpm_wait_BndCommS3D_       cpm_wait_bndcomms3d_
  #define cpm_wait_BndCommV3D_       cpm_wait_bndcommv3d_
  #define cpm_wait_BndCommS4D_       cpm_wait_bndcomms4d_
  #define cpm_BndCommV3DEx_          cpm_bndcommv3dex_
  #define cpm_BndCommS4DEx_          cpm_bndcomms4dex_
  #define cpm_BndCommV3DEx_nowait_   cpm_bndcommv3dex_nowait_
  #define cpm_BndCommS4DEx_nowait_   cpm_bndcomms4dex_nowait_
  #define cpm_wait_BndCommV3DEx_     cpm_wait_bndcommv3dex_
  #define cpm_wait_BndCommS4DEx_     cpm_wait_bndcomms4dex_
#else
  #define cpm_Initialize_            CPM_INITIALIZE
  #define cpm_VoxelInit_             CPM_VOXELINIT
  #define cpm_VoxelInit_nodiv_       CPM_VOXELINIT_NODIV
  #define cpm_IsParallel_            CPM_ISPARALLEL
  #define cpm_GetDivNum_             CPM_GETDIVNUM
  #define cpm_GetPitch_              CPM_GETPITCH
  #define cpm_GetGlobalVoxelSize_    CPM_GETGLOBALVOXELSIZE
  #define cpm_GetGlobalOrigin_       CPM_GETGLOBALORIGIN
  #define cpm_GetGlobalRegion_       CPM_GETGLOBALREGION
  #define cpm_GetLocalVoxelSize_     CPM_GETLOCALVOXELSIZE
  #define cpm_GetLocalOrigin_        CPM_GETLOCALORIGIN
  #define cpm_GetLocalRegion_        CPM_GETLOCALREGION
  #define cpm_GetDivPos_             CPM_GETDIVPOS
  #define cpm_GetVoxelHeadIndex_     CPM_GETVOXELHEADINDEX
  #define cpm_GetVoxelTailIndex_     CPM_GETVOXELTAILINDEX
  #define cpm_GetNeighborRankID_     CPM_GETNEIGHBORRANKID
  #define cpm_GetPeriodicRankID_     CPM_GETPERIODICRANKID
  #define cpm_GetMyRankID_           CPM_GETMYRANKID
  #define cpm_GetNumRank_            CPM_GETNUMRANK
  #define cpm_Abort_                 CPM_ABORT
  #define cpm_Barrier_               CPM_BARRIER
  #define cpm_Wait_                  CPM_WAIT
  #define cpm_Waitall_               CPM_WAITALL
  #define cpm_Bcast_                 CPM_BCAST
  #define cpm_Send_                  CPM_SEND
  #define cpm_Recv_                  CPM_RECV
  #define cpm_Isend_                 CPM_ISEND
  #define cpm_Irecv_                 CPM_IRECV
  #define cpm_Allreduce_             CPM_ALLREDUCE
  #define cpm_Gather_                CPM_GATHER
  #define cpm_Allgather_             CPM_ALLGATHER
  #define cpm_Gatherv_               CPM_GATHERV
  #define cpm_Allgatherv_            CPM_ALLGATHERV
  #define cpm_SetBndCommBuffer_      CPM_SETBNDCOMMBUFFER
  #define cpm_BndCommS3D_            CPM_BNDCOMMS3D
  #define cpm_BndCommV3D_            CPM_BNDCOMMV3D
  #define cpm_BndCommS4D_            CPM_BNDCOMMS4D
  #define cpm_BndCommS3D_nowait_     CPM_BNDCOMMS3D_NOWAIT
  #define cpm_BndCommV3D_nowait_     CPM_BNDCOMMV3D_NOWAIT
  #define cpm_BndCommS4D_nowait_     CPM_BNDCOMMS4D_NOWAIT
  #define cpm_wait_BndCommS3D_       CPM_WAIT_BNDCOMMS3D
  #define cpm_wait_BndCommV3D_       CPM_WAIT_BNDCOMMV3D
  #define cpm_wait_BndCommS4D_       CPM_WAIT_BNDCOMMS4D
  #define cpm_BndCommV3DEx_          CPM_BNDCOMMV3DEX
  #define cpm_BndCommS4DEx_          CPM_BNDCOMMS4DEX
  #define cpm_BndCommV3DEx_nowait_   CPM_BNDCOMMV3DEX_NOWAIT
  #define cpm_BndCommS4DEx_nowait_   CPM_BNDCOMMS4DEX_NOWAIT
  #define cpm_wait_BndCommV3DEx_     CPM_WAIT_BNDCOMMV3DEX
  #define cpm_wait_BndCommS4DEx_     CPM_WAIT_BNDCOMMS4DEX
#endif

////////////////////////////////////////////////////////////////////////////////
/** Abort
 *  - AbortのFortranインターフェイス関数
 *  @param[in]  errorcode MPI_Abortに渡すエラーコード
 */
CPM_EXTERN
void
cpm_Abort_( int *errorcode )
{
  int err = 0;
  if( errorcode )
  {
    err = *errorcode;
  }

  // インスタンス取得
  cpm_ParaManager *paraMngr = cpm_ParaManager::get_instance();
  if( !paraMngr )
  {
    MPI_Abort( MPI_COMM_WORLD, err );
    exit(err);
    return;
  }

  paraMngr->Abort( err );
}

////////////////////////////////////////////////////////////////////////////////
/** Barrier
 *  - BarrierのFortranインターフェイス関数
 *  @param[in]  procGrpNo プロセスグループ番号
 *  @param[out] ierr      終了コード(0=正常終了、0以外=cpm_ErrorCodeの値)
 */
CPM_EXTERN
void
cpm_Barrier_( int *procGrpNo, int *ierr )
{
  if( !procGrpNo || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  // インスタンス取得
  cpm_ParaManager *paraMngr = cpm_ParaManager::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  // Barrier
  *ierr = paraMngr->Barrier( *procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
/** Wait
 *  - WaitのFortranインターフェイス関数
 *  @param[in]  reqNo     リクエスト番号(0以上の整数)
 *  @param[out] ierr      終了コード(0=正常終了、0以外=cpm_ErrorCodeの値)
 */
CPM_EXTERN
void
cpm_Wait_( int *reqNo, int *ierr )
{
  if( !reqNo || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  // インスタンス取得
  cpm_ParaManager *paraMngr = cpm_ParaManager::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  // cpm_Wait
  *ierr = paraMngr->cpm_Wait( *reqNo );
}

////////////////////////////////////////////////////////////////////////////////
/** Waitall
 *  - WaitallのFortranインターフェイス関数
 *  @param[in]  count     リクエストの数
 *  @param[in]  reqlist   リクエスト番号のリスト(0以上の整数)
 *  @param[out] ierr      終了コード(0=正常終了、0以外=cpm_ErrorCodeの値)
 */
CPM_EXTERN
void
cpm_Waitall_( int *count, int *reqlist, int *ierr )
{
  if( !count || !reqlist || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  // インスタンス取得
  cpm_ParaManager *paraMngr = cpm_ParaManager::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  // cpm_Waitall
  *ierr = paraMngr->cpm_Waitall( *count, reqlist );
}

////////////////////////////////////////////////////////////////////////////////
/** Bcast
 *  - BcastのFortranインターフェイス関数
 *  @param[inout] buf       送受信バッファ
 *  @param[in]    count     送信バッファのサイズ(ワード数)
 *  @param[in]    datatype  データタイプ(fparam.fiを参照)
 *  @param[in]    root      送信元のランク番号(procGrpNo内でのランク番号)
 *  @param[in]    procGrpNo プロセスグループ番号
 *  @param[out]   ierr      終了コード(0=正常終了、0以外=cpm_ErrorCodeの値)
 */
CPM_EXTERN
void
cpm_Bcast_( void *buf, int *count, int *datatype, int *root, int *procGrpNo, int *ierr )
{
  if( !buf || !count || !datatype || !root || !procGrpNo || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  // インスタンス取得
  cpm_ParaManager *paraMngr = cpm_ParaManager::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  // MPI_Datatype
  MPI_Datatype dtype = cpm_ParaManager::GetMPI_Datatype( *datatype );
  if( dtype == MPI_DATATYPE_NULL )
  {
    *ierr = CPM_ERROR_MPI_INVALID_DATATYPE;
    return;
  }

  // Bcast
  *ierr = paraMngr->Bcast( dtype, buf, *count, *root, *procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
/** Send
 *  - SendのFortranインターフェイス関数
 *  @param[inout] buf       送信バッファ
 *  @param[in]    count     送信バッファのサイズ(ワード数)
 *  @param[in]    datatype  データタイプ(fparam.fiを参照)
 *  @param[in]    dest      送信先のランク番号(procGrpNo内でのランク番号)
 *  @param[in]    procGrpNo プロセスグループ番号
 *  @param[out]   ierr      終了コード(0=正常終了、0以外=cpm_ErrorCodeの値)
 */
CPM_EXTERN
void
cpm_Send_( void *buf, int *count, int *datatype, int *dest, int *procGrpNo, int *ierr )
{
  if( !buf || !count || !datatype || !dest || !procGrpNo || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  // インスタンス取得
  cpm_ParaManager *paraMngr = cpm_ParaManager::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  // MPI_Datatype
  MPI_Datatype dtype = cpm_ParaManager::GetMPI_Datatype( *datatype );
  if( dtype == MPI_DATATYPE_NULL )
  {
    *ierr = CPM_ERROR_MPI_INVALID_DATATYPE;
    return;
  }

  // Send
  *ierr = paraMngr->Send( dtype, buf, *count, *dest, *procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
/** Recv
 *  - RecvのFortranインターフェイス関数
 *  @param[inout] buf       受信バッファ
 *  @param[in]    count     受信バッファのサイズ(ワード数)
 *  @param[in]    datatype  データタイプ(fparam.fiを参照)
 *  @param[in]    source    送信元のランク番号(procGrpNo内でのランク番号)
 *  @param[in]    procGrpNo プロセスグループ番号
 *  @param[out]   ierr      終了コード(0=正常終了、0以外=cpm_ErrorCodeの値)
 */
CPM_EXTERN
void
cpm_Recv_( void *buf, int *count, int *datatype, int *source, int *procGrpNo, int *ierr )
{
  if( !buf || !count || !datatype || !source || !procGrpNo || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  // インスタンス取得
  cpm_ParaManager *paraMngr = cpm_ParaManager::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  // MPI_Datatype
  MPI_Datatype dtype = cpm_ParaManager::GetMPI_Datatype( *datatype );
  if( dtype == MPI_DATATYPE_NULL )
  {
    *ierr = CPM_ERROR_MPI_INVALID_DATATYPE;
    return;
  }

  // Recv
  *ierr = paraMngr->Recv( dtype, buf, *count, *source, *procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
/** Isend
 *  - IsendのFortranインターフェイス関数
 *  @param[inout] buf       送信バッファ
 *  @param[in]    count     送信バッファのサイズ(ワード数)
 *  @param[in]    datatype  データタイプ(fparam.fiを参照)
 *  @param[in]    dest      送信先のランク番号(procGrpNo内でのランク番号)
 *  @param[in]    procGrpNo プロセスグループ番号
 *  @param[out]   reqNo     リクエスト番号
 *  @param[out]   ierr      終了コード(0=正常終了、0以外=cpm_ErrorCodeの値)
 */
CPM_EXTERN
void
cpm_Isend_( void *buf, int *count, int *datatype, int *dest, int *procGrpNo, int *reqNo, int *ierr )
{
  if( !buf || !count || !datatype || !dest || !procGrpNo || !reqNo || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  // インスタンス取得
  cpm_ParaManager *paraMngr = cpm_ParaManager::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  // cpm_Isend
  if( (*ierr = paraMngr->cpm_Isend( buf, *count, *datatype, *dest, reqNo, *procGrpNo ) ) != CPM_SUCCESS )
  {
    return;
  }

  *ierr = CPM_SUCCESS;
  return;
}

////////////////////////////////////////////////////////////////////////////////
/** Irecv
 *  - IrecvのFortranインターフェイス関数
 *  @param[inout] buf       受信バッファ
 *  @param[in]    count     受信バッファのサイズ(ワード数)
 *  @param[in]    datatype  データタイプ(fparam.fiを参照)
 *  @param[in]    source    送信元先のランク番号(procGrpNo内でのランク番号)
 *  @param[in]    procGrpNo プロセスグループ番号
 *  @param[in]    reqNo     リクエスト番号
 *  @param[out]   ierr      終了コード(0=正常終了、0以外=cpm_ErrorCodeの値)
 */
CPM_EXTERN
void
cpm_Irecv_( void *buf, int *count, int *datatype, int *source, int *procGrpNo, int *reqNo, int *ierr )
{
  if( !buf || !count || !datatype || !source || !procGrpNo || !reqNo || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  // インスタンス取得
  cpm_ParaManager *paraMngr = cpm_ParaManager::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  // cpm_Irecv
  if( (*ierr = paraMngr->cpm_Irecv( buf, *count, *datatype, *source, reqNo, *procGrpNo ) ) != CPM_SUCCESS )
  {
    return;
  }

  *ierr = CPM_SUCCESS;
  return;
}

////////////////////////////////////////////////////////////////////////////////
// cpm_Wait
cpm_ErrorCode
cpm_ParaManager::cpm_Wait( int reqNo )
{
  // MPI_Request
  MPI_Request *req = m_reqList.Get(reqNo);
  if( !req )
  {
    return CPM_ERROR_INVALID_OBJKEY;
  }

  // MPI_Wait
  MPI_Status status;
  if( MPI_Wait( req, &status ) != MPI_SUCCESS )
  {
    return CPM_ERROR_MPI_WAIT;
  }

  // 削除
  return m_reqList.Delete(reqNo);
}

////////////////////////////////////////////////////////////////////////////////
// cpm_Waitall
cpm_ErrorCode
cpm_ParaManager::cpm_Waitall( int count, int reqNoList[] )
{
  // MPI_Request
  MPI_Status* stat = new MPI_Status [count];
  MPI_Request* req = new MPI_Request[count];
  int cnt = 0;
  for( int i=0;i<count;i++ )
  {
    MPI_Request *r = m_reqList.Get( reqNoList[i] );
    if( r )
    {
      r[cnt++] = *req;
    }
  }
  if( cnt == 0 )
  {
    delete [] stat;
    delete [] req;
    return CPM_ERROR_INVALID_OBJKEY;
  }

  // MPI_Waitall
  if( MPI_Waitall( cnt, req, stat ) != MPI_SUCCESS )
  {
    delete [] stat;
    delete [] req;
    return CPM_ERROR_MPI_WAITALL;
  }

  // 削除
  for( int i=0;i<count;i++ )
  {
    m_reqList.Delete( reqNoList[i] );
  }

  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// cpm_Isend
cpm_ErrorCode
cpm_ParaManager::cpm_Isend( void *buf, int count, int datatype, int dest, int *reqNo, int procGrpNo )
{
  // MPI_Datatype
  MPI_Datatype dtype = cpm_ParaManager::GetMPI_Datatype( datatype );
  if( dtype == MPI_DATATYPE_NULL )
  {
    return CPM_ERROR_MPI_INVALID_DATATYPE;
  }

  // MPI_Request
  MPI_Request *req = m_reqList.Create();
  if( !req )
  {
    return CPM_ERROR_INVALID_PTR;
  }

  // Isend
  cpm_ErrorCode ret = Isend( dtype, buf, count, dest, req, procGrpNo );
  if( ret != MPI_SUCCESS )
  {
    delete req;
    return ret;
  }

  // MPI_Requestを登録
  if( (*reqNo = m_reqList.Add(req) ) < 0 )
  {
    delete req;
    return CPM_ERROR_REGIST_OBJKEY;
  }

  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// cpm_Irecv
cpm_ErrorCode
cpm_ParaManager::cpm_Irecv( void *buf, int count, int datatype, int source, int *reqNo, int procGrpNo )
{
  // MPI_Datatype
  MPI_Datatype dtype = cpm_ParaManager::GetMPI_Datatype( datatype );
  if( dtype == MPI_DATATYPE_NULL )
  {
    return CPM_ERROR_MPI_INVALID_DATATYPE;
  }

  // Irecv
  MPI_Request req;
  cpm_ErrorCode ret = Irecv( dtype, buf, count, source, &req, procGrpNo );
  if( ret != MPI_SUCCESS )
  {
    return ret;
  }

  // MPI_Requestを登録
  MPI_Request *r = m_reqList.Create();
  *r = req;
  if( (*reqNo = m_reqList.Add(r) ) < 0 )
  {
    delete r;
    return CPM_ERROR_REGIST_OBJKEY;
  }

  return CPM_SUCCESS;
}

