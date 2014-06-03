!********************************************************************
!
!   FFV : Frontflow / violet
!
! Copyright (c) 2012-2014 Institute of Industrial Science, The University of Tokyo.
! All right reserved.
!
! Copyright (c) 2012-2014 Advanced Institute for Computational Science, RIKEN.
! All right reserved.
!
!********************************************************************

!> @file   FB_util.f90
!! @brief  FlowBase utilities
!! @author kero
!<

!> ********************************************************************
!! @brief スカラー値の書き出し
!! @param [in] s        スカラー値
!! @param [in] sz       配列長
!! @param [in] g        ガイドセル長
!! @param [in] fname    ファイル名
!! @param [in] step     ステップ数
!! @param [in] time     時刻
!! @param [in] org      起点座標
!! @param [in] pit      格子幅
!! @param [in] d_type   (1-float, 2-double)
!! @param [in] gs       guide cell (0-without, others-with)
!! @param [in] avs      平均値識別子 (0-average, 1-instantaneous) 
!! @param [in] step_avr 平均ステップ
!! @param [in] time_avr 平均時刻
!<
  subroutine fb_write_sph_s(s, sz, g, fname, step, time, org, pit, d_type, gs, avs, step_avr, time_avr)
  implicit none
  integer                                                   ::  i, j, k, ix, jx, kx, g, step, gs, step_avr, avs
  integer                                                   ::  sv_type, d_type, imax, jmax, kmax
  integer, dimension(3)                                     ::  sz
  real                                                      ::  time, time_avr
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  s
  real, dimension(3)                                        ::  org, pit
  character(64)                                             ::  fname

  ix = sz(1)
  jx = sz(2)
  kx = sz(3)

  sv_type = 1 ! scalar

  if ( gs == 0 ) then
    imax = ix
    jmax = jx
    kmax = kx
  else
    imax = ix + 2*g
    jmax = jx + 2*g
    kmax = kx + 2*g
  end if

  open(16, file=fname, form='unformatted')
  write(16) sv_type, d_type
  write(16) imax, jmax, kmax
  write(16) org(1), org(2), org(3)
  write(16) pit(1), pit(2), pit(3)
  write(16) step, time

  if ( gs /= 0 ) then
    write(16) (((s(i,j,k),i=1-g,ix+g),j=1-g,jx+g),k=1-g,kx+g)
  else
    write(16) (((s(i,j,k),i=1,ix),j=1,jx),k=1,kx)
  end if
  
  if ( avs == 0 ) then
    write(16) step_avr, time_avr
  end if
  
  close (unit=16)

  return
  end subroutine fb_write_sph_s


!> ********************************************************************
!! @brief ベクトル値の書き出し
!! @param [in] v        ベクトル値
!! @param [in] sz       配列長
!! @param [in] g        ガイドセル長
!! @param [in] fname    ファイル名
!! @param [in] step     ステップ数
!! @param [in] time     時刻
!! @param [in] org      起点座標
!! @param [in] pit      格子幅
!! @param [in] d_type   (1-float, 2-double)
!! @param [in] gs       guide cell (0-without, others-with)
!! @param [in] avs      平均値識別子 (0-average, 1-instantaneous)
!! @param [in] step_avr 平均ステップ
!! @param [in] time_avr 平均時刻
!<
  subroutine fb_write_sph_v(v, sz, g, fname, step, time, org, pit, d_type, gs, avs, step_avr, time_avr)
  implicit none
  integer                                                   ::  i, j, k, l, ix, jx, kx, g, step, gs, step_avr, avs
  integer                                                   ::  sv_type, d_type, imax, jmax, kmax
  integer, dimension(3)                                     ::  sz
  real                                                      ::  time, time_avr
  real, dimension(3, 1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  v
  real, dimension(3)                                        ::  org, pit
  character(64)                                             ::  fname

  ix = sz(1)
  jx = sz(2)
  kx = sz(3)

  sv_type = 2 ! vector

  if ( gs == 0 ) then
    imax = ix
    jmax = jx
    kmax = kx
  else
    imax = ix + 2*g
    jmax = jx + 2*g
    kmax = kx + 2*g
  end if

  open(16, file=fname, form='unformatted')
  write(16) sv_type, d_type
  write(16) imax, jmax, kmax
  write(16) org(1), org(2), org(3)
  write(16) pit(1), pit(2), pit(3)
  write(16) step, time

  if ( gs /= 0 ) then
    write(16) ((((v(l,i,j,k),l=1,3),i=1-g,ix+g),j=1-g,jx+g),k=1-g,kx+g)
  else
    write(16) ((((v(l,i,j,k),l=1,3),i=1,ix),j=1,jx),k=1,kx)
  end if
  
  if ( avs == 0 ) then
    write(16) step_avr, time_avr
  end if
  
  close (unit=16)

  return
  end subroutine fb_write_sph_v

  
!> ********************************************************************
!! @brief 有効セルに対する，1タイムステップ進行時のベクトルの絶対値の変化量の和と平均値
!! @param [out] d    戻り値（変化量の2乗和と平均値）
!! @param [in]  sz   配列長
!! @param [in]  g    ガイドセル長
!! @param [in]  vn   ベクトル値 n+1 step
!! @param [in]  vo   ベクトル値 n step
!! @param [in]  bx   BCindex
!! @param [out] flop 浮動小数演算数
!<
  subroutine fb_delta_v (d, sz, g, vn, vo, bx, flop)
  implicit none
  include 'ffv_f_params.h'
  integer                                                   ::  i, j, k, ix, jx, kx, g
  integer, dimension(3)                                     ::  sz
  double precision                                          ::  u, v, w, x, y, z, actv
  double precision                                          ::  flop, av, rm
  real, dimension(3, 1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  vn, vo
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bx
  double precision, dimension(2)                            ::  d

  ix = sz(1)
  jx = sz(2)
  kx = sz(3)

  flop = flop + dble(ix)*dble(jx)*dble(kx)*58.0d0
  ! sqrt float->10, double->20

  av = 0.0
  rm = 0.0

!$OMP PARALLEL &
!$OMP REDUCTION(+:av) &
!$OMP REDUCTION(+:rm) &
!$OMP PRIVATE(actv, u, v, w, x, y, z) &
!$OMP FIRSTPRIVATE(ix, jx, kx)

!$OMP DO SCHEDULE(static)

  do k=1,kx
  do j=1,jx
  do i=1,ix

    actv = real(ibits(bx(i,j,k), State, 1))
    
    u = dble(vn(1,i,j,k))
    v = dble(vn(2,i,j,k))
    w = dble(vn(3,i,j,k))
    av = av + sqrt(u*u + v*v + w*w)*actv
    
    x = u - dble(vo(1,i,j,k))
    y = v - dble(vo(2,i,j,k))
    z = w - dble(vo(3,i,j,k))
    rm = rm + sqrt(x*x + y*y + z*z)*actv

  end do
  end do
  end do
!$OMP END DO
!$OMP END PARALLEL

  d(1) = rm
  d(2) = av

  return
  end subroutine fb_delta_v

!> ********************************************************************
!! @brief 有効セルに対する，1タイムステップ進行時の変化量の2乗和と平均値(RootMean)
!! @param [out] d    戻り値（変化量の2乗和と平均値）
!! @param [in]  sz   配列長
!! @param [in]  g    ガイドセル長
!! @param [in]  sn   スカラー値 n+1 step
!! @param [in]  so   スカラー値 n step
!! @param [in]  bx   BCindex
!! @param [out] flop 浮動小数演算数
!<
  subroutine fb_delta_s (d, sz, g, sn, so, bx, flop)
  implicit none
  include 'ffv_f_params.h'
  integer                                                   ::  i, j, k, ix, jx, kx, g
  integer, dimension(3)                                     ::  sz
  double precision                                          ::  a, s, av, rm, actv, flop
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  sn, so
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bx
  double precision, dimension(2)                            ::  d

  ix = sz(1)
  jx = sz(2)
  kx = sz(3)

  flop = flop + dble(ix)*dble(jx)*dble(kx)*7.0d0

  av = 0.0
  rm = 0.0

!$OMP PARALLEL &
!$OMP REDUCTION(+:av) &
!$OMP REDUCTION(+:rm) &
!$OMP PRIVATE(actv, s, a) &
!$OMP FIRSTPRIVATE(ix, jx, kx)

!$OMP DO SCHEDULE(static)

  do k=1,kx
  do j=1,jx
  do i=1,ix
    actv = dble(ibits(bx(i,j,k), State,  1))
    
    s = dble(sn(i,j,k))
    av = av + s * actv
    
    a = ( s - dble(so(i,j,k)) )*actv
    rm = rm + a*a
  end do
  end do
  end do
!$OMP END DO
!$OMP END PARALLEL
  
  d(1) = rm
  d(2) = av

  return
  end subroutine fb_delta_s

!> ********************************************************************
!! @brief ベクトル値を設定する
!! @param [out] var ベクトル配列
!! @param [in]  sz  配列長
!! @param [in]  g   ガイドセル長
!! @param [in]  val ベクトル値
!! @param [in]  bv  BCindex V
!<
    subroutine fb_set_vector (var, sz, g, val, bv)
    implicit none
    include '../FB/ffv_f_params.h'
    integer                                                   ::  i, j, k, ix, jx, kx, g, bvx
    integer, dimension(3)                                     ::  sz
    real                                                      ::  u1, u2, u3
    real, dimension(3)                                        ::  val
    real, dimension(3, 1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  var
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bv

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    u1 = val(1)
    u2 = val(2)
    u3 = val(3)

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(u1, u2, u3) &
!$OMP FIRSTPRIVATE(ix, jx, kx, g) &
!$OMP PRIVATE(bvx)

!$OMP DO SCHEDULE(static)

    do k=1-g, kx+g
    do j=1-g, jx+g
    do i=1-g, ix+g
        bvx = ibits(bv(i,j,k), State, 1)
        var(1,i,j,k) = u1 * real(bvx)
        var(2,i,j,k) = u2 * real(bvx)
        var(3,i,j,k) = u3 * real(bvx)
    end do
    end do
    end do
!$OMP END DO
!$OMP END PARALLEL

    return
    end subroutine fb_set_vector

