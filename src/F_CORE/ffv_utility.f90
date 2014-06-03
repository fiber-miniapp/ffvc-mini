!********************************************************************
!
!   FFV : Frontflow / violet Cartesian
!
! Copyright (c) 2012-2014 Institute of Industrial Science, The University of Tokyo.
! All right reserved.
!
! Copyright (c) 2012-2014 Advanced Institute for Computational Science, RIKEN.
! All right reserved.
!
!********************************************************************

!> @file   ffv_utility.f90
!! @brief  Utilitiy functions
!! @author kero
!<


!> ********************************************************************
!! @brief 速度成分の最大値を計算する
!! @param [out] ds   最大値
!! @param [in]  sz   配列長
!! @param [in]  g    ガイドセル長
!! @param [in]  div  速度の発散
!! @param [in]  coef 係数
!! @param [in]  bp   BCindex P
!! @param [out] flop flop count
!<
    subroutine norm_v_div_max (ds, sz, g, div, coef, bp, flop)
    implicit none
    include 'ffv_f_params.h'
    integer                                                   ::  i, j, k, ix, jx, kx, g
    integer, dimension(3)                                     ::  sz
    double precision                                          ::  flop, ds, r
    real                                                      ::  coef
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  div
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bp

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    ds = 0.0

!   flop = flop + dble(ix)*dble(jx)*dble(kx)*6.0d0
    flop = flop + dble(ix)*dble(jx)*dble(kx)*4.0d0 ! exclude cast ops

!$OMP PARALLEL &
!$OMP REDUCTION(max:ds) &
!$OMP PRIVATE(r) &
!$OMP FIRSTPRIVATE(ix, jx, kx, coef)

#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif

    do k=1,kx
    do j=1,jx
    do i=1,ix
      r = dble(div(i,j,k) * coef) * dble(ibits(bp(i,j,k), vld_cnvg, 1)) ! 有効セルの場合 1.0
      ds = max(ds, abs(r) )
    end do
    end do
    end do
!$OMP END DO
!$OMP END PARALLEL

    return
    end subroutine norm_v_div_max


!> ********************************************************************
!! @brief 速度成分の最大値を計算する
!! @param [out] v_max 最大値
!! @param [in]  sz    配列長
!! @param [in]  g     ガイドセル長
!! @param [in]  v00   参照速度
!! @param [in]  v     速度ベクトル
!! @param [out] flop  flop count
!<
    subroutine find_vmax (v_max, sz, g, v00, v, flop)
    implicit none
    integer                                                   ::  i, j, k, ix, jx, kx, g
    integer, dimension(3)                                     ::  sz
    double precision                                          ::  flop
    real                                                      ::  vm1, vm2, vm3, v_max, vx, vy, vz
    real, dimension(3, 1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  v
    real, dimension(0:3)                                      ::  v00

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    vm1 = 0.0
    vm2 = 0.0
    vm3 = 0.0
    vx = v00(1)
    vy = v00(2)
    vz = v00(3)
    flop = flop + dble(ix)*dble(jx)*dble(kx)*9.0d0 + 2.0d0

!$OMP PARALLEL &
!$OMP REDUCTION(max:vm1) &
!$OMP REDUCTION(max:vm2) &
!$OMP REDUCTION(max:vm3) &
!$OMP FIRSTPRIVATE(ix, jx, kx, vx, vy, vz)

#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif

    do k=1,kx
    do j=1,jx
    do i=1,ix
      vm1 = max(vm1, abs(v(1,i,j,k)-vx ) )
      vm2 = max(vm2, abs(v(2,i,j,k)-vy ) )
      vm3 = max(vm3, abs(v(3,i,j,k)-vz ) )
    end do
    end do
    end do
!$OMP END DO
!$OMP END PARALLEL
    
    v_max = max(vm1, vm2, vm3) ! maxss %xmm0, %xmm1, x 2 times > 2 flop

    return
    end subroutine find_vmax
