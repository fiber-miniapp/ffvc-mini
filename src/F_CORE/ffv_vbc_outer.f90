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

!> @file   ffv_bc_outer.f90
!! @brief  外部境界条件
!! @author kero
!<

!> ********************************************************************
!! @brief 外部速度境界条件による対流項と粘性項の流束の修正
!! @param[out] wv 疑似ベクトルの空間項の評価値
!! @param sz 配列長
!! @param g ガイドセル長
!! @param dh 格子幅
!! @param v00 参照速度
!! @param rei Reynolds数の逆数
!! @param v0 速度ベクトル（n-step）
!! @param bv BCindex V
!! @param vec 指定する速度ベクトル
!! @param face 外部境界処理のときの面番号
!! @param[out] flop
!! @note vecには，流入条件のとき指定速度，流出境界の流束はローカルのセルフェイス速度を使うこと
!! @todo 内部と外部の分離 do loopの内側に条件分岐を入れているので修正
!<
    subroutine pvec_vobc_wall  (wv, sz, g, dh, v00, rei, v0, bv, vec, face, flop)
    implicit none
    include 'ffv_f_params.h'
    integer                                                   ::  i, j, k, g, bvx, face
    integer                                                   ::  ix, jx, kx
    integer, dimension(3)                                     ::  sz
    double precision                                          ::  flop
    real                                                      ::  Up0, Ue1, Uw1, Us1, Un1, Ub1, Ut1
    real                                                      ::  Vp0, Ve1, Vw1, Vs1, Vn1, Vb1, Vt1
    real                                                      ::  Wp0, We1, Ww1, Ws1, Wn1, Wb1, Wt1
    real                                                      ::  dh, dh1, dh2, EX, EY, EZ, rei
    real                                                      ::  u_ref, v_ref, w_ref, m
    real                                                      ::  u_bc, v_bc, w_bc, u_bc_ref, v_bc_ref, w_bc_ref
    real                                                      ::  u_bc_ref2, v_bc_ref2, w_bc_ref2
    real, dimension(3, 1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  v0, wv
    real, dimension(0:3)                                      ::  v00
    real, dimension(3)                                        ::  vec
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bv
    
    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    
    dh1= 1.0/dh
    dh2= rei*dh1*dh1

    ! 参照座標系の速度
    u_ref = v00(1)
    v_ref = v00(2)
    w_ref = v00(3)
    
    ! u_bcは境界速度
    u_bc = vec(1)
    v_bc = vec(2)
    w_bc = vec(3)
    
    ! u_bc_refは参照座標系での境界速度
    u_bc_ref = u_bc + u_ref
    v_bc_ref = v_bc + v_ref
    w_bc_ref = w_bc + w_ref
    
    u_bc_ref2 = 2.0*u_bc_ref
    v_bc_ref2 = 2.0*v_bc_ref
    w_bc_ref2 = 2.0*w_bc_ref
    
    flop = flop + 16.0d0 ! DP 21 flops

    m = 0.0

!$OMP PARALLEL REDUCTION(+:m) &
!$OMP FIRSTPRIVATE(ix, jx, kx, u_bc_ref2, v_bc_ref2, w_bc_ref2, dh2, face) &
!$OMP PRIVATE(i, j, k, bvx, Up0, Vp0, Wp0, EX, EY, EZ) &
!$OMP PRIVATE(Ue1, Uw1, Us1, Un1, Ub1, Ut1) &
!$OMP PRIVATE(Ve1, Vw1, Vs1, Vn1, Vb1, Vt1) &
!$OMP PRIVATE(We1, Ww1, Ws1, Wn1, Wb1, Wt1)

    FACES : select case (face)
    
    case (X_minus)
      i = 1
      
#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif
      do k=1,kx
      do j=1,jx
        bvx = bv(i,j,k)
        
        if ( ibits(bvx, bc_face_W, bitw_5) == obc_mask ) then
          Up0 = v0(1,i,j,k)
          Vp0 = v0(2,i,j,k)
          Wp0 = v0(3,i,j,k)
          
          Uw1  = u_bc_ref2 - Up0
          Vw1  = v_bc_ref2 - Vp0
          Ww1  = w_bc_ref2 - Wp0
          
          EX = Uw1 - Up0
          EY = Vw1 - Vp0
          EZ = Ww1 - Wp0

          wv(1,i,j,k) = wv(1,i,j,k) + EX*dh2
          wv(2,i,j,k) = wv(2,i,j,k) + EY*dh2
          wv(3,i,j,k) = wv(3,i,j,k) + EZ*dh2
          m = m + 1.0
        endif
      end do
      end do
!$OMP END DO
      
      
    case (X_plus)
      i = ix
      
#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif
      do k=1,kx
      do j=1,jx
        bvx = bv(i,j,k)
        
        if ( ibits(bvx, bc_face_E, bitw_5) == obc_mask ) then
          Up0 = v0(1,i,j,k)
          Vp0 = v0(2,i,j,k)
          Wp0 = v0(3,i,j,k)

          Ue1  = u_bc_ref2 - Up0
          Ve1  = v_bc_ref2 - Vp0
          We1  = w_bc_ref2 - Wp0

          EX = Ue1 - Up0
          EY = Ve1 - Vp0
          EZ = We1 - Wp0

          wv(1,i,j,k) = wv(1,i,j,k) + EX*dh2
          wv(2,i,j,k) = wv(2,i,j,k) + EY*dh2
          wv(3,i,j,k) = wv(3,i,j,k) + EZ*dh2
          m = m + 1.0
        endif
      end do
      end do
!$OMP END DO
      
      
    case (Y_minus)
      j = 1
      
#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif
      do k=1,kx
      do i=1,ix
        bvx = bv(i,j,k)
        
        if ( ibits(bvx, bc_face_S, bitw_5) == obc_mask ) then
          Up0 = v0(1,i,j,k)
          Vp0 = v0(2,i,j,k)
          Wp0 = v0(3,i,j,k)

          Us1  = u_bc_ref2 - Up0
          Vs1  = v_bc_ref2 - Vp0
          Ws1  = w_bc_ref2 - Wp0
        
          EX = Us1 - Up0
          EY = Vs1 - Vp0
          EZ = Ws1 - Wp0

          wv(1,i,j,k) = wv(1,i,j,k) + EX*dh2
          wv(2,i,j,k) = wv(2,i,j,k) + EY*dh2
          wv(3,i,j,k) = wv(3,i,j,k) + EZ*dh2
          m = m + 1.0
        endif
      end do
      end do
!$OMP END DO
      
      
    case (Y_plus)
      j = jx
      
#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif
      do k=1,kx
      do i=1,ix
        bvx = bv(i,j,k)
        
        if ( ibits(bvx, bc_face_N, bitw_5) == obc_mask ) then
          Up0 = v0(1,i,j,k)
          Vp0 = v0(2,i,j,k)
          Wp0 = v0(3,i,j,k)

          Un1  = u_bc_ref2 - Up0
          Vn1  = v_bc_ref2 - Vp0
          Wn1  = w_bc_ref2 - Wp0

          EX = Un1 - Up0
          EY = Vn1 - Vp0
          EZ = Wn1 - Wp0
          
          wv(1,i,j,k) = wv(1,i,j,k) + EX*dh2
          wv(2,i,j,k) = wv(2,i,j,k) + EY*dh2
          wv(3,i,j,k) = wv(3,i,j,k) + EZ*dh2
          m = m + 1.0
        endif
      end do
      end do
!$OMP END DO
      
      
    case (Z_minus)
      k = 1
      
#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif
      do j=1,jx
      do i=1,ix
        bvx = bv(i,j,k)
        
        if ( ibits(bvx, bc_face_B, bitw_5) == obc_mask ) then
          Up0 = v0(1,i,j,k)
          Vp0 = v0(2,i,j,k)
          Wp0 = v0(3,i,j,k)

          Ub1  = u_bc_ref2 - Up0
          Vb1  = v_bc_ref2 - Vp0
          Wb1  = w_bc_ref2 - Wp0
          
          EX = Ub1 - Up0
          EY = Vb1 - Vp0
          EZ = Wb1 - Wp0
          
          wv(1,i,j,k) = wv(1,i,j,k) + EX*dh2
          wv(2,i,j,k) = wv(2,i,j,k) + EY*dh2
          wv(3,i,j,k) = wv(3,i,j,k) + EZ*dh2
          m = m + 1.0
        endif
      end do
      end do
!$OMP END DO
      
      
    case (Z_plus)
      k = kx
      
#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif
      do j=1,jx
      do i=1,ix
        bvx = bv(i,j,k)
        
        if ( ibits(bvx, bc_face_T, bitw_5) == obc_mask ) then
          Up0 = v0(1,i,j,k)
          Vp0 = v0(2,i,j,k)
          Wp0 = v0(3,i,j,k)

          Ut1  = u_bc_ref2 - Up0
          Vt1  = v_bc_ref2 - Vp0
          Wt1  = w_bc_ref2 - Wp0
          
          EX = Ut1 - Up0
          EY = Vt1 - Vp0
          EZ = Wt1 - Wp0
          
          wv(1,i,j,k) = wv(1,i,j,k) + EX*dh2
          wv(2,i,j,k) = wv(2,i,j,k) + EY*dh2
          wv(3,i,j,k) = wv(3,i,j,k) + EZ*dh2
          m = m + 1.0
        endif
      end do
      end do
!$OMP END DO
  
    case default
    end select FACES
    
!$OMP END PARALLEL

    flop = flop + dble(m)*12.0d0
    
    return
    end subroutine pvec_vobc_wall


!> ********************************************************************
!! @brief 外部指定境界条件による速度の発散の修正
!! @param [in,out] div   速度の発散
!! @param [in]     sz    配列長
!! @param [in]     g     ガイドセル長
!! @param [in]     face  面番号
!! @param [in]     v00   参照速度
!! @param [in]     bv    BCindex V
!! @param [in]     vec   指定する速度ベクトル
!! @param [out]    flop flop count
!! @note 指定面でも固体部分は対象外とするのでループ中に判定あり
!<
    subroutine div_obc_drchlt (div, sz, g, face, v00, bv, vec, flop)
    implicit none
    include 'ffv_f_params.h'
    integer                                                   ::  i, j, k, g, ix, jx, kx, face, bvx
    integer, dimension(3)                                     ::  sz
    double precision                                          ::  flop, rix, rjx, rkx
    real                                                      ::  u_bc_ref, v_bc_ref, w_bc_ref
    real, dimension(0:3)                                      ::  v00
    real, dimension(3)                                        ::  vec
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  div
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bv

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)

    rix = dble(jx)*dble(kx)
    rjx = dble(ix)*dble(kx)
    rkx = dble(ix)*dble(jx)
    
    ! 参照座標系の速度に係数をかけておく
    u_bc_ref = (vec(1) + v00(1))
    v_bc_ref = (vec(2) + v00(2))
    w_bc_ref = (vec(3) + v00(3))
    
    flop = flop + 3.0d0

!$OMP PARALLEL REDUCTION(+:flop) &
!$OMP FIRSTPRIVATE(ix, jx, kx, u_bc_ref, v_bc_ref, w_bc_ref, face) &
!$OMP PRIVATE(i, j, k, bvx)

    FACES : select case (face)
    case (X_minus)
      i = 1
      
#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif
      do k=1,kx
      do j=1,jx
        bvx = bv(i,j,k)
        if ( ibits(bvx, bc_face_W, bitw_5) == obc_mask ) then
          div(i,j,k) = div(i,j,k) - u_bc_ref * real(ibits(bvx, State, 1))
        endif
      end do
      end do
!$OMP END DO

!     flop = flop + rix*3.0d0 ! 2+ real*1
      flop = flop + rix*2.0d0 ! exclude cast op
      
      
    case (X_plus)
      i = ix
      
#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif
      do k=1,kx
      do j=1,jx
        bvx = bv(i,j,k)
        if ( ibits(bvx, bc_face_E, bitw_5) == obc_mask ) then
          div(i,j,k) = div(i,j,k) + u_bc_ref * real(ibits(bvx, State, 1))
        endif
      end do
      end do
!$OMP END DO
      
!     flop = flop + rix*3.0d0 ! 2+ real*1
      flop = flop + rix*2.0d0 ! exclude cast op
      
      
    case (Y_minus)
      j = 1
      
#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif
      do k=1,kx
      do i=1,ix
        bvx = bv(i,j,k)
        if ( ibits(bvx, bc_face_S, bitw_5) == obc_mask ) then
          div(i,j,k) = div(i,j,k) - v_bc_ref * real(ibits(bvx, State, 1))
        endif
      end do
      end do
!$OMP END DO
      
!     flop = flop + rjx*3.0d0 ! 2+ real*1
      flop = flop + rjx*2.0d0 ! exclude cast op
      
      
    case (Y_plus)
      j = jx
      
#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif
      do k=1,kx
      do i=1,ix
        bvx = bv(i,j,k)
        if ( ibits(bvx, bc_face_N, bitw_5) == obc_mask ) then
          div(i,j,k) = div(i,j,k) + v_bc_ref * real(ibits(bvx, State, 1))
        endif
      end do
      end do
!$OMP END DO
      
!     flop = flop + rjx*3.0d0 ! 2+ real*1
      flop = flop + rjx*2.0d0 ! exclude cast op
    
    
    case (Z_minus)
      k = 1
      
#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif
      do j=1,jx
      do i=1,ix
        bvx = bv(i,j,k)
        if ( ibits(bvx, bc_face_B, bitw_5) == obc_mask ) then
          div(i,j,k) = div(i,j,k) - w_bc_ref * real(ibits(bvx, State, 1))
        endif
      end do
      end do
!$OMP END DO
      
!     flop = flop + rkx*3.0d0 ! 2+ real*1
      flop = flop + rkx*2.0d0 ! exclude cast op
      
      
    case (Z_plus)
      k = kx
      
#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif
      do j=1,jx
      do i=1,ix
        bvx = bv(i,j,k)
        if ( ibits(bvx, bc_face_T, bitw_5) == obc_mask ) then
          div(i,j,k) = div(i,j,k) + w_bc_ref * real(ibits(bvx, State, 1))
        endif
      end do
      end do
!$OMP END DO
      
!     flop = flop + rkx*3.0d0 ! 2+ real*1
      flop = flop + rkx*2.0d0 ! exclude cast op
    
    case default
    end select FACES

!$OMP END PARALLEL

    return
    end subroutine div_obc_drchlt

