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

!> @file   ffv_velocity_binary.f90
!! @brief  速度計算のルーチン群（バイナリモデル）
!! @author kero
!<

!> ********************************************************************
!! @brief 対流項と粘性項の計算
!! @param[out] wv 疑似ベクトルの空間項の評価値
!! @param sz 配列長
!! @param g ガイドセル長
!! @param dh 格子幅
!! @param c_scheme 対流項スキームのモード（1-UWD, 2-center, 3-MUSCL）
!! @param v00 参照速度
!! @param rei レイノルズ数の逆数
!! @param v 速度ベクトル（n-step, collocated）
!! @param bv BCindex V
!! @param bp BCindex P
!! @param v_mode 粘性項のモード (0=対流項のみ, 1=対流項と粘性項，2=粘性項は壁法則)
!! @param ut 摩擦速度
!! @param wall_type 壁面条件 (0=no_slip, 1=slip)
!! @param bd BCindex ID
!! @param [in]  vcs_coef 粘性項の係数
!! @param [out] flop
!<
    subroutine pvec_muscl (wv, sz, g, dh, c_scheme, v00, rei, v, bv, bp, v_mode, ut, wall_type, bd, vcs_coef, flop)
    implicit none
    include 'ffv_f_params.h'
    integer ::  i, j, k, ix, jx, kx, g, c_scheme, bvx, v_mode, bpx, wall_type, bdx
    integer ::  b_e1, b_w1, b_n1, b_s1, b_t1, b_b1, b_e2, b_w2, b_n2, b_s2, b_t2, b_b2, b_p
    integer, dimension(3) ::  sz
    double precision ::  flop, vflop
    real ::  UPe, UPw, VPn, VPs, WPt, WPb, u1, u2, u3, ug, e1, e2, e3, u_tau
    real ::  Up0, Ue1, Ue2, Uw1, Uw2, Us1, Us2, Un1, Un2, Ub1, Ub2, Ut1, Ut2
    real ::  Vp0, Ve1, Ve2, Vw1, Vw2, Vs1, Vs2, Vn1, Vn2, Vb1, Vb2, Vt1, Vt2
    real ::  Wp0, We1, We2, Ww1, Ww2, Ws1, Ws2, Wn1, Wn2, Wb1, Wb2, Wt1, Wt2
    real ::  ck, dh, dh1, dh2, vcs, vcs_coef
    real ::  u_ref, v_ref, w_ref, u_ref2, v_ref2, w_ref2
    real ::  c_e, c_w, c_n, c_s, c_t, c_b, wls, wm1, wm2, cm1, cm2, ss_4, tmp1, tmp2
    real ::  w_e, w_w, w_n, w_s, w_t, w_b
    real ::  uu_e, uu_w, uu_s, uu_n, uu_b, uu_t
    real ::  vv_e, vv_w, vv_s, vv_n, vv_b, vv_t
    real ::  ww_e, ww_w, ww_s, ww_n, ww_b, ww_t
    real ::  dv1, dv2, dv3, dv4, g1, g2, g3, g4, g5, g6, s1, s2, s3, s4, b
    real ::  Urr, Url, Ulr, Ull, Vrr, Vrl, Vlr, Vll, Wrr, Wrl, Wlr, Wll
    real ::  cr, cl, acr, acl, cnv_u, cnv_v, cnv_w, EX, EY, EZ, rei, beta, qtz
    real ::  fu_r, fu_l, fv_r, fv_l, fw_r, fw_l, uq, vq, wq, ss
    real ::  lmt_w, lmt_e, lmt_s, lmt_n, lmt_b, lmt_t
    real, dimension(3, 1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  v, wv
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  ut
    real, dimension(0:3)                                      ::  v00
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bv, bp, bd
    
    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    
    dh1= 1.0/dh
    dh2= rei*dh1
    qtz = 1.0/255.0
    
    ! 粘性項のマスク
    if ( v_mode == 0 ) then ! 粘性項は計算しない
      vcs = 0.0
    else
      vcs = 1.0
    endif

    ! vcs = 1.0 (Euler Explicit) / 0.5 (CN)
    vcs = vcs * vcs_coef
    
    ! 参照座標系の速度
    u_ref = v00(1)
    v_ref = v00(2)
    w_ref = v00(3)
    u_ref2 = 2.0*u_ref
    v_ref2 = 2.0*v_ref
    w_ref2 = 2.0*w_ref
    
    ck = 0.0
    b  = 0.0
    ss = 1.0
            
    if ( c_scheme == 1 ) then      !     1st order upwind
      ss = 0.0
    else if ( c_scheme == 2 ) then !     2nd order central 
      ck = 1.0     
    else if ( c_scheme == 3 ) then !     3rd order MUSCL
      ck = 1.0/3.0
      b  = (3.0-ck)/(1.0-ck)
    else
      write(*,*) 'out of scheme selection'
      stop
    endif
    
    ss_4 = 0.25*ss
    
    cm1 = 1.0 - ck
    cm2 = 1.0 + ck
    
    ! 壁面条件
    if ( wall_type == 0 ) then ! no slip
      wls = 0.0
    else if ( wall_type == 1 ) then ! slip
      wls = 1.0
    endif
    
    wm1 = 2.0*(1.0-wls)
    wm2 = 2.0*wls-1.0
    
! /4 + 12 = 8*4+13  = 45 flops
! /4 + 12 = 13*4+13 = 65 flops ! DP
! ループ内セットアップ   6 + 36 + 2 = 44 flops
! 対流項の各方向　( 7+6+7 + 78*3 + 12 ) * 3dir = 798
! 粘性項 36 + 36 + 12 = 84
! 粘性項置換部　67 ! DP 92
! Total : 45 + 44 + 798 + 84 = 971 ! DP 65 + 44 + 798 + 84 = 991

!   flop = flop + dble(ix)*dble(jx)*dble(kx)*798.0d0
    flop = flop + dble(ix)*dble(jx)*dble(kx)*792.0d0 ! exclude cast ops
    ! flop = flop + dble(ix)*dble(jx)*dble(kx)*991.0d0 ! DP
    
    vflop = 0.0d0

!$OMP PARALLEL &
!$OMP REDUCTION(+:vflop) &
!$OMP FIRSTPRIVATE(ix, jx, kx, dh1, dh2, qtz, vcs, b, ck, ss_4, ss, cm1, cm2, wls) &
!$OMP FIRSTPRIVATE(u_ref, v_ref, w_ref, u_ref2, v_ref2, w_ref2, wm1, wm2, rei, v_mode) &
!$OMP PRIVATE(cnv_u, cnv_v, cnv_w, bvx, bpx, uq, vq, wq, tmp1, tmp2, bdx) &
!$OMP PRIVATE(Up0, Ue1, Ue2, Uw1, Uw2, Us1, Us2, Un1, Un2, Ub1, Ub2, Ut1, Ut2) &
!$OMP PRIVATE(Vp0, Ve1, Ve2, Vw1, Vw2, Vs1, Vs2, Vn1, Vn2, Vb1, Vb2, Vt1, Vt2) &
!$OMP PRIVATE(Wp0, We1, We2, Ww1, Ww2, Ws1, Ws2, Wn1, Wn2, Wb1, Wb2, Wt1, Wt2) &
!$OMP PRIVATE(b_e1, b_w1, b_n1, b_s1, b_t1, b_b1, b_e2, b_w2, b_n2, b_s2, b_t2, b_b2, b_p) &
!$OMP PRIVATE(w_e, w_w, w_n, w_s, w_t, w_b) &
!$OMP PRIVATE(c_e, c_w, c_n, c_s, c_t, c_b) &
!$OMP PRIVATE(UPe, UPw, VPn, VPs, WPt, WPb) &
!$OMP PRIVATE(cr, cl, acr, acl, beta) &
!$OMP PRIVATE(dv1, dv2, dv3, dv4, g1, g2, g3, g4, g5, g6, s1, s2, s3, s4) &
!$OMP PRIVATE(Urr, Url, Ulr, Ull, Vrr, Vrl, Vlr, Vll, Wrr, Wrl, Wlr, Wll) &
!$OMP PRIVATE(fu_r, fu_l, fv_r, fv_l, fw_r, fw_l) &
!$OMP PRIVATE(uu_e, uu_w, uu_s, uu_n, uu_b, uu_t) &
!$OMP PRIVATE(vv_e, vv_w, vv_s, vv_n, vv_b, vv_t) &
!$OMP PRIVATE(ww_e, ww_w, ww_s, ww_n, ww_b, ww_t) &
!$OMP PRIVATE(u1, u2, u3, ug, e1, e2, e3, u_tau) &
!$OMP PRIVATE(lmt_w, lmt_e, lmt_s, lmt_n, lmt_b, lmt_t)

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
      cnv_u = 0.0
      cnv_v = 0.0
      cnv_w = 0.0

      ! 変数のロード
      include 'load_var_stencil5.h'
      
      bvx = bv(i,j,k)
      bpx = bp(i,j,k)
      
      ! セル状態 (0-solid / 1-fluid)
      b_p = ibits(bvx,             State, 1)
      b_w2= ibits(bv(i-2,j  ,k  ), State, 1)
      b_w1= ibits(bv(i-1,j  ,k  ), State, 1)
      b_e1= ibits(bv(i+1,j  ,k  ), State, 1)
      b_e2= ibits(bv(i+2,j  ,k  ), State, 1)
      b_s2= ibits(bv(i  ,j-2,k  ), State, 1)
      b_s1= ibits(bv(i  ,j-1,k  ), State, 1)
      b_n1= ibits(bv(i  ,j+1,k  ), State, 1)
      b_n2= ibits(bv(i  ,j+2,k  ), State, 1)
      b_b2= ibits(bv(i  ,j  ,k-2), State, 1)
      b_b1= ibits(bv(i  ,j  ,k-1), State, 1)
      b_t1= ibits(bv(i  ,j  ,k+1), State, 1)
      b_t2= ibits(bv(i  ,j  ,k+2), State, 1)
      
      ! セル界面のフラグ (0-wall face / 1-fluid) > real*6= 6 flops
      w_e = real(b_e1 * b_p)
      w_w = real(b_w1 * b_p)
      w_n = real(b_n1 * b_p)
      w_s = real(b_s1 * b_p)
      w_t = real(b_t1 * b_p)
      w_b = real(b_b1 * b_p)
      
      ! 各面のVBCフラグ ibits() = 0(Normal) / others(BC) >> c_e = 1.0(Normal) / 0.0(BC) 
      c_e = 1.0
      c_w = 1.0
      c_n = 1.0
      c_s = 1.0
      c_t = 1.0
      c_b = 1.0
      if ( ibits(bvx, bc_face_E, bitw_5) /= 0 ) c_e = 0.0
      if ( ibits(bvx, bc_face_W, bitw_5) /= 0 ) c_w = 0.0
      if ( ibits(bvx, bc_face_N, bitw_5) /= 0 ) c_n = 0.0
      if ( ibits(bvx, bc_face_S, bitw_5) /= 0 ) c_s = 0.0
      if ( ibits(bvx, bc_face_T, bitw_5) /= 0 ) c_t = 0.0
      if ( ibits(bvx, bc_face_B, bitw_5) /= 0 ) c_b = 0.0
      
      ! 延びたステンシルの参照先がvspec, outflowである場合のスキームの破綻を回避，１次精度におとす
      lmt_w = 1.0
      lmt_e = 1.0
      lmt_s = 1.0
      lmt_n = 1.0
      lmt_b = 1.0
      lmt_t = 1.0

      if ( (ibits(bv(i-1, j  , k  ), bc_face_W, bitw_5) /= 0) .and. (ibits(bp(i-1, j  , k  ), vbc_uwd, 1) == 1) ) lmt_w = 0.0
      if ( (ibits(bv(i+1, j  , k  ), bc_face_E, bitw_5) /= 0) .and. (ibits(bp(i+1, j  , k  ), vbc_uwd, 1) == 1) ) lmt_e = 0.0
      if ( (ibits(bv(i  , j-1, k  ), bc_face_S, bitw_5) /= 0) .and. (ibits(bp(i  , j-1, k  ), vbc_uwd, 1) == 1) ) lmt_s = 0.0
      if ( (ibits(bv(i  , j+1, k  ), bc_face_N, bitw_5) /= 0) .and. (ibits(bp(i  , j+1, k  ), vbc_uwd, 1) == 1) ) lmt_n = 0.0
      if ( (ibits(bv(i  , j  , k-1), bc_face_B, bitw_5) /= 0) .and. (ibits(bp(i  , j  , k-1), vbc_uwd, 1) == 1) ) lmt_b = 0.0
      if ( (ibits(bv(i  , j  , k+1), bc_face_T, bitw_5) /= 0) .and. (ibits(bp(i  , j  , k+1), vbc_uwd, 1) == 1) ) lmt_t = 0.0
      
      ! 外部境界条件の場合
      if ( (i == 1)  .and. (ibits(bp(0   , j   , k   ), vbc_uwd, 1) == 1) ) lmt_w = 0.0
      if ( (i == ix) .and. (ibits(bp(ix+1, j   , k   ), vbc_uwd, 1) == 1) ) lmt_e = 0.0
      if ( (j == 1)  .and. (ibits(bp(i   , 0   , k   ), vbc_uwd, 1) == 1) ) lmt_s = 0.0
      if ( (j == jx) .and. (ibits(bp(i   , jx+1, k   ), vbc_uwd, 1) == 1) ) lmt_n = 0.0
      if ( (k == 1)  .and. (ibits(bp(i   , j   , 0   ), vbc_uwd, 1) == 1) ) lmt_b = 0.0
      if ( (k == kx) .and. (ibits(bp(i   , j   , kx+1), vbc_uwd, 1) == 1) ) lmt_t = 0.0
      
      ! 界面速度（スタガード位置） > 36 flops
      UPe = 0.5*(Up0+Ue1)*w_e + u_ref*(1.0-w_e)
      UPw = 0.5*(Up0+Uw1)*w_w + u_ref*(1.0-w_w)
      VPn = 0.5*(Vp0+Vn1)*w_n + v_ref*(1.0-w_n)
      VPs = 0.5*(Vp0+Vs1)*w_s + v_ref*(1.0-w_s)
      WPt = 0.5*(Wp0+Wt1)*w_t + w_ref*(1.0-w_t)
      WPb = 0.5*(Wp0+Wb1)*w_b + w_ref*(1.0-w_b)
      
      ! セルセンターからの壁面修正速度 > 2 flops
      uq = u_ref2 - Up0
      vq = v_ref2 - Vp0
      wq = w_ref2 - Wp0
			
      ! X方向 ---------------------------------------
      
      ! 速度指定の場合にMUSCLスキームの参照先として，固体内にテンポラリに与えた値を使う
      if ( (b_e2 == 0)  ) then  ! 7 flops
        Ue2 = u_ref2    - v(1, i+1,j  ,k  )
        Ve2 = wm1*v_ref + v(2, i+1,j  ,k  )*wm2
        We2 = wm1*w_ref + v(3, i+1,j  ,k  )*wm2
        !Ve2 = v_ref2 - v(2, i+1,j  ,k  )
        !We2 = w_ref2 - v(3, i+1,j  ,k  )
      endif
      
      ! 壁面の場合の参照速度の修正
      tmp1 = wm1*v_ref + Vp0*wm2 ! 6 flops
      tmp2 = wm1*w_ref + Wp0*wm2
      
      if ( b_e1 == 0 ) then  
        Ue1 = uq
        Ve1 = tmp1
        We1 = tmp2
        !Ve1 = vq
        !We1 = wq
      endif
      
      if ( b_w1 == 0 ) then
        Uw1 = uq
        Vw1 = tmp1
        Ww1 = tmp2
        !Vw1 = vq
        !Ww1 = wq
      end if
      
      if ( (b_w2 == 0)  ) then ! 7 flops
        Uw2 = u_ref2    - v(1, i-1,j  ,k  )
        Vw2 = wm1*v_ref + v(2, i-1,j  ,k  )*wm2
        Ww2 = wm1*w_ref + v(3, i-1,j  ,k  )*wm2
        !Vw2 = v_ref2 - v(2, i-1,j  ,k  )
        !Ww2 = w_ref2 - v(3, i-1,j  ,k  )
      end if
      
      ! 流束　流体のみと固体壁の影響を含み，隣接セルが固体の場合にはマスクする
      cr  = UPe - u_ref
      cl  = UPw - u_ref
      acr = abs(cr)
      acl = abs(cl)
      
      dv4 = Ue2-Ue1
      dv3 = Ue1-Up0
      dv2 = Up0-Uw1
      dv1 = Uw1-Uw2
      
      include 'muscl.h'
      
      Urr = Ue1 - (cm1*g6+cm2*g5)*ss_4 * lmt_e
      Url = Up0 + (cm1*g3+cm2*g4)*ss_4
      Ulr = Up0 - (cm1*g4+cm2*g3)*ss_4
      Ull = Uw1 + (cm1*g1+cm2*g2)*ss_4 * lmt_w
      fu_r = 0.5*(cr*(Urr+Url) - acr*(Urr-Url)) * w_e
      fu_l = 0.5*(cl*(Ulr+Ull) - acl*(Ulr-Ull)) * w_w ! > 4 + 4 + 36 + 5*4+7*2 = 78 flops

      dv4 = Ve2-Ve1
      dv3 = Ve1-Vp0
      dv2 = Vp0-Vw1
      dv1 = Vw1-Vw2
      
      include 'muscl.h'
      
      Vrr = Ve1 - (cm1*g6+cm2*g5)*ss_4 * lmt_e
      Vrl = Vp0 + (cm1*g3+cm2*g4)*ss_4
      Vlr = Vp0 - (cm1*g4+cm2*g3)*ss_4
      Vll = Vw1 + (cm1*g1+cm2*g2)*ss_4 * lmt_w
      fv_r = 0.5*(cr*(Vrr+Vrl) - acr*(Vrr-Vrl)) * w_e
      fv_l = 0.5*(cl*(Vlr+Vll) - acl*(Vlr-Vll)) * w_w

      dv4 = We2-We1
      dv3 = We1-Wp0
      dv2 = Wp0-Ww1
      dv1 = Ww1-Ww2
      
      include 'muscl.h'
      
      Wrr = We1 - (cm1*g6+cm2*g5)*ss_4 * lmt_e
      Wrl = Wp0 + (cm1*g3+cm2*g4)*ss_4
      Wlr = Wp0 - (cm1*g4+cm2*g3)*ss_4
      Wll = Ww1 + (cm1*g1+cm2*g2)*ss_4 * lmt_w
      fw_r = 0.5*(cr*(Wrr+Wrl) - acr*(Wrr-Wrl)) * w_e
      fw_l = 0.5*(cl*(Wlr+Wll) - acl*(Wlr-Wll)) * w_w
      
      ! 流束の加算　VBCでない面の寄与のみを評価する
      cnv_u = cnv_u + fu_r*c_e - fu_l*c_w
      cnv_v = cnv_v + fv_r*c_e - fv_l*c_w
      cnv_w = cnv_w + fw_r*c_e - fw_l*c_w ! > 4*3 = 12 flops
			
      ! Y方向 ---------------------------------------
      
      if ( (b_n2 == 0)  ) then
        Un2 = wm1*u_ref + v(1, i  ,j+1,k  )*wm2
        Vn2 = v_ref2    - v(2, i  ,j+1,k  )
        Wn2 = wm1*w_ref + v(3, i  ,j+1,k  )*wm2
        !Un2 = u_ref2 - v(1, i  ,j+1,k  )
        !Wn2 = w_ref2 - v(3, i  ,j+1,k  )
      endif
      
      tmp1 = wm1*u_ref + Up0*wm2
      tmp2 = wm1*w_ref + Wp0*wm2
      
      if ( b_n1 == 0 ) then
        Un1 = tmp1
        Vn1 = vq
        Wn1 = tmp2
        !Un1 = uq
        !Wn1 = wq
      endif
      
      if ( b_s1 == 0 ) then
        Us1 = tmp1
        Vs1 = vq
        Ws1 = tmp2
        !Us1 = uq
        !Ws1 = wq
      endif
      
      if ( (b_s2 == 0)  ) then
        Us2 = wm1*u_ref + v(1, i  ,j-1,k  )*wm2
        Vs2 = v_ref2    - v(2, i  ,j-1,k  )
        Ws2 = wm1*w_ref + v(3, i  ,j-1,k  )*wm2
        !Us2 = u_ref2 - v(1, i  ,j-1,k  )
        !Ws2 = w_ref2 - v(3, i  ,j-1,k  )
      endif
      
      ! 流束　流体のみと固体壁の影響を含み，隣接セルが固体の場合にはマスクする
      cr  = VPn - v_ref
      cl  = VPs - v_ref
      acr = abs(cr)
      acl = abs(cl)
      
      dv4 = Un2-Un1
      dv3 = Un1-Up0
      dv2 = Up0-Us1
      dv1 = Us1-Us2
      
      include 'muscl.h'
      
      Urr = Un1 - (cm1*g6+cm2*g5)*ss_4 * lmt_n
      Url = Up0 + (cm1*g3+cm2*g4)*ss_4
      Ulr = Up0 - (cm1*g4+cm2*g3)*ss_4
      Ull = Us1 + (cm1*g1+cm2*g2)*ss_4 * lmt_s
      fu_r = 0.5*(cr*(Urr+Url) - acr*(Urr-Url)) * w_n
      fu_l = 0.5*(cl*(Ulr+Ull) - acl*(Ulr-Ull)) * w_s

      dv4 = Vn2-Vn1
      dv3 = Vn1-Vp0
      dv2 = Vp0-Vs1
      dv1 = Vs1-Vs2
      
      include 'muscl.h'
      
      Vrr = Vn1 - (cm1*g6+cm2*g5)*ss_4 * lmt_n
      Vrl = Vp0 + (cm1*g3+cm2*g4)*ss_4
      Vlr = Vp0 - (cm1*g4+cm2*g3)*ss_4
      Vll = Vs1 + (cm1*g1+cm2*g2)*ss_4 * lmt_s
      fv_r = 0.5*(cr*(Vrr+Vrl) - acr*(Vrr-Vrl)) * w_n
      fv_l = 0.5*(cl*(Vlr+Vll) - acl*(Vlr-Vll)) * w_s

      dv4 = Wn2-Wn1
      dv3 = Wn1-Wp0
      dv2 = Wp0-Ws1
      dv1 = Ws1-Ws2
      
      include 'muscl.h'
      
      Wrr = Wn1 - (cm1*g6+cm2*g5)*ss_4 * lmt_n
      Wrl = Wp0 + (cm1*g3+cm2*g4)*ss_4
      Wlr = Wp0 - (cm1*g4+cm2*g3)*ss_4
      Wll = Ws1 + (cm1*g1+cm2*g2)*ss_4 * lmt_s
      fw_r = 0.5*(cr*(Wrr+Wrl) - acr*(Wrr-Wrl)) * w_n
      fw_l = 0.5*(cl*(Wlr+Wll) - acl*(Wlr-Wll)) * w_s
      
      ! 流束の加算　VBCでない面の寄与のみを評価する
      cnv_u = cnv_u + fu_r*c_n - fu_l*c_s
      cnv_v = cnv_v + fv_r*c_n - fv_l*c_s
      cnv_w = cnv_w + fw_r*c_n - fw_l*c_s
			
      ! Z方向 ---------------------------------------
      
      ! 壁面の場合の参照速度の修正
      if ( (b_t2 == 0)  ) then
        Ut2 = wm1*u_ref + v(1, i  ,j  ,k+1)*wm2
        Vt2 = wm1*v_ref + v(2, i  ,j  ,k+1)*wm2
        Wt2 = w_ref2    - v(3, i  ,j  ,k+1)
        !Ut2 = u_ref2 - v(1, i  ,j  ,k+1)
        !Vt2 = v_ref2 - v(2, i  ,j  ,k+1)
      end if
      
      tmp1 = wm1*u_ref + Up0*wm2
      tmp2 = wm1*v_ref + Vp0*wm2
      
      if ( b_t1 == 0 ) then
        Ut1 = tmp1
        Vt1 = tmp2
        Wt1 = wq
        !Ut1 = uq
        !Vt1 = vq
      end if
      
      if ( b_b1 == 0 ) then
        Ub1 = tmp1
        Vb1 = tmp2
        Wb1 = wq
        !Ub1 = uq
        !Vb1 = vq
      end if
      
      if ( (b_b2 == 0)  ) then
        Ub2 = wm1*u_ref + v(1, i  ,j  ,k-1)*wm2
        Vb2 = wm1*v_ref + v(2, i  ,j  ,k-1)*wm2
        Wb2 = w_ref2    - v(3, i  ,j  ,k-1)
        !Ub2 = u_ref2 - v(1, i  ,j  ,k-1)
        !Vb2 = v_ref2 - v(2, i  ,j  ,k-1)
      end if
      
      ! 流束　流体のみと固体壁の影響を含み，隣接セルが固体の場合にはマスクする
      cr  = WPt - w_ref
      cl  = WPb - w_ref
      acr = abs(cr)
      acl = abs(cl)

      dv4 = Ut2-Ut1
      dv3 = Ut1-Up0
      dv2 = Up0-Ub1
      dv1 = Ub1-Ub2
      
      include 'muscl.h'
      
      Urr = Ut1 - (cm1*g6+cm2*g5)*ss_4 * lmt_t
      Url = Up0 + (cm1*g3+cm2*g4)*ss_4
      Ulr = Up0 - (cm1*g4+cm2*g3)*ss_4
      Ull = Ub1 + (cm1*g1+cm2*g2)*ss_4 * lmt_b
      fu_r = 0.5*(cr*(Urr+Url) - acr*(Urr-Url)) * w_t
      fu_l = 0.5*(cl*(Ulr+Ull) - acl*(Ulr-Ull)) * w_b

      dv4 = Vt2-Vt1
      dv3 = Vt1-Vp0
      dv2 = Vp0-Vb1
      dv1 = Vb1-Vb2
      
      include 'muscl.h'
            
      Vrr = Vt1 - (cm1*g6+cm2*g5)*ss_4 * lmt_t
      Vrl = Vp0 + (cm1*g3+cm2*g4)*ss_4
      Vlr = Vp0 - (cm1*g4+cm2*g3)*ss_4
      Vll = Vb1 + (cm1*g1+cm2*g2)*ss_4 * lmt_b
      fv_r = 0.5*(cr*(Vrr+Vrl) - acr*(Vrr-Vrl)) * w_t
      fv_l = 0.5*(cl*(Vlr+Vll) - acl*(Vlr-Vll)) * w_b

      dv4 = Wt2-Wt1
      dv3 = Wt1-Wp0
      dv2 = Wp0-Wb1
      dv1 = Wb1-Wb2
      
      include 'muscl.h'
            
      Wrr = Wt1 - (cm1*g6+cm2*g5)*ss_4 * lmt_t
      Wrl = Wp0 + (cm1*g3+cm2*g4)*ss_4
      Wlr = Wp0 - (cm1*g4+cm2*g3)*ss_4
      Wll = Wb1 + (cm1*g1+cm2*g2)*ss_4 * lmt_b
      fw_r = 0.5*(cr*(Wrr+Wrl) - acr*(Wrr-Wrl)) * w_t
      fw_l = 0.5*(cl*(Wlr+Wll) - acl*(Wlr-Wll)) * w_b
      
      ! 流束の加算　VBCでない面の寄与のみを評価する
      cnv_u = cnv_u + fu_r*c_t - fu_l*c_b
      cnv_v = cnv_v + fv_r*c_t - fv_l*c_b
      cnv_w = cnv_w + fw_r*c_t - fw_l*c_b

      
      ! 粘性項の計算　セル界面の剪断力を計算し，必要に応じて置換する
      uu_e = ( Ue1 - Up0 ) * dh2 ! 1/Re du/dx, e
      uu_w = ( Up0 - Uw1 ) * dh2
      uu_n = ( Un1 - Up0 ) * dh2
      uu_s = ( Up0 - Us1 ) * dh2
      uu_t = ( Ut1 - Up0 ) * dh2
      uu_b = ( Up0 - Ub1 ) * dh2
      
      vv_e = ( Ve1 - Vp0 ) * dh2
      vv_w = ( Vp0 - Vw1 ) * dh2
      vv_n = ( Vn1 - Vp0 ) * dh2
      vv_s = ( Vp0 - Vs1 ) * dh2
      vv_t = ( Vt1 - Vp0 ) * dh2
      vv_b = ( Vp0 - Vb1 ) * dh2
      
      ww_e = ( We1 - Wp0 ) * dh2
      ww_w = ( Wp0 - Ww1 ) * dh2
      ww_n = ( Wn1 - Wp0 ) * dh2
      ww_s = ( Wp0 - Ws1 ) * dh2
      ww_t = ( Wt1 - Wp0 ) * dh2
      ww_b = ( Wp0 - Wb1 ) * dh2 ! > (2*6)*3 = 36 flops
      
      if ( v_mode == 2 ) then ! 壁法則の場合の壁面摩擦による剪断応力の置換
        if ( ibits(bpx, facing_W, 6) /= 0 ) then ! 6面のうちのどれか方向フラグが立っている場合，つまり壁面隣接セル
          u1 = Up0 - u_ref
          u2 = Vp0 - v_ref
          u3 = Wp0 - w_ref
          ug = sqrt(u1*u1 + u2*u2 + u3*u3)
          
          if ( ug > 0.0 ) then
            e1 = abs(u1/ug)
            e2 = abs(u2/ug)
            e3 = abs(u3/ug)
          else
            e1 = 0.0
            e2 = 0.0
            e3 = 0.0
          endif
          
          u_tau = ut(i,j,k)*ut(i,j,k)
          
          if ( ibits(bpx, facing_E, 1) == 1 ) then ! ビットON -> 壁面
          ! uu_eは法線方向で，既に壁面条件で計算されている
            vv_e = u_tau * e2*e2
            ww_e = u_tau * e3*e3
          endif
          
          if ( ibits(bpx, facing_W, 1) == 1 ) then
            vv_w = u_tau * e2*e2
            ww_w = u_tau * e3*e3
          endif
          
          if ( ibits(bpx, facing_N, 1) == 1 ) then
            uu_n = u_tau * e1*e1
            ww_n = u_tau * e3*e3
          endif
          
          if ( ibits(bpx, facing_S, 1) == 1 ) then
            uu_s = u_tau * e1*e1
            ww_s = u_tau * e3*e3
          endif
          
          if ( ibits(bpx, facing_T, 1) == 1 ) then
            uu_t = u_tau * e1*e1
            vv_t = u_tau * e2*e2
          endif
          
          if ( ibits(bpx, facing_B, 1) == 1 ) then
            uu_b = u_tau * e1*e1
            vv_b = u_tau * e2*e2
          endif
          
          vflop = vflop + 67.0d0
          ! > 9 + sqrt*1 + /3 + 4*6 = 9+10+8*3+24 = 67 ! DP 9+20+13*3+24 = 92
        endif
      endif 
      
      beta = 1.0
      bdx = bd(i,j,k)
      
      EX = ( uu_e * c_e &
           - uu_w * c_w &
           + uu_n * c_n &
           - uu_s * c_s &
           + uu_t * c_t &
           - uu_b * c_b ) * dh1
           
      EY = ( vv_e * c_e &
           - vv_w * c_w &
           + vv_n * c_n &
           - vv_s * c_s &
           + vv_t * c_t &
           - vv_b * c_b ) * dh1
           
      EZ = ( ww_e * c_e &
           - ww_w * c_w &
           + ww_n * c_n &
           - ww_s * c_s &
           + ww_t * c_t &
           - ww_b * c_b ) * dh1 ! > 2*6*3 = 36 flops
			
      ! 対流項と粘性項の和 > 4*3 = 12 flops
      wv(1,i,j,k) = -cnv_u*dh1 + beta*EX*vcs
      wv(2,i,j,k) = -cnv_v*dh1 + beta*EY*vcs
      wv(3,i,j,k) = -cnv_w*dh1 + beta*EZ*vcs
    end do
    end do
    end do
!$OMP END DO
!$OMP END PARALLEL

    flop = flop + vflop

    return
    end subroutine pvec_muscl

!> ********************************************************************
!! @brief 次ステップのセルセンターの速度を更新
!! @param [out] v    n+1時刻の速度ベクトル
!! @param [out] div  \sum {u^{n+1}}
!! @param [in]  sz   配列長
!! @param [in]  g    ガイドセル長
!! @param [in]  dt   時間積分幅
!! @param [in]  dh   格子幅
!! @param [in]  vc   疑似速度ベクトル
!! @param [in]  p    圧力
!! @param [in]  bp   BCindex P
!! @param [in]  bv   BCindex V
!! @param [in]  v00  参照速度
!! @param [out] flop flop count
!! @note 
!!    - actvのマスクはSPEC_VEL/OUTFLOWの参照セルをマスクしないようにbvを使う
!!    - VBC(OUTFLOW, SPEC_VEL)の参照セルでは不定値となるが，InnerVBC()で上書き
!<
    subroutine update_vec (v, div, sz, g, dt, dh, vc, p, bp, bv, v00, flop)
    implicit none
    include 'ffv_f_params.h'
    integer                                                   ::  i, j, k, ix, jx, kx, g, bpx, bvx
    integer, dimension(3)                                     ::  sz
    double precision                                          ::  flop
    real                                                      ::  dh, dt, dd, actv, r_actv
    real                                                      ::  pc, px, py, pz, pxw, pxe, pys, pyn, pzb, pzt
    real                                                      ::  u_ref, v_ref, w_ref
    real                                                      ::  Ue0, Uw0, Vn0, Vs0, Wt0, Wb0, Up0, Vp0, Wp0
    real                                                      ::  Ue, Uw, Vn, Vs, Wt, Wb
    real                                                      ::  c1, c2, c3, c4, c5, c6
    real                                                      ::  N_e, N_w, N_n, N_s, N_t, N_b
    real                                                      ::  b_w, b_e, b_s, b_n, b_b, b_t
    real, dimension(3, 1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  v, vc
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  div, p
    real, dimension(0:3)                                      ::  v00
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bp, bv
    
    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    dd = dt/dh
    u_ref = v00(1)
    v_ref = v00(2)
    w_ref = v00(3)
    
!   flop = flop + dble(ix)*dble(jx)*dble(kx)*107.0 + 16.0d0
    flop = flop + dble(ix)*dble(jx)*dble(kx)*95.0 + 8.0d0 ! exclude cast ops
    ! flop = flop + dble(ix)*dble(jx)*dble(kx)*108.0 + 26.0d0 ! DP


!$OMP PARALLEL &
!$OMP PRIVATE(bpx, bvx, actv, r_actv) &
!$OMP PRIVATE(c1, c2, c3, c4, c5, c6) &
!$OMP PRIVATE(N_e, N_w, N_n, N_s, N_t, N_b) &
!$OMP PRIVATE(b_w, b_e, b_s, b_n, b_b, b_t) &
!$OMP PRIVATE(Ue0, Uw0, Vn0, Vs0, Wt0, Wb0, Up0, Vp0, Wp0) &
!$OMP PRIVATE(Ue, Uw, Vn, Vs, Wt, Wb) &
!$OMP PRIVATE(pc, px, py, pz, pxw, pxe, pys, pyn, pzb, pzt) &
!$OMP FIRSTPRIVATE(ix, jx, kx, dd, u_ref, v_ref, w_ref)

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
      bpx = bp(i,j,k)
      bvx = bv(i,j,k)
      actv = real(ibits(bvx, State,  1))
      r_actv = 1.0 - actv
      
      c1 = 1.0
      c2 = 1.0
      c3 = 1.0
      c4 = 1.0
      c5 = 1.0
      c6 = 1.0
      
      ! Neumann条件のとき，0.0 > 6 flop
      N_w = real(ibits(bpx, bc_n_W, 1))  ! w
      N_e = real(ibits(bpx, bc_n_E, 1))  ! e
      N_s = real(ibits(bpx, bc_n_S, 1))  ! s
      N_n = real(ibits(bpx, bc_n_N, 1))  ! n
      N_b = real(ibits(bpx, bc_n_B, 1))  ! b
      N_t = real(ibits(bpx, bc_n_T, 1))  ! t
      
      ! 疑似ベクトル
      Up0 = vc(1, i  ,j  ,k  )
      Vp0 = vc(2, i  ,j  ,k  )
      Wp0 = vc(3, i  ,j  ,k  )
      Uw0 = vc(1, i-1,j  ,k  )
      Ue0 = vc(1, i+1,j  ,k  )
      Vs0 = vc(2, i  ,j-1,k  )
      Vn0 = vc(2, i  ,j+1,k  )
      Wb0 = vc(3, i  ,j  ,k-1)
      Wt0 = vc(3, i  ,j  ,k+1)      
      
      ! 壁面による修正 (0-solid / 1-fluid) > 6 flop
      b_w = real( ibits(bv(i-1,j  ,k  ), State, 1) )
      b_e = real( ibits(bv(i+1,j  ,k  ), State, 1) )
      b_s = real( ibits(bv(i  ,j-1,k  ), State, 1) )
      b_n = real( ibits(bv(i  ,j+1,k  ), State, 1) )
      b_b = real( ibits(bv(i  ,j  ,k-1), State, 1) )
      b_t = real( ibits(bv(i  ,j  ,k+1), State, 1) )

      Uw = 0.5*( Up0 + Uw0 )*b_w + (1.0-b_w)*u_ref ! 36 flop
      Ue = 0.5*( Up0 + Ue0 )*b_e + (1.0-b_e)*u_ref
      Vs = 0.5*( Vp0 + Vs0 )*b_s + (1.0-b_s)*v_ref
      Vn = 0.5*( Vp0 + Vn0 )*b_n + (1.0-b_n)*v_ref
      Wb = 0.5*( Wp0 + Wb0 )*b_b + (1.0-b_b)*w_ref
      Wt = 0.5*( Wp0 + Wt0 )*b_t + (1.0-b_t)*w_ref
      
      ! c=0.0(VBC), 1.0(Fluid); VBCは内部と外部の両方
      if ( ibits(bvx, bc_face_W, bitw_5) /= 0 ) c1 = 0.0
      if ( ibits(bvx, bc_face_E, bitw_5) /= 0 ) c2 = 0.0
      if ( ibits(bvx, bc_face_S, bitw_5) /= 0 ) c3 = 0.0
      if ( ibits(bvx, bc_face_N, bitw_5) /= 0 ) c4 = 0.0
      if ( ibits(bvx, bc_face_B, bitw_5) /= 0 ) c5 = 0.0
      if ( ibits(bvx, bc_face_T, bitw_5) /= 0 ) c6 = 0.0
      
      ! 圧力勾配 18flop
      pc  = p(i,  j,  k  )
      pxw = (pc - p(i-1,j  ,k  )) * N_w
      pxe = (p(i+1,j  ,k  ) - pc) * N_e
      pys = (pc - p(i  ,j-1,k  )) * N_s
      pyn = (p(i  ,j+1,k  ) - pc) * N_n
      pzb = (pc - p(i  ,j  ,k-1)) * N_b
      pzt = (p(i  ,j  ,k+1) - pc) * N_t
      px = 0.5*(pxe + pxw)
      py = 0.5*(pyn + pys)
      pz = 0.5*(pzt + pzb)

      ! 発散値 VBCの寄与と壁面の影響は除外 25flop
      div(i,j,k) = ((Ue - dd * pxe) * c2 &
                  - (Uw - dd * pxw) * c1 &
                  + (Vn - dd * pyn) * c4 &
                  - (Vs - dd * pys) * c3 &
                  + (Wt - dd * pzt) * c6 &
                  - (Wb - dd * pzb) * c5 ) * actv
      
      ! セルセンタの速度更新 15flop
      v(1,i,j,k) = ( Up0-px*dd )*actv + r_actv*u_ref
      v(2,i,j,k) = ( Vp0-py*dd )*actv + r_actv*v_ref
      v(3,i,j,k) = ( Wp0-pz*dd )*actv + r_actv*w_ref

    end do
    end do
    end do
!$OMP END DO
!$OMP END PARALLEL
    
    return
    end subroutine update_vec

!> ********************************************************************
!! @brief 速度の発散に使う\sum{u^*}を計算する
!! @param [out] div  速度の発散値
!! @param [in]  sz   配列長
!! @param [in]  g    ガイドセル長
!! @param [in]  v0   疑似ベクトル
!! @param [in]  bv   BCindex V
!! @param [in]  v00  参照速度
!! @param [out] flop flop count
!<
    subroutine divergence (div, sz, g, v0, bv, v00, flop)
    implicit none
    include 'ffv_f_params.h'
    integer                                                   ::  i, j, k, ix, jx, kx, g, bvx
    integer, dimension(3)                                     ::  sz
    double precision                                          ::  flop
    real                                                      ::  Ue, Uw, Vn, Vs, Wt, Wb, actv
    real                                                      ::  c_e, c_w, c_n, c_s, c_t, c_b
    real                                                      ::  u_ref, v_ref, w_ref
    real                                                      ::  Ue0, Uw0, Vn0, Vs0, Wt0, Wb0, Up0, Vp0, Wp0
    real                                                      ::  b_w, b_e, b_s, b_n, b_b, b_t
    real, dimension(3, 1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  v0
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  div
    real, dimension(0:3)                                      ::  v00
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bv

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    u_ref = v00(1)
    v_ref = v00(2)
    w_ref = v00(3)
    
    ! loop : 1 + 42 + 13
!   flop  = flop + dble(ix)*dble(jx)*dble(kx)*55.0d0
    flop  = flop + dble(ix)*dble(jx)*dble(kx)*48.0d0 ! exclude cast ops

!$OMP PARALLEL &
!$OMP PRIVATE(bvx, actv) &
!$OMP PRIVATE(Ue0, Uw0, Vn0, Vs0, Wt0, Wb0, Up0, Vp0, Wp0) &
!$OMP PRIVATE(b_w, b_e, b_s, b_n, b_b, b_t) &
!$OMP PRIVATE(Ue, Uw, Vn, Vs, Wt, Wb) &
!$OMP PRIVATE(c_e, c_w, c_n, c_s, c_t, c_b) &
!$OMP FIRSTPRIVATE(ix, jx, kx, u_ref, v_ref, w_ref)

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
      bvx = bv(i,j,k)
      actv= real(ibits(bvx, State, 1))
      
      include 'd_o_o_p.h' ! 42flops
      
      ! 各面のVBCフラグ ibits() = 0(Normal) / others(BC) >> c_e = 1.0(Normal) / 0.0(BC)
      c_e = 1.0
      c_w = 1.0
      c_n = 1.0
      c_s = 1.0
      c_t = 1.0
      c_b = 1.0
      if ( ibits(bvx, bc_face_E, bitw_5) /= 0 ) c_e = 0.0
      if ( ibits(bvx, bc_face_W, bitw_5) /= 0 ) c_w = 0.0
      if ( ibits(bvx, bc_face_N, bitw_5) /= 0 ) c_n = 0.0
      if ( ibits(bvx, bc_face_S, bitw_5) /= 0 ) c_s = 0.0
      if ( ibits(bvx, bc_face_T, bitw_5) /= 0 ) c_t = 0.0
      if ( ibits(bvx, bc_face_B, bitw_5) /= 0 ) c_b = 0.0
      
      ! VBC面の影響をフラグで無効化 >> OBC_SPEC_VEL, OBC_WALL, OBC_SYMMETRIC, OBC_OUTFLOW  13flops
      div(i,j,k) = ( Ue*c_e - Uw*c_w + Vn*c_n - Vs*c_s + Wt*c_t - Wb*c_b ) * actv
    end do
    end do
    end do
!$OMP END DO
!$OMP END PARALLEL

    return
    end subroutine divergence
    
!> ********************************************************************
!! @brief 疑似ベクトルの時間積分（Euler陽解法）
!! @param[in,out] vc 疑似ベクトル
!! @param sz 配列長
!! @param g ガイドセル長
!! @param dt 時間積分幅
!! @param v 速度ベクトル（n-step, collocated）
!! @param bd BCindex ID
!! @param[out] flop
!! @note ここのマスクはIDのこと，VSPEC, OUTFLOWの増分をキャンセルするため
!<
    subroutine euler_explicit (vc, sz, g, dt, v, bd, flop)
    implicit none
    include 'ffv_f_params.h'
    integer                                                   ::  i, j, k, ix, jx, kx, g
    integer, dimension(3)                                     ::  sz
    double precision                                          ::  flop
    real                                                      ::  actv, dt
    real, dimension(3, 1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  vc, v
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bd
    
    ix = sz(1)
    jx = sz(2)
    kx = sz(3)

!   flop = flop + dble(ix)*dble(jx)*dble(kx)*8.0d0
    flop = flop + dble(ix)*dble(jx)*dble(kx)*7.0d0 ! exclude cast ops

!$OMP PARALLEL &
!$OMP PRIVATE(actv) &
!$OMP FIRSTPRIVATE(ix, jx, kx, dt)

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
      actv = dt * real(ibits(bd(i,j,k), State, 1))

      vc(1,i,j,k) = v(1,i,j,k) + vc(1,i,j,k)* actv
      vc(2,i,j,k) = v(2,i,j,k) + vc(2,i,j,k)* actv
      vc(3,i,j,k) = v(3,i,j,k) + vc(3,i,j,k)* actv
    end do
    end do
    end do
!$OMP END DO
!$OMP END PARALLEL
    
    return
    end subroutine euler_explicit

