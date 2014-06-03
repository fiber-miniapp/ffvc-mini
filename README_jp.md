FFVC-MINI
=========

* version: 1.0.0 (based on FFVC 0.9.2)
* date: 2014/06/03
* contact: miniapp@riken.jp


FFVCについて
-------------

本ミニアプリは，東京大学生産技術研究所で開発が進められている三次元非定常非圧縮性熱流体シミュレーションプログラムFFVC(FrontFlow/Violet Cartesian)をベースにしている．
FFVCは，三次元非圧縮Navier-Stokes方程式を，直交等間隔格子上での有限体積法により離散化して計算する．
詳細についてはオリジナルFFVCのユーザガイド(ffvc_ug.pdf)を参照のこと．

オリジナルFFVC連絡先: 小野謙二<keno@riken.jp>


インストール
------------

本プログラムのコンパイル・実行には，OpenMPをサポートしたC++およびFortran90コンパイラ，MPIライブラリ，GNU Makeが必要である．

1. 本パッケージの入手とファイルの展開

2. srcディレクトリに移動し、"make_setting"ファイルの内容を実行環境に合わせて編集する．
   同ディレクトリには，以下の設定例が含まれる:
    * make_setting.intel : Intelコンパイラ
    * make_setting.fx10  : 富士通コンパイ(京/FX10)
    * make_setting.gcc   : GCCコンパイラ
    * make_setting.pgi   : PGIコンパイラ

3. srcディレクトリでmakeを実行．
    正しくmakeが完了すると、binディレクトリ下に実行プログラム`ffvc_mini`が作成される．


テスト
------

インタラクティブにMPIジョブが実行できる環境用に，簡単なテストスクリプトをtestディレクトリに用意してある．
テストを実行するには， testディレクトリでシェルスクリプトgo.shを実行するか，srcディレクトリで「make test」を実行する．
このテストスクリプトでは，8プロセス，2スレッドでの計算を行い，見本計算結果との間で流速データの比較を行い，計算結果の妥当性を確認する．


プログラム計算内容
------------------

本ミニアプリは，オリジナルFFVCのPMT(performance monitor test)組み込み例題を計算することに特化されている．
PMT組み込み例題では，三次元直方体領域でのキャビティフローを計算する．
計算は性能測定モードで行われ，反復計算での収束判定は行わす，反復回数は固定となる．
PMT組み込み例題については，オリジナルFFVCのユーザガイド(ffvc_ug.pdf)18頁を参照のこと．

オリジナルFFVCでは，設定可能な境界条件は，計算領域の周囲に位置する外部境界条件と，計算領域内に設置可能な内部境界条件からなる．内部境界条件には， 物質境界，湧き出し/吸い込み口などの他に，熱交換器のモデル化に対応可能な圧力損失境界条件も設定可能である．
本ミニアプリでは，計算対象を三次元キャビティフローに限定しているため，外部境界条件のみを使用している．そのため，未使用な境界条件に対応した関数・サブルーチンの多くは削除してある．
ただし，プログラムの制御構造は変更せず，オリジナルコードのV-P反復とPoisson反復からなる二重ループ構造(下図)を継承している．


    時間ステップループ {
        - 対流項，拡散項計算
        - Poissonソース項計算
  
        V-P反復ループ {
            Poisson反復ループ {
                - SOR計算
            }
  
            - 速度更新
            /*
            - 圧力損失境界条件によるPoissonソース項の更新
              (ミニアプリでは削除)
            */
        }
    }

Poisson反復ループでは，圧力Poisson方程式を反復解法で解く．
一方，V-P反復ループは，速度の関数として与えられた圧力損失境界条件を， Newton-Raphson法的に扱うためのループである．
そのため，圧力損失境界条件を使用しない本ミニアプリでは，実際にはV-P反復の中身は1回実行すれば充分である．
ただし，本ミニアプリでは，ベンチマークプログラムとしての使用を考慮して，V-P反復，Poisson反復とも，任意の反復回数を指定して実行可能である．

オリジナルFFVCがターゲットとしている計算では，典型的な反復回数は,
V-P反復は5〜100回，Poisson反復は20〜1000回程度である．

なお，ミニアプリ化にともない，使用アルゴリズムを，時間積分は1次精度Euler陽解法，対流項は3次精度MUSCLスキーム，反復解法はストライドメモリアクセス型の2色SOR法に固定している．


プログラム実行方法
------------------

### コマンドライン引数 ####

本プログラムには，入力ファイルはなく，実行時に必要なパラメータは全てコマンド
ラインで指定する．指定可能なコマンドラインオプションは以下のとおり．

    --scale=str     サイズ指定モードをstrongまたはweakで指定 (必須)
    --size=int      計算領域サイズを整数(一辺のセル数)で指定 (必須)
    --division=str  領域分割をLxMxN形式で指定 (ディフォルト 1x1x1)
    --step=int      計算ステップ数を整数で指定 (ディフォルト 20)
    --dt=float      時間刻み幅を実数(CFL数単位で)で指定 (ディフォルト 0.2)
    --comm=str      通信方法をsync(同期)またはasync(非同期)で指定
                    (ディフォルト async)
    --p_itr=int     Poisson反復の最大反復回数 (ディフォルト 30)
    --vp_itr=int    V-P反復の最大反復回数 (ディフォルト 20)
    --practical     practical計算フラグ (ディフォルト 未指定)
    --output_interval=int
                    出力ステップ間隔 (ディフォルト 0(出力なし))
    --help          ヘルプメッセージを表示して終了


#### サイズ指定モード ####

scaleパラメータに指定する文字列(strongまたはweak)は，それぞれ，strong scaling, weak scalingでの性能測定に対応する．
指定された文字列により，以下の例のように，sizeパラメータの値に対して異なった解釈がなされる．

以下のコマンドライン指定は，storong scalingでの性能計測に対応し，128x128x128セルを全計算領域として，それを2x2x2に分割して，8プロセスで計算する．

    $ ffvc_mini --scale=strong --size=128 --division=2x2x2

一方，以下のコマンドライン指定では，weak scalingでの性能計測に対応し，128x128x128セルの計算領域を基本単位として，それを2x2x2に積み上げた領域を，8プロセスで計算する．

    $ ffvc_mini --scale=weak --size=128 --division=2x2x2


#### practical計算フラグ ####

ディフォルト(practicalフラグ指定なし)では，性能測定モードで計算が行われる．
性能測定モードでは，毎時間ステップにおいて，Poisson反復，V-P反復とも，収束状況によらず常に最大反復回数の反復がなされる．
practical計算フラグが指定されると，実際の流体シミュレーションに対応した計算が行われる．その場合，毎時間ステップの反復計算において，収束条件が満たされた時点で反復を終了する．


### 出力ファイル ####

各時間ステップでの，反復収束状況，物理量の統計値，計算時間が，ファイルhistory_base.txtに出力される．

また，output_intervalオプションに1以上の整数が指定されると，そのステップ間隔で，圧力と速度の瞬時値が，出力される．
出力は，各時間ステップ毎，計算ノード毎に，別ファイルに出力される．
この時，出力ファイル名は以下のとおりとなる．

  * prs_.dfi   ... 圧力値分散出力でのインデックスファイル(テキスト)
  * vel_.dfi   ... 速度値分散出力でのインデックスファイル(テキスト)
  * prs_TTTTTTTTTT_idNNNNNN.sph  ... 圧力値(バイナリ)
  * vel_TTTTTTTTTT_idNNNNNN.sph  ... 圧力値(バイナリ)

ここで，TTTTTTTTTTは10桁整数で表した時間ステップ番号，NNNNNNは6桁整数で表したノード番号である．


エクサスケールでの想定計算規模
------------------------------

一例として，1mm幅の格子，1000億セル、実時間3秒の計算を想定している．
この系に対して，ディフォルトの反復回数の下で，100万ステップの単精度計算を
3時間程度で完了したい．


オリジナルFFVCからの変更点
--------------------------

* 機能を限定することにより，コードを削減

    - PMT組み込み例題(3次元キャビティフローでの性能測定)の実行に限定

    - 流体解法アルゴリズムをFS_EE_EE(Fractional Step, Euler Explicit)に限定

    - 反復解法アルゴリズムをSOR2SMA(ストライドメモリアクセスの2色SOR)に限定

    - モニタリング，サンプリング機能は，history_base.txt相当を残して削除

    - その他の機能限定，パラメータのハードコード

* 実行時パラメータ指定方法を，設定ファイルからコマンドライン引数に変更

* 未使用なライブラリ，コード，変数の削除

* その他のコード整備


詳細情報
--------

### ホットスポット ###

Poisson反復内でのSOR計算がホットスポットとなる．

  * psor2sma_core_     SOR計算コア部 (N_vp*N_p*2回)
  * poi_residual_      残差計算 (N_vp*N_p回)
  * sma_comm_, sma_comm_wait_    袖通信 (N_vp*N_p*2回)

カッコ内は，V-P反復，Poisson反復，それぞれの反復回数をN_vp, N_pとした時の，毎時間ステップでの実行回数である．

また，実際の流体計算では，反復が早く収束して最大反復回数よりかなり少ない回数で反復を抜けるような状況も発生する．そのような場合には，Poisson反復の外で実行される

  * update_vec_    圧力勾配計算と速度の更新 (N_vp回)

や，V-P反復の外で実行される

  * pvec_muscl_    対流項と粘性項の計算 (1回)

などのサブルーチンもホットスポットとなる．


### ビット情報配列 ###

FFVCでは，各セルでの，媒質情報，境界条件情報，隣接セル情報などは，
ビット単位のフラグとして32ビット整数内にエンコードし，これらの配列として保持
している．プログラム中では，以下の3つのビット情報配列が使用される．

 * d_bcd  媒質情報
 * d_bcp  圧力計算に必要な情報
 * d_bcv  速度計算に必要な情報

連立一次方程式の係数行列などは，実際の計算時にこれらのビット情報配列よりデコードされる．


### 単精度 vs. 倍精度 ###

本ミニアプリでは，浮動小数点演算は単精度で行っている．
非圧縮流体計算では，精度的に単精度計算で充分と考えられる場合がある．

本ミニアプリを倍精度計算で評価するには，srcディレクトリのmake_settingにおいて，マクロ`CXXFLAGS`に`-D_REAL_IS_DOUBLE_`を，`F90FLAGS`に基本実数型を倍精度実数型とするコンパイラオプション(Intelコンパイラなら`-r8`，京/FX10では`-CcdRR8`)を
追加し，実行プログラムをビルドする．