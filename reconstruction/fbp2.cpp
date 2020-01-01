#include<bits/stdc++.h>

using namespace std;

const char writeFileName_map[] = "projection_test\\testmpap2e5";
const char writeFileName_xy[] = "projection_test\\testxyo2e5";
const char writeFileName_project[] = "projection_test\\outproj2e5";

char wname_map[27];
char wname_xy[27];
char wname_project[27];

void writeRawFile(const char fname[], const size_t size, const size_t num, void* image);

int main(void){

  int H=65,W=65;
  int outH =256, outW=256;
  
  float*	projection=(float*)calloc(H, sizeof(float));	/** 投影画像用配列 **/
  float*	projection3d=(float*)calloc(H*360, sizeof(float));
  float*	projection_w=(float*)calloc(H*360, sizeof(float));	/** 重み付け画像用配列 **/
  float*	projection_out=(float*)calloc(H * 360, sizeof(float));	/** 出力画像用配列 **/
  FILE	*fpi;	/** ファイルポインタ **/

  //3あり画像か無しかは要確認
  fpi = fopen("map5_20_2e5.raw" , "rb");
  fread(projection3d, sizeof(float), 65*360, fpi);

  cout<<endl<<endl;
  
  float* image_xy=(float*)calloc(outH * outW, sizeof(float));	/* 出力画像用配列 */

//step1 重み付け
#if 1
for(int b = 0; b<360;b++){
  for(int zeta = 0; zeta < H; zeta++){
    projection_w[b*W + zeta] = projection3d[b*W + zeta]*(60./(sqrt(pow(60.,2)+pow(-16.25+0.5*zeta,2))));
  }
}

#endif

  //step2 Ramp filterによる補正
  //rampf filter配列作成
#if 1
  //doubleに変更//////////////////////////////////////////////
  double* ramp =(double*)calloc(W*2 - 1, sizeof(double));
  ramp[H-1] = 0.25; //1/4を代入すると0が入る事に注意

  for(int n = 1; n<W; n++){
    if(n%2 == 0){
      ramp[H-1 + n] = 0;
      ramp[H-1 - n] = 0;
    }           
    else{
      ramp[H-1 + n] = -1./pow(n*M_PI,2);
      ramp[H-1 - n] = -1./pow(n*M_PI,2);
    }
  }

  for(int a = 0;a<H;a++){
    //cout<<a<<" "<<ramp[a]<<" "<<128-a<<" "<<ramp[128-a]<<endl;
  }

  //重畳積分
  float tmp1 = 0;
  for(int d=0;d<360;d++){//縦移動
    for(int b = 0; b<W; b++){//出力の注目画素
      double tmp = 0.;
      for(int c = 0; c<W; c++){//横にスライドしてそれぞれ掛けて足す
        tmp += projection_w[d * H + c]*0.5*ramp[H-1-b+c];
        //if(a == 0 && c == 32){cout<< c <<" "<<64-b+c<<endl;}
      }
      projection_out[d*H + b] = tmp;
       //map[a*H + b]+=tmp;
      if(b == 32){
        tmp1 += projection_out[d*H + b]; 
      }
    }
  }

  cout<<tmp1<<" ";

  cout<<"filtered"<<endl;


#endif

#if 1
  //step3 逆投影
  int beta_max = 360;
  float beta_span = 1;
for(int beta = 1; beta < beta_max; beta+=beta_span){//角度
    //光源をxz平面に関して1度ごと360方向回転させる
    float start_x = -160;
    float start_y = 0;//?

    double primary_x = start_x*cos(M_PI*beta/180) - start_y*sin(M_PI*beta/180);
    double primary_y = start_x*sin(M_PI*beta/180) + start_y*cos(M_PI*beta/180);

    //25~29,34
      for(int t = 0; t < outH; t++){//高さ↑ y 15 to H-15
        for(int s = 0; s < outW; s++){//幅→ 右半分だけ再構成 x 32 to W

          double photon_vec[2];//再考空間における光子の座標が入るベクトル
          photon_vec[0] = (-12.8 + s*0.1) - primary_x;//x
          photon_vec[1] = (12.8 - t*0.1) - primary_y;//y

          double tmp_x = photon_vec[0];

          //逆回転させて0度の位置に持ってくる
          photon_vec[0] = photon_vec[0]*cos(-1*M_PI*beta/180) - photon_vec[1]*sin(-1*M_PI*beta/180);
          photon_vec[1] = tmp_x*sin(-1*M_PI*beta/180) + photon_vec[1]*cos(-1*M_PI*beta/180);

          //x座標がditectorの位置になるようにベクトルを拡大し光子を動かす
          double to_ditector = 220 / photon_vec[0];
          for(int pv_index = 0; pv_index < 2; pv_index++){
            photon_vec[pv_index] *= to_ditector;
          }

          //y,zが検出器外ならこの先処理しない
          //z=32の時の球領域より外側のマイナス部分の値が入っていない．
          //ここが怪しそう
          //16.25, 23.25を小数点無しにずらした影響もありそう
          if(abs(photon_vec[1]) > 16.25){
            //cout<<photon_vec[1]<<" ";
            continue;
          }

          //検出器座標
          int index_y = -2 * (((photon_vec[1]) - 16.25));

          //////////////////////ここから下をcheck///////////////////////////////////////

          //再構成空間における各点をbetaに従い回転させる
          //角度マイナスか?
          double tmp_s =(-12.8+s*0.1)*cos(M_PI*beta/180) - (12.8-0.1*t)*sin(M_PI*beta/180);
          double tmp_t = (-12.8+s*0.1)*sin(M_PI*beta/180) + (12.8-0.1*t)*cos(M_PI*beta/180);
          //s軸を表す方程式は y = tan(beta) * x なので点と直線の距離公式から
          double d = abs(-1*tan(beta)*tmp_t + 1*tmp_s)/(sqrt(1 + pow(tan(beta),2)));

          //s軸よりも右側の場合はマイナスにする
          if(tmp_s<0){
            d *= -1; 
          }
          //検出器より後ろの場合は再構成領域に加算しない

          //if((t-128)*(t-128)+(s-128)*(s-128)<=118*118){//380
            //球領域にのみ逆投影100
            //dbetaをかける
            //Dはオブジェクト中心 (pow(22.5,2)/pow(22.5 - d,2)) * 
            double tmp_output = (pow(60,2)/pow(60 - d,2)) * projection_out[(beta)*H + index_y] * beta_span * 2 * M_PI/360;
            //double((pow(22.5,2)/pow(22.5 - d,2)) * 
            image_xy[outH*t + s] += tmp_output*1.7;//(t-1)->t
          //}
       }
     }
   }
//}
#endif

  sprintf(wname_project,"%s.raw",writeFileName_project);
  writeRawFile(wname_project, sizeof(float), H*360, projection_out);
  FILE* fp = fopen(wname_project, "wb");
  fwrite(projection_out , sizeof(float), H*360, fp);

  sprintf(wname_xy,"%s.raw",writeFileName_xy);
  writeRawFile(wname_xy, sizeof(float), outH * outW, image_xy);
  FILE* fp_xy = fopen(wname_xy, "wb");
  fwrite(image_xy , sizeof(float), outH * outW, fp_xy);
}

void writeRawFile(const char fname[], const size_t size, const size_t num, void* image)
{
    // ファイルを開く
    FILE* fp = fopen(fname, "wb");

    // ファイルを開くことができなかった場合のエラー処理
    if ( NULL == fp )
    {
        printf("failed to open %s\n", fname);
        exit(-1);
    }

    // データの書き出し
    size_t ret = fwrite(image, size, num, fp);

    // データを書き込むことができなかった場合のエラー処理
    if ( num != ret )
    {
        printf("failed to write %s\n", fname);
        fclose(fp);
        exit(-1);
    }

    // ファイルを閉じる
    fclose(fp);
}