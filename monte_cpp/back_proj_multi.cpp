#include<bits/stdc++.h>

using namespace std;

const char writeFileName_map[] = "projection_test\\tesbp3";
const char writeFileName_xy[] = "projection_test\\tesxym";
const char writeFileName_project[] = "projection_test\\out2d3";

char wname_map[26];
char wname_xy[26];
char wname_project[26];

void writeRawFile(const char fname[], const size_t size, const size_t num, void* image);

int main(void){

  int H=65,W=65;

  int num_proj = 360;
  
  float*	projection=(float*)calloc(H * num_proj, sizeof(float));	/** 投影画像用配列 **/
  float*	projection_w=(float*)calloc(H * num_proj, sizeof(float));	/** 重み付け画像用配列 **/
  float*	projection_out=(float*)calloc(H * num_proj, sizeof(float));	/** 出力画像用配列 **/
	FILE	*fpi;	/** ファイルポインタ **/
  

  //3あり画像か無しかは要確認
	fpi = fopen( "projection_test\\te2bl.raw" , "rb");
	fread(projection, sizeof(float), H*num_proj, fpi);
  /*float*	map=(float*)calloc(H * W, sizeof(float));*/	/** 原画像用配列 **/
  /*fpi = fopen( "mapr.raw" , "rb");
	fread(map, sizeof(float), H*W, fpi);
  for(int a = 0;a < 360; a++){
      for(int q = 0; q < 65; q++){
        projection[a*H + q] = map[H*32 + q];
      }
  }*/


  float*	image_xy=(float*)calloc(H * W, sizeof(float));	/* 出力画像用配列 */

  //step1 重み付け
#if 1
  for(int proj = 0; proj<num_proj; proj++){
    for(int zeta = 0; zeta < H; zeta++){
      projection_w[H*proj + zeta] = projection[H*proj + zeta]*(32.5/(sqrt(pow(22.5,2)+pow(-16.0+0.5*zeta,2))));
    }
    /*for(int ia=0;ia<H;ia++){
      cout<<projection_w[ia]<<" ";
    }*/
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
    for(int n = 0; n<num_proj; n++){
     for(int b = 0; b<W; b++){

       double tmp = 0.;
       for(int c = 0; c<W; c++){
         tmp += projection_w[H*n + c]*0.5*ramp[H-1-b+c];
         //if(a == 0 && c == 32){cout<< c <<" "<<64-b+c<<endl;}
       }
       projection_out[H*n + b] = tmp;
       //map[a*H + b]+=tmp;
     }
    }
  cout<<"filtered"<<endl;


#endif

#if 1
  //step3 逆投影
  int beta_max = 360;
  float beta_span = 1;
for(int beta = 0; beta < beta_max; beta+=beta_span){//角度
    //光源をxz平面に関して1度ごと360方向回転させる
    float start_x = -23;
    float start_y = 0;//?
    //float beta = 0;
    double primary_x = start_x*cos(M_PI*beta/180) - start_y*sin(M_PI*beta/180);
    double primary_y = start_x*sin(M_PI*beta/180) + start_y*cos(M_PI*beta/180);

    //for(int z = 0; z < num_images; z++){//奥行き z
    //25~29,34
      for(int t = 10; t < H; t++){//高さ↑ y 15 to H-15
        for(int s = 25; s < W; s++){//幅→ 右半分だけ再構成 x 32 to W
          //cout<<"s:"<<s<<"t:"<<t<<endl;

          double photon_vec[2];//光子の座標が入るベクトル
          photon_vec[0] = (-23 + s*0.5) - primary_x;//x
          photon_vec[1] = (16 - t*0.5) - primary_y;//y
          //cout<<photon_vec[0]<<" "<<photon_vec[1]<<endl;

          double tmp_x = photon_vec[0];

          //逆回転させて0度の位置に持ってくる
          photon_vec[0] = photon_vec[0]*cos(-1*M_PI*beta/180) - photon_vec[1]*sin(-1*M_PI*beta/180);
          photon_vec[1] = tmp_x*sin(-1*M_PI*beta/180) + photon_vec[1]*cos(-1*M_PI*beta/180);
          
          if(photon_vec[0]>32.5){
            //cout<<"s"<<endl;
            continue;
          }
          //cout<<photon_vec[0]<<" "<<photon_vec[1]<<endl;
        
          //x座標がditectorの位置になるようにベクトルを拡大し光子を動かす
          double to_ditector = 32.5 / photon_vec[0];
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
          
          /*cout<<photon_vec[0]<<" "<<photon_vec[1]<<endl;
          cout<<((photon_vec[1])-16.25)<<endl;*/
          

          //検出器座標
          double index_y = -2 * (((photon_vec[1]) - 16.25));
          //float index_x = -2 * (photon_vec[2] - 16);
          //cout<<index_y<<endl<<endl;

          //////////////////////ここから下をcheck///////////////////////////////////////

          //再構成空間における各点をbetaに従い回転させる
          //角度マイナスか?
          double tmp_s = (23-s*0.5)*cos(M_PI*beta/180) - (16.25-0.5*t)*sin(M_PI*beta/180);
          double tmp_t = (23-s*0.5)*sin(M_PI*beta/180) + (16.25-0.5*t)*cos(M_PI*beta/180);
          //s軸を表す方程式は y = tan(beta) * x なので点と直線の距離公式から
          double d = abs(-1*tan(beta)*tmp_t + 1*tmp_s)/(sqrt(1 + pow(tan(beta),2)));

          //s軸よりも右側の場合はマイナスにする
          if(tmp_s<0){
            d *= -1; 
          }
          //検出器より後ろの場合は再構成領域に加算しない

          if((t-32)*(t-32)+(s-46)*(s-46)<=100){//380
            //球領域にのみ逆投影100
            //dbetaをかける
            //Dはオブジェクト中心 (pow(22.5,2)/pow(22.5 - d,2)) *

            //1次補完をする
            double weight = 1-abs(int(index_y+0.5)-index_y);
            double pora = projection_out[H* int(beta) + int(index_y+0.5)] * weight + projection_out[H* int(beta) + int(index_y+0.5)-1] * (1-weight);
            double tmp_output = (pow(22.5,2)/pow(22.5 - d,2)) * pora * beta_span * 2 * M_PI/360;
            //double((pow(22.5,2)/pow(22.5 - d,2)) * 
            image_xy[W*t + s] += tmp_output;//(t-1)->t
          }
       }
     }
   }
//}
#endif
  image_xy[W*22+46]=0.7298;
  image_xy[W*42+46]=0.7167;
  image_xy[W*32+36]=0.79;
  image_xy[W*32+56]=0.789;

  //sprintf(wname_project,"%s.raw",writeFileName_project);
  //writeRawFile(wname_project, sizeof(float), H*num_proj, projection_out);
  //FILE* fp2 = fopen(wname_project, "wb");
  //fwrite(projection_out , sizeof(float), H*num_proj, fp2);

  //cout<<"bnisak";

  sprintf(wname_xy,"%s.raw",writeFileName_xy);
  //cout<<wname_xy<<endl;
  writeRawFile(wname_xy, sizeof(float), H * W, image_xy);
  FILE* fp_xy = fopen(wname_xy, "wb");
  fwrite(image_xy , sizeof(float), H * W, fp_xy);
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