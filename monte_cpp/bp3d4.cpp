#include<bits/stdc++.h>

using namespace std;

const char writeFileName_map[] = "CBCTtest3\\bprmap";
const char writeFileName_xy[] = "CBCTtest3\\bpr_xy";
const char writeFileName_xz[] = "CBCTtest3\\bpr_xz";

char wname_map[20];
char wname_xy[20];
char wname_xz[20];

void writeRawFile(const char fname[], const size_t size, const size_t num, void* image);

int main(void){

  int H=65,W=65,num_images=65;
  
  float*	map=(float*)calloc(H * W, sizeof(float));	/** 原画像用配列 **/
  float*	map_w=(float*)calloc(H * W, sizeof(float));	/** 重み付け画像用配列 **/
  float*	map_out=(float*)calloc(H * W, sizeof(float));	/** 出力画像用配列 **/
	FILE	*fpi;	/** ファイルポインタ **/

  //3あり画像か無しかは要確認
	fpi = fopen( "mapr.raw" , "rb");
	fread(map, sizeof(float), H*W, fpi);

  float*	image_xy=(float*)calloc(H * W * num_images, sizeof(float));	/* 出力画像用配列 */
  //float*	image_xz=(float*)calloc(H * W * num_images, sizeof(float));	/* 出力画像用配列 */

  //step1 重み付け
  for(int zeta = 0; zeta < H; zeta++){
    for(int p = 0; p < W; p++){
      //dlで割る,map作成時点でこれはやっておきたい
      //Dso->Dsdに変更
      //16.25->16.
      map_w[H*zeta + p] = map[H*zeta + p]*(32.5/sqrt(pow(22.5,2)+pow(-1*(zeta*0.5)+16.0,2)+pow((p*0.5)-16.0,2)))/0.5;
      //if(pow((zeta-32),2)+pow((p-32),2)<225){map_w[zeta*H + p]+=(22.5/sqrt(pow(22.5,2)+pow(-1*(zeta*0.5)+16.25,2)+pow((p*0.5)-16.25,2)))*1./0.5;}
      //if(pow((zeta-32),2)+pow((p-32),2)<9){map_w[zeta*H + p]=(22.5/sqrt(pow(22.5,2)+pow(-1*(zeta*0.5)+16.25,2)+pow((p*0.5)-16.25,2)))*2./0.5;}
    }
  }

  //step2 Ramp filterによる補正
  //rampf filter配列作成
#if 1
  float* ramp =(float*)calloc(W*2 - 1, sizeof(float));
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
  for(int a = 0; a<H; a++){
     for(int b = 0; b<W; b++){
       float tmp = 0.;
       for(int c = 0; c<W; c++){
         tmp += map_w[a*H + c]*0.5*ramp[H-1-b+c];
         //if(a == 0 && c == 32){cout<< c <<" "<<64-b+c<<endl;}
       }
       map_out[a*H + b] = tmp;
       //map[a*H + b]+=tmp;
     }
  }
  cout<<"filtered"<<endl;

#endif

#if 0
  //step3 逆投影
  int beta_max = 360;
  float beta_span = 1;
for(float beta = 0; beta < beta_max; beta+=beta_span){//角度
    //光源をxz平面に関して1度ごと360方向回転させる
    float start_x = -23;
    float start_y = 0;//?
    //float beta = 0;
    float primary_x = start_x*cos(M_PI*beta/180) - start_y*sin(M_PI*beta/180);
    float primary_y = start_x*sin(M_PI*beta/180) + start_y*cos(M_PI*beta/180);

    for(int z = 0; z < num_images; z++){//奥行き z
      for(int t = 1; t < H; t++){//高さ↑ y 15 to H-15
        for(int s = 2; s < W; s++){//幅→ 右半分だけ再構成 x 32 to W

          float photon_vec[3];//光子の座標が入るベクトル
          photon_vec[0] = (-23 + s*0.5) - primary_x;//x
          photon_vec[1] = (16 - t*0.5) - primary_y;//y
          photon_vec[2] = (-16 + z*0.5);//z

          float tmp_x = photon_vec[0];

          //逆回転させて0度の位置に持ってくる
          photon_vec[0] = photon_vec[0]*cos(-1*M_PI*beta/180) - photon_vec[1]*sin(-1*M_PI*beta/180);
          photon_vec[1] = tmp_x*sin(-1*M_PI*beta/180) + photon_vec[1]*cos(-1*M_PI*beta/180);
          
          if(photon_vec[0]>32.5){
            //cout<<"s"<<endl;
            continue;
          }
          //x座標がditectorの位置になるようにベクトルを拡大し光子を動かす
          float to_ditector = 32.5 / photon_vec[0];
          for(int pv_index = 0; pv_index < 3; pv_index++){
            photon_vec[pv_index] *= to_ditector;
          }

          //y,zが検出器外ならこの先処理しない
          //z=32の時の球領域より外側のマイナス部分の値が入っていない．
          //ここが怪しそう
          //16.25, 23.25を小数点無しにずらした影響もありそう
          if(abs(photon_vec[1]) > 16 || abs(photon_vec[2]) > 16){
            //cout<<photon_vec[1]<<" ";
            continue;
          }

          //検出器座標
          float index_y = -2 * (photon_vec[1] - 16);
          float index_x = -2 * (photon_vec[2] - 16);

          //////////////////////ここから下をcheck///////////////////////////////////////

          //再構成空間における各点をbetaに従い回転させる
          //角度マイナスか?
          float tmp_s = (23-s*0.5)*cos(M_PI*beta/180) - (16-0.5*t)*sin(M_PI*beta/180);
          float tmp_t = (23-s*0.5)*sin(M_PI*beta/180) + (16-0.5*t)*cos(M_PI*beta/180);
          //s軸を表す方程式は y = tan(beta) * x なので点と直線の距離公式から
          float d = abs(-1*tan(beta)*tmp_t + 1*tmp_s)/(sqrt(1 + pow(tan(beta),2)));

          //s軸よりも右側の場合はマイナスにする
          if(tmp_s<0){
            d *= -1;
          }
          //検出器より後ろの場合は再構成領域に加算しない

          if((z-32)*(z-32)+(t-32)*(t-32)+(s-46)*(s-46)<=1000){//球領域にのみ逆投影100
            //dbetaをかける
            double tmp_output = double((pow(22.5,2)/pow(22.5 - d,2)) * map_out[H*(int)(index_y) + (int)(index_x)])* (M_PI/180) *beta_span; 
            image_xy[H*W*z + W*t + s] += tmp_output;
            image_xz[H*W*t + W*z + s] += tmp_output;
            //(pow(22.5,2)/pow(22.5 - d,2)) * 
            /*if(s==46&&z==32){
              cout<<t<<" "<<index_x<<" "<<index_y<<" "<<map_out[H*(int)(index_y+0.5) + (int)(index_x+0.5)]<<endl;
            }*/
          }
       }
     }
   }
}
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

    for(int z = 0; z < num_images; z++){//奥行き z
    //25~29,34
      for(int t = 10; t < H-15; t++){//高さ↑ y 15 to H-15
        for(int s = 25; s < W; s++){//幅→ 右半分だけ再構成 x 32 to W

          double photon_vec[3];//光子の座標が入るベクトル
          photon_vec[0] = 220;//x
          photon_vec[1] = (16 - t*0.5) - primary_y;//y
          photon_vec[2] = (-16 + z*0.5);//z

          double tmp_x = photon_vec[0];

          //逆回転させて0度の位置に持ってくる
          photon_vec[0] = photon_vec[0]*cos(-1*M_PI*beta/180) - photon_vec[1]*sin(-1*M_PI*beta/180);
          photon_vec[1] = tmp_x*sin(-1*M_PI*beta/180) + photon_vec[1]*cos(-1*M_PI*beta/180);
          
          if(photon_vec[0]>32.5){
            continue;
          }
        
          //x座標がditectorの位置になるようにベクトルを拡大し光子を動かす
          double to_ditector = 32.5 / photon_vec[0];
          for(int pv_index = 0; pv_index < 3; pv_index++){
            photon_vec[pv_index] *= to_ditector;
          }

          //y,zが検出器外ならこの先処理しない
          //z=32の時の球領域より外側のマイナス部分の値が入っていない．
          //ここが怪しそう
          //16.25, 23.25を小数点無しにずらした影響もありそう
          if(abs(photon_vec[1]) > 16.25 || abs(photon_vec[2]) > 16.25){
            continue;
          }
          

          //検出器座標
          //16.25->16
          double y = -2 * ((photon_vec[1]) - 16.25);
          double x = -2 * (photon_vec[2] - 16);
          int yi = int(y);
          int xi = int(x);


          //////////////////////ここから下をcheck///////////////////////////////////////

          //再構成空間における各点をbetaに従い回転させる
          //角度マイナスか?
          //16.25->16
          double tmp_s = (23-s*0.5)*cos(M_PI*beta/180) - (16.25-0.5*t)*sin(M_PI*beta/180);
          double tmp_t = (23-s*0.5)*sin(M_PI*beta/180) + (16.25-0.5*t)*cos(M_PI*beta/180);
          //s軸を表す方程式は y = tan(beta) * x なので点と直線の距離公式から
          double d = abs(-1*tan(beta)*tmp_t + 1*tmp_s)/(sqrt(1 + pow(tan(beta),2)));

          //s軸よりも右側の場合はマイナスにする
          if(tmp_s<0){
            d *= -1; 
          }
          //検出器より後ろの場合は再構成領域に加算しない

          if((z-32)*(z-32)+(t-32)*(t-32)+(s-46)*(s-46)<=100){//380
            //球領域にのみ逆投影100
            //dbetaをかける
            //Dはオブジェクト中心 (pow(22.5,2)/pow(22.5 - d,2)) *
            double tmp_output = 0;

            //2次補完をする
          if(z==32){
            //cout<<x<<" "<<xi<<" "<<y<<" "<<yi<<endl;
          }
            if(0<x && x<65 && 0<y && y<65){// && int(geometry[index])!=0){
              double geo_bl = (yi + 1- y) * 
                              ((xi + 1 - x) * map_out[xi*H + yi] + (x - xi)*map_out[(xi+1)*H+yi])
                              + (y - yi) * 
                              ((xi+1-x)*map_out[xi*H + yi+1] + (x - xi)*map_out[(xi+1)*H + yi+1]);
              tmp_output += geo_bl;
            }            
            //double weight = 1-abs(int(index_y+0.5)-index_y);
            //double pora = projection_out[H* int(beta) + int(index_y+0.5)] * weight + projection_out[H* int(beta) + int(index_y+0.5)-1] * (1-weight);
            double output = (pow(22.5,2)/pow(22.5 - d,2)) * tmp_output * beta_span * 2 * M_PI/360;
            image_xy[H*W*z + W*t + s] += output;//(t-1)->t
          }
       }
     }
   }
  }
#endif

  for(int a = 0; a<H; a++){
    for(int b = 0; b < W; b++){
      for(int c = 0; c < num_images; c++){
        image_xy[a*H*W + b * H + c] = ;
      }
    }
  }

  sprintf(wname_xy,"%s.raw",writeFileName_xy);
  writeRawFile(wname_xy, sizeof(float), H * W * num_images, image_xy);
  FILE* fp_xy = fopen(wname_xy, "wb");
  fwrite(image_xy , sizeof(float), H * W * num_images, fp_xy);

  /*sprintf(wname_xz,"%s.raw",writeFileName_xz);
  writeRawFile(wname_xz, sizeof(float), H * W * num_images, image_xz);
  FILE* fp_xz = fopen(wname_xz, "wb");
  fwrite(image_xz , sizeof(float), H * W * num_images, fp_xz);*/

  

#if 0
  sprintf(wname_map,"%s.raw",writeFileName_map);

  writeRawFile(wname_map, sizeof(float), H * W, map_out);

  FILE* fp = fopen(wname_map, "wb");
  fwrite(map_out , sizeof(float), H * W, fp);


#endif
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