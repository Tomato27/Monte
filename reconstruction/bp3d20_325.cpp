#include<bits/stdc++.h>

using namespace std;

const char writeFileName_map[] = "CBCTrecon\\mapca0_2t";
const char writeFileName_xy[] = "CBCTrecon\\xyca0_2t";
const char writeFileName_zy[] = "CBCTrecon\\zyca0_2t";

char wname_map[20];
char wname_xy[20];
char wname_zy[20];

void writeRawFile(const char fname[], const size_t size, const size_t num, void* image);

int main(void){

  int H=65*5,W=65*5,num_images=65;
  //int H_out =37*5, W_out = 37*5 ,num_out = 65*5;
  int H_out = 256, W_out = 256 ,num_out = 256;
  int num_proj  = 360;
  
  float*	map=(float*)calloc(H * W * num_proj, sizeof(float));	/** 原画像用配列 **/
  float*	map_w=(float*)calloc(H * W * num_proj, sizeof(float));	/** 重み付け画像用配列 **/
  float*	map_out=(float*)calloc(H * W * num_proj, sizeof(float));	/** 出力画像用配列 **/
	FILE	*fpi;	/** ファイルポインタ **/

  //3あり画像か無しかは要確認
	//fpi = fopen( "mapgpu0_20kyu.raw" , "rb");
  fpi = fopen( "mapg325_0_2tca.raw" , "rb");
	fread(map, sizeof(float), H*W * num_proj, fpi);


  float*	image_xy=(float*)calloc(H_out * W_out * num_out, sizeof(float));	/* 出力画像用配列 */
  float*	image_zy=(float*)calloc(H_out * W_out * num_out, sizeof(float));	/* 出力画像用配列 */
  //34,29,202
  //map[202*H*W + H*29 + 34] = -log(1) + log(20000);

  //step1 重み付け
  for(int num = 0; num < num_proj; num++){  
    for(int zeta = 0; zeta < H; zeta++){
      for(int p = 0; p < W; p++){
        //dlで割る,map作成時点でこれはやっておきたい
        //Dso->Dsdに変更
        //16.25->16.
        map_w[H*W*num+ H*zeta + p] = map[H*W*num + H*zeta + p]*(60./sqrt(pow(60.,2)+pow(-1*(zeta*0.1)+16.25,2)+pow(((p*0.1)-16.25),2)));
        //if(pow((zeta-32),2)+pow((p-32),2)<225){map_w[zeta*H + p]+=(22.5/sqrt(pow(22.5,2)+pow(-1*(zeta*0.5)+16.25,2)+pow((p*0.5)-16.25,2)))*1./0.5;}
        //if(pow((zeta-32),2)+pow((p-32),2)<9){map_w[zeta*H + p]=(22.5/sqrt(pow(22.5,2)+pow(-1*(zeta*0.5)+16.25,2)+pow((p*0.5)-16.25,2)))*2./0.5;}
      }
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

  float tmp1 = 0;

  //重畳積分
  for(int num = 0; num<num_proj;num++){
    for(int d = 0; d<H; d++){
      for(int b = 0; b<W; b++){
        float tmp = 0.;
        for(int c = 0; c<W; c++){
          if(isnan(map_w[H-1-b+c])){
            cout<<(H-1-b+c)<<" ";
          }
          tmp += map_w[num*H*W + d*H + c]*0.5*ramp[H-1-b+c];
          //if(a == 0 && c == 32){cout<< c <<" "<<64-b+c<<endl;}
        }
        map_out[num*H*W + d*H + b] = tmp;
        //map[a*H + b]+=tmp;
      /*if(b == 32 && a==32){
        tmp1 += map_out[num*H*W + a*H + b]; 
      }*/
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
for(int beta = 0; beta < beta_max; beta+=beta_span){//角度
    //光源をxz平面に関して1度ごと360方向回転させる
    float start_x = -160;
    float start_y = 0;//?
    //float beta = 0;
    double primary_x = start_x*cos(M_PI*beta/180) - start_y*sin(M_PI*beta/180);
    double primary_y = start_x*sin(M_PI*beta/180) + start_y*cos(M_PI*beta/180);

    for(int z = 0; z < num_out; z++){//奥行き z 0~num_out約60分 20~num_out-20 
    //25~29,34 (22~233でok)
      for(int t = 0; t < H_out; t++){//高さ↑ y 20 to H-20
        for(int s = 125; s < 130; s++){//幅→ 右半分だけ再構成 x 20 to W-20

          double photon_vec[3];//光子の座標が入るベクトル
          photon_vec[0] = (-12.8 + s*0.1) - primary_x;//x
          photon_vec[1] = (12.8 - t*0.1) - primary_y;//y
          photon_vec[2] = (12.8 - z*0.1);//z,マイナスに変更

          double tmp_x = photon_vec[0];

          //逆回転させて0度の位置に持ってくる
          photon_vec[0] = photon_vec[0]*cos(-1*M_PI*beta/180) - photon_vec[1]*sin(-1*M_PI*beta/180);
          photon_vec[1] = tmp_x*sin(-1*M_PI*beta/180) + photon_vec[1]*cos(-1*M_PI*beta/180);
        
          //x座標がditectorの位置になるようにベクトルを拡大し光子を動かす
          double to_ditector = 220 / photon_vec[0];
          for(int pv_index = 0; pv_index < 3; pv_index++){
            photon_vec[pv_index] *= to_ditector;
          }

          //y,zが検出器外ならこの先処理しない
          //z=32の時の球領域より外側のマイナス部分の値が入っていない．
          //ここが怪しそう
          //16.25, 23.25を小数点無しにずらした影響もありそう
          if(abs(photon_vec[1]) > 16.25 || abs(photon_vec[2]) > 16.25){
            //cout<<"al";
            continue;
          }
          

          //検出器座標
          //16.25->16
          double y = -1 * ((photon_vec[1]) * 10. - 325./2.);
          double x = -1 * (photon_vec[2] * 10. - 325./2.);
          int yi = int(y);
          int xi = int(x);


          //////////////////////ここから下をcheck///////////////////////////////////////

          //再構成空間における各点をbetaに従い回転させる
          //角度マイナスか?
          //16.25->16
          double tmp_s = (-12.8 + s*0.1)*cos(M_PI*beta/180) - (12.8 - t*0.1)*sin(M_PI*beta/180);
          double tmp_t = (-12.8 + s*0.1)*sin(M_PI*beta/180) + (12.8 - t*0.1)*cos(M_PI*beta/180);
          //s軸を表す方程式は y = tan(beta) * x なので点と直線の距離公式から
          double d = abs(-1*tan(beta)*tmp_t + 1*tmp_s)/(sqrt(1 + pow(tan(beta),2)));

          //s軸よりも右側の場合はマイナスにする
          if(tmp_s<0){
            d *= -1; 
          }
          //検出器より後ろの場合は再構成領域に加算しない

          //if((z-128)*(z-128)+(t-128)*(t-128)+(s-128)*(s-128)<=118*118){//pixel単位
            //cout<<" ;ipahlo";
            //球領域にのみ逆投影100
            //dbetaをかける
            //Dはオブジェクト中心 (pow(22.5,2)/pow(22.5 - d,2)) *
            double tmp_output = 0;

            //2次補完をする
            if(0<=x && x<=65*5 && 0<=y && y<=65*5){// && int(geometry[index])!=0){
              double geo_bl = (yi + 1- y) * 
                              ((xi + 1 - x) * map_out[beta*H*W + xi*H + yi] + (x - xi)*map_out[beta*H*W + (xi+1)*H+yi])
                              + (y - yi) * 
                              ((xi+1-x)*map_out[beta*H*W + xi*H + yi+1] + (x - xi)*map_out[beta*H*W + (xi+1)*H + yi+1]);
              tmp_output = geo_bl;
            }  
            double output = (pow(60,2)/pow(60. - d,2)) * tmp_output * beta_span * 2 * M_PI/360;
            image_xy[H_out*W_out*z + W_out*t + s] += output*2.7*5;//(t-1)->t
            image_zy[H_out*W_out*s + W_out*t + z] += output*2.7*5;
          //}
        }
      }
    }
}
#endif

  
  sprintf(wname_xy,"%s.raw",writeFileName_xy);
  writeRawFile(wname_xy, sizeof(float), H_out * W_out * num_out, image_xy);
  FILE* fp_xy = fopen(wname_xy, "wb");
  fwrite(image_xy , sizeof(float), H_out * W_out * num_out, fp_xy);

  
  sprintf(wname_zy,"%s.raw",writeFileName_zy);
  writeRawFile(wname_zy, sizeof(float), H_out * W_out * num_out, image_zy);
  FILE* fp_zy = fopen(wname_zy, "wb");
  fwrite(image_zy , sizeof(float), H_out * W_out * num_out, fp_zy);
  

  

#if 1
  sprintf(wname_map,"%s.raw",writeFileName_map);

  writeRawFile(wname_map, sizeof(float), H * W *num_proj, map_out);

  FILE* fp = fopen(wname_map, "wb");
  fwrite(map_out , sizeof(float), H * W *num_proj, fp);


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