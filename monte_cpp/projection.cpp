#include<bits/stdc++.h>

using namespace std;

const char writeFileName_project[] = "projection_test\\te2bl";

char wname_project[26];

void writeRawFile(const char fname[], const size_t size, const size_t num, void* image);

int main(void){

  int H=65,W=65;
  
  unsigned char*	geometry=(unsigned char*)calloc(H * W, sizeof(unsigned char));	/** 原画像用配列 **/
  float*	projection=(float*)calloc(H * W, sizeof(float));	/** 出力画像用配列 **/
	FILE	*fpi;	/** ファイルポインタ **/

	fpi = fopen( "ball_fan.raw" , "rb");
	fread(geometry, sizeof(unsigned char), H*W, fpi);
  for(int j=0;j<H;j++){
    for(int k=0;k<W;k++){
      //if(int(geometry[H*j+k])!=0)cout<<(int)geometry[H*j+k]<<" "<<j<<" "<<k<<endl;
    }
  }

  //投影
  int num_pr = 360;//投影数
  float ditector[65 * num_pr]={0};//検出器の輝度値が入る配列
  //double photon[2]= {};
  int itr_num[65] = {};//ditector到達のイテレーション回数
  float per_step[65] = {};

  //0度の時
  for(int t = 0; t < H; t++){//検出器の幅 
    //int t=31;

    //double start_x = -23;
    //double start_y = 0;
    //photon[0]=-23;
    //photon[1]=0;

    double photon_vec[2];//光子の座標が入るベクトル
    photon_vec[0] = 32.5;//x
    photon_vec[1] = -16 + t*0.5;//y

    int tmp = 0;
    int index = 0;

    //cout<<(int)geometry[31*H+50]<<endl;

    double vec_length = sqrt(pow(photon_vec[0],2)+pow(abs(photon_vec[1]),2));
    double theta0 = atan(abs(photon_vec[1])/photon_vec[0]);
    per_step[t] = 0.125*cos(theta0);
    //cout<<"theta:"<<theta<<endl;
    int count = 0;
    for(float i=17;i<32.5;i+=per_step[t]){//+ (photon_vec[1]*2*(i/32.5)))
      double x = 2*photon_vec[0]*i/32.5;
      double y = 2*(16.25+(photon_vec[1]*i/32.5));//16.25が良いが四捨五入されて下ピクセルに入ってしまうので

      index = (H*int(y)) + int(x);//int(H*y)とかやると死ぬ
      if(int(geometry[index])!=0){
        tmp += (int)geometry[index];
      }
      count++;
    }
    itr_num[t] = count;
    //cout<<t<<":"<<itr_num[t]<<" ";

    //ditector[t] = tmp * 0.25;
  }

  //0度以外の時
  for(int theta = 0; theta<num_pr; theta++){
    for(int t2 = 0; t2 < H; t2++){//検出器の幅 
      //int t=31;

      //0ではしっかり動くのに，sinのファクターがダメっぽい

      float start_x = -23;
      float start_y = 0;
      double primary_x = start_x*cos(M_PI*theta/180) - start_y*sin(M_PI*theta/180);
      double primary_y = start_x*sin(M_PI*theta/180) + start_y*cos(M_PI*theta/180);
      //cout<<"p:"<<primary_x<<" "<<primary_y<<" ";


      double photon_vec[2];//光子の座標が入るベクトル
      //photon_vec[0] = 32.5;//x
      //photon_vec[1] = -16 + t*0.5;//y
      
      photon_vec[0] = (32.5 - 23)*cos(M_PI*theta/180) - (- 16 + t2*0.5)*sin(M_PI*theta/180) - primary_x;//32.5*cos(theta*M_PI/180);//x
      //回転させた検出器座標-光源
      photon_vec[1] = (32.5 - 23)*sin(M_PI*theta/180) + (- 16 + t2*0.5)*cos(M_PI*theta/180) - primary_y;//y
      //cout<<"t2:"<<t2<<": "<<photon_vec[0]<<" "<<photon_vec[1]<<endl;

      double tmp = 0;
      int index = 0;

      //theta_違いそう,最悪これも配列にして渡せばよい
      //double theta_ = atan(abs(photon_vec[1])/photon_vec[0]);

      int count = 0;
      int num = 1;
      //cout<<itr_num[t2]<<" ";
      for(float i2=17; num<itr_num[t2]; i2+=per_step[t2]){//+ (photon_vec[1]*2*(i/32.5)))
        num++;
        //ここから先，書き換え必要?
        double x = (2 * (primary_x+23) + 2 * photon_vec[0]*i2/32.5);
        double y = 2*(-primary_y + 16.25 - (photon_vec[1]*i2/32.5));
        //if(t2==32)cout<<x<<" "<<y<<endl;
        int xi = int(x);
        int yi = int(y);

        index = (H*int(y)) + int(x);//int(H*y)とかやると死ぬ
        if(0<=x && x<65 && 0<=y && y<65){// && int(geometry[index])!=0){
          double geo_bl = (xi + 1- x) * 
                          ((yi + 1 - y) * geometry[yi*H + xi] + (y - yi)*geometry[(yi+1)*H+xi])
                          + (x - xi) * 
                          ((yi+1-y)*geometry[yi*H + xi+1] + (y - yi)*geometry[(yi+1)*H + xi+1]);
          tmp += geo_bl;
        }
      }
      ditector[theta * H + t2] = tmp * 0.25;
    }
  }

  

  sprintf(wname_project,"%s.raw",writeFileName_project);

  writeRawFile(wname_project, sizeof(float), H * num_pr, ditector);

  FILE* fp = fopen(wname_project, "wb");
  fwrite(ditector , sizeof(float), H * num_pr, fp);
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