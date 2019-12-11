#include<bits/stdc++.h>
#include <windows.h>
#include "Mersenne_twister.h"
#include "Mersenne_twister.cpp"
#pragma comment(lib, "winmm.lib");


#define NUM 10000000   //発生光子数
#define ScatterNUM  2   //散乱回数


using namespace std;

const char writeFileName[] = "circle3\\monte_c_";

char wname[22];

class csv{
public:
  string fname;
  bool csv_get=false;
  csv(string filename, int sizex, int sizey);
  /*
  string getfname(){
    return fname;
  }
  bool getflag(){
    return csv_get;
  }
  */
};
csv::csv(string filename, int sizex, int sizey){
  fname=filename;
  int sx=sizex;
  int sy=sizey;
}

class photon{
public:
  double x=0,y=0,z=0;
  double x_p=0, y_p=0, z_p=0;
  double before_vec[3]={0,0,0};
  double Energy;
  double theta=0;
  double theta_r=0;
  double phi=0;
  double theta_sum=0;
  double length = 0;
};

void readcsv(csv csvfele, int sizex, int sizey, double** csvarray);
void writeRawFile(const char fname[], const size_t size, const size_t num, void* image);
void writecsv(int height, int width, int* image, int* es, int num);
void delta_sampling(photon* p, double mu_H2O, double mu_Ca, int* geometry, double sin_theta_a, double cos_theta_a, double sin_phi_a, double cos_phi_a);
void add_result(photon* p, double phi, int* image, int* es, int count, int ditector_index,int Energy, int a, int q);

int main(void){

  int H=65,W=65;
  
  int*	geometry=(int*)calloc(H * W, sizeof(int));	/** 原画像用配列 **/
  unsigned char* image_i = (unsigned char*)calloc(H * W, sizeof(unsigned char));
	FILE	*fpi;	/** ファイルポインタ **/
	fpi = fopen( "circle_.raw" , "rb");
	fread(geometry, sizeof(int), H*W, fpi);
	fclose(fpi);
  int count_o9=0;
  for(int check_geometry=0;check_geometry<H*W;check_geometry++){   
  if(geometry[check_geometry]!=0&&geometry[check_geometry]!=1) cout<<"geometry"<<check_geometry<<":"<<geometry[check_geometry]<<" ";
  }

  DWORD	start, end, start_read , end_read, start_calc, end_calc;



  /*** 時間計測開始 ***/
  //start = timeGetTime();//コンパイル時にオプション-lwinmmが必要

  int i, count=0, sizex, sizey;
  int start_keV=140, end_keV=140;
	double sum_length = 0.0;
  double dens_H2O = 1.0; //水の密度i
  double dens_Ca = 1.550;

  float min_e;
	/*時定数によりSEED値を初期化*/
  init_genrand((unsigned long)time(NULL));

	/*線減衰係数の値をmuに代入*/
  sizex=4; //線減衰係数配列のインデックス
  sizey=200;

  double* output_H2O = new double[sizex];
  double** csv_H2O = new double*[sizex];
  double* output_Ca = new double[sizex];
  double** csv_Ca = new double*[sizex];
  for(int k=0;k<sizex;k++){
    csv_H2O[k]=new double[sizey+5];
    csv_Ca[k]=new double[sizey+5];
  }

int es[start_keV+1]={0};
int theta_s[181]={0};


  /*csv読み込み*/
  csv c_H2O("xcom2.csv", 4, 200);
  readcsv(c_H2O, 4, 200, csv_H2O);

  csv c_Ca("Ca.csv", 4, 200);
  readcsv(c_Ca, 4, 200, csv_Ca);


  //int keV=140;
  min_e=start_keV;
  
  for(int q=ScatterNUM;q<=ScatterNUM;q++){

    for(int keV=start_keV;keV<=end_keV;keV++){     
      double Energy = start_keV;
      double mu_H2O = csv_H2O[3][(int)(Energy+0.5)]*dens_H2O;
      double mu_Ca = csv_Ca[3][(int)(Energy+0.5)]*dens_Ca;
      double mu_Max = (mu_H2O <= mu_Ca) ? mu_Ca : mu_H2O ;

      sum_length=0;
      cout<<"--------------------"<<keV<<"keV"<<"--------------------\n";
      cout<<"(mu_H2O, mu_Ca, mu_Max) = ("<<mu_H2O<<","<<mu_Ca<<","<<mu_Max<<")"<<endl;

      /*(0,1]の一様乱数を生成して光子NUM個の光路長の和を計算*/
      count=0;
      /*calc start*/
      start_calc=timeGetTime();


      double theta=0., phi=0.;
      double min_e=140.;
      int primary_num=0;

      int minus_num=0;
      int coh_num=0;

      int* image = (int*)calloc(W*180, sizeof(int));

      for(i = 0; i <= NUM; i++){
        double cos_theta_a_new=1., cos_phi_a_new=0., sin_theta_a_new=0., sin_phi_a_new=1.;
        double cos_theta_a=0., cos_phi_a, sin_theta_a=1., sin_phi_a;
        //double theta_a=0., phi_a=0.;
        double phi_a_result;

        photon *p;
        p = new photon;
        double theta_sum = 0.;
        p -> x = 0;//0.0002*genrand_real3() - 0.0001;


        Energy = start_keV;
        mu_H2O=csv_H2O[3][(int)(Energy+0.5)]*dens_H2O;
        mu_Ca=csv_Ca[3][(int)(Energy+0.5)]*dens_Ca;


        //double r= genrand_real3();
        //double length = -1*log(1-r)/mu;
        double phia = genrand_real3() * 2 * M_PI;
        //phia  =M_PI/6;
        int ditector_index = 180*phia/(2*M_PI);//検出器のインデックス
        sin_phi_a = sin(phia);
        cos_phi_a = cos(phia);
        //cout<<"before ds"<<endl; //print debug
        delta_sampling(p , mu_H2O, mu_Ca, geometry, sin_theta_a, cos_theta_a, sin_phi_a, cos_phi_a);
        //cout<<p->x*p->x+p->y*p->y<<endl;

        if(p->x*p->x+p->y*p->y>=10*10){//1回目で検出器到着---------------primary-----------------change10->100
          primary_num++;
          if(q==0){
            count++;
            //image[H*32+32]++;
            //cout<<(int)(1+360*(phia/(2*M_PI))/2)<<" ";

            image[65*ditector_index+32]++;//32bit signedで見ような

            //image[32]++;
            //(int)(361*(phia/(2*M_PI))/2)
            es[140]++;
          }

          //continue;//lengthの加算まで飛ばしていた
          //break;
        }
        else{//散乱
          //cout<<"Scatterring \n";

          //p->z=length;
          p->before_vec[2]=p->z/sqrt(pow(p->x,2)+pow(p->y,2)+pow(p->z,2));//z方向に移動

          bool coh_flag=false;
          bool com_flag=false;

          for(int a=0; a<q; a++){//相互作用開始
            //mu = csvarray[3][(int)(Energy+0.5)]*dens_H2O;
            p->length = 0;
            double dens, ab, coh, com, mu;

            //H2O内かCa内かの判別
            //if(p->x<0){
              dens = dens_H2O;
              ab = csv_H2O[2][(int)(Energy+0.5)];
              coh = csv_H2O[0][(int)(Energy+0.5)];
              com = csv_H2O[1][(int)(Energy+0.5)];
              mu =  csv_H2O[3][(int)(Energy+0.5)];
            /*}
            else{
              dens = dens_Ca;
              ab = csv_Ca[2][(int)(Energy+0.5)];
              coh = csv_Ca[0][(int)(Energy+0.5)];
              com = csv_Ca[1][(int)(Energy+0.5)];
              mu =  csv_Ca[3][(int)(Energy+0.5)];;
            }*/

            double sc_rand = genrand_real3();

            if(sc_rand <=  ab / mu){//光電効果
              //cout<<"光電 \n";
              break;
            }
          else if(ab / mu < sc_rand  && sc_rand <= (ab + coh)/mu){

            //コヒーレント散乱----------------------------------------------------------------------------
              //cout<<"coh \n";
              coh_num++;
              //cout<<coh;

              float ditector_x=0,ditector_y=0;          //------------------------------------------------comと2重定義   
              double x_rotate, y_rotate, x_p_rotate, y_p_rotate;       
              double cos_check;      
              double p_r;

              if(com_flag==false){
                delta_sampling(p, mu_H2O, mu_Ca, geometry, sin_theta_a, cos_theta_a, sin_phi_a, cos_phi_a);
                p_r = phia;
                coh_flag = true;
              }
              else{
                delta_sampling(p, mu_H2O, mu_Ca, geometry, sin_theta_a_new, cos_theta_a_new, sin_phi_a_new, cos_phi_a_new);
                p_r = phi_a_result;           
              }

            add_result(p, p_r, image, es, count, ditector_index,Energy, a, q);
          }

            else{//コンプトン散乱
              //cout<<"com \n";
              

              //1:散乱角，エネルギー計算
              double lambda = 511.0 / Energy;
              //double lambda = 9.109383*pow(2.99792,2)*1E-15/(1.602*1E-16*Energy);

              double lambda_d=0.;
              bool track_flag=true;

              while(track_flag){
                double r1=genrand_real3();
                //r1=0.1;
                if(r1 < (lambda+2.0)/(9.0*lambda+2.0)){//track1 <=or<?
                  double r2 = genrand_real3();
                  //r2=0.2;
                  double ro = 1.0+(2.0/lambda)*r2;
                  double r3 = genrand_real3();
                  //r3=0.3;

                  if(r3<=4.0*((1./ro)-(1./(ro*ro)))){
                    lambda_d = ro*lambda;
                    track_flag=false;
                  }
                }
                else{//track2
                  double r2 = genrand_real3();
                  double ro = (lambda+2.)/(lambda+2.*(1.-r2));
                  double r3 = genrand_real3();
                  if(r3<=0.5*(pow((lambda-ro*lambda+1.),2)+(1./ro))){
                    lambda_d = ro*lambda;
                    track_flag=false;
                  }
                }
              }

              double theta = acos(1.-(lambda_d - lambda));///----------------------何かまずいかも
              double cos_theta = (1-(lambda_d - lambda));
              double sin_theta = sqrt(1.-pow((cos_theta),2));

              Energy = 511. / lambda_d;
              mu_H2O = csv_H2O[3][(int)(Energy+0.5)]*dens_H2O;
              mu_Ca =  csv_Ca[3][(int)(Energy+0.5)]*dens_Ca;

              //2：方位角，光路長計算
              double phi = genrand_real3()*2.*M_PI;

              //3:相対座標→絶対座標

              if(com_flag == true){
                sin_theta_a=sin_theta_a_new;
                cos_theta_a=cos_theta_a_new;
                sin_phi_a=sin_phi_a_new;
                cos_phi_a=cos_phi_a_new;
              }

              com_flag=true;

              cos_theta_a_new = cos_theta_a*cos_theta-sin_theta_a*sin_theta*cos(phi);
              sin_theta_a_new = sqrt(1.-pow(cos_theta_a_new,2));//絶対座標系の新しい角度


              cos_phi_a_new = (cos_theta_a*cos_phi_a*sin_theta*cos(phi) + sin_theta_a*cos_phi_a*cos_theta - sin_phi_a*sin_theta*sin(phi))/sin_theta_a_new;
              sin_phi_a_new = (cos_theta_a*sin_phi_a*sin_theta*cos(phi) + sin_theta_a*sin_phi_a*cos_theta + cos_phi_a*sin_theta*sin(phi))/sin_theta_a_new;


              
              /*if(ditector_index<0){
              ditector_index = -ditector_index;
              ditector_index = 90+(90-ditector_index);
              }*/
              
              //cout<<ditector_index<<endl;

              delta_sampling(p, mu_H2O, mu_Ca,geometry, sin_theta_a_new, cos_theta_a_new, sin_phi_a_new, cos_phi_a_new);

              float ditector_x=0, ditector_y=0;

              //if(p->x*p->x + p->y*p->y >= 100){//&& p->theta_r <= 90){
                //散乱で検出器到着(z座標，角度，散乱回数check)


              if(asin(sin_phi_a_new)>0){
                ditector_index = 180*asin(sin_phi_a_new)/(2*M_PI)+0.5;
              }
              else{
                ditector_index = 180*asin(sin_phi_a_new)/(2*M_PI)-0.5;
              }
              phi_a_result = asin(sin_phi_a_new);


              if(cos_phi_a_new<0){//asinの定義域は-pi/2 ~ pi/2なのでそれ以外の範囲の時は変換が必要
              //cout<<" "<<"transrated ";

                if(ditector_index>0){
                  ditector_index = 90-ditector_index;

                  phi_a_result  = M_PI - phi_a_result;
                  }
                else{
                  ditector_index = -ditector_index;
                  ditector_index = 90+ditector_index;
                  
                  phi_a_result = -M_PI - phi_a_result;
                }
              }

              if(ditector_index<0){
                ditector_index = 180+ditector_index;
                phi_a_result = 2*M_PI+phi_a_result;
              }
              //cout<<ditector_index<<" "<<phi_a_result<<endl;

              add_result(p, phi_a_result, image, es, count, ditector_index, Energy, a, q);

 /*                 if(Energy < min_e && a == q-1){
                    //最小エネルギーは検出されたものの内での最小エネルギーを指す
                    min_e = Energy;
                    cout<<"min_e:"<< min_e <<" sum_theta:"<<theta_sum*180/M_PI <<" theta "<< theta*180/M_PI <<endl;
                    //cout<<"thet_a:"<<theta_a*180/M_PI<<endl;
                    if(Energy>139.5){
                      //cout<<"com"<<endl;
                    }
                    //break;
                  }*/
                  //break;
                }
              //}
            }
          }
        }
        /*sum_length += p->length;
        p->length = 0;*/
      }

      end_calc=timeGetTime();
      cout<<"simulation_time = "<<end_calc-start_calc<<"[ms]"<<"\n\n";
      cout<<"-----result-----\n";
      //cout<<"primary_num:"<<primary_num<<endl;

      /*平均光路長の出力*/
      //cout<<"average length = " << sum_length/NUM << endl;
      cout<<"count = " << count <<" / " << NUM << "\n\n";
      //cout<<"min Energy("<<q<<"):"<<min_e<<endl;
      //cout<<"coh:"<<coh_num<<endl;

      //cout<<"writeFileName"+to_string(q)<<endl;
      //string a=to_string(q);

      sprintf(wname,"%s%d.raw",writeFileName, q);

      writeRawFile(wname, sizeof(int), 180 * 65, image);

      FILE* fp = fopen(wname, "wb");
      fwrite(image, sizeof(int), 180 * 65, fp);
      fclose(fp);
      writecsv(180, W, image, es,q);


   }

   
  }
    //return 0;
//}


void readcsv(csv csvfile, int sizex,int sizey ,double** csvarray){
  string fname=(string)csvfile.fname;
  ifstream ifs(fname);
  if(!ifs.is_open()){cout<<"ERROR : file didn't open"<<endl;}

  csvfile.csv_get=true;

  double com[sizey+5]={},coh[sizey+5]={}, photoel[sizey+5]={},mu[sizey+5]={};
  string tmp1[sizey+5];
  string tmp2, tmp3, tmp4;
  int i=0; //csv取得用配列のindex

  while(ifs.good()){
    i++; //indexは0start故keVより１先んじるので最初にindex++で帳尻合わせ
    std::getline(ifs, tmp1[i],',');
    std::getline(ifs, tmp2,',');
    std::getline(ifs, tmp3,',');
    std::getline(ifs, tmp4,'\n');
    com[i]=std::stod(tmp2);
    photoel[i]=std::stod(tmp3);
    mu[i]=std::stod(tmp4);

    csvarray[1][i]=com[i];
    csvarray[2][i]=photoel[i];
    csvarray[3][i]=mu[i];
  }
  ifs.close();

  string tmp;

  csvarray[0][1]=1.372;
  for(int a=2;a<sizey+1;a++){//謎------------------------------------------?
    tmp=tmp1[a];
    csvarray[0][a]=stod(tmp);
  }
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

void writecsv(int height, int width, int* image, int* es, int num){
  // ファイル入力ストリームの初期化
    ofstream ofs("circle3\\"+to_string(num)+"monte_c.csv");
    //cout<<"write csv in"<<endl;

    for(int i=0;i<width;i++){
      if(image[height*0+i] == 0){
        ofs << i << "," << 1 <<endl;
      }
      else{
        ofs << i << "," << image[height*0+i] <<endl;
      }
      cout<<image[30*65+i]<<" ";
    }
    ofs.close();
    //cout<<"write end"<<endl;

    ofstream ofs2("circle3\\"+to_string(num)+"es_c.csv");
    cout<<"es \n";
    for(int i=0;i<=140;i++){
      if(es[i] == 0){
        ofs2 << i << "," << 1 <<endl;
      }
      else{
        ofs2 << i <<","<< es[i] <<endl;
      }
    }
    ofs2.close();
    cout<<"es end";
    cout<<endl;
}

void delta_sampling(photon* p, double mu_H2O, double mu_Ca,int* geometry, double sin_theta_a, double cos_theta_a, double sin_phi_a, double cos_phi_a){
  double mu_max;
  //mu_max = max(mu_H2O, mu_Ca);
  mu_max=mu_H2O;//1媒質の場合
  bool loop_flag = true;
  double x1 = p->x,y1= p->y,z1 = p->z, length1 = p->length;
  double x = p->x,y= p->y,z = p->z, length = 0;
  p->x_p = x1, p->y_p = y1, p->z_p = z1;

  while(loop_flag){
    double beta = genrand_real3();
    double r = -log(beta)/mu_max;
    x +=  r * sin_theta_a * cos_phi_a;
    y +=  r * sin_theta_a * sin_phi_a;
    z +=  r * cos_theta_a;
    length += r;

    //cout<<mu_Ca/mu_max;
    float nu = genrand_real3();
    //cout<<"yo";
    //cout<<geometry[65*(int)(y*2+32.5)+(int)(x*2+32.5)]<<" ";
    if(geometry[65*(int)(y*2+32.5)+(int)(x*2+32.5)]==2){//Ca領域の時
      if(nu <=  mu_Ca/mu_max){
        loop_flag = false;
        //break;
      }
    }
    if(geometry[65*(int)(y*2+32.5)+(int)(x*2+32.5)]==0){//空気の時，直進し続ける
      if(sqrt(x*x+y*y)>10+0.5){//sin1度の分軽く足しておく 5*sqrt(5)+0.1
        loop_flag = false;
      }
    }
    
    else{
      if(nu <=  mu_H2O/mu_max){
        loop_flag = false;
        //break;
      }
    }
  }

  p->x = x;
  p->y = y;
  p->z = z;
  p->length = length;
  //cout<<p->z_p<<" "<<p->z<<endl;
  //cout<<length<<" ";
}

void add_result(photon* p, double phi, int* image, int* es, int count, int ditector_index,int Energy ,int a ,int q){
  double x_rotate = p->x*cos(-phi) - p->y*sin(-phi);
  double y_rotate = p->x*sin(-phi) + p->y*cos(-phi);            
  double x_p_rotate = p->x_p*cos(-phi) - p->y_p*sin(-phi);
  double y_p_rotate = p->x_p*sin(-phi) + p->y_p*cos(-phi);

  double ditector_y=(10-x_p_rotate)*(y_rotate-y_p_rotate)/(x_rotate-x_p_rotate)+y_p_rotate;

  if(x_rotate>=10){
    if(a != q-1){//特定の散乱回数の場合のみ検出
      return;
    }
    if(abs(ditector_y)>10){
      //cout<<ditector_index<<" "<<ditector_y<<endl;
    }       
    image[ditector_index*65+(int)(ditector_y*2+32.5)]++;
    count++;
    es[(int)Energy]++;
    return;
  }
}