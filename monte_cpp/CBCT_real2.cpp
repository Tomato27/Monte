#include<bits/stdc++.h>
#include <windows.h>
#include "Mersenne_twister.h"
#include "Mersenne_twister.cpp"
#pragma comment(lib, "winmm.lib");


#define NUM 10000000   //発生光子数
#define ScatterNUM  5  //散乱回数


using namespace std;

const char writeFileName[] = "CBCTtest3\\monte00";

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
void delta_sampling(photon* p, double mu_H2O, double mu_Ca, unsigned char* geometry, double sin_theta_a, double cos_theta_a, double sin_phi_a, double cos_phi_a);
void add_result(photon* p, double phi, int* image, int* es, int count, int detector_index,int Energy, int a, int q);

int main(void){
  //cout<<"joxl";

  int H=65,W=65,num_images=65;
  int num_nd = 0;
  int num_scatter =0;

   int err = 0;

   double r_xp;
  
  unsigned char*	geometry=(unsigned char*)calloc(37*5 * 37*5 *num_images*5, sizeof(unsigned char));	/** 原画像用配列 **/
  //unsigned char* image_i = (unsigned char*)calloc(H * W* num_images, sizeof(unsigned char));
	FILE	*fpi;	/** ファイルポインタ **/
	fpi = fopen( "spher01.raw" , "rb");
	fread(geometry, sizeof(unsigned char), 37*5*37*5*num_images*5, fpi);//phantomは0.1cm間隔
  
	/*for(int a=0 ; a<H ; a++){
		for(int b=0 ; b<H ; b++){
      for(int c=0 ; c<H ; c++){  
        if(geometry[a*H*H+b*H+c]==1){cout<<a<<","<<b<<","<<c<<","<<(int)geometry[a*H*H+b*H+c]<<","<<endl;  
        }
      }
    }
  }*/
	fclose(fpi);
  int count_o9=0;

  int* image = (int*)calloc(H*W*1, sizeof(int));//出力画像

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

  double** csv_al2 = new double*[sizex];
  double** csv_al10 = new double*[sizex];
  for(int k=0;k<sizex;k++){
    csv_al2[k]=new double[250+5];
    csv_al10[k]=new double[250+5];
  }

int es[start_keV+1]={0};
int theta_s[181]={0};


  /*csv読み込み*/
  csv c_H2O("xcom2.csv", 4, 200);
  readcsv(c_H2O, 4, 200, csv_H2O);

  csv c_Ca("Ca.csv", 4, 200);
  readcsv(c_Ca, 4, 200, csv_Ca);

  csv c_al2("125kv_al2mm.csv", 4, 200);
  readcsv(c_al2, 4, 250, csv_al2);

  csv c_al10("125kv_al10mm.csv", 4, 200);
  readcsv(c_al10, 4, 250, csv_al10);


  //int keV=140;
  min_e=start_keV;

  int per = 2000*5000;
  //int nlo=
  for(int q=ScatterNUM;q<=ScatterNUM;q++){
    //for(int keV=start_keV;keV<=end_keV;keV++){

      sum_length=0;

      /*(0,1]の一様乱数を生成して光子NUM個の光路長の和を計算*/
      count=0;
      /*calc start*/
      //start_calc=timeGetTime();

      /*for(int t = 0; t<H*W; t++){
        image[t] = 0;
      }*/

      double theta=0., phi=0.;
      double min_e=140.;
      int primary_num=0;

      int minus_num=0;
      int coh_num=0;

      int projection_step = 90;
      for(int num_p = 0; num_p < 360 ; num_p += projection_step*10){//投影数(xy平面で回転)
        for(i = 32; i < H-32; i++){//theta横
          for(int j = 32;j < W-32; j++){//phi縦
            for(int k  = 0;k < per;k++){
      
      int s_flag = -1;//エネルギーの探査方向
      
      double random_num = genrand_real3();
      
      if(random_num >= csv_al2[3][125]){
        s_flag = 1;
      }
      int Energy;//start_keV;
      
      for(int en = 125; en<=1||en>=250; en+=s_flag){
        if(csv_al2[3][en-1] <= random_num && random_num < csv_al2[3][en]){
          Energy = en*0.5;
        }
      }
      Energy = 140;

      double mu_H2O = csv_H2O[3][(int)(Energy+0.5)]*dens_H2O;
      double mu_Ca = csv_Ca[3][(int)(Energy+0.5)]*dens_Ca;
      double mu_Max = (mu_H2O <= mu_Ca) ? mu_Ca : mu_H2O ;
      double cos_theta_a_new=1., cos_phi_a_new=0., sin_theta_a_new=0., sin_phi_a_new=1.;

        //double cos_theta_a=0., cos_phi_a, sin_theta_a=1., sin_phi_a;
        if(i==0&&j==0&&k==0){
          cout<<"--------------------"<<"125keV2mm"<<"--------------------\n";
          cout<<"(mu_H2O, mu_Ca, mu_Max) = ("<<mu_H2O<<","<<mu_Ca<<","<<mu_Max<<")"<<endl;
          cout<<"-----"<<num_p<<"[°]-----"<<endl<<endl;
        }

        photon *p;
        p = new photon;
        double theta_sum = 0.;
        p -> x = -160;//-16.25*sqrt(3)/2;
        p -> y = 0;
        p -> z = 0;

        //zlをマイナスに
        double  xl = 60, yl = 16 - 0.5 * i, zl = 16 - 0.5 * j;

        double cos_phi_a, sin_phi_a;
        double theta_a = 0.5*M_PI - atan(zl/220);

        double cos_theta_a = cos(theta_a), sin_theta_a = sin(theta_a);

        double phi_a_result;

        Energy = start_keV;
        mu_H2O=csv_H2O[3][(int)(Energy+0.5)]*dens_H2O;
        mu_Ca=csv_Ca[3][(int)(Energy+0.5)]*dens_Ca;


        //double phia = (29.95-(59.9*j/64.))*M_PI/180.;///////////////////////////////////////////////////////////////////////////???
        //pixel中心から
        double phia = atan(yl/220);
        phia += M_PI*num_p/180.;
        //cout<<phia<<" ";


        int detector_index = 180*phia/(2*M_PI);//検出器の検出画素のインデックス
        sin_phi_a = sin(phia);
        cos_phi_a = cos(phia);

        //double primary_x = p->x;
        //double primary_y = p->y;

        //線源の回転 
        double primary_x = p->x*cos(M_PI*num_p/180) - p->y*sin(M_PI*num_p/180);
        double primary_y = p->x*sin(M_PI*num_p/180) + p->y*cos(M_PI*num_p/180);

        p->x = primary_x;
        p->y = primary_y;

        //x=-6の位置まで飛ばす
        double photon_vec[3];//光子の座標が入るベクトル
        photon_vec[0] = 60*cos(M_PI*num_p/180) - yl*sin(M_PI*num_p/180) - primary_x;//60 - primary_x;//x
        photon_vec[1] = 60*sin(M_PI*num_p/180) + yl*cos(M_PI*num_p/180) - primary_y;//yl - primary_y;//y
        photon_vec[2] = zl;//(-16 + z*0.5);//z

        //cout<<photon_vec[0]<<" "<<photon_vec[1]<<endl;

        //ここ書き換え必要
        double to = double(160-6)/220.;//double(210./220.);//
        //cout<<to<<endl;]

        double to_detector_x = p->x + photon_vec[0];
        double to_detector_y = p->y + photon_vec[1];

        for(int to_n = 0; to_n<3;to_n++){
          photon_vec[to_n]*= to;
        }

        p->x += photon_vec[0];
        p->y += photon_vec[1];
        p->z += photon_vec[2];


        //cout<< p->x <<" "<< p->y <<" "<< p->z <<endl;
        if(i==0&&j==0&&k<2){
          cout<< primary_x <<" "<< primary_y <<endl;
          cout<< photon_vec[0] <<" "<< photon_vec[1] <<" "<< photon_vec[2] <<endl;
          cout<< p->x <<" "<< p->y <<" "<< p->z <<endl<<endl;
        }

        delta_sampling(p , mu_H2O, mu_Ca, geometry, sin_theta_a, cos_theta_a, sin_phi_a, cos_phi_a);
        //p->x=50.1;

        double tmp_x = p->x;
        double tmp_y = p->y;
        double tmp_z = p->z;

        double tmp_x2 =  tmp_x*cos(M_PI*-num_p/180) - tmp_y*sin(M_PI*-num_p/180);    

        //if(q!=0)cout<<k<<" ";    

        if(tmp_x2>60){//1回目で検出器到着---------------primary-----------------&&yの範囲
        //abs(p->x) >= abs(to_detector_x) && abs(p->y) >= abs(to_detector_y)
          primary_num++;
          if(q==0){
            count++;
            //image[H*32+32]++;

            //image[65*i+j]++;//32bit signedで見ような
            double p_ = 60;
            //double d_z = (p->z/(p->x+p_))*p_*2;//*cos(phia);
            //double d_y = (p->y/(p->x+p_))*p_*2;//ok?
            
            //線源を元の座標に戻すために逆回転
            //double tmp_x = p->x;
            //double tmp_y = p->y;
            //double tmp_z = p->z;

            p->x =  tmp_x*cos(M_PI*-num_p/180) - tmp_y*sin(M_PI*-num_p/180);
            p->y =  tmp_x*sin(M_PI*-num_p/180) + tmp_y*cos(M_PI*-num_p/180);

            double d_z = ((p->z - 0)/(p->x-(-160)))*60 + 160*p->z/(p->x+160);
            double d_y = ((p->y - 0)/(p->x-(-160)))*60 + 160*p->y/(p->x+160);

            
            int result_y = -1*(int(d_y*2-32.5));
            int result_x = -1*(int(d_z*2-32.5));
            if(result_y<0||result_y>64||result_x<0||result_x>64){
              err++;
            }

          /*if(i==0&&j==0&&k<100){
            cout<<k<<": ";
            cout<< d_z <<" "<< d_y <<"  ";//" "<< p->z <<"  ";
            cout<<result_y<<" "<<result_x<<endl<<endl;
          }*/
            
            image[result_y*65+result_x]++;
            //continue;

            es[140]++;
          }
         //continue;//lengthの加算まで飛ばしていた
          //break;
        }
        else{//散乱
          num_scatter++;
          bool coh_flag=false;
          bool com_flag=false;
          double phi_a_result=0;

          for(int a=0; a<q; a++){//相互作用開始
            //cout<<k<<" ";
            //mu = csvarray[3][(int)(Energy+0.5)]*dens_H2O;
            p->length = 0;
            double dens, ab, coh, com, mu;

            double x_rotate_c = p->x*cos(M_PI*-phi/180) - p->y*sin(M_PI*-phi/180);
            double y_rotate_c = p->x*cos(M_PI*-phi/180) - p->y*sin(M_PI*-phi/180);
            if(x_rotate_c>=60||abs(y_rotate_c)>=16.25||abs(p->z)>16.25){
              //cout<<k<<" "<<"break\n";
              break;//a+=10;
            }            

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
            //else 
            else if(ab / mu < sc_rand  && sc_rand <= (ab + coh)/mu){
            //コヒーレント散乱----------------------------------------------------------------------------
              //cout<<a<<"coh ";//  \n";
              coh_num++;
              num_scatter--;
              //cout<<coh;

              //float ditector_x=0,ditector_y=0;             
              double x_rotate, y_rotate, x_p_rotate, y_p_rotate;       
              double cos_check;      
              double p_r;

              if(com_flag==false){
                delta_sampling(p, mu_H2O, mu_Ca, geometry, sin_theta_a, cos_theta_a, sin_phi_a, cos_phi_a);
                p_r = phia;
                coh_flag = true;
              }
              else{//過去にコンプトンあるなら，その時の角度を参照
                delta_sampling(p, mu_H2O, mu_Ca, geometry, sin_theta_a_new, cos_theta_a_new, sin_phi_a_new, cos_phi_a_new);
                p_r = phi_a_result;   
              }
              //if(a==q-1){
                //add_result(p, num_p, image, es, count, detector_index,Energy, a, q);//p_r->num_pに変更
              //}
          }

            else{//コンプトン散乱
              //cout<<a<<"com ";//\n";
              com_flag =true;

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
              //lambda_d=lambda+0.1;

              double theta = acos(1.-(lambda_d - lambda));///----------------------何かまずいかも

              double cos_theta = (1.-(lambda_d - lambda));//cos(theta);//cos(0.5*M_PI + atan(50./220.));
              //0.5*_PIだとnanに
              double sin_theta = sqrt(1.-pow((cos_theta),2));

              //double cos_theta = cos_tmp;
              //double sin_theta = sqrt(1.-pow((cos_theta),2));

              Energy = 511. / lambda_d;
              mu_H2O = csv_H2O[3][(int)(Energy+0.5)]*dens_H2O;
              mu_Ca =  csv_Ca[3][(int)(Energy+0.5)]*dens_Ca;

              //2：方位角，光路長計算
              double phi = -genrand_real3()*2.*M_PI;
              //乗算では.を忘れるな
              //cout<<phi<<" ";

              //phi = M_PI/6;
              //r = genrand_real3();
              //length = -1.*log(1.-r)/mu;

              //3:相対座標→絶対座標
              sin_theta_a=sin_theta_a_new;
              cos_theta_a=cos_theta_a_new;
              sin_phi_a=sin_phi_a_new;
              cos_phi_a=cos_phi_a_new;

              cos_theta_a_new = cos_theta_a*cos_theta-sin_theta_a*sin_theta*cos(phi);//cos(0.5*M_PI) - 
              sin_theta_a_new = sqrt(1.-pow(cos_theta_a_new,2));//絶対座標系の新しい角度


              cos_phi_a_new = (cos_theta_a*cos_phi_a*sin_theta*cos(phi) + sin_theta_a*cos_phi_a*cos_theta - sin_phi_a*sin_theta*sin(phi))/sin_theta_a_new;
              sin_phi_a_new = (cos_theta_a*sin_phi_a*sin_theta*cos(phi) + sin_theta_a*sin_phi_a*cos_theta + cos_phi_a*sin_theta*sin(phi))/sin_theta_a_new;
            
              /*if(ditector_index<0){
              ditector_index = -ditector_index;
              ditector_index = 90+(90-ditector_index);
              }*/
              
              //cout<<ditector_index<<endl;

              delta_sampling(p, mu_H2O, mu_Ca, geometry, sin_theta_a_new, cos_theta_a_new, sin_phi_a_new, cos_phi_a_new);
              //if(p->x>60)cout<<p->x<<" "<<p->y<<" "<<p->z<<endl;

              p->theta=theta;
              theta_sum+=theta;

              double v_length = sqrt(pow(p->length*sin_theta_a_new*cos_phi_a_new,2) + pow(p->length*sin_theta_a_new*sin_phi_a_new,2) +pow(p->length*cos_theta_a_new,2));
              p->before_vec[0] = p->length*sin_theta_a_new*cos_phi_a_new / v_length;
              p->before_vec[1] = p->length*sin_theta_a_new*sin_phi_a_new / v_length;
              p->before_vec[2] = p->length*cos_theta_a_new / v_length;
              //vectorは長さ1に正規化

              p->theta_r=(int)(1+theta_sum*180/M_PI)%180;

              float detector_x=0, detector_y=0;

              //if(p->x*p->x + p->y*p->y >= 100){//&& p->theta_r <= 90){
                //散乱で検出器到着(z座標，角度，散乱回数check)


                if(asin(sin_phi_a_new)>0){
                  detector_index = 180*asin(sin_phi_a_new)/(2*M_PI)+0.5;
                }
                else{
                  detector_index = 180*asin(sin_phi_a_new)/(2*M_PI)-0.5;
                }
                phi_a_result = asin(sin_phi_a_new);

                if(cos_phi_a_new<0){//asinの定義域は-pi/2 ~ pi/2なのでそれ以外の範囲の時は変換が必要
                
                if(detector_index>0){
                  detector_index = 90-detector_index;

                  phi_a_result  = M_PI - phi_a_result;
                  }
                else{
                  detector_index = -detector_index;
                  detector_index = 90+detector_index;
                  
                  phi_a_result = -M_PI - phi_a_result;
                }
              }

              if(detector_index<0){
                detector_index = 180+detector_index;
                phi_a_result = 2*M_PI+phi_a_result;
              }
              //cout<<ditector_index<<" "<<phi_a_result<<endl;

                double x_rotate=0, y_rotate=0, x_p_rotate=0, y_p_rotate=0;
                //座標の回転は今ついてる角度と逆方向に回さねば
                phi_a_result = num_p;//M_PI*210./180.;
                x_rotate = p->x*cos(-phi_a_result) - p->y*sin(-phi_a_result);
                y_rotate = p->x*sin(-phi_a_result) + p->y*cos(-phi_a_result);
                x_p_rotate = p->x_p*cos(-phi_a_result) - p->y_p*sin(-phi_a_result);
                y_p_rotate = p->x_p*sin(-phi_a_result) + p->y_p*cos(-phi_a_result);

                //cout<<p->x<<" "<<p->y<<" "<<p->z<<endl;

                //detector_y=(10-x_p_rotate)*(y_rotate-y_p_rotate)/(x_rotate-x_p_rotate)+y_p_rotate;//x==10の時のy座標

                double d_z=(60-p->x_p)*(p->z-p->z_p)/(p->x-p->x_p)+p->z_p;
                double d_y=(60-p->x_p)*(p->y-p->y_p)/(p->x-p->x_p)+p->y_p;


                if(x_rotate>=60&&abs(d_z)<16.25&&abs(d_y)<16.25){//&&a==q-1){// && abs(p->z)<16.25 && abs(y_rotate)<16.25){// && x_rotate>10
                //検出器を通過したかcheck,今回は必ず検出器まで到達する
                  /*if(a!=q-1){
                    break;
                  }*/
                  //image[detector_index*65+(int)(detector_y*2+32.5)]++;
                  
                  //double d_z = (p->z - p->z_p)/(x_rotate-x_p_rotate)*60 + (x_rotate*p->z_p - x_p_rotate*p->z)/(x_rotate - x_p_rotate);
                  //double d_y = (y_rotate - y_p_rotate)/(x_rotate-x_p_rotate)*60 + (x_rotate*y_p_rotate - x_p_rotate*y_rotate)/(x_rotate - x_p_rotate);

                  //if(abs(d_z)<16.25&&abs(d_y)<16.25){

                  int result_y = -1*(int(d_y*2-32.5));
                  int result_x = -1*(int(d_z*2-32.5));
                  //cout<<"s:"<<d_z<<" "<<d_y<<endl;


                  count++;
                  es[(int)(Energy+0.5)]++;

                  image[result_y*65+result_x]++;
                  /*}
                  else{
                    //cout<<d_z<<" ";
                  }*/
                  //cout<<"yo";

                  /*if(Energy < min_e && a == q-1){
                    //最小エネルギーは検出されたものの内での最小エネルギーを指す
                    min_e = Energy;
                    cout<<"min_e:"<< min_e <<" sum_theta:"<<theta_sum*180/M_PI <<" theta "<< theta*180/M_PI <<endl;
                    //cout<<"thet_a:"<<theta_a*180/M_PI<<endl;
                    break;
                  }*/
                  //break;
                }
                else if(p->x>60||abs(p->y)>16.25||abs(p->z)>16.25){
                  //cout<<p->x<<" "<<p->y<<" "<<p->z<<" "<<num <<endl;
                  num_nd++;
                }
              //}
            }
           //sum_length += p->length;  //primaryの光路長を平均光路長に含めぬ場合
           //p->length = 0;
          }

        }//散乱

        }//k
        }//j
      } //for i
      }//num_p
    //}//kev

      //end_calc=timeGetTime();
      //cout<<"simulation_time = "<<end_calc-start_calc<<"[ms]"<<"\n\n";
      cout<<"-----result-----\n";
      cout<<"num_scatter "<<num_scatter<<endl<<
      "num not ditected "<<num_nd<<endl;
      //cout<<"primary_num:"<<primary_num<<endl;

      /*平均光路長の出力*/
      //cout<<"average length = " << sum_length/NUM << endl;
      cout<<"count = " << count <<" / " << H*W*per << "\n\n";
      //cout<<"min Energy("<<q<<"):"<<min_e<<endl;
      //cout<<"coh:"<<coh_num<<endl;

      //cout<<"writeFileName"+to_string(q)<<endl;
      //string a=to_string(q);

      sprintf(wname,"%s%d.raw",writeFileName, q);

      writeRawFile(wname, sizeof(int), H * W *1, image);

      FILE* fp = fopen(wname, "wb");
      fwrite(image, sizeof(int), H * W *1, fp);
      fclose(fp);
      //writecsv(180, W, image, es,q);
      cout<<"err;"<<err;


   }//ScatterNum

   
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
    ofstream ofs("CBCTtest3\\"+to_string(num)+"monte_c.csv");
    //cout<<"write csv in"<<endl;
    int sum = 0;

    for(int i=0;i<width;i++){
      if(image[height*0+i] == 0){
        ofs << i << "," << 1 <<endl;
      }
      else{
        ofs << i << "," << image[height*0+i] <<endl;
      }
      cout<<image[32*65+i]<<" ";
      sum+=image[32*65+i];
    }
    ofs.close();
    //cout<<"write end"<<endl;

    ofstream ofs2("CBCTtest3\\"+to_string(num)+"es_c.csv");
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
    cout<<sum<<endl;
    cout<<endl;
}


void delta_sampling(photon* p, double mu_H2O, double mu_Ca,unsigned char* geometry, double sin_theta_a, double cos_theta_a, double sin_phi_a, double cos_phi_a){
  double mu_max;
  //mu_max = max(mu_H2O, mu_Ca);
  mu_max=mu_H2O;//1媒質の場合
  bool loop_flag = true;
  bool air_flag = true;
  double x2 = p->x,y2= p->y,z2 = p->z;
  double x1 = p->x,y1= p->y,z1 = p->z, length1 = p->length;
  double x = p->x,y= p->y,z = p->z, length = 0;
  p->x_p = x1, p->y_p = y1, p->z_p = z1;

  int l = 0;
  //int geo_index;
  int check;


  while(loop_flag){
    double beta = genrand_real3();
    double r = -log(beta)/mu_max;
    
    x +=  r * sin_theta_a * cos_phi_a;
    y +=  r * sin_theta_a * sin_phi_a;
    z +=  r * cos_theta_a;
    length += r;

    check=0;

    //cout<<mu_Ca/mu_max;
    float nu = genrand_real3();

    /*if(abs(x)>=61||abs(y)>=16.25||abs(z)>=16.25){
      check = 0;
    }
    else{*/
      //geo_index = (int)(rint(rint(x)*2+rint(32.5))+rint(rint(y)*2+rint(32.5))*65+rint(rint(z)*2+rint(32.5))*65*65);
      //check = (int)geometry[geo_index];
      //if(50<x&&x<60){//((x+160)-160)*((x+160)-160)+((y+16.25)-16.25)*((y+16.25)-16.25)+((z+16.25)-16.25)*((z+16.25)-16.25)<=25){
      if(-6<x&&x<10&&-6<y&&y<6&&-10<z&&z<10){
        if(geometry[int(185*185*(rint(z*10+160)) + 185*(rint(y*10)+90) + (rint(x*10)+90))]==1){ 
          check = 1;
          //air_flag=false;
        }
      }

      //check=0;//-----------------------------------------
    //}

    if(check==0){//空気の時
      //空気の時，直進し続ける
      if(abs(x)>=62||abs(y)>=62||abs(z)>=17){
        loop_flag = false;
      }
    }
    
    else if(check==1){//H2Oの時
      if(nu <=  mu_H2O/mu_max){
        loop_flag = false;
        //break;
      }
    }
#if 0
    else{//Ca領域の時
      if(air_flag){
        air_flag = false;
        continue;
      }
      if(nu <=  mu_Ca/mu_max){
        loop_flag = false;
        //break;
      }
    }
#endif
    l++;
    //cout<<x<<" "<<y<<" "<<z<<endl;
  }

  p->x = x;
  p->y = y;
  p->z = z;
  //cout<<p->x<<" "<<p->y<<" "<<p->z<<endl;

  p->length = length;
}

void add_result(photon* p, double phi, int* image, int* es, int count, int detector_index,int Energy ,int a ,int q){
  double x_rotate = p->x*cos(M_PI*-phi/180.) - p->y*sin(M_PI*-phi/180.);
  double y_rotate = p->x*sin(M_PI*-phi/180.) + p->y*cos(M_PI*-phi/180.);            
  double x_p_rotate = p->x_p*cos(M_PI*-phi/180.) - p->y_p*sin(M_PI*-phi/180.);
  double y_p_rotate = p->x_p*sin(M_PI*-phi/180.) + p->y_p*cos(M_PI*-phi/180.);
  //z, z_pはそのままでおｋ

  //double detector_y = (60-x_p_rotate)*(y_rotate-y_p_rotate)/(x_rotate-x_p_rotate)+y_p_rotate;
  //double detector_z = (60-x_p_rotate)*(p->z-p->z_p)/(x_rotate-x_p_rotate)+p->z_p;
  //double d_y = (60-x_p_rotate)*(y_rotate-y_p_rotate)/(x_rotate-x_p_rotate)+y_p_rotate;
  //double d_z = (60-x_p_rotate)*(p->z-p->z_p)/(x_rotate-x_p_rotate)+p->z_p;

  double d_z = (p->z - p->z_p)/(x_rotate-x_p_rotate)*60 + (x_rotate*p->z_p - x_p_rotate*p->z)/(x_rotate - x_p_rotate);
  double d_y = (y_rotate - y_p_rotate)/(x_rotate-x_p_rotate)*60 + (x_rotate*y_p_rotate - x_p_rotate*y_rotate)/(x_rotate - x_p_rotate);

  //double d_z = ((p->z - 0)/(p->x-(-160)))*60 + 160*p->z/(p->x+160);
  //double d_y = ((p->y - 0)/(p->x-(-160)))*60 + 160*p->y/(p->x+160);


  if(x_p_rotate<60&&x_rotate>=60&&abs(d_y)<16.25&&abs(d_z)<16.25){
    if(a != q-1){//特定の散乱回数の場合のみ検出
      return;
    }
  
    int result_y = -1*(int(d_y*2-32.5));
    int result_x = -1*(int(d_z*2-32.5));

  /*if(result_y<0||result_y>64||result_x<0||result_x>64){
    err++;
  }*/
    
    //image[result_y*65+result_x]++;

    //image[ditector_index*65+(int)(ditector_y*2+32.5)]++;
    a+=10;
    count++;
    //es[(int)Energy]++;
    //return;
  }
}