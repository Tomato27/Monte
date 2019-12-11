#include<bits/stdc++.h>
#include <windows.h>
#include "Mersenne_twister.h"
#include "Mersenne_twister.cpp"
#pragma comment(lib, "winmm.lib");
#define num 512
#define num_images 145

using namespace std;

class csv{
public:
  string fname;
  bool csv_get=false;
  csv(string filename, int sizex, int sizey);
};
csv::csv(string filename, int sizex, int sizey){
  fname=filename;
  int sx=sizex;
  int sy=sizey;
}

void readcsv(csv csvfele, int sizex, int sizey, double** csvarray);

int main(void){
  cout<<"start";

  int i, count=0, sizex, sizey;
  int start_keV=140, end_keV=140;
  //Energy Spectrumはデータが取れていないので一様とする
  float Energy = start_keV; 

  double dens_H2O = 1.0; //水の密度i

  float min_e;
	//時定数によりSEED値を初期化
  init_genrand((unsigned long)time(NULL));

	//線減衰係数の値をmuに代入
  sizex=4; //線減衰係数配列のインデックス
  sizey=200;
  cout<<"a";

  double** csv_H2O = new double*[sizex];

  for(int k=0;k<sizex;k++){
    csv_H2O[k]=new double[sizey+5];
  }

  //csv読み込み
  csv c_H2O("xcom2.csv", 4, 200);
  readcsv(c_H2O, 4, 200, csv_H2O);
  
  
  float mu_H2O = csv_H2O[3][(int)(Energy+0.5)]*dens_H2O;

  float *f, *g;
  //int	f[num*num*num_images/49];	// 原画像用配列 
  f = (float*)malloc(num*num*num_images);

	FILE	*fp1, *fp2;	//ファイルポインタ
	fp1 = fopen( "tmp.raw" , "wb");
  cout<<"jj";
	fread( f, sizeof(float), num*num*num_images, fp1 );
	fclose(fp1);
  printf("\n%d bit exe",sizeof(void*)*8);


  g = (float*)malloc((420-90+1)*(336-115+1)*num_images);
  int index_crop = 0;
  for(int z = 0;z < num_images ;z++){
    for(int y = 105;y < 337 ;y++){
      for(int x = 90;x < 421 ;x++){
        //if(index_crop>2000000&&index_crop%100==0){cout<<index_crop<<" "<<f[x + y*(num) +z*(num)*(num)]<<" ";}
        g[index_crop] = f[x + y*(num) +z*(num)*(num)];//zがまずい
        index_crop++;
      }
    }
  }
  fp2 = fopen( "crop_3d_.raw" , "wb" );
  fwrite( g , sizeof(float) ,  (420-90+1)*(336-105+1)*num_images, fp2 );
  fclose( fp2 );

  free(f);
  free(g);
}


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