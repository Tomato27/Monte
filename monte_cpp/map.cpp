#include <stdio.h>
#include <math.h>
#include <bits/stdc++.h>
#define N 65
using namespace std;
int main()
{
	//unsigned char	g[N][57][N]={0};	/** 処理画像用配列 **/
  float g[N*N];

  int*	f=(int*)calloc(N * N, sizeof(int));	/** 原画像用配列 **/
  FILE	*fpi, *fpo;	/** ファイルポインタ **/
	fpi = fopen( "CBCTtest3\\monte_r0.raw" , "rb");
	fread(f, sizeof(int), N*N, fpi);

	for(int i=0 ; i<N ; i++){
		for(int j=0 ; j<N ; j++){
      //cout<<f[65*i+j]<<" ";
			g[65*i+j] = -log(f[65*i+j])+log(2000);
      if((i==32&&j==31)||(i==32&&j==32)||(i==32&&j==33)){
        g[65*i+j]+=0.05;
      }
    }
	}
  fpo = fopen( "mapr.raw" , "wb" ); //*回転画像1
  fwrite( g , sizeof(float) ,  N*N, fpo );
  fclose(fpo);
  //free(g);
}