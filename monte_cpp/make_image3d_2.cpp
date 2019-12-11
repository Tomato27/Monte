#include <stdio.h>
#include <math.h>
#include <bits/stdc++.h>
#define L 37
#define M 37
#define N 65
using namespace std;
int main()
{
	//unsigned char	g[N][57][N]={0};	/** 処理画像用配列 **/
  unsigned char	g[L*M*N]={0};
	FILE	*fp1;

	for(int k=0 ; k<N ; k++){
		for(int j=0 ; j<M ; j++){
      for(int i=0 ; i<L ; i++){
			  //g[i][j][k]=0;
			  if((i-18)*(i-18)+(j-18)*(j-18)+(k-32)*(k-32)<=100){
				  g[k*L*M+j*M+i]++;
			  }
        //if(1){
          //g[i][j][k]++;
        //}
			  //if(g[i*N*57+j*57+k]>1)cout<<"index"<<"("<<i<<","<<j<<","<<k<<")"<<":"<<g[i*N*57+j*57+k]<<" ";
      }
    }
	}
  fp1 = fopen( "spher_e3.raw" , "wb" ); //*回転画像1
  fwrite( g , sizeof(unsigned char) ,  L*M*N, fp1 );
  fclose(fp1);
  //free(g);
}