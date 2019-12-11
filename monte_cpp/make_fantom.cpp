#include <bits/stdc++.h>
#define N 65
using namespace std;
int main()
{
  /** 処理画像用配列 **/
unsigned char	g[N*N]={0};
  FILE	*fp1;

		for(int j=0 ; j<N ; j++){
      for(int k=0 ; k<N ; k++){
			  if((j-32)*(j-32)+(k-46)*(k-46)<=100){
				  g[j*N+k]++;
          //cout<<"nhai"<<g[j*N+k];
        }
      }
    }
  fp1 = fopen( "ball_fan.raw" , "wb" ); //*出力
  fwrite( g , sizeof(unsigned char) ,  N*N, fp1 );
  fclose(fp1);
}