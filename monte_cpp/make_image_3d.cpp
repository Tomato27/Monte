#include <stdio.h>
#include <math.h>
#include <bits/stdc++.h>
#define N 65
using namespace std;
int main()
{
	//unsigned char	g[N][57][N]={0};	/** 処理画像用配列 **/
  unsigned char	g[N*N*N]={0};
  unsigned char go[36*36*N]={0};
	FILE	*fp1;

	for(int i=0 ; i<N ; i++){
		for(int j=0 ; j<N ; j++){
      for(int k=0 ; k<N ; k++){
			  //g[i][j][k]=0;
			  if((i-32)*(i-32)+(j-32)*(j-32)+(k-46)*(k-46)<=100){
				  g[i*N*N+j*N+k]++;
          if(14<j&&j<50&&28<k&&k<64){go[i*N*N + (j-14) * N + (k-28)]++;}
          //g[i][j][k]++;
				  if((i-32)*(i-32)+(j-32)*(j-32)+(k-46)*(k-46)<=4){
					  //g[i*N*N+j*N+k]++;
            //g[i][j][k]++;
				  }
			  }
        //if(1){
          //g[i][j][k]++;
        //}
			  //if(g[i*N*57+j*57+k]>1)cout<<"index"<<"("<<i<<","<<j<<","<<k<<")"<<":"<<g[i*N*57+j*57+k]<<" ";
      }
    }
	}
  //unsigned char go[36*36*N]={0};
  for(int a = 0; a<N; a++){
    for(int b = 0; b <N ; b++){//14~N^14
      for(int c = 0; c < N; c++){//N-36~N
      if(b>=14&&b<=50&&c>=28&&c<=64){
        //go[a*N*N + (b-14) * N + (c-28)] = g[a*N*N+N*b+c];
      }
        //go[a*N*N+(b-14)*N+(c-(N-36))]++;
    }
  }
}

  fp1 = fopen( "spher_e3.raw" , "wb" );
  fwrite( go , sizeof(unsigned char) ,  36*36*N, fp1 );
  fclose(fp1);
  //free(g);
}