#include<bits/stdc++.h>

using namespace std;
void writeRawFile(const char fname[], const size_t size, const size_t num, void* image);

int main(void){

  int H=65,W=65;
  
  unsigned char*	im=(unsigned char*)calloc(H*W, sizeof(unsigned char));	/** 投影画像用配列 **/
  float*	im_out=(float*)calloc(H * W, sizeof(float));	/** 出力画像用配列 **/
	FILE	*fpi;	/** ファイルポインタ **/

  //3あり画像か無しかは要確認
	fpi = fopen("fbp2d_non(pora).raw" , "rb");
	fread(im, sizeof(unsigned char), H*W, fpi);

  for(int a = 0;a<H;a++){
    for(int b=0;b<W;b++){
      im_out[a*H+b]=float(im[a*H+b])/205.; 
    }
  }


  FILE* fp_xy = fopen("fbp2d_non(pora)f.raw", "wb");
  fwrite(im_out , sizeof(float), H * W, fp_xy);
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