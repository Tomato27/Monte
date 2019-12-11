/*
 *
 * generate_uni_rand_wMersennetwister.c -- メルセンヌツイスタを用いて一様乱数を生成
 *
 */

#include<stdio.h>
#include<time.h>
#include<math.h>
#include "Mersenne_twister.h"
#include "Mersenne_twister.cpp"
#define NUM 1E5

int main(void){

	unsigned long seed;
	double r;
	FILE *fp;
	char filename[64];

	//時定数で乱数生成器のseed値を初期化
	seed = (unsigned long)time(NULL);
	init_genrand(seed);

	//生成した乱数のヒストグラム用のcsvファイル
	snprintf(filename, sizeof(filename),  "uniform_rand_hist.csv");
	fp = fopen(filename, "w");

	//区間(0,1)の乱数をNUM個生成して書き出し
	for(int i = 0; i < NUM; i++){
		r = genrand_real3();				//ここに関数を追加
		fprintf(fp, "%lf\n", r);
	}

	fclose(fp);

	printf("%s is generated!\n", filename);

	return 0;
}
