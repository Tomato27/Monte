#include <bits/stdc++.h>
#include <curand.h>
#include <curand_kernel.h>
#pragma comment(lib, "winmm.lib");

//とりあえず単一エネルギー

#define NUM 10000000 //発生光子数
#define ScatterNUM 5 //散乱回数

#define detector_x 65*5
#define detector_y 65*5

#define H 65*5
#define W 65*5
#define num_images 65
#define num_proj 360 

#define num_photon 10000

using namespace std; 

//energy125kv,al10mm

const char writeFileName0[] = "proj325_teth0cyu8e.raw";
const char writeFileName_m0[] = "map325_teth0cyu8e.raw";

const char writeFileName5[] = "proj325_teth5cyu8e.raw";
const char writeFileName_m5[] = "map325_teth5cyu8e.raw";

const char *H2O_c = "xcom2.txt";
const char *Ca_c = "Ca.txt";

const char geo_name[] = "spher01.raw";

char wname[11];

class csv
{
public:
  string fname;
  bool csv_get = false;
  csv(string filename, int sizex, int sizey);
};
csv::csv(string filename, int sizex, int sizey)
{
  fname = filename;
  //int sx=sizex;
  //int sy=sizey;
}

class photon
{
public:
  float x = 0, y = 0, z = 0;
  float x_p = 0, y_p = 0, z_p = 0;
  float before_vec0 = 0, before_vec1 = 0, before_vec2 = 0; //配列から変更
  float Energy;
  float theta = 0;
  float phi = 0;
  float length = 0;
  __device__ void delta_sampling(float mu_H2O, float mu_Ca, unsigned char *geometry, float sin_theta_a, float cos_theta_a, float sin_phi_a, float cos_phi_a, curandStateMRG32k3a *st);
};

//void readXcom(const char* csv_c, float* ab, float* coh, float* com, float* mu);
void readXcom();
void readxray();
void writeRawFile(const char fname[], const size_t size, const size_t num, void *image);
void writecsv(int height, int width, int *image, int *es, int num);

__global__ void projection(int per, float mu_H2O, float mu_Ca, unsigned char *geometry, int *image0, int *image5, curandStateMRG32k3a *state_gpu, float *ab_H2O, float *coh_H2O, float *com_H2O, float *mua_H2O,
                           float *ab_Ca, float *coh_Ca, float *com_Ca, float *mua_Ca, float *al10mm);
//後で散乱を足して，csv_H2Oを渡すように変える
__global__ void LaunchPhoton(curandStateMRG32k3a *state, int seed);
//__device__ void delta_sampling(photon p_gpu, float mu_H2O, float mu_Ca, unsigned char* geometry, float sin_theta_a, float cos_theta_a, float sin_phi_a, float cos_phi_a,curandStateMRG32k3a* state_gpu);
__global__ void RandStateGenerator(curandStateMRG32k3a *state_gpu);
void add_result(photon *p, float phi, int *image, int *es, int count, int ditector_index, int Energy, int a, int q);

float ab_H2O[201], coh_H2O[201], com_H2O[201], mua_H2O[201];
float ab_Ca[201], coh_Ca[201], com_Ca[201], mua_Ca[201];
float al2mm[251];
float al10mm[251];

int main(void)
{
  //cout<<"joxl";

  //int H=65,W=65,num_images=65 ,num_proj = 360;

  unsigned char *geometry = (unsigned char *)calloc(37 * 5 * 37 * 5 * num_images * 5, sizeof(unsigned char)); /** 原画像用配列 **/
  unsigned char *geometry_gpu;

  FILE *fpi; /** ファイルポインタ **/

  fpi = fopen(geo_name, "rb");
  fread(geometry, sizeof(unsigned char), 37 * 5 * 37 * 5 * num_images * 5, fpi); //phantomは0.1cm間隔

  cudaMalloc((void **)&geometry_gpu, sizeof(unsigned char) * 37 * 5 * 37 * 5 * num_images * 5);
  cudaMemcpy(geometry_gpu, geometry, sizeof(unsigned char) * 37 * 5 * 37 * 5 * num_images * 5, cudaMemcpyHostToDevice);

  fclose(fpi);
  //int count_o9=0;

  int *image0 = (int *)calloc(detector_x * detector_y * 1 * num_proj, sizeof(int));     //処理用画像65*65
  int *image_out0 = (int *)calloc(detector_x * detector_y * 1 * num_proj, sizeof(int)); //処理用画像65*65

  int *image_gpu0;
  cudaMalloc((void **)&image_gpu0, sizeof(int) * detector_y * detector_x * num_proj);
  cudaMemcpy(image_gpu0, image0, sizeof(int) * detector_y * detector_x * num_proj, cudaMemcpyHostToDevice);

  int *image5 = (int *)calloc(detector_x * detector_y * 1 * num_proj, sizeof(int));     //処理用画像65*65
  int *image_out5 = (int *)calloc(detector_x * detector_y * 1 * num_proj, sizeof(int)); //処理用画像65*65

  int *image_gpu5;
  cudaMalloc((void **)&image_gpu5, sizeof(int) * detector_y * detector_x * num_proj);
  cudaMemcpy(image_gpu5, image5, sizeof(int) * detector_y * detector_x * num_proj, cudaMemcpyHostToDevice);

  //int i,
  int count = 0, sizex, sizey;
  int start_keV = 140, end_keV = 140;
  float sum_length = 0.0;
  float dens_H2O = 1.0; //水の密度i
  float dens_Ca = 1.550;

  /*線減衰係数の値をmuに代入*/
  sizex = 4; //線減衰係数配列のインデックス
  sizey = 200;

  readXcom();
  readxray();

  float *ab_H2O_gpu;
  float *coh_H2O_gpu;
  float *com_H2O_gpu;
  float *mua_H2O_gpu;

  float *ab_Ca_gpu;
  float *coh_Ca_gpu;
  float *com_Ca_gpu;
  float *mua_Ca_gpu;

  float *al2mm_gpu;
  float *al10mm_gpu;

  cudaMalloc((void **)&ab_H2O_gpu, sizeof(float) * 201);
  cudaMemcpy(ab_H2O_gpu, ab_H2O, sizeof(float) * 201, cudaMemcpyHostToDevice);

  cudaMalloc((void **)&coh_H2O_gpu, sizeof(float) * 201);
  cudaMemcpy(coh_H2O_gpu, coh_H2O, sizeof(float) * 201, cudaMemcpyHostToDevice);

  cudaMalloc((void **)&com_H2O_gpu, sizeof(float) * 201);
  cudaMemcpy(com_H2O_gpu, com_H2O, sizeof(float) * 201, cudaMemcpyHostToDevice);

  cudaMalloc((void **)&mua_H2O_gpu, sizeof(float) * 201);
  cudaMemcpy(mua_H2O_gpu, mua_H2O, sizeof(float) * 201, cudaMemcpyHostToDevice);

  cudaMalloc((void **)&ab_Ca_gpu, sizeof(float) * 201);
  cudaMemcpy(ab_Ca_gpu, ab_Ca, sizeof(float) * 201, cudaMemcpyHostToDevice);

  cudaMalloc((void **)&coh_Ca_gpu, sizeof(float) * 201);
  cudaMemcpy(coh_Ca_gpu, coh_Ca, sizeof(float) * 201, cudaMemcpyHostToDevice);

  cudaMalloc((void **)&com_Ca_gpu, sizeof(float) * 201);
  cudaMemcpy(com_Ca_gpu, com_Ca, sizeof(float) * 201, cudaMemcpyHostToDevice);

  cudaMalloc((void **)&mua_Ca_gpu, sizeof(float) * 201);
  cudaMemcpy(mua_Ca_gpu, mua_Ca, sizeof(float) * 201, cudaMemcpyHostToDevice);

  cudaMalloc((void **)&al2mm_gpu, sizeof(float) * 251);
  cudaMemcpy(al2mm_gpu, al2mm, sizeof(float) * 251, cudaMemcpyHostToDevice);

  cudaMalloc((void **)&al10mm_gpu, sizeof(float) * 251);
  cudaMemcpy(al10mm_gpu, al10mm, sizeof(float) * 251, cudaMemcpyHostToDevice);

  //int keV=140;
  //min_e=start_keV;

  int per = num_photon;

  int block_x_n = 1;
  int block_y_n = 1;

  //一番近い2のべき乗を求める
  while (true)
  {
    block_x_n *= 2;
    if (block_x_n >= detector_x)
    {
      break;
    }
  }

  while (true)
  {
    block_y_n *= 2;
    if (block_y_n >= detector_y)
    {
      break;
    }
  }
  cout << "\n blockx, block_y = (" << block_x_n << ", " << block_y_n << ")";

  // dim3変数の宣言
  dim3 blocks(block_x_n / 8, block_y_n / 8);
  dim3 threads(8, 8);


  for (int q = ScatterNUM; q <= ScatterNUM; q++)
  {
    for (int keV = start_keV; keV <= end_keV; keV++)
    {
      float Energy = start_keV;
      float mu_H2O = mua_H2O[(int)(Energy + 0.5)] * dens_H2O; //csv_H2O[3][(int)(Energy+0.5)]*dens_H2O;
      float mu_Ca = mua_Ca[(int)(Energy + 0.5)] * dens_Ca;    //csv_Ca[3][(int)(Energy+0.5)]*dens_Ca;
      float mu_Max = (mu_H2O <= mu_Ca) ? mu_Ca : mu_H2O;

      sum_length = 0;
      cout << "--------------------" << keV << "keV"
           << "--------------------\n";
      cout << "(mu_H2O, mu_Ca, mu_Max) = (" << mu_H2O << "," << mu_Ca << "," << mu_Max << ")" << endl;

      count = 0;

      //float theta=0., phi=0.;
      //float min_e=140.;
      //int primary_num=0;

      //int minus_num=0;
      //int coh_num=0;
      //int projection_step = 90;

      /*photon p;
      //p = new photon;
      photon p_gpu;
      //p_gpu = ;
      cudaMalloc((void**)&p_gpu, sizeof(photon) * 1);
      cudaMemcpy(p_gpu , &p ,sizeof(photon) * 1, cudaMemcpyHostToDevice);*/
      //配列無い時は&が必要

      //シード値
      curandStateMRG32k3a *state_gpu;
      cudaMalloc((void **)&state_gpu, sizeof(curandStateMRG32k3a) * detector_x * detector_y);

      RandStateGenerator<<<blocks, threads>>>(state_gpu);
      cudaThreadSynchronize();

      projection<<<blocks, threads>>>(per, mu_H2O, mu_Ca, geometry_gpu, image_gpu0, image_gpu5, state_gpu,
                                      ab_H2O_gpu, coh_H2O_gpu, com_H2O_gpu, mua_H2O_gpu, ab_Ca_gpu, coh_Ca_gpu, com_Ca_gpu, mua_Ca_gpu, al10mm_gpu);

      cudaThreadSynchronize();

      cudaMemcpy(image_out0, image_gpu0, sizeof(int) * detector_x * detector_y * num_proj, cudaMemcpyDeviceToHost);

      cudaMemcpy(image_out5, image_gpu5, sizeof(int) * detector_x * detector_y * num_proj, cudaMemcpyDeviceToHost);

      cout << "-----result-----\n";

      /*平均光路長の出力*/
      cout << "count = " << count << " / " << H * W * per << "\n\n";

      //sprintf(wname,"%s.raw",writeFileName);
      //cout<<"spok"<<endl;

      float *map_out0 = (float *)calloc(H * W * 1 * num_proj, sizeof(float)); //処理用画像65*65*360(real)
      float *map_out5 = (float *)calloc(H * W * 1 * num_proj, sizeof(float)); //処理用画像65*65*360(real)

      for (int a = 0; a < num_proj; a++)
      {
        for (int b = 0; b < H; b++)
        {
          for (int c = 0; c < W; c++)
          {
            if (image_out0[a * H * W + b * H + c] > per)
            { //||image_out[a*H*W + b*H + c]>2000){
              image_out0[a * H * W + b * H + c] = per;
            }
            if (image_out5[a * H * W + b * H + c] > per)
            { //||image_out[a*H*W + b*H + c]>2000){
              image_out5[a * H * W + b * H + c] = per;
            }
            if (image_out0[a * H * W + b * H + c] == 0)
            {
              image_out0[a * H * W + b * H + c] = 1;
            }
            if (image_out5[a * H * W + b * H + c] == 0)
            {
              image_out5[a * H * W + b * H + c] = 1;
            }
            //cout<<image_out[a*65*65 + b*65 + c]<<" "<<-log(image_out[a*65*65 + b*65 + c])+log(2000)<<endl;
            map_out0[a * H * W + b * H + c] = -log(image_out0[a * H * W + b * H + c]) + log(float(per));
            map_out5[a * H * W + b * H + c] = -log(image_out5[a * H * W + b * H + c]) + log(float(per));
          }
        }
      }

      writeRawFile(writeFileName0, sizeof(int), H * W * 1 * num_proj, image_out0);
      writeRawFile(writeFileName_m0, sizeof(float), H * W * 1 * num_proj, map_out0);

      writeRawFile(writeFileName5, sizeof(int), H * W * 1 * num_proj, image_out5);
      writeRawFile(writeFileName_m5, sizeof(float), H * W * 1 * num_proj, map_out5);
    }
  }
}

void readXcom()
{
  FILE *fp_h;
  fp_h = fopen(H2O_c, "r");

  if (fp_h == NULL)
  {
    printf("failed to open %s\n", H2O_c);
    exit(-1);
  }
  else
  {
    printf("%s \n", H2O_c);
    float buf[4];

    for (int i = 0; i <= 200; i++)
    {
      fscanf(fp_h, "%f\t%f\t%f\t%f", &coh_H2O[i], &com_H2O[i], &ab_H2O[i], &mua_H2O[i]);
      //printf("%f\t%f\t%f\t%f", coh_H2O[i], com_H2O[i], ab_H2O[i], mua_H2O[i]);
    }
  }
  fclose(fp_h);

  FILE *fp_c;
  int ret2 = 1;

  fp_c = fopen(Ca_c, "r");
  if (fp_c == NULL)
  {
    printf("failed to open %s\n", Ca_c);
    exit(-1);
  }
  printf("%s \n", Ca_c);

  for (int j = 0; j <= 200; j++)
  {
    //cout<<fp<<" ";
    fscanf(fp_c, "%f\t%f\t%f\t%f", &coh_Ca[j], &com_Ca[j], &ab_Ca[j], &mua_Ca[j]);
    //cout<<ab_Ca[j]<<" "<<coh_Ca[j]<<" "<<com_Ca[j]<<" "<<mua_Ca[j]<<endl;
  }

  fclose(fp_c);
}

void readxray()
{
  FILE *fp_h;
  int ret1 = 1;

  fp_h = fopen("125kv_al2mm.txt", "r");
  //printf("%d \n",fp_h[0]);
  if (fp_h == NULL)
  {
    printf("failed to open %s\n", "125kv_2mm.txt");
    exit(-1);
  }
  printf("125kv_al2mm.txt \n");

  for (int i = 0; i <= 250; i++)
  {
    fscanf(fp_h, "%f", &al2mm[i]);
    //printf("%f ",al2mm[i]);
  }

  fclose(fp_h);

  FILE *fp_c;
  int ret2 = 1;

  fp_c = fopen("125kv_al10mm.txt", "r");
  if (fp_c == NULL)
  {
    printf("failed to open 125kv_al10mm.txt\n");
    exit(-1);
  }
  printf("125kv_al10mm.txt ");

  for (int j = 0; j <= 250; j++)
  {
    fscanf(fp_h, "%f", &al10mm[j]);
    //printf("%f ",al10mm[j]);
  }

  fclose(fp_c);
}

void writeRawFile(const char fname[], const size_t size, const size_t num, void *image)
{
  // ファイルを開く
  FILE *fp = fopen(fname, "wb");

  // ファイルを開くことができなかった場合のエラー処理
  if (NULL == fp)
  {
    printf("failed to open %s\n", fname);
    exit(-1);
  }

  // データの書き出し
  //cout<<"loadhvoval";
  size_t ret = fwrite(image, size, num, fp);

  // データを書き込むことができなかった場合のエラー処理
  if (num != ret)
  {
    printf("failed to write %s\n", fname);
    fclose(fp);
    exit(-1);
  }

  // ファイルを閉じる
  fclose(fp);
}

__global__ void projection(int per, float mu_H2O, float mu_Ca, unsigned char *geometry, int *image0, int *image5, curandStateMRG32k3a *state_gpu,
                           float *ab_H2O, float *coh_H2O, float *com_H2O, float *mua_H2O,
                           float *ab_Ca, float *coh_Ca, float *com_Ca, float *mua_Ca, float *al10mm)
{
  //関数にditector_x, ditector_y渡す
  int i = blockIdx.y * blockDim.y + threadIdx.y;
  int j = blockIdx.x * blockDim.x + threadIdx.x;

  int countp[325] = {0};

  float dens_H2O = 1.0; //水の密度i
  float dens_Ca = 1.550;

  int num_scatter = 0;

  if (j >= detector_x || i >= detector_y)
  {
    return;
  }
  int err = 0;

  int s_index = detector_y * i + j;

  curandStateMRG32k3a st = state_gpu[s_index];

  int projection_step = 1;
  int num_add = 0;

  const float detector_height = 16.25, Dso = 160., Dod = 60., Dsd = 220., pixel_size_d = 0.1;
  const float start_fantom = 10.1;

  for (int num_p = 0; num_p < num_proj; num_p += projection_step)
  { //投影数(xy平面で回転)for(int num_proj = 0; num_proj)
    for (int num_ph = 0; num_ph < per; num_ph++)
    {
      float cos_theta_a_new = 1., cos_phi_a_new = 0., sin_theta_a_new = 0., sin_phi_a_new = 1.;
      float cos_theta_a = 0., sin_theta_a = 1.;
      photon p1 = photon();
      //p1 = new photon;//p_gpu[0];

      p1.x = -Dso; //-16.25*sqrt(3)/2;
      p1.y = 0;
      p1.z = 0;

      //zlをマイナスに
      float xl = Dod,
            yl = detector_height - 0.5 * pixel_size_d /*ピクセル中心*/ - pixel_size_d * i,
            zl = detector_height - 0.5 * pixel_size_d - pixel_size_d * j;
            /*if(num_p == 0 && num_ph == 0 && j == 0){
              //atomicAdd(&countp[i],1);
              printf("%f, %f \n",yl,zl);
            }*/

      float cos_phi_a, sin_phi_a;
      float theta_a = 0.5 * M_PI - atan(zl / Dsd);

      float cos_theta_a2 = cos(theta_a), sin_theta_a2 = sin(theta_a);
      cos_theta_a = cos_theta_a2, sin_theta_a = sin_theta_a2;


      float Energy = 140; //start_keV
      float e_rand = curand_uniform(&st);
      for(int dec_e = 0; dec_e <250; dec_e++){
        if( al10mm[dec_e] <= e_rand  && e_rand <= al10mm[dec_e+1]){
          Energy = (dec_e+1)*0.5;
          break;
        }
      }
      //printf("Energy: %f", Energy);

      mu_H2O = mua_H2O[(int)(Energy + 0.5)] * dens_H2O;
      mu_Ca = mua_Ca[(int)(Energy + 0.5)] * dens_Ca;

      //printf("%f\n",mua_H2O[(int)(Energy + 0.5)]);
      //mu_H2O=csv_H2O[3][(int)(Energy+0.5)]*dens_H2O;-----------------------?
      //mu_Ca=csv_Ca[3][(int)(Energy+0.5)]*dens_Ca;

      float phia; // = (29.95-(59.9*j/64.))*M_PI/180.;
      //pixel中心から
      phia = atan(yl / Dsd);
      phia += M_PI * num_p / 180.;

      //int ditector_index = 180*phia/(2*M_PI);//検出器の検出画素のインデックス
      sin_phi_a = sin(phia);
      cos_phi_a = cos(phia);
      //cos_phi_a_new = cos_phi_a, sin_phi_a_new = sin_phi_a;

      //線源の回転
      float primary_x = p1.x * cos(M_PI * num_p / 180.) - p1.y * sin(M_PI * num_p / 180.);
      float primary_y = p1.x * sin(M_PI * num_p / 180.) + p1.y * cos(M_PI * num_p / 180.);

      p1.x = primary_x;
      p1.y = primary_y;

      //x=-start_fantomの位置まで飛ばす
      float photon_vec[3];                                                                      //光子の座標が入るベクトル
      photon_vec[0] = Dod * cos(M_PI * num_p / 180.) - yl * sin(M_PI * num_p / 180.) - primary_x; //60 - primary_x;//x
      photon_vec[1] = Dod * sin(M_PI * num_p / 180.) + yl * cos(M_PI * num_p / 180.) - primary_y; //yl - primary_y;//y
      photon_vec[2] = zl;

      //ここ書き換え必要
      float to = float(Dso - start_fantom) / Dsd;

      float to_ditector_x = p1.x + photon_vec[0];
      float to_ditector_y = p1.y + photon_vec[1];

      for (int to_n = 0; to_n < 3; to_n++){
        photon_vec[to_n] *= to;
      }

      p1.x += photon_vec[0];
      p1.y += photon_vec[1];
      p1.z += photon_vec[2];
      if(isnan(p1.x)||isnan(p1.y)||isnan(p1.z)){
        printf("%f,%f,%f\n",p1.x,p1.y,p1.z);
      }

      if(isnan(sin_theta_a2)||isnan(cos_theta_a2)||isnan(sin_phi_a)||isnan(cos_phi_a)){
        printf("%f,%f,%f, %f\n",sin_theta_a2, cos_theta_a2, sin_phi_a, cos_phi_a);
      }

      p1.delta_sampling(mu_H2O, mu_Ca, geometry, sin_theta_a2, cos_theta_a2, sin_phi_a, cos_phi_a, &st);
      //printf("%f\n\n",p1.x);
      //線源を元の座標に戻すために逆回転
      float tmp_x = p1.x;
      float tmp_y = p1.y;
      float tmp_z = p1.z;

      float x_r = tmp_x * cosf((-num_p * M_PI) / 180.) - tmp_y * sinf((-num_p * M_PI) / 180.);
      float y_r = tmp_x * sinf((-num_p * M_PI) / 180.) + tmp_y * cosf((-num_p * M_PI) / 180.);

      //p1.x =  tmp_x*cos(M_PI*-num_p/180) - tmp_y*sin(M_PI*-num_p/180);
      //p1.y =  tmp_x*sin(M_PI*-num_p/180) + tmp_y*cos(M_PI*-num_p/180);

      //if(num_ph==0&&num_p == 0){printf("%d %d %f %f %f %f %f\n", i, j , p1.x, p1.y, p1.z, to_ditector_y, to_ditector_x);}

      /*if(i==32&&j==32){
        printf(" %f %f ", x_r,y_r);
      }*/
      if (x_r >= Dod){ //1回目で検出器到着---------------primary-----------------&&yの範囲

        p1.x = x_r;
        p1.y = y_r;
        float d_z = ((p1.z - 0) / (p1.x - (-160))) * 60 + 160 * p1.z / (p1.x + 160);
        float d_y = ((p1.y - 0) / (p1.x - (-160))) * 60 + 160 * p1.y / (p1.x + 160);

        int result_y = -1 * (int(d_y * 10 - 325./2.)); //x,y変更
        int result_x = -1 * (int(d_z * 10 - 325./2.));
        if (result_y < 0 || result_y > 325 || result_x < 0 || result_x > 325)
        {
          err++;
        }

        //if(ScatterNUM == 0)
        //{
        num_add++;
        atomicAdd(&image0[num_p * detector_y * detector_x + result_y * detector_y + result_x], 1);
        atomicAdd(&image5[num_p * detector_y * detector_x + result_y * detector_y + result_x], 1);
        //if(i==32&&j==32)printf("%d add\n",num_add);
        //}

        //break;
      }
      else
      { //散乱
        num_scatter++;
        //if(i==32&&j==32)
        //printf("%d 散乱  ",num_ph);
        bool coh_flag = false;
        bool com_flag = false;

        for (int a = 0; a < ScatterNUM; a++)
        {
          //state_gpu[s_index] = st;///////////?
          float dens, ab, coh, com, mu;

          if (isnan(cos_theta_a) || isnan(cos_phi_a) || isnan(sin_theta_a) || isnan(sin_phi_a) ||
              isnan(cos_theta_a_new) || isnan(cos_phi_a_new) || isnan(sin_theta_a_new) || isnan(sin_phi_a_new))
          {
            printf("a %d %d %d \t",a, num_p,num_ph);
            printf("%f %f %f %f %f %f %f %f\n",cos_theta_a,cos_phi_a,sin_theta_a,sin_phi_a
            ,cos_theta_a_new, cos_phi_a_new, sin_theta_a_new, sin_phi_a_new);
            break;
          }

          float x_rotate_c = p1.x * cosf(-num_p * M_PI / 180) - p1.y * sinf(-num_p * M_PI / 180);
          float y_rotate_c = p1.x * sinf(-num_p * M_PI / 180) + p1.y * cosf(-num_p * M_PI / 180);
          if (x_rotate_c >= 60 || abs(y_rotate_c) >= 16.25 || abs(p1.z) >= 16.25)
          {
            //printf("break");
            break;
          }
          
          //if(pow((p1.y-5),2)+pow(p1.z,2) <= 9||pow((p1.y+3),2)+pow((p1.z+3),2) <= 9)

          if(pow((p1.y-5),2)+pow(p1.z,2) <= pow(1.5,2) || pow((p1.y+5),2)+pow(p1.z,2) <= pow(1.5,2)
	          || pow(p1.y,2)+pow((p1.z-5),2) <= pow(1.5,2) || pow(p1.y,2)+pow((p1.z+5),2) <= pow(1.5,2)
	          || pow((p1.y-5*sin(M_PI/4.)),2)+pow((p1.z-5*cos(M_PI/4.)),2) <= pow(1.5,2) || pow((p1.y+5*sin(M_PI/4.)),2)+pow((p1.z+5*cos(M_PI/4.)),2) <= pow(1.5,2)
	          || pow((p1.y-5*sin(M_PI*3./4.)),2)+pow((p1.z-5*cos(M_PI*3./4.)),2) <= pow(1.5,2) || pow((p1.y+5*sin(M_PI*3./4.)),2)+pow((p1.z+5*cos(M_PI*3./4.)),2) <= pow(1.5,2)
	        )
          {
            dens = dens_Ca;
            ab = ab_Ca[(int)(Energy+0.5)];
            coh = coh_Ca[(int)(Energy+0.5)];
            com = com_Ca[(int)(Energy+0.5)];
            mu =  mua_Ca[(int)(Energy+0.5)];
          }
          else
          {
            dens = dens_H2O;
            ab = ab_H2O[(int)(Energy + 0.5)];
            coh = coh_H2O[(int)(Energy + 0.5)];
            com = com_H2O[(int)(Energy + 0.5)];
            mu = mua_H2O[(int)(Energy + 0.5)];
          }

          float sc_rand = curand_uniform(&st); //////////?

          if (sc_rand <= ab / mu)
          { //光電効果
            //if(i==32&&j==32)printf("光電 ");
            break;
          }
          if (ab / mu < sc_rand && sc_rand <= (ab + coh) / mu)
          { //コヒーレント散乱----------------------------------------------------------------------------
            //if(i==32&&j==32)printf("coh ");
            //coh_num++;
            //num_scatter--;

            if (com_flag == false)
            {
              p1.delta_sampling(mu_H2O, mu_Ca, geometry, sin_theta_a2, cos_theta_a2, sin_phi_a, cos_phi_a, &st);
              coh_flag = true;
            }
            else
            { //過去にコンプトンあるなら，その時の角度を参照
              p1.delta_sampling(mu_H2O, mu_Ca, geometry, sin_theta_a_new, cos_theta_a_new, sin_phi_a_new, cos_phi_a_new, &st);
            }

            float x_rotate_ch = p1.x * cosf((-num_p * M_PI) / 180.) - p1.y * sinf((-num_p * M_PI) / 180.);
            float y_rotate_ch = p1.x * sinf((-num_p * M_PI) / 180.) + p1.y * cosf((-num_p * M_PI) / 180.);
            //printf("x_rotate_ch");

            float x_p_rotate_ch = p1.x_p * cosf((-num_p * M_PI) / 180.) - p1.y_p * sinf((-num_p * M_PI) / 180.);
            float y_p_rotate_ch = p1.x_p * sinf((-num_p * M_PI) / 180.) + p1.y_p * cosf((-num_p * M_PI) / 180.);
            //float x_p_rotate_ch = p1.x_p*cosf(M_PI*-num_p/180.) - p1.y_p*sin(M_PI*-num_p/180.);
            //float y_p_rotate_ch = p1.x_p*sin(M_PI*-num_p/180.) + p1.y_p*cos(M_PI*-num_p/180.);

            //float ditector_y=(10-x_p_rotate)*(y_rotate-y_p_rotate)/(x_rotate-x_p_rotate)+y_p_rotate;
            //float d_z_ch = ((p1.z - p1.z_p)/(x_rotate_ch-x_p_rotate_ch))*60 + (x_rotate_ch*p1.z_p - x_p_rotate_ch*p1.z)/(x_rotate_ch - x_p_rotate_ch);
            //float d_y_ch = ((y_rotate_ch - y_p_rotate_ch)/(x_rotate_ch-x_p_rotate_ch))*60 + (x_rotate_ch*y_p_rotate_ch - x_p_rotate_ch*y_rotate_ch)/(x_rotate_ch - x_p_rotate_ch);
            float d_z_ch = ((p1.z - p1.z_p) / (x_rotate_ch - x_p_rotate_ch)) * 60 + (x_rotate_ch * p1.z_p - x_p_rotate_ch * p1.z) / (x_rotate_ch - x_p_rotate_ch);
            float d_y_ch = ((y_rotate_ch - y_p_rotate_ch) / (x_rotate_ch - x_p_rotate_ch)) * 60 + (x_rotate_ch * y_p_rotate_ch - x_p_rotate_ch * y_rotate_ch) / (x_rotate_ch - x_p_rotate_ch);

            int result_x_ch = -1 * (int(d_z_ch * 10 - 325./2.)); //x,y変更
            int result_y_ch = -1 * (int(d_y_ch * 10 - 325./2.));

            if (x_rotate_ch >= 60 && result_y_ch <= 325 && result_x_ch <= 325)// && com_flag == true
            {
              atomicAdd(&image5[num_p * detector_y * detector_x + result_y_ch * detector_y + result_x_ch], 1);
              break;
            }
          }
          else
          { //コンプトン散乱
            //if(i==32&&j==32)printf("com ");//\n";

            //1:散乱角，エネルギー計算
            float lambda = 511.0 / Energy;

            float lambda_d = 0.;
            bool track_flag = true;

            while (track_flag)
            {
              //double r1=genrand_real3();
              float r1 = curand_uniform(&st);
              //r1=0.1;
              if (r1 < (lambda + 2.0) / (9.0 * lambda + 2.0))
              { //track1 <=or<?
                float r2 = curand_uniform(&st);
                //r2=0.2;
                float ro = 1.0 + (2.0 / lambda) * r2;
                float r3 = curand_uniform(&st);
                //r3=0.3;

                if (r3 <= 4.0 * ((1. / ro) - (1. / (ro * ro))))
                {
                  lambda_d = ro * lambda;
                  track_flag = false;
                }
              }
              else
              { //track2
                float r2 = curand_uniform(&st);
                float ro = (lambda + 2.) / (lambda + 2. * (1. - r2));
                float r3 = curand_uniform(&st);
                if (r3 <= 0.5 * (pow((lambda - ro * lambda + 1.), 2) + (1. / ro)))
                {
                  lambda_d = ro * lambda;
                  track_flag = false;
                }
              }
            }
            //lambda_d=lambda+0.1;

            float theta = acos(1. - (lambda_d - lambda)); ///----------------------何かまずいかも
            if (isnan(theta))
            {
              //printf("mz");
            }
            //float theta = atan(((lambda_d - lambda)*(lambda_d - lambda))/(lambda_d - lambda))
            float cos_theta = (1 - (lambda_d - lambda)); //cos(theta);//cos(0.5*M_PI + atan(50./220.));
            if (cos_theta < -1)
              cos_theta = -1;
            //0.5*_PIだとnanに
            float sin_theta = sqrt(1. - pow((cos_theta), 2));

            if (abs(cos_theta) > 1.0)
            {
              //cos_theta+=0.01;
              //printf(" %f %f\n", cos_theta,sin_theta);
            }

            Energy = 511. / lambda_d;
            //mu_H2O = csv_H2O[3][(int)(Energy+0.5)]*dens_H2O;
            //mu_Ca =  csv_Ca[3][(int)(Energy+0.5)]*dens_Ca;
            mu_H2O = mua_H2O[(int)(Energy + 0.5)] * dens_H2O;
            mu_Ca = mua_Ca[(int)(Energy + 0.5)] * dens_Ca;

            //2：方位角，光路長計算
            float phi = curand_uniform(&st) * 2. * M_PI;
            //乗算では.を忘れるな

            //3:相対座標→絶対座標
            if (com_flag == true)
            {
              sin_theta_a = sin_theta_a_new;
              cos_theta_a = cos_theta_a_new;
              cos_phi_a = cos_phi_a_new;
            }

            cos_theta_a_new = cos_theta_a * cos_theta - sin_theta_a * sin_theta * cos(phi); //cos(0.5*M_PI) -
            if (cos_theta_a_new < -1)cos_theta_a_new = -1;
            sin_theta_a_new = sqrt(1. - pow(cos_theta_a_new, 2)); //絶対座標系の新しい角度

            cos_phi_a_new = (cos_theta_a * cos_phi_a * sin_theta * cos(phi) + sin_theta_a * cos_phi_a * cos_theta - sin_phi_a * sin_theta * sin(phi)) / sin_theta_a_new;
            sin_phi_a_new = (cos_theta_a * sin_phi_a * sin_theta * cos(phi) + sin_theta_a * sin_phi_a * cos_theta + cos_phi_a * sin_theta * sin(phi)) / sin_theta_a_new;

            if (isnan(cos_phi_a_new) || isnan(sin_phi_a_new) || isnan(cos_theta_a_new) || isnan(sin_theta_a_new))
            {
              printf("%f, %f\n", cos_theta_a, sin_theta_a);
              printf("a: %d theta: %f, %f  phi; %f , %f\n", a, cos_theta_a_new, sin_theta_a_new, cos_phi_a_new, sin_phi_a_new);
              if (com_flag)
              {
                printf("com\n");
              }
              break;
            }

            com_flag = true;

            p1.delta_sampling(mu_H2O, mu_Ca, geometry, sin_theta_a_new, cos_theta_a_new, sin_phi_a_new, cos_phi_a_new, &st);
            //delta_sampling(p, mu_H2O, mu_Ca, geometry, sin_theta_a_new, cos_theta_a_new, sin_phi_a_new, cos_phi_a_new);

            p1.theta = theta;

            float v_length = sqrt(pow(p1.length * sin_theta_a_new * cos_phi_a_new, 2) + pow(p1.length * sin_theta_a_new * sin_phi_a_new, 2) + pow(p1.length * cos_theta_a_new, 2));
            p1.before_vec0 = p1.length * sin_theta_a_new * cos_phi_a_new / v_length;
            p1.before_vec1 = p1.length * sin_theta_a_new * sin_phi_a_new / v_length;
            p1.before_vec2 = p1.length * cos_theta_a_new / v_length;
            //vectorは長さ1に正規化

            //cout<<ditector_index<<" "<<phi_a_result<<endl;

            //float x_rotate=0, y_rotate=0, x_p_rotate=0, y_p_rotate=0;
            //座標の回転は今ついてる角度と逆方向に回さねば
            /*phi_a_result = num_p;//M_PI*210./180.;
              x_rotate = p1.x*cos(-phi_a_result) - p1.y*sin(-phi_a_result);
              y_rotate = p1.x*sin(-phi_a_result) + p1.y*cos(-phi_a_result);
              x_p_rotate = p1.x_p*cos(-phi_a_result) - p1.y_p*sin(-phi_a_result);
              y_p_rotate = p1.x_p*sin(-phi_a_result) + p1.y_p*cos(-phi_a_result);*/

            //cout<<p->x<<" "<<p->y<<" "<<p->z<<endl;

            //detector_y=(10-x_p_rotate)*(y_rotate-y_p_rotate)/(x_rotate-x_p_rotate)+y_p_rotate;//x==10の時のy座標

            //float d_z=(60-p1.x_p)*(p1.z-p1.z_p)/(p1.x-p1.x_p)+p1.z_p;
            //float d_y=(60-p1.x_p)*(p1.y-p1.y_p)/(p1.x-p1.x_p)+p1.y_p;

            float x_rotate = p1.x * cosf((-num_p * M_PI) / 180.) - p1.y * sinf((-num_p * M_PI) / 180.); //phi_a_resultから変更
            float y_rotate = p1.x * sinf((-num_p * M_PI) / 180.) + p1.y * cosf((-num_p * M_PI) / 180.);
            float x_p_rotate = p1.x_p * cosf((-num_p * M_PI) / 180.) - p1.y_p * sinf((-num_p * M_PI) / 180.);
            float y_p_rotate = p1.x_p * sinf((-num_p * M_PI) / 180.) + p1.y_p * cosf((-num_p * M_PI) / 180.);

            float d_z = ((p1.z - p1.z_p) / (x_rotate - x_p_rotate)) * 60 + (x_rotate * p1.z_p - x_p_rotate * p1.z) / (x_rotate - x_p_rotate);
            float d_y = ((y_rotate - y_p_rotate) / (x_rotate - x_p_rotate)) * 60 + (x_rotate * y_p_rotate - x_p_rotate * y_rotate) / (x_rotate - x_p_rotate);

            if (x_rotate >= 60 && abs(d_z) <= 16.25 && abs(d_y) <= 16.25)
            { //&&a==q-1){// && abs(p->z)<16.25 && abs(y_rotate)<16.25){// && x_rotate>10
              //検出器を通過したかcheck,今回は必ず検出器まで到達する

              int result_y = -1 * (int(d_y * 10 - 325./2.));
              int result_x = -1 * (int(d_z * 10- 325./2.));

              //count++;

              //image[+result_y*65+result_x]++;
              atomicAdd(&image5[num_p * detector_y * detector_x + result_y * detector_y + result_x], 1);
              break;
            }
          }
        }
      }
    }
  }
  state_gpu[s_index] = st;
  
  if (err != 0)
  {
    printf("err: %d", err);
  }
  //printf(" %d ",num_add);

  /*for(int out = 0; out<325;out++){
    printf("%f\n", countp[out]);
  }*/
}

__device__ void photon::delta_sampling(float mu_H2O, float mu_Ca, unsigned char *geometry, float sin_theta_a, float cos_theta_a, float sin_phi_a, float cos_phi_a, curandStateMRG32k3a *st)
{
  //printf("%f\n",x);
  float mu_max;
  mu_max = max(mu_H2O, mu_Ca); //2媒質の時
  //mu_max = mu_Ca; //1媒質の場合
  //printf("%f ",mu_max);
  bool loop_flag = true;
  bool air_flag = true;
  float x1 = x, y1 = y, z1 = z, length1 = length;
  float x2 = x, y2 = y, z2 = z, length3 = 0;
  x_p = x1, y_p = y1, z_p = z1;

  int num_itr = 0;
  //int geo_index;
  int check;
  if (isnan(sin_theta_a) || isnan(cos_theta_a) || isnan(sin_phi_a) || isnan(cos_phi_a))
  {
    printf("%f, %f , %f, %f\n", sin_theta_a, cos_theta_a, sin_phi_a, cos_theta_a);
    return;
  }

  while (loop_flag)
  {
    float beta = curand_uniform(st);
    float r = -log(beta) / mu_max;
    //printf("%f",r);

    x2 += r * sin_theta_a * cos_phi_a;
    y2 += r * sin_theta_a * sin_phi_a;
    z2 += r * cos_theta_a;
    length3 += r;
    //printf("%f, %f ,%f \n", x2,y2,z2);

    check = 0;

    //cout<<mu_Ca/mu_max;
    float nu = curand_uniform(st); //genrand_real3();

    /*if (x2 * x2 + y2 * y2 + z2 * z2 <= 100)//半径10cm球
    {
      check = 1;
    }*/
    /*if(x2 * x2 + y2 * y2 + z2 * z2 <= 100){//半径10球右Ca
      if(z2>=0){//Ca
        check = 2;
      }
      else{//h2o
        check = 2;
      }
    }*/
    /*if(y2 * y2 + z2 * z2 <= 100 && -10 <= x2 && x2 <= 10){//半径10円柱in 3cm,2cm Ca
      if(pow(y2-5,2)+pow(z2,2) <= 9||pow(y2+5,2)+pow(z2,2) <= 9
      || pow(y2,2)+pow(z2-5,2) <= 9||pow(y2,2)+pow(z2+5,2) <= 9){
        check = 2;
      }
      else{
        check = 1;
      }
    }*/
    /*if(x2 * x2 + y2 * y2 + z2 * z2 <= 16){//半径4球in 1.5cm Ca球
      if(pow((x2 - 2),2) + pow((y2),2) + pow(z2,2) <= 2.25){
        check = 2;
      }
      else{
        check = 1;
      }
    }*/
  if(y2 * y2 + z2 * z2 <= 100 && -10 <= x2 && x2 <= 10){//半径10円柱in 3cm Ca
      if(pow((y2-5),2)+pow(z2,2) <= pow(1.5,2) || pow((y2+5),2)+pow(z2,2) <= pow(1.5,2)
	|| pow(y2,2)+pow((z2-5),2) <= pow(1.5,2) || pow(y2,2)+pow((z2+5),2) <= pow(1.5,2)
	|| pow((y2-5*sin(M_PI/4.)),2)+pow((z2-5*cos(M_PI/4.)),2) <= pow(1.5,2) || pow((y2+5*sin(M_PI/4.)),2)+pow((z2+5*cos(M_PI/4.)),2) <= pow(1.5,2)
	|| pow((y2-5*sin(M_PI*3./4.)),2)+pow((z2-5*cos(M_PI*3./4.)),2) <= pow(1.5,2) || pow((y2+5*sin(M_PI*3./4.)),2)+pow((z2+5*cos(M_PI*3./4.)),2) <= pow(1.5,2)
	){

        check = 2;
      }
      else{
        check = 1;
      }
    }



    if (check == 0)
    { //空気の時
      //空気の時，直進し続ける
      if (abs(x2) >= 62 || abs(y2) >= 62 || abs(z2) >= 17)
      {
        loop_flag = false;
      }
      else
      {
        num_itr++;
      }
    }

    else if (check == 1)
    { //H2Oの時
      if (nu <= mu_H2O / mu_max)
      {
        loop_flag = false;
        num_itr++;
        //break;
      }
    }
//#if 1
    else{//Ca領域の時
      if(nu <=  mu_Ca/mu_max){
        loop_flag = false;
        //break;
      }
    }
//#endif

    if (num_itr > 100)
    {
      //printf("%d %f, %f, %f\n",num_itr, x2,y2,z2);
    }
  }

  x = x2;
  y = y2;
  z = z2;

  length = length3;
  //printf("%f\n",x);
}

__global__ void RandStateGenerator(curandStateMRG32k3a *state_gpu)
{
  int a = blockIdx.y * blockDim.y + threadIdx.y;
  int b = blockIdx.x * blockDim.x + threadIdx.x;


  if (a >= detector_x || b >= detector_y)
  {
    return;
  }
  int index = detector_x * a + b;

  curand_init(0, index, 0, &state_gpu[index]);
}

__global__ void LaunchPhoton(curandStateMRG32k3a *state, int seed)
{
  int a = blockIdx.y * blockDim.y + threadIdx.y;
  int b = blockIdx.x * blockDim.x + threadIdx.x;

  int index = detector_x * a + b;
  curand_init(seed, index, 0, &state[index]);
}

void add_result(photon *p, float phi, int *image, int *es, int count, int ditector_index, int Energy, int a, int q)
{
  float x_rotate = p->x * cos(-phi) - p->y * sin(-phi);
  float y_rotate = p->x * sin(-phi) + p->y * cos(-phi);
  float x_p_rotate = p->x_p * cos(-phi) - p->y_p * sin(-phi);
  float y_p_rotate = p->x_p * sin(-phi) + p->y_p * cos(-phi);

  float ditector_y = (10 - x_p_rotate) * (y_rotate - y_p_rotate) / (x_rotate - x_p_rotate) + y_p_rotate;

  if (x_rotate >= 16.25)
  {
    if (a != q - 1)
    { //特定の散乱回数の場合のみ検出
      return;
    }
    if (abs(ditector_y) > 10)
    {
      //cout<<ditector_index<<" "<<ditector_y<<endl;
    }
    image[ditector_index * 65 + (int)(ditector_y * 2 + 32.5)]++;
    count++;
    es[(int)Energy]++;
    return;
  }
}