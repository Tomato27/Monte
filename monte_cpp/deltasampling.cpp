<謎の画像が出てしまう>
void delta_sampling(photon* p, double mu_H2O, double mu_Ca,unsigned char* geometry, double sin_theta_a, double cos_theta_a, double sin_phi_a, double cos_phi_a){
  double mu_max;
  //mu_max = max(mu_H2O, mu_Ca);
  mu_max=mu_H2O;//1媒質の場合
  bool loop_flag = true;
  bool air_flag = true;
  double x1 = p->x,y1= p->y,z1 = p->z, length1 = p->length;
  double x = p->x,y= p->y,z = p->z, length = 0;
  p->x_p = x1, p->y_p = y1, p->z_p = z1;
  //cout<<x<<" ";

  int l = 0;
  int geo_index;

  int check = 10;
  
  double to = 17 / (sin_theta_a * cos_phi_a);
  x +=  to * sin_theta_a * cos_phi_a;
  y +=  to * sin_theta_a * sin_phi_a;
  z +=  to * cos_theta_a;
  length += to;
  //cout<<x<<" ";

  while(loop_flag){
    double beta = genrand_real3();
    double r = -log(beta)/mu_max;
    /*if(air_flag){
      if(x<-6){
        r = 1;
      }
      else{
        r = 0.1;
      }
    }
    10.25*/
    x +=  r * sin_theta_a * cos_phi_a;
    y +=  r * sin_theta_a * sin_phi_a;
    z +=  r * cos_theta_a;
    length += r;

    //cout<<mu_Ca/mu_max;
    float nu = genrand_real3();

    if(x+23>=32.5||abs(y)>=16.25||abs(z)>=16.25){
      check = 0;
    }
    else{
      //geo_index = (int)(rint(rint(x)*2+rint(32.5))+rint(rint(y)*2+rint(32.5))*65+rint(rint(z)*2+rint(32.5))*65*65);
      //check = (int)geometry[geo_index];
      if(((x+16.25)*2-46)*((x+16.25)*2-46)+((y+16.25)*2-32.5)*((y+16.25)*2-32.5)+((z+16.25)*2-32.5)*((z+16.25)*2-32.5)<=100){
        check = 1;
      }
      //check=0;//-----------------------------------------
    }

    if(check==0){//空気の時
      if(x+23>25){
        air_flag=false;
      }
      //空気の時，直進し続ける
      if(x+23>=32.5||abs(y)>=25||abs(z)>=25){
        loop_flag = false;
      }
    }
    
    else if(check==1){//H2O領域の時
      if(air_flag){
        air_flag = false;
        continue;
      }
      //cout<<"h2o";
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
  }

  p->x = x;
  p->y = y;
  p->z = z;

  p->length = length;
}


<３つの輪っかが出来てしまう>
void delta_sampling(photon* p, double mu_H2O, double mu_Ca,unsigned char* geometry, double sin_theta_a, double cos_theta_a, double sin_phi_a, double cos_phi_a){
  double mu_max;
  //mu_max = max(mu_H2O, mu_Ca);
  mu_max=mu_H2O;//1媒質の場合
  bool loop_flag = true;
  bool air_flag = true;
  double x1 = p->x,y1= p->y,z1 = p->z, length1 = p->length;
  double x = p->x,y= p->y,z = p->z, length = 0;
  p->x_p = x1, p->y_p = y1, p->z_p = z1;

  int l = 0;
  int geo_index;

  while(loop_flag){
    double beta = genrand_real3();
    double r = -log(beta)/mu_max;

    if(air_flag){
      if(x<20){
        r = 1;
      }
      else{
        r = 0.05;
      }
    }
    x +=  r * sin_theta_a * cos_phi_a;
    y +=  r * sin_theta_a * sin_phi_a;
    z +=  r * cos_theta_a;
    length += r;

    int check=0;

    //cout<<mu_Ca/mu_max;
    float nu = genrand_real3();

    if(abs(x)>=16.25||abs(y)>=16.25||abs(z)>=16.25){
      check = 0;
    }
    else{
      geo_index = (int)(rint(rint(x)*2+rint(32.5))+rint(rint(y)*2+rint(32.5))*65+rint(rint(z)*2+rint(32.5))*65*65);
      //check = (int)geometry[geo_index];
      if(((x+16.25)*2-46)*((x+16.25)*2-46)+((y+16.25)*2-32.5)*((y+16.25)*2-32.5)+((z+16.25)*2-32.5)*((z+16.25)*2-32.5)<=100){
        check = 1;
      }
      //check=0;//-----------------------------------------
    }

    if(check==0){//空気の時
      if(x>25){
        air_flag=false;
      }
      //空気の時，直進し続ける
      if(abs(x)>=25||abs(y)>=25||abs(z)>=25){
        loop_flag = false;
      }
    }
    
    else if(check==1){
      if(air_flag){
        air_flag = false;
        continue;
      }
      //cout<<"h2o";
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
  }

  p->x = x;
  p->y = y;
  p->z = z;

  p->length = length;
}