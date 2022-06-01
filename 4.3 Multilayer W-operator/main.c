///////////////////////////////////////////////////////////////////////////////////////////////
//////////A data-driven systematic, consistent and feasible approach to Model Selection////////
//////////4.3 Multilayer W-operator                                                    ////////
//////////PhD Thesis, Diego Marcondes                                                  ////////
//////////Universisty of São Paulo, 2022                                               ////////
///////////////////////////////////////////////////////////////////////////////////////////////

//gcc main.c -o mlwo -fopenmp

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <string.h>
#include "input.h"
#include "utils.h"
#include "matrix_utils.h"
#include "window_utils.h"
#include <omp.h>

int main(){
    //Declare objects
    int ***train,***val,***W,***currentW,****Wtrain,****Wval,****currentWtrain,****currentWval,**ytrain,**yval,*wlength,*currentwlength,i,j,k,l,nlayer,keep,increase,digit,ntest;
    float *error_step,*error_nei;
    long start, time_read_data, time_Ucurve_start, time_digit_start, time_digit_end;
    char name_tmp[100],****joint,****currentjoint,binary[2][2],char_tmp[100],**buffer,tname[64];
    struct timeval timecheck;
    struct tm *timenow;
    time_t now = time(NULL);

    //Allocate vector for error and for the name of folder to save trace of algorithm
    error_nei = alloc_vector_float(2);
    buffer = alloc_vector_char(1,64);

    //Start random and create file to store trace of algorithm
    srand(time(NULL));
    FILE *fp, *Wfile;

    //Erase files of broken runs
    system("rm ./trace/*.txt");

    //Welcome message
    printf("\nWelcome to multilayerWOperators v0.1 Beta (November 2021).\n");
    printf("Developed by Diego Marcondes (dmarcondes@ime.usp.br).\n\n");

    //Start measuring time
    gettimeofday(&timecheck, NULL);
    start = (long)timecheck.tv_sec * 1000 + (long)timecheck.tv_usec / 1000;

    //Number of pixels to increase image so it is divisble by window size
    increase = (wsize - dim % (2 * (wsize/2))) % (2 * (wsize/2));

    //Increase dimension of images
    dim = dim + increase;

    //Number of layers of W-operator
    nlayer = dim/(2 * (wsize/2));

    //Start binary vector to create domain of each W-operator
    strcpy(binary[0],"0");
    strcpy(binary[1],"1");

    //File to save results
    fp = fopen(file_name, "w");

    //Read data
    printf("Reading data...\n");

    //Alocate array to save input and output of train and validation samples
    train = (int***) malloc(train_size*sizeof(int**));
    val = (int***) malloc(val_size*sizeof(int**));
    ytrain = alloc_matrix_int(train_size,1);
    yval = alloc_matrix_int(val_size,1);

    //Read each training image
    for(i = 0;i < train_size;i++){
      train[i] = alloc_matrix_int(dim,dim);
      strcpy(name_tmp,"");
      sprintf(name_tmp,"./mnist/train%06d.txt",i+1);
      readmatrix(dim,dim,train[i],name_tmp,increase);
    }

    //Read training label
    strcpy(name_tmp,"");
    sprintf(name_tmp,"./mnist/ytrain.txt");
    readmatrix(train_size,1,ytrain,name_tmp,0);

    //Read each validation image
    for(i = 0;i < val_size;i++){
      val[i] = alloc_matrix_int(dim,dim);
      strcpy(name_tmp,"");
      sprintf(name_tmp,"./mnist/val%06d.txt",i+1);
      readmatrix(dim, dim, val[i],name_tmp,increase);
    }

    //Read validation label
    strcpy(name_tmp,"");
    sprintf(name_tmp,"./mnist/yval.txt");
    readmatrix(val_size,1,yval,name_tmp,0);

    //Calculate time it took to read data
    gettimeofday(&timecheck, NULL);
    time_read_data = (long)timecheck.tv_sec * 1000 + (long)timecheck.tv_usec / 1000;
    printf("Data was read in %ld seconds.\n\n",(time_read_data-start)/1000);

    //Initialize window, current window, joint, current joint, window length and current window length
    W =  (int***) malloc(nlayer*sizeof(int**));
    currentW = (int***) malloc(nlayer*sizeof(int**));
    joint = (char****) malloc(nlayer*sizeof(char***));
    currentjoint = (char****) malloc(nlayer*sizeof(char***));
    wlength = alloc_vector_int(nlayer);
    currentwlength = alloc_vector_int(nlayer);

    if(random_init){
       //Create window and joint of each layer
      for(k = 0;k < nlayer;k++){
        //Each window has size 5
        wlength[k] = 5;
        currentwlength[k] = 5;

        //Alloc window of layer k
        W[k] = alloc_matrix_int(wsize*wsize,2);
        currentW[k] = alloc_matrix_int(wsize*wsize,2);

        //Atribute values
        W[k][0][0] = 0;
        W[k][0][1] = 0;
        W[k][1][0] = 0;
        W[k][1][1] = 1;
        W[k][2][0] = 0;
        W[k][2][1] = -1;
        W[k][3][0] = 1;
        W[k][3][1] = 0;
        W[k][4][0] = -1;
        W[k][4][1] = 0;

        currentW[k][0][0] = 0;
        currentW[k][0][1] = 0;
        currentW[k][1][0] = 0;
        currentW[k][1][1] = 1;
        currentW[k][2][0] = 0;
        currentW[k][2][1] = -1;
        currentW[k][3][0] = 1;
        currentW[k][3][1] = 0;
        currentW[k][4][0] = -1;
        currentW[k][4][1] = 0;

        //Alloc joint distribution of W operator in layer k
        joint[k] = alloc_matrix_char(2,(1 << wlength[k]),wsize*wsize);
        currentjoint[k] = alloc_matrix_char(2,(1 << wlength[k]),wsize*wsize);

        //Start with 0 and 1
        strcpy(joint[k][0][0],"0");
        strcpy(joint[k][0][1],"1");
        strcpy(currentjoint[k][0][0],"0");
        strcpy(currentjoint[k][0][1],"1");

        //Cartesian product with zero and one wlength[k] times to obtain operator domain
        for(i = 1;i < wlength[k];i++){
          cartesian(joint[k][0],binary,(1 << i),2,wsize*wsize);
          cartesian(currentjoint[k][0],binary,(1 << i),2,wsize*wsize);
        }

        //Randomly start w-operator
        for(i = 0;i < 32;i++){
          j = rand() % 2;
          if(j == 0){
            strcpy(joint[k][1][i],"0");
            strcpy(currentjoint[k][1][i],"0");
          }
          else{
            strcpy(joint[k][1][i],"1");
            strcpy(currentjoint[k][1][i],"1");
          }
        }
      }
    }
    else{
      //Calculate size of each window
      for(k = 0;k < nlayer;k++){
        wlength[k] = 0;
        currentwlength[k] = 0;
      }
      Wfile = fopen(init_file_W, "r");
      ntest = 0;
      while(ntest != -9){
        fscanf(Wfile,"%d %d %d",&ntest,&i,&j);
        if(ntest != -9){
          wlength[ntest] = wlength[ntest] + 1;
          currentwlength[ntest] = currentwlength[ntest] + 1;
        }
      }
      fclose(Wfile);
      Wfile = NULL;

      //Allocate memory
      for(k = 0;k < nlayer;k++){
        W[k] = alloc_matrix_int(wsize*wsize,2);
        currentW[k] = alloc_matrix_int(wsize*wsize,2);
        joint[k] = alloc_matrix_char(2,(1 << wlength[k]),wsize*wsize);
        currentjoint[k] = alloc_matrix_char(2,(1 << wlength[k]),wsize*wsize);
      }

      //Reading W
      Wfile = fopen(init_file_W, "r");
      fscanf(Wfile,"%d", &k);
      while(k != -9){
        for(j = 0;j < wlength[k];j++){
          if(k != -9){
            fscanf(Wfile,"%d", &W[k][j][0]);
            fscanf(Wfile,"%d", &W[k][j][1]);
            currentW[k][j][0] = W[k][j][0];
            currentW[k][j][1] = W[k][j][1];
            fscanf(Wfile,"%d", &k);
          }
        }
      }
      fclose(Wfile);
      Wfile = NULL;

      //Reading joint
      Wfile = fopen(init_file_joint, "r");
      for(k = 0;k < nlayer;k++)
        for(j = 0;j < (1 << wlength[k]);j++){
          fscanf(Wfile,"%d %s %s",&i,joint[k][0][j],joint[k][1][j]);
          strcpy(currentjoint[k][0][j],joint[k][0][j]);
          strcpy(currentjoint[k][1][j],joint[k][1][j]);
        }
      fclose(Wfile);
      Wfile = NULL;
    }

    //for(k = 0;k < nlayer;k++)
      //for(j = 0;j < (1 << wlength[k]);j++)
        //printf("%d %s %s\n",k,joint[k][0][j],joint[k][1][j]);

    //for(k = 0;k < nlayer;k++)
      //for(j = 0;j < wlength[k];j++)
        //printf("%d %d %d\n",k,W[k][j][0],W[k][j][1]);

    //Allocate array with current W operator layers applied to each data point
    Wtrain =  (int****) malloc(train_size*sizeof(int***));
    Wval =  (int****) malloc(val_size*sizeof(int***));
    currentWtrain =  (int****) malloc(train_size*sizeof(int***));
    currentWval =  (int****) malloc(val_size*sizeof(int***));

    //Apply sequence of w operators to each training image
    #pragma omp parallel for num_threads(30) private(i)
    for(k = 0;k < train_size;k++){
      //Alocate array with result of k-th training image
      Wtrain[k] = (int***) malloc(nlayer*sizeof(int**));
      currentWtrain[k] = (int***) malloc(nlayer*sizeof(int**));
      for(i = 0;i < nlayer;i++){
        //For each layer, apply its operator to the output of previous layer
        if(i == 0){
          Wtrain[k][i] = apply_window(train[k],W[i],joint[i],dim - 2*i*(wsize/2),wlength[i],wsize);
          currentWtrain[k][i] = apply_window(train[k],W[i],joint[i],dim - 2*i*(wsize/2),wlength[i],wsize);
        }
        else{
          Wtrain[k][i] = apply_window(Wtrain[k][i-1],W[i],joint[i],dim - 2*i*(wsize/2),wlength[i],wsize);
          currentWtrain[k][i] = apply_window(currentWtrain[k][i-1],W[i],joint[i],dim - 2*i*(wsize/2),wlength[i],wsize);
        }
      }
    }

    //Apply sequence of w operators to each validation image
    #pragma omp parallel for num_threads(30) private(i)
    for(k = 0;k < val_size;k++){
      //Alocate array with result of k-th Validation image
      Wval[k] = (int***) malloc(nlayer*sizeof(int**));
      currentWval[k] = (int***) malloc(nlayer*sizeof(int**));
      for(i = 0;i < nlayer;i++){
        //For each layer, apply its operator to the output of previous layer
        if(i == 0){
          Wval[k][i] = apply_window(val[k],W[i],joint[i],dim - 2*i*(wsize/2),wlength[i],wsize);
          currentWval[k][i] = apply_window(val[k],W[i],joint[i],dim - 2*i*(wsize/2),wlength[i],wsize);
        }
        else{
          Wval[k][i] = apply_window(Wval[k][i-1],W[i],joint[i],dim - 2*i*(wsize/2),wlength[i],wsize);
          currentWval[k][i] = apply_window(currentWval[k][i-1],W[i],joint[i],dim - 2*i*(wsize/2),wlength[i],wsize);
        }
      }
    }

    //Get time start of U-curve algorithm
    gettimeofday(&timecheck, NULL);
    time_Ucurve_start = (long)timecheck.tv_sec * 1000 + (long)timecheck.tv_usec / 1000;

    //Start U-curve algorithm
    printf("Initializing U-curve Algorithm...\n\n");

    //For each digit apply U-curve algorithm
    for(digit = 0;digit < 10;digit++){
      //Start keep to keep searching for windows
      keep = 1;

      //Get time start digit learning
      gettimeofday(&timecheck, NULL);
      time_digit_start = (long)timecheck.tv_sec * 1000 + (long)timecheck.tv_usec / 1000;

      printf("Learning digit %d\n",digit);

      //Error of first window
      error_step = get_error_node(W,Wtrain,Wval,train,ytrain,val,yval,joint,dim,train_size,val_size,wlength,nlayer,wsize,digit,fp);
      printf("\nCurrent Validation Error = %.3f\n",error_step[1]);

      //Save current node
      save_window(W,joint,wlength,nlayer,wsize,digit,number);
      number = number + 1;

      //Find a strong local minimum window
      while(keep){
        //Start error as error of current node
        error_nei[0] = error_step[0];
        error_nei[1] = error_step[1];

        //Check greater neighboors
        check_great_neighboors(W,Wtrain,Wval,currentWtrain,currentWval,currentW,train,ytrain,val,yval,joint,currentjoint,dim,train_size,val_size,wlength,currentwlength,nlayer,wsize,digit,fp,error_nei);

        //Check lesser neighboors
        check_lesser_neighboors(W,Wtrain,Wval,currentWtrain,currentWval,currentW,train,ytrain,val,yval,joint,currentjoint,dim,train_size,val_size,wlength,currentwlength,nlayer,wsize,digit,fp,error_nei);

        printf("\n\nFinished Exhausting Neighboors Current Validation Error = %.3f\n",error_nei[1]);

        //If found a neighboor with lesser error
        if(error_nei[1] < error_step[1]){
          //Free joint
          for(i = 0;i < nlayer;i++)
            free_matrix_char(joint[i],2,(1 << wlength[i]));
          //Store new node
          for(k = 0;k < nlayer;k++){
            //Store window length
            wlength[k] = currentwlength[k];

            //Store new window
            equal_matrix_int(W[k],currentW[k],wlength[k],2);

            //REcreate joint
            joint[k] = alloc_matrix_char(2,(1 << wlength[k]),wsize*wsize);
            strcpy(joint[k][0][0],"0");
            strcpy(joint[k][0][1],"1");
            for(i = 1;i < wlength[k];i++)
              cartesian(joint[k][0],binary,(1 << i),2,wsize*wsize);

            //Copy joint
            for(i = 0;i < (1 << wlength[k]);i++)
              strcpy(joint[k][1][i],currentjoint[k][1][i]);
          }

          //Save new window
          save_window(W,joint,wlength,nlayer,wsize,digit,number);
          number = number + 1;

          //Store new errors
          error_step[0] = error_nei[0];
          error_step[1] = error_nei[1];

          //Copy wtrain
          #pragma omp parallel for num_threads(30) private(i)
          for(k = 0;k < train_size;k++)
            for(i = 0;i < nlayer;i++)
              equal_matrix_int(Wtrain[k][i],currentWtrain[k][i],dim - 2*(i+1)*(wsize/2),dim - 2*(i+1)*(wsize/2));

          //Copy wval
          #pragma omp parallel for num_threads(30) private(i)
          for(k = 0;k < val_size;k++)
            for(i = 0;i < nlayer;i++)
              equal_matrix_int(Wval[k][i],currentWval[k][i],dim - 2*(i+1)*(wsize/2),dim - 2*(i+1)*(wsize/2));
        }
        else
          keep = 0; //Found strong local minimum
        }

      //End time digit learning
      gettimeofday(&timecheck, NULL);
      time_digit_end = (long)timecheck.tv_sec * 1000 + (long)timecheck.tv_usec / 1000;
      printf("Classifier of digit %d was learned in %f hours.\n\n",digit,(float) ((time_digit_end-time_digit_start)/1000)/360.0);
    }

    //Close files of results
    fclose(fp);
    fp = NULL;

    //Free memory
    printf("Freeing memory...\n\n");
    for(i = 0;i < train_size;i++){
      free_matrix_int(train[i],dim);
      train[i] = NULL;
    }
    free(train);

    for(i = 0;i < val_size;i++){
      free_matrix_int(val[i],dim);
      val[i] = NULL;
    }
    free(val);

    for(i = 0;i < nlayer;i++){
      free_matrix_int(W[i],wsize*wsize);
      W[i] = NULL;
    }
    free(W);
    W = NULL;

    free_matrix_int(yval,val_size);
    yval = NULL;

    free_matrix_int(ytrain,train_size);
    ytrain = NULL;

    for(k = 0;k < nlayer;k++)
      free_matrix_char(joint[k],2,(1 << wlength[k]));
    free(joint);
    joint = NULL;

    free(wlength);
    wlength = NULL;

    free(error_nei);
    error_nei = NULL;

    #pragma omp parallel for num_threads(30) private(i)
    for(k = 0;k < val_size;k++){
      for(i = 0;i < nlayer;i++){
        free_matrix_int(Wval[k][i],dim - 2*(i+1)*(wsize/2));
        free_matrix_int(currentWval[k][i],dim - 2*(i+1)*(wsize/2));
      }
      free(Wval[k]);
      free(currentWval[k]);
    }

    #pragma omp parallel for num_threads(30) private(i)
    for(k = 0;k < train_size;k++){
      for(i = 0;i < nlayer;i++){
        free_matrix_int(Wtrain[k][i],dim - 2*(i+1)*(wsize/2));
        free_matrix_int(currentWtrain[k][i],dim - 2*(i+1)*(wsize/2));
      }
      free(Wtrain[k]);
      free(currentWtrain[k]);
    }

    //Organize files saved
    timenow = gmtime(&now);
    strftime (buffer[0],64,"%Y-%m-%d_%H:%M:%S",timenow);
    sprintf(tname,"mkdir ./trace/%s",buffer[0]);
    system(tname);
    sprintf(tname,"mv ./trace/*.txt ./trace/%s/",buffer[0]);
    system(tname);
    free(buffer);

    return 0;
}
