///////////////////////////////////////////////////////////////////////////////////////////////
//////////A data-driven systematic, consistent and feasible approach to Model Selection////////
//////////4.3 Multilayer W-operator                                                    ////////
//////////PhD Thesis, Diego Marcondes                                                  ////////
//////////Universisty of São Paulo, 2022                                               ////////
///////////////////////////////////////////////////////////////////////////////////////////////

#ifndef WINDOW_UTILS_H_INCLUDED
#define WINDOW_UTILS_H_INCLUDED

//Apply window to matrix
int **apply_window(int **M,int **W,char*** joint,int dim,int wlength,int wsize){
  int **WM,i,j,k,pos;
  char res[100],a[100];

  //Allocate resulting matrix
  WM = alloc_matrix_int(dim - 2*(wsize/2),dim - 2*(wsize/2));

  //For each valida position in image find the input of Woperator and determine the output
  for(i = wsize/2;i < dim - (wsize/2);i++){
    for(j = wsize/2;j < dim - (wsize/2);j++){
      strcpy(res,"");
      pos = 0;
      //Apply window to find input
      for(k = 0;k < wlength;k++){
        if(M[i + W[k][0]][j + W[k][1]] == 0)
          strcat(res,"0");
        else{
          strcat(res,"1");
          //If equal to one refresh position of the point at the joint distribution
          pos = pos + (1 << (wlength-k-1));
        }
      }

      //Test if position is right
      if(strcmp(&joint[0][pos][0],res) != 0)
        printf("We have a problem!\n");

      //Calculate the image of W-operator
      WM[i - wsize/2][j - wsize/2] = joint[1][pos][0] - '0';
    }
  }

  return WM;
}

float get_error_window_preprocessed(int ****Wtrain,int nlayer,int **ytrain,int sample_size,int digit){
  int k;
  float *e,e_final = 0.0;

  //Alocate vector to represent simple loss function
  e = alloc_vector_float(sample_size);

  //For each point in the training sample get error
  #pragma omp parallel for num_threads(30)
  for(k = 0;k < sample_size;k++){
    e[k] = 0;
    if(ytrain[k][0] == digit & Wtrain[k][nlayer-1][0][0] == 0)
      e[k] = 1.0;
    else if(ytrain[k][0] != digit & Wtrain[k][nlayer-1][0][0] == 1)
      e[k] = 1.0;
    //printf("Digit = %d Train = %d Output = %d Erro = %f\n",digit,ytrain[k][0],Wtrain[k][nlayer-1][0][0],e[k]);
  }

  for(k = 0;k < sample_size;k++)
    e_final = e_final + e[k];
  free(e);

  return e_final/sample_size;
}

float get_error_window_changed(int ***W,int ****Wtrain,char ****joint,int *wlength,int nlayer,int wsize,int ***train,int **ytrain,int sample_size,int digit,int layer,int pos,int dim){
  int k,i,j,**M, **M_tmp;;
  float *e,e_final = 0.0;
  char ****joint_tmp;

  //Alloc vector to represent simple loss function
  e = alloc_vector_float(sample_size);

  //Alloc array to represent change joint distribution
  joint_tmp = (char****) malloc(nlayer*sizeof(char***));

  //For each layer copy joint distribution but change given layer and position
  for(i = 0;i < nlayer;i++){
    //Allocate new joint distribution of layer i
    joint_tmp[i] = alloc_matrix_char(2,(1 << wlength[i]),wsize*wsize);
    //For each position, copy domain value and copy image only if not layer pos; otherwise, change value
    for(j = 0;j < (1 << wlength[i]);j++){
      strcpy(joint_tmp[i][0][j],joint[i][0][j]);
      if(i != layer){
        strcpy(joint_tmp[i][1][j],joint[i][1][j]);
      }
      else{
        if(j != pos)
          strcpy(joint_tmp[i][1][j],joint[i][1][j]);
        else{
          if(strcmp(joint[i][1][j],"0") == 0)
            strcpy(joint_tmp[i][1][j],"1");
          else
            strcpy(joint_tmp[i][1][j],"0");
        }
      }
    }
  }

  //For each training point, calculate new w operator and get its erro
  #pragma omp parallel for num_threads(30) private(i,M,M_tmp)
  for(k = 0;k < sample_size;k++){
    //Get preprocessed calculation until layer prior to the changed
    if(layer > 0)
      M = Wtrain[k][layer-1];
    else
      M = train[k];
    //Start recalculation from the changed layer
    for(i = layer;i < nlayer;i++){
      M_tmp = apply_window(M,W[i],joint_tmp[i],dim - 2*i*(wsize/2),wlength[i],wsize);
      if(i > layer)
        free_matrix_int(M,dim - 2*i*(wsize/2));
      M = M_tmp;
    }

    //Calculate error
    e[k] = 0;
    if(ytrain[k][0] == digit & M[0][0] == 0)
      e[k] = 1.0;
    else if(ytrain[k][0] != digit & M[0][0] == 1)
      e[k] = 1.0;
    //printf("Digit = %d Train = %d Output = %d Erro = %f\n",digit,ytrain[k][0],M[0][0],e[k]);
    free_matrix_int(M,1);
  }

  //Sum error
  for(k = 0;k < sample_size;k++)
    e_final = e_final + e[k];
  free(e);

  //Free arrays
  for(k = 0;k < nlayer;k++)
    free_matrix_char(joint_tmp[k],2,(1 << wlength[k]));
  free(joint_tmp);
  joint_tmp = NULL;

  return e_final/sample_size;
}

int *search_hood(int ***W,int ****Wtrain,int ***train,int **ytrain,char ****joint,int dim,int train_size,int *wlength,int nlayer,int wsize,int digit,float w_error){
  int i,k,*res,*point;
  float e_tmp,*error_layer;

  //Allocate vector to store the neighboor and its error
  point = alloc_vector_int(nlayer);
  error_layer = alloc_vector_float(nlayer);

  //For each layer get neighboor with least error
  for(k = 0;k < nlayer;k++){
    //Start with current error
    error_layer[k] = w_error;
    //For each point in the domain of the layer change its value and get error
    for(i = 0;i < (1 << wlength[k]);i++){
      e_tmp = get_error_window_changed(W,Wtrain,joint,wlength,nlayer,wsize,train,ytrain,train_size,digit,k,i,dim);

      //If error lesser, store point
      if(e_tmp < error_layer[k]){
        error_layer[k] = e_tmp;
        point[k] = i;
      }
    }
  }

  //Allocate vector to return result
  res = alloc_vector_int(2);
  res[0] = -1;
  res[1] = -1;
  e_tmp = w_error;

  //Test which layer has lesser error
  for(k = nlayer-1;k >= 0;k--){
    if(error_layer[k] < e_tmp){
      e_tmp = error_layer[k];
      res[0] = k;
      res[1] = point[k];
    }
  }

  //Free vectors
  free(error_layer);
  free(point);

  return res;
}

float *get_error_node(int ***W,int ****Wtrain,int ****Wval,int ***train,int **ytrain,int ***val,int **yval,char ****joint,int dim,int train_size,int val_size,int *wlength,int nlayer,int wsize,int digit,FILE *fp){
  int keep = 1, *nei,i,k;
  float error_val_step = 1.0, error_train_step, *error;

  //Get error of current joint distribution
  error_train_step = get_error_window_preprocessed(Wtrain,nlayer,ytrain,train_size,digit);

  while(keep){
    //Search neighboor with lesser error
    nei = search_hood(W,Wtrain,train,ytrain,joint,dim,train_size,wlength,nlayer,wsize,digit,error_train_step);

    //If there is a neighboor with less error, get its error
    if(nei[0] != -1){
      error_train_step = get_error_window_changed(W,Wtrain,joint,wlength,nlayer,wsize,train,ytrain,train_size,digit,nei[0],nei[1],dim);
      error_val_step = get_error_window_changed(W,Wval,joint,wlength,nlayer,wsize,val,yval,val_size,digit,nei[0],nei[1],dim);

      //Change joint to this neighboor
      if(strcmp(joint[nei[0]][1][nei[1]],"0") == 0)
        strcpy(joint[nei[0]][1][nei[1]],"1");
      else
        strcpy(joint[nei[0]][1][nei[1]],"0");

      //Recalculate preprocessed data to new joint
      #pragma omp parallel for num_threads(30) private(i)
      for(k = 0;k < train_size;k++){
        for(i = nei[0];i < nlayer;i++){
          //Free current preprocessed array
          free_matrix_int(Wtrain[k][i],dim - 2*(i+1)*(wsize/2));
          //Calculate new preprocessed data
          if(i == 0)
            Wtrain[k][i] = apply_window(train[k],W[i],joint[i],dim - 2*i*(wsize/2),wlength[i],wsize);
          else
            Wtrain[k][i] = apply_window(Wtrain[k][i-1],W[i],joint[i],dim - 2*i*(wsize/2),wlength[i],wsize);
        }
      }

      //Recalculate preprocessed data to new joint
      #pragma omp parallel for num_threads(30) private(i)
      for(k = 0;k < val_size;k++){
        for(i = nei[0];i < nlayer;i++){
          //Free current preprocessed array
          free_matrix_int(Wval[k][i],dim - 2*(i+1)*(wsize/2));
          //Calculate new preprocessed data
          if(i == 0)
            Wval[k][i] = apply_window(val[k],W[i],joint[i],dim - 2*i*(wsize/2),wlength[i],wsize);
          else
            Wval[k][i] = apply_window(Wval[k][i-1],W[i],joint[i],dim - 2*i*(wsize/2),wlength[i],wsize);
        }
      }
    }
    else
      keep = 0; //Stop because there is no more neighboor with lesser training sample

    //Free vector
    free(nei);
  }

  //Allocate vector and get error of strong local minimum within node W
  error = alloc_vector_float(2);
  error[0] = error_train_step;
  error[1] = error_val_step;

  return error;
}

void check_great_neighboors(int ***W,int ****Wtrain,int ****Wval,int ****currentWtrain,int ****currentWval,int ***currentW,int ***train,int **ytrain,int ***val,int **yval,char ****joint,char ****currentjoint,int dim,int train_size,int val_size,int *wlength,int *currentwlength,int nlayer,int wsize,int digit,FILE *fp,float *min_error){
  int i,j,k,l,i1,k1,j1,l1,r,stop,c1,c2,**visited,vis,*Nwlength,***NW,progress = 0,max_progress = 0,****NWtrain,****NWval;
  char ****Njoint,rchar,binary[2][2];
  float *error_nei;

  //Allocate newwindow, new preprocessed data, new joint and new window length
  NW = (int***) malloc(nlayer*sizeof(int**));
  Nwlength = alloc_vector_int(nlayer);
  Njoint = (char****) malloc(nlayer*sizeof(char***));
  NWtrain =  (int****) malloc(train_size*sizeof(int***));
  NWval =  (int****) malloc(val_size*sizeof(int***));
  for(k = 0;k < nlayer;k++)
    NW[k] = alloc_matrix_int(wsize*wsize,2);

  //Star binary vector used on cartesian product to generate domain of windows
  strcpy(binary[0],"0");
  strcpy(binary[1],"1");

  //Initialize progress bar
  for(k = 0;k < nlayer;k++)
    max_progress = max_progress + wlength[k];
  max_progress = 4*max_progress;

  //For each layer, search its greater neighboors
  for(k = 0;k < nlayer;k++){
    //Store visited neighboors of the layer
    visited = alloc_matrix_int(4*wlength[k],2);
    for(i = 0;i < 4*wlength[k];i++)
      visited[i][0] = visited[i][1] = 0;
    vis = 0; //Number of visited points

    //For each point in the layer find its neighboors
    for(j = 0;j < wlength[k];j++){
      //Possible nighboors of each point
      for(i = 0;i < 4;i++){
        //Start control if neighboor is valid
        stop = 0;
        //Updtae progress bar
        progress = progress + 1;
        printf("\rExhausting greater neighboors: %d out of %d (%.2f %%) Current Validation Error = %.3f",progress,max_progress,100.0 * ((float) progress)/ ((float) max_progress),min_error[1]);
        fflush(stdout);

        //Choose the neighboor of point
        if(i == 0){
          c1 = W[k][j][0] - 1;
          c2 = W[k][j][1];
        }
        else if(i == 1){
          c1 = W[k][j][0] + 1;
          c2 = W[k][j][1];
        }
        else if(i == 2){
          c1 = W[k][j][0];
          c2 = W[k][j][1] - 1;
        }
        else if(i == 3){
          c1 = W[k][j][0];
          c2 = W[k][j][1] + 1;
        }

        //If outside window, not a valid neighboor
        if(c1 > wsize/2 || c1 < -wsize/2 || c2 > wsize/2 || c2 < -wsize/2)
          continue;

        //Test if new point not already at window
        for(l = 0;l < wlength[k];l++)
          if(W[k][l][0] == c1 && W[k][l][1] == c2)
            stop = 1;

        //Test if new point not already visited
        for(r = 0;r < vis;r++)
          if(visited[r][0] == c1 && visited[r][1] == c2)
            stop = 1;

        //Next iteration if not valid neighboor
        if(stop == 1)
          continue;

        //Store as visisted
        visited[vis][0] = c1;
        visited[vis][1] = c2;
        vis = vis + 1;

        //Copy W and wlength
        for(k1 = 0;k1 < nlayer;k1++){
          Nwlength[k1] = wlength[k1];
          equal_matrix_int(NW[k1],W[k1],wlength[k1],2);
        }

        //Add new point to NW
        Nwlength[k] = Nwlength[k] + 1;
        NW[k][Nwlength[k]-1][0] = c1;
        NW[k][Nwlength[k]-1][1] = c2;

        //Create Njoint
        for(k1 = 0;k1 < nlayer;k1++){
          //Create domain
          Njoint[k1] = alloc_matrix_char(2,(1 << Nwlength[k1]),wsize*wsize);
          strcpy(Njoint[k1][0][0],"0");
          strcpy(Njoint[k1][0][1],"1");
          for(i1 = 1;i1 < Nwlength[k1];i1++)
            cartesian(Njoint[k1][0],binary,(1 << i1),2,wsize*wsize);

          //Copy layers which have not changed and sample image of new layer
          for(i1 = 0;i1 < (1 << Nwlength[k1]);i1++){
            if(k1 != k)
              strcpy(Njoint[k1][1][i1],joint[k1][1][i1]);
            else{
              rchar = rand() % 2;
              if(rchar == 0)
                strcpy(Njoint[k1][1][i1],"0");
              else
                strcpy(Njoint[k1][1][i1],"1");
            }
          }
        }

        //Apply new W operator to train sample
        #pragma omp parallel for num_threads(30) private(i1)
        for(k1 = 0;k1 < train_size;k1++){
          NWtrain[k1] = (int***) malloc(nlayer*sizeof(int**));
          for(i1 = 0;i1 < nlayer;i1++){
            //If layer greater than changed, apply window
            if(i1 >= k){
              if(i1 == 0)
                NWtrain[k1][i1] = apply_window(train[k1],NW[i1],Njoint[i1],dim - 2*i1*(wsize/2),Nwlength[i1],wsize);
              else
                NWtrain[k1][i1] = apply_window(NWtrain[k1][i1-1],NW[i1],Njoint[i1],dim - 2*i1*(wsize/2),Nwlength[i1],wsize);
            }
            //If layer lesser than changed, copy wtrain
            else{
              NWtrain[k1][i1] = alloc_matrix_int(dim - 2*(i1+1)*(wsize/2),dim - 2*(i1+1)*(wsize/2));
              equal_matrix_int(NWtrain[k1][i1],Wtrain[k1][i1],dim - 2*(i1+1)*(wsize/2),dim - 2*(i1+1)*(wsize/2));
            }
          }
        }

        //Apply new W operator to val sample
        #pragma omp parallel for num_threads(30) private(i1)
        for(k1 = 0;k1 < val_size;k1++){
          NWval[k1] = (int***) malloc(nlayer*sizeof(int**));
          for(i1 = 0;i1 < nlayer;i1++){
            //If layer greater than changed, apply window
            if(i1 >= k){
              if(i1 == 0)
                NWval[k1][i1] = apply_window(val[k1],NW[i1],Njoint[i1],dim - 2*i1*(wsize/2),Nwlength[i1],wsize);
              else
                NWval[k1][i1] = apply_window(NWval[k1][i1-1],NW[i1],Njoint[i1],dim - 2*i1*(wsize/2),Nwlength[i1],wsize);
            }
            //If layer lesser than changed, copy wtrain
            else{
              NWval[k1][i1] = alloc_matrix_int(dim - 2*(i1+1)*(wsize/2),dim - 2*(i1+1)*(wsize/2));
              equal_matrix_int(NWval[k1][i1],Wval[k1][i1],dim - 2*(i1+1)*(wsize/2),dim - 2*(i1+1)*(wsize/2));
            }
          }
        }

        //Get error of neighboor
        error_nei = get_error_node(NW,NWtrain,NWval,train,ytrain,val,yval,Njoint,dim,train_size,val_size,Nwlength,nlayer,wsize,digit,fp);

        //If error lesser, store new minimum
        if(error_nei[1] < min_error[1]){
          //Store its error
          min_error[1] = error_nei[1];
          min_error[0] = error_nei[0];

          //Free currentjoint
          for(i1 = 0;i1 < nlayer;i1++)
            free_matrix_char(currentjoint[i1],2,(1 << currentwlength[i1]));
          //Update currentwlength
          for(i1 = 0;i1 < nlayer;i1++)
              currentwlength[i1] = Nwlength[i1];

          //For each layer, store NW and Njoint as current
          for(k1 = 0;k1 < nlayer;k1++){
            //Equal currentW to NW
            equal_matrix_int(currentW[k1],NW[k1],Nwlength[k1],2);

            //Recreate currentjoint
            currentjoint[k1] = alloc_matrix_char(2,(1 << currentwlength[k1]),wsize*wsize);
            strcpy(currentjoint[k1][0][0],"0");
            strcpy(currentjoint[k1][0][1],"1");
            for(i1 = 1;i1 < currentwlength[k1];i1++)
              cartesian(currentjoint[k1][0],binary,(1 << i1),2,wsize*wsize);

            //Copy Njoint to currentjoint
            for(i1 = 0;i1 < (1 << currentwlength[k1]);i1++)
              strcpy(currentjoint[k1][1][i1],Njoint[k1][1][i1]);
          }

          //Copy NWtrain to currentWtrain
          #pragma omp parallel for num_threads(30) private(i1)
          for(k1 = 0;k1 < train_size;k1++)
            for(i1 = 0;i1 < nlayer;i1++)
              equal_matrix_int(currentWtrain[k1][i1],NWtrain[k1][i1],dim - 2*(i1+1)*(wsize/2),dim - 2*(i1+1)*(wsize/2));

          //Copy NWval to currentWval
          #pragma omp parallel for num_threads(30) private(i1)
          for(k1 = 0;k1 < val_size;k1++)
            for(i1 = 0;i1 < nlayer;i1++)
              equal_matrix_int(currentWval[k1][i1],NWval[k1][i1],dim - 2*(i1+1)*(wsize/2),dim - 2*(i1+1)*(wsize/2));

          //free Njoint, Nwtrain and NWval
          for(i1 = 0;i1 < nlayer;i1++)
            free_matrix_char(Njoint[i1],2,(1 << Nwlength[i1]));

          #pragma omp parallel for num_threads(30) private(i1)
          for(k1 = 0;k1 < val_size;k1++){
            for(i1 = 0;i1 < nlayer;i1++)
              free_matrix_int(NWval[k1][i1],dim - 2*(i1+1)*(wsize/2));
            free(NWval[k1]);
          }

          #pragma omp parallel for num_threads(30) private(i1)
          for(k1 = 0;k1 < train_size;k1++){
            for(i1 = 0;i1 < nlayer;i1++)
              free_matrix_int(NWtrain[k1][i1],dim - 2*(i1+1)*(wsize/2));
            free(NWtrain[k1]);
          }
        }
        else{
          //Free Njoint, Nwtrain and NWval
          for(i1 = 0;i1 < nlayer;i1++)
            free_matrix_char(Njoint[i1],2,(1 << Nwlength[i1]));

          #pragma omp parallel for num_threads(30) private(i1)
          for(k1 = 0;k1 < val_size;k1++){
            for(i1 = 0;i1 < nlayer;i1++)
              free_matrix_int(NWval[k1][i1],dim - 2*(i1+1)*(wsize/2));
            free(NWval[k1]);
          }

          #pragma omp parallel for num_threads(30) private(i1)
          for(k1 = 0;k1 < train_size;k1++){
            for(i1 = 0;i1 < nlayer;i1++)
              free_matrix_int(NWtrain[k1][i1],dim - 2*(i1+1)*(wsize/2));
            free(NWtrain[k1]);
          }
        }
      }
    }

    //Free visited
    free_matrix_int(visited,4*wlength[k]);
    visited = NULL;
  }

  //free
  free(Njoint);
  Njoint = NULL;
  free(NWval);
  NWval = NULL;
  free(NWtrain);
  NWtrain == NULL;
  for(k = 0;k < nlayer;k++)
    free_matrix_int(NW[k],wsize*wsize);
  free(NW);
  NW = NULL;
  free(Nwlength);
  Nwlength = NULL;
}

void check_lesser_neighboors(int ***W,int ****Wtrain,int ****Wval,int ****currentWtrain,int ****currentWval,int ***currentW,int ***train,int **ytrain,int ***val,int **yval,char ****joint,char ****currentjoint,int dim,int train_size,int val_size,int *wlength,int *currentwlength,int nlayer,int wsize,int digit,FILE *fp,float* min_error){
  int j,k,k1,j1,i1,*Nwlength,***NW,****NWtrain,****NWval,progress = 0,max_progress = 0;
  char ****Njoint,rchar,binary[2][2];
  float *error_nei;

  //Alloc new window, new joint, new wlength, new Wtrain e new Wval
  NW = (int***) malloc(nlayer*sizeof(int**));
  Njoint = (char****) malloc(nlayer*sizeof(char***));
  Nwlength = alloc_vector_int(nlayer);
  NWtrain =  (int****) malloc(train_size*sizeof(int***));
  NWval =  (int****) malloc(val_size*sizeof(int***));
  for(k1 = 0;k1 < nlayer;k1++)
    NW[k1] = alloc_matrix_int(wsize*wsize,2);

  //Initialize vector for cartesian product to generate joint
  strcpy(binary[0],"0");
  strcpy(binary[1],"1");

  //Initialize progress bar
  for(k = 0;k < nlayer;k++)
    max_progress = max_progress + wlength[k];
  max_progress = max_progress;

  //For each layer
  for(k = 0;k < nlayer;k++){
    //Erase each point of the layer
    for(j = 0;j < wlength[k];j++){
      //Update progress bar
      progress = progress + 1;
      printf("\rExhausting lesser neighboors: %d out of %d (%.2f %%) Current Validation Error = %.3f",progress,max_progress,100.0 * ((float) progress)/ ((float) max_progress),min_error[1]);
      fflush(stdout);

      //Copy W and wlength
      for(k1 = 0;k1 < nlayer;k1++){
        //If other layer, copy it
        if(k1 != k){
          Nwlength[k1] = wlength[k1];
          equal_matrix_int(NW[k1],W[k1],Nwlength[k1],2);
        }
        else{
          //If changed layer, copy all but dropped point
          Nwlength[k1] = wlength[k1] - 1;
          if(j > 0)
            for(j1 = 0;j1 < j;j1++){
              NW[k1][j1][0] = W[k1][j1][0];
              NW[k1][j1][1] = W[k1][j1][1];
            }
          for(j1 = j;j1 < Nwlength[k1];j1++){
            NW[k1][j1][0] = W[k1][j1+1][0];
            NW[k1][j1][1] = W[k1][j1+1][1];
          }
        }
      }

      //Create Njoint
      for(k1 = 0;k1 < nlayer;k1++){
        //Create domain
        Njoint[k1] = alloc_matrix_char(2,(1 << Nwlength[k1]),wsize*wsize);
        strcpy(Njoint[k1][0][0],"0");
        strcpy(Njoint[k1][0][1],"1");
        for(i1 = 1;i1 < Nwlength[k1];i1++)
          cartesian(Njoint[k1][0],binary,(1 << i1),2,wsize*wsize);

        //Copy layers which have not changed and sample image of new layer
        for(i1 = 0;i1 < (1 << Nwlength[k1]);i1++){
          if(k1 != k)
            strcpy(Njoint[k1][1][i1],joint[k1][1][i1]);
          else{
            rchar = rand() % 2;
            if(rchar == 0)
              strcpy(Njoint[k1][1][i1],"0");
            else
              strcpy(Njoint[k1][1][i1],"1");
          }
        }
      }

      //Apply neighboor to train sample
      #pragma omp parallel for num_threads(30) private(i1)
      for(k1 = 0;k1 < train_size;k1++){
        //Alloc array
        NWtrain[k1] = (int***) malloc(nlayer*sizeof(int**));
        for(i1 = 0;i1 < nlayer;i1++){
          //If layer greater than changed, apply windows
          if(i1 >= k){
            if(i1 == 0)
              NWtrain[k1][i1] = apply_window(train[k1],NW[i1],Njoint[i1],dim - 2*i1*(wsize/2),Nwlength[i1],wsize);
            else
              NWtrain[k1][i1] = apply_window(NWtrain[k1][i1-1],NW[i1],Njoint[i1],dim - 2*i1*(wsize/2),Nwlength[i1],wsize);
          }
          else{
            //If layer lesser than changed, just copy Wtrain
            NWtrain[k1][i1] = alloc_matrix_int(dim - 2*(i1+1)*(wsize/2),dim - 2*(i1+1)*(wsize/2));
            equal_matrix_int(NWtrain[k1][i1],Wtrain[k1][i1],dim - 2*(i1+1)*(wsize/2),dim - 2*(i1+1)*(wsize/2));
          }
        }
      }

      //Apply neighboor to val sample
      #pragma omp parallel for num_threads(30) private(i1)
      for(k1 = 0;k1 < val_size;k1++){
        //Alloc array
        NWval[k1] = (int***) malloc(nlayer*sizeof(int**));
        for(i1 = 0;i1 < nlayer;i1++){
          //If layer greater than changed, apply windows
          if(i1 >= k){
            if(i1 == 0)
              NWval[k1][i1] = apply_window(val[k1],NW[i1],Njoint[i1],dim - 2*i1*(wsize/2),Nwlength[i1],wsize);
            else
              NWval[k1][i1] = apply_window(NWval[k1][i1-1],NW[i1],Njoint[i1],dim - 2*i1*(wsize/2),Nwlength[i1],wsize);
          }
          else{
            //If layer lesser than changed, just copy Wtrain
            NWval[k1][i1] = alloc_matrix_int(dim - 2*(i1+1)*(wsize/2),dim - 2*(i1+1)*(wsize/2));
            equal_matrix_int(NWval[k1][i1],Wval[k1][i1],dim - 2*(i1+1)*(wsize/2),dim - 2*(i1+1)*(wsize/2));
          }
        }
      }

      //Search neihboor
      error_nei = get_error_node(NW,NWtrain,NWval,train,ytrain,val,yval,Njoint,dim,train_size,val_size,Nwlength,nlayer,wsize,digit,fp);

      //If error lesser, store new minimum
      if(error_nei[1] < min_error[1]){
        //Store error of new minimum
        min_error[1] = error_nei[1];
        min_error[0] = error_nei[0];
        printf("\nCurrent Validation Error = %f\n",min_error[1]);

        //Free currentjoint
        for(i1 = 0;i1 < nlayer;i1++)
          free_matrix_char(currentjoint[i1],2,(1 << currentwlength[i1]));

        //Update currentwlength
        for(i1 = 0;i1 < nlayer;i1++)
            currentwlength[i1] = Nwlength[i1];

        //Store new W, new wlength, new Wtrain and new Wval as current
        for(k1 = 0;k1 < nlayer;k1++){
          //Store NW
          equal_matrix_int(currentW[k1],NW[k1],Nwlength[k1],2);

          //Recreate joint domain
          currentjoint[k1] = alloc_matrix_char(2,(1 << currentwlength[k1]),wsize*wsize);
          strcpy(currentjoint[k1][0][0],"0");
          strcpy(currentjoint[k1][0][1],"1");
          for(i1 = 1;i1 < currentwlength[k1];i1++)
            cartesian(currentjoint[k1][0],binary,(1 << i1),2,wsize*wsize);

          //Copy new image
          for(i1 = 0;i1 < (1 << currentwlength[k1]);i1++)
            strcpy(currentjoint[k1][1][i1],Njoint[k1][1][i1]);
        }

        //Copy new wtrain to current
        #pragma omp parallel for num_threads(30) private(i1)
        for(k1 = 0;k1 < train_size;k1++)
          for(i1 = 0;i1 < nlayer;i1++)
            equal_matrix_int(currentWtrain[k1][i1],NWtrain[k1][i1],dim - 2*(i1+1)*(wsize/2),dim - 2*(i1+1)*(wsize/2));

        //Copy new wval to current
        #pragma omp parallel for num_threads(30) private(i1)
        for(k1 = 0;k1 < val_size;k1++)
          for(i1 = 0;i1 < nlayer;i1++)
            equal_matrix_int(currentWval[k1][i1],NWval[k1][i1],dim - 2*(i1+1)*(wsize/2),dim - 2*(i1+1)*(wsize/2));

        //Free new objects
        for(i1 = 0;i1 < nlayer;i1++)
          free_matrix_char(Njoint[i1],2,(1 << Nwlength[i1]));

        #pragma omp parallel for num_threads(30) private(i1)
        for(k1 = 0;k1 < val_size;k1++){
          for(i1 = 0;i1 < nlayer;i1++)
            free_matrix_int(NWval[k1][i1],dim - 2*(i1+1)*(wsize/2));
          free(NWval[k1]);
        }

        #pragma omp parallel for num_threads(30) private(i1)
        for(k1 = 0;k1 < train_size;k1++){
          for(i1 = 0;i1 < nlayer;i1++)
            free_matrix_int(NWtrain[k1][i1],dim - 2*(i1+1)*(wsize/2));
          free(NWtrain[k1]);
        }
      }
      else{
        //Free new objects
        for(i1 = 0;i1 < nlayer;i1++)
          free_matrix_char(Njoint[i1],2,(1 << Nwlength[i1]));

        #pragma omp parallel for num_threads(30) private(i1)
        for(k1 = 0;k1 < val_size;k1++){
          for(i1 = 0;i1 < nlayer;i1++)
            free_matrix_int(NWval[k1][i1],dim - 2*(i1+1)*(wsize/2));
          free(NWval[k1]);
        }

        #pragma omp parallel for num_threads(30) private(i1)
        for(k1 = 0;k1 < train_size;k1++){
          for(i1 = 0;i1 < nlayer;i1++)
            free_matrix_int(NWtrain[k1][i1],dim - 2*(i1+1)*(wsize/2));
          free(NWtrain[k1]);
        }
      }
    }
  }

  //free
  free(Njoint);
  Njoint = NULL;
  for(k = 0;k < nlayer;k++)
    free_matrix_int(NW[k],wsize*wsize);
  free(NW);
  NW = NULL;
  free(Nwlength);
  Nwlength = NULL;
}

void save_window(int ***W,char ****joint,int *wlength,int nlayer,int wsize,int digit,int number){
  char namew[100],namej[100];
  int k,j;
  FILE *fw,*fj;

  //Name of files joint and window
  strcpy(namew,"");
  strcpy(namej,"");
  sprintf(namew,"./trace/W_%d_%05d.txt",digit,number);
  sprintf(namej,"./trace/Joint_%d_%05d.txt",digit,number);

  //Open files
  fw = fopen(namew, "w");
  fj = fopen(namej, "w");

  //Print to window file
  for(k = 0;k < nlayer;k++)
    for(j = 0;j < wlength[k];j++)
      fprintf(fw,"%d %d %d\n",k,W[k][j][0],W[k][j][1]);
  fclose(fw);

  //Print to joint file
  for(k = 0;k < nlayer;k++)
    for(j = 0;j < (1 << wlength[k]);j++)
      fprintf(fj,"%d %s %s\n",k,joint[k][0][j],joint[k][1][j]);
  fclose(fj);
}

#endif // WINDOW_UTILS_H_INCLUDED
