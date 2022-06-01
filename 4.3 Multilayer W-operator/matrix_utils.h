///////////////////////////////////////////////////////////////////////////////////////////////
//////////A data-driven systematic, consistent and feasible approach to Model Selection////////
//////////4.3 Multilayer W-operator                                                    ////////
//////////PhD Thesis, Diego Marcondes                                                  ////////
//////////Universisty of São Paulo, 2022                                               ////////
///////////////////////////////////////////////////////////////////////////////////////////////

#ifndef MATRIX_UTILS_H_INCLUDED
#define MATRIX_UTILS_H_INCLUDED

//Allocate matrix float
float **alloc_matrix_float(int row,int col){
    //Allocate a matrix with 'row' rows and 'col' collumns
    float **M;
    int i;

    M = (float **) malloc(row * sizeof(float *));

    for(i = 0;i < row;i++)
        M[i] = (float *) malloc(col * sizeof(float));

    return M;
}

//Allocate matrix int
int **alloc_matrix_int(int row,int col){
    //Allocate a matrix with 'row' rows and 'col' collumns
    int **M;
    int i;

    M = (int **) malloc(row * sizeof(int *));

    for(i = 0;i < row;i++)
        M[i] = (int *) malloc(col * sizeof(int));

    return M;
}

//Allocate matrix char
char ***alloc_matrix_char(int row,int col,int charsize){
    //Allocate a matrix with 'row' rows and 'col' collumns
    char ***M;
    int i,j;

    M = malloc(row * sizeof(char **));

    for(i = 0;i < row;i++){
      M[i] = malloc(col * sizeof(char*));
      for(j = 0;j < col;j++)
        M[i][j] = malloc((charsize+1) * sizeof(char));
    }

    return M;
}

//Free matrix float
void free_matrix_float(float **M,int row){
    //Free the matrix 'M', with 'row' rows
    int i;

    for(i = 0;i < row;i++)
        free(M[i]);

    free(M);
}

//Free matrix int
void free_matrix_int(int **M,int row){
    //Free the matrix 'M', with 'row' rows
    int i;

    for(i = 0;i < row;i++)
        free(M[i]);

    free(M);
}

//Free matrix char
void free_matrix_char(char ***M,int row,int col){
    //Free the matrix 'M', with 'row' rows and 'col' collums
    int i,j;

    for(i = 0;i < row;i++){
      for(j = 0;j < col;j++){
        free(M[i][j]);
      }
      free(M[i]);
    }

    free(M);
}

//Allocate vector float
float *alloc_vector_float(int col){
    //Allocate a vector with 'col' coordinates
    float *M;

    M = (float *) malloc(col * sizeof(float));

    return M;
}

//Allocate vector int
int *alloc_vector_int(int col){
    //Allocate a vector with 'col' coordinates
    int *M;

    M = (int *) malloc(col * sizeof(int));

    return M;
}

//Allocate vector int
char **alloc_vector_char(int col,int charsize){
    //Allocate a vector with 'col' coordinates
    char **M;
    int i;

    M = (char **) malloc(col * sizeof(char*));

    for(i = 0;i < col;i++)
      M[i] = (char*) malloc(charsize * sizeof(char));

    return M;
}


//Print matrix
void print_matrix(float **M,int row,int col){
    //Print the matrix 'M', with 'row' rows and 'col' collumns
    int i,j;

    for(i = 0;i < row;i++){
        for(j = 0;j < col;j++)
            printf("%0.lf ",M[i][j]);
        printf("\n");
    }
}

//Print matrix int
void print_matrix_int(int **M,int row,int col){
    //Print the matrix 'M', with 'row' rows and 'col' collumns
    int i,j;

    for(i = 0;i < row;i++){
        for(j = 0;j < col;j++)
            printf("%d ",M[i][j]);
        printf("\n");
    }
}

//Print matrix char
void print_matrix_char(char ***M,int row,int col){
    //Print the matrix 'M', with 'row' rows and 'col' collumns
    int i,j;

    for(i = 0;i < row;i++){
        for(j = 0;j < col;j++)
            printf("%s ",&M[i][j][0]);
        printf("\n");
    }
}

//Delete collumn of matrix
void delete_column(float **M,int row,int col,int del){
    int i,j;

    if(del < col - 1){
        for(j = del + 1;j < col;j++)
            for(i = 0;i < row;i++)
                M[i][j-1] = M[i][j];
    }
    for(i = 0;i < row;i++)
        M[i] = realloc(M[i], sizeof(float)*(col-1));
}

//Print vector
void print_vector(float *M,int col){
    //Print the vector 'M' with 'col' coordinates
    int i;

    for(i = 0;i < col;i++)
        printf("%.6lf ",M[i]);
    printf("\n");

}

//Print vector int
void print_vector_int(int *M,int col){
    //Print the vector 'M' with 'col' coordinates
    int i;

    for(i = 0;i < col;i++)
        printf("%d ",M[i]);
    printf("\n");

}

//Print vector int
void print_vector_char(char **M,int col){
    //Print the vector 'M' with 'col' coordinates
    int i;

    for(i = 0;i < col;i++)
        printf("%s ",&M[0][i]);
    printf("\n");

}

//Print matrix on file
void fprint_matrix(float **M,int row,int col,FILE *fp){
    //Print the matrix 'M', with 'row' rows and 'col' collumns, on the file 'fp'
    int i,j;

    for(i = 0;i < row;i++){
        for(j = 0;j < col-1;j++)
            fprintf(fp,"%.16lf;",M[i][j]); //The values are separated by ';' (a .csv file format)
        fprintf(fp,"%.16lf\n",M[i][col-1]);
    }
}

//Make matrices equal
void equal_matrix(float **M1,float **M2,int row,int col){
    //Make matrix 'M1' equal to matrix 'M2', both having 'row' rows and 'col' collumns
    int i, j;

    for(i = 0;i < row;i++)
        for(j = 0;j < col;j++)
            M1[i][j] = M2[i][j];
}

//Make matrices equal int
void equal_matrix_int(int **M1,int **M2,int row,int col){
    //Make matrix 'M1' equal to matrix 'M2', both having 'row' rows and 'col' collumns
    int i, j;

    for(i = 0;i < row;i++)
        for(j = 0;j < col;j++)
          M1[i][j] = M2[i][j];
}

//Make vectors equal
void equal_vector(float *M1,float *M2,int col){
    //Make vector 'M1' equal to vector 'M2', both having 'col' coordinates
    int j;

    for(j = 0;j < col;j++)
        M1[j] = M2[j];
}

//Sum matrices
float **sum_matrix(float **M1,float **M2,int row,int col){
    //Sum matrices 'M1' and 'M1', which have 'row' rows and 'col' collumns
    float **s;
    int i,j;

    s = alloc_matrix_float(row,col);

    for(i = 0;i < row;i++)
        for(j = 0;j < col;j++)
            s[i][j] = M1[i][j] + M2[i][j];

    return s;
}

//Subtract matrices
float **subtract_matrix(float **M1,float **M2,int row,int col){
    //Subtract matrices 'M1' and 'M1', which have 'row' rows and 'col' collumns
    float **s;
    int i,j;

    s = alloc_matrix_float(row,col);

    for(i = 0;i < row;i++)
        for(j = 0;j < col;j++)
            s[i][j] = M1[i][j] - M2[i][j];

    return s;
}

//Sum vectors
float *sum_vector(float *M1,float *M2,int col){
    //Sum vectors 'M1' and 'M2', which have 'col' coordinates
    float *s;
    int i;

    s = alloc_vector_float(col);

    for(i = 0;i < col;i++)
        s[i] = M1[i] + M2[i];

    return s;
}

//Subtract vectors
float *subtract_vector(float *M1,float *M2,int col){
    //Subtract vectors 'M1' and 'M2', which have 'col' coordinates
    float *s;
    int i;

    s = alloc_vector_float(col);

    for(i = 0;i < col;i++)
        s[i] = M1[i] - M2[i];

    return s;
}

//Multiply matrices
float **multiply_matrix(float **M1,float **M2,int row1,int col1,int col2){
    //Multiply matrix 'M2', with 'col1' rows and 'col2' collumns, by matrix 'M1', with 'row1' rows and 'col1' collumns
    float **M;
    int i,j,k;

    M = alloc_matrix_float(row1,col2);

    for(i = 0;i < row1;i++)
        for(j = 0;j < col2;j++){
            M[i][j] = 0;
            for(k = 0;k < col1;k++)
                M[i][j] = M[i][j] + M1[i][k] * M2[k][j];
        }
    return M;
}

//Transpose matrix
float **transpose_matrix(float **M,int row,int col){
    //Transpose matrix 'M', with 'row' rows and 'col' collumns
    float **t;
    int i,j;

    t = alloc_matrix_float(col,row);

    for(i = 0;i < col;i++)
        for(j = 0;j < row;j++)
            t[i][j] = M[j][i];

    return t;
}

//Multiply a matrix by a scalar
float **smultiply_matrix(float **M,int row,int col,float a){
    //Multiply matrix 'M', with 'row' rows and 'col' collumns, by the scalar 'a'
    int i,j;
    float **M1;

    M1 = alloc_matrix_float(row,col);

    for(i = 0;i < row;i++)
        for(j = 0;j < col;j++)
            M1[i][j] = a * M[i][j];

    return M1;
}

//Multiply a vector by a scalar
float *smultiply_vector(float *M,int col,float a){
    //Multiply vector 'M', with 'col' coordinates, by the scalar 'a'
    int i;
    float *M1;

    M1 = alloc_vector_float(col);

    for(i = 0;i < col;i++)
        M1[i] = a * M[i];

    return M1;
}

//Sum a constant to a matrix
float **csum_matrix(float **M,int row,int col,float a){
    //Sum the constant 'a' to matrix 'M', with 'row' rows and 'col' collumns
    int i,j;
    float **M1;

    M1 = alloc_matrix_float(row,col);

    for(i = 0;i < row;i++)
        for(j = 0;j < col;j++)
            M1[i][j] = M[i][j] + a;

    return M1;
}

//Sum a constant to a vector
float *csum_vector(float *M,int col,float a){
    //Sum constant 'a' to vector 'M', with 'col' coordinates
    int i;
    float *M1;

    M1 = alloc_vector_float(col);

    for(i = 0;i < col;i++)
        M1[i] = M[i] + a;

    return M1;
}

//Read matrix from file
int readmatrix(size_t rows, size_t cols, int** a, const char* filename,int increase){

    FILE *pf;
    pf = fopen(filename, "r");

    if(pf == NULL)
      return 0;

    for(size_t i = 0; i < rows; ++i){
      for(size_t j = 0; j < cols; ++j){
        if(i >= increase & j >= increase)
          fscanf(pf, "%d", &a[i][j]);
        else
          a[i][j] = 0;
      }
    }


    fclose (pf);
    return 1;
}

#endif // MATRIX_UTILS_H_INCLUDED
