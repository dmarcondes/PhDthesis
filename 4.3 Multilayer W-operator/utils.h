///////////////////////////////////////////////////////////////////////////////////////////////
//////////A data-driven systematic, consistent and feasible approach to Model Selection////////
//////////4.3 Multilayer W-operator                                                    ////////
//////////PhD Thesis, Diego Marcondes                                                  ////////
//////////Universisty of São Paulo, 2022                                               ////////
///////////////////////////////////////////////////////////////////////////////////////////////

#ifndef UTILS_H_INCLUDED
#define UTILS_H_INCLUDED
#include "matrix_utils.h"
#include <math.h>
#define pi 3.14159265358979323846

//Sampling from uniform distribution
int sample_discrete(int n,float *prob){
    double sample,s = 0;
    int i,u;

    sample = ((double) random())/((double) RAND_MAX);
    for(i = 0;i < n;i++){
        s = s + prob[i];
        if(s >= sample)
            return i;
    }
}

//Minimum between two numbers
double min(double x, double y){
    return (x <= y) * x + (x > y) * y;
}

//Maximum between two int
double max_int(int x, int y){
    return (x <= y) * y + (x > y) * x;
}

//Factorial of a number
int factorial(int n){
    if (n == 0)
        return 1;
    else
        return n * factorial(n - 1);
}

//Order the indexes of a vector by its values
int *order(double *dist,int n){
    //Order the indexes of a vector in ascending order of its values
    int i, j, a, *index;

    index = (int *) malloc(n * sizeof(int));

    for(i = 0;i < n;i++)
        index[i] = i;//Initialize the index vector

    for(i=0;i < n;i++)
        for(j = i+1;j < n;j++)
            if(dist[index[i]] > dist[index[j]]){//If the value at the index in position i is greater than the index at position j (i < j)
                a = index[i]; //Swap the index at position i and the index at position j
                index[i] = index[j];
                index[j] = a;
            }

    return index;
}

//Order in ascending order a vector of m integers
void order_int(int *index,int m){
    //Order the a vector of integers in ascending order of its values
    int i, j, a;

    for(i=0;i < m;i++)
        for(j = i+1;j < m;j++)
            if(index[i] > index[j]){//If the value in position i is greater than the value at position j (i < j)
                a = index[i]; //Swap the value at position i and the value at position j
                index[i] = index[j];
                index[j] = a;
            }
}

int ipow(int base, int exp)
{
    int result = 1;
    for (;;)
    {
        if (exp & 1)
            result *= base;
        exp >>= 1;
        if (!exp)
            break;
        base *= base;
    }

    return result;
}

//Cartesian product between string vectors
void cartesian(char** vector,char s[2][2],int len_vector,int len_s,int charsize){
    char** cpvec, tmp[100];
    int i,j,k = 0;

    //Alloc vector to copy vector
    cpvec = alloc_vector_char(len_vector,charsize);

    //Copy vector
    for(i = 0;i < len_vector;i++)
      strcpy(cpvec[i],vector[i]);

    //Cartesian product between vector and s saved on vector
    for(i = 0;i < len_vector;i++)
      for(j = 0;j < len_s;j++){
        sprintf(tmp,"%s%s",cpvec[i],s[j]);
        strcpy(vector[k],tmp);
        k = k + 1;
      }

    free(cpvec);
    cpvec = NULL;
}

#endif // UTILS_H_INCLUDED
