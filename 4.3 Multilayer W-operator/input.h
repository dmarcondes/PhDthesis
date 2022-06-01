///////////////////////////////////////////////////////////////////////////////////////////////
//////////A data-driven systematic, consistent and feasible approach to Model Selection////////
//////////4.3 Multilayer W-operator                                                    ////////
//////////PhD Thesis, Diego Marcondes                                                  ////////
//////////Universisty of São Paulo, 2022                                               ////////
///////////////////////////////////////////////////////////////////////////////////////////////

//ifndef INPUT_H_INCLUDED
//define INPUT_H_INCLUDED
////include <math.h>

//Size of window
int wsize = 5;

//File name to save results
char file_name[100] = "Wresults.txt";

//Size of training and validation sample
int train_size = 45000;
int val_size = 15000;

//Dimension of images
int dim = 28;

//Random initialize (files should end with -9)
int random_init = 1;
char init_file_W[100] = ".txt";
char init_file_joint[100] = ".txt";

//Number of next node visited
int number = 16;

//Subdirectory of files
char fdir[8] = "";

//endif // INPUT_H_INCLUDED
