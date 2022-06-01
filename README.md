# A data-driven systematic, consistent and feasible approach to Model Selection
Code of applications in my PhD thesis

**Diego Marcondes. A data-driven systematic, consistent and feasible approach to Model Selection. Universidade de São Paulo, 2022**

The codes refer to three sections in the application section of the thesis:

## 4.1 Learning via the Partition Lattice Learning Space

Simulation studies to illustrate the learning on the Partition Lattice Learning Space. The simulations were performed in **R** with the **partitionUcurve** package, which has its folder in this repository. The folder of this application contains the **R** code of the simulations, the results obtained, and the fst files of the Partition Lattices of up to 10 points in the classifier domain, which are necessary to run optimal U-curve algorithms with the **partitionUcurve** package.

## 4.2 Forecasting variable order Markov chains

Illustration of how learning via Learning Spaces may be applied to forecast a binary sequence by predicting the next value based on the last k under a variable order Markov chain model learned via a subset of the respective Partition Lattice Learning Space. The learning was performed in **R** with the **MarkovLS** package, which has its folder in this repository. The algorithm is applied to forecast the daily variation (positive or negative) of bitcoin in order to develop an investment strategy to it. The folder of this application contains the **R** code used, the data set of bitcoin history and the results of the learned models.

## 4.3 Multilayer W-operator

This application aimed to learn a classifier for digit zero in the MNIST data set via multilayer W-operators. The algorithm was implemented in **C**. This folder of this application contains the **C** code, the **R** code which analysed the results, a folder with each image in the MNIST data set in a txt file, a folder with windows W and operators (joint) at each step of the algorithm, and files with the learned windows W and operators (joint) and the prediction of the learned model on the validation and test samples.


