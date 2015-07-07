##################################################################################################################
## Autoencoder.R
## Autoencoder is an unsupervised learning appraoch to get the feature of the data. We consider it as a non-linear
## PCA. Our goal is to do clustering based on deep learning.
## Author: Meilei
##################################################################################################################
library(dplyr)
library(ggplot2)
library(autoencoder)


# load the data in the package --------------------------------------------
data('training_matrix_N=5e3_Ninput=100')
## look at a smaller training set
trainM = training.matrix[1:500,]


# set up the autoencoder architecture  ------------------------------------

nl=3                          ## number of layers (default is 3: input, hidden, output)
unit.type = "logistic"        ## specify the network unit type, i.e., the unit's 
## activation function ("logistic" or "tanh")
Nx.patch=10                   ## width of training image patches, in pixels
Ny.patch=10                   ## height of training image patches, in pixels
N.input = Nx.patch*Ny.patch   ## number of units (neurons) in the input layer (one unit per pixel)
N.hidden = 10*10                ## number of units in the hidden layer
lambda = 0.0002               ## weight decay parameter     
beta = 6                      ## weight of sparsity penalty term 
rho = 0.01                    ## desired sparsity parameter
epsilon <- 0.001              ## a small parameter for initialization of weights 
## as small gaussian random numbers sampled from N(0,epsilon^2)
max.iterations = 2000         ## number of iterations in optimizer


# Train the auto encoder on training.matrix -------------------------------

autoencoder.object <- autoencode(X.train=trainM,nl=nl,N.hidden=N.hidden,
                                 unit.type=unit.type,lambda=lambda,beta=beta,rho=rho,epsilon=epsilon,
                                 optim.method="BFGS",max.iterations=max.iterations,
                                 rescale.flag=TRUE,rescaling.offset=0.001)
