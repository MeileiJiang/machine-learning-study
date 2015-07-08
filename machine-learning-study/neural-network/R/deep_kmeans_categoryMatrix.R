##############################################################################
## deep_kmeans_categoryMatrix.R
## This is a test on categoryMatrix data under the framework of deep k-means: 
## Do k-means on the weighted matrix between output layer and hidden layer.
## Author: Meilei
##############################################################################
library(dplyr)
library(autoencoder)
library(fpc)


# data processing ------------------------------------------------

load("neural-network/Data/categoryMatrix.RData")
dim(category_matrix)
# [1] 98 2301
## each row is a data object: a 177 * 11 0-1 matrix 

# set up the autoencoder architecture  ------------------------------------

nl=3                          ## number of layers (default is 3: input, hidden, output)
unit.type = "logistic"        ## specify the network unit type, i.e., the unit's 
## activation function ("logistic" or "tanh")
Nx.patch=13                   ## width of training image patches, in pixels
Ny.patch=177                   ## height of training image patches, in pixels
N.input = Nx.patch*Ny.patch   ## number of units (neurons) in the input layer (one unit per pixel)
N.hidden = 177              ## number of units in the hidden layer
lambda = 0.0002               ## weight decay parameter     
beta = 6                      ## weight of sparsity penalty term 
rho = 0.01                    ## desired sparsity parameter
epsilon <- 0.001              ## a small parameter for initialization of weights 
## as small gaussian random numbers sampled from N(0,epsilon^2)
max.iterations = 1000         ## number of iterations in optimizer

# Train the auto encoder on training.matrix -------------------------------

autoencoder.object <- autoencode(X.train=t(category_matrix),nl=nl,N.hidden=N.hidden,
                                 unit.type=unit.type,lambda=lambda,beta=beta,rho=rho,epsilon=epsilon,
                                 optim.method="BFGS",max.iterations=max.iterations,
                                 rescale.flag=TRUE,rescaling.offset=0.001)

## Extract weights W and biases b from autoencoder.object:
W <- autoencoder.object$W
b <- autoencoder.object$b
WM = W[[2]]

km = kmeansruns(WM, krange = 2:10, criterion = 'ch', critout=T)
cl = km$cluster
table(cl)
visualize.hidden.units(autoencoder.object,Nx.patch,Ny.patch)
save(autoencoder.object, file = "neural-network/Data/autoencoderResult.RData")
