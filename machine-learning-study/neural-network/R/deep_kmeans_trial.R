##############################################################################
## deep_kmeans_trial.R
## This is a simulation study on the idea of deep k-means: 
## Do k-means on the weighted matrix between output layer and hidden layer.
## Author: Meilei
##############################################################################
library(dplyr)
library(autoencoder)
library(fpc)


# make the simulation data ------------------------------------------------

r1 = c(rep(0, 20), rep(1, 80))
r2 = c(rep(1, 20), rep(0, 30), rep(1, 50))
r3 = c(rep(1, 50), rep(0, 50))

trainM = matrix(c(rep(r1, 60), rep(r2, 30), rep(r3, 80)), ncol = 100, byrow =T )
noise = matrix(rnorm(170*100, 0 ,1), ncol = 100)
train = trainM + noise

km0 = kmeansruns(train, krange = 2:10, criterion = 'asw', critout=T)
cl0 = km0$cluster
table(cl0)

# set up the autoencoder architecture  ------------------------------------

nl=3                          ## number of layers (default is 3: input, hidden, output)
unit.type = "logistic"        ## specify the network unit type, i.e., the unit's 
## activation function ("logistic" or "tanh")
Nx.patch=100                   ## width of training image patches, in pixels
Ny.patch=1                   ## height of training image patches, in pixels
N.input = Nx.patch*Ny.patch   ## number of units (neurons) in the input layer (one unit per pixel)
N.hidden = 10              ## number of units in the hidden layer
lambda = 0.0002               ## weight decay parameter     
beta = 6                      ## weight of sparsity penalty term 
rho = 0.01                    ## desired sparsity parameter
epsilon <- 0.001              ## a small parameter for initialization of weights 
## as small gaussian random numbers sampled from N(0,epsilon^2)
max.iterations = 1000         ## number of iterations in optimizer

# Train the auto encoder on training.matrix -------------------------------

autoencoder.object <- autoencode(X.train=t(trainM),nl=nl,N.hidden=N.hidden,
                                 unit.type=unit.type,lambda=lambda,beta=beta,rho=rho,epsilon=epsilon,
                                 optim.method="BFGS",max.iterations=max.iterations,
                                 rescale.flag=TRUE,rescaling.offset=0.001)

## Extract weights W and biases b from autoencoder.object:
W <- autoencoder.object$W
b <- autoencoder.object$b
WM = W[[2]]

km = kmeansruns(WM, krange = 2:10, criterion = 'asw', critout=T)
cl = km$cluster
table(cl)
