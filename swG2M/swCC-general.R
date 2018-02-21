path <- "C:/Users/Rosa/Downloads"

source(paste(path, "swCC-models2.r", sep = "/"))

#model
mdl <- "phogkp"
pmdl <- "NCCchk"

#saving file
sfile <- "phogkp-ch.Rdata"

#graphical conditions
mVariables <- c("y0","q0", "p0")#"z0","s0","r0","y0","q0", "p0",
cols <- c("z0"="#B5231B", "y0"="#01386F", "r0"="#73D8D8", "q0"="#0B7ABF", "s0"="#F09837", "p0"="#E35D01")


#bifurcation analysis
mol <- "z0" #name of molecule
vn1 <- "ch" #names of parameter to change
vv1 <- seq(0,35,0.01) #define parameter values to check
jac <- T #calculate the Jacobian?
st <- seq(0.0,3,0.1) #initial condidions of state variables
tt <- c(0,200,0.1) #time frame

#initial conditions
molec <- c("z", "y","q", "p", "r", "s") #c("kn","ph")#
names(molec)<- c("Cdk1", "PP2A", "PP1", "Gwl", "Cdc25", "Wee1")
# names(molec)<- c("Kin0", "Pho0", "Pho1", "Kin2", "Kin1", "Pho2")
nics <- apply(expand.grid(molec,c(0,2)), 1, paste, collapse="")

ics <- matrix(data = NA, nrow = length(st), ncol = length(nics), dimnames = list(NULL,nics))