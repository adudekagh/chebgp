start_time <- Sys.time()

library(INLA)

# Generated data (from python)
#x <- c(-0.21558522, -0.1329833, 0.08068001, 0.86035093, 0.76285724, -0.64821812, 0.27034352, 0.6716487, -0.19474706, 0.59998129, 0.86436037, -0.30727194, -0.41886138, 0.21547172, -0.84083455, -0.62773623, 0.93725498, -0.56009593, 0.94870193, 0.02186769)
#y <- c(-0.39898039, -0.34918657, -0.11447868, -0.20011555, 0.44504202, -1.64524772, 0.09360612, 0.33980915, -0.42562159, 0.2208891, 0.15916043, -0.65115996, -0.80375508, 0.15148805, -2.57718936, -1.83607748, 0.13775842, -1.47140522, 0.18098801, 0.21166059)

x <- c(0.32080647, 0.26587014, 0.72187268, -0.30351804, 0.2354536, -0.83286119,
       -0.12406531, 0.75046623, 0.17505236, -0.43120703, -0.26736266, 0.41574915,
       -0.13673504, -0.53063669, -0.58944929, 0.71710496, 0.96346998, -0.60412132,
       -0.27278486, -0.94576153)

y <- c(2.76114155, 2.49891946, 3.49267932, 2.19077668, 2.25096868, 2.3522071,
       1.84837552, 3.94204235, 2.39402669, 1.86377928, 1.36257218, 2.76642226,
       1.76145556, 2.32005902, 2.16071458, 4.01261904, 5.19621366, 1.92165487,
       1.70098854, 2.6992103)



name <- "rbf"
#n - nodes value
n <- 200
k <- 1:n
x_nodes <- cos(pi * (2 * k - 1) / (2 * n))
x_nodes <- rev(x_nodes)


mesh <- inla.mesh.1d(x_nodes)
# Define the INLA Model
#formula <- y ~ f(x, model = "rw2")

A = inla.spde.make.A(mesh, loc=x)
######
Xcov = data.frame(value=seq(1,1,length.out=length(y)))
Xcov = as.matrix(Xcov)
colnames(Xcov)
######???
print("Creating stack")
stack <- inla.stack(tag='est',
                    data=list(y=y),
                    
                    effects=list(
                      s=1:mesh$n,
                      Xcov=Xcov),
                    
                    A = list(A, 1)
)


print("Creating model")
prior.median.sd = 75; prior.median.range = 10;
spde = inla.spde2.pcmatern(mesh, prior.range = c(prior.median.range, .5), prior.sigma = c(prior.median.sd, .5))
######
formula = y ~ -1 + f(s, model=spde)
######
prior.median.gaus.sd = 75
family = 'gaussian'
control.family = list(hyper = list(prec = list(prior = "pc.prec", fixed = FALSE, param = c(prior.median.gaus.sd,1))))
print("Calculating result")
res <- inla(formula, data=inla.stack.data(stack),
            control.predictor=list(A = inla.stack.A(stack), compute=T),
            # compute=T to get posterior for fitted values
            family = family,
            # control.family = control.family,
            
            # #control.compute = list(config=T, dic=T, cpo=T, waic=T),
            # # - Model comparisons
            # control.inla = list(int.strategy='eb'),
            # # - faster computation
            # #control.inla = list(int.strategy='grid'),
            # # - More accurate integration over hyper-parameters
            # verbose=F
)

summary(res)




# Fit the Model
#model <- inla(formula, data = data.frame(x = x, y = y), family = "gaussian")

#predictions <- inla.mesh.predict(model, mesh)
#plot(x,y)
#plot(x_nodes, res[["summary.random"]][["s"]][["mean"]], main = "Simple Linear Regression with INLA")


x <- x_nodes
avg <- res[["summary.random"]][["s"]][["mean"]]
sdev <- 2*res[["summary.random"]][["s"]][["sd"]]
plot(x, avg,
     ylim=range(c(avg-sdev, avg+sdev)),
     pch=19, col = 'red', xlab="x", ylab="f(x)",
     main="Plot with 95% confidence interval"
)
arrows(x, avg-sdev, x, avg+sdev, length=0.05, angle=90, code=0, col = 'blue')
points(x, avg, col = 'red',pch=19)
points(x, avg, col = 'black',pch=21)

#model_summary <- summary(model)



coef_matrix <- res[["summary.random"]][["s"]][["mean"]]


value <- coef_matrix
  
extended_vector <- c(rev(value), value[2:(length(value) - 1)])
fourier_coeffs <- Re(fft(extended_vector))
cheb_coef <- fourier_coeffs[1:n] / n
cheb_coef[1] <- cheb_coef[1] / 2
cheb_coef[n] <- cheb_coef[n] / 2
  
coef_matrix <- cheb_coef

test <- abs(coef_matrix)
x_plt <- 1:n
plot(x_plt, test, type = "l", log = "y", main = paste("Order N =", n),
     ylab = "Value", xlab = "Coefficient order", col = "blue", lwd = 2)
legend("topright", legend = "Coefficient values", col = "blue", lwd = 2)

# Assuming x_plt, test, x, avg, sdev are your variables

# Create a data frame with the variables
df <- data.frame(x_plt, test, x, avg, sdev)

# Specify the path where you want to save the CSV file
csv_file_path <- paste("r", name, n, ".csv", sep="")

# Write the data frame to a CSV file
write.csv(df, file = csv_file_path, row.names = FALSE)


end_time <- Sys.time()
elapsed_time <- end_time - start_time

# Print the elapsed time
print(elapsed_time)