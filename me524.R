y = c(10.07, 14.73, 17.94, 23.93, 29.61,
                            35.18, 40.02, 44.82, 50.76, 55.05,
                            61.01, 66.4, 75.47, 81.78)
x = c(77.6, 114.9, 141.1, 190.8, 239.9,
                            289, 332.8, 378.4, 434.8, 477.3,
                            536.8, 593.1, 689.1, 760)
Misrala <- data.frame(y,x)

Loss_f1 <- function(b){
  b1 <- b[1]
  b2 <- b[2]
  Q <- sum( (y - b1 * (1 - exp(-x * b2)))^2 )
  return(Q)
}

Grad_f1 <- function(b){
  b1 <- b[1]
  b2 <- b[2]
  del_b1 <- sum(2 * ( exp(-b2 * x) - 1) * (y - (1 - exp(-b2 * x)) * b1))
  del_b2 <- sum(-2 * b1 * x * (y - b1 * (1 - exp(-x * b2))) * exp(-x * b2))
  return(c(del_b1,del_b2))
}

Hess_f1 <- function(b){
  b1 <- b[1]
  b2 <- b[2]
  del_b1_2 <- sum( 2 *(exp( -x * b2) - 1) ^2 )
  del_b1_b2 <- sum( -2 * x * exp(-2 * x * b2) * ((y - 2 * b1) * exp(x * b2) + 2 * b1))
  del_b2_2 <- sum((2 * b1 * x^2 * exp(-2 * x * b2) * ((y - b1) * exp(x * b2) + 2 * b1)))
  return(matrix(c(del_b1_2, del_b1_b2, del_b1_b2, del_b2_2),nrow = 2))
}

# Loss_f1 <- function(b1,b2){
#   return(sum( (y - b1 * (1 - exp(-x * b2)))^2 ))
# }
# 
# x_graph <- seq(0,2000, length.out = 1000)
# y_graph <- seq(0,0.005, length.out = 100)
# 
# z <- matrix(0, nrow = length(x_graph), ncol = length(y_graph))
# 
# # Loop through all combinations of b1 and b2
# for (i in 1:length(x_graph)) {
#   for (j in 1:length(y_graph)) {
#     z[i, j] <- Loss_f1(x_graph[i], y_graph[j])
#   }
# }
# contour(x_graph, y_graph, z, main = "Contour Plot of Loss Function", xlab = "b1", ylab = "b2",nlevels = 30)
# point_b1 <- 2.3894212918E+02
# point_b2 <- 5.5015643181E-04
# 
# # Add the point to the plot using points()
# points(point_b1, point_b2, col = "red", pch = 19, cex = 1.5) 
# points(500,0.0001,col = "blue", pch = 19, cex = 1)
# points(250,0.0005,col = "blue", pch = 19, cex = 1)
# 
# arrows(500, 0.0001, 500 + 9.620245e+02 , 0.0001 + 1.873369e-04,
#        col = "blue", length = 0.1, lwd = 2)
# arrows(250, 0.0005, 250 + 5.058065e+00 , 0.0005  + 1.058997e-05,
#        col = "blue", length = 0.1, lwd = 2)
find_step <- function(hess,grad,b, g = (25:1)/25){
  for(i in g){
    if(rcond(hess(b -  i*solve(hess(b),grad(b) ))) > 1e-15){
      print(i)
      return(i*solve(hess(b),grad(b)))
      break
    }
  }
}

Newton_Raph <- function(b, hess, grad) {
  i <- 0
  b0 <- b
  while (i <= 1000) {
    if (rcond(hess(b0)) <= 1e-15) {
      print("Non-invertible Hessian at the current point")
      return(b0)
    }
    step <- solve(hess(b0), grad(b0))
    gamma <- 1
    b <- b0 - step
    while (rcond(hess(b)) < 1e-16 && gamma > 0) {
      gamma <- gamma - 0.01
      b <- b0 - gamma * step
    }
    if (gamma <= 0) {
      print("Failed to find a valid gamma")
      return(b0)
    }
    b0 <- b
    if (max(abs(step)) <= 1e-16) {
      print("Converged")
      return(b0)
    }
    i <- i + 1
  }
  print("Reached maximum iterations")
  return(b0)
}

Newton_Raph(c(250, 0.0005),hess = Hess_f1,grad = Grad_f1)
Newton_Raph(c(500, 0.0001),hess = Hess_f1,grad = Grad_f1)

targetStep <- function(b, g = (1:50)/50,func,hess,grad){
  pk <- solve(hess(b),grad(b))
  Loss <- numeric(length(g))
  for(i in 1:length(g)){
    Loss[i] <- func(b-g[i] * pk)
  }
  return(data.frame(g = g, Loss = Loss))
}

Newton_Raph_Line <- function(b,func,hess,grad){
  i <- 0
  while(i <= 1000) {
    b0 <- b
    eta <- targetStep(b, func = func,grad = grad,hess = hess)[which.min(targetStep(b, func = func,grad = grad,hess = hess)[,2]),1]
    b <- b0 - eta * solve(hess(b0),grad(b0))
    gamma <- 1
    while (rcond(hess(b)) < 1e-15 && gamma > 0) {
      gamma <- gamma - 0.04
      b <- b0 - gamma * eta * solve(hess(b0),grad(b0))
    }
    print(grad(b))
    if (gamma <= 0) {
      print("Failed to find a valid gamma")
      return(b0)
    }
    if (max(abs(gamma * eta * solve(hess(b0),grad(b0)))) <= 0.0001) {
      print("Converged")
      return(b0)
    }
    i <- i + 1
  }
  return(b0)
}

Newton_Raph_Line(c(250, 0.0005),func = Loss_f1,grad = Grad_f1,hess = Hess_f1)
Newton_Raph_Line(c(500, 0.0001),func = Loss_f1,grad = Grad_f1,hess = Hess_f1)

targetStep_grad <- function(b, g = (1:25)/25,func,grad){
  pk <- grad(b)
  Loss <- numeric(length(g))
  for(i in 1:length(g)){
    Loss[i] <- func(b-g[i]*pk)
  }
  return(data.frame(g = g, Loss = Loss))
}

Grad_Desc <- function(b,func,grad){
  i <- 0
  b0 <- c(rep(0,length(b)))
  while (i <= 1000) {
    b0 <- b
    eta <- targetStep_grad(b,func = func,grad = grad)[which.min(targetStep_grad(b,func = func,grad = grad)[,2]),1]
    print(eta)
    b <- b0 - eta * (grad(b0))
    if(max(abs(b0 - b)) <= 0.0001){
      return(b0)
    }
    i <- i +1
  }
  return(b)
}

Grad_Desc(c(250, 0.0005),grad = Grad_f1,func = Loss_f1)
Grad_Desc(c(500, 0.0001),grad = Grad_f1,func = Loss_f1)

Loss_f2 <- function(b){
  for (i in seq_along(b)) {
    assign(paste0("b_", i), b[i])
  }
  return( sum( (y - (b_1 + b_2 * x + b_3 * x^2 + b_4 * x^3)/(1 + b_5 * x + b_6 * x^2 + b_7 * x^3 ))^2 ))
}

Grad_f2 <- function(b){
  for (i in seq_along(b)) {
    assign(paste0("b_", i), b[i])
  }
  
  del_b1 <- sum( -(2 * (y - (b_1 + b_4 * x^3 + b_3 * x^2 + b_2 * x) / (b_7 * x^3 + b_6 * x^2 + b_5 * x + 1))) / (b_7 * x^3 + b_6 * x^2 + b_5 * x + 1))
  del_b2 <- sum( -(2 * x * (y - (x * b_2 + b_4 * x^3 + b_3 * x^2 + b_1) / (b_7 * x^3 + b_6 * x^2 + b_5 * x + 1))) / (b_7 * x^3 + b_6 * x^2 + b_5 * x + 1))
  del_b3 <- sum( -(2 * x^2 * (y - (x^2 * b_3 + b_4 * x^3 + b_2 * x + b_1) / (b_7 * x^3 + b_6 * x^2 + b_5 * x + 1))) / (b_7 * x^3 + b_6 * x^2 + b_5 * x + 1))
  del_b4 <- sum( -(2 * x^3 * (y - (x^3 * b_4 + b_3 * x^2 + b_2 * x + b_1) / (b_7 * x^3 + b_6 * x^2 + b_5 * x + 1))) / (b_7 * x^3 + b_6 * x^2 + b_5 * x + 1))
  del_b5 <- sum( (2 * x * (b_4 * x^3 + b_3 * x^2 + b_2 * x + b_1) * (y - (b_4 * x^3 + b_3 * x^2 + b_2 * x + b_1) / (x * b_5 + b_7 * x^3 + b_6 * x^2 + 1))) / (x * b_5 + b_7 * x^3 + b_6 * x^2 + 1)^2)
  del_b6 <- sum( (2 * x^2 * (b_4 * x^3 + b_3 * x^2 + b_2 * x + b_1) * (y - (b_4 * x^3 + b_3 * x^2 + b_2 * x + b_1) / (x^2 * b_6 + b_7 * x^3 + b_5 * x + 1))) / (x^2 * b_6 + b_7 * x^3 + b_5 * x + 1)^2)
  del_b7 <- sum( (2 * x^3 * (b_4 * x^3 + b_3 * x^2 + b_2 * x + b_1) * (y - (b_4 * x^3 + b_3 * x^2 + b_2 * x + b_1) / (x^3 * b_7 + b_6 * x^2 + b_5 * x + 1))) / (x^3 * b_7 + b_6 * x^2 + b_5 * x + 1)^2)
  return(c(del_b1,del_b2,del_b3,del_b4,del_b5,del_b6,del_b7))
}

Hess_f2 <- function(b){
  for (i in seq_along(b)) {
    assign(paste0("b_", i), b[i])
  }
  del_b1_2  <- sum( 2 / (b_7 * x^3 + b_6 * x^2 + b_5 * x + 1)^2)
  del_b1_b2 <- sum( (2 * x) / (b_7 * x^3 + b_6 * x^2 + b_5 * x + 1)^2)
  del_b1_b3 <- sum( (2 * x^2) / (b_7 * x^3 + b_6 * x^2 + b_5 * x + 1)^2)
  del_b1_b4 <- sum( (2 * x^3) / (b_7 * x^3 + b_6 * x^2 + b_5 * x + 1)^2)
  del_b1_b5 <- sum( (2 * x * (x * y * b_5 + (b_7 * x^3 + b_6 * x^2 + 1) * y - 2 * b_4 * x^3 - 2 * b_3 * x^2 - 2 * b_2 * x - 2 * b_1)) / (x * b_5 + b_7 * x^3 + b_6 * x^2 + 1)^3)
  del_b1_b6 <- sum( (2 * x^2 * (x^2 * y * b_6 + (b_7 * x^3 + b_5 * x + 1) * y - 2 * b_4 * x^3 - 2 * b_3 * x^2 - 2 * b_2 * x - 2 * b_1)) / (x^2 * b_6 + b_7 * x^3 + b_5 * x + 1)^3)
  del_b1_b7 <- sum( (2 * x^3 * (x^3 * y * b_7 + (b_6 * x^2 + b_5 * x + 1) * y - 2 * b_4 * x^3 - 2 * b_3 * x^2 - 2 * b_2 * x - 2 * b_1)) / (x^3 * b_7 + b_6 * x^2 + b_5 * x + 1)^3)
  del_b2_2  <- sum( (2 * x^2) / (b_7 * x^3 + b_6 * x^2 + b_5 * x + 1)^2)
  del_b2_b3 <- sum( (2 * x^3) / (b_7 * x^3 + b_6 * x^2 + b_5 * x + 1)^2)
  del_b2_b4 <- sum( (2 * x^4) / (b_7 * x^3 + b_6 * x^2 + b_5 * x + 1)^2)
  del_b2_b5 <- sum( (2 * x^2 * (x * y * b_5 + (b_7 * x^3 + b_6 * x^2 + 1) * y - 2 * b_4 * x^3 - 2 * b_3 * x^2 - 2 * b_2 * x - 2 * b_1)) / (x * b_5 + b_7 * x^3 + b_6 * x^2 + 1)^3)
  del_b2_b6 <- sum( (2 * x^3 * (x^2 * y * b_6 + (b_7 * x^3 + b_5 * x + 1) * y - 2 * b_4 * x^3 - 2 * b_3 * x^2 - 2 * b_2 * x - 2 * b_1)) / (x^2 * b_6 + b_7 * x^3 + b_5 * x + 1)^3)
  del_b2_b7 <- sum( (2 * x^4 * (x^3 * y * b_7 + (b_6 * x^2 + b_5 * x + 1) * y - 2 * b_4 * x^3 - 2 * b_3 * x^2 - 2 * b_2 * x - 2 * b_1)) / (x^3 * b_7 + b_6 * x^2 + b_5 * x + 1)^3)
  del_b3_2  <- sum( (2 * x^4) / (b_7 * x^3 + b_6 * x^2 + b_5 * x + 1)^2)
  del_b3_b4 <- sum( (2 * x^5) / (b_7 * x^3 + b_6 * x^2 + b_5 * x + 1)^2)
  del_b3_b5 <- sum( (2 * x^3 * (x * y * b_5 + (b_7 * x^3 + b_6 * x^2 + 1) * y - 2 * b_4 * x^3 - 2 * b_3 * x^2 - 2 * b_2 * x - 2 * b_1)) / (x * b_5 + b_7 * x^3 + b_6 * x^2 + 1)^3)
  del_b3_b6 <- sum( (2 * x^4 * (x^2 * y * b_6 + (b_7 * x^3 + b_5 * x + 1) * y - 2 * b_4 * x^3 - 2 * b_3 * x^2 - 2 * b_2 * x - 2 * b_1)) / (x^2 * b_6 + b_7 * x^3 + b_5 * x + 1)^3)
  del_b3_b7 <- sum( (2 * x^5 * (x^3 * y * b_7 + (b_6 * x^2 + b_5 * x + 1) * y - 2 * b_4 * x^3 - 2 * b_3 * x^2 - 2 * b_2 * x - 2 * b_1)) / (x^3 * b_7 + b_6 * x^2 + b_5 * x + 1)^3)
  del_b4_2  <- sum( (2 * x^6) / (b_7 * x^3 + b_6 * x^2 + b_5 * x + 1)^2)
  del_b4_b5 <- sum( (2 * x^4 * (x * y * b_5 + (b_7 * x^3 + b_6 * x^2 + 1) * y - 2 * b_4 * x^3 - 2 * b_3 * x^2 - 2 * b_2 * x - 2 * b_1)) / (x * b_5 + b_7 * x^3 + b_6 * x^2 + 1)^3)
  del_b4_b6 <- sum( (2 * x^5 * (x^2 * y * b_6 + (b_7 * x^3 + b_5 * x + 1) * y - 2 * b_4 * x^3 - 2 * b_3 * x^2 - 2 * b_2 * x - 2 * b_1)) / (x^2 * b_6 + b_7 * x^3 + b_5 * x + 1)^3)
  del_b4_b7 <- sum( (2 * x^6 * (x^3 * y * b_7 + (b_6 * x^2 + b_5 * x + 1) * y - 2 * b_4 * x^3 - 2 * b_3 * x^2 - 2 * b_2 * x - 2 * b_1)) / (x^3 * b_7 + b_6 * x^2 + b_5 * x + 1)^3)
  del_b5_2 <-  sum( -(2 * x^2 * (b_4 * x^3 + b_3 * x^2 + b_2 * x + b_1) * (2 * x * y * b_5 + (2 * b_7 * x^3 + 2 * b_6 * x^2 + 2) * y - 3 * b_4 * x^3 - 3 * b_3 * x^2 - 3 * b_2 * x - 3 * b_1)) / (x * b_5 + b_7 * x^3 + b_6 * x^2 + 1)^4)
  del_b5_b6 <- sum( -(2 * x^3 * (b_4 * x^3 + b_3 * x^2 + b_2 * x + b_1) * (2 * x^2 * y * b_6 + (2 * b_7 * x^3 + 2 * b_5 * x + 2) * y - 3 * b_4 * x^3 - 3 * b_3 * x^2 - 3 * b_2 * x - 3 * b_1)) / (x^2 * b_6 + b_7 * x^3 + b_5 * x + 1)^4)
  del_b5_b7 <- sum( -(2 * x^4 * (b_4 * x^3 + b_3 * x^2 + b_2 * x + b_1) * (2 * x^3 * y * b_7 + (2 * b_6 * x^2 + 2 * b_5 * x + 2) * y - 3 * b_4 * x^3 - 3 * b_3 * x^2 - 3 * b_2 * x - 3 * b_1)) / (x^3 * b_7 + b_6 * x^2 + b_5 * x + 1)^4) 
  del_b6_2 <-  sum( -(2 * x^4 * (b_4 * x^3 + b_3 * x^2 + b_2 * x + b_1) * (2 * x^2 * y * b_6 + (2 * b_7 * x^3 + 2 * b_5 * x + 2) * y - 3 * b_4 * x^3 - 3 * b_3 * x^2 - 3 * b_2 * x - 3 * b_1)) / (x^2 * b_6 + b_7 * x^3 + b_5 * x + 1)^4)
  del_b6_b7 <- sum( -(2 * x^5 * (b_4 * x^3 + b_3 * x^2 + b_2 * x + b_1) * (2 * x^3 * y * b_7 + (2 * b_6 * x^2 + 2 * b_5 * x + 2) * y - 3 * b_4 * x^3 - 3 * b_3 * x^2 - 3 * b_2 * x - 3 * b_1)) / (x^3 * b_7 + b_6 * x^2 + b_5 * x + 1)^4)
  del_b7_2 <-  sum( -(2 * x^6 * (b_4 * x^3 + b_3 * x^2 + b_2 * x + b_1) * (2 * x^3 * y * b_7 + (2 * b_6 * x^2 + 2 * b_5 * x + 2) * y - 3 * b_4 * x^3 - 3 * b_3 * x^2 - 3 * b_2 * x - 3 * b_1)) / (x^3 * b_7 + b_6 * x^2 + b_5 * x + 1)^4)
  
  hessian <- matrix(c(
    del_b1_2,  del_b1_b2, del_b1_b3, del_b1_b4, del_b1_b5, del_b1_b6, del_b1_b7,
    del_b1_b2, del_b2_2,  del_b2_b3, del_b2_b4, del_b2_b5, del_b2_b6, del_b2_b7,
    del_b1_b3, del_b2_b3, del_b3_2,  del_b3_b4, del_b3_b5, del_b3_b6, del_b3_b7,
    del_b1_b4, del_b2_b4, del_b3_b4, del_b4_2,  del_b4_b5, del_b4_b6, del_b4_b7,
    del_b1_b5, del_b2_b5, del_b3_b5, del_b4_b5, del_b5_2,  del_b5_b6, del_b5_b7,
    del_b1_b6, del_b2_b6, del_b3_b6, del_b4_b6, del_b5_b6, del_b6_2,  del_b6_b7,
    del_b1_b7, del_b2_b7, del_b3_b7, del_b4_b7, del_b5_b7, del_b6_b7, del_b7_2
  ), nrow = 7, ncol = 7, byrow = TRUE)
  
return(hessian)
}


y = c(80.574, 84.248, 87.264, 87.195, 89.076, 89.608,
      89.868, 90.101, 92.405, 95.854, 100.696, 101.06,
      401.672, 390.724, 567.534, 635.316, 733.054,
      759.087, 894.206, 990.785, 1090.109, 1080.914,
      1122.643, 1178.351, 1260.531, 1273.514, 1288.339,
      1327.543, 1353.863, 1414.509, 1425.208, 1421.384,
      1442.962, 1464.35, 1468.705, 1447.894, 1457.628)
x = c(-3.067, -2.981, -2.921, -2.912, -2.84, -2.797,
      -2.702, -2.699, -2.633, -2.481, -2.363, -2.322,
      -1.501, -1.46, -1.274, -1.212, -1.1, -1.046,
      -0.915, -0.714, -0.566, -0.545, -0.4, -0.309,
      -0.109, -0.103, 0.01, 0.119, 0.377, 0.79, 0.963,
      1.006, 1.115, 1.572, 1.841, 2.047, 2.2)
Thurber <- data.frame(y,x)

Res_1_N_r <- Newton_Raph_Line(c(1000, 1000, 400, 40, 0.7,0.3,0.03),func = Loss_f2,grad = Grad_f2,hess = Hess_f2)
Res_2_N_r <- Newton_Raph_Line(c(1300, 1500, 500, 75, 1, 0.4, 0.05),func = Loss_f2,grad = Grad_f2,hess = Hess_f2)

Res_1_G_d <- Grad_Desc(c(1000, 1000, 400, 40, 0.7,0.3,0.03),func = Loss_f2,grad = Grad_f2)
Res_2_G_d <- Grad_Desc(c(1300, 1500, 500, 75, 1, 0.4, 0.05),func = Loss_f2,grad = Grad_f2)
