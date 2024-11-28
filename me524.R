Misrala <- data.frame(y = c(10.07, 14.73, 17.94, 23.93, 29.61,
                            35.18, 40.02, 44.82, 50.76, 55.05,
                            61.01, 66.4, 75.47, 81.78),
                      x = c(77.6, 114.9, 141.1, 190.8, 239.9,
                            289, 332.8, 378.4, 434.8, 477.3,
                            536.8, 593.1, 689.1, 760))

Thurber <- data.frame(y = c(80.574, 84.248, 87.264, 87.195, 89.076, 89.608,
                            89.868, 90.101, 92.405, 95.854, 100.696, 101.06,
                            401.672, 390.724, 567.534, 635.316, 733.054,
                            759.087, 894.206, 990.785, 1090.109, 1080.914,
                            1122.643, 1178.351, 1260.531, 1273.514, 1288.339,
                            1327.543, 1353.863, 1414.509, 1425.208, 1421.384,
                            1442.962, 1464.35, 1468.705, 1447.894, 1457.628),
                      x = c(-3.067, -2.981, -2.921, -2.912, -2.84, -2.797,
                            -2.702, -2.699, -2.633, -2.481, -2.363, -2.322,
                            -1.501, -1.46, -1.274, -1.212, -1.1, -1.046,
                            -0.915, -0.714, -0.566, -0.545, -0.4, -0.309,
                            -0.109, -0.103, 0.01, 0.119, 0.377, 0.79, 0.963,
                            1.006, 1.115, 1.572, 1.841, 2.047, 2.2))

Loss_f1 <- function(b){
  Q <- sum( (Misrala$y - b[1] * (1 - exp(-Misrala$x * b[2])))^2 )
  return(Q)
}

Grad_f1 <- function(b){
  d_b1 <- sum(2 * ( exp(-b[2] * Misrala$x) - 1) * (Misrala$y - (1 - exp(-b[2] * Misrala$x)) * b[1]))
  d_b2 <- sum(-2 * b[1] * Misrala$x * (Misrala$y - b[1] * (1 - exp(-Misrala$x * b[2]))) * exp(-Misrala$x * b[2]))
  return(c(d_b1,d_b2))
}

Hess_f1 <- function(b){
  d_b1_2 <- sum( 2 *(exp( -Misrala$x * b[2]) - 1) ^2 )
  d_b1_b2 <- sum( 2 * Misrala$x * exp(-2 * b[2] * Misrala$x) * ((2 * exp(b[2] * Misrala$x) - 2) * b[1] - exp(b[2] * Misrala$x) * Misrala$y))
  d_b2_2 <- sum((2 * b[1] * Misrala$x^2 * exp(-2 * Misrala$x * b[2]) * ((Misrala$y - b[1]) * exp(Misrala$x * b[2]) + 2 * b[1])))
  return(matrix(c(d_b1_2, d_b1_b2, d_b1_b2, d_b2_2),nrow = 2))
}

Loss_f2 <- function(b){
  return( sum( (Thurber$y - (b[1] + b[2] * Thurber$x + b[3] * Thurber$x^2 + b[4] * Thurber$x^3)/(1 + b[5] * Thurber$x + b[6] * Thurber$x^2 + b[7] * Thurber$x^3 ))^2 ))
}

Grad_f2 <- function(b){
  x <- Thurber$x
  y <- Thurber$y
  denom <- (b[7] * x^3 + b[6] * x^2 + b[5] * x + 1)
  num <- (b[4] * x^3 + b[3] * x^2 + b[2] * x + b[1])
  d_b1 <- sum( -(2 * (y - (b[1] + b[4] * x^3 + b[3] * x^2 + b[2] * x) / denom)) / denom)
  d_b2 <- sum( -(2 * x * (y - (x * b[2] + b[4] * x^3 + b[3] * x^2 + b[1]) / denom)) / denom)
  d_b3 <- sum( -(2 * x^2 * (y - (x^2 * b[3] + b[4] * x^3 + b[2] * x + b[1]) / denom)) / denom)
  d_b4 <- sum( -(2 * x^3 * (y - (x^3 * b[4] + b[3] * x^2 + b[2] * x + b[1]) / denom)) / denom)
  d_b5 <- sum( (2 * x * num * (y - num / (x * b[5] + b[7] * x^3 + b[6] * x^2 + 1))) / (x * b[5] + b[7] * x^3 + b[6] * x^2 + 1)^2)
  d_b6 <- sum( (2 * x^2 * num * (y - num / (x^2 * b[6] + b[7] * x^3 + b[5] * x + 1))) / (x^2 * b[6] + b[7] * x^3 + b[5] * x + 1)^2)
  d_b7 <- sum( (2 * x^3 * num * (y - num / denom)) / denom^2)
  return(c(d_b1,d_b2,d_b3,d_b4,d_b5,d_b6,d_b7))
}

Hess_f2 <- function(b){
  x <- Thurber$x
  y <- Thurber$y
  denom <- c((b[7] * x^3 + b[6] * x^2 + b[5] * x + 1),(2 * x^3 * y * b[7] + (2 * b[6] * x^2 + 2 * b[5] * x + 2)))
  num <- c((2 * b[4] * x^3 - 2 * b[3] * x^2 - 2 * b[2] * x - 2 * b[1]),(b[4] * x^3 + b[3] * x^2 + b[2] * x + b[1]))
  d_b1_2  <- sum( 2 / denom[1]^2)
  d_b1_b2 <- sum( (2 * x) / denom[1]^2)
  d_b1_b3 <- sum( (2 * x^2) / denom[1]^2)
  d_b1_b4 <- sum( (2 * x^3) / denom[1]^2)
  d_b1_b5 <- sum( (2 * x * (x * y * b[5] + (b[7] * x^3 + b[6] * x^2 + 1) * y - num[1])) / denom[1]^3)
  d_b1_b6 <- sum( (2 * x^2 * (x^2 * y * b[6] + (b[7] * x^3 + b[5] * x + 1) * y - num[1])) / denom[1]^3)
  d_b1_b7 <- sum( (2 * x^3 * (x^3 * y * b[7] + (b[6] * x^2 + b[5] * x + 1) * y - num[1])) / denom[1]^3)
  d_b2_2  <- sum( (2 * x^2) / denom[1]^2)
  d_b2_b3 <- sum( (2 * x^3) / denom[1]^2)
  d_b2_b4 <- sum( (2 * x^4) / denom[1]^2)
  d_b2_b5 <- sum( (2 * x^2 * (x * y * b[5] + (b[7] * x^3 + b[6] * x^2 + 1) * y - num[1])) / denom[1]^3)
  d_b2_b6 <- sum( (2 * x^3 * (x^2 * y * b[6] + (b[7] * x^3 + b[5] * x + 1) * y - num[1])) / denom[1]^3)
  d_b2_b7 <- sum( (2 * x^4 * (x^3 * y * b[7] + (b[6] * x^2 + b[5] * x + 1) * y - num[1])) / denom[1]^3)
  d_b3_2  <- sum( (2 * x^4) / denom[1]^2)
  d_b3_b4 <- sum( (2 * x^5) / denom[1]^2)
  d_b3_b5 <- sum( (2 * x^3 * (x * y * b[5] + (b[7] * x^3 + b[6] * x^2 + 1) * y - num[1])) / denom[1]^3)
  d_b3_b6 <- sum( (2 * x^4 * (x^2 * y * b[6] + (b[7] * x^3 + b[5] * x + 1) * y - num[1])) / denom[1]^3)
  d_b3_b7 <- sum( (2 * x^5 * (x^3 * y * b[7] + (b[6] * x^2 + b[5] * x + 1) * y - num[1])) / denom[1]^3)
  d_b4_2  <- sum( (2 * x^6) / denom[1]^2)
  d_b4_b5 <- sum( (2 * x^4 * (x * y * b[5] + (b[7] * x^3 + b[6] * x^2 + 1) * y - num[1])) / denom[1]^3)
  d_b4_b6 <- sum( (2 * x^5 * (x^2 * y * b[6] + (b[7] * x^3 + b[5] * x + 1) * y - num[1])) / denom[1]^3)
  d_b4_b7 <- sum( (2 * x^6 * (x^3 * y * b[7] + (b[6] * x^2 + b[5] * x + 1) * y - num[1])) / denom[1]^3)
  d_b5_2 <-  sum( -(2 * x^2 * num[2]* (2 * x * y * b[5] + (2 * b[7] * x^3 + 2 * b[6] * x^2 + 2) * y - 3 * b[4] * x^3 - 3 * b[3] * x^2 - 3 * b[2] * x - 3 * b[1])) / denom[1]^4)
  d_b5_b6 <- sum( -(2 * x^3 * num[2]* (2 * x^2 * y * b[6] + (2 * b[7] * x^3 + 2 * b[5] * x + 2) * y - 3 * b[4] * x^3 - 3 * b[3] * x^2 - 3 * b[2] * x - 3 * b[1])) / denom[1]^4)
  d_b5_b7 <- sum( -(2 * x^4 * num[2]* (2 * x^3 * y * b[7] + (2 * b[6] * x^2 + 2 * b[5] * x + 2) * y - 3 * b[4] * x^3 - 3 * b[3] * x^2 - 3 * b[2] * x - 3 * b[1])) / denom[1]^4) 
  d_b6_2 <-  sum( -(2 * x^4 * num[2]* (2 * x^2 * y * b[6] + (2 * b[7] * x^3 + 2 * b[5] * x + 2) * y - 3 * b[4] * x^3 - 3 * b[3] * x^2 - 3 * b[2] * x - 3 * b[1])) / denom[1]^4)
  d_b6_b7 <- sum( -(2 * x^5 * num[2]* (2 * x^3 * y * b[7] + (2 * b[6] * x^2 + 2 * b[5] * x + 2) * y - 3 * b[4] * x^3 - 3 * b[3] * x^2 - 3 * b[2] * x - 3 * b[1])) / denom[1]^4)
  d_b7_2 <-  sum( -(2 * x^6 * num[2]* (2 * x^3 * y * b[7] + (2 * b[6] * x^2 + 2 * b[5] * x + 2) * y - 3 * b[4] * x^3 - 3 * b[3] * x^2 - 3 * b[2] * x - 3 * b[1])) / denom[1]^4)
  
  hessian <- matrix(c(
    d_b1_2,  d_b1_b2, d_b1_b3, d_b1_b4, d_b1_b5, d_b1_b6, d_b1_b7,
    d_b1_b2, d_b2_2,  d_b2_b3, d_b2_b4, d_b2_b5, d_b2_b6, d_b2_b7,
    d_b1_b3, d_b2_b3, d_b3_2,  d_b3_b4, d_b3_b5, d_b3_b6, d_b3_b7,
    d_b1_b4, d_b2_b4, d_b3_b4, d_b4_2,  d_b4_b5, d_b4_b6, d_b4_b7,
    d_b1_b5, d_b2_b5, d_b3_b5, d_b4_b5, d_b5_2,  d_b5_b6, d_b5_b7,
    d_b1_b6, d_b2_b6, d_b3_b6, d_b4_b6, d_b5_b6, d_b6_2,  d_b6_b7,
    d_b1_b7, d_b2_b7, d_b3_b7, d_b4_b7, d_b5_b7, d_b6_b7, d_b7_2
  ), nrow = 7, ncol = 7, byrow = TRUE)
  
  return(hessian)
}

is.pos.def <- function(hess) all(eigen(hess, only.values = TRUE)$values > 0)
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

targetStep <- function(b, g = (1:50)/50,func,step){
  Loss <- numeric(length(g))
  for(i in 1:length(g)){
    Loss[i] <- func(b-g[i] * step)
  }
  return(data.frame(g = g, Loss = Loss))
}

Newton_Raph <- function(b, hess, grad,func, line.search) {
  i <- c(1,as.numeric(Sys.time()))
  b0 <- b
  while (i[1] <= 1000) {
    if (rcond(hess(b0)) <= 1e-15) {
      print("Hessiana não invertível")
      return(list(b = b0,Iter = i[1],grad = grad(b0),pos_Def = is.pos.def(hess(b0)), tempo = paste(as.numeric(Sys.time()) - i[2 ], "segundos"),perda = func(b0)))
    }
    step <- solve(hess(b0),grad(b0))
    
    if(line.search){
      eta <- targetStep(b, func = func,step = step)[which.min(targetStep(b, func = func,step = step)[,2]),1]
      
      b <- b0 - eta * step
    }
    if(!line.search){
      b <- b0 - solve(hess(b0)) %*% grad(b0)
    }
    if (norm(b-b0,type = "2") <= 1e-6) {
      print("Convergiu")
      return(list(b = b,Iter = i[1],grad = grad(b),pos_Def = is.pos.def(hess(b)),tempo = paste(as.numeric(Sys.time()) - i[2 ],"segundos"),perda = func(b0)))
    }
    b0 <- b
    i[1] <- i[1] + 1
  }
  print("Número máximo de iterações")
  return(list(b = b,Iter = i[1],grad = grad(b),pos_Def = is.pos.def(hess(b)), tempo = paste(as.numeric(Sys.time()) - i[2], "segundos"),perda = func(b0)))
}

targetStep_grad <- function(b, g = (1:25)/25,func,grad){
  pk <- grad(b)
  Loss <- numeric(length(g))
  for(i in 1:length(g)){
    Loss[i] <- func(b-g[i]*pk)
  }
  return(data.frame(g = g, Loss = Loss))
}

Grad_Desc <- function(b,func,grad,hess){
  i <- c(1,as.numeric(Sys.time()))
  b0 <- c(rep(0,length(b)))
  while (i[1] <= 1000) {
    b0 <- b
    eta <- targetStep_grad(b,func = func,grad = grad)[which.min(targetStep_grad(b,func = func,grad = grad)[,2]),1]
    b <- b0 - eta * (grad(b0))
    if(max(abs(b0 - b)) <= 0.0001){
      print("Convergiu")
      return(list(b = b0,Iter = i[1],grad = grad(b0), tempo = paste(as.numeric(Sys.time()) - i[2], "segundos"),perda = func(b0)))
    }
    i[1] <- i[1] +1
  }
  print("Número máximo de iterações")
  return(list(b = b0,Iter = i[1],grad = grad(b0), tempo = paste(as.numeric(Sys.time()) - i[2], "segundos"),perda = func(b0)))
}

Newton_Raph(c(250, 0.0005),hess = Hess_f1,grad = Grad_f1,func = Loss_f1,line.search = F)
Newton_Raph(c(500, 0.0001),hess = Hess_f1,grad = Grad_f1,func = Loss_f1,line.search = F)

Newton_Raph(c(500, 0.0001),func = Loss_f1,grad = Grad_f1,hess = Hess_f1,line.search = T)
Newton_Raph(c(250, 0.0005),func = Loss_f1,grad = Grad_f1,hess = Hess_f1,line.search = T)

Grad_Desc(c(250, 0.0005),grad = Grad_f1,func = Loss_f1)
Grad_Desc(c(500, 0.0001),grad = Grad_f1,func = Loss_f1)

Newton_Raph(c(1000, 1000, 400, 40, 0.7,0.3,0.03),func = Loss_f2,grad = Grad_f2,hess = Hess_f2,line.search = F)
Newton_Raph(c(1300, 1500, 500, 75, 1, 0.4, 0.05),func = Loss_f2,grad = Grad_f2,hess = Hess_f2,line.search = F)

Newton_Raph(c(1000, 1000, 400, 40, 0.7,0.3,0.03),func = Loss_f2,grad = Grad_f2,hess = Hess_f2,line.search = T)
Newton_Raph(c(1300, 1500, 500, 75, 1, 0.4, 0.05),func = Loss_f2,grad = Grad_f2,hess = Hess_f2,line.search = T)

Grad_Desc(c(1000, 1000, 400, 40, 0.7,0.3,0.03),func = Loss_f2,grad = Grad_f2)
Grad_Desc(c(1300, 1500, 500, 75, 1, 0.4, 0.05),func = Loss_f2,grad = Grad_f2)


reg_hessian <- function(H, epsilon = 1e-3) {
  eig <- eigen(H, symmetric = TRUE)
  eig$values <- pmax(eig$values, epsilon)
  H_regularized <- eig$vectors %*% diag(eig$values) %*% t(eig$vectors)
  return(H_regularized)
}


Newton_Raph_Reg <- function(b, hess, grad,func, line.search) {
  i <- c(1,as.numeric(Sys.time()))
  b0 <- b
  while (i[1] <= 1000) {
    if (rcond(hess(b0)) <= 1e-15) {
      step <- solve(reg_hessian(hess(b0)),grad(b0))
    }
    else step <- solve(hess(b0),grad(b0))
    if(line.search){
      eta <- targetStep(b, func = func,step = step)[which.min(targetStep(b, func = func,step = step)[,2]),1]
      b <- b0 - eta * step
    }
    if(!line.search){
      b <- b0 - step
    }
    if (norm(b-b0,type = "2") <= 1e-6) {
      print("Convergiu")
      return(list(b = b,Iter = i[1],grad = grad(b),pos_Def = is.pos.def(hess(b)),tempo = paste(as.numeric(Sys.time()) - i[2 ],"segundos"),perda = func(b0)))
    }
    b0 <- b
    i[1] <- i[1] + 1
  }
  print("Número máximo de iterações")
  return(list(b = b,Iter = i[1],grad = grad(b),pos_Def = is.pos.def(hess(b)), tempo = paste(as.numeric(Sys.time()) - i[2], "segundos"),perda = func(b0)))
}

Newton_Raph_Reg(c(250, 0.0005),hess = Hess_f1,grad = Grad_f1,func = Loss_f1,line.search = F)
Newton_Raph_Reg(c(500, 0.0001),hess = Hess_f1,grad = Grad_f1,func = Loss_f1,line.search = F)

Newton_Raph_Reg(c(500, 0.0001),func = Loss_f1,grad = Grad_f1,hess = Hess_f1,line.search = T)
Newton_Raph_Reg(c(250, 0.0005),func = Loss_f1,grad = Grad_f1,hess = Hess_f1,line.search = T)

Newton_Raph_Reg(c(1000, 1000, 400, 40, 0.7,0.3,0.03),func = Loss_f2,grad = Grad_f2,hess = Hess_f2,line.search = F)
Newton_Raph_Reg(c(1300, 1500, 500, 75, 1, 0.4, 0.05),func = Loss_f2,grad = Grad_f2,hess = Hess_f2,line.search = F)

Newton_Raph_Reg(c(1000, 1000, 400, 40, 0.7,0.3,0.03),func = Loss_f2,grad = Grad_f2,hess = Hess_f2,line.search = T)
Newton_Raph_Reg(c(1300, 1500, 500, 75, 1, 0.4, 0.05),func = Loss_f2,grad = Grad_f2,hess = Hess_f2,line.search = T)
