# The Pleiades problem is a celestial mechanics problem of seven stars in the plane of coordinates xi,
# yi and masses $m_i = i (i = 1, \ldots 7)$.
# If we consider the star i. According to the second law of Newton this
# star is subjected to the action $F_i = m_i p_i''$, where $p_i = (x_i, y_i)^T$. 
# The law of gravity states that the force working on body i implied by body j, denoted by $F_{ij}$, is:
# $F_{ij} = g \frac{mi \dot mj}{||p_i-p_j||^2_2}d_{ij}$
# Here, $F_i$, $F_{ij}$ $\in R^2$, $g$ is the gravitational constant, 
# which is assumed to be one here, and $dij = \frac{p_j-p_i}{||p_j-p_i||_2}$
# represents the direction of the distance between the two stars. 
# According to the principle of superposition of actions, 
# $F_i$ will be the sum of the interactions between body i and all the others,
# $F_i =\sum_{i\neq j} F_{ij}$
# It is easily checked that this problem can be rewritten as:
# $z'' = f(z), z(0)=z_0, z'(0) = z0'$, where $z = (x,y)^T$ and
# $x_i'' = \sum_{j \neq i} m_j(x_j-x_i)/r_{ij}^{3/2}, 
# $y_i'' = \sum_{j \neq i} m_j(y_j-y_i)/r_{ij}^{3/2},
# $r_{ij} = (x_i-x_j)^2 + (y_i-y_j)^2$ 
#During the movement of the 7 bodies several quasi-collisions occur.

require(deTestSet)

pleiade <- function (t, Y, pars) {
   x <- Y[1:7]
   y <- Y[8:14]
   u <- Y[15:21]
   v <- Y[22:28]
   rij <- as.matrix(dist(cbind(x, y)))^3
   dx <- outer(x, x, FUN = function(x, y) x - y) * (1:7)
   dy <- outer(y, y, FUN = function(x, y) x - y) * (1:7)
   fx <- dx / rij ; diag(fx) <- 0
   fy <- dy / rij ; diag(fy) <- 0
   list(c(u, v, colSums(fx), colSums(fy)))
}

yini <- c(x1 = 3, x2 = 3, x3 =-1, x4 =-3,    x5 = 2, x6 =-2,   x7 = 2,
          y1 = 3, y2 =-3, y3 = 2, y4 = 0,    y5 = 0, y6 =-4,   y7 = 4,
          u1 = 0, u2 = 0, u3 = 0, u4 = 0,    u5 = 0, u6 =1.75, u7 = -1.5,
          v1 = 0, v2 = 0, v3 = 0, v4 =-1.25, v5 = 1, v6 = 0,   v7 = 0)

print(system.time(
out <- ode(func = pleiade, parms = NULL, y = yini, times = seq(0, 50, 0.01), 
             atol = 1e-10, rtol = 1e-10)
))
                                                                                                
par(mfrow = c(3,3))
for (i in 1:7) plot(out[,i+1], out[,i+8], type = "l", main = paste("body ",i),
     xlab = "x", ylab = "y")


plot(0, 0 , type = "n", main = "ALL",
     xlab = "x", ylab = "y", xlim = c(-3, 4), ylim = c(-4, 5))
for (i in 1:7) lines(out[,i+1], out[,i+8], col = i, lwd = 2)

matplot(out[,"time"], out[, c("u1", "u7")], type = "l", lwd = 2, 
        col = c("black", "grey"), lty = 1, xlab = "time", 
        ylab = "velocity", main = "stars 1 and 7")


