p <- c(0,
       0.05,
       0.10,
       0.15,
       0.20,
       0.25,
       0.30,
       0.35,
       0.40,
       0.45,
       0.50,
       0.55,
       0.60)

##Lower bound
h1 <- c(0.396,
        0.355,
        0.320,
        0.291,
        0.290,
        0.275,
        0.275,
        0.300,
        0.270,
        0.250
)

h2 <- c(0.857,
        0.805,
        0.767,
        0.722,
        0.671,
        0.647,
        0.634,
        0.642,
        0.639,
        0.626
)

h3 <- c(0.989,
        0.980,
        0.973,
        0.964,
        0.950,
        0.936,
        0.928,
        0.913,
        0.905,
        0.911
)

pdf("power_at_lower_bound.pdf")
plot(p[1:length(h1)], h1, main = expression(paste("Power at ", psi, " = ", psi[l]- h/sqrt(n))), 
     type = "o", col = "red", xlim = c(0, 0.5), ylim = c(0,1), ylab = "Power", xlab = "p",
     cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5)
par(new = TRUE)
plot(p[1:length(h2)], h2, type = "o", col = "green", xlim = c(0, 0.5), ylim = c(0,1), ylab = NA, xlab = NA, cex.axis = 1.5)
par(new = TRUE)
plot(p[1:length(h3)], h3, type = "o", col = "blue", xlim = c(0, 0.5), ylim = c(0,1), ylab = NA, xlab = NA, cex.axis = 1.5)

legend("bottomleft", legend = c("h = 1", "h = 2", "h = 3"), lty = 1, lwd = 2, col = c("red", "green", "blue"))
dev.off()

##Upper bound
h1 <- c(0.409,
        0.404,
        0.403,
        0.411,
        0.402,
        0.400,
        0.397,
        0.398,
        0.397,
        0.381,
        0.590,
        0.361,
        0.389
)

h2 <- c(0.881,
        0.876,
        0.877,
        0.878,
        0.878,
        0.878,
        0.878,
        0.876,
        0.870,
        0.848,
        0.984,
        0.852,
        0.891
)

h3 <- c(0.996,
        0.997,
        0.996,
        0.996,
        0.996,
        0.996,
        0.997,
        0.996,
        0.996,
        0.994,
        1,
        0.993,
        0.990
)

pdf("power_at_upper_bound.pdf")
plot(p, h1, main = expression(paste("Power at ", psi, " = ", psi[u] + h/sqrt(n))), 
     type = "o", col = "red", xlim = c(0, 0.6), ylim = c(0,1), ylab = "Power", xlab = "p",
     cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5)
par(new = TRUE)
plot(p, h2, type = "o", col = "green", xlim = c(0, 0.6), ylim = c(0,1), ylab = NA, xlab = NA, cex.axis = 1.5)
par(new = TRUE)
plot(p, h3, type = "o", col = "blue", xlim = c(0, 0.6), ylim = c(0,1), ylab = NA, xlab = NA, cex.axis = 1.5)

legend("bottomleft", legend = c("h = 1", "h = 2", "h = 3"), lty = 1, lwd = 2, col = c("red", "green", "blue"))
dev.off()