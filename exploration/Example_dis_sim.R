p <- 5
D <- matrix(c(0, 1, 2, 4, 4,
              1, 0, 2, 4, 4,
              2, 2, 0, 4, 4,
              4, 4, 4, 0, 2,
              4, 4, 4, 2, 0), 5, 5)

# An real example distance matrix


S <- 1 - D/5
eigen <- eigen(S)
V <- eigen$vectors
J <- diag(rep(1, p)) - 1/p * matrix(1, p, p)
V_sub <- J %*% V
qr(V_sub)$rank

S_new <- J %*% S %*% J
eigen_new <- eigen(S_new)
V_new <- eigen_new$vectors

