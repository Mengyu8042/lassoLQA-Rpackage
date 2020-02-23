#' Select variables via lasso in linear regression.
#'
#' @param x A matrix.
#' @param y A vector.
#' @return The result of selected variables.
#' @examples
#' lassoLQA(matrix(rnorm(90),30,3), matrix(rnorm(90),30,3)%*%c(2,1,0)+rnorm(30,0,4))
#' lassoLQA(matrix(rnorm(1200),200,6), matrix(rnorm(1200),200,6)%*%c(3,2,1,0,0,0)+rnorm(200,0,4))

#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

lassoLQA <- function(x, y) {
  n = length(y)
  p = length(x[1,])
  K = 1000  #遍历lambda的次数
  b = matrix(0, p, K)

  #确定选择变量的准则：AIC、BIC、GCV
  SSe = c()
  for (i in 1:K) SSe[i] = 0
  AIC = c()
  BIC = c()
  GCV = c()

  #使用牛顿迭代法求解
  bb = matrix(NA, p, 10^5)
  v = rnorm(p, 0, 0.1)
  bb[,1] = solve(t(x) %*% x + diag(v)) %*% t(x) %*% y

  for(k in 1:K){
    lambda = 0.01 * k
    j = 2
    repeat {
      d2Q = t(x) %*% x
      dQ = -2 * t(x) %*% (y - x %*% bb[,j-1])
      sigma = matrix(NA, p, p)
      for (s1 in 1:p) {
        for (s2 in 1:p) {
          if (s1 == s2) {
            if (bb[s1,j-1] != 0)
              sigma[s1,s2] = lambda / abs(bb[s1,j-1])
            else
              sigma[s1,s2] = lambda * 10000
          }
          else sigma[s1,s2] = 0
        }
      }
      U = sigma %*% bb[,j-1]
      bb[,j] = bb[,j-1] - solve(d2Q + n * sigma) %*% (dQ + n * U)
      for (i in 1:p) {
        if (bb[i,j] < 1e-2) {
          bb[i,j] = 0
        }
      }
      if (sqrt(sum((bb[,j]-bb[,j-1])^2)) < 1e-5) {
        s = j
        break
      }
      j = j+1
    }
    b[,k] = bb[,s]   #每个lambda下估出的参数

    #AIC & BIC
    for (i in 1:n) {
      SSe[k] = SSe[k] + (y[i] - t(x[i,]) %*% b[,k])^2
    }
    q = 0
    for (l in 1:p) {
      if (b[l,k] != 0) q = q + 1
    }

    AIC[k] = n * log(SSe[k]) + 2 * q
    BIC[k] = n * log(SSe[k]) + log(n) * q

    #GCV
    sigma2 = matrix(NA, p, p)
    for (s1 in 1:p) {
      for (s2 in 1:p) {
        if (s1 == s2) {
          if (b[s1,k] != 0)
            sigma2[s1,s2] = lambda / abs(b[s1,k])
          else
            sigma2[s1,s2] = lambda * 10000
        }
        else sigma2[s1,s2] = 0
      }
    }
    P = x %*% solve(t(x) %*% x + n * sigma2) %*% t(x)
    e = sum(diag(P))
    GCV[k] = SSe[k] / ((1-e/n)^2)

  }

  #输出变量选择的结果
  k1 = which.min(AIC)
  sel_AIC = b[,k1]

  k2 = which.min(BIC)
  sel_BIC = b[,k2]

  k3 = which.min(GCV)
  sel_GCV = b[,k3]

  lam = c(0.01*k1, 0.01*k2, 0.01*k3)  #lambda_AIC, lambda_BIC, lambda_GCV
  b_lasso = list(sel_AIC, sel_BIC, sel_GCV)  #b_AIC, b_BIC, b_GCV

  result <- list(criteria = c("AIC", "BIC", "GCV"),
                lambda = c(0.01*k1, 0.01*k2, 0.01*k3),
                b_lasso = matrix(c(sel_AIC, sel_BIC, sel_GCV), nrow = p))
  return(result)

}
