# BDA Practical #1 Solutions
# J. Gair and R. Amaros-Salvador, 25 January 2019

# ------------------------------------------------------------------------
#  Problem 1 - analysis of binomial data
# prior is Be(9.2, 13.8); data are 15 successes in 20 trials
alpha.prior <- 9.2
beta.prior  <- 13.8
n           <- 20
num.succ    <- 15
p.hat       <- num.succ/n

# (a) The posterior dist'n is alpha' = alpha.prior + x
#                         and beta'  = beta.prior + n-x
alpha.post  <- alpha.prior+num.succ
beta.post   <- beta.prior +n-num.succ
cat("Posterior alpha=",alpha.post,"and beta=",beta.post,"\n")
#   Posterior alpha= 24.2 and beta= 18.8 

# Q ii. posterior mean=a'/(a'+b') 
cat("Posterior mean=",alpha.post/(alpha.post+beta.post),"\n")
# Posterior mean= 0.5627907 

x <- seq(0.1,0.99,by=0.01)
plot(x,dbeta(x,shape1=alpha.post,shape2=beta.post),type="l")

# Highest 95% Interval: need code for this
require(TeachingDemos)
cat("The 95% HPD for probability of success is",
    hpd(qbeta, shape1=alpha.post, shape2=beta.post),"\n")
# The 95% HPD for probability of success is 0.4162469 0.7077353 

#Q iii. not the same as above:
cat("0.025 and 0.975=",qbeta(p=c(0.025,0.975),shape1=alpha.post,
                             shape2=beta.post))
# 0.414, 0.706

#Q iv. Pr(p >0.6)
cat("Pr(p>-.6)=",1-pbeta(0.6,shape1=alpha.post,shape2=beta.post),"\n")
# Pr(p>-.6)= 0.3156323 

#Q v. with Uniform prior
alpha.post2 <- 1+num.succ
beta.post2  <- 1+n-num.succ
cat("Pr(p>0.6=",1-pbeta(0.6,shape1=alpha.post2,shape2=beta.post2),"\n")
# Pr(p>0.6= 0.9042598 

#Q v. with Jeffreys prior
alpha.post2 <- 0.5+num.succ
beta.post2  <- 0.5+n-num.succ
cat("Pr(p>0.6=",1-pbeta(0.6,shape1=alpha.post2,shape2=beta.post2),"\n")
# Pr(p>0.6= 0.9176342 


#Q vi. Posterior marginal for event (posterior predictive dist'n):
# Pr(x>=25|n=40); need the beta-binomial distribution

#
incomplete.beta.fun <- function(alpha,beta) {
  out <- gamma(alpha)*gamma(beta)/gamma(alpha+beta)
  return(out)
}

dbetabinomial <- function(x,alpha,beta,n) {
  icb.1 <- incomplete.beta.fun(alpha,beta)
  icb.2 <- incomplete.beta.fun(x+alpha,n-x+beta)
  bin.coef <- gamma(n+1)/(gamma(x+1)*gamma(n-x+1))
  out <- (1/icb.1)*bin.coef*icb.2
  return(out)
}

pbetabinomial <- function(X,alpha,beta,n) {
  out <- sum(dbetabinomial(x=0:X,alpha,beta,n))
  return(out)
}

cat("Posterior Pr(X>=25|n=40) based on Beta-Binom(alpha.post, beta.post, n) \n")
cat(1-pbetabinomial(X=24,alpha=alpha.post,beta=beta.post,n=40),"\n")
# 0.3290134 

#Q vii. Bayesian P-value based on the Prior
#    Probability of at least 15 successes in n=20 trials. 
cat("Prior Pr(X>=15|n=20) based on Beta-Binom(alpha, beta, n) \n")
cat(1-pbetabinomial(X=14,alpha=alpha.prior,beta=beta.prior,n=20),"\n")
# 0.01525992

#Q viii.  pictures of densities and obs'n, and likelihood- 10-15% overlap?
theta <- seq(0.01,0.99,by=0.01)
binomial.likelihood.norm.constant <- 
  integrate(f=function(theta) { dbinom(x=num.succ,size=n,prob=theta)},
            lower=0.01,upper=0.99)$value

prior.dens <- dbeta(x=theta,shape1=alpha.prior,shape2=beta.prior)
post.dens  <- dbeta(x=theta,shape1=alpha.post, shape2=beta.post)
likelihood <- dbinom(x=num.succ,size=n,prob=theta)/
  binomial.likelihood.norm.constant
my.ylim <- range(c(prior.dens,post.dens,likelihood))
plot(theta,prior.dens,type='l',col='blue',ylim=my.ylim,xlab='Pr(Success)',
     ylab='',main='Prior, Likelihood, Posterior for Pr(Success)')
lines(theta,likelihood,col='green',lty=2)
lines(theta,post.dens,col='red',lty=3)
legend('topleft',legend=c('Prior','Likelihood','Posterior'),
       col=c('blue','green','red'),lty=1:3 )

# ------------------------------------------------------------------------
# Question 2 - Mixture Prior: 95% for Be(9.2, 13.8) and 5% for Be such that E=0.8, SD=0.1

#Q i. parameters of Beta such that E=0.8 and SD=0.1
# solving for Beta shape parameters given mu and sigma
beta.param.calc <- function(mu,sigma) {
  alpha <- (mu^2*(1-mu)-mu*sigma^2)/sigma^2
  beta  <- alpha*(1-mu)/mu
  return(list(alpha=alpha,beta=beta))
}

winner.beta.par <- unlist(beta.param.calc(mu=0.8,sigma=0.1))
alpha.winner <- winner.beta.par[1]
beta.winner  <- winner.beta.par[2]
cat("Beta dist'n for winner, shape1=",alpha.winner,"shape2=",
    beta.winner ,"\n")
# Beta dist'n for winner, shape1= 12 shape2= 3 

#Q ii. Draw picture of the mixture prior- seems sensible
pi.1  <- 0.95
theta <- seq(0.01,0.99,by=0.01)
prior.dens.winner <- dbeta(x=theta,shape1=alpha.winner,shape2=beta.winner)
prior.mix.density <- pi.1*prior.dens+(1-pi.1)*prior.dens.winner
plot(theta,prior.mix.density,xlab="Prob. of Success",ylab="",
     main="Mixture Prior for Binomial p",type="l")

#Q iii. Pr(p>0.6)? Need posterior distribution of the mixture
y <- 15; n<- 20
marginal.y <- pi.1*dbetabinomial(x=y,alpha=alpha.prior,beta=beta.prior,n=n) +
  (1-pi.1)*dbetabinomial(x=y,alpha=alpha.winner,beta=beta.winner,n=n)

w1 <- pi.1*choose(n=n,k=y)*beta(alpha.prior+y,beta.prior+n-y)/
  (beta(alpha.prior,beta.prior))
w2 <- (1-pi.1)*choose(n=n,k=y)*beta(alpha.winner+y,beta.winner+n-y)/
  (beta(alpha.winner,beta.winner))
w1.scale <- w1/marginal.y
w2.scale <- w2/marginal.y
cat("wt1=",w1.scale,"wt2=",w2.scale,"sum=",w1.scale+w2.scale,"\n")


alpha.winner.post=alpha.winner+y
beta.winner.post=beta.winner+n-y

cat("Pr(p> .6)=",w1.scale*(1-pbeta(0.6,shape1=alpha.post,shape2=beta.post)) +
      w2.scale*(1-pbeta(0.6,shape1=alpha.winner.post,shape2=beta.winner.post)),
    "\n")
# Pr(p> .6)= 0.5806201, 

#Q iv. Bayesian P-value calculation, based on the Mixture Prior
# marginal dist'n is a mixture of Beta-Binom
#    Probability of at least 15 successes in n=20 trials. 
cat("Prior Pr(X>=15|n=20) based on mixture of Beta-Binom \n")
piece1 <- 1-pbetabinomial(X=14,alpha=alpha.prior,beta=beta.prior,n=20)
piece2 <- 1-pbetabinomial(X=14,alpha=alpha.winner,
                          beta=beta.winner,n=20)
cat(pi.1*piece1+(1-pi.1)*piece2,"\n")
# 0.05144791  # 5.1% vs 1.5%

#Q v. prior/likelihood/posterior plot
post.mix.density  <- w1.scale*dbeta(x=theta,shape1=alpha.post, shape2=beta.post) +
  w2.scale*dbeta(x=theta,shape1=alpha.winner.post,shape2=beta.winner.post)
likelihood <- dbinom(x=num.succ,size=n,prob=theta)/
  binomial.likelihood.norm.constant
my.ylim <- range(c(prior.mix.density,post.mix.density,likelihood))
plot(theta,prior.mix.density,type='l',col='blue',ylim=my.ylim,xlab='Pr(Success)',
     ylab='',main='Mixed Prior, Likelihood, Posterior for Pr(Success)')
lines(theta,likelihood,col='green',lty=2)
lines(theta,post.mix.density,col='red',lty=3)
legend('topleft',legend=c('Prior','Likelihood','Posterior'),
       col=c('blue','green','red'),lty=1:3 )

# ------------------------------------------------------------------------
# 3. Analysis of normal data: systolic blood pressure, sigma=5, mu unknown
obs <- c(127,133)

# Q i: prior for mu: Normal(mean=120,sd=10)

# Q ii: posterior given 2 measurements of 127 and 133
normal.mu.posterior.mean <- function(mu.prior,mu.sd,sigma,ybar,n) {
  wt.prior <- sigma^2/(sigma^2+n*mu.sd^2)
  out <- wt.prior*mu.prior + (1-wt.prior)*ybar
  return(out)
}

normal.mu.posterior.sd <- function(mu.sd,sigma,n) {
  out <- sqrt(mu.sd^2*sigma^2/(sigma^2+n*mu.sd^2))
  return(out)
}

post.mean <- normal.mu.posterior.mean(mu.prior=120,
                                      mu.sd=10,sigma=5,ybar=mean(obs),n=2)
post.sd   <- normal.mu.posterior.sd(mu.sd=10,sigma=5,n=2)
cat("posterior mean for mu=",post.mean ,"\n")
# posterior mean for mu= 128.8889 

cat("posterior sd for mu=",post.sd ,"\n")
# posterior sd for mu= 3.33 

# 95% credible interval
cat("95% CI for mu=",qnorm(p=c(0.025,0.975),mean=post.mean,sd=post.sd),"\n")
# 95% CI for mu= 122.3557 135.4221  

#Q iii: effect of 2 more readings of 130 on the posterior
obs <- c(127,133,130,130)
post.mean <- normal.mu.posterior.mean(mu.prior=120,
                                      mu.sd=10,sigma=5,ybar=mean(obs),n=4)
post.sd   <- normal.mu.posterior.sd(mu.sd=10,sigma=5,n=4)
cat("posterior mean for mu=",post.mean ,"\n")
# posterior mean for mu= 129.4118 

cat("posterior sd for mu=",post.sd ,"\n")
# posterior sd for mu= 2.425356  

# 95% credible interval
cat("95% CI for mu=",qnorm(p=c(0.025,0.975),mean=post.mean,sd=post.sd),"\n")
# 95% CI for mu= 124.6582 134.1654  

