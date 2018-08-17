from numpy import interp,load # directly takes a query point
from scipy.interpolate import pchip # creates an interpolation function that can then be given a query point
from scipy.optimize import fminbound

npzfile = load('BTM.npz')
beta = npzfile['beta']; theta = npzfile['theta']; M2 = npzfile['M2'];

inds = M2.argsort()
M2 = M2[inds]; beta = beta[inds]; theta = theta[inds];
beta_star = interp(1.0,M2,beta)

inds = beta.argsort()
M2 = M2[inds]; beta = beta[inds]; theta = theta[inds];
theta_star = interp(beta_star,beta,theta)

# find maximum deflection estimate
theta_max_guess = max(theta)

inds = theta.argsort()
M2 = M2[inds]; beta = beta[inds]; theta = theta[inds];
beta_max_guess = interp(theta_max_guess,theta,beta)


inds = beta.argsort()
M2 = M2[inds]; beta = beta[inds]; theta = theta[inds];
fun = pchip(beta,theta)
beta_max = fminbound(lambda x: -fun(x), .9*beta_max_guess, 1.1*beta_max_guess)
theta_max = fun(beta_max)


