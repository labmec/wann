import torch
from torch import nn
import matplotlib.pyplot as plt
import sklearn
from sklearn.model_selection import train_test_split
import requests
from pathlib import Path
import sys
import numpy as np
import torch.distributed as dist
import os

data = np.loadtxt('../../kpts.txt', delimiter='\t',dtype=np.float32)
print(f"data: {data}")
x = data[:,0]
y = data[:,1]
print(f"x: {x}")
print(f"y: {y}")
mydeg = 30

coefficients = np.polynomial.legendre.legfit(x, y, deg=mydeg)
print(f"coefficients: {coefficients}")
# Compute the Legendre polynomials with the fitted coefficients
fitted_polynomials = [coef * np.polynomial.legendre.Legendre.basis(deg=i)(x) for i,coef in enumerate(coefficients)]
fitted_polynomial = np.sum(fitted_polynomials, axis=0)
# Plot fitted_polynomial and sin(x)
plt.plot(x, y, 'o', label='Kreal')

plt.plot(x, fitted_polynomial, label='Fitted polynomial Legendre')
plt.legend()
# plt.show()

# coefficients = np.polyfit(x, y, deg=mydeg)
# print(f"coefficients: {coefficients}")
# # Compute the polynomials with the fitted coefficients

# fitted_polynomial = np.poly1d(coefficients)
# fitted_polynomial = fitted_polynomial(x)
# Plot fitted_polynomial and sin(x)
# plt.plot(x, sinx, label='sin(x)', linewidth=4.0)

# plt.plot(x, fitted_polynomial, label='Fitted polynomial polyfit', linestyle='--')
plt.legend()
plt.show()
