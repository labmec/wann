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

x = np.linspace(-1, 1, 100)  # Generate x values
mydeg = 4
# Compute Legendre polynomials up to order 6
# legendre_polynomials = [np.polynomial.legendre.Legendre.basis(deg=i)(x) for i in range(7)]

# Plot Legendre polynomials
# for i, polynomial in enumerate(legendre_polynomials):
#     plt.plot(x, polynomial, label=f'Order {i}')

# plt.legend()  # Show legend
# plt.xlabel('x')
# plt.ylabel('P(x)')
# plt.title('Legendre Polynomials')
# # plt.show()

# Fit Legendre polynomials to sin(x)

sinx = np.sin(x)
coefficients = np.polynomial.legendre.legfit(x, sinx, deg=mydeg)
print(f"coefficients: {coefficients}")
# Compute the Legendre polynomials with the fitted coefficients
fitted_polynomials = [coef * np.polynomial.legendre.Legendre.basis(deg=i)(x) for i,coef in enumerate(coefficients)]
fitted_polynomial = np.sum(fitted_polynomials, axis=0)
# Plot fitted_polynomial and sin(x)
plt.plot(x, sinx, label='sin(x)', linewidth=4.0)

plt.plot(x, fitted_polynomial, label='Fitted polynomial')
plt.legend()
# plt.show()

sinx = np.sin(x)
coefficients = np.polyfit(x, sinx, deg=mydeg)
print(f"coefficients: {coefficients}")
# Compute the polynomials with the fitted coefficients

fitted_polynomial = np.poly1d(coefficients)
fitted_polynomial = fitted_polynomial(x)
# Plot fitted_polynomial and sin(x)
# plt.plot(x, sinx, label='sin(x)', linewidth=4.0)

plt.plot(x, fitted_polynomial, label='Fitted polynomial polyfit', linestyle='--')
plt.legend()
plt.show()
