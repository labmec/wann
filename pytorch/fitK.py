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
from knn import KNN
# Download helper functions (if it is not already there)
if Path("helper_functions.py").is_file():
  print("helper_functions.py already exists, skipping download")
else:
  print("Downloading helper_functions.py...")
  request = requests.get("https://raw.githubusercontent.com/mrdbourke/pytorch-deep-learning/main/helper_functions.py")
  with open("helper_functions.py","wb") as f:
    f.write(request.content)

from helper_functions import plot_predictions, plot_decision_boundary

torch.set_num_threads(8)
torch.set_num_interop_threads(8)
# Create directory
MODEL_PATH = Path("models")
# Create model save path
MODEL_NAME = "fitK-one-curve.pth"
saveModel = False
isUseBCasExtraCost = False

data = torch.from_numpy(np.loadtxt('../../kpts.txt', delimiter='\t',dtype=np.float32))
X = data[:, :-1]
y = data[:, -1].unsqueeze(dim=1)

if isUseBCasExtraCost:
  databc = torch.from_numpy(np.loadtxt('../../kpts_bc.txt', delimiter='\t',dtype=np.float32))
  Xbc = databc[:, :-1]
  ybc = databc[:, -1].unsqueeze(dim=1)

# Split data into train and test
X_train, X_test, y_train, y_test = train_test_split(X,y,
                                                    test_size=0.2, # 20% of data will be test
                                                    random_state=42)

X_train = X
y_train = y

t_boundaryA = torch.tensor(-1.).view(-1,1).requires_grad_(True)
t_boundaryB = torch.tensor(1.).view(-1,1).requires_grad_(True)

# Create model with input features
torch.manual_seed(42)
  
model0 = KNN(input_features=2,output_features=1,hidden_units=40,n_hiddenlayers=10)
# model0 = torch.nn.parallel.DistributedDataParallel(
#     model,
#     device_ids=[local_rank],
#     output_device=local_rank,
# )

# loss_fn = nn.L1Loss() # sigmoiod activation function built in. It uses log tricks to be more numerically stable
loss_fn = nn.MSELoss()
# loss_fn = nn.BCEWithLogitsLoss()

# optimizer = torch.optim.SGD(params=model0.parameters(),lr=0.002)
optimizer = torch.optim.AdamW(params=model0.parameters(),lr=0.00009)

# ---------------- Train the model -----------------
torch.manual_seed(42)
# epochs = 200
epochs = 15000
for epoch in range(epochs):
  model0.train() # sets requires_grad = True

  if epoch == 10000:
    for g in optimizer.param_groups:
      g['lr'] = 0.00006
  elif epoch == 15000:
    for g in optimizer.param_groups:
      g['lr'] = 0.00003
  elif epoch == 18000:
    for g in optimizer.param_groups:
      g['lr'] = 0.000015
  elif epoch == 22000:
    for g in optimizer.param_groups:
      g['lr'] = 0.000005   

  # 1. Forward pass
  y_pred = model0(X_train)

  # 2. Calculate loss/accuracy
  lossav = loss_fn(y_pred,
                 y_train) 
  # print(f"lossav: {lossav}")

  # Calculate the error in the first and last point
  if isUseBCasExtraCost:
    weight = 0.0025
    y_pred_bc = model0(Xbc)
    bc_error = 0
    for i in range(len(ybc)):
      bc_error += weight * torch.abs(y_pred_bc[i] - ybc[i])
    # bc_error = bc_error / len(ybc)
    loss = lossav + torch.squeeze(bc_error)
  else:
    loss = lossav

  # 3. Optimizer zero grad
  optimizer.zero_grad()

  # 4. Loss backward
  loss.backward()

  # 5. Optimizer step (gradient descent)
  optimizer.step()

  # 6. Testing
  model0.eval()
  with torch.inference_mode():
    test_pred = model0(X_test)
    test_loss = loss_fn(test_pred,y_test)


  # Print out what is happenning
  if epoch % 100 == 0:
    print(f"Epoch: {epoch} | Loss: {loss:.7f} | Test loss: {test_loss:.7f}")

# Lets see the results

# Plot decision boundary of the model
Xtensor = torch.tensor([[-1 + i * 0.0005, 0.0005] for i in range(4001)], dtype=torch.float32)
test_pred = model0(Xtensor)
plot_predictions(train_data=X_train[:, 0].detach().numpy(),
                train_labels=y_train.detach().numpy(),
                test_data=Xtensor[:, 0].detach().numpy(),
                test_labels=test_pred.detach().numpy(),
                predictions=test_pred.detach().numpy(),
                all_data=Xtensor[:, 0].detach().numpy())
plt.show()

# Save model
MODEL_SAVE_PATH = MODEL_PATH / MODEL_NAME
if saveModel:
  MODEL_PATH.mkdir(parents=True,exist_ok=True)  
  print(f"Saving model to: {MODEL_SAVE_PATH}")
  torch.save(obj=model0.state_dict(),f=MODEL_SAVE_PATH)