import torch
from torch import nn
import matplotlib.pyplot as plt
import sklearn
from sklearn.model_selection import train_test_split
import requests
from pathlib import Path
import sys
import numpy as np
# Download helper functions (if it is not already there)
if Path("helper_functions.py").is_file():
  print("helper_functions.py already exists, skipping download")
else:
  print("Downloading helper_functions.py...")
  request = requests.get("https://raw.githubusercontent.com/mrdbourke/pytorch-deep-learning/main/helper_functions.py")
  with open("helper_functions.py","wb") as f:
    f.write(request.content)

from helper_functions import plot_predictions, plot_decision_boundary

weight = 0.7
bias = 0.3
start = 0
end = 1
step = 0.1

# X = torch.arange(start,end,step).unsqueeze(dim=1)
# y = weight * X + bias


data = torch.from_numpy(np.loadtxt('kpts.txt', delimiter='\t',dtype=np.float32))
X = data[:, :-1]
y = data[:, -1].unsqueeze(dim=1)

# Split data into train and test
X_train, X_test, y_train, y_test = train_test_split(X,y,
                                                    test_size=0.01, # 20% of data will be test
                                                    random_state=42)

# train_split = int(0.8 * len(X))
# X_train, y_train = X[:train_split], y[:train_split]
# X_test, y_test = X[train_split:], y[train_split:]

# plot_predictions(train_data=X_train,train_labels=y_train,test_data=X_test,test_labels=y_test)
# plt.show()

# Construct model
torch.manual_seed(42)
model0 = nn.Sequential(
  nn.Linear(in_features=1,out_features=32),
  nn.Tanh(),
  nn.Linear(in_features=32,out_features=32),
  nn.Tanh(),
  nn.Linear(in_features=32,out_features=32),
  nn.Tanh(), 
  nn.Linear(in_features=32,out_features=1),
  nn.Tanh(),
)

loss_fn = nn.L1Loss() # sigmoiod activation function built in. It uses log tricks to be more numerically stable
# loss_fn = nn.MSELoss()

optimizer = torch.optim.SGD(params=model0.parameters(),lr=0.001)

# ---------------- Train the model -----------------
torch.manual_seed(42)
epochs = 15000
for epoch in range(epochs):
  model0.train() # sets requires_grad = True

  # 1. Forward pass
  y_pred = model0(X_train)

  # 2. Calculate loss/accuracy
  loss = loss_fn(y_pred, # We are passing logits bcz BCEWithLogits expects logits and does sigmoid by itself
                 y_train) 

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
    print(f"Epoch: {epoch} | Loss: {loss:.5f} | Test loss: {test_loss:.5f}")

# Lets see the results


# Plot decision boundary of the model
test_pred = model0(X)
plot_predictions(train_data=X_train.detach().numpy(),
                 train_labels=y_train.detach().numpy(),
                 test_data=X_test.detach().numpy(),
                 test_labels=y_test.detach().numpy(),
                 predictions=test_pred.detach().numpy(),
                 all_data=X.detach().numpy())
plt.show()