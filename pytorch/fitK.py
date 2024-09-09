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

data = torch.from_numpy(np.loadtxt('kpts.txt', delimiter='\t',dtype=np.float32))
X = data[:, :-1]
y = data[:, -1].unsqueeze(dim=1)

# Split data into train and test
X_train, X_test, y_train, y_test = train_test_split(X,y,
                                                    test_size=0.01, # 20% of data will be test
                                                    random_state=42)


# Create model with input features
torch.manual_seed(42)
class MyModel(nn.Module):
  def __init__(self,input_features: int,output_features: int,hidden_units=8,n_hiddenlayers=1):    
    """Initializes multi-class classification model.

    Args:
        input_features (int): Number of inputs to the model
        output_features (int): 
        hidden_units (int, optional): Number of hidden units between layers. Defaults to 8.
    """
    super().__init__()
    if n_hiddenlayers < 1:
      raise ValueError("n_hiddenlayers must be at least 1")
    activation = nn.ReLU
    self.entry =  nn.Sequential(*[
                        nn.Linear(in_features=input_features, out_features=hidden_units),
                        activation()])
    self.hidden = nn.Sequential(*[
                        nn.Sequential(*[
                            nn.Linear(in_features=hidden_units, out_features=hidden_units),
                            activation()]) for _ in range(n_hiddenlayers)])
    self.end = nn.Linear(in_features=hidden_units, out_features=output_features)


  def forward(self,x):
    x = self.entry(x)
    x = self.hidden(x)
    x = self.end(x)
    return x

  
model0 = MyModel(input_features=1,output_features=1,hidden_units=16,n_hiddenlayers=2)

loss_fn = nn.L1Loss() # sigmoiod activation function built in. It uses log tricks to be more numerically stable
# loss_fn = nn.MSELoss()
# loss_fn = nn.BCEWithLogitsLoss()

optimizer = torch.optim.SGD(params=model0.parameters(),lr=0.002)

# ---------------- Train the model -----------------
torch.manual_seed(42)
epochs = 20000
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