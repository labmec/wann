import torch
from torch import nn

class KNN(nn.Module):
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
    activation = nn.Tanh
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