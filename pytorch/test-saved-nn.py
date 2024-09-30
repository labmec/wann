import torch
from pathlib import Path
import matplotlib.pyplot as plt
from knn import KNN
 
MODEL_PATH = Path("models")
MODEL_NAME = "fitK-one-curve.pth"  
MODEL_SAVE_PATH = MODEL_PATH / MODEL_NAME  
model = KNN(input_features=1,output_features=1,hidden_units=40,n_hiddenlayers=10) # These numbers have to match the model that was saved
model.load_state_dict(torch.load(f=MODEL_SAVE_PATH))

# ---------- Testing the model ----------
isTestModel = True

if isTestModel:
    X_test = torch.linspace(-1,1,200).view(-1,1)
    y_test = model(X_test)

    # Plot the X_test and y_test using matplotlib
    plt.plot(X_test.numpy(), y_test.detach().numpy())
    plt.xlabel('Position')
    plt.ylabel('K')
    plt.grid(True)
    plt.show()


