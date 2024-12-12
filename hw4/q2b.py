import numpy as np
import pandas as pd
import os
import torch
import torch.nn as nn
import torch.optim as optim
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.model_selection import train_test_split

"""
MAX_JOBS=4
for i in {1..10}; do
  echo "Starting simulation iteration $i"
  python simulation.py
  
  # Limit the number of background jobs
  if (( $(jobs -r | wc -l) >= MAX_JOBS )); then
    wait -n
  fi
done

# Wait for all remaining background jobs to finish
wait
echo "All simulations completed."
"""


def mu0_function(X, s):
    if s == 0:
        return np.zeros(X.shape[0])
    return np.sum(np.sin(X[:, :s]), axis=1)

def g_function(X, s):
    if s == 0:
        return np.zeros(X.shape[0])
    g = np.zeros(X.shape[0])
    if s >= 1:
        g += 0.5 * np.sin(X[:, 0])
    if s >= 2:
        g += (X[:, 1] ** 2) / 4
    if s >= 3:
        g -= np.cos(X[:, 2])
    for j in range(3, s):
        g += 0.3 * (X[:, j] ** 3)
    return g

# Generate synthetic data
def generate_data(n, d, s, tau=1, seed=None):
    if seed is not None:
        np.random.seed(seed)
    X = np.random.normal(0, 1, (n, d))
    mu0_vals = mu0_function(X, s)
    g_vals = g_function(X, s)
    p_x = 1 / (1 + np.exp(-g_vals))
    T = np.random.binomial(1, p_x)
    Y = T * (mu0_vals + tau) + (1 - T) * mu0_vals + np.random.normal(0, 1, n)
    return X, Y, T, p_x, mu0_vals

# Define PyTorch models
class MLPRegressor(nn.Module):
    def __init__(self, input_dim):
        super(MLPRegressor, self).__init__()
        self.model = nn.Sequential(
            nn.Linear(input_dim, 64),
            nn.ReLU(),
            nn.Linear(64, 32),
            nn.ReLU(),
            nn.Linear(32, 1)
        )

    def forward(self, x):
        return self.model(x)

class MLPClassifier(nn.Module):
    def __init__(self, input_dim):
        super(MLPClassifier, self).__init__()
        self.model = nn.Sequential(
            nn.Linear(input_dim, 64),
            nn.ReLU(),
            nn.Linear(64, 32),
            nn.ReLU(),
            nn.Linear(32, 1),
            nn.Sigmoid()
        )

    def forward(self, x):
        return self.model(x)

# Train a PyTorch model
def train_pytorch_model(model, X, y, epochs=10, batch_size=32):
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model.to(device)
    criterion = nn.MSELoss() if isinstance(model, MLPRegressor) else nn.BCELoss()
    optimizer = optim.Adam(model.parameters(), lr=0.001)
    X_tensor = torch.tensor(X, dtype=torch.float32).to(device)
    y_tensor = torch.tensor(y, dtype=torch.float32).to(device)

    dataset = torch.utils.data.TensorDataset(X_tensor, y_tensor)
    loader = torch.utils.data.DataLoader(dataset, batch_size=batch_size, shuffle=True)

    model.train()
    for epoch in range(epochs):
        for batch_X, batch_y in loader:
            optimizer.zero_grad()
            outputs = model(batch_X).squeeze()
            loss = criterion(outputs, batch_y)
            loss.backward()
            optimizer.step()

    model.eval()
    return model

# Fit propensity score and mu0 using PyTorch models
def fit_p_deepnet(X, T, epochs=10, batch_size=32):
    model = MLPClassifier(X.shape[1])
    trained_model = train_pytorch_model(model, X, T, epochs, batch_size)
    return lambda newX: trained_model(torch.tensor(newX, dtype=torch.float32)).detach().numpy().squeeze()

def fit_mu0_deepnet(X, Y, T, epochs=10, batch_size=32):
    X0 = X[T == 0]
    Y0 = Y[T == 0]
    model = MLPRegressor(X.shape[1])
    trained_model = train_pytorch_model(model, X0, Y0, epochs, batch_size)
    return lambda newX: trained_model(torch.tensor(newX, dtype=torch.float32)).detach().numpy().squeeze()

# Plugin estimator
def compute_plugin_estimator(Y, T, p_hat_func, X):
    p_hat_vals = p_hat_func(X)
    p_hat_vals = np.clip(p_hat_vals, 1e-3, 1 - 1e-3)
    p_bar = np.mean(p_hat_vals)
    mu1_hat = np.mean((T * Y) / p_bar)
    mu0_hat = np.mean(((1 - T) * p_hat_vals * Y) / (1 - p_hat_vals))
    return mu1_hat - mu0_hat

# DR estimator
def compute_DR_estimator(Y, T, p_hat_func, mu0_hat_func, X):
    p_hat_vals = p_hat_func(X)
    p_hat_vals = np.clip(p_hat_vals, 1e-3, 1 - 1e-3)
    p_bar = np.mean(p_hat_vals)
    mu0_hat_vals = mu0_hat_func(X)
    mu1_hat = np.mean((T * Y) / p_bar)
    term = T * mu0_hat_vals + ((1 - T) * p_hat_vals * Y) / (1 - p_hat_vals)
    mu0_hat = np.mean(term) / p_bar
    return mu1_hat - mu0_hat

# Simulation setup
def run_simulation(ns, ds, sparsities, methods, reps=50, tau=1, seed=123):
    np.random.seed(seed)
    results = []

    for n in ns:
        for d in ds:
            for s in [sp for sp in sparsities if sp <= d]:
                for method in methods:
                    tau_pi_vals = []
                    tau_dr_vals = []

                    for _ in range(reps):
                        X, Y, T, _, _ = generate_data(n, d, s, tau)
                        if method == "rf":
                            p_model = RandomForestClassifier().fit(X, T)
                            p_hat_func = lambda newX: p_model.predict_proba(newX)[:, 1]
                            mu0_model = RandomForestRegressor().fit(X[T == 0], Y[T == 0])
                            mu0_hat_func = lambda newX: mu0_model.predict(newX)
                        elif method == "deepnet":
                            p_hat_func = fit_p_deepnet(X, T)
                            mu0_hat_func = fit_mu0_deepnet(X, Y, T)

                        tau_pi_vals.append(compute_plugin_estimator(Y, T, p_hat_func, X))
                        tau_dr_vals.append(compute_DR_estimator(Y, T, p_hat_func, mu0_hat_func, X))

                    new_result = {
                        'n': n,
                        'd': d,
                        's': s,
                        'method': method,
                        'plugin_mean': np.mean(tau_pi_vals),
                        'plugin_sd': np.std(tau_pi_vals),
                        'dr_mean': np.mean(tau_dr_vals),
                        'dr_sd': np.std(tau_dr_vals)
                    }
                    print(new_result)

                    results.append(new_result)

    return results

seed_opt = np.random.randint(0, 1000)

results = run_simulation(
    ns=[1000, 5000, 15000],
    ds=[1, 3, 5, 10],
    sparsities=[1, 3, 5, 10],
    methods=["rf", "deepnet"],
    reps=1,
    tau=1,
    seed=seed_opt
)

results_df = pd.DataFrame(results)
timestamp = pd.Timestamp.now().strftime("%Y%m%d_%H%M%S")
# if "hw4" in the path name, the file will be saved in the hw4 folder
try:
    if "hw4" in os.getcwd():
        results_df.to_csv(f"q2b-data/q2b_all_simulation_results_{timestamp}_{seed_opt:.0f}.csv", index=False)
    else:
        results_df.to_csv(f"hw4/q2b-data/q2b_all_simulation_results_{timestamp}_{seed_opt:.0f}.csv", index=False)
except Exception:
    results_df.to_csv(f"q2b_all_simulation_results_{timestamp}_{seed_opt:.0f}.csv", index=False)