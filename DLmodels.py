# Import necessary libraries
import os
import itertools
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, roc_auc_score
from sklearn.preprocessing import StandardScaler
import torch
import torch.nn as nn
from torch.utils.data import DataLoader, Dataset
from sklearn.feature_selection import RFE
from sklearn.ensemble import RandomForestClassifier
from tqdm import tqdm  # For progress bar
from sklearn.model_selection import train_test_split, StratifiedKFold, cross_val_score
from sklearn.utils import resample
from imblearn.over_sampling import SMOTE
from sklearn.metrics import precision_recall_curve

# Set random seed for reproducibility
RANDOM_SEED = 42
np.random.seed(RANDOM_SEED)

# Define Transformer Model
class TransformerModel(nn.Module):
    def __init__(self, input_dim, output_dim, nhead=8, num_layers=4, hidden_dim=256):
        super(TransformerModel, self).__init__()
        self.embedding = nn.Linear(input_dim, hidden_dim)
        self.dropout = nn.Dropout(0.2)
        self.layer_norm = nn.LayerNorm(hidden_dim)
        self.transformer = nn.Transformer(hidden_dim, nhead, num_layers, batch_first=True)
        self.fc = nn.Linear(hidden_dim, output_dim)

    def forward(self, x):
        x = self.embedding(x)
        x = self.dropout(x)
        x = self.layer_norm(x)
        x = x.unsqueeze(1)
        x = self.transformer(x, x)
        x = x.mean(dim=1)
        x = self.fc(x)
        return torch.sigmoid(x)

# Define MLP Model
class MLPModel(nn.Module):
    def __init__(self, input_dim, output_dim):
        super(MLPModel, self).__init__()
        self.fc1 = nn.Linear(input_dim, 256)
        self.fc2 = nn.Linear(256, 128)
        self.fc3 = nn.Linear(128, 64)
        self.fc4 = nn.Linear(64, output_dim)
        self.dropout = nn.Dropout(0.2)

    def forward(self, x):
        x = torch.relu(self.fc1(x))
        x = self.dropout(x)
        x = torch.relu(self.fc2(x))
        x = self.dropout(x)
        x = torch.relu(self.fc3(x))
        x = self.fc4(x)
        return torch.sigmoid(x)

# Define RNN Model
class RNNModel(nn.Module):
    def __init__(self, input_dim, output_dim):
        super(RNNModel, self).__init__()
        self.rnn = nn.LSTM(input_dim, 256, num_layers=3, batch_first=True, dropout=0.3, bidirectional=True)
        self.fc = nn.Linear(256 * 2, output_dim)

    def forward(self, x):
        x = x.unsqueeze(1)
        _, (hidden, _) = self.rnn(x)
        x = torch.cat((hidden[-2], hidden[-1]), dim=1)
        x = self.fc(x)
        return torch.sigmoid(x)

# Define CNN Model
class CNNModel(nn.Module):
    def __init__(self, input_dim, output_dim):
        super(CNNModel, self).__init__()
        self.conv1 = nn.Conv1d(in_channels=1, out_channels=128, kernel_size=5, padding=2)
        self.conv2 = nn.Conv1d(in_channels=128, out_channels=256, kernel_size=3, padding=1)
        self.pool = nn.MaxPool1d(2)  # Default stride is 2
        self.dropout = nn.Dropout(0.3)
        self.bn1 = nn.BatchNorm1d(128)
        self.bn2 = nn.BatchNorm1d(256)
        # Dynamically adjust the input features of the linear layer
        self.fc = nn.Linear(self._calculate_fc_input(input_dim), output_dim)

    def _calculate_fc_input(self, input_dim):
        """
        Calculate the input features for the fully connected layer.
        Dynamically adjust pooling to ensure no invalid sizes.
        """
        x = torch.zeros((1, 1, input_dim))  # Assume one sample
        x = self.pool(self.bn1(torch.relu(self.conv1(x))))
        if x.shape[2] < 2:  # Skip second pooling if size is too small
            return x.size(1) * x.size(2)
        x = self.pool(self.bn2(torch.relu(self.conv2(x))))
        return x.size(1) * x.size(2)

    def forward(self, x):
        x = x.unsqueeze(1)  # Add channel dimension (batch_size, 1, input_dim)
        x = self.pool(self.bn1(torch.relu(self.conv1(x))))
        if x.shape[2] >= 2:  # Only apply second pooling if feature dimension is large enough
            x = self.pool(self.bn2(torch.relu(self.conv2(x))))
        x = x.view(x.size(0), -1)  # Flatten to (batch_size, num_filters * feature_dim)
        x = self.fc(x)
        return torch.sigmoid(x)

# Custom Dataset Class
class GeneDataset(Dataset):
    def __init__(self, X, y):
        self.X = X
        self.y = y

    def __len__(self):
        return len(self.y)

    def __getitem__(self, idx):
        return torch.tensor(self.X[idx], dtype=torch.float32), torch.tensor(self.y[idx], dtype=torch.float32)

# Custom output file prefix
output_prefix = "GSE51799-HDHC"

# Read all CSV files in the folder
folder_path = "/home/xxu/AD/PD/GEOdata/deg0.1/module/GSE51799-HDHC/"
files = [os.path.join(folder_path, f) for f in os.listdir(folder_path) if f.endswith('.csv')]

# Load and combine data
dataframes = []
for idx, file in enumerate(files):
    df = pd.read_csv(file)
    df['Sample_ID'] = range(idx * len(df), (idx + 1) * len(df))  # Add unique row IDs for each matrix
    dataframes.append(df)

file_names = [os.path.basename(file) for file in files]
# Generate all possible combinations of dataframes
combinations = []
for r in range(1, len(dataframes) + 1):
    combinations += list(itertools.combinations(range(len(dataframes)), r))

print("\n### Prediction Combinations ###")
for idx, comb in enumerate(combinations):
    print(f"Combination {idx + 1}: {[file_names[i] for i in comb]}")

# Store results
all_results = []

# Iterate over each combination
for comb in tqdm(combinations, desc="Processing combinations"):
    # Merge data
    combined_df = pd.concat([dataframes[i] for i in comb], axis=1)
    combined_df = combined_df.loc[:, ~combined_df.columns.duplicated()]  # Remove duplicate columns
    print(f"Combined DataFrame shape: {combined_df.shape}")
    print(f"Combined DataFrame columns: {combined_df.columns}")
    
    # Extract target and features
    y = combined_df.iloc[:, 0].values
    X = combined_df.iloc[:, 1:-1].values  # Remove the last Sample_ID column
    print(f"X shape before removing duplicates: {X.shape}")

    # Remove duplicate rows
    df_deduplicated = pd.DataFrame(X).drop_duplicates()
    X = df_deduplicated.values
    y = y[:len(X)]
    print(f"X shape after removing duplicates: {X.shape}")

    # Check if any features are lost
    if X.shape[1] == 0:
        raise ValueError("Feature matrix X has 0 features. Check your data processing pipeline.")

    # Check class distribution
    unique, counts = np.unique(y, return_counts=True)
    print(f"Class distribution: {dict(zip(unique, counts))}")

    # Stratified train-test split
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=RANDOM_SEED, stratify=y)
    print(f"X_train shape: {X_train.shape}, y_train shape: {y_train.shape}")

    # Apply SMOTE for oversampling
    if X_train.shape[1] > 0:
        smote = SMOTE(random_state=RANDOM_SEED)
        X_train, y_train = smote.fit_resample(X_train, y_train)
        print(f"X_train shape after SMOTE: {X_train.shape}, y_train shape: {y_train.shape}")
    else:
        raise ValueError("SMOTE cannot be applied because X_train has 0 features.")

    # Standardize data
    scaler = StandardScaler()
    X_train = scaler.fit_transform(X_train)
    X_test = scaler.transform(X_test)

    # Convert to PyTorch datasets
    train_dataset = GeneDataset(X_train, y_train)
    test_dataset = GeneDataset(X_test, y_test)
    train_loader = DataLoader(train_dataset, batch_size=32, shuffle=True)
    test_loader = DataLoader(test_dataset, batch_size=32, shuffle=False)

    # Check class distribution
    unique_train, counts_train = np.unique(y_train, return_counts=True)
    unique_test, counts_test = np.unique(y_test, return_counts=True)
    print(f"Train class distribution: {dict(zip(unique_train, counts_train))}")
    print(f"Test class distribution: {dict(zip(unique_test, counts_test))}")

    # List of models
    models = {
        'Transformer': TransformerModel(X_train.shape[1], 1),
        'MLP': MLPModel(X_train.shape[1], 1),
        'RNN': RNNModel(X_train.shape[1], 1),
        'CNN': CNNModel(X_train.shape[1], 1)
    }

    # Record performance for each model
    for model_name, model in models.items():
        print(f"\n### Training {model_name} on combination: {[file_names[i] for i in comb]} ###")
        model = model.to('cuda' if torch.cuda.is_available() else 'cpu')
        criterion = nn.BCELoss()
        optimizer = torch.optim.Adam(model.parameters(), lr=0.0001)
        scheduler = torch.optim.lr_scheduler.StepLR(optimizer, step_size=5, gamma=0.5)

        accuracy_list, precision_list, recall_list, f1_list, roc_auc_list = [], [], [], [], []

        best_roc_auc = 0
        patience, max_patience = 0, 3

        for epoch in range(20):  # Increase number of epochs
            model.train()
            train_loss = 0
            for batch_idx, (batch_X, batch_y) in enumerate(train_loader):
                batch_X, batch_y = batch_X.to('cuda'), batch_y.to('cuda')
                optimizer.zero_grad()
                outputs = model(batch_X).squeeze()
                loss = criterion(outputs, batch_y)
                loss.backward()
                optimizer.step()
                train_loss += loss.item()

                if (batch_idx + 1) % 10 == 0:
                    print(f"Epoch {epoch + 1}, Batch {batch_idx + 1}, Loss: {loss.item():.4f}")
            
            scheduler.step()

            model.eval()
            y_true, y_pred, y_scores = [], [], []
            with torch.no_grad():
                for batch_X, batch_y in test_loader:
                    batch_X, batch_y = batch_X.to('cuda'), batch_y.to('cuda')
                    outputs = model(batch_X).squeeze()
                    y_true.extend(batch_y.cpu().numpy())
                    y_scores.extend(outputs.cpu().numpy())

                # Dynamic threshold selection
                precisions, recalls, thresholds = precision_recall_curve(y_true, y_scores)
                f1_scores = 2 * (precisions * recalls) / (precisions + recalls + 1e-6)
                best_threshold = thresholds[np.argmax(f1_scores)]
                y_pred = (np.array(y_scores) > best_threshold).astype(int)

            roc_auc = roc_auc_score(y_true, y_scores)
            accuracy = accuracy_score(y_true, y_pred)
            precision = precision_score(y_true, y_pred, zero_division=1)
            recall = recall_score(y_true, y_pred, zero_division=1)
            f1 = f1_score(y_true, y_pred, zero_division=1)

            print(f"Epoch {epoch + 1}, ROC AUC: {roc_auc:.4f}, Accuracy: {accuracy:.4f}, Precision: {precision:.4f}, Recall: {recall:.4f}, F1: {f1:.4f}")

            accuracy_list.append(accuracy)
            precision_list.append(precision)
            recall_list.append(recall)
            f1_list.append(f1)
            roc_auc_list.append(roc_auc)

            if roc_auc > best_roc_auc:
                best_roc_auc = roc_auc
                patience = 0
            else:
                patience += 1

            if patience > max_patience:
                print(f"Early stopping at epoch {epoch} for {model_name}")
                break

        # Calculate mean performance metrics
        accuracy_mean = np.mean(accuracy_list)
        precision_mean = np.mean(precision_list)
        recall_mean = np.mean(recall_list)
        f1_mean = np.mean(f1_list)
        roc_auc_mean = np.mean(roc_auc_list)

        # Cross-validation evaluation
        skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
        cross_val_auc = []
        for train_idx, val_idx in skf.split(X, y):
            X_train_cv, X_val_cv = X[train_idx], X[val_idx]
            y_train_cv, y_val_cv = y[train_idx], y[val_idx]
            rf = RandomForestClassifier(n_estimators=100, max_depth=10, random_state=42)
            rf.fit(X_train_cv, y_train_cv)
            y_val_scores = rf.predict_proba(X_val_cv)[:, 1]
            cross_val_auc.append(roc_auc_score(y_val_cv, y_val_scores))
        print(f"Cross-validated AUC for {model_name}: {np.mean(cross_val_auc):.4f}")

        # Save mean performance metrics
        metrics = {
            'Combination': [file_names[i] for i in comb],
            'Model': model_name,
            'Accuracy (Mean)': accuracy_mean,
            'Precision (Mean)': precision_mean,
            'Recall (Mean)': recall_mean,
            'F1 Score (Mean)': f1_mean,
            'ROC AUC (Mean)': roc_auc_mean,
            'Cross-validated AUC': np.mean(cross_val_auc),
            'Final Score (Mean of AUCs)': (best_roc_auc + np.mean(cross_val_auc)) / 2
        }

        # Save results
        all_results.append(metrics)

        print(f"\n### Results for {model_name} ###")
        print(metrics)

# Save all results to a CSV file
results_df = pd.DataFrame(all_results)
results_df = results_df.sort_values(by='Final Score (Mean of AUCs)', ascending=False)
results_df.to_csv(f"{output_prefix}_performance.csv", index=False)

# Extract the best combination for each model
best_combinations = results_df.groupby('Model').first().reset_index()

# Save gene importance
for _, row in best_combinations.iterrows():
    model_name = row['Model']
    best_comb = row['Combination']
    comb_indices = [file_names.index(file) for file in best_comb]
    comb_df = pd.concat([dataframes[i] for i in comb_indices], axis=1)
    comb_df = comb_df.loc[:, ~comb_df.columns.duplicated()]

    # Extract gene names (starting from the second column)
    gene_names = comb_df.columns[1:]  # Exclude the first target column
    if len(gene_names) == 0:
        print(f"No genes found for combination: {best_comb}")
        continue

    # Extract module names
    module_names = []
    for i in comb_indices:
        module_name = file_names[i]
        module_gene_count = len(dataframes[i].columns[1:])  # Exclude the target column
        module_names.extend([module_name] * module_gene_count)
    module_names = module_names[:len(gene_names)]  # Truncate or extend to match gene count

    # Train-test split
    X = comb_df.iloc[:, 1:].values
    y = comb_df.iloc[:, 0].values
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=RANDOM_SEED)

    # Use Random Forest as base model for RFE
    rf = RandomForestClassifier(n_estimators=100, random_state=42)
    rfe = RFE(estimator=rf, n_features_to_select=len(gene_names))  # Select all features
    rfe.fit(X_train, y_train)

    # Get selected features
    selected_features = np.where(rfe.support_)[0]
    selected_gene_names = [gene_names[i] for i in selected_features]

    # If RFE selects all features, use all gene names; otherwise, use selected gene names
    gene_names = selected_gene_names if len(selected_gene_names) != 0 else gene_names

    # Explicitly fit Random Forest to get feature importance
    rf.fit(X_train[:, selected_features], y_train)  # Fit on selected features
    importances = rf.feature_importances_  # Get feature importance

    # Normalize importance
    normalized_importances = (importances - importances.min()) / (importances.max() - importances.min())
    
    # Save gene, module, and importance
    importance_data = pd.DataFrame({
        'Gene': gene_names,
        'Module': module_names[:len(gene_names)],  # Ensure module names match gene names
        'Importance': normalized_importances
    })

    # Write to CSV file
    importance_file = f"{output_prefix}_{model_name}_best_combination_genes.csv"
    importance_data.to_csv(importance_file, index=False)
    print(f"Saved importance data for {model_name} to {importance_file}")

print("\n### Final Results Summary ###")
print(results_df)