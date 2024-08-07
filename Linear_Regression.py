import pandas as pd
from sklearn.linear_model import LinearRegression

# Load the data from an Excel file
file_path = 'data.xlsx'  # Replace with your file path
data = pd.read_excel(file_path)

# Extract the shift (x) and tensor (y) values
X = data[['shift']].values  # Feature matrix
y = data['tensor'].values   # Target vector

# Initialize and train the linear regression model
model = LinearRegression()
model.fit(X, y)

# Get the slope (coefficient) and intercept of the model
slope = model.coef_[0]
intercept = model.intercept_

print(f'Slope: {slope}')
print(f'Intercept: {intercept}')

# Create a DataFrame to save the slope and intercept
results = pd.DataFrame({
    'Parameter': ['Slope', 'Intercept'],
    'Value': [slope, intercept]
})

# Load existing Excel file and write the new data to a new sheet
with pd.ExcelWriter(file_path, mode='a', engine='openpyxl', if_sheet_exists='new') as writer:
    results.to_excel(writer, sheet_name='Regression_Results', index=False)

