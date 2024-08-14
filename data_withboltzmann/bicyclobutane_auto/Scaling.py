import numpy as np
import pandas as pd
import argparse

def update_chemical_shifts_from_excel(file_path, slope, intercept):
    """
    Update chemical shifts using a linear regression model from tensor values in an Excel sheet.

    Args:
        file_path (str): Path to the Excel file containing tensor values.
        slope (float): Slope of the linear regression model.
        intercept (float): Intercept of the linear regression model.

    Returns:
        dict: Updated chemical shifts.
    """
    # Read the Excel file
    df = pd.read_excel(file_path)

    # Ensure the required columns are present
    if 'Atoms' not in df.columns or 'Average Tensor' not in df.columns:
        raise ValueError("The Excel file must contain 'Atoms' and 'Average Tensor' columns.")

    # Convert the DataFrame to a dictionary
    tensor_values = dict(zip(df['Atoms'], df['Average Tensor']))

    # Calculate updated chemical shifts
    updated_shifts = {}
    for atom, tensor in tensor_values.items():
        updated_shifts[atom] = (tensor - intercept) / slope

    return updated_shifts
  
def main():
 
    parser = argparse.ArgumentParser(description="Scale the Shielding Tensor Values.")
    parser.add_argument("--slope", type=float, help="Slope value of linear regression model")
    parser.add_argument("--intercept", type=float, help="Intercept Value of Linear Regression Model")
    parser.add_argument("--output_excel", type=str, default="chemical_shifts.xlsx", help="Output Excel file (default: chemical_shifts.xlsx).")
    args = parser.parse_args()
    
    updated_shifts = update_chemical_shifts_from_excel(args.output_excel,args.slope,args.intercept)
    print(updated_shifts)
    
    # Load the existing data from the input Excel file
    df = pd.read_excel(args.output_excel)
    
    # Convert the updated shifts dictionary to a DataFrame
    updated_shifts_df = pd.DataFrame(list(updated_shifts.items()), columns=['Atoms', 'Shift Values'])
    
    
    # Merge the existing DataFrame with the updated shifts DataFrame
    merged_df = pd.merge(df, updated_shifts_df, on='Atoms', how='left')
    
    # Save the updated DataFrame back to the same Excel file
    merged_df.to_excel(args.output_excel, index=False)
    print(f"Updated chemical shifts have been added to {args.output_excel}")
    
if __name__ == "__main__":
    main()


 
