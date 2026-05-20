import argparse
import pandas as pd
import matplotlib.pyplot as plt
import glob
import os
import sys

def main():
    parser = argparse.ArgumentParser(description="Plot CFD relaxation data from CSV files.")
    parser.add_argument("--case", nargs=2, action='append', metavar=('LABEL', 'PATH'),
                        help="Case label and directory/file path (can be repeated)")
    parser.add_argument("--field", type=str, required=True, 
                        help="Field name to search for in filenames (e.g., TTR, TVib)")
    parser.add_argument("--title", type=str, default="Temperature Relaxation", 
                        help="Title of the plot")
    parser.add_argument("--output", type=str, default="plot.png", 
                        help="Output image filename")
    parser.add_argument("--export-csv", action='store_true', 
                        help="Export combined data to 'combined_data.csv'")

    args = parser.parse_args()

    if not args.case:
        print("Error: No cases provided. Use --case 'Label' 'Path'")
        sys.exit(1)

    plt.style.use('seaborn-v0_8-paper') # Scientific style
    plt.figure(figsize=(8, 6))
    
    combined_dfs = []

    for label, path in args.case:
        # Search for the file if path is a directory
        if os.path.isdir(path):
            search_pattern = os.path.join(path, f"*{args.field}*.csv")
            files = glob.glob(search_pattern)
            if not files:
                print(f"Warning: No files matching *{args.field}*.csv found in {path}")
                continue
            file_path = files[0]
        else:
            file_path = path

        # Load data
        try:
            df = pd.read_csv(file_path)
            # Find columns (assuming 'Time' is in the first column)
            time_col = [c for c in df.columns if 'Time' in c][0]
            temp_col = [c for c in df.columns if 'Temp' in c][0]
            
            # Sort by time for plotting consistency
            df = df.sort_values(by=time_col)
            
            # Plot data (skipping t=0 for log scale if necessary)
            plot_df = df[df[time_col] > 0]
            plt.plot(plot_df[time_col], plot_df[temp_col], label=label, linewidth=1.5)
            
            if args.export_csv:
                df_to_add = df[[time_col, temp_col]].copy()
                df_to_add.columns = ['Time (s)', f'{label}_{args.field}']
                combined_dfs.append(df_to_add)
                
        except Exception as e:
            print(f"Error processing {file_path}: {e}")

    # Formatting
    plt.xscale('log')
    plt.xlabel(r'Time ($s$)', fontsize=12)
    plt.ylabel(r'Temperature ($K$)', fontsize=12)
    plt.title(args.title, fontsize=14)
    plt.grid(True, which="both", ls="-", alpha=0.5)
    plt.legend()
    plt.tight_layout()
    
    # Save plot
    plt.savefig(args.output, dpi=300)
    print(f"Plot saved to {args.output}")

    # Export CSV if requested
    if args.export_csv and combined_dfs:
        from functools import reduce
        final_df = reduce(lambda left, right: pd.merge(left, right, on='Time (s)', how='outer'), combined_dfs)
        final_df.to_csv("combined_data.csv", index=False)
        print("Combined data exported to combined_data.csv")

if __name__ == "__main__":
    main()
