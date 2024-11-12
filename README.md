# HD_volcano_woods
Volcano and woods plots for HDX data
Peptide Data Analysis Streamlit App
This Streamlit app allows users to analyze peptide uptake data, generate volcano plots and Woods plots, and visualize differences in protein uptake between two states. It supports interactive data visualization, enabling easy exploration of statistical significance and uptake differences.

Features
GST Calculation: Computes the Global Significance Threshold (GST) to assess differences in peptide uptake between two states.
Volcano Plot: Displays differences in uptake versus p-value, highlighting statistically significant differences.
Woods Plots: Generates heatmap-based Woods plots for each exposure time, with the option to filter by significant differences.
Customizable Parameters: Users can select the significance threshold (alpha), states for comparison, and decide whether to display only significant differences.

Usage
Upload Data: Upload a CSV file containing peptide uptake data with columns such as State, Sequence, Start, End, Exposure, Uptake, and Uptake SD.
Set Parameters: Choose the significance level (alpha), and select two states for comparison.
Generate Plots: Click the “Calculate” button to compute the GST, generate a volcano plot, and Woods plots for each exposure time.
View Results: Interactive volcano plots will be displayed on the right panel. Woods plots will be saved to the specified output directory.
Example
Upload a CSV file with peptide data (for example, peptide_data.csv).
Select State 1 and State 2 from the dropdowns.
Adjust the Alpha Level (e.g., 0.05) and check the option to Show Only Significant Differences if desired.
Click Calculate to view the results and download the Woods plots as HTML files.


Contributions
Contributions are welcome! If you find bugs or have ideas for new features, feel free to open an issue or submit a pull request.

License
