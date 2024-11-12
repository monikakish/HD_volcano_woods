
import streamlit as st
st.set_page_config(layout="wide")  # Set the page to wide layout

import pandas as pd
import numpy as np
import plotly.graph_objs as go
import os
from scipy.stats import ttest_ind_from_stats
from matplotlib.colors import Normalize
import streamlit as st
from matplotlib.colors import LinearSegmentedColormap



def calculate_gst(data, state1, state2, alpha, n_replicates=3):
    """
    Calculate the GST (global significance threshold) value for differences in protein uptake 
    between two states, using the specified alpha level for significance.

    Parameters:
        data (DataFrame): The dataset containing protein uptake data, including 'State' and 'Uptake SD'.
        state1 (str): Name of the first state for comparison.
        state2 (str): Name of the second state for comparison.
        alpha (float): Significance level (e.g., 0.1, 0.05, 0.01) for calculating the t-value.
        n_replicates (int): Number of replicates per state (default is 3).

    Returns:
        float: The GST value rounded to 2 decimal places.
    """

    # Filter data for the selected states
    state1_data = data[data['State'] == state1]
    state2_data = data[data['State'] == state2]

    # Calculate the pooled standard deviation (SD)
    sd_pool_numerator = np.sum((data['Uptake SD'] ** 2) * (n_replicates - 1))  # Sum of squared SDs across replicates
    sd_pool_denominator = ((len(state1_data) - 1) + (len(state2_data) - 1)) * 2  # Degrees of freedom
    sd_pool = np.sqrt(sd_pool_numerator / sd_pool_denominator)  # Pooled SD

    # Calculate the standard error of the mean (SEM)
    SEM = np.sqrt((sd_pool ** 2 / n_replicates) * 2)

    # Lookup table for t-values corresponding to different alpha levels
    t_value = {0.1: 1.3, 0.05: 2.1, 0.01: 3.7}.get(alpha, None)  # Lookup based on alpha level

    # Calculate the GST as the product of the t-value and SEM
    CI_result = t_value * SEM

    return round(CI_result, 2)  # Return the GST, rounded to 2 decimal places


def prepare_merged_data(data, state1, state2, alpha, gst):
    """
    Prepare and merge data from two specified states, then calculate the significance 
    of differences in protein uptake between these states, along with p-values.

    Parameters:
        data (DataFrame): Dataset containing protein uptake data, with columns 'State', 'Uptake', and 'Uptake SD'.
        state1 (str): Name of the first state for comparison.
        state2 (str): Name of the second state for comparison.
        alpha (float): Significance level for p-value threshold (e.g., 0.05).
        gst (float): Global significance threshold (GST) for assessing the magnitude of differences.

    Returns:
        DataFrame: Merged dataset with additional columns for difference, p-value, and significance.
    """

    # Filter data for the selected states and rename columns for clarity
    state1_data = data[data['State'] == state1].copy()
    state2_data = data[data['State'] == state2].copy()
    state2_data = state2_data.rename(columns={'Uptake': 'Uptake_State_2', 'Uptake SD': 'UptakeSD_State_2'})

    # Merge datasets based on shared columns
    merged_data = pd.merge(state1_data, state2_data, on=['Sequence', 'Start', 'End', 'Exposure'], how='inner')

    # Calculate the difference in uptake between the two states
    merged_data['Difference'] = merged_data['Uptake'] - merged_data['Uptake_State_2']

    # Calculate p-values using t-test for each row and store them in a list
    p_values = []
    for _, row in merged_data.iterrows():
        res = ttest_ind_from_stats(
            mean1=row['Uptake'], std1=row['Uptake SD'], nobs1=3,
            mean2=row['Uptake_State_2'], std2=row['UptakeSD_State_2'], nobs2=3
        )
        p_values.append(res.pvalue)
    
    # Add p-values to the merged data
    merged_data['p_value'] = p_values

    # Determine significance based on alpha level and GST threshold
    merged_data['Significance'] = (merged_data['p_value'] < alpha) & \
                                  ((merged_data['Difference'] > gst) | (merged_data['Difference'] < -gst))

    return merged_data


def create_volcano_plot(merged_data, gst, alpha, output_dir):
    """
    Create and display a volcano plot to show differences in protein uptake and significance levels.
    
    Parameters:
        merged_data (DataFrame): Data containing 'Difference' and 'p_value' columns for plotting.
        gst (float): Global significance threshold for assessing the magnitude of differences.
        alpha (float): Significance level for determining p-value cutoff (e.g., 0.05).
        output_dir (str): Directory to save the volcano plot as an HTML file.

    Returns:
        None: Displays the plot in Streamlit and saves it as an HTML file.
    """
    
    # Determine maximum values for the plot range, with buffer for better visibility
    valid_p_values = merged_data['p_value'][np.isfinite(merged_data['p_value'])]
    max_y_value = max(-np.log10(valid_p_values))
    max_x_value = max(abs(merged_data['Difference'].min()), abs(merged_data['Difference'].max()))
    buffer1 = max_y_value / 10
    buffer2 = max_x_value / 10

    # Normalize color mapping from -max_x_value to max_x_value
    norm = Normalize(vmin=-max_x_value, vmax=max_x_value)
    volcano_fig = go.Figure()

    # Define a blue-to-white-to-red colormap
    custom_cmap = LinearSegmentedColormap.from_list("custom_cmap", ["blue", "white", "red"])

    # Generate consistent marker colors based on significance and difference values
    marker_colors = []
    for _, row in merged_data.iterrows():
        if row['Significance']:
            color = custom_cmap(norm(row['Difference']))
            hex_color = '#%02x%02x%02x' % (int(color[0] * 255), int(color[1] * 255), int(color[2] * 255))
        else:
            hex_color = '#E5E5E5'  # Gray color for non-significant points
        marker_colors.append(hex_color)

    # Scatter plot with markers for each point in merged_data
    volcano_fig.add_trace(go.Scatter(
        x=merged_data['Difference'],
        y=-np.log10(merged_data['p_value']),
        mode='markers',
        marker=dict(color=marker_colors, size=10, line=dict(color='black', width=0.5)),
        opacity=0.9,
        text=merged_data['Start-End_y'],
    ))

    # Add threshold and cutoff lines for significance and GST thresholds
    volcano_fig.add_shape(type='line', x0=-max_x_value - buffer2, x1=max_x_value + buffer2,
                          y0=-np.log10(alpha), y1=-np.log10(alpha),
                          line=dict(color='black', dash='dash', width=0.5))
    volcano_fig.add_shape(type='line', x0=gst, x1=gst, y0=0, y1=max_y_value + buffer1,
                          line=dict(color='black', dash='dash', width=0.5))
    volcano_fig.add_shape(type='line', x0=-gst, x1=-gst, y0=0, y1=max_y_value + buffer1,
                          line=dict(color='black', dash='dash', width=0.5))

    # Update layout properties for plot appearance
    volcano_fig.update_layout(
        title="Volcano Plot",
        xaxis_title="ΔD (Da)",
        yaxis_title="-log10(p-value)",
        margin=dict(l=40, r=40, t=40, b=40),
        plot_bgcolor='white',
        paper_bgcolor='white',
        width=600,
        height=400,
    )

    # Add shaded rectangles to highlight significance regions
    volcano_fig.add_shape(type='rect',
                          x0=-max_x_value - buffer2, x1=-gst,
                          y0=-np.log10(alpha), y1=max_y_value + buffer1,
                          fillcolor='rgba(173, 216, 230, 0.05)', line=dict(width=0))
    volcano_fig.add_shape(type='rect',
                          x0=gst, x1=max_x_value + buffer2,
                          y0=-np.log10(alpha), y1=max_y_value + buffer1,
                          fillcolor='rgba(255, 182, 193, 0.05)', line=dict(width=0))

    # Ensure symmetrical x-axis and adjust axis appearance
    volcano_fig.update_xaxes(range=[-max_x_value-buffer2, max_x_value+buffer2], showline=True, linecolor='black', linewidth=1, mirror=True)
    volcano_fig.update_yaxes(range=[-buffer1, max_y_value+buffer1], showline=True, linecolor='black', linewidth=1, mirror=True)

    # Define color scale for the heatmap, centered at 0
    color_scale = [[0, "blue"], [0.5, 'white'], [1, "red"]]
    z_data = np.array([merged_data['Difference'].values, merged_data['Difference'].values])

    # Add heatmap trace for color legend and visual effect
    volcano_fig.add_trace(go.Heatmap(
        z=z_data,
        x=merged_data['Start-End_y'],
        colorscale=color_scale,
        showscale=True,
        zmin=-max_x_value,
        zmax=max_x_value,
        colorbar=dict(title="ΔD (Da)", titleside='right', len=0.5, thickness=15),
        visible=True
    ))

    # Display the plot in Streamlit
    st.plotly_chart(volcano_fig)

    # Save the plot as an HTML file in the specified output directory
    fig_path = os.path.join(output_dir, "volcano_plot.html")
    volcano_fig.write_html(fig_path)
    

def create_woods_plots(merged_data, gst, output_dir="woods_plots", show_significant_only=False):
    """
    Create and save Woods plots for each time point, showing protein uptake differences.
    
    Parameters:
        merged_data (DataFrame): Data containing 'Difference' and 'Significance' columns for plotting.
        gst (float): Global significance threshold for evaluating the magnitude of differences.
        output_dir (str): Directory to save the Woods plot HTML files.
    
    Returns:
        fig_paths (list): List of file paths to saved Woods plots.
    """
    
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Set normalization based on the range of differences
    max_diff = max(abs(merged_data['Difference'].min()), abs(merged_data['Difference'].max()))
    norm = Normalize(vmin=-max_diff, vmax=max_diff)

    # Define a blue-to-white-to-red colormap for consistent coloring
    custom_cmap = LinearSegmentedColormap.from_list("custom_cmap", ["blue", "white", "red"])
    color_scale = [[0, "blue"], [0.5, "white"], [1, "red"]]

    fig_paths = []
    for time_point, group in merged_data.groupby('Exposure'):
        # Skip plots for time point 0 if present
        if time_point == 0:
            continue

        # Initialize figure for each time point
        woods_fig = go.Figure()

        # Reshape data for the heatmap
        z_data = group['Difference'].values.reshape(1, -1)

        # Add bars for each protein segment’s uptake difference
        for _, row in group.iterrows():
            # If 'Show only significant differences' is checked, skip non-significant rows
            if show_significant_only and not row['Significance']:
                continue
        
            # Determine color based on significance and difference value
            color = custom_cmap(norm(row['Difference'])) if row['Significance'] else (0.9, 0.9, 0.9)
            hex_color = '#%02x%02x%02x' % (int(color[0] * 255), int(color[1] * 255), int(color[2] * 255))
            # Set segment details for the hover template
            segment_start = row['Start']
            segment_end = row['End']
            segment_length = segment_end - segment_start
            segment_name = row['Start-End_y']
            
            
            
            woods_fig.add_trace(go.Bar(
                x=[row['End'] - row['Start']],  # Protein segment length
                y=[row['Difference']],          # Difference in uptake
                orientation='h',
                marker=dict(color=hex_color, line=dict(color='black', width=0.5)),
                opacity=0.7,
                base=row['Start'],              # Start point of the protein segment
                width=0.03,                     # Bar width for visual clarity
                name=row['Start-End_y'],         # Label each bar with protein segment information
               customdata=[[segment_start, segment_end]],  # Pass Start and End as custom data
               hovertemplate="<b>Segment:</b> %{customdata[0]}-%{customdata[1]}<br>" +  # Use customdata for Start and End
                      "<b>Difference:</b> %{y}<br><extra></extra>"  # Display the Difference in uptake
                          ))
        # Add heatmap for color scaling
        woods_fig.add_trace(go.Heatmap(
            z=z_data,
            x=group['Start-End_y'],
            colorscale=color_scale,
            showscale=True,
            colorbar=dict(title="ΔD (Da)", titleside='right', len=0.5, thickness=15),
            zmin=-max_diff,
            zmax=max_diff
        ))

        # Define y-axis range with a buffer for clarity in visualization
        buffer = 0.05
        absolute_max = group['Difference'].abs().max()

        # Configure y-axis with symmetrical range and buffer
        woods_fig.update_yaxes(range=[-max_diff-buffer, max_diff+buffer], showline=True, linecolor='black', linewidth=2, mirror=True, title_standoff=15)

        # Add threshold lines for gst regions to highlight significance boundaries
        woods_fig.add_hline(y=-gst, line=dict(color="rgba(169, 169, 169, 0.2)", width=1))
        woods_fig.add_hline(y=gst, line=dict(color="rgba(169, 169, 169, 0.2)", width=1))
        woods_fig.add_hrect(y0=-gst, y1=gst, line_width=0, fillcolor="lightgrey", opacity=0.1)

        # Final layout adjustments and axis configuration
        woods_fig.update_xaxes(showline=True, linecolor='black', linewidth=2, mirror=True, zeroline=False, zerolinecolor='black', zerolinewidth=2)
        woods_fig.update_layout(
            title=f"Woods Plot for Time Point: {time_point}",
            xaxis_title="Protein Sequence",
            yaxis_title="ΔD (Da)",
            plot_bgcolor='white',
            paper_bgcolor='white',
            margin=dict(l=40, r=40, t=40, b=40),
            showlegend=False,
            width=1000,  # Set plot width
            height=400,
            # Set plot height
        )

        # Save the figure and record file path
        fig_path = os.path.join(output_dir, f"woods_plot_{time_point}.html")
        woods_fig.write_html(fig_path)
        fig_paths.append(fig_path)

    return fig_paths  # Return the list of saved file paths


def display_woods_plot(plot_path):
    """Display a Woods plot HTML file."""
    try:
        with open(plot_path, 'r', encoding='utf-8') as f:
            html_content = f.read()
        st.components.v1.html(html_content, height=600, scrolling=True)
    except FileNotFoundError:
        st.error(f"File {plot_path} not found.")


def main():
    """
    Streamlit app for peptide data analysis. Allows users to upload data, select parameters,
    generate volcano and Woods plots, and view results interactively.
    """
    
    # Set the title of the app
    st.title("Peptide Data Analysis")

    # Apply custom CSS to adjust font size for input fields in the left column
    st.markdown(
        """
        <style>
        .css-1r7m3o2 { font-size: 10px; }
        .css-1ht1p4w { font-size: 10px; }
        .stSelectbox label, .stNumberInput label, .stTextInput label {
            font-size: 12px !important;
        }
        </style>
        """, 
        unsafe_allow_html=True
    )

    # Create two columns with adjusted widths (left for inputs, right for plots)
    left_column, right_column = st.columns([0.3, 1.7])

    # Left column for input fields
    with left_column:
        st.subheader("Input Parameters")

        # File uploader for CSV input data
        input_file = st.file_uploader("Upload Data CSV", type=["csv"], key="input_file_uploader")
        if input_file is not None:
            # Load data and format 'Start-End' column
            data = pd.read_csv(input_file)
            data["Start-End"] = data["Start"].astype(str) + "-" + data["End"].astype(str)

            # Alpha level selection
            alpha = st.number_input(
                "Select Alpha Level (e.g., 0.1, 0.05, 0.01):",
                min_value=0.01, max_value=1.0, value=0.05, step=0.01
            )
            
            # State selection for comparison
            unique_states = data['State'].unique()
            state1 = st.selectbox("Select State 1:", unique_states)
            state2 = st.selectbox("Select State 2:", unique_states)
            
            
            # Add a checkbox for showing only significant differences
            show_significant_only = st.checkbox("Show only significant differences", value=False)
            
            
            # Output directory for saving Woods plots
            output_dir = st.text_input("Enter directory to save Woods plots (leave empty for current directory):", "")

            # Button to initiate calculations and plot generation
            if st.button("Calculate"):
                # Calculate GST for the selected states
                gst = calculate_gst(data, state1, state2, alpha)
                st.success(f"Calculated GST for {state1} vs {state2}: {gst}")

                # Prepare merged data for plotting
                merged_data = prepare_merged_data(data, state1, state2, alpha, gst)
                
                # Generate Woods plots and store file paths in session state
                woods_fig_paths = create_woods_plots(merged_data, gst, output_dir, show_significant_only)
                st.session_state['woods_fig_paths'] = woods_fig_paths

                # Generate and display volcano plot in the right column
                with right_column:
                    st.subheader("Volcano Plot")
                    create_volcano_plot(merged_data, gst, alpha, output_dir)

    # Right column for displaying volcano plot (top) and Woods plot viewer (bottom)
    with right_column:
        # Check if Woods plot paths are stored in session state
        if 'woods_fig_paths' in st.session_state:
            st.subheader("Woods Plot Viewer")

            woods_fig_paths = st.session_state['woods_fig_paths']
            time_points = range(len(woods_fig_paths))

            # Slider to select time point for Woods plot display
            selected_time = st.slider(
                "Select Time Point for Woods Plot",
                min_value=min(time_points),
                max_value=max(time_points),
                step=1
            )
            
            # Display the selected Woods plot
            display_woods_plot(woods_fig_paths[selected_time])

# Run the app
if __name__ == "__main__":
    main()



