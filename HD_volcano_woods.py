
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
from scipy.stats import t
import zipfile
import io





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
    t_value = t.ppf(1 - alpha / 2, df=(2 * n_replicates - 2))  # Two-tailed

    # Calculate the GST as the product of the t-value and SEM
    CI_result = t_value * SEM

    return round(CI_result, 2)  # Return the GST, rounded to 2 decimal places



def plot_sd_histogram_plotly(data, state1, state2):
    """
    Plot an interactive histogram of the Uptake SD values for the specified states using Plotly.

    Parameters:
        data (DataFrame): The dataset containing protein uptake data, including 'Uptake SD' and 'UptakeSD_State_2' columns.
        state1 (str): Name of the first state for comparison.
        state2 (str): Name of the second state for comparison.
    """
    # Extract Uptake SD for both states from merged_data
    state1_data = data['Uptake SD']
    state2_data = data['UptakeSD_State_2']

    # Create the histogram
    fig = go.Figure()

    # Add histogram for state1 data
    fig.add_trace(go.Histogram(
        x=state1_data,
        nbinsx=20,
        name=state1,
        marker=dict(color='skyblue', line=dict(color='black', width=1)),
        opacity=0.7
    ))

    # Add histogram for state2 data
    fig.add_trace(go.Histogram(
        x=state2_data,
        nbinsx=20,
        name=state2,
        marker=dict(color='lightgreen', line=dict(color='black', width=1)),
        opacity=0.7
    ))

    # Update layout for better appearance
    fig.update_layout(
        title="Distribution of Uptake SDs",
        xaxis_title="Uptake SD",
        yaxis_title="Frequency",
        template="simple_white",
        bargap=0.2,  # Space between bars
        margin=dict(l=40, r=40, t=40, b=40),
        xaxis=dict(
            showgrid=True, 
            gridcolor='lightgrey',
            showline=True,  # Show x-axis line
            linecolor='black',  # Set color for the x-axis line
            linewidth=1,  # Set thickness for the x-axis line
            mirror=True,  # Show axis line on the other side
        ),
        yaxis=dict(
            showgrid=True, 
            gridcolor='lightgrey',
            showline=True,  # Show y-axis line
            linecolor='black',  # Set color for the y-axis line
            linewidth=1,  # Set thickness for the y-axis line
            mirror=True,  # Show axis line on the other side
            title_standoff=15  # Space between title and axis line
        ),
        barmode='overlay',  # Overlay the histograms for comparison
        
        # Set the background color to white
        plot_bgcolor='white',  # Background color of the plot area
        paper_bgcolor='white', # Background color of the entire figure (including outside the plot area)
        
        width=700,
        height=500,
    )
 
    
    # Show the interactive plot in the Streamlit app
    st.plotly_chart(fig)

    # Prepare the HTML content for download
    html_data = fig.to_html(full_html=False)
    st.download_button(
        label="Download Histogram Plot (HTML)",
        data=html_data,
        file_name="histogram_plot.html",
        mime="text/html"
    )


def prepare_merged_data(data, state1, state2, alpha, gst):
    """
    Prepare and merge data from two specified states, calculate sum uptake and SD per peptide,
    and assess the significance of differences in uptake between the two states. 
    Add summary rows for summed values for each peptide.

    Parameters:
        data (DataFrame): Dataset containing protein uptake data, with columns 'State', 'Uptake', and 'Uptake SD'.
        state1 (str): Name of the first state for comparison.
        state2 (str): Name of the second state for comparison.
        alpha (float): Significance level for p-value threshold (e.g., 0.05).
        gst (float): Global significance threshold (GST) for assessing the magnitude of differences.

    Returns:
        DataFrame: Merged dataset with additional columns for difference, p-value, and significance.
    """

    # Calculate sum of uptake and uptake SD per peptide
    sum_stats = data.groupby(['State', 'Sequence']).agg(
        Sum_Uptake=('Uptake', 'sum'),
        Sum_Uptake_SD=('Uptake SD', 'sum')
    ).reset_index()

    # Filter data for each state
    state1_data = data[data['State'] == state1].copy()
    state2_data = data[data['State'] == state2].copy()

    # Initialize final datasets with original rows
    final_state1 = []
    final_state2 = []

    # Add sum rows for each peptide in State 1
    for peptide, group in state1_data.groupby('Sequence', sort=False):
        final_state1.append(group)
        sum_values = sum_stats[(sum_stats['State'] == state1) & (sum_stats['Sequence'] == peptide)]
        sum_row = pd.DataFrame([{
            'Sequence': peptide,
            'Start': group['Start'].iloc[0],
            'End': group['End'].iloc[0],
            'State': state1,
            'Exposure': '100000000',
            'Uptake': sum_values['Sum_Uptake'].values[0],
            'Uptake SD': sum_values['Sum_Uptake_SD'].values[0],
            'Start-End': group['Start-End'].iloc[0]
        }])
        final_state1.append(sum_row)

    # Add sum rows for each peptide in State 2
    for peptide, group in state2_data.groupby('Sequence', sort=False):
        final_state2.append(group)
        sum_values = sum_stats[(sum_stats['State'] == state2) & (sum_stats['Sequence'] == peptide)]
        sum_row = pd.DataFrame([{
            'Sequence': peptide,
            'Start': group['Start'].iloc[0],
            'End': group['End'].iloc[0],
            'State': state2,
            'Exposure': '100000000',
            'Uptake': sum_values['Sum_Uptake'].values[0],
            'Uptake SD': sum_values['Sum_Uptake_SD'].values[0],
            'Start-End': group['Start-End'].iloc[0]
        }])
        final_state2.append(sum_row)

    # Combine rows back into full DataFrame for each state
    state1_data = pd.concat(final_state1, ignore_index=True)
    state2_data = pd.concat(final_state2, ignore_index=True)

    # Rename summed columns for clarity
    state2_data = state2_data.rename(columns={'Uptake': 'Uptake_State_2', 'Uptake SD': 'UptakeSD_State_2'})

    # Merge datasets based on shared columns
    merged_data = pd.merge(
        state1_data, state2_data, 
        on=['Sequence', 'Start', 'End', 'Exposure'], 
        how='inner'
    )

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
    valid_p_values = merged_data['p_value'][np.isfinite(merged_data['p_value'])]

    merged_data['-np.log10(valid_p_values)'] = -np.log10(valid_p_values)

    # Determine significance based on alpha level and GST threshold
    merged_data['Significance'] = (merged_data['p_value'] < alpha) & \
                                  ((merged_data['Difference'] > gst) | (merged_data['Difference'] < -gst))

    # Save merged data to session state
    st.session_state['merged_data'] = merged_data

    # Add a button for downloading merged data as CSV
    if 'merged_data' in st.session_state:
        merged_data = st.session_state['merged_data']
        
        # Convert merged data to CSV format
        csv_data = merged_data.to_csv(index=False)
        st.download_button(
            label="Download Merged Data (CSV)",
            data=csv_data,
            file_name="merged_data.csv",
            mime="text/csv"
        )

    return merged_data




def create_volcano_plot(merged_data, gst, alpha):
    """
    Create and display a volcano plot to show differences in protein uptake and significance levels.
    
    Parameters:
        merged_data (DataFrame): Data containing 'Difference' and 'p_value' columns for plotting.
        gst (float): Global significance threshold for assessing the magnitude of differences.
        alpha (float): Significance level for determining p-value cutoff (e.g., 0.05).
        

    Returns:
        None: Displays the plot in Streamlit and saves it as an HTML file.
    """
    
    # Clip p-values to avoid extreme values and ensure realistic y-axis scaling
    merged_data['p_value'] = merged_data['p_value'].clip(lower=1e-100)
    
    # Debugging: Check p-value distribution
    print("P-value distribution:")
    print(merged_data['p_value'].describe())
    
    # Calculate the maximum y-axis value based on adjusted p-values
    valid_p_values = merged_data['p_value'][np.isfinite(merged_data['p_value'])]
    max_y_value = max(-np.log10(valid_p_values))

    # Calculate x-axis range based on the absolute differences
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
            hex_color = '#E5E5E5'  # Gray for non-significant points
        marker_colors.append(hex_color)

    # Scatter plot with markers for each point
    volcano_fig.add_trace(go.Scatter(
        x=merged_data['Difference'],
        y=-np.log10(merged_data['p_value']),
        mode='markers',
        marker=dict(color=marker_colors, size=10, line=dict(color='black', width=0.5)),
        opacity=0.9,
        text=merged_data['Start-End_y'],
    ))

    # Add significance and GST threshold lines
    volcano_fig.add_shape(type='line', x0=-max_x_value - buffer2, x1=max_x_value + buffer2,
                          y0=-np.log10(alpha), y1=-np.log10(alpha),
                          line=dict(color='black', dash='dash', width=0.5))
    volcano_fig.add_shape(type='line', x0=gst, x1=gst, y0=0, y1=max_y_value + buffer1,
                          line=dict(color='black', dash='dash', width=0.5))
    volcano_fig.add_shape(type='line', x0=-gst, x1=-gst, y0=0, y1=max_y_value + buffer1,
                          line=dict(color='black', dash='dash', width=0.5))

    # Update plot layout for a clean appearance
    volcano_fig.update_layout(
        title="Volcano Plot",
        xaxis_title="ΔD (Da)",
        yaxis_title="-log10(p-value)",
        margin=dict(l=40, r=40, t=40, b=40),
        plot_bgcolor='white',
        paper_bgcolor='white',
        width=700,
        height=500,
    )

    # Add shaded rectangles for significant regions
    volcano_fig.add_shape(type='rect',
                          x0=-max_x_value - buffer2, x1=-gst,
                          y0=-np.log10(alpha), y1=max_y_value + buffer1,
                          fillcolor='rgba(173, 216, 230, 0.05)', line=dict(width=0))
    volcano_fig.add_shape(type='rect',
                          x0=gst, x1=max_x_value + buffer2,
                          y0=-np.log10(alpha), y1=max_y_value + buffer1,
                          fillcolor='rgba(255, 182, 193, 0.05)', line=dict(width=0))

    # Ensure symmetrical x-axis and adjust axis appearance
    volcano_fig.update_xaxes(range=[-max_x_value - buffer2, max_x_value + buffer2], showline=True,
                             linecolor='black', linewidth=1, mirror=True)
    volcano_fig.update_yaxes(range=[-buffer1, max_y_value + buffer1], showline=True,
                             linecolor='black', linewidth=1, mirror=True)

    
    # Display the plot in Streamlit
    st.plotly_chart(volcano_fig)
# Convert the figure to HTML
    html_data = volcano_fig.to_html(full_html=False)

    # Add a download button for the volcano plot (HTML)
    st.download_button(
        label="Download Volcano Plot (HTML)",
        data=html_data,
        file_name="volcano_plot.html",
        mime="text/html")


def create_woods_plots(merged_data, gst, show_significant_only=False):
    """
    Create and display Woods plots for each time point, showing protein uptake differences.
    
    Parameters:
        merged_data (DataFrame): Data containing 'Difference' and 'Significance' columns for plotting.
        gst (float): Global significance threshold for evaluating the magnitude of differences.
        show_significant_only (bool): If True, only significant differences will be shown.
    
    Returns:
        fig_paths (list): List of file paths for the Woods plots.
    """
  
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

        # Filter out non-significant rows
        significant_group = group[group['Significance']]
    
        # Skip if no significant data
        if significant_group.empty:
            continue

        # Initialize figure for each time point
        woods_fig = go.Figure()

        # Reshape data for the heatmap
        z_data = group['Difference'].values.reshape(1, -1)
        max_end = group['End'].max()  # Maximum endpoint for x-axis range

        for _, row in group.iterrows():
            # If 'Show only significant differences' is checked, skip non-significant rows
            if show_significant_only and not row['Significance']:
                continue
        
            # Determine color based on significance and difference value
            color = custom_cmap(norm(row['Difference'])) if row['Significance'] else (0.9, 0.9, 0.9)
            hex_color = '#%02x%02x%02x' % (int(color[0] * 255), int(color[1] * 255), int(color[2] * 255))
            segment_start = row['Start']
            segment_end = row['End']
            
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
               hovertemplate="<b>Segment:</b> %{customdata[0]}-%{customdata[1]}<br>" +
                             "<b>Difference:</b> %{y}<br><extra></extra>"
            ))

        # Add heatmap for color scaling
        woods_fig.add_trace(go.Heatmap(
            z=z_data,
            x=group['Start-End_y'],
            colorscale=color_scale,
            showscale=True,
            colorbar=dict(title="ΔD (Da)", len=0.5, thickness=15),
            zmin=-max_diff,
            zmax=max_diff
        ))

        # Configure y-axis and add threshold lines
        woods_fig.update_yaxes(range=[-max_diff-0.05, max_diff+0.05], showline=True, linecolor='black', linewidth=2, mirror=True)
        woods_fig.add_hline(y=-gst, line=dict(color="rgba(169, 169, 169, 0.2)", width=1))
        woods_fig.add_hline(y=gst, line=dict(color="rgba(169, 169, 169, 0.2)", width=1))
        woods_fig.add_hrect(y0=-gst, y1=gst, line_width=0, fillcolor="lightgrey", opacity=0.1)

        # Final layout adjustments
        woods_fig.update_xaxes(range=[0, max_end], showline=True, linecolor='black', linewidth=2, mirror=True)
        woods_fig.update_layout(
            title=f"Woods Plot for Time Point: {time_point}",
            xaxis_title="Protein Sequence",
            yaxis_title="ΔD (Da)",
            plot_bgcolor='white',
            paper_bgcolor='white',
            margin=dict(l=40, r=40, t=40, b=40),
            showlegend=False,
            width=1500,
            height=500,
        )
        # Save the plot to a buffer in memory
        buffer = io.StringIO()
        woods_fig.write_html(buffer)
        buffer.seek(0)
        fig_paths.append((time_point, buffer))

    # Create a zip file containing all the HTML plots
    zip_buffer = io.BytesIO()
    with zipfile.ZipFile(zip_buffer, 'w', zipfile.ZIP_DEFLATED) as zip_file:
        for time_point, fig_buffer in fig_paths:
            zip_file.writestr(f"woods_plot_{time_point}.html", fig_buffer.getvalue())
    
    zip_buffer.seek(0)

    # Add a download button for the zip file
    st.download_button(
        label="Download All Woods Plots (ZIP)",
        data=zip_buffer,
        file_name="woods_plots.zip",
        mime="application/zip"
    )

    return fig_paths

def display_woods_plot(plot_tuple):
    """Display a Woods plot from a buffer."""
    try:
        # Extract the buffer from the tuple
        _, plot_buffer = plot_tuple
        
        # Get the HTML content from the buffer
        html_content = plot_buffer.getvalue()
        
        # Display the HTML content
        st.components.v1.html(html_content, height=600, scrolling=True)
    except Exception as e:
        st.error(f"Error displaying plot: {e}")




# Set the title of the app
st.title("Peptide Data Analysis")

# Custom CSS for styling
st.markdown("""
    <style>
        /* Smaller font size for the entire app */
        body {
            font-size: 12px;
        }
        
        /* Smaller font size for the title */
        h1, h2, h3, h4, h5, h6 {
            font-size: 20px;
        }
        
        /* Smaller font size for subheaders */
        .streamlit-expanderHeader {
            font-size: 14px;
        }
        
        /* Adjusting button size */
        .css-1emrehy.edgvbvh3 {
            font-size: 12px;  /* Decrease font size */
            padding: 5px 10px;  /* Reduce button padding */
        }
        
        /* Styling all columns to have white background and white outline */
        div[data-testid="stHorizontalBlock"] > div {
            background-color: white;  /* Set all columns to white */
            border: 1px solid #cccccc;  /* Change the border color to white */
            padding: 10px;  /* Add some padding inside the columns */
        }
    </style>
""", unsafe_allow_html=True)

# Create two columns with adjusted widths (left for inputs, right for plots)
left_column, right_column = st.columns([0.3, 1.7])

with left_column:
    st.markdown('<div class="stVerticalBlock">', unsafe_allow_html=True)
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
        st.session_state['alpha'] = alpha  # Save alpha to session state
        
        # State selection for comparison
        unique_states = data['State'].unique()
        state1 = st.selectbox("Select State 1:", unique_states)
        state2 = st.selectbox("Select State 2:", unique_states)
        
        # Add a checkbox for showing only significant differences
        show_significant_only = st.checkbox("Show only significant differences", value=False)

        # Input for specific peptides to skip (semicolon-separated Start-End pairs)
        peptides_to_skip = st.text_input("Enter peptide Start-End pairs to skip (e.g., 9-38;111-114):", "")
        if peptides_to_skip:
            # Parse the input into a list of tuples (Start, End), split by semicolon
            skip_peptides = [
                tuple(map(int, pair.split('-')))
                for pair in peptides_to_skip.split(';') if '-' in pair
            ]
            st.session_state['skip_peptides'] = skip_peptides  # Save to session state
        
        
        
        if st.button("Calculate"):
            # Filtering peptides based on exact Start-End pairs
            if peptides_to_skip:
                for start, end in st.session_state['skip_peptides']:
                    
                    
                    # Exclude data where both 'Start' and 'End' exactly match any specified pair
                    data = data[~((data['Start'] == start) & (data['End'] == end))]
        
            # Save the filtered data to session state
            st.session_state['filtered_data'] = data
           

            
            
            # Calculate GST and prepare merged data
            gst = calculate_gst(data, state1, state2, alpha)
            st.session_state['gst'] = gst
            st.success(f"Calculated GST for {state1} vs {state2}: {gst}")
        
            merged_data = prepare_merged_data(data, state1, state2, alpha, gst)
            st.session_state['merged_data'] = merged_data
            st.session_state['alpha'] = alpha  # Store alpha level


            


            # Generate Woods plots and store file paths in session state

            woods_fig_paths = create_woods_plots(merged_data, gst, show_significant_only=show_significant_only)
            st.session_state['woods_fig_paths'] = woods_fig_paths
    
    st.markdown('</div>', unsafe_allow_html=True)


# Top row: Volcano Plot and Uptake SD Histogram
volcano_col, histogram_col = st.columns([1, 1])

with volcano_col:
    st.subheader("Volcano Plot")
    if 'merged_data' in st.session_state and 'gst' in st.session_state and 'alpha' in st.session_state:
        create_volcano_plot(
            merged_data=st.session_state['merged_data'], 
            gst=st.session_state['gst'], 
            alpha=st.session_state['alpha'] 
            )
        
    else:
        st.info("Upload data and calculate GST to view the Volcano Plot.")

with histogram_col:
    st.subheader("Uptake SD Histogram")
    if 'merged_data' in st.session_state:
        plot_sd_histogram_plotly(
            data=st.session_state['merged_data'], 
            state1=state1, 
            state2=state2
        )
    else:
        st.info("Upload data and calculate GST to view the Histogram.")

# Second row: Woods Plot Viewer in the right column
with right_column:
    st.subheader("Woods Plot Viewer")
    st.markdown('<div class="stVerticalBlock">', unsafe_allow_html=True)
    if 'woods_fig_paths' in st.session_state:
        woods_fig_paths = st.session_state['woods_fig_paths']
        time_points = range(len(woods_fig_paths))

        selected_time = st.slider(
            "Select Time Point for Woods Plot",
            min_value=min(time_points),
            max_value=max(time_points),
            step=1,
            key="woods_plot_slider"
        )

        display_woods_plot(woods_fig_paths[selected_time])
    else:
        st.info("Upload data and generate Woods plots to view them.")
    st.markdown('</div>', unsafe_allow_html=True)
