
import streamlit as st
st.set_page_config(layout="wide")  # Set the page to wide layout

import pandas as pd
import numpy as np
import plotly.graph_objs as go
import os
from scipy.stats import ttest_ind_from_stats
import matplotlib.cm as cm
from matplotlib.colors import Normalize
import os
import numpy as np
import plotly.graph_objects as go
from matplotlib import cm
from matplotlib.colors import Normalize
import streamlit as st
from matplotlib.colors import LinearSegmentedColormap




def calculate_gst(data, state1, state2, alpha, n_replicates=3):
    """Calculate GST value based on the selected states and alpha level."""
    state1_data = data[data['State'] == state1]
    state2_data = data[data['State'] == state2]

    sd_pool_numerator = np.sum((data['Uptake SD']**2) * (n_replicates - 1))
    sd_pool_denominator = ((len(state1_data) - 1) + (len(state2_data) - 1)) * 2
    sd_pool = np.sqrt(sd_pool_numerator / sd_pool_denominator)
    SEM = np.sqrt((sd_pool**2 / n_replicates) * 2)

    t_value = {0.1: 1.3, 0.05: 2.1, 0.01: 3.7}.get(alpha, None)
    CI_result = t_value * SEM
    return round(CI_result, 2)

def prepare_merged_data(data, state1, state2, alpha, gst):
    """Prepare and merge data from the two states, calculating significance and p-values."""
    state1_data = data[data['State'] == state1].copy()
    state2_data = data[data['State'] == state2].copy()
    state2_data = state2_data.rename(columns={'Uptake': 'Uptake_State_2', 'Uptake SD': 'UptakeSD_State_2'})
    merged_data = pd.merge(state1_data, state2_data, on=['Sequence', 'Start', 'End', 'Exposure'], how='inner')
    merged_data['Difference'] = merged_data['Uptake'] - merged_data['Uptake_State_2']

    p_values = []
    for _, row in merged_data.iterrows():
        res = ttest_ind_from_stats(row['Uptake'], row['Uptake SD'], 3, row['Uptake_State_2'], row['UptakeSD_State_2'], 3)
        p_values.append(res.pvalue)
    merged_data['p_value'] = p_values
    merged_data['Significance'] = (merged_data['p_value'] < alpha) & \
                                  ((merged_data['Difference'] > gst) | (merged_data['Difference'] < -gst))

    return merged_data



def create_volcano_plot(merged_data, gst, alpha, output_dir):
    import os
    import numpy as np
    import plotly.graph_objs as go
    from matplotlib.colors import Normalize, LinearSegmentedColormap

    valid_p_values = merged_data['p_value'][np.isfinite(merged_data['p_value'])]
    max_y_value = max(-np.log10(valid_p_values))
    max_x_value = max(abs(merged_data['Difference'].min()), abs(merged_data['Difference'].max()))
    buffer1 = max_y_value / 10
    buffer2 = max_x_value / 10

    # Set normalization from -max_x_value to max_x_value
    norm = Normalize(vmin=-max_x_value, vmax=max_x_value)
    volcano_fig = go.Figure()

    # Define a blue to white to red colormap
    custom_cmap = LinearSegmentedColormap.from_list("custom_cmap", ["blue", "white", "red"])

    # Create consistent marker colors based on the colormap
    marker_colors = []
    for _, row in merged_data.iterrows():
        color = custom_cmap(norm(row['Difference'])) if row['Significance'] else (0.9, 0.9, 0.9)
        hex_color = '#%02x%02x%02x' % (int(color[0] * 255), int(color[1] * 255), int(color[2] * 255))
        marker_colors.append(hex_color)

    # Scatter plot for markers
    volcano_fig.add_trace(go.Scatter(
        x=merged_data['Difference'],
        y=-np.log10(merged_data['p_value']),
        mode='markers',
        marker=dict(color=marker_colors, size=10, line=dict(color='black', width=0.5)),
        opacity=0.9,
        text=merged_data['Start-End_y'],
    ))

    # Add threshold and cutoff lines
    volcano_fig.add_shape(type='line', x0=-max_x_value - buffer2, x1=max_x_value + buffer2,
                          y0=-np.log10(alpha), y1=-np.log10(alpha),
                          line=dict(color='black', dash='dash',width=0.5))
    volcano_fig.add_shape(type='line', x0=gst, x1=gst,
                          y0=0, y1=max_y_value + buffer1,
                          line=dict(color='black', dash='dash',width=0.5))
    volcano_fig.add_shape(type='line', x0=-gst, x1=-gst,
                          y0=0, y1=max_y_value + buffer1,
                          line=dict(color='black', dash='dash',width=0.5))

    # Update layout
    volcano_fig.update_layout(
        title="Volcano Plot",
        xaxis_title="ΔD (Da)",
        yaxis_title="-log10(p-value)",
        margin=dict(l=40, r=40, t=40, b=40),
        plot_bgcolor='white',
        paper_bgcolor='white',
        width=600,  # Set the width of the plot (adjust as needed)
        height=400,  # Set the height of the plot (adjust as needed)
    )
    

    # Shaded rectangles between lines
    volcano_fig.add_shape(type='rect',
                          x0=-max_x_value - buffer2, x1=-gst,
                          y0=-np.log10(alpha), y1=max_y_value + buffer1,
                          fillcolor='rgba(173, 216, 230, 0.05)', line=dict(width=0))
    volcano_fig.add_shape(type='rect',
                          x0=gst, x1=max_x_value + buffer2,
                          y0=-np.log10(alpha), y1=max_y_value + buffer1,
                          fillcolor='rgba(255, 182, 193, 0.05)', line=dict(width=0))

    # Ensure symmetrical x-axis
    volcano_fig.update_xaxes(range=[-max_x_value-buffer2, max_x_value+buffer2], showline=True, linecolor='black', linewidth=1, mirror=True)
    
    
    
    volcano_fig.update_yaxes(range=[-buffer1, max_y_value+buffer1],showline=True, linecolor='black', linewidth=1, mirror=True)

    # Define color scale for the heatmap with center at 0
    color_scale = [[0, "blue"], [0.5, 'white'], [1, "red"]]
    z_data = np.array([merged_data['Difference'].values, merged_data['Difference'].values])

    # Create heatmap with the blue-to-red colormap
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

    # Save the figure
    fig_path = os.path.join(output_dir, "volcano_plot.html")
    volcano_fig.write_html(fig_path)
    print(f"Volcano plot saved to {fig_path}")




def create_woods_plots(merged_data, gst, output_dir="woods_plots"):
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Set normalization based on the range of differences
    max_diff = max(abs(merged_data['Difference'].min()), abs(merged_data['Difference'].max()))
    norm = Normalize(vmin=-max_diff, vmax=max_diff)

    # Define a blue to white to red colormap for consistent coloring
    custom_cmap = LinearSegmentedColormap.from_list("custom_cmap", ["blue", "white", "red"])
    color_scale = [[0, "blue"], [0.5, "white"], [1, "red"]]

    fig_paths = []
    for time_point, group in merged_data.groupby('Exposure'):
        # Skip plots for time point 0
        if time_point == 0:
            continue

        woods_fig = go.Figure()

        # Reshape data for the heatmap
        z_data = group['Difference'].values.reshape(1, -1)

        # Add bars for each protein segment’s uptake difference
        for _, row in group.iterrows():
            color = custom_cmap(norm(row['Difference'])) if row['Significance'] else (0.9, 0.9, 0.9)
            hex_color = '#%02x%02x%02x' % (int(color[0] * 255), int(color[1] * 255), int(color[2] * 255))

            woods_fig.add_trace(go.Bar(
                x=[row['End'] - row['Start']],
                y=[row['Difference']],
                orientation='h',
                marker=dict(color=hex_color, line=dict(color='black', width=0.5)),
                opacity=0.7,
                base=row['Start'],
                width=0.03,
                name=row['Start-End_y']
            ))

        # Add heatmap trace using the consistent color scale
        woods_fig.add_trace(go.Heatmap(
            z=z_data,
            x=group['Start-End_y'],
            colorscale=color_scale,
            showscale=True,
            colorbar=dict(title="ΔD (Da)", titleside='right', len=0.5, thickness=15),
            zmin=-max_diff,
            zmax=max_diff
        ))

        # Define y-axis range with a buffer for clarity
        buffer = 0.2
        y_min = group['Difference'].min() - buffer
        y_max = group['Difference'].max() + buffer
        absolute_max = group['Difference'].abs().max()

        woods_fig.update_yaxes(range=[-absolute_max, absolute_max], showline=True, linecolor='black', linewidth=2, mirror=True, title_standoff=15)
        # Add threshold lines for gst regions
        woods_fig.add_hline(y=-gst, line=dict(color="rgba(186, 85, 211, 0.2)", width=1))
        woods_fig.add_hline(y=gst, line=dict(color="rgba(186, 85, 211, 0.2)", width=1))
        woods_fig.add_hrect(y0=-gst, y1=gst, line_width=0, fillcolor="lightgrey", opacity=0.1)

        # Configure final layout and axis settings
        woods_fig.update_xaxes(showline=True, linecolor='black', linewidth=2, mirror=True, zeroline=False, zerolinecolor='black', zerolinewidth=2)
        woods_fig.update_layout(
            title=f"Woods Plot for Time Point: {time_point}",
            xaxis_title="Protein Sequence",
            yaxis_title="ΔD (Da)",
            plot_bgcolor='white',
            paper_bgcolor='white',
            margin=dict(l=40, r=40, t=40, b=40),
            showlegend=False,
            width=1000,  # Set the width of the plot (adjust as needed)
            height=400,  # Set the height of the plot (adjust as needed)
        )
        

        # Save the figure and store file path
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



# Define other functions (calculate_gst, prepare_merged_data, create_volcano_plot, create_woods_plots, display_woods_plot)

def main():
    # Set the title of the app
    st.title("Peptide Data Analysis")

    # Apply custom CSS to adjust the font size for input fields in the left column
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

    # Create two columns with adjusted widths (left is smaller for inputs)
    left_column, right_column = st.columns([0.3, 1.7])

    # Left column for input fields
    with left_column:
        st.subheader("Input Parameters")
        
        # Upload file
        input_file = st.file_uploader("Upload Data CSV", type=["csv"], key="input_file_uploader")
        if input_file is not None:
            data = pd.read_csv(input_file)
            data["Start-End"] = data["Start"].astype(str) + "-" + data["End"].astype(str)

            # Parameter inputs
            alpha = st.number_input("Select Alpha Level (e.g., 0.1, 0.05, 0.01):", min_value=0.01, max_value=1.0, value=0.05, step=0.01)
            unique_states = data['State'].unique()
            state1 = st.selectbox("Select State 1:", unique_states)
            state2 = st.selectbox("Select State 2:", unique_states)
            
            
            
            
            output_dir = st.text_input("Enter directory to save Woods plots (leave empty for current directory):", "")

            # Calculate GST and generate plots
            if st.button("Calculate"):
                gst = calculate_gst(data, state1, state2, alpha)
                st.success(f"Calculated GST for {state1} vs {state2}: {gst}")

                # Prepare merged data
                merged_data = prepare_merged_data(data, state1, state2, alpha, gst)
                
                # Generate Woods plots and store paths in session state
                woods_fig_paths = create_woods_plots(merged_data, gst, output_dir)
                st.session_state['woods_fig_paths'] = woods_fig_paths

                # Generate volcano plot and display in the right column
                with right_column:
                    st.subheader("Volcano Plot")
                    create_volcano_plot(merged_data, gst, alpha, output_dir)

    # Right column for displaying volcano plot (top) and Woods plot viewer (bottom)
    with right_column:
        if 'woods_fig_paths' in st.session_state:
            st.subheader("Woods Plot Viewer")

            woods_fig_paths = st.session_state['woods_fig_paths']
            time_points = range(len(woods_fig_paths))

            # Slider to select time point for the Woods plot
            selected_time = st.slider("Select Time Point for Woods Plot", min_value=min(time_points), max_value=max(time_points), step=1)
            
            # Display the selected Woods plot below the volcano plot
            display_woods_plot(woods_fig_paths[selected_time])

if __name__ == "__main__":
    main()








