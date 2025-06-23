import matplotlib.pyplot as plt
import numpy as np
import ipywidgets as widgets
from IPython.display import display, clear_output
import pandas as pd
import os
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

np.float_ = np.float64 # Solve the issue. Also possible to pip install "numpy<2"
from obspy.taup import velocity_model
from obspy.taup import TauPyModel



def getDfFromND(model):
    # Extract relevant data from model.layers
    depths = []
    vp_values = []
    vs_values = []
    qp_values = []
    qs_values = []
    density_values=[]
    for layer in model.layers:
        top_depth = layer[0]
        bot_depth = layer[1]
        depths.extend([top_depth, bot_depth])
        
        vp_values.extend([layer[2], layer[3]])
        vs_values.extend([layer[4], layer[5]])
        density_values.extend([layer[6], layer[7]])
        qp_values.extend([layer[8], layer[9]])
        qs_values.extend([layer[10], layer[11]])

    # Create DataFrame
    df = pd.DataFrame({
        "Depth": depths,
        "Vp": vp_values,
        "Vs": vs_values,
        'Density' : density_values,
        "Qkappa": qp_values,
        "Qmu": qs_values
    })
    
    return df

def interactive_point_selector(x, model_values, npts=1):
    """
    Function to interactively select points from a model plot, with discontinuity flag.
    Returns the selected points as a list of (index, value, discontinuity_flag).
    """
    selected_points = []
    max_points = 0
    click_counter = 0
    pending_selection = None

    # Widgets
    point_number_selector = widgets.IntText(
        value=npts, description='Number of points:', disabled=False
    )
    start_button = widgets.Button(description="Start Selection")
    output_area = widgets.Textarea(
        value='', placeholder='Selected points will appear here...',
        description='Selected:', layout=widgets.Layout(width='50%', height='120px')
    )
    confirm_button = widgets.Button(description="Confirm Selection", button_style='success')
    reject_button = widgets.Button(description="Reject Selection", button_style='danger')

    # 💡 Here's your missing discontinuity selector
    discontinuity_selector = widgets.ToggleButton(
        value=False,
        description='Discontinuity',
        button_style='warning'
    )

    selection_controls = widgets.HBox([confirm_button, reject_button, discontinuity_selector])

    final_output = widgets.Output()

    # Display widgets
    display(point_number_selector, start_button, output_area)
    display(selection_controls)
    display(final_output)

    # Create the plot
    fig, ax = plt.subplots(figsize=(4, 8))
    line, = ax.plot(x, model_values, '-o', picker=5)
    ax.set_title("Click on points to select them")
    selection_controls.layout.visibility = 'hidden'

    def start_selection(b):
        nonlocal max_points, click_counter, pending_selection
        selected_points.clear()
        max_points = point_number_selector.value
        click_counter = 0
        pending_selection = None
        final_output.clear_output()
        output_area.value = ''
        selection_controls.layout.visibility = 'hidden'
        ax.set_title(f"Click to select {max_points} points")
        fig.canvas.draw()

    start_button.on_click(start_selection)

    def on_click(event):
        nonlocal click_counter, pending_selection
        if event.inaxes != ax or click_counter >= max_points:
            return
        click_x, click_y = event.xdata, event.ydata
        distances = np.sqrt((x - click_x)**2 + (model_values - click_y)**2)
        ix = distances.argmin()
        pending_selection = (ix, model_values[ix])
        ax.plot(x[ix], model_values[ix], 'yo')
        fig.canvas.draw()
        selection_controls.layout.visibility = 'visible'

    def confirm_selection(b):
        nonlocal click_counter, pending_selection
        if pending_selection is None:
            return
        ix, val = pending_selection
        disc_flag = discontinuity_selector.value
        selected_points.append((ix, val, disc_flag))  # ✅ include the discontinuity flag
        ax.plot(x[ix], model_values[ix], 'ro' if disc_flag else 'go')  # red for discontinuity, green for regular
        fig.canvas.draw()
        click_counter += 1
        output_area.value = '\n'.join([
            f"Index: {i}, Value: {v:.2f}, Disc: {d}" for i, v, d in selected_points
        ])

        pending_selection = None
        selection_controls.layout.visibility = 'hidden'

        if click_counter >= max_points:
            ax.set_title("Selection complete.")
            fig.canvas.draw()
            with final_output:
                clear_output()
                print("Final selected points:")
                for i, v, d in selected_points:
                    print(f"Index: {i}, Value: {v:.2f}, Discontinuity: {d}")

    def reject_selection(b):
        nonlocal pending_selection
        pending_selection = None
        fig.canvas.draw()
        selection_controls.layout.visibility = 'hidden'

    fig.canvas.mpl_connect('button_press_event', on_click)
    confirm_button.on_click(confirm_selection)
    reject_button.on_click(reject_selection)

    return selected_points
def AddScatteringLayerTopCrust(Hscatt, dF, epsilon, kappa, a, nu):
    """
    Add a scattering layer to the top of the crustal model.
    Adjust the underlying crust layer so that total model thickness is preserved.
    """
    # Define new column names
    newParamNames = ["Depth", "Vp", "Vs", "Density", "Qkappa", "Qmu", "epsilon", "kappa", "a", "nu"]
    NoEpsilon, NoKappa, NoA, NoNu = 1E-3, 0.8, 1E4, 0.2

    # Fill missing values
    dF['epsilon'] = dF['epsilon'].fillna(NoEpsilon)
    dF['kappa']   = dF['kappa'].fillna(NoKappa)
    dF['a']       = dF['a'].fillna(NoA)
    dF['nu']      = dF['nu'].fillna(NoNu)

    # Extract crustal layer properties (from 2nd row = index 1)
    Vp = dF.loc[1, 'Vp']
    Vs = dF.loc[1, 'Vs']
    rho = dF.loc[1, 'Density']
    Qmu = dF.loc[1, 'Qmu']
    Qkappa = dF.loc[1, 'Qkappa']

    # Original depth of the first crustal layer
    crust_thickness = dF.loc[1, 'Depth']

    # Check if Hscatt is less than crust thickness
    if Hscatt >= crust_thickness:
        raise ValueError("Hscatt must be smaller than the original first crustal layer thickness")

    # Adjust first crustal layer depth by subtracting Hscatt
    
    #dF.loc[1, 'Depth'] = crust_thickness + (crust_thickness - Hscatt)
    #dF.loc[2, 'Depth'] = crust_thickness + (crust_thickness - Hscatt)

    # Create scattering layer: two points — top and bottom
    scattering_layer = pd.DataFrame({
        'Depth': [0, Hscatt],
        'Vp': [Vp/2, Vp/2],
        'Vs': [Vs/2, Vs/2],
        'Density': [2*rho/3, 2*rho/3],
        'Qkappa': [Qkappa, Qkappa],
        'Qmu': [Qmu, Qmu],
        'epsilon': [epsilon, epsilon],
        'kappa': [kappa, kappa],
        'a': [a, a],
        'nu': [nu, nu]
    })

    # Concatenate scattering layer and adjusted original model
    new_dF = pd.concat([scattering_layer, dF[1:]], ignore_index=True)

    # Ensure columns are ordered
    new_dF = new_dF[newParamNames]

    return new_dF


def GetLayerDepth(model,layer):
    """
    Function to get the depth of a layer
    """
    if layer == 'Mantle':
        return model.cmb_depth
    elif layer == 'OuterCore':
        try:
            return model.icb_iocb_depth        
        except AttributeError:
            print("Inner Core boundary depth not found in model.")
            return None
    else:
        raise ValueError("Unknown layer")
    
def AddScatteringLayerTop(model,layerDepth, Hscatt, dF, epsilon, kappa, a, nu):
    """
    Add a scattering layer to the top of the crustal model.
    :param layer: Layer number (0 for crustal layer)
    :param Hscatt: Scattering layer thickness in km
    :param dF: DataFrame containing the crustal model
    :param epsilon: Epsilon parameter for the scattering layer
    :param kappa: Kappa parameter for the scattering layer
    :param a: a parameter for the scattering layer
    :param nu: nu parameter for the scattering layer
    :return: DataFrame with the new model including the scattering layer
    """
    NoEpsilon, NoKappa, NoA, NoNu = 1E-3,0.8, 1E4, 0.2
    # Parse the NoValues in the original dataframe
    dF['epsilon'] = dF['epsilon'].replace(0, NoEpsilon)
    dF['kappa'] = dF['kappa'].replace(0, NoKappa)
    dF['a'] = dF['a'].replace(0, NoA)
    dF['nu'] = dF['nu'].replace(0, NoNu)
    # parameters
    newParamNames = ["Depth", "Vp", "Vs", "Density", "Qkappa", "Qmu", "epsilon", "kappa", "a", "nu"]

    new_dF = dF.copy()
    # identify mantle layer index
    index = np.where(new_dF['Depth'] == layerDepth)[0][0]

    index_previous_layer = index - 1
    # Copy parts before, create scattering layer, then copy parts after
    before_df = new_dF.iloc[:index].copy()

    # two new rows for scattering layer
    scattering_layer = new_dF.iloc[index:index+2].copy()
    
    # If the depth match with an existing depth, just use the existing depths
    if scattering_layer['Depth'].iloc[0] - Hscatt == new_dF['Depth'].iloc[index_previous_layer]:
        print("Layer already exists at the specified depth. Adding scattering properties")
        new_dF['epsilon'][index_previous_layer -1:index_previous_layer+1] = epsilon
        new_dF['kappa'][index_previous_layer-1:index_previous_layer+1] = kappa
        new_dF['a'][index_previous_layer-1:index_previous_layer+1] = a
        new_dF['nu'][index_previous_layer-1:index_previous_layer+1] = nu
        return new_dF
    
    else:
            
        scattering_layer['Depth'] -= Hscatt
        scattering_layer['epsilon'] = epsilon
        scattering_layer['kappa'] = kappa
        scattering_layer['a'] = a
        scattering_layer['nu'] = nu

        # after the scattering layer
        after_df = new_dF.iloc[index:].copy()

        # concatenate them together
        new_dF2 = pd.concat([before_df, scattering_layer, after_df], ignore_index=True)

        new_dF2 = new_dF2[newParamNames] # Ensure columns are ordered

        return new_dF2
    
    
    




def getQmuQkappa_PREM_AK(PREM_FilePath, AK_FilePath):
    """
    Reads the PREM and AK files and returns the Qmu and Qkappa values.
    """
        

    colNamesPREM = ['Radius','Depth','Density','Vpv','Vph','Vsv','Vsh','eta','Qmu','Qkappa']
    colNamesAK = ['Depth','Density','Vp','Vs','Qkappa','Qmu']

    colors = ['-or','-ob']
    labels = ['AK','PREM']
    paths = [AK_FilePath, PREM_FilePath]
    types = [r'\mu',r'\kappa']

    # Depths for the color plots
    depth_layers = {
        "Crust": (0, 35),  # Crust extends from surface to ~35 km
        "Mantle": (35, 2891),  # Mantle extends from 35 km to ~2891 km
        "Outer Core": (2891, 5153),  # OC from 2891 to ~5153 km
        "Inner Core": (5153, 6371)  # IC from 5153 to Earth's center (6371 km)
    }

    # Colors
    layer_colors = {
        "Crust": "lightgray",
        "Mantle": "lightcoral",
        "Outer Core": "gold",
        "Inner Core": "lightblue"
    }


    mean_Qmu_values = {}
    mean_Qkappa_values = {}
    fig,(ax1,ax2)=plt.subplots(1,2,figsize=(10,8))

    # Read parse and compute mean Q values
    for nb, cols in enumerate([colNamesAK, colNamesPREM]):
        mParams = pd.read_csv(paths[nb], delimiter=',', header=None, names=cols)

        # Initialize model dictionary
        ModelDict = {key: [] for key in cols}

        # Extract depth, Qmu, Qkappa
        depth_values = mParams['Depth'].to_numpy()
        Qmu = mParams['Qmu'].to_numpy()
        Qkappa = mParams['Qkappa'].to_numpy()

        # Compute mean Q values per layer
        model_Qmu_means = {}
        model_Qkappa_means = {}

        for layer, (top, bottom) in depth_layers.items():
            mask = (depth_values >= top) & (depth_values < bottom)
            if np.any(mask):  # Ensure there's data in the range
                model_Qmu_means[layer] = np.mean(Qmu[mask])
                model_Qkappa_means[layer] = np.mean(Qkappa[mask])

        mean_Qmu_values[labels[nb]] = model_Qmu_means
        mean_Qkappa_values[labels[nb]] = model_Qkappa_means

        # Append the first value
        for name in cols:
            mParam = mParams[name].to_numpy()
            ModelDict[name].append(mParam[0])

        # Initialize node counter
        nNodes = 0

        # Detect discontinuities
        for i in range(1, len(depth_values)):
            if depth_values[i] == depth_values[i - 1]:  # Discontinuity detected
                for name in cols:
                    mParam = mParams[name].to_numpy()
                    ModelDict[name].append(mParam[i - 1])
                    ModelDict[name].append(mParam[i])
                    if name == 'Depth':
                        nNodes += 1

        # Append the last value
        for name in cols:
            mParam = mParams[name].to_numpy()
            ModelDict[name].append(mParam[-1])

        ax1.plot(ModelDict['Qmu'], ModelDict['Depth'], colors[nb], ms=5, label=labels[nb])
        ax2.plot(ModelDict['Qkappa'], ModelDict['Depth'], colors[nb], ms=5, label=labels[nb])

    # Background colors and mean values
    for nb, ax in enumerate([ax1, ax2]):
        for layer, (top, bottom) in depth_layers.items():
            ax.axhspan(top, bottom, color=layer_colors[layer], alpha=0.3)

        # Print mean Q values on the figure
        for i, (model, valuesQmu) in enumerate(mean_Qmu_values.items()):
            valuesQkappa = mean_Qkappa_values[model]
            for j, (layer, mean_Qmu) in enumerate(valuesQmu.items()):
                mean_Qkappa = valuesQkappa[layer]
                layer_mid = (depth_layers[layer][0] + depth_layers[layer][1]) / 2

                if nb == 0:
                    ax.text(400, layer_mid + 250 * i, 
                            f"{layer} ({model}): $Q_{{\\mu}}={mean_Qmu:.1f}$",
                            fontsize=10, ha='center', bbox=dict(facecolor='white', alpha=0.7))
                else:
                    ax.text(30000, layer_mid + 250 * i, 
                            f"{layer} ({model}): $Q_{{\\kappa}}={mean_Qkappa:.1f}$",
                            fontsize=10, ha='center', bbox=dict(facecolor='white', alpha=0.7))

        ax.set_title(f'$Q_{types[nb]}$')
        ax.set_ylim([6371 + 500, -200])  # Earth's depth limit
        ax.set_ylabel('Depth (km)')
        ax.set_xlabel(f'$Q_{types[nb]}$')
        ax.legend()
    ax1.set_xlim(-50, 800)
    plt.subplots_adjust()
    plt.show()
    which_layers = {
        'Crust': 'PREM',
        'Mantle': 'PREM',
        'Outer Core': 'AK',
        'Inner Core': 'PREM'
    }
    Qmu = {}
    Qkappa = {}

    for layer, model in which_layers.items():
        Qmu[layer] = mean_Qmu_values[model][layer]
        Qkappa[layer] = mean_Qkappa_values[model][layer]
    
    
    return mean_Qmu_values, mean_Qkappa_values
