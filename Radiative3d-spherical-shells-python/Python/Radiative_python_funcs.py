import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import cartopy.crs as ccrs
import cartopy.feature as cfeature
#from geopy.distance import geodesic
np.float_ = np.float64 # Solve the issue. Also possible to pip install "numpy<2"
from obspy.taup import plot_ray_paths
from obspy.taup import TauPyModel
from obspy import Stream, Trace, UTCDateTime
from obspy import read, Stream
import math
from obspy.taup import plot_travel_times
import os
import importlib.util
import subprocess
import sys
import pickle
import json
import multiprocessing as mp
from multiprocessing import cpu_count
from matplotlib.colors import PowerNorm


# Browse '###' to navigate through the code
# The code is structured as follows:
# 1. Envelope / Map Plotting
# 2. Helper Functions
# 3. Travel Time curves
# 4. User security check ( Imported libraries and dependancies)


### Envelope / Map Plotting ##

def plot_taup(metadata,depth,earth_radius_km=6371.0,plot=True):
    'Plots the taup raypath for a given epicentral distance in the soruce-receiver plane'

    #Get the carthesian coordinates
    SeisLoc = metadata['SeisLocation']
    EventLoc=metadata['EventLocation']

    #Compute the range
    _,_,Range = distance_to_coords(SeisLoc,EventLoc)
    #Convert range to epicentral distance in degree

    dist = (Range / earth_radius_km) * (180 / np.pi)

    # taup
    model = TauPyModel(model="prem")
    range_seis = metadata['SeisLocation'][0] # km
    if plot==True:
        fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
        arrivals = model.get_ray_paths(depth, dist,phase_list=["P","S"])
        arrivals.plot_rays(ax=ax)
        ax.set_title(f'Ray path geometry for epicentral distance {np.round(dist,3)} °')
        return arrivals
    else:
        arrivals = model.get_ray_paths(depth, dist,phase_list=["P","S"])
        return arrivals
    
    return arrivals


def plot_geometry_2D(metadata):
    """
    Visualizes the source-receiver geometry on a 2D map using Cartopy.
    
    Parameters:
    - metadata: dict, contains the source and receiver location data.
    """
    # Unpack metadata
    EventLoc = metadata['EventLocation']  # [range, azimuth, elevation]
    SeisLoc = metadata['SeisLocation']  # Receiver location in Cartesian
    # Convert the source location to lat/lon using range and azimuth
    source_lat,source_lon=[0,0]
    receiver_lat,receiver_lon,_ = distance_to_coords(SeisLoc,EventLoc)  # range and azimuth
    # Map projection
    projection = ccrs.PlateCarree()
    fig, ax = plt.subplots(figsize=(10, 6), subplot_kw={'projection': projection})
    
    # Add features to the map
    ax.add_feature(cfeature.LAND, edgecolor='black', color='lightgray')
    ax.add_feature(cfeature.OCEAN, color='lightblue')
    ax.add_feature(cfeature.COASTLINE)
    ax.add_feature(cfeature.BORDERS, linestyle=':')
    ax.gridlines(draw_labels=True)
    
    # Plot the source and receiver
    ax.scatter(source_lon, source_lat, color='red', label='Source', zorder=5, s=100)
    ax.scatter(receiver_lon, receiver_lat, color='blue', label='Receiver', zorder=5, s=100)
    
    # Plot the great-circle path
    ax.plot(
        [source_lon, receiver_lon],
        [source_lat, receiver_lat],
        color='orange',
        linestyle='--',
        label='Great Circle Path',
        transform=ccrs.Geodetic()
    )
    ax.set_extent([-120, 120, -60,60])
    # Add legend and title
    ax.legend(loc='upper left')
    ax.set_title("Source-Receiver Geometry on a 2D Map", fontsize=14)
    
    plt.tight_layout()
    plt.show()

def seisplot(metadata, gpow=1):
    """
    Makes an envelope plot from the trace data in `tracefile`.
    Plots X, Y, Z, P, and S components with optional gamma scaling, shaded area under the envelope,
    and a vertical line at the time of maximum summed energy.

    Parameters:
    - tracefile: str, path to the trace data file.
    - gpow: float, gamma scaling factor : 1 for linear, 2 for square root... The unit is square amplitude ~ Energy density, so gamma=1 is default, and gamma=0.5 gives the ray amplitudes.
    """
    # Load the seismic data
    fig,ax=plt.subplots(figsize=(11,8))


    TraceXYZ = metadata['EnergyXYZ']
    TracePS = metadata['EnergyPS']
    CountPS = metadata['CountPS']
    AxisDesc = metadata['AxisDesc']
    Axis1 = AxisDesc[0]
    Axis2 = AxisDesc[1]
    Axis3 = AxisDesc[2]

    # Compute the time axis
    dt = metadata['binSize']
    time_axis = metadata['TimeAxis']
    NumPhon = metadata['NumPhonon']
    freq = metadata['Frequency']
    x_min, x_max = time_axis[0], time_axis[-1]

    # Gamma scaling
    TraceXYZ = TraceXYZ ** (gpow / 2)
    TracePS = TracePS ** (gpow / 2)

    # Normalization
    norm_detected = max(np.max(TraceXYZ), np.max(TracePS)) ** (gpow / 2)
    norm_goal = 190
    norm_factor = norm_goal / norm_detected

    TraceXYZ *= norm_factor
    TracePS *= norm_factor
    EEX, EEY, EEZ = TraceXYZ[0], TraceXYZ[1], TraceXYZ[2]
    EEP, EES = TracePS[0], TracePS[1]

    # Plot setup
    staging = 200
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(0.05 * staging, 6.35 * staging)
    ax.set_xlabel("Time (s)")
    ax.set_yticks(np.array([5, 4, 3, 2, 1]) * staging)
    ax.set_yticklabels([Axis1, Axis2, Axis3, "P", "S"])

    # Plot traces
    baseline = np.zeros_like(time_axis)
    colors = {
        Axis1: [0.85, 0.00, 0.00],
        Axis2: [0.00, 0.60, 0.00],
        Axis3: [0.00, 0.10, 0.75],
        "P": [0.90, 0.00, 0.60],
        "S": [0.00, 0.70, 0.90],
    }
    envelopes = {Axis1: EEX, Axis2: EEY, Axis3: EEZ, "P": EEP, "S": EES}
    baselines = [5, 4, 3, 2, 1]

    for label, baseline_level in zip(envelopes.keys(), baselines):
        baseline = baseline + baseline_level * staging
        ax.plot(time_axis, baseline + envelopes[label], linewidth=0.5, color='k')
        sumEdt = np.sum(envelopes[label]*dt) # Sum of Energy
        
        #print()
        ax.text(x=400,y=baseline[0] + 100 ,s=r"$\int_{t_0}^{t_f} E_{\mathrm{"+label+"}}(t) dt$ = " f' {sumEdt:.2f}',fontsize=13) # Add the sum of energy to the plot
        # Shading the area under the envelope
        ax.fill_between(time_axis, baseline, baseline + envelopes[label], color=colors[label], alpha=0.3)

        baseline = np.zeros_like(time_axis)

    # Compute the summed energy at each time step (sum of squares of each component)
    energy_sum = np.sum(np.array([EEX, EEY, EEZ, EEP, EES]) ** 2, axis=0)

    # Find the time index where the energy sum is maximized
    max_energy_index = np.argmax(energy_sum)

    # Add a vertical line at the time of maximum summed energy
    max_energy_time = time_axis[max_energy_index]
    ax.axvline(max_energy_time, color='r', linestyle='--', linewidth=0.8,)
    ax.text(max_energy_time +10, 6*staging , f'Max energy time : {max_energy_time:.1f} seconds ', color='k', fontsize=11, rotation=0)

    SeisLoc = metadata['SeisLocation']
    EventLoc = metadata['EventLocation']

    azimuth, Range, distance_degrees = carthesian_2_azi(SeisLoc, EventLoc)

    # Add titles and annotations
    if azimuth < 0:
        azimuth += 360
    range_km = np.sqrt(np.sum((metadata['SeisLocation'] - metadata['EventLocation']) ** 2))
    title = f"Location: {range_km:.0f} km at {azimuth:.1f}° Azimuth from source located at {-metadata['EventLocation'][2]:.0f} km depth"
    ax.set_title(title)

    # print informations about the simulation : Frequency, Number of Phonons, % of event isotropy, Total phonon caught
    IsoRatio = metadata['IsoPercentage']
    
    ax.text(x=300,y=150,s=f'Frequency: {freq} Hz',fontsize=11)
    ax.text(x=300,y=120 ,s=f'Number of Phonons: {format_large_number(int(NumPhon))}',fontsize=11)
    ax.text(x=300,y=90 ,s=f'Isotropy Ratio: {IsoRatio:.2f} %',fontsize=11)
    ax.text(x=300,y=60 ,s=f'Total Phonon Caught: {format_large_number(np.sum(CountPS))}',fontsize=11)
    plt.tight_layout()
    #Make the directory if not created
    os.makedirs('Figures',exist_ok=True)
    
    return


import os
import numpy as np
import matplotlib.pyplot as plt
from obspy import read
from obspy.core.utcdatetime import UTCDateTime

def seisplot(streamfile, gpow=1):
    """
    Plot seismic data stored in a pickled ObsPy Stream.
    Assumes traces are labeled by channels 'ER', 'ET', 'EZ', 'EP', 'ES', 'CP', 'CS'.

    Parameters:
    - streamfile: str, path to the stream pickle file.
    - gpow: float, gamma scaling factor.
    """

    # Load the seismic Stream from pickle
    st = read(streamfile)

    # Extract traces by channel
    data = {}
    for ch in ['ER', 'ET', 'EZ', 'EP', 'ES', 'CP', 'CS']:
        data[ch] = st.select(channel=ch)[0].data

    # Metadata — same for all traces
    stats = st[0].stats
    dt = stats.delta
    npts = stats.npts
    time_axis = np.arange(0, npts * dt, dt)
    NumPhon = int(stats.NumPhonon)
    freq = stats.Frequency
    AxisDesc = stats.AxisDesc
    Axis1, Axis2, Axis3 = AxisDesc[0], AxisDesc[1], AxisDesc[2]

    # Gamma scaling
    for key in ['ER', 'ET', 'EZ', 'EP', 'ES']:
        data[key] = data[key] ** (gpow / 2)

    # Normalization
    norm_detected = max([np.max(data[key]) for key in ['ER', 'ET', 'EZ', 'EP', 'ES']])
    norm_goal = 190
    norm_factor = norm_goal / norm_detected
    for key in ['ER', 'ET', 'EZ', 'EP', 'ES']:
        data[key] *= norm_factor

    EEX, EEY, EEZ = data['ER'], data['ET'], data['EZ']
    EEP, EES = data['EP'], data['ES']
    CountPS = data['CP'] + data['CS']

    # Plot setup
    staging = 200
    fig, ax = plt.subplots(figsize=(11, 8))
    ax.set_xlim(time_axis[0], time_axis[-1])
    ax.set_ylim(0.05 * staging, 6.35 * staging)
    ax.set_xlabel("Time (s)")
    ax.set_yticks(np.array([5, 4, 3, 2, 1]) * staging)
    ax.set_yticklabels([Axis1, Axis2, Axis3, "P", "S"])

    colors = {
        Axis1: [0.85, 0.00, 0.00],
        Axis2: [0.00, 0.60, 0.00],
        Axis3: [0.00, 0.10, 0.75],
        "P": [0.90, 0.00, 0.60],
        "S": [0.00, 0.70, 0.90],
    }
    envelopes = {Axis1: EEX, Axis2: EEY, Axis3: EEZ, "P": EEP, "S": EES}
    baselines = [5, 4, 3, 2, 1]

    # Plot envelopes and integrals
    for label, baseline_level in zip(envelopes.keys(), baselines):
        baseline = np.zeros_like(time_axis) + baseline_level * staging
        ax.plot(time_axis, baseline + envelopes[label], linewidth=0.5, color='k')
        sumEdt = np.sum(envelopes[label] * dt)
        ax.text(x=400, y=baseline[0] + 100, s=r"$\int E_{\mathrm{" + label + "}} dt$ = " f' {sumEdt:.2f}', fontsize=13)
        ax.fill_between(time_axis, baseline, baseline + envelopes[label], color=colors[label], alpha=0.3)

    # Summed energy
    energy_sum = np.sum(np.array([EEX, EEY, EEZ, EEP, EES]) ** 2, axis=0)
    max_energy_time = time_axis[np.argmax(energy_sum)]
    ax.axvline(max_energy_time, color='r', linestyle='--', linewidth=0.8)
    ax.text(max_energy_time + 10, 6 * staging, f'Max energy time : {max_energy_time:.1f} s', color='k', fontsize=11, rotation=0)

    # Compute azimuth etc.
    azimuth, Range, distance_degrees = carthesian_2_azi(stats.SeisLocation, stats.EventLocation)
    if azimuth < 0:
        azimuth += 360
    range_km = np.sqrt(np.sum((stats.SeisLocation - stats.EventLocation) ** 2))
    title = f"Location: {range_km:.0f} km at {azimuth:.1f}° Azimuth from source at {-stats.EventLocation[2]:.0f} km depth"
    ax.set_title(title)

    # Add annotations
    IsoRatio = stats.IsoPercentage
    ax.text(x=300, y=150, s=f'Frequency: {freq} Hz', fontsize=11)
    ax.text(x=300, y=120, s=f'Number of Phonons: {format_large_number(NumPhon)}', fontsize=11)
    ax.text(x=300, y=90, s=f'Isotropy Ratio: {IsoRatio:.2f} %', fontsize=11)
    ax.text(x=300, y=60, s=f'Total Phonon Caught: {format_large_number(int(np.sum(CountPS)))}', fontsize=11)

    plt.tight_layout()
    os.makedirs('Figures', exist_ok=True)

    return

### Helper Functions ###
def distance_to_coords(SeisLoc,EventLoc):
    # Convert range to angular distance (in degrees)
    azimuth,Range,distance_deg = carthesian_2_azi(SeisLoc,EventLoc)
    
    # Latitude is based on the distance traveled along the great circle
    latitude = distance_deg * np.cos(np.deg2rad(azimuth))
    
    # Longitude is based on the azimuth and the distance traveled along the great circle
    longitude = distance_deg * np.sin(np.deg2rad(azimuth))
    
    return latitude, longitude,Range

def read_asc_dat(file):
    parent_dir=os.path.dirname(os.path.dirname(file)) ## Go two directory before
    out_mparam_path = f'{parent_dir}/stdout.txt'

    metadata = {}
    trace_data = []
    begin_trace_idx = None  # Initialize to None


    with open(out_mparam_path,'r') as f:
        lines=f.readlines()
        for idx,line in enumerate(lines):
            if line.startswith('Number of Phonons:'): # Get the number of Phonon
                metadata['NumPhonon'] = line.split()[3]
    
            elif line.startswith('TOA_Degree:'): # TOA degree
                metadata['TOA_Degree'] = line.split()[1]

            elif line.startswith('Time to Live:'): # Phonon Time To live
                metadata['TimeToLive'] = line.split()[3]

            elif 'Seismometers Requested.' in line: # nSeismometers
                metadata['nSeismometer'] = line.split()[0]

    with open(file, 'r') as f:
        lines = f.readlines()

    for idx, line in enumerate(lines):
        line = line.strip()
        if not line:
            continue
        # Get the radius size for surface normalization
        if "Gather Radius (Outer)" in line:
            rad_P = float(line.split()[5])
            rad_S = float(line.split()[8])
            metadata['rad_P'] = rad_P
            metadata['rad_S'] = rad_S

        elif "Bin-size is" in line:
            
            metadata["binSize"] = float(line.split()[3])
            metadata["numBin"] = int(line.split()[5])
        elif line.startswith('SEIS: AxisDesc:'):
            metadata["AxisDesc"] = line.split()[2]

        elif line.startswith('SEIS: EventLocation'):
            metadata["EventLocation"] = np.array([float(line.split()[2]),float(line.split()[3]),float(line.split()[4])])

        elif line.startswith('SEIS: Seismometer Location:'):
            metadata["SeisLocation"] = np.array([float(line.split()[6].rstrip(',')),float(line.split()[7].rstrip(',')),float(line.split()[8].rstrip(','))])
        elif line.startswith('SEIS: MomentTensor: xx'):
            metadata['MTx']= (float(line.split()[3]), #xx
                                float(line.split()[4]), #xy
                                float(line.split()[5])) #xz
        elif line.startswith('SEIS: yy'):
            metadata['MTy']= (float(line.split()[2]), #yx
                                float(line.split()[3]), #yy
                                float(line.split()[4])) #yz

        elif line.startswith('SEIS: zz'):
            metadata['MTz']= (float(line.split()[2]), #zx
                                float(line.split()[3]), #zy
                                float(line.split()[4])) #zz
        elif line.startswith('SEIS: Frequency'):
            metadata['Frequency'] = line.split()[2]
        

        ## Start Parsing the Data

        elif line.startswith("SEIS: #### BEGIN TRACE ####"):
            idx+=1 # Remove the SEIS : ## BEGIN TRACE LINE 
            nlinebrowsed=0
            EnergyX =[]
            EnergyY =[]
            EnergyZ =[]
            EnergyP =[]
            EnergyS =[]
            CountP =[]
            CountS =[]
            for line_idx in range(0,metadata['numBin']): # Browse all the bins individually
                data = lines[line_idx + idx].split() # Get the whole data
                EnergyX.append(float(data[0]))
                EnergyY.append(float(data[1]))
                EnergyZ.append(float(data[2]))
                EnergyP.append(float(data[3]))
                EnergyS.append(float(data[4]))
                CountP.append(float(data[5]))
                CountS.append(float(data[6]))

                    
            
    # Convert trace data to numpy array    
    metadata["EnergyXYZ"] = np.array([EnergyX,EnergyY,EnergyZ])
    metadata["EnergyPS"] = np.array([EnergyP,EnergyS])
    metadata["CountPS"]=np.array([CountP,CountS])

    metadata['TimeAxis'] = np.arange(0,metadata['binSize']*metadata['numBin'],metadata['binSize'])
    metadata['EventMT'] = np.array([metadata['MTx'],metadata['MTy'],metadata['MTz']])

    # Extracting the moment tensor components from the metadata
    MTx = np.array(metadata['MTx'])
    MTy = np.array(metadata['MTy'])
    MTz = np.array(metadata['MTz'])
    
    # Constructing the full moment tensor (3x3 matrix)
    MT = np.array([
        [MTx[0], MTx[1], MTx[2]],
        [MTy[0], MTy[1], MTy[2]],
        [MTz[0], MTz[1], MTz[2]]
    ])
    
    # Calculate the trace (sum of diagonal elements)
    tr = np.trace(MT)
    # Explosive (isotropic) part of the moment tensor (only the trace part contributes to isotropy)
    explpart = (tr ** 2) / 3
    # Calculate the total magnitude squared of the whole tensor
    total = np.sum(MT * MT)
    # Isotropic fraction
    expl = explpart / total
    # Adjust sign of expl if trace is negative
    if tr < 0:
        expl *= -1
    # Calculate Iso% (percentage of isotropy)
    metadata["IsoPercentage"] = expl * 100  # Convert to percentage

    return metadata

def merge_traces(Directory, SeisName, NCore):
    """
    Directory : Base directory of the batch simulation
    SeisName : Name of the seismic MSEED file (e.g., 'seis_0.mseed')
    NCore : Number of cores in the simulation
    """
    merged_stream = None
    total_phonons = 0

    for Core in range(NCore):
        seis_path = f"{Directory}/process_{Core}/seisfiles/{SeisName}"
        if not os.path.exists(seis_path):
            continue
        
        st = read(seis_path)

        # Initialize with first stream
        if merged_stream is None:
            merged_stream = st.copy()
            # Initialize phonon count
            total_phonons += int(st[0].stats.NumPhonon)
        else:
            # Add data from each trace
            for tr_new in st:
                tr_existing = merged_stream.select(channel=tr_new.stats.channel)[0]
                tr_existing.data += tr_new.data
            # Update phonon count
            total_phonons += int(st[0].stats.NumPhonon)

        # Optionally remove individual file after merging
        os.remove(seis_path)

    if merged_stream is None:
        print("No traces found to merge.")
        return 1

    # Update phonon count metadata in all traces
    for tr in merged_stream:
        tr.stats.NumPhonon = total_phonons

    # Write final merged MSEED
    output_dir = f"{Directory}/seisfiles"
    os.makedirs(output_dir, exist_ok=True)
    merged_stream.write(f"{output_dir}/{SeisName}", format="PICKLE")

    return 0

def dat_to_pickle(file):
    parent_dir = os.path.dirname(os.path.dirname(file))
    out_mparam_path = f'{parent_dir}/stdout.txt'

    station_number = file.split('/')[-1].split('_')[1]  # Extract station number from filename
    metadata = {}

    # Read global metadata
    with open(out_mparam_path, 'r') as f:
        lines = f.readlines()
        for line in lines:
            if line.startswith('Number of Phonons:'):
                metadata['NumPhonon'] = line.split()[3]
            elif line.startswith('TOA_Degree:'):
                metadata['TOA_Degree'] = line.split()[1]
            elif line.startswith('Time to Live:'):
                metadata['TimeToLive'] = line.split()[3]
            elif 'Seismometers Requested.' in line:
                metadata['nSeismometer'] = line.split()[0]

    # Read data and trace-specific metadata
    with open(file, 'r') as f:
        lines = f.readlines()

    for idx, line in enumerate(lines):
        line = line.strip()
        if not line:
            continue
        if "Gather Radius (Outer)" in line:
            metadata['rad_P'] = float(line.split()[5])
            metadata['rad_S'] = float(line.split()[8])

        elif "Bin-size is" in line:
            metadata["binSize"] = float(line.split()[3])
            metadata["numBin"] = int(line.split()[5])

        elif line.startswith('SEIS: AxisDesc:'):
            metadata["AxisDesc"] = line.split()[2]

        elif line.startswith('SEIS: EventLocation'):
            metadata["EventLocation"] = np.array([float(line.split()[2]), float(line.split()[3]), float(line.split()[4])])

        elif line.startswith('SEIS: Seismometer Location:'):
            metadata["SeisLocation"] = np.array([float(line.split()[6].rstrip(',')),
                                                 float(line.split()[7].rstrip(',')),
                                                 float(line.split()[8].rstrip(','))])

        elif line.startswith('SEIS: MomentTensor: xx'):
            metadata['MTx'] = tuple(map(float, line.split()[3:6]))

        elif line.startswith('SEIS: yy'):
            metadata['MTy'] = tuple(map(float, line.split()[2:5]))

        elif line.startswith('SEIS: zz'):
            metadata['MTz'] = tuple(map(float, line.split()[2:5]))

        elif line.startswith('SEIS: Frequency'):
            metadata['Frequency'] = line.split()[2]

        elif line.startswith("SEIS: #### BEGIN TRACE ####"):
            axis1 = metadata["AxisDesc"][0]
            axis2 = metadata["AxisDesc"][1]
            axis3 = metadata["AxisDesc"][2]
            idx += 1
            numBin = metadata['numBin']
            trace_data = {k: [] for k in [f'E{axis1}',f'E{axis2}',f'E{axis3}','EP','ES','CP','CS']}
            for line_idx in range(numBin):
                data = lines[line_idx + idx].split()
                trace_data[f'E{axis1}'].append(float(data[0]))
                trace_data[f'E{axis2}'].append(float(data[1]))
                trace_data[f'E{axis3}'].append(float(data[2]))
                trace_data['EP'].append(float(data[3]))
                trace_data['ES'].append(float(data[4]))
                trace_data['CP'].append(float(data[5]))
                trace_data['CS'].append(float(data[6]))
            break

    # Construct moment tensor and isotropy
    MT = np.array([metadata['MTx'], metadata['MTy'], metadata['MTz']])
    tr_sum = np.trace(MT)
    explpart = (tr_sum ** 2) / 3
    total = np.sum(MT * MT)
    expl = (explpart / total) * (-1 if tr_sum < 0 else 1)
    metadata["IsoPercentage"] = expl * 100

    # Build Stream
    st = Stream()
    npts = metadata['numBin']
    sampling_rate = 1.0 / metadata['binSize']
    starttime = UTCDateTime(0)

    for channel, data in trace_data.items():
        tr = Trace(data=np.array(data, dtype=np.float32))
        tr.stats.station = station_number
        tr.stats.network = "R3D"
        tr.stats.channel = channel
        tr.stats.starttime = starttime
        tr.stats.sampling_rate = sampling_rate

        # Attach relevant metadata to this trace's stats
        tr.stats.rad_P = metadata['rad_P']
        tr.stats.rad_S = metadata['rad_S']
        tr.stats.EventLocation = metadata['EventLocation']
        tr.stats.SeisLocation = metadata['SeisLocation']
        tr.stats.NumPhonon = metadata['NumPhonon']
        tr.stats.IsoPercentage = metadata["IsoPercentage"]
        tr.stats.MomentTensor = MT.tolist()
        tr.stats.Frequency = metadata['Frequency']
        tr.stats.AxisDesc = metadata['AxisDesc']
        tr.stats.TimeToLive = metadata['TimeToLive']
        tr.stats.TOA_Degree = metadata['TOA_Degree']

        st.append(tr)

    directory = os.path.dirname(file)
    with open(f'{directory}/seis_{station_number}.pkl', 'wb') as f:
        pickle.dump(st, f)

    os.remove(file)
    return st





def get_azi_range_km(metadata):
    SeisLoc = metadata['SeisLocation'][0:2] # Only XY
    SourceLoc = metadata['EventLocation'][0:2] # Only XY
    delta = SeisLoc - SourceLoc
    range_km = np.round(np.sqrt(delta.dot(delta)),2) # km
    azimuth = math.atan2(delta[0],delta[1]) * (180/np.pi) # FROM SOURCE TO STATION
    if azimuth <0:azimuth+=360 # degree
    return range_km,azimuth

def format_large_number(number):
    if number >= 1e9:  # 1 billion
        return f"{number / 1e9:0.1f} B"
    else:
        return f"{number / 1e6:0.1f} M"

def carthesian_2_azi(Seisloc,EventLoc, earth_radius_km=6371.0):
    """
    Converts range (in kilometers) to angular distance (in degrees) on the sphere.
    
    Parameters:
    - range_km: float, the range in kilometers.
    - earth_radius_km: float, the radius of the sphere (default: 6371 km, Earth's average radius).
    
    Returns:
    - float: Angular distance in degrees.
    """
    # Calculate Relative distance

    RSeisLoc = Seisloc - EventLoc
    azimuth = math.atan2(RSeisLoc[0],RSeisLoc[1]) * 180/np.pi
    if azimuth<0:
        azimuth+=360

    Range= np.sqrt(RSeisLoc[0]**2+RSeisLoc[1]**2)
    distance_degrees = (Range / earth_radius_km) * (180 / np.pi)
    return azimuth,Range,distance_degrees

### Travel Time curves ##
def moving_average(arr, window_size):
    """
    Applies a simple moving average to a 1D array.
    
    Parameters:
    - arr: 1D numpy array of data to be smoothed.
    - window_size: The size of the moving window to average over.
    
    Returns:
    - Smoothed 1D numpy array.
    """
    if window_size <= 1:
        return arr  # No smoothing if window size is 1 or less
    
    # Create a moving average filter (a simple average kernel)
    kernel = np.ones(window_size) / window_size
    
    # Apply the moving average filter using convolution
    smoothed = np.convolve(arr, kernel, mode='same')
    
    return smoothed
def overlayplot(Y0, dY, X, Y, R, caption, tsize=7, lwidth=2, ax=None):
    """
    Draws an overlay on the given axis.
    
    Parameters:
      Y0      : Fraction of the vertical view where the overlay’s baseline is placed 
                (e.g. 0.8 means 80% up from the bottom).
      dY      : Fraction of the view height that defines the overlay’s peak-to-baseline scaling.
      X       : 1D array of x-values.
      Y       : 1D array of the data series (e.g. energy) to be overlaid.
                The series will be scaled so that its peak fits within [baseline-dY, baseline+dY].
      R       : (Optional) 1D array for a reference curve (pass None or an empty array if not used).
      caption : Either a single caption string or a list/tuple of up to three captions:
                [main caption, upper-right caption, lower-right caption].
      tsize   : Text size for the captions.
      lwidth  : Base line width for plotting the overlay curves.
      ax      : Matplotlib axis to plot on. If None, uses the current axis.
    """
    if ax is None:
        ax = plt.gca()

    # Process captions: if multiple captions are provided, ensure we have three.
    if isinstance(caption, (list, tuple)):
        caption1 = caption[0] if len(caption) >= 1 else ""
        caption2 = caption[1] if len(caption) >= 2 else ""
        caption3 = caption[2] if len(caption) >= 3 else ""
    else:
        caption1 = caption
        caption2 = ""
        caption3 = ""

    # Get current axis limits and compute view dimensions.
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    viewwidth = xmax - xmin
    viewheight = ymax - ymin

    # Determine the baseline (plot_y0) and scaling (plot_dy) for the overlay.
    plot_y0 = ymin + Y0 * viewheight
    plot_dy = dY * viewheight

    # For caption positioning.
    rightmarg = xmin + 0.96 * viewwidth
    topmarg = ymin + 0.96 * viewheight
    botmarg = ymin + 0.05 * viewheight

    # Scale the overlay trace.
    # Find the maximum absolute value of Y to determine the scaling factor.
    absY = np.abs(Y)
    PeakY = np.max(absY)
    iPY = np.argmax(absY)
    
    # Flat baseline at plot_y0.
    Yb = np.full_like(Y, plot_y0)
    # Scale Y so that its values vary by plot_dy about the baseline.
    YY = plot_y0 + plot_dy * (Y / PeakY)
    
    # Similarly, if a reference curve R is provided, scale it.
    if R is not None and len(R) > 0:
        R = np.asarray(R)
        RR = plot_y0 + plot_dy * (R / PeakY)
    else:
        RR = None

    # Use a slightly thinner line for auxiliary elements.
    lwidth_th = lwidth / 1.5
    tag = "OLCurve"

    # Plot the reference curve as a dotted line, if available.
    if RR is not None:
        ax.plot(X, RR, color="black", linestyle=":", linewidth=lwidth_th, label=tag)
    # Plot the baseline.
    ax.plot(X, Yb, color="black", linewidth=lwidth, label=tag)
    # Plot the scaled data series.
    ax.plot(X, YY, color="black", linewidth=lwidth, label=tag)
    # Mark the peak with a vertical dotted line.
    ax.plot([X[iPY], X[iPY]], [plot_y0, YY[iPY]], color="black", linestyle=":", 
            linewidth=lwidth_th, label=tag)

    # Add the captions.
    ax.text(X[0], plot_y0, caption1, fontsize=tsize,
            horizontalalignment="left", verticalalignment="top", color="black")
    ax.text(rightmarg, topmarg, caption2, fontsize=tsize,
            horizontalalignment="right", verticalalignment="top", color="black")
    ax.text(rightmarg, botmarg, caption3, fontsize=tsize,
            horizontalalignment="right", verticalalignment="bottom", color="black")    
    
    ax.set_xlim([xmin, xmax])
    ax.set_ylim([ymin, ymax])

    return


def plot_ttcurve(base_dir,array, ax, i, title, gamma=2, norm=0.3, theoretical=True,R=6371):
    """
    Plots a Time-Distance (TD) curve from a seismic array.
    Normalizes the energy traces with gamma scaling and a weighted area/peak scheme.
    Optionally overlays theoretical travel-time curves and an energy summary overlay.
    
    Parameters:
      array       : Dictionary with array data and metadata.
      ax          : Matplotlib axis to plot on.
      i           : Plot index (for naming, etc.).
      title       : Title for the plot.
      gamma       : Gamma scaling factor.
      norm        : Normalization parameter (0.0 => area-based; 1.0 => peak-based).
      theoretical : If True, overlays theoretical travel-time curves.
    """
    # Reshape the array to its matrix form.
    ArraySummed = np.array(array['EnergyXYZ'])  # (num_seismometers x num_samples)
    # Energy sum by seismometer (with bin size in seconds).
    EE = np.sum(ArraySummed,axis=0) * array['BinSize']

    # Gamma scaling.
    # Initialize distance array (convert from radians to degrees using Earth's radius).
    X0 = array['Distances'][0] / ((np.pi / 180) * R)  # in degrees
    X1 = array['Distances'][-1] / ((np.pi / 180) * R)  # in degrees
    X = np.linspace(X0, X1, ArraySummed.shape[1])
    # Initialize time array (convert seconds to minutes).
    T = array['TimeAxis'] / 60
    T0 = T[0]
    T1 = T[-1]

    # Plot the image.
    im = ax.imshow(
        ArraySummed.T,  # Transposed to have time on vertical axis
        cmap='hot_r',
        extent=[X0, X1, T0, T1],
        aspect='auto',
        origin='lower',
        norm=PowerNorm(gamma=gamma)
    )
    cbar = plt.colorbar(im, ax=ax)
    # Set the color limits.
    im.set_clim(0, 1)
    #cbar.set_label('Energy (Normalized)', fontsize=12)
    #cbar.ax.tick_params()
    #remove tick labels
    #cbar.set_ticks([0,0.1,0.2,0.3])  # This removes the ticks
    #cbar.ax.set_yticklabels([''] * len(cbar.get_ticks()))

    cbar.ax.tick_params(
        direction='in',  # Place the ticks inside the colorbar
        length=5,        # Length of the ticks
        width=1,         # Width of the ticks
        colors='black',  # Color of the ticks
        grid_color='black', # Color of grid lines, if any
        grid_alpha=0.5,  # Transparency of grid lines
    )

    

    # Overlay theoretical travel-time curves if requested.
    if theoretical:
#        with open(f"{base_dir}/travel_time_curves.json", "r") as f:
        # Plot the theoretical travel-time curves.
        plot_theoretical_curves(base_dir,ax)
    
    # Set the theoretical travel-time curves.
    PP_curve = get_travel_time_curve('PP')
    
    Xmax = np.max(array['Distances']) / ((np.pi / 180) * R)  # in degrees
    # Set plot labels and limits.
    ax.set_title(title, fontsize=14)
    ax.set_xlabel(r'$\Delta$ (°)', fontsize=14)
    ax.set_ylabel("Time (mins)", fontsize=14)
    ax.set_ylim([T0, T1])
    ax.set_xlim([X0, Xmax])

    # --- Overlay Plot Section ---
    # Select a subset of the energy series for the overlay.
    CutIn = max(int(np.ceil(len(EE) / 10)), 2)
    # Python indices: subtract one to match 0-indexing.
    EEsubset = EE[CutIn - 1:]
    Xsubset = X[CutIn - 1:]
    PeakE = np.max(EEsubset)
    iPE = np.argmax(EEsubset)
    TermE = EEsubset[-1]
    # Placeholder for event isotropy fraction; replace with getexplfrac if available.
    explfrac = array["IsoPercentage"] / 100
    NPhonCast = array.get('NumPhonon', -1)
    # Optionally, extra reference caption if norm is a structure (here, legacy norm so left empty).
    RefCaption = ""

    # If these variables are strings, try to convert them to float (if applicable):
    PeakE = float(PeakE) if isinstance(PeakE, str) else PeakE
    TermE = float(TermE) if isinstance(TermE, str) else TermE
    NPhonCast = float(NPhonCast) if isinstance(NPhonCast, str) else NPhonCast
    explfrac = float(explfrac) if isinstance(explfrac, str) else explfrac
    frequency = float(array.get('Frequency', 0)) if isinstance(array.get('Frequency', 0), str) else array.get('Frequency', 0)


    caption_str = (
        f'Sum(E*dt); Max:{PeakE:8.1e} at {Xsubset[iPE]:.0f}' r'$\degree$'
        #f"\n           Term: {TermE:8.2e} at {Xsubset[-1]:1.0f}  degree;"
        f"\n\nPhonons Cast: {format_large_number(NPhonCast)}"
        f"\nEvent Isotropy: {explfrac * 100:0.0f}%"
        f"\nFrequency: {frequency:0.2f} Hz"
        f"{RefCaption}"
    )
    # Overlay parameters.
    OLTextSize = 8.5
    OLY0 = 0.85   # Vertical position (fraction of view height)
    OLdY =0.15   # Vertical scaling (fraction of view height)
    CaptUR = ""
    CaptLR = ""
    OLCscale = 0.5 # Exponent for scaling (e.g. 0.5 for square-root scaling)
    # Apply scaling to the energy subset.
    EEoverlay = EEsubset ** OLCscale

    EEoverlay = moving_average(EEoverlay,800)

    # Call the overlay plotting function.
    overlayplot(OLY0, OLdY, Xsubset, EEoverlay, None,
                caption=[caption_str, CaptUR, CaptLR],
                tsize=OLTextSize, lwidth=1, ax=ax)

    return
def get_azi_range_km(metadata):
    SeisLoc = metadata['SeisLocation'][0:2] # Only XY
    SourceLoc = metadata['EventLocation'][0:2] # Only XY
    delta = SeisLoc - SourceLoc
    range_km = np.round(np.sqrt(delta.dot(delta)),0) # km
    azimuth = math.atan2(delta[0],delta[1]) * (180/np.pi) # FROM SOURCE TO STATION
    if azimuth <0:azimuth+=360 # degree
    return range_km,azimuth

def plot_theoretical_curves(base_dir,ax):
    """
    Plots the theoretical travel-time curves for a given seismic phase.
    
    Parameters:
      phase          : The seismic phase to plot (e.g., 'P', 'S', etc.).
    """
    
    #Load travel_time_curves.json
    with open(f"{base_dir}/travel_time_curves.json", "r") as f:
        travel_time_data = json.load(f)
        
    for phase in travel_time_data:
        if phase not in travel_time_data:
            print(f"Phase {phase} not found in travel_time_data.")
            continue
        if phase == 'SKSSKS':continue
        if len(travel_time_data[phase]) == 0:
            print(f"No data for phase {phase}.")
            continue
        xdata, ydata = zip(*travel_time_data[phase])
        ax.plot(xdata, ydata,'k-',linewidth=0.1,alpha=0.8)

    return

def get_travel_time_curve(phase):
    fig_func, ax_func = plt.subplots(figsize=(9, 9))  # Creates a figure

    plot_travel_times(
        source_depth=10,
        phase_list=[phase],
        ax=ax_func,  # Use existing axis
        fig=fig_func,  # Use created figure
        legend=False,
        verbose=True,
        plot_all=False,
        show=False  # Prevents display
    )

    travel_time_data = []
    for line in ax_func.get_lines():
        x_data = line.get_xdata()
        y_data = line.get_ydata()
        travel_time_data.append(np.column_stack((x_data, y_data)))

    travel_time_array = np.vstack(travel_time_data)

    plt.close(fig_func)  # Close the figure to prevent extra plots
    return travel_time_array
def GetDistanceFromTravelTime(PP_curve,time):
    " Get the distance from the travel time curve"
    " PP_curve : Travel time curve"
    " time : Travel time in minutes"

    " Return : Distance in degrees"
    times = PP_curve[:,1] # Travel time in minutes
    idx = (np.abs(times - time)).argmin()
    return PP_curve[idx,0]




def sprintf(nseis,pattern='%03d'):
    filename = f'seis_{pattern}' % (nseis)
    return filename

def get_azi_range_km_from_coords(seisloc, eventloc):
    """Return range in km and azimuth from seismometer to event."""
    delta_xyz = np.array(seisloc) - np.array(eventloc)
    range_km = np.linalg.norm(delta_xyz)
    azimuth = np.degrees(np.arctan2(delta_xyz[1], delta_xyz[0]))
    return range_km, azimuth

def assembleArray(fpath, ibegin, iend):
    array = {
        'NumSeismometer': 0,
        'Distances': [],
        'Azimuths': [],
        'EnergyXYZ': []
    }

    for i in range(ibegin, iend+1):
        filename = f"{fpath}/seis_{i:03d}.pkl"
        if not os.path.exists(filename):
            continue

        with open(filename, 'rb') as f:
            st = pickle.load(f)

        # Get distances & azimuth
        tr0 = st[0]  # one trace for metadata
        seisloc = tr0.stats.SeisLocation
        eventloc = tr0.stats.EventLocation
        IsoPercentage = tr0.stats.IsoPercentage
        array['IsoPercentage'] = IsoPercentage
        range_km,azimuth = get_azi_range_km_from_coords(seisloc, eventloc)
        array['NumSeismometer'] += 1


        # Sum EnergyXYZ components
        ex = st.select(channel="ER")[0].data
        ey = st.select(channel="ET")[0].data
        ez = st.select(channel="EZ")[0].data
        total_energy_xyz = ex + ey + ez
        range_km = np.round(range_km, 0)  # Here I round the range so that values from different seismometers can be merged. This is important !
        if range_km in array['Distances']: # Check if the distance is already in the array

            # If the distance is already in the array, add the energy to the existing seismometer
            # Find the index of the duplicate distance
        
            idx = array['Distances'].index(range_km)
            # Add the energy to the existing seismometer
            array['EnergyXYZ'][idx] += total_energy_xyz
        else:
            # Add the new distance and energy
            array['Distances'].append(range_km)
            array['EnergyXYZ'].append(total_energy_xyz)
            array['Azimuths'].append(azimuth)
            # Add the energy to the existing seismometer

        if i == ibegin:
            array['TimeAxis'] = np.arange(len(ex)) * (1.0 / st[0].stats.sampling_rate)
            array['NumBin'] = len(ex)
            array['Frequency'] = tr0.stats.Frequency
            array['EventMT'] = tr0.stats.MomentTensor
            array['TOA_degree'] = tr0.stats.TOA_Degree
            array['NumPhonon'] = tr0.stats.NumPhonon
            array['TimeToLive'] = tr0.stats.TimeToLive
            array['BinSize'] = 1.0 / st[0].stats.sampling_rate
            array['MaxTime'] = array['TimeAxis'][-1]

    return array

### User security check ##

def check_and_install(package):
    """Check if a package is installed, and prompt the user to install it or create a new virtual environment."""
    if importlib.util.find_spec(package) is None:
        print(f"\n  The package '{package}' is missing.")

        # Ask the user what they want to do
        choice = input(f"Do you want to: \n"
                       f"Install '{package}' in the current environment \n"
                       f"2Create a new virtual environment and install it there \n"
                       f"Skip installation \n"
                       f"Enter 1, 2, or 3: ").strip()

        if choice == "1":
            print(f"Installing {package} in the current environment...")
            subprocess.check_call([sys.executable, "-m", "pip", "install", package])

        elif choice == "2":
            venv_name = input("Enter a name for the new virtual environment: ").strip()
            if not venv_name:
                venv_name = "new_env"

            print(f"\nCreating a virtual environment '{venv_name}' and installing {package}...\n")

            # Create a virtual environment
            subprocess.check_call([sys.executable, "-m", "venv", venv_name])

            # Determine the correct pip path inside the venv
            if sys.platform == "win32":
                pip_path = os.path.join(venv_name, "Scripts", "pip.exe")
                activate_cmd = f"{venv_name}\\Scripts\\activate"
            else:
                pip_path = os.path.join(venv_name, "bin", "pip")
                activate_cmd = f"source {venv_name}/bin/activate"

            # Install the package in the new environment
            subprocess.check_call([pip_path, "install", package])

            print(f"\n Virtual environment '{venv_name}' created.")
            print(f"⚡ To activate it, run: \n  `{activate_cmd}`\n")

        else:
            print(f" Skipping installation of '{package}'. Some functionality may be missing.")

    else:
        print(f" {package} is already installed.")



###GRID READING ##

def plot_earth_layers(grid,ax,Radius_km,colormap,**kwargs):
    """
    Plots the radial layers of the Earth model in a planar XY view.

    Parameters:
        grid (ndarray): The Earth model grid with shape (1, 1, n_layers, 3),
                        where each layer is defined as [x, y, z].
        colormap (ndarray): An optional Nx3 matrix to color the layers.
    """
    # Extract the radial coordinates (z-values)
    radii = grid[0, 0, :, 2]
    # Prepare the plot
    ax.set_aspect('equal')
    ax.set_xlim([-Radius_km, Radius_km])
    ax.set_ylim([-Radius_km, Radius_km])

    # Plot each layer as a circle
    for i, radius in enumerate(radii):
        # Define the color if colormap is provided
        if colormap is not None and i < len(colormap):
            color = colormap[i]
        else:
            color = 'r'  # Default color

        circle = plt.Circle((0, 0), radius, edgecolor='black', facecolor='white', alpha=0.1,lw=1,**kwargs)
        ax.add_artist(circle)

    # Label the plot
    ax.set_title("Grid used in the simulation", fontsize=14)
    ax.set_xlabel("X (km)")
    ax.set_ylabel("Y (km)")
    
    return ax
    
def read_gridgeom(gridfile):
    # Load the raw grid data from the file
    GG_Raw = np.loadtxt(gridfile, skiprows=8)

    # Initialize base index (assuming it is 0, as in the original code)
    ibase = 0

    # Extract grid dimensions (NX, NY, NZ) based on the maximum values in each coordinate
    NX = int(np.max(GG_Raw[:, 0]) - ibase + 1)
    NY = int(np.max(GG_Raw[:, 1]) - ibase + 1)
    NZ = int(np.max(GG_Raw[:, 2]) - ibase + 1)
    # Initialize the GG matrix (shape: NX x NY x NZ x 3)
    GG = np.zeros((NX, NY, NZ, 3))

    # Loop through each line of the raw grid and pack the values into the GG matrix
    for j in range(GG_Raw.shape[0]):
        ix = int(GG_Raw[j, 0] - ibase)
        iy = int(GG_Raw[j, 1] - ibase)
        iz = int(GG_Raw[j, 2] - ibase)
        
        # Assign the X, Y, Z coordinates to the GG matrix
        GG[ix, iy, iz, 0] = GG_Raw[j, 3]  # X
        GG[ix, iy, iz, 1] = GG_Raw[j, 4]  # Y
        GG[ix, iy, iz, 2] = GG_Raw[j, 5]  # Z

    return GG

def read_params(file_path):
    params = {}  # Dictionary to store the parameters
    
    # Open the file and read line by line
    with open(file_path, 'r') as file:
        for line in file:
            # Remove leading/trailing whitespaces
            line = line.strip()
            
            # Skip empty lines or lines starting with a comment
            if not line or line.startswith('#'):
                continue
            
            # Split the line into key and value based on '='
            key, value = line.split('=')
            
            # Convert the value to the appropriate type (int, float, or string)
            if value.isdigit():  # If it's a digit, convert to int
                value = int(value)
            elif value.replace('.', '', 1).isdigit():  # If it's a float
                value = float(value)
            
            # Add to the dictionary
            params[key] = value
    
    return params
