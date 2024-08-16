import numpy as np


def save_beam_data_to_csv(beam, filename):
    """
    Save beam data to a CSV file with predefined column headers.

    Parameters:
    - beam : ShadowLib.Beam
        The Shadow Beam object containing the data to be saved.
    - filename : str
        The filename (including path) to save the data.

    Returns:
    None
    """
    # Fixed columns list
    cols = [11, 23, 24, 25, 1, 3, 2, 4, 6, 5]
    
    # Extract data from Beam object
    data = np.asarray(beam._beam.getshcol(cols)).T
    
    # Define column headers
    headers = [
        "energy",
        "intensity",
        "intensity_s-pol",
        "intensity_p-pol",
        "X",
        "Y",
        "Z",
        "Xp",
        "Yp",
        "Zp",
    ]
    
    # Save data to file with headers
    np.savetxt(filename, data, header=",".join(headers), fmt='%1.6e', delimiter=',', comments='')

    print(f"{filename} saved!")

# Example usage:
# Save data to 'beam_data.csv'
beam = in_object_1
save_beam_data_to_csv(beam, 'beam_data.csv')