#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 31 14:22:52 2024

@author: idariash
"""

import pyart
import numpy as np

def find_closest_gate(radar, target_lat, target_lon, target_alt):
    """Find the closest radar gate to a given latitude, longitude, and altitude."""
    radar_lat = radar.latitude['data'][0]
    radar_lon = radar.longitude['data'][0]
    radar_alt = radar.altitude['data'][0]

    # Convert latitude, longitude, and altitude to Cartesian coordinates
    target_x, target_y = pyart.core.geographic_to_cartesian_aeqd(
        target_lon, target_lat, radar_lon, radar_lat)
    target_z = target_alt

    # Compute distance to each gate
    X = radar.gate_x['data']
    Y = radar.gate_y['data']
    Z = radar.gate_z['data'] + radar_alt

    distances = np.sqrt((X - target_x)**2 + (Y - target_y)**2 + (Z - target_z)**2)

    # Find the index of the closest gate
    closest_gate_index = np.unravel_index(np.argmin(distances, axis=None), distances.shape)

    return closest_gate_index

def main():
    # Replace 'your_radar_file.nc' with the path to your radar file
    radar_file = '/net/k2/storage/people/idariash/home/CSU/Scattering/data/Chill/CHL20000623_001514.uf'
    radar = pyart.io.read(radar_file)

    # Replace target_lat, target_lon, and target_alt with the desired coordinates
    target_lat = 40
    target_lon = -102.32
    target_alt = 5000.0  # Altitude in meters

    # Find the closest gate
    closest_gate_index = find_closest_gate(radar, target_lat, target_lon, target_alt)

    # Print information about the closest gate
    print(f"Closest gate index: {closest_gate_index}")
    print(f"Reflectivity: {radar.fields['reflectivity']['data'][closest_gate_index]}")
    # print(f"Latitude: {radar.latitude['data'][closest_gate_index]}")
    # print(f"Longitude: {radar.longitude['data'][closest_gate_index]}")
    # print(f"Altitude: {radar.altitude['data'][closest_gate_index]} meters")

if __name__ == "__main__":
    main()
