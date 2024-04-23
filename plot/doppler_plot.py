#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  3 18:09:09 2021

@author: idariash
"""
import scipy.ndimage as spyi
import pyart
import matplotlib.pyplot as plt
import numpy as np
import netCDF4
import matplotlib.colors as colors
from scipy.interpolate import RegularGridInterpolator


def cappi(filename_1b, filename_dd, field, ind, x0=None, y0=None, x1=None, y1=None):

    pyart_grid_1b = pyart.io.read_grid(filename_1b)
    pyart_grid_dd = pyart.io.read_grid(filename_dd)
    fig = plt.figure(figsize=(6, 5))
    
    try:
        W = pyart_grid_dd.fields['w']['data']
    except:
        W = pyart_grid_dd.fields['upward_air_velocity']['data']
    
    W = np.ma.masked_where( W < -50, W)
    Z = pyart_grid_1b.fields['corrected_reflectivity']['data']
    W = np.ma.masked_where( Z < 15, W)
    W = W[ind]
    W = spyi.gaussian_filter(W, sigma = 1)
    x, y = np.meshgrid(0.001*pyart_grid_1b.x['data'], 0.001*pyart_grid_1b.y['data']) 
    W = np.ma.masked_where( x > -20, W)
    W = np.ma.masked_where( y < -20, W)
    
    if ('REF' in field.upper()) and (not 'DIF' in field.upper()):
        field = 'Reflectivity'
        units = 'dBZ'        
        field_var = pyart_grid_1b.fields['corrected_reflectivity']['data']
        field_var = np.ma.masked_where( Z < 10, field_var)        
        cs = plt.pcolormesh(x, y, field_var[ind], vmin=0, vmax=70, cmap='pyart_NWSRef')
          
    elif 'REF' in field.upper() and 'DIF' in field.upper():
        field = 'Diff. Reflectivity'
        units = 'dB'        
        field_var = pyart_grid_1b.fields['corrected_differential_reflectivity']['data']       
        field_var = np.ma.masked_where( Z < 10, field_var)        
        cs = plt.pcolormesh(x, y, field_var[ind], vmin=-2, vmax=6, cmap=pyart.graph.cm.RefDiff)
        
    elif 'SPECIFIC' in field.upper():
        field = 'Specific Diff. Phase'
        units = 'deg/km'        
        field_var = pyart_grid_1b.fields['corrected_specific_differential_phase']['data']      
        field_var = np.ma.masked_where( Z < 10, field_var)        
        cs = plt.pcolormesh(x, y, field_var[ind], vmin=0, vmax=6, cmap=pyart.graph.cm.Theodore16)
        
    elif 'CORRELATION' in field.upper():
        field = 'Copolar Corr. Ratio'
        units = ''        
        field_var = pyart_grid_1b.fields['corrected_copol_correlation_ratio']['data']       
        field_var = np.ma.masked_where( Z < 10, field_var)        
        cs = plt.pcolormesh(x, y, field_var[ind], vmin=0.5, vmax=1, cmap='jet')
        
    plt.contour(x, y, W, levels=[1, 5, 10], colors=['k', 'k', 'k'], antialised = True)
    #plt.contour(x, y, W[ind], levels=[-100, -5], colors=['w', 'w', 'w'], antialised = True)
    #CHIVO mark
    plt.plot(0, 0, 'ok', markersize=5)
    plt.text(0, 0, ' C', fontsize=12)
    # CSAPR mark
    plt.plot(-54, -54, 'ok', markersize=5)
    plt.text(-54, -54, ' Y', fontsize=12)
    # RMA mark
    plt.plot(-10, 20, 'ok', markersize=5)
    plt.text(-10, 20, ' R', fontsize=12)
    
    if x0 is not None:
        
        plt.plot(x0, y0, 'ok', markersize=5)
        plt.plot(x1, y1, 'ok', markersize=5)
        
        plt.plot([x0, x1],[y0, y1], 'k--',)
    
    plt.xlim(-90, 40)
    plt.ylim(-90, 40)
    
    
    time_start = netCDF4.num2date(pyart_grid_1b.time['data'][0], pyart_grid_1b.time['units'])
    time_text = time_start.strftime('%Y-%m-%dT%H:%M UTC')
    plt.colorbar(cs, label= field + ' (' + units + ')')
    plt.title('CSU-CHIVO '+ field + ' at '+ str(ind + 1) +' km \n'+ time_text)
    plt.xlabel('Distance east of CHIVO (km)')
    plt.ylabel('Distance north of CHIVO (km)')
    return fig




def cross_sec(filename_1b, filename_dd, field, x0, y0, x1, y1):
    
    pyart_grid_1b = pyart.io.read_grid(filename_1b)
    pyart_grid_dd = pyart.io.read_grid(filename_dd)
    
    x = 0.001*pyart_grid_dd.y['data']
    y = 0.001*pyart_grid_dd.y['data']
    z = 0.001*pyart_grid_dd.z['data']
    
    distance = np.sqrt((x1- x0)**2 + (y1 - y0)**2)
    xi = np.linspace(x0, x1, int(distance))
    yi = np.linspace(y0, y1, int(distance))

    
    try:
        U = pyart_grid_dd.fields['u']['data']
        V = pyart_grid_dd.fields['v']['data']
        W = pyart_grid_dd.fields['w']['data']
        Z = pyart_grid_dd.fields['DT']['data']
    except:
        U = pyart_grid_dd.fields['eastward_wind']['data']
        V = pyart_grid_dd.fields['northward_wind']['data']
        W = pyart_grid_dd.fields['upward_air_velocity']['data']
        Z = pyart_grid_dd.fields['reflectivity']['data']
        
    #Z = pyart_grid_1b.fields['corrected_reflectivity']['data']
    
    
    
    U = np.ma.masked_where( Z < 0, U)
    V = np.ma.masked_where( Z < 0, V)
    
    # Substract strom motion per layer
    for i in range(len(z)):
        Um = np.ma.mean(U[i])
        Vm = np.ma.mean(V[i])
        
        U[i] = U[i] - Um
        V[i] = V[i] - Vm
    
    U = np.swapaxes(U, 0, 2)
    V = np.swapaxes(V, 0, 2)
    W = np.swapaxes(W, 0, 2)
    Z = np.swapaxes(Z, 0, 2)
    
    fn_u = RegularGridInterpolator((x, y, z), U, method = 'nearest')
    fn_v = RegularGridInterpolator((x, y, z), V, method = 'nearest')
    fn_w = RegularGridInterpolator((x, y, z), W, method = 'nearest')
    fn_z = RegularGridInterpolator((x, y, z), Z, method = 'nearest')
    
    Ui_pj = np.zeros((len(yi), len(z))) # Projection horizontal winds along the direction of the cross section
    Wi = np.zeros((len(yi), len(z)))
    Zi = np.zeros((len(yi), len(z)))
    field_var_i = np.zeros((len(yi), len(z)))
    
    for i in range(len(yi),):
        for j in range(len(z)):
            u = fn_u(np.array([xi[i], yi[i], z[j]]))
            v = fn_v(np.array([xi[i], yi[i], z[j]]))
            Ui_pj[i, j] = ((x1-x0)*u + (y1-y0)*v)/distance # Projection of a vector on another vector
            Wi[i, j] = fn_w(np.array([xi[i], yi[i], z[j]]))
            Zi[i, j] = fn_z(np.array([xi[i], yi[i], z[j]]))
            
    Ui_pj = np.ma.masked_where( Ui_pj < -100, Ui_pj)

    
    fig = plt.figure(figsize=(8, 5))
    
    if ('REF' in field.upper()) and (not 'DIF' in field.upper()):
        field = 'Reflectivity'
        units = 'dBZ'
        field_var = pyart_grid_1b.fields['corrected_reflectivity']['data']
        field_var = np.swapaxes(field_var, 0, 2)
        fn_field_var = RegularGridInterpolator((x, y, z), field_var, method = 'nearest')
        for i in range(len(yi),):
            for j in range(len(z)):
                field_var_i[i, j] = fn_field_var(np.array([xi[i], yi[i], z[j]]))
        
        Zi = np.ma.masked_where( Zi < 10, field_var_i)
        cs = plt.pcolormesh(xi, z, Zi.T, vmin=0, vmax=70, cmap='pyart_NWSRef')
    
    elif 'REF' in field.upper() and 'DIF' in field.upper():
        field = 'Diff. Reflectivity'
        units = 'dB'        
        field_var = pyart_grid_1b.fields['corrected_differential_reflectivity']['data']
        field_var = np.swapaxes(field_var, 0, 2)
        fn_field_var = RegularGridInterpolator((x, y, z), field_var, method = 'nearest')
        for i in range(len(yi),):
            for j in range(len(z)):
                field_var_i[i, j] = fn_field_var(np.array([xi[i], yi[i], z[j]]))
        
        field_var_i = np.ma.masked_where( Zi < 10, field_var_i)
        cs = plt.pcolormesh(xi, z, field_var_i.T, vmin=-2, vmax=6, cmap=pyart.graph.cm.RefDiff)
        
    elif 'SPECIFIC' in field.upper():
        field = 'Specific Diff. Phase'
        units = 'deg/km'        
        field_var = pyart_grid_1b.fields['corrected_specific_differential_phase']['data']
        field_var = np.swapaxes(field_var, 0, 2)
        fn_field_var = RegularGridInterpolator((x, y, z), field_var, method = 'nearest')
        for i in range(len(yi),):
            for j in range(len(z)):
                field_var_i[i, j] = fn_field_var(np.array([xi[i], yi[i], z[j]]))
        
        field_var_i = np.ma.masked_where( Zi < 10, field_var_i)
        cs = plt.pcolormesh(xi, z, field_var_i.T, vmin=0, vmax=6, cmap=pyart.graph.cm.Theodore16)
        
    elif 'CORRELATION' in field.upper():
        field = 'Copolar Corr. Ratio'
        units = ''        
        field_var = pyart_grid_1b.fields['corrected_copol_correlation_ratio']['data']   
        field_var = np.swapaxes(field_var, 0, 2)
        fn_field_var = RegularGridInterpolator((x, y, z), field_var, method = 'nearest')
        for i in range(len(yi),):
            for j in range(len(z)):
                field_var_i[i, j] = fn_field_var(np.array([xi[i], yi[i], z[j]]))
        
        field_var_i = np.ma.masked_where( Zi < 10, field_var_i)
        cs = plt.pcolormesh(xi, z, field_var_i.T,  vmin=0.7, vmax=1, cmap='jet')
        
    elif 'HYDROCLASS' in field.upper() or 'HYC' in field.upper():
        field = 'Hydrometeor Classification'
        units = ''        
        field_var = pyart_grid_1b.fields['HydroClass']['data']   
        field_var = np.swapaxes(field_var, 0, 2)
        fn_field_var = RegularGridInterpolator((x, y, z), field_var, method = 'nearest')
        for i in range(len(yi),):
            for j in range(len(z)):
                field_var_i[i, j] = fn_field_var(np.array([xi[i], yi[i], z[j]]))
        field_var_i = np.ma.masked_where( Zi < 10, field_var_i)
        hid_colors = ['White', 'Black','LightBlue', 'MediumBlue', 'DarkOrange', 'LightPink',
                  'Cyan', 'DarkGray', 'Lime', 'Yellow', 'Red', 'Fuchsia']
        cmaphid = colors.ListedColormap(hid_colors)
        cs = plt.pcolormesh(xi, z, field_var_i.T,  vmin=5, vmax=16, cmap=cmaphid)
    
    q = plt.quiver(xi, z, Ui_pj.T, Wi.T,  scale = 400, color='black',)

    plt.quiverkey(q, X=0.3, Y=1.1, U=10, label='Vertical Motion 10(m/s)', labelpos='E')
    
    plt.xlim(-60, 0)
    plt.ylim(0, 20)
    
    time_start = netCDF4.num2date(pyart_grid_1b.time['data'][0], pyart_grid_1b.time['units'])
    time_text = time_start.strftime('%Y-%m-%dT%H:%M UTC')
    rotated_angle = round(np.arctan((y1-y0)/(x1-x0))*180/np.pi)
    plt.colorbar(cs, label= field + ' (' + units + ')')
    plt.title('CSAPR '+ field + ' rotated '+ str(rotated_angle) +' deg '+ time_text + '\n \n \n')
    plt.xlabel('East West distance from CHIVO (km)')
    plt.ylabel('Altitude (km)')
    
    return fig

def cross_sec_NorthSouth(filename_1b, filename_dd, field, x0, y0, x1, y1):
    
    pyart_grid_1b = pyart.io.read_grid(filename_1b)
    pyart_grid_dd = pyart.io.read_grid(filename_dd)
    
    x = 0.001*pyart_grid_dd.y['data']
    y = 0.001*pyart_grid_dd.y['data']
    z = 0.001*pyart_grid_dd.z['data']
    
    distance = np.sqrt((x1- x0)**2 + (y1 - y0)**2)
    xi = np.linspace(x0, x1, int(distance))
    yi = np.linspace(y0, y1, int(distance))
    
    
    try:
        U = pyart_grid_dd.fields['u']['data']
        V = pyart_grid_dd.fields['v']['data']
        W = pyart_grid_dd.fields['w']['data']
        Z = pyart_grid_dd.fields['DT']['data']
    except:
        U = pyart_grid_dd.fields['eastward_wind']['data']
        V = pyart_grid_dd.fields['northward_wind']['data']
        W = pyart_grid_dd.fields['upward_air_velocity']['data']
        Z = pyart_grid_dd.fields['reflectivity']['data']
    
    
    U = np.ma.masked_where( Z < 0, U)
    V = np.ma.masked_where( Z < 0, U)
    W = np.ma.masked_where( Z < 0, W)
    
    #Substract strom motion per layer
    
    # for i in range(len(z)):
    #     Um = np.ma.mean(U[i])
    #     Vm = np.ma.mean(V[i])
        
    #     U[i] = U[i] - Um
    #     V[i] = V[i] - Vm
    
    U = np.swapaxes(U, 0, 2)
    V = np.swapaxes(V, 0, 2)
    W = np.swapaxes(W, 0, 2)
    Z = np.swapaxes(Z, 0, 2)
    
    
    fn_u = RegularGridInterpolator((x, y, z), U, method = 'nearest')
    fn_v = RegularGridInterpolator((x, y, z), V, method = 'nearest')
    fn_w = RegularGridInterpolator((x, y, z), W, method = 'nearest')
    fn_z = RegularGridInterpolator((x, y, z), Z, method = 'nearest')
    
    Ui_pj = np.zeros((len(yi), len(z))) # Projection horizontal winds along the direction of the cross section
    Wi = np.zeros((len(yi), len(z)))
    Zi = np.zeros((len(yi), len(z)))
    field_var_i = np.zeros((len(yi), len(z)))
    
    for i in range(len(yi),):
        for j in range(len(z)):
            u = fn_u(np.array([xi[i], yi[i], z[j]]))
            v = fn_v(np.array([xi[i], yi[i], z[j]]))
            Ui_pj[i, j] = ((x1-x0)*u + (y1-y0)*v)/distance # Projection of a vector on another vector
            Wi[i, j] = fn_w(np.array([xi[i], yi[i], z[j]]))
            Zi[i, j] = fn_z(np.array([xi[i], yi[i], z[j]]))
            
    Ui_pj = np.ma.masked_where( Ui_pj < -100, Ui_pj)

    
    fig = plt.figure(figsize=(8, 5))
    
    if ('REF' in field.upper()) and (not 'DIF' in field.upper()):
        field = 'Reflectivity'
        units = 'dBZ'
        try:
            field_var = pyart_grid_1b.fields['corrected_reflectivity']['data']
        except:
            field_var = pyart_grid_1b.fields['reflectivity']['data']
        field_var = np.swapaxes(field_var, 0, 2)
        fn_field_var = RegularGridInterpolator((x, y, z), field_var, method = 'nearest')
        for i in range(len(yi),):
            for j in range(len(z)):
                field_var_i[i, j] = fn_field_var(np.array([xi[i], yi[i], z[j]]))
        
        Zi = np.ma.masked_where( Zi < 0, field_var_i)
        cs = plt.pcolormesh(yi, z, Zi.T, vmin=0, vmax=70, cmap='pyart_NWSRef')
    
    elif 'REF' in field.upper() and 'DIF' in field.upper():
        field = 'Diff. Reflectivity'
        units = 'dB'        
        field_var = pyart_grid_1b.fields['corrected_differential_reflectivity']['data']
        field_var = np.swapaxes(field_var, 0, 2)
        fn_field_var = RegularGridInterpolator((x, y, z), field_var, method = 'nearest')
        for i in range(len(yi),):
            for j in range(len(z)):
                field_var_i[i, j] = fn_field_var(np.array([xi[i], yi[i], z[j]]))
        
        field_var_i = np.ma.masked_where( Zi < 10, field_var_i)
        cs = plt.pcolormesh(yi, z, field_var_i.T, vmin=-2, vmax=6, cmap=pyart.graph.cm.RefDiff)
        
    elif 'SPECIFIC' in field.upper():
        field = 'Specific Diff. Phase'
        units = 'deg/km'        
        field_var = pyart_grid_1b.fields['corrected_specific_differential_phase']['data']
        field_var = np.swapaxes(field_var, 0, 2)
        fn_field_var = RegularGridInterpolator((x, y, z), field_var, method = 'nearest')
        for i in range(len(yi),):
            for j in range(len(z)):
                field_var_i[i, j] = fn_field_var(np.array([xi[i], yi[i], z[j]]))
        
        field_var_i = np.ma.masked_where( Zi < 10, field_var_i)
        cs = plt.pcolormesh(yi, z, field_var_i.T, vmin=0, vmax=6, cmap=pyart.graph.cm.Theodore16)
        
    elif 'CORRELATION' in field.upper():
        field = 'Copolar Corr. Ratio'
        units = ''        
        field_var = pyart_grid_1b.fields['corrected_copol_correlation_ratio']['data']   
        field_var = np.swapaxes(field_var, 0, 2)
        fn_field_var = RegularGridInterpolator((x, y, z), field_var, method = 'nearest')
        for i in range(len(yi),):
            for j in range(len(z)):
                field_var_i[i, j] = fn_field_var(np.array([xi[i], yi[i], z[j]]))
        
        field_var_i = np.ma.masked_where( Zi < 10, field_var_i)
        cs = plt.pcolormesh(yi, z, field_var_i.T,  vmin=0.7, vmax=1, cmap='jet')
        
    elif 'HYDROCLASS' in field.upper() or 'HYC' in field.upper():
        field = 'Hydrometeor Classification'
        units = ''        
        field_var = pyart_grid_1b.fields['HydroClass']['data']   
        field_var = np.swapaxes(field_var, 0, 2)
        fn_field_var = RegularGridInterpolator((x, y, z), field_var, method = 'nearest')
        for i in range(len(yi),):
            for j in range(len(z)):
                field_var_i[i, j] = fn_field_var(np.array([xi[i], yi[i], z[j]]))
        field_var_i = np.ma.masked_where( Zi < 10, field_var_i)
        hid_colors = ['White', 'Black','LightBlue', 'MediumBlue', 'DarkOrange', 'LightPink',
                  'Cyan', 'DarkGray', 'Lime', 'Yellow', 'Red', 'Fuchsia']
        cmaphid = colors.ListedColormap(hid_colors)
        cs = plt.pcolormesh(yi, z, field_var_i.T,  vmin=5, vmax=16, cmap=cmaphid)
    
    q = plt.quiver(yi, z, Ui_pj.T, Wi.T,  scale = 400, color='black',)

    plt.quiverkey(q, X=0.3, Y=1.1, U=10, label='Vertical Motion 10(m/s)', labelpos='E')
    
    plt.xlim(-40, 20)
    plt.ylim(0, 20)
    
    time_start = netCDF4.num2date(pyart_grid_1b.time['data'][0], pyart_grid_1b.time['units'])
    time_text = time_start.strftime('%Y-%m-%dT%H:%M UTC')
    rotated_angle = 90 # round(np.arctan((y1-y0)/(x1-x0))*180/np.pi)
    plt.colorbar(cs, label= field + ' (' + units + ')')
    plt.title('CSAPR '+ field + ' rotated '+ str(rotated_angle) +' deg '+ time_text + '\n \n \n')
    plt.xlabel('East West distance from CHIVO (km)')
    plt.ylabel('Altitude (km)')
    
    return fig






def cross_sec_east(filename_1b, filename_dd, field, distance_north):
    
    ind = 151 - 50 + distance_north
    thin = 1
    pyart_grid_1b = pyart.io.read_grid(filename_1b)
    pyart_grid_dd = pyart.io.read_grid(filename_dd)
    fig = plt.figure(figsize=(8, 5))
    W = pyart_grid_dd.fields['w']['data']
    U = pyart_grid_dd.fields['u']['data']
    Z = pyart_grid_1b.fields['corrected_reflectivity']['data']
    x, z = np.meshgrid(0.001*pyart_grid_1b.x['data'], 0.001*pyart_grid_1b.z['data']) 
    
    
    U = np.swapaxes(U,0,1)
    W= np.swapaxes(W,0,1)
    Z = np.swapaxes(Z,0,1)
    
    if ('REF' in field.upper()) and (not 'DIF' in field.upper()):
        field = 'Reflectivity'
        units = 'dBZ'
        field_var = pyart_grid_1b.fields['corrected_reflectivity']['data']
        field_var = np.swapaxes(field_var,0,1)
        field_var = np.ma.masked_where( Z < 10, field_var)
        cs = plt.pcolormesh(x, z, field_var[ind], vmin=0, vmax=70, cmap='pyart_NWSRef')
    
    elif 'REF' in field.upper() and 'DIF' in field.upper():
        field = 'Diff. Reflectivity'
        units = 'dB'        
        field_var = pyart_grid_1b.fields['corrected_differential_reflectivity']['data']
        field_var = np.swapaxes(field_var,0,1)
        field_var = np.ma.masked_where( Z < 10, field_var)
        cs = plt.pcolormesh(x, z, field_var[ind], vmin=-2, vmax=6, cmap=pyart.graph.cm.RefDiff)
        
    elif 'SPECIFIC' in field.upper():
        field = 'Specific Diff. Phase'
        units = 'deg/km'        
        field_var = pyart_grid_1b.fields['corrected_specific_differential_phase']['data']
        field_var = np.swapaxes(field_var,0,1)
        field_var = np.ma.masked_where( Z < 10, field_var)
        cs = plt.pcolormesh(x, z, field_var[ind], vmin=0, vmax=6, cmap=pyart.graph.cm.Theodore16)
        
    elif 'CORRELATION' in field.upper():
        field = 'Copolar Corr. Ratio'
        units = ''        
        field_var = pyart_grid_1b.fields['corrected_copol_correlation_ratio']['data']   
        field_var = np.swapaxes(field_var,0,1)
        field_var = np.ma.masked_where( Z < 10, field_var)
        cs = plt.pcolormesh(x, z, field_var[ind],  vmin=0.7, vmax=1, cmap='jet')
        
    elif 'HYDROCLASS' in field.upper() or 'HYC' in field.upper():
        field = 'Hydrometeor Classification'
        units = ''        
        field_var = pyart_grid_1b.fields['HydroClass']['data']   
        field_var = np.swapaxes(field_var,0,1)
        field_var = np.ma.masked_where( Z < 10, field_var)
        
        hid_colors = ['White', 'Black','LightBlue', 'MediumBlue', 'DarkOrange', 'LightPink',
                  'Cyan', 'DarkGray', 'Lime', 'Yellow', 'Red', 'Fuchsia']
        cmaphid = colors.ListedColormap(hid_colors)
        
        cs = plt.pcolormesh(x, z, field_var[ind],  vmin=5, vmax=16, cmap=cmaphid)
    
    q = plt.quiver(x[::thin, ::thin], z[::thin, ::thin], 
            U[ind][::thin, ::thin], W[ind][::thin, ::thin],
            scale = 400, color='black',)

    plt.quiverkey(q, X=0.3, Y=1.1, U=10, label='Vertical Motion 10(m/s)', labelpos='E')
    
    plt.xlim(-90, -30)
    plt.ylim(0, 20)
    
    time_start = netCDF4.num2date(pyart_grid_1b.time['data'][0], pyart_grid_1b.time['units'])
    time_text = time_start.strftime('%Y-%m-%dT%H:%M UTC')
    plt.colorbar(cs, label= field + ' (' + units + ')')
    plt.title('CSU-CHIVO '+ field + ' at '+ str(distance_north) +' km '+ time_text + '\n \n \n')
    plt.xlabel('East West distance from CHIVO (km)')
    plt.ylabel('Altitude (km)')
    
    return fig

def cross_sec_north(filename_1b, filename_dd, field, distance_east):
    
    ind = 151 - 50 - distance_east
    thin = 1
    pyart_grid_1b = pyart.io.read_grid(filename_1b)
    pyart_grid_dd = pyart.io.read_grid(filename_dd)
    fig = plt.figure(figsize=(8, 5))
    W = pyart_grid_dd.fields['w']['data']
    V = pyart_grid_dd.fields['v']['data']
    Z = pyart_grid_1b.fields['corrected_reflectivity']['data']
    y, z = np.meshgrid(0.001*pyart_grid_1b.y['data'], 0.001*pyart_grid_1b.z['data']) 
    
    V = np.swapaxes(V,1,2)
    W= np.swapaxes(W,1,2)
    Z = np.swapaxes(Z,1,2)
    
    V = np.swapaxes(V,0,1)
    W= np.swapaxes(W,0,1)
    Z = np.swapaxes(Z,0,1)
    
    if ('REF' in field.upper()) and (not 'DIF' in field.upper()):
        field = 'Reflectivity'
        units = 'dBZ'
        field_var = pyart_grid_1b.fields['corrected_reflectivity']['data']
        field_var = np.swapaxes(field_var,1,2)
        field_var = np.swapaxes(field_var,0,1)
        field_var = np.ma.masked_where( Z < 10, field_var)
        cs = plt.pcolormesh(y, z, field_var[ind], vmin=0, vmax=70, cmap='pyart_NWSRef')
    
    elif 'REF' in field.upper() and 'DIF' in field.upper():
        field = 'Diff. Reflectivity'
        units = 'dB'        
        field_var = pyart_grid_1b.fields['corrected_differential_reflectivity']['data']
        field_var = np.swapaxes(field_var,1,2)
        field_var = np.swapaxes(field_var,0,1)
        field_var = np.ma.masked_where( Z < 10, field_var)
        cs = plt.pcolormesh(y, z, field_var[ind], vmin=-2, vmax=6, cmap=pyart.graph.cm.RefDiff)
        
    elif 'SPECIFIC' in field.upper():
        field = 'Specific Diff. Phase'
        units = 'deg/km'        
        field_var = pyart_grid_1b.fields['corrected_specific_differential_phase']['data']
        field_var = np.swapaxes(field_var,1,2)
        field_var = np.swapaxes(field_var,0,1)
        field_var = np.ma.masked_where( Z < 10, field_var)
        cs = plt.pcolormesh(y, z, field_var[ind], vmin=0, vmax=6, cmap=pyart.graph.cm.Theodore16)
        
    elif 'CORRELATION' in field.upper():
        field = 'Copolar Corr. Ratio'
        units = ''        
        field_var = pyart_grid_1b.fields['corrected_copol_correlation_ratio']['data']   
        field_var = np.swapaxes(field_var,1,2)
        field_var = np.swapaxes(field_var,0,1)
        field_var = np.ma.masked_where( Z < 10, field_var)
        cs = plt.pcolormesh(y, z, field_var[ind],  vmin=0.7, vmax=1, cmap='jet')
    
    q = plt.quiver(y[::thin, ::thin], z[::thin, ::thin], 
            V[ind][::thin, ::thin], W[ind][::thin, ::thin],
            scale = 400, color='black',)

    plt.quiverkey(q, X=0.3, Y=1.1, U=10, label='Vertical Motion 10(m/s)', labelpos='E')
    
    plt.xlim(0, 60)
    plt.ylim(0, 20)
    
    time_start = netCDF4.num2date(pyart_grid_1b.time['data'][0], pyart_grid_1b.time['units'])
    time_text = time_start.strftime('%Y-%m-%dT%H:%M UTC')
    plt.colorbar(cs, label= field + ' (' + units + ')')
    plt.title('CSU-CHIVO '+ field + ' at '+ str(distance_east) +' km '+ time_text + '\n \n \n')
    plt.xlabel('North South distance from CHIVO (km)')
    plt.ylabel('Altitude (km)')
    
    return fig
    

    




def cappi_ref(filename_1b, filename_dd, field, ind):
    
    pyart_grid_1b = pyart.io.read_grid(filename_1b)
    pyart_grid_dd = pyart.io.read_grid(filename_dd)
    fig = plt.figure(figsize=(6, 5))
    W = pyart_grid_dd.fields['w']['data']
    W = np.ma.masked_where( W < -50, W)
    W = W[ind]
    W = spyi.gaussian_filter(W, sigma = 1)
    Z = pyart_grid_1b.fields['corrected_reflectivity']['data']
    x, y = np.meshgrid(0.001*pyart_grid_1b.x['data'], 0.001*pyart_grid_1b.y['data']) 
    W = np.ma.masked_where( x > -20, W)
    W = np.ma.masked_where( y < -20, W)
    
    if ('REF' in field.upper()) and (not 'DIF' in field.upper()):
        field = 'Reflectivity'
        units = 'dBZ'        
        field_var = pyart_grid_1b.fields['corrected_reflectivity']['data']
        field_var = np.ma.masked_where( Z < 10, field_var)        
        cs = plt.pcolormesh(x, y, field_var[ind], vmin=0, vmax=70, cmap='pyart_NWSRef')
          
    elif 'REF' in field.upper() and 'DIF' in field.upper():
        field = 'Diff. Reflectivity'
        units = 'dB'        
        field_var = pyart_grid_1b.fields['corrected_differential_reflectivity']['data']       
        field_var = np.ma.masked_where( Z < 10, field_var)        
        cs = plt.pcolormesh(x, y, field_var[ind], vmin=-2, vmax=6, cmap=pyart.graph.cm.RefDiff)
        
    elif 'SPECIFIC' in field.upper():
        field = 'Specific Diff. Phase'
        units = 'deg/km'        
        field_var = pyart_grid_1b.fields['corrected_specific_differential_phase']['data']      
        field_var = np.ma.masked_where( Z < 10, field_var)        
        cs = plt.pcolormesh(x, y, field_var[ind], vmin=0, vmax=6, cmap=pyart.graph.cm.Theodore16)
        
    elif 'CORRELATION' in field.upper():
        field = 'Copolar Corr. Ratio'
        units = ''        
        field_var = pyart_grid_1b.fields['corrected_copol_correlation_ratio']['data']       
        field_var = np.ma.masked_where( Z < 10, field_var)        
        cs = plt.pcolormesh(x, y, field_var[ind], vmin=0.5, vmax=1, cmap='jet')
        
    plt.contour(x, y, W, levels=[1, 10,20], colors=['k', 'k', 'k'], antialised = True)
    #plt.contour(x, y, W[ind], levels=[-100, -5], colors=['w', 'w', 'w'], antialised = True)
    #CHIVO mark
    plt.plot(0, 0, 'ok', markersize=5)
    plt.text(0, 0, ' C', fontsize=12)
    # CSAPR mark
    plt.plot(-54, -54, 'ok', markersize=5)
    plt.text(-54, -54, ' Y', fontsize=12)
    # CSAPR mark
    plt.plot(-10, 20, 'ok', markersize=5)
    plt.text(-10, 20, ' R', fontsize=12)
    
    plt.xlim(-90, 40)
    plt.ylim(-90, 40)
    
    time_start = netCDF4.num2date(pyart_grid_1b.time['data'][0], pyart_grid_1b.time['units'])
    time_text = time_start.strftime('%Y-%m-%dT%H:%M UTC')
    plt.colorbar(cs, label= field + ' (' + units + ')')
    plt.title('CSU-CHIVO '+ field + ' at '+ str(ind + 1) +' km \n'+ time_text)
    plt.xlabel('Distance east of CHIVO (km)')
    plt.ylabel('Distance north of CHIVO (km)')
    return fig