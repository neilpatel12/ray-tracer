#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  4 11:03:48 2020

@author: neilpatel
"""
# from RayTracerProjectModular 
import ray_tracer_tools as tool
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

class Ray:
    '''
    A class used to represent a ray that can be propagated through optical 
    elements. 

    The positions and directions of the ray at each interaction with an optical 
    element is recorded in hidden attributes that can be accessed using access 
    methods. The best desrciption of the ray's motionis returned by access 
    method get_path() as both position and direction history is recorded.
        
    Methods
    -------
    
    get_path()
        Returns 3D array, containing position and direction history as ray 
        propagates through a system of optical elements.
        
    vertices()
        Returns 2D numpy array of ray's positions.
        array([p0, p1, p2, ..., pm])
        
    directions()
        Returns 2D numpy array of ray's directions.
        array([k0, k1, k2, ..., km])
        
    p()
        Returns 1D numpy array of length 3 representing current ray position in
        3D space
        
    k()
        Returns normalised 1D numpy array of length 3 representing current ray 
        direction in 3D space
        
    is_alive()
        Returns boolean value representing 
        
    set_p(new_p)
        overwrites ray's last known position, _p, with new_p
    
    set_k(new_k)
        overwrites ray's last known direction, _k, with new_k
        
    append(new_p, new_k)
        appends new position and direction to 3D array containing point and 
        direction history.
         
        
    
    Notes
    -----
        Hidden Attributes
        -----------------
        _p : 1D numpy array of length 3
            Current position of ray in 3D space
            
        _k : 1D numpy array of length 3
            Current direction of ray in 3D space
            
        _path : 3D array, length depends on how many objects ray propagates through
                Contains position and direction history as ray propagates through a 
                system of optical elements
                
                array([[p0, k0], [p1, k1], [p2, k2], ..., [pm, km]]) 
                , where pn and kn are numpy arrays containing point position and 
                direction.
                
         _alive : bool
                 False if ray has been terminated (by _terminate method) due to 
                 poor propagation through optical element, e.g. due to Total 
                 Internal Reflection or ray missing optical element.
                 True if ray propagates optical element successfully.
                 Value is used to filter ray_list in RayBundle instance to consider 
                 behaviour of remaining alive rays in the bundle.

    '''

    
    def __init__(self, p = [0, 0, 0], k = [0, 0, 1]):   

        '''    
        Parameters
        ----------
        p : list object with length 3, optional
            Represents initial ray position in 3D space
           
        k : list object with length 3, optional
            Represents initial ray direction in 3D space
    
        '''                                             
        # EXCEPTIONS AND ASSIGNMENTS
# =============================================================================
#         if type(p) == list and len(p) == 3 :
#             self._p = np.array(p, dtype = float)
#         else:
#             raise Exception('p must be type list with length 3')
#              
#         if type(k) == list and len(k) == 3 :   
#             self._k = tool.normalise_vector(k)
#         else:
#             raise Exception('k must be type list with length 3')
# =============================================================================
            
            
        self._p = np.array(p, dtype = float)
        self._k = tool.normalise_vector(k)
        self._path = np.array([[self._p, self._k]])
        self._alive = True
            
        
    def __repr__(self):
        return f'Ray(p = {self._p}, k = {self._k}), with ray path {self._path}'
    
    def __str__(self):
        return f'Ray(p = {self._p}, k = {self._k})'

# ACCESS METHODS
    def get_path(self):
        return self._path
    
    def vertices(self):
        return self._path[:,0,:]
    
    def directions(self):
        return self._path[:,1,:]
    
    def p(self):
        return self._path[:,0,:][-1] # or self._p
    
    def k(self):
        return self._path[:,1,:][-1] # or self._k
    
    def is_alive(self):
        return self._alive
    
# MODIFIER METHODS

    def set_p(self, new_p):
        self._p = np.array(new_p, dtype = float) 
     
    def set_k(self, new_k):
        self._k = tool.normalise_vector(new_k)
     
# UPDATING METHODS (for when ray hits surface causing position or direction changes)  
    def append(self, new_p, new_k):
        self.set_p(new_p)
        self.set_k(new_k)
        self._path = np.append(self._path, [[self._p, self._k]], axis = 0) # updates old path by appending most recent point an direction

# TERMINATE METHOD
    def _terminate(self):
        print(f'{str(self)} has been terminated')
        self._alive = False  # can be used to filter when we have list of rays
        
        
    
class RayBundle:
    '''
    A class used to represent a ray bundle. A ray bundle is a collection of 
    rays that are evenly spaced to achieve uniform density within the bundle.
    The rays are evenly spaced by setting the initial positions of rays in the 
    bundle to be fixed on a grid, centred on a specific position. The rays 
    travel in the same direction.
     
    The rays that form the bundle propagate independently of each other and 
    do not interfere, but the overall behaviour of the bundle provides insight
    into how a beam of light may act when passed through an optical system.
 
   
    Attributes
    ----------
    p_centre : list object with length 3, optional
        Central ray's initial position in 3D space   
                
    k : list object with length 3, optional
        Initial direction in 3D space of all rays in the bundle
        
    bundle_radius : float, optional
        Radius of ray bundle
        
    xy_spacing : float, optional
        Distance between horizontally and vertically adjacent rays in ray 
        bundle
        
        
    Methods
    -------
    get_paths_list()
        Returns list of 3D arrays which contain position and direction
        history of a ray in the ray bundle
        
    get_vertices_list()
        Returns list of 2D arrays which contain position
        history of a ray in the ray bundle
        
    get_xs()
        Returns list of 1D arrays which contain x position
        history of a ray in the ray bundle
        
    get_ys()
        Returns list of 1D arrays which contain y position
        history of a ray in the ray bundle
        
    get_zs()
        Returns list of 1D arrays which contain z position
        history of a ray in the ray bundle
        
    get_directions_list()
        Returns list of 2D arrays which contain direction
        history of a ray in the ray bundle
        
    get_ps_list()
        Returns list of 1D arrays which contain current position in 3D space 
        of a ray in the ray bundle
        
    get_ks_list()
        Returns list of 1D arrays which contain current direction in 3D space 
        of a ray in the ray bundle
    
    RMS_dev()
        Returns Root Mean Square deviation of ray bundle using the deviation of
        the current position of each ray in the bundle from the optical axis
        
    num_alive_rays()
        Returns number of rays in the bundle that have succesfully propagated 
        through all of the optical elements in the system so far
        
    spot_diagram(color = 'b', Mean = False)
        Produces 2D spot diagram of current ray positions.
        If Mean is True then mean average ray position is plotted on diagram
           
    plot_paths_3D()
        Produces 3D plot of the rays that have successfully propagated through
        the entire optical system. Rays that do not successfully propagate 
        through any one of the optical elements in the optical system are not 
        plotted.
        
    plot_2D_projection()
        Produces 2D plot of the rays that have successfully propagated through
        the entire optical system from side view.  Rays that do not 
        successfully propagate through any one of the optical elements in the 
        optical system are not plotted.

    '''
    
    def __init__(self, p_centre = [0, 0, 0], k = [0, 0, 1], bundle_radius = 25, xy_spacing = 1):
        self.p_centre = np.array(p_centre, dtype = float)
        self.k = tool.normalise_vector(k)
        self.bundle_radius = bundle_radius
        self.xy_spacing = xy_spacing
        self.ray_list = []
        for x in np.arange(- (self.bundle_radius + self.xy_spacing) , self.bundle_radius + 2 * self.xy_spacing, self.xy_spacing):
            for y in np.arange(- (self.bundle_radius + self.xy_spacing) , self.bundle_radius + 2 * self.xy_spacing, self.xy_spacing):
                if x ** 2 + y ** 2 <= self.bundle_radius ** 2:
                    grid_coord = np.array([p_centre[0] + x, p_centre[1] + y, 0])
                    ray = Ray(p = grid_coord, k = self.k)
                    (self.ray_list).append(ray)
                    
       
# REPRESENTATION     
                    
    def __repr__(self):
        return f'RayBundle(p_centre = {self.p_centre}, k = {self.k} , bundle_radius = {self.bundle_radius} , xy_spacing = {self.xy_spacing})'
    
    def __str__(self):
        return
    
# ACCESS METHODS
        
    def get_paths_list(self):   
        return np.array([ray._path for ray in self.ray_list])
            
    def get_vertices_list(self):
        return np.array([ray.vertices() for ray in self.ray_list])
    
    def get_xs(self):
        return self.get_vertices_list()[:,:, 0]
     
    def get_ys(self):
        return self.get_vertices_list()[:,:, 1]
    
    def get_zs(self):
        return self.get_vertices_list()[:,:, 2]
            
            
    def get_directions_list(self):
        return np.array([ray.get_directions() for ray in self.ray_list])
    
    def get_ps_list(self):
        return np.array([ray.p() for ray in self.ray_list])
    
    def get_ks_list(self):
        return np.array([ray.k() for ray in self.ray_list])
        
    def RMS_dev(self):
        S_list = []
        for ray in self.ray_list:
            if ray._alive:
                x = ray.p()[0] 
                y = ray.p()[1]
                S = x ** 2 + y ** 2
                S_list.append(S)
        RMS = np.sqrt(sum(S_list)/self.num_alive_rays())
        return RMS        
        
    def num_alive_rays(self):
        num_alive_rays_in_list = 0
        for ray in self.ray_list:
            if ray._alive:
                num_alive_rays_in_list += 1
        return num_alive_rays_in_list
                
    def spot_diagram(self, color = 'b', mean = False): # plot spot diagram of last interaction
        fig = plt.figure()
        ax = fig.add_subplot(111)
        x_list = []
        y_list = []
        for ray in self.ray_list:
            if ray._alive:
                x = ray.p()[0] 
                y = ray.p()[1]
                ax.scatter(x,y, color = color)
                x_list.append(x)
                y_list.append(y)
        if mean:
            x_bar = sum(x_list)/self.num_alive_rays()
            y_bar = sum(y_list)/self.num_alive_rays()
            plt.scatter(x_bar,y_bar, color = 'r')
            print(f'mean ray position is ({x_bar}, {y_bar})')
        plt.title(f'Spot Diagram of ray bundle, {self.num_alive_rays()} alive rays')
        plt.xlabel('X (mm)')
        plt.ylabel('Y (mm)')
        plt.show()
        
    def plot_paths_3D(self): #? might need to sort out axis labels
        fig = plt.figure(figsize = plt.figaspect(1)*1.5)
        ax = fig.add_subplot(111, projection='3d') 
        ax.set_xlabel('X (mm)')
        ax.set_ylabel('Z (mm)')
        ax.set_zlabel('Y (mm)')
        for ray in self.ray_list:
            if ray._alive:
                x = ray.vertices()[:, 0]
                y = ray.vertices()[:, 1]
                z = ray.vertices()[:, 2] 
                ax.plot3D(x, y, z, zdir = 'y', alpha = 0.5)
        plt.title(f'Rays with direction {self.k}')      
        plt.show() 
       
    def plot_2D_projection(self): #? might need to sort out axis labels
        fig = plt.figure(figsize = plt.figaspect(1)*1.5)
        ax = fig.add_subplot(111) 
        ax.set_xlabel('Z (mm)')
        ax.set_ylabel('Y (mm)')
        for ray in self.ray_list:
            if ray._alive:
                y = ray.vertices()[:, 1]
                z = ray.vertices()[:, 2] 
                ax.plot(z, y, '-')
        plt.title(f'Rays with direction {self.k}')      
        plt.show() 

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
        
