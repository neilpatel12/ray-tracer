#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  4 11:48:02 2020

@author: neilpatel
"""

# from RayTracerProjectModular 
import ray_tracer_tools as tool
# from RayTracerProjectModular 
import ray_and_ray_bundle as rrb
import numpy as np


class OpticalElement:
    '''
    Base class to SphericalRefraction and OpticalElements classes.
    
    Forms basic structure of an optical element i.e. the position of the 
    element along the optical axis and the aperture radius of the element.
    

        
    Methods
    -------
    get_z_0()
        Returns position of optical element along optical axis
    
    get_radius()
        Returns optical element aperture radius
    
    set_z_0()
        Changes position of optical element along optical axis
        
    set_radius() 
        Changes optical element aperture radius
    
    propagate_ray()
        raises NotImplementedError if derived class does not define 
        propagate_ray method
        
    intercept()
        raises NotImplementedError if derived class does not define intercept
        method
        
    propagate_ray_bundle(ray_bundle)
        Iterates through rays in ray_bundle ray_list, applying the optical 
        elements propagate_ray method to each ray. Includes counter to keep 
        track of ray number and boolean to see when a ray successfully 
        or unsuccessfully propagates through an optical element.
    
    Notes
    -----
        Hidden Attributes
        -----------------
        _z_0 : float
            z position of optical element
            
        _r : float
            aperture radius of optical element
            

    '''
    def __init__(self, z_0, r):  
        '''
        Parameters
        ----------
        p : list object with length 3, optional
            Represents initial ray position in 3D space
           
        k : list object with length 3, optional
            Represents initial ray direction in 3D space
        '''
        self._z_0 = z_0
        self._r = r
    
# ACCESS METHODS
    
    def get_z_0(self):
        return self._z_0
    
    def get_radius(self):
        return self._r
    
# MODIFIER METHODS

    def set_z_0(self, new_z_0):
        self._z_0 = new_z_0 
     
    def set_radius(self, new_rad):
        self._r = new_rad
    
# METHODS THAT MUST BE IMPLEMENTED
    def propagate_ray(self,ray):
        raise NotImplementedError('propagate a ray through the optical element')
   
    def intercept(self, ray):
        raise NotImplementedError('optical element must have intercept method')

# LOOPS PROPAGATE_RAY()
    def propagate_ray_bundle(self, ray_bundle):
        i = 0
        for ray in ray_bundle.ray_list:
            i += 1
            print(i, f'{self}')
            print(ray)
            print('ray._alive is', ray._alive)
            self.propagate_ray(ray)
            print('ray._alive is', ray._alive)              
            print('')    
            
class SphericalRefraction(OpticalElement):     
    '''
     Derived from OpticalElement class.
     
     Represents spherically curved surface at position z_0 along the optical 
     axis with aperture radius r and a specific curvature and refractive index.
   
    Methods
    -------
    
    get_n1()
        Returns refractive index of area at z < z_0
    
    get_n2()
        Returns refractive index of area at z > z_0
    
    get_curvature()
        Returns curvature of spherical surface
    
    intercept(ray)
        Returns numpy array describing position of intercept of ray and 
        surface.
        Takes account of intercept position to ensure correct position is
        calculated for a positive lens or negative lens. Extra care is taken 
        for rays that could intercept the lens twice.
    
    snell(k_hat, n_hat, n1, n2)
        Returns numpy array describing new refracted direction of ray.
        Or returns None if there is Total Internal Reflection.
    
    propagate_ray(ray)
        Appends new position and direction of ray to 3D array containing ray's
        position and direction history. Calculates intercept position and
        new direction of ray using intercept() and snell() 
    
    parax_focus()
        Returns float describing focal length of the refracting surface
        with positive curvature, value is calculated by propagating paraxial 
        ray with very small perpendicular displacement from optical axis and
        calculating where the ray intercepts the optical axis.
        
    focus_from_eq()
        Returns float describing focal length of the refracting surface
        with positive curvature. 
        Uses approximate calculation: f = R * n2 / (n2 - n1), where R is radius 
        of surface curvature
        
    Notes
    -----
        Hidden Attributes
        -----------------
        _z_0 : float
            z position optical element
            
        _r : float
            aperture radius of optical element
            
        _curvature : float
            curvature of surface of optical element
            
        _R : float
            radius of curvature of surface
        
        _origin : float
            position of centre curvature of surface
      
        _n1 : float
            refractive index of area at z < _z_0
        
        _n2 : float
            refractive index of area at z > _z_0
    '''
    def __init__(self, z_0, curv, n1 = 1, n2 = 1, r = 100):
        OpticalElement.__init__(self, z_0, r)
        self._curvature = curv
        if self._curvature != 0: 
            self._origin = z_0 + 1/self._curvature
            self._R = 1/self._curvature
        elif self._curvature == 0:
            self._origin = self._z_0
            self._R = None
        self._n1 = n1
        self._n2 = n2
        
    def __repr__(self):
        return f'SphericalRefraction(z_0 = {self._z_0}, curv = {self._curvature}, n1 = {self._n1}, n2 = {self._n2}, r = {self._r})'
    
# ACCESS METHODS       
        
    def get_n1(self):
        return self._n1
    
    def get_n2(self):
        return self._n2
    
    def get_curvature(self):
        return self._curvature
    
# PROPAGATION TOOLS AND METHODS   
    def intercept(self, ray):   
        r_vect = ray._p - np.array([0, 0, self._origin], dtype = float) # r_vect is vector from origin of curvature to point ray._p
        if self._curvature == 0 :   # for lens with zero curvature, on z axis at z
            l = (self._origin - ray._p[2])/ ray._k[2]
        elif (np.dot(r_vect, ray._k) * np.dot(r_vect, ray._k)) - (np.dot(r_vect,r_vect) - self._R * self._R ) < 0:
            print(f'{str(ray)} does not intercept with lens surface, ray does not even intercept sphere of lens curvature radius')
            return None
        
        else:
            # Assuming rays come in from left, i.e. z=0
            l1 = - np.dot(r_vect, ray._k) + np.sqrt( (np.dot(r_vect, ray._k) * np.dot(r_vect, ray._k)) - (np.dot(r_vect, r_vect) - self._R * self._R )) 
            l2 = - np.dot(r_vect, ray._k) - np.sqrt( (np.dot(r_vect, ray._k) * np.dot(r_vect, ray._k)) - (np.dot(r_vect, r_vect) - self._R * self._R ))         
            
            if self._curvature < 0:  # prevents rays that intersect positive lens sphere twice on same side from claiming to be negative intersection 
                l_poss = max(l1,l2)
                possible_int_posn = ray._p + l_poss * ray._k
                if possible_int_posn[2] > self._origin:         # ray must intercept lens at z position greater than origin of curvature to intercept negative lens
                    l = l_poss
                else:
                    print(f'{str(ray)} intercepts with incorrect lens, i.e. intercepts with positive lens not the negative lens as required')
                    return None
            elif self._curvature > 0: # prevents rays that intersect negative lens sphere twice on same side from claiming to be positive intersection 
                l_poss = min(l1,l2) 
                possible_int_posn = ray._p + l_poss * ray._k
                if possible_int_posn[2] < self._origin:         # ray must intercept lens at z position less than origin of curvature to intercept positive lens
                    l = l_poss
                else:
                    print(f'{str(ray)} intercepts with incorrect lens, i.e. intercepts with negative lens not the positive lens as required')
                    return None
        if (ray._p + l * ray._k)[0] ** 2 + (ray._p + l * ray._k)[1] ** 2 > self._r **2 :# if resulting position is outside aperture radius then position must be excluded
            print(f'{str(ray)} does not intercept with lens, it is outside the aperture radius of the lens')  
            return None
        else:
            return ray._p + l * ray._k
        
    
    def snell(self, k_hat, n_hat, n1, n2):  
        c = - np.dot(tool.normalise_vector(k_hat), tool.normalise_vector(n_hat)) # c = cos(theta1), must be postive
        r = n1/n2
        if c < 0:
            raise TypeError('Cos(theta1) must be positive')
        if np.sqrt(1 - c ** 2) > 1/r and r > 1: #n2/n1: # Total internal refraction
            return None
        else:
            return tool.normalise_vector(r * k_hat + (r * c - np.sqrt(1 - (r ** 2) * (1 - c ** 2))) * n_hat)
     
        
    def propagate_ray(self,ray):
        
        intercept_posn = self.intercept(ray)
        if isinstance(intercept_posn, type(None)):
            ray._terminate()
        else:
            if self._curvature == 0:
                n_surf = np.array([0, 0, -1], dtype = float)
            else:
                n_surf = tool.normalise_vector((intercept_posn - np.array([0, 0, self._origin], dtype = float)) * self._curvature/abs(self._curvature)) # multiply surface vector by curvature sign to get normal to point towards light source      
            
            new_k = self.snell(k_hat = ray._k, n_hat = n_surf, n1 = self._n1, n2 = self._n2)
            
            if new_k is None:
                print(f'{str(ray)} was reflected by TIR')
                ray._terminate()
            else:
                print(f'{str(ray)} is refracted')
                ray.append(intercept_posn, new_k)   
                print(f'ray is now {str(ray)}')
       
    def parax_focus(self):
        
        if self._curvature > 0:
            ray = rrb.Ray(p = [0, 0.001, 0], k = [0,0,1])
            self.propagate_ray(ray)
            output_plane = OutputPlane(z_0 = 10000, r = 10000 ) #z_0 large to allow for very long focal lengths, r is large to allow for very short focal lengths so ray still intercepts output plane
            output_plane.propagate_ray(ray)
            focal_length = - ray.vertices()[1][1] * (ray.directions()[2][2]/ray.directions()[2][1])
            print(f'From simulation, for spherical surface with curvature = {self._curvature} mm^-1, focal length = {focal_length} mm')
            return float(focal_length)
        else:
            return None
    
    def focus_from_eq(self): # equation: f = (R * n2)/(n2 - n1)
        if self._curvature > 0:
            focal_length = ((1/self._curvature) * self._n2)/ (self._n2 - self._n1)
            print(f'From spherical surface refraction equation, for spherical surface with curvature = {self._curvature} mm^-1, focal length = {focal_length} mm')
            return float(focal_length)
        else:
            return None
        
class OutputPlane(OpticalElement):
    
    '''
     Derived from OpticalElement class.
     
     Represents planar screen at position z_0 that has radius r and is 
     perpendicular to optical axis.
     
     Rays intercept this screen and the position is recorded, but ray direction
     is not changed.
     
    Methods
    -------
    
    get_n1()
        Returns refractive index of area at z < z_0
    
    get_n2()
        Returns refractive index of area at z > z_0
    
    get_curvature()
        Returns curvature of spherical surface
    
    intercept(ray)
        Returns numpy array describing position of intercept of ray and 
        screen.
    
    propagate_ray(ray)
        Appends new position and unchanged direction of ray to 3D array 
        containing ray's position and direction history. 
        Calculates intercept position using intercept().
    

        
    Notes
    -----
        Hidden Attributes
        -----------------
        _z_0 : float
            z position optical element
            
        _r : float
            aperture radius of optical element
    '''
    
    def __init__(self, z_0, r): # r is aperture radius
        OpticalElement.__init__(self, z_0, r)

    def __repr__(self):
        return f'OpticalElement(z_0 = {self._z_0}, r = {self._r})'
    def intercept(self, ray):
        l = (self._z_0 - ray._p[2])/ ray._k[2]
        if (ray._p + l * ray._k)[0] ** 2 + (ray._p + l * ray._k)[1] ** 2 > self._r **2 :# if resulting position is outside aperture radius then position must be excluded
            print(f'{str(ray)} does not intercept with plane, it is outside the aperture radius of the plane')  
            return None
        else:
            print(f'{str(ray)} intercepts plane')
            return ray._p + l * ray._k

    def propagate_ray(self, ray):
        intercept_posn = self.intercept(ray)
        if isinstance(intercept_posn, type(None)):
            ray._terminate()
        else:
            ray.append(intercept_posn, ray._k) # direction of ray remains unchanged
            print(f'ray intercepts plane, is now {str(ray)}')