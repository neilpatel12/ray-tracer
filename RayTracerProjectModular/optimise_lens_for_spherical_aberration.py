#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  5 10:36:03 2020

@author: neilpatel
"""
'''
optimise_lens_for_spherical_aberration.py covers optimising biconvex lens 

'''
#%%
import numpy as np
# from RayTracerProjectModular 
import ray_and_ray_bundle as rrb
# from RayTracerProjectModular 
import optical_elements as oe
# from RayTracerProjectModular
import ray_tracer_tools as tool
import matplotlib.pyplot as plt
import scipy.optimize as opt

lens_p = oe.SphericalRefraction(z_0 = 100, curv = 0.02, n1 = 1, n2 = 1.5168, r = 100)
lens_f = oe.SphericalRefraction(z_0 = lens_p.get_z_0() + 5, curv = 0, n1 = lens_p.get_n2(), n2 = 1, r = 100)
plane = oe.OutputPlane(z_0 = lens_p.get_z_0() + 100, r = 1000)
bundle_5mm = rrb.RayBundle(p_centre = [0,0,0], k = [0,0,1], bundle_radius = 5, xy_spacing = 1)

def lens_tester_constrict_fl_opt(lens_p_curvaturelens_f_curvature, focal_length = 100): #focal length is fixed in function
    bundle_5mm = rrb.RayBundle(p_centre = [0,0,0], k = [0,0,1], bundle_radius = 5, xy_spacing = 1)
    lens_p_curvature = lens_p_curvaturelens_f_curvature[0]
    lens_f_curvature = lens_p_curvaturelens_f_curvature[1]
    lens_p = oe.SphericalRefraction(z_0 = 100, curv = lens_p_curvature, n1 = 1, n2 = 1.5168, r = 100)
    lens_f = oe.SphericalRefraction(z_0 = lens_p.get_z_0() + 5, curv = lens_f_curvature, n1 = lens_p.get_n2(), n2 = 1, r = 100)
    focal_plane = oe.OutputPlane(z_0 = lens_f.get_z_0() + focal_length, r = 1000)   #treat focal length as 100
    
    opt_elements_foc = [lens_p, lens_f, focal_plane]
    
    i=0
    for element in opt_elements_foc:
        element.propagate_ray_bundle(bundle_5mm)
        i+=1

    return bundle_5mm.RMS_dev()

#%%
bundle_5mm = rrb.RayBundle(p_centre = [0,0,0], k = [0,0,1], bundle_radius = 5, xy_spacing = 1)

curv1_guess = 0.001
curv2_guess = 0.001
#curvs 0.01269652, -0.00633395
curv1_guess = 0.02 
curv2_guess = 0.02
#curvs -0.57361683, -0.65209902 # only one ray remaining, other diverge
curv1_guess = 0.01269652
curv2_guess = -0.00633395
# gives itself
curv1_guess = 0.014
curv2_guess = 0
#curvs 0.01624465, -0.00261178 
curv1_guess = 0.01624465
curv2_guess = -0.00261178
#gives itself, message: 'Desired error not necessarily achieved due to precision loss.'
curv1_guess = 0.01
curv2_guess = 0
#curvs 0.01450175, -0.00444576
curv1_guess = 0.02
curv2_guess = 0
#curvs -0.71851219,  0.68898457, 'Optimization terminated successfully.'
opt.minimize(lens_tester_constrict_fl_opt, (curv1_guess, curv2_guess), args = (100)) 

#%%



def lens_tester_constrict_fl(lens_p_curvaturelens_f_curvature, focal_length = 100): #focal length is fixed in function
    bundle_5mm = rrb.RayBundle(p_centre = [0,0,0], k = [0,0,1], bundle_radius = 5, xy_spacing = 1)
    lens_p_curvature = lens_p_curvaturelens_f_curvature[0]
    lens_f_curvature = lens_p_curvaturelens_f_curvature[1]
    lens_p = oe.SphericalRefraction(z_0 = 100, curv = lens_p_curvature, n1 = 1, n2 = 1.5168, r = 100)
    lens_f = oe.SphericalRefraction(z_0 = lens_p.get_z_0() + 5, curv = lens_f_curvature, n1 = lens_p.get_n2(), n2 = 1, r = 100)
    focal_plane = oe.OutputPlane(z_0 = lens_f.get_z_0() + focal_length, r = 1000)   #treat focal length as 100
    
    opt_elements_foc = [lens_p, lens_f, focal_plane]
    
    i=0
    for element in opt_elements_foc:
        element.propagate_ray_bundle(bundle_5mm)
        i+=1
    bundle_5mm.spot_diagram() #######FLAG THIS TO MARKER IN INTERVIEW
    plt.title(f'Spot diagram in focal plane\n of optimised lens with curvature {lens_p_curvaturelens_f_curvature}')
    bundle_5mm.plot_paths_3D()
    return bundle_5mm.RMS_dev()

#%%
bundle_5mm = rrb.RayBundle(p_centre = [0,0,0], k = [0,0,1], bundle_radius = 5, xy_spacing = 1)
lens_tester_constrict_fl((0.01269652, -0.00633395), focal_length = 100) #curvs 0.01269652, -0.00633395

#%%
bundle_5mm = rrb.RayBundle(p_centre = [0,0,0], k = [0,0,1], bundle_radius = 5, xy_spacing = 1)
lens_tester_constrict_fl((-0.57361683, -0.65209902), focal_length = 100) #curvs 0.01269652, -0.00633395

#%% Plot in report
bundle_5mm = rrb.RayBundle(p_centre = [0,0,0], k = [0,0,1], bundle_radius = 5, xy_spacing = 1)
lens_tester_constrict_fl((0.01624465, -0.00261178), focal_length = 100) #curvs 0.01269652, -0.00633395

#%% Plot in report
bundle_5mm = rrb.RayBundle(p_centre = [0,0,0], k = [0,0,1], bundle_radius = 5, xy_spacing = 1)
lens_tester_constrict_fl((0.01450175, -0.00444576), focal_length = 100)

#%%
bundle_5mm = rrb.RayBundle(p_centre = [0,0,0], k = [0,0,1], bundle_radius = 5, xy_spacing = 1)
lens_tester_constrict_fl((-0.71851219,  0.68898457), focal_length = 100)

#%%





bundle_5mm.plot_paths_3D()