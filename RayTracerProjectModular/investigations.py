#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  4 15:52:55 2020

@author: neilpatel
"""
'''
investigations.py runs Tasks 12 and 13

produces Fig.4, Fig.5, Fig.6 from report
'''
#%%
import numpy as np
from RayTracerProjectModular import ray_and_ray_bundle as rrb
from RayTracerProjectModular import optical_elements as oe
from RayTracerProjectModular import ray_tracer_tools as tool
import matplotlib.pyplot as plt

#%% TASK 12 and 13
'''
Task 12: Trace a bundle of rays for a uniform collimated beam for a larger 
diameter, e.g., 5mm to the paraxial focal plane.
'''
lens = oe.SphericalRefraction(z_0 = 100, curv = 0.03, n1 = 1, n2 = 1.5,\
                              r = 1/0.03)

focal_plane = oe.OutputPlane(z_0 = 200, r = 100)

bundle_5mm = rrb.RayBundle(p_centre = [0, 0, 0], k = [0, 0, 1], \
                           bundle_radius = 5, xy_spacing = 0.5)
# RayBundle generates rays with same direction and with initial positions
# centred on a grid on the xy plane, so rays are evenly spaced, giving a 
# uniform density of rays in the bundle.



bundle_5mm.spot_diagram('b')
plt.title(f'Instantiated Ray Bundle')

lens.propagate_ray_bundle(bundle_5mm)
bundle_5mm.spot_diagram('b')
plt.title(f'lens surface, {bundle_5mm.num_alive_rays()} rays remaining')

focal_plane.propagate_ray_bundle(bundle_5mm)
bundle_5mm.spot_diagram('g')

plt.title('Task 13\n Spot Diagram at focal plane')
print('RMS deviation from paraxial axis is', bundle_5mm.RMS_dev(), 'mm')

bundle_5mm.plot_paths_3D()

bundle_5mm.plot_2D_projection()
plt.title('Task 12 \n 10mm diameter ray bundle propagates refracting surface')