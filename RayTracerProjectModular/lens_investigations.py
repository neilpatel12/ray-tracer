#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  4 22:51:29 2020

@author: neilpatel
"""

'''
lens_investigations runs Task 15
'''

#%% THIS IS PLOTTED, RUN THIS FOR PLOT
#convex then flat (|
'''
focal length = 93.45281250408226 mm
lens_maker focal length = 96.74922600619198
RMS spot radius at surface 3: 0.06521823592293968 mm
'''
import numpy as np
# from RayTracerProjectModular 
import ray_and_ray_bundle as rrb
# from RayTracerProjectModular 
import optical_elements as oe
# from RayTracerProjectModular 
import ray_tracer_tools as tool
import matplotlib.pyplot as plt


lens_p = oe.SphericalRefraction(z_0 = 100, curv = 0.02, n1 = 1, n2 = 1.5168, r = 10)
lens_f = oe.SphericalRefraction(z_0 = lens_p.get_z_0() + 5, curv = 0, n1 = lens_p.get_n2(), n2 = 1, r = 100)
plane = oe.OutputPlane(z_0 = lens_p.get_z_0() + 100, r = 1000)

bundle = rrb.RayBundle(p_centre = [0,0,0], k = [0,0,1], bundle_radius = 5, xy_spacing = 1)

thin_bundle = rrb.RayBundle(p_centre = [0,0,0], k = [0,0,1], bundle_radius = 0.001, xy_spacing = 0.001)

opt_elements = [lens_p, lens_f, plane ]

i = 0
for element in opt_elements:
    i+=1
    element.propagate_ray_bundle(bundle)
    #bundle.spot_diagram()                         ####### FLAG TO MARKER THAT THIS IS NOT HASHED OUT
    #plt.title(f'ray on surface {i}')              #######
    print(f'RMS spot radius on surface {i}:', bundle.RMS_dev(), 'mm')

# estimating paraxial focus

test_ray = rrb.Ray(p = [0,0.000001,0], k = [0,0,1])

for element in opt_elements:
    element.propagate_ray(test_ray)

focal_length = - test_ray.vertices()[-2][1] * (test_ray.directions()[-2][2]/test_ray.directions()[-2][1])
print(f'focal length = {focal_length} mm')


# lensmakers equation

focal_length_lme = 1/((lens_p.get_n2() - 1) * (lens_p.get_curvature() - lens_f.get_curvature() + (lens_p.get_n2() - 1) * ((lens_f.get_z_0() - lens_p.get_z_0()) * lens_p.get_curvature() * lens_f.get_curvature()) /(lens_p.get_n2())))
print(f'lens_maker focal length = {focal_length_lme}')
focal_length_hand = 1/(0.5168 * (0.02))

print(f'RMS spot radius at surface {i}:', bundle.RMS_dev(), 'mm')

bundle.spot_diagram()
plt.title('Task 15 \n plano-convex') ####FLAG TO MARKER
bundle.plot_paths_3D()
plt.title('Task 15\n plano-convex, convex first')
#%%
#propagate through focal plane, ^^^above was not through focal plane
bundle = rrb.RayBundle(p_centre = [0,0,0], k = [0,0,1], bundle_radius = 5, xy_spacing = 1)
focal_plane = oe.OutputPlane(z_0 = lens_f.get_z_0() + focal_length, r = 1000)
opt_elements_foc = [lens_p, lens_f, focal_plane]

i=0
for element in opt_elements_foc:
    i+=1
    element.propagate_ray_bundle(bundle)
    bundle.spot_diagram()
    plt.title(f'ray on {element} surface {i}')
    print(f'RMS spot radius at {element}:', bundle.RMS_dev(), 'mm')
bundle.plot_paths_3D()
#RMS 0.00781537033105665 mm
#%% THIS IS PLOTTED, RUN THIS FOR PLOT
# flat then convex |)
'''
focal length = 96.74922600619195 mm
lens_maker focal length = 96.74922600619198
RMS spot radius at surface 3: 0.03783731067228958 mm
'''
import numpy as np
from RayTracerProjectModular import ray_and_ray_bundle as rrb
from RayTracerProjectModular import optical_elements as oe
from RayTracerProjectModular import ray_tracer_tools as tool
import matplotlib.pyplot as plt

lens_p = oe.SphericalRefraction(z_0 = 100, curv = 0, n1 = 1, n2 = 1.5168, r = 100)
lens_f = oe.SphericalRefraction(z_0 = lens_p.get_z_0() + 5, curv = -0.02, n1 = lens_p.get_n2(), n2 = 1, r = 100)
plane = oe.OutputPlane(z_0 = lens_p.get_z_0() + 100, r = 1000)

bundle = rrb.RayBundle(p_centre = [0,0,0], k = [0,0,1], bundle_radius = 5, xy_spacing = 1)

thin_bundle = rrb.RayBundle(p_centre = [0,0,0], k = [0,0,1], bundle_radius = 0.001, xy_spacing = 0.001)

opt_elements = [lens_p, lens_f, plane ]


i = 0
for element in opt_elements:
    i+=1
    element.propagate_ray_bundle(bundle)
    bundle.spot_diagram()                  ####### FLAG TO MARKER THAT THIS IS NOT HASHED OUT 
    plt.title(f'ray on surface {i}')       #######
    print(f'RMS spot radius on surface {i}:', bundle.RMS_dev(), 'mm')

# estimating paraxial focus

test_ray = rrb.Ray(p = [0,0.000001,0], k = [0,0,1])

for element in opt_elements:
    element.propagate_ray(test_ray)

focal_length = - test_ray.vertices()[-2][1] * (test_ray.directions()[-2][2]/test_ray.directions()[-2][1])
print(f'focal length = {focal_length} mm')


# lensmakers equation

focal_length_lme = 1/((lens_p.get_n2() - 1) * (lens_p.get_curvature() - lens_f.get_curvature() + (lens_p.get_n2() - 1) * ((lens_f.get_z_0() - lens_p.get_z_0()) * lens_p.get_curvature() * lens_f.get_curvature()) /(lens_p.get_n2())))
print(f'lens_maker focal length = {focal_length_lme}')
focal_length_hand = 1/(0.5168 * (0.02))

print(f'RMS spot radius at surface {i}:', bundle.RMS_dev(), 'mm')

bundle.spot_diagram()
plt.title('Task 15 \n plano-convex')  ####FLAG TO MARKER
bundle.plot_paths_3D()
plt.title('Task 15\n plano-convex, plane first')
#%%
#propagate through focal plane, ^^^above was not through focal plane
bundle = rrb.RayBundle(p_centre = [0,0,0], k = [0,0,1], bundle_radius = 5, xy_spacing = 1)
focal_plane = oe.OutputPlane(z_0 = lens_f.get_z_0() + focal_length, r = 1000)
opt_elements_foc = [lens_p, lens_f, focal_plane]

i=0
for element in opt_elements_foc:
    i+=1
    element.propagate_ray_bundle(bundle)
    bundle.spot_diagram()
    plt.title(f'ray on {element} surface {i}')
    print(f'RMS spot radius at {element}:', bundle.RMS_dev(), 'mm')
bundle.plot_paths_3D()

#RMS 0.031129492486365772 mm
#%%
# flat then concave |(
'''
focal length = -96.74922600619195 mm
lens_maker focal length = -96.74922600619198
RMS spot radius at surface 3: 7.158232195770108 mm
'''
import numpy as np
from RayTracerProjectModular import ray_and_ray_bundle as rrb
from RayTracerProjectModular import optical_elements as oe
from RayTracerProjectModular import ray_tracer_tools as tool
import matplotlib.pyplot as plt


lens_p = oe.SphericalRefraction(z_0 = 100, curv = 0, n1 = 1, n2 = 1.5168, r = 100)
lens_f = oe.SphericalRefraction(z_0 = lens_p.get_z_0() + 5, curv = 0.02, n1 = lens_p.get_n2(), n2 = 1, r = 100)
plane = oe.OutputPlane(z_0 = lens_p.get_z_0() + 100, r = 1000)

bundle = rrb.RayBundle(p_centre = [0,0,0], k = [0,0,1], bundle_radius = 5, xy_spacing = 1)

opt_elements = [lens_p, lens_f, plane ]


i = 0
for element in opt_elements:
    i+=1
    element.propagate_ray_bundle(bundle)
    bundle.spot_diagram()
    plt.title(f'ray on {element} surface {i}')
    print(f'RMS spot radius at {element}:', bundle.RMS_dev(), 'mm')
bundle.plot_paths_3D()

    
# estimating paraxial focus

test_ray = rrb.Ray(p = [0,0.000001,0], k = [0,0,1])

for element in opt_elements:
    element.propagate_ray(test_ray)

focal_length = - test_ray.vertices()[-2][1] * (test_ray.directions()[-2][2]/test_ray.directions()[-2][1])
print(f'focal length = {focal_length} mm')


# lensmakers equation

focal_length_lme = 1/((lens_p.get_n2() - 1) * (lens_p.get_curvature() - lens_f.get_curvature() + (lens_p.get_n2() - 1) * ((lens_f.get_z_0() - lens_p.get_z_0()) * lens_p.get_curvature() * lens_f.get_curvature()) /(lens_p.get_n2())))
print(f'lens_maker focal length = {focal_length_lme}')

focal_length_hand = 1/(0.5168 * (0.02))

print(f'RMS spot radius at surface {i}:', bundle.RMS_dev(), 'mm')
#%%
'''
focal length = -100.04563950830166 mm
lens_maker focal length = -96.74922600619198
RMS spot radius at surface 3: 7.281450032291046 mm
'''

# concave then flat )|
import numpy as np
from RayTracerProjectModular import ray_and_ray_bundle as rrb
from RayTracerProjectModular import optical_elements as oe
from RayTracerProjectModular import ray_tracer_tools as tool
import matplotlib.pyplot as plt

lens_p = oe.SphericalRefraction(z_0 = 100, curv = -0.02, n1 = 1, n2 = 1.5168, r = 100)
lens_f = oe.SphericalRefraction(z_0 = lens_p.get_z_0() + 5, curv = 0, n1 = lens_p.get_n2(), n2 = 1, r = 100)
plane = oe.OutputPlane(z_0 = lens_p.get_z_0() + 100, r = 1000)

bundle = rrb.RayBundle(p_centre = [0,0,0], k = [0,0,1], bundle_radius = 5, xy_spacing = 1)

opt_elements = [lens_p, lens_f, plane ]


i = 0
for element in opt_elements:
    i+=1
    element.propagate_ray_bundle(bundle)
    bundle.spot_diagram()
    plt.title(f'ray on {element} surface {i}')
    print(f'RMS spot radius at {element}:', bundle.RMS_dev(), 'mm')
bundle.plot_paths_3D()
    
# estimating paraxial focus

test_ray = rrb.Ray(p = [0,0.000001,0], k = [0,0,1])

for element in opt_elements:
    element.propagate_ray(test_ray)

focal_length = - test_ray.vertices()[-2][1] * (test_ray.directions()[-2][2]/test_ray.directions()[-2][1])
print(f'focal length = {focal_length} mm')


# lensmakers equation

focal_length_lme = 1/((lens_p.get_n2() - 1) * (lens_p.get_curvature() - lens_f.get_curvature() + (lens_p.get_n2() - 1) * ((lens_f.get_z_0() - lens_p.get_z_0()) * lens_p.get_curvature() * lens_f.get_curvature()) /(lens_p.get_n2())))
print(f'lens_maker focal length = {focal_length_lme}')

focal_length_hand = 1/(0.5168 * (0.02))


print(f'RMS spot radius at surface {i}:', bundle.RMS_dev(), 'mm')


