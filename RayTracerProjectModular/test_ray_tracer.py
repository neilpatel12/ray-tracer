#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  4 12:27:13 2020

@author: neilpatel
"""
'''
test_ray_tracer goes over Tasks 9,10,11

produces Fig.2, Fig.3 from report

Do not need to restart kernal in each code block, but do restart kernel when new task is started.
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

#%% TASK 9
'''
Task 9: Start with a simple system comprising a single spherical refracting 
surface and an output plane. Consider a spherical surface at z=100 mm with 
curvature 0.03mm−1 and refractive indices n1=1.0 and n2=1.5 and an output plane
 at z=250mm .
Try tracing a few example rays from different points in the input plane (z=0mm)
 to the output plane.

When you have traced rays through the optical system, you can plot their paths 
using matplotlib. Remember; the method vertices should return the points along 
the ray.
'''
lens = oe.SphericalRefraction(z_0 = 100, curv = 0.03, n1 = 1, n2 = 1.5,\
                              r = 1/0.03)
plane = oe.OutputPlane(z_0 = 250, r = 100)

#%% Task 9
## 2D PLOT, paraxial

k = [0,0,1]

bundle_rad = 10 
step = 2

#fig1 = plt.figure()
fig1, (ax01,ax02) = plt.subplots(2, 1, constrained_layout= False)

fig1.suptitle('Task 9 \n\n')
# generating rays with initial positions distributed in the vertical axis
for a in range(-bundle_rad, bundle_rad + step, step):
    ray = rrb.Ray(p = [0, a, 0], k = k)
    lens.propagate_ray(ray)
    plane.propagate_ray(ray)    
    x = ray.vertices()[:, 0] # finding coordinates to plot rays
    y = ray.vertices()[:, 1]
    z = ray.vertices()[:, 2] 
    ax01.plot(z, y, '-')    # plot rays, viewing system facing in the +x dir

ax01.set_title(f'Rays with direction {k}, in plane x = 0')    
ax01.set_ylabel('Y(mm)')   
ax01.set_xlabel('Z(mm)') 
fig1.show()

## 2D PLOT, not paraxial

k = [0,0.1,1]

#fig2 = plt.figure()
# generating rays with initial positions distributed in the vertical axis
for a in range(-bundle_rad, bundle_rad + step, step):
    ray = rrb.Ray(p = [0, a, 0], k = k)
    lens.propagate_ray(ray)
    plane.propagate_ray(ray)    
    x = ray.vertices()[:, 0] # finding coordinates to plot rays
    y = ray.vertices()[:, 1]
    z = ray.vertices()[:, 2] 
    ax02.plot(z, y, '-')    # plot rays, viewing system facing in the +x dir

ax02.set_title(f'Rays with direction {k}, in plane x = 0')    
ax02.set_ylabel('Y(mm)')   
ax02.set_xlabel('Z(mm)') 
plt.show()




#%% Task 9
## 3D PLOT, arb angle
k = [0,0.1,1]

bundle_rad = 10 # radius of the bundle
step = 5        # step in radius of concentric circles

# creating figure and axis labels
fig3 = plt.figure(figsize = plt.figaspect(1)*1.5)
ax = fig3.add_subplot(111, projection='3d')
ax.set_xlabel('X')
ax.set_ylabel('Z')
ax.set_zlabel('Y')
ax.set_xlim(-25,25)
ax.set_zlim(-25,25)
ax.set_ylim(0,260)


for a in range(-bundle_rad, bundle_rad + step, step): # Generates rays evenly 
                                                      #spaced around concentric
                                                      #circles 
    for theta in np.linspace(0, np.pi, 15):
        ray = rrb.Ray(p = [a * np.cos(theta), a * np.sin(theta), 0], k = k)
        lens.propagate_ray(ray)
        plane.propagate_ray(ray)
        x = ray.vertices()[:, 0] # finding coordinates to plot rays
        y = ray.vertices()[:, 1]
        z = ray.vertices()[:, 2] 
        ax.plot3D(x,y,z, zdir = 'y', alpha = 0.5) # plot ray path in each iter
plt.title(f'Rays with direction {k}')         

#%% TASK 10
'''
Task 10: Trace a few rays parallel to the optical axis through your system. 
Estimate the position of the paraxial focus using a ray close to the optical 
axis (say 0.1 mm).
How does this compare to the expected paraxial focal point for a spherical 
surface?
'''
lens = oe.SphericalRefraction(z_0 = 100, curv = 0.03, n1 = 1, n2 = 1.5,\
                              r = 1/0.03)
plane = oe.OutputPlane(z_0 = 250, r = 100)

fig5 = plt.figure()
for a in [-0.1, 0.1]:
    #generates two paraxial rays with small perpendicular displacement
    # suitable for obtaining paraxial focus
    ray = rrb.Ray(p = [0, a, 0], k = [0,0,1])
    lens.propagate_ray(ray)
    plane.propagate_ray(ray)
    x = ray.vertices()[:, 0] # finding coordinates to plot rays
    y = ray.vertices()[:, 1]
    z = ray.vertices()[:, 2] 
    plt.plot(z, y, '-')       # 2D plot, facing in +x direction
    
plt.title(f'Rays with direction {k}, in plane x = 0') 
plt.ylabel('Y(mm)')   
plt.xlabel('Z(mm)') 
fig5.show()



focal = - ray.vertices()[1][1] * \
(ray.directions()[2][2]/ray.directions()[2][1]) 
# finds focal length by calculating where the ray intercepts optical axis
# simple y = mx +c

print('')

print('focal length calculated by SphericalRefraction method is'\
      , lens.parax_focus())

print('paraxial focal length is ', focal, 'mm') # 99.99964999931262 mm

f_theory = (lens._R * lens._n2)/(lens._n2 - lens._n1)
print('theoretical focal lenth is', f_theory) #100 mm

perc_diff = 100 * (focal - f_theory)/f_theory
print('percentage difference in simulation determined paraxial focal length \
      and theoretical approximate focal length is ', perc_diff)


# percentage difference = -0.0003500006873764505%
# could measure how perc diff in f varies with R.

#%% TASK 11
'''
Task 11:  Can you think of any other simple test cases, given that the 
spherical surface is supposed to be capable of forming an image?
Does it work as you expect?
'''

# Capable of imaging:
# - Try arrow, see if arrow is inverted on other side of lens, when ray is
#   propagated paraxially
# - Try 

#generate rays with inital positions that form upward facing arrow

lens = oe.SphericalRefraction(z_0 = 100, curv = 0.03, n1 = 1, n2 = 1.5, \
                              r = 1/0.03)
plane = oe.OutputPlane(z_0 = 250, r = 100)

# instantiating rays with positions chosen to form upward facing arrow
ray0 = rrb.Ray(p = [0, 0, 0], k = [0, 0, 1])
ray1 = rrb.Ray(p = [0, -1, 0], k = [0, 0, 1])
ray2 = rrb.Ray(p = [0, -2, 0], k = [0, 0, 1])
ray3 = rrb.Ray(p = [0, -3, 0], k = [0, 0, 1])
ray4 = rrb.Ray(p = [0, 1, 0], k = [0, 0, 1])
ray5 = rrb.Ray(p = [0, 2, 0], k = [0, 0, 1])
ray6 = rrb.Ray(p = [0, 3, 0], k = [0, 0, 1])   
ray7 = rrb.Ray(p = [0, 4, 0], k = [0, 0, 1])   
ray8 = rrb.Ray(p = [-2, 2, 0], k = [0, 0, 1])
ray9 = rrb.Ray(p = [2, 2, 0], k = [0, 0, 1])
ray10 = rrb.Ray(p = [-1, 3, 0], k = [0, 0, 1])
ray11 = rrb.Ray(p = [1, 3, 0], k = [0, 0, 1])
ray12 = rrb.Ray(p = [-1, 2, 0], k = [0, 0, 1])
ray13 = rrb.Ray(p = [1, 2, 0], k = [0, 0, 1])


ray_list = [ray0, ray1, ray2, ray3, ray4, ray5, ray6, ray7, ray8, ray9, \
            ray10, ray11, ray12, ray13]

# creating seperate figures

fig6 = plt.figure(figsize = plt.figaspect(1)*1.5)
ax1 = fig6.add_subplot(111, projection='3d')
ax1.set_xlabel('X')
ax1.set_ylabel('Z')
ax1.set_zlabel('Y')



fig7, (ax2, ax3) = plt.subplots(1, 2, constrained_layout= False)

fig7.suptitle('Task 11')
ax2.set_title('Spot diagram of arrow before propagation')
ax2.set_xlabel('X (mm)')
ax2.set_ylabel('Y (mm)')


ax3.set_title('Spot diagram of arrow image after propagation')
ax3.set_xlabel('X (mm)')
ax3.set_ylabel('Y (mm)')

ax1.set_title('Arrow imaging')
for ray in ray_list:
        x_before = ray.p()[0]         # Coordinates of ray's position before
        y_before = ray.p()[1]         # propagation through system.
        lens.propagate_ray(ray)
        plane.propagate_ray(ray)
        x = ray.vertices()[:, 0] # Coordinates of ray's position at each stage 
        y = ray.vertices()[:, 1] # of propagation through the system.
        z = ray.vertices()[:, 2] #
        x_after = ray.p()[0]          # Coordinates of ray's position at end of
        y_after = ray.p()[1]          # propagating through system.
        
        ax1.plot3D(x,y,z, zdir = 'y', alpha = 0.5) # Plot ray path in each iter
        ax2.scatter(x_before, y_before, color = 'g') # Spot diagram before
        ax3.scatter(x_after, y_after, color = 'b')   # Spot diagram after


# OUTPUT FROM THIS CODE BLOCK
## 3D plot of arrow imaging
## spot diagram before and after spherical surface imaging

#%%

