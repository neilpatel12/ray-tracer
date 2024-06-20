#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  4 11:48:41 2020

@author: neilpatel
"""
'''
ray_tracer_tools.py contains functions necessary to run optical_elements.py and 
ray_and_ray_bundle.py
'''
import numpy as np

def normalise_vector(vector_list):
    if (vector_list[0], vector_list[1], vector_list[2]) == (0, 0, 0): #to prevent division by zero error
        return np.array([0,0,0], dtype = float)
    else:
        return np.array(vector_list, dtype = float)/ np.sqrt(vector_list[0] ** 2 + vector_list[1] ** 2 + vector_list[2] ** 2)
    