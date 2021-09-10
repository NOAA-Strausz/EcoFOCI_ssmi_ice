#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 28 16:30:45 2019

@author: strausz
"""

from haversine import haversine

point1=(0, -140)
point2=(1, -140)

test=haversine(point1, point2, unit="nmi" )

print(test)