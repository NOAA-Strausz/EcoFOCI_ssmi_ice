#!/usr/bin/env python3
# -*- coding: utf-8 -*-


lat=65.75
lon=-169

import pygmt 

fig = pygmt.figure
region = [-169.5, -168.5, 65.5, 66]

fig.basemap(region=region, projection="M6i", Frame=True)

fig.plot(x=[-169.1,168.9], y=[65.65, 65.85], style='r', pen='2p, red')

fig.show()


