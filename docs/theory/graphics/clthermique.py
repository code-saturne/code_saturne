#!/usr/bin/env python3

import math

for i in range(300):
    x = (i+1)*0.1
    print(' %17.9e %17.9e' % (x, x*x*x/1000))

# To use the curve,
# two points must be placed for (a+at)/nu = 1/Pr = 1 
# of which one at x = 0,
# and two points must be placed for (a+at)/nu = 0.42 x / Prt = 0.42 x  
# of which one at x = 30.
