#!/usr/bin/env python3

import math

for i in range(2000):
    x = (i+1)*0.01
    print(' %17.9e %17.9e %17.9e %17.9e %17.9e %17.9e'
          % (x, 1/0.42*math.log(x)+5.2,
             x, 2/0.42*math.log(x)+5.2,
             2*x, 2/0.42))
