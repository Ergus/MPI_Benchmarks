#!/usr/bin/env python3

'''This is a test to check that the matrix multiplication was properly performed'''

#!/usr/bin/env python3

import sys
import numpy as np

def main(args):
    if (len(args)<5):
        print("Usage: executable a_value initial_x initial_y final_y")
        sys.exit(1)

    a=np.float(args[1])
    x=np.genfromtxt(args[2])
    y1=np.genfromtxt(args[3])
    yc=np.genfromtxt(args[4])

    ypy=y1+a*x
    L1=np.amax(np.abs(ypy-yc))

    if (L1>1E-5):
        print("Arrays don't match")
        fmt="%3d\t y_py: %+.6f\t y_c: %+.6f = %d * %.6f %+.6f \t Diff: %.6f"
        arrsize=len(ypy)
        for var in zip(np.arange(arrsize),ypy,yc,np.repeat(a,repeats=arrsize),x,y1,np.abs(yc-ypy)):
            print(fmt % var)
        return 1
    
    print("Arrays match")
    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv))
