#!/usr/bin/env python3

'''This is a test to check that the matrix multiplication was properly performed'''

import sys
import numpy as np

def main(args):
    if (len(args)<3):
        print("Usage: executable A.mat B.mat C.mat")
        sys.exit(1)

    A=np.genfromtxt(args[1])
    B=np.genfromtxt(args[2])
    C=np.genfromtxt(args[3])

    C2=A.dot(B)
    L1=np.amax(np.abs(C-C2))

    if (L1>1E-5):        
        print("Arrays don't match by: ", L1)
        print("A: \n", A)
        print("B: \n", B)
        print("C: \n", C)
        print("C2: \n", C2)
        return 1
    
    print("Arrays match")
    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv))
