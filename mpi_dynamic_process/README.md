README
======

This is a testing example that shows how to implement an MPI application that uses dynamically a number of processes.

The project contains 2 versions. 

1. dynamic.x is a C like implementation using function oriented programming. The only need for C++-11 in the compilation is to use threads. But this can be implemented using pthreads. This is a one file programm.

2. dynamicOOP.x is an OOP version for the previous program. This one uses clases and typical OOP tools. It needs many other files.

3. The reduction (process elimination) is implemented now from the tail always. It performs an split and then the extra processes goes to finalize.

Author: Jimmy Aguilar Mena
