Parallelization thesis
======================

Parallelization is especially useful and worthwhile in large-scale numerical computations in which a vast amount of independent, relatively simple, equations need to be solved, for example in estimating required pension reserves. In this thesis we will consider such problems drawn from pension and life insurance mathematics, and in particular how to estimate required pension reserves and cashflows, within the framework of the Actulus research project. Such work in particular can make use of the power of GPGPU's and the CUDA C programming languages. 

CUDA C is, however, based on the somewhat archaic C language. For safety and flexibility industry wants to automatically generate GPGPU code and to control its execution via the .NET platform. 

This project aims to allow for such calculations to be written in the F\# language on the .NET platform using Alea.cuBase which in turn is translated into CUDA C (or PTX, the CUDA C bytecode). Furthermore, I hope to implement a parser of the Actulus CalcSpec language files to be translated to F\# kernel-representations. 
This will allow for business to stick to a .NET workflow while using a modern language (F\#) and avoid having to work directly with (CUDA) C. 

The generated F\# CUDA code will be evaluated to to both an existing single-threaded C\# application provided by Peter Sestoft, manually written F\# cuBase code and manually written CUDA code. 

Through this project I hope to learn how to write high-performance GPGPU software for numerical problems, to demonstrate that I can design, implement, describe and critically evaluate the performance and other quality aspects of such software. 
\keywords{CUDA, parallelization, F\# Alea.cuBase, pension reserve estimation, Actulus CalcSpec}