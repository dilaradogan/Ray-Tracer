# Ray-Tracer
This is a simple ray tracer project for Computer Graphics class. All methods and source code used for the Ray-Tracer are included in RayTracer.cpp file. While the program takes an input file which contains parameters for ray tracer, it gives a ppm file as output.

- Includes sphere objects and lights.
- Supports shadows, ambient, diffuse and specular shading. 
- Refraction and Phong reflection have also implemented. 
- For execution, the first command line argument should be file name of input.

  *RayTracer.exe test1.in*

  **Notes**: There is an extra 2 fields for the refraction in the sample input file. If you are going to use such a file as input, the second command line argument should be expressed as letter of y.
(y -> refraction yes)

  *RayTracer.exe test2refractive.in y*
- Sample input and output files are also avaliable in the project structure.
