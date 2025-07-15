Bader Paralellized

The function GetCPCL_Spatial is a paralellized version of GetCPCL and can be enabled/disabled by setting thread_count in options_mod.f90 to 1 or greater than 1. 

GetCPCL_Spatial splits the 3 dimensional grid into spatial regions for each thread, and each thread checks for candidates in its region, ensuring that there is only one candidate in a 3 point radius within each spatial region. These points are saved to a local thread array and merged at the end of the subroutine, with another check removing candidates in close proximity across the thread regions. Then, the cpcl array is resized to fit the exact amount of elements, and then the cpl array is resized to the same size.

Please email agar5333@gmail.com with any questions.
