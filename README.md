# Mesh-Construction-from-Point-Cloud-in-Parallel

The following instructions explain how to reproduce the results shown in our
final report. A link to our github repo is:
https://github.com/aditn/Mesh-Construction-from-Point-Cloud-in-Parallel

1. Unzip the inputs/ folder. These consist of the point clouds use for surface
   reconstruction:
   `unzip inputs.zip`

2. Download Eigen libary 3.3.1 from http://eigen.tuxfamily.org/index.php?title=Main_Page
   Extract the library to the `util/` folder. Enter the newly extracted folder and copy
   the `Eigen/` folder to `util/`. The extracted folder may now be deleted.

3. If `compile.bash` is not an executable run:
   `chmod a+x compile.bash`

4. Compile code. This produces the executable `runnable`.
   `./compile.bash`

5. Run our exectuable from its directory on an input in the inputs directory. It is best to use a rho value
   specified in our writeup, because too small a rho value will crash our program (we are unsure if this is a malloc issue or an algorithm issue). Fixing this bug is another future direction.
   `./runnable -f inputs/[filename] -t [number of processors] -p [rho-value specified in report]`

6. A detailed timing results will be printed to the console upon completion.

7. Outputs are located in `outputs/` in a `.obj` file.

8. If you would like to view debug images shown in report add a `-d` flag to
   Step 4. The debug files are placed in the home directory.


We generated our own cube and sphere point clouds through the `makeMesh.py`
script located in `inputs/`. To generate a cube or sphere with a user-defined
number of points perform the follwing steps.

1. Start the python interpreter from the commandline.
   `python`

2. Import the script file to use.
   `from makeMesh.py import *`

3. For cube run the follwing:
   `makeCubeMesh([(x,y,z) Center], [half side length],"[filename]",[numPoints=10000])`

4. For sphere run the following:
   `makeSphereMesh([(x,y,z) Center], [radius],"[filename]",[numPoints=10000])`

5. The point cloud will be in the inputs directory.
