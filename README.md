# Main Features
**runPFEM** is an open-source application for the simulation of nonlinear solid mechanics, free surface flow and fluid-structure interaction problems using a positional-based Particle Finite Element Method (PFEM).

It was designed and implemented in C++ aiming high performance, thus it uses the PETSc interface of MPI protocol for running in parallel under Unix systems.

# Dependencies
In order to build **runPFEM**, the user must install previously the following libraries:

## CMake
CMake can be installed directly from unix repositories with the command:
```bash
sudo apt-get install cmake
```
## MPICH
PETSc uses MPICH to deal with parallelism, so you can decide either install the MPICH and give its path to PETSc or let PETSc download it during the configuration process. The first option is recommended because the path to the *mpiexec* becomes fixed and independent of the PETSc build configuration.
You can download the file [mpich-3.4.2.tar.gz](https://github.com/pmodels/mpich/releases/tag/v3.4.2) and follow the instructions on [github](https://github.com/pmodels/mpich):
- Unpack the tar file and go to the top level directory:
```bash
tar xzf mpich-3.4.2.tar.gz
cd mpich-3.4.2
```
- Configure MPICH specifying the installation directory (could be either an empty directory or an non existent directory) and device:
```bash
./configure --prefix=/path/to/mpi/installation/directory --with-device=ch4:ofi 2>&1 | tee c.txt
```
- Build MPICH:
```bash
make 2>&1 | tee m.txt
```
- Install the MPICH commands:
```bash
make install 2>&1 | tee mi.txt
```
- Add the bin subdirectory of the installation directory to your PATH by adding the following line to the file ~/.bashrc:
```bash
export PATH="/path/to/mpi/installation/directory/bin:$PATH"
```
Check that everything is in order at this point by doing:
```bash
which mpicc
which mpicxx
which mpiexec
```
These commands should print the path to the bin subdirectory of the MPICH installation directory.

## METIS
Metis can be downloaded and installed during the configuration process of PETSc.

## PETSc
The best way to download PETSc is cloning it from gitlab. This way, you can get any update or new releases by just using the command *git pull*.
- Clone PETSc repository and go to the top level directory:
```bash
git clone -b release https://gitlab.com/petsc/petsc.git /path/to/petsc/directory
cd /path/to/petsc/directory
```
- Configure
If you want to build a Release (fast) and a Debug (slow) version of **runPFEM**, it is recommended to have two different configurations of PETSc, namely "architectures", builded in your system. By doing that, you can simply switch from one version to another by properly setting the environmental variable PETSC_ARCH. The Release build can be configured by doing:
```bash
./configure PETSC_ARCH=arch-linux2-c-opt --with-mpi-dir=/path/to/mpi/installation/directory --with-cxx-dialect=C++11 --with-debugging=0 --with-X=1 COPTFLAGS='-O3 -march=native -mtune=native' CXXOPTFLAGS='-O3 -march=native -mtune=native' FOPTFLAGS='-O3 -march=native -mtune=native' --download-metis --download-parmetis --download-mumps --download-scalapack --download-ptscotch --download-fblaslapack --download-hdf5
```
Here, */path/to/mpi/installation/directory* corresponds to the directory where MPICH was installed.
- Build
After the configuration, it will appear on the terminal the *make* command to build the respectively PETSc architecture. Just copy it and paste it in terminal.
- Test
Once the build process is finalized, the *make* command to execute the tests will appear on the terminal.
- Add the PETSc directory to your PATH by adding the following line to the file ~/.bashrc:
```bash
export PETSC_DIR="/path/to/petsc/directory"
```
The same process can be applied to install a Debug version of PETSc, you only need to change the arguments of the configure command to this one:
```bash
./configure PETSC_ARCH=arch-linux2-c-debug --with-mpi-dir=/path/to/mpi/installation/directory --with-cxx-dialect=C++11 --with-debugging=1 --with-X=1 --download-metis --download-parmetis --download-mumps --download-scalapack --download-ptscotch --download-fblaslapack --download-hdf5
```

## Gmsh
Gmsh 3.0.6 is available in Linux repository and has all the functionality needed in **runPFEM**. You can get it by the command:
```bash
sudo apt-get install gmsh
```
If you decide to download its latest version from the Gmsh [website](https://gmsh.info/#Download), make sure that the binary directory is included in the global PATH. Otherwise, you must specify its absolute path in the argument *gmshPath* when calling the function *generateMesh()*.

## Lapacke
Lapacke library can be installed directly from the linux repository through the command:
```bash
sudo apt-get install liblapacke-dev
```

# Configuration and Build
The process of configuring the build is performed only once for each build type. The user can configure and build the code in three different ways:

## Using the extension CMake Tools inside VSCode (Preferred)

Configuring and building the code inside VSCode can be easily done using the extension CMake Tools.
- Firstly, go to the *Extension* tab on the left panel, search and install CMake Tools. 
- Open the file *settings.json* located inside the .vscode directory and change the variable *PETSC_ARCH* to the name used during PETSc installation for both build types. Note: if you used the same command to configure PETSc as the ones suggested here, please go to the next step.
- The CMake commands can be accessed opening the Command Palette (Ctrl+Shift+P) and typing CMake:command, or through the buttons at the blue bar at the lower left corner. Use the command *CMake:Scan for Kits*, and then select the available compilers to be used.
- Select the build type through the command *CMake:Select Variant*. The available options are Debug and Release. The active variant will also appear in the blue lower bar as *CMake: \[variant\]: Ready*. By clicking on this icon you can change the variant as well.
- To configure the build type, just type *CMake:Configure*. After this step, CMake will generate a makefile to actually build the code, located in the build/BuildType directory.
- To compile, simply type *Cmake:Build* or click on the Build button of the blue lower bar.
Bear in mind that you always can navigate from one build type to another by simply changing the active variant with *CMake:Select Variant*.

## Using the terminal interface CCMake

- Install the CMake GUI with the command:
```bash
sudo apt-get install cmake-curses-gui
```
- Create the build directory and go to the top directory:
```bash
mkdir build/BUILD_TYPE
cd build/BUILD_TYPE
```
- Set the PETSC_ARCH environmental variable according to the building type:
```bash
export PETSC_ARCH="arch-type"
```
- Run the CCMake on the top directory informing the root directory:
```bash
ccmake ../../.
```
and press C to configure.
- After the process is done, the environmental variables will be displayed on the terminal. Move the cursor to CMAKE_BUILD_TYPE, press ENTER to edit it according to the desired variant and press ENTER again. Check if PETSC_ARCH is set properly. Press C to configure again and G to terminate CCMake.
- To compile, simply type:
```bash
make
```
or to compile in parallel, simply add the option -j followed by the number of processors. For example, to compile using 4 threads:
```bash
make -j4
```

## Manually using CMake: user must specify the building type (Debug or Release)

- Create the build directory and go to the top directory:
```bash
mkdir build/BUILD_TYPE
cd build/BUILD_TYPE
```
- Set the PETSC_ARCH environmental variable according to the building type:
```bash
export PETSC_ARCH="arch-type"
```
- Use the CMake command specifying the root directory (-H), the build directory (-B) and the build type (-DCMAKE_BUILD_TYPE) as follows:
```bash
cmake -H"../../." -B"." -DCMAKE_BUILD_TYPE=BUILD_TYPE
```
- To compile, simply type:
```bash
make
```
or to compile in parallel, simply add the option -j followed by the number of processors. For example, to compile using 4 threads:
```bash
make -j4
```

# Running the application

## From Terminal
Once compiled, a executable binary should be created in the build directory. It can be executed from the terminal as:
```bash
mpiexec -n n_process path_to_executable
```
where n_process is the number of processes to be used in the application and path_to_executable is the absolute path to the binary, located on path-to-runPFEM/build/BUILD_TYPE/runPFEM.

## From VSCode
There is also the option to run the application inside the VSCode. This is particularly useful when debugging the code because VSCode has a built-in interface which makes it easier to add breakpoints and to check memory issues. To do so, just make sure that PETSC_ARCH is correctly set in *settings.json*.

- Debugger: To start the debugger, simply press F5 or go to the Debug icon on the left bar and click on the green play button. When you add a breakpoint, the execution should stop at that point and you can have access to all the needed information.

- Release: In order to run in Release mode and also in parallel, click on Terminal/Run Task and select the task named Run. Choose Release from the popup window and next inform the number of processes to be used.

