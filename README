# Argo Version 1.3. #

## Overview ##

The following is an implementation of a three-axis correction routine for images displaying distortion in Scanning Probe Microscopy.
 The distortion correction uses a grid-search routine for thermal drift correction followed by a Nelder-Mead Simplex routine for
  simultaneous correction of thermal drift and z-axis piezo creep.

## Compiling ##

Compile using Visual Studio 2013.
Requires the latest version of Eigen, Intel Thread Building Blocks (TBB), FreeImage and an OpenCL SDK.

The latest version of Eigen can be retreived from : http://eigen.tuxfamily.org/.
The latest version of TBB can be retreived from : https://www.threadingbuildingblocks.org/.
The latest version of FreeImage can be retreived from : http://freeimage.sourceforge.net/.
For OpenCL, you need to use the appropriate SDK associated with your hardware. This can be NVIDIA, Intel, or AMD.


To simplify, retreive the libraries and unzip them to folders named TBB, FreeImage and Eigen respectively. Place these
folders in the parent folder relative to where your solution is placed (e.g. in C:\ if you have the solution in C:\Project).
This simplifies the next few steps significantly such that you only need to configure the solution for your appropriate OpenCL
SDK. The steps have been outlined to include information for how to configure with respect to the libraries as well
if they were stored in a location different or named differently than the one suggested here.

The following needs to be configured properly in order for the Solution to compile successfully :

1)	Add in the full include path (usually in a folder called include) for the OpenCL SDK (TBB, Eigen, FreeImage as
	well if the folders are located differently than suggested above) to :

	Project | Properties | Configuration Properties | C/C++ | General | Additional Include Directories


2)	Add the full library path for the 32-bit/64-bit library files (usually in a folder called lib) for
	the OpenCL SDK (TBB and FreeImage as well if the folders are located differently than suggested above) to :

	Project | Properties | Configuration Properties | Linker | General | Additional Library Directories


3)	Only follow if the TBB or FreeImage library is stored differently than suggested above.
	The configuration has a custom build step which can be modified to your purpose or removed
	such that you may implement it as you see fit. The dynamically linkable libraries for TBB and FreeImage,
	FreeImage.dll and tbb.dll (or tbb_debug.dll), need to be made available to the binaries compiled with this project.
	This can be done using custom paths or manually moving the files. At present there is a custom build-step that achieves this at :

	Project | Properties | Configuration Properties | Custom Build Step | General


The OpenCL SDK currently in use is AMD's SDK. Thus there are some AMD dependent compiler flags in some of the kernels
that may need to be removed or configured appropriately. The OpenCL kernels will require specific tuning in order to
execute properly. Please refer to your hardware for specifics. The kernels may need to be modified slightly to make
them work. The kernels have been tested with Intel CPU's and with AMD's R9 290 GPU. The kernels have been modified in
the past to work with Nvidia GPU's with success.

Please refer to documentation provided by FreeImage, Intel's TBB, Eigen and the appropriate OpenCL SDK for general questions about
those libraries.


## Related Publications ##

Nathan D. Follin, Keefer D. Taylor, Christopher J. Musalo, and Matthew L. Trawick, “Three-axis correction of distortion due to 
positional drift in scanning probe microscopy,” Review of Scientific Instruments (in press), 2012

Brian S. Salmons, Daniel R. Katz, and Matthew L. Trawick, "Correction of distortion due
 to thermal drift in scanning probe microscopy," Ultramicroscopy, 110, No. 4, 339, 2010.