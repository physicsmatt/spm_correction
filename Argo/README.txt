Compile using Visual Studio 2013.
To include updated Eigen source code from http://eigen.tuxfamily.org modify as such :

ArgoThreaded -> Properties -> Configuration Properties -> C/C++ -> General -> Additional Include Directories -> $PATH$

Where if your Eigen code was in H:\..\eigen_3.0.2\Eigen\
$PATH$ would correspond to H:\..\eigen_3.0.2\

This informs the compiler of the proper location of all Eigen dependency header files.