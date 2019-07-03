
# Setting up OpenCV in Visual Studio 2015

- Download latest OpenCV for Windows from https://opencv.org/releases/. I downloaded v4.1.0 at the time of writing this document.
- Extract the package to a folder of your choice (E.g. C:\Project\Softwares\OpenCV\opencv)
- Download opencv_contrib from https://github.com/opencv/opencv_contrib. Copy the file xfeatures2d.hpp and folder xfeatures2d from opencv_contrib-master\modules\xfeatures2d\include\opencv2 to C:\Project\Softwares\OpenCV\opencv\build\include\opencv2. This is required for using extra features and/or modules like BRISK desciptor 
- Add the path "user-path\opencv\build\x64\vc14\bin" to the environment variable. Log Off and Log On.
- Open VS2015. Create a new empty VC++ project. Set the platform to x64. Create an empty main.cpp.
- Copy a sample code(for e.g. user-path\opencv\sources\samples\cpp\tutorial_code\ImgProc\Morphology_1.cpp) to the main.cpp
- The OpenCV library has to be linked to the VS solution in order to build successfully.
  - Right click "Solution->Properties->C/C++->General" and add the include path "user-path\opencv\build\include" to "Additional include Directories"
  - Right click "Solution->Properties->Linker->General" and add the library path "user-path\opencv\build\x64\vc14\lib" to "Additional Library Directories"
  - Depending on the build configuration, either the opencv_world410d.lib(for Debug) or opencv_world410.lib(for Release) shoould be added to "Solution->Properties->Linker->Input->Additional Dependencies"
  - Create a folder called "data" in the solution directory and place an image file. Modify the code accordingly to read in the proper image file.
  - Build and run the code.
