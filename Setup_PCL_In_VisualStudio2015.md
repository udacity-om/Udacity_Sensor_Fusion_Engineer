# Setting up Point Cloud Library(PCL) in Visual Studio 2015
1. Download the PCL 1.8.0 All-In-One Installer for Microsoft Visual Studio 2015 from https://1drv.ms/u/s!ApoY_0Ymu57sg5QiqO2sR1k-zcSi_w
2. Below are the third party softwares included in the installer
  ![alt text](ReadMe_Images/Third_Party_Softwares_Included_In_The_PCL_Installer.JPG)
  Install the PCL-1.8.0-AllInOne-msvc2015-win64.exe.
3. Also install OpenCV and Freeglut(https://www.transmissionzero.co.uk/files/software/development/GLUT/freeglut-MSVC.zip)
4. Once installed, make sure you have the following environment variables:
    * PCL_ROOT → The install path you chose to install PCL. It should have folders like 3rdParty under Eg: C:\Program Files\PCL 1.8.0
    * OPENNI2_LIB64 → The Lib folder under the place where your OpenNI2 is installed. Eg: C:\Program Files\OpenNI2\Lib
    * OPENNI2_INCLUDE64 → The Include folder under the place where your OpenNI2 is installed. Eg: C:\Program Files\OpenNI2\Include
    * OPENNI2_REDIST64 → The Redist folder under the place where your OpenNI2 is installed. Eg: C:\Program Files\OpenNI2\Redist
    * OPENCV_DIR → Go to your opencv install/unpack directory. Navigate to build → x64 (or 32 if you opted for 32 bits) → vc14 (in case you are using VS 2015). Set the value of the environment value to this folder. Eg: C:\Project\Softwares\OpenCV\opencv\build\x64\vc14
    * FREEGLUT_DIR → Set this to the location where you installed/extracted freeGLUT. Directories like bin and include should be UNDER the pointed folder. Eg: C:\Project\Softwares\freeglut
    Note: Log Off and Log On once the environment variables are set.
5. Launch VS2015 and create an empty project. Add an empty file main.cpp. 
6. Add the required include directories to "Solution->Properties->C/C++->General->Additional Include Directories"
    * $(OPENCV_DIR)\..\..\include
    * $(PCL_ROOT)\include\pcl-1.8\
    * $(PCL_ROOT)\3rdParty\VTK\include\vtk-7.0
    * $(PCL_ROOT)\3rdParty\Boost\include\boost-1_61
    * $(PCL_ROOT)\3rdParty\Qhull\include
    * $(PCL_ROOT)\3rdParty\FLANN\include
    * $(PCL_ROOT)\3rdParty\Eigen\eigen3
    * $(OPENNI2_INCLUDE64)
    * $(FREEGLUT_DIR)\include
7. Add the required library directories to "Solution->Properties->Linker->General->Additional Library Directories"
    * $(OPENCV_DIR)\lib
    * $(PCL_ROOT)\lib
    * $(PCL_ROOT)\3rdParty\VTK\lib
    * $(PCL_ROOT)\3rdParty\Boost\lib
    * $(PCL_ROOT)\3rdParty\Qhull\lib
    * $(PCL_ROOT)\3rdParty\FLANN\lib
    * $(OPENNI2_LIB64)
    * $(FREEGLUT_DIR)\lib\x64
8. Add the required .lib files to "Solution->Properties->Linker->Input->Additional Dependencies". 
   (Tip: Open up a command prompt and for each of the library locations set up above do the following: 
   Navigate to the library location in the command prompt. Enter "dir \*.lib /B > C:\Users\dz31jl\Desktop\store.txt", without the quotes. Copy the contents of the store.txt file and add it to the Additional Dependencies. Make sure to add the right .lib to the corresponding Solution Configuration(For e.g. opencv_world410d.lib should be added to Debug and opencv_world410.lib to Release). The contents of store.txt has to be copied before executing the command again otherwise store.txt will be overwritten.)
9. Write a simple program to check if the settings are correct
```c++
  #include <iostream>

  #include "pcl\conversions.h"

  int main()
  {
    std::cout << "Hello PCL\n";

    std::cin.get();
  }
```
10. Save these settings. File->Export template->Project template. Give a name for the template, for e.g. PCL_Poject_template. Make sure that the "Automatically import the template into Visual Studio" is checked.
11. Create a new project by selecting "PCL_Project_Template" and name it "Lidat_Obstacle_Detection". Get the files from https://github.com/udacity/SFND_Lidar_Obstacle_Detection.
  Below is the file structure of my project:
  ![alt text](ReadMe_Images/PCL_Solution_File_Structure.JPG)
  Make sure to build in x64 platform.
    * Add "\_CRT_SECURE_NO_WARNINGS" and "\_SCL_SECURE_NO_WARNINGS" to "Solution->Properties->C/C++->Preprocessor->Preprocessor Definitions if required for successful build.
    * Place OpenNI2.dll(from C:\Program Files\OpenNI2\Redist) in the project output folder(where the project .exe is generated, user-path\Lidar_Obstacle_Detection\x64\Debug)

