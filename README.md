# SCALE_DTW
## Environment and tools
+ Software environment  
  windows 11
+ Development and debugging tools  
  MinGW;Clion 2023.2.2;Cmake
## DEPENDENCY
+ FFTW library
  Considering the need for convolutional computation, we therefore need to install the FFTW library.
  Download the installation package of FFTW library [here](https://fftw.org/install/windows.html)

  then open the "x64 Native Tools Command Prompt for VS 2022",Go to the decompressed folder and run the following commands:
  
  ```
  lib /machine:x64 /def:libfftw3f-3.def  
  lib /machine:x64 /def:libfftw3-3.def  
  lib /machine:x64 /def:libfftw3l-3.def
   ```
  Add the **fftw3.h** file to the clion project folder;  
  Add the **libfftw3-3.dll**, **libfftw3f-3.dll**,**libfftw3l-3.dll** files to the Executable directory;  
  Add **libfftw3-3.lib**, **libfftw3f-3.lib**, **libfftw3l-3.lib** to the lib folder of the clion project.

  Then add the bin folder under the clion project folder to the environment variables.
  In CMakeList.txt,add the following lines:
  ```
  set(LIB_DIR "${CMAKE_CURRENT_SOURCE_DIR}/lib")
  link_directories(${LIB_DIR})
  target_link_libraries(your_project_name libfftw3-3.lib libfftw3f-3.lib  libfftw3l-3.lib)
  ```
  After completing the above steps, you are ready to use the fftw library in this clion project.
  If there is an issue with linking the dynamic libraries, you can add the three DLL libraries to the system's PATH environment variable.
## INSTUCTIONS OF THE CODE
+ Both the query and the data to be searched should be formatted,where the numbers should be separated by spaces(' ') or line breaks('\n').  
  Before the start of the searching,we should turn the data to be searched into **binary** format;
  These two formats can be converted using two.exe programs in the **binaryfile** folder.  
  the usage is as follows:
  ```
  file2binary.exe normal-file-path binary-file-path
  readbinary.exe binary-file-path normal-file-path
  ```
  the normal-file-path points to an ASCII or Unicode file.
  the binary-file-path points to the path where you store the binary file.

+ Configuration  
  you should set the Path to the query file and the data to be searched file in **def.h** ( which defines all file paths).

  + For Query: 
  ```
  #define Q "path-to-your-query-file"
  ```
  Q is an ASCII or Unicode Query sequence file,so map your path to it.
  
  + For the long series to be searched:
  ```
  #define BI_T "path-to-your-queried-binaryfile"
  ```
  Note the long series must be in the binary format, which can be obtained through file2binary.exe as mentioned above.

+ The executable file needs 3 parameters:
  1. argv[1] is the length of query subsequence, which corresponds to $m$ in the paper.
  2. argv[2] is the parameter controlling the width of warping path for Sakoe-Chiba band,
where w is set to  $\lceil m * argv[2] \rceil$.
  3. argv[3] is the distance threshold for the subsequence matching problem, corresponds to $\epsilon$ in the paper.
+ For Example, 
  input:
  ```
  scale_dtw_optimize.exe 1024 0.02 0.5
  ```
  is to perform query with $\epsilon = 0.5$, $\lceil w= 1024* 0.02 \rceil =21$, the length of query $m=1024$ and the files storing the query and the series to be queryied are setted in "./def.h".
+ Some examples of queries and data series  are in "./Data/" (from the ECG dataset).
  modify  "./def.h"  according to your needs.
  

  
  

  
