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
  Add the **libfftw3-3.dll**, **libfftw3f-3.dll**,**libfftw3l-3.dll** files to the bin folder of the clion project;  
  Add **libfftw3-3.lib**, **libfftw3f-3.lib**, **libfftw3l-3.lib** to the lib folder of the clion project.

  Then add the bin folder under the clion project folder to the environment variables.
  In CMakeList.txt,add the following lines:
  ```
  set(MY_BIN_DIR "${CMAKE_BINARY_DIR}/bin")
  set(LIB_DIR "${CMAKE_CURRENT_SOURCE_DIR}/lib")
  link_directories(${LIB_DIR})
  target_link_libraries(your_project_name libfftw3-3.lib libfftw3f-3.lib  libfftw3l-3.lib)
  ```
  After completing the above steps, you are ready to use the fftw library in this clion project.
  If it does not work,you need to configure the FFTW library yourself.
## Usage Method 
+ File format
  Sequence data should be separated by Spaces or line breaks and stored in.txt files  
  The query sequence is stored as a text type, while the queried sequence is stored as a **binary** format
  These two formats can be converted using two.exe programs in the **binaryfile** folder.

  the usage is as follows:
  ```
  file2binary.exe origin-file-path binary-file-path
  readbinary.exe binary-file-path normal-file-path
  ```
+ File reading
  you have to change the Path to the query sequence file and queried sequence file in **def.h**
  + Query sequence
  ```
  #define Q0 
  #define Q1 
  #define Q2 
  #define Q3 
  #define Q4 
  #define Q5 
  #define Q6 
  #define Q7 
  #define Q8 
  #define Q9  
  ```
  Q0~Q9 are normal .txt Query sequence files,so map your paths to them.
  
  + Queried sequence file
  ```
  #define BI_T
  ```
  you MUST put the path to zhe **Binary-read-and-write** file here as the Queried sequence file.

+ Input parameter
  there are totally 4 parameters you need to input.
  1. argv[1] input the length of query subsequence.
  2. argv[2] input the value of the relative width of warping path for Sakoe-Chiba band.
  3. argv[2] input the threshold of DTW values.
  4. argv[4] input the number of query subsequence,deciding which of the Q0-9 options to choose.

+ For Instance
  input:
  ```
  scale_dtw_optimize.exe 1024 0.02 0.5 0
  ```
  as an example,we use the data in the Data folder as the default query and the queried sequence, in this case the ECG data. You can change the data set in the def.h file according to your needs.
  
+ Supplementary statement
  + All functions are defined in fun.h.
  + All file paths are in def.h.
  
  

  
