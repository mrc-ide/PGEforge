## Installation

The C source code must be compiled before use, and requires no special libraries. If you do not already have a C compiler installed, you will have to install one. (Note if you are looking for a compiler in your environment: common names you might see for a compiler include CC and GCC.) The installation procedure varies with operation system:
<details>
  <summary>Mac</summary>
  
 1. Open a Terminal window. You can do this by searching (keyboard shortcut command+space) for "Terminal".
  2. Inside the terminal window, type the command `xcode-select --install`. This will prompt you to install Xcode command line tools.
  3. Once you have Xcode command line tools installed you should be able to compile.
</details>

<details>
  <summary>Linux</summary>
    Follow these directions: https://phoenixnap.com/kb/install-gcc-ubuntu.

</details>
<details>
  <summary>Windows</summary>
    Installation is complicated on Windows. These two guides may help: https://www.guru99.com/c-gcc-install.html and https://dev.to/gamegods3/how-to-install-gcc-in-windows-10-the-easier-way-422j.

</details>
If your compiler is named 'CC', the following command 

```
cc -o hmmIBD -O3 -lm -Wall hmmIBD.c
```

should yield the executable file *hmmIBD*.

If you wish to use them, the Python scripts for extracting VCF data and thinning sites require that Python be installed, version 3.6 or later. Mac and most versions of Linux already have Python3 installed. For Windows, directions can be found [here](https://www.codecademy.com/article/install-python3).
