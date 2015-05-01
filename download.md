# Downloads #

The Google Code platform does not provide anymore the ability to offer direct downloads.

| Please, find the latest version of the source package or the binary installers for Windows on [Pypi](https://pypi.python.org/pypi/pydsm/) |
|:------------------------------------------------------------------------------------------------------------------------------------------|

Past versions of the source package and binary installers are available on [Google Drive](https://drive.google.com/folderview?id=0B0XINMKmLTjHWWFqUE5sRHFCZ0E&usp=sharing).

As far as you can please prefer [Pypi](https://pypi.python.org/pypi/pydsm/) for your downloads. Google enforces some restrictions to public downloads from Google Drive. Excessive traffic may end up in downloads being refused for some time.


## Files for installing on Linux ##

To install on Linux, there is no need to explicitly download the package. The preferred way of installing is now via `pip`.
```
   pip install pydsm
```
This downloads the source package, compiles its binary parts and installs on the system.


## Files for installing on MacOs ##

To install on MacOs, there is no need to explicitly download the package. The preferred way of installing is now via `pip`.
```
   pip install pydsm
```
This downloads the source package, compiles its binary parts and installs on the system.  For the procedure to work Xcode must be present on your system.


## Files for installing on Windows ##

To install on Windows, the `pip` installer works as in Linux (see above) _provided that you have the C compiler from Microsoft Visual Studio 2008_. The _Express_ Edition of Visual Studio 2008 is a free download from Microsoft. Note that in any case there are some complications in building from source in Windows (refer to the documentation for details).

If you prefer not to build from sources, we provide binary installers for 32 and 64 bit versions of Windows as a convenience. These are the `.exe` files found on [Pypi](https://pypi.python.org/pypi/pydsm/) (or on
[Google Drive](https://drive.google.com/folderview?id=0B0XINMKmLTjHWWFqUE5sRHFCZ0E&usp=sharing) for less recent versions).
Just download the installer and execute it (or drop it on the _Control Panel_ if you are using [WinPython](http://winpython.sourceforge.net/)).

Note that **your mileage with the binary installers may vary**. Pre-compiled code is very sensitive to the rest of the software that is present on your system and that it needs to work with. The code in the binary installers may work with a specific Python distribution and not with another, or may work if your system uses a specific `numpy` version and not with another.  We test the installers with recent versions of [WinPython](http://winpython.sourceforge.net/).
In any case, consider that we develop on Linux and that our ability to test the code in Windows is somehow limited. Be so kind to report any issues that you may encounter.


## Precompiled documentation ##

Online documentation is available for consultation on [Pypi](http://pythonhosted.org//pydsm).

Pre-built `html` documentation can also be found packaged on a `zip` file on  [Google Drive](https://drive.google.com/folderview?id=0B0XINMKmLTjHWWFqUE5sRHFCZ0E&usp=sharing).


## Pre-requisites necessary to use pydsm versions before and including 0.8.x ##

Very old versions of PyDSM need a Python package that has been retired. Please avoid these versions. If you really need to use them, the missing pre-requisite can be found on [Google Drive](https://drive.google.com/folderview?id=0B0XINMKmLTjHWWFqUE5sRHFCZ0E&usp=sharing).