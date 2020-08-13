# Core_BEC_Analysis
**[Bryce M. Henson](https://github.com/brycehenson), [Jacob A. Ross](https://github.com/GroundhogState),[Kieran F. Thomas](https://github.com/KF-Thomas),[David Shin](https://github.com/spicydonkey)**  
A reasonably comprehensive toolset for analysis of data generated in the He* BEC group.  
**Status:** This Code is **ready for use in other projects**. Unit Testing is implemented for **most** functions. Integration/system testing is **not** implemented.

## Install 

As a submodule, most common usage. Usually installed in the lib folder of a project.
```
git submodule add -b dev  https://github.com/brycehenson/Core_BEC_Analysis.git lib/Core_BEC_Analysis
```

Solo install. Not recomended for including Core_BEC_Analysis in a project. It makes install for new users harder and does not point to a particular version of this code making debuging harder.
``` 
git clone --recursive https://github.com/brycehenson/Core_BEC_Analysis.git
```

to update and initalize all the sub-sub-modules
```
git submodule update --init --recursive --remote --merge
```

## To Do
Contributions from students in the He* BEC group or other experienced persons are highly encouraged. Drop me an [email](mailto:bryce.m.henson+github.Core_BEC_Analysis@gmail.com?subject=I%20would%20Like%20to%20Contribute[github][Core_BEC_Analysis]) to gain acess.
- [X] standard project structure
- [X] Combine functions that are found in multiple repos
- [ ] More documentation of the functionality provided


## Contributions  
This project would not have been possible without the many open source tools that it is based on. In no particular order: 

* ***James Conder*** [gaussfilt](https://au.mathworks.com/matlabcentral/fileexchange/43182-gaussfilt-t-z-sigma)
* ***Ander Biguri*** [Perceptually uniform colormaps](https://au.mathworks.com/matlabcentral/fileexchange/51986-perceptually-uniform-colormaps)
* ***Jan*** [FileTime](https://au.mathworks.com/matlabcentral/fileexchange/24671-filetime)
* ***Benjamin Kraus*** [nanconv](https://au.mathworks.com/matlabcentral/fileexchange/41961-nanconv)
* ***M. A. Hopcroft**** [allan](https://au.mathworks.com/matlabcentral/fileexchange/13246-allan)
* ***Daniel Eaton***  [sfigure](https://au.mathworks.com/matlabcentral/fileexchange/8919-smart-silent-figure)
* ***Denis Gilbert***  [M-file Header Template](https://au.mathworks.com/matlabcentral/fileexchange/4908-m-file-header-template)
* ***DrosteEffect***  [CIECAM02](https://github.com/DrosteEffect/CIECAM02)
* ***Bruno Luong*** [histcn](https://au.mathworks.com/matlabcentral/fileexchange/23897-n-dimensional-histogram)
* ***Kenneth Johnson*** [elementwise power](https://au.mathworks.com/matlabcentral/fileexchange/44574-elementwise-power)
* ***James*** [isalmost](https://au.mathworks.com/matlabcentral/fileexchange/15816-isalmost)
* ***Yair Altman*** [export_fig](https://github.com/altmany/export_fig)
* ***Christopher Thissen*** [errorbarxy](https://github.com/cthissen/errorbarxy)

