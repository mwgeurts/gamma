## CalcGamma

by Mark Geurts <mark.w.geurts@gmail.com>
<br>Copyright &copy; 2014, University of Wisconsin Board of Regents

`CalcGamma()` computes 1-D, 2-D, or 3-D global or absolute gamma between two datasets (reference and target) given a defined coordinate space. The datasets must have the same number of dimensions, although they can be different sizes. Gamma will be computed for each target dose point by shifting the reference image (using linear interpolation) and determining the minimum Gamma index across all shifts.

This function optionally uses the Parallel Computing Toolbox GPU interp functions to increase computation speed. A try-catch statement is used to test for GPU support. In addition, for memory management, the meshgrid and data arrays are converted to single precision during interpolation.

## Contents

* [MATLAB Function Use](README.md#matlab-function-use)
* [Example](README.md#example)
* [Gamma Computation Methods](README.md#gamma-computation-methods)
* [Compatibility and Requirements](README.md#compatibility-and-requirements)
* [Event Calling](README.md#event-calling)
* [License](README.md#license)

## MATLAB Function Use

The following variables are required for proper execution: 

* varargin{1}: structure containing the reference data, where the field start is an array containing the coordinates along each dimension of the first voxel, width is an array containing the width of each voxel along each dimension, and data is an n-dimensional array
* varargin{2}: structure containing the target data, where the field start is an array containing the coordinates along each dimension of the first voxel, width is an array containing the width of each voxel along each dimension, and data is an n-dimensional array
* varargin{3}: Gamma absolute criterion percentage
* varargin{4}: Gamma Distance To Agreement (DTA) criterion, in the same units as the reference and target width structure fields * varargin{5} (optional): boolean, indicates whether to perform a local (1) or global (0) Gamma computation.  If not present, the function will assume a global Gamma computation.
* varargin{6} (optional): reference value for the global absolute criterion.  Is used with the percentage from varargin{3} to compute absolute value.  If not present, the maximum value in the reference data is used.
* varargin{7} (optional): restricted search flag. If 1, only the gamma values along the X/Y/Z axes are computed during 3D comptation. If 0 or not provided, the entire rectangular search space is computed.

The following variables are returned upon succesful completion:

* gamma: array of the same dimensions as varargin{2}.data containing the computed gamma values

## Example

```matlab
reference.start = [-10 -10]; % mm
reference.width = [0.1 0.1]; % mm
reference.data = rand(200);

target.start = [-10 -10]; % mm
target.width = [0.1 0.1]; % mm
target.data = rand(200);

percent = 3;
dta = 0.5; % mm
local = 0; % Perform global gamma
   
gamma = CalcGamma(reference, target, percent, dta, local);
```

## Gamma Computation Methods

The Gamma analysis is performed based on the formalism presented by D. A. Low et. al., [A technique for the quantitative evaluation of dose distributions.](http://www.ncbi.nlm.nih.gov/pubmed/9608475), Med Phys. 1998 May; 25(5): 656-61.  In this formalism, the Gamma quality index *&gamma;* is defined as follows for each point along the measured profile *Rm* given the reference profile *Rc*:

*&gamma; = min{&Gamma;(Rm,Rc}&forall;{Rc}*

where:

*&Gamma; = &radic; (r^2(Rm,Rc)/&Delta;dM^2 + &delta;^2(Rm,Rc)/&Delta;DM^2)*,

*r(Rm,Rc) = | Rc - Rm |*,

*&delta;(Rm,Rc) = Dc(Rc) - Dm(Rm)*,

*Dc(Rc)* and *Dm(Rm)* represent the reference and measured signal at each *Rc* and *Rm*, respectively, and

*&Delta;dM* and *&Delta;DM* represent the absolute and Distance To Agreement Gamma criterion (by default 2%/1mm), respectively.  

The absolute criterion is typically given in percent and can refer to a percent of the maximum dose (commonly called the global method) or a percentage of the voxel *Rm* being evaluated (commonly called the local method).  The application is capable of computing gamma using either approach, and can be set when calling CalcGamma.m by passing a boolean value of 1 (for local) or 0 (for global).  By default, the global method (0) is used.

The computation applied in the tool is the 1D algorithm, in that the distance to agreement criterion is evaluated only along the dimension of the reference profile when determining *min{&Gamma;(Rm,Rc}&forall;{Rc}*. To accomplish this, the reference profile is shifted relative to the measured profile using linear 1D CUDA (when available) interpolation.  For each shift, *&Gamma;(Rm,Rc}* is computed, and the minimum value *&gamma;* is determined.  To improve computation efficiency, the computation space *&forall;{Rc}* is limited to twice the distance to agreement parameter.  Thus, the maximum "real" Gamma index returned by the application is 2.

## Compatibility and Requirements

This tool has been tested with MATLAB 8.3 and 8.4.  The Parallel Computing toolbox (versions 6.4 and 6.5 tested and a CUDA-compatible GPU are required to run GPU based interpolation (CPU interpolation is automatically supported if not present).

## Event Calling

These functions optionally return execution status and error information to an `Event()` function. If available in the MATLAB path, `Event()` will be called with one or two variables: the first variable is a string containing the status information, while the second is the status classification (WARN or ERROR). If the status information is only informative, the second argument is not included.  Finally, if no `Event()` function is available errors will still be thrown via the standard `error()` MATLAB function.


## License

This program is free software: you can redistribute it and/or modify it 
under the terms of the GNU General Public License as published by the  
Free Software Foundation, either version 3 of the License, or (at your 
option) any later version.

This program is distributed in the hope that it will be useful, but 
WITHOUT ANY WARRANTY; without even the implied warranty of 
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General 
Public License for more details.

You should have received a copy of the GNU General Public License along 
with this program. If not, see http://www.gnu.org/licenses/.
