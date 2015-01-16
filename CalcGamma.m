function gamma = CalcGamma(varargin)
% CalcGamma computes 1-D, 2-D, or 3-D global or absolute gamma between two
% datasets (reference and target) given a defined coordinate space. The 
% datasets must have the same number of dimensions, although they can be 
% different sizes. Gamma will be computed for each target dose point by
% shifting the reference image (using linear interpolation) and determining
% the minimum Gamma index across all shifts.
%
% This function optionally uses the Parallel Computing Toolbox GPU interp
% functions to increase computation speed. A try-catch statement is used
% to test for GPU support. In addition, for memory management, the
% meshgrid and data arrays are converted to single precision during
% interpolation. This function calls Event.m to log execution status, if 
% available.
%
% For more information on the Gamma evaluation function, see D. A. Low et 
% al., "A technique for the quantitative evaluation of dose distributions", 
% Med Phys. 1998 May; 25(5): 656-61.
%
% The following variables are required for proper execution: 
%   varargin{1}: structure containing the reference data, where the field
%       start is an array containing the coordinates along each dimension
%       of the first voxel, width is an array containing the width of each
%       voxel along each dimension, and data is an n-dimensional array
%   varargin{2}: structure containing the target data, where the field
%       start is an array containing the coordinates along each dimension
%       of the first voxel, width is an array containing the width of each
%       voxel along each dimension, and data is an n-dimensional array
%   varargin{3}: Gamma absolute criterion percentage
%   varargin{4}: Gamma Distance To Agreement (DTA) criterion, in the same
%       units as the reference and target width structure fields  
%   varargin{5} (optional): boolean, indicates whether to perform a local 
%       (1) or global (0) Gamma computation.  If not present, the function
%       will assume a global Gamma computation.
%   varargin{6} (optional): reference value for the global absolute 
%       criterion.  Is used with the percentage from varargin{3} to compute
%       absolute value.  If not present, the maximum value in the reference
%       data is used.
%   varargin{7} (optional): restricted search flag. If 1, only the gamma 
%       values along the X/Y/Z axes are computed during 3D comptation. If 
%       0 or not provided, the entire rectangular search space is computed.
%
% The following variables are returned upon succesful completion:
%   gamma: array of the same dimensions as varargin{2}.data containing the
%       computed gamma values
%
% Below is an example of how the function is used:
%
%   reference.start = [-10 -10]; % mm
%   reference.width = [0.1 0.1]; % mm
%   reference.data = rand(200);
%
%   target.start = [-10 -10]; % mm
%   target.width = [0.1 0.1]; % mm
%   target.data = rand(200);
%
%   percent = 3;
%   dta = 0.5; % mm
%   local = 0; % Perform global gamma
%   
%   gamma = CalcGamma(reference, target, percent, dta, local);
%
% Author: Mark Geurts, mark.w.geurts@gmail.com
% Copyright (C) 2014 University of Wisconsin Board of Regents
%
% This program is free software: you can redistribute it and/or modify it 
% under the terms of the GNU General Public License as published by the  
% Free Software Foundation, either version 3 of the License, or (at your 
% option) any later version.
%
% This program is distributed in the hope that it will be useful, but 
% WITHOUT ANY WARRANTY; without even the implied warranty of 
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General 
% Public License for more details.
% 
% You should have received a copy of the GNU General Public License along 
% with this program. If not, see http://www.gnu.org/licenses/.

% Log initialization and start timer
if exist('Event', 'file') == 2
    Event('Beginning Gamma calculation');
    tic;
end


% Check if the reference structure contains width, start, and data fields,
% and if the size of the width and start vectors are equal
if ~isfield(varargin{1}, 'width') || ~isfield(varargin{1}, 'start') || ...
        ~isfield(varargin{1}, 'data') || ~isequal(size(varargin{1}.width), ...
        size(varargin{1}.start))
    
    % If not, throw an error and stop execution
    if exist('Event', 'file') == 2
        Event(['Incorrect reference data format.  Must contain width, ', ...
            'start, and data fields and be of equal dimensions'], 'ERROR');
    else
        error(['Incorrect reference data format.  Must contain width, ', ...
            'start, and data fields and be of equal dimensions']);
    end
    
% Check if the target structure contains width, start, and data fields,
% and if the size of the width and start vectors are equal
elseif ~isfield(varargin{2}, 'width') || ~isfield(varargin{2}, 'start') || ...
        ~isfield(varargin{2}, 'data') || ~isequal(size(varargin{2}.width), ...
        size(varargin{2}.start))
    
    % If not, throw an error and stop execution
    if exist('Event', 'file') == 2
        Event(['Incorrect target data format.  Must contain width, ', ...
            'start, and data fields and be of equal dimensions'], 'ERROR');
    else
        error(['Incorrect target data format.  Must contain width, ', ...
            'start, and data fields and be of equal dimensions']);
    end
    
% Check if the reference and target data arrays are the same number of
% dimensions.  Calculating the gamma from a lower dimensional dataset to a
% higher dimensional reference is currently not supported
elseif ~isequal(size(size(varargin{1}.data)), size(size(varargin{2}.data)))
    
    % If not, throw an error and stop execution
    if exist('Event', 'file') == 2
        Event(['The fixed and target data arrays must be of the same', ...
            ' dimensions'], 'ERROR');
    else
        error(['The fixed and target data arrays must be of the same', ...
            ' dimensions']);
    end
end

% Log validation completed
if exist('Event', 'file') == 2
    Event('Data validation completed');
end

% If a local/global Gamma flag was not provided
if nargin < 4
    % Assume the computation is global
    varargin{5} = 0;
    
    % Log type
    if exist('Event', 'file') == 2
        Event('Gamma calculation assumed to global');
    end
elseif varargin{5} == 0  
    % Log type
    if exist('Event', 'file') == 2
        Event('Gamma calculation set to global');
    end
elseif exist('Event', 'file') == 2
    % Log type
    Event('Gamma calculation set to local');  
end

% If a reference absolute value was not provided
if nargin < 6
    % Assume the reference value is the maximum value in the reference data
    % array (ie, Gamma % criterion is % of max value)
    varargin{6} = max(max(max(varargin{1}.data)));
    if exist('Event', 'file') == 2
        Event(['No reference value was provided, maximum value in ', ...
            'dataset used']);
    end
end

% If a restrict search flag was not provided
if nargin < 7
    % Disable restricted search
    varargin{7} = 0;
    
    % Log result
    if exist('Event', 'file') == 2
        Event('Restricted search disabled by default');
    end
end

% If the reference dataset is 1-D
if size(varargin{1}.width,2) == 1
    if exist('Event', 'file') == 2
        Event('Reference dataset is 1-D');
    end
    
    % Check if the data is in rows or columns (this is only needed for 1-D)
    if size(varargin{1}.data,1) > size(varargin{1}.data,2)
        
        % If in rows, transpose
        varargin{1}.data = varargin{1}.data';
    end
    
    % Compute the reference X coordinates using the start and width values
    refX = single(varargin{1}.start(1):varargin{1}.width(1):varargin{1}.start(1) ...
        + varargin{1}.width(1) * (size(varargin{1}.data,2) - 1));
    
% Otherwise, if the reference dataset is 2-D
elseif size(varargin{1}.width,2) == 2
    if exist('Event', 'file') == 2
        Event('Reference dataset is 2-D');
    end
    
    % Compute X and Y meshgrids for the reference dataset positions using 
    % the start and width values
    [refX, refY] = meshgrid(single(varargin{1}.start(1):varargin{1}.width(1): ...
        varargin{1}.start(1) + varargin{1}.width(1) * ...
        (size(varargin{1}.data,1) - 1)), single(varargin{1}.start(2): ...
        varargin{1}.width(2):varargin{1}.start(2)...
        + varargin{1}.width(2) * (size(varargin{1}.data,2) - 1)));
    
% Otherwise, if the reference dataset is 3-D
elseif size(varargin{1}.width,2) == 3
    if exist('Event', 'file') == 2
        Event('Reference dataset is 3-D');
    end
    
    % Compute X, Y, and Z meshgrids for the reference dataset positions
    % using the start and width values, permuting X/Y
    [refX, refY, refZ] = meshgrid(single(varargin{1}.start(2): ...
        varargin{1}.width(2):varargin{1}.start(2) + varargin{1}.width(2) * ...
        (size(varargin{1}.data,2) - 1)), single(varargin{1}.start(1): ...
        varargin{1}.width(1):varargin{1}.start(1) + varargin{1}.width(1)...
        * (size(varargin{1}.data,1) - 1)), single(varargin{1}.start(3):...
        varargin{1}.width(3):varargin{1}.start(3) + varargin{1}.width(3)...
        * (size(varargin{1}.data,3) - 1)));

% Otherwise, if the reference data is of higher dimension
else
    % Throw an error and stop execution
    if exist('Event', 'file') == 2
        Event('The fixed data structure contains too many dimensions', ...
            'ERROR');
    else
        error('The fixed data structure contains too many dimensions');
    end
end

% If the target dataset is 1-D
if size(varargin{2}.width,2) == 1
    if exist('Event', 'file') == 2
        Event('Target dataset is 1-D');
    end
    
    % Check if the data is in rows or columns (this is only needed for 1-D)
    if size(varargin{2}.data,1) > size(varargin{2}.data,2)
        
        % If in rows, transpose
        varargin{2}.data = varargin{2}.data';
    end
    
    % Compute the target X coordinates using the start and width values
    tarX = single(varargin{2}.start(1):varargin{2}.width(1):varargin{2}.start(1) ...
        + varargin{2}.width(1) * (size(varargin{2}.data,2) - 1));
    
% Otherwise, if the target dataset is 2-D
elseif size(varargin{2}.width,2) == 2
    if exist('Event', 'file') == 2
        Event('Target dataset is 2-D');
    end
    
    % Compute X and Y meshgrids for the target dataset positions using the
    % start and width values
    [tarX, tarY] = meshgrid(single(varargin{2}.start(1):varargin{2}.width(1): ...
        varargin{2}.start(1) + varargin{2}.width(1) * ...
        (size(varargin{2}.data,1) - 1)), single(varargin{2}.start(2): ...
        varargin{2}.width(2):varargin{2}.start(2)...
        + varargin{2}.width(2) * (size(varargin{2}.data,2) - 1)));
    
% Otherwise, if the target dataset is 3-D
elseif size(varargin{2}.width,2) == 3
    if exist('Event', 'file') == 2
        Event('Target dataset is 3-D');
    end
    
    % Compute X, Y, and Z meshgrids for the target dataset positions using
    % the start and width values, permuting X/Y
    [tarX, tarY, tarZ] = meshgrid(single(varargin{2}.start(2):...
        varargin{2}.width(2):varargin{2}.start(2) + varargin{2}.width(2) * ...
        (size(varargin{2}.data,2) - 1)), single(varargin{2}.start(1): ...
        varargin{2}.width(1):varargin{2}.start(1) + varargin{2}.width(1) ...
        * (size(varargin{2}.data,1) - 1)), single(varargin{2}.start(3):...
        varargin{2}.width(3):varargin{2}.start(3) + varargin{2}.width(3) ...
        * (size(varargin{2}.data,3) - 1)));
    
% Otherwise, if the reference data is of higher dimension
else
    % Throw an error and stop execution
    if exist('Event', 'file') == 2
        Event('The target data structure contains too many dimensions', ...
            'ERROR');
    else
        error('The target data structure contains too many dimensions');
    end
end

% The resolution parameter determines the number of steps (relative to 
% the distance to agreement) that each reference voxel will be
% interpolated to and gamma calculated.  A value of 5 with a DTA of 3
% mm means that gamma will be calculated at intervals of 3/5 = 0.6 mm.
% Different resolutions can be set for different dimensions of data. 
if size(varargin{2}.width,2) == 1
    % Set 1-D resolution
    res = 100;
elseif size(varargin{2}.width,2) == 2
    % Set 2-D resolution
    res = 100;
elseif size(varargin{2}.width,2) == 3
    % Set 3-D resolution
    res = 20;
end

% Log resolution
if exist('Event', 'file') == 2
    Event(sprintf('Interpolation resolution set to %i', res));
end

% Generate an initial gamma volume with values of 2 (this is the maximum
% reliable value of gamma).
gamma = ones(size(varargin{2}.data)) * 2;

% Log number of gamma calculations (for status updates on 3D calcs)
if varargin{7} == 1
    % Compute number of restricted search calcs
    num = res * 4 * size(varargin{2}.width,2);
else
    % Compute total number of calcs
    num = res * 4 ^ size(varargin{2}.width,2);
end

% num is the number of iterations, num * numel the total number of
% interpolations being performed
if exist('Event', 'file') == 2
    Event(sprintf('Number of gamma calculations = %g', num * numel(gamma)));
end

% Initialize counter (for progress indicator)
n = 0;

% Start try-catch block to safely test for CUDA functionality
try
    % Clear and initialize GPU memory.  If CUDA is not enabled, or if the
    % Parallel Computing Toolbox is not installed, this will error, and the
    % function will automatically rever to CPU computation via the catch
    % statement
    gpuDevice(1);
    
    % Start a for loop to interpolate the dose array along the x-direction.  
    % Note to support parfor loops indices must be integers, so x varies 
    % from -2 to +2 multiplied by the number of interpolation steps.  
    % Effectively, this evaluates gamma from -2 * DTA to +2 * DTA.
    for x = -2*res:2*res
        
        % i is the x axis step value
        i = x/res * varargin{4};
        
        % Initialize j and k as zero (they will be updated if the data is
        % of higher dimension)
        j = 0;
        k = 0;
        
        % If the data contains a second dimension
        if size(varargin{1}.width,2) > 1
   
            % Start a for loop to interpolate the dose array along the
            % y-direction.  Note to support parfor loops indices must be
            % integers, so y varies from -2 to +2 multiplied by the number
            % of interpolation steps.  Effectively, this evaluates gamma
            % from -2 * DTA to +2 * DTA.
            for y = -2*res:2*res
                
                % j is the y axis step value
                j = y/res * varargin{4};
                
                % Initialize k as zero (it will be updated if the data is
                % of higher dimension)
                k = 0;
                
                % If the data contains a third dimension
                if size(varargin{1}.width,2) > 2
                    
                    % Start a for loop to interpolate the dose array along 
                    % the z-direction.  Note to support parfor loops 
                    % indices must be integers, so z varies from -2 to +2 
                    % multiplied by the number of interpolation steps.
                    % Effectively, this evaluates gamma from -2 * DTA to 
                    % +2 * DTA.
                    for z = -2*res:2*res
                        
                        % k is the z axis step value
                        k = z/res * varargin{4};

                        % Check restricted search flag
                        if varargin{7} == 0 || sum(abs([x y z]) > 0) == 1
                            
                            % Run GPU interp3 function to compute the reference
                            % values at the specified target coordinate points
                            interp = gather(interp3(gpuArray(refX), gpuArray(refY), ...
                                gpuArray(refZ), gpuArray(single(varargin{1}.data)), ...
                                gpuArray(tarX + i), gpuArray(tarY + j), ...
                                gpuArray(tarZ + k), 'linear', 0));

                            % Update the gamma array by returning the minimum
                            % of the existing value or the new value
                            gamma = min(gamma, GammaEquation(interp, ...
                                varargin{2}.data, i, j, k, varargin{3}, varargin{4}, ...
                                varargin{6}, varargin{5}));
                            
                            % Update counter 
                            n = n + 1;
                            
                            % If counter is at an even %, display progress
                            if mod((n-1)/num, 0.01) > 0.005 && ...
                                    mod(n/num, 0.01) < 0.005
                                fprintf('%0.1f%%\n', n/num*100);
                            end
                        end
                    end
                    
                % Otherwise, the data is 2-D
                else
                    % Run GPU interp2 function to compute the reference
                    % values at the specified target coordinate points
                    interp = gather(interp2(gpuArray(refX), gpuArray(refY), ...
                        gpuArray(single(varargin{1}.data)), gpuArray(tarX + i), ...
                        gpuArray(tarY + j), 'linear', 0));
                    
                    % Update the gamma array by returning the minimum
                    % of the existing value or the new value
                    gamma = min(gamma, GammaEquation(interp, varargin{2}.data, ...
                        i, j, k, varargin{3}, varargin{4}, ...
                        varargin{6}, varargin{5}));
                end
            end
            
        % Otherwise, the data is 1-D
        else
            % Run GPU interp function to compute the reference values at 
            % the specified target coordinate points
            interp = gather(interp1(gpuArray(refX), ...
                gpuArray(single(varargin{1}.data)), gpuArray(tarX + i), ...
                'linear', 0));
            
            % Update the gamma array by returning the minimum of the 
            % existing value or the new value
            gamma = min(gamma, GammaEquation(interp, varargin{2}.data, ...
                i, j, k, varargin{3}, varargin{4}, ...
                varargin{6}, varargin{5}));
        end
    end
   
% If GPU fails, revert to CPU computation
catch
    % Log GPU failure
    if exist('Event', 'file') == 2
        Event('GPU failed, reverting to CPU method', 'WARN'); 
    end
    
    % Start a for loop to interpolate the dose array along the x-direction.  
    % Note to support parfor loops indices must be integers, so x varies 
    % from -2 to +2 multiplied by the number of interpolation steps.  
    % Effectively, this evaluates gamma from -2 * DTA to +2 * DTA.
    for x = -2*res:2*res
        % i is the x axis step value
        i = x/res * varargin{4};
        
        % Initialize j and k as zero (they will be updated if the data is
        % of higher dimension)
        j = 0;
        k = 0;
        
        % If the data contains a second dimension
        if size(varargin{1}.width,2) > 1
            
            % Start a for loop to interpolate the dose array along the
            % y-direction.  Note to support parfor loops indices must be
            % integers, so y varies from -2 to +2 multiplied by the number
            % of interpolation steps.  Effectively, this evaluates gamma
            % from -2 * DTA to +2 * DTA.
            for y = -2*res:2*res
                
                % j is the y axis step value
                j = y/res * varargin{4};
                
                % Initialize k as zero (it will be updated if the data is
                % of higher dimension)
                k = 0;
                
                % If the data contains a third dimension
                if size(varargin{1}.width,2) > 2
                    
                    % Start a for loop to interpolate the dose array along 
                    % the z-direction.  Note to support parfor loops 
                    % indices must be integers, so z varies from -2 to +2 
                    % multiplied by the number of interpolation steps.
                    % Effectively, this evaluates gamma from -2 * DTA to 
                    % +2 * DTA.
                    for z = -2*res:2*res
                        
                        % k is the z axis step value
                        k = z/res * varargin{4};

                        % Check restricted search flag
                        if varargin{7} == 0 || sum(abs([x y z]) > 0) == 1
                            
                            % Run CPU interp3 function to compute the reference
                            % values at the specified target coordinate points
                            interp = interp3(refX, refY, refZ, ...
                                single(varargin{1}.data), tarX + i, ...
                                tarY + j, tarZ + k, '*linear', 0);

                            % Update the gamma array by returning the minimum
                            % of the existing value or the new value
                            gamma = min(gamma, GammaEquation(interp, ...
                                varargin{2}.data, i, j, k, varargin{3}, ...
                                varargin{4}, varargin{6}, varargin{5}));
                            
                            % Update counter 
                            n = n + 1;
                            
                            % If counter is at an even %, display progress
                            if mod((n-1)/num, 0.01) > 0.005 && ...
                                    mod(n/num, 0.01) < 0.005
                                fprintf('%0.1f%%\n', n/num*100);
                            end
                        end
                    end
                    
                % Otherwise, the data is 2-D
                else
                    % Run CPU interp2 function to compute the reference
                    % values at the specified target coordinate points
                    interp = interp2(refX, refY, single(varargin{1}.data), ...
                        tarX + i, tarY + j, '*linear', 0);
                    
                    % Update the gamma array by returning the minimum
                    % of the existing value or the new value
                    gamma = min(gamma, GammaEquation(interp, ...
                        varargin{2}.data, i, j, k, varargin{3}, ...
                        varargin{4}, varargin{6}, varargin{5}));
                end
            end
            
        % Otherwise, the data is 1-D
        else
            % Run CPU interp function to compute the reference values at 
            % the specified target coordinate points
            interp = interp1(refX, single(varargin{1}.data), tarX + i, ...
                '*linear', 0);
            
            % Update the gamma array by returning the minimum of the 
            % existing value or the new value
            gamma = min(gamma, GammaEquation(interp, varargin{2}.data, ...
                i, j, k, varargin{3}, varargin{4}, ...
                varargin{6}, varargin{5}));
        end
    end
end
    
% Log completion
if exist('Event', 'file') == 2
    Event(sprintf(['Gamma calculation completed successfully in ', ...
        '%0.3f seconds'], toc));
end

% Clear temporary variables
clear refX refY refZ tarX tarY tarZ interp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function gamma = GammaEquation(ref, tar, i, j, k, perc, dta, refval, local)
% GammaEquation is the programmatic form of the Gamma definition as given
% by Low et al in matrix form.  This function computes both local and
% global Gamma, and is a subfunction for CalcGamma.
%
% The following inputs are used for computation and are required:
%   ref: the reference 3D array.  Must be the same size as tar
%   tar: the target 3D array.  Must be the same size as ref
%   i: magnitude of x position offset of tar to ref, relative to dta
%   j: magnitude of y position offset of tar to ref, relative to dta
%   k: magnitude of z position offset of tar to ref, relative to dta
%   perc: the percent Gamma criterion, given in % (i.e. 3 for 3%)
%   dta: the distance to agreement Gamma criterion, unitless but relative
%       to i, j, and k
%   refval: if global, the reference value to base the % criterion from 
%   local: boolean, indicates whether to perform a local (1) or global (0)
%       Gamma computation
%
% The following variables are returned:
%   gamma: a 3D array of the same dimensions as ref and interp of the
%       computed gamma value for each voxel based on interp and i,j,k
%

% If local is set to 1, perform a local Gamma computation
if local == 1
    % Gamma is defined as the sqrt((abs difference/relative tolerance)^2 +
    % sum((voxel offset/dta)^2))
    gamma = sqrt(((tar-ref)./(ref*perc/100)).^2 + ...
        (i/dta)^2 + (j/dta)^2 + (k/dta)^2);
else
    % Gamma is defined as the sqrt((abs difference/absolute  tolerance)^2 +
    % sum((voxel offset/dta)^2))
    gamma = sqrt(((tar-ref)/(refval*perc/100)).^2 + ...
        (i/dta)^2 + (j/dta)^2 + (k/dta)^2);
end
