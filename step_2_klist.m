function exitcode = step_2_klist(count_matrix_no_ext, klist)

% This file is part of GPCCA.
%
% Copyright (c) 2018, 2017 Bernhard Reuter
%
% If you use this code or parts of it, cite the following reference:
%
% Reuter, B., Weber, M., Fackeldey, K., Röblitz, S., & Garcia, M. E. (2018). Generalized
% Markov State Modeling Method for Nonequilibrium Biomolecular Dynamics: Exemplified on
% Amyloid β Conformational Dynamics Driven by an Oscillating Electric Field. Journal of
% Chemical Theory and Computation, 14(7), 3579–3594. https://doi.org/10.1021/acs.jctc.8b00079
%
% GPCCA is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU Lesser General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
% -------------------------------------------------------------------------
% main file for clustering of a row-stochastic matrix by GPCCA
% Written by Bernhard Reuter, Theoretical Physics II,
% University of Kassel, 2017
% -------------------------------------------------------------------------

%   set precision with miltiprecision toolbox for certain numerics
%   uncomment if you want to use multiprecision:
%mp.Digits(50) ;

%   uncomment if you want to use multiprecision:
%disp (['number of digits used in multiprecision numerics (mp): ' ...
    %int2str(mp.Digits)])

%   global variable to define precision to use
global class_t ;

%   set precision of variables or expressions wrapped by numeric_t(),
%   i.e. 'double' or 'mp' for multipresicion
class_t = 'double' ;
disp (' ')
disp (['precision to use in sensitive numerics ' ...
    '(i.e. Eigenvalue and Schur decomposition): ' class_t])

% -------------------------------------------------------------------------

%   Parameters for gpcca
%klist                          % list of number of clusters for which to optimize

wk.schur = 1 ;                  % calculate Schurvectors (schur=1) 
                                % or use existing from file (schur=0)
wk.b = 0 ;                      % if b < 0 then -b blocks will be sorted,
                                % if b > 0 then  b or b+1 eigenvalues will 
                                % be sorted, depending on the sizes of 
                                % the blocks,
                                % if b = 0 then the whole Schur form will
                                % be sorted.
wk.init = 1 ;                   % if 1 use A=inv(EVS(index,:)) as starting
                                % guess, 
                                % if =0 read A from file.
wk.solver = 'nelder-mead' ;     % solver for unconstrained optimization 
                                % problem, either 'nelder-mead',
                                % 'levenberg-marquardt', 'gauss-newton'
wk.maxiter = -1 ;
wk.parallel = 0 ;
wk.tolx = 1e-8 ;
wk.tolfun = 1e-8 ;
iopt.init = 2 ;                 % If =1 use A=inv(EVS(index,:)) as starting 
                                % guess, 
                                % if =0 read A from file with identifier 
                                % interactively passed from the command 
                                % window,
                                % if =2 use the the optimized A matrices 
                                % from the  first optimization loop as 
                                % input for the final optimization.
iopt.solver = 'gauss-newton' ;   % solver for optional final optimization
iopt.maxiter = 10 ;
iopt.parallel = 0 ;

% -------------------------------------------------------------------------

%   read the count matrix from file, calculate the stochastic matrix,
%   and call gpcca rotine gpcca.m

%   load the count matrix from file
%disp (' ')
%COUNTMATRIX = input('Enter the name of the matrix file (IN QUOTES): ') ;
%wk.id = input('Enter the id of this simulation (to be used in output file names) (IN QUOTES): ') ;

COUNTMATRIX = strcat(count_matrix_no_ext, '.txt')
wk.id = count_matrix_no_ext

%   perform optimization 
gpcca_step_2_klist(klist, wk, iopt) ;

exitcode = 0

end
