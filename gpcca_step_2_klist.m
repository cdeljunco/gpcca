function [ Pc, chi, A, wk, iopt ] = gpcca_step_2_klist(klist, wk, iopt)
% This file is part of GPCCA.
%
% Copyright (c) 2020, 2018, 2017 Bernhard Reuter
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
% This algorithm generates a fuzzy clustering such that the resulting
% membership functions are as crisp (characteristic) as possible.
%
%
% Input:
%   P                   (N,N)-matrix to be clustered (row-stochastic)
%   sd                  "initial distribution"
%   klist		list of cluster number for which to do the optimization
%   wk                  (possibly empty) structure with the following 
%                       fields,
%                       that are initialized here, if they arent passed.
%	wk.b		Serves to tell SRSchur_num_t.m how many eigenvalues
%			or blocks (and Schurvectors) of the Schurmatrix 
%			(and the Schurvector matrix) shall be sorted:              
% 			If b < 0 then -b blocks will be sorted,
% 			if b > 0 then  b or b+1 eigenvalues will be sorted, 
%			depending on the sizes of the blocks,
% 			if b = 0 then the whole Schur form will be sorted.
%			WARNING: The number of sorted eigenvalues needs to 
%			be larger than the maximal number of clusters!
%			WARNING: If you choose b < 0 you should REALLY
%			know what you are doing!
%       wk.display      If 1, iterative output of the optimization progress
%                       shown. CAUTION: this slows the process
%                       significantly down.
%                       if 0, then not.
%       wk.id           Identification string for naming of the output
%                       files.
%       wk.init         If =1 use A=inv(EVS(index,:)) as starting guess,
%                       if =0 read A from file
%                       with identifier interactively passed from the
%                       command window.
%       wk.maxiter      Maximum number of optimization iterations:
%                       typically 2000-100000 for 'nelder-mead', 
%                       typically ~500 for 'levenberg-marquardt',
%                       typically 100-1000 for 'gauss-newton'.
%                       -1 means off (internal maximum number of iterations
%                       of the algorithm is used).
%                       -2 means to use 50*(k-1)^2 iterations (only in case
%                       'nelder-mead').
%       wk.parallel     Selection of parallel execution of the first 
%                       optimization loop by choosing wk.parallel=1,
%                       else if wk.parallel=0 the program will be executed 
%                       serially (Standard mode).
%       wk.schur        Calculate Schurvectors (schur=1) or use existing
%                       from file (schur=0)
%                       with identifier interactively passed from the 
%                       command window.
%       wk.solver       Solver for unconstrained optimization problem.
%                       One of the following can be chosen:
%                       'melder-mead',
%                       'levenberg-marquardt',
%                       'gauss-newton'.
%       wk.tolfun       Optimization termination tolerance on the function 
%                       value, a positive scalar:
%                       typically 1e-4 for 'nelder-mead',
%                       typically 1e-8 for 'levenberg-marquardt'.
%       wk.tolx         Optimization termination tolerance on x, a positive 
%                       scalar: 
%                       typically 1e-4 for 'nelder-mead',
%                       typically 1e-10 for 'levenberg-marquardt'.
%       wk.xscale       Gauss-Newton optimization factor for the scaling  
%                       vector xscale = xscale * ones(size(xguess)).
%       wk.xtol         Gauss-Newton optimization relative (!) tolerance in 
%                       x (real tolerance related to xscale*xtol).
%   iopt                (possibly empty) structure with the following
%                       optional fields,
%                       that are only necessary, if one wants to perform a
%                       optional final optimization.
%                       If one passes at least iopt.solver, the other
%                       fields are initialized here.
%       iopt.display    If 1, iterative optimization output is shown;
%                       if 0, then not.
%       iopt.init       If =1 use A=inv(EVS(index,:)) as starting guess,
%                       if =0 read A from file
%                       with identifier interactively passed from the
%                       command window,
%                       if =2 use the the optimized A matrices from the 
%                       first optimization loop as input for the optional 
%                       final optimization.
%       iopt.maxiter    Maximum number of optimization iterations:
%                       -1 means off (internal maximum number of iterations
%                       of the algorithm is used),
%                       typically 2000-100000 for 'nelder-mead',
%                       typically ~500 for 'levenberg-marquardt',
%                       typically 100-1000 for 'gauss-newton'.
%       iopt.parallel   Selection of parallel execution of the optional 
%                       optimization loop by choosing iopt.parallel=1,
%                       else if iopt.parallel=0 the program will be  
%                       executed serially (Standard mode).
%       iopt.solver     Solver for unconstrained optimization problem.
%                       One of the following can be chosen:
%                       'melder-mead',
%                       'levenberg-marquardt',
%                       'gauss-newton'.
%       iopt.tolfun     Optimization termination tolerance on the function 
%                       value, a positive scalar:
%                       typically 1e-4 for 'nelder-mead',
%                       typically 1e-8 for 'levenberg-marquardt'.
%       iopt.tolx       Optimization termination tolerance on x, a positive 
%                       scalar:
%                       typically 1e-4 for 'nelder-mead',
%                       typically 1e-10 for 'levenberg-marquardt'.
%       iopt.xscale     Gauss-Newton optimization factor for the scaling 
%                       vector xscale = xscale * ones(size(xguess)).
%       iopt.xtol       Gauss-Newton optimization relative (!) tolerance in 
%                       x (real tolerance related to xscale*xtol).
%
% Output:
%   Pc                  Coarse grained stochastic transition matrix as
%                       returned from the last optimization.
%   chi                 (N,k)-matrix with membership values.
%   A                   (k,k)-tranformation matrix A, chi = evs * A
%   wk                  See above. All values are given back. Might vary
%                       from input values, if those were unfeasible. If no 
%                       input was passed, the here initialized values
%                       are given back.
%   iopt                See above. Values are only given back, if at least
%                       iopt.solver was passed as input.
%                       Might vary from input values, if those were 
%                       unfeasible. If only iopt.solver was passed, the 
%                       here initialized values are given back.
% Cite:
%
% [1] P.Deuflhard, M.Weber: Robust Perron Cluster Analysis in Conformation
%     Dynamics,
%     Lin. Alg. App. 2005, 398c, 161-184.
% [2] S. Roeblitz and M. Weber: Fuzzy Spectral Clustering by PCCA+.
%     In: Classification and Clustering: Models, Software and Applications,
%     WIAS Report No. 26, Berlin 2009.
%     http://www.wias-berlin.de/publications/wias-publ/run.jsp?template=...
%     abstract&type=Report&year=2009&number=26
%
% Written by by Bernhard Reuter, Theoretical Physics II, 
% University of Kassel, 2017,
% based on algorithms and code from Susanna Roeblitz and Marcus Weber, 
% Zuse Institute Berlin, 2012
%--------------------------------------------------------------------------

%   parse command-line klist argument to get double values for cluster numbers 
    %disp(klist)

    if (isa(klist, 'char')) 
        char = strsplit(klist);
        klist = zeros(1, length(char) - 2);
	for i = 2:length(char)-1
	    klist(i-1) = str2double(char(i));
        end
     end

    %disp(klist)    

    kmax = max(klist);
    kmin = min(klist);
%	make sure that one isnt messing with wk.b
    if (wk.b > 0)
        assert(wk.b >= kmax, 'gpcca:k_InputError', ['wk.b !>= kmax: ', ...
            'The number of sorted eigenvalues needs to be larger ', ...
            'than the maximal number of clusters!']) ;
    elseif (wk.b < 0)
        disp(['You chose wk.b < 0: You should REALLY know what you ', ...
            'are doing!'])
    end
%--------------------------------------------------------------------------
    %  initialize the workspace
    [wk, fileid] = initialize_workspace(wk)
    assert(~(wk.init==0 & wk.parallel==1), ...
        'gpcca:IncompatibleOptionsError', ...
        ['You choose parallel execution of the first optimization loop ' ...
        'combined with initialization of A from file. These two ' ...
        'options are mutually exclusive!']) ;
   
%   load the stochastic matrix P and the initial distribution sd
    name = strcat(fileid,'-','P','.txt') ;
    P = load_t(name, '-ascii', 'double');
    name=strcat(fileid,'-','sd','.txt') ;
    sd = load_t(name, '-ascii', 'double');

    class_t1 = class(P) ;
    class_t2 = class(sd) ;
    if (strcmpi(class_t1,class_t2))
        class_t = class_t1 ;
    else
        error('gpcca:P_sd_DataTypeError', ...
            ['class(P) is not equal class(sd)! This will lead to ', ...
            'numeric precision errors!'])
    end

    function r = num_t(expression)
        if (nargin > 0)
            if(strcmpi(class_t,'mp')), r = mp(expression) ;
            else
                if isnumeric(expression)
                    r = expression ;
                else
                    r = eval(expression) ;
                end ;
            end ;
        else
            r = class_t ;
        end ;
    end  % num_t
%   -----------------------------------------------------------------------
%       Global variable 'figureVisible' defines, if figures are visible.
%       Even if figures are invisible (figureVisible=false), the figures
%       are still generated and printed to files.
    global figureVisible
    if figureVisible == false
        % all subsequent figures invisible
        set(0,'DefaultFigureVisible','off');  
    end
%   -----------------------------------------------------------------------
    
    %       initialize recording to diary
    name = strcat(fileid,'-','diary','.txt') ;
    diary(name) ;
   
%   -----------------------------------------------------------------------
%   If you want to use gpcca on eigen- instead of Schur-vectors, replace the
%   calculation and sd-orthonormalization of the Schur-vectors by the
%   calculation and sd-orthonormalization of the eigenvectors. The
%   sd-orthonormalization is imperative for the correct function of the 
%   following implementation of GPCCA!
     	%Read in Schur vectors, etc.
	inputname = strcat(fileid, '-X.txt')
     	X = load_t(inputname,'-ascii',num_t) ;
        dummy = (X'*diag(sd)*X - num_t(eye(size(X,2))) > (num_t('10000')*eps(num_t))) ;
        assert(~any(dummy(:)),'gpcca:X_OrthogonalityError', ...
            'Schurvectors appear to not be orthogonal')
        clearvars dummy
        dummy = ( abs(X(:,1) - 1) < ( numeric_t('1000') * eps(numeric_t) ) ) ;
        assert(all(dummy(:)), 'gpcca:FirstColumnError', ...
            'X(:,1) isnt equal 1!')
       
        inputname = strcat(fileid, '-RR.txt')
	RR = load_t(inputname,'-ascii','double') ;
%           check RR for correctness
        assert(size(RR,1)==size(RR,2),'gpcca:MatrixShapeError1', ...
        'RR Matrix is not quadratic but %d x %d',size(RR,1),size(RR,2))
%           check if the number of rows of X matches with the shape of RR
        assert( size(X,1)==size(RR,1), 'gpcca:MatrixShapeError2', ...
        ['The number of rows size(X,1)=%d of X doesnt match with the ', ...
        'shape (%d,%d) of RR!'],...
        size(X,1), size(RR,1), size(RR,2) )

    [ badblocks, badblockrows ] = find_twoblocks( double(X), double(RR), fileid ) ;
    clearvars RR
    disp(' ')
    disp('The following cluster numbers will be excluded in all further analysis ')
    disp('since using them would lead to splitting of 2x2-blocks of eigenvalues!')
    disp(badblockrows)
    
%   -----------------------------------------------------------------------

    EVS = X ;

%       convert to double for futher calculations (but not EVS!) 
    sd = double(sd) ;
    %pi = double(pi) ;
    P = double(P) ;


%       decide for "best" number of clusters (other criteria are possible):
%       the optimal value for trace(S) can be at most kt
%       therefore crispness=trace(S)/kt<1, but should be as large as 
%       possible

    disp (' ')
    disp ('Optimize crispness vs. number of clusters')
    disp ('=============================================')
    disp (' ')

    assert(isa(EVS,num_t),'gpcca:EVS_DataTypeError', ...
        'Variable is type %s not %s',class(EVS),num_t)

    optimization_ks = length(klist);
%       initialize A
    disp(klist)
    disp('max(klist)')
    disp(max(klist))
    A_cell = cell([1,max(klist)]) ;
    for kt = klist
        if badblocks(kt) == false
            evs = EVS(:,1:kt) ;
            A_cell{kt} = initialize_A(evs, wk.init) ;
	else
	    optimization_ks = optimization_ks - 1;
        end
    end
    
    if optimization_ks > 0
	    [ Pc, A_cell, chi, val_vec, opt_vec ] = optimization_loop_klist( P, sd, ...
	        double(EVS), A_cell, klist, wk, fileid ) ;
	    A = A_cell{kmax} ;
	    
	%   -----------------------------------------------------------------------
	
	%       save the 'val'-vector and the crispness-vector
	    name=strcat(fileid,'-','opt_vec','-','n=',int2str(kmin),'-', ...
	        int2str(kmax),'.txt') ;
	    save_t(name,opt_vec,'-ascii')
	    name=strcat(fileid,'-','val_vec','-','n=',int2str(kmin),'-', ...
	        int2str(kmax),'.txt') ;
	    save_t(name,val_vec,'-ascii')
	
	%       plot the crispness over the number of clusters    
	    fig5 = figure ;
	    plot(klist, opt_vec(:,2),'-x')
	    xlabel('number of clusters')
	    ylabel('crispness')
	    name=strcat(fileid,'-','crispness-figure','-','n=',int2str(kmin), ...
	        '-',int2str(kmax)) ;
	    h = gcf ; % Current figure handle
	    set(h,'PaperPositionMode','manual') ;
	    set(h,'PaperUnits','inches') ;
	    set(h,'PaperPosition',[0 0 3.33 2.63]) ;
	    set(h,'PaperUnits','inches') ;
	    set(h,'PaperSize',[3.33 2.63]) ;
	    savefig(fig5,strcat(name,'.fig'),'compact')
	    print(fig5,strcat(name,'.pdf'),'-dpdf')
	    print(fig5,strcat(name,'.png'),'-dpng','-r600')
	 
	%   -----------------------------------------------------------------------
	
	%       initialize optional parameters for optional final optimization 
	    [iopt, fileid, final_opt] = initialize_optional(iopt, wk.id)
	    
	%       test, if final optimization shall be performed
	    if final_opt == true
	        
	        assert(~(iopt.init==0 & iopt.parallel==1), ...
	        'gpcca:IncompatibleOptionsError', ...
	        ['You choose parallel execution of the optional optimization loop ' ...
	        'combined with initialization of A from file. These two ' ...
	        'options are mutually exclusive!']) ;
	        
	%           change the diary
	        name = strcat(fileid,'-','diary','.txt') ;
	        diary(name) ;
	    
	%       -------------------------------------------------------------------
	%           optional final optimization for beforehand determined optimal
	%           cluster number
	        
	        disp (' ')
	        disp ('Perform final optimization')
	        disp ('==========================')
	        disp (' ')
	    
	        assert(isa(EVS,num_t),'gpcca:EVS_DataTypeError', ...
	            'Variable is type %s not %s', class(EVS), num_t)
	       
	%           initialize A_cell, if iopt.init~=2, else use A_cell from the
	%           previous optimization
	        if ( iopt.init == 0 || iopt.init == 1 )  
	            A_cell = cell([1,kmax]) ;
	            for kt = klist
	                if badblocks(kt) == false
	                    evs = EVS(:,1:kt) ;
	                    A_cell{kt} = initialize_A(evs, wk.init) ;
	                end
	            end
	        elseif iopt.init == 2
	        else 
	            error('gpcca:IoptInitError',['Couldnt '...
	            'initialize A for the final optimization since iopt.init ' ...
	            'was neither =0, =1 nor =2!'])
	        end
	        
	        [ Pc, A_cell, chi, val_vec, opt_vec ] = optimization_loop_klist( P, sd,...
	            double(EVS), A_cell, klist, iopt, fileid ) ;
	        A = A_cell{kmax} ;
	        
	%       -------------------------------------------------------------------
	
	%           plot the crispness over the number of clusters, if more than
	%           one final optimization was performed
	        if optimization_ks > 1
	   
	            fig6 = figure ;
	            plot(klist, opt_vec(:,2),'-x')
	            xlabel('number of clusters')
	            ylabel('crispness')
	            name=strcat(fileid,'-','crispness-figure','-','n=', ...
	                int2str(kmin),'-',int2str(kmax)) ;
	            h = gcf ; % Current figure handle
	            set(h,'PaperPositionMode','manual') ;
	            set(h,'PaperUnits','inches') ;
	            set(h,'PaperPosition',[0 0 3.33 2.63]) ;
	            set(h,'PaperUnits','inches') ;
	            set(h,'PaperSize',[3.33 2.63]) ;
	            savefig(fig6,strcat(name,'.fig'),'compact')
	            print(fig6,strcat(name,'.pdf'),'-dpdf')
	            print(fig6,strcat(name,'.png'),'-dpng','-r600')
	            
	%               save the 'val'-vector and the crispness-vector
	            name=strcat(fileid,'-','opt_vec','-','n=',int2str(kmin),'-', ...
	                int2str(kmax),'.txt') ;
	            save_t(name,opt_vec,'-ascii')
	            name=strcat(fileid,'-','val_vec','-','n=',int2str(kmin),'-', ...
	                int2str(kmax),'.txt') ;
	            save_t(name,val_vec,'-ascii')
	        
	        end
	        
	%       -------------------------------------------------------------------
	
	    end
	
	%   -----------------------------------------------------------------------
	    set(0,'DefaultFigureVisible','on');  % all subsequent figures visible
	    diary OFF ;
    end
end
