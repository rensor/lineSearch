classdef lineSearch
  
  properties(Constant)
    name = 'lineSearch';
  end

  properties
    options = [];
  end
  
  properties(SetAccess = public, Hidden = true)

    % Switch
    initialized = false;
    
    % Inputs
    fun = [];
    x0 = [];
    d = [];
    method = [];
    f0 = [];
    df0 = [];
    
  end
  
  methods
    
    % Construct
    function this = lineSearch(fun,x0,d,method,f0,df0,varargin)
      
      % initialize data structures
      if nargin < 1 || ~isa(fun,'function_handle')
        error([this.name,': fun is required to be a function handle'])
      else
        this.fun = fun;
      end
      
      if nargin < 2 || isempty(x0)
        this.x0 = [];
      else
        if isnumeric(x0)
          this.x0 = x0(:);
        else
          error([this.name,': x0 is required to be of type numeric'])
        end
      end
      
      if nargin < 3 || isempty(d)
         error([this.name,': direction vector d is required'])
      else
        if isnumeric(d)
          this.d = d(:);
        else
          error([this.name,' direction is required to be of type numeric'])
        end
      end
      
     if nargin < 4 || isempty(method)
         error([this.name,': Please specify which line search method to apply'])
      else
        if ischar(method)
          this.method = method;
        else
          error([this.name,': Select method by string input'])
        end
      end
      
      if nargin < 5 || isempty(f0)
         this.f0 = [];
      else
        if isnumeric(f0)
          this.f0 = f0;
        else
          error([this.name,': Initial function value, f0, is required to be of type numeric'])
        end
      end
      
      if nargin < 6 || isempty(df0)
         this.df0 = [];
      else
        if isnumeric(df0)
          this.df0 = df0(:);
        else
          error([this.name,': Initial gradient value, df0, is required to be of type numeric'])
        end
      end
      % initialize options structure
      this.options = lineSearch.setOptions(varargin);
      
      % We made it this far
      this.initialized = true;
    end
    
    %% LineSearch algorithm
    function [xout,fval,stepSize,relStepSize,nFeval,nGradEval,exitflag,message] = solve(this)
      
      if ~this.initialized
        error([this.name,': Please run construct to initialize data structures'])
      end
      
      % Initialize variables
      nGradEval = 0;
      alpha = inf;
      fval = inf;
      nf1 = 0;
      nf2 = 0;
      xout = [];
      % calculate norm of step
      dNorm = norm(this.d);
      dUnit = this.d./dNorm;
          
      if ~isempty(this.f0) && ~isempty(this.df0)
        dphi0 = this.df0'*dUnit;
        phi0 = this.f0;
      else
        dphi0 = [];
        phi0 = [];
      end
      
      % Define linesearch function
      phiFun = @(alpha) lineSearchObj(this.fun,this.x0,dUnit,alpha);
      
      switch this.method
        case 'none'
          alpha = 1;
          exitflag = 1;
          message = 'lineSearch: Unit step applied';
          xout = this.x0 + this.d; % New point
          
        case 'wolfe'
          % bracket linesearch function
          alpha0 = 0;
          stepIncrement = max(dNorm*1e-4,this.options.StepTolerance/2);
          [alphaU, alphaL, ~, ~, nf1, exitflag,message] = lineSearch.bracket(phiFun,alpha0,phi0,stepIncrement,...
            'Display',this.options.Display,'MaxFunctionEvaluations',this.options.MaxFunctionEvaluations);
          % If the bracket function worked
          if exitflag > 0
            [alpha,fval,nf2,nGradEval,exitflag,message] = lineSearch.wolfe(phiFun,alphaL,alphaU,phi0,dphi0,...
              'c1',this.options.c1,'c2',this.options.c2,'Display',this.options.Display,'MaxFunctionEvaluations',this.options.MaxFunctionEvaluations,'StepTolerance',this.options.StepTolerance,'FunctionTolerance',this.options.FunctionTolerance);
            if exitflag > 0
              xout = this.x0 + dUnit*alpha; % New point
            end
          end
          
        case 'golden'
          % bracket linesearch function
          alpha0 = 0;
          stepIncrement = max(dNorm*1e-4,this.options.StepTolerance/2);
          [alphaU, alphaL, ~, ~, nf1, exitflag,message] = lineSearch.bracket(phiFun,alpha0,phi0,stepIncrement,...
            'Display',this.options.Display,'MaxFunctionEvaluations',this.options.MaxFunctionEvaluations);
          % If the bracket function worked
          if exitflag > 0
            [alpha,fval,nf2,exitflag,message] = lineSearch.goldenSection(phiFun,alphaL,alphaU,...
              'Display',this.options.Display,'MaxFunctionEvaluations',this.options.MaxFunctionEvaluations,'StepTolerance',this.options.StepTolerance,'FunctionTolerance',this.options.FunctionTolerance);
            if exitflag > 0
              xout = this.x0 + dUnit*alpha; % New point
            end
          end
          
        case 'backtrack'
          % Set upper bound for alpha
          alphaU = dNorm*2;
          [alpha,fval,nf2,nGradEval,exitflag,message] = lineSearch.backtracking(phiFun,alphaU,phi0,dphi0,...
             'Display',this.options.Display,'MaxFunctionEvaluations',this.options.MaxFunctionEvaluations,'StepTolerance',this.options.StepTolerance);
           if exitflag > 0
            xout = this.x0 + dUnit*alpha; % New point
           end
        otherwise
          exitflag = -3;
          message = sprintf('lineSearch: Unknown algorithm, %s',this.method);
      end
      
      nFeval = nf1 + nf2;
      relStepSize = alpha;
      stepSize = dNorm*relStepSize;
      
      function [phi,dphi] = lineSearchObj(fun,x,d,alpha)
        
        if nargout > 1
          [phi, df] = fun(x+d*alpha);
          dphi = df'*d;
        else
          phi = fun(x+d*alpha);
        end
      end
      
    end
  end
  
  methods(Static = true)
  
    function [alphaU, alphaL, phiU, phiL, nFeval, exitflag, message, dsign] = bracket(phi,alpha0,f0,delta,varargin)
      % Help for bracket
      % This function brackets a local minimum for a 1D line search.
      % Problem is formulated as: phi(alpha) = fun(x+alpha*d), where the goal is to
      % determine upper and lower bounds of alpha which bounds a local minimum
      % of the function fun. 
      % The output is the upper and lower bounds of the minimum location.
      %
      % Inputs:
      %         phi:    Function handle for the 1D line search, scalar output. phi(alpha) = fun(x+alpha*d)
      % Optional inputs: 
      %         alpha0:  default = 0. Initial lower bound for the
      %                 minimum location.
      %         f0:   phi(alpha0). function value at input alpha
      %         delta:  default = 1e-3, scalar > 0. Controls the
      %                 step size used for bracketing. alpha = alpha + delta*(1.618)^IterNo;
      %
      % Parsed inputs: 
      %         Display: string, default ='off', set to 'on' for iteration print-out to commandwindow
      %         MaxFunctionEvaluations: integer > 0, maximum allowed function evaluations before termination of algorithm. 
      %         StepTolerance: scalar > 0, algorithm terminates if "distance" between lower and upper bound is below threashold. 
      %         directionSwitch, logical, default=false, if set to true, the sign of the direction/decent vector can be changed if the algorithm sees divergence.
      %
      % Outputs:
      %         alphaU  : Upperbound on alpha
      %         alphaL  : Lowerbound on alpha
      %         phiU    : phi(alphaU)
      %         phiL    : phi(alphaL)
      %         nFeval  : Number of used function evaluations
      %         exitflag: 
      %           0 (failed to bracket, something unforseen went wrong)
      %           1 (sucessfully bracket)
      %           2 (direction switch used, sucessfully bracket)
      %          -1 (maximum number of function evaluations reached)
      %          -2 (direction switch used, failed)
      %         message: output message from algorithm
      %         dsign: output sign multiplied on the input decent/direction vector used for bracket
      
      p = inputParser;
      p.CaseSensitive = false;
      % Helper functions for input parser
      checkEmpetyOrChar = @(x) (isempty(x) || ischar(x));
      checkEmptyOrNumericPositive = @(x) (isempty(x) || (isnumeric(x) && all(x > 0)));
      % Set parameters;
      p.addParameter('Display','off',  @(x) checkEmpetyOrChar(x));
      p.addParameter('MaxFunctionEvaluations',1000,  @(x) checkEmptyOrNumericPositive(x));
      p.addParameter('directionSwitch',false,  @(x) islogical(x));

      parse(p,varargin{:});
      options = p.Results;
      
      % Initialize function evaluation counter 
      nFeval = 0;
      
      if nargin < 1 || isempty(phi)
        error('bracket: Missing input \n Please supply linesearch function \n The input is either missing or provided empty')
      end
      
      if nargin < 2 || isempty(alpha0)
        alpha0 = 0;
      end
      
      if nargin < 3 || isempty(f0)
        f0 = phi(alpha0);
        nFeval = nFeval + 1;
      end
      
      if nargin < 4 || isempty(delta)
        delta = 1e-3;
      else
        if delta < 0
          error('bracket: Missing input \n Incorrect input \n Please supply delta > 0 \n The supplied delta must be larger than zero')
        end
      end
      
      alphaU = [];
      alphaL = [];
      phiU = [];
      phiL = [];
      
      % If decent direction switch is allowed, we initialize a parameter
      % telling us that this is the first time we switch direction.
      dSwitchUsed = false;
            
      % Initialize step counter
      kk = 0;
      % Store initial values
      deltaInit = delta;
      % Initialize first "previous" value
      prevPhi = f0;
      % Assume this will fail
      exitflag = 0;
      % Initial sign multiplied on the decent direction. No change
      dsign = 1;
      % First alpha value
      alpha = alpha0;
      % Run algorithm
      while true
        alpha = alpha + delta*(1.618)^kk;
        [curPhi] = phi(alpha*dsign);
        nFeval = nFeval + 1;
        if prevPhi <= curPhi && kk > 0
          % Success, problem has been bracked/bound
          alphaU = alpha;
          if kk > 1
            alphaL = alpha - delta*(1.618)^kk - delta*(1.618)^(kk-1);
          else
            alphaL = alpha;
            prevPrevPhi = f0;
          end
          
          % Set exitflag
          if options.directionSwitch
            if dSwitchUsed
              exitflag = 2;
              message = sprintf('MathUtilities.bracket: \n Direction switch used, Sucessfully bracket');
            else
              exitflag = 1;
              message = sprintf('MathUtilities.bracket: \n Sucessfully bracket');
            end
          else
            exitflag = 1;
            message = sprintf('MathUtilities.bracket: \n Sucessfully bracket');
          end
          % Exit while loop
          break
        else
          if kk == 0 && (curPhi > prevPhi)
             % If here, we failed in the first attempt to find an alpha value
             % which could minimize the line search function. We now have two
             % options. The first is to reduce the initial delta value. If this
             % fails, we check if the user has allowed for a decent
             % direction switch. If so, we try to reduce the delta value again.
             % If the user has not given this option, the algorithm terminates.
             
             % initialize "second" trial counter
            jj = 0;
            
            % Check the input alpha values.
            if alpha0 > 0
              trialAlpha = alpha0;
            else
              trialAlpha = 1;
            end
            % Try to reduce delta to see if we are able to minimize the linear
            % search function.
            while true
              delta = delta / 10;
              [curPhi] = phi(delta*trialAlpha*dsign);
              jj = jj + 1;
              nFeval = nFeval + 1;
              if curPhi <= f0
                prevPhi = curPhi;
                alpha = alpha0;
                break
              end
              % Sanity check, have we tried this too many times 
              if jj == 10
                if options.directionSwitch
                  if ~dSwitchUsed
                    delta = deltaInit;
                    dSwitchUsed = true;
                    dsign = -1;
                    jj = 0;
                  else
                    exitflag = -2;
                    message = sprintf('MathUtilities.bracket: \n Applied decent direction seems wrong, tried to switch direction, did not help');
                    return
                  end
                else
                  exitflag = -1;
                  message = sprintf('MathUtilities.bracket: \n Applied decent direction seems wrong, \n check supplied decent direction, initial alpha, and delta values');
                  return
                end
              end
                
              if nFeval >= options.MaxFunctionEvaluations 
                exitflag = -1;
                message = sprintf('MathUtilities.bracket: \n unable to determine delta parameter, this indicates the either the applied decent direction is wrong, \ n or the supplied function (phi) is faulty');
                break
              end
            end
          else
            % Everything looks good, update previous values
            prevPrevPhi = prevPhi;
            prevPhi = curPhi;
            kk = kk + 1;
          end
        end
        if nFeval >= options.MaxFunctionEvaluations 
          message = sprintf('MathUtilities.bracket: \n Maximum number of allowed function iterations reached: %i >= %i',nFeval,options.MaxFunctionEvaluations );
          exitflag = -1;
          break
        end
      end
      
      phiU = curPhi;
      phiL = prevPrevPhi;
        
      if dSwitchUsed
        phiU = prevPrevPhi;
        phiL = curPhi;
        temp = alphaL;
        alphaL = alphaU;
        alphaU = temp;
      end
      
    end
    
    function [alpha,fval,nFeval,nGrad,exitflag,message]= wolfe(phi,alphaL,alphaU,phi0,dphi0,varargin)
      % Help for wolfe
      % This function determine a local minimum for a 1D linesearch problem.
      % Problem is formulated as: phi(alpha) = fun(x+alpha*d), where the goal is to
      % determine optimum value of alpha which satisfies the Strong Wolfe Conditions
      % Based on the text book: Numerical Optimization by Nocedal and Wright, second edition, page 60, algorithm 3.5
      %
      % The output is the optimum step size (alpha) which satisfies the Strong Wolfe Conditions
      %
      % Inputs:
      %         phi     : Function handle for the 1D line search, scalar output. phi(alpha) = fun(x+alpha*d)
      %         alphaL  : Lower bound of alpha
      %         alphaU  : Upper bound of alpha
      % Optional inputs: 
      %         f0    : phi(0). function value at input alpha
      %         dphi0   : Derivative of phi at alpha=0 i.e., dphi/dalpha|alpha=0
      %
      % Parsed inputs: 
      %         Display: string, (default ='off'), set to 'on' for iteration print-out to commandwindow
      %         MaxIterations: integer > 0,(default=1000), maximum allowed iteration by the algorithm
      %         MaxFunctionEvaluations: integer > 0,(default=1000), maximum allowed function evaluations before termination of algorithm. 
      %         StepTolerance: scalar > 0,(default=1e-10), algorithm terminates if "distance" between lower and upper bound is below threashold. 
      %         FunctionTolerance: scalar > 1 (default=1e-6), algorithm terminates if the difference in objective value between two sucessive iterations is below threashold%         FunctionTolerance: scalar > 1 (default=1e-6), algorithm terminates if the difference in objective value between two sucessive iterations is below threashold
      %         c1: scalar < 1, (default=1e-4), parameter used in the strong wolfe convergence criteria
      %         c2: scalar > c1 < 1 , (default=0.1), parameter used in the strong wolfe convergence criteria
      %
      % Outputs:
      %         alpha   : Optimum step size which satisfies the Strong Wolfe Conditions
      %         fval    : phi(alpha)
      %         nFeval  :     Number of used function evaluations
      %         nGrad   :     Number of used gradient evaluations
      %         exitflag: 
      %           0 (failed, something unforseen went wrong)
      %           1 (sucessfully determine optimum step size)
      %          -1 (maximum number of function evaluations reached)
      %          -2 (step size reached upper bound, check bounds)
      %         message: output message from algorithm
      
      p = inputParser;
      p.CaseSensitive = false;
      % Helper functions for input parser
      checkEmpetyOrChar = @(x) (isempty(x) || ischar(x));
      checkEmptyOrNumericPositive = @(x) (isempty(x) || (isnumeric(x) && all(x > 0)));
      % Set parameters;
      p.addParameter('Display','off',  @(x) checkEmpetyOrChar(x));
      p.addParameter('MaxIterations',1000,  @(x) checkEmptyOrNumericPositive(x));
      p.addParameter('MaxFunctionEvaluations',1000,  @(x) checkEmptyOrNumericPositive(x));
      p.addParameter('StepTolerance',1e-10,  @(x) checkEmptyOrNumericPositive(x));
      p.addParameter('FunctionTolerance',1e-6,  @(x) checkEmptyOrNumericPositive(x));
      p.addParameter('c1',1e-4,  @(x) checkEmptyOrNumericPositive(x));
      p.addParameter('c2',0.1,  @(x) checkEmptyOrNumericPositive(x));
      parse(p,varargin{:});
      
      options = p.Results;
      
      c1 = options.c1;
      c2 = options.c2;
      nFeval = 1;
      nGrad = 1;
      exitflag = 0;
      
      if isempty(dphi0) || isempty(f0)
        [phi0,dphi0] = phi(0);
      end

      % The algorim starts with a low value of alpha, and increases it until
      % convergence
      alpha = alphaL;
      % alpha > alpha_m1
      
      % Initialize previous values
      alpha_m1 = 0;
      phi_m1 = phi0;
      dphi_m1 = dphi0;
      
      ii = 0;
      while true
        ii = ii + 1;
        nFeval = nFeval + 1;
        nGrad = nGrad + 1;
        [curPhi,curDphi] = phi(alpha);
        if (curPhi > phi0 + c1*alpha*dphi0) || (ii > 1) && (curPhi >= phi_m1)
          [alpha,fval,nFevalZ,nGradZ,exitflag,message] = wolfeZoom(phi,alpha_m1,alpha,phi_m1,curPhi,phi0,dphi0);
          nFeval = nFeval + nFevalZ;
          nGrad = nGrad + nGradZ;
          return
        else
          if abs(curDphi) <= -c2*dphi0
            fval = curDphi;
            exitflag = 1;
            message = sprintf('wolfe: Sucessfully determine optimum step size based on gradient requirement.\n abs(dphi/dalpha) <= -c2*dphi0 \n %0.5e <= %0.5e',abs(curDphi),-c2*dphi0);
            return
          elseif ii > 1 && curDphi >=0
            [alpha,fval,nFevalZ,nGradZ,exitflag,message] = wolfeZoom(phi,alpha,alpha_m1,curPhi,phi_m1,phi0,dphi0);
            nFeval = nFeval + nFevalZ;
            nGrad = nGrad + nGradZ;
            return
          end
        end
        
        if ii >= options.MaxIterations
          message = sprintf('wolfe: Maximum number of allowed iterations reached: %i >= i%', ii,options.MaxIterations);
          exitflag = -1;
          alpha = inf;
          fval = inf;
          return
        elseif alpha <= options.StepTolerance
          message = sprintf('wolfe: Step size tolerance reached, alpha=%0.5e <= %0.5e',alpha,options.StepTolerance);
          exitflag = 1;
          nFeval = nFeval + 1;
          fval = phi(alpha);
          return
        elseif alpha == alphaU
          exitflag = -2;
          message = sprintf('wolfe: Step size reached upper bound, check bounds.\n alpha=alphaU \n %0.5e=%0.5e',alpha,alphaU);
          nFeval = nFeval + 1;
          [fval] = phi(alpha);
          return
        elseif abs(phi_m1-curPhi) <= options.FunctionTolerance
          message = sprintf('wolfe: Function Tolerance reached, f(alpha_(i-1))-f(alpha_i)%0.5e <= %0.5e',abs(phi_m1-curPhi),options.FunctionTolerance);
          exitflag = 1;
          fval = curPhi;
          return
        end
        
        % update previous values
        alpha_m1 = alpha;
        phi_m1 = curPhi;
        % Set new alpha
        alpha = alpha+(alphaU-alpha)*0.1;
             
      end
      
      function [alpha, fval, nFeval, nGrad, exitflag, message] = wolfeZoom(phi,alphaHi,alphaLo,phiHi,phiLo,phi0,dphi0)
        % Helper function for the wolfe algorithm. 
        % Based on the text book: Numerical Optimization by Nocedal and Wright, second edition, page 61, algorithm 3.6
        jj = 0;
        nFeval = 0;
        nGrad = 0;
        alpha = 0;
        fval = inf;
        exitflag = 0;
        message = '';
        while true
          L = (alphaHi-alphaLo);
          alphaj = L/2+alphaLo;
          phij = phi(alphaj);
          nFeval = nFeval + 1;
          jj = jj + 1;
          if phij > phi0+c1*alphaj*dphi0 || phij >= phiLo
            alphaHi = alphaj;
            phiHi = phij;
          else
            [phij,dphij] = phi(alphaj);
            nFeval = nFeval + 1;
            nGrad = nGrad + 1;
            if abs(dphij) <= -c2*dphi0
              alpha = alphaj;
              exitflag = 1;
              message = sprintf('wolfe.wolfeZoom: Sucessfully determine optimum step size based on gradient requirement.\n abs(dphi/dalpha) <= -c2*dphi0 \n %0.5e <= %0.5e',abs(dphij),-c2*dphi0);
              return
            elseif dphij*(phiHi-phiLo) >= 0
              alphaHi = alphaLo;
              phiHi = phiLo;
            end
            alphaLo = alphaj;
            phiLo = phij;
          end
          
          if jj >= options.MaxIterations 
            exitflag = -1;
            message = sprintf('wolfe.wolfeZoom: Exceeded maximum number of iterations %i',options.MaxIterations);
            alpha = alphaj;
            fval = phij;
            return
          elseif alphaj <= options.StepTolerance
            exitflag = 1;
            message = sprintf('wolfe.wolfeZoom: Step size tolerance reached, alpha=%0.5e <= %0.5e',alphaj,options.StepTolerance);
            alpha = alphaj;
            fval = phij;
            return
          elseif abs(phiHi-phiLo) <= options.FunctionTolerance 
            message = sprintf('wolfe.wolfeZoom: Function Tolerance reached, f(alphaU))-f(alphaL)%0.5e <= %0.5e',abs(phiHi-phiLo),options.FunctionTolerance);
            exitflag = 1;
            fval = curPhi;
            return
          end
        end
      end
    end
    
    function [alpha,phiAlpha,nFeval,nGrad,exitflag,message] = backtracking(phi,alphaU,f0,dphi0,varargin)
      p = inputParser;
      p.CaseSensitive = false;
      % Helper functions for input parser
      checkEmpetyOrChar = @(x) (isempty(x) || ischar(x));
      checkEmptyOrNumericPositive = @(x) (isempty(x) || (isnumeric(x) && all(x > 0)));
      % Set parameters;
      p.addParameter('Display','off',  @(x) checkEmpetyOrChar(x));
      p.addParameter('MaxIterations',1000,  @(x) checkEmptyOrNumericPositive(x));
      p.addParameter('MaxFunctionEvaluations',1000,  @(x) checkEmptyOrNumericPositive(x));
      p.addParameter('StepTolerance',1e-10,  @(x) checkEmptyOrNumericPositive(x));
      p.addParameter('c1',1e-4,  @(x) checkEmptyOrNumericPositive(x));
      parse(p,varargin{:});
      
      options = p.Results;
      
      nFeval = 0;
      nGrad = 0;
      if isempty(dphi0) || isempty(f0)
        [f0,dphi0] = phi(0);
        nFeval = nFeval + 1;
        nGrad = nGrad + 1;
      end
      
      alpha = alphaU;
      while true
        curPhi = phi(alpha);
        nFeval = nFeval + 1;
        if curPhi <= f0 + options.c1*alpha*dphi0
          exitflag = 1;
          message = sprintf('backtracking: Sucessfully determine step size based on the Armijo-Goldstein condition');
          phiAlpha = curPhi;
          return
        elseif nFeval >= options.MaxFunctionEvaluations
          message = sprintf('backtracking: Maximum number of allowed function evaluations reached: %i >= %i', nFeval,options.MaxFunctionEvaluations);
          exitflag = -1;
          alpha = alpha;
          phiAlpha = curPhi;
          return
        elseif alpha <= options.StepTolerance
          exitflag = 1;
          message = sprintf('backtracking: Step size tolerance reached, alpha=%0.5e <= %0.5e',alpha,options.StepTolerance);
          phiAlpha = curPhi;
          return
        end
        alpha = alpha*0.5;
      end
    end
    
    function [alpha,phiAlpha,nFeval,exitflag,message] = goldenSection(phi,alphaL,alphaU,varargin)
      p = inputParser;
      p.CaseSensitive = false;
      % Helper functions for input parser
      checkEmpetyOrChar = @(x) (isempty(x) || ischar(x));
      checkEmptyOrNumericPositive = @(x) (isempty(x) || (isnumeric(x) && all(x > 0)));
      % Set parameters;
      p.addParameter('Display','off',  @(x) checkEmpetyOrChar(x));
      p.addParameter('MaxFunctionEvaluations',1000,  @(x) checkEmptyOrNumericPositive(x));
      p.addParameter('StepTolerance',1e-10,  @(x) checkEmptyOrNumericPositive(x));
      p.addParameter('FunctionTolerance',1e-6,  @(x) checkEmptyOrNumericPositive(x));
      parse(p,varargin{:});
      options = p.Results;
      
      % Define "golden" constants
      gC = (sqrt(5)+1)/2;
      invGc = 1/gC;
      invGc2 = 1/gC^2;
      
      % Set function evaluation counter
      nFeval = 0;
      funcTol = false;
      % Interval reduction
      h = alphaU-alphaL;
      
      if h <= options.StepTolerance
        alpha = (alphaL+alphaU)/2;
        nFeval = nFeval + 1;
        [phiAlpha] = phi(alpha);
        exitflag = 1;
        message = sprintf('goldenSection: Initial upper and lower bounds within step size tolerance, alpha=%0.5e <= %0.5e',alpha,options.StepTolerance);
        return
      end
      
      % required steps to reach tolerance
      n = ceil(log(options.StepTolerance/h)/log(invGc));
      
      % Define inner interval
      c = alphaL + invGc2*h;
      d = alphaL + invGc*h;
      % Evaluate point c
      nFeval = nFeval + 1;
      yc = phi(c);
      % Evaluate point d
      nFeval = nFeval + 1;
      yd = phi(d);
      
      for ii = 1:n
        if yc < yd
          alphaU = d;
          d = c;
          yd = yc;
          h = invGc*h;
          c = alphaL + invGc2*h;
          % Evaluate point c
          nFeval = nFeval + 1;
          yc = phi(c);
        else
          alphaL = c;
          c = d;
          yc = yd;
          h = invGc*h; 
          d = alphaL + invGc*h;
          % Evaluate point d
          nFeval = nFeval + 1;
          yd = phi(d);
        end
        if abs(yc-yd) <= options.FunctionTolerance
          funcTol = true;
          break
        end
      end
      
      if yc < yd
        c = alphaL;
        alpha = (c+d)/2;
      else
        d = alphaU;
        alpha = (c+d)/2;
      end
      
      % Evaluate final point
      nFeval = nFeval + 1;
      [phiAlpha] = phi(alpha);
      if funcTol
        exitflag = 1;
        message = sprintf('goldenSection: Function tolerance reached, f(alphaU)-f(alphaL) %0.3e <= %0.3e',abs(yc-yd),options.FunctionTolerance);
      else
        exitflag = 1;
        message = sprintf('goldenSection: Step size tolerance reached, (alphaU-alphaL) = %0.3e <= %0.3e',alphaU-alphaL,options.StepTolerance);
      end
    end
    
  end
  
  methods(Static = true, Hidden = true)
    % options initialization
    function options = setOptions(input)
      % Here you can add new options if needed
      p = inputParser;
      p.CaseSensitive = false;
      % Helper functions for input parser
      checkEmpetyOrChar = @(x) (isempty(x) || ischar(x));
      checkEmptyOrNumericPositive = @(x) (isempty(x) || (isnumeric(x) && all(x > 0)));

      % Set parameters;
      p.addParameter('Display','off',  @(x) checkEmpetyOrChar(x));
      p.addParameter('MaxGradientEvaluations',1000,  @(x) checkEmptyOrNumericPositive(x));
      p.addParameter('MaxFunctionEvaluations',1000,  @(x) checkEmptyOrNumericPositive(x));
      p.addParameter('StepTolerance',1e-10,  @(x) checkEmptyOrNumericPositive(x));
      p.addParameter('FunctionTolerance',1e-6,  @(x) checkEmptyOrNumericPositive(x));
      p.addParameter('c1',1e-4,  @(x) checkEmptyOrNumericPositive(x));
      p.addParameter('c2',0.1,  @(x) checkEmptyOrNumericPositive(x));
      
      % pars input
      if nargin < 1 || isempty(input)
        parse(p);
      else
        parse(p,input{:});
      end
      
      % Output results to options structure
      options = p.Results;
      
    end
  end
end
