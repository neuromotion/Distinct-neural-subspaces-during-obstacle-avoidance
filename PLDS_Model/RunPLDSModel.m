function [TrainTraj,TestTraj,LL]=RunPLDSModel(TrainSpikeCount,TestSpikeCount,Dim,MaxIter,MaxTime)
    % [TrainTraj,TestTraj,LL]=RunPLDSModel(TrainSpikeCount,TestSpikeCount,Dim,MaxIter,MaxTime)    
    
    % This function is a concise adaptation of the general code provided by 
    % Macke JH, Buesing L, Sahani M (2015) 
    % Estimating state and parameters in state space models of spike trains
    % in Advanced state space methods for neural and clinical data, volume 137
    % (Chen Z, editor), Cambridge University Press, Cambridge, UK   
        
    seq=[];
    for k=1:length(TrainSpikeCount)
        seq(k).y=TrainSpikeCount{k};
        seq(k).T=size(TrainSpikeCount{k},2);
    end
    params=PLDSInitialize(seq,Dim,'NucNormMin',[]);
    %params=PLDSInitialize(seq,Dim,'ExpFamPCA',[]);    
    params.model.inferenceHandle=@PLDSLaplaceInference;
    %params.model.inferenceHandle=@PLDSVariationalInference;
    params.opts.algorithmic.EMIterations.maxCPUTime=MaxTime;
    params.opts.algorithmic.EMIterations.maxIter=MaxIter;
    [params,seq,LL]=PopSpikeEM(params,seq);
    for k=1:length(seq),TrainTraj{k}=seq(k).posterior.xsm;end
    
    seq_=[];
    for k=1:length(TestSpikeCount)
        seq_(k).y=TestSpikeCount{k};
        seq_(k).T=size(TestSpikeCount{k},2);
    end
    
    if isempty(seq_) 
        TestTraj =[];
    else
        InferenceMethod=params.model.inferenceHandle;
        seq_=InferenceMethod(params,seq_);
        for k=1:length(seq_),TestTraj{k}=seq_(k).posterior.xsm;end
    end
end

function [C, X, d] = ExpFamPCA(Y,xDim,varargin)
    %
    % [C, X, d] = ExpFamPCA(Y,xDim)
    %
    % exponential family PCA, currently only implemented for
    % exponential-Poisson observation model, i.e. learns C, X and d for model
    % Y ~ Poisson(exp(Cx+d+s));
    %
    % inputs:
    % Y:     matrix of size yDim x T
    % s:     additive , observed input, same size as Y or scalar (optional)
    % xDim:  scalar, dimensionality of latent space
    %
    % output:
    % C:      loading matrix, of size yDim x xDim
    % X:      recovered latent factors, of size xDim x T
    % d:      mean offset
    %
    %
    % (c) Lars Buesing 2014
    
    
    s                   = 0;
    dt                  = 10;    % rebin factor %!!! this is very much dependent on the firing rates
    lam                 = 1e-1;  % penalizer
    CposConstrain       = false; % constrain C to be have pos elements
    
    options.display     = 'none';
    options.MaxIter     = 10000;
    options.maxFunEvals = 50000;
    options.Method      = 'lbfgs';
    options.progTol     = 1e-9;
    options.optTol      = 1e-5;
    
    assignopts(who,varargin);
    
    
    seqDum.y = Y;
    seqDum = rebinRaster(seqDum,dt);
    Y = seqDum.y;
    
    if numel(s)>1.5; s = subsampleSignal(s,dt);end
    
    [yDim T] = size(Y);
    
    %% rough initialization for ExpFamPCA based on SVD
    my = mean(Y-s,2);
    [Uy Sy Vy] = svd(bsxfun(@minus,Y-s,my),0);
    my = max(my,0.1);
    
    Cinit = Uy(:,1:xDim);
    if CposConstrain
        Cinit = 0.1*randn(yDim,xDim);
    end
    
    Xinit = 0.01*randn(xDim,T);
    dinit = log(my);
    CXdinit = [vec([Cinit; Xinit']); dinit];
    
    %run ExpFamCPA
    CXdOpt  = minFunc(@ExpFamPCACost,CXdinit,options,Y,xDim,lam,s,CposConstrain);
    
    % Function returns all parameters lumped together as one vector, so need to disentangle:
    d  = CXdOpt(end-yDim+1:end);
    CX = reshape(CXdOpt(1:end-yDim),yDim+T,xDim);
    C  = CX(1:yDim,:);
    if CposConstrain; C = exp(C); end
    X  = CX(yDim+1:end,:)';
    d  = d-log(dt);
    Xm = mean(X,2);
    X  = bsxfun(@minus,X,Xm);
    d  = d+C*Xm;
    
    if ~CposConstrain
        % transform C to have orthonormal columns
        [UC SC VC] = svd(C);
        M = SC(1:xDim,1:xDim)*VC(:,1:xDim)';
        C = C/M;
        X = M*X;
        
        % transform for X to have orthogonal rows
        [MU MS MV] = svd((X*X')./T);
        M = MU';
        C = C/M;
        X = M*X;
    end
end

function [f df] = ExpFamPCACost(CXd,Y,xDim,lambda,s,CposConstrain)
    %
    % [f df] = ExpFamPCACost(CXd,Y,xDim,lambda)
    %
    % (c) L Buesing 2014
    
    [yDim T] = size(Y);
    
    
    d  = CXd(end-yDim+1:end);
    CX = reshape(CXd(1:end-yDim),yDim+T,xDim);
    C  = CX(1:yDim,:);
    if CposConstrain; C = exp(C); end;
    X  = CX(yDim+1:end,:)';
    
    nu = bsxfun(@plus,C*X+s,d);
    Yhat = exp(nu);
    
    f = sum(vec(-Y.*nu+Yhat))+lambda/2*(norm(C,'fro')^2+norm(X,'fro'));
    
    YhatmY = Yhat-Y;
    
    gX = C'*YhatmY+lambda*X;
    gC = YhatmY*X'+lambda*C;
    if CposConstrain; gC = gC.*C; end
    gd = sum(YhatmY,2);
    
    df = [vec([gC;gX']);gd];
    
end

function [Y,Xu,Xs,Xv,d,cost] = MODnucnrmminWithd( S, opts, varargin )
    %
    % [Y,Xu,Xs,Xv,d,cost] = MODnucnrmminWithd( S, opts, varargin )
    %
    %  obj function is -log(S|x) + lambda||X||_* + Tr[Z'(x-X)] + rho/2||x-X||^2_F
    %
    % opts -
    %   rho - dual gradient ascent rate for ADMM
    %   eps_abs - absolute threshold for ADMM
    %   eps_rel - relative threshold for ADMM
    %   lambda - strength of the nuclear norm penalty. multiplied by the square
    %       root of the number of elements in y to make sure it properly scales
    %       as the size of the problem increases
    %   maxIter - maximum number of iterations to run if the stopping threshold
    %       is not reached
    %   nlin - the nonlinearity to be used in the Poisson log likelihood:
    %       exp - exponential
    %       soft - exponential for negative x, linear otherwise
    %       logexp - log(1+exp(x)), a smooth version of "soft"
    %   verbose - verbosity level.
    %       0 - print nothing
    %       1 - print every outer loop of ADMM
    %       2 - print every inner loop of Newton's method and ADMM
    %
    % David Pfau, 2012-2013
    % minor modifications by Lars Buesing, 2014
    %
    
    Yinit = [];
    dinit = [];
    Yext  = [];
    
    f   = [];
    df  = [];
    d2f = [];
    
    
    assignopts(who,varargin);
    
    % set default values
    rho     = 1.3;
    eps_abs = 1e-10;%1e-6;
    eps_rel = 1e-5;%1e-3;
    maxIter = 250;
    nlin    = 'exp';
    lambda  = 1;
    verbose = 0;
    
    if nargin > 1
        for field = {'rho','eps_abs','eps_rel','maxIter','nlin','lambda','verbose'}
            if isfield(opts, field)
                eval([field{1} ' = opts.' field{1} ';'])
            end
        end
    end
    
    [N,T] = size(S);
    if isempty(Yext)
        Yext = zeros(size(S));
    end
    
    switch nlin
        case 'exp'
            f = @exp;
        case 'soft'
            f   = @(x) exp(x).*(x<0) + (1+x).*(x>=0);
            df  = @(x) exp(x).*(x<0) + (x>=0);
            d2f = @(x) exp(x).*(x<0);
        otherwise
            disp('NucNormMin: using user-defined link function')
    end
    
    if ~isempty(Yinit)
        x = Yinit;
    else
        switch nlin  % crude ADMM initialization
            case 'exp'
                x = log(max(S,1))-Yext;
            otherwise
                x = max(S-1,0)-Yext;
        end
    end
    
    if ~isempty(dinit)
        d = dinit;
        if isempty(Yinit)
            x = zeros(size(S));
        end
    else
        d = mean(x,2);
        x = bsxfun(@minus,x,d);
    end
    
    X = zeros(N,T);
    Z = zeros(N,T);
    
    
    nr_p = Inf; nr_d= Inf;
    e_p = 0;    e_d = 0;
    iter = 0;
    if verbose>0;
        fprintf('Iter:\t Nuc nrm:\t Loglik:\t Objective:\t dX:\t\t r_p:\t\t e_p:\t\t r_d:\t\t e_d:\n')
    end
    
    while ( nr_p > e_p || nr_d > e_d ) && iter < maxIter % Outer loop of ADMM
        stopping = Inf; % stopping criterion
        x_old = x;
        if verbose == 2, fprintf('\tNewton:\t Obj\t\t Stopping\t\n'); end
        while stopping/norm(x,'fro') > 1e-6 % Outer loop of Newton's method
            xEval = bsxfun(@plus,x+Yext,d);
            switch nlin
                case 'exp'
                    h =  exp( xEval ); % diagonal of Hessian
                    g =  exp( xEval ) - S; % gradient
                otherwise
                    h =  d2f( xEval ) - S .* ( d2f( xEval ) .* f( xEval ) - df( xEval ).^2 ) ./ f( xEval ).^2;
                    g =  df ( xEval ) - S .* df( xEval ) ./ f( xEval );
                    g(isnan(g)) = 0;
                    h(isnan(h)) = 0;
            end
            
            grad = g + rho*(x-X) + Z;
            dx = -inv_hess_mult(h,grad);
            
            %!!! fix here
            gd = sum(g,2);%+lambda*d;
            hd = sum(h,2);%+lambda;
            dd = -gd./hd;
            
            % upadate
            x = x + dx;
            d = d + dd;
            
            stopping = abs(grad(:)'*dx(:)+dd'*gd);
            if verbose == 2 % verbosity level
                fprintf('\t\t %1.2e \t %1.2e\n', obj(x), stopping)
            end
        end
        dx = norm(x_old-x,'fro')/norm(x,'fro');
        
        
        if T > N
            [v,s,u] = svd( x' + Z'/rho, 0 );
        else
            [u,s,v] = svd( x + Z/rho, 0);
        end
        
        Xs = max( s - eye(min(N,T))*lambda/rho, 0 );
        Xu = u;
        Xv = v;
        
        X_ = u*max( s - eye(min(N,T))*lambda/rho, 0 )*v';
        
        Z_ = Z + rho * ( x - X_ );
        
        % compute residuals and thresholds
        r_p = x - X_;
        r_d = rho * ( X - X_ );
        e_p = sqrt(N*T) * eps_abs + eps_rel * max( norm( x, 'fro' ), norm( X_, 'fro' ) );
        e_d = sqrt(N*T) * eps_abs + eps_rel * norm( Z , 'fro' );
        nr_p = norm( r_p, 'fro' );
        nr_d = norm( r_d, 'fro' );
        
        % heuristics to adjust dual gradient ascent rate to balance primal and
        % dual convergence. Seems to work pretty well: we almost always
        % converge in under 100 iterations.
        if nr_p > 10*nr_d
            rho = 2*rho;
        elseif nr_d > 10*nr_p
            rho = rho/2;
        end
        
        % update
        X = X_;
        Z = Z_;
        
        % !!! mod here, report stuff from low rank rates
        fval = sum( sum( f( bsxfun(@plus,X+Yext,d) ) - S .* log( f( bsxfun(@plus,X+Yext,d) ) ) ) );
        nn = sum( svd( X ) );
        
        cost.loglike   = -fval;
        cost.nucnorm   = nn;
        cost.total     = -cost.loglike+lambda*cost.nucnorm;
        
        % print
        iter = iter + 1;
        if verbose>0
            fprintf('%i\t %1.4e\t %1.4e\t %1.4e\t %1.4e\t %1.4e\t %1.4e\t %1.4e\t %1.4e\n',iter, nn, fval, lambda*nn+fval, dx, nr_p, e_p, nr_d, e_d);
        end
        
    end
    
    Y = X;   %!!! returns low rank rates
    sqrtXs = diag(sqrt(diag(Xs)));
    cost.penaltyC = norm(Xu*sqrtXs,'fro')^2;
    cost.penaltyX = norm(Xv*sqrtXs,'fro')^2;
    
    
    function y = inv_hess_mult(H,x)
        y = x./(H + rho);
    end
    function y = obj(x)
        foo = S.*log(f( bsxfun(@plus,x+Yext,d) ));
        foo(isnan(foo)) = 0; % 0*log(0) = 0 by convention
        y = sum(sum(f( bsxfun(@plus,x+Yext,d)  ) - foo)) + sum(sum(Z.*x)) + rho/2*norm(x-X,'fro')^2;
    end
    
end

function [Y,D] = nucnrmmin( S, opts )
    % x = nucnrmmin( Y, opts )
    %   Find the minimum of lambda*||A(Y-DH)||_* + gamma*||D||_1 + f(Y)
    %   where A(Y) is the operator that mean-centers Y (that is, we are looking
    %   for a low-dimensional *affine* space rather than a low dimensional
    %   subspace, as a way of accounting for mean firing rates.) and H is a
    %   matrix of spike history terms (ignored if k = 0). This focuses strictly
    %   on the Poisson case for simplicity. The algorithm follows the ADMM
    %   method introduced by Liu, Hansson and Vandenberghe.
    %
    % opts -
    %   rho - dual gradient ascent rate for ADMM
    %   eps_abs - absolute threshold for ADMM
    %   eps_rel - relative threshold for ADMM
    %   lambda - strength of the nuclear norm penalty. multiplied by the square
    %       root of the number of elements in y to make sure it properly scales
    %       as the size of the problem increases
    %   maxIter - maximum number of iterations to run if the stopping threshold
    %       is not reached
    %   nlin - the nonlinearity to be used in the Poisson log likelihood:
    %       exp - exponential
    %       soft - exponential for negative x, linear otherwise
    %       logexp - log(1+exp(x)), a smooth version of "soft"
    %   center - if true, miniminze the nuclear norm of the mean-centered data.
    %       if false, just do vanilla nuclear norm minimization (A(x) = x)
    %   q - the number of history steps to fit (if zero, ignore history)
    %   gamma - the weight on the \ell_1 term, if k > 0
    %   verbose - verbosity level.
    %       0 - print nothing
    %       1 - print every outer loop of ADMM
    %       2 - print every inner loop of Newton's method and ADMM
    %
    % David Pfau, 2012-2013
    
    % set default values
    rho     = 1.3;
    eps_abs = 1e-6;
    eps_rel = 1e-3;
    maxIter = 250;
    nlin    = 'logexp';
    lambda  = 1;
    center  = 1;
    verbose = 1;
    q       = 0;
    gamma   = 1;
    if nargin > 1
        for field = {'rho','eps_abs','eps_rel','maxIter','nlin','lambda','center','k','gamma','verbose'}
            if isfield(opts, field)
                eval([field{1} ' = opts.' field{1} ';'])
            end
        end
    end
    
    nz = logical(sum(S,2));
    S = S(nz,:); % remove rows with no spikes
    [N,T] = size(S);
    T = T-q;
    H = zeros(N*q,T); % history term
    for i = 1:q
        H((i-1)*N+(1:N),:) = S(:,i:end-q+i-1);
    end
    S = S(:,q+1:end);
    
    if center && strcmpi(nlin,'soft')
        % There's probably some alternative to the standard Newton step that
        % can take advantage of what we know about the nonlinearity, but I'll
        % save that for another day.
        error('Newton''s method is ill-conditioned when combining mean-centering with locally-flat nonlinearities!')
    end
    lambda = lambda * sqrt(N*T);
    gamma  = gamma  * T/N;
    
    switch nlin
        case 'exp'
            f = @exp;
        case 'soft'
            f   = @(x) exp(x).*(x<0) + (1+x).*(x>=0);
            df  = @(x) exp(x).*(x<0) + (x>=0);
            d2f = @(x) exp(x).*(x<0);
        case 'logexp'
            f   = @(x) log(1+exp(x));
            df  = @(x) 1./(1+exp(-x));
            d2f = @(x) exp(-x)./(1+exp(-x)).^2;
    end
    
    if center
        A = @(x) bsxfun(@minus,x,mean(x,2));
    else
        A = @(x) x;
    end
    A_adj = A;
    
    switch nlin  % crude ADMM initialization
        case 'exp'
            x = log(max(S,1));
        case {'soft','logexp'}
            x = max(S-1,0);
    end
    
    X = zeros(N,T);
    Z = zeros(N,T);
    D = zeros(N,N*q);
    G = zeros(N,N*q); % lagrange multiplier for inner loop of ADMM, mostly used for initialization near optimal point
    
    nr_p = Inf; nr_d = Inf;
    e_p = 0;    e_d = 0;
    iter = 0;
    if q == 0
        fprintf('Iter:\t Nuc nrm:\t Loglik:\t Objective:\t dX:\t\t r_p:\t\t e_p:\t\t r_d:\t\t e_d:\n')
    else
        fprintf('Iter:\t Nuc nrm:\t Loglik:\t ||D||_1:\t Objective:\t dX:\t\t r_p:\t\t e_p:\t\t r_d:\t\t e_d:\n')
    end
    while ( nr_p > e_p || nr_d > e_d ) && iter < maxIter % Outer loop of ADMM
        stopping = Inf; % stopping criterion
        x_old = x;
        if verbose == 2, fprintf('\tNewton:\t Obj\t\t Stopping\t\n'); end
        while stopping/norm(x,'fro') > 1e-6 % Outer loop of Newton's method
            switch nlin
                case 'exp'
                    h =  exp( x ); % diagonal of Hessian
                    g =  exp( x ) - S; % gradient
                otherwise
                    h =  d2f( x ) - S .* ( d2f( x ) .* f( x ) - df( x ).^2 ) ./ f( x ).^2;
                    g =  df ( x ) - S .* df( x ) ./ f( x );
                    g(isnan(g)) = 0;
                    h(isnan(h)) = 0;
            end
            
            grad = g + rho*A( x ) - A_adj( rho * X + rho * A(D*H) - Z ); % update this to include spike history term.
            dx = -inv_hess_mult(h,grad);
            x = x + dx;
            stopping = abs(grad(:)'*dx(:));
            if verbose == 2 % verbosity level
                fprintf('\t\t %1.2e \t %1.2e\n', obj(x), stopping)
            end
        end
        dx = norm(x_old-x,'fro')/norm(x,'fro');
        
        if q > 0 % run ADMM for LASSO
            [D,G] = admm_lasso_mat( H - mean(H,2)*ones(1,T), A_adj( A(x) - X + Z/rho ), gamma/rho, D, G, verbose == 2 );
        end
        
        Ax_ = A( x - D*H );
        
        if T > N
            [v,s,u] = svd( Ax_' + Z'/rho, 0 );
        else
            [u,s,v] = svd( Ax_ + Z/rho, 0);
        end
        X_ = u*max( s - eye(min(N,T))*lambda/rho, 0 )*v';
        
        Z_ = Z + rho * ( Ax_ - X_ );
        
        % compute residuals and thresholds
        r_p = Ax_ - X_;
        r_d = rho * A_adj( X - X_ );
        e_p = sqrt(N*T) * eps_abs + eps_rel * max( norm( Ax_, 'fro' ), norm( X_, 'fro' ) );
        e_d = sqrt(N*T) * eps_abs + eps_rel * norm( A_adj( Z ), 'fro' );
        nr_p = norm( r_p, 'fro' );
        nr_d = norm( r_d, 'fro' );
        
        % heuristics to adjust dual gradient ascent rate to balance primal and
        % dual convergence. Seems to work pretty well: we almost always
        % converge in under 100 iterations.
        if nr_p > 10*nr_d
            rho = 2*rho;
        elseif nr_d > 10*nr_p
            rho = rho/2;
        end
        
        % update
        X = X_;
        Z = Z_;
        
        fval = sum( sum( f( x ) - S .* log( f( x ) ) ) );
        nn = sum( svd( A( x ) ) );
        
        % print
        iter = iter + 1;
        if q == 0
            fprintf('%i\t %1.4e\t %1.4e\t %1.4e\t %1.4e\t %1.4e\t %1.4e\t %1.4e\t %1.4e\n', ...
                iter, nn, fval, lambda*nn+fval, dx, nr_p, e_p, nr_d, e_d);
        else
            fprintf('%i\t %1.4e\t %1.4e\t %1.4e\t %1.4e\t %1.4e\t %1.4e\t %1.4e\t %1.4e\t %1.4e\n', ...
                iter, nn, fval, sum(abs(D(:))), lambda*nn+gamma*sum(abs(D(:)))+fval, dx, nr_p, e_p, nr_d, e_d);
        end
    end
    Y = zeros(length(nz),T);
    Y( nz,:) = x;
    Y(~nz,:) = -Inf; % set firing rates to zero for rows with no data. Used to make sure the returned value is aligned with the input
    
    function y = inv_hess_mult(H,x)
        % For the Newton step of ADMM, we must multiply a vector by the inverse of
        % H + rho*A'*A, where H is the Hessian of the smooth penalty term (which in
        % our case is diagonal and therefore simple) and A is the linear operator
        % in the objective ||A(x)||_* + f(x). In our problem the linear operator is
        % mean-centering the matrix x, which is a symmetric and idempotent operator
        % that can be written out in matrix form for m-by-n matrices as:
        %
        % A = eye(m*n) - 1/n*kron(ones(n,1),eye(m))*kron(ones(n,1),eye(m))'
        %
        % and so (H + rho*A*A')^-1*x can be efficiently computed by taking
        % advantage of Woodbury's lemma.
        
        if center % update this to include spike history term
            Hi = 1./(H + rho);
            foo = (sum(Hi.*x,2)./(ones(N,1)-sum(Hi,2)*rho/T));
            y = Hi.*x + rho/T*Hi.*(foo*ones(1,T));
        else
            y = x./(H + rho);
        end
    end
    
    function y = obj(x)
        
        foo = S.*log(f(x));
        foo(isnan(foo)) = 0; % 0*log(0) = 0 by convention
        y = sum(sum(f(x) - foo)) + sum(sum(Z.*A(x-D*H))) + rho/2*norm(A(x-D*H)-X,'fro')^2;
        
    end
end

function penalty = PLDSemptyParamPenalizerHandle(params)
    %
    % function penalty = PLDSemptyParamPenalizerHandle(params)
    %
    
    penalty = 0;
end

function params = PLDSgenerateExample(varargin)
    %
    % trueparams = PLDSgenerateExample(varargin)
    %
    % generate a random PLDS model based on some inputs,
    % using Poisson or Bernoulli observations
    
    params   = [];
    useR     = false;
    uDim     = 0;
    xDim     = 10;
    yDim     = 100;
    
    Aspec    = 0.99;
    Arand    = 0.03;
    Q0max    = 0.3;
    BernFlag = false;
    doff     = -1.9;
    statFlag = false;
    Bscale   = 1.0;
    
    assignopts(who,varargin);
    
    
    %%%%%%%%%  generate parameters %%%%%%%%%
    
    A  = eye(xDim)+Arand*randn(xDim);
    A  = A./max(abs(eig(A)))*Aspec;
    Q  = diag(rand(xDim,1));
    Q0 = dlyap(A,Q);
    M  = diag(1./sqrt(diag(Q0)));
    A  = M*A*pinv(M);
    Q  = M*Q*M'; Q=(Q+Q)/2;
    
    O  = orth(randn(xDim));
    Q0 = O*diag(rand(xDim,1)*Q0max)*O'/3;
    x0 = randn(xDim,1)/3;
    
    C  = randn(yDim,xDim)./sqrt(3*xDim);
    d  = 0.3*randn(yDim,1)+doff;
    
    params.model.A    = A;
    params.model.Q    = Q;
    params.model.Q0   = Q0;
    params.model.x0   = x0;
    params.model.C    = C;
    params.model.d    = d;
    
    
    if BernFlag
        params.model.dualHandle = @LogisticBernoulliDualHandle;
        params.model.likeHandle = @LogisticBernoulliHandle;
    end
    
    if uDim>0
        cQ = max(abs(diag(chol(params.model.Q))));
        params.model.B = cQ*(rand(xDim,uDim)+0.5)/(uDim)*Bscale;
        params.model.notes.useB = true;
    end
    
    if statFlag
        params.model.x0 = zeros(xDim,1);
        params.model.Q0 = dlyap(params.model.A,params.model.Q);
    end
    
    params = PLDSsetDefaultParameters(params,xDim,yDim);
end

function params = PLDSInitialize(seq,xDim,initMethod,params)
    %
    % function params = PLDSInitialize(params,seq)
    %
    % inititalize parameters of Population-LDS model with different methods. At
    % the moment, focusses on exponential link function and Poisson
    % observations.
    %
    %input:
    % seq:      standard data struct
    % xdim:     desired latent dimensionality
    % initMethods:
    %
    % - params											% just initialize minimal undefiened fields with standard values
    % - PLDSID											% moment conversion + Ho-Kalman SSID
    % - ExpFamPCA											% exponential family PCA
    % - NucNormMin											% nuclear norm minimization, see [Robust learning of low-dimensional dynamics from large neural ensembles David Pfau, Eftychios A. Pnevmatikakis, Liam Paninski. NIPS2013]
    % params: if initialization method 'params' is chosen, the params-struct
    % that one should use
    %
    % (c) L Buesing 2014
    
    
    yDim       = size(seq(1).y,1);
    Trials     = numel(seq);
    params.opts.initMethod = initMethod;
    params     = PLDSsetDefaultParameters(params,xDim,yDim);					% set standard parameter values
    
    
    switch initMethod
        
        case 'params'
            % do nothing
            disp('Initializing PLDS parameters with given parameters')
            
            
        case 'PLDSID'
            % !!! debugg SSID stuff separately & change to params.model convention
            disp('Initializing PLDS parameters using PLDSID')
            if params.model.notes.useB
                warning('SSID initialization with external input: not implemented yet!!!')
            end
            PLDSIDoptions = struct2arglist(params.opts.algorithmic.PLDSID);
            params.model = FitPLDSParamsSSID(seq,xDim,'params',params.model,PLDSIDoptions{:});
            
            
        case 'ExpFamPCA'
            % this replaces the FA initializer from the previous verions...
            disp('Initializing PLDS parameters using exponential family PCA')
            
            dt = params.opts.algorithmic.ExpFamPCA.dt;
            Y  = [seq.y];
            if params.model.notes.useS; s = [seq.s];
            else s=0;end
            [Cpca, Xpca, dpca] = ExpFamPCA(Y,xDim,'dt',dt,'lam',params.opts.algorithmic.ExpFamPCA.lam,'options',params.opts.algorithmic.ExpFamPCA.options,'s',s);
            params.model.Xpca = Xpca;
            params.model.C = Cpca;
            params.model.d = dpca;
            
            if params.model.notes.useB; u = [seq.u];else;u = [];end
            params.model = LDSObservedEstimation(Xpca,params.model,dt,u);
            
            
        case 'NucNormMin'
            disp('Initializing PLDS parameters using Nuclear Norm Minimization')
            
            dt = params.opts.algorithmic.NucNormMin.dt;
            seqRebin.y = [seq.y]; seqRebin = rebinRaster(seqRebin,dt);
            Y  = [seqRebin.y];
            options = params.opts.algorithmic.NucNormMin.options;
            options.lambda = options.lambda*sqrt(size(Y,1)*size(Y,2));
            if params.model.notes.useS
                Yext = subsampleSignal([seq.s],dt);
            else
                Yext = [];
            end
            [Y,Xu,Xs,Xv,d] = MODnucnrmminWithd( Y, options , 'Yext', Yext );
            params.model.d = d-log(dt);
            
            if ~params.opts.algorithmic.NucNormMin.fixedxDim
                disp('Variable dimension; still to implement!')
            else
                params.model.C = Xu(:,1:xDim)*Xs(1:xDim,1:xDim);
                if params.model.notes.useB; u = [seq.u];else;u = [];end
                params.model = LDSObservedEstimation(Xv(:,1:xDim)',params.model,dt,u);
                params.model.Xpca = Xv(:,1:xDim)';
                params.model.Xs   = diag(Xs(1:xDim,1:xDim));
            end
            
            
        otherwise
            warning('Unknown PLDS initialization method')
            
    end
    
    if params.model.notes.useB && (numel(params.model.B)<1)
        params.model.B = zeros(xDim,size(seq(1).u,1));
    end
    
    params = LDSTransformParams(params,'TransformType',params.opts.algorithmic.TransformType);	% clean up parameters
    params.modelInit = params.model;
end

function [f df ddf] = PLDSLaplaceCost(x,y,Lambda,mu,W,d)
    %
    %
    %
    over_m = W*x+d;
    yhat   = exp(over_m);
    xmu    = x-mu;
    Lamxmu = Lambda*mu;
    yhaty  = yhat-y;
    
    f   = 0.5*xmu'*Lamxmu + y'*over_m-sum(yhat);
    
    df  = Lamxmu + W'*(yhat-y);
    
    ddf = Lambda + W'*diag(yhat-y)*W;
end

function [seq varBound] = PLDSLaplaceInference(params,seq)
    %
    % function seq = PLDSlpinf(params,seq)
    %
    
    computeVarBound = true; %!!! make this nice
    
    Trials      = numel(seq);
    [yDim xDim] = size(params.model.C);
    varBound    = 0;
    
    mps = params.model;
    mps.E = [];
    mps.D = [];
    mps.initV = mps.Q0;
    mps.initx = mps.x0;
    
    runinfo.nStateDim = xDim;
    runinfo.nObsDim   = yDim;
    
    infparams.initparams = mps;
    infparams.runinfo    = runinfo;
    infparams.notes      = params.model.notes;
    
    Tmax = max([seq.T]);
    for tr=1:Trials
        T = size(seq(tr).y,2);
        indat = seq(tr);
        if isfield(seq(tr),'posterior') && isfield(seq(tr).posterior,'xsm')
            indat.xxInit = seq(tr).posterior.xsm;
        else
            indat.xxInit = getPriorMeanLDS(params,T,'seq',seq(tr));
        end
        seqRet = PLDSLaplaceInferenceCore(indat, infparams);
        seq(tr).posterior.xsm      = seqRet.x;
        seq(tr).posterior.Vsm      = reshape(seqRet.V,xDim,xDim*T)';
        seq(tr).posterior.VVsm     = reshape(permute(seqRet.VV(:,:,2:end),[2 1 3]),xDim,xDim*(T-1))';
        seq(tr).posterior.lamOpt   = exp(vec(seqRet.Ypred));
    end
    
    
    if computeVarBound
        
        VarInfparams    = params.model;
        VarInfparams.CC = zeros(xDim,xDim,yDim);
        for yy=1:yDim
            VarInfparams.CC(:,:,yy) = params.model.C(yy,:)'*params.model.C(yy,:);
        end
        VarInfparams.CC = reshape(VarInfparams.CC,xDim^2,yDim);
        
        Cl = {}; for t=1:Tmax; Cl = {Cl{:} params.model.C}; end
        Wmax = sparse(blkdiag(Cl{:}));
        
        % iterate over trials
        
        for tr = 1:Trials
            
            T = size(seq(tr).y,2);
            VarInfparams.d = repmat(params.model.d,T,1); if params.model.notes.useS; VarInfparams.d = VarInfparams.d + vec(seq(tr).s); end
            VarInfparams.mu         = vec(getPriorMeanLDS(params,T,'seq',seq(tr)));
            VarInfparams.W          = Wmax(1:yDim*T,1:xDim*T);
            VarInfparams.y          = seq(tr).y;
            VarInfparams.Lambda     = buildPriorPrecisionMatrixFromLDS(params,T);
            VarInfparams.WlamW      = sparse(zeros(xDim*T));
            VarInfparams.dualParams = [];
            
            if isfield(params.model,'baseMeasureHandle')
                VarInfparams.DataBaseMeasure = feval(params.model.baseMeasureHandle,seq(tr).y,params);
                seq(tr).posterior.DataBaseMeasure = VarInfparams.DataBaseMeasure;
            end
            
            lamOpt = seq(tr).posterior.lamOpt;
            [DualCost, ~, varBound] = VariationalInferenceDualCost(lamOpt,VarInfparams);
            seq(tr).posterior.varBound = varBound;
            
        end
        varBound = 0;
        for tr=1:Trials; varBound = varBound + seq(tr).posterior.varBound; end;
    end
    
end

function smooth = PLDSLaplaceInferenceCore(data, params)
    % Compute laplace-approximated posterior means and covariances for PLDS
    % (c) Evan Archer, 2014
    %
    runinfo = params.runinfo;
    mps = params.initparams;
    
    % initialize latent states to 0
    %xx = randn(runinfo.nStateDim, size(data.y,2)); % for subsequent EM iterations in the poisson case we'll want to initialize with previous x's.
    
    xx = data.xxInit;
    
    if params.notes.useB
        H = mps.B*data.u;% zeros(size(xx)); %mps.nlin.f(mps.nlin, data.h); % compute nlin first.
    else
        H = zeros(size(xx));
    end
    
    Qinv = pinv(mps.Q);
    Q0inv = pinv(mps.initV);
    AQiA = mps.A'*Qinv*mps.A;
    
    LL = -inf;
    LLprev = -inf;
    while 1
        
        XX = xx(:,2:end) - mps.A*xx(:,1:end-1) - H(:,2:end);
        QiH = zeros(size(H));
        QiH(:,2:end) = - Qinv*XX;
        QiH(:,1) = - Q0inv*( xx(:,1) - H(:,1) - mps.initx );
        
        T = size(data.y,2);
        
        % sum over s and u
        if params.notes.useS
            ds = data.s;
        else
            ds = sparse(size(data.y,1), size(data.y,2));
        end
        
        
        Ypred = bsxfun(@plus, mps.C*xx + ds, mps.d);
        
        %% Latent-state grad & hessian (this should be the same for both poisson and gaussian likelihoods)
        
        lat_grad =  [ QiH(:,1:end-1) - mps.A'*QiH(:,2:end),   QiH(:,end) ];
        
        II = speye(size(data.y,2)-2);
        lat_hess_diag = -blkdiag(sparse(Q0inv+AQiA), kron( II, Qinv + AQiA), sparse(Qinv));
        
        II_c = circshift(speye(size(data.y,2)), [0 1]); II_c(end,1) = 0;
        lat_hess_off_diag = kron(II_c, mps.A'*Qinv); lat_hess_off_diag = lat_hess_off_diag + lat_hess_off_diag';
        
        lat_hess = lat_hess_diag  + lat_hess_off_diag;
        
        %% Poisson Observation gradient and hessian
        
        Lambda = exp(Ypred);
        YL = data.y-Lambda;
        
        YC = zeros(runinfo.nObsDim, size(data.y,2), runinfo.nStateDim);
        poiss_hess_blk = zeros(runinfo.nStateDim, runinfo.nStateDim, size(data.y,2));
        for idx = 1:runinfo.nStateDim
            YC(:,:,idx) = bsxfun(@times, YL, mps.C(:,idx));
            poiss_hess_blk(idx,:,:) = -mps.C'*bsxfun(@times, Lambda, mps.C(:,idx));
        end
        
        poiss_grad = sparse(squeeze(sum(YC,1))');
        if runinfo.nStateDim==1
            poiss_grad = poiss_grad';
        end
        poiss_hess = spblkdiag(poiss_hess_blk);
        
        %% Add the latent and observation hessian & gradient
        
        Hess = poiss_hess + lat_hess;
        
        Grad = poiss_grad + lat_grad;
        
        %% Compute newton step, perform line search
        
        if LL - LLprev < 1e-10 % TODO: put tolerance in setup file
            break
        end
        
        xold = xx;
        UPDATE  = reshape(Hess \ Grad(:), size(lat_grad));
        
        dx = 0;
        LLprev = LL;
        ii = 0;
        
        LOGDETS = logdet(mps.initV) + (T-1)*logdet(mps.Q);
        while 1
            dx = 2^(-ii); ii = ii + .1;
            xx = xold - dx*UPDATE;
            
            %% Recompute just the likelihood @ dx step
            Ypred = bsxfun(@plus, mps.C*xx + ds, mps.d);
            Lambda = exp(Ypred);
            
            if dx < .001
                break
            end
            
            % Let's see if we can't compute the likelihood.
            % fprintf('\ncomputing likelihood:\n')
            
            XX = xx(:,2:end) - mps.A*xx(:,1:end-1) - H(:,2:end);
            X0 = (xx(:,1) - mps.initx - H(:,1));
            
            % Dropping constant terms that depend on dimensionality of states and observations
            % (might want to change that for model comparison purposes later)
            PLL = sum(sum(data.y.*Ypred-Lambda));
            GnLL = LOGDETS + sum(sum(XX.*(Qinv*XX))) + sum(sum(X0.*(Q0inv*X0)));
            
            LL = PLL - GnLL/2;
            %% Finish recomputing the likelihood
            
            if LL > LLprev
                break
            end
        end
        
        %fprintf('%0.2f --> %0.3f', dx, LL);
        
    end
    
    %%
    
    AA0 = -(Q0inv + AQiA);
    AAE = -Qinv;
    AA = -(Qinv + AQiA);
    BB = (mps.A'*Qinv);
    smooth = [];
    smooth.VV = zeros(size(AA,1),size(AA,1),T);
    AAz = poiss_hess_blk; AAz(:,:,1) = AA0 + AAz(:,:,1); AAz(:,:,2:end-1) = bsxfun(@plus, AA, poiss_hess_blk(:,:,2:end-1)); AAz(:,:,end) = AAE + AAz(:,:,end);
    
    AAz = -AAz; BB = -BB;
    
    [smooth.V,smooth.VV(:,:,2:end)]=sym_blk_tridiag_inv_v1(AAz, BB, [1:T]', ones(T-1,1));
    
    smooth.x = xx;
    smooth.loglik = LL;
    smooth.Ypred = Ypred;
    
    %% With minfunc
    %
    % opt = struct('Method','csd', 'maxFunEvals', 500, 'Display', 'on'); % TODO: Make this an option in opt struc
    % xxmin = minFunc( @(x) poiss_obj_fun(x, data, params), xx(:), opt);
    
    
end

function [params seq] = PLDSMStep(params,seq)
    %
    % params = PLDSMStep(params,seq)
    %
    
    
    params = LDSMStepLDS(params,seq);
    params = PLDSMStepObservation(params,seq);
    
    params = LDSTransformParams(params,'TransformType',params.opts.algorithmic.TransformType);
end

function params = PLDSMStepObservation(params,seq)
    %
    % function params = PLDSMStepObservation(params,seq)
    %
    
    
    minFuncOptions = params.opts.algorithmic.MStepObservation.minFuncOptions;
    
    [yDim xDim] = size(params.model.C);
    
    if params.model.notes.useCMask; params.model.C = params.model.C.*params.model.CMask; end
    CdInit = vec([params.model.C params.model.d]); % warm start at current parameter values
    MStepCostHandle = @PLDSMStepObservationCost;
    
    %%% optimization %%%
    
    CdOpt = minFunc(MStepCostHandle,CdInit,minFuncOptions,seq,params);
    CdOpt = reshape(CdOpt,yDim,xDim+1);
    
    params.model.C = CdOpt(:,1:xDim);
    if params.model.notes.useCMask; params.model.C = params.model.C.*params.model.CMask; end
    params.model.d = CdOpt(:,end);
    
    
end

function [f, df] = PLDSMStepObservationCost(vecCd,seq,params)
    %
    % function [f, df] = PLDSMStepObservationCost(vecCd,seq,params)
    %
    % Mstep for observation parameters C,d for standard PLDS with exp-Poisson observations
    %
    % Input:
    %	- convention Cd = [C d]  and vecCd = vec(Cd)
    %
    % to do:
    %
    %       0) analyze run time
    %
    %
    % (c) L Buesing 2014
    
    
    Trials  = numel(seq);
    yDim    = size(seq(1).y,1);
    xDim    = size(params.model.A,1);
    
    
    CdMat   = reshape(vecCd,yDim,xDim+1);
    C       = CdMat(:,1:xDim);
    if params.model.notes.useCMask; C = C.*params.model.CMask; end
    d       = CdMat(:,end);
    
    CC      = zeros(yDim,xDim^2);
    for yd=1:yDim
        CC(yd,:) = vec(C(yd,:)'*C(yd,:));
    end
    
    
    f   = 0;				% current value of the cost function = marginal likelihood
    df  = zeros(size(C));			% derviative wrt C
    dfd = zeros(yDim,1);			% derivative wrt d
    
    for tr=1:Trials
        
        T    = size(seq(tr).y,2);
        y    = seq(tr).y;
        m    = seq(tr).posterior.xsm;
        Vsm  = reshape(seq(tr).posterior.Vsm',xDim.^2,T);
        
        h    = bsxfun(@plus,C*m,d);
        if params.model.notes.useS; h = h+seq(tr).s; end
        rho  = CC*Vsm;
        
        
        yhat = exp(h+rho/2);
        f    = f+sum(vec(y.*h-yhat));
        
        TT   = yhat*Vsm';
        TT   = reshape(TT,yDim*xDim,xDim);
        TT   = squeeze(sum(reshape(bsxfun(@times,TT,vec(C)),yDim,xDim,xDim),2));
        
        df   = df  + (y-yhat)*m'-TT;
        dfd  = dfd + sum((y-yhat),2);
        
    end
    
    f  = -f;
    if params.model.notes.useCMask; df = df.*params.model.CMask; end
    df = -vec([df dfd]);
    
end

function [ypred xpred xpredCov seqInf] = PLDSPredictRange(params,y,condRange,predRange,varargin);
    %
    % function [ypred xpred xpredCov] = PLDSPredictRange(params,y);
    %
    %
    % Lars Buesing, 2014
    %
    
    s       = [];
    u       = [];
    lamInit = [];
    
    assignopts(who,varargin);
    
    [yDim xDim] = size(params.model.C);
    
    Tcond = numel(condRange);
    Tpred = numel(predRange);
    tcRlo = min(condRange);
    tcRhi = max(condRange);
    tpRlo = min(predRange);
    tpRhi = max(predRange);
    
    if size(y,2)<Tcond
        error('Conditioning range larger than data, aborting')
    elseif size(y,2)>Tcond
        y = y(:,condRange);
        if params.model.notes.useB
            ucond = u(:,condRange);
        end
        if params.model.notes.useS
            scond = s(:,condRange);
        end
    end
    
    
    paramsInf = params;
    if paramsInf.model.notes.useB; paramsInf.model.x0 = paramsInf.model.x0+paramsInf.model.B*u(:,1);end
    for t=2:tcRlo                                                   % get the starting distribution right
        paramsInf.model.x0 = paramsInf.model.A*paramsInf.model.x0;
        if paramsInf.model.notes.useB&&(t<tcRlo); paramsInf.model.x0 = paramsInf.model.x0+paramsInf.model.B*u(:,t);end
        paramsInf.model.Q0 = paramsInf.model.A*paramsInf.model.Q0*paramsInf.model.A'+paramsInf.model.Q;
    end
    
    seqInf.y  = y;
    seqInf.T  = size(seqInf.y,2);
    if numel(lamInit)==(yDim*Tcond)
        disp('Warm-starting predction inference')
        seqInf.posterior.lamOpt = lamInit;
    end
    if paramsInf.model.notes.useB; seqInf.u = ucond; end;
    if paramsInf.model.notes.useS; seqInf.s = scond; end;
    
    seqInf   = params.model.inferenceHandle(paramsInf,seqInf);
    xpred    = zeros(xDim,Tpred);
    xpredCov = zeros(xDim,xDim,Tpred);
    ypred    = zeros(yDim,Tpred);
    
    xNow     = seqInf(1).posterior.xsm(:,end);
    xCovNow  = seqInf(1).posterior.Vsm(end+1-xDim:end,:);
    
    for t = (tcRhi+1):tpRhi    % progagate prediction
        xNow    = paramsInf.model.A*xNow;
        if params.model.notes.useB; xNow = xNow+params.model.B*u(:,t); end;
        xCovNow = paramsInf.model.A*xCovNow*paramsInf.model.A'+paramsInf.model.Q;
        if t>=tpRlo
            xpred(:,t-tpRlo+1) = xNow;
            xpredCov(:,:,t-tpRlo+1) = xCovNow;
            yr = params.model.C*xNow+params.model.d+0.5*diag(params.model.C*xCovNow*params.model.C');
            if params.model.notes.useS
                yr = yr+s(:,t);
            end
            ypred(:,t-tpRlo+1) = exp(yr);
        end
    end
end

function seq = PLDSsample(params,T,Trials,varargin)
    %
    % seq = PLDSsample(params,T,Trials)
    %
    % sample from a PLDS model with exponential nonlinearity (or
    % user defined link function); uses LDSsample
    %
    % (c) L Buesing 2014
    %
    
    
    
    yMax = inf;
    assignopts(who,varargin);
    
    if isfield(params.model,'linkFunc')
        linkFunc = params.model.linkFunc;
        %disp('Using non-exp link function for sampling')
    else
        linkFunc = @exp;
    end
    
    seq = LDSsample(params,T,Trials,varargin{:});
    
    
    for tr=1:Trials
        seq(tr).yr   = real(seq(tr).y);
        seq(tr).y    = poissrnd(linkFunc(seq(tr).yr));
        seq(tr).y(:) = min(yMax,seq(tr).y(:));
    end
    
    
end

function params = PLDSsetDefaultParameters(params,xDim,yDim)
    %
    % params = PLDSsetDefaultParameters(params,xDim,yDim)
    %
    %
    % Lars Buesing, Jakob H Macke
    %
    
    
    %%%%%%%%%%% set standard parameters %%%%%%%%%%%%%%%%%%%%
    % these standard settings make sense for data binned at 10ms with average rates of roughly 10Hz
    params = touchField(params,'model');
    params.model = touchField(params.model,'A',0.9*eye(xDim));    % dynamics matrix A
    params.model = touchField(params.model,'Q',(1-0.9.^2)*eye(xDim)); %innovation covariance Q
    params.model = touchField(params.model,'Q0',eye(xDim)); %initial state covariance Q0
    params.model = touchField(params.model,'x0',zeros(xDim,1)); %initial mean x0
    params.model = touchField(params.model,'C',randn(yDim,xDim)./sqrt(xDim)); %couplings from latent to observed C
    params.model = touchField(params.model,'d',zeros(yDim,1)-2.0); %mean-controlling offset d for each neuron
    params.model = touchField(params.model,'B',zeros(xDim,0));
    
    params.model = touchField(params.model,'notes');
    params.model.notes = touchField(params.model.notes,'learnx0', true);
    params.model.notes = touchField(params.model.notes,'learnQ0', true);
    params.model.notes = touchField(params.model.notes,'learnA',  true);
    params.model.notes = touchField(params.model.notes,'learnR',  false);
    params.model.notes = touchField(params.model.notes,'useR',    false);
    params.model.notes = touchField(params.model.notes,'useB',    false);
    params.model.notes = touchField(params.model.notes,'useS',    false);
    params.model.notes = touchField(params.model.notes,'useCMask',false);
    
    
    %%%%%%%%%%% set standard observation model handles for variational inference %%%%%%%%%%%%%%%%%%%%
    %note: these functions all have to be consistent with each other, i.e.
    %the likeHandle, dualHandle, domainHandle, baseMeasureHandle, MstepHandle all
    %have to be corresponding to the same likelihood-model and inference
    %procedure, otherwise funny things will happen.
    %at the moment, only Poisson model with exponential nonlinarity is
    %implemented
    params.model = touchField(params.model,'likeHandle',       @ExpPoissonHandle); %use exponential Poisson likelihood
    params.model = touchField(params.model,'dualHandle',       @ExpPoissonDualHandle); %and its dual
    params.model = touchField(params.model,'domainHandle',     @ExpPoissonDomain); %specify the domain of addmissable parameters
    params.model = touchField(params.model,'baseMeasureHandle',@PoissonBaseMeasure); %base measure, i.e. constant part which does not need to be evaluated at each step
    params.model = touchField(params.model,'inferenceHandle',  @PLDSVariationalInference); % function that does the actual inference
    params.model = touchField(params.model,'MStepHandle',      @PLDSMStep); %handle to function that does the M-step
    params.model = touchField(params.model,'ParamPenalizerHandle',@PLDSemptyParamPenalizerHandle);
    
    
    %%%%%%%%%%% set standard algorithmic parameters %%%%%%%%%%%%%%%%%%%%
    
    params = touchField(params,'opts');
    params.opts = touchField(params.opts,'algorithmic');
    
    
    %%%% set parameters for Variational Inference %%%%
    %these parameters are handed over to the function 'minFunc' that is used
    %for optimization, so see the documentation of minFunc for what the
    %parameters mean and do
    params.opts.algorithmic = touchField(params.opts.algorithmic,'VarInfX');
    params.opts.algorithmic.VarInfX = touchField(params.opts.algorithmic.VarInfX,'minFuncOptions');
    
    params.opts.algorithmic.VarInfX.minFuncOptions = touchField(params.opts.algorithmic.VarInfX.minFuncOptions,'display',	'none');
    params.opts.algorithmic.VarInfX.minFuncOptions = touchField(params.opts.algorithmic.VarInfX.minFuncOptions,'maxFunEvals',50000);
    params.opts.algorithmic.VarInfX.minFuncOptions = touchField(params.opts.algorithmic.VarInfX.minFuncOptions,'MaxIter',	5000);
    params.opts.algorithmic.VarInfX.minFuncOptions = touchField(params.opts.algorithmic.VarInfX.minFuncOptions,'progTol',	1e-6); % this might be too agressive, maybe 1e-9 is the better option
    params.opts.algorithmic.VarInfX.minFuncOptions = touchField(params.opts.algorithmic.VarInfX.minFuncOptions,'optTol',	1e-5);
    params.opts.algorithmic.VarInfX.minFuncOptions = touchField(params.opts.algorithmic.VarInfX.minFuncOptions,'Method',	'lbfgs');
    
    
    %%%% set parameters for MStep of observation model %%%%%%%%
    %these parameters are handed over to the function 'minFunc' that is used
    %for optimization, so see the documentation of minFunc for what the
    %parameters mean and do
    params.opts.algorithmic = touchField(params.opts.algorithmic,'MStepObservation');
    params.opts.algorithmic.MStepObservation = touchField(params.opts.algorithmic.MStepObservation,'minFuncOptions');
    
    params.opts.algorithmic.MStepObservation.minFuncOptions = touchField(params.opts.algorithmic.MStepObservation.minFuncOptions,'maxFunEvals', 5000);
    params.opts.algorithmic.MStepObservation.minFuncOptions = touchField(params.opts.algorithmic.MStepObservation.minFuncOptions,'MaxIter',	  500);
    params.opts.algorithmic.MStepObservation.minFuncOptions = touchField(params.opts.algorithmic.MStepObservation.minFuncOptions,'Method',	  'lbfgs');
    params.opts.algorithmic.MStepObservation.minFuncOptions = touchField(params.opts.algorithmic.MStepObservation.minFuncOptions,'progTol',     1e-9);
    params.opts.algorithmic.MStepObservation.minFuncOptions = touchField(params.opts.algorithmic.MStepObservation.minFuncOptions,'optTol',      1e-5);
    params.opts.algorithmic.MStepObservation.minFuncOptions = touchField(params.opts.algorithmic.MStepObservation.minFuncOptions,'display',	  'none');
    
    
    %%%% set parameters for EM iterations %%%%%%%%
    
    params.opts.algorithmic = touchField(params.opts.algorithmic,'TransformType','0');         % transform LDS parameters after each MStep to canonical form?
    params.opts.algorithmic = touchField(params.opts.algorithmic,'EMIterations');
    params.opts.algorithmic.EMIterations = touchField(params.opts.algorithmic.EMIterations,'maxIter',100);			% max no of EM iterations
    params.opts.algorithmic.EMIterations = touchField(params.opts.algorithmic.EMIterations,'maxCPUTime',inf);		% max CPU time for EM
    params.opts.algorithmic.EMIterations = touchField(params.opts.algorithmic.EMIterations,'progTolvarBound',1e-6);     	% progress tolerance on var bound per data time bin
    params.opts.algorithmic.EMIterations = touchField(params.opts.algorithmic.EMIterations,'abortDecresingVarBound',true);
    
    
    %%%% set parameters for initialization methods %%%%
    
    params.opts = touchField(params.opts,'initMethod','params');
    
    switch params.opts.initMethod
        
        case 'params'
            % do nothing
            
        case 'PLDSID'
            %use Poisson-Linear-Dynamics System Identification Method, see
            %documentation of 'FitPLDSParamsSSID' for details
            params.opts.algorithmic = touchField(params.opts.algorithmic,'PLDSID');
            params.opts.algorithmic.PLDSID = touchField(params.opts.algorithmic.PLDSID,'algo','SVD');
            params.opts.algorithmic.PLDSID = touchField(params.opts.algorithmic.PLDSID,'hS',xDim);
            params.opts.algorithmic.PLDSID = touchField(params.opts.algorithmic.PLDSID,'minFanoFactor',1.01);
            params.opts.algorithmic.PLDSID = touchField(params.opts.algorithmic.PLDSID,'minEig',1e-4);
            params.opts.algorithmic.PLDSID = touchField(params.opts.algorithmic.PLDSID,'useB',0);
            params.opts.algorithmic.PLDSID = touchField(params.opts.algorithmic.PLDSID,'doNonlinTransform',1);
            
        case 'ExpFamPCA'
            %use Exponential Family PCA, see function ExpFamPCA for details
            params.opts.algorithmic = touchField(params.opts.algorithmic,'ExpFamPCA');
            params.opts.algorithmic.ExpFamPCA = touchField(params.opts.algorithmic.ExpFamPCA,'dt',10);				% rebinning factor, choose such that roughly E[y_{kt}] = 1 forall k,t
            params.opts.algorithmic.ExpFamPCA = touchField(params.opts.algorithmic.ExpFamPCA,'lam',1);		  		% regularization coeff for ExpFamPCA
            params.opts.algorithmic.ExpFamPCA = touchField(params.opts.algorithmic.ExpFamPCA,'options');				% options for minFunc
            params.opts.algorithmic.ExpFamPCA.options = touchField(params.opts.algorithmic.ExpFamPCA.options,'display','none');
            params.opts.algorithmic.ExpFamPCA.options = touchField(params.opts.algorithmic.ExpFamPCA.options,'MaxIter',10000);
            params.opts.algorithmic.ExpFamPCA.options = touchField(params.opts.algorithmic.ExpFamPCA.options,'maxFunEvals',50000);
            params.opts.algorithmic.ExpFamPCA.options = touchField(params.opts.algorithmic.ExpFamPCA.options,'Method','lbfgs');
            params.opts.algorithmic.ExpFamPCA.options = touchField(params.opts.algorithmic.ExpFamPCA.options,'progTol',1e-9);
            params.opts.algorithmic.ExpFamPCA.options = touchField(params.opts.algorithmic.ExpFamPCA.options,'optTol',1e-5);
            
            
        case 'NucNormMin'
            %use Exponential Family PCA, see function MODnucnrmminWithd for details
            params.opts.algorithmic = touchField(params.opts.algorithmic,'NucNormMin');
            params.opts.algorithmic.NucNormMin = touchField(params.opts.algorithmic.NucNormMin,'dt',10);
            params.opts.algorithmic.NucNormMin = touchField(params.opts.algorithmic.NucNormMin,'fixedxDim',true);
            params.opts.algorithmic.NucNormMin = touchField(params.opts.algorithmic.NucNormMin,'options');
            params.opts.algorithmic.NucNormMin.options = touchField(params.opts.algorithmic.NucNormMin.options,'rho',	1.3);
            params.opts.algorithmic.NucNormMin.options = touchField(params.opts.algorithmic.NucNormMin.options,'eps_abs',	1e-6);
            params.opts.algorithmic.NucNormMin.options = touchField(params.opts.algorithmic.NucNormMin.options,'eps_rel',	1e-3);
            params.opts.algorithmic.NucNormMin.options = touchField(params.opts.algorithmic.NucNormMin.options,'maxIter',  	250);
            params.opts.algorithmic.NucNormMin.options = touchField(params.opts.algorithmic.NucNormMin.options,'nlin',     	'exp');
            params.opts.algorithmic.NucNormMin.options = touchField(params.opts.algorithmic.NucNormMin.options,'lambda',	0.03);
            params.opts.algorithmic.NucNormMin.options = touchField(params.opts.algorithmic.NucNormMin.options,'verbose',	0);
            
        otherwise
            
            warning('Unknown PLDS initialization method, cannot set parameters')
            
    end
    
end

function [seq varBound] = PLDSVariationalInference(params,seq)
    %
    % [seq] = PLDSVariationalInferenceDualLDS(params,seq)
    %
    
    Trials = numel(seq);
    
    for tr=1:Trials
        optparams.dualParams{tr} = [];
    end
    optparams.minFuncOptions = params.opts.algorithmic.VarInfX.minFuncOptions;
    
    [seq] = VariationalInferenceDualLDS(params,seq,optparams);
    
    varBound = 0;
    for tr=1:Trials; varBound = varBound + seq(tr).posterior.varBound; end;
end

function Lambda = buildPriorPrecisionMatrixFromLDS(params,T)
    %
    % Lambda = buildPrecisionMatrixFromLDS(params,T)
    %
    % construct the precision matrix of the prior across all time points and
    % dimensions, as described in Paninski et al, A new look at state-space
    % models for neural data, 2009
    %
    % c/o L Buesing and J Macke, 01/2014
    
    
    xDim   = size(params.model.A,1);
    invQ   = pinv(params.model.Q);
    invQ0  = pinv(params.model.Q0);
    AinvQ  = params.model.A'*invQ;
    AinvQA = AinvQ*params.model.A;
    
    
    Lambda = sparse(T*xDim,T*xDim);
    Lambda(1:xDim,1:xDim) = invQ0;
    
    for t=1:T-1
        xidx = ((t-1)*xDim+1):(t*xDim);
        Lambda(xidx,xidx) = Lambda(xidx,xidx)+AinvQA;
        Lambda(xidx,xidx+xDim) = -AinvQ;
        Lambda(xidx+xDim,xidx) = -AinvQ';
        Lambda(xidx+xDim,xidx+xDim) = Lambda(xidx+xDim,xidx+xDim)+invQ;
    end
    Lambda = sparse((Lambda+Lambda')/2);
end

function Mu = getPriorMeanLDS(params,T,varargin)
    %
    %
    %
    
    seq  = [];
    A    = params.model.A;
    x0   = params.model.x0;
    xAdd = [];
    
    assignopts(who,varargin);
    
    xDim = size(params.model.A,1);
    
    Mu = zeros(xDim,T);
    Mu(:,1) = x0;
    
    if params.model.notes.useB
        if ~isempty(seq)
            Mu(:,1) = Mu(:,1)+params.model.B*seq.u(:,1);
        else
            error('params.model.notes.useB == true   but no seq given!')
        end
    end
    
    if ~isempty(xAdd)
        Mu(:,1) = Mu(:,1)+xAdd(:,1);
    end
    
    
    for t=2:T
        Mu(:,t) = A*Mu(:,t-1);
        
        if ~isempty(seq) && params.model.notes.useB
            Mu(:,t) = Mu(:,t)+params.model.B*seq.u(:,t);
        end
        if ~isempty(xAdd)
            Mu(:,t) = Mu(:,t)+xAdd(:,t);
        end
        
    end
end

function [params, seq] = LDSApplyParamsTransformation(M,params,varargin)
    %
    % [params seq] = LDSApplyParamsTransformation(M,params,varargin)
    %
    % Applies M from left and inv(M)/M' from the right
    %
    % to A,Q,Q0,x0,B,C
    %
    %
    % L Buesing 2014
    
    seq = [];
    
    assignopts(who,varargin);
    
    xDim = size(params.model.A,1);
    
    if cond(M)>1e3
        warning('Attempting LDSApplyParamsTransformation with ill-conditioned transformation')
    end
    
    params.model.C  =     params.model.C  / M;
    params.model.A  = M * params.model.A  / M;
    if ~iscell(params.model.Q)
        params.model.Q  = M * params.model.Q  * M';
    else
        for mm=1:numel(params.model.Q)
            params.model.Q{mm}  = M * params.model.Q{mm}  * M';
        end
    end
    
    if ~iscell(params.model.Q0)
        params.model.Q0 = M * params.model.Q0 * M';
    else
        for mm=1:numel(params.model.Q0)
            params.model.Q0{mm} = M * params.model.Q0{mm} * M';
        end
    end
    params.model.x0 = M * params.model.x0;
    
    if isfield(params.model,'B')
        params.model.B = M*params.model.B;
    end
    
    if isfield(params.model,'Pi')
        if ~iscell(params.model.Q)
            params.model.Pi = dlyap(params.model.A,params.model.Q);
        else
            params.model.Pi = dlyap(params.model.A,params.model.Q{1});
        end
    end
    
    if ~isempty(seq)
        for tr=1:numel(seq)
            if isfield(seq,'posterior')
                seq(tr).posterior.xsm  = M*seq(tr).posterior.xsm;
                for t = 1:size(seq(tr).y,2);
                    xidx = ((t-1)*xDim+1):(xDim*t);
                    seq(tr).posterior.Vsm(xidx,:)  = M*seq(tr).posterior.Vsm(xidx,:)*M';
                    if t>1;
                        seq(tr).posterior.VVsm(xidx-xDim,:) = M*seq(tr).posterior.VVsm(xidx-xDim,:)*M';
                    end
                end
            end
            if isfield(seq,'x')
                seq(tr).x = M*seq(tr).x;
            end
        end
    end
end

function penalty = LDSemptyParamPenalizerHandle(params)
    %
    % function penalty = LDSemptyParamPenalizerHandle(params)
    %
    
    penalty = 0;
end

function params = LDSgenerateExample(varargin)
    %
    % params = generateLDS(varargin)
    %
    % make a random LDS given some parameters
    
    %NOT DOCUMENTED YET%
    
    uDim     = 0;
    xDim     = 10;
    yDim     = 100;
    
    Arot     = 0.1;
    Aspec    = 0.99;
    Arand    = 0.03;
    Q0max    = 0.3;
    Rmin     = 0.1;
    Rmax     = 0.1;
    
    assignopts(who,varargin);
    
    
    %%%%%%%%%  generate parameters %%%%%%%%%
    
    A  = eye(xDim)+Arand*randn(xDim);
    A  = A./max(abs(eig(A)))*Aspec;
    MAS = randn(xDim); MAS = (MAS-MAS')/2;A  = expm(Arot.*(MAS))*A;
    Q  = diag(rand(xDim,1));
    Q0 = dlyap(A,Q);
    M  = diag(1./sqrt(diag(Q0)));
    A  = M*A*pinv(M);
    Q  = M*Q*M'; Q=(Q+Q)/2;
    
    O  = orth(randn(xDim));
    Q0 = O*diag(rand(xDim,1)*Q0max)*O'/3;
    x0 = randn(xDim,1)/3;
    
    C  = randn(yDim,xDim)./sqrt(3*xDim);
    R  = diag(rand(yDim,1)*Rmax+Rmin);
    d  = 0.3*randn(yDim,1);
    
    params.model.A    = A;
    params.model.Q    = Q;
    params.model.Q0   = Q0;
    params.model.x0   = x0;
    params.model.C    = C;
    params.model.d    = d;
    params.model.R    = R;
    params.model.Pi   = dlyap(params.model.A,params.model.Q);
    params.model.notes.useR = true;
    params.model.notes.useS = false;
    
    if uDim>0
        cQ = max(abs(diag(chol(params.model.Q))));
        params.model.B    = cQ*(rand(xDim,uDim)+0.5)/(uDim);
        params.model.notes.useB = true;
    else
        params.model.notes.useB = false;
    end
    
    params = LDSsetDefaultParameters(params,xDim,yDim);
    
end

function [seq varBound Lambda LambdaPost] = LDSInference(params,seq)
    %
    % simplest Kalman smoother in O(T)
    %
    % assume all trials are of the same length T
    %
    
    
    
    [yDim,xDim] = size(params.model.C);
    
    Trials      = numel(seq);
    varBound    = 0;
    
    for tr=1:Trials
        
        
        T           = size(seq(tr).y,2);
        
        
        %%%%%%%%%%%%% covariances %%%%%%%%%%%%
        
        Lambda     = buildPriorPrecisionMatrixFromLDS(params,T); % prior precision
        LambdaPost = Lambda;                                 % posterior precision
        
        CRC    = zeros(xDim*T,xDim);
        CRinvC = params.model.C'*pinv(params.model.R)*params.model.C;
        for t=1:T
            xidx = ((t-1)*xDim+1):(t*xDim);
            CRC(xidx,:) = CRinvC;
            LambdaPost(xidx,xidx)  = LambdaPost(xidx,xidx) + CRinvC;
        end
        
        [Vsm, VVsm] = smoothedKalmanMatrices(params.model,CRC);
        
        %%%%%%%%%%%%%%%%% means %%%%%%%%%%%%%%%%
        
        varBound = varBound + -0.5*T*yDim;
        
        Mu = getPriorMeanLDS(params,T);
        LamMu = Lambda*vec(Mu);
        
        
        
        Y = bsxfun(@minus,seq(tr).y,params.model.d);
        if params.model.notes.useS
            Y = Y-seq.s;
        end
        Yraw = Y;
        Y  = params.model.R\Y;
        Y  = params.model.C'*Y;
        seq(tr).posterior.xsm  = reshape(LambdaPost\(LamMu+vec(Y)),xDim,T);
        seq(tr).posterior.Vsm  = Vsm;
        seq(tr).posterior.VVsm = VVsm;
        
        Yraw = Yraw-params.model.C*Mu;
        varBound = varBound - 0.5*trace(params.model.R\(Yraw*Yraw'));
        
        Yraw = vec(params.model.C'*(params.model.R\Yraw));
        varBound = varBound + 0.5*Yraw'*(LambdaPost\Yraw);
        
    end
    
    varBound = varBound - 0.5*Trials*(logdet(LambdaPost,'chol')-logdet(Lambda,'chol')+T*logdet(params.model.R));
end

function params = LDSInitialize(seq,xDim,initMethod,params)
    %
    % function params = LDSInitialize(params,seq)
    %
    % inititalize parameters of LDS model
    % the moment, focusses on exponential link function and Poisson
    % observations.
    %
    %input:
    % seq:      standard data struct
    % xdim:     desired latent dimensionality
    %
    % initMethods:
    %
    % - params											% just initialize minimal undefiened fields with standard values
    % - SSID											% moment conversion + Ho-Kalman SSID
    % - PCA
    %
    % (c) L Buesing 2014
    
    
    yDim       = size(seq(1).y,1);
    Trials     = numel(seq);
    params.opts.initMethod = initMethod;
    params     = LDSsetDefaultParameters(params,xDim,yDim);					% set standard parameter values
    
    
    switch initMethod
        
        case 'params'
            % do nothing
            disp('Initializing PLDS parameters with given parameters')
            
            
        case 'SSID'
            % !!! debugg SSID stuff separately & change to params.model convention
            disp('Initializing LDS parameters using SSID')
            if params.model.notes.useB
                warning('SSID initialization with external input: not implemented yet!!!')
            end
            SSIDoptions  = struct2arglist(params.opts.algorithmic.SSID);
            params.model = FitPLDSParamsSSID(seq,xDim,'params',params.model,SSIDoptions{:});
            
            
        otherwise
            warning('Unknown PLDS initialization method')
            
    end
    
    if params.model.notes.useB && (numel(params.model.B)<1)
        params.model.B = zeros(xDim,size(seq(1).u,1));
    end
    
    params = LDSTransformParams(params,'TransformType',params.opts.algorithmic.TransformType);	% clean up parameters
    params.modelInit = params.model;
    
end

function [params seq] = LDSMStep(params,seq)
    %
    % params = LDSMStep(params,seq)
    %
    
    
    params = LDSMStepLDS(params,seq);
    params = LDSMStepObservation(params,seq);
    
    params = LDSTransformParams(params,'TransformType',params.opts.algorithmic.TransformType);
end

function params = LDSMStepLDS(params,seq)
    %
    % function params = MStepLDS(params,seq)
    %
    % Parameters to update: A,Q,Q0,x0,B
    %
    %
    % (c) L Buesing 2014
    %
    
    xDim    = size(params.model.A,1);
    Trials  = numel(seq);
    
    
    %% compute posterior statistics
    
    if params.model.notes.useB
        uDim = size(seq(1).u,1);
    else
        uDim = 0;
    end
    
    S11 = zeros(xDim,xDim);
    S01 = zeros(xDim+uDim,xDim);
    S00 = zeros(xDim+uDim,xDim+uDim);
    
    x0 = zeros(xDim,Trials);
    Q0 = zeros(xDim,xDim);
    
    Tall = [];
    
    for tr = 1:Trials
        
        T = size(seq(tr).y,2);
        Tall  = [Tall T];
        
        if isfield(seq(tr).posterior,'Vsm')
            Vsm   = reshape(seq(tr).posterior.Vsm' ,xDim,xDim,T);
            VVsm  = reshape(seq(tr).posterior.VVsm',xDim,xDim,T-1);
        else
            Vsm   = reshape(seq(1).posterior.Vsm' ,xDim,xDim,T);
            VVsm  = reshape(seq(1).posterior.VVsm',xDim,xDim,T-1);
        end
        
        MUsm0 = seq(tr).posterior.xsm(:,1:T-1);
        MUsm1 = seq(tr).posterior.xsm(:,2:T);
        
        S11                = S11                + sum(Vsm(:,:,2:T),3)  + MUsm1*MUsm1';
        S01(1:xDim,:)      = S01(1:xDim,:)      + sum(VVsm(:,:,1:T-1),3) + MUsm0*MUsm1';
        S00(1:xDim,1:xDim) = S00(1:xDim,1:xDim) + sum(Vsm(:,:,1:T-1),3)  + MUsm0*MUsm0';
        
        if params.model.notes.useB
            u = seq(tr).u(:,1:T-1);
            S01(1+xDim:end,:)          = S01(1+xDim:end,:)          + u*MUsm1';
            S00(1+xDim:end,1:xDim)     = S00(1+xDim:end,1:xDim)     + u*MUsm0';
            S00(1:xDim,1+xDim:end)     = S00(1:xDim,1+xDim:end)     + MUsm0*u';
            S00(1+xDim:end,1+xDim:end) = S00(1+xDim:end,1+xDim:end) + u*u';
        end
        
        x0(:,tr) = MUsm0(:,1);
        Q0 = Q0 + Vsm(:,:,1);
        
    end
    
    S00 = (S00+S00')/2;
    S11 = (S11+S11')/2;
    
    if params.model.notes.learnA
        params.model.A  = S01'/S00;
    end
    params.model.Q  = (S11+params.model.A*S00*params.model.A'-S01'*params.model.A'-params.model.A*S01)./(sum(Tall)-Trials);
    %params.model.Q  = S11-S01'/S00*S01;
    params.model.Q  = (params.model.Q+params.model.Q')/2;
    
    %[aQ bQ] = eig(params.model.Q);
    %params.model.Q = aQ*diag(max(diag(bQ),0))*aQ';
    
    if params.model.notes.useB
        params.model.B = params.model.A(:,1+xDim:end);
        params.model.A = params.model.A(:,1:xDim);
    end
    
    if params.model.notes.learnx0
        params.model.x0 = mean(x0,2);
    end
    
    x0dev = bsxfun(@minus,x0,params.model.x0);
    
    if params.model.notes.learnQ0
        params.model.Q0 = (Q0 + x0dev*x0dev')./Trials;
    else
        params.model.Q0 = dlyap(params.model.A,params.model.Q);
    end
    
    if (min(eig(params.model.Q))<0) || (min(eig(params.model.Q0))<0)
        keyboard
        params.model.Q=params.model.Q+1e-9*eye(xDim);
    end
end


function params = LDSMStepObservation(params,seq)
    %
    % function params = LDSMStepObservation(params,seq)
    %
    
    [yDim xDim] = size(params.model.C);
    Trials = numel(seq);
    Tall   = sum([seq.T]);
    
    if params.model.notes.useCMask;
        warning('params.model.notes.useCMask == true: not implemented for LDS yet')
    end
    
    %%% optimization %%%
    
    Psi  = []; for tr=1:Trials; Psi = [Psi seq(tr).posterior.xsm]; end
    Psi  = [Psi;ones(1,Tall)];
    
    Yall = [seq.y];
    if params.model.notes.useS; Yall = Yall-[seq.s]; end
    YPsi = Yall*Psi';
    
    PsiPsi = Psi*Psi';
    VsmAll = zeros(xDim);
    for tr=1:Trials;  VsmAll = VsmAll + sum(reshape(seq(tr).posterior.Vsm',xDim,xDim,seq(tr).T),3);end
    PsiPsi(1:end-1,1:end-1) = PsiPsi(1:end-1,1:end-1) + VsmAll;
    
    CdOpt = YPsi/(PsiPsi+eye(xDim+1)*1e-10);
    params.model.C = CdOpt(:,1:end-1);
    params.model.d = CdOpt(:,end);
    params.model.R = cov((Yall - CdOpt*Psi)',1)+params.model.C*VsmAll*params.model.C'./sum(Tall);
    
    if params.model.notes.diagR
        params.model.R = diag(diag(params.model.R));
    end
end

function model = LDSObservedEstimation(X,model,dt,u)
    %
    % estimate paramaters of a fully observed LDS
    %
    % Input:
    %
    %   - X: observations   xDim x T
    %   - model: struct to save parameters to
    %   - dt: sub-sampling factor
    %   - u: external input (optional, only used if model.notes.useB == true)
    %
    % Ouput:
    %
    %   - model: standard LDS model structure
    %
    %
    % Lars Buesing, 2014
    
    xDim = size(X,1);
    T    = size(X,2);
    Pi   = X*X'./T;
    
    if ~model.notes.useB
        
        A  = X(:,2:end)/X(:,1:end-1);         % find A via regression
        if dt>1
            %A = diag(min(max((diag(abs(A))).^(1/dt),0),1));
            % clearly, this needs fixing !!!
            A = diag(min(max(diag(abs(A)),0.1),1));
        end
        Q  = Pi-A*Pi*A';                      % residual covariance
        
    else
        
        uDim = size(u,1);
        if size(u,2)>size(X,2)
            uSub = subsampleSignal(u,dt);
        else
            uSub = u;
        end
        AB = X(:,2:end)/[X(:,1:end-1);uSub(:,2:end)];
        A  = AB(1:xDim,1:xDim);
        B  = AB(1:xDim,1+xDim:end);
        
        if dt>1
            A = diag(min(max((diag(abs(A))).^(1/dt),0),1));
            
            uauto = zeros(uDim,dt);
            for ud=1:uDim
                uauto(ud,:) = autocorr(u(ud,:),dt-1);
            end
            Aauto = zeros(xDim,dt);
            for tt=1:dt
                Aauto(:,tt) = diag(A).^(tt-1);
            end
            M = Aauto*uauto'*diag(1./uauto(:,1));
            B = B./M;
        end
        
        xuCov = A*X(:,1:end-1)*uSub(:,2:end)'/(size(uSub,2)-1)*B';
        Q = Pi-A*Pi*A'-B*cov(uSub')*B'-xuCov-xuCov';
        model.B = B;
        
    end
    
    [Uq Sq Vq] = svd(Q);                   % ensure that Q is pos def
    Q  = Uq*diag(max(diag(Sq),0))*Uq';
    x0 = zeros(xDim,1);
    Q0 = dlyap(A,Q);
    
    
    model.A  = A;
    model.Pi = Pi;
    model.Q  = Q;
    model.Q0 = Q0;
    model.x0 = x0;
    
end

function seq = LDSsample(params,T,Trials,varargin)
    %
    % sampleLDS(params,T,Trials)
    %
    % sample from linear dynamical system model
    %
    %
    % (c) L Buesing 2014
    %
    
    u = [];
    s = [];
    assignopts(who,varargin);
    
    if numel(T)==1
        T = ones(Trials,1)*T;
    end
    
    Trials = numel(T);
    
    [yDim xDim] = size(params.model.C);
    CQ          = chol(params.model.Q);
    CQ0         = chol(params.model.Q0);
    
    if isfield(params.model,'R') && params.model.notes.useR
        R = params.model.R;
        if size(R,2)==1
            R = diag(R);
        end
        CR = chol(R);
    else
        CR = zeros(yDim);
    end
    
    
    for tr=1:Trials
        
        if params.model.notes.useB
            % !!! take this out
            if isempty(u)
                %      warning('You are using an extremely dirty convenience function. It will disappear')
                uDim = size(params.model.B,2);
                gpsamp = sampleGPPrior(1,T(tr),uDim-1,'tau',10);
                tpsamp = (vec(repmat(rand(1,floor(T(tr)/10))>0.5,10,1))-0.5); tpsamp = [tpsamp' zeros(1,T(tr)-floor(T(tr)/10)*10)];
                seq(tr).u = [gpsamp{1}/3;tpsamp];
            else
                seq(tr).u = u{tr};
            end
        end
        
        seq(tr).x = zeros(xDim,T(tr));
        seq(tr).x(:,1) = params.model.x0+CQ0'*randn(xDim,1);
        if params.model.notes.useB; seq(tr).x(:,1) = seq(tr).x(:,1)+params.model.B*seq(tr).u(:,1);end;
        for t=2:T(tr)
            seq(tr).x(:,t) = params.model.A*seq(tr).x(:,t-1)+CQ'*randn(xDim,1);
            if params.model.notes.useB; seq(tr).x(:,t) = seq(tr).x(:,t)+params.model.B*seq(tr).u(:,t);end;
        end
        seq(tr).y = bsxfun(@plus,params.model.C*seq(tr).x,params.model.d)+CR'*randn(yDim,T(tr));
        seq(tr).T = T(tr);
        
        if params.model.notes.useS
            seq(tr).y = seq(tr).y+s{tr};
            seq(tr).s = s{tr};
        end
        
    end
    
end

function params = LDSsetDefaultParameters(params,xDim,yDim)
    %
    % params = LDSsetDefaultParameters(params,xDim,yDim)
    %
    %
    % Lars Buesing, Jakob H Macke
    %
    
    
    %%%%%%%%%%% set standard parameters %%%%%%%%%%%%%%%%%%%%
    
    params = touchField(params,'model');
    params.model = touchField(params.model,'A',0.9*eye(xDim));        % dynamics matrix A
    params.model = touchField(params.model,'Q',(1-0.9.^2)*eye(xDim)); %innovation covariance Q
    params.model = touchField(params.model,'Q0',eye(xDim));           % initial state covariance Q0
    params.model = touchField(params.model,'x0',zeros(xDim,1));       %initial mean x0
    params.model = touchField(params.model,'C',randn(yDim,xDim)./sqrt(xDim)); %couplings from latent to observed C
    params.model = touchField(params.model,'d',zeros(yDim,1)-2.0);    %mean-controlling offset d for each neuron
    params.model = touchField(params.model,'B',zeros(xDim,0));
    params.model = touchField(params.model,'R',0.1*eye(yDim));        % private noise covariance matrix
    
    params.model = touchField(params.model,'notes');
    params.model.notes = touchField(params.model.notes,'learnx0', true);
    params.model.notes = touchField(params.model.notes,'learnQ0', true);
    params.model.notes = touchField(params.model.notes,'learnA',  true);
    params.model.notes = touchField(params.model.notes,'learnR',  true);
    params.model.notes = touchField(params.model.notes,'diagR',   true);  % learn diagonal private variances?
    params.model.notes = touchField(params.model.notes,'useR',    true);
    params.model.notes = touchField(params.model.notes,'useB',    false);
    params.model.notes = touchField(params.model.notes,'useS',    false);
    params.model.notes = touchField(params.model.notes,'useCMask',false);
    
    
    %%%%%%%%%%% set standard observation model handles for variational inference %%%%%%%%%%%%%%%%%%%%
    
    params.model = touchField(params.model,'inferenceHandle',     @LDSInference);
    params.model = touchField(params.model,'MStepHandle',         @LDSMStep);
    params.model = touchField(params.model,'ParamPenalizerHandle',@LDSemptyParamPenalizerHandle);
    
    
    %%%%%%%%%%% set standard algorithmic parameters %%%%%%%%%%%%%%%%%%%%
    
    params = touchField(params,'opts');
    params.opts = touchField(params.opts,'algorithmic');
    
    
    %%%% set parameters for MStep of observation model %%%%%%%%
    %these parameters are handed over to the function 'minFunc' that is used
    %for optimization, so see the documentation of minFunc for what the parameters mean and do
    params.opts.algorithmic = touchField(params.opts.algorithmic,'MStepObservation');
    
    
    %%%% set parameters for EM iterations %%%%%%%%
    
    params.opts.algorithmic = touchField(params.opts.algorithmic,'TransformType','0');         % transform LDS parameters after each MStep to canonical form?
    params.opts.algorithmic = touchField(params.opts.algorithmic,'EMIterations');
    params.opts.algorithmic.EMIterations = touchField(params.opts.algorithmic.EMIterations,'maxIter',100);			% max no of EM iterations
    params.opts.algorithmic.EMIterations = touchField(params.opts.algorithmic.EMIterations,'maxCPUTime',inf);		% max CPU time for EM
    params.opts.algorithmic.EMIterations = touchField(params.opts.algorithmic.EMIterations,'progTolvarBound',1e-6);     	% progress tolerance on var bound per data time bin
    params.opts.algorithmic.EMIterations = touchField(params.opts.algorithmic.EMIterations,'abortDecresingVarBound',true);
    
    
    %%%% set parameters for initialization methods %%%%
    
    params.opts = touchField(params.opts,'initMethod','params');
    
    switch params.opts.initMethod
        
        case 'params'
            % do nothing
            
        case 'SSID'
            %Subspace System Identification Method
            params.opts.algorithmic = touchField(params.opts.algorithmic,'SSID');
            params.opts.algorithmic.SSID = touchField(params.opts.algorithmic.SSID,'algo','SVD');
            params.opts.algorithmic.SSID = touchField(params.opts.algorithmic.SSID,'hS',xDim);
            params.opts.algorithmic.SSID = touchField(params.opts.algorithmic.SSID,'useB',0);
            params.opts.algorithmic.SSID = touchField(params.opts.algorithmic.SSID,'doNonlinTransform',0);
            
        otherwise
            warning('Unknown LDS initialization method, cannot set parameters')
            
    end
end

function [params, seq] = LDSTransformParams(params,varargin)
    %
    % function [params, seq] = LDSTransformParams(params,varargin)
    %
    %
    % transform parameters of LDS by imposing constraints on C and the
    % stationary distribution Pi.
    %
    %  Pi := dlyap(params.model.A,params.model.Q)   stationary distribution
    %
    % TransformType:
    %
    % 0: do nothing
    % 1: C'*C = eye		&&		Pi = diag    [default]
    % 2: C'*C = diag 	&& 		Pi = eye
    % 3: ||C(:,k)||_2 = 1
    % 4: P_ii = 1
    % 5: A = blk-diag       &&              ||C(:,k)||_2 = 1
    %
    % (c) 2014 Lars Busing   lars@stat.columbia.edu
    %
    %
    % also see: LDSApplyParamsTransformation(M,params,varargin)
    %
    
    
    seq = [];
    TransformType   = '1';
    
    assignopts(who,varargin);
    
    xDim = size(params.model.A,1);
    
    switch TransformType
        
        case '0'
            
            % do nothing
            
        case '1'
            
            [UC,SC,VC]    = svd(params.model.C,0);
            [params seq]  = LDSApplyParamsTransformation(SC*VC',params,'seq',seq);
            
            params.model.Pi     = dlyap(params.model.A,params.model.Q);
            if min(eig(params.model.Pi))<0
                params.model.Pi = params.model.Q;
            end
            [UPi SPi VPi] = svd(params.model.Pi);
            [params seq]  = LDSApplyParamsTransformation(UPi',params,'seq',seq);
            
        case '2'
            
            params.model.Pi = dlyap(params.model.A,params.model.Q);
            if min(eig(params.model.Pi))<0
                params.model.Pi = params.model.Q;
            end
            [UPi SPi VPi] = svd(params.model.Pi);
            M    	      = diag(1./sqrt(diag(SPi)))*UPi';
            [params seq]  = LDSApplyParamsTransformation(M,params,'seq',seq);
            
            [UC,SC,VC]    = svd(params.model.C,0);
            [params seq]  = LDSApplyParamsTransformation(VC',params,'seq',seq);
            
        case '3'
            
            [params seq] = LDSApplyParamsTransformation(diag(sqrt(sum(params.model.C.^2,1))),params,'seq',seq);
            
        case '4'
            
            Pi = dlyap(params.model.A,params.model.Q);
            if min(eig(Pi))<0
                Pi = params.model.Q;
            end
            M = diag(1./sqrt(diag(Pi)));
            [params seq] = LDSApplyParamsTransformation(M,params,'seq',seq);
            
        case '5'
            
            [T B] = bdschur(params.model.A,inf);
            [params seq] = LDSApplyParamsTransformation(pinv(T),params);
            [params seq] = LDSApplyParamsTransformation(diag(sqrt(sum(params.model.C.^2,1))),params,'seq',seq);
            
            
        otherwise
            
            warning('Unknow paramter transformation type')
            
    end
end

function [Vsm VVsm F0 F1] = smoothedKalmanMatrices(model,CRinvC)
    %
    % [Vsm VVsm F0 F1] = smoothedKalmanMatrices(model,CRinvC)
    %
    %
    % computes posterior covariances by Kalman smoothing
    %
    % INPUT:
    %
    %  - CRinvC is a array of dimension (xDim*T) x (xDim) where
    %    CRinvC((t-1)*xDim+1:t*xDim) = C'/Rt*C where Rt is observation noise
    %    covariance matrix at time t and C is loading matrix
    %
    %  - model with fields, also works if first argument is params
    %    - A
    %    - Q
    %    - Q0
    %
    %
    % OUTPUT: all matrices but VVsm are of dimension (xDim*T) x (xDim)
    %
    % the t-th block of dimension (xDim) x (xDim) is defined as:
    % F0   = predictive cov Cov(x_t|y_{1:t-1})
    % F1   = filtering cov  Cov(x_t|y_{1:t})
    % Vsm  = smoothed cov   Cov(x_t|y_{1:T})
    %
    % VVsm is of dimension (xDim*(T-1)) x (xDim)
    %
    % VVsm = smoothed cov   Cov(x_{t+1},x_t}|y_{1:T})
    %
    %
    %
    % (c) Lars Buesing 2013,2014
    %
    
    
    if isfield(model,'model');   % check if first input is params or model
        model = model.model;
    end
    
    xDim = size(model.A,1);
    T    = round(size(CRinvC,1)/xDim);
    Vsm  = zeros(xDim*T,xDim);
    VVsm = zeros(xDim*(T-1),xDim);
    F0   = zeros(xDim*T,xDim);
    F1   = zeros(xDim*T,xDim);
    
    
    % forward pass
    
    F0(1:xDim,1:xDim) = model.Q0;
    for t=1:(T-1)
        xidx = ((t-1)*xDim+1):(t*xDim);
        %F1(xidx,:) = pinv(eye(xDim)+F0(xidx,:)*CRinvC(xidx,:))*F0(xidx,:);  % debug line, do not use
        F1(xidx,:) = (eye(xDim)+F0(xidx,:)*CRinvC(xidx,:))\F0(xidx,:);
        F0(xidx+xDim,:) = model.A*F1(xidx,:)*model.A'+model.Q;
    end
    t=T;xidx = ((t-1)*xDim+1):(t*xDim);
    F1(xidx,:) = eye(xDim)/(eye(xDim)/(F0(xidx,:))+CRinvC(xidx,:));
    
    
    % backward pass using Rauch???Tung???Striebel smoother
    
    t = T; xidx = ((t-1)*xDim+1):(t*xDim);
    Vsm(xidx,:) = F1(xidx,:);
    for t=(T-1):(-1):1
        xidx = ((t-1)*xDim+1):(t*xDim);
        Ck = F1(xidx,:)*model.A'/F0(xidx+xDim,:);
        Vsm(xidx,:)  = F1(xidx,:)+Ck*(Vsm(xidx+xDim,:)-F0(xidx+xDim,:))*Ck';
        VVsm(xidx,:) = (Ck*Vsm(xidx+xDim,:))';
    end
end

function domFlag = ExpPoissonDomain(lam);
    %
    
    if (min(lam)<0)||any(isnan(lam))||any(isinf(lam))
        domFlag = false;
    else
        domFlag = true;
    end
end

function [f df] = ExpPoissonDualHandle(lam,varargin);
    %
    % dual of standard exp-poisson likelihood
    %
    
    if min(lam)<0
        f  = inf;
        df = nan(size(lam));
    else
        loglam = log(lam);
        f  = lam'*(loglam-1);
        df = loglam;
    end
end

function [f] = ExpPoissonHandle(y,OverM_ast,OverV_ast,varargin);
    %
    % NEGATIVE exp-poisson likelihood without base measure
    %
    
    f = -y.*OverM_ast+exp(OverM_ast+OverV_ast/2);
end

function [f df] = ExpPoissonMixDualHandle(lam,logPi);
    %
    % dual of standard exp-poisson likelihood
    %
    
    if min(lam)<0
        f  = inf;
        df = nan(size(lam));
    else
        loglam = log(lam);
        f  = lam'*(loglam-1-logPi);
        df = loglam-logPi;
    end
end

function [f] = ExpPoissonMixHandle(y,OverM_ast,OverV_ast,logPi);
    %
    % NEGATIVE exp-poisson likelihood
    %
    
    [yDim T] = size(y);
    
    Pi = reshape(exp(logPi),yDim,T);
    f  = -y.*OverM_ast+exp(OverM_ast+OverV_ast/2).*Pi;
end

function [f df] = LogisticBernoulliDualHandle(lam,varargin);
    %
    % dual of standard exp-poisson likelihood
    %
    
    if (min(lam)<0)||(max(lam)>1)
        f  = inf;
        df = nan(size(lam));
    else
        loglam = log(lam);
        f  = lam'*loglam+(1-lam)'*(log(1-lam));
        df = log(lam./(1-lam));
    end
    
end

function [f] = PoissonBaseMeasure(y,params);
    %
    % log y!
    %
    
    f = -sum(log(gamma(vec(y)+1)));
    
end

function [f, df, varBound, m_ast, invV_ast, Vsm, VVsm, over_m, over_v] = VariationalInferenceDualCost(lam,VarInfparams);
    %
    % [f, df, varBound, m_ast, invV_ast, Vsm, VVsm, over_m, over_v] = VariationalInferenceDualCost(lam,VarInfparams)
    %
    % Cost function for variational inference via dual optimization for
    % Gaussian LDS with exponential family observations
    %
    % see [M. E. Khan, A. Aravkin, M. Friedlander, and M. Seeger. Fast Dual Variational Inference for Non-Conjugate Latent Gaussian Models. In JMLR W&CP, volume 28, pages 951-959, 2013]
    %
    % VarInfparams.Lambda
    % VarInfparams.y
    % VarInfparams.mu
    % VarInfparams.W
    % VarInfparams.A
    % VarInfparams.WlamW
    % VarInfparams.d
    % VarInfparams.CC
    %
    % OUTPUT:
    % f        = dual cost
    % df       = gradient of dual cost wrt lam
    % varBound = variational lower bound to marignal log-likelihood log p(y)
    % m_ast    = variational posterior mean xDim x T
    % invV_ast = variational posterior precision
    % Vsm	   = var smoothed cov Cov[x(t),x(t)|y_{1:T}]
    % VVsm	   = var smoothed cov Cov[x(t+1),x(t)|y_{1:T}]
    % over_m   = W*m_ast+d
    % over_v   = diag(W/invV_ast*W)
    %
    % (c) Lars Buesing, 2013
    %
    
    xDim     = size(VarInfparams.A,1);
    y        = VarInfparams.y; % observed data
    [yDim T] = size(y);
    
    % check domain of lam
    if ~feval(VarInfparams.domainHandle,lam); f = inf; df = nan(size(lam)); return; end
    
    W  = VarInfparams.W;          		     	     % loading matrix, sparse, only blk-diag
    mu = VarInfparams.mu;         			     % prior mean
    if isfield(VarInfparams,'d')
        d = VarInfparams.d;		                     % bias for likelihood
    else
        d = zeros(yDim*T,1);
    end
    Lambda = VarInfparams.Lambda;			     % precision matrix of LDS prior, assumed to be tri-diagonal
    
    % compute quadratric term of cost function  (lam-y)'*W*Sig*W'*(lam-y)
    
    lamY          = lam - vec(VarInfparams.y);
    WlamY         = W'*lamY;
    SigWlamY      = Lambda\WlamY;
    WlamYSigWlamY = WlamY'*SigWlamY;
    
    
    %VarInfparams.WlamW  = sparse(zeros(xDim*T)); %debug-line, do note use, already pre-allocated
    for t=1:T
        xidx = ((t-1)*xDim+1):(t*xDim);
        yidx = ((t-1)*yDim+1):(t*yDim);
        %VarInfparams.WlamW(xidx,xidx) = VarInfparams.C'*diag(lam(yidx))*VarInfparams.C; %debug-line, use below
        VarInfparams.WlamW(xidx,xidx) = reshape(VarInfparams.CC*lam(yidx),xDim,xDim);
    end
    Alam       = Lambda+VarInfparams.WlamW;  % precision matrix of current variational approximation
    logdetAlam = logdet(Alam,'chol');
    
    % function value
    
    [like_f, like_df] = feval(VarInfparams.dualHandle,lam,VarInfparams.dualParams);
    
    f = 0.5*WlamYSigWlamY-mu'*WlamY-d'*lamY-0.5*logdetAlam+like_f;
    
    % catch infeasible lambdas
    if isinf(f)
        f = inf; df = nan(size(lam));
        varBound = -inf;
        m_ast = nan(xDim*T,1); invV_ast = nan;
        Vsm = nan(xDim*T,xDim); VVsm = nan(xDim*T,xDim);
        over_m = nan(yDim*T,1); over_v = nan(yDim*T,1);
        return
    end
    
    % gradient
    
    CRinvC = zeros(xDim*T,xDim);
    for t=1:T
        xidx = ((t-1)*xDim+1):(t*xDim);
        CRinvC(xidx,:) = VarInfparams.WlamW(xidx,xidx);
    end
    
    [Vsm, VVsm] = smoothedKalmanMatrices(VarInfparams,CRinvC);
    lam_con = zeros(yDim*T,1); % equals diag(W/Alam*W')
    for t=1:T
        xidx = ((t-1)*xDim+1):(t*xDim);
        yidx = ((t-1)*yDim+1):(t*yDim);
        %lam_con(yidx) = diag(VarInfparams.C*P(xidx,:)*VarInfparams.C');  %debug-line, use below
        lam_con(yidx) = VarInfparams.CC'*vec(Vsm(xidx,:));
    end
    
    %df = W*SigWlamY-W*mu+loglam-0.5*diag(W/Alam*W'); %debug-line, use below
    df  = W*SigWlamY-W*mu-d-lam_con/2+like_df;
    
    
    %%%%% compute variational lower bound
    
    m_ast    = mu-SigWlamY;  % optimal mean
    invV_ast = Alam;         % optimal inverse covariance
    
    over_m 	 = W*m_ast+d;
    over_v	 = lam_con;
    over_m 	 = reshape(over_m,yDim,T);
    over_v   = reshape(over_v,yDim,T);
    
    varBound = -0.5*logdetAlam-0.5*WlamYSigWlamY+0.5*lam'*lam_con;                                            %prior contribution
    varBound = varBound - 0.5*logdet(VarInfparams.Q)*(T-1)-0.5*logdet(VarInfparams.Q0);
    
    varBound = varBound-sum(vec(feval(VarInfparams.likeHandle,y,over_m,over_v,VarInfparams.dualParams)));     %likelihood contribution
    if isfield(VarInfparams,'DataBaseMeasure');   % !!! put this into PLDSInitialize or smth
        varBound = varBound+VarInfparams.DataBaseMeasure;
    end
    
end

function [seq] = VariationalInferenceDualLDS(params,seq,optparams)
    %
    % [seq] = VariationalInferenceDualLDS(params,seq);
    %
    %
    % L Buesing 2014
    %
    
    
    Trials      = numel(seq);
    [yDim xDim] = size(params.model.C);
    Tmax        = max([seq.T]);
    
    % set up parameters for variational inference
    
    VarInfparams    = params.model;
    VarInfparams.CC = zeros(xDim,xDim,yDim);
    for yy=1:yDim
        VarInfparams.CC(:,:,yy) = params.model.C(yy,:)'*params.model.C(yy,:);
    end
    VarInfparams.CC = reshape(VarInfparams.CC,xDim^2,yDim);
    
    Cl = {}; for t=1:Tmax; Cl = {Cl{:} params.model.C}; end
    Wmax = sparse(blkdiag(Cl{:}));
    
    % iterate over trials
    
    for tr = 1:Trials
        
        T = size(seq(tr).y,2);
        
        VarInfparams.d          = repmat(params.model.d,T,1); if params.model.notes.useS; VarInfparams.d = VarInfparams.d + vec(seq(tr).s); end
        VarInfparams.mu         = vec(getPriorMeanLDS(params,T,'seq',seq(tr)));
        VarInfparams.W          = Wmax(1:yDim*T,1:xDim*T);
        VarInfparams.y          = seq(tr).y;
        VarInfparams.Lambda     = buildPriorPrecisionMatrixFromLDS(params,T);  % generate prior precision matrix
        VarInfparams.WlamW      = sparse(zeros(xDim*T)); %allocate sparse observation matrix
        VarInfparams.dualParams = optparams.dualParams{tr};
        
        if isfield(params.model,'baseMeasureHandle')
            VarInfparams.DataBaseMeasure = feval(params.model.baseMeasureHandle,seq(tr).y,params);
            seq(tr).posterior.DataBaseMeasure = VarInfparams.DataBaseMeasure;
        end
        
        % init value
        if isfield(seq(tr),'posterior')&&isfield(seq(tr).posterior,'lamInit')
            lamInit = seq(tr).posterior.lamInit;
        else
            lamInit = zeros(yDim*T,1)+mean(vec(seq(tr).y))+1e-3;
        end
        % warm start inference if possible
        if isfield(seq(tr),'posterior')&&isfield(seq(tr).posterior,'lamOpt')
            lamInit = seq(tr).posterior.lamOpt;
        end
        
        
        lamOpt = minFunc(@VariationalInferenceDualCost,lamInit,optparams.minFuncOptions,VarInfparams);
        
        [DualCost, ~, varBound, m_ast, invV_ast, Vsm, VVsm, over_m, over_v] = VariationalInferenceDualCost(lamOpt,VarInfparams);
        
        
        seq(tr).posterior.xsm        = reshape(m_ast,xDim,T);	      % posterior mean   E[x(t)|y(1:T)]
        seq(tr).posterior.Vsm        = Vsm;			      % posterior covariances Cov[x(t),x(t)|y(1:T)]
        seq(tr).posterior.VVsm       = VVsm;			      % posterior covariances Cov[x(t+1),x(t)|y(1:T)]
        seq(tr).posterior.lamOpt     = lamOpt;		      % optimal value of dual variable
        seq(tr).posterior.lamInit    = lamInit;
        seq(tr).posterior.varBound   = varBound;		      % variational lower bound for trial
        seq(tr).posterior.DualCost   = DualCost;
        seq(tr).posterior.over_m     = over_m;		      % C*xsm+d
        seq(tr).posterior.over_v     = over_v;		      % diag(C*Vsm*C')
        
    end
end

function [NOWparams seq varBound EStepTimes MStepTimes] = PopSpikeEM(params,seq,Kin)
    %
    % [NOWparams seq varBound EStepTimes MStepTimes] = PopSpikeEM(params,seq)
    %
    % Expectation maximization algorithm for learning parameters of population model with spikes
    %
    % input:
    % params:       struct,  see PopSikeEngine.m for a definition and description
    % seq:          struct with multiple elements, see PopSikeEngine.m for a defintion and description
    %
    % output:
    % NOWparams:    struct, same as input-params but with updated and added fields
    % seq:          struct, same as input-struct but with added field 'posterior'
    % varBound:     vector, variational bound (or other cost function) for each  iteration of EM
    % EStepTimes, MStepTimes: vector, cpu-time taken by each iteration
    %
    % (c) L Buesing 01/2014
    
    
    Trials          = numel(seq);
    maxIter         = params.opts.algorithmic.EMIterations.maxIter;
    progTolvarBound = params.opts.algorithmic.EMIterations.progTolvarBound;
    maxCPUTime      = params.opts.algorithmic.EMIterations.maxCPUTime;
    
    ParamPenalizerHandle = params.model.ParamPenalizerHandle;
    InferenceMethod      = params.model.inferenceHandle;
    MstepMethod          = params.model.MStepHandle;
    
    
    EStepTimes      = nan(maxIter,1);
    MStepTimes      = nan(maxIter+1,1);
    varBound        = nan(maxIter,1);
    PREVparams      = params;			     % params for backtracking!
    NOWparams       = params;
    varBoundMax     = -inf;
    
    
    disp(['Starting PopSpikeEM using InferenceMethod  >>' char(InferenceMethod) '<<    and MStepMethod  >>' char(MstepMethod) '<<'])
    disp('----------------------------------------------------------------------------------------------------------------------------')
    
    
    Tall        = sum([seq.T]);
    EMbeginTime = cputime;
    
    %%%%%%%%%%% outer EM loop
    for ii=1:maxIter
        %try
        NOWparams.state.EMiter = ii;        
        %%%%%%% E-step: inference
        
        % do inference
        infTime = tic;
        NOWparams.opts.EMiter = ii;
        %try
        [seq varBound(ii)] = InferenceMethod(NOWparams,seq);            %For variational method, varBound for each trials is saved in seq.posterior... ?
        % catch
        %  disp('Error in inference, aborting EM iterations')
        %  break
        %end
        EStepTimes(ii) = toc(infTime);
        
        % add regularizer costs to varBound !!!
        varBound(ii) = varBound(ii) - ParamPenalizerHandle(NOWparams);
        
        fprintf('\rIteration: %i     Elapsed time (EStep): %d     Elapsed time (MStep): %d     Variational Bound: %d',ii,EStepTimes(ii),MStepTimes(ii),varBound(ii))
        
        % check termination criteria
        if params.opts.algorithmic.EMIterations.abortDecresingVarBound && (varBound(ii)<varBoundMax)    % check if varBound is increasing!
            NOWparams = PREVparams;	   % parameter backtracking
            fprintf('\n ');
            warning('Variational lower bound is decreasing, aborting EM & backtracking');
            break;
        end
        
        if params.opts.algorithmic.EMIterations.abortDecresingVarBound && ((abs(varBound(ii)-varBoundMax)/Tall)<progTolvarBound)
            fprintf('\nReached progTolvarBound for EM, aborting')
            break
        end
        
        if (cputime-EMbeginTime)>maxCPUTime            
            fprintf('\nReached maxCPUTime for EM, aborting')
            break
        end
        
        varBoundMax = varBound(ii);
        PREVparams  = NOWparams;
        
        %%%%%%% M-step
        
        mstepTime = tic;
        [NOWparams seq] = MstepMethod(NOWparams,seq);
        %NOWparams = PLDSMStep(NOWparams,seq);
        MStepTimes(ii+1) = toc(mstepTime);
        %{
catch
    NOWparams = PREVparams;     % parameter backtracking
    fprintf('\n ');
    warning('Aborting EM & backtracking');
    disp('Error in PopSpikeEM')
    break
  end
        %}
        if nargin>2
            for k=1:length(Kin),Proj{k}=seq(k).posterior.xsm;end
            [~,~,NOWparams.R2(ii,:),NOWparams.CC(ii,:)]=MyFunctionsVer1('LeaveOneOutKalmanDecoding',Kin,Proj,1);
        end
    end
    
    NOWparams.opts = rmfield(NOWparams.opts,'EMiter');
    
    fprintf('\n----------------------------------------------------------------------------------------------------------------------------\n')
    disp('EM iterations done')
end

function remain = assignopts (opts, varargin)
    % assignopts - assign optional arguments (matlab 5 or higher)
    %
    %   REM = ASSIGNOPTS(OPTLIST, 'VAR1', VAL1, 'VAR2', VAL2, ...)
    %   assigns, in the caller's workspace, the values VAL1,VAL2,... to
    %   the variables that appear in the cell array OPTLIST and that match
    %   the strings 'VAR1','VAR2',... .  Any VAR-VAL pairs that do not
    %   match a variable in OPTLIST are returned in the cell array REM.
    %   The VAR-VAL pairs can also be passed to ASSIGNOPTS in a cell
    %   array: REM = ASSIGNOPTS(OPTLIST, {'VAR1', VAL1, ...});
    %
    %   By default ASSIGNOPTS matches option names using the strmatch
    %   defaults: matches are case sensitive, but a (unique) prefix is
    %   sufficient.  If a 'VAR' string is a prefix for more than one
    %   option in OPTLIST, and does not match any of them exactly, no
    %   assignment occurs and the VAR-VAL pair is returned in REM.
    %
    %   This behaviour can be modified by preceding OPTLIST with one or
    %   both of the following flags:
    %      'ignorecase' implies case-insensitive matches.
    %      'exact'      implies exact string matches.
    %   Both together imply case-insensitive, but otherwise exact, matches.
    %
    %   ASSIGNOPTS useful for processing optional arguments to a function.
    %   Thus in a function which starts:
    %		function foo(x,y,varargin)
    %		z = 0;
    %		assignopts({'z'}, varargin{:});
    %   the variable z can be given a non-default value by calling the
    %   function thus: foo(x,y,'z',4);  When used in this way, a list
    %   of currently defined variables can easily be obtained using
    %   WHO.  Thus if we define:
    %		function foo(x,y,varargin)
    %		opt1 = 1;
    %               opt2 = 2;
    %		rem = assignopts('ignorecase', who, varargin);
    %   and call foo(x, y, 'OPT1', 10, 'opt', 20); the variable opt1
    %   will have the value 10, the variable opt2 will have the
    %   (default) value 2 and the list rem will have the value {'opt',
    %   20}.
    %
    %   See also WARNOPTS, WHO.
    %
    % Copyright (C) by Maneesh Sahani
    
    ignorecase = 0;
    exact = 0;
    
    % check for flags at the beginning
    while (~iscell(opts))
        switch(lower(opts))
            case 'ignorecase',
                ignorecase = 1;
            case 'exact',
                exact = 1;
            otherwise,
                error(['unrecognized flag :', opts]);
        end
        
        opts = varargin{1};
        varargin = varargin{2:end};
    end
    
    % if passed cell array instead of list, deal
    if length(varargin) == 1 & iscell(varargin{1})
        varargin = varargin{1};
    end
    
    if rem(length(varargin),2)~=0,
        error('Optional arguments and values must come in pairs')
    end
    
    done = zeros(1, length(varargin));
    
    origopts = opts;
    if ignorecase
        opts = lower(opts);
    end
    
    for i = 1:2:length(varargin)
        
        opt = varargin{i};
        if ignorecase
            opt = lower(opt);
        end
        
        % look for matches
        
        if exact
            match = strmatch(opt, opts, 'exact');
        else
            match = strmatch(opt, opts);
        end
        
        % if more than one matched, try for an exact match ... if this
        % fails we'll ignore this option.
        
        if (length(match) > 1)
            match = strmatch(opt, opts, 'exact');
        end
        
        % if we found a unique match, assign in the corresponding value,
        % using the *original* option name
        
        if length(match) == 1
            assignin('caller', origopts{match}, varargin{i+1});
            done(i:i+1) = 1;
        end
    end
    
    varargin(find(done)) = [];
    remain = varargin;
end

function varargout = autocorr(y,numLags,numMA,numSTD)
    %AUTOCORR Sample autocorrelation
    %
    % Syntax:
    %
    %   [acf,lags,bounds] = autocorr(y)
    %   [acf,lags,bounds] = autocorr(y,numLags,numMA,numSTD)
    %   autocorr(...)
    %
    % Description:
    %
    %   Compute the sample autocorrelation function (ACF) of a univariate,
    %   stochastic time series y. When called with no output arguments,
    %   AUTOCORR plots the ACF sequence with confidence bounds.
    %
    % Input Arguments:
    %
    %   y - Vector of observations of a univariate time series for which the
    %     sample ACF is computed or plotted. The last element of y contains the
    %     most recent observation.
    %
    % Optional Input Arguments:
    %
    %   numLags - Positive integer indicating the number of lags of the ACF
    %     to compute. If empty or missing, the default is to compute the ACF at
    %     lags 0,1,2, ... T = min[20,length(y)-1]. Since ACF is symmetric
    %     about lag zero, negative lags are ignored.
    %
    %   numMA - Nonnegative integer indicating the number of lags beyond which
    %     the theoretical ACF is deemed to have died out. Under the hypothesis
    %     that the underlying y is really an MA(numMA) process, the large-lag
    %     standard error is computed via Bartlett's approximation for lags >
    %     numMA as an indication of whether the ACF is effectively zero beyond
    %     lag numMA. If numMA is empty or missing, the default is numMA = 0, in
    %     which case y is assumed to be Gaussian white noise. If y is a
    %     Gaussian white noise process of length N, the standard error will be
    %     approximately 1/sqrt(N). numMA must be less than numLags.
    %
    %   numSTD - Positive scalar indicating the number of standard deviations
    %     of the sample ACF estimation error to compute, assuming the
    %     theoretical ACF of y is zero beyond lag numMA. When numMA = 0 and y
    %     is a Gaussian white noise process of length numMA, specifying numSTD
    %     will result in confidence bounds at +/-(numSTD/sqrt(numMA)). If empty
    %     or missing, the default is numSTD = 2 (approximate 95% confidence).
    %
    % Output Arguments:
    %
    %   acf - Sample autocorrelation function of y. acf is a vector of
    %     length numLags+1 corresponding to lags 0,1,2,...,numLags. The first
    %     element of acf is unity (i.e., acf(1) = 1 at lag 0).
    %
    %   lags - Vector of lags corresponding to acf (0,1,2,...,numLags).
    %
    %   bounds - Two-element vector indicating the approximate upper and lower
    %     confidence bounds, assuming that y is an MA(numMA) process. Note that
    %     bounds is approximate for lags > numMA only.
    %
    % Example:
    %
    %   % Create an MA(2) process from a sequence of 1000 Gaussian deviates,
    %   % and assess whether the ACF is effectively zero for lags > 2:
    %
    %     x = randn(1000,1);         % 1000 Gaussian deviates ~ N(0,1)
    %     y = filter([1 -1 1],1,x);  % Create an MA(2) process
    %     autocorr(y,[],2)           % Inspect the ACF with 95% confidence
    %
    % Reference:
    %
    %   [1] Box, G. E. P., G. M. Jenkins, and G. C. Reinsel. Time Series
    %       Analysis: Forecasting and Control. 3rd edition. Upper Saddle River,
    %       NJ: Prentice-Hall, 1994.
    %
    % See also CROSSCORR, PARCORR, FILTER.
    
    % Copyright 1999-2010 The MathWorks, Inc.
    % $Revision: 1.1.8.4 $  $Date: 2010/10/08 16:41:02 $
    
    % Ensure the sample data is a vector:
    
    [rows,columns] = size(y);
    
    if (rows ~= 1) && (columns ~= 1)
        
        error(message('econ:autocorr:NonVectorInput'))
        
    end
    
    rowSeries = (size(y,1) == 1);
    
    y = y(:);         % Ensure a column vector
    N = length(y);    % Sample size
    defaultLags = 20; % Recommendation of [1]
    
    % Ensure numLags is a positive integer or set default:
    
    if (nargin >= 2) && ~isempty(numLags)
        
        if numel(numLags) > 1
            
            error(message('econ:autocorr:NonScalarLags'))
            
        end
        
        if (round(numLags) ~= numLags) || (numLags <= 0)
            
            error(message('econ:autocorr:NonPositiveInteger'))
            
        end
        
        if numLags > (N-1)
            
            error(message('econ:autocorr:LagsTooLarge'))
            
        end
        
    else
        
        numLags = min(defaultLags,N-1); % Default
        
    end
    
    
    % Ensure numMA is a nonnegative integer or set default:
    
    if (nargin >= 3) && ~isempty(numMA)
        
        if numel(numMA) > 1
            
            error(message('econ:autocorr:NonScalarNMA'))
            
        end
        
        if (round(numMA) ~= numMA) || (numMA < 0)
            
            error(message('econ:autocorr:NegativeIntegerNMA'))
            
        end
        
        if numMA >= numLags
            
            error(message('econ:autocorr:NMATooLarge'))
            
        end
        
    else
        
        numMA = 0; % Default
        
    end
    
    % Ensure numSTD is a positive scalar or set default:
    
    if (nargin >= 4) && ~isempty(numSTD)
        
        if numel(numSTD) > 1
            
            error(message('econ:autocorr:NonScalarSTDs'))
            
        end
        
        if numSTD < 0
            
            error(message('econ:autocorr:NegativeSTDs'))
            
        end
        
    else
        
        numSTD = 2; % Default
        
    end
    
    % Convolution, polynomial multiplication, and FIR digital filtering are all
    % the same operation. The FILTER command could be used to compute the ACF
    % (by convolving the de-meaned y with a flipped version of itself), but
    % FFT-based computation is significantly faster for large data sets.
    
    % The ACF computation is based on [1], pages 30-34, 188:
    
    nFFT = 2^(nextpow2(length(y))+1);
    F = fft(y-mean(y),nFFT);
    F = F.*conj(F);
    acf = ifft(F);
    acf = acf(1:(numLags+1)); % Retain non-negative lags
    acf = acf./acf(1); % Normalize
    acf = real(acf);
    
    % Compute approximate confidence bounds using the approach in [1],
    % equations 2.1.13 and 6.2.2, pp. 33 and 188, respectively:
    
    sigmaNMA = sqrt((1+2*(acf(2:numMA+1)'*acf(2:numMA+1)))/N);
    bounds = sigmaNMA*[numSTD;-numSTD];
    lags = (0:numLags)';
    
    if nargout == 0
        
        %  Plot the sample ACF:
        
        lineHandles = stem(lags,acf,'filled','r-o');
        set(lineHandles(1),'MarkerSize',4)
        grid('on')
        xlabel('Lag')
        ylabel('Sample Autocorrelation')
        title('Sample Autocorrelation Function')
        hold('on')
        
        %  Plot confidence bounds (horizontal lines) under the hypothesis that the
        %  underlying y is really an MA(numMA) process. Bartlett's approximation
        %  gives an indication of whether the ACF is effectively zero beyond lag
        %  numMA. For this reason, the confidence bounds appear over the ACF only
        %  for lags greater than numMA (i.e., numMA+1, numMA+2, ... numLags). In
        %  other words, the confidence bounds enclose only those lags for which the
        %  null hypothesis is assumed to hold.
        
        plot([numMA+0.5 numMA+0.5; numLags numLags],[bounds([1 1]) bounds([2 2])],'-b');
        plot([0 numLags],[0 0],'-k');
        hold('off')
        a = axis;
        axis([a(1:3) 1]);
        
    else
        
        %  Re-format outputs for compatibility with the y input. When y is input as
        %  a row vector, then pass the outputs as a row vectors; when y is a column
        %  vector, then pass the outputs as a column vectors.
        
        if rowSeries
            
            acf = acf';
            lags = lags';
            bounds = bounds';
            
        end
        
        varargout = {acf,lags,bounds};
        
    end
end

function seq = binSpikeTimes(spikeTimes,dt,varargin);
    %
    % seq = binSpikeTimes(spikeTimes,dt,varargin)
    %
    % convert data in noise-workshop data format to seq struct
    %
    
    Tmax = [];
    
    assignopts(who,varargin);
    
    yDim = numel(spikeTimes);
    seq = [];
    
    
    Trials = numel(spikeTimes{1});
    
    for tr=1:Trials
        
        if isempty(Tmax)
            Tall = [];
            for yd=1:yDim
                Tall = [Tall vec(spikeTimes{yd}{tr})'];
            end
            T = ceil(max(Tall));
        else
            T = Tmax;
        end
        
        bins = 0:dt:T;
        
        seq(tr).y = zeros(yDim,numel(bins));
        
        for yd=1:yDim
            if ~isempty(spikeTimes{yd}{tr})
                seq(tr).y(yd,:) = histc(spikeTimes{yd}{tr},bins)';
            end
        end
        seq(tr).y = seq(tr).y(:,1:end-1);
        seq(tr).T = size(seq(tr).y,2);
        
        
    end
end

function H=circle(center,radius,NOP,style)
    %---------------------------------------------------------------------------------------------
    % H=CIRCLE(CENTER,RADIUS,NOP,STYLE)
    % This routine draws a circle with center defined as
    % a vector CENTER, radius as a scaler RADIS. NOP is
    % the number of points on the circle. As to STYLE,
    % use it the same way as you use the rountine PLOT.
    % Since the handle of the object is returned, you
    % use routine SET to get the best result.
    %
    %   Usage Examples,
    %
    %   circle([1,3],3,1000,':');
    %   circle([2,4],2,1000,'--');
    %
    %   Zhenhai Wang <zhenhai@ieee.org>
    %   Version 1.00
    %   December, 2002
    %---------------------------------------------------------------------------------------------
    
    if (nargin <3),
        error('Please see help for INPUT DATA.');
    elseif (nargin==3)
        style='b-';
    end;
    THETA=linspace(0,2*pi,NOP);
    RHO=ones(1,NOP)*radius;
    [X,Y] = pol2cart(THETA,RHO);
    X=X+center(1);
    Y=Y+center(2);
    H=plot(X,Y,style);
    axis square;
    
end

function v = logdet(A, op)
    %LOGDET Computation of logarithm of determinant of a matrix
    %
    %   v = logdet(A);
    %       computes the logarithm of determinant of A.
    %
    %       Here, A should be a square matrix of double or single class.
    %       If A is singular, it will returns -inf.
    %
    %       Theoretically, this function should be functionally
    %       equivalent to log(det(A)). However, it avoids the
    %       overflow/underflow problems that are likely to
    %       happen when applying det to large matrices.
    %
    %       The key idea is based on the mathematical fact that
    %       the determinant of a triangular matrix equals the
    %       product of its diagonal elements. Hence, the matrix's
    %       log-determinant is equal to the sum of their logarithm
    %       values. By keeping all computations in log-scale, the
    %       problem of underflow/overflow caused by product of
    %       many numbers can be effectively circumvented.
    %
    %       The implementation is based on LU factorization.
    %
    %   v = logdet(A, 'chol');
    %       If A is positive definite, you can tell the function
    %       to use Cholesky factorization to accomplish the task
    %       using this syntax, which is substantially more efficient
    %       for positive definite matrix.
    %
    %   Remarks
    %   -------
    %       logarithm of determinant of a matrix widely occurs in the
    %       context of multivariate statistics. The log-pdf, entropy,
    %       and divergence of Gaussian distribution typically comprises
    %       a term in form of log-determinant. This function might be
    %       useful there, especially in a high-dimensional space.
    %
    %       Theoretially, LU, QR can both do the job. However, LU
    %       factorization is substantially faster. So, for generic
    %       matrix, LU factorization is adopted.
    %
    %       For positive definite matrices, such as covariance matrices,
    %       Cholesky factorization is typically more efficient. And it
    %       is STRONGLY RECOMMENDED that you use the chol (2nd syntax above)
    %       when you are sure that you are dealing with a positive definite
    %       matrix.
    %
    %   Examples
    %   --------
    %       % compute the log-determinant of a generic matrix
    %       A = rand(1000);
    %       v = logdet(A);
    %
    %       % compute the log-determinant of a positive-definite matrix
    %       A = rand(1000);
    %       C = A * A';     % this makes C positive definite
    %       v = logdet(C, 'chol');
    %
    
    %   Copyright 2008, Dahua Lin, MIT
    %   Email: dhlin@mit.edu
    %
    %   This file can be freely modified or distributed for any kind of
    %   purposes.
    %
    
    %% argument checking
    
    assert(isfloat(A) && ndims(A) == 2 && size(A,1) == size(A,2), ...
        'logdet:invalidarg', ...
        'A should be a square matrix of double or single class.');
    
    if nargin < 2
        use_chol = 0;
    else
        assert(strcmpi(op, 'chol'), ...
            'logdet:invalidarg', ...
            'The second argument can only be a string ''chol'' if it is specified.');
        use_chol = 1;
    end
    
    %% computation
    
    if use_chol
        v = 2 * sum(log(diag(chol(A))));
    else
        [L, U, P] = lu(A);
        du = diag(U);
        c = det(P) * prod(sign(du));
        v = log(c) + sum(log(abs(du)));
    end
    
end

function plotMatrixSpectrum(A,varargin)
    %
    % plot matrix spectra
    %
    % (c) Lars Buesing, 2014
    
    figh = -1;
    col  = 'r';
    linw  = 2.0;
    linwc = 2.0;
    assignopts(who,varargin);
    
    if figh<0
        figure
    else
        figure(figh)
    end
    
    hold on;
    p=circle([0,0],1,1000,'k');
    set(p,'LineWidth',linwc)
    
    if (size(A,1)==1) || (size(A,2)==1)
        EigA = A;
    else
        EigA = eig(A);
    end
    
    plot(real(EigA),imag(EigA),['x' col],'Linewidth',linw)
    axis off;
end

function plotPosterior(seq,trId,params)
    %
    %
    %
    
    [xDim T] = size(seq(trId).posterior.xsm);
    
    try
        Pidx = [1 1+params.dE];
    catch
        Pidx = [1 2];
    end
    
    if nargin<1.5
        trId = 1;
    end
    
    xsm  = seq(trId).posterior.xsm;
    xerr = zeros(size(xsm));
    for t=1:T
        xidx = (t-1)*xDim+1:t*xDim;
        xerr(:,t) = sqrt(diag(seq(trId).posterior.Vsm(xidx,:)));
    end
    
    figure; hold on; title('posterior')
    
    for i=1:numel(Pidx)
        subplot(numel(Pidx),1,i); hold on;
        pidx = Pidx(i);
        errorbar(1:T,xsm(pidx,:),xerr(pidx,:),'r')
        try;plot(1:T,seq(trId).x(pidx,:),'linewidth',2);end;
        plot(1:T,xsm(pidx,:),'r','linewidth',2)
        ylabel('x(t)');
        if i==numel(Pidx);xlabel('t');end
        %figSize = {14,10};
        %figuresize(figSize{:},'centimeters')
    end
    
    if isfield(seq(1).posterior,'phi')
        try
            figure
            plot(seq(1).posterior.phi,'x','MarkerSize',10,'linewidth',2)
            xlabel('neuron no');ylabel('phi')
            %figSize = {14,10};
            %figuresize(figSize{:},'centimeters')
            
            %{
  figure
  piNow = exp(seq(1).posterior.phi);
  piNow = bsxfun(@times,piNow,1./sum(piNow,2));
  plot(piNow,'x','MarkerSize',10,'linewidth',2)
  %plot(seq(1).posterior.phi,'rx','MarkerSize',10,'linewidth',2)
  xlabel('neuron no');ylabel('pi')
  %figSize = {14,10};
  %figuresize(figSize{:},'centimeters')
            %}
        end
    end
end

function D = PoissonRegression(Y,U,varargin)
    %
    % D = PoissonRegression(Y,U,varargin)
    %
    % Y  ~  Poisson(D*U+x)
    %
    % where x ~ Normal(over_m,V) with over_v = diag(V);
    %
    %
    % INPUT:
    %
    % - data Y
    % - observed variates U
    % - lam:     penalizer for L2 regularization of D  [optional, default = 0.1]
    % - over_m:  mean of x   [optional, default 0]
    % - over_v:  diagonal of covariance of x  [optional, default 0]
    % - Dinit:   inital value of D  [optinal, default 0]
    % - options: minFunc optimzation options [optinal]
    %
    % L Buesing, 2014
    
    
    yDim = size(Y,1);
    uDim = size(U,1);
    
    lam    = 0.1;
    over_m = [];
    over_v = [];
    Dinit  = zeros(yDim,uDim);
    
    options.Display     = 'iter';
    options.Method      = 'lbfgs';
    options.MaxIter     = 5000;
    options.maxFunEvals = 50000;
    options.progTol     = 1e-9;
    options.optTol      = 1e-5;
    
    assignopts(who,varargin);
    
    
    D = minFunc(@PoissonRegressionCost,vec(Dinit),options,Y,U,lam,over_m,over_v);
    D = reshape(D,yDim,uDim);
    
    
end

function [f,df] =  PoissonRegressionCost(vecD,y,u,lam,over_m,over_v)
    %
    % [f,df] =  PoissonRegressionCost(vecD,y,u,lam,over_m,over_v)
    %
    % Poisson regression cost funtion
    %
    %
    
    [yDim T] = size(y);
    uDim     = size(u,1);
    
    D    = reshape(vecD,yDim,uDim);
    Du   = D*u;
    nu   = Du;
    if ~isempty(over_m); nu = nu + over_m; end
    yhat = nu;
    if ~isempty(over_m); yhat = yhat + 0.5*over_v;end
    yhat = exp(yhat);
    
    f  = sum(vec(-y.*nu+yhat));
    df = (yhat-y)*u';
    
    % L2 regularization
    f  = f  + 0.5*lam*norm(D,'fro').^2;
    df = df + lam*D;
    
    df = vec(df);
end

function seqnew = rebinRaster(seq,dt)
    %
    % function seq = rebinRaster(seq,dt)
    %
    % rebin seq by a factor of dt
    %
    
    
    Trials = numel(seq);
    yDim   = size(seq(1).y,1);
    
    if isfield(seq,'x')
        xDim = size(seq(1).x,1);
    end
    
    seqnew = seq;
    
    for tr=1:Trials
        
        T    = size(seq(tr).y,2);
        Tnew = floor(T/dt);
        
        yold = reshape(seq(tr).y(:,1:Tnew*dt),yDim,dt,Tnew);
        %     ynew = squeeze(sum(yold,2)); % wrong result if yDim==1
        ynew = reshape(sum(yold,2), yDim, Tnew); %
        
        seqnew(tr).y = ynew;
        seqnew(tr).T = Tnew;
        
        if isfield(seq,'yr')
            yrold = reshape(seq(tr).yr(:,1:Tnew*dt),yDim,dt,Tnew);
            %        yrnew = squeeze(sum(yrold,2));
            yrnew = reshape(sum(yrold,2), yDim, Tnew);
            seqnew(tr).yr = yrnew;
        end
        
        if isfield(seq,'x')
            xold = reshape(seq(tr).x(:,1:Tnew*dt),xDim,dt,Tnew);
            %        xnew = squeeze(sum(xold,2));
            xnew = reshape(sum(xold,2), xDim, Tnew);
            seqnew(tr).x = xnew;
        end
        
    end
end

function samples = sampleGPPrior(N,T,Bdim,varargin)
    %
    % samples = sampleGPPrior(N,T,Bdim,varargin)
    %
    %
    
    
    xpos = 1:T; % sampling locations
    tau  = 10; % sqaured-exp length scale
    sig  = 0.01; % uncorrelated noise
    
    assignopts(who, varargin);
    
    sqx  = xpos.^2;
    Dist = repmat(sqx,T,1)+repmat(sqx,T,1)'-2*xpos'*xpos;
    K    = exp(-Dist./(tau.^2))*(1-sig.^2)+sig.^2*eye(T);
    C    = chol(K)';
    
    for n=1:N
        samples{n} = (C*randn(T,Bdim))';
    end
end

function M = spblkdiag(blocks)
    %function M = spblkdiag(blocks)
    %
    % function generating efficiently a sparse matrix containing
    % subblocks blocks(:,:,i) as block i along the diagonal.
    % this is 1000 times faster than blkdiag!!!!
    %
    % Bernard Haasdonk 31.8.2009
    
    % This program is open source.  For license terms, see the COPYING file.
    %
    % --------------------------------------------------------------------
    % ATTRIBUTION NOTICE:
    % This product includes software developed for the RBmatlab project at
    % (C) Universities of Stuttgart and MÃ?nster, Germany.
    %
    % RBmatlab is a MATLAB software package for model reduction with an
    % emphasis on Reduced Basis Methods. The project is maintained by
    % M. Dihlmann, M. Drohmann, B. Haasdonk, M. Ohlberger and M. Schaefer.
    % For Online Documentation and Download we refer to www.morepas.org.
    % --------------------------------------------------------------------
    
    
    [n,m,k] = size(blocks);
    
    row_ind = (1:n)'*ones(1,m*k);
    row_offs =(ones(n*m,1)*(0:(k-1)))*n;
    row_ind = row_ind(:)+row_offs(:);
    
    col_ind = ones(n,1)*(1:(m*k));
    col_ind = col_ind(:);
    
    M = sparse(row_ind,col_ind,blocks(:));%| \docupdate
end

function c = struct2arglist(s)
    %STRUCT2ARGLIST  Convert structure to cell array of fields/values.
    % STRUCT2ARGLIST(S) returns a cell array {'field1',value1,'field2',value2,...}
    % It is the opposite of MAKESTRUCT.
    %
    % Example:
    %   function f(varargin)
    %   opt.FontSize = 10;
    %   opt = setfields(opt,makestruct(varargin),'ignore');
    %   varargin = struct2arglist(opt);
    %   g(varargin{:});
    %
    % See also MAKESTRUCT.
    
    % Written by Tom Minka
    % (c) Microsoft Corporation. All rights reserved.
    
    f = fieldnames(s);
    c = cell(1,2*length(f));
    for i = 1:length(f)
        c{2*i-1} = f{i};
        c{2*i} = s.(f{i});
    end
end

function uSub = subsampleSignal(u,dt);
    %
    % function uSub = subsampleSignal(u,dt);
    %
    
    uSub = [u];
    T    = size(uSub,2);
    Tf   = floor(T/dt)*dt;
    uSub = reshape(uSub(:,1:Tf),size(uSub,1),dt,Tf/dt);
    uSub = squeeze(uSub(:,1,:));
end

function  [D, OD] = sym_blk_tridiag_inv_v1(AA,BB,adx,bdx)
    % Compute block tridiagonal terms of the inverse of a *symmetric* block
    % tridiagonal matrix.
    %
    % Note: Could be generalized to non-symmetric matrices, but it is not currently implemented.
    % Note: Could be generalized to compute *all* blocks of the inverse, but
    % it's not necessary for my current application and hence not implemented.
    %
    % Note: AA and BB could be combined into a single variable, but I think that might get confusing.
    %
    % Input:
    %   AA - (n x n x Ka) unique diagonal blocks
    %   BB - (n x n x Kb) unique off-diagonal blocks
    %  adx - (T x 1) vector such that (i,i)th block of A is
    %                   A_{ii} = AA(:,:,adx(ii))
    %  bdx - (T-1 x 1) vector such that (i,i+1) block of A is
    %                   A_{i,i+1} = AA(:,:,bdx(ii))
    %
    % Output:
    %   D  - (n x n x T) diagonal blocks of the inverse
    %  OD  - (n x n x T-1) off-diagonal blocks of the inverse
    %
    % From:
    % Jain et al, 2006
    % "Numerically Stable Algorithms for Inversion of Block Tridiagonal and Banded Matrices"
    % (c) Evan Archer, 2014
    %
    assert(numel(size(AA)) == 3, 'Always a 3d-array. For scalar blocks, make A of size [1 x 1 x Ka].')
    assert(length(adx) == length(bdx)+1, 'Should be one less upper diagonal term than diagonal term.')
    assert(length(size(AA)) == 3, 'Always expect the block index on the 3rd dimension');
    assert(size(AA,1) == size(BB,1) && size(AA,2) == size(BB,2) );
    assert(size(AA,1) == size(AA,2), 'Input matrix must be square.')
    
    % We only need R when our matrix is non-symmetric. for us, it always will be
    % R = zeros(runinfo.nStateDim, runinfo.nStateDim, T);
    BB = -BB; % we gotta make them the negative of the blocks
    T = numel(adx);
    n = size(AA,1);
    S = zeros(n, n, T-1);
    D = zeros(n, n, T); % diagonal
    OD = zeros(n, n, T-1); % off diagonal
    III = eye(n);
    % R(:,:,1) = AA0\BB;
    S(:,:,end) = BB(:,:,bdx(end))/ AA(:,:,adx(end));
    for idx = (T-2):-1:1
        %    R(:,:,idx) = (AA - BB'*R(:,:,idx-1))\BB;
        S(:,:,idx) = BB(:,:,bdx(idx)) / (AA(:,:,adx(idx+1)) - S(:,:,idx+1)*BB(:,:,bdx(idx+1))');
    end
    % Compute diagonal and off-diagonal blocks
    D(:,:,1) = pinv(AA(:,:,adx(1)) - BB(:,:,bdx(1))*S(:,:,1)');
    OD(:,:,1) = S(:,:,1)'*D(:,:,1);
    for idx = 2:T-1
        D(:,:,idx) = (AA(:,:,adx(idx)) - BB(:,:,bdx(idx)) * S(:,:,idx)')\(III + BB(:,:,bdx(idx-1))'*D(:,:,idx-1)*S(:,:,idx-1));
        OD(:,:,idx) = S(:,:,idx)'*D(:,:,idx);
    end
    D(:,:,end) = AA(:,:,adx(end)) \ (III + BB(:,:,bdx(T-1))'*D(:,:,end-1)*S(:,:,end));
end

function params = touchField(params,ftouch,fval);
    %
    % params = touchField(params,ftouch,fval);
    %
    
    
    if nargin<2.5
        fval = [];
    end
    
    if ~isfield(params,ftouch)
        params = setfield(params,ftouch,fval);
    end
end

function v=vec(v)
    
    v = v(:);
end