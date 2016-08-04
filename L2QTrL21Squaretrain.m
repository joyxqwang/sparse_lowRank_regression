%% function L2QTrL21Squaretrain for Bioinformatics
% Author: Xiaoqian Wang
% \min_{W,Qi} ||X'*W-Y||_F^2 + gamma1*(\sum_i^numG||W*Qi||_Sp^p)^k + gamma2*||W||_{2,q}
function [outW, outQ, outObj, outNumIter] = L2QTrL21Squaretrain(X, Y, numG, k, p, q, gamma1, gamma2, maxIter, inThresh)
% input:
    % n: number of samples
    % dim: number of SNPs
    % c: number of QTs
    % X: dim*n SNP data
    % Y: n*(cT) phenotype matrix
    % numG: number of groups
% output:
    % outW: dim*cT weight matrix
    % outQ: cT*numG group indicator
    % outObj: objective function value

%% Initialization
%
if nargin < 10
    inThresh = 10^-8;
end
if nargin < 9
    maxIter = 10^4;
end
if nargin < 8
    gamma2 = 1;
end
if nargin < 7
    gamma1 = 1;
end
if nargin < 6
    q = 0.15;
end
if nargin < 5
    p = 0.15;
end
if nargin < 4
    k = 2.5;
end
delta = 10^-8;

%
cT = size(Y, 2);
dim = size(X, 1);

%
s = RandStream.create('mt19937ar','seed',7);  %% seed, "457" can be changed
RandStream.setGlobalStream(s);
A = rand(cT, numG);
Q = A./(sum(A,2)*ones(1,numG));

%
XY = X*Y;
XX = X*X';
W = (XX + gamma2*eye(dim))\XY;

%
D = getQtraceSquare(W, Q, delta, p, k);

%% Main code
obj = zeros(maxIter, 1);

for iter = 1: maxIter
    
    % fix D, Q, W, update B    
    b22 = sum(W.*W, 2) + delta;
    b2 = (p/2)*b22.^(p/2-1);
    B = diag(b2);

    % fix D, W, update Q
    for g = 1: numG
        A(:, g) = diag(W'*D{g}*W);
    end
    for t = 1: cT
        a = A(t, :);
        Q(t,:) = 1/sum(1./a)*1./a;
    end

    % fix B, Q, D, update W    
    for t = 1: cT
        tmpD = zeros(dim, dim);
        for g = 1: numG
            tmpD = tmpD + D{g}*(Q(t, g)^2);
        end
        W(:,t) = (XX + gamma1*tmpD + gamma2*B)\XY(:,t);
    end
    
    % fix W, Q, update D
    [D, Qtr] = getQtraceSquare(W, Q, delta, p, k);

    % calculate obj
    Loss  =  norm(Y - X'*W, 'fro')^2;
    tmpW = W.*W;   
    obj(iter) = Loss + gamma1*Qtr + gamma2*sum(sum(tmpW,2).^(q/2));

    if(iter > 1)
        if((obj(iter-1) - obj(iter))/obj(iter-1) < inThresh)
            break;
        end
    end

    if mod(iter, 10) == 0
        fprintf('process iteration %d, the obj is %d ...\n', iter, obj(iter));
    end
    
end

%% Outputs
%
[min_val, qMat] = max(Q, [], 2);
outQ = zeros(cT,numG);
for t = 1:cT     
    outQ(t, qMat(t)) = 1;
end

%
outNumIter = iter;

%
outObj = obj(1:iter);

%
outW = W;

end



% Di = k*p/2 * (||W*Qi||_Sp^p)^(k-1) * (W*Qi*W' + delta*eye(dim))^(0.5*p-1)
function [D, Qtr] = getQtraceSquare(W, Q, delta, p, k)
% W: dim*cT weight matrix
% Q: cT*numG group indicator

%% Initialization
%
numG = size(Q, 2);
[dim, cT] = size(W);
Qtr = 0;

if dim > cT
    z = zeros(dim-cT,1);
else
    z = [];
end       

%% Main code
for g = 1:numG

    WQ = W*diag(Q(:,g));
    [U,S,V] = svd(WQ);
    s = [diag(S); z]; 
    sv(:,g)=[sum(s); sum(s.^p)^k; s];
    s0 = sqrt(s.^2+delta).^p;
    Qtrg = sum(s0);
    Qtr = Qtr + Qtrg^k;
    d = p/2*(sqrt(s.^2+delta).^(p-2));
    D{g} = k*Qtrg^(k-1)*U*diag(d)*U';
    
end

end
