%6/23/17
%Get distribution of T for psi = 0, 0.01, 0.25, 0.49, 0.5

clear;


L = 2;

%1-theta
pC = 0.5;
%theta
pT = 0.5;

%vector of probability of failure and being assigned to control and probability of success and being assigned to control
pCvec = [0.5,0.5]*pC;
%vector of probability of failure and being assigned to treatment and probability of success and being assigned to treatment
pTvec = [0.5,0.5]*pT;



%Options for quadprog
options = optimoptions('quadprog','Algorithm','interior-point-convex','Display','off','TolFun',1e-10,'TolX',1e-10);





numGam = 2*L; %these numbers will be used often
numPI = L^2;
numBoth = numGam + numPI;

j = (1:numPI)+numGam;
Amat1 = sparse(1:numPI,j,-1,numPI,numBoth);
bvec1 = zeros(numPI,1);
dvec1 = -2 * [pCvec, pTvec, zeros(1,numPI)];
v = [pC * ones(1, L),pT* ones(1, L)]*2;
Dmat1 = sparse(1:numGam,1:numGam,v,numBoth,numBoth);

mat = transpose(reshape(1:numPI,[L,L]));
benPI = transpose(nonzeros(triu(mat,1)));
i1 = ones(1,numPI);
j1 = (1:numPI) + numGam;
i2 = 2*ones(1,L*(L-1)/2);
j2 = numGam + benPI;
i3 = sort(repmat(1:numGam,1,L))+2;
j3 = [1:numPI,reshape(mat,[1,numPI])]+numGam;
i4 = (1:numGam)+2;
j4 = 1:numGam;
i = [i4,i1,i2,i3];
j = [j4,j1,j2,j3];
v = ones(1,numel(i));
v(1:numel(i4)) = -1;
Aeq1 = sparse(i,j,v,2+numGam,numBoth);

clear i j i1 i2 i3 i4 j1 j2 j3 j4 

bvec2 = zeros(numGam,1);
Amat2 = sparse(1:numGam,1:numGam,-1,numGam,numGam);
beq2 = [1,1]';
Aeq2 = [ones(1,L),zeros(1,L); zeros(1,L),ones(1,L)];
dvec2 = -2*[pCvec, pTvec];
v = [pC * ones(1, L),pT* ones(1, L)]*2;
Dmat2 = sparse(1:numGam,1:numGam,v,numGam,numGam);
[gammaHat,val2] = quadprog(Dmat2,dvec2,Amat2,bvec2,Aeq2,beq2,[],[],[],options);
val2 = val2 + 1;


sigma = zeros(numGam,numGam);
v = [pC * ones(1, L),pT* ones(1, L)];
j = [pCvec, pTvec];
for r = 1:numGam
    for c = 1:numGam
      if r == c
          sigma(r,r) = j(r) - 2*gammaHat(r)*j(r)+gammaHat(r)^2*v(r);
      elseif r > L && c <= L
          sigma(r,c) = 0;
      elseif c > L && r <= L
          sigma(r,c) = 0;
      else
          sigma(r,c) = -gammaHat(r)*j(c)-gammaHat(c)*j(r) + gammaHat(c)*gammaHat(r)*v(r); 
      end
    end
end
sigma = sigma * 4; 


bvec34 = zeros(numPI,1);
beq4 = zeros(numGam * 2,1);
beq3 = zeros(numGam * 2 + 1,1);
    
Amat34 = sparse(1:numPI,1:numPI,-1,numPI,numPI+numGam*2); 
Dmat34 = sparse((1:numGam)+numPI,(1:numGam)+numPI,1,numPI+numGam*2,numPI+numGam*2);

i1 = sort(repmat(1:numGam,1,L));
j1 = [1:numPI,reshape(mat,1,numPI)];
v1 = ones(1,numel(i1));
i2 = 1:numGam*2;
j2 = [(1:numGam)+numPI + numGam,(1:numGam)+numPI + numGam];
v2 = -1 * ones(1,numel(i2));
i3 = (numGam+1):(2*numGam);
v3 = ones(1,numel(i3));
j3 = (1:numGam)+numPI;
i4 = repmat(i3,1,numPI);
j4 = sort(repmat(1:numPI,1,numGam));
v4 = repmat(transpose(gammaHat),1,numPI);
i = [i1,i2,i3,i4];
j = [j1,j2,j3,j4];
v = [v1,v2,v3,v4];
Aeq4 = sparse(i,j,v,numGam*2,numPI+numGam*2);
clear i1 i2 i3 i4 j1 j2 j3 j4 v1 v2 v3 v4

psi = [0,0.01,0.25,0.49,0.5];
ndraw  = 100000;
tdraws = zeros(length(psi),ndraw);
    
%iPsi,jPsi are used to make Aeq3 inside the loop
jPsi = [j,1:numPI,benPI];
iPsi = [i,(1+numGam*2)*ones(1,numPI+numel(benPI))];

clear i j

for p = 1:length(psi)
    rng('default');
    rng(p);
    
    vPsi = [v,psi(p)*ones(1,numel(1:numPI)),-1*ones(1,numel(benPI))];
    Aeq3 = sparse(iPsi,jPsi,vPsi,1+numGam*2,numPI+numGam*2);
        
    for d = 1:ndraw
        Z = mvnrnd(zeros(numGam,1),sigma);
        dvec34 = [zeros(1,numPI),Z,zeros(1,numGam)];
        [~,val4] = quadprog(Dmat34,dvec34,Amat34,bvec34,Aeq4,beq4,[],[],[],options);
        [~,val3] = quadprog(Dmat34,dvec34,Amat34,bvec34,Aeq3,beq3,[],[],[],options);
        tdraws(p,d) = val3 - val4; 
    end
end

FileName = ['result',datestr(now, 'dd-mmm-yyyy'),'.mat'];
save(FileName, 'tdraws', 'psi')

exit
