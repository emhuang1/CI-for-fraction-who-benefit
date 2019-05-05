clear;

job_id_string = getenv('SGE_TASK_ID');
val = str2double(job_id_string);

%addpath('/users/mrosen/CPLEX_STUDIO/cplex/matlab')
addpath('/users/emhuang/CIpropBenefit/binary_noHarm_5050/testDatasets')

n = 500;

options = optimoptions('quadprog','Algorithm','interior-point-convex','Display','off','TolFun',1e-10,'TolX',1e-10);
%msg = 'Problem with solving quadratic program.';

data = importdata('N500.txt');
pC = data(val,1);
pT = data(val,2);
pC1 = data(val,3);
pC2 = data(val,4);
pT1 = data(val,5);
pT2 = data(val,6);

psi = 0:0.01:1;

Amat1 = [eye(2);[1,1];[-1,-1]];
Dmat1 = 2 * eye(2);

bvec2 = [0; 0; 0; 1; -1];
Amat2 = [eye(3); [1,1,1]; [-1,-1,-1]];
dvec2 = 2 * [pC1 + pT1, pC1 + pT2, pC2 + pT2];
Dmat2 = eye(3);
Dmat2(1,2) = pC;
Dmat2(2,1) = pC;
Dmat2(2,3) = pT;
Dmat2(3,2) = pT;
Dmat2 = Dmat2 * 2;
%[x2, val2] = cplexqp(Dmat2,-dvec2,-Amat2,-bvec2);
[x2, val2] = quadprog(Dmat2,-dvec2,-Amat2,-bvec2,[],[],[],[],[],options);
%if exitflag ~= 1
%    error(msg)
%end
val2 = val2 + 1;
gammaHat = [x2(1)+x2(2); x2(3); x2(1); x2(2)+x2(3)];

sigma = zeros(4);
sigma(1,1) = pC1 - 2*gammaHat(1)*pC1 + gammaHat(1)^2*pC;
sigma(2,2) = pC2 - 2*gammaHat(2)*pC2 + gammaHat(2)^2*pC;
sigma(3,3) = pT1 - 2*gammaHat(3)*pT1 + gammaHat(3)^2*pT;
sigma(4,4) = pT2 - 2*gammaHat(4)*pT2 + gammaHat(4)^2*pT;
sigma(1,2) = -gammaHat(2)*pC1 - gammaHat(1)*pC2 + gammaHat(1)*gammaHat(2)*pC;
sigma(2,1) = sigma(1,2);
sigma(3,4) = -gammaHat(4)*pT1 - gammaHat(3)*pT2 + gammaHat(3)*gammaHat(4)*pT;
sigma(4,3) = sigma(3,4);
sigma = sigma * 4; 
    
    
bvec3 = [zeros(19,1);-10*sqrt(n)/log(n);zeros(2,1)];
bvec4 = [zeros(19,1);-10*sqrt(n)/log(n)];

    
Amat4 = zeros(20,11);
    
Amat4(1:3,1:3) = eye(3);

Amat4(4,1:2)   = 1;
Amat4(4,8)   = -1;
Amat4(5,1:2)   = -1;
Amat4(5,8)   = 1;
Amat4(6,1) = 1;
Amat4(6,10) = -1;
Amat4(7,1) = -1;
Amat4(7,10) = 1;
Amat4(8,3) = 1;
Amat4(8,9) = -1;
Amat4(9,3) = -1;
Amat4(9,9) = 1;
Amat4(10,2:3) = 1;
Amat4(10,11) = -1;
Amat4(11,2:3) = -1;
Amat4(11,11) = 1;
Amat4(12,1:3) = -gammaHat(1);
Amat4(13,1:3) = gammaHat(1);
Amat4(14,1:3) = -gammaHat(2);
Amat4(15,1:3) = gammaHat(2);
Amat4(16,1:3) = -gammaHat(3);
Amat4(17,1:3) = gammaHat(3);
Amat4(18,1:3) = -gammaHat(4);
Amat4(19,1:3) = gammaHat(4);
Amat4(12,4) = -1;
Amat4(12,8) = 1;
Amat4(13,4) = 1;
Amat4(13,8) = -1;
Amat4(14,5) = -1;
Amat4(14,9) = 1;
Amat4(15,5) = 1;
Amat4(15,9) = -1;
Amat4(16,6) = -1;
Amat4(16,10) = 1;
Amat4(17,6) = 1;
Amat4(17,10) = -1;
Amat4(18,7) = -1;
Amat4(18,11) = 1;
Amat4(19,7) = 1;
Amat4(19,11) = -1;

Amat4(20,1:3) = -1; 



Dmat34 = zeros(11);
Dmat34(4:7,4:7) = eye(4);
    
ndraw  = 1000;
T95 = zeros(1,length(psi));
CI = ones(1,length(psi))*99;
    
for j = 1:length(psi)
    rng('default');
    rng(j);
    temp  = [psi(j), psi(j)-1, psi(j),zeros(1,8)];
    Amat3 = [Amat4;temp;-temp];
    
    tdraws = zeros(1,ndraw);
    
    for i = 1:ndraw
        Z = mvnrnd(zeros(4,1),sigma);
        dvec34 = [zeros(1,3),Z,zeros(1,4)];
        
        %[~,val4] = cplexqp(Dmat34,dvec34,-Amat4,-bvec4);
        [~,val4] = quadprog(Dmat34,dvec34,-Amat4,-bvec4,[],[],[],[],[],options);
        %if exitflag ~= 1
        %    [~,val4,exitflag] = quadprog(Dmat34,dvec34,-Amat4,-bvec4,[],[],[],[],[],options);
        %    if exitflag ~= 1
        %        error(msg)
        %    end
        %end
        %[~,val3] = cplexqp(Dmat34,dvec34,-Amat3,-bvec3);
        [~,val3] = quadprog(Dmat34,dvec34,-Amat3,-bvec3,[],[],[],[],[],options);
        %if exitflag ~= 1
        %    [~,val3,exitflag] = quadprog(Dmat34,dvec34,-Amat3,-bvec3,[],[],[],[],[],options);
        %    if exitflag ~= 1
        %        error(msg)
        %    end
        %end
        tdraws(i) = val3 - val4;
      
        clear val3 val4 x3 x4 Z dvec34 
    end
    tdraws = sort(tdraws);
    T95(j) = tdraws(0.95*ndraw);
    %T95(j) = quantile(tdraws,0.95); %Segmentation (core dumped)
    %error on cluster
    
    bvec1 = [0; 0; 1-psi(j); psi(j)-1];
    dvec1 = 2 * [pC1 + pT1 - pC*psi(j), pC2 + pT2 - pT*psi(j)];
    %[~, val1] = cplexqp(Dmat1,-dvec1,-Amat1,-bvec1);
    [~, val1] = quadprog(Dmat1,-dvec1,-Amat1,-bvec1,[],[],[],[],[],options);
    %if exitflag ~= 1
    %    error(msg)
    %end
    val1 = val1 + 1 + psi(j)^2 - 2*(pC1 + pT2)*psi(j);
    
    Tn = n * (val1 - val2);

    if Tn <= T95(j)+1e-10
        break
    end
    
end

leftLim = psi(j);


for j = length(psi):-1:1
    rng('default');
    rng(j);
    temp  = [psi(j), psi(j)-1, psi(j),zeros(1,8)];
    Amat3 = [Amat4;temp;-temp];
    
    tdraws = zeros(1,ndraw);
    
    for i = 1:ndraw
        Z = mvnrnd(zeros(4,1),sigma);
        dvec34 = [zeros(1,3),Z,zeros(1,4)];
        
        %[~,val4] = cplexqp(Dmat34,dvec34,-Amat4,-bvec4);
        [~,val4] = quadprog(Dmat34,dvec34,-Amat4,-bvec4,[],[],[],[],[],options);
        %if exitflag ~= 1
        %    [~,val4,exitflag] = quadprog(Dmat34,dvec34,-Amat4,-bvec4,[],[],[],[],[],options);
        %    if exitflag ~= 1
        %        error(msg)
        %    end
        %end
        %[~,val3] = cplexqp(Dmat34,dvec34,-Amat3,-bvec3);
        [~,val3] = quadprog(Dmat34,dvec34,-Amat3,-bvec3,[],[],[],[],[],options);
        %if exitflag ~= 1
        %    [~,val3,exitflag] = quadprog(Dmat34,dvec34,-Amat3,-bvec3,[],[],[],[],[],options);
        %    if exitflag ~= 1
        %        error(msg)
        %    end
        %end
        tdraws(i) = val3 - val4;
      
        clear val3 val4 x3 x4 Z dvec34 
    end
    tdraws = sort(tdraws);
    T95(j) = tdraws(0.95*ndraw);
    %T95(j) = quantile(tdraws,0.95); %Segmentation (core dumped)
    %error on cluster
    
    bvec1 = [0; 0; 1-psi(j); psi(j)-1];
    dvec1 = 2 * [pC1 + pT1 - pC*psi(j), pC2 + pT2 - pT*psi(j)];
    %[~, val1] = cplexqp(Dmat1,-dvec1,-Amat1,-bvec1);
    [~, val1] = quadprog(Dmat1,-dvec1,-Amat1,-bvec1,[],[],[],[],[],options);
    %if exitflag ~= 1
    %    error(msg)
    %end
    val1 = val1 + 1 + psi(j)^2 - 2*(pC1 + pT2)*psi(j);
    
    Tn = n * (val1 - val2);

    if Tn <= T95(j)+1e-10
        break
    end
    
end

rightLim = psi(j);

CI = [leftLim, rightLim];

FileName = ['result',datestr(now, 'dd-mmm-yyyy'),'-',num2str(val),'.mat'];
save(FileName, 'CI')

exit


