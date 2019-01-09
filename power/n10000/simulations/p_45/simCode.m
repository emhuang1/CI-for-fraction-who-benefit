clear;

job_id_string = getenv('SGE_TASK_ID');
val = str2double(job_id_string);

%addpath('/users/mrosen/CPLEX_STUDIO/cplex/matlab')
addpath('/users/emhuang/CIpropBenefit/power/binary_varyMarginals/n10000/testDatasets')

p = 0.45;
n = 10000;

options = optimoptions('quadprog','Algorithm','interior-point-convex','Display','off','TolFun',1e-10,'TolX',1e-10);
%msg = 'Problem with solving quadratic program.';

data = importdata('p_45.txt');
pC = data(val,1);
pT = data(val,2);
pC1 = data(val,3);
pC2 = data(val,4);
pT1 = data(val,5);
pT2 = data(val,6);

lb = max([0, 1 - p - 0.5]);
ub = min([0.5, 1 - p]);

psi = [lb - 1/sqrt(n), lb - 2/sqrt(n), lb - 3/sqrt(n), ub + 1/sqrt(n), ub + 2/sqrt(n), ub + 3/sqrt(n)];


Amat1 = [eye(3); [1,1,1]; [-1,-1,-1]];
Dmat1 = [[1,pT,0]; [pT,1,pC]; [0,pC,1]];
Dmat1 = Dmat1 * 2;


bvec2 = [0,0,0,0,1,-1,1,-1]';
Amat2 = [eye(4);[1,1,0,0];[-1,-1,0,0];[0,0,1,1];[0,0,-1,-1]];
dvec2 = 2*[pC1, pC2, pT1, pT2];
Dmat2 = zeros(4);
Dmat2(1,1) = 2*pC;
Dmat2(2,2) = 2*pC;
Dmat2(3,3) = 2*pT;
Dmat2(4,4) = 2*pT;    
%[gammaHat,val2] = cplexqp(Dmat2,-dvec2,-Amat2,-bvec2);
[gammaHat,val2] = quadprog(Dmat2,-dvec2,-Amat2,-bvec2,[],[],[],[],[],options);
%if exitflag ~= 1
%    error(msg)
%end
val2 = val2 + 1;


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
    
    
bvec3 = zeros(22,1);
bvec4 = zeros(20,1);
    
Amat4 = zeros(20,12);
Amat4(1:4,1:4) = eye(4);
Amat4(5,1:2)   = 1;
Amat4(6,1:2)   = -1;
Amat4(7,3:4)   = 1;
Amat4(8,3:4)   = -1;
Amat4(9,[1,3]) = 1;
Amat4(10,[1,3])= -1;
Amat4(11,[2,4])   = 1;
Amat4(12,[2,4])   = -1;    
Amat4(5,9)   = -1;
Amat4(7,10)  = -1;
Amat4(6,9)   = 1;
Amat4(8,10)  = 1;
Amat4(9,11)  = -1;
Amat4(11,12) = -1;
Amat4(10,11) = 1;
Amat4(12,12) = 1;
Amat4(13,1:4) = -gammaHat(1);
Amat4(14,1:4) = gammaHat(1);
Amat4(15,1:4) = -gammaHat(2);
Amat4(16,1:4) = gammaHat(2);
Amat4(17,1:4) = -gammaHat(3);
Amat4(18,1:4) = gammaHat(3);
Amat4(19,1:4) = -gammaHat(4);
Amat4(20,1:4) = gammaHat(4);    
Amat4(13,5)  = -1;
Amat4(14,9)  = -1;
Amat4(13,9)  = 1;
Amat4(14,5)  = 1;
Amat4(15,6)  = -1;
Amat4(16,10) = -1;
Amat4(15,10) = 1;
Amat4(16,6)  = 1;
Amat4(17,7)  = -1;
Amat4(18,11) = -1;
Amat4(17,11) = 1;
Amat4(18,7)  = 1;
Amat4(19,8)  = -1;
Amat4(20,12) = -1;
Amat4(19,12) = 1;
Amat4(20,8)  = 1;   
Dmat34 = zeros(12);
Dmat34(5:8,5:8) = eye(4);
    
ndraw  = 1000;
T95 = zeros(1,length(psi));
CI = ones(1,length(psi))*99;

for j = 1:length(psi)
    rng('default');
    rng(j);
    temp  = [psi(j), psi(j)-1, psi(j), psi(j),zeros(1,8)];
    Amat3 = [Amat4;temp;-temp];
    
    tdraws = zeros(1,ndraw);
    
    for i = 1:ndraw
        Z = mvnrnd(zeros(4,1),sigma);
        dvec34 = [zeros(1,4),-Z,zeros(1,4)];
        
        %[~,val4] = cplexqp(Dmat34,-dvec34,-Amat4,-bvec4);
        [~,val4] = quadprog(Dmat34,-dvec34,-Amat4,-bvec4,[],[],[],[],[],options);
        %if exitflag ~= 1
        %    [~,val4,exitflag] = quadprog(Dmat34,-dvec34,-Amat4,-bvec4,[],[],[],[],[],options);
        %    if exitflag ~= 1
        %        error(msg)
        %    end
        %end
        %[~,val3] = cplexqp(Dmat34,-dvec34,-Amat3,-bvec3);
        [~,val3] = quadprog(Dmat34,-dvec34,-Amat3,-bvec3,[],[],[],[],[],options);
        %if exitflag ~= 1
        %    [~,val3,exitflag] = quadprog(Dmat34,-dvec34,-Amat3,-bvec3,[],[],[],[],[],options);
        %    if exitflag ~= 1
        %        error(msg)
        %    end
        %end
        tdraws(i) = val3 - val4;
        clear val3 val4 x3 x4 Z dvec34 
    end
    tdraws = sort(tdraws);
    T95(j) = tdraws(0.95*ndraw);
    
    bvec1 = [0;0;0;(1-psi(j));(psi(j)-1)];
    dvec1 = 2 * [pC1 + pT1 - pC * psi(j), pC2 + pT1, pC2 + pT2 - pT * psi(j)];
    %[~,val1] = cplexqp(Dmat1, -dvec1, -Amat1, -bvec1);
    [~,val1] = quadprog(Dmat1, -dvec1, -Amat1, -bvec1,[],[],[],[],[],options);
    %if exitflag ~= 1
    %    error(msg)
    %end
    val1 = val1 + 1 + psi(j)^2 - 2 * (pC1 + pT2)*psi(j);
    
    Tn = n * (val1 - val2);
    
    CI(j) = (Tn <= T95(j)+1e-10); 
    clear tdraws Amat3 temp Tn
end
FileName = ['result',datestr(now, 'dd-mmm-yyyy'),'-',num2str(val),'.mat'];
save(FileName, 'CI', 'T95')

exit

