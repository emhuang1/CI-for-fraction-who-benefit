%Code to compute one-sided CI for the upper bound

clear;

job_id_string = getenv('SGE_TASK_ID');
jobval = str2double(job_id_string);

addpath('/users/emhuang/cvx')

cvx_setup

addpath('/users/emhuang/CIpropBenefit/binary_noRes_5050/testDatasets')

data = importdata('N200.txt');

pc1 = data(jobval,3)/data(jobval,1);
pt1 = data(jobval,5)/data(jobval,2);


n  = 200;
theta = 0.5;

p = 0.975;

R = 1000;

V = [1,1,1;1,0,1;0,1,1;0,0,1];
b = [0;1;0;0];

%%%%%%%%%%%%%%%%Get upper limit of CI%%%%%%%%%%%%%%%%%%%%%%%%

theta1 = [pc1,pt1,1];

%Step 1
gamma = 1 - 0.1/log(n);
Zr = mvnrnd([0,0], eye(2), R);

%Step 2
Omega  = [pc1*(1-pc1)/(1-theta),0;0,pt1*(1-pt1)/theta];
Omegat = [Omega,zeros(2,1);zeros(1,3)];
Omegas = sqrt(Omegat);

%Step 4: Finding kappa
kappaV = chi2inv(gamma,2);
kappaV = sqrt(kappaV);

%Step 4: Finding Vn

cvx_begin
variable v(3)
minimize theta1*v + kappaV*norm(Omegas*v)/sqrt(n)
subject to 
     V*v >= b;
cvx_end

m = cvx_optval;

angle = pi/180 * (1:360);
x = cos(angle);
y = sin(angle);
angleLen = length(angle);

marker = zeros(1,angleLen); 

for i = 1:angleLen
	u1 = x(i);
	u2 = y(i);
	cons = pc1 * u1 + pt1 * u2 - 2 * kappaV/sqrt(n) * sqrt(pc1*(1-pc1)/(1-theta) * u1^2 + pt1*(1-pt1)/theta * u2^2);	
	
	val = 1;
	
	if u1 > 0 
		a = 1/u1;
		val2 = max([cons * a, cons * a - u2 * a, cons * a + 1 - u1 * a, cons * a - u1 * a - u2 * a]);
		if val2 < val
			val = val2;
		end
	end
	
	if u2 < 0
		a = -1/u2;
		val2 = max([cons * a, cons * a - u2 * a, cons * a + 1 - u1 * a, cons * a - u1 * a - u2 * a]);
		if val2 < val
			val = val2;
		end
	end
	
	if u1 > u2
		a = 1/(u1-u2);
		val2 = max([cons * a, cons * a - u2 * a, cons * a + 1 - u1 * a, cons * a - u1 * a - u2 * a]);
		if val2 < val
			val = val2;
		end
	end	
	
	limit_slopes = [cons, cons - u2, cons - u1, cons - u1 - u2];
	
	if sum(limit_slopes < 0) == 4
		marker(i) = 1;
	elseif sum(limit_slopes > 0) > 0
		val = val;
	else
		val2 = 0;
		if val2 < val
			val = val2;
		end
	end
	
	if val <= m
		marker(i) = 1;
	end
end


nV = sum(marker == 1);
Vn = zeros(nV, 2);

i = 1;

for j = 1:angleLen
	if marker(j) == 1
		Vn(i,:) = [x(j), y(j)];
		i = i + 1;
	end
end

%Step 5: Finding kappa

val = zeros(R,1);

for r = 1:R
	val2 = zeros(nV,1);
	for j = 1:nV
		v = Vn(j,:);
		gv = v*sqrt(Omega);
		val2(j) = gv*transpose(Zr(r,:))/sqrt(gv(1)^2 + gv(2)^2);
	end
	val(r) = max(val2);
end

kappaVn = quantile(val, p);

%Step 5: Finding theta_n0

cvx_begin
variable v(3)
minimize theta1*v + kappaVn*norm(Omegas*v)/sqrt(n)
subject to 
     V*v >= b;
cvx_end

theta_n0 = cvx_optval;


FileName = ['result',datestr(now, 'dd-mmm-yyyy'),'-',num2str(jobval),'.mat'];
save(FileName, 'theta_n0')

exit
