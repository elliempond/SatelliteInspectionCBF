clc; clear; 
 

%% Rendevzvous barrier

tic;

cone = 'SOS';
chase = 1;

mcmax = 2;
mcmin = mcmax;
mc = mcmin + (mcmax-mcmin).*rand(1,chase);

n = .0010;
Rt = .5;
acc = .5/1000;



% Part 1
[sol,tout] = rendezvousCBFPart1(chase,mc,n,Rt,acc,cone);
compile(1,1) = sol.info.solverInfo.pinf;
compile(2,1) = sol.info.solverInfo.feasratio;
compile(3,1) = sol.info.solverInfo.numerr;
compile(4,1) = tout;

% Part 2
if chase >1
[sol2,tout2] = rendezvousCBFPart2(chase,mc,n,Rt,acc,cone);
compile(1,2) = sol2.info.solverInfo.pinf;
compile(2,2) = sol2.info.solverInfo.feasratio;
compile(3,2) = sol2.info.solverInfo.numerr;
compile(4,2) = tout2;
end


