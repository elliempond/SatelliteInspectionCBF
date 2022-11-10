function [sol,tout] = rendezvousCBFPart1(chase,mc,n,Rt,acc,cone)
tic;
% Create prog
prog = spotsosprog;
z = msspoly('z',6*chase);
prog = prog.withIndeterminate(z); 


% Target Barrier
for L = 1:chase
mcCur = mc(1,L);
zc = z(6*L-5:6*L);
mx = monomials(zc(1:3),0:4);
    
% DSOS
switch cone
    case 'DSOS'
        [prog,coef_s1] = prog.newDD(length(mx));
        s1 = [mx]'*coef_s1*[mx];
        [prog,coef_s2] = prog.newDD(length(mx));
        s2 = [mx]'*coef_s2*[mx];
    case 'SOS'
        [prog,coef_s1] = prog.newPSD(length(mx));
        s1 = [mx]'*coef_s1*[mx];
        [prog,coef_s2] = prog.newPSD(length(mx));
        s2 = [mx]'*coef_s2*[mx];
end

% Polys
[prog,coef_p10] = prog.newFree(length(mx)); 
p10 = coef_p10'*mx;  
[prog,coef_p20] = prog.newFree(length(mx)); 
p20 = coef_p20'*mx;  
[prog,coef_p1a] = prog.newFree(length(mx)); 
p1a = coef_p1a'*mx;
[prog,coef_p1b] = prog.newFree(length(mx)); 
p1b = coef_p1b'*mx;
[prog,coef_p1c] = prog.newFree(length(mx)); 
p1c = coef_p1c'*mx;
p1 = [p1a; p1b; p1c];
[prog,coef_p2a] = prog.newFree(length(mx)); 
p2a = coef_p2a'*mx;
[prog,coef_p2b] = prog.newFree(length(mx)); 
p2b = coef_p2b'*mx;
[prog,coef_p2c] = prog.newFree(length(mx)); 
p2c = coef_p2c'*mx;
p2 = [p2a; p2b; p2c];


% Expressions
a = 1;
b = zc(1:3)'*zc(1:3) + (mcCur/acc)*zc(4:6)'*zc(4:6) - Rt^2;
f = [0 0 0 1 0 0; 0 0 0 0 1 0; 0 0 0 0 0 1; 0 0 0 0 2*n 0; 0 0 0 -2*n 0 0; 0 0 -n^2 0 0 0]*[zc] +[3*n^2; 0; 0; 0; 0; 0];
g = [0 0 0; 0 0 0; 0 0 0; 1/mcCur 0 0; 0 1/mcCur 0; 0 0 1/mcCur];
delbdelz = [2*zc(1),2*zc(2),2*zc(3),2*(mcCur/acc)*zc(4),2*(mcCur/acc)*zc(5),2*(mcCur/acc)*zc(6)];
Lfb = delbdelz*f;
Lgb = delbdelz*g;

h1 = s1 + b*p10 + Lgb*p1;
h2 = s2 + b*p20 + Lgb*p2;

expr = h1*Lfb - h2 - (Lfb)^(2*a);
prog = prog.withPolyEqs(expr);
end

% Solve
sol = prog.minimize(0);
tout = toc;

end
