function [sol,tout] = rendezvousCBFPart2(chase,mc,n,Rt,acc,cone)
tic;
prog = spotsosprog;
z = msspoly('z',6*chase);
prog = prog.withIndeterminate(z); 

mx = [];
for L = 1:chase
    zc = z(6*L-5:6*L);
    if L ==1
        temp1 = monomials(zc(1:3),0:2);
        temp2 = monomials(zc(4:6),1:2);
        temp2 =[];
    else
        temp1 = monomials(zc(1:3),1:2);
        temp2 = monomials(zc(4:6),1:1);
        temp2 =[];
    end 
    mx = [mx; temp1; temp2];
end

switch cone
    case 'DSOS'
        [prog,coef_s0] = prog.newDD(length(mx));
        s0 = [mx]'*coef_s0*[mx];
        for L = 1:chase
            [prog,coef_s] = prog.newDD(length(mx));
            s(L,1) = [mx]'*coef_s*[mx];
        end
    case 'SOS'
        [prog,coef_s0] = prog.newPSD(length(mx));
        s0 = [mx]'*coef_s0*[mx];
        for L = 1:chase
            [prog,coef_s] = prog.newPSD(length(mx));
            s(L,1) = [mx]'*coef_s*[mx];
        end
end

% Expressions
for L = 1:chase
zc = z(6*L-5:6*L);
mcCur = mc(1,L);
b(1,L) = zc(1:3)'*zc(1:3) + (mcCur/acc)*zc(4:6)'*zc(4:6) - Rt^2;
end

expr = 1 + s0 + b*s;
prog = prog.withPolyEqs(expr);


% Solve
sol = prog.minimize(0);
tout = toc;
end