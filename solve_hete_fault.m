function [flag, Mb, Lvsb, T, U] = solve_hete_fault(Bfb, Cfb)
%SOLVE_HETE_FAULT Solve for heterogeneous fault
% flag: 1: solve, 0: non

global A Ab B Bb C Cb Kb Lob Hb Pb Nab
global N % number of followers
global n m r



% define gains
Nab = zeros(m*N);
Pb = zeros(r*N);
Mb = zeros(m*N,n*N);
Lvsb = zeros(n*N,r*N);
for i = 1:N
    Bfi = Bfb(n*i-n+1:n*i, m*i-m+1:m*i);
    Nai = pinv(Bfi)*B;
    Nab(m*i-m+1:m*i, m*i-m+1:m*i) = Nai;
    Mi = place(A, Bfi, linspace(-10,-5,n));
    Mb(m*i-m+1:m*i, n*i-n+1:n*i) = Mi;


    Cfi = Cfb(r*i-r+1:r*i, n*i-n+1:n*i);
    Pi = C*pinv(Cfi);
    Pb(r*i-r+1:r*i, r*i-r+1:r*i) = Pi;
    Lvsi = place(A',Cfi',linspace(-10,-7,n))';
    Lvsb(n*i-n+1:n*i, r*i-r+1:r*i) = Lvsi;
end

Q1 = sdpvar(n);
Q2 = sdpvar(n);
Q3 = sdpvar(n);
Q4 = sdpvar(n);
Q14 = sdpvar(n, n, 'full');
Q24 = sdpvar(n, n, 'full');
Q34 = sdpvar(n, n, 'full');
T = sdpvar(n);
U = sdpvar(n);

Q1b = kron(eye(N), Q1);
Q2b = kron(eye(N), Q2);
Q3b = kron(eye(N), Q3);
Q4b = kron(eye(N), Q4);
Q14b = kron(eye(N), Q14);
Q24b = kron(eye(N), Q24);
Q34b = kron(eye(N), Q34);
Tb = kron(eye(N),T);
Ub = kron(eye(N),U);

S11 = Q1b*(Ab-Hb*Bb*Kb) + (Ab-Hb*Bb*Kb)'*Q1b ...
    -Q14b*(Bb-Bfb*Nab)*Kb*Hb-(Q14b*(Bb-Bfb*Nab)*Kb*Hb)' ...
    +Hb*Hb*Tb;
S12 = -Q1b*Hb*Bb*Kb ...
    -Q14b*(Bb-Bfb*Nab)*Kb*Hb-(Q24b*(Bb-Bfb*Nab)*Kb*Hb)' ...
    +Hb*Hb*Tb;
S13 = -(Q3b*(Bb-Bfb)*Nab*Kb*Hb)' ...
    -(Q34b*(Bb-Bfb*Nab)*Kb*Hb)';
S14 = -(Q4b*(Bb-Bfb*Nab)*Kb*Hb)' ...
    +Q14b*(Ab-Bfb*Mb) + Ab'*Q14b - Hb*Kb'*Bb'*Q14b - Kb'*Nab'*(Bb-Bfb)'*Hb*Q34b;
S15 = -Q1b*Hb*Bb*Kb ...
    -Q14b*(Bb-Bfb*Nab)*Kb*Hb ...
    +Hb*Hb*Tb;

S22 = Q2b*(Ab+Lob*Cb) + (Ab+Lob*Cb)'*Q2b ...
    -Q24b*(Bb-Bfb*Nab)*Kb*Hb - (Q24b*(Bb-Bfb*Nab)*Kb*Hb)'...
    +Hb*Hb*Tb;
S23 = -Q2b*Lob*(Cb-Pb*Cfb) - (Q3b*(Bb-Bfb)*Nab*Kb*Hb)' ...
    -(Q34b*(Bb-Bfb*Nab)*Kb*Hb)';
S24 = Q2b*Lob*Cb - (Q4b*(Bb-Bfb*Nab)*Kb*Hb)' ...
    -Hb*Kb'*Bb'*Q14b + Q24b*(Ab-Bfb*Mb) + (Ab+Lob*Cb)'*Q24b - Kb'*Nab'*(Bb-Bfb)'*Hb*Q34b;
S25 = zeros(N*n) ...
    -Q24b*(Bb-Bfb*Nab)*Kb*Hb ...
    +Hb*Hb*Tb;

S33 = Q3b*(Ab-Lvsb*Cfb) + (Ab-Lvsb*Cfb)'*Q3b;
S34 = Q3b*(Bb-Bfb)*Mb ...
    -(Cb-Pb*Cfb)'*Lob'*Q24b + Q34b*(Ab-Bfb*Mb) + (Ab-Lvsb*Cfb)'*Q34b;
S35 = -Q3b*(Bb-Bfb)*Nab*Kb*Hb ...
    -Q34b*(Bb-Bfb*Nab)*Kb*Hb;

S44 = Q4b*(Ab-Bfb*Mb) + (Ab-Bfb*Mb)'*Q4b ...
    -Cb'*Lob'*Q24b - (Cb'*Lob'*Q24b)' + Mb'*(Bb-Bfb)'*Q34b + (Mb'*(Bb-Bfb)'*Q34b)';
S45 = -Q4b*(Bb-Bfb*Nab)*Kb*Hb ...
    -Q14b'*Bb*Kb*Hb - Q34b'*(Bb-Bfb)*Nab*Kb*Hb;

S55 = Hb*Hb*Tb-Ub;

S = [S11  S12  S13  S14  S15;
     S12' S22  S23  S24  S25;
     S13' S23' S33  S34  S35;
     S14' S24' S34' S44  S45;
     S15' S25' S35' S45' S55];

PP = [Q1b                 zeros(N*n)    zeros(N*n)  Q14b;
      zeros(N*n)          Q2b           zeros(N*n)  Q24b;
      zeros(N*n)          zeros(N*n)    Q3b         Q34b;
      Q14b'               Q24b'         Q34b'       Q4b];

constraints = [S<=0, PP>=0, T>=0, U>=0];
options = sdpsettings('solver','sdpt3','verbose',1);
optimize(constraints,[],options)
[pri, ~] = check(constraints)

T = value(T);
U = value(U);

if (min(pri) > 0)
    flag = 1;
else
    flag = 0;
end

end

