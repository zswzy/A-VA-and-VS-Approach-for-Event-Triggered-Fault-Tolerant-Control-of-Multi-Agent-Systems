function DYN_FAULT = get_fault_dynamic()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

global A Ab B Bb C Cb Kb Lob Hb Bfb Cfb Pb Nab K D L
global N % number of followers
global n m r
global Mb Lvsb T U

D11 = Ab-Hb*Bb*Kb;
D12 = -Hb*Bb*Kb;
D13 = zeros(N*n);
D14 = zeros(N*n);
D21 = zeros(N*n);
D22 = Ab+Lob*Cb;
D23 = -Lob*(Cb-Pb*Cfb);
D24 = Lob*Cb;
D31 = -(Bb-Bfb)*Nab*Kb*Hb;
D32 = -(Bb-Bfb)*Nab*Kb*Hb;
D33 = Ab-Lvsb*Cfb;
D34 = (Bb-Bfb)*Mb;
D41 = -(Bb-Bfb*Nab)*Kb*Hb;
D42 = -(Bb-Bfb*Nab)*Kb*Hb;
D43 = zeros(N*n);
D44 = Ab-Bfb*Mb;

D15 = -Hb*Bb*Kb;
D25 = zeros(N*n);
D35 = -(Bb-Bfb)*Nab*Kb*Hb;
D45 = -(Bb-Bfb*Nab)*Kb*Hb;
D51 = Hb*Bb*Kb;
D52 = Hb*Bb*Kb - Lob*Cb;
D53 = Lob*(Cb-Pb*Cfb);
D54 = -Lob*Cb;
D55 = Ab+Hb*Bb*Kb;

DYN = [D11  D12  D13  D14;
       D21  D22  D23  D24;
       D31  D32  D33  D34;
       D41  D42  D43  D44];


DYN_FAULT = [ [DYN [D15;D25;D35;D45]]; [D51 D52 D53 D54 D55] ];
end

