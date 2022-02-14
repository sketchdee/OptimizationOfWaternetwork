%处理约束

m=12;
n=15;

FL=[53.3;2176;0;0;0;0;0;59;385;45;102;0];

%简化模型，不自身回用，两个末端单元，通过等式约束来控制。
Aeq1=[];
beq1=[];
for i=1:1:m
    M1=zeros(m,n);
    M1(i,i)=1;
    Aeq1=[Aeq1;Tool.M2V(M1)];
    beq1=[beq1;0];
end

%终端单元
for i=1:1:n
    M1=zeros(m,n);
    M1(10,i)=1;
    Aeq1=[Aeq1;Tool.M2V(M1)];
    beq1=[beq1;0];
end
for i=1:1:n
    M1=zeros(m,n);
    M1(11,i)=1;
    Aeq1=[Aeq1;Tool.M2V(M1)];
    beq1=[beq1;0];
end

%两个处理单元
for i=1:1:m
    if(i~=4)
        M1=zeros(m,n);
        M1(i,5)=1;
        Aeq1=[Aeq1;Tool.M2V(M1)];
        beq1=[beq1;0];
    end
    if(i==4)
        M1=zeros(m,n);
        M1(4,5)=-1/0.432;
        M1(1:9,4)=ones(9,1);
        Aeq1=[Aeq1;Tool.M2V(M1)];
        beq1=[beq1;0];
    end
end

for i=1:1:m
    if(i~=4)
        M1=zeros(m,n);
        M1(i,7)=1;
        Aeq1=[Aeq1;Tool.M2V(M1)];
        beq1=[beq1;0];
    end
    if(i==4)
        M1=zeros(m,n);
        M1(6,7)=-1/0.008;
        M1(1:9,6)=ones(9,1);
        Aeq1=[Aeq1;Tool.M2V(M1)];
        beq1=[beq1;0];
    end
end

Cw=[18.9,0,133];
%新鲜水负荷
M1(12,13)=1;
Aeq1=[Aeq1;Tool.M2V(M1)];
beq1=[beq1;Cw(1)];
M1(12,14)=1;
Aeq1=[Aeq1;Tool.M2V(M1)];
beq1=[beq1;Cw(2)];
M1(12,15)=1;
Aeq1=[Aeq1;Tool.M2V(M1)];
beq1=[beq1;Cw(3)];


%流量守恒约束
Aeq2=[];
beq2=[];
for i=1:1:m
    M1=zeros(m,n);
    M1(i,:)=ones(1,n).*(-1);
    M1(:,i)=ones(m,1);
    M1(i,i)=0;
    Aeq2=[Aeq2;Tool.M2V(M1)];
    beq2=[beq2;FL(i,1)]; 
end


%不等式约束
CoutMax=[300	30	1000
1600	100	28000
1600	100	28000
200	1000	0
300	0	28000
0.1	0	500
1000	0	500
10	0	2000
20000	20000	0
];
ub=ones(m,n).*500;
ub(1:9,13:15)=CoutMax;
ub=Tool.M2V(ub);
lb=Tool.M2V(zeros(m,n));

Aeq=[Aeq1;Aeq2];
beq=[beq1;beq2];

x0=ones(1,m*n);
options = optimoptions('fmincon');
options = optimoptions(options,'Display', 'iter',...
    'FunctionTolerance',1e-100,...
    'MaxIterations', 100000,...
    'Algorithm', 'active-set ',...
    'FiniteDifferenceType', 'central',...
    'StepTolerance',1e-1000);
options.MaxFunctionEvaluations = 100000;

[x,fval,exitflag,output] = fmincon(@ObjectFunction,x0,[],[],Aeq,beq,lb,ub,@nonlcon);
for i=1:1:100
    x0=x;
    [x,fval,exitflag,output] = fmincon(@ObjectFunction,x0,[],[],Aeq,beq,lb,ub,@nonlcon);

end

%x = ga(@ObjectFunction,180,[],[],Aeq,beq,lb,ub,nonlcon)
%result=reshape(x,[m,n])

%杂质守恒约束
%由于初始的出口浓度为未知的，所以只能通过联立方程组来求解
%非线性等式约束
function [c,ceq]=nonlcon(x)
m=12;
n=15;
%M=ones(12,15);
%M=Tool.V2M(x',m,n);
M=reshape(x,[12,15]);
MF=[
314.28	237.6	0
76750.32	20220	265165
15445.5	0	89355
22880	110000	-216480
22880	-2850	383040
-33737.5	0	33000
2149.9	0	264
244	0	78336
285263.3	65877	-283077
];
CinMax=[
100	0	300
300	10	2000
300	10	2000
200	500	2000
200	300	0
300	0	300
0.1	0	300
0.1	0	500
1000	1000	2000
];
for i=1:1:9
     ceq(i)=M(1:9,i)'*diag(ones(1,9)*(-1))*M(1:9,13) + M(i,1:12)*ones(12,1)*M(i,13)+M(12,i)*M(12,13)-MF(i,1);
% ceq(1+3*(i-1))=M(1:9,i)'*diag(ones(1,9)*(-1))*M(1:9,13) + M(i,1:12)*ones(12,1)*M(i,13)+M(12,i)*M(12,13)-MF(i,1);
% ceq(2+3*(i-1))=M(1:9,i)'*diag(ones(1,9)*(-1))*M(1:9,14) + M(i,1:12)*ones(12,1)*M(i,14)+M(12,i)*M(12,14)-MF(i,2);
% ceq(3+3*(i-1))=M(1:9,i)'*diag(ones(1,9)*(-1))*M(1:9,15) + M(i,1:12)*ones(12,1)*M(i,15)+M(12,i)*M(12,15)-MF(i,3);
end

for i=1:1:9
    c(i)= M(1:9,i)'*(M(1:9,13)-CinMax(1:9,1)) +M(12,i)*M(12,13);
%c(1+3*(i-1))= M(1:9,i)'*(M(1:9,13)-CinMax(1:9,1)) +M(12,i)*M(12,13);
%c(2+3*(i-1))= M(1:9,i)'*(M(1:9,13)-CinMax(1:9,2)) +M(12,i)*M(12,14);
%c(3+3*(i-1))= M(1:9,i)'*(M(1:9,13)-CinMax(1:9,3)) +M(12,i)*M(12,15);
end
end

function [Z]=ObjectFunction(x)
    Z=x(166)+x(167)+x(168)+x(169)+x(170)+x(171)+x(172)+x(173)+x(174)+x(175)+x(176);

end
