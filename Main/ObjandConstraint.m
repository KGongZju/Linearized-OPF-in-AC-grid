function [Constraints,Obj1,Obj2,Obj3,PostiveBranchFlow,ReactiveBranchFlow]...
         =ObjandConstraint(PostiveBranchFlow,ReactiveBranchFlow,NodeVoltageS,NodeTheta,...
                         UnitP,GenOutP,OperationIF1,UnitQ,GenOutQ,...
                         InstallLocation,InstallCapacity,InstallPowerrating,...
                         BatteryDc,BatteryCh,BatteryIF,...
                         AGV1,AGV2,NodeIniVoltage,NodeIniTheta)
                                           
[Busdata,Gendata,branchdata,Gencostdata]=Data(); %读取电网数据


% 线路信息
LineI=branchdata(:,1);
LineJ=branchdata(:,2);
LineR=branchdata(:,3);
LineX=branchdata(:,4);
LineB=branchdata(:,5);
LineNum=length(LineI);

% 节点信息
NodeI=Busdata(:,1);
NodePl=Busdata(:,3)/100;
NodeQl=Busdata(:,4)/100;
NodeNum=length(NodeI);
% 发电机信息
GenI=Gendata(:,1);
GenNum=length(GenI);
GenPmin=Gendata(:,10)/100;
GenPmax=Gendata(:,9)/100;
GenQmin=Gendata(:,5)/100;
GenQmax=Gendata(:,4)/100;
GenC=Gencostdata(:,7);
GenB=Gencostdata(:,6);
GenA=Gencostdata(:,5);



%% OPF基础声明变量
% 支路有功变量声明，采用sparse稀疏矩阵形式，方便coding
% PositiveBFvarible=sdpvar(LineNum,1);
% PostiveBranchFlow=sparse([LineI;LineJ],[LineJ,LineI],[PositiveBFvarible;-PositiveBFvarible],NodeNum,NodeNum);
% IniPostiveBranchFlow=sparse([LineI;LineJ],[LineJ,LineI],[ones(LineNum,1);-ones(LineNum,1)],NodeNum,NodeNum);
% 支路无功变量声明，采用sparse稀疏矩阵形式，方便coding
% ReactiveBFvarible=sdpvar(LineNum,1);
% ReactiveBranchFlow=sparse([LineI;LineJ],[LineJ,LineI],[ReactiveBFvarible;-ReactiveBFvarible],NodeNum,NodeNum);
% IniReactiveBranchFlow=sparse([LineI;LineJ],[LineJ,LineI],[zeros(LineNum,1);zeros(LineNum,1)],NodeNum,NodeNum);
% 节点电压平方变量声明
% NodeVoltageS=sdpvar(NodeNum,1);
% NodeIniVoltage=ones(NodeNum,1);
% 节点相角变量声明
% NodeTheta=sdpvar(NodeNum,1);
% NodeIniTheta=zeros(NodeNum,1); %随迭代进行更改
% 声明发电机有功出力变量
% UnitP=sdpvar(GenNum,1);
% OperationIF=binvar(GenNum,1);
% GenOutP=sparse(GenI,ones(1,length(GenI)),UnitP,NodeNum);
% 声明发电机无功出力变量
% UnitQ=sdpvar(GenNum,1);
% GenOutQ=sparse(GenI,ones(1,length(GenI)),UnitQ,NodeNum);
%% 基于传统OPF基础上，对储能选址定容变量进行声明
% InstallLocation=binvar(NodeNum,1);
% InstallCapacity=sdpvar(NodeNum,1);
% InstallPowerrating=sdpvar(NodeNum,1);
% BatteryDc=sdpvar(NodeNum,1); %电池放电变量
% BatteryCh=sdpvar(NodeNum,1); %电池充电变量
% BatterySoc=sdpvar(NodeNum,1); %电池充电状态
% BatteryIF=binvar(NodeNum,1); %并网电池的充电状态约束
MaxCapacity=0; % 最大新建电池容量
MaxPowerrating=0; %最大新建电池额定有功功率
%% 形成节点导纳矩阵Ymatrix

% 形成节点导纳矩阵Ymatrix的对角阵
LineY=1./(LineR+1i*LineX); % 计算每条线路的导纳       
Y1=sparse([LineI;LineJ],[LineJ;LineI],[-LineY;-LineY],NodeNum,NodeNum);    
% 形成节点导纳矩阵Ymatrix的非对角阵
Y2=sparse([LineI;LineJ],[LineI;LineJ],[LineY+1i*LineB;LineY+1i*LineB],NodeNum,NodeNum);
AdmittanceMatrix=Y1+Y2;
% 形成每条线路的电导g和电纳b
Conductanceij=sparse([LineI;LineJ],[LineJ;LineI],[real(LineY);real(LineY)],NodeNum,NodeNum);   %电导g
Susceptanceij=sparse([LineI;LineJ],[LineJ;LineI],[imag(LineY);imag(LineY)],NodeNum,NodeNum);   %电纳b
% 形成相角差系数矩阵
DeffTheta=sparse([LineI;LineJ],[LineJ;LineI],[NodeTheta(LineI)-NodeTheta(LineJ);NodeTheta(LineJ)-NodeTheta(LineI)],NodeNum,NodeNum);
DeffIniTheta=sparse([LineI;LineJ],[LineJ;LineI],[NodeIniTheta(LineI)-NodeIniTheta(LineJ);NodeIniTheta(LineJ)-NodeIniTheta(LineI)],NodeNum,NodeNum);
% 形成初始节点电压系数矩阵及迭代电压系数矩阵
MultiVoltageij=sparse([LineI;LineJ],[LineJ;LineI],[NodeIniVoltage(LineI).*NodeIniVoltage(LineJ);NodeIniVoltage(LineJ).*NodeIniVoltage(LineI)],NodeNum,NodeNum);
VoltageijS=sparse([LineI;LineJ],[LineJ;LineI],...
           [2*(NodeIniVoltage(LineI)-NodeIniVoltage(LineJ))./(NodeIniVoltage(LineI)+NodeIniVoltage(LineJ)).*(NodeVoltageS(LineI)-NodeVoltageS(LineJ))-(NodeIniVoltage(LineI)-NodeIniVoltage(LineJ)).^2;...
            2*(NodeIniVoltage(LineJ)-NodeIniVoltage(LineI))./(NodeIniVoltage(LineI)+NodeIniVoltage(LineJ)).*(NodeVoltageS(LineJ)-NodeVoltageS(LineI))-(NodeIniVoltage(LineJ)-NodeIniVoltage(LineI)).^2],NodeNum,NodeNum);
% 形成每条线路的有功等效电导和等效电纳
EPConductanceij=Conductanceij.*cos(DeffIniTheta)+Susceptanceij.*sin(DeffIniTheta);
EPSusceptanceij=-Conductanceij.*MultiVoltageij.*sin(DeffIniTheta)+Susceptanceij.*MultiVoltageij.*cos(DeffIniTheta);
% 形成每条线路的无功等效电导和等效电纳
EQConductanceij=Conductanceij.*MultiVoltageij.*cos(DeffIniTheta)+Susceptanceij.*MultiVoltageij.*sin(DeffIniTheta);
EQSusceptanceij=-Conductanceij.*sin(DeffIniTheta)+Susceptanceij.*cos(DeffIniTheta);

%% *********构造约束*********
Constraints=[];
%% *********正常运行状态*********
%% 储能电池充放电约束（考虑选址）
M=100;
Constraints=[Constraints, 0<=InstallPowerrating(NodeI)<=InstallLocation(NodeI)*MaxPowerrating];  %是否配置电池，及配置电池的额定有功上限约束
Constraints=[Constraints, 0<=InstallCapacity(NodeI)<=InstallLocation(NodeI)*MaxCapacity];        %是否配置电池，及配置电池的容量上限约束
Constraints=[Constraints, 0<=BatteryDc(NodeI)<=InstallPowerrating(NodeI)];                       %电池充电有功约束
Constraints=[Constraints, 0<=BatteryCh(NodeI)<=InstallPowerrating(NodeI)];                       %电池放电有功约束
Constraints=[Constraints, 0<=BatteryDc(NodeI)<=M*BatteryIF(NodeI)];        %充放电约束          
Constraints=[Constraints, 0<=BatteryCh(NodeI)<=M*(1-BatteryIF(NodeI))];    %充放电约束
Kcharge=0.9;
Kdischarge=1.1;
Constraints=[Constraints,  Kcharge*BatteryCh(NodeI)-Kdischarge*BatteryDc(NodeI)<=InstallCapacity(NodeI)];
Constraints=[Constraints,  Kcharge*BatteryCh(NodeI)-Kdischarge*BatteryDc(NodeI)>=0];
%% 构造交流线路的Input node 和 Output node，方便引用

GetACBranchI=branchdata(:,1);
GetACBranchJ=branchdata(:,2);
EAV=ones(length(GetACBranchI),1); %扩展辅助矩阵，方便后续向量程序编写
% 构造交流线路有功无功潮流

PostiveBranchFlow(GetACBranchI,GetACBranchJ) =Conductanceij(GetACBranchI,GetACBranchJ).*NodeVoltageS(GetACBranchI,EAV)-...
                                              EPConductanceij(GetACBranchI,GetACBranchJ).*(NodeVoltageS(GetACBranchI,EAV)+NodeVoltageS(GetACBranchJ,EAV))/2-...
                                              EPSusceptanceij(GetACBranchI,GetACBranchJ).*(DeffTheta(GetACBranchI,GetACBranchJ)-DeffIniTheta(GetACBranchI,GetACBranchJ))+...
                                              EPConductanceij(GetACBranchI,GetACBranchJ).*VoltageijS(GetACBranchI,GetACBranchJ)/2;
PostiveBranchFlow(GetACBranchJ,GetACBranchI) =Conductanceij(GetACBranchJ,GetACBranchI).*NodeVoltageS(GetACBranchJ,EAV)-...
                                              EPConductanceij(GetACBranchJ,GetACBranchI).*(NodeVoltageS(GetACBranchI,EAV)+NodeVoltageS(GetACBranchJ,EAV))/2-...
                                              EPSusceptanceij(GetACBranchJ,GetACBranchI).*(DeffTheta(GetACBranchJ,GetACBranchI)-DeffIniTheta(GetACBranchJ,GetACBranchI))+...
                                              EPConductanceij(GetACBranchJ,GetACBranchI).*VoltageijS(GetACBranchJ,GetACBranchI)/2;

ReactiveBranchFlow(GetACBranchI,GetACBranchJ)=-Susceptanceij(GetACBranchI,GetACBranchJ).*NodeVoltageS(GetACBranchI,EAV)+...
                                              EQSusceptanceij(GetACBranchI,GetACBranchJ).*(NodeVoltageS(GetACBranchI,EAV)+NodeVoltageS(GetACBranchJ,EAV))/2-...
                                              EQConductanceij(GetACBranchI,GetACBranchJ).*(DeffTheta(GetACBranchI,GetACBranchJ)-DeffIniTheta(GetACBranchI,GetACBranchJ))-...
                                              EQSusceptanceij(GetACBranchI,GetACBranchJ).*VoltageijS(GetACBranchI,GetACBranchJ)/2;
ReactiveBranchFlow(GetACBranchJ,GetACBranchI)=-Susceptanceij(GetACBranchJ,GetACBranchI).*NodeVoltageS(GetACBranchJ,EAV)+...
                                              EQSusceptanceij(GetACBranchJ,GetACBranchI).*(NodeVoltageS(GetACBranchI,EAV)+NodeVoltageS(GetACBranchJ,EAV))/2-...
                                              EQConductanceij(GetACBranchJ,GetACBranchI).*(DeffTheta(GetACBranchJ,GetACBranchI)-DeffIniTheta(GetACBranchJ,GetACBranchI))-...
                                              EQSusceptanceij(GetACBranchJ,GetACBranchI).*VoltageijS(GetACBranchJ,GetACBranchI)/2;




%% 构造节点平衡方程-交流节点

GetACNodeI=Busdata(:,1);
for i=1:length(GetACNodeI)
    SumGij(GetACNodeI(i))=sum(real(AdmittanceMatrix((GetACNodeI(i)),:)));
    SumBij(GetACNodeI(i))=sum(imag(AdmittanceMatrix((GetACNodeI(i)),:)));
end
for i=1:length(GetACNodeI)
    if GetACNodeI(i)==1  %平衡节点独立约束
      Constraints=[Constraints,sum(UnitP)+sum(BatteryDc)-sum(BatteryCh)-SumGij(GetACNodeI)*NodeVoltageS(GetACNodeI)-sum(sum(PostiveBranchFlow))==sum(NodePl)];
      Constraints=[Constraints,sum(UnitQ)+SumBij(GetACNodeI)*NodeVoltageS(GetACNodeI)-sum(sum(ReactiveBranchFlow))==sum(NodeQl)];
%                 +ReactiveBranchFlow(ACDCconnectiondata(1,1),ACDCconnectiondata(1,2))+ReactiveBranchFlow(ACDCconnectiondata(2,1),ACDCconnectiondata(2,2))
%                 %********************************************************上面添加项存疑，需要讨论********************************************************%
      Constraints=[Constraints,  NodeVoltageS(GetACNodeI(i))<=1.05^2 , 0.95^2<=NodeVoltageS(GetACNodeI(i))];  
    else
    Corrlbranchij=SearchNodeConnection(LineI,LineJ,GetACNodeI(i)); %获取每个节点对应的支路信息
    InjectionACNodeP(GetACNodeI(i))=GenOutP(GetACNodeI(i))-NodePl(GetACNodeI(i))+BatteryDc(GetACNodeI(i))-BatteryCh(GetACNodeI(i));
    InjectionACNodeQ(GetACNodeI(i))=GenOutQ(GetACNodeI(i))-NodeQl(GetACNodeI(i));   
    SumCorrlBranchACP(GetACNodeI(i))=sum(PostiveBranchFlow(GetACNodeI(i),Corrlbranchij(:,2))); 
    SumCorrlBranchACQ(GetACNodeI(i))=sum(ReactiveBranchFlow(GetACNodeI(i),Corrlbranchij(:,2))); 
    % 节点有功平衡
    Constraints=[Constraints,...
        InjectionACNodeP(GetACNodeI(i))==SumCorrlBranchACP(GetACNodeI(i))+NodeVoltageS(GetACNodeI(i))*SumGij(GetACNodeI(i))
    ];

    % 节点无功平衡
    Constraints=[Constraints,...
        InjectionACNodeQ(GetACNodeI(i))==SumCorrlBranchACQ(GetACNodeI(i))-NodeVoltageS(GetACNodeI(i))*SumBij(GetACNodeI(i))
    ];

    % 节点电压约束
    Constraints=[Constraints,  NodeVoltageS(GetACNodeI(i))<=1.05^2 , 0.95^2<=NodeVoltageS(GetACNodeI(i))];  
    end
end



%% 发电机组约束集
% AGV1=sdpvar(GenNum,1); % anxillary generator varible +
% AGV2=sdpvar(GenNum,1); % anxillary generator varible -
Constraints=[Constraints, GenPmin.*OperationIF1<=UnitP<=GenPmax.*OperationIF1, GenQmin.*OperationIF1<=UnitQ<=GenQmax.*OperationIF1];
Constraints=[Constraints, UnitQ./(GenQmax-GenQmin)+AGV1-AGV2==0,  AGV1>=0, AGV2>=0];
% Constraints=[Constraints, -1<=PostiveBranchFlow(GetACBranchI,GetACBranchJ)<=1];
% Constraints=[Constraints, -0.5<=ReactiveBranchFlow(GetACBranchI,GetACBranchJ)<=0.5];
%% 交流线路潮流上下限约束集
Npart=20; %将上、下半圆截为20份线段
alpha=pi/6;
beta=(pi-2*alpha)/Npart;
M=Npart;  %上半圆份数
N=Npart;  %下半圆份数
Smax=2;
%**********上半圆**********%
KAU=zeros(Npart,1);KBU=zeros(Npart,1);XPAU=zeros(Npart,1);
YPAU=zeros(Npart,1);XPBU=zeros(Npart,1);YPBU=zeros(Npart,1);
%**********下半圆**********%
KAD=zeros(Npart,1); KBD=zeros(Npart,1);XPAD=zeros(Npart,1);
YPAD=zeros(Npart,1);XPBD=zeros(Npart,1);YPBD=zeros(Npart,1);

%**********上半圆**********%  
for i=1:M
    KAU(i)=tan(alpha+(i-1)*beta);
    KBU(i)=tan(alpha+i*beta);
    XPAU(i)=1/( sqrt( 1+( 1/(KAU(i))^2 ) )*KAU(i) ) * Smax;
    YPAU(i)=1/( sqrt( 1+( 1/(KAU(i))^2 ) ) ) * Smax;
    XPBU(i)=1/( sqrt( 1+( 1/(KBU(i))^2 ) )*KBU(i) ) * Smax;
    YPBU(i)=1/( sqrt( 1+( 1/(KBU(i))^2 ) ) ) * Smax;
    Constraints=[Constraints, ( ( YPBU(i)-YPAU(i) )/( XPBU(i)-XPAU(i) )*( ReactiveBranchFlow(GetACBranchI,GetACBranchJ) - XPAU(i) )+YPAU(i)-PostiveBranchFlow(GetACBranchI,GetACBranchJ) )>=0,... %交流支路
    ];  
end

%**********下半圆**********%

for i=1:N
    KAD(i)=tan(-alpha-(i-1)*beta);
    KBD(i)=tan(-alpha-i*beta);
    XPAD(i)=-1/( sqrt( 1+( 1/(KAD(i))^2 ) )*KAD(i) ) * Smax;
    YPAD(i)=-1/( sqrt( 1+( 1/(KAD(i))^2 ) ) ) * Smax;
    XPBD(i)=-1/( sqrt( 1+( 1/(KBD(i))^2 ) )*KBD(i) ) * Smax;
    YPBD(i)=-1/( sqrt( 1+( 1/(KBD(i))^2 ) ) ) * Smax;
    Constraints=[Constraints, ( ( YPBD(i)-YPAD(i) )/( XPBD(i)-XPAD(i) )*( ReactiveBranchFlow(GetACBranchI,GetACBranchJ) - XPAD(i) )+YPAD(i)-PostiveBranchFlow(GetACBranchI,GetACBranchJ) )<=0,...%交流支路
    ]; %交流支路
end

%% 目标成本函数
CPV=18.95;   %有功成本 百$/kW
CEV=9.01;    %容量成本 百$/kWh
OMC=0.132;   %运行成本 百$/kW
InstallCost=sum(CPV*InstallPowerrating)+sum(CEV*InstallCapacity);
OPM=sum(OMC*(BatteryDc)+OMC*(BatteryCh));
Obj1=InstallCost;
Obj2=GenA'*UnitP+sum(AGV1+AGV2)+OPM;
Obj3=0;
end
