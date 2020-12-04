function [PostiveBranchFlow,ReactiveBranchFlow,NodeVoltageS,DNodeVoltageS,NodeTheta,BatteryCh,BatteryDc,DBatteryCh,DBatteryDc,VIMultiDVI]=...
          FormulationDeclear(IniPostiveBranchFlow,IniReactiveBranchFlow,NodeIniVoltage,NodeIniTheta,DNodeIniVoltage)

      
DNodeIniVoltage=0; %拓展数据接口，暂不用
DNodeVoltageS=0;   %拓展数据接口，暂不用
DBatteryCh=0;      %拓展数据接口，暂不用
DBatteryDc=0;      %拓展数据接口，暂不用
VIMultiDVI=0;      %拓展数据接口，暂不用

[Busdata,Gendata,branchdata,Gencostdata]=Data(); %读取电网数据
LineI=branchdata(:,1);
LineJ=branchdata(:,2);
LineNum=length(LineI);
% 节点信息
NodeI=Busdata(:,1);
NodeNum=length(NodeI);
% 发电机信息
GenI=Gendata(:,1);
GenNum=length(GenI);


%% OPF基础声明变量
% 支路有功变量声明，采用sparse稀疏矩阵形式，方便coding
PositiveBFvarible=sdpvar(LineNum,1);
PostiveBranchFlow=sparse([LineI;LineJ],[LineJ,LineI],[PositiveBFvarible;-PositiveBFvarible],NodeNum,NodeNum);

% 支路无功变量声明，采用sparse稀疏矩阵形式，方便coding
ReactiveBFvarible=sdpvar(LineNum,1);
ReactiveBranchFlow=sparse([LineI;LineJ],[LineJ,LineI],[ReactiveBFvarible;-ReactiveBFvarible],NodeNum,NodeNum);

% 节点电压平方变量声明
NodeVoltageS=sdpvar(NodeNum,1);

% 节点相角变量声明
NodeTheta=sdpvar(NodeNum,1);

% 声明发电机有功出力变量
UnitP=sdpvar(GenNum,1);
OperationIF1=binvar(GenNum,1);
OperationIF2=binvar(GenNum,1);
GenOutP=sparse(GenI,ones(1,length(GenI)),UnitP,NodeNum);

% 声明发电机无功出力变量
UnitQ=sdpvar(GenNum,1);
GenOutQ=sparse(GenI,ones(1,length(GenI)),UnitQ,NodeNum);

%% 基于传统OPF基础上，对储能选址定容变量进行声明
InstallLocation=binvar(NodeNum,1);
InstallCapacity=sdpvar(NodeNum,1);
InstallPowerrating=sdpvar(NodeNum,1);

BatteryDc=sdpvar(NodeNum,1); %电池放电变量
BatteryCh=sdpvar(NodeNum,1); %电池充电变量
BatteryIF=binvar(NodeNum,1); %并网电池的充电状态约束

%% 发电机组约束集
AGV1=sdpvar(GenNum,1); % anxillary generator varible +
AGV2=sdpvar(GenNum,1); % anxillary generator varible -

% %% *********故障运行状态*********
% DPositiveBFvarible=sdpvar(LineNum,1);
% DPostiveBranchFlow=sparse([LineI;LineJ],[LineJ,LineI],[DPositiveBFvarible;-DPositiveBFvarible],NodeNum,NodeNum);
% DReactiveBFvarible=sdpvar(LineNum,1);
% DReactiveBranchFlow=sparse([LineI;LineJ],[LineJ,LineI],[DReactiveBFvarible;-DReactiveBFvarible],NodeNum,NodeNum);

% % 系统频率变化
% DeltaF=sdpvar(1,1);

% % 故障状态节点电压平方变量声明 首字母D表示delta，表示变化量
% DNodeVoltageS=sdpvar(NodeNum,1);

% % 故障状态节点相角变量声明
% DNodeTheta=sdpvar(NodeNum,1);

% % 有功甩负荷
% CurtailP=sdpvar(NodeNum,1);

% % 无功甩负荷
% CurtailQ=sdpvar(NodeNum,1);

% % 声明发电机故障有功出力变化量
% DUnitP=sdpvar(GenNum,1);
% DGenOutP=sparse(GenI,ones(1,length(GenI)),DUnitP,NodeNum);

% % 电池充放电变化量声明
% DBatteryDc=sdpvar(NodeNum,1); %电池放电变量
% DBatteryCh=sdpvar(NodeNum,1); %电池充电变量

% DNodePl=sdpvar(NodeNum,1);
% DNodeQl=sdpvar(NodeNum,1);

% %电压变化量
% VIMultiDVI=zeros(NodeNum,1);
%% 调用构造函数
% IniPostiveBranchFlow=sparse([LineI;LineJ],[LineJ,LineI],[ones(LineNum,1);-ones(LineNum,1)],NodeNum,NodeNum);
% IniReactiveBranchFlow=sparse([LineI;LineJ],[LineJ,LineI],[zeros(LineNum,1);zeros(LineNum,1)],NodeNum,NodeNum);
% NodeIniVoltage=ones(NodeNum,1);
% NodeIniTheta=zeros(NodeNum,1);
[Constraints,Obj1,Obj2,Obj3,PostiveBranchFlow,ReactiveBranchFlow]...
         =ObjandConstraint(PostiveBranchFlow,ReactiveBranchFlow,NodeVoltageS,NodeTheta,...
                         UnitP,GenOutP,OperationIF1,UnitQ,GenOutQ,...
                         InstallLocation,InstallCapacity,InstallPowerrating,...
                         BatteryDc,BatteryCh,BatteryIF,...
                         AGV1,AGV2,NodeIniVoltage,NodeIniTheta);



OBJ=Obj1+Obj2+Obj3;
optimize(Constraints,OBJ);

NodeTheta=double(NodeTheta);
NodeVoltageS=double(NodeVoltageS);
PostiveBranchFlow=double(PostiveBranchFlow); ReactiveBranchFlow=double(ReactiveBranchFlow);


BatteryDc=double(BatteryDc);
BatteryCh=double(BatteryCh);


AGV1=double(AGV1); AGV2=double(AGV2);
NodeTheta=double(NodeTheta);
GenOutP=double(GenOutP); GenOutQ=double(GenOutQ);
NodeVoltageS=double(NodeVoltageS); UnitP=double(UnitP); UnitQ=double(UnitQ); OperationIF1=double(OperationIF1); OperationIF2=double(OperationIF2);
PostiveBranchFlow=double(PostiveBranchFlow); ReactiveBranchFlow=double(ReactiveBranchFlow);
BatteryCh=double(BatteryCh); BatteryDc=double(BatteryDc); InstallPowerrating=double(InstallPowerrating); InstallCapacity=double(InstallCapacity);
InstallLocation=double(InstallLocation);BatteryIF=double(BatteryIF);

Obj1=double(Obj1); Obj2=double(Obj2); Obj3=double(Obj3);

end
