%% 读取线路信息
clear;clc;
[Busdata,Gendata,branchdata,Gencostdata]=Data();
% 线路信息
LineI=branchdata(:,1);
LineJ=branchdata(:,2);
LineR=branchdata(:,3);
LineX=branchdata(:,4);
LineB=branchdata(:,5);
LineNum=length(LineI);
% 节点信息
NodeI=Busdata(:,1);
NodeNum=length(NodeI);
% 发电机信息
GenI=Gendata(:,1);

%% 初始化
IniPostiveBranchFlow=sparse([LineI;LineJ],[LineJ,LineI],[ones(LineNum,1);-ones(LineNum,1)],NodeNum,NodeNum);
IniReactiveBranchFlow=sparse([LineI;LineJ],[LineJ,LineI],[zeros(LineNum,1);zeros(LineNum,1)],NodeNum,NodeNum);
NodeIniVoltage=ones(NodeNum,1);
NodeIniTheta=zeros(NodeNum,1);
DNodeIniVoltage=zeros(NodeNum,1);


%% 首次带入
[PostiveBranchFlow,ReactiveBranchFlow,NodeVoltageS,DNodeVoltageS,NodeTheta,BatteryCh,BatteryDc,DBatteryCh,DBatteryDc,VIMultiDVI]=...
          FormulationDeclear(IniPostiveBranchFlow,IniReactiveBranchFlow,NodeIniVoltage,NodeIniTheta,DNodeIniVoltage);
        
DNodeIniVoltage=sqrt(DNodeVoltageS+NodeVoltageS+2*VIMultiDVI)-sqrt(NodeVoltageS);

%% 带入第二次解
IniPostiveBranchFlow=PostiveBranchFlow;
IniReactiveBranchFlow=ReactiveBranchFlow;
NodeIniVoltage=sqrt(NodeVoltageS);
NodeIniTheta=NodeTheta;


[PostiveBranchFlow,ReactiveBranchFlow,NodeVoltageS,DNodeVoltageS,NodeTheta,BatteryCh,BatteryDc,DBatteryCh,DBatteryDc,VIMultiDVI]=...
          FormulationDeclear(IniPostiveBranchFlow,IniReactiveBranchFlow,NodeIniVoltage,NodeIniTheta,DNodeIniVoltage);



%% 计算偏差
ACN=[branchdata(:,1) branchdata(:,2)]; %读取线路节点信息
for iteration=1:100
for i=1:length(ACN(:,1))
    PQ(i)=PQS(PostiveBranchFlow(ACN(i,1),ACN(i,2)),ReactiveBranchFlow(ACN(i,1),ACN(i,2)));
end
    MAXPQ=max(PQ);

for i=1:length(ACN(:,1))
    DeltaPQ(i,iteration)=abs(PQ(i)-PQS(IniPostiveBranchFlow(ACN(i,1),ACN(i,2)),IniReactiveBranchFlow(ACN(i,1),ACN(i,2))))/MAXPQ;
end   

% 判断收敛条件
   if DeltaPQ(:,iteration)<=0.01
      fprintf('迭代结束，收敛完毕\n'); break
   else
   IniPostiveBranchFlow=PostiveBranchFlow;   %更新初始迭代点
   IniReactiveBranchFlow=ReactiveBranchFlow; %更新初始迭代点
% NodeIniVoltage=sqrt(NodeVoltageS);
% NodeIniTheta=NodeTheta;
   DNodeIniVoltage=sqrt(DNodeVoltageS+NodeVoltageS+2*VIMultiDVI)-sqrt(NodeVoltageS);


   [PostiveBranchFlow,ReactiveBranchFlow,NodeVoltageS,DNodeVoltageS,NodeTheta,BatteryCh,BatteryDc,DBatteryCh,DBatteryDc,VIMultiDVI]=...
          FormulationDeclear(IniPostiveBranchFlow,IniReactiveBranchFlow,NodeIniVoltage,NodeIniTheta,DNodeIniVoltage);
   end
end
IniPostiveBranchFlow=PostiveBranchFlow;
IniReactiveBranchFlow=ReactiveBranchFlow;
% NodeIniVoltage=sqrt(NodeVoltageS);
% NodeIniTheta=NodeTheta;



DNodeIniVoltage=sqrt(DNodeVoltageS+NodeVoltageS+2*VIMultiDVI)-sqrt(NodeVoltageS);


[PostiveBranchFlow,ReactiveBranchFlow,NodeVoltageS,DNodeVoltageS,NodeTheta,BatteryCh,BatteryDc,DBatteryCh,DBatteryDc,VIMultiDVI]=...
          FormulationDeclear(IniPostiveBranchFlow,IniReactiveBranchFlow,NodeIniVoltage,NodeIniTheta,DNodeIniVoltage);
