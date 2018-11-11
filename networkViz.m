function network2 = networkViz(testCase);
%% Network visualisation
%
OriginBusLoc = testCase.bus(:,1); % number of buses 
OriginLinesFrom = testCase.branch(:,1);
OriginLinesTo = testCase.branch(:,2);
[tf, linesFrom]=ismember(OriginLinesFrom,OriginBusLoc,'rows');
[tf, linesTo]=ismember(OriginLinesTo,OriginBusLoc,'rows');
%
R = testCase.branch(:,3); % resistance
j = sqrt(-1);
%delta_R = 1e-4; 
%R(R < delta_R) = delta_R; % adding small resistance to every tansformer with zero resistance
X = testCase.branch(:,4); % reactance
Z = R + 1j * X;
%
network2 = digraph(linesFrom,linesTo, real(Z));
network2.Edges; % shows the number of Edges and Nodes 
%A = adjacency(network); % shows all nodes that llnes are connected to 
figure(2) 
postNetwork.Edges.LWidths = 3*abs(network2.Edges.Weight)/max(abs(network2.Edges.Weight))
if exist('testCase.bus_name','var') == 1
    netgraph = plot(network2,'Layout','force','Linewidth',2,'NodeLabel',BusName,'ArrowSize',12,'EdgeLabel',round(network2.Edges.Weight,2));    
else
    netgraph = plot(network2,'Layout','force','Linewidth',2,'ArrowSize',12,'EdgeLabel',round(network2.Edges.Weight,2));        
end
highlight(netgraph,find(testCase.bus(:,2)==1),'NodeColor','black', 'MarkerSize',6); % loads 
highlight(netgraph,find(testCase.bus(:,2)==3),'NodeColor','blue', 'MarkerSize',10); % gen 
highlight(netgraph,find(testCase.bus(:,2)==2),'NodeColor','red', 'MarkerSize',6); % gen 
set(gcf,'color','w')
end