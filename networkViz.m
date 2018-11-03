function network2 = networkViz(mpc,linesFrom, linesTo, Z);
%% Network visualisation

network2 = digraph(linesFrom,linesTo, real(Z));  
network2.Edges; % shows the number of Edges and Nodes 
%A = adjacency(network); % shows all nodes that llnes are connected to 
figure(2) 
postNetwork.Edges.LWidths = 3*abs(network2.Edges.Weight)/max(abs(network2.Edges.Weight))
if exist('mpc.bus_name','var') == 1
    netgraph = plot(network2,'Layout','force','Linewidth',2,'NodeLabel',BusName,'ArrowSize',12,'EdgeLabel',round(network2.Edges.Weight,2));    
else
    netgraph = plot(network2,'Layout','force','Linewidth',2,'ArrowSize',12,'EdgeLabel',round(network2.Edges.Weight,2));        
end
highlight(netgraph,find(mpc.bus(:,2)==1),'NodeColor','black', 'MarkerSize',6); % loads 
highlight(netgraph,find(mpc.bus(:,2)==3),'NodeColor','blue', 'MarkerSize',10); % gen 
highlight(netgraph,find(mpc.bus(:,2)==2),'NodeColor','red', 'MarkerSize',6); % gen 
set(gcf,'color','w')
end