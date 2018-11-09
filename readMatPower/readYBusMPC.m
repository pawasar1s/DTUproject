function [ZBus, YBus, YBusSlack, Ymn] = readYBusMPC(linesFrom, linesTo, Z, Ybus, nLines, nBuses)
%% Bus Admittance Matrix
%ZBusSlack = zeros(nBuses,nBuses); % initiating the impedance matrix Zbus
YBusSlack = zeros(nBuses,nBuses); % initiating the admittance matrix Ybus
%Formation of the off diagonal elements
for m = 1 : nLines
    %ZBusSlack(linesFrom(m),linesTo(m)) = ZBusSlack(linesFrom(m),linesTo(m)) + Z(m); % setting values for the right side
    %ZBusSlack(linesTo(m),linesFrom(m)) = ZBusSlack(linesFrom(m),linesTo(m)); % mirroring the right side
    YBusSlack(linesFrom(m),linesTo(m)) = YBusSlack(linesFrom(m),linesTo(m)) + Ybus(m); % setting values for the right side
    YBusSlack(linesTo(m),linesFrom(m)) = YBusSlack(linesFrom(m),linesTo(m)); % mirroring the right side
end
% Formation of the diagonal elements
 for n = 1:nBuses
     for m = 1:nLines
         if linesFrom(m) == n || linesTo(m) == n
             %ZBusSlack(n,n) = ZBusSlack(n,n) + Z(m);
             YBusSlack(n,n) = YBusSlack(n,n) + Ybus(m);
         end
     end
 end
%ZBusSlack;
%
%rr = real(ZBusSlack);
%xx = imag(ZBusSlack);
Ymn = 2*[real(conj(YBusSlack)) zeros(nBuses,nBuses); zeros(nBuses,nBuses) real(conj(YBusSlack))];
%ZBus = ZBusSlack(2:end,2:end);
YBus = YBusSlack(2:end,2:end);
ZBus = inv(YBus);
end