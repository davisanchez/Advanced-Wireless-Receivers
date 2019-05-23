function h = channel(RX,NumberOfSymbols,NumberOfCarriers,CoherenceTime,CoherenceBandwidth)
tfactor = 2;
tborder = max(100,2*CoherenceTime*tfactor);
fborder = max(100,2*CoherenceBandwidth);

grid = 1/sqrt(2) * (randn(NumberOfCarriers+2*fborder,NumberOfSymbols+2*tborder,RX) + 1i * randn(NumberOfCarriers+2*fborder,NumberOfSymbols+2*tborder,RX));

% 1D filtering

wnt = 2/(CoherenceTime*tfactor);

if wnt == 1
    wnt = wnt - 0.00001;
end

fct = 1 / CoherenceTime;
lowpass = fir1(2*tborder+1,wnt);

for i=1:RX
    lraw = conv(grid(fborder+1,:,i),lowpass,'same');
    vraw = 1/length(lraw)*sum(abs(lraw).^2);
    grid(fborder+1,:,i) = 1/sqrt(vraw)*lraw;
end

h = grid(fborder+1:end-fborder,tborder+1:end-tborder,:);