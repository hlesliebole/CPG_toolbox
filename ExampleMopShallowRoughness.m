RS=GetMopShoalsRoughnessV2(666,11,5);
idx=find(RS(2).Z < -1);
figure;ColorScatter(RS(2).Xutm(idx),RS(2).Yutm(idx),RS(2).Sigma(idx));view(2)
figure;ColorScatter(RS(2).Xutm(idx),RS(2).Yutm(idx),RS(2).Z(idx));view(0,0)