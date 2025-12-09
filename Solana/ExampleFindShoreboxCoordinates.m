% Example to find shorebox coordiates of a point

% The Solana Shorebox coordinate frame is based on Mops 635 to Mop 665

% A guess at the camera location is 32.9907  -117.2742

% use deg2ShoreBox.m to convert lat lon to shorebox x,y
[Xsb,Ysb]=deg2ShoreBox(32.9907,-117.2742,635,665);
fprintf('Shorebox x,y coordinates: %10.2f %10.2f\n',Xsb,Ysb)

% There is also a utm2ShoreBox.m to convert Eutm Nutm to shorebox x,y