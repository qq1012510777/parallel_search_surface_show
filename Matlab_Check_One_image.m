clc
clear all
close all
currentPath = fileparts(mfilename('fullpath'));

A1 = (h5read([currentPath, '/inp/cut1.h5'], '/data'));
A2 = (h5read([currentPath, '/bin/distance_of_pore_to_bone.h5'], '/data'));
A3 = (h5read([currentPath, '/bin/size_of_pore.h5'], '/data'));

A1([2:end], :, :) = [];
A2([2:end], :, :) = [];
A3([2:end], :, :) = [];

A1 = squeeze(A1);
A2 = squeeze(A2);
A3 = squeeze(A3);

[xx, yy] = meshgrid(1:size(A2, 1),1:size(A2, 2));

figure(1)
subplot(1, 3, 1)
pbaspect([1, 1, 1]); hold on
title('Digital sample'); hold on
contourf(xx, yy, A1)
colorbar

figure(1)
subplot(1, 3, 2)
pbaspect([1, 1, 1]); hold on
title('Distance between a pore to its closest coral bone'); hold on
contourf(xx, yy, A2)
colorbar

figure(1)
subplot(1, 3, 3)
pbaspect([1, 1, 1]); hold on
title('Visualization of pore sizes'); hold on
contourf(xx, yy, A3)
colorbar

AS = find(A2 <= 1);
JK = A3(AS);

a = find(JK > 1)