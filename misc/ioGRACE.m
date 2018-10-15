%% IO of Rory's GRACE data

% Read in dat
fid1=fopen("coding/GSFC.glb.200301_201607_v02.3b_1d_trend_2005_2015.dat");
fid2=fopen("coding/GSFC.glb.200301_201607_v02.3b-ICE6G_1d_trend_2005_2015.dat");
trnd1=fread(fid1,'float32');
trnd2=fread(fid2,'float32');
fclose(fid1);
fclose(fid2);


% Output as txt
dlmwrite("coding/grace_grid.txt", trnd1(2:64801))
dlmwrite("coding/grace-gia_grid.txt", trnd2(2:64801))