X=[]; Y=[]; FRS=[]; DENS = [];
load('/home/zach/Documents/Sunny Data/asdf_networks_2013_jan/asdf_data/1.mat')
load('/home/zach/Documents/Sunny Data/asdf_networks_2013_jan/XY locations Jan 2013/xy_02_0.mat')
DENS = [DENS findDense(x, y)]; X = [X x]; Y = [Y y]; FRS = [FRS; cellfun(@(x) length(x)./asdf_raw{end}(2), asdf_raw(1:end-2))];
load('/home/zach/Documents/Sunny Data/asdf_networks_2013_jan/asdf_data/2.mat')
load('/home/zach/Documents/Sunny Data/asdf_networks_2013_jan/XY locations Jan 2013/xy_03_0.mat')
DENS = [DENS findDense(x, y)];X = [X x]; Y = [Y y]; FRS = [FRS; cellfun(@(x) length(x)./asdf_raw{end}(2), asdf_raw(1:end-2))];
load('/home/zach/Documents/Sunny Data/asdf_networks_2013_jan/asdf_data/3.mat')
load('/home/zach/Documents/Sunny Data/asdf_networks_2013_jan/XY locations Jan 2013/xy_03_1.mat')
DENS = [DENS findDense(x, y)];X = [X x]; Y = [Y y]; FRS = [FRS; cellfun(@(x) length(x)./asdf_raw{end}(2), asdf_raw(1:end-2))];
load('/home/zach/Documents/Sunny Data/asdf_networks_2013_jan/asdf_data/4.mat')
load('/home/zach/Documents/Sunny Data/asdf_networks_2013_jan/XY locations Jan 2013/xy_04_0.mat')
DENS = [DENS findDense(x, y)];X = [X x]; Y = [Y y]; FRS = [FRS; cellfun(@(x) length(x)./asdf_raw{end}(2), asdf_raw(1:end-2))];
load('/home/zach/Documents/Sunny Data/asdf_networks_2013_jan/asdf_data/5.mat')
load('/home/zach/Documents/Sunny Data/asdf_networks_2013_jan/XY locations Jan 2013/xy_04_1.mat')
DENS = [DENS findDense(x, y)];X = [X x]; Y = [Y y]; FRS = [FRS; cellfun(@(x) length(x)./asdf_raw{end}(2), asdf_raw(1:end-2))];
load('/home/zach/Documents/Sunny Data/asdf_networks_2013_jan/asdf_data/6.mat')
load('/home/zach/Documents/Sunny Data/asdf_networks_2013_jan/XY locations Jan 2013/xy_05_0.mat')
DENS = [DENS findDense(x, y)];X = [X x]; Y = [Y y]; FRS = [FRS; cellfun(@(x) length(x)./asdf_raw{end}(2), asdf_raw(1:end-2))];
load('/home/zach/Documents/Sunny Data/asdf_networks_2013_jan/asdf_data/7.mat')
load('/home/zach/Documents/Sunny Data/asdf_networks_2013_jan/XY locations Jan 2013/xy_05_1.mat')
DENS = [DENS findDense(x, y)];X = [X x]; Y = [Y y]; FRS = [FRS; cellfun(@(x) length(x)./asdf_raw{end}(2), asdf_raw(1:end-2))];
load('/home/zach/Documents/Sunny Data/asdf_networks_2013_jan/asdf_data/8.mat')
load('/home/zach/Documents/Sunny Data/asdf_networks_2013_jan/XY locations Jan 2013/xy_06_0.mat')
DENS = [DENS findDense(x, y)];X = [X x]; Y = [Y y]; FRS = [FRS; cellfun(@(x) length(x)./asdf_raw{end}(2), asdf_raw(1:end-2))];
load('/home/zach/Documents/Sunny Data/asdf_networks_2013_jan/asdf_data/9.mat')
load('/home/zach/Documents/Sunny Data/asdf_networks_2013_jan/XY locations Jan 2013/xy_06_1.mat')
DENS = [DENS findDense(x, y)];X = [X x]; Y = [Y y]; FRS = [FRS; cellfun(@(x) length(x)./asdf_raw{end}(2), asdf_raw(1:end-2))];
load('/home/zach/Documents/Sunny Data/asdf_networks_2013_jan/asdf_data/10.mat')
load('/home/zach/Documents/Sunny Data/asdf_networks_2013_jan/XY locations Jan 2013/xy_07_0.mat')
DENS = [DENS findDense(x, y)];X = [X x]; Y = [Y y]; FRS = [FRS; cellfun(@(x) length(x)./asdf_raw{end}(2), asdf_raw(1:end-2))];
load('/home/zach/Documents/Sunny Data/asdf_networks_2013_jan/asdf_data/11.mat')
load('/home/zach/Documents/Sunny Data/asdf_networks_2013_jan/XY locations Jan 2013/xy_07_1.mat')
DENS = [DENS findDense(x, y)];X = [X x]; Y = [Y y]; FRS = [FRS; cellfun(@(x) length(x)./asdf_raw{end}(2), asdf_raw(1:end-2))];
load('/home/zach/Documents/Sunny Data/asdf_networks_2013_jan/asdf_data/12.mat')
load('/home/zach/Documents/Sunny Data/asdf_networks_2013_jan/XY locations Jan 2013/xy_08_0.mat')
DENS = [DENS findDense(x, y)];X = [X x]; Y = [Y y]; FRS = [FRS; cellfun(@(x) length(x)./asdf_raw{end}(2), asdf_raw(1:end-2))];
load('/home/zach/Documents/Sunny Data/asdf_networks_2013_jan/asdf_data/13.mat')
load('/home/zach/Documents/Sunny Data/asdf_networks_2013_jan/XY locations Jan 2013/xy_08_1.mat')
DENS = [DENS findDense(x, y)];X = [X x]; Y = [Y y]; FRS = [FRS; cellfun(@(x) length(x)./asdf_raw{end}(2), asdf_raw(1:end-2))];
load('/home/zach/Documents/Sunny Data/asdf_networks_2013_jan/asdf_data/14.mat')
load('/home/zach/Documents/Sunny Data/asdf_networks_2013_jan/XY locations Jan 2013/xy_09_0.mat')
DENS = [DENS findDense(x, y)];X = [X x]; Y = [Y y]; FRS = [FRS; cellfun(@(x) length(x)./asdf_raw{end}(2), asdf_raw(1:end-2))];
load('/home/zach/Documents/Sunny Data/asdf_networks_2013_jan/XY locations Jan 2013/xy_09_1.mat')
load('/home/zach/Documents/Sunny Data/asdf_networks_2013_jan/asdf_data/15.mat')
DENS = [DENS findDense(x, y)];X = [X x]; Y = [Y y]; FRS = [FRS; cellfun(@(x) length(x)./asdf_raw{end}(2), asdf_raw(1:end-2))];

[rho, p] = corr([DENS' FRS])
[FO G] = fit(DENS',FRS,  'poly1')
