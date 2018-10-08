
% fully regular sampled data case:
 reg1longr = load ('sampled_lon_1.txt');
 reg1latgr = load ('sampled_lat_1.txt');
 reg2datgr = load ('sampled_data_1.txt');

figure(1);
plotit_reg(reg1longr, reg1latgr, reg2datgr);

% fully regular sampled data case:
 reg1longr = load ('sampled_lon_2.txt');
 reg1latgr = load ('sampled_lat_2.txt');
 reg2datgr = load ('sampled_data_2.txt');

figure(2);
plotit_reg(reg1longr, reg1latgr, reg2datgr);

% rotated sampled data case:
 reg1longr = load ('sampled_lon_3.txt');
 reg1latgr = load ('sampled_lat_3.txt');
 reg2datgr = load ('sampled_data_3.txt');

figure(3);
plotit_reg(reg1longr, reg1latgr, reg2datgr);

% rotated, deformed sampled data case:
 reg1longr = load ('sampled_lon_4.txt');
 reg1latgr = load ('sampled_lat_4.txt');
 reg2datgr = load ('sampled_data_4.txt');

figure(4);
plotit_reg(reg1longr, reg1latgr, reg2datgr);

%@% % regular case:
%@% 
%@% % data grid - fully regular, values output only for plotting
%@%  reg1longr = load ('data_lons_1d_reg_test.txt');
%@%  reg1latgr = load ('data_lats_1d_reg_test.txt');
%@%  reg2datgr = load ('data_data_2d_reg_test.txt');
%@% 
%@% % sampled data - regular grid rotated by some angle
%@% %                so considered fully irregular
%@%  reg2lonsp = load ('sample_lons_2d_reg_test.txt');
%@%  reg2latsp = load ('sample_lats_2d_reg_test.txt');
%@%  reg2datsp = load ('sample_data_2d_reg_test.txt');
%@% 
%@% 
%@% % irregular case:
%@% 
%@% % data grid - fully irregular quads
%@%  irr2longr = load ('data_lons_2d_irreg_test.txt');
%@%  irr2latgr = load ('data_lats_2d_irreg_test.txt');
%@%  irr2datgr = load ('data_data_2d_irreg_test.txt');
%@% 
%@% % sampled data - fully regular, values output only for plotting
%@%  irr1lonsp = load ('sample_lons_1d_irreg_test.txt');
%@%  irr1latsp = load ('sample_lats_1d_irreg_test.txt');
%@%  irr2datsp = load ('sample_data_2d_irreg_test.txt');
%@% 
%@% 
%@% figure(1);
%@% plotit_reg(reg1longr, reg1latgr, reg2datgr);
%@% figure(2);
%@% plotit_irr(reg2lonsp, reg2latsp, reg2datsp, 1);
%@% 
%@% figure(3);
%@% plotit_irr(irr2longr, irr2latgr, irr2datgr, 40);
%@% figure(4);
%@% plotit_reg(irr1lonsp, irr1latsp, irr2datsp);
