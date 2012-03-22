function corr = jeff_correl(base_var, state_var)
%% jeff_correl  Computes the temporal evolution of the correlation of ...

%% DART software - Copyright 2004 - 2011 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% <next few lines under version control, do not edit>
% $URL$
% $Id$
% $Revision$
% $Date$

% create correlations at every timestep
figure(9); clf; hold on;

num_times = size(state_var, 1);
for i = 1:num_times
   x = corrcoef(base_var(i, :), state_var(i, :));
   corr(i) = x(1, 2);
end 

plot(base_var(:, 1), state_var(:, 1));
title('Correlation through time')

