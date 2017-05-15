function copy_index = get_copy_index(fname, copystring)
%% GET_COPY_INDEX  Gets an index corresponding to copy meta_data string
% Retrieves index associated with a given copy meta_data string in
% file fname. If string is not found in meta_data list, a -1 is returned.
%
% Example:
% fname = 'obs_diag_output.nc';
% copystring = 'N_DARTqc_5';
% copy_index = get_copy_index(fname, copystring);

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

if ( exist(fname,'file') ~= 2 ), error('%s does not exist.',fname); end

% Matlab seems to always need to transpose character variables.
copy_meta_data  = ncread(fname,'CopyMetaData')';
[num_copies, ~] = nc_dim_info(fname,'copy');
[metalen, ~]    = nc_dim_info(fname,'stringlength');

if( size(copy_meta_data,1) ~= num_copies || size(copy_meta_data,2) ~= metalen)
    error('%s from %s does not have the shape expected',copystring,fname)
end

nowhitecs = dewhite(copystring);

% Figure out which copy is the matching one
copy_index = -1;
for i = 1:num_copies,

   % for matching -- we want to ignore whitespace -- find it & remove it
   nowhitemd = dewhite(copy_meta_data(i,:));

   if strcmp(nowhitemd , nowhitecs) == 1
      copy_index = i;
   end
end

% Provide modest error support

if (copy_index < 0)
   fprintf('ERROR: %s is not a valid metadata string for file %s\n', ...
                strtrim(copystring), fname)
   disp('valid metadata strings are: ')
   for i = 1:num_copies,
      fprintf('%s\n',deblank(copy_meta_data(i,:)))
   end
   error('Thats all. Start over')
end


function str2 = dewhite(str1)
%  function to remove ALL whitespace from a character string

str2 = str1(~isspace(str1));


% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
