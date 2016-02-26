function oplus = plot_interp_diffs(fname)
%%
% 
%    fname = '../data/dart.oplus2.bin';
%    plot_interp_diffs(fname)
%

%% DART software - Copyright 2004 - 2013 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

ggcm = fopen(fname,'r');
byte = fread(ggcm,1,'int32');
dim1 = fread(ggcm,1,'int32');
dim2 = fread(ggcm,1,'int32');
dim3 = fread(ggcm,1,'int32');
byte = fread(ggcm,1,'int32');

byte = fread(ggcm,1,'int32');
gdat = fread(ggcm,dim1*dim2*dim3,'float64');
byte = fread(ggcm,1,'int32');
stat = fclose(ggcm);

datmat = reshape(gdat,dim1,dim2,dim3);

f2 = '../data/DATA.ionos2.nc';
oplus = nc_varget(f2,'oplus');

size(datmat)
size(oplus)

for ilon = 1:dim3

   fprintf('longitude %d\n',ilon)

   figure(1); orient tall;
   subplot(2,1,1)
   level = 1;
   imagesc(squeeze(datmat(:,:,ilon)));
   
   h= title({sprintf('%s longitude %d',fname, ilon), sprintf('%d -x- %d -x- %d',dim1,dim2,dim3)});
   set(gca,'YDir','normal');
   ylabel(sprintf('dimension %d (length=%d)',1,dim1))
   xlabel(sprintf('dimension %d (length=%d)',2,dim2))
   axis image
   wysiwyg
   set(gca,'FontSize',20)
   colorbar('eastoutside')
   
   subplot(2,1,2)
   imagesc(squeeze(oplus(ilon,:,:))');
   set(gca,'Ydir','normal')
   axis image
   set(gca,'FontSize',20)
   colorbar('eastoutside')

   disp('Pausing, hit any key to continue ...')
   pause

end



% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$

