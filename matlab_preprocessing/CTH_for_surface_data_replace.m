clear all;clc;
% standard = imread('F:\cth_with_pft\1.tif');
nc_canopy_top_old = ncread('H:\CLM_input\surfdata_1.9x2.5_2005.nc','MONTHLY_HEIGHT_TOP');
CTH_2020 = nc_canopy_top_old;
MAP = nc_canopy_top_old(:,:,:,1);
MAPNEW = MAP;
temp = zeros(144,96);

pos = zeros(144,96);
for i=1:8
    filepath1 = ["H:\CTH_with_pft_max\MAX_cth_pft_"+i+".tiff"];
    %filepath2 = ["H:\CTH_with_pft_median\MEDIAN_percentage_pft_"+i+".tiff"];
    
    MAX_cth_pft = imread(filepath1);
    %PFT = imread(filepath2);

    for j = 1:144
        line = j-71;
        if(line <= 0)
        line = line + 144;
        end
        % temp_PFT(line,:) = PFT(j,:);
    end

    for j = 1:144
        line = j-71;
        if(line <= 0)
        line = line + 144;
        end
        temp_CTH(line,:) = MAX_cth_pft(j,:);
    end
    for a = 1:144
        for b=1:96
            if(MAPNEW(a,b,i+1)>0 && temp_CTH(a,b)>0)
            MAPNEW(a,b,i+1) = temp_CTH(a,b);
            pos(a,b) = 1;
            end
        end
    end

    for m = 1:12
        CTH_2020(:,:,i+1,m) = MAPNEW(:,:,i+1);
    end
%    figure();subplot(1,2,2);imshow(MAPNEW(:,:,2),[]);subplot(1,2,1);imshow(MAP(:,:,2),[]);
end
imwrite(pos,"H:\CTH_with_pft_median\position.tiff");
new_file = 'H:\CLM_input\cth0.nc';
ncdisp(new_file);
ncid_new = netcdf.open(new_file,'WRITE'); % open nc file as 'read-write' mode
read_CFT_id = netcdf.inqVarID(ncid_new,'MONTHLY_HEIGHT_TOP'); % read var id
read_CFT_org = netcdf.getVar(ncid_new,read_CFT_id); % get variable
read_CFT_org = CTH_2020(:,:,:,:);
netcdf.putVar(ncid_new,read_CFT_id,read_CFT_org); % write new CFT
netcdf.close(ncid_new); % close nc file
disp('finish all for 2020');

pft = ncread('H:\CLM_input\surfdata_1.9x2.5_2005.nc','PCT_NAT_PFT');
for i=1:15
    meanpft(i) = nanmean(nanmean(pft(:,:,i)));
end
sumpft = sum(sum(sum(meanpft)))