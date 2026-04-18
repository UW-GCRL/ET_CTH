clc;
disp("start!");
new_file = 'H:\CLM_input\mean.nc';

nc_canopy_top_old = ncread('H:\CLM_input\surfdata_1.9x2.5_2005.nc','MONTHLY_HEIGHT_TOP');
pft = ncread('H:\CLM_input\surfdata_1.9x2.5_2005.nc','PCT_NAT_PFT');
old_data = nc_canopy_top_old(:,:,:,1);
nc_canopy_top_new_mean = ncread('H:\CLM_input\mean.nc','MONTHLY_HEIGHT_TOP');
mean_data = nc_canopy_top_new_mean(:,:,:,1);
nc_canopy_top_new_max = ncread('H:\CLM_input\max.nc','MONTHLY_HEIGHT_TOP');
max_data = nc_canopy_top_new_max(:,:,:,1);
nc_canopy_top_new_median = ncread('H:\CLM_input\median.nc','MONTHLY_HEIGHT_TOP');
median_data = nc_canopy_top_new_median(:,:,:,1);



for i=2:9
    % diffmean(:,:,i-1) = mean_data(:,:,i) - old_data(:,:,i);
    % diffmedian(:,:,i-1) = median_data(:,:,i) - old_data(:,:,i);
    % diffmax(:,:,i-1) = max_data(:,:,i) - old_data(:,:,i);
    diffmean(:,:,i-1) = mean_data(:,:,i);
    diffmedian(:,:,i-1) = median_data(:,:,i);
    diffmax(:,:,i-1) = max_data(:,:,i);

    a(i) = nanmean(nanmean(pft(:,:,i-1).*diffmean(:,:,i-1)));
    b(i) = nanmean(nanmean(pft(:,:,i-1).*diffmedian(:,:,i-1)));
    c(i) = nanmean(nanmean(pft(:,:,i-1).*diffmax(:,:,i-1)));
end

valmax(1:8) = 0;
count(1:8) = 0;
for i=1:8
    for m=1:144
        for n=1:96
            if(diffmax(m,n,i)~=0)
                count(i) = count(i)+1;
                valmax(i) = valmax(i)+diffmax(m,n,i);
            end
        end
    end
end
averagemax = valmax./count;

valmean(1:8) = 0;
count(1:8) = 0;
for i=1:8
    for m=1:144
        for n=1:96
            if(diffmean(m,n,i)~=0)
                count(i) = count(i)+1;
                valmean(i) = valmean(i)+diffmean(m,n,i);
            end
        end
    end
end
averagemean = valmean./count;

valmedian(1:8) = 0;
count(1:8) = 0;
for i=1:8
    for m=1:144
        for n=1:96
            if(diffmedian(m,n,i)~=0)
                count(i) = count(i)+1;
                valmedian(i) = valmedian(i)+diffmedian(m,n,i);
            end
        end
    end
end
averagemedian = valmedian./count;

for i=1:8
    %figure();subplot(1,3,1);imshow(diffmax(:,:,i),[-40,40]);subplot(1,3,2);imshow(diffmedian(:,:,i),[-40,40]);subplot(1,3,3);imshow(diffmean(:,:,i),[-40,40]);
    MEANofDIFFmax(i) = nanmean(nanmean(diffmax(:,:,i)));
    MEANofDIFFmean(i) = nanmean(nanmean(diffmean(:,:,i)));
    MEANofDIFFmedian(i) = nanmean(nanmean(diffmedian(:,:,i)));
end


figure('DefaultAxesFontSize',30);
subplot(3,1,1);bar(1:8,MEANofDIFFmax);title("Mean of MAX");xlabel('PFT');ylabel('Average height (m)');
subplot(3,1,2);bar(1:8,MEANofDIFFmedian);title("Mean of MEDIAN");xlabel('PFT');ylabel('Average height (m)');
subplot(3,1,3);bar(1:8,MEANofDIFFmean);title("Mean of MEAN");xlabel('PFT');ylabel('Average height (m)');

figure('DefaultAxesFontSize',30);
subplot(3,1,1);bar(1:8,averagemax);title("Mean of MAX vs Mean in 2020");xlabel('PFT');ylabel('Difference (m)');
subplot(3,1,2);bar(1:8,averagemedian);title("Mean of MEDIAN vs Mean in 2020");xlabel('PFT');ylabel('Difference (m)');
subplot(3,1,3);bar(1:8,averagemean);title("Mean of MEAN vs Mean in 2020");xlabel('PFT');ylabel('Difference (m)');
% 
% MAXpercentagepft1 = load("H:\CTH_with_pft_max\MAX_percentage_pft_1.txt");
% MAXpercentagepft2 = load("H:\CTH_with_pft_max\MAX_percentage_pft_2.txt");
% MAXpercentagepft3 = load("H:\CTH_with_pft_max\MAX_percentage_pft_3.txt");
% MAXpercentagepft4 = load("H:\CTH_with_pft_max\MAX_percentage_pft_4.txt");
% MAXpercentagepft5 = load("H:\CTH_with_pft_max\MAX_percentage_pft_5.txt");
% MAXpercentagepft6 = load("H:\CTH_with_pft_max\MAX_percentage_pft_6.txt");
% MAXpercentagepft7 = load("H:\CTH_with_pft_max\MAX_percentage_pft_7.txt");
% MAXpercentagepft8 = load("H:\CTH_with_pft_max\MAX_percentage_pft_8.txt");
% 
% pftpercent(1) = sum(sum(MAXpercentagepft1));
% pftpercent(2) = sum(sum(MAXpercentagepft2));
% pftpercent(3) = sum(sum(MAXpercentagepft3));
% pftpercent(4) = sum(sum(MAXpercentagepft4));
% pftpercent(5) = sum(sum(MAXpercentagepft5));
% pftpercent(6) = sum(sum(MAXpercentagepft6));
% pftpercent(7) = sum(sum(MAXpercentagepft7));
% pftpercent(8) = sum(sum(MAXpercentagepft8));
% A = sum(pftpercent);
% pftpercent = pftpercent/A*100;
% 
% figure('DefaultAxesFontSize',30);bar(1:8,pftpercent);title("Tree percentage");xlabel('PFT');ylabel('Difference (m)');
