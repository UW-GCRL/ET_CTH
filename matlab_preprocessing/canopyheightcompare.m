clc;
clear all;
close all;
CTH_max_2020 = imread("H:\123\final_map\Aggrega_Fina11.tif");
CTH_max_2005 = imread("H:\123\final_map\Forest_Canopy_HeightWGS84_Pr1.tif");
[wid,len] = size(CTH_max_2020);
NUM = wid/10;
MAP = single(nan(wid,len));
increase = single(nan(wid,len));
decrease = single(nan(wid,len));
for i=1:wid
    for j=1:len
        if(CTH_max_2020(i,j)>0 && CTH_max_2005(i,j)>0)
            MAP(i,j) = double(CTH_max_2005(i,j)) - double(CTH_max_2020(i,j));
            if(MAP(i,j)>5)
                decrease(i,j) = MAP(i,j);
            elseif(MAP(i,j)<5)
                increase(i,j) = MAP(i,j);
            end
        end
    end

    if(i==NUM)
        disp("10%");
    elseif(i==NUM*2)
        disp("20%");
    elseif(i==NUM*3)
        disp("30%");
    elseif(i==NUM*4)
        disp("40%");
    elseif(i==NUM*5)
        disp("50%");
    elseif(i==NUM*6)
        disp("60%");
    elseif(i==NUM*7)
        disp("70%");
    elseif(i==NUM*8)
        disp("80%");
    elseif(i==NUM*9)
        disp("90%");
    elseif(i==NUM*10)
        disp("100%");
    end

end
figure('DefaultAxesFontSize',25);imshow(-MAP,[-20,20]);colormap("turbo");colorbar;
figure('DefaultAxesFontSize',18);imshow(abs(decrease),[-20,20]);colormap("parula");colorbar;title(["Decrease magnitude.";"Starting from 0, the minimum colorbar of 20 in the figure is to distinguish the background color."]);
figure('DefaultAxesFontSize',18);imshow(abs(increase),[-20,20]);colormap("parula");colorbar;title(["Increase magnitude.";"Starting from 0, the minimum colorbar of 20 in the figure is to distinguish the background color."]);

mean = nanmean(nanmean(MAP));
mean_abs = nanmean(nanmean(abs(MAP)));
MAP_data = reshape(MAP,[1,wid*len]);
MAP_data = MAP_data(MAP_data>-1000);
[h,p,ci]=ttest(MAP_data);
