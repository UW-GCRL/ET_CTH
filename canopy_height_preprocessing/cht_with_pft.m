CTH_ori = imread("H:\ETH_GlobalCanopyHeight_10m_2020_version1\ETH_GlobalCanopyHeight_10m_2020_version1\Global Canopy Height\cth2020_GEE_0.tif");
CTH_ori = CTH_ori';
disp("CTH_ori read finished!");

PFT_ori = imread("H:\PFT\PFT_500m_2020.tif");
PFT_ori = PFT_ori';
disp("PFT_ori read finished!");

for k =1:8
    clear CTH;
    clear MAP1;
    clear MAP2;
    clear MAPcth;
    clear PFT_per;

CTH = single(single(CTH_ori).*single((PFT_ori==k)));

MAP1 = zeros(144,96);
MAP2 = MAP1;
MAPcth = MAP1;

for i = 1:80150
    deltax = ceil((i/80150)*144);
    for j = 1:40076
        if((CTH(i,j)>0&&(CTH(i,j)<90)))
            deltay = ceil((j/40076)*96);
            MAP1(deltax,deltay) = MAP1(deltax,deltay)+1;
            MAP2(deltax,deltay) = MAP2(deltax,deltay)+CTH(i,j)*0.001;
        end
    end
    
    if(i==8000)
        disp("10%");
    elseif(i==8000*2)
        disp("20%");
    elseif(i==8000*3)
        disp("30%");
    elseif(i==8000*4)
        disp("40%");
    elseif(i==8000*5)
        disp("50%");
    elseif(i==8000*6)
        disp("60%");
    elseif(i==8000*7)
        disp("70%");
    elseif(i==8000*8)
        disp("80%");
    elseif(i==8000*9)
        disp("90%");
    elseif(i==8000*10)
        disp("100%");
    end
end

MAP1 = fliplr(MAP1);
MAP2 = fliplr(MAP2);
for i=1:144
    for j=1:96
        if(MAP1(i,j) > 0)
            MAPcth(i,j) = 1000.0 * MAP2(i,j) / MAP1(i,j);
        else
            MAPcth(i,j) = 0;
        end
    end
end
PFT_per = zeros(144,96);
num = 40076*80150/(144*96);

PFT_per = MAP1/num;

title = ["Figure of "+k];
figure('Name' , title);subplot(1,2,1);imshow(MAP1);subplot(1,2,2);imshow(PFT_per);
filepath1 = ["H:\cth_with_pft\cth_pft_"+k+".tiff"];
filepath2 = ["H:\cth_with_pft\count_pft_"+k+".tiff"];
filepath3 = ["H:\cth_with_pft\oriheight_pft_"+k+".tiff"];
filepath4 = ["H:\cth_with_pft\percentage_pft_"+k+".tiff"];
imwrite(MAPcth,filepath1);
imwrite(MAP1,filepath2);
imwrite(MAP2,filepath3);
imwrite(PFT_per,filepath3);

disp("PFT of "+k+" finished!");
end