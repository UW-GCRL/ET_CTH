ori_cht = imread('H:/ETH_GlobalCanopyHeight_10m_2020_version1/ETH_GlobalCanopyHeight_10m_2020_version1/Global Canopy Height/CanopyTopHeight2020_500m.tif');
ori_pft = imread('H:/CLM input/PFT_500m_2020.tif');

i=1;
    PFT = ori_pft==i;
    disp("PFT "+i+" mapping finished.");

    for a = 1:width(ori_cht)
        if (a == floor(length(ori_cht)/10))
            disp("10");
        end
        if (a == floor(length(ori_cht)/5))
            disp("20");
        end
        if (a == floor(length(ori_cht)/2))
            disp("50");
        end
        if (a == 8*floor(length(ori_cht)/10))
            disp("80");
        end
        if (a == 9*floor(length(ori_cht)/10))
            disp("90");
        end
        for b = 1:length(ori_cht)
            cth(a,b) = int8(int8(PFT(a,b))*ori_cht(a,b));
        end
    end
    disp("CTH wih PFT "+i+" mapping finished.");

    char = "H:\cthwithpft\CTH_with_PFT"+i+".png";
    imwrite(cth,char);
    clear PFT;
    clear cth;
    disp("CTH wih PFT "+i+" output finished.");

