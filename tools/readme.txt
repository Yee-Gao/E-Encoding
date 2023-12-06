The function to do FDR correction can be downloaded from 
https://www.mathworks.com/matlabcentral/fileexchange/27418-fdr_bh

[h, crit_p, adj_p]=fdr_bh(pvals,q,method,report);



The colormap used in our figures can be downloaded from
https://www.mathworks.com/matlabcentral/fileexchange/30564-othercolor

colormap(flipud(othercolor('RdBu4')))



The function to compute modulation spectrum can be downloaded from
https://www.sciencedirect.com/science/article/pii/S0149763416305668?via%3Dihub

[ms,f]=Modulation_Spectrum(v,fs,'log');