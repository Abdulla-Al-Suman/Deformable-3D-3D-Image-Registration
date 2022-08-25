function DSC_value = DSC(A_A,A_B)

A_AiA_B = and(A_A,A_B);A_A = double(A_A);A_B = double(A_B);A_AiA_B = double(A_AiA_B);
NA_A= sum(A_A(:)); NA_B= sum(A_B(:));NA_AiA_B=sum(A_AiA_B(:));
DSC_value = ((2*NA_AiA_B)/(NA_A+NA_B));
clear A_AiA_B NA_A NA_B NA_AiA_B;   

