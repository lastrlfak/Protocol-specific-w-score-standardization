%% Protocol-specific w-score standardization

% Code from Jinyong Chung (06/06/2017)
% It is tested using MATLAB R2016b.
% You can contact me through e-mail (lastrlfak@kaist.ac.kr).

% %% Input data after freesurfer preprocessing 
% 
% AMC_NC_ct_ws: Cortical thickness data with Protocol 1 (third dataset), Matrix Dimension = Number of vertices X Number of subjects.
% SMC_NC_ct_ws: Cortical thickness data with Protocol 2 (third dataset), Matrix Dimension = Number of vertices X Number of subjects.
% 
% AMC_NC_ct: Cortical thickness data of NC with Protocol 1 (first dataset), Matrix Dimension = Number of vertices X Number of subjects.
% SMC_NC_ct: Cortical thickness data of NC with Protocol 2 (first dataset), Matrix Dimension = Number of vertices X Number of subjects.
% AMC_AD_ct: Cortical thickness data of AD with Protocol 1 (first dataset), Matrix Dimension = Number of vertices X Number of subjects.
% SMC_AD_ct: Cortical thickness data of AD with Protocol 2 (first dataset), Matrix Dimension = Number of vertices X Number of subjects.
% 
% "~_age": Age information (year), D = Number of subjects X 1.
% "~_ICV": Intracranial volume information (10^5 mm3), D = Number of subjects X 1.

%% Vertex-wise multiple linear regression using third dataset

% Outlier exclusion
mdl = fitlm([AMC_NC_age_ws AMC_NC_ICV_ws],mean(AMC_NC_ct_ws));
OLs_AMC_ws = find((mdl.Diagnostics.CooksDistance)<=3*mean(mdl.Diagnostics.CooksDistance));

mdl = fitlm([SMC_NC_age_ws SMC_NC_ICV_ws],mean(SMC_NC_ct_ws));
OLs_SMC_ws = find((mdl.Diagnostics.CooksDistance)<=3*mean(mdl.Diagnostics.CooksDistance));

% Vertex-wise multiple linear regression using third dataset (each protocol)

X = [ones(size(AMC_NC_age_ws(OLs_AMC_ws),1),1), AMC_NC_age_ws(OLs_AMC_ws), AMC_NC_ICV_ws(OLs_AMC_ws)];
X\AMC_NC_ct_ws(:,OLs_AMC_ws)';

%Beta weights
AMC_NC_ws_b0 = ans(1,:);
AMC_NC_ws_b1 = ans(2,:);
AMC_NC_ws_b2 = ans(3,:);

X = [ones(size(SMC_NC_age_ws(OLs_SMC_ws),1),1), SMC_NC_age_ws(OLs_SMC_ws), SMC_NC_ICV_ws(OLs_SMC_ws)];
X\SMC_NC_ct_ws(:,OLs_SMC_ws)';

%Beta weights
SMC_NC_ws_b0 = ans(1,:);
SMC_NC_ws_b1 = ans(2,:);
SMC_NC_ws_b2 = ans(3,:);

% Residual calculation for w-score standardization (each protocol)
SMC_NC_Residual_ws = (SMC_NC_ct_ws(:,OLs_SMC_ws)'-(SMC_NC_ws_b1.*SMC_NC_age_ws(OLs_SMC_ws)+SMC_NC_ws_b2.*SMC_NC_ICV_ws(OLs_SMC_ws)+repmat(SMC_NC_ws_b0,[size(OLs_SMC_ws,1) 1])));
AMC_NC_Residual_ws = (AMC_NC_ct_ws(:,OLs_AMC_ws)'-(AMC_NC_ws_b1.*AMC_NC_age_ws(OLs_AMC_ws)+AMC_NC_ws_b2.*AMC_NC_ICV_ws(OLs_AMC_ws)+repmat(AMC_NC_ws_b0,[size(OLs_AMC_ws,1) 1])));

%% Protocol-specific w-score standardization applying to first dataset

% Protocol-specific w-score standardization (each protocol)
AMC_NC_Wscore = (AMC_NC_ct' - (AMC_NC_ws_b1.*AMC_NC_age+AMC_NC_ws_b2.*AMC_NC_ICV+AMC_NC_ws_b0))./std(AMC_NC_Residual_ws);
AMC_AD_Wscore = (AMC_AD_ct' - (AMC_NC_ws_b1.*AMC_AD_age+AMC_NC_ws_b2.*AMC_AD_ICV+AMC_NC_ws_b0))./std(AMC_NC_Residual_ws);

SMC_NC_Wscore = (SMC_NC_ct' - (SMC_NC_ws_b1.*SMC_NC_age+SMC_NC_ws_b2.*SMC_NC_ICV+SMC_NC_ws_b0))./std(SMC_NC_Residual_ws);
SMC_AD_Wscore = (SMC_AD_ct' - (SMC_NC_ws_b1.*SMC_AD_age+SMC_NC_ws_b2.*SMC_AD_ICV+SMC_NC_ws_b0))./std(SMC_NC_Residual_ws);
