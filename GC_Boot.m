function [in_sample,out_sample,index,N_oob]=GC_Boot(Data,mode)
%% Gyeongcheol Cho's Bootstrapping function
    % input : 
    %   - Data: a (N x 1) vector or (N x J) matrix
    %           for the latter, bootstrapping will be implemented not for each column of A.
    %           Each row is considered a set of data and they are bootsampled for the entire data.
    %   - mode =0(default) if out_samples are not necessary
    %          =1          if out_samples are necessary
    % Output :
    %   - in_sample: Bootstrapping sample.
    %                (N by 1) if Data is a vector
    %                (N by J) if Data is a matrix.
    %   - out_sample: Sample not including in_sample
    %   - index : id for in_sample
    if nargin == 1
        mode=0;
    end
    N_oob=[]; out_sample=[];
    N=size(Data,1); %%% after 1901
    index=ceil(N*rand(N,1)); % Bootstrapping in terms of index
    in_sample=Data(index,:); 
    if mode ==1
        index_oob=(1:N)'; index_oob(index)=[];
        out_sample=Data(index_oob,:);
        N_oob=length(index_oob);
    end
end