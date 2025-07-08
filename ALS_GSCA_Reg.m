function [W,C,B,vec_err]=ALS_GSCA_Reg(Z,W,W0,C0,B0,ind_Adep,lambda_b,Min_limit,Max_iter,N,J,P,T,Jy,Py,loc_Cdep,loc_Bdep)
% to set initial value, and check the characteristics of Data
    iter=0;
    improve=100000000;
    dif_c=100000000;   
    eye_J=eye(J);
% Normalize & Sample Covariance    
    Z=zscore(Z,1)/sqrt(N);    
    cov_Z=Z'*Z;
    for p=1:P
        w_p=W(W0(:,p),p);
        cov_zp=cov_Z(W0(:,p),W0(:,p));
        norm_cv=w_p'*cov_zp*w_p;
        W(W0(:,p),p)=W(W0(:,p),p)/sqrt(norm_cv);
    end
% Reduce the computational cost
    if N>J
        [R,flag]=chol(cov_Z);
        if flag==0; Z=R; end
    end
    CV = Z*W;

    C=double(C0); B=double(B0);
    while improve > Min_limit && iter < Max_iter
        iter=iter+1;
        V=[eye_J, W]; 
    % (2-1) to estimate A given W ... SS(ZV-ZWA)=SS(Psi-GamA)
        %{ 
            oneshot
            % Psi = Z * V;
            % CV = Z * W; % N by P 
            % Phi=kron(eye(T),CV);
            % Phi=Phi(:,ind_A);
            % a=(Phi'*Phi)\Phi'*Psi(:);
            % A(ind_A)=a;
        %}
        if Jy>0
            for j=loc_Cdep
                CVx = CV(:,C0(:,j));
                C(C0(:,j),j)=(CVx'*CVx)\(CVx'*Z(:,j));
            end
        end
        if Py>0
            penalty_b=0;
            for q=loc_Bdep
                CVx = CV(:,B0(:,q));             
                bq_est=(CVx'*CVx+lambda_b*eye(sum(B0(:,q),1)))\(CVx'*CV(:,q));
                B(B0(:,q),q)=bq_est;
                penalty_b = penalty_b + bq_est'*bq_est;
            end
        end
        A=[C,B];
    % (2-2) to estimate W given A    
        for p=1:P
            dif_p=dif_c;
            t=J+p;

            V_m_t=V; V_m_t(:,t)=0; % _m_ : minus
            A_m_p=A; A_m_p(p,:)=0; % _p_ : existent
            a_p=A(p,:);

            W_m_p=W; W_m_p(:,p)=0;

            Lam_m_p=W_m_p*A_m_p;

            i_t=zeros(1,T);
            i_t(1,t)=1;

            beta=i_t - a_p;
            delta=Lam_m_p - V_m_t;

            Ksi=kron(beta',Z);
            Ksi=Ksi(:,W0(:,p));
            vec_y=Z*delta;
            vec_y=vec_y(:);
            theta=((Ksi'*Ksi)\Ksi')*vec_y;

            norm_cv=theta'*cov_Z(W0(:,p),W0(:,p))*theta;
            theta=theta/sqrt(norm_cv);
            W(W0(:,p),p)=theta;
        end
        CV=Z*W;
        ERROR=[Z CV] - CV*A;
        vec_err=sum(ERROR(:,ind_Adep).^2,1);
        dif_c=sum(vec_err,2)+lambda_b*penalty_b;   
        improve=dif_p-dif_c; % sum(error_past^2) - sum(error_current^2) 
     end  
     C=A(:,1:J);
     B=A(:,(J+1):T);
     if Jy>0
         for j=loc_Cdep
             CVx = CV(:,C0(:,j));
             C(C0(:,j),j)=(CVx'*CVx)\(CVx'*Z(:,j));
         end
     end
     if Py>0
         for q=loc_Bdep
             CVx = CV(:,B0(:,q));                
             B(B0(:,q),q)=(CVx'*CVx+lambda_b*eye(sum(B0(:,q),1)))\(CVx'*CV(:,q)); 
         end
     end
end