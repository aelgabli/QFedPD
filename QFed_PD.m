function [logs, itr, number_of_bits_toSend, Energy]=QFed_PD(para,data,gra,fnc,init_x, bitsToSend, P, rate)
    logs=zeros(para.iter,2);
    itr=zeros(para.iter,1);
    X_0=init_x;
    X=X_0;
    lam = X_0-X;
    Q_X_0=zeros(para.node,size(X_0,2));
    prev_Q_X_0=zeros(para.node,size(X_0,2));
    QResolution=2*ones(para.node,1);
    Range=ones(para.node,1);
    flag=zeros(para.node,1);
    %number_of_bits_toSend=zeros(para.iter,para.node);
    number_of_bits_toSend=zeros(para.iter,1);
    Energy = zeros(para.iter,1);
    for i=1:para.iter
        i
        if (i) > 1
            number_of_bits_toSend(i)=number_of_bits_toSend(i-1);
            Energy(i) = Energy(i-1);
        end
        Dataset=transform(data,para.node,para.bs,1);
        logs(i,:)=ComputeLog(para.W*X_0,Dataset,gra,fnc,1);
        if (para.VR==1 && mod(i,para.I)==1)
            g_ex=getgrad(X,Dataset,gra,para.stepsize);
        end
        for j=1:para.localiter
            Dataset=transform(data,para.node,para.bs,para.mini_bs);
            if (para.VR==0)
                g_ex=getgrad(X,Dataset,gra,para.stepsize);
                if (para.mini_bs<1)
                    eta = 1.0/sqrt(para.localiter*(i-1)/5+j);
                else
                    eta=1;
                end
            else 
                eta=1;
            end
            X_old=X;
            X = X - eta*((X-X_0)+g_ex+para.stepsize*lam);
            if (para.VR==1)
                g_ex=g_ex+getgrad(X,Dataset,gra,para.stepsize)-getgrad(X_old,Dataset,gra,para.stepsize);
            end
        end
        lam = lam + (X-X_0)/para.stepsize;
        X_0=X_0 + lam*para.stepsize;
        if (para.R==-1)
            p = 1-i/para.iter;
            %p = i/para.iter;
        else
            p = 1-1.0/para.R;
        end
        if (rand()>p)
            
         for ii=1:para.node
            %[Q_X_0(ii,:),number_of_bits_toSend]=stochasticQ(Q_X_0(ii,:),X_0(ii,:),prev_Q_X_0(ii,:),bitsToSend);

            [Q_X_0(ii,:),number_of_bits_toSend(i), Range(ii), QResolution(ii)]=stochasticQ_DR(Q_X_0(ii,:),X_0(ii,:),prev_Q_X_0(ii,:),bitsToSend, Range(ii), QResolution(ii), flag(ii), number_of_bits_toSend(i));
            
            prev_Q_X_0(ii,:)=Q_X_0(ii,:);
            flag(ii) =1;
         end
            X_0=para.W*Q_X_0;
            itr(i) = 1;
            Energy(i) = Energy(i) + P*number_of_bits_toSend(i)/rate;
            %tot_number_of_bits_toSend(i)=sum(number_of_bits_toSend,2);
        end
    end