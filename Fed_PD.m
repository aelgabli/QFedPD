function [logs, itr, number_of_bits_toSend, Energy]=Fed_PD(para,data,gra,fnc,init_x, P, rate)
    logs=zeros(para.iter,2);
    itr=zeros(para.iter,1);
    X_0=init_x;
    X=X_0;
    lam = X_0-X;
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
        else
            p = 1-1.0/para.R;
        end
        if (rand()>p)
            X_0=para.W*X_0;
            itr(i) = 1;
            number_of_bits_toSend(i)=number_of_bits_toSend(i)+32*para.dimx*para.node;
            Energy(i) = Energy(i) + P*number_of_bits_toSend(i)/rate;
        end
    end