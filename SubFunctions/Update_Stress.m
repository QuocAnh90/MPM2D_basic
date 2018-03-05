function [F_sp,V_sp,s_sp,p_sp] = Update_Stress(dt,cellCount,mspoints,CONNECT,nvelo_si,dG,Lambda,Mu,F_sp,V_spo,m_sp,s_sp,p_sp,V_sp)

for c=1:cellCount
    mpts = mspoints{c};
    
    for sp = 1:length(mpts)
        spid = mpts(sp);
        L_sp = zeros(2,2);
    
        for j=1:4
            L_sp = L_sp + (nvelo_si(CONNECT(spid,j),:)'*dG((spid-1)*4+j,:));
        end       
        F_sp{spid} = (eye(2,2)+L_sp*dt)*F_sp{spid};
        J = det(F_sp{spid});
        V_sp(spid)=V_spo(spid)*J;
        SSP = Lambda*log(J)/J*eye(2,2) + Mu/J*(F_sp{spid}*F_sp{spid}'-eye(2,2));
        s_sp(spid,1) = SSP(1,1);
        s_sp(spid,2) = SSP(2,2);
        s_sp(spid,3) = SSP(1,2);
        
        p_sp(spid) = m_sp(spid)/V_sp(spid);
    end      
end