function [W,A] = ntcalc(NX,N,X,HK,TH,UE,VE,NW)
% C------------------------------------------------------------------
% C     Calculates range of frequencies which span the
% C     critical frequency.  Also calculates the amplitude
% C     distribution for each frequency.
% C
% C    Input:  NX     array physical dimension
% C            N      number of streamwise points i
% C                    (i = N  point is assumed turbulent)
% C            X (i)  streamwise coordinate array for integrating A(x)
% C            HK(i)  kinematic shape parameter
% C            TH(i)  momentum thickness
% C            UE(i)  edge velocity
% C            VE(i)  edge kinematic viscosity (in same units as UE*TH)
% C            NW     number of frequencies to be set
% C
% C    Output: W(k)   radian frequencies in same units as UE/TH
% C            A(i,k) amplitude distribution for frequency W(k)
% C------------------------------------------------------------------



DW = -1.50/(NW-1);

W = zeros(1,NX);
A = zeros(NW,NX);

for ii = 1:N-1
    RDL = log10( HK(ii)*TH(ii)*UE(ii)/VE(ii) );
    
    HKB = 1.0 / (HK(ii) - 1.0);
    RDLC = 2.23 + 1.35*HKB + 0.85*tanh(10.4*HKB - 7.07) - 0.1;

    if RDL >= RDLC
        ISTART = ii;
        UOT = UE(ii)/TH(ii);
        for IW = 1:NW
            WLOG = -1 + DW*(IW-1);
            W(IW) = (10^WLOG)*UOT;
        end
        for jj = ISTART+1:N
            IM = jj-1;
            UA = (UE(IM) + UE(jj))/2;
            VA = (VE(IM) + VE(jj))/2;
            TA = (TH(IM) + TH(jj))/2;
            HA = (HK(IM) + HK(jj))/2;
            if jj == N
                UA = 1.5*UE(IM) - 0.5*UE(IM-1);
                VA = 1.5*VE(IM) - 0.5*VE(IM-1);
                TA = 1.5*TH(IM) - 0.5*TH(IM-1);
                HA = 1.5*HK(IM) - 0.5*HK(IM-1);
            end
            RSP = UA*TA/VA;
            HSP = HA;
            HSP = min(HSP , 19.999 );
        
        for IW = 1:NW
            WSP = W(IW)*TA/UA;
            osmap(RSP,WSP,HSP,AR,AR_R,AR_W,AR_H,ARW_R,ARW_W,ARW_H,AI,AI_R,AI_W,AI_H,AIW_R,AIW_W,AIW_H,OK)
            DX = X(ii) - X(IM);
            A(IW,ii) = A(IW,IM) - AI * DX/TA;
            A(IW,ii) = max( A(IW,ii) , 0.0 );
        end
        end
    else
        disp('RCrit not exceeded')
    end
end
