% C***********************************************************************
% C    Module:  naca.f
% C
% C    Copyright (C) 2000 Mark Drela
% C
% C    This program is free software; you can redistribute it and/or modify
% C    it under the terms of the GNU General Public License as published by
% C    the Free Software Foundation; either version 2 of the License, or
% C    (at your option) any later version.
% C
% C    This program is distributed in the hope that it will be useful,
% C    but WITHOUT ANY WARRANTY; without even the implied warranty of
% C    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% C    GNU General Public License for more details.
% C
% C    You should have received a copy of the GNU General Public License
% C    along with this program; if not, write to the Free Software
% C    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
% C***********************************************************************
function [XX,YT,YC,XB,YB,NB,NAME] = naca5(IDES,NSIDE)
XX = zeros(1,NSIDE);
YT = zeros(1,NSIDE);
YC = zeros(1,NSIDE);
XB = zeros(1,2*NSIDE);
YB = zeros(1,2*NSIDE);


AN = 1.5;

N5 =  IDES                                        / 10000;
N4 = (IDES - N5*10000                           ) / 1000;
N3 = (IDES - N5*10000 - N4*1000                 ) / 100;
N2 = (IDES - N5*10000 - N4*1000 - N3*100        ) / 10;
N1 = (IDES - N5*10000 - N4*1000 - N3*100 - N2*10);

N543 = 100*N5 + 10*N4 + N3;

if N543 == 210 
    P = 0.05;
    M = 0.0580;
    C = 361.4;
elseif N543 == 220
    P = 0.10;
    M = 0.1260;
    C = 51.64;
elseif N543 == 230
    P = 0.15;
    M = 0.2025;
    C = 15.957;
elseif N543 == 240
    P = 0.2;
    M = 0.2900;
    C = 6.643;
elseif N543 == 250
    P = 0.25;
    M = 0.391;
    C = 3.230;
else
    disp('Illegal 5-digit designation')
    disp('First three digits must be 210, 220, ... 250')
    IDES = 0;
    return
end

T = (N2*10 + N1) / 100.0;

ANP = AN + 1.0;
for ii = 1:NSIDE
    FRAC = (I-1)/(NSIDE-1);
    if ii == NSIDE
        XX(ii) = 1;
    else
        XX(ii) = 1.0 - ANP*FRAC*(1.0-FRAC)^AN - (1.0-FRAC)^ANP;
    end
    YT(ii) = ( 0.29690*sqrt(XX(ii)) - 0.12600*XX(ii) - 0.35160*XX(ii)^2 + 0.28430*XX(ii)^3 - 0.10150*XX(ii)^4) * T / 0.20;
    if XX(ii) < M
        YC(ii) = (C/6.0) * (XX(ii)^3 - 3.0*M*XX(ii)^2 + M*M*(3.0-M)*XX(ii));
    
    else
        YC(ii) = (C/6.0) * M^3 * (1.0 - XX(ii));
    end
end

IB = 0;
for ii = NSIDE:-1:1
    IB = IB+1;
    XB(IB) = XX(ii);
    YB(IB) = YC(ii) + YT(ii);
end
for ii = 2:NSIDE
    IB = IB + 1;
    XB(IB) = XX(ii);
    YB(IB) = YC(ii) - YT(ii);
end


NB = IB;
NAME = strcat('NACA',N5+1,N4+1,N3+1,N2+1,N1+1);

